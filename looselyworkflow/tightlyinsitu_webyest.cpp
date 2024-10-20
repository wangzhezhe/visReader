#include <iostream>
#include <fstream>
#include <cstdlib>

#include <thallium/serialization/stl/string.hpp>
#include <thallium/serialization/stl/vector.hpp>
#include <thallium.hpp>
#include <mpi.h>

// vtkm header
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>

// padvection filter
#include "../AssignStrategy.hpp"
#include "../LoadData.hpp"
#include "../filter.hpp"
#include "../workloadEstimation/Block.h"
#include "../workloadEstimation/SLTool.h"

#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/filter/flow/Tracer.h>

#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
// communicator marshal
#include <vtkCommunicator.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>

namespace tl = thallium;

int globalRank = 0;
int totalRanks = 0;
std::vector<std::string> globalAddrList;

vtkSmartPointer<vtkDataSet> LoadDataIntoVTK(const std::string &visitfileName)
{

    // open .visit file
    std::string buff;
    std::ifstream is;
    is.open(visitfileName);
    if (globalRank == 0)
    {
        std::cout << "Opening: " << visitfileName << std::endl;
    }
    if (!is)
    {
        std::cout << "File not found! : " << visitfileName << std::endl;
        throw "unknown file: " + visitfileName;
    }

    auto p0 = visitfileName.rfind(".visit");
    if (p0 == std::string::npos)
        throw "Only .visit files are supported.";
    auto tmp = visitfileName.substr(0, p0);
    auto p1 = tmp.rfind("/");
    auto dir = tmp.substr(0, p1);

    std::getline(is, buff);
    auto numBlocks = std::stoi(buff.substr(buff.find("!NBLOCKS ") + 9, buff.size()));
    if (globalRank == 0)
        std::cout << "numBlocks= " << numBlocks << std::endl;

    // get block id list
    if (numBlocks != totalRanks)
    {
        throw std::runtime_error("numblocks is supposed to equal to totalRanks");
        exit(0);
    }

    // load each vtk file
    vtkSmartPointer<vtkDataSet> ds;
    for (int i = 0; i < numBlocks; i++)
    {
        // get a new vtk entry with the associated block id in visit file
        std::getline(is, buff);
        if (i == globalRank)
        {
            // std::cout << "rank " << globalRank << " assign blockid " << i << std::endl;
            std::string vtkFilePath = dir + "/" + buff;
            // std::cout << "rank " << globalRank << " open file " << vtkFilePath << std::endl;
            vtkSmartPointer<vtkDataSetReader> reader =
                vtkSmartPointer<vtkDataSetReader>::New();
            reader->SetFileName(vtkFilePath.c_str());
            reader->Update();
            ds = reader->GetOutput();
            break;
        }
    }
    return ds;
}

std::string loadMasterAddr(std::string masterConfigFile)
{

    std::ifstream infile(masterConfigFile);
    std::string content = "";
    std::getline(infile, content);
    // spdlog::debug("load master server conf {}, content -{}-", masterConfigFile,content);
    if (content.compare("") == 0)
    {
        std::getline(infile, content);
        if (content.compare("") == 0)
        {
            throw std::runtime_error("failed to load the master server\n");
        }
    }
    return content;
}

bool dirExists(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
        return false; // Cannot access path
    }
    else if (info.st_mode & S_IFDIR)
    {
        return true; // Path is a directory
    }
    return false; // Path exists but is not a directory
}

vtkm::cont::DataSet ConvertToVtkmDS(vtkStructuredPoints *outDataSet, std::string filedName)
{
    int dimensions[3];
    outDataSet->GetDimensions(dimensions);
    double spacijk[3];
    outDataSet->GetSpacing(spacijk);
    double origin[3];
    outDataSet->GetOrigin(origin);

    int arrayNum = outDataSet->GetPointData()->GetNumberOfArrays();
    vtkDataArray *dataArray = outDataSet->GetPointData()->GetArray(filedName.c_str());

    vtkIdType numTuples = dataArray->GetNumberOfTuples();

    vtkm::cont::DataSetBuilderUniform dsb;
    vtkm::Id3 vtkmDimensions(dimensions[0], dimensions[1], dimensions[2]);
    vtkm::Vec<double, 3> vtkmSpacing(spacijk[0], spacijk[1], spacijk[2]);
    vtkm::Vec<double, 3> vtkmOrigin(origin[0], origin[1], origin[2]);

    vtkm::cont::DataSet ds = dsb.Create(vtkmDimensions, vtkmOrigin, vtkmSpacing);

    // adding vector field into vtkm data set from vtkDoubleArray
    vtkm::cont::ArrayHandle<vtkm::Vec3f> vecField;
    vecField.Allocate(numTuples);
    auto wPortal = vecField.WritePortal();

    std::array<double, 3> tuple;
    for (vtkIdType t = 0; t < numTuples; ++t)
    {
        dataArray->GetTuple(t, tuple.data());
        wPortal.Set(vtkm::Id(t), vtkm::Vec<double, 3>(tuple[0], tuple[1], tuple[2]));
    }
    //it is still mesh_mesh/velocity here
    ds.AddPointField(filedName.c_str(), vecField);
    return ds;
}

std::vector<int> GetServerIdListByEst(int c, int numServers, std::string assignFileName, vtkSmartPointer<vtkDataSet> inData)
{
    std::vector<int> serverIDList;
    // For advection case, we do not need care about the first round
    // there is no cold start issue

    // converting vtk data sets to vtkm data sets for each rank
    std::vector<vtkm::cont::DataSet> vtkmDataSets;
    vtkStructuredPoints *outDataSet = vtkStructuredPoints::SafeDownCast(inData);

    std::string fieldNm = "mesh_mesh/velocity";
    vtkm::cont::DataSet vtkmInDS = ConvertToVtkmDS(outDataSet, fieldNm);
    vtkmDataSets.push_back(vtkmInDS);

    // make sure all data sets are ready before moving to the next step
    MPI_Barrier(MPI_COMM_WORLD);

    vtkm::filter::flow::internal::BoundsMap boundsMap(vtkmDataSets);
    std::vector<DomainBlock *> blockInfo;
    /*
    NUM_TEST_POINTS=50
    NXYZ=2
    WIDTH_PCT=0.1
    */
    // the advection filter call the tracer anyway
    vtkm::filter::flow::GetTracer().Get()->Init(globalRank);
    vtkm::filter::flow::GetTracer().Get()->StartTimer();
    // reset cycle before running the filter
    vtkm::filter::flow::GetTracer().Get()->ResetIterationStep(c);
    
    //Be carefule, the effects of the estimation is highly related to the initial parameters
    //and associated data sets
    vtkm::Id Nxyz = 4;
    int NX = Nxyz, NY = Nxyz, NZ = Nxyz;
    bool subdivUniform = false;
    int totNumBlocks = totalRanks;
    double pctWidth = 0.1;
    double stepSize = 0.005;
    int maxSteps = 2000;
    int numFacePts = 20;
    int numTestPts = 1000;

    // if(globalRank==0){
    //     std::cout << "---Debug CreateBlockInfo " << totNumBlocks << " " << subdivUniform << " "
    //      << NX << " " << NY << " " << NZ << " " << pctWidth << std::endl;
    // }
    DomainBlock::CreateBlockInfo(blockInfo, totNumBlocks, boundsMap, subdivUniform, NX, NY, NZ, pctWidth);

    std::vector<int> blockPopularity(totNumBlocks, 0), particlesIn(totNumBlocks, 0), particlesOut(totNumBlocks, 0), cycleCnt(totNumBlocks, 0);

    std::vector<int> blockAdvMultiStages(totNumBlocks * 3, 0);

    int blockID = boundsMap.GetLocalBlockId(0);
    const auto &ds = vtkmDataSets[0];
    auto block = blockInfo[blockID];


    // std::cout << "debug CreateFlowMap rank " << globalRank << " Nxyz " << Nxyz << " blockID " << blockID
    //           << " fieldNm " << fieldNm << " stepSize " << stepSize << " maxSteps " << maxSteps << " numFacePts " << numFacePts << std::endl;
    auto flowMap = CreateFlowMap(Nxyz, blockInfo, ds, blockID, fieldNm, stepSize, maxSteps, numFacePts, boundsMap);

    CalcBlockPopularityMultiStages(blockInfo, ds, blockID, flowMap, boundsMap, blockAdvMultiStages, particlesIn, particlesOut, cycleCnt, numTestPts, fieldNm, stepSize, maxSteps);

    if (globalRank == 0)
    {
        std::cout << std::endl;
        int sum = std::accumulate(blockAdvMultiStages.begin(), blockAdvMultiStages.end(), 0);
        std::vector<float> blockPopNorm;
        for (const auto &v : blockAdvMultiStages)
        {
            // std::cout << "debug blockAdvMultiStages " << v << std::endl;
            blockPopNorm.push_back((float)(v) / sum);
        }

        // output blockPopNorm into file
        std::ofstream logfile;
        logfile.open("adv_step_stages_list.json");
        logfile << "[";
        for (int i = 0; i < blockPopNorm.size() / 3; i++)
        {
            int startPos = i * 3;
            if (i > 0)
            {
                logfile << ",";
            }
            logfile << "[" << blockPopNorm[startPos] << "," << blockPopNorm[startPos + 1] << "," << blockPopNorm[startPos + 2] << "]";
        }
        logfile << "]" << std::endl;
        logfile.close();

        // ok for write out the workload estimation results
        // create assignment plan based on python call (todo, move to cpp in future)

        // create assignment plan
        // parse the file
        // for tightly in situ case, #block equals to #ranks
        int totalBlocks = totalRanks;

        // parser the log of the execution results
        // make sure the output name of log file is tightlyinsitu_webyest_est.log
        // std::string command1 = "python3 parser_estimation_run.py tightlyinsitu_webyest_execution.log " + std::to_string(totalBlocks) + " 3 ./adv_step_stages_list.json";
        // std::cout << "execute command1: " << command1 << std::endl;
        // int result = system(command1.c_str());
        // if (result != 0)
        //{
        //    throw std::runtime_error("Failed to execute command1 and generate assignment plan");
        //}

        std::string command2 = "python3 generate_assignment_actual_bpacking_dup_capacity_vector.py " + std::to_string(totalBlocks) + " " + std::to_string(numServers) + " ./adv_step_stages_list.json";
        std::cout << "execute command2: " << command2 << std::endl;
        int result = system(command2.c_str());
        if (result != 0)
        {
            throw std::runtime_error("Failed to execute command2 and generate assignment plan");
        }
    }

    // For each rank, pasrse the block assignment file
    MPI_Barrier(MPI_COMM_WORLD);
    // load results from assignFileName to know
    // which client need to send data to which server
    std::string line;
    //make sure each rank knows the assignFileName
    assignFileName = "assign_options.config";
    std::ifstream infile(assignFileName.c_str());
    int serverIDTemp = 0;
    // assuming rankid is block id
    std::string global_rank_str = std::to_string(globalRank);
    while (std::getline(infile, line))
    {
        // if (globalRank == 0)
        // {
        //     std::cout << "get line " << line << std::endl;
        // }

        // split the line based on space
        // put results into a vector
        std::vector<int> blockids;
        int leftIdx = 0;
        int rightIdx = 0;
        while (rightIdx < line.size())
        {
            while (rightIdx < line.size() && line[rightIdx] != ' ')
                rightIdx++;
            blockids.push_back(std::stoi(line.substr(leftIdx, rightIdx - leftIdx)));
            leftIdx = rightIdx + 1;
            rightIdx = rightIdx + 1;
        }

        // for (int k = 0; k < blockids.size(); k++)
        // {
        //     if (globalRank == 0)
        //     {
        //         std::cout << "debug block id " << blockids[k] << std::endl;
        //     }
        // }
        for (int k = 0; k < blockids.size(); k++)
        {
            if (blockids[k] == globalRank)
            {
                // find element
                // std::cout << "debug serverIDList, globalRank " << globalRank << " push server id " << serverIDTemp << std::endl;
                serverIDList.push_back(serverIDTemp);
                break;
            }
        }
        // move to the next line
        serverIDTemp++;
    }

    return serverIDList;
}

// The tightly part should be like a client and send data into the server(loosely part)
int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);

    // declared in SLTool.h
    SLTOOL_Rank = globalRank;
    PRINT_RANK = -1;
    PRINT_DETAILS = -1;
    SLTOOL_TOTALPROCS = totalRanks;

    // set necessary vtkm arguments and timer information
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    // vtkm::cont::SetStderrLogLevel(vtkm::cont::LogLevel::Perf);
    vtkm::cont::Timer timer{initResult.Device};
    timer.Start();

    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0] << " <protocol> <masterconf> <dataset_prefix> <dataset_suffix> <#cycle> <simTime>" << std::endl;
        exit(0);
    }

    // let rank 0 to detect the address of vis server
    std::string protocol = argv[1];
    std::string masterconf = argv[2];
    std::string dataset_dir_prefix = argv[3];
    std::string dataset_dir_suffix = argv[4];
    int numCycles = std::stoi(argv[5]);
    int simTime = std::stoi(argv[6]);

    // typical suffix can be ".2_2_2.128_128_128.visit"

    std::string masterAddr = loadMasterAddr(masterconf);

    std::cout << "detect server addr: " << masterAddr << std::endl;

    // broad cast store all server in list
    tl::engine myEngine(protocol, THALLIUM_CLIENT_MODE);
    // tl::remote_procedure sum   = myEngine.define("sum");
    // tl::remote_procedure prod  = myEngine.define("prod");
    tl::remote_procedure hello = myEngine.define("hello").disable_response();
    tl::remote_procedure print = myEngine.define("print").disable_response();
    tl::remote_procedure getaddrlist = myEngine.define("getaddrlist");
    tl::remote_procedure stagetest = myEngine.define("stagetest");
    tl::remote_procedure runfilter = myEngine.define("runfilter");
    tl::remote_procedure stage = myEngine.define("stage");
    tl::remote_procedure finalize = myEngine.define("finalize");

    // loosely coupled in situ mode
    // load the master server addr

    tl::endpoint masterep = myEngine.lookup(masterAddr);
    // just fix it
    uint16_t provider_id = 22;

    // create provider handler
    tl::provider_handle ph(masterep, provider_id);
    // std::cout << "identity: " << ph.get_identity() << std::endl;

    if (globalRank == 0)
    {
        std::string testStr = "tightly in situ client";
        hello.on(ph)(testStr);
        std::cout << "Done sending hello RPC, no response expected" << std::endl;
    }

    // get addrlist
    std::vector<std::string> globalAddrList = getaddrlist.on(ph)();

    if (globalRank == 0)
    {
        std::cout << "addr list len:" << globalAddrList.size() << std::endl;
        for (auto addr : globalAddrList)
        {
            std::cout << "get server addr: " << addr << std::endl;
        }
    }

    std::vector<tl::async_response> reqlist;
    if (globalRank == 0)
    {
        double initTime = timer.GetElapsedTime();
        std::cout << "Time init ok is: " << initTime << std::endl;
    }

    for (int c = 0; c < numCycles; c++)
    {
        // load the data
        // make sure #data blocks equals to #ranks
        if (globalRank == 0)
        {
            std::cout << "---start to process cycle " << c << " ---" << std::endl;
        }
        std::vector<vtkm::cont::DataSet> vtkmDataSets;
        std::vector<int> blockIDList;
        std::string visitfileName = dataset_dir_prefix + std::to_string(c) + dataset_dir_suffix;
        std::string assignFileName = ""; // do not use it for round roubin

        // load the data to vtk file
        vtkSmartPointer<vtkDataSet> inData = LoadDataIntoVTK(visitfileName);
        // marshal the vtk data (data transfer)
        vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
        bool oktoMarshal = vtkCommunicator::MarshalDataObject(inData.GetPointer(), buffer);
        //buffer->Print(std::cout);
        if (oktoMarshal == false)
        {
            throw std::runtime_error("failed to marshal vtk data");
        }

        if (globalRank == 0)
        {
            std::cout << "ok to load the data, start sim" << std::endl;
        }

        // staging request is ok, start to sleep
        // if there is existing staging, we are overlapping the
        // staging and the simulation here
        sleep(simTime);

        if (globalRank == 0)
        {
            double loadAndSimOk = timer.GetElapsedTime();
            std::cout << "Time load and sim ok for cycle " << c << " is: " << loadAndSimOk << std::endl;
        }

        // we can do some data processing opeartions when waiting the server
        // such as processing the log etc.
        //  consider sending the vtk data to the remote server
        int numServers = globalAddrList.size();
        // This is using the previous log to find server id
        // use new server id when this part is ready
        std::vector<int> serverIDList = GetServerIdListByEst(c, numServers, assignFileName, inData);
        MPI_Barrier(MPI_COMM_WORLD);
        if (globalRank == 0)
        {
            double processLog = timer.GetElapsedTime();
            std::cout << "Time processing log ok is: " << processLog << std::endl;

            // make sure all response for runfilter is ok, this is last steps thing for exeecuting the run filter command

            for (int i = 0; i < reqlist.size(); i++)
            {
                // std::cout << "global rank " << globalRank << " req list " << reqlist.size() << std::endl;
                int status = reqlist[i].wait();
                if (status != 0)
                {
                    std::cout << "failed to call runfilter for some request " << i << std::endl;
                }
            }
        }

        // All processes need to wait for the completion of last step's run
        // we can start to send the data to server
        // after this point, the vis opertaion is ok
        // and the data buffer at server end can be reused.
        MPI_Barrier(MPI_COMM_WORLD);
        if (globalRank == 0)
        {
            if (c > 0)
            {
                double runfilterok = timer.GetElapsedTime();
                std::cout << "Time runfilter ok is: " << runfilterok << std::endl;
            }
        }

        // start to send data to corresponding servers
        // only give to first one
        std::vector<tl::async_response> stgReqlist;
        // we may send data to multiple servers
        for (int s = 0; s < serverIDList.size(); s++)
        {
            //std::cout << "rank " << globalRank << "  send data to server with id " << serverIDList[s] << std::endl;
            tl::endpoint serverEndpoint = myEngine.lookup(globalAddrList[serverIDList[s]]);
            tl::provider_handle stage_ph(serverEndpoint, provider_id);

            std::vector<std::pair<void *, std::size_t>> stgsegments(1);
            stgsegments[0].first = buffer->GetPointer(0);
            std::size_t sizeofstgArray = buffer->GetNumberOfTuples() * sizeof(char);
            stgsegments[0].second = sizeofstgArray;

            tl::bulk stgBulk = myEngine.expose(stgsegments, tl::bulk_mode::read_only);
            // using list here, it might send to multiple servers
            auto stgResponse = stage.on(stage_ph).async(sizeofstgArray, stgBulk, globalRank);
            // wait until the server pull the data from client
            // we can then do the next step
            // so we can not use the async now, we need to wait until pull completes
            stgResponse.wait();
        }

        // Data staging is ok for this point
        // we start to call each server to start to run the filter
        MPI_Barrier(MPI_COMM_WORLD);
        if (globalRank == 0)
        {
            double stageOk = timer.GetElapsedTime();
            std::cout << "Time stage ok is: " << stageOk << std::endl;
        }

        if (globalRank == 0)
        {
            reqlist.clear();
            for (auto addr : globalAddrList)
            {
                std::cout << "sent run filter api to " << addr << std::endl;
                tl::endpoint addrEndPoint = myEngine.lookup(addr);
                tl::provider_handle phVisServer(addrEndPoint, provider_id);
                // use async call
                std::string field = "velocity";
                auto response = runfilter.on(phVisServer).async(c, field);
                reqlist.push_back(std::move(response));
            }
        }
    }

    // when all steps complete
    std::cout << "client " << globalRank << " close" << std::endl;

    // finalize all servers
    if (globalRank == 0)
    {

        // for last round, we still need another workload estimation process
        // TODO

        for (int i = 0; i < reqlist.size(); i++)
        {
            // std::cout << "global rank " << globalRank << " req list " << reqlist.size() << std::endl;
            int status = reqlist[i].wait();
            if (status != 0)
            {
                std::cout << "failed to call runfilter for some request " << i << std::endl;
            }
        }

        double runfilterok = timer.GetElapsedTime();
        std::cout << "Time runfilter ok is: " << runfilterok << std::endl;
        for (auto addr : globalAddrList)
        {
            std::cout << "sent finalize api to " << addr << std::endl;
            tl::endpoint addrEndPoint = myEngine.lookup(addr);
            tl::provider_handle phVisServer(addrEndPoint, provider_id);
            finalize.on(phVisServer)();
        }
    }

    timer.Stop();
    if (globalRank == 0)
    {
        std::cout << "Time when workflow exit is " << timer.GetElapsedTime() << " seconds" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
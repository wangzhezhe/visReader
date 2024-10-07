#include <iostream>
#include <fstream>
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

// The tightly part should be like a client and send data into the server(loosely part)
int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);

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
        AssignStrategy assignStrategy = AssignStrategy::ROUNDROUBIN;
        std::string assignFileName = ""; // do not use it for round roubin
        // load the data to vtk file
        vtkSmartPointer<vtkDataSet> inData = LoadDataIntoVTK(visitfileName);
        // there is only one block per rank
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
            std::cout << "ok for sim, send stage api" << std::endl;
        }
        // make sure all response is ok
        for (int i = 0; i < reqlist.size(); i++)
        {
            // std::cout << "global rank " << globalRank << " req list " << reqlist.size() << std::endl;
            int status = reqlist[i].wait();
            if (status != 0)
            {
                std::cout << "failed to call runfilter for some request " << i << std::endl;
            }
        }
        // after this point, the vis opertaion is ok
        // and the data buffer at server end can be reused.

        // marshal the vtk data (data transfer)
        vtkSmartPointer<vtkCharArray> buffer = vtkSmartPointer<vtkCharArray>::New();
        bool oktoMarshal = vtkCommunicator::MarshalDataObject(inData.GetPointer(), buffer);

        if (oktoMarshal == false)
        {
            throw std::runtime_error("failed to marshal vtk data");
        }

        // buffer->Print(std::cout);
        //  consider sending the vtk data to the remote server
        int numServers = globalAddrList.size();
        int serverID = globalRank % numServers;

        // std::cout << "rank " << globalRank << " send data to server with id " << serverID << std::endl;

        // send the marshaled object to remote server
        // using bulk transfer

        tl::endpoint serverEndpoint = myEngine.lookup(globalAddrList[serverID]);
        tl::provider_handle stage_ph(serverEndpoint, provider_id);

        std::vector<std::pair<void *, std::size_t>> stgsegments(1);
        stgsegments[0].first = buffer->GetPointer(0);
        std::size_t sizeofstgArray = buffer->GetNumberOfTuples() * sizeof(char);
        stgsegments[0].second = sizeofstgArray;

        tl::bulk stgBulk = myEngine.expose(stgsegments, tl::bulk_mode::read_only);
        auto stgResponse = stage.on(stage_ph).async(sizeofstgArray, stgBulk);
        int status = stgResponse.wait();
        if (status != 0)
        {
            std::cout << "failed to stage the data for rank " << globalRank << std::endl;
        }

        // when sending ok for all ranks

        // TODO call run filter for testing
        // maybe we should call the master, and master will call run for all associated ranks
        MPI_Barrier(MPI_COMM_WORLD);
        if (globalRank == 0)
        {
            // // for controller rank
            // // send run api for all staging service
            // for (auto addr : globalAddrList)
            // {
            //     std::cout << "sent stagetest api to " << addr << std::endl;
            //     tl::endpoint addrEndPoint = myEngine.lookup(addr);
            //     tl::provider_handle phVisServer(addrEndPoint, provider_id);
            //     // use async call
            //     auto response = stagetest.on(phVisServer).async();
            //     reqlist.push_back(std::move(response));
            // }

            // // make sure all response is ok
            // for (int i = 0; i < reqlist.size(); i++)
            // {
            //     int status = reqlist[i].wait();
            //     if (status != 0)
            //     {
            //         std::cout << "failed to stagetest the data for some request " << i << std::endl;
            //     }
            //     // TODO, hangs at the reqlist here
            // }
            // std::cout << "Controller ok to call the stagetest" << std::endl;

            // clean the reqlist
            // call the runfilter on remote server
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

    std::cout << "client " << globalRank << " close" << std::endl;
    // TODO, notify serer to close the connection
    timer.Stop();
    if (globalRank == 0)
    {
        std::cout << "whole workflow exec time " << timer.GetElapsedTime() << " seconds" << std::endl;
    }

    // finalize all servers
    if (globalRank == 0)
    {
        for (auto addr : globalAddrList)
        {
            std::cout << "sent finalize api to " << addr << std::endl;
            tl::endpoint addrEndPoint = myEngine.lookup(addr);
            tl::provider_handle phVisServer(addrEndPoint, provider_id);
            finalize.on(phVisServer)();
        }
    }

    MPI_Finalize();
    return 0;
}
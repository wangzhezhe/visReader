#include <iostream>
#include <thallium.hpp>
#include <thallium/serialization/stl/string.hpp>
#include <thallium/serialization/stl/vector.hpp>
#include <spdlog/spdlog.h>
#include <mpi.h>
#include <fstream>

// vtkm header
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

// padvection filter
#include "../AssignStrategy.hpp"
#include "../LoadData.hpp"
#include "../filter.hpp"
#include <vtkm/filter/flow/Tracer.h>

// vtk data
#include <vtkCommunicator.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>

// The loosely part can be a visualization service
// it need to be started when all the tightly insitu (client) is well started
namespace tl = thallium;

int globalRank = 0;
int totalRanks = 0;

std::vector<std::string> globalAddrList;
// TODO there should be a switch betweem vtkm data sets 1 and 2
// duing filtering operation, the vtkmDataSets is occupied by filter
// and can not be occupied by another set of data staging operation.
bool occupyingOne = false;
tl::mutex vtkmDataSetsMutex;
std::vector<vtkm::cont::DataSet> vtkmDataSets;
tl::engine *globalServerEnginePtr = nullptr;

std::vector<vtkm::Id> local_blockids;

enum VisOpEnum
{
    UNDEFINEDOP,
    ISO,
    VOLUME,
    ADVECT,
    CLOVERADVECT
};

VisOpEnum myVisualizationOperation = ADVECT;
// serial is 1
vtkm::cont::DeviceAdapterId GlobalDeviceID = vtkm::cont::make_DeviceAdapterId(1);

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
    ds.AddPointField("velocity", vecField);
    return ds;
}

class my_provider : public tl::provider<my_provider>
{

private:
    // for testing
    tl::auto_remote_procedure m_prod_rpc;
    tl::auto_remote_procedure m_sum_rpc;
    tl::auto_remote_procedure m_hello_rpc;
    tl::auto_remote_procedure m_print_rpc;

    // get all addresses
    tl::auto_remote_procedure m_getaddrlist_rpc;

    // call the vtkm particle advection filter
    tl::auto_remote_procedure m_stagetest_rpc;
    tl::auto_remote_procedure m_runfilter_rpc;

    tl::auto_remote_procedure m_stage_rpc;
    tl::auto_remote_procedure m_finalize_rpc;

    void prod(const tl::request &req, int x, int y)
    {
        std::cout << "Computing " << x << "*" << y << std::endl;
        req.respond(x * y);
    }

    int sum(int x, int y) const
    {
        std::cout << "Computing " << x << "+" << y << std::endl;
        return x + y;
    }

    void hello(const std::string &name)
    {
        std::cout << "Hello, " << name << ", from " << identity() << std::endl;
    }

    int print(const std::string &word)
    {
        std::cout << "Printing " << word << std::endl;
        return word.size();
    }

    void getaddrlist(const tl::request &req)
    {
        req.respond(globalAddrList);
    }

    // TODO, udpate it to the mode that recieve data from the api call
    void stagetest(const tl::request &req)
    {
        // load the data, this will be replaced by data transfer function
        // the code here is for temporary testing
        std::vector<int> blockIDList;
        std::string visitfileName = "/home/ubuntu/dataset/astro.2_2_2.visit";
        AssignStrategy assignStrategy = AssignStrategy::ROUNDROUBIN;
        std::string assignFileName = ""; // do not use it for round roubin
        // load the data
        LoadData(visitfileName, assignStrategy, assignFileName, vtkmDataSets, blockIDList, globalRank, totalRanks);

        MPI_Barrier(MPI_COMM_WORLD);
        spdlog::debug("rank {} load number of data blocks {}", globalRank, vtkmDataSets.size());

        // ok to stage data, run particle advection algorithm
        req.respond(0);
    }

    void stage(const tl::request &req, std::size_t &stgDataSize, tl::bulk &dataBulk, int blockid)
    {
        // this will be called multiple times
        // pull the buffer
        // unmarshal back to vtk data, convert to vtkm data
        // add atomic operation for vtkmDataSets
        size_t transferSize = stgDataSize;

        auto buffer = vtkSmartPointer<vtkCharArray>::New();
        buffer->SetNumberOfComponents(1);
        buffer->SetNumberOfTuples(transferSize);

        // get the transfered data by exposing the mem space
        std::vector<std::pair<void *, std::size_t>> segments(1);
        segments[0].first = (void *)(buffer->GetPointer(0));
        segments[0].second = transferSize;
        spdlog::debug("rank {} recv size {}", globalRank, transferSize);

        // transfer data
        tl::bulk currentBulk = globalServerEnginePtr->expose(segments, tl::bulk_mode::write_only);

        tl::endpoint ep = req.get_endpoint();
        // pull the data onto the server
        dataBulk.on(ep) >> currentBulk;

        // unmarshal back to vtk data
        vtkSmartPointer<vtkDataObject> recvbj = vtkCommunicator::UnMarshalDataObject(buffer);

        // Assuming we know the dedicated results
        vtkStructuredPoints *outDataSet = vtkStructuredPoints::SafeDownCast(recvbj);

        // std::cout << "\n---Unmarshal object:" << std::endl;
        // outDataSet->PrintSelf(std::cout, vtkIndent(5));
        // converting to vtkm data set
        vtkm::cont::DataSet vtkmInDs = ConvertToVtkmDS(outDataSet, "mesh_mesh/velocity");
        //vtkmInDs.PrintSummary(std::cout);

        // using atomic operation to add the data into the vector
        vtkmDataSetsMutex.lock();
        vtkmDataSets.push_back(vtkmInDs);
        vtkmDataSetsMutex.unlock();

        if (globalRank == 0)
        {
            spdlog::debug("size of vtkmDataSets is  {}", vtkmDataSets.size());
        }

        // insert corresponding blockid
        local_blockids.push_back(blockid);
        req.respond(0);
    }

    void runfilter(const tl::request &req, int cycle, std::string fieldToOperateOn)
    {
        // run the particle advection filter
        // this is called by controller at client

        // create the vtkm data set
        auto partitionedDataSet = vtkm::cont::PartitionedDataSet(vtkmDataSets);

        vtkm::filter::flow::GetTracer().Get()->Init(globalRank);
        vtkm::filter::flow::GetTracer().Get()->StartTimer();
        // reset cycle before running the filter
        vtkm::filter::flow::GetTracer().Get()->ResetIterationStep(cycle);
        vtkm::filter::flow::GetTracer().Get()->TimeTraceToBuffer("FilterStart");

        // Only testing particle advection now
        std::string seedMethod = "domainrandom";
        bool use_basic = true;
        bool recordTrajectories = false;
        bool outputResults = false;
        bool outputSeeds = false;

        MPI_Barrier(MPI_COMM_WORLD);
        FILTER::CommStrategy == "async_probe";
        FILTER::GLOBAL_ADVECT_NUM_STEPS = 2000;
        FILTER::GLOBAL_ADVECT_STEP_SIZE = 0.001;
        FILTER::GLOBAL_ADVECT_NUM_SEEDS = 1000; // per domain
        FILTER::G_xMin = 0.01;
        FILTER::G_xMax = 0.99;
        FILTER::G_yMin = 0.01;
        FILTER::G_yMax = 0.99;
        FILTER::G_zMin = 0.01;
        FILTER::G_zMax = 0.99;
        // using manual blockids
        FILTER::LOCAL_BLOCKIDS = local_blockids;
        FILTER::GLOBAL_BLOCK_MANUALID = true;

        // for (int bid = 0; bid < local_blockids.size(); bid++)
        // {
        //     std::cout << "rank " << globalRank << " local_blockid " << local_blockids[bid] << std::endl;
        // }

        FILTER::runAdvection(partitionedDataSet, globalRank, totalRanks, cycle, seedMethod, fieldToOperateOn, true, recordTrajectories, outputResults, outputSeeds, GlobalDeviceID);
        vtkm::filter::flow::GetTracer().Get()->TimeTraceToBuffer("FilterEnd");

        // ouptut trace
        std::string traceDirName = "tracelog_cycle" + std::to_string(cycle);

        if (globalRank == 0)
        {
            vtkm::filter::flow::GetTracer().Get()->CreateDir(traceDirName);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        vtkm::filter::flow::GetTracer().Get()->OutputBuffer(globalRank, traceDirName);
        vtkm::filter::flow::GetTracer().Get()->StopTimer();

        // clear the data set after each filter running
        // otherwise, the dataset size will increase forever
        vtkmDataSets.clear();
        // clear blockids
        local_blockids.clear();
        req.respond(0);
    }

    void finalize(const tl::request &req)
    {
        std::cout << "recv finalize api call" << std::endl;
        globalServerEnginePtr->finalize();
        req.respond(0);
    }

public:
    my_provider(const tl::engine &e, uint16_t provider_id)
        : tl::provider<my_provider>(e, provider_id, "myprovider"),
          m_prod_rpc{define("prod", &my_provider::prod)},
          m_sum_rpc{define("sum", &my_provider::sum)},
          m_hello_rpc{define("hello", &my_provider::hello).disable_response()},
          m_print_rpc{define("print", &my_provider::print, tl::ignore_return_value())},
          m_getaddrlist_rpc{define("getaddrlist", &my_provider::getaddrlist)},
          m_stagetest_rpc{define("stagetest", &my_provider::stagetest)},
          m_runfilter_rpc{define("runfilter", &my_provider::runfilter)},
          m_stage_rpc{define("stage", &my_provider::stage)},
          m_finalize_rpc{define("finalize", &my_provider::finalize)}
    {
    }
};

std::vector<std::string> AddrSplit(const char *s, int size, char seperatorH, char seperatorE)
{
    std::vector<std::string> result;
    typedef std::string::size_type string_size;
    string_size i = 0;
    int flag = 0;

    while (i != size)
    {
        // if flag =0 , estra str
        // if flag =1 , real content
        // int flag = 0;
        while (i != size && flag == 0)
        {
            // caculate start position

            if (s[i] == seperatorH)
            {
                flag = 1;
                ++i;
                break;
            }
            else
            {
                ++i;
            }
        }

        // caculate end position
        string_size j = i;

        while (j != size && flag == 1)
        {

            if (s[j] == seperatorE)
            {
                flag = 0;
                break;
            }

            if (flag == 1)
                ++j;
        }

        if (i != j)
        {
            char substr[100];
            std::memcpy(substr, s + i, j - i);
            substr[j - i] = '\0';
            result.push_back(std::string(substr));
            i = j;
        }
    }
    return result;
}

void gatherAddr(std::string endpoint)
{
    int msgPaddingLen = 200;
    if (endpoint.size() > msgPaddingLen)
    {
        throw std::runtime_error("current addr is longer than msgPaddingLen, reset the addr buffer to "
                                 "make it larger");
        return;
    }

    int sendLen = msgPaddingLen;
    int sendSize = sendLen * sizeof(char);
    char *sendipStr = (char *)malloc(sendSize);
    sprintf(sendipStr, "H%sE", endpoint.c_str());

    std::cout << "check send ip: " << std::string(sendipStr) << std::endl;
    int rcvLen = sendLen;
    char *rcvString = NULL;
    int rcvSize = msgPaddingLen * totalRanks * sizeof(char);
    rcvString = (char *)malloc(rcvSize);
    {
        if (rcvString == NULL)
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }
    }

    int error_code =
        MPI_Gather(sendipStr, sendLen, MPI_CHAR, rcvString, rcvLen, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (error_code != MPI_SUCCESS)
    {
        std::cout << "error for rank " << globalRank << " get MPI_GatherError: " << error_code
                  << std::endl;
    }

    if (globalRank == 0)
    {
        std::vector<std::string> ipList = AddrSplit(rcvString, msgPaddingLen * totalRanks, 'H', 'E');

        for (int i = 0; i < ipList.size(); i++)
        {
            // store all server addrs
            spdlog::debug("rank {} add raw data server {}", globalRank, ipList[i]);
            globalAddrList.push_back(ipList[i]);
        }
    }

    free(sendipStr);
    free(rcvString);
}

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);

    // vtkm init
    //  set necessary vtkm arguments and timer information
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    // vtkm::cont::SetStderrLogLevel(vtkm::cont::LogLevel::Perf);
    vtkm::cont::Timer timer{initResult.Device};
    GlobalDeviceID = initResult.Device;
    timer.Start();
    if (globalRank == 0)
    {
        std::cout << "vtkm device is " << initResult.Device.GetName() << std::endl;
    }
    if (argc != 3)
    {
        std::cerr << "looselyinsitu <protocol> <log level>" << std::endl;
        exit(0);
    }

    uint16_t provider_id = 22;
    std::string protocol = argv[1];
    std::string logLevel = argv[2];

    if (logLevel == "info")
    {
        spdlog::set_level(spdlog::level::info);
    }
    else
    {
        spdlog::set_level(spdlog::level::debug);
    }

    tl::engine myEngine(protocol, THALLIUM_SERVER_MODE);
    globalServerEnginePtr = &myEngine;

    std::string selfAddr = myEngine.self();
    std::cout << "Server running at address " << selfAddr
              << " with provider id " << provider_id << std::endl;

    tl::abt scope;

    if (globalRank == 0)
    {
        std::cout << "total ranks num: " << totalRanks << std::endl;
    }

    // start provider
    auto provider = new my_provider(myEngine, provider_id);
    myEngine.push_finalize_callback(provider, [provider]()
                                    { delete provider; });

    // collecting all addresses
    gatherAddr(selfAddr);

    if (globalRank == 0)
    {
        std::string confile = "masterinfo.conf";
        std::ofstream confFile;
        spdlog::debug("master info file: {}", confile);
        confFile.open(confile);

        if (!confFile.is_open())
        {
            spdlog::info("Could not open file: {}", confile);
            exit(-1);
        }
        confFile << selfAddr << "\n";
        confFile.close();
    }

    myEngine.wait_for_finalize();

    // we only need to finalize the tracer at the end of the program
    vtkm::filter::flow::GetTracer().Get()->Finalize();
    std::cout << "server close" << std::endl;
    timer.Stop();
    if (globalRank == 0)
    {
        std::cout << "server running time is:" << timer.GetElapsedTime() << std::endl;
    }
    MPI_Finalize();
    return 0;
}
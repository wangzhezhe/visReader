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

namespace tl = thallium;

int globalRank = 0;
int totalRanks = 0;
std::vector<std::string> globalAddrList;

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

    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <dataset dir> <sleep time for sim> <#cycle>" << std::endl;
        exit(0);
    }

    // let rank 0 to detect the address of vis server
    std::string dataset_dir_suffix = argv[1];
    int sim_sleep_time=std::stoi(argv[2]);
    int total_cycle=std::stoi(argv[3]);

    for (int cycle = 0; cycle < total_cycle; cycle++)
    {
        //TODO add time tracer       
        if(globalRank==0){
            std::cout << "---start to process data generated at cycle " << cycle  << " synthetic sim"<< std::endl;
        }
        //synthetic sim

        sleep(sim_sleep_time);

        // load the vtk data set for current cycle
        std::vector<vtkm::cont::DataSet> vtkmDataSets;
        std::vector<int> blockIDList;
        std::string visitfileName = dataset_dir_suffix+std::to_string(cycle)+".2_2_2.128_128_128.visit";
        AssignStrategy assignStrategy = AssignStrategy::ROUNDROUBIN;
        std::string assignFileName = ""; // do not use it for round roubin
        // load the data from file
        LoadData(visitfileName, assignStrategy, assignFileName, vtkmDataSets, blockIDList, globalRank, totalRanks);
        
        if(globalRank==0){
            std::cout << "---load number of data block " << vtkmDataSets.size() << std::endl;
            //checking details of vtkmDataSets
            vtkmDataSets[0].PrintSummary(std::cout);
        }

        // create the vtkm data set
        auto partitionedDataSet = vtkm::cont::PartitionedDataSet(vtkmDataSets);

        // call the tightly in situ workflow
        MPI_Barrier(MPI_COMM_WORLD);

        //setting parameters and run filter
        vtkm::filter::flow::GetTracer().Get()->Init(globalRank);
        vtkm::filter::flow::GetTracer().Get()->ResetIterationStep(cycle);
        vtkm::filter::flow::GetTracer().Get()->StartTimer();
        vtkm::filter::flow::GetTracer().Get()->TimeTraceToBuffer("FilterStart");

        // Only testing particle advection now
        std::string seedMethod = "domainrandom";
        bool use_basic = true;
        bool recordTrajectories = false;
        bool outputResults = false;
        bool outputSeeds = false;

        // TODO, adding parameters to support more config
        // use hard coded config now
        FILTER::CommStrategy == "async_probe";
        FILTER::GLOBAL_ADVECT_NUM_STEPS = 2000;
        std::string fieldToOperateOn = "mesh_mesh/velocity";
        FILTER::GLOBAL_ADVECT_STEP_SIZE=0.001;
        FILTER::GLOBAL_ADVECT_NUM_SEEDS=1000; //per domain
        FILTER::G_xMin = 0.01;
        FILTER::G_xMax = 0.99;
        FILTER::G_yMin = 0.01;
        FILTER::G_yMax = 0.99;
        FILTER::G_zMin = 0.01;
        FILTER::G_zMax = 0.99;
        // other parameters need to be set from outside or client
        // numsteps, stepsize, number of seeds

        FILTER::runAdvection(partitionedDataSet, globalRank, totalRanks, cycle, seedMethod, fieldToOperateOn, true, recordTrajectories, outputResults, outputSeeds, initResult.Device);
        vtkm::filter::flow::GetTracer().Get()->TimeTraceToBuffer("FilterEnd");

        // ouptut trace 
        // there are some errors for output resutls with multiple cycle
        //vtkm::filter::flow::GetTracer().Get()->OutputBuffer(globalRank);
        //vtkm::filter::flow::GetTracer().Get()->Finalize();
    }

    if (globalRank == 0)
    {
        std::cout << "client completing all tasks" << std::endl;
    }
    timer.Stop();
    if(globalRank==0){
        std::cout << "whole workflow time " << timer.GetElapsedTime() << " seconds" << std::endl;
    }
    MPI_Finalize();
    return 0;
}
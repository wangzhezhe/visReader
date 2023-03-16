#include <iostream>
#include <chrono>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include <random>
#include <string>
#include <vector>
#include <fstream>

// only load visit data generated by one step
// and then execute it by vtkh filter

#include "vtkhfilter.hpp"
#include "utilities.h"

enum VisOpEnum
{
  UNDEFINEDOP,
  ISO,
  VOLUME,
  ADVECT,
  CLOVERADVECT
};

enum AssignStrategy
{
  CONTINUOUS,
  ASSIGNFILE,
  ROUNDROUBIN
};

VisOpEnum myVisualizationOperation = ADVECT;

std::string visitfileName = "";
std::string jsonFile = "";
std::string readMethod = "";
std::string fieldToOperateOn = "velocity";
std::string seedMethod = "box";
int unknownArg = 0;
std::ofstream *timingInfo = NULL;
// use default fileName for current solution
std::string assignFileName = "assign_options.config";

// continuous round
AssignStrategy assignStrategy = AssignStrategy::CONTINUOUS;

bool cloverleaf = true;
// this is used to decide if run the streamline or only particle advection
bool recordTrajectories = false;
bool outputResults = false;

std::vector<int> assignBlocksToRank(int totalBlocks, int nRanks, int rank)
{

  std::vector<int> assignedBlocks;

  int nPer = totalBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;

  if (rank == (nRanks - 1))
    b1 = totalBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (i >= b0 && i < b1)
    {
      assignedBlocks.push_back(i);
    }
  }

  return assignedBlocks;
}

std::vector<int> assignBlocksToRankOffset(int totalBlocks, int nRanks, int rank, int offset)
{

  std::vector<int> assignedBlocks;

  int nPer = totalBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;

  // TODO, add better mechanism for offset setting
  // exp
  // only work for 2 ranks case
  if (rank == (nRanks - 1))
    b1 = totalBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (rank == 0)
    {
      if (i >= 2 && i <= 5)
      {
        assignedBlocks.push_back(i);
      }
    }
    else
    {
      if (i >= 2 && i <= 5)
      {
        continue;
      }
      else
      {
        assignedBlocks.push_back(i);
      }
    }
  }

  return assignedBlocks;
}

// round roubin / clock wise assign
std::vector<int> assignBlocksToRankRR(int totalBlocks, int nRanks, int rank)
{
  std::vector<int> assignedBlocks;

  for (int i = 0; i < totalBlocks; i++)
  {
    if (i % nRanks == rank)
    {
      assignedBlocks.push_back(i);
    }
  }

  return assignedBlocks;
}

vector<int> getIntList(std::string str)
{
  vector<int> intList;
  int leftIdx = 0;
  int rightIdx = 0;
  while (rightIdx < str.size())
  {
    while (rightIdx < str.size() && str[rightIdx] != ' ')
      rightIdx++;
    intList.push_back(std::stoi(str.substr(leftIdx, rightIdx - leftIdx)));
    leftIdx = rightIdx + 1;
    rightIdx = rightIdx + 1;
  }
  return intList;
}

// assign by file
std::vector<int> assignByFile(int totalBlocks, int nRanks, int rank)
{
  // use the default name
  std::ifstream infile(assignFileName);
  std::vector<std::vector<int>> allBlockList;
  // load the file
  std::string line;
  std::vector<int> idlist;
  // assume the input file is correct
  int blockNum = 0;
  while (std::getline(infile, line))
  {

    // split by space
    idlist = getIntList(line);
    // std::cout << idlist.size() << std::endl;
    allBlockList.push_back(idlist);
    blockNum += idlist.size();
  }

  // check the file content
  if (allBlockList.size() != nRanks)
  {
    throw std::runtime_error("allBlockList size not equals to nRanks");
  }
  if (blockNum != totalBlocks)
  {
    throw std::runtime_error("blockNum size not equals to totalBlocks");
  }

  return allBlockList[rank];
}

void LoadData(std::vector<vtkm::cont::DataSet> &dataSets, std::vector<int> &blockIDList, int rank, int nRanks)
{
  std::string buff;
  std::ifstream is;
  is.open(visitfileName);
  std::cout << "Opening: " << visitfileName << std::endl;
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
  if (rank == 0)
    std::cout << "numBlocks= " << numBlocks << std::endl;

  if (assignStrategy == AssignStrategy::CONTINUOUS)
  {
    blockIDList = assignBlocksToRank(numBlocks, nRanks, rank);
  }
  else if (assignStrategy == AssignStrategy::ROUNDROUBIN)
  {
    blockIDList = assignBlocksToRankRR(numBlocks, nRanks, rank);
  }
  else if (assignStrategy == AssignStrategy::ASSIGNFILE)
  {
    blockIDList = assignByFile(numBlocks, nRanks, rank);
  }
  else
  {
    throw std::runtime_error("unsupported assignStrategy");
  }

  for (int i = 0; i < numBlocks; i++)
  {
    std::getline(is, buff);

    auto it = std::find(blockIDList.begin(), blockIDList.end(), i);
    if (it != blockIDList.end())
    {
      std::cout << "rank " << rank << " assign blockid " << i << std::endl;
      vtkm::cont::DataSet ds;
      std::string vtkFile = dir + "/" + buff;
      vtkm::io::VTKDataSetReader reader(vtkFile);
      ds = reader.ReadDataSet();

      // wrap filed by the arrayhandle
      /*
      auto f = ds.GetField(fieldToOperateOn).GetData();
      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>> fieldArray;
      f.AsArrayHandle(fieldArray);
      int n = fieldArray.GetNumberOfValues();

      // TODO, how to set data as an approporiate type?
      auto portal = fieldArray.WritePortal();
      for (int ii = 0; ii < n; ii++)
        portal.Set(ii, vtkm::Vec<vtkm::Float32, 3>(1, 0, 0));
      */
      dataSets.push_back(ds);
    }
  }
}

void initBackend(int rank)
{
  // init the vtkh device
  char const *tmp = getenv("ASCENT_VTKM_BACKEND");
  std::string backend;
  if (tmp == nullptr)
  {
    std::cout << "no ASCENT_VTKM_BACKEND env, use openmp" << std::endl;
    backend = "openmp";
  }
  else
  {
    backend = std::string(tmp);
  }

  if (rank == 0)
  {
    std::cout << "vtkm backend is:" << backend << std::endl;
  }

  if (backend == "serial")
  {
    vtkh::ForceSerial();
  }
  else if (backend == "openmp")
  {
    vtkh::ForceOpenMP();
  }
  else if (backend == "cuda")
  {
    vtkh::ForceCUDA();
  }
  else
  {
    std::cerr << " unrecognized backend " << backend << std::endl;
  }
  return;
}

void runCoordinator(VisOpEnum myVisualizationOperation, vtkh::DataSet *data_set, int rank, int numRanks, int step)
{

  initBackend(rank);
  // Use vtkh logger for stuff
  std::stringstream ss;
  ss << "cycle_" << step + 1;
  vtkh::DataLogger::GetInstance()->OpenLogEntry(ss.str());
  vtkh::DataLogger::GetInstance()->AddLogData("cycle", step + 1);

  // select visualization operation to run
  if (myVisualizationOperation == VisOpEnum::ADVECT)
  {
    if (cloverleaf)
    {
      // vtkh::DataSet *ghosts = FILTER::runGhostStripper(data_set, rank, numRanks, step, "ascent_ghosts");
      // void runAdvection(vtkh::DataSet *data_set, int rank, int numRanks, int step, std::string seedMethod, std::string fieldToOperateOn, bool cloverleaf, bool recordTrajectories, bool outputResults);
      FILTER::runAdvection(data_set, rank, numRanks, step, seedMethod, fieldToOperateOn, cloverleaf, recordTrajectories, outputResults, false);
      // FILTER::runAdvection(ghosts, rank, numRanks, step, seedMethod, fieldToOperateOn, cloverleaf, recordTrajectories, outputResults, false);

      // runAdvection(data_set, rank, numRanks, step);
    }
    else
    {

      // rename ghost cells because VTKm currently intrinsically only looks for one specific name
      vtkm::Id numDS = data_set->GetNumberOfDomains();
      for (vtkm::Id i = 0; i < numDS; i++)
      {
        if (data_set->GetDomain(i).HasField("topo_ghosts"))
        {
          auto temp = data_set->GetDomain(i).GetField("topo_ghosts");

          if (temp.GetNumberOfValues() >= 1)
          {
            auto ghostArr = temp.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
            const vtkm::FloatDefault *buff = ghostArr.GetReadPointer();
            vtkm::cont::ArrayHandle<vtkm::UInt8> ghosts;
            ghosts.Allocate(temp.GetNumberOfValues());
            for (vtkm::Id z = 0; z < temp.GetNumberOfValues(); z++)
            {
              ghosts.WritePortal().Set(z, static_cast<vtkm::UInt8>(buff[z]));
            }
            data_set->GetDomain(i).AddCellField("vtkmGhostCells", ghosts);
            // data.GetDomain(i).RemoveField("topo_ghosts");
          }
        }
        else
        {
          throw std::logic_error("topo_ghosts does not exist in data_set");
        }
      }

      // create velocity field from x y z velocity
      for (int currentPartition = 0;
           currentPartition < data_set->GetNumberOfDomains();
           currentPartition++)
      {
        auto xarrayV = data_set->GetDomain(currentPartition).GetField("velocityx").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
        auto yarrayV = data_set->GetDomain(currentPartition).GetField("velocityy").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
        auto zarrayV = data_set->GetDomain(currentPartition).GetField("velocityz").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
        data_set->GetDomain(currentPartition).AddCellField("velocityVector", vtkm::cont::make_ArrayHandleSOA<vtkm::Vec3f_64>({xarrayV, yarrayV, zarrayV}));
      }

      vtkh::EVENT_BEGIN("recenter");
      std::chrono::steady_clock::time_point start = startTimer();
      vtkh::Recenter recenter;
      recenter.SetInput(data_set);
      recenter.SetField("velocityVector");
      recenter.SetResultAssoc(vtkm::cont::Field::Association::POINTS);
      recenter.Update();
      vtkh::DataSet *recenteredData = recenter.GetOutput();
      // recenteredData->PrintSummary(cerr);
      float totalTime = endTimer(start);
      fprintf(stderr, "\n%i, VISapp_%i_%i, recenter, %f", step, rank, numRanks, totalTime);
      (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                    << ", recenter, " << totalTime << std::endl;
      vtkh::EVENT_END("recenter");

      // runAdvection(recenteredData, rank, numRanks, step);

      FILTER::runAdvection(recenteredData, rank, numRanks, step, seedMethod, fieldToOperateOn, cloverleaf, recordTrajectories, outputResults, false);
    }
  }
  else
  {
    std::cout << "unsupported operation" << std::endl;
    return;
  }

  // Close logger for this step
  vtkh::DataLogger::GetInstance()->CloseLogEntry();
}

void runTest(int totalRanks, int myRank)
{

  std::vector<vtkm::cont::DataSet> vtkmDataSets;
  std::vector<int> blockIDList;
  // load the data
  LoadData(vtkmDataSets, blockIDList, myRank, totalRanks);

  if (vtkmDataSets.size() != blockIDList.size())
  {
    throw std::runtime_error("vtkmDataSets is supposed to equal blockIDList");
  }

  // create the vtkh data set
  vtkh::DataSet vtkhDataSets;
  for (unsigned int i = 0; i < vtkmDataSets.size(); i++)
  {
    // vtkmDataSets[i].PrintSummary(std::cout);
    // std::cout << "domain id is " << blockIDList[i] << std::endl;
    vtkhDataSets.AddDomain(vtkmDataSets[i], blockIDList[i]);
  }
  // it is always 0 since we only look at one step here
  vtkhDataSets.SetCycle(0);

  // vtkhDataSets.PrintSummary(std::cout);

  // Init this ranks timing file
  if (timingInfo == NULL)
  {
    timingInfo = new std::ofstream;
    char nm[32];
    sprintf(nm, "timing.vis.%d.out", myRank);
    timingInfo->open(nm, ofstream::out);
    FILTER::timingInfo = timingInfo;
  }

  // make sure all reader goes to same step
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < 2; i++)
  {
    runCoordinator(myVisualizationOperation, &vtkhDataSets, myRank, totalRanks, i);
  }

  return;
}

void printUsage(int argc, char **argv, int rank, int numTasks)
{
  if (rank == 0)
  {
    fprintf(stderr,
            "\n\tUSAGE: %s \n"
            "\t\tRequired Arguments:\n"
            "\t\t\t--file=<path+visit file (one step)>\n"
            "\t\t\t--field-name=<string> (default 'pressure')\n"
            "\t\t\t--advect-num-steps=<int> (default 10)\n"
            "\t\t\t--advect-num-seeds=<int> (default 100)\n"
            "\t\t\t--advect-step-size=<float> (default .1)\n"
            "\t\t\t--advect-seed-box-extents=<float> (default 0,0,0,0,0,0)\n"
            "\t\t\t--sim-code=<cloverleaf/wind> (default cloverleaf)\n"
            "\t\t\t--record-trajectories=<true/false> (default false)\n"
            "\t\t\t--output-results=<true/false> (default false)\n"
            "\t\t\t--help || --h || -h\n"
            "\n\n",
            argv[0]);
    printf("\n\t--Ran with :: Number of tasks=%d--\n\n", numTasks);
  }
  exit(13);
} // END printUsage

// method to check user input args and set up the program
void checkArgs(int argc, char **argv, int rank, int numTasks)
{
  char repeatargs[2048];
  sprintf(repeatargs, "\n\tRunning %s with:\n", argv[0]);

  char unrecognizedArgs[2048];
  sprintf(unrecognizedArgs, "\n\t!!WARNING!! Passed unrecognized argument:\n");

  for (int i = 1; i < argc; i++)
  {
    std::string longvarNm(argv[i]);
    std::string optionName = longvarNm.substr(0, longvarNm.find("=", 1, 1) + 1);
    std::string optionValue =
        longvarNm.substr(longvarNm.find("=", 1, 1) + 1, longvarNm.length());

    if (optionName == "")
    {
      optionName = longvarNm;
    }

    if (optionName == "--help" || optionName == "-h" || optionName == "--h")
    {
      printUsage(argc, argv, rank, numTasks);
    }
    else if (optionName == "--file=")
    {
      visitfileName = optionValue;
      // set args to repeat to user
      char str[1024];
      sprintf(str, "\t\t--file=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--sim-code=")
    {
      if (optionValue == "cloverleaf")
      {
        strcat(repeatargs, "\t\t--sim-code=cloverleaf\n");
        cloverleaf = true;
      }
      else if (optionValue == "wind")
      {
        strcat(repeatargs, "\t\t--sim-code=wind\n");
        cloverleaf = false;
      }
    }
    else if (optionName == "--field-name=")
    {
      char str[1024];
      sprintf(str, "\t\t--field-name=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
      fieldToOperateOn = optionValue.c_str();
    }
    else if (optionName == "--advect-num-steps=")
    {
      FILTER::GLOBAL_ADVECT_NUM_STEPS = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-num-steps=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--advect-num-seeds=")
    {
      FILTER::GLOBAL_ADVECT_NUM_SEEDS = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-num-seeds=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--advect-step-size=")
    {
      FILTER::GLOBAL_ADVECT_STEP_SIZE = atof(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-step-size=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--record-trajectories=")
    {
      if (optionValue == "true")
      {
        strcat(repeatargs, "\t\t--record-trajectories=true\n");
        recordTrajectories = true;
      }
      else if (optionValue == "false")
      {
        strcat(repeatargs, "\t\t--record-trajectories=false\n");
        recordTrajectories = false;
      }
    }
    else if (optionName == "--output-results=")
    {
      if (optionValue == "true")
      {
        strcat(repeatargs, "\t\t--output-results=true\n");
        outputResults = true;
      }
      else if (optionValue == "false")
      {
        strcat(repeatargs, "\t\t--output-results=false\n");
        outputResults = false;
      }
    }
    else if (optionName == "--seeding-method=")
    {
      if (optionValue == "box")
      {
        strcat(repeatargs, "\t\t--seeding-method=box\n");
        seedMethod = "box";
      }
      else if (optionValue == "boxsample")
      {
        strcat(repeatargs, "\t\t--seeding-method=boxsample\n");
        seedMethod = "boxsample";
      }
    }
    else if (optionName == "--seeding-sample=")
    {
      // this will be useful when the seeding method is the boxsample
      int position = 0;
      char *token;
      char *rest = new char[optionValue.length() + 1];
      strcpy(rest, optionValue.c_str());
      int sampleVlaue[3];
      while ((token = strtok_r(rest, ",", &rest)))
      {
        sampleVlaue[position] = atof(token);
        // printf("%s\n", token);
        position++;
      }
      FILTER::G_SampleX = sampleVlaue[0];
      FILTER::G_SampleY = sampleVlaue[1];
      FILTER::G_SampleZ = sampleVlaue[2];

      char str[1024];
      sprintf(str, "\t\t--seeding-sample==%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--advect-seed-box-extents=")
    {
      int position = 0;
      char *token;
      char *rest = new char[optionValue.length() + 1];
      strcpy(rest, optionValue.c_str());
      vtkm::FloatDefault minMax[6];
      while ((token = strtok_r(rest, ",", &rest)))
      {
        minMax[position] = atof(token);
        // printf("%s\n", token);
        position++;
      }
      FILTER::G_xMin = minMax[0];
      FILTER::G_xMax = minMax[1];
      FILTER::G_yMin = minMax[2];
      FILTER::G_yMax = minMax[3];
      FILTER::G_zMin = minMax[4];
      FILTER::G_zMax = minMax[5];

      // std::cout << "debug extends " << FILTER::G_xMin << " " << FILTER::G_xMax << " "
      //<< FILTER::G_yMin << " " << FILTER::G_yMax << " " << FILTER::G_zMin << " "
      //<< FILTER::G_zMax << std::endl;

      char str[1024];
      sprintf(str, "\t\t--advect-seed-box-extents=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--assign-strategy=")
    {
      if (optionValue == "continuous")
      {
        strcat(repeatargs, "\t\t--assign-strategy=continuous\n");
        assignStrategy = AssignStrategy::CONTINUOUS;
      }
      else if (optionValue == "roundroubin")
      {
        strcat(repeatargs, "\t\t--assign-strategy=roundroubin\n");
        assignStrategy = AssignStrategy::ROUNDROUBIN;
      }
      else if (optionValue == "file")
      {
        strcat(repeatargs, "\t\t--assign-strategy=file\n");
        assignStrategy = AssignStrategy::ASSIGNFILE;
      }
      else
      {
        throw std::runtime_error("--assign-strategy=continuous/roundroubin/file");
      }
    }
    else
    {
      unknownArg = 1;
      char str[1024];
      sprintf(str, "\t\t%s\n", longvarNm.c_str());
      strcat(unrecognizedArgs, str);
    }
  }

  // test for required args
  if (visitfileName == "")
  {
    if (rank == 0)
    {
      if (unknownArg == 1)
        printf("%s\n", unrecognizedArgs);
      printf("\n\n\t-*-*-ERROR-*-*- \t%s\n", repeatargs);
      printUsage(argc, argv, rank, numTasks);
    }
    exit(13);
  }
  else
  {
    if (rank == 0)
    {
      if (unknownArg == 1)
        printf("%s\n", unrecognizedArgs);
      printf("%s\n", repeatargs);
    }
  }
} // END checkAndSetProgramArgs

int main(int argc, char **argv)
{
  int rc, mpiRank, numTasks, len;
  char my_hostname[MPI_MAX_PROCESSOR_NAME];

  //----Init MPI args
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
  {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Get_processor_name(my_hostname, &len);
  //--

  //----create timers/counters
  vtkh::ADD_TIMER("total");
  vtkh::ADD_TIMER("ts_time");
  vtkh::ADD_EVENT("main");
  vtkh::ADD_EVENT("render");
  vtkh::ADD_EVENT("read");
  vtkh::ADD_EVENT("read_metadata");
  vtkh::ADD_EVENT("createRendering");
  vtkh::ADD_EVENT("contour");
  vtkh::ADD_EVENT("recenter");
  vtkh::ADD_EVENT("contour_setup");
  vtkh::ADD_EVENT("contour_getOutput");
  vtkh::ADD_EVENT("ghost_stripper");
  vtkh::ADD_EVENT("streamline/particleadv_total");
  vtkh::ADD_EVENT("streamline_generate_seeds");
  vtkh::ADD_EVENT("streamline_setup");
  vtkh::ADD_EVENT("streamline_update");
  vtkh::ADD_EVENT("streamline_get_output");
  vtkh::ADD_EVENT("streamline_save_advection_files");
  vtkh::ADD_EVENT("particle_advection_setup");
  vtkh::ADD_EVENT("particle_advection_update");
  vtkh::ADD_EVENT("particle_advection_get_output");
  vtkh::ADD_EVENT("particle_advection_save_advection_files");

  vtkh::SET_EVENT_T0();
  //--

  //----start total run time clock
  vtkh::TIMER_START("total");
  vtkh::EVENT_BEGIN("main");
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  //--

  //----Test and then set args
  checkArgs(argc, argv, mpiRank, numTasks);
  //--

  //----Enable vtkh mpi
  vtkh::SetMPICommHandle(MPI_Comm_c2f(MPI_COMM_WORLD));
  //--

  //----Print run setup info
  if (mpiRank == 0)
  {
    // #ifdef USE_READ_OMP
    //  printf("\t - Max number of openmp threads: %i\n", omp_get_max_threads());
    // #endif
    printf("\t - Number of tasks=%d My rank=%d Running on %s\n", numTasks,
           mpiRank, my_hostname);
  }
  //--
  // Create then run vis
  try
  {
    runTest(numTasks, mpiRank);
  }
  catch (std::exception &e)
  {
    std::cerr << __LINE__ << std::endl;
    std::cerr << "ERROR :: Exception from main loop" << std::endl;
    std::cout << e.what() << std::endl;
  }

  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;
  printf("\nRank - %i - Total time for program = %f min\n", mpiRank,
         std::chrono::duration<double>(diff).count() / 60);

  vtkh::TIMER_STOP("total");
  vtkh::EVENT_END("main");
  vtkh::DUMP_STATS("run.stats.txt");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
} // END main

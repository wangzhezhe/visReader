#include "utilities.h"

//--define vis op enum
enum VisOpEnum
{
  UNDEFINEDOP,
  ISO,
  VOLUME,
  ADVECT,
  CLOVERADVECT
};
VisOpEnum myVisualizationOperation = UNDEFINEDOP;
//--

string fileName = "";
string jsonFile = "";
string readMethod = "";
string fieldToOperateOn = "pressure";
string seedMethod = "";
long unsigned int totalSteps = 0;
int unknownArg = 0;
int GLOBAL_NUM_LEVELS = 1;
static ofstream *timingInfo = NULL;
int GLOBAL_ADVECT_NUM_STEPS = 10;
int GLOBAL_ADVECT_NUM_SEEDS = 100;
vtkm::FloatDefault GLOBAL_ADVECT_STEP_SIZE = 0.1;
vtkm::FloatDefault G_xMin = 0, G_xMax = 0, G_yMin = 0, G_yMax = 0, G_zMin = 0, G_zMax = 0;
bool cloverleaf = true;
bool recordTrajectories = true;
// print program usage message to user
void printUsage(int argc, char **argv, int rank, int numTasks)
{
  if (rank == 0)
  {
    fprintf(stderr,
            "\n\tUSAGE: %s \n"
            "\t\tRequired Arguments:\n"
            "\t\t\t--file=<path+file>\n"
            "\t\t\t--read-method=<SST/BP4>\n"
            "\t\t\t--visualization-op=<iso/volume/advect>\n"
            "\t\tAdditional Arguments:\n"
            "\t\t\t--sst-json-file=<string> \n"
            "\t\t\t--field-name=<string> (default 'pressure')\n"
            "\t\t\t--advect-num-steps=<int> (default 10)\n"
            "\t\t\t--advect-num-seeds=<int> (default 100)\n"
            "\t\t\t--advect-step-size=<float> (default .1)\n"
            "\t\t\t--advect-seed-box-extents=<float> (default 0,0,0,0,0,0)\n"
            "\t\t\t--iso-levels=<int> (default 1)\n"
            "\t\t\t--seed-method=<string> (cell/domain/box) \n"
            "\t\t\t--num-steps=<int> (default infinity for staging, all "
            "available for mpi)\n"
            "\t\t\t--sim-code=<cloverleaf/wind> (default cloverleaf)\n"
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
    string longvarNm(argv[i]);
    string optionName = longvarNm.substr(0, longvarNm.find("=", 1, 1) + 1);
    string optionValue =
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
      fileName = optionValue;

      // set args to repeat to user
      char str[1024];
      sprintf(str, "\t\t--file=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--read-method=")
    {
      if (optionValue == "SST" || optionValue == "sst")
      {
        strcat(repeatargs, "\t\t--read-method=SST\n");
        readMethod = "SST";
      }
      else if (optionValue == "BP4" || optionValue == "bp4")
      {
        strcat(repeatargs, "\t\t--read-method=BP4\n");
        readMethod = "BP4";
      }
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
    else if (optionName == "--seed-method=")
    {
      if (optionValue == "box")
      {
        strcat(repeatargs, "\t\t--seed-method=box\n");
        seedMethod = "box";
      }
      else if (optionValue == "cell")
      {
        strcat(repeatargs, "\t\t--seed-method=cell\n");
        seedMethod = "cell";
      }
      else if (optionValue == "domain")
      {
        strcat(repeatargs, "\t\t--seed-method=domain\n");
        seedMethod = "domain";
      }
    }
    else if (optionName == "--visualization-op=")
    {
      if (optionValue == "ISO" || optionValue == "iso")
      {
        strcat(repeatargs, "\t\t--visualization-op=iso\n");
        myVisualizationOperation = ISO;
      }
      else if (optionValue == "VOLUME" || optionValue == "volume")
      {
        strcat(repeatargs, "\t\t--visualization-op=volume\n");
        myVisualizationOperation = VOLUME;
      }
      else if (optionValue == "ADVECT" || optionValue == "advect")
      {
        strcat(repeatargs, "\t\t--visualization-op=ADVECT\n");
        myVisualizationOperation = ADVECT;
      }
      else if (optionValue == "CLOVERADVECT" ||
               optionValue == "cloverAdvect")
      {
        strcat(repeatargs, "\t\t--visualization-op=CLOVERADVECT\n");
        myVisualizationOperation = CLOVERADVECT;
      }
    }
    else if (optionName == "--sst-json-file=")
    {
      char str[1024];
      sprintf(str, "\t\t--sst-json-file=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
      jsonFile = optionValue.c_str();
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
      GLOBAL_ADVECT_NUM_STEPS = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-num-steps=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--advect-num-seeds=")
    {
      GLOBAL_ADVECT_NUM_SEEDS = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-num-seeds=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--advect-step-size=")
    {
      GLOBAL_ADVECT_STEP_SIZE = atof(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--advect-step-size=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
    }
    else if (optionName == "--num-steps=")
    {
      totalSteps = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--num-steps=%s\n", optionValue.c_str());
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
    else if (optionName == "--iso-levels=")
    {
      GLOBAL_NUM_LEVELS = atoi(optionValue.c_str());
      char str[1024];
      sprintf(str, "\t\t--iso-levels=%s\n", optionValue.c_str());
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
      G_xMin = minMax[0];
      G_xMax = minMax[1];
      G_yMin = minMax[2];
      G_yMax = minMax[3];
      G_zMin = minMax[4];
      G_zMax = minMax[5];

      char str[1024];
      sprintf(str, "\t\t--advect-seed-box-extents=%s\n", optionValue.c_str());
      strcat(repeatargs, str);
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
  if (fileName == "" || readMethod == "" ||
      myVisualizationOperation == UNDEFINEDOP)
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

// start a new timer
std::chrono::steady_clock::time_point startTimer()
{
  return std::chrono::steady_clock::now();
}

// stop an existing timer and print timing message
float endTimer(std::chrono::steady_clock::time_point start)
{
  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;
  return std::chrono::duration<double, std::milli>(diff).count();
}

vtkm::rendering::Camera createWindCamera(vtkm::Bounds bounds, float zoom)
{
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);

  vtkm::FloatDefault myElevation = -75;
  camera.Elevation(myElevation);
  //-make the camera equal to ascent
  /*vtkm::Vec<vtkm::Float32, 3> look = {5,5,5};
  vtkm::Vec<vtkm::Float32, 3> pos = {7.60472, 13.6603, 19.7721};
  camera.SetLookAt(look);
  camera.SetPosition(pos);
  camera.SetClippingRange(1.73205, 173.205);
  camera.Zoom(zoom);*/
  //
  return camera;
}

vtkm::rendering::Camera createCamera(vtkm::Bounds bounds, float zoom)
{
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);

  vtkm::FloatDefault myAzimuth = 65;
  vtkm::FloatDefault myElevation = 30;
  camera.Azimuth(myAzimuth);
  camera.Elevation(myElevation);
  //-make the camera equal to ascent
  /*vtkm::Vec<vtkm::Float32, 3> look = {5,5,5};
  vtkm::Vec<vtkm::Float32, 3> pos = {7.60472, 13.6603, 19.7721};
  camera.SetLookAt(look);
  camera.SetPosition(pos);
  camera.SetClippingRange(1.73205, 173.205);
  camera.Zoom(zoom);*/
  //
  return camera;
}

static vtkm::FloatDefault random01()
{
  return (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
}

void createLineOfSeeds(std::vector<vtkm::Particle> &seeds, vtkm::Particle startPoint,
                       vtkm::Particle endPoint, int numSeeds, int rank)
{
  vtkm::Particle v;
  v.Pos = startPoint.Pos - endPoint.Pos;
  // cerr << "Vector from points = " << v.Pos <<endl;

  vtkm::FloatDefault t, dt;
  if (numSeeds == 1)
  {
    t = 0.5;
    dt = 0.5;
  }
  else
  {
    t = 0.0;
    dt = 1.0 / (vtkm::FloatDefault)(numSeeds - 1);
  }

  for (int i = 0; i < numSeeds; i++)
  {
    vtkm::Particle p;
    p.Pos = endPoint.Pos + t * v.Pos;
    p.ID = static_cast<vtkm::Id>(i);
    seeds.push_back(p);
    t = t + dt;
  }

  if (verbose)
    printLineOhSeeds(seeds, startPoint, endPoint, rank);
}

void createBoxOfSeeds(vtkh::DataSet *data,
                      std::vector<vtkm::Particle> &seeds,
                      vtkm::FloatDefault xMin,
                      vtkm::FloatDefault xMax,
                      vtkm::FloatDefault yMin,
                      vtkm::FloatDefault yMax,
                      vtkm::FloatDefault zMin,
                      vtkm::FloatDefault zMax,
                      int numSeeds, int rank, int numRanks, int step)
{

  // Dave begin changes
  // Set seeds for BBox
  std::vector<vtkm::Vec3f> allSeeds;
  for (int i = 0; i < numSeeds; i++)
  {
    float x = xMin + (xMax - xMin) * random01();
    float y = yMin + (yMax - yMin) * random01();
    float z = zMin + (zMax - zMin) * random01();
    allSeeds.push_back({x, y, z});
  }

  // auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
  vtkm::Id numDomains = data->GetNumberOfDomains();
  std::vector<vtkm::cont::DataSet> dataSetVec;
  for (vtkm::Id i = 0; i < numDomains; i++)
    dataSetVec.push_back(data->GetDomain(i));

  vtkm::filter::particleadvection::BoundsMap boundsMap(dataSetVec);

  // select seeds that belongs to current domain that the rank owns
  for (int i = 0; i < numSeeds; i++)
  {
    auto blockIds = boundsMap.FindBlocks(allSeeds[i]);
    if (!blockIds.empty() && boundsMap.FindRank(blockIds[0]) == rank)
      seeds.push_back({allSeeds[i], i});
  }

  std::vector<int> seedCounts(numRanks, 0);
  seedCounts[rank] = seeds.size();
  MPI_Allreduce(MPI_IN_PLACE, seedCounts.data(), numRanks, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
  int totNum = 0;
  for (int i = 0; i < numRanks; i++)
  {
    totNum += seedCounts[i];
  }
  if (totNum != numSeeds)
  {
    throw std::runtime_error("totNum is supposed to equal numSeeds");
  }

  /*original way to create seeds
  vtkm::Particle diff, startPoint, endPoint;

  startPoint.Pos = {xMin, yMin, zMin};
  endPoint.Pos = {xMax, yMax, zMax};
  diff.Pos = endPoint.Pos - startPoint.Pos;

  for (int i = 0; i < numSeeds; i++)
  {
      vtkm::Particle p;
      p.Pos = {startPoint.Pos[0] + (diff.Pos[0] * random01()),
               startPoint.Pos[1] + (diff.Pos[1] * random01()),
               startPoint.Pos[2] + (diff.Pos[2] * random01())};
      p.ID = static_cast<vtkm::Id>(i);
      seeds.push_back(p);
  }

  //Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = numSeeds;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0)
  {
    for (auto& p : seeds)
      p.ID += offset;
  }
  */

  if (verbose)
    printBoxOhSeeds(seeds, rank, step);
}

void createSeedInEveryCellCloverleaf(std::vector<vtkm::Particle> &seeds,
                                     vtkh::DataSet *data_set,
                                     int rank,
                                     int numRanks,
                                     int step)
{
  using AxisArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
  using CartesianProduct = vtkm::cont::ArrayHandleCartesianProduct<AxisArrayType, AxisArrayType, AxisArrayType>;

  int particleCount = 0;
  vtkm::Id numDomains = data_set->GetNumberOfDomains();
  for (vtkm::Id i = 0; i < numDomains; i++)
  {
    auto tempDS = data_set->GetDomain(i);
    vtkm::cont::CellSetStructured<3> tempCS =
        tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
    // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
    auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<CartesianProduct>();

    vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
    vtkm::cont::Invoker invoke;
    invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

    auto buff = cellCenters.ReadPortal();

    for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
    {
      vtkm::Particle p;
      p.Pos = buff.Get(z);
      p.ID = static_cast<vtkm::Id>(particleCount);
      seeds.push_back(p);
      particleCount++;
    }
  }

  // Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = particleCount;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0)
  {
    for (auto &p : seeds)
      p.ID += offset;
  }

  if (verbose)
    printAllOhSeeds(seeds, rank, step);
}

void createSeedInEveryCell(std::vector<vtkm::Particle> &seeds,
                           vtkh::DataSet *data_set,
                           int rank,
                           int numRanks,
                           int step)
{
  using UniformCoordType = vtkm::cont::ArrayHandleUniformPointCoordinates;

  int particleCount = 0;
  vtkm::Id numDomains = data_set->GetNumberOfDomains();
  for (vtkm::Id i = 0; i < numDomains; i++)
  {
    auto tempDS = data_set->GetDomain(i);
    vtkm::cont::CellSetStructured<3> tempCS =
        tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
    // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
    auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<UniformCoordType>();

    vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
    vtkm::cont::Invoker invoke;
    invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

    auto tempGhosts = tempDS.GetField("topo_ghosts");
    auto ghostArr = tempGhosts.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
    const vtkm::FloatDefault *ghostBuff = ghostArr.GetReadPointer();
    auto buff = cellCenters.ReadPortal();

    for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
    {
      if (ghostBuff[z] == 0)
      {
        vtkm::Particle p;
        p.Pos = buff.Get(z);
        p.ID = static_cast<vtkm::Id>(particleCount);
        seeds.push_back(p);
        particleCount++;
      }
    }
  }

  // Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = particleCount;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0)
  {
    for (auto &p : seeds)
      p.ID += offset;
  }

  if (verbose)
    printAllOhSeeds(seeds, rank, step);
}

// Method to create a seed in every domain in a data set.
// The seed is randomely placed within the bounds of the domain.
void createSeedInEveryDomainCloverleaf(std::vector<vtkm::Particle> &seeds,
                                       vtkh::DataSet *data_set,
                                       int rank,
                                       int numRanks,
                                       int step)
{
  using AxisArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
  using CartesianProduct = vtkm::cont::ArrayHandleCartesianProduct<AxisArrayType, AxisArrayType, AxisArrayType>;

  int particleCount = 0;
  vtkm::Id numDomains = data_set->GetNumberOfDomains();
  for (vtkm::Id i = 0; i < numDomains; i++)
  {
    auto tempDS = data_set->GetDomain(i);
    vtkm::cont::CellSetStructured<3> tempCS =
        tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
    // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
    auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<CartesianProduct>();

    vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
    vtkm::cont::Invoker invoke;
    invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);
    auto buff = cellCenters.ReadPortal();

    for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
    {
      vtkm::Particle p;
      p.Pos = buff.Get(z);
      p.ID = static_cast<vtkm::Id>(particleCount);
      seeds.push_back(p);
      particleCount++;
      break;
    }
  }

  // Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = particleCount;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0)
  {
    for (auto &p : seeds)
      p.ID += offset;
  }

  if (verbose)
    printAllOhSeeds(seeds, rank, step);
}

// Method to create a seed in every domain in a data set.
// The seed is randomely placed within the bounds of the domain.
void createSeedInEveryDomain(std::vector<vtkm::Particle> &seeds,
                             vtkh::DataSet *data_set,
                             int rank,
                             int numRanks,
                             int step)
{
  using UniformCoordType = vtkm::cont::ArrayHandleUniformPointCoordinates;

  int particleCount = 0;
  vtkm::Id numDomains = data_set->GetNumberOfDomains();
  for (vtkm::Id i = 0; i < numDomains; i++)
  {
    auto tempDS = data_set->GetDomain(i);
    vtkm::cont::CellSetStructured<3> tempCS =
        tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
    // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
    auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<UniformCoordType>();

    vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
    vtkm::cont::Invoker invoke;
    invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

    auto tempGhosts = tempDS.GetField("topo_ghosts");
    auto ghostArr = tempGhosts.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
    const vtkm::FloatDefault *ghostBuff = ghostArr.GetReadPointer();
    auto buff = cellCenters.ReadPortal();

    for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
    {
      if (ghostBuff[z] == 0)
      {
        vtkm::Particle p;
        p.Pos = buff.Get(z);
        p.ID = static_cast<vtkm::Id>(particleCount);
        seeds.push_back(p);
        particleCount++;
        break;
      }
    }
  }

  // Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = particleCount;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0)
  {
    for (auto &p : seeds)
      p.ID += offset;
  }

  if (verbose)
    printAllOhSeeds(seeds, rank, step);
}

void runAdvection(vtkh::DataSet *data_set, int rank, int numRanks, int step)
{
  //---- Time the render
  std::chrono::steady_clock::time_point start = startTimer();
  vtkh::EVENT_BEGIN("streamline/particleadv_total");

  vtkh::EVENT_BEGIN("streamline_generate_seeds");
  std::vector<vtkm::Particle> seeds;
  // createLineOfSeeds(seeds, startPoint, endPoint, GLOBAL_ADVECT_NUM_SEEDS, rank);
  if (seedMethod == "box")
  {
    createBoxOfSeeds(data_set, seeds, G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax, GLOBAL_ADVECT_NUM_SEEDS, rank, numRanks, step);
  }
  else if (seedMethod == "cell")
  {
    if (cloverleaf)
      createSeedInEveryCellCloverleaf(seeds, data_set, rank, numRanks, step);
    else
      createSeedInEveryCell(seeds, data_set, rank, numRanks, step);
  }
  else if (seedMethod == "domain")
  {
    if (cloverleaf)
      createSeedInEveryDomainCloverleaf(seeds, data_set, rank, numRanks, step);
    else
      createSeedInEveryDomain(seeds, data_set, rank, numRanks, step);
  }
  vtkh::EVENT_END("streamline_generate_seeds");

  if (rank == 0)
  {
    std::cerr << "\nstep size " << GLOBAL_ADVECT_STEP_SIZE << std::endl;
    std::cerr << "maxSteps " << GLOBAL_ADVECT_NUM_STEPS << std::endl;
  }

  if (recordTrajectories)
  {
    vtkh::Streamline streamline;
    vtkh::EVENT_BEGIN("streamline_setup");
    streamline.SetSeeds(seeds);
    streamline.SetInput(data_set);
    streamline.SetField(fieldToOperateOn);
    streamline.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
    streamline.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
    streamline.SetIterationStep(step);
    vtkh::EVENT_END("streamline_setup");
    vtkh::DataSet *streamline_output = NULL;

    vtkh::EVENT_BEGIN("streamline_update");
    streamline.Update();
    vtkh::EVENT_END("streamline_update");

    vtkh::EVENT_BEGIN("streamline_get_output");
    streamline_output = streamline.GetOutput();
    vtkh::EVENT_END("streamline_get_output");

    if (writeStreamlines)
    {
      vtkh::EVENT_BEGIN("streamline_save_advection_files");
      writeDataSet(streamline_output, "advection_streamlinesOutput", rank, step);
      vtkh::EVENT_END("streamline_save_advection_files");
    }

    delete streamline_output;
  }
  else
  {
    vtkh::ParticleAdvection pa;
    vtkh::EVENT_BEGIN("particle_advection_setup");
    pa.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
    pa.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
    pa.SetSeeds(seeds);
    pa.SetField(fieldToOperateOn);
    pa.SetInput(data_set);
    pa.SetIterationStep(step);
    vtkh::EVENT_END("particle_advection_setup");

    vtkh::EVENT_BEGIN("particle_advection_update");
    pa.Update();
    vtkh::EVENT_END("particle_advection_update");

    vtkh::DataSet *particleAdvectOutput = NULL;

    vtkh::EVENT_BEGIN("particle_advection_get_output");
    particleAdvectOutput = pa.GetOutput();
    vtkh::EVENT_END("particle_advection_get_output");
    delete particleAdvectOutput;
  }

  vtkh::EVENT_END("streamline/particleadv_total");
  float totalTime = endTimer(start);
  if (recordTrajectories)
  {
    fprintf(stderr, "\n%i, VISapp_%i_%i, advect(streamline), %f", step, rank, numRanks,
            totalTime);
    (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                  << ", advect, " << totalTime << endl;
  }
  else
  {
    fprintf(stderr, "\n%i, VISapp_%i_%i, advect(particleadev), %f", step, rank, numRanks,
            totalTime);
    (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                  << ", advect, " << totalTime << endl;
  }

  //--
}

void runVolumeRender(vtkh::DataSet *data_set, int rank, int numRanks,
                     int step)
{
  int height = 1024;
  int width = 1024;

  //---- Time the render
  std::chrono::steady_clock::time_point start = startTimer();
  vtkm::Bounds bounds = data_set->GetGlobalBounds();

  vtkh::EVENT_BEGIN("createRendering");
  vtkh::Scene scene;
  // scene.SetRenderBatchSize(1);
  //-- use this to run multiple renders of the same dataset
  for (int j = 0; j < 1; j++)
  {
    char imgNm[128];
    if (cloverleaf)
      sprintf(imgNm, "clover_out.step-%05i_render#%02i", step, j);
    else
      sprintf(imgNm, "wind_out.step-%05i_render#%02i", step, j);
    vtkm::rendering::Camera cam;

    if (cloverleaf)
      cam = createCamera(bounds, .2);
    else
      cam = createWindCamera(bounds, .2);

    float bg_color[4] = {0.f, 0.f, 0.f, 1.f};
    vtkh::Render render =
        vtkh::MakeRender(height, width, cam, *data_set, imgNm, bg_color);
    render.DoRenderAnnotations(false);
    // vtkh::Scene scene;
    scene.AddRender(render);

    /*vtkh::VolumeRenderer volume;
    volume.SetInput(&data_set);
    volume.SetField("pressure");
    volume.SetNumberOfSamples(1000);
    scene.AddRenderer(&volume);
    scene.Render();   */
  }

  vtkh::VolumeRenderer volume;
  volume.SetInput(data_set);
  volume.SetField(fieldToOperateOn);
  volume.SetNumberOfSamples(100);
  scene.AddRenderer(&volume);
  vtkh::EVENT_END("createRendering");

  vtkh::EVENT_BEGIN("render");
  scene.Render();
  vtkh::EVENT_END("render");

  float totalTime = endTimer(start);
  fprintf(stderr, "\n%i, VISapp_%i_%i, render, %f", step, rank, numRanks,
          totalTime);
  (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                << ", render, " << totalTime << endl;
  //--
}

// ray trace an image
void runRender(vtkh::DataSet iso_output, int rank, int numRanks, int step)
{
  int height = 1024;
  int width = 1024;

  char imgNm[128];
  if (cloverleaf)
    sprintf(imgNm, "clover_out.step-%05d", step);
  else
    sprintf(imgNm, "wind_out.step-%05d", step);

  //---- Time the render
  std::chrono::steady_clock::time_point start = startTimer();
  vtkm::Bounds bounds = iso_output.GetGlobalBounds();

  vtkh::EVENT_BEGIN("createRendering");
  auto cam6 = createCamera(bounds, 0.0);

  float bg_color[4] = {0.f, 0.f, 0.f, 1.f};
  vtkh::Render render =
      vtkh::MakeRender(height, width, cam6, iso_output, imgNm, bg_color);
  render.DoRenderAnnotations(false);
  vtkh::Scene scene;
  vtkh::RayTracer tracer;
  tracer.SetInput(&iso_output);
  tracer.SetField(fieldToOperateOn);
  scene.AddRenderer(&tracer);
  scene.AddRender(render);
  vtkh::EVENT_END("createRendering");

  vtkh::EVENT_BEGIN("render");
  scene.Render();
  vtkh::EVENT_END("render");

  float totalTime = endTimer(start);
  fprintf(stderr, "\n%i, VISapp_%i_%i, render, %f", step, rank, numRanks,
          totalTime);
  (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                << ", render, " << totalTime << endl;
  //--
}

vtkh::DataSet *runContour(vtkh::DataSet *data_set, int rank, int numRanks,
                          int step)
{
  //---- Time the contour operation
  vtkh::EVENT_BEGIN("contour_setup");
  std::chrono::steady_clock::time_point start = startTimer();
  vtkh::MarchingCubes marcher;
  marcher.SetInput(data_set);
  marcher.SetField(fieldToOperateOn);

  const int num_levels = GLOBAL_NUM_LEVELS;
  marcher.SetLevels(num_levels);
  vtkh::EVENT_END("contour_setup");

  vtkh::EVENT_BEGIN("contour");
  marcher.Update();
  vtkh::EVENT_END("contour");

  vtkh::EVENT_BEGIN("contour_getOutput");
  vtkh::DataSet *iso_output = marcher.GetOutput();
  float totalTime = endTimer(start);
  fprintf(stderr, "\n%i, VISapp_%i_%i, contour, %f", step, rank, numRanks,
          totalTime);
  (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                << ", contour, " << totalTime << endl;
  //--
  vtkh::EVENT_END("contour_getOutput");
  return iso_output;
}

vtkh::DataSet *runGhostStripper(vtkh::DataSet *data_set, int rank, int numRanks,
                                int step, string ghostFieldName)
{
  //---- Time the ghosts
  vtkh::EVENT_BEGIN("ghost_stripper");
  std::chrono::steady_clock::time_point start = startTimer();
  vtkh::GhostStripper ghost;
  ghost.SetInput(data_set);
  ghost.SetField(ghostFieldName);

  ghost.SetMaxValue(0);
  ghost.SetMinValue(0);
  ghost.Update();

  vtkh::DataSet *ghosts = ghost.GetOutput();
  float totalTime = endTimer(start);
  fprintf(stderr, "\n%i, VISapp_%i_%i, ghost_stripper, %f", step, rank,
          numRanks, totalTime);
  (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                << ", ghost_stripper, " << totalTime << endl;
  vtkh::EVENT_END("ghost_stripper");
  //--
  return ghosts;
}

void runCoordinator(vtkh::DataSet *data_set, int rank, int numRanks, int step)
{
  // Use vtkh logger for stuff
  std::stringstream ss;
  ss << "cycle_" << step + 1;
  vtkh::DataLogger::GetInstance()->OpenLogEntry(ss.str());
  vtkh::DataLogger::GetInstance()->AddLogData("cycle", step + 1);

  // select visualization operation to run
  if (myVisualizationOperation == ISO)
  {
    vtkh::DataSet *iso_output = runContour(data_set, rank, numRanks, step);
    runRender(*iso_output, rank, numRanks, step);
    delete iso_output;
  }
  else if (myVisualizationOperation == VOLUME)
  {
    if (cloverleaf)
    {
      vtkh::DataSet *ghosts = runGhostStripper(data_set, rank, numRanks, step, "ascent_ghosts");
      runVolumeRender(ghosts, rank, numRanks, step);
    }
    else
    {
      vtkh::DataSet *ghosts = runGhostStripper(data_set, rank, numRanks, step, "topo_ghosts");
      runVolumeRender(ghosts, rank, numRanks, step);
    }
  }
  else if (myVisualizationOperation == ADVECT)
  {
    if (cloverleaf)
    {
      vtkh::DataSet *ghosts = runGhostStripper(data_set, rank, numRanks, step, "ascent_ghosts");
      runAdvection(ghosts, rank, numRanks, step);
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
          throw logic_error("topo_ghosts does not exist in data_set");
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
                    << ", recenter, " << totalTime << endl;
      vtkh::EVENT_END("recenter");

      runAdvection(recenteredData, rank, numRanks, step);
    }
  }

  // Close logger for this step
  vtkh::DataLogger::GetInstance()->CloseLogEntry();
}

/**
 * Do all of the reading and run orchestration
 */
void runTest(int rank, int numRanks)
{
  // Print basic vtkh configuration information
  if (rank == 0)
    std::cout << vtkh::AboutVTKH() << std::endl;

  // Declare timers
  std::chrono::steady_clock::time_point start;
  float totalTime, totalTimeInternal;

  // Init this ranks timing file
  if (timingInfo == NULL)
  {
    timingInfo = new ofstream;
    char nm[32];
    sprintf(nm, "timing.vis.%d.out", rank);
    timingInfo->open(nm, ofstream::out);
  }

  // initialze fides reader
  fides::io::DataSetReader reader(jsonFile);
  std::unordered_map<std::string, std::string> paths;
  paths["source"] = fileName;
  fides::DataSourceParams params;
  if (readMethod == "SST")
    params["engine_type"] = "SST";
  reader.SetDataSourceParameters("source", params);

  // Variables needed for read loop
  fides::metadata::Size nBlocks = 0;
  fides::metadata::Size nSteps = 0;
  bool metadataRead = false;
  long unsigned int currentStep = 0;

  // Loop until we are out of steps to read or reach a predefined maximum
  while (true)
  {
    vtkh::TIMER_START("ts_time");
    start = startTimer();

    if (rank == 0)
      cerr << "\n\nWorking on step -> " << currentStep << endl;

    auto status = reader.PrepareNextStep(paths);
    if (status == fides::StepStatus::NotReady)
    {
      continue;
    }
    else if (status == fides::StepStatus::EndOfStream)
    {
      // TODO, how we process this for the last step?
      // currently, the write is closed before the reader for the
      // last step, how to writer wait the reader read, then close?
      break;
    }

    // When using SST, a BeginStep call must happen before any attempt to read
    // data. For XGC, ReadMetaData requires reading the nphi variable in order
    // to determine the number of blocks, this means for XGC ReadMetaData must
    // be called after PrepareNextStep.
    fides::metadata::MetaData metaData;
    if (!metadataRead || cloverleaf == false)
    {
      vtkh::EVENT_BEGIN("read_metadata");
      metaData = reader.ReadMetaData(paths);
      metadataRead = true;
      nBlocks = 0;
      nSteps = 0;

      nBlocks =
          metaData.Get<fides::metadata::Size>(fides::keys::NUMBER_OF_BLOCKS());
      if (rank == 0)
        std::cout << "num blocks " << nBlocks.NumberOfItems << std::endl;

      // if running sst skip this check and set nsteps to infinity
      if (readMethod != "SST")
      {
        nSteps =
            metaData.Get<fides::metadata::Size>(fides::keys::NUMBER_OF_STEPS());
        if (rank == 0)
          std::cout << "num steps " << nSteps.NumberOfItems << std::endl;

        //--Check to see if we only need to process a limited number of steps
        if (totalSteps == 0)
        {
          totalSteps = nSteps.NumberOfItems;
        }
      }
      else
      {
        totalSteps = ULONG_MAX;
      }
      //--
      vtkh::EVENT_END("read_metadata");
    }

    //----figure out how much work to do
    fides::metadata::Vector<std::size_t> workBlocks =
        assignWorkBlocks(nBlocks.NumberOfItems, rank, numRanks);
    if (currentStep == 0 || cloverleaf == false) // check how much work I am doing
    {
      cerr << "\nRank[" << rank << "] is working on blocks ";
      for (unsigned int z = 0; z < workBlocks.Data.size(); z++)
        cerr << workBlocks.Data[z] << " ";
      cerr << endl;
    }
    //--

    //----read my blocks
    vtkh::EVENT_BEGIN("read");
    std::chrono::steady_clock::time_point startInternal = startTimer();
    metaData.Set(fides::keys::BLOCK_SELECTION(), workBlocks);
    // fides::metadata::Index idx(currentStep);
    // metaData.Set(fides::keys::STEP_SELECTION(), idx);
    auto output = reader.ReadStep(paths, metaData);
    if (output.GetNumberOfPartitions() == 0)
    {
      throw logic_error("Reader failed to generate the expected partitions");
    }
    /*
        //If we have ghost data we need to covert this to use for streamlines
         for (long unsigned int currentPartition = 0;
                 currentPartition < workBlocks.Data.size(); currentPartition++) {
            if(output.GetPartition(currentPartition).HasField("topo_ghosts"))
            {
                auto temp = output.GetPartition(currentPartition).GetField("topo_ghosts");

                if(temp.GetNumberOfValues() >= 1)
                {
                    auto ghostArr = temp.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
                    const vtkm::FloatDefault* buff = ghostArr.GetReadPointer();
                    vtkm::cont::ArrayHandle<vtkm::UInt8> ghosts;
                    ghosts.Allocate(temp.GetNumberOfValues());
                    for (vtkm::Id z = 0; z < temp.GetNumberOfValues(); z++)
                    {
                        ghosts.WritePortal().Set(z, static_cast<vtkm::UInt8>(buff[z]));
                    }
                    output.GetPartition(currentPartition).AddCellField("vtkmGhostCells", ghosts);
                }
            }
        }
    */

    totalTimeInternal = endTimer(startInternal);
    fprintf(stderr, "\n%li, VISapp_%i_%i, read, %f", currentStep, rank,
            numRanks, totalTimeInternal);
    (*timingInfo) << currentStep << ", VISapp_" << rank << "_" << numRanks
                  << ", read, " << totalTimeInternal << endl;
    vtkh::EVENT_END("read");
    //--

    //----dump data to vtk for verification
    if (debugVTK)
    {
      // output.PrintSummary(cerr);
      std::chrono::steady_clock::time_point startInternal = startTimer();
      char fileNm[128];
      for (long unsigned int currentPartition = 0;
           currentPartition < workBlocks.Data.size(); currentPartition++)
      {
        sprintf(fileNm, "fileReader_output_step%03ld_partition%03ld.vtk",
                currentStep, workBlocks.Data[currentPartition]);
        vtkm::io::VTKDataSetWriter writer(fileNm);
        writer.WriteDataSet(output.GetPartition(currentPartition));
      }
      totalTimeInternal = endTimer(startInternal);
      fprintf(stderr, "\n%li, VISapp_%i_%i, save_VTK, %f", currentStep, rank,
              numRanks, totalTimeInternal);
      (*timingInfo) << currentStep << ", VISapp_" << rank << "_" << numRanks
                    << ", save_VTK, " << totalTimeInternal << endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
    //--

    //---- Get all the work for this rank
    startInternal = startTimer();
    vtkh::DataSet data_set;
    for (unsigned int myBlocks = 0; myBlocks < workBlocks.Data.size();
         myBlocks++)
    {
      data_set.AddDomain(output.GetPartition(myBlocks),
                         workBlocks.Data[myBlocks]);
    }
    data_set.SetCycle(currentStep);
    totalTimeInternal = endTimer(startInternal);
    fprintf(stderr, "\n%li, VISapp_%i_%i, vtkm->vtkh, %f", currentStep, rank,
            numRanks, totalTimeInternal);
    (*timingInfo) << currentStep << ", VISapp_" << rank << "_" << numRanks
                  << ", vtkm->vtkh, " << totalTimeInternal << endl;
    //--

    // cerr << "data in the vtkh data set" << endl;
    // data_set.PrintSummary(cerr);

    //----do all of the real work
    runCoordinator(&data_set, rank, numRanks, currentStep);
    totalTime = endTimer(start);
    fprintf(stderr, "\n%li, VISapp_%i_%i, ts_time, %f", currentStep, rank,
            numRanks, totalTime);
    (*timingInfo) << currentStep << ", VISapp_" << rank << "_" << numRanks
                  << ", ts_time, " << totalTime << endl;
    vtkh::TIMER_STOP("ts_time");
    //--

    MPI_Barrier(MPI_COMM_WORLD);

    currentStep++;

    // Check if we are done based on how many steps the user asked for
    if (currentStep >= totalSteps)
    {
      break;
    }
  }

  if (timingInfo != NULL)
    timingInfo->close();
}

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
  vtkh::SET_EVENT_T0();
  //--

  //----start total run time clock
  vtkh::TIMER_START("total");
  vtkh::EVENT_BEGIN("main");
  std::chrono::steady_clock::time_point start = startTimer();
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
    //#ifdef USE_READ_OMP
    // printf("\t - Max number of openmp threads: %i\n", omp_get_max_threads());
    //#endif
    printf("\t - Number of tasks=%d My rank=%d Running on %s\n", numTasks,
           mpiRank, my_hostname);
  }
  //--

  // Create then run vis
  try
  {
    runTest(mpiRank, numTasks);
  }
  catch (std::exception &e)
  {
    cerr << __LINE__ << endl;
    cerr << "ERROR :: Exception from main loop" << endl;
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

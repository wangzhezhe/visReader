//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/Particle.h>
#include <vtkm/cont/AssignerPartitionedDataSet.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/EnvironmentTracker.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>

#include <vtkm/filter/flow/ParticleAdvection.h>
#include <vtkm/filter/flow/worklet/RK4Integrator.h>
#include <vtkm/filter/flow/worklet/Stepper.h>
#include <vtkm/filter/flow/worklet/ParticleAdvection.h>
#include <vtkm/filter/flow/worklet/Particles.h>

#include <vtkm/thirdparty/diy/diy.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>

#include <mpi.h>
#include "SLTool.h"
#include "../LoadData.hpp"
#include "../AssignStrategy.hpp"

#include <limits>

int Rank, Size;
// using long int to prevent the `larger than max_size()` error for some configurations
long int NUMVALS = -1;
int SEED_PING_PONG = 0;
bool UseMultiStages = true;

// Example computing streamlines.
// An example vector field is available in the vtk-m data directory: magField.vtk
// Example usage:
//   this will advect 200 particles 50 steps using a step size of 0.01
//
// Particle_Advection <path-to-data-dir>/magField.vtk vec 200 50 0.01 output.vtk
//

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
{
  out << "[";
  for (int i = 0; i < v.size(); i++)
  {
    out << v[i];
    if (i < v.size() - 1)
      out << ",";
  }
  out << "]";

  return out;
}

template <typename T, typename U>
std::ostream &operator<<(std::ostream &out, const std::map<T, U> &m)
{
  out << "{";
  for (auto it = m.begin(); it != m.end(); it++)
    out << "(" << it->first << " : " << it->second << "), ";
  out << "}";

  return out;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &s)
{
  out << "(";
  for (const auto &v : s)
    out << v << ", ";
  out << ")";

  return out;
}

std::map<int, std::vector<int>>
DetectCycles(const std::vector<int> &blockPath, int len,
             std::vector<int> &blockCycles)
{
  std::map<int, std::vector<int>> ret;

  int N = blockPath.size();

  // Create a sequence of ints of length len
  std::map<int, std::vector<int>> sequences;
  for (int i = 0; i < N - len + 1; i++)
  {
    std::vector<int> seq;
    for (int j = 0; j < len; j++)
      seq.push_back(blockPath[i + j]);
    sequences[i] = seq;
  }

  // Put them in a set to remove duplicates.
  std::set<std::vector<int>> candidates;
  for (auto it = sequences.begin(); it != sequences.end(); it++)
    candidates.insert(it->second);

  // std::cout<<"Sequences : "<<sequences<<std::endl;
  // std::cout<<"Candidates: "<<candidates<<std::endl;

  // Duplicates if the sizes are different.
  if (sequences.size() > candidates.size())
  {
    std::map<std::vector<int>, int> counter;
    for (auto it = sequences.begin(); it != sequences.end(); it++)
    {
      int n = 1;
      if (counter.find(it->second) != counter.end())
        n = counter[it->second] + 1;
      counter[it->second] = n;
    }

    // std::cout<<"Cycles: "<<counter<<std::endl;

    // update the cycles data.
    for (auto it = counter.begin(); it != counter.end(); it++)
    {
      int val = it->second;
      for (int j = 0; j < it->first.size(); j++)
        blockCycles[it->first[j]] += val;
    }
  }

  return ret;
}

void GenerateFaceSeeds2(const vtkm::cont::DataSet &ds,
                        std::vector<vtkm::Particle> &seeds)
{
  auto coords = ds.GetCoordinateSystem().GetData().template AsArrayHandle<vtkm::cont::ArrayHandleUniformPointCoordinates>();
  auto coordsPortal = coords.ReadPortal();
  auto cellSet = ds.GetCellSet().template AsCellSet<vtkm::cont::CellSetStructured<3>>();
  auto cellDims = cellSet.GetCellDimensions();

  // Cell center for every boundary face.
  // Dumb way to do this... Use a WorkletPointNeighborhood worklet....
  vtkm::Id idx = 0, pid = 0;
  for (vtkm::Id i = 0; i < cellDims[0]; i++)
  {
    for (vtkm::Id j = 0; j < cellDims[1]; j++)
    {
      for (vtkm::Id k = 0; k < cellDims[2]; k++)
      {
        if (i == 0 || i == cellDims[0] - 1 ||
            j == 0 || j == cellDims[1] - 1 ||
            k == 0 || k == cellDims[2] - 1)
        {
          // pick every 'nth' cell.
          if (pid < 100)
          {
            vtkm::Id ptIds[8];
            cellSet.GetCellPointIds(idx, ptIds);

            vtkm::Vec3f pt = coordsPortal.Get(ptIds[0]);
            for (vtkm::Id ii = 1; ii < 8; ii++)
            {
              pt += coordsPortal.Get(ptIds[ii]);
            }
            pt[0] /= 8.0;
            pt[1] /= 8.0;
            pt[2] /= 8.0;
            seeds.push_back(vtkm::Particle(pt, pid));
          }
          pid++;
        }
        idx++;
      }
    }
  }

  // Save out for validation.
  /*
  char fname[128];
  sprintf(fname, "output_%d.txt", Rank);
  std::ofstream of(fname, std::ofstream::out);

  idx = 0;
  for (const auto& p : seeds)
  {
    const auto pt = p.GetPosition();
    of<<idx<<", "<<pt[0]<<", "<<pt[1]<<", "<<pt[2]<<std::endl;
    idx++;
  }
  */
}







/*
std::map<int, std::vector<FlowStat>>
CreateFlowMap(vtkm::Id Nxyz,
              const std::vector<DomainBlock *> &blockInfo,
              const vtkm::cont::DataSet &ds,
              const vtkm::Id blockID,
              const std::string &fieldNm,
              vtkm::FloatDefault stepSize,
              vtkm::Id maxSteps,
              int numFacePts,
              const vtkm::filter::flow::internal::BoundsMap &boundsMap)
{

  vtkm::cont::Timer timer;
  timer.Start();
  // Go through all subdomains and compute the number of leaves
  std::vector<vtkm::Particle> seeds;
  // Each face is broken up into Nxyz x Nxyz pieces.
  // Althgough the total number of seeds in face is numFacePts * Nxyz * Nxyz
  // When we build the index, we divide each face into Nxyz*Nxyz partitions
  GenerateFaceSeeds1(ds, seeds, numFacePts * Nxyz * Nxyz);

  auto res = AdvectFaceSeeds(ds, fieldNm, stepSize, maxSteps, seeds);

  timer.Stop();
  double executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to AdvectFaceSeeds in CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  // create the flow map for the test seeds.
  //  {src : {dst, totNumSteps, totNumP}}
  std::map<int, std::vector<FlowEntry>> flowMap;
  int numLeafs = blockInfo[blockID]->NumLeafs();
  std::map<int, int> countFromSource;
  BuildFlowMap(res.Particles.ReadPortal(), seeds, blockID, blockInfo, boundsMap, flowMap, countFromSource);

  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to BuildFlowMap in CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  auto globalFlowMap = ComputeGlobalFlowMap(Nxyz, flowMap, blockID, blockInfo, countFromSource);

  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to ComputeGlobalFlowMap in CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  return globalFlowMap;
}
*/



void CalcBlockPopularity(std::vector<DomainBlock *> blockInfo,
                         const vtkm::cont::DataSet &ds,
                         int blockID,
                         std::map<int, std::vector<FlowStat>> &flowMap,
                         const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                         std::vector<int> &blockPopularity,
                         std::vector<int> &particlesIn,
                         std::vector<int> &particlesOut,
                         std::vector<int> &cycleCnt,
                         int numPts,
                         const std::string &fieldNm,
                         vtkm::FloatDefault stepSize,
                         int maxSteps)
{

  std::vector<vtkm::Particle> seeds;
  GenerateInteriorSeeds(ds, seeds, numPts);
  auto res = AdvectFaceSeeds(ds, fieldNm, stepSize, maxSteps, seeds);
  auto portal = res.Particles.ReadPortal();
  vtkm::Id n = portal.GetNumberOfValues();
  vtkm::FloatDefault maxStepsFloat = static_cast<vtkm::FloatDefault>(maxSteps);

  for (vtkm::Id i = 0; i < n; i++)
  {
    std::vector<int> domPath;

    vtkm::Particle seed = portal.Get(i);

    blockPopularity[blockID] += seed.GetNumberOfSteps();
    domPath.push_back(blockID);

    auto bids = boundsMap.FindBlocks(seed.GetPosition());
    if (bids.empty())
      continue;
    int bid = bids[0];

    // seed terminated inside the block.
    if (seed.GetStatus().CheckTerminate() || seed.GetNumberOfSteps() >= maxSteps || bid == blockID)
      continue;

    vtkm::FloatDefault numSteps = static_cast<vtkm::FloatDefault>(seed.GetNumberOfSteps());
    domPath.push_back(bid);
    int gid = blockInfo[bid]->GetLeaf(seed.GetPosition())->gid;

    bool inPingPong = false;
    while (true)
    {

      // Determine the destination based on a random percentage.
      auto it = flowMap.find(gid);
      if (it == flowMap.end())
      {
        auto p0 = seeds[i];
        auto p1 = seed;
        auto l = blockInfo[bid]->GetLeaf(seed.GetPosition());
        std::cout << p0 << " --> " << p1 << " into " << l->nm << std::endl;
        throw std::runtime_error("Source not found! " + std::to_string(gid));
      }
      FlowStat dstEntry = PickRandomWeightedDst(it->second);

      numSteps += dstEntry.avgSteps;
      blockPopularity[bid] += static_cast<int>(dstEntry.avgSteps + .5);
      particlesOut[bid]++;

      // Terminate or take max number of steps, we are done.
      int nextGID = dstEntry.dst;
      auto nextLeaf = DomainBlock::GetBlockFromGID(blockInfo, nextGID);

      if (nextGID == -1 || numSteps >= maxStepsFloat)
        break;

      particlesIn[nextLeaf->dom]++;
      bid = nextLeaf->dom;
      gid = nextGID;
      domPath.push_back(bid);
    }
    if (inPingPong)
    {
      std::cout << "PingPongCycle: " << domPath << std::endl;
    }

    // Detect any cycles.
    // DetectCycles(domPath, 2, cycleCnt);
    // DetectCycles(domPath, 3, cycleCnt);
    // DetectCycles(domPath, 4, cycleCnt);
    // DetectCycles(domPath, 5, cycleCnt);
    // DetectCycles(domPath, 6, cycleCnt);
  }
  ReduceToRoot(blockPopularity);
  ReduceToRoot(particlesIn);
  ReduceToRoot(particlesOut);
  ReduceToRoot(cycleCnt);
}


int main(int argc, char **argv)
{
  if (Rank == 0)
  {
    if (argc != 9)
    {
      std::cout << "<executable> <visitfileName> <fieldNm> <stepSize> <maxSteps> <numFacePts> <numTestPts> <Nxyz> <UseMultiStages>" << std::endl;
      exit(0);
    }
  }

  MPI_Init(&argc, &argv);
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  Rank = comm.rank();
  Size = comm.size();
  
  //declared in SLTool.h
  SLTOOL_Rank=Rank;
  PRINT_RANK = -1;
  PRINT_DETAILS = -1;
  SLTOOL_TOTALPROCS = Size;

  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);

  vtkm::filter::flow::GetTracer().Get()->Init(Rank);
  vtkm::filter::flow::GetTracer().Get()->ResetIterationStep(0);
  vtkm::filter::flow::GetTracer().Get()->StartTimer();
  vtkm::cont::Timer timer;

  auto MPIComm = vtkmdiy::mpi::mpi_cast(comm.handle());

  std::string visitfileName = argv[1];
  std::string fieldNm = argv[2];
  vtkm::FloatDefault stepSize = std::atof(argv[3]);
  vtkm::Id maxSteps = std::atoi(argv[4]);
  vtkm::Id numFacePts = std::atoi(argv[5]);
  vtkm::Id numTestPts = std::atoi(argv[6]);
  vtkm::Id Nxyz = std::atoi(argv[7]);
  std::string useMultiStagesStr = (argv[8]);
  
  if(useMultiStagesStr=="true"){
    UseMultiStages=true;
  }else{
    UseMultiStages=false;
  }
  vtkm::FloatDefault pctWidth = 0.10;

  if (Rank == 0)
  {
    std::cout << "Checking input parameters:" << std::endl;
    std::cout << " visitfileName: " << visitfileName << std::endl;
    std::cout << " fieldNm: " << fieldNm << std::endl;
    std::cout << " stepSize: " << stepSize << std::endl;
    std::cout << " maxSteps: " << maxSteps << std::endl;
    std::cout << " numFacePts: " << numFacePts << std::endl;
    std::cout << " numTestPts: " << numTestPts << std::endl;
    std::cout << " Nxyz: " << Nxyz << std::endl;
    std::cout << " UseMultiStages: " << UseMultiStages << std::endl;

  }

  std::vector<vtkm::cont::DataSet> dataSets;
  std::vector<int> blockIDList;
  AssignStrategy assignStrategy = AssignStrategy::CONTINUOUS;

  LoadData(visitfileName, assignStrategy, "", dataSets, blockIDList, Rank, Size);

  timer.Start();
  int numLocalBlocks = dataSets.size();
  int totNumBlocks = numLocalBlocks;
  MPI_Allreduce(MPI_IN_PLACE, &totNumBlocks, 1, MPI_INT, MPI_SUM, MPIComm);
  if (totNumBlocks != Size)
    throw std::runtime_error("Num blocks must be the same as number of ranks");
  
  std::cout << "debug rank " << Rank << " numLocalBlocks " << numLocalBlocks << std::endl;
  if(Rank==0){
    dataSets[0].PrintSummary(std::cout);
  }

  vtkm::filter::flow::internal::BoundsMap boundsMap(dataSets);
  std::vector<DomainBlock *> blockInfo;
  int NX = Nxyz, NY = Nxyz, NZ = Nxyz;
  bool subdivUniform = false;

    if(Rank==0){
        std::cout << "Debug CreateBlockInfo " << totNumBlocks << " " << subdivUniform << " "
         << NX << " " << NY << " " << NZ << " " << pctWidth << std::endl;
    }

  DomainBlock::CreateBlockInfo(blockInfo, totNumBlocks, boundsMap, subdivUniform, NX, NY, NZ, pctWidth);

  timer.Stop();
  double executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to CreateBlockInfo is: " << executionTime << std::endl;
  }
  timer.Start();

  std::vector<int> blockPopularity(Size, 0), particlesIn(Size, 0), particlesOut(Size, 0), cycleCnt(Size, 0);
  std::vector<int> blockAdvMultiStages(Size * 3, 0);

  int blockID = boundsMap.GetLocalBlockId(0);
  const auto &ds = dataSets[0];
  auto block = blockInfo[blockID];
  auto flowMap = CreateFlowMap(Nxyz, blockInfo, ds, blockID, fieldNm, stepSize, maxSteps, numFacePts, boundsMap);

  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  if (UseMultiStages)
  {
    CalcBlockPopularityMultiStages(blockInfo, ds, blockID, flowMap, boundsMap, blockAdvMultiStages, particlesIn, particlesOut, cycleCnt, numTestPts, fieldNm, stepSize, maxSteps);
  }
  else
  {
    CalcBlockPopularity(blockInfo, ds, blockID, flowMap, boundsMap, blockPopularity, particlesIn, particlesOut, cycleCnt, numTestPts, fieldNm, stepSize, maxSteps);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to CalcBlockPopularity is: " << executionTime << std::endl;
  }
  timer.Start();

  if (Rank == 0)
  {
    std::cout << std::endl
              << std::endl;
    if (UseMultiStages)
    {
      //std::cout << "blockAdvMultiStages: " << blockAdvMultiStages << std::endl;
      //TODO, write out the data as the format that can be recognized by the parser
      //std::cout << "blockAdvMultiStages: " << blockAdvMultiStages << std::endl;
      // normalize..
      int sum = std::accumulate(blockAdvMultiStages.begin(), blockAdvMultiStages.end(), 0);
      std::vector<float> blockPopNorm;
      for (const auto &v : blockAdvMultiStages){
          blockPopNorm.push_back((float)(v) / sum);
      }
      //Make sure the output can be parsed by json.
      std::cout <<"MultiStagesBlockPopularityNorm:[";
      for (int i=0;i<blockPopNorm.size()/3;i++){
        int startPos = i*3;
        if(i>0){
          std::cout <<",";
        }
        std::cout <<"[" << blockPopNorm[startPos]<<","<<blockPopNorm[startPos+1]<<","<<blockPopNorm[startPos+2]<<"]";
      }
      std::cout<<"]"<<std::endl;
    }
    else
    {
      std::cout << "BlockPopularity: " << blockPopularity << std::endl;
      // normalize..
      int sum = std::accumulate(blockPopularity.begin(), blockPopularity.end(), 0);
      std::vector<float> blockPopNorm;
      for (const auto &v : blockPopularity)
        blockPopNorm.push_back((float)(v) / sum);
      std::cout << "NormBlockPopularity: " << blockPopNorm << std::endl;

      std::vector<float> particlesInNorm;
      sum = std::accumulate(particlesIn.begin(), particlesIn.end(), 0);
      for (const auto &v : particlesIn)
        particlesInNorm.push_back((float)(v) / sum);

      std::vector<float> particlesOutNorm;
      sum = std::accumulate(particlesOut.begin(), particlesOut.end(), 0);
      for (const auto &v : particlesOut)
        particlesOutNorm.push_back((float)(v) / sum);

      std::vector<float> particlesInOutNorm;
      std::vector<int> particlesInOut;
      sum = std::accumulate(particlesIn.begin(), particlesIn.end(), 0) + std::accumulate(particlesOut.begin(), particlesOut.end(), 0);
      for (int i = 0; i < particlesIn.size(); i++)
      {
        int v = particlesIn[i] + particlesOut[i];
        particlesInOut.push_back(v);
        particlesInOutNorm.push_back((float)(v) / sum);
      }

      std::vector<float> cycleCntNorm;
      sum = std::accumulate(cycleCnt.begin(), cycleCnt.end(), 0);
      for (const auto &v : cycleCnt)
        cycleCntNorm.push_back((float)(v) / (float)sum);

      std::cout << "ParticlesIn:  " << particlesIn << std::endl;
      std::cout << "NormParticlesIn:  " << particlesInNorm << std::endl;

      std::cout << "ParticlesOut:  " << particlesIn << std::endl;
      std::cout << "NormParticlesOut: " << particlesOutNorm << std::endl;

      std::cout << "ParticlesInOut: " << particlesInOut << std::endl;
      std::cout << "NormParticlesInOut: " << particlesInOutNorm << std::endl;

      std::cout << "CycleCnt:     " << cycleCnt << std::endl;
      std::cout << "NormCycleCnt: " << cycleCntNorm << std::endl;

      // output to a single file.
      //  std::ofstream outF("table.txt", std::ofstream::out);
      //  outF<<"Block, popN, pinN, poutN, pinoutN"<<std::endl;
      //  for (int i = 0; i < Size; i++)
      //  {
      //    outF<<i<<", "<<blockPopNorm[i]<<", "<<particlesInNorm[i]<<", "<<particlesOutNorm[i]<<", "<<particlesInOutNorm[i]<<std::endl;
      //  }
    }
  }

  MPI_Finalize();
  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (Rank == 0)
  {
    std::cout << "Execution time to Output log is: " << executionTime << std::endl;
  }
}

// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecXY >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecX >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit velocity >&out

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

#include <mpi.h>
#include <vtkm/thirdparty/diy/diy.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>

#include "Block.h"
#include "../LoadData.hpp"
#include "../AssignStrategy.hpp"

int Rank, Size;
int NUMVALS = -1;

const int PRINT_RANK = 0;


static vtkm::FloatDefault random_1()
{
  return (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
}

/*
void LoadData(std::string &fname, std::string &fieldNm, std::vector<vtkm::cont::DataSet> &dataSets, int rank, int nRanks)
{
  std::string buff;
  std::ifstream is;
  is.open(fname);
  std::cout << "Opening: " << fname << std::endl;
  if (!is)
  {
    std::cout << "File not found! : " << fname << std::endl;
    throw "unknown file: " + fname;
  }

  auto p0 = fname.rfind(".visit");
  if (p0 == std::string::npos)
    throw "Only .visit files are supported.";
  auto tmp = fname.substr(0, p0);
  auto p1 = tmp.rfind("/");
  auto dir = tmp.substr(0, p1);

  std::getline(is, buff);
  auto numBlocks = std::stoi(buff.substr(buff.find("!NBLOCKS ") + 9, buff.size()));
  if (rank == 0)
    std::cout << "numBlocks= " << numBlocks << std::endl;

  int nPer = numBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;
  if (rank == (nRanks - 1))
    b1 = numBlocks;

  for (int i = 0; i < numBlocks; i++)
  {
    std::getline(is, buff);
    if (i >= b0 && i < b1)
    {
      vtkm::cont::DataSet ds;
      std::string vtkFile = dir + "/" + buff;
      vtkm::io::VTKDataSetReader reader(vtkFile);
      ds = reader.ReadDataSet();
      auto f = ds.GetField(fieldNm).GetData();
      vtkm::cont::ArrayHandle<vtkm::Vec<double, 3>> fieldArray;
      f.AsArrayHandle(fieldArray);
      int n = fieldArray.GetNumberOfValues();
      auto portal = fieldArray.WritePortal();
      for (int ii = 0; ii < n; ii++)
        portal.Set(ii, vtkm::Vec<double, 3>(1, 0, 0));

      dataSets.push_back(ds);
    }
  }
}
*/

// Example computing streamlines.
// An example vector field is available in the vtk-m data directory: magField.vtk
// Example usage:
//   this will advect 200 particles 50 steps using a step size of 0.01
//
// Particle_Advection <path-to-data-dir>/magField.vtk vec 200 50 0.01 output.vtk
//

void GenerateTestPts(DomainBlock *leaf, int numPts, std::vector<vtkm::Particle> &pts)
{
  vtkm::FloatDefault bb[6];
  leaf->GetBBox(bb);
  vtkm::FloatDefault dx = bb[1] - bb[0], dy = bb[3] - bb[2], dz = bb[5] - bb[4];

  for (int i = 0; i < numPts; i++)
  {
    vtkm::Vec3f pt(bb[0] + random_1() * dx,
                   bb[2] + random_1() * dy,
                   bb[4] + random_1() * dz);
    // Doesn't matter what the ID is....
    pts.push_back(vtkm::Particle(pt, i));
  }
}

class nextBlock
{
public:
  nextBlock()
  {
    cnt = 0;
    numIters = 0;
  }

  void visit(int i)
  {
    cnt++;
    numIters += i;
  }

  int cnt, numIters;
};

void RunTestPts(std::vector<DomainBlock *> blockInfo,
                vtkm::Id blkId,
                DomainBlock *leaf,
                const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                std::vector<vtkm::Particle> &seeds,
                const vtkm::cont::DataSet &ds,
                const std::string &fieldNm,
                std::vector<int> &leafData,
                vtkm::FloatDefault stepSize,
                vtkm::Id maxSteps)
{
  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::flow::VelocityField<ArrayType>;
  using GridEvalType = vtkm::worklet::flow::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::flow::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::flow::Stepper<RK4Type, GridEvalType>;

  auto seedsCopy = seeds;

  auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::Off);

  ArrayType arr;
  vtkm::cont::ArrayCopyShallowIfPossible(ds.GetField(fieldNm).GetData(), arr);
  FieldType velocities(arr, ds.GetField(fieldNm).GetAssociation());
  GridEvalType gridEval(ds, velocities);
  Stepper rk4(gridEval, stepSize);

  vtkm::worklet::flow::ParticleAdvection pa;
  auto res = pa.Run(rk4, seedArray, maxSteps);

  std::map<int, nextBlock *> nextBlocks;

  auto ptsPortal = res.Particles.ReadPortal();
  vtkm::Id totNumSteps = 0, totNumPts = 0;

  for (vtkm::Id i = 0; i < ptsPortal.GetNumberOfValues(); i++)
  {
    auto pt = ptsPortal.Get(i);
    auto pt0 = seedsCopy[i];
    auto destinations = boundsMap.FindBlocks(pt.GetPosition(), {blkId});
    vtkm::Id dst = -1;
    DomainBlock *dstLeaf = nullptr;

    if (destinations.empty())
    {
      auto ptr = leaf->parent;
      while (ptr->parent != nullptr)
        ptr = ptr->parent;
      dst = ptr->gid;
      dstLeaf = ptr;
    }
    else
    {
      dst = destinations[0];
      dstLeaf = blockInfo[dst]->GetLeaf(pt.GetPosition());
    }

    nextBlock *n;
    auto it = nextBlocks.find(dstLeaf->gid);
    if (it == nextBlocks.end())
    {
      n = new nextBlock();
      nextBlocks[dstLeaf->gid] = n;
    }
    else
      n = it->second;

    n->visit(pt.GetNumberOfSteps());
    totNumSteps += pt.GetNumberOfSteps();
    totNumPts++;

    if (Rank == PRINT_RANK)
    {
      if (i == 0) std::cout<<" *** ID, pos0, posn, #steps, id0-->idn, leaf0-->leafn"<<std::endl;
      std::cout << pt.GetID() << " : " << pt0.GetPosition() << " --> " << pt.GetPosition() << ", " << pt.GetNumberOfSteps();
      std::cout << " :: (" << blkId << " --> " << dst << ")";
      std::cout << " (" << leaf->gid << " --> ";
      if (dstLeaf)
        std::cout << dstLeaf->gid << ")" << std::endl;
      else
        std::cout << "-1)" << std::endl;
    }
  }

  // Put the nextBlocks into an array.
  int index = 0;
  leafData.resize(nextBlocks.size() * 4, -1);
  for (auto it = nextBlocks.begin(); it != nextBlocks.end(); it++)
  {
    auto n = it->second;
    leafData[index++] = it->first;
    leafData[index++] = n->cnt;
    leafData[index++] = n->numIters;
    leafData[index++] = totNumPts;

    delete n;
  }
  if (index >= NUMVALS)
  {
    throw "memory overflow in leafData. NUMVALS is too small!";
  }

  nextBlocks.clear();
}

vtkm::cont::DataSet
CreateUniformDataSet(const vtkm::Bounds &bounds,
                     const vtkm::Id3 &dims)
{
  vtkm::Vec3f origin(static_cast<vtkm::FloatDefault>(bounds.X.Min),
                     static_cast<vtkm::FloatDefault>(bounds.Y.Min),
                     static_cast<vtkm::FloatDefault>(bounds.Z.Min));
  vtkm::Vec3f spacing(static_cast<vtkm::FloatDefault>(bounds.X.Length()) /
                          static_cast<vtkm::FloatDefault>((dims[0] - 1)),
                      static_cast<vtkm::FloatDefault>(bounds.Y.Length()) /
                          static_cast<vtkm::FloatDefault>((dims[1] - 1)),
                      static_cast<vtkm::FloatDefault>(bounds.Z.Length()) /
                          static_cast<vtkm::FloatDefault>((dims[2] - 1)));

  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet ds = dataSetBuilder.Create(dims, origin, spacing);
  return ds;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  Rank = comm.rank();
  Size = comm.size();

  vtkm::filter::flow::GetTracer().Get()->Init(Rank);
  vtkm::filter::flow::GetTracer().Get()->ResetIterationStep(0);
  vtkm::filter::flow::GetTracer().Get()->StartTimer();
  auto MPIComm = vtkmdiy::mpi::mpi_cast(comm.handle());

  vtkm::FloatDefault stepSize = 0.01;
  vtkm::Id maxSteps = 1000;

  std::string visitfileName = argv[1];
  std::string fieldNm = argv[2];
  std::vector<vtkm::cont::DataSet> dataSets;
  std::vector<int> blockIDList;
  AssignStrategy assignStrategy = AssignStrategy::CONTINUOUS;

  LoadData(visitfileName, assignStrategy, "", dataSets, blockIDList, Rank, Size);

  int numLocalBlocks = dataSets.size();
  int totNumBlocks = numLocalBlocks;
  MPI_Allreduce(MPI_IN_PLACE, &totNumBlocks, 1, MPI_INT, MPI_SUM, MPIComm);
  if (Rank == 0) std::cout<<" Tot numBlocks= "<<totNumBlocks<<std::endl;

  vtkm::filter::flow::internal::BoundsMap boundsMap(dataSets);
  std::vector<DomainBlock *> blockInfo;
  int NX = 2, NY = 2, NZ = 2;
  bool subdivUniform = false;
  DomainBlock::CreateBlockInfo(blockInfo, totNumBlocks, boundsMap, subdivUniform, NX, NY, NZ, 0.10f);

  int nVals = 4;      //(dst, numPts, numIters, totalPtsFromSrc)
  int nNeighbors = 7; // 6 neighboring blocks, plus self.
  int nSubdiv = std::max(std::max(NX, NY), NZ);
  if (!subdivUniform)
    nSubdiv++;
  NUMVALS = nVals * (nNeighbors * nSubdiv);

  //if (Rank == 0) DomainBlock::Dump(blockInfo, std::cout, 0);
  //MPI_Barrier(MPI_COMM_WORLD);

  int totNumLeafs = DomainBlock::TotalNumLeaves(blockInfo);
  int numLeafs = blockInfo[0]->NumLeafs();
  // std::vector<std::vector<int>> blockData(numLeafs);
  // for (auto& bd : blockData)
  // bd.resize(NUMVALS, -1);
  /*
  int **blockData = new int*[totNumLeafs];
  for (int i = 0; i < totNumLeafs; i++)
  {
    blockData[i] = new int[NUMVALS];
    for (int j = 0; j < NUMVALS; j++)
      blockData[i][j] = -1;
  }
  */

  std::vector<int> allLeafData(totNumLeafs * NUMVALS, -1);

  // generate test points.
  const int NUM_PTS = 10;
  int totalNumPts = 0;
  for (int i = 0; i < numLocalBlocks; i++)
  {
    //int blockID = boundsMap.GetBlockIdFromLocalIndex(i);
    int blockID = boundsMap.GetLocalBlockId(i);

    const auto &ds = dataSets[i];
    auto block = blockInfo[blockID];

    for (int j = 0; j < numLeafs; j++)
    {
      std::vector<vtkm::Particle> seeds;
      auto leaf = block->GetLeafFromIndex(j);
      if (leaf->leafBlockType == DomainBlock::INTERNAL)
      {
        // for internal points, we really want to add the actual seeds.
        GenerateTestPts(leaf, NUM_PTS, seeds);
        totalNumPts += NUM_PTS;
      }
      else
      {
        GenerateTestPts(leaf, NUM_PTS, seeds);
        totalNumPts += NUM_PTS;
      }
      std::vector<int> leafData;
      RunTestPts(blockInfo, blockID, leaf, boundsMap, seeds, dataSets[i], fieldNm, leafData, stepSize, maxSteps);

      // print out leafData
      int sz = leafData.size();
      if (Rank == PRINT_RANK)
      {
        std::cout << "leafData[" << j << "]= (dstLeaf, #pts, #steps, totPtsFromSrc)" << std::endl;
        for (int k = 0; k < sz; k += 4)
        {
          std::cout << "  ";
          std::cout << (k / 4) << " : ";
          std::cout << leafData[k + 0] << " " << leafData[k + 1] << " " << leafData[k + 2] << " " << leafData[k + 3];
          std::cout << std::endl;
        }
      }

      // Update the global block info
      int idx = blockID * numLeafs * NUMVALS + (j * NUMVALS);
      for (int k = 0; k < sz; k++)
        allLeafData[idx + k] = leafData[k];
      if (Rank == PRINT_RANK) std::cout << std::endl;
    }
  }

  if (Size > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, allLeafData.data(), allLeafData.size(), MPI_INT, MPI_MAX, MPIComm);
  }

  std::vector<DomainBlock*> allDomainBlocks;
  if (Rank == PRINT_RANK)
    std::cout<<"******* Printing ALL Leaves *******: tot= "<<totNumLeafs<<std::endl;
  for (int i = 0; i < totNumLeafs; i++)
  {
    auto leaf = DomainBlock::GetBlockFromGID(blockInfo, i);
    DomainBlock tmp = *leaf;

    if (Rank == PRINT_RANK)
      std::cout << "leafData[" << i << "]("<<leaf->nm<<")= (idx : dstLeaf, #pts, #steps, totPtsFromSrc, dstLeafName)"<<std::endl;
    int idx = i * NUMVALS;
    for (int j = 0; j < NUMVALS; j += 4)
    {
      if (allLeafData[idx + j + 1] == -1)
        break;

      auto dstLeaf = DomainBlock::GetBlockFromGID(blockInfo, allLeafData[idx + j + 0]);
      int dstGID = allLeafData[idx+j+0];
      int numICs = allLeafData[idx+j+1];
      int numIters = allLeafData[idx+j+2];
      int totNumICs = allLeafData[idx+j+3];
      leaf->AddBlockData(dstLeaf, 1,2,3); //numICs, numIters, totNumICs);
      if (Rank == PRINT_RANK)
      {
        std::cout << "  ";
        std::cout << (j / 4) << " : ";
        std::cout<<dstGID<<" "<<numICs<<" "<<numIters<<" "<<totNumICs;
        if (dstLeaf == nullptr)
          std::cout << " INTERNAL ";
        else
          std::cout << "  " << dstLeaf->nm;
        std::cout << std::endl;
      }
      if (dstGID == -1)
        dstLeaf = leaf->GetInternal();
      leaf->AddBlockData(dstLeaf, numICs, numIters, totNumICs);
    }
    if (Rank == PRINT_RANK)
      std::cout<<"UNIFIED DATA "<<std::endl;
    leaf->UnifyData();
    if (Rank == PRINT_RANK)
    {
      DomainBlock::Dump(leaf, std::cout, 1);
      std::cout<<std::endl<<std::endl;
    }
  }

  //Put allLeafData back into blockInfo.
  for (int i = 0; i < totNumLeafs; i++)
  {
    auto leaf = DomainBlock::GetBlockFromGID(blockInfo, i);
    if (leaf->sub == 0)
      allDomainBlocks.push_back(leaf);
  }

  if (Rank == PRINT_RANK)
  {
    std::cout<<std::endl<<std::endl;
    std::cout<<"********************  blockInfo ********************"<<std::endl;
    DomainBlock::Dump(allDomainBlocks, std::cout, 2);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return(0);



  /*
  //print out blockData
  for (int i = 0; i < totNumLeafs; i++)
  {
    std::cout<<"leafData["<<i<<"]= (dstLeaf, #pts, #steps, totPtsFromSrc)"<<std::endl;
    for (int j = 0; j < NUMVALS; j+=4)
    {
      std::cout<<"  ";
      std::cout<<(j/4)<<" : ";
      std::cout<<leafData[i][j+0]<<" "<<leafData[i][j+1]<<" "<<leafData[i][j+2]<<" "<<leafData[i][j+3];
      std::cout<<std::endl;
      int numLeafs = 25;
      int idx = i*numLeafs*NUMVALS;
      std::cout<<"      "<<allLeafData[idx+j+0]<<" "<<allLeafData[idx+j+1]<<" "<<allLeafData[idx+j+2]<<" "<<allLeafData[idx+j+3]<<std::endl;
    }
  */

#if 0
  vtkm::filter::flow::ParticleAdvection pa;

  vtkm::cont::ArrayHandle<vtkm::Particle> seedArray;
  seedArray = vtkm::cont::make_ArrayHandle({ vtkm::Particle(vtkm::Vec3f(.1f, .1f, .9f), 0),
                                             vtkm::Particle(vtkm::Vec3f(.1f, .6f, .6f), 1),
                                             vtkm::Particle(vtkm::Vec3f(.1f, .9f, .1f), 2) });
  pa.SetStepSize(0.001f);
  pa.SetNumberOfSteps(10000);
  pa.SetSeeds(seedArray);
  pa.SetActiveField(fieldNm);

  vtkm::cont::PartitionedDataSet pds(dataSets);
  auto output = pa.Execute(pds);
  output.PrintSummary(std::cout);
#endif

  MPI_Finalize();
  return 0;
}


//mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecXY >&out
//mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecX >&out
//mpirun -np 8 ./StreamlineMPI ./data/clover.visit velocity >&out

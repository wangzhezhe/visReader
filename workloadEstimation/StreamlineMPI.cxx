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
//using long int to prevent the `larger than max_size()` error for some configurations
long int NUMVALS = -1;

const int PRINT_RANK = -1;

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
  // get results information after particle advection
  auto ptsPortal = res.Particles.ReadPortal();
  vtkm::Id totNumSteps = 0, totNumPts = 0;
  // go through each particle
  for (vtkm::Id i = 0; i < ptsPortal.GetNumberOfValues(); i++)
  {
    auto pt = ptsPortal.Get(i);
    auto pt0 = seedsCopy[i];
    // find the destination of current particle according to particle's position
    // ignore the current block id
    auto destinations = boundsMap.FindBlocks(pt.GetPosition(), blkId);
    vtkm::Id dst = -1;
    DomainBlock *dstLeaf = nullptr;

    if (destinations.empty())
    {
      // if the destination is empty, maybe the particle still locate in the same block
      // get the parent node recursively until the one at the top level
      auto ptr = leaf->parent;
      while (ptr->parent != nullptr)
        ptr = ptr->parent;
      // the destination is the global id of that parent
      // the destination leaf is still the same ptr node
      dst = ptr->gid;
      dstLeaf = ptr;
    }
    else
    {
      // if the destination is not empty, the dst leaf is the leaf
      // of associated destination block
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
    // update the the next block value
    n->visit(pt.GetNumberOfSteps());
    totNumSteps += pt.GetNumberOfSteps();
    totNumPts++;

    if (Rank == PRINT_RANK)
    {
      // id of the particle, start position(seeds position), end position, #of integration steps
      // block id of start position, block id of end position, id of start leaf, id of end leaf
      if (i == 0)
        std::cout << " *** ID, pos0-->posn, #steps, id0-->idn (large block id), leaf0-->leafn (global subdomain leaf id)" << std::endl;
      std::cout << pt.GetID() << "," << pt0.GetPosition() << "-->" << pt.GetPosition() << "," << pt.GetNumberOfSteps();
      std::cout << ",(" << blkId << "-->" << dst << ")";
      std::cout << ",(" << leaf->gid << "(" << leaf->nm << ")-->" << dstLeaf->gid << "(";
      if (dstLeaf->gid == -1)
        std::cout << "NULL";
      else
        std::cout << dstLeaf->nm;
      std::cout << "))" << std::endl;
    }
  }

  // Put the nextBlocks into an array.
  // This leafData is empty when calling the RunTestPts function
  // It stors next blocks information for current subdomain
  // for the nextBlocks, the key is the global id of block, the value is the nextBlock entity
  int index = 0;
  leafData.resize(nextBlocks.size() * 4, -1);
  for (auto it = nextBlocks.begin(); it != nextBlocks.end(); it++)
  {
    auto n = it->second;
    // global id of next one
    leafData[index++] = it->first;
    // the number of particles in the next block, this increase one when there is a new particle go to next block
    leafData[index++] = n->cnt;
    // sum of steps for all particles to go to the next block
    // all steps for particles goes into the next block are added into this value
    leafData[index++] = n->numIters;
    // total number of points
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

void CalcuateBlockPopularity(std::vector<DomainBlock *> &blockInfo,
                             const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                             std::vector<int> &blockPopularity,
                             std::vector<int> &particlesIn,
                             std::vector<int> &particlesOut,
                             std::vector<int> &cycleCnt,
                             int numPts, int maxSteps)
{
  int blockID = Rank;
  auto block = blockInfo[blockID];

  std::vector<vtkm::Particle> seeds;
  GenerateTestPts(block, numPts, seeds);

  int seedCnt = 0;
  for (const auto &seed : seeds)
  {
    auto destinations = boundsMap.FindBlocks(seed.GetPosition());
    if (destinations.empty())
      continue;

    int dst = destinations[0];
    DomainBlock *dstLeaf = blockInfo[dst]->GetLeaf(seed.GetPosition());

    vtkm::FloatDefault numSteps = static_cast<vtkm::FloatDefault>(maxSteps);

    const int winSZ = 2;
    std::vector<int> window(winSZ, -1);
    std::vector<int> domPath;
    while (true)
    {
      domPath.push_back(dstLeaf->dom);
      // Determine the destination based on a random percentage.
      vtkm::FloatDefault pct = random_1();

      // make sure this seed goes to which direaction
      int idx = dstLeaf->GetDataIdxFromPct(pct, Rank);
      // if (Rank == PRINT_RANK && seedId >= 0)
      // {
      //   std::cout << "debug pct " << pct << " destidx " << idx << std::endl;
      // }
      if (idx == -1)
        break;
      const DomainBlock::blockData &data = dstLeaf->data[idx];

      // Take avgIt steps in this block.
      // dom represent the original block id before subdivision
      blockPopularity[dstLeaf->dom] += static_cast<int>(data.avgIt + 0.5);
      // blockPopularity[dstLeaf->dom] += static_cast<int>(data.avgIt );
      numSteps -= data.avgIt;
      particlesOut[dstLeaf->dom]++;

      // Next destination.
      auto nextLeaf = data.blk;

      // If we took the max number of steps, or we terminated (ended up in the same block), we are done.
      if (numSteps <= 0 || nextLeaf->leafBlockType == DomainBlock::INTERNAL)
        break;

      // continue in the next leaf.
      particlesIn[nextLeaf->dom]++;
      dstLeaf = nextLeaf;
    }

    // if (Rank == 0 && seedCnt == 0)
    {
      // std::cout<<"DOMPath: "<<domPath<<std::endl;
      DetectCycles(domPath, 2, cycleCnt);
      // std::cout<<"Cycles: "<<cycleCnt<<std::endl;
    }
    seedCnt++;
  }

  // Calculate the block popularity over all ranks.
  MPI_Allreduce(MPI_IN_PLACE, blockPopularity.data(), blockPopularity.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, particlesIn.data(), particlesIn.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, particlesOut.data(), particlesOut.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, cycleCnt.data(), cycleCnt.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
  if (Rank == 0)
  {
    if (argc != 9)
    {
      std::cout << "<executable> <visitfileName> <fieldNm> <stepSize> <maxSteps> <NUM_TEST_POINTS> <NUM_SIM_POINTS_PER_DOM> <Nxyz> <pctWidth>" << std::endl;
      exit(0);
    }
  }

  MPI_Init(&argc, &argv);
  auto comm = vtkm::cont::EnvironmentTracker::GetCommunicator();
  Rank = comm.rank();
  Size = comm.size();

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
  vtkm::Id NUM_TEST_POINTS = std::atoi(argv[5]);
  vtkm::Id NUM_SIM_POINTS_PER_DOM = std::atoi(argv[6]);
  vtkm::Id Nxyz = std::atoi(argv[7]);
  vtkm::FloatDefault pctWidth = std::stod(argv[8]);

  if (pctWidth <= 0 || pctWidth >= 0.5)
  {
    std::cout << "pctWidth shoule be a value between 0 to 0.5" << std::endl;
    exit(0);
  }

  if (Rank == 0)
  {
    std::cout << "Checking input parameters:"
              << "\n visitfileName: " << visitfileName << "\n fieldNm: " << fieldNm << "\n stepSize: " << stepSize << "\n maxSteps: " << maxSteps << "\n NUM_TEST_POINTS: " << NUM_TEST_POINTS << "\n NUM_SIM_POINTS_PER_DOM: " << NUM_SIM_POINTS_PER_DOM << "\n Nxyz: " << Nxyz << "\n pctWidth: " << pctWidth << std::endl;
  }

  std::vector<vtkm::cont::DataSet> dataSets;
  std::vector<int> blockIDList;
  AssignStrategy assignStrategy = AssignStrategy::CONTINUOUS;

  LoadData(visitfileName, assignStrategy, "", dataSets, blockIDList, Rank, Size);
  
  timer.Start();
  int numLocalBlocks = dataSets.size();
  int totNumBlocks = numLocalBlocks;
  MPI_Allreduce(MPI_IN_PLACE, &totNumBlocks, 1, MPI_INT, MPI_SUM, MPIComm);
  if (Rank == 0)
    std::cout << " Tot numBlocks= " << totNumBlocks << std::endl;
  // index the block in distributed way
  vtkm::filter::flow::internal::BoundsMap boundsMap(dataSets);
  std::vector<DomainBlock *> blockInfo;
  int NX = Nxyz, NY = Nxyz, NZ = Nxyz;
  bool subdivUniform = false;
  DomainBlock::CreateBlockInfo(blockInfo, totNumBlocks, boundsMap, subdivUniform, NX, NY, NZ, pctWidth);

  int nVals = 4;      //(dst, numPts, numIters, totalPtsFromSrc)
  int nNeighbors = 7; // 6 neighboring blocks, plus self.
  // all possible destinations
  int nSubdiv = std::max(std::max(NX, NY), NZ);
  // square the value since for each boundary region, we divide it through two dimentions
  // for each slab such as x, it is divided into ny*nz small blocks, for each small block
  // there are 7 adjacent neighbors(including itsself)
  // we accume ny=nx=nz here
  nSubdiv = nSubdiv * nSubdiv;

  int numSmallBlocksPerBlock = 6 * nSubdiv;

  if (!subdivUniform)
  {
    // there is 1 sub block for internal region
    numSmallBlocksPerBlock++;
  }
  else
  {
    numSmallBlocksPerBlock += nSubdiv;
  }

  // nSubdiv represents small blocks in each large block
  if (Rank == 0)
  {
    std::cout << "Checking numSmallBlocksPerBlock: " << numSmallBlocksPerBlock << std::endl;
  }
  // Total values for metadata for all subdomains
  // This value shows how many values are stored for each leaf node
  NUMVALS = nVals * totNumBlocks * numSmallBlocksPerBlock;

  // if (Rank == 0) DomainBlock::Dump(blockInfo, std::cout, 0);
  // MPI_Barrier(MPI_COMM_WORLD);
  // Go through all subdomains and compute the number of leaves
  int totNumLeafs = DomainBlock::TotalNumLeaves(blockInfo);
  // The number of leaf for the first block
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
  // assign vector that store all metadata for leaves
  std::vector<int> allLeafData(totNumLeafs * NUMVALS, -1);

  // generate test points.
  int totalNumPts = 0;
  for (int i = 0; i < numLocalBlocks; i++)
  {
    // int blockID = boundsMap.GetBlockIdFromLocalIndex(i);
    int blockID = boundsMap.GetLocalBlockId(i);

    const auto &ds = dataSets[i];
    auto block = blockInfo[blockID];

    // For the leaf associated with each block
    // generate seeds, which are NUM_TEST_POINTS
    // this number of seeds will be placed in each leaf no matter it is bounds or internal
    for (int j = 0; j < numLeafs; j++)
    {
      std::vector<vtkm::Particle> seeds;
      auto leaf = block->GetLeafFromIndex(j);
      if (leaf->leafBlockType == DomainBlock::INTERNAL)
      {
        // for internal points, we really want to add the actual seeds.
        // make the seeds number in interal region match with it region
        // we use the smallest region as the unit the complete equattion is (1-2p)^3/[(1-2p)^2*p/Nxyz]
        vtkm::Id scale = Nxyz*((1.0 - 2.0 * pctWidth) / pctWidth);
        scale = vtkm::Max(1.0, scale);
        vtkm::Id NumPtsInternal = scale * NUM_TEST_POINTS;
        GenerateTestPts(leaf, NumPtsInternal, seeds);
        totalNumPts += NumPtsInternal;
        if (Rank == PRINT_RANK)
        {
          std::cout << "Checking seeds, scale is "<< scale << " #internal block seeds: " << NumPtsInternal << std::endl;
        }
      }
      else
      {
        GenerateTestPts(leaf, NUM_TEST_POINTS, seeds);
        totalNumPts += NUM_TEST_POINTS;
        if (Rank == PRINT_RANK)
        {
          std::cout << "Checking seeds, leaf block: " << leaf->nm << " #seeds " << NUM_TEST_POINTS << std::endl;
        }
      }
      std::vector<int> leafData;
      // Run the particle advection based on current seeds and blocks
      // The execution here is serialized way, we go through each leaf region
      // and run the particle advection within each leaf region
      RunTestPts(blockInfo, blockID, leaf, boundsMap, seeds, dataSets[i], fieldNm, leafData, stepSize, maxSteps);

      // print out leafData (next block data)
      int sz = leafData.size();
      if (Rank == PRINT_RANK)
      {
        // j represent the id of leaf
        std::cout << "block id " << blockID << " leafName " << leaf->nm << " leafData[" << j << "]= (dstLeafGlobalId, #pts, #steps, totPtsFromSrc)" << std::endl;
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
      if (Rank == PRINT_RANK)
        std::cout << std::endl;
    }
  }

  if (Size > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, allLeafData.data(), allLeafData.size(), MPI_INT, MPI_MAX, MPIComm);
  }

  std::vector<DomainBlock *> allDomainBlocks;
  if (Rank == PRINT_RANK)
    std::cout << "******* Printing ALL Leaves *******: tot= " << totNumLeafs << std::endl;
  for (int i = 0; i < totNumLeafs; i++)
  {
    auto leaf = DomainBlock::GetBlockFromGID(blockInfo, i);
    DomainBlock tmp = *leaf;

    if (Rank == PRINT_RANK)
      std::cout << "leafData[" << i << "](" << leaf->nm << ")= (idx : dstLeaf, #pts, #steps, totPtsFromSrc, dstLeafName)" << std::endl;
    int idx = i * NUMVALS;
    for (int j = 0; j < NUMVALS; j += 4)
    {
      if (allLeafData[idx + j + 1] == -1)
        break;

      auto dstLeaf = DomainBlock::GetBlockFromGID(blockInfo, allLeafData[idx + j + 0]);
      int dstGID = allLeafData[idx + j + 0];
      int numICs = allLeafData[idx + j + 1];
      int numIters = allLeafData[idx + j + 2];
      int totNumICs = allLeafData[idx + j + 3];
      leaf->AddBlockData(dstLeaf, numICs, numIters, totNumICs);
      if (Rank == PRINT_RANK)
      {
        std::cout << "  ";
        std::cout << (j / 4) << " : ";
        std::cout << dstGID << " " << numICs << " " << numIters << " " << totNumICs;
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
      std::cout << "UNIFIED DATA " << std::endl;
    // compute the percentage here
    leaf->UnifyData();
    if (Rank == PRINT_RANK)
    {
      DomainBlock::Dump(leaf, std::cout, 1);
      std::cout << std::endl
                << std::endl;
    }
  }

  // Put allLeafData back into blockInfo.
  for (int i = 0; i < totNumLeafs; i++)
  {
    auto leaf = DomainBlock::GetBlockFromGID(blockInfo, i);
    if (leaf->sub == 0)
      allDomainBlocks.push_back(leaf);
  }

  if (Rank == PRINT_RANK)
  {
    // check output information for all blocks, how much possibility to go to the dest region
    std::cout << std::endl
              << std::endl;
    // The three digit at the end of the each boundary block represent the index of the subblock
    // if the dx is 2 for x, the first digit is 0 or 1, if the dx is 3, the first digit is 0 1 2
    std::cout << "********************  blockInfo ********************" << std::endl;
    DomainBlock::Dump(allDomainBlocks, std::cout, 2);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<int> blockPopularity(Size, 0), particlesIn(Size, 0), particlesOut(Size, 0), cycleCnt(Size, 0);
  CalcuateBlockPopularity(allDomainBlocks, boundsMap,
                          blockPopularity, particlesIn, particlesOut, cycleCnt,
                          NUM_SIM_POINTS_PER_DOM, maxSteps);

  if (Rank == 0)
  {
    std::cout << std::endl
              << std::endl;
    std::cout << "BlockPopularity: " << blockPopularity << std::endl;

    // normalize..
    int sum = std::accumulate(blockPopularity.begin(), blockPopularity.end(), 0);
    std::vector<float> blockPopNorm;
    for (const auto &v : blockPopularity)
      blockPopNorm.push_back((float)(v) / sum);
    std::cout << "NormBlockPopularity: " << blockPopNorm << std::endl;

    std::cout << "ParticlesIn:  " << particlesIn << std::endl;
    std::cout << "ParticlesOut: " << particlesOut << std::endl;
    std::cout << "CycleCnt:     " << cycleCnt << std::endl;
  }

  /*
    if (Rank == PRINT_RANK)
    {

      std::cout<<"************** WALK PARTICLES ******************"<<std::endl;
      vtkm::Vec3f pt(0.122048,0.835823,1.01106);
      pt = vtkm::Vec3f(1.71681,1.01744,1.54469);

      auto destinations = boundsMap.FindBlocks(pt);
      vtkm::Id dst = destinations[0];
      vtkm::FloatDefault maxSteps = 1000;

      DomainBlock* dstLeaf = allDomainBlocks[dst]->GetLeaf(pt);
      std::cout<<"PT= "<<pt<<" :: "<<dstLeaf->nm<<std::endl;

      std::vector<int> blockPopularity(Size, 0);
      int cnt = 0;
      while (true)
      {
        vtkm::FloatDefault pct = random_1();
        int idx = dstLeaf->GetDataIdxFromPct(pct);
        if (idx == -1)
          break;

        std::cout<<" DstLeaves:  ";
        DomainBlock::Dump(dstLeaf, std::cout, 10); std::cout<<std::endl;
        std::cout<<"  Pct= "<<pct<<" idx= "<<idx<<std::endl;
        const DomainBlock::blockData& data = dstLeaf->data[idx];
        std::cout<<"    Take steps: "<<data.avgIt<<" --> "<<data.blk->nm<<std::endl;

        blockPopularity[dstLeaf->dom] += static_cast<int>(data.avgIt + 0.5);
        maxSteps -= data.avgIt;
        //DomainBlock::Dump(dstLeaf, std::cout, 10); std::cout<<std::endl;
        auto nextLeaf = data.blk;


        //std::cout<<pct<<" ;; "<<cnt<<": "<<dstLeaf->nm<<" "<<data.avgIt<<" ms: "<<maxSteps<<" --> "<<nextLeaf->nm<<" "<<nextLeaf->dom<<" "<<nextLeaf->gid<<std::endl;

        if (maxSteps <= 0 || nextLeaf->leafBlockType == DomainBlock::INTERNAL)
          break;
        cnt++;
        dstLeaf = nextLeaf;
        std::cout<<"******************************************"<<std::endl<<std::endl;
      }

      std::cout<<std::endl;
      std::cout<<"BlockPopularity: [";
      for (auto& p : blockPopularity)
        std::cout<<p<<" ";
      std::cout<<"]"<<std::endl;
    }
  */

  MPI_Finalize();
  timer.Stop();
  double executionTime = timer.GetElapsedTime();
  if(Rank==0){
    std::cout << "Execution time (exclude data loading) is: " << executionTime << std::endl;
  }
  return (0);
}

// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecXY >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecX >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit velocity >&out

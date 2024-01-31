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

#include <limits>

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
    pt[1] = 0.25;
    pt[2] = 0.25;
    // Doesn't matter what the ID is....
    pts.push_back(vtkm::Particle(pt, i));
  }
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

class FlowEntry
{
public:
  FlowEntry() {}
  FlowEntry(int d, int n, int np, const std::string& _nm)
    : dst(d), totNumSteps(n), totNumP(np), nm(_nm) {}
  FlowEntry(int d, int n, const std::string& _nm)
    : dst(d), totNumSteps(n), totNumP(1), nm(_nm) {}

  void Update(int n)
  {
    this->totNumSteps += n;
    this->totNumP++;
  }

  std::string nm = "--";
  int dst = 0;
  int totNumSteps = 0;
  int totNumP = 0;
};

class FlowStat
{
public:
  FlowStat() {}
  FlowStat(int d, int n, int nTot, int nSteps, const std::string& _nm)
  {
    this->dst = d;
    this->pct = static_cast<float>(n) / static_cast<float>(nTot);
    this->avgSteps = static_cast<float>(nSteps) / static_cast<float>(n);
    this->nm = _nm;
  }

  std::string nm = "--";
  int dst = -1;
  float pct = 0.f;
  float avgSteps = 0.f;

  static inline bool rev_sorter(const FlowStat& s1, const FlowStat& s2)
  {
    return s2.pct < s1.pct;
  }
};

void
BoxOfSeeds(const vtkm::Bounds& bb, int numPts, std::vector<vtkm::Particle>& seeds)
{
  vtkm::FloatDefault x = bb.X.Min, y = bb.Y.Min, z = bb.Z.Min;
  vtkm::FloatDefault dx = bb.X.Length(), dy = bb.Y.Length(), dz = bb.Z.Length();
  for (int i = 0; i < numPts; i++)
  {
    vtkm::Vec3f pt(x + random_1() * dx,
                   y + random_1() * dy,
                   z + random_1() * dz);
    // Doesn't matter what the ID is....
    seeds.push_back(vtkm::Particle(pt, i));
  }
}

void
GenerateFaceSeeds(DomainBlock* block,
                  std::vector<vtkm::Particle>& seeds,
                  int numSeedsPerFace)
{
  int numLeafs = block->NumLeafs();

  for (int i = 0; i < numLeafs; i++)
  {
    auto leaf = block->GetLeafFromIndex(i);
    if (leaf->leafBlockType == DomainBlock::INTERNAL)
    {
      continue;
    }
    else
    {
      GenerateTestPts(leaf, numSeedsPerFace, seeds);
      break;
    }
  }
}


void
GenerateFaceSeeds1(const vtkm::cont::DataSet& ds,
                   std::vector<vtkm::Particle>& seeds,
                   int numSeedsPerFace)
{
  auto bounds = ds.GetCoordinateSystem().GetBounds();
  auto x0 = bounds.X.Min, x1 = bounds.X.Max;
  auto y0 = bounds.Y.Min, y1 = bounds.Y.Max;
  auto z0 = bounds.Z.Min, z1 = bounds.Z.Max;
  double delta = 0.01;

  vtkm::Bounds xbb(x0+delta,x0+delta, y0+delta,y1-delta, z0+delta,z1-delta);
  vtkm::Bounds Xbb(x1-delta,x1-delta, y0+delta,y1-delta, z0+delta,z1-delta);
  vtkm::Bounds ybb(x0+delta,x1-delta, y0+delta,y0+delta, z0+delta,z1-delta);
  vtkm::Bounds Ybb(x0+delta,x1-delta, y1-delta,y1-delta, z0+delta,z1-delta);
  vtkm::Bounds zbb(x0+delta,x1-delta, y0+delta,y1-delta, z0+delta,z0+delta);
  vtkm::Bounds Zbb(x0+delta,x1-delta, y0+delta,y1-delta, z1-delta,z1-delta);

  BoxOfSeeds(xbb, numSeedsPerFace, seeds);
  BoxOfSeeds(Xbb, numSeedsPerFace, seeds);

  BoxOfSeeds(ybb, numSeedsPerFace, seeds);
  BoxOfSeeds(Ybb, numSeedsPerFace, seeds);

  BoxOfSeeds(zbb, numSeedsPerFace, seeds);
  BoxOfSeeds(Zbb, numSeedsPerFace, seeds);
}

void
GenerateInteriorSeeds(const vtkm::cont::DataSet& ds,
                      std::vector<vtkm::Particle>& seeds,
                      int numSeeds,
                      vtkm::FloatDefault delta = 0.0)
{
  auto bounds = ds.GetCoordinateSystem().GetBounds();
  auto x0 = bounds.X.Min, x1 = bounds.X.Max;
  auto y0 = bounds.Y.Min, y1 = bounds.Y.Max;
  auto z0 = bounds.Z.Min, z1 = bounds.Z.Max;

  x0 += delta;
  y0 += delta;
  z0 += delta;

  x1 -= delta;
  y1 -= delta;
  z1 -= delta;

  vtkm::Bounds bbox(x0+delta, x1-delta,
                    y0+delta, y1-delta,
                    z0+delta, z1-delta);

  vtkm::FloatDefault dx = x1-x0, dy = y1-y0, dz = z1-z0;
  BoxOfSeeds(bbox, numSeeds, seeds);

}

void
GenerateFaceSeeds2(const vtkm::cont::DataSet& ds,
                   std::vector<vtkm::Particle>& seeds)
{
  auto coords = ds.GetCoordinateSystem().GetData().template AsArrayHandle<vtkm::cont::ArrayHandleUniformPointCoordinates>();
  auto coordsPortal = coords.ReadPortal();
  auto cellSet = ds.GetCellSet().template AsCellSet<vtkm::cont::CellSetStructured<3>>();
  auto cellDims = cellSet.GetCellDimensions();

  // Cell center for every boundary face.
  //Dumb way to do this... Use a WorkletPointNeighborhood worklet....
  vtkm::Id idx = 0, pid = 0;
  for (vtkm::Id i = 0; i < cellDims[0]; i++)
  {
    for (vtkm::Id j = 0; j < cellDims[1]; j++)
    {
      for (vtkm::Id k = 0; k < cellDims[2]; k++)
      {
        if (i == 0 || i == cellDims[0]-1 ||
            j == 0 || j == cellDims[1]-1 ||
            k == 0 || k == cellDims[2]-1)
        {
          //pick every 'nth' cell.
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

  //Save out for validation.
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

vtkm::worklet::flow::ParticleAdvectionResult<vtkm::Particle>
AdvectFaceSeeds(const vtkm::cont::DataSet& ds,
                const std::string& fieldNm,
                vtkm::FloatDefault stepSize,
                vtkm::Id maxSteps,
                std::vector<vtkm::Particle>& seeds)
{
  //Advect all these points and see where they go.
  auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
  using FieldType = vtkm::worklet::flow::VelocityField<ArrayType>;
  using GridEvalType = vtkm::worklet::flow::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::flow::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::flow::Stepper<RK4Type, GridEvalType>;

  ArrayType arr;
  vtkm::cont::ArrayCopyShallowIfPossible(ds.GetField(fieldNm).GetData(), arr);
  auto arrP = arr.WritePortal();

  /*
  vtkm::Vec3f vecX(1,0,0);
  for (vtkm::Id i = 0; i < arr.GetNumberOfValues(); i++)
  {
    arrP.Set(i, vecX);
  }
  if (Rank == 0) std::cout<<"**********************************  using VecX!!!"<<std::endl;
  */


  FieldType velocities(arr, ds.GetField(fieldNm).GetAssociation());
  GridEvalType gridEval(ds, velocities);
  Stepper rk4(gridEval, stepSize);

  vtkm::worklet::flow::ParticleAdvection pa;
  auto res = pa.Run(rk4, seedArray, maxSteps);

  //Now, figure out where each point goes.
  if (res.Particles.ReadPortal().GetNumberOfValues() != seeds.size())
    throw "Array sizes do not match!";

  return res;
}

template <typename T>
void
BuildFlowMap(const T& endPtsPortal,
             const std::vector<vtkm::Particle>& seeds,
             const vtkm::Id blockID,
             const std::vector<DomainBlock *>& blockInfo,
             const vtkm::filter::flow::internal::BoundsMap& boundsMap,
             std::map<int, std::vector<FlowEntry>>& flowMap,
             std::map<int,int>& countFromSource)

{
  int totNumLeafs = DomainBlock::TotalNumLeaves(blockInfo);
  int numLeafs = blockInfo[blockID]->NumLeafs();
  vtkm::Id numPts = endPtsPortal.GetNumberOfValues();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    auto const p = endPtsPortal.Get(i);
    auto p0 = seeds[i].GetPosition();
    auto p1 = p.GetPosition();
    auto numSteps = p.GetNumberOfSteps();
    auto status = p.GetStatus();
    // p0 goes to p1 in numSteps.

    DomainBlock *srcLeaf = blockInfo[blockID]->GetLeaf(p0);
    int dst = -1;

    if (status.CheckTerminate())
    {
      dst = -1;
    }
    else if (status.CheckSpatialBounds())
    {
      auto destinations = boundsMap.FindBlocks(p1, blockID);
      if (destinations.size() > 1)
        throw "Should be a single destination.";

      if (destinations.empty())
        dst = -1;
      else
      {
        auto dstBlock = blockInfo[destinations[0]];
        auto dstLeaf = blockInfo[destinations[0]]->GetLeaf(p1);
        dst = dstLeaf->gid;
        //Randomize the destination.
        //dst = random_1() * totNumLeafs;
        //dstLeaf = DomainBlock::GetBlockFromGID(blockInfo, dst);
      }
    }
    else
    {
      throw "Particle has unknown status.";
    }

    std::string dstNm = "TERM";
    if (dst != -1) dstNm = DomainBlock::GetBlockFromGID(blockInfo, dst)->nm;


    if (Rank == 0)
    {
      /*
      if (i == 0)
      {
        std::vector<DomainBlock*> leaves;
        blockInfo[1]->GetLeaves(leaves);
        std::cout<<"Domain 1 leaves"<<std::endl;
        vtkm::FloatDefault bb[6];
        for (int ii = 0; ii < leaves.size(); ii++)
        {
          leaves[ii]->GetBBox(bb);
          std::cout<<"  "<<leaves[ii]->nm<<": ("<<bb[0]<<" "<<bb[1]<<") ("<<bb[2]<<" "<<bb[3]<<") ("<<bb[4]<<" "<<bb[5]<<")"<<std::endl;
          std::cout<<"    ptin= "<<leaves[ii]->InBBox(p1)<<std::endl;
        }
      }
      */

      //std::cout<<i<<": "<<p0<<" "<<srcLeaf->gid<<" "<<srcLeaf->nm<<" -- "<<numSteps<<" --> "<<p1<<" "<<dst<<" "<<dstNm<<" "<<p.GetStatus()<<std::endl;
    }

    //We have src/dst. Add it to the map.
    if (countFromSource.find(srcLeaf->gid) == countFromSource.end())
      countFromSource[srcLeaf->gid] = 0;
    countFromSource[srcLeaf->gid]++;

    auto it = flowMap.find(srcLeaf->gid);
    if (it == flowMap.end()) // new src. Add the dst.
      flowMap[srcLeaf->gid].push_back(FlowEntry(dst, numSteps, dstNm));
    else
    {
      //find the dst
      bool found = false;
      for (auto& entry : it->second)
      {
        if (entry.dst == dst) //update a found dst
        {
          entry.Update(numSteps);
          found = true;
          break;
        }
      }
      if (!found) // not found, add a new dst.
        it->second.push_back({dst, numSteps, dstNm});
    }
  }

  if (Rank == PRINT_RANK)
  {
    std::cout<<"*******************************  Local flow map"<<std::endl;
    for (auto it : flowMap)
    {
      std::cout<<"Src: "<<it.first<<" totPts= "<<countFromSource[it.first]<<"  (dst, totNumSteps, TotNumPts)"<<std::endl;
      for (auto entry : it.second)
      {
        std::cout<<"  "<<entry.dst<<" "<<entry.totNumSteps<<" "<<entry.totNumP<<std::endl;
      }
    }
    std::cout<<"*******************************"<<std::endl;
  }
}

std::map<int, std::vector<FlowStat>>
ComputeGlobalFlowMap(vtkm::Id Nxyz,
                     const std::map<int, std::vector<FlowEntry>>& flowMap,
                     vtkm::Id blockID,
                     const std::vector<DomainBlock *>& blockInfo,
                     const std::map<int,int>& countFromSource)
{
  //determine largest size of an entry for a flow map.
  //maxSz[0]: max number of keys in flowMap
  //maxSz[1]: max number of array size in values for each key.

  int maxSz[2] = {(int)flowMap.size(), -1};
  for (auto it = flowMap.begin(); it != flowMap.end(); it++)
    maxSz[1] = std::max(maxSz[1], (int)it->second.size());

  MPI_Allreduce(MPI_IN_PLACE, &maxSz, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int totNumBlocks = Size;
  int nVals = 5;  //(src, dst, numPts, numIters, totalPtsFromSrc)

  int szPerBlock = maxSz[0] * maxSz[1] * nVals;

  int EMPTY_VALUE = std::numeric_limits<int>::min();

  std::vector<int> flowMapData(totNumBlocks * szPerBlock, EMPTY_VALUE);

  int offset = Rank * szPerBlock;
  for (auto it = flowMap.begin(); it != flowMap.end(); it++)
  {
    int src = it->first;
    int n = it->second.size();
    for (int i = 0; i < n; i++)
    {
      flowMapData[offset + 0] = src;
      flowMapData[offset + 1] = it->second[i].dst;
      flowMapData[offset + 2] = it->second[i].totNumP;
      flowMapData[offset + 3] = it->second[i].totNumSteps;
      flowMapData[offset + 4] = countFromSource.find(src)->second;
      offset += nVals;
    }
  }

  //Now everyone has the array...
  MPI_Allreduce(MPI_IN_PLACE, flowMapData.data(), flowMapData.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  /*
  if (Rank == 0)
  {
    for (int i = 0; i < flowMapData.size(); i+= nVals)
    {
      std::cout<<i<<": ";
      for (int j = 0; j < nVals; j++) std::cout<<flowMapData[i+j]<<" ";
      std::cout<<std::endl;
    }
  }
  */

  //Build the global flow map
  std::map<int, std::vector<FlowStat>> allFlowMap;
  for (int i = 0; i < flowMapData.size(); i += nVals)
  {
    int src = flowMapData[i];
    if (src == EMPTY_VALUE)
      continue;

    int dst = flowMapData[i+1];
    int numPts = flowMapData[i+2];
    int totNumSteps = flowMapData[i+3];
    int totPtsFromSrc = flowMapData[i+4];

    std::string dstNm = "TERM";
    if (dst != -1) dstNm = DomainBlock::GetBlockFromGID(blockInfo, dst)->nm;

    FlowStat fs(dst, numPts, totPtsFromSrc, totNumSteps, dstNm);
    auto it = allFlowMap.find(src);
    if (it == allFlowMap.end())
      allFlowMap[src].clear();
    allFlowMap[src].push_back(fs);
  }

  //sort the entries.  should be sorted before communication...
  for (auto it = allFlowMap.begin(); it != allFlowMap.end(); it++)
  {
    if (it->second.size() > 1)
      std::sort(it->second.begin(), it->second.end(), FlowStat::rev_sorter);
  }

  if (Rank == PRINT_RANK)
  {
    std::cout<<"********* Global Flow Map "<<std::endl;
    for (auto it = allFlowMap.begin(); it != allFlowMap.end(); it++)
    {
      auto leaf = DomainBlock::GetBlockFromGID(blockInfo, it->first);
      std::cout<<"Src: "<<it->first<<" "<<leaf->nm<<std::endl;
      for (int i = 0; i < it->second.size(); i++)
      {
        const auto fs = it->second[i];
        std::cout<<"   "<<fs.dst<<" "<<fs.nm<<" : "<<fs.pct<<" "<<fs.avgSteps<<std::endl;
      }
      std::cout<<std::endl;
    }
    std::cout<<"*********"<<std::endl;
  }

  return allFlowMap;
}


std::map<int, std::vector<FlowStat>>
CreateFlowMap(vtkm::Id Nxyz,
              const std::vector<DomainBlock *>& blockInfo,
              const vtkm::cont::DataSet& ds,
              const vtkm::Id blockID,
              const std::string& fieldNm,
              vtkm::FloatDefault stepSize,
              vtkm::Id maxSteps,
              int numFacePts,
              int numInteriorPts,
              const vtkm::filter::flow::internal::BoundsMap& boundsMap)
{
  // Go through all subdomains and compute the number of leaves
  std::vector<vtkm::Particle> seeds;
  //Each face is broken up into Nxyz x Nxyz pieces.
  GenerateFaceSeeds1(ds, seeds, numFacePts * Nxyz * Nxyz);
  //GenerateInteriorSeeds(ds, seeds, numInteriorPts, 0.01);

  auto res = AdvectFaceSeeds(ds, fieldNm, stepSize, maxSteps, seeds);

  //create the flow map for the test seeds.

  // {src : {dst, totNumSteps, totNumP}}
  std::map<int, std::vector<FlowEntry>> flowMap;
  int numLeafs = blockInfo[blockID]->NumLeafs();
  std::map<int,int> countFromSource;
  BuildFlowMap(res.Particles.ReadPortal(), seeds, blockID, blockInfo, boundsMap, flowMap, countFromSource);

  auto globalFlowMap = ComputeGlobalFlowMap(Nxyz, flowMap, blockID, blockInfo, countFromSource);
  return globalFlowMap;
}

FlowStat
PickRandomWeightedDst(const std::vector<FlowStat>& entries)
{
  if (entries.empty())
    throw "Empty FlowStat";

  if (entries.size() == 1)
    return entries[0];

  vtkm::FloatDefault pct = random_1();
  for (const auto& e : entries)
  {
    if (pct < e.pct)
      return e;
    pct -= e.pct;
  }

  throw "Error! pct not found";
  return entries[entries.size()-1];

}

void
ReduceToRoot(std::vector<int>& data)
{
  std::vector<int> tmp;
  if (Rank == 0)
    tmp.resize(data.size(), 0);

  MPI_Reduce(data.data(), tmp.data(), data.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Rank == 0)
    data = tmp;
}

void CalcBlockPopularity(std::vector<DomainBlock *> blockInfo,
                         const vtkm::cont::DataSet& ds,
                         int blockID,
                         std::map<int, std::vector<FlowStat>>& flowMap,
                         const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                         std::vector<int> &blockPopularity,
                         std::vector<int> &particlesIn,
                         std::vector<int> &particlesOut,
                         std::vector<int> &cycleCnt,
                         int numPts,
                         const std::string& fieldNm,
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

    //seed terminated inside the block.
    if (seed.GetStatus().CheckTerminate() || seed.GetNumberOfSteps() >= maxSteps || bid == blockID)
      continue;

    vtkm::FloatDefault numSteps = static_cast<vtkm::FloatDefault>(seed.GetNumberOfSteps());
    domPath.push_back(bid);
    int gid = blockInfo[bid]->GetLeaf(seed.GetPosition())->gid;

    while (true)
    {
      // Determine the destination based on a random percentage.
      auto it = flowMap.find(gid);
      if (it == flowMap.end())
      {
        throw "Source not found! " + std::to_string(gid);
      }
      FlowStat dstEntry = PickRandomWeightedDst(it->second);

      numSteps += dstEntry.avgSteps;
      blockPopularity[bid] += static_cast<int>(dstEntry.avgSteps + .5);
      particlesOut[bid]++;

      //Terminate or take max number of steps, we are done.
      int nextGID = dstEntry.dst;
      if (nextGID == -1 || numSteps >= maxStepsFloat)
        break;

      auto nextLeaf = DomainBlock::GetBlockFromGID(blockInfo, nextGID);
      particlesIn[nextLeaf->dom]++;
      bid = nextLeaf->dom;
      gid = nextGID;
    }

    //Detect any cycles.
    DetectCycles(domPath, 2, cycleCnt);
    DetectCycles(domPath, 3, cycleCnt);
    DetectCycles(domPath, 4, cycleCnt);
  }

  // Calculate the block popularity over all ranks.
  ReduceToRoot(blockPopularity);
  ReduceToRoot(particlesIn);
  ReduceToRoot(particlesOut);
  ReduceToRoot(cycleCnt);
}


int main(int argc, char **argv)
{
  if (Rank == 0)
  {
    if (argc != 10)
    {
      std::cout << "<executable> <visitfileName> <fieldNm> <stepSize> <maxSteps> <numFacePts> <numInteriorPts> <numTestPts> <Nxyz> <pctWidth>" << std::endl;
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
  vtkm::Id numFacePts = std::atoi(argv[5]);
  vtkm::Id numInteriorPts = std::atoi(argv[6]);
  vtkm::Id numTestPts = std::atoi(argv[7]);
  vtkm::Id Nxyz = std::atoi(argv[8]);
  vtkm::FloatDefault pctWidth = std::stod(argv[9]);

  if (pctWidth <= 0 || pctWidth >= 0.5)
  {
    std::cout << "pctWidth shoule be a value between 0 to 0.5" << std::endl;
    exit(0);
  }

  if (Rank == 0)
  {
    std::cout<<"Checking input parameters:" <<std::endl;
    std::cout<<" visitfileName: "<<visitfileName<<std::endl;
    std::cout<<" fieldNm: "<<fieldNm<<std::endl;
    std::cout<<" stepSize: "<<stepSize<<std::endl;
    std::cout<<" maxSteps: "<<maxSteps<<std::endl;
    std::cout<<" numFacePts: "<<numFacePts<<std::endl;
    std::cout<<" numInteriorPts: "<<numInteriorPts<<std::endl;
    std::cout<<" numTestPts: "<<numTestPts<<std::endl;
    std::cout<<" Nxyz: "<<Nxyz<<std::endl;
    std::cout<<" pctWidth: "<<pctWidth<<std::endl;
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
    throw "Num blocks must be the same as number of ranks";

  vtkm::filter::flow::internal::BoundsMap boundsMap(dataSets);
  std::vector<DomainBlock *> blockInfo;
  int NX = Nxyz, NY = Nxyz, NZ = Nxyz;
  bool subdivUniform = false;
  DomainBlock::CreateBlockInfo(blockInfo, totNumBlocks, boundsMap, subdivUniform, NX, NY, NZ, pctWidth);


  std::vector<int> blockPopularity(Size, 0), particlesIn(Size, 0), particlesOut(Size, 0), cycleCnt(Size, 0);

  int blockID = boundsMap.GetLocalBlockId(0);
  const auto &ds = dataSets[0];
  auto block = blockInfo[blockID];
  auto flowMap = CreateFlowMap(Nxyz, blockInfo, ds, blockID, fieldNm, stepSize, maxSteps, numFacePts, numInteriorPts, boundsMap);

  CalcBlockPopularity(blockInfo, ds, blockID, flowMap, boundsMap, blockPopularity, particlesIn, particlesOut, cycleCnt, numTestPts, fieldNm, stepSize, maxSteps);

  MPI_Barrier(MPI_COMM_WORLD);
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

  MPI_Finalize();
}

// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecXY >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit vecX >&out
// mpirun -np 8 ./StreamlineMPI ./data/clover.visit velocity >&out

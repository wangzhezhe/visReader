#include "SLTool.h"
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSet.h>
#include <mpi.h>


int SLTOOL_Rank=0;
int SLTOOL_TOTALPROCS=0;
int PRINT_RANK=-1;
int PRINT_DETAILS=-1;

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

FlowStat
PickRandomWeightedDst(const std::vector<FlowStat> &entries)
{
  if (entries.empty())
    throw std::runtime_error("Empty FlowStat");

  if (entries.size() == 1)
    return entries[0];

  vtkm::FloatDefault pct = random_1();
  for (const auto &e : entries)
  {
    if (pct < e.pct)
      return e;
    pct -= e.pct;
  }

  throw std::runtime_error("Error! pct not found");
  return entries[entries.size() - 1];
}


void BoxOfSeeds(const vtkm::Bounds &bb, int numPts, std::vector<vtkm::Particle> &seeds)
{
  vtkm::FloatDefault x = bb.X.Min, y = bb.Y.Min, z = bb.Z.Min;
  vtkm::FloatDefault dx = bb.X.Length(), dy = bb.Y.Length(), dz = bb.Z.Length();
  for (int i = 0; i < numPts; i++)
  {
    vtkm::Vec3f pt(x + random_1() * dx,
                   y + random_1() * dy,
                   z + random_1() * dz);
    // if (SEED_PING_PONG)
    // {
    //   pt[0] = .32;
    //   pt[1] = .495;
    //   pt[2] = .47;
    // }

    // Doesn't matter what the ID is....
    seeds.push_back(vtkm::Particle(pt, i));
  }
}

void GenerateFaceSeeds1(const vtkm::cont::DataSet &ds,
                        std::vector<vtkm::Particle> &seeds,
                        int numSeedsPerFace)
{
  auto bounds = ds.GetCoordinateSystem().GetBounds();
  auto x0 = bounds.X.Min, x1 = bounds.X.Max;
  auto y0 = bounds.Y.Min, y1 = bounds.Y.Max;
  auto z0 = bounds.Z.Min, z1 = bounds.Z.Max;
  double delta = 0.01;

  vtkm::Bounds xbb(x0 + delta, x0 + delta, y0 + delta, y1 - delta, z0 + delta, z1 - delta);
  vtkm::Bounds Xbb(x1 - delta, x1 - delta, y0 + delta, y1 - delta, z0 + delta, z1 - delta);
  vtkm::Bounds ybb(x0 + delta, x1 - delta, y0 + delta, y0 + delta, z0 + delta, z1 - delta);
  vtkm::Bounds Ybb(x0 + delta, x1 - delta, y1 - delta, y1 - delta, z0 + delta, z1 - delta);
  vtkm::Bounds zbb(x0 + delta, x1 - delta, y0 + delta, y1 - delta, z0 + delta, z0 + delta);
  vtkm::Bounds Zbb(x0 + delta, x1 - delta, y0 + delta, y1 - delta, z1 - delta, z1 - delta);

  BoxOfSeeds(xbb, numSeedsPerFace, seeds);
  BoxOfSeeds(Xbb, numSeedsPerFace, seeds);

  BoxOfSeeds(ybb, numSeedsPerFace, seeds);
  // if (Rank == 5) SEED_PING_PONG = 1;
  BoxOfSeeds(Ybb, numSeedsPerFace, seeds);
  // if (Rank == 5) SEED_PING_PONG = 0;

  BoxOfSeeds(zbb, numSeedsPerFace, seeds);
  BoxOfSeeds(Zbb, numSeedsPerFace, seeds);
}

template <typename T>
void BuildFlowMap(const T &endPtsPortal,
                  const std::vector<vtkm::Particle> &seeds,
                  const vtkm::Id blockID,
                  const std::vector<DomainBlock *> &blockInfo,
                  const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                  std::map<int, std::vector<FlowEntry>> &flowMap,
                  std::map<int, int> &countFromSource)

{
  int totNumLeafs = DomainBlock::TotalNumLeaves(blockInfo);
  int numLeafs = blockInfo[blockID]->NumLeafs();
  vtkm::Id numPts = endPtsPortal.GetNumberOfValues();

  if (SLTOOL_Rank == PRINT_DETAILS)
    std::cout << "************************  BuildFlowMap" << std::endl;

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

    if (SLTOOL_Rank == PRINT_DETAILS)
    {
      vtkm::Vec3f pingPongPt(.32, .495, .47);
      auto pingLeaf = blockInfo[blockID]->GetLeaf(pingPongPt);
      if (srcLeaf == nullptr)
      {
        srcLeaf = blockInfo[blockID]->GetLeaf(p0);
      }
      std::cout << i << ": " << p0 << " " << srcLeaf->nm << " " << srcLeaf->gid << " ---> " << p1 << " ";
    }

    if (status.CheckTerminate())
    {
      dst = -1;
      if (SLTOOL_Rank == PRINT_DETAILS)
        std::cout << "TERMINATE" << std::endl;
    }
    else if (status.CheckSpatialBounds())
    {
      auto destinations = boundsMap.FindBlocks(p1, blockID);
      if (destinations.size() > 1)
      {
        std::cout << "warning, multiple destination.";
        for (int k = 0; k < destinations.size(); k++)
        {
          std::cout << destinations[i] << ",";
        }
        std::cout << std::endl;
      }

      if (destinations.empty())
      {
        dst = -1;
        if (SLTOOL_Rank == PRINT_DETAILS)
          std::cout << "TERMINATE" << std::endl;
      }
      else
      {
        auto dstBlock = blockInfo[destinations[0]];
        auto dstLeaf = dstBlock->GetLeaf(p1);
        dst = dstLeaf->gid;
        if (SLTOOL_Rank == PRINT_DETAILS)
          std::cout << " DST= " << dstLeaf->nm << std::endl;
        // Randomize the destination.
        // dst = random_1() * totNumLeafs;
        // dstLeaf = DomainBlock::GetBlockFromGID(blockInfo, dst);
      }
    }
    else
    {
      throw std::runtime_error("Particle has unknown status.");
    }

    std::string dstNm = "TERM";
    if (dst != -1)
      dstNm = DomainBlock::GetBlockFromGID(blockInfo, dst)->nm;

    if (PRINT_RANK == 0)
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

      // std::cout<<i<<": "<<p0<<" "<<srcLeaf->gid<<" "<<srcLeaf->nm<<" -- "<<numSteps<<" --> "<<p1<<" "<<dst<<" "<<dstNm<<" "<<p.GetStatus()<<std::endl;
    }

    // We have src/dst. Add it to the map.
    if (countFromSource.find(srcLeaf->gid) == countFromSource.end())
      countFromSource[srcLeaf->gid] = 0;
    countFromSource[srcLeaf->gid]++;

    auto it = flowMap.find(srcLeaf->gid);
    if (it == flowMap.end()) // new src. Add the dst. dstName is only used for debugging
      flowMap[srcLeaf->gid].push_back(FlowEntry(dst, numSteps, dstNm));
    else
    {
      // find the dst
      bool found = false;
      for (auto &entry : it->second)
      {
        if (entry.dst == dst) // update a found dst
        {
          // update both numsteps and number of particles
          entry.Update(numSteps);
          found = true;
          break;
        }
      }
      if (!found) // not found, add a new dst.
        it->second.push_back({dst, numSteps, dstNm});
    }
  }

  if (SLTOOL_Rank == PRINT_RANK)
  {
    std::cout << "*******************************  Local flow map" << std::endl;
    for (auto it : flowMap)
    {
      std::cout << "Src: " << it.first << " totPts= " << countFromSource[it.first] << "  (dst, totNumSteps, TotNumPts)" << std::endl;
      for (auto entry : it.second)
      {
        std::cout << "  " << entry.dst << " " << entry.totNumSteps << " " << entry.totNumP << std::endl;
      }
    }
    std::cout << "*******************************" << std::endl;
  }
}


std::map<int, std::vector<FlowStat>>
ComputeGlobalFlowMap(vtkm::Id Nxyz,
                     const std::map<int, std::vector<FlowEntry>> &flowMap,
                     vtkm::Id blockID,
                     const std::vector<DomainBlock *> &blockInfo,
                     const std::map<int, int> &countFromSource)
{
  // determine largest size of an entry for a flow map.
  // maxSz[0]: max number of keys in flowMap
  // maxSz[1]: max number of array size in values for each key.

  int maxSz[2] = {(int)flowMap.size(), -1};
  for (auto it = flowMap.begin(); it != flowMap.end(); it++)
    maxSz[1] = std::max(maxSz[1], (int)it->second.size());

  MPI_Allreduce(MPI_IN_PLACE, &maxSz, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int totNumBlocks = SLTOOL_TOTALPROCS;
  int nVals = 5; //(src, dst, numPts, numIters, totalPtsFromSrc)

  int szPerBlock = maxSz[0] * maxSz[1] * nVals;

  int EMPTY_VALUE = std::numeric_limits<int>::min();

  std::vector<int> flowMapData(totNumBlocks * szPerBlock, EMPTY_VALUE);

  int offset = SLTOOL_Rank * szPerBlock;
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

  // Now everyone has the array...
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

  // Build the global flow map
  std::map<int, std::vector<FlowStat>> allFlowMap;
  for (int i = 0; i < flowMapData.size(); i += nVals)
  {
    int src = flowMapData[i];
    if (src == EMPTY_VALUE)
      continue;

    int dst = flowMapData[i + 1];
    int numPts = flowMapData[i + 2];
    int totNumSteps = flowMapData[i + 3];
    int totPtsFromSrc = flowMapData[i + 4];

    std::string dstNm = "TERM";
    if (dst != -1)
      dstNm = DomainBlock::GetBlockFromGID(blockInfo, dst)->nm;

    FlowStat fs(dst, numPts, totPtsFromSrc, totNumSteps, dstNm);
    auto it = allFlowMap.find(src);
    if (it == allFlowMap.end())
      allFlowMap[src].clear();
    allFlowMap[src].push_back(fs);
  }

  // sort the entries.  should be sorted before communication...
  for (auto it = allFlowMap.begin(); it != allFlowMap.end(); it++)
  {
    if (it->second.size() > 1)
      std::sort(it->second.begin(), it->second.end(), FlowStat::rev_sorter);
  }

  if (SLTOOL_Rank == PRINT_RANK)
  {
    std::cout << "********* Global Flow Map " << std::endl;
    for (auto it = allFlowMap.begin(); it != allFlowMap.end(); it++)
    {
      auto leaf = DomainBlock::GetBlockFromGID(blockInfo, it->first);
      std::cout << "Src: " << it->first << " " << leaf->nm << std::endl;
      for (int i = 0; i < it->second.size(); i++)
      {
        const auto fs = it->second[i];
        std::cout << "   " << fs.dst << " " << fs.nm << " : " << fs.pct << " " << fs.avgSteps << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << "*********" << std::endl;
  }

  return allFlowMap;
}


std::map<int, std::vector<FlowStat>> CreateFlowMap(vtkm::Id Nxyz,
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
  if (SLTOOL_Rank == 0)
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
  if (SLTOOL_Rank == 0)
  {
   std::cout << "Execution time to BuildFlowMap in CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  auto globalFlowMap = ComputeGlobalFlowMap(Nxyz, flowMap, blockID, blockInfo, countFromSource);
  //TODO, there is runtime error at this step
  timer.Stop();
  executionTime = timer.GetElapsedTime();
  if (SLTOOL_Rank == 0)
  {
   std::cout << "Execution time to ComputeGlobalFlowMap in CreateFlowMap is: " << executionTime << std::endl;
  }
  timer.Start();

  return globalFlowMap;
}


vtkm::worklet::flow::ParticleAdvectionResult<vtkm::Particle>
AdvectFaceSeeds(const vtkm::cont::DataSet &ds,
                const std::string &fieldNm,
                vtkm::FloatDefault stepSize,
                vtkm::Id maxSteps,
                std::vector<vtkm::Particle> &seeds)
{
  // Advect all these points and see where they go.
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

  // After this step, we could figure out where each point goes.
  if (res.Particles.ReadPortal().GetNumberOfValues() != seeds.size())
    throw std::runtime_error("Array sizes do not match!");

  return res;
}

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

void GenerateFaceSeeds(DomainBlock *block,
                       std::vector<vtkm::Particle> &seeds,
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

void GenerateInteriorSeeds(const vtkm::cont::DataSet &ds,
                           std::vector<vtkm::Particle> &seeds,
                           int numSeeds,
                           vtkm::FloatDefault delta)
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

  vtkm::Bounds bbox(x0 + delta, x1 - delta,
                    y0 + delta, y1 - delta,
                    z0 + delta, z1 - delta);

  vtkm::FloatDefault dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
  BoxOfSeeds(bbox, numSeeds, seeds);
}

void ReduceToRoot(std::vector<int> &data)
{
  std::vector<int> tmp;
  if (SLTOOL_Rank == 0)
    tmp.resize(data.size(), 0);

  MPI_Reduce(data.data(), tmp.data(), data.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (SLTOOL_Rank == 0)
    data = tmp;
}

void CalcBlockPopularityMultiStages(std::vector<DomainBlock *> blockInfo,
                                    const vtkm::cont::DataSet &ds,
                                    int blockID,
                                    std::map<int, std::vector<FlowStat>> &flowMap,
                                    const vtkm::filter::flow::internal::BoundsMap &boundsMap,
                                    std::vector<int> &BlockStepsMultiStages,
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

  // for each block, there should be a list
  // size of stepRecords should be numBlocks
  std::vector<std::vector<std::pair<int, int>>> stepRecords(particlesIn.size());

  int localMaxCycle = 0;
  for (vtkm::Id i = 0; i < n; i++)
  {
    int cycle = 0;
    std::vector<int> domPath;

    vtkm::Particle seed = portal.Get(i);

    stepRecords[blockID].push_back({cycle, seed.GetNumberOfSteps()});

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

      cycle = cycle + 1;

      numSteps += dstEntry.avgSteps;

      particlesOut[bid]++;
      stepRecords[bid].push_back({cycle, dstEntry.avgSteps});
      localMaxCycle = std::max(localMaxCycle, cycle);

      // Terminate or take max number of steps, we are done.
      int nextGID = dstEntry.dst;
      auto nextLeaf = DomainBlock::GetBlockFromGID(blockInfo, nextGID);

      if (nextGID == -1 || numSteps >= maxStepsFloat)
        break;

      particlesIn[nextLeaf->dom]++;
      bid = nextLeaf->dom;
      gid = nextGID;
      domPath.push_back(bid);
      cycle += 1;
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

  // Calculate the block popularity over all ranks.
  // check the maximal steps
  // according to max cycle and number of stages
  // split the current results into two parts
  // do a mpi reduce to get globalMaxCycle
  int globalMaxCycle = 0;
  MPI_Allreduce(&localMaxCycle, &globalMaxCycle, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // if (Rank==0){
  // print out the block id results

  // }
  // std::cout << "rank " << Rank << " max cycle is " << localMaxCycle << " globalMaxCycle is " << globalMaxCycle << std::endl;

  // if there are multiple stages use another block steps datastructure
  // the index of the vector is block id, and the element is staged block workload for each block
  int numStages = 3;
  int keyCycle0 = 0;
  int keyCycle1 = globalMaxCycle / 3;
  int keyCycle2 = (globalMaxCycle / 3) * 2;
  int keyCycle3 = globalMaxCycle;

  //std::cout << "keycycle " << keyCycle1 << "," << keyCycle2 << "," << keyCycle3 << std::endl;

  for (int i = 0; i < stepRecords.size(); i++)
  {
    // for each recorded block
    int bid = i;
    // for each advected recrod
    for (int j = 0; j < stepRecords[i].size(); j++)
    {
      // put results into the BlockStepsMultiStages
      // put into BlockStepsMultiStages 0 1 2 ?
      int advCycle = stepRecords[i][j].first;
      int advSteps = stepRecords[i][j].second;

      if (advCycle >= keyCycle0 && advCycle < keyCycle1)
      {
        BlockStepsMultiStages[bid * numStages + 0] += advSteps;
      }
      else if (advCycle >= keyCycle1 && advCycle < keyCycle2)
      {
        BlockStepsMultiStages[bid * numStages + 1] += advSteps;
      }
      else
      {
        BlockStepsMultiStages[bid * numStages + 2] += advSteps;
      }
    }
  }

  //std::cout << " Rank " << Rank << " blockAdvMultiStages: " << BlockStepsMultiStages << std::endl;
  ReduceToRoot(BlockStepsMultiStages);
}

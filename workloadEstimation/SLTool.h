#ifndef SLTOOL_H
#define SLTOOL_H

#include "Block.h"

#include <vtkm/cont/DataSet.h>
#include <vtkm/Particle.h>

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>

#include <vtkm/filter/flow/ParticleAdvection.h>
#include <vtkm/filter/flow/worklet/RK4Integrator.h>
#include <vtkm/filter/flow/worklet/Stepper.h>
#include <vtkm/filter/flow/worklet/ParticleAdvection.h>
#include <vtkm/filter/flow/worklet/Particles.h>

#include <map>
#include <vector>
#include <string>

class FlowEntry
{
public:
  FlowEntry() {}
  FlowEntry(int d, int n, int np, const std::string &_nm)
      : dst(d), totNumSteps(n), totNumP(np), nm(_nm) {}
  FlowEntry(int d, int n, const std::string &_nm)
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
  FlowStat(int d, int n, int nTot, int nSteps, const std::string &_nm)
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

  static inline bool rev_sorter(const FlowStat &s1, const FlowStat &s2)
  {
    return s2.pct < s1.pct;
  }
};


std::map<int, std::vector<FlowStat>> CreateFlowMap(vtkm::Id Nxyz,
              const std::vector<DomainBlock *> &blockInfo,
              const vtkm::cont::DataSet &ds,
              const vtkm::Id blockID,
              const std::string &fieldNm,
              vtkm::FloatDefault stepSize,
              vtkm::Id maxSteps,
              int numFacePts,
              const vtkm::filter::flow::internal::BoundsMap &boundsMap);


vtkm::worklet::flow::ParticleAdvectionResult<vtkm::Particle>
AdvectFaceSeeds(const vtkm::cont::DataSet &ds,
                const std::string &fieldNm,
                vtkm::FloatDefault stepSize,
                vtkm::Id maxSteps,
                std::vector<vtkm::Particle> &seeds);

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
                                    int maxSteps);


FlowStat PickRandomWeightedDst(const std::vector<FlowStat> &entries);


void ReduceToRoot(std::vector<int> &data);

void GenerateInteriorSeeds(const vtkm::cont::DataSet &ds,
                           std::vector<vtkm::Particle> &seeds,
                           int numSeeds,
                           vtkm::FloatDefault delta = 0.0);

extern int SLTOOL_Rank;
extern int SLTOOL_TOTALPROCS;
extern int PRINT_RANK;
extern int PRINT_DETAILS;

#endif
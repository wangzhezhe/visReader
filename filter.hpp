#ifndef __FILTER_H
#define __FILTER_H

//--vtkm
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/cont/DeviceAdapterTag.h>
namespace FILTER
{
    extern vtkm::FloatDefault G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax;
    extern vtkm::FloatDefault GLOBAL_ADVECT_STEP_SIZE;
    extern int G_SampleX,G_SampleY,G_SampleZ;
    extern int GLOBAL_ADVECT_NUM_STEPS;
    extern int GLOBAL_ADVECT_NUM_SEEDS;
    extern int GLOBAL_NUM_LEVELS;
    extern std::ofstream *timingInfo;
    extern std::string CommStrategy;
    extern int GLOBAL_NUM_RECIEVERS;
    extern int GLOBAL_NUM_PARTICLE_PER_PACKET;
    extern bool GLOBAL_BLOCK_DUPLICATE;

    void runAdvection(const vtkm::cont::PartitionedDataSet &pds, int rank, int numRanks, int step,
                      std::string seedMethod, std::string fieldToOperateOn, bool cloverleaf,
                      bool recordTrajectories, bool outputResults, bool outputseeds, vtkm::cont::DeviceAdapterId& deviceID);
}

#endif
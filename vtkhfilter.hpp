#ifndef __VTKH_FILTER_H
#define __VTKH_FILTER_H

//--VTKh
#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/MarchingCubes.hpp>
#include <vtkh/filters/PointAverage.hpp>
#include <vtkh/rendering/RayTracer.hpp>
#include <vtkh/rendering/VolumeRenderer.hpp>
#include <vtkh/rendering/MeshRenderer.hpp>
#include <vtkh/rendering/Scene.hpp>
#include <vtkh/Logger.hpp>
#include <vtkh/filters/GhostStripper.hpp>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/filters/Streamline.hpp>
#include <vtkh/filters/Recenter.hpp>

//--vtkm
#include <vtkm/filter/CleanGrid.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/ArrayHandleCartesianProduct.h>
#include <vtkm/cont/ArrayHandleUniformPointCoordinates.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/exec/ParametricCoordinates.h>
#include <vtkm/exec/CellInterpolate.h>
#include <vtkm/filter/particleadvection/BoundsMap.h>
//--statistics
#include <vtkh/StatisticsDB.hpp>

namespace VTKH_FILTER
{

    extern vtkm::FloatDefault G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax;
    extern vtkm::FloatDefault GLOBAL_ADVECT_STEP_SIZE;
    extern int G_SampleX,G_SampleY,G_SampleZ;
    extern int GLOBAL_ADVECT_NUM_STEPS;
    extern int GLOBAL_ADVECT_NUM_SEEDS;
    extern int GLOBAL_NUM_LEVELS;
    extern std::ofstream *timingInfo;

    void runAdvection(vtkh::DataSet *data_set, int rank, int numRanks, int step,
                      std::string seedMethod, std::string fieldToOperateOn, bool cloverleaf,
                      bool recordTrajectories, bool outputResults, bool outputseeds);
    vtkh::DataSet *runGhostStripper(vtkh::DataSet *data_set, int rank, int numRanks,
                                    int step, std::string ghostFieldName);
}

#endif
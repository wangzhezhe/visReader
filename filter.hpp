#ifndef __FILTER_H
#define __FILTER_H

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

namespace FILTER
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
}

#endif
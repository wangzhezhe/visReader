#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <iostream>
#include <chrono>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#include <random>
#include <string>
#include <vector>

//--Fides
#include <fides/DataSetReader.h>
#include <fides/DataSetWriter.h>

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
#include <vtkm/exec/ParametricCoordinates.h>
#include <vtkm/exec/CellInterpolate.h>
#include <vtkm/filter/particleadvection/BoundsMap.h>

//--statistics
#include <vtkh/StatisticsDB.hpp>


using namespace std;

//debug value
bool verbose = false;
bool debugVTK = false;
bool writeStreamlines = true;


namespace vtkm
{
namespace worklet
{
struct CellCenter: public vtkm::worklet::WorkletVisitCellsWithPoints
{
    public:
    using ControlSignature = void (CellSetIn cellSet, FieldInPoint inputPointField, FieldOut outputCellField);
    using ExecutionSignature = void (_1, PointCount, _2, _3);
    using InputDomain = _1;
    template <typename CellShape, typename InputPointFieldType, typename OutputType>
    VTKM_EXEC void operator()(CellShape shape, vtkm::IdComponent numPoints,
    const InputPointFieldType &inputPointField, OutputType &centerOut) const
    {
        vtkm::Vec3f parametricCenter;
        vtkm::exec::ParametricCoordinatesCenter(numPoints, shape, parametricCenter);
        vtkm::exec::CellInterpolate(inputPointField, parametricCenter, shape, centerOut);
    }
};
}
}

void writeDataSet(vtkh::DataSet *data, std::string fName, int rank, int step)
{
  int numDomains = data->GetNumberOfDomains();
  std::cerr << "num domains " << numDomains << std::endl;
  for(int i = 0; i < numDomains; i++)
  {
    char fileNm[128];
    sprintf(fileNm, "%s.step%d.rank%d.domain%d.vtk", fName.c_str(), step, rank, i);
    vtkm::io::VTKDataSetWriter write(fileNm);
    write.WriteDataSet(data->GetDomain(i));
  }
}

void printLineOhSeeds(std::vector<vtkm::Particle> &seeds,
                vtkm::Particle startPoint, vtkm::Particle endPoint, int rank)
{
    if(rank == 0)
    {
        ofstream *seedFile = new ofstream;
        char nm[32];
        sprintf(nm, "generatedSeeds.out");
        seedFile->open(nm, ofstream::out);


        (*seedFile) << "x " << "y " << "z " << "value" << endl;

        /*(*seedFile) << startPoint.Pos[0] << " "
                    << startPoint.Pos[1] << " "
                    << startPoint.Pos[2] << " "
                    << "-1" << endl;

        (*seedFile) << endPoint.Pos[0] << " "
                    << endPoint.Pos[1] << " "
                    << endPoint.Pos[2] << " "
                    << "1" << endl;
*/
        for(long unsigned int i = 0; i < seeds.size(); i++)
        {
           (*seedFile) << seeds[i].Pos[0] << " "
                    << seeds[i].Pos[1] << " "
                    << seeds[i].Pos[2] << " "
                    << "0" << endl;
        }

        seedFile->close();
    }
}

void printBoxOhSeeds(std::vector<vtkm::Particle> &seeds,int rank, int step)
{
    if(rank == 0)
    {
        ofstream *seedFile = new ofstream;
        char nm[32];
        sprintf(nm, "generatedBoxOfSeeds_step%d.out", step);
        seedFile->open(nm, ofstream::out);


        (*seedFile) << "x " << "y " << "z " << "value" << endl;

        /*(*seedFile) << startPoint.Pos[0] << " "
                    << startPoint.Pos[1] << " "
                    << startPoint.Pos[2] << " "
                    << "-1" << endl;

        (*seedFile) << endPoint.Pos[0] << " "
                    << endPoint.Pos[1] << " "
                    << endPoint.Pos[2] << " "
                    << "1" << endl;
*/
        for(long unsigned int i = 0; i < seeds.size(); i++)
        {
           (*seedFile) << seeds[i].Pos[0] << " "
                    << seeds[i].Pos[1] << " "
                    << seeds[i].Pos[2] << " "
                    << "0" << endl;
        }

        seedFile->close();
    }
}

void printAllOhSeeds(std::vector<vtkm::Particle> &seeds, int rank, int step)
{
    ofstream *seedFile = new ofstream;
    char nm[64];
    sprintf(nm, "generatedSeeds_step%d_rank%d.out", step, rank);
    seedFile->open(nm, ofstream::out);


    (*seedFile) << "x " << "y " << "z " << "value" << endl;

    for(long unsigned int i = 0; i < seeds.size(); i++)
    {
       (*seedFile) << seeds[i].Pos[0] << " "
                << seeds[i].Pos[1] << " "
                << seeds[i].Pos[2] << " "
                << "0" << endl;
    }

    seedFile->close();
}


fides::metadata::Vector<std::size_t> assignWorkBlocks(int totalBlocks, int rank, int numRanks)
{
    fides::metadata::Vector<std::size_t> blockSelection;

    if(totalBlocks <= rank)
    {
        printf("Warning! :: Process[%i] is wasting cycles :: too many processors for data blocks\n", rank);
        return blockSelection;
    }

    //----Set indexes for each reader if there is more than one
    int endIndex;
    int startIndex = (totalBlocks / numRanks) * rank;
    if (totalBlocks % numRanks > rank)
    {
        startIndex += rank;
        endIndex = startIndex + (totalBlocks / numRanks) + 1;
    }
    else
    {
        startIndex += totalBlocks % numRanks;
        endIndex = startIndex + (totalBlocks / numRanks);
    }

    for(int i = startIndex; i < endIndex; i++)
    {
        blockSelection.Data.push_back(static_cast<std::size_t>(i));
    }

    return blockSelection;
}

vtkh::DataSet *runGhostStripper(vtkh::DataSet *data_set, int rank, int numRanks,
                                int step, string ghostFieldName);
#endif /* __UTILITIES_H */

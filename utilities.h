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

#ifdef USE_FIDES
#include <fides/DataSetReader.h>
#include <fides/DataSetWriter.h>
#endif

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

// debug value
static bool verbose = false;

namespace vtkm
{
    namespace worklet
    {
        struct CellCenter : public vtkm::worklet::WorkletVisitCellsWithPoints
        {
        public:
            using ControlSignature = void(CellSetIn cellSet, FieldInPoint inputPointField, FieldOut outputCellField);
            using ExecutionSignature = void(_1, PointCount, _2, _3);
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

inline std::chrono::steady_clock::time_point startTimer()
{
    return std::chrono::steady_clock::now();
}

// stop an existing timer and print timing message
inline float endTimer(std::chrono::steady_clock::time_point start)
{
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    return std::chrono::duration<double, std::milli>(diff).count();
}

inline void writeDataSet(vtkh::DataSet *data, std::string fName, int rank, int step)
{
    int numDomains = data->GetNumberOfDomains();
    std::cerr << "rank " << rank << " num domains " << numDomains << std::endl;
    for (int i = 0; i < numDomains; i++)
    {
        char fileNm[128];
        sprintf(fileNm, "%s.step%d.rank%d.domain%d.vtk", fName.c_str(), step, rank, i);
        vtkm::io::VTKDataSetWriter write(fileNm);
        write.WriteDataSet(data->GetDomain(i));
    }
}

inline void printLineOhSeeds(std::vector<vtkm::Particle> &seeds,
                             vtkm::Particle startPoint, vtkm::Particle endPoint, int rank)
{
    if (rank == 0)
    {
        ofstream *seedFile = new ofstream;
        char nm[32];
        sprintf(nm, "generatedSeeds.out");
        seedFile->open(nm, ofstream::out);

        (*seedFile) << "x "
                    << "y "
                    << "z "
                    << "value" << endl;

        /*(*seedFile) << startPoint.Pos[0] << " "
                    << startPoint.Pos[1] << " "
                    << startPoint.Pos[2] << " "
                    << "-1" << endl;

        (*seedFile) << endPoint.Pos[0] << " "
                    << endPoint.Pos[1] << " "
                    << endPoint.Pos[2] << " "
                    << "1" << endl;
*/
        for (long unsigned int i = 0; i < seeds.size(); i++)
        {
            (*seedFile) << seeds[i].Pos[0] << " "
                        << seeds[i].Pos[1] << " "
                        << seeds[i].Pos[2] << " "
                        << "0" << endl;
        }

        seedFile->close();
    }
}

inline void printBoxOhSeeds(std::vector<vtkm::Particle> &seeds, int rank, int step)
{
    if (rank == 0)
    {
        ofstream *seedFile = new ofstream;
        char nm[32];
        sprintf(nm, "generatedBoxOfSeeds_step%d.out", step);
        seedFile->open(nm, ofstream::out);

        (*seedFile) << "x "
                    << "y "
                    << "z "
                    << "value" << endl;

        /*(*seedFile) << startPoint.Pos[0] << " "
                    << startPoint.Pos[1] << " "
                    << startPoint.Pos[2] << " "
                    << "-1" << endl;

        (*seedFile) << endPoint.Pos[0] << " "
                    << endPoint.Pos[1] << " "
                    << endPoint.Pos[2] << " "
                    << "1" << endl;
*/
        for (long unsigned int i = 0; i < seeds.size(); i++)
        {
            (*seedFile) << seeds[i].Pos[0] << " "
                        << seeds[i].Pos[1] << " "
                        << seeds[i].Pos[2] << " "
                        << "0" << endl;
        }

        seedFile->close();
    }
}

inline void printAllOhSeeds(std::vector<vtkm::Particle> &seeds, int rank, int step)
{
    ofstream *seedFile = new ofstream;
    char nm[64];
    sprintf(nm, "generatedSeeds_step%d_rank%d.out", step, rank);
    seedFile->open(nm, ofstream::out);

    (*seedFile) << "x "
                << "y "
                << "z "
                << "value" << endl;

    for (long unsigned int i = 0; i < seeds.size(); i++)
    {
        (*seedFile) << seeds[i].Pos[0] << " "
                    << seeds[i].Pos[1] << " "
                    << seeds[i].Pos[2] << " "
                    << "0" << endl;
    }

    seedFile->close();
}

// test two blocks case, horizental division
#ifdef USE_FIDES
inline fides::metadata::Vector<std::size_t> assignWorkBlocks4(int totalBlocks, int rank, int numRanks)
{
    fides::metadata::Vector<std::size_t> blockSelection;
    for (int i = 0; i < totalBlocks; i++)
    {
        if (i == 2 || i == 3 || i == 6 || i == 7 || i == 10 || i == 11)
        {
            if (rank == 0)
            {
                blockSelection.Data.push_back(static_cast<std::size_t>(i));
                std::cout << "rank:" << rank << " id " << i << std::endl;
            }
        }
        else
        {
            if (rank != 0)
            {
                blockSelection.Data.push_back(static_cast<std::size_t>(i));
                std::cout << "rank:" << rank << " id " << i << std::endl;
            }
        }
    }

    return blockSelection;
}

// testing for 12 blocks 2 readers
inline fides::metadata::Vector<std::size_t> assignWorkBlocks3(int totalBlocks, int rank, int numRanks)
{
    fides::metadata::Vector<std::size_t> blockSelection;
    for (int i = 0; i < totalBlocks; i++)
    {
        if (rank == 0)
        {
            if (i == 10)
            {
                blockSelection.Data.push_back(static_cast<std::size_t>(i));
                std::cout << "rank:" << rank << " id " << i << std::endl;
            }
        }
        else
        {
            if (i != 10)
            {
                blockSelection.Data.push_back(static_cast<std::size_t>(i));
                std::cout << "rank:" << rank << " id " << i << std::endl;
            }
        }
    }

    return blockSelection;
}

inline fides::metadata::Vector<std::size_t> assignWorkBlocks2(int totalBlocks, int rank, int numRanks)
{
    fides::metadata::Vector<std::size_t> blockSelection;
    int blocksPerRank = totalBlocks / numRanks;
    std::cout << "blocksPerRank: " << blocksPerRank << std::endl;
    for (int i = 0; i < totalBlocks; i++)
    {
        if (i % numRanks == rank)
        {
            blockSelection.Data.push_back(static_cast<std::size_t>(i));
            std::cout << "rank:" << rank << " id " << i << std::endl;
        }
    }

    return blockSelection;
}

// maybe add more here
// one is continuous assigining, one is discontinuous assigning
// one is adaptive assigning
inline fides::metadata::Vector<std::size_t> assignWorkBlocks(int totalBlocks, int rank, int numRanks)
{
    fides::metadata::Vector<std::size_t> blockSelection;

    if (totalBlocks <= rank)
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

    std::cout << "rank " << rank << " allocated block index start:" << startIndex << " end: " << endIndex << std::endl;

    for (int i = startIndex; i < endIndex; i++)
    {
        blockSelection.Data.push_back(static_cast<std::size_t>(i));
    }

    return blockSelection;
}
#endif

#endif /* __UTILITIES_H */

#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <iostream>
#include <chrono>
#include <sys/time.h>
#include <mpi.h>
#include <random>
#include <string>
#include <vector>

#include <vtkm/Particle.h>
#include <vtkm/io/VTKDataSetWriter.h>

using namespace std;

// debug value
static bool verbose = false;

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


inline void writeDataSet(const vtkm::cont::PartitionedDataSet &pds, std::string fName, int rank, int step)
{
    int numDomains = pds.GetNumberOfPartitions();
    for (int i = 0; i < numDomains; i++)
    {
        char fileNm[128];
        sprintf(fileNm, "%s.step%d.rank%d.domain%d.vtk", fName.c_str(), step, rank, i);
        vtkm::io::VTKDataSetWriter write(fileNm);
        write.WriteDataSet(pds.GetPartition(i));
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
            (*seedFile) << seeds[i] << endl;
        }

        seedFile->close();
    }
}

inline void printBoxOhSeeds(std::vector<vtkm::Particle> &seeds, int rank, int step)
{
    //if (rank == 0)
    //{
        ofstream *seedFile = new ofstream;
        char nm[256];
        sprintf(nm, "generatedBoxOfSeeds_step%d_rank%d.out", step, rank);
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
            (*seedFile) << seeds[i] << endl;
        }

        seedFile->close();
    //}
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
        (*seedFile) << seeds[i] << endl;
    }

    seedFile->close();
}

#endif /* __UTILITIES_H */

#include "filter.hpp"
#include "utilities.h"
#include <vector>
#include <vtkm/filter/flow/internal/BoundsMap.h>
#include <vtkm/filter/flow/ParticleAdvection.h>
#include <vtkm/filter/flow/Streamline.h>

namespace FILTER
{
    vtkm::FloatDefault G_xMin = 0, G_xMax = 0, G_yMin = 0, G_yMax = 0, G_zMin = 0, G_zMax = 0;
    vtkm::FloatDefault GLOBAL_ADVECT_STEP_SIZE = 0.1;
    int G_SampleX = 10, G_SampleY = 10, G_SampleZ = 10;
    int GLOBAL_ADVECT_NUM_STEPS = 10;
    int GLOBAL_ADVECT_NUM_SEEDS = 100;
    int GLOBAL_NUM_LEVELS = 1;

    // the random number between 0 and 1
    static vtkm::FloatDefault random01()
    {
        return (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
    }

    // fixed seeds position for debug
    void createPointsOfSeeds(std::vector<vtkm::Particle> &seeds)
    {
        vtkm::Particle p = vtkm::Particle(vtkm::Vec3f(1, 1, 6), 1);
        seeds.push_back(p);
    }

    // specify how many sample in each dimension
    // instead of using the random seeds position
    void createBoxOfSeedsFixed(const vtkm::cont::PartitionedDataSet &pds,
                               std::vector<vtkm::Particle> &seeds,
                               vtkm::FloatDefault xMin,
                               vtkm::FloatDefault xMax,
                               vtkm::FloatDefault yMin,
                               vtkm::FloatDefault yMax,
                               vtkm::FloatDefault zMin,
                               vtkm::FloatDefault zMax,
                               int sampleNumX, int sampleNumY, int sampleNumZ,
                               int rank, int numRanks, int step)
    {

        float deltax = (xMax - xMin) / (sampleNumX + 1);
        float deltay = (yMax - yMin) / (sampleNumY + 1);
        float deltaz = (zMax - zMin) / (sampleNumZ + 1);

        int totalNumSeeds = sampleNumX * sampleNumY * sampleNumZ;

        std::vector<vtkm::Vec3f> allSeeds;

        for (int i = 0; i < sampleNumX; i++)
        {
            for (int j = 0; j < sampleNumY; j++)
            {
                for (int k = 0; k < sampleNumZ; k++)
                {
                    float x = xMin + (i + 1) * deltax;
                    float y = yMin + (j + 1) * deltay;
                    float z = zMin + (k + 1) * deltaz;
                    allSeeds.push_back({x, y, z});
                }
            }
        }
        // auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
        vtkm::Id numDomains = pds.GetNumberOfPartitions();
        std::vector<vtkm::cont::DataSet> dataSetVec;
        for (vtkm::Id i = 0; i < numDomains; i++)
            dataSetVec.push_back(pds.GetPartition(i));

        vtkm::filter::flow::internal::BoundsMap boundsMap(dataSetVec);
        // this is used for debug
        // if (rank == 0 && step == 0)
        //{
        //     boundsMap.BoundsInfo();
        // }
        //  select seeds that belongs to current domain that the rank owns
        if (allSeeds.size() != totalNumSeeds)
        {
            throw std::runtime_error("allSeeds.size()!=numSeeds");
        }
        for (int i = 0; i < totalNumSeeds; i++)
        {
            auto blockIds = boundsMap.FindBlocks(allSeeds[i]);
            if (!blockIds.empty() && boundsMap.FindRank(blockIds[0]) == rank)
            {
                seeds.push_back({allSeeds[i], i});
            }
        }

        // Counting the seeds number
        std::vector<int> seedCounts(numRanks, 0);
        seedCounts[rank] = seeds.size();
        // std::cout << "debug seed count rank " << rank << " " << seedCounts[rank] << std::endl;
        MPI_Allreduce(MPI_IN_PLACE, seedCounts.data(), numRanks, MPI_INT,
                      MPI_SUM, MPI_COMM_WORLD);
        int totNum = 0;
        for (int i = 0; i < numRanks; i++)
        {
            totNum += seedCounts[i];
        }
        if (totNum != totalNumSeeds)
        {
            std::cout << "Warn: valid num " << totNum << " total numSeeds " << totalNumSeeds << std::endl;
            // set the extends a little bit smaller to the actual one can avoid this issue
            // throw std::runtime_error("totNum is supposed to equal numSeeds");
        }
    }

    void createBoxOfSeedsRandom(const vtkm::cont::PartitionedDataSet &pds,
                                std::vector<vtkm::Particle> &seeds,
                                vtkm::FloatDefault xMin,
                                vtkm::FloatDefault xMax,
                                vtkm::FloatDefault yMin,
                                vtkm::FloatDefault yMax,
                                vtkm::FloatDefault zMin,
                                vtkm::FloatDefault zMax,
                                int numSeeds, int rank, int numRanks, int step, bool output)
    {

        // std::cout << xMin << " " << xMax << std::endl;
        // std::cout << rank << " " << numRanks << std::endl;
        //  Dave begin changes
        //  Set seeds for BBox
        std::vector<vtkm::Vec3f> allSeeds;

        // make sure all ranks use the same time and have the same seeds
        MPI_Barrier(MPI_COMM_WORLD);
        srand(time(NULL));
        for (int i = 0; i < numSeeds; i++)
        {
            float x = xMin + (xMax - xMin) * random01();
            float y = yMin + (yMax - yMin) * random01();
            float z = zMin + (zMax - zMin) * random01();
            allSeeds.push_back({x, y, z});
            // std::cout << "push seed " << x << " " << y << " " << z << std::endl;
        }

        // auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
        vtkm::Id numPartitions = pds.GetNumberOfPartitions();
        std::vector<vtkm::cont::DataSet> dataSetVec;
        for (vtkm::Id i = 0; i < numPartitions; i++)
            dataSetVec.push_back(pds.GetPartition(i));

        vtkm::filter::flow::internal::BoundsMap boundsMap(dataSetVec);
        // this is used for debug
        // if (rank == 0 && step == 0)
        //{
        //     boundsMap.BoundsInfo();
        // }
        //  select seeds that belongs to current domain that the rank owns
        //  std::cout << "numSeeds " << numSeeds << std::endl;
        if (allSeeds.size() != numSeeds)
        {
            throw std::runtime_error("allSeeds.size()!=numSeeds");
        }
        for (int i = 0; i < numSeeds; i++)
        {
            auto blockIds = boundsMap.FindBlocks(allSeeds[i]);
            if (!blockIds.empty() && boundsMap.FindRank(blockIds[0]) == rank)
            {
                seeds.push_back({allSeeds[i], i});
            }
            // else
            //{
            //     std::cout << "debug seed no bounds " << allSeeds[i][0] << " " << allSeeds[i][1] << " " << allSeeds[i][2] << std::endl;
            //     std::cout << blockIds.empty() << " " << boundsMap.FindRank(blockIds[0]) << " " << rank << std::endl;
            //     boundsMap.DebugFindBlocks(allSeeds[i]);
            // }
        }

        // Counting the seeds number
        std::vector<int> seedCounts(numRanks, 0);
        seedCounts[rank] = seeds.size();
        // std::cout << "debug seed count rank " << rank << " " << seedCounts[rank] << std::endl;
        MPI_Allreduce(MPI_IN_PLACE, seedCounts.data(), numRanks, MPI_INT,
                      MPI_SUM, MPI_COMM_WORLD);
        int totNum = 0;
        for (int i = 0; i < numRanks; i++)
        {
            totNum += seedCounts[i];
        }
        if (totNum != numSeeds)
        {
            std::cout << "Warn: totNum " << totNum << " actual numSeeds " << numSeeds << std::endl;
            // set the extends a little bit smaller to the actual one can avoid this issue
            // throw std::runtime_error("totNum is supposed to equal numSeeds");
        }

        if (output)
            printBoxOhSeeds(seeds, rank, step);
    }

    void runAdvection(const vtkm::cont::PartitionedDataSet &pds, int rank, int numRanks, int step, std::string seedMethod, std::string fieldToOperateOn, bool cloverleaf, bool recordTrajectories, bool outputResults, bool outputseeds)
    {

        std::vector<vtkm::Particle> seeds;

        // createLineOfSeeds(seeds, startPoint, endPoint, GLOBAL_ADVECT_NUM_SEEDS, rank);
        if (seedMethod == "point")
        {
            createPointsOfSeeds(seeds);
        }
        else if (seedMethod == "boxrandom")
        {
            // random seeding in box
            createBoxOfSeedsRandom(pds, seeds, G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax, GLOBAL_ADVECT_NUM_SEEDS, rank, numRanks, step, outputseeds);
        }
        else if (seedMethod == "boxfixed")
        {
            // fixed position in box
            createBoxOfSeedsFixed(pds, seeds, G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax, G_SampleX, G_SampleY, G_SampleZ, rank, numRanks, step);
        }
        else
        {
            throw std::runtime_error("unsupported seed method");
        }
        /*TODO
          adding these cases as needed
        else if (seedMethod == "cell")
        {
            if (cloverleaf)
                createSeedInEveryCellCloverleaf(seeds, pds, rank, numRanks, step, outputseeds);
            else
                createSeedInEveryCell(seeds, pds, rank, numRanks, step, outputseeds);
        }
        else if (seedMethod == "domain")
        {
            if (cloverleaf)
                createSeedInEveryDomainCloverleaf(seeds, pds, rank, numRanks, step, outputseeds);
            else
                createSeedInEveryDomain(seeds, pds, rank, numRanks, step, outputseeds);
        }
        */

        if (rank == 0)
        {
            
            std::cout << "advect maxSteps " << GLOBAL_ADVECT_NUM_STEPS << std::endl;
            std::cout << "advect seeds number " << GLOBAL_ADVECT_NUM_SEEDS << std::endl;
            std::cout << "advect step size " << GLOBAL_ADVECT_STEP_SIZE << std::endl;
        }

        std::chrono::steady_clock::time_point filterStart = std::chrono::steady_clock::now();

        if (recordTrajectories)
        {
            vtkm::filter::flow::Streamline streamline;
            streamline.SetSeeds(seeds);
            streamline.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
            streamline.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
            streamline.SetActiveField(fieldToOperateOn);
            auto streamlineOutput = streamline.Execute(pds);

            // TODO, maybe setting the timestep as a separate parameter
            // outside of the logging class

            if (outputResults)
            {
                writeDataSet(streamlineOutput, "advection_streamlinesOutput", rank, step);
            }

            // delete streamline_output;
        }
        else
        {
            vtkm::filter::flow::ParticleAdvection pa;
            pa.SetSeeds(seeds);
            pa.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
            pa.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
            pa.SetActiveField(fieldToOperateOn);
            auto paOutput = pa.Execute(pds);
        }

        float filterTime = endTimer(filterStart);

        // output execution time
        if (recordTrajectories)
        {
            if (rank == 0)
            {
                fprintf(stderr, "%i, VISapp_%i_%i, advect(streamline), %f\n\n", step, rank, numRanks,
                        filterTime);
            }
        }
        else
        {
            if (rank == 0)
            {
                fprintf(stderr, "%i, VISapp_%i_%i, advect(particleadev), %f\n\n", step, rank, numRanks,
                        filterTime);
            }
        }
    }

}

#include "vtkhfilter.hpp"
#include "utilities.h"

namespace VTKH_FILTER
{
    vtkm::FloatDefault G_xMin = 0, G_xMax = 0, G_yMin = 0, G_yMax = 0, G_zMin = 0, G_zMax = 0;
    vtkm::FloatDefault GLOBAL_ADVECT_STEP_SIZE = 0.1;
    int GLOBAL_ADVECT_NUM_STEPS = 10;
    int GLOBAL_ADVECT_NUM_SEEDS = 100;
    int GLOBAL_NUM_LEVELS = 1;
    std::ofstream *timingInfo = NULL;

    vtkm::rendering::Camera createWindCamera(vtkm::Bounds bounds, float zoom)
    {
        vtkm::rendering::Camera camera;
        camera.ResetToBounds(bounds);

        vtkm::FloatDefault myElevation = -75;
        camera.Elevation(myElevation);
        //-make the camera equal to ascent
        /*vtkm::Vec<vtkm::Float32, 3> look = {5,5,5};
        vtkm::Vec<vtkm::Float32, 3> pos = {7.60472, 13.6603, 19.7721};
        camera.SetLookAt(look);
        camera.SetPosition(pos);
        camera.SetClippingRange(1.73205, 173.205);
        camera.Zoom(zoom);*/
        //
        return camera;
    }

    vtkm::rendering::Camera createCamera(vtkm::Bounds bounds, float zoom)
    {
        vtkm::rendering::Camera camera;
        camera.ResetToBounds(bounds);

        vtkm::FloatDefault myAzimuth = 65;
        vtkm::FloatDefault myElevation = 30;
        camera.Azimuth(myAzimuth);
        camera.Elevation(myElevation);
        //-make the camera equal to ascent
        /*vtkm::Vec<vtkm::Float32, 3> look = {5,5,5};
        vtkm::Vec<vtkm::Float32, 3> pos = {7.60472, 13.6603, 19.7721};
        camera.SetLookAt(look);
        camera.SetPosition(pos);
        camera.SetClippingRange(1.73205, 173.205);
        camera.Zoom(zoom);*/
        //
        return camera;
    }
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

    void createLineOfSeeds(std::vector<vtkm::Particle> &seeds, vtkm::Particle startPoint,
                           vtkm::Particle endPoint, int numSeeds, int rank, bool output)
    {
        vtkm::Particle v;
        v.Pos = startPoint.Pos - endPoint.Pos;
        // cerr << "Vector from points = " << v.Pos <<endl;

        vtkm::FloatDefault t, dt;
        if (numSeeds == 1)
        {
            t = 0.5;
            dt = 0.5;
        }
        else
        {
            t = 0.0;
            dt = 1.0 / (vtkm::FloatDefault)(numSeeds - 1);
        }

        for (int i = 0; i < numSeeds; i++)
        {
            vtkm::Particle p;
            p.Pos = endPoint.Pos + t * v.Pos;
            p.ID = static_cast<vtkm::Id>(i);
            seeds.push_back(p);
            t = t + dt;
        }

        if (output)
            printLineOhSeeds(seeds, startPoint, endPoint, rank);
    }

    void createBoxOfSeeds(vtkh::DataSet *data,
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
            //std::cout << "push seed " << x << " " << y << " " << z << std::endl;
        }

        // auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
        vtkm::Id numDomains = data->GetNumberOfDomains();
        std::vector<vtkm::cont::DataSet> dataSetVec;
        for (vtkm::Id i = 0; i < numDomains; i++)
            dataSetVec.push_back(data->GetDomain(i));

        vtkm::filter::particleadvection::BoundsMap boundsMap(dataSetVec);
        if (rank == 0 && step==0)
        {
            boundsMap.BoundsInfo();
        }
        // select seeds that belongs to current domain that the rank owns
        // std::cout << "numSeeds " << numSeeds << std::endl;
        if(allSeeds.size()!=numSeeds){
            throw std::runtime_error("allSeeds.size()!=numSeeds");
        }
        for (int i = 0; i < numSeeds; i++)
        {
            auto blockIds = boundsMap.FindBlocks(allSeeds[i]);
            if (!blockIds.empty() && boundsMap.FindRank(blockIds[0]) == rank)
            {
                seeds.push_back({allSeeds[i], i});
            }
            //else
            //{
            //    std::cout << "debug seed no bounds " << allSeeds[i][0] << " " << allSeeds[i][1] << " " << allSeeds[i][2] << std::endl;
            //    std::cout << blockIds.empty() << " " << boundsMap.FindRank(blockIds[0]) << " " << rank << std::endl;
            //    boundsMap.DebugFindBlocks(allSeeds[i]);
            //}
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

        /*
        original way to create seeds
        vtkm::Particle diff, startPoint, endPoint;

        startPoint.Pos = {xMin, yMin, zMin};
        endPoint.Pos = {xMax, yMax, zMax};
        diff.Pos = endPoint.Pos - startPoint.Pos;

        for (int i = 0; i < numSeeds; i++)
        {
            vtkm::Particle p;
            p.Pos = {startPoint.Pos[0] + (diff.Pos[0] * random01()),
                     startPoint.Pos[1] + (diff.Pos[1] * random01()),
                     startPoint.Pos[2] + (diff.Pos[2] * random01())};
            p.ID = static_cast<vtkm::Id>(i);
            seeds.push_back(p);
        }

        //Make the paricle ID's unique
        std::vector<int> particlesPerRank(numRanks, 0);
        particlesPerRank[rank] = numSeeds;
        MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int offset = 0;
        for (int i = 0; i < rank; i++)
          offset += particlesPerRank[i];

        if (offset > 0)
        {
          for (auto& p : seeds)
            p.ID += offset;
        }
        */

        if (output)
            printBoxOhSeeds(seeds, rank, step);
    }

    void createSeedInEveryCellCloverleaf(std::vector<vtkm::Particle> &seeds,
                                         vtkh::DataSet *data_set,
                                         int rank,
                                         int numRanks,
                                         int step, bool output)
    {
        using AxisArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
        using CartesianProduct = vtkm::cont::ArrayHandleCartesianProduct<AxisArrayType, AxisArrayType, AxisArrayType>;

        int particleCount = 0;
        vtkm::Id numDomains = data_set->GetNumberOfDomains();
        for (vtkm::Id i = 0; i < numDomains; i++)
        {
            auto tempDS = data_set->GetDomain(i);
            vtkm::cont::CellSetStructured<3> tempCS =
                tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
            // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
            auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<CartesianProduct>();

            vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
            vtkm::cont::Invoker invoke;
            invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

            auto buff = cellCenters.ReadPortal();

            for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
            {
                vtkm::Particle p;
                p.Pos = buff.Get(z);
                p.ID = static_cast<vtkm::Id>(particleCount);
                seeds.push_back(p);
                particleCount++;
            }
        }

        // Make the paricle ID's unique
        std::vector<int> particlesPerRank(numRanks, 0);
        particlesPerRank[rank] = particleCount;
        MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int offset = 0;
        for (int i = 0; i < rank; i++)
            offset += particlesPerRank[i];

        if (offset > 0)
        {
            for (auto &p : seeds)
                p.ID += offset;
        }

        if (output)
            printAllOhSeeds(seeds, rank, step);
    }

    void createSeedInEveryCell(std::vector<vtkm::Particle> &seeds,
                               vtkh::DataSet *data_set,
                               int rank,
                               int numRanks,
                               int step, bool output)
    {
        using UniformCoordType = vtkm::cont::ArrayHandleUniformPointCoordinates;

        int particleCount = 0;
        vtkm::Id numDomains = data_set->GetNumberOfDomains();
        for (vtkm::Id i = 0; i < numDomains; i++)
        {
            auto tempDS = data_set->GetDomain(i);
            vtkm::cont::CellSetStructured<3> tempCS =
                tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
            // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
            auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<UniformCoordType>();

            vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
            vtkm::cont::Invoker invoke;
            invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

            auto tempGhosts = tempDS.GetField("topo_ghosts");
            auto ghostArr = tempGhosts.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
            const vtkm::FloatDefault *ghostBuff = ghostArr.GetReadPointer();
            auto buff = cellCenters.ReadPortal();

            for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
            {
                if (ghostBuff[z] == 0)
                {
                    vtkm::Particle p;
                    p.Pos = buff.Get(z);
                    p.ID = static_cast<vtkm::Id>(particleCount);
                    seeds.push_back(p);
                    particleCount++;
                }
            }
        }

        // Make the paricle ID's unique
        std::vector<int> particlesPerRank(numRanks, 0);
        particlesPerRank[rank] = particleCount;
        MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int offset = 0;
        for (int i = 0; i < rank; i++)
            offset += particlesPerRank[i];

        if (offset > 0)
        {
            for (auto &p : seeds)
                p.ID += offset;
        }

        if (output)
            printAllOhSeeds(seeds, rank, step);
    }

    // Method to create a seed in every domain in a data set.
    // The seed is randomely placed within the bounds of the domain.
    void createSeedInEveryDomainCloverleaf(std::vector<vtkm::Particle> &seeds,
                                           vtkh::DataSet *data_set,
                                           int rank,
                                           int numRanks,
                                           int step, bool output)
    {
        using AxisArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
        using CartesianProduct = vtkm::cont::ArrayHandleCartesianProduct<AxisArrayType, AxisArrayType, AxisArrayType>;

        int particleCount = 0;
        vtkm::Id numDomains = data_set->GetNumberOfDomains();
        for (vtkm::Id i = 0; i < numDomains; i++)
        {
            auto tempDS = data_set->GetDomain(i);
            vtkm::cont::CellSetStructured<3> tempCS =
                tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
            // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
            auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<CartesianProduct>();

            vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
            vtkm::cont::Invoker invoke;
            invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);
            auto buff = cellCenters.ReadPortal();

            for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
            {
                vtkm::Particle p;
                p.Pos = buff.Get(z);
                p.ID = static_cast<vtkm::Id>(particleCount);
                seeds.push_back(p);
                particleCount++;
                break;
            }
        }

        // Make the paricle ID's unique
        std::vector<int> particlesPerRank(numRanks, 0);
        particlesPerRank[rank] = particleCount;
        MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int offset = 0;
        for (int i = 0; i < rank; i++)
            offset += particlesPerRank[i];

        if (offset > 0)
        {
            for (auto &p : seeds)
                p.ID += offset;
        }

        if (output)
            printAllOhSeeds(seeds, rank, step);
    }

    // Method to create a seed in every domain in a data set.
    // The seed is randomely placed within the bounds of the domain.
    void createSeedInEveryDomain(std::vector<vtkm::Particle> &seeds,
                                 vtkh::DataSet *data_set,
                                 int rank,
                                 int numRanks,
                                 int step, bool output)
    {
        using UniformCoordType = vtkm::cont::ArrayHandleUniformPointCoordinates;

        int particleCount = 0;
        vtkm::Id numDomains = data_set->GetNumberOfDomains();
        for (vtkm::Id i = 0; i < numDomains; i++)
        {
            auto tempDS = data_set->GetDomain(i);
            vtkm::cont::CellSetStructured<3> tempCS =
                tempDS.GetCellSet().Cast<vtkm::cont::CellSetStructured<3>>();
            // auto t = tempDS.GetCoordinateSystem(0).GetDataAsMultiplexer();
            auto t = tempDS.GetCoordinateSystem().GetData().AsArrayHandle<UniformCoordType>();

            vtkm::cont::ArrayHandle<vtkm::Vec3f> cellCenters;
            vtkm::cont::Invoker invoke;
            invoke(vtkm::worklet::CellCenter{}, tempCS, t, cellCenters);

            auto tempGhosts = tempDS.GetField("topo_ghosts");
            auto ghostArr = tempGhosts.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
            const vtkm::FloatDefault *ghostBuff = ghostArr.GetReadPointer();
            auto buff = cellCenters.ReadPortal();

            for (int z = 0; z < cellCenters.GetNumberOfValues(); z++)
            {
                if (ghostBuff[z] == 0)
                {
                    vtkm::Particle p;
                    p.Pos = buff.Get(z);
                    p.ID = static_cast<vtkm::Id>(particleCount);
                    seeds.push_back(p);
                    particleCount++;
                    break;
                }
            }
        }

        // Make the paricle ID's unique
        std::vector<int> particlesPerRank(numRanks, 0);
        particlesPerRank[rank] = particleCount;
        MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int offset = 0;
        for (int i = 0; i < rank; i++)
            offset += particlesPerRank[i];

        if (offset > 0)
        {
            for (auto &p : seeds)
                p.ID += offset;
        }

        if (output)
            printAllOhSeeds(seeds, rank, step);
    }

    void runAdvection(vtkh::DataSet *data_set, int rank, int numRanks, int step, std::string seedMethod, std::string fieldToOperateOn, bool cloverleaf, bool recordTrajectories, bool outputResults, bool outputseeds)
    {

        vtkh::EVENT_BEGIN("streamline/particleadv_total");
        vtkh::EVENT_BEGIN("streamline_generate_seeds");
        std::vector<vtkm::Particle> seeds;

        // createLineOfSeeds(seeds, startPoint, endPoint, GLOBAL_ADVECT_NUM_SEEDS, rank);
        if (seedMethod == "point")
        {
            createPointsOfSeeds(seeds);
        }
        else if (seedMethod == "box")
        {
            /*
                void createBoxOfSeeds(vtkh::DataSet *data,
                          std::vector<vtkm::Particle> &seeds,
                          vtkm::FloatDefault xMin,
                          vtkm::FloatDefault xMax,
                          vtkm::FloatDefault yMin,
                          vtkm::FloatDefault yMax,
                          vtkm::FloatDefault zMin,
                          vtkm::FloatDefault zMax,
                          int numSeeds, int rank, int numRanks, int step, bool output)
    {
            */
            createBoxOfSeeds(data_set, seeds, G_xMin, G_xMax, G_yMin, G_yMax, G_zMin, G_zMax, GLOBAL_ADVECT_NUM_SEEDS, rank, numRanks, step, outputseeds);
        }
        else if (seedMethod == "cell")
        {
            if (cloverleaf)
                createSeedInEveryCellCloverleaf(seeds, data_set, rank, numRanks, step, outputseeds);
            else
                createSeedInEveryCell(seeds, data_set, rank, numRanks, step, outputseeds);
        }
        else if (seedMethod == "domain")
        {
            if (cloverleaf)
                createSeedInEveryDomainCloverleaf(seeds, data_set, rank, numRanks, step, outputseeds);
            else
                createSeedInEveryDomain(seeds, data_set, rank, numRanks, step, outputseeds);
        }
        vtkh::EVENT_END("streamline_generate_seeds");

        if (rank == 0)
        {
            std::cerr << "advect maxSteps " << GLOBAL_ADVECT_NUM_STEPS << std::endl;
            std::cerr << "advect seeds number " << GLOBAL_ADVECT_NUM_SEEDS << std::endl;
            std::cerr << "advect step size " << GLOBAL_ADVECT_STEP_SIZE << std::endl;
        }

        std::chrono::steady_clock::time_point filterStart = std::chrono::steady_clock::now();

        if (recordTrajectories)
        {
            vtkh::Streamline streamline;
            vtkh::EVENT_BEGIN("streamline_setup");
            streamline.SetSeeds(seeds);
            streamline.SetInput(data_set);
            streamline.SetField(fieldToOperateOn);
            streamline.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
            streamline.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
            streamline.SetIterationStep(step);
            vtkh::EVENT_END("streamline_setup");
            vtkh::DataSet *streamline_output = NULL;

            vtkh::EVENT_BEGIN("streamline_update");
            streamline.Update();
            vtkh::EVENT_END("streamline_update");

            vtkh::EVENT_BEGIN("streamline_get_output");
            streamline_output = streamline.GetOutput();
            vtkh::EVENT_END("streamline_get_output");

            if (outputResults)
            {
                vtkh::EVENT_BEGIN("streamline_save_advection_files");
                writeDataSet(streamline_output, "advection_streamlinesOutput", rank, step);
                vtkh::EVENT_END("streamline_save_advection_files");
            }

            delete streamline_output;
        }
        else
        {
            vtkh::ParticleAdvection pa;
            vtkh::EVENT_BEGIN("particle_advection_setup");
            pa.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
            pa.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);
            pa.SetSeeds(seeds);
            pa.SetField(fieldToOperateOn);
            pa.SetInput(data_set);
            pa.SetIterationStep(step);
            vtkh::EVENT_END("particle_advection_setup");
            vtkh::EVENT_BEGIN("particle_advection_update");
            pa.Update();
            vtkh::EVENT_END("particle_advection_update");

            vtkh::DataSet *particleAdvectOutput = NULL;

            vtkh::EVENT_BEGIN("particle_advection_get_output");
            particleAdvectOutput = pa.GetOutput();
            vtkh::EVENT_END("particle_advection_get_output");

            // there is no streamline for the particle advection filter
            delete particleAdvectOutput;
        }

        vtkh::EVENT_END("streamline/particleadv_total");
        float filterTime = endTimer(filterStart);
        if (recordTrajectories)
        {
            if (rank == 0)
            {
                fprintf(stderr, "\n%i, VISapp_%i_%i, advect(streamline), %f\n", step, rank, numRanks,
                        filterTime);
                (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                              << ", advect, " << filterTime << endl;
            }
        }
        else
        {
            if (rank == 0)
            {
                fprintf(stderr, "\n%i, VISapp_%i_%i, advect(particleadev), %f\n", step, rank, numRanks,
                        filterTime);
                (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                              << ", advect, " << filterTime << endl;
            }
        }

        //--
    }

    void runVolumeRender(vtkh::DataSet *data_set, int rank, int numRanks,
                         int step, bool cloverleaf, std::string fieldToOperateOn)
    {
        int height = 1024;
        int width = 1024;

        //---- Time the render
        std::chrono::steady_clock::time_point start = startTimer();
        vtkm::Bounds bounds = data_set->GetGlobalBounds();

        vtkh::EVENT_BEGIN("createRendering");
        vtkh::Scene scene;
        // scene.SetRenderBatchSize(1);
        //-- use this to run multiple renders of the same dataset
        for (int j = 0; j < 1; j++)
        {
            char imgNm[128];
            if (cloverleaf)
                sprintf(imgNm, "clover_out.step-%05i_render#%02i", step, j);
            else
                sprintf(imgNm, "wind_out.step-%05i_render#%02i", step, j);
            vtkm::rendering::Camera cam;

            if (cloverleaf)
                cam = createCamera(bounds, .2);
            else
                cam = createWindCamera(bounds, .2);

            float bg_color[4] = {0.f, 0.f, 0.f, 1.f};
            vtkh::Render render =
                vtkh::MakeRender(height, width, cam, *data_set, imgNm, bg_color);
            render.DoRenderAnnotations(false);
            // vtkh::Scene scene;
            scene.AddRender(render);

            /*vtkh::VolumeRenderer volume;
            volume.SetInput(&data_set);
            volume.SetField("pressure");
            volume.SetNumberOfSamples(1000);
            scene.AddRenderer(&volume);
            scene.Render();   */
        }

        vtkh::VolumeRenderer volume;
        volume.SetInput(data_set);
        volume.SetField(fieldToOperateOn);
        volume.SetNumberOfSamples(100);
        scene.AddRenderer(&volume);
        vtkh::EVENT_END("createRendering");

        vtkh::EVENT_BEGIN("render");
        scene.Render();
        vtkh::EVENT_END("render");

        float totalTime = endTimer(start);
        if (rank == 0)
        {
            fprintf(stderr, "\n%i, VISapp_%i_%i, render, %f", step, rank, numRanks,
                    totalTime);
            (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                          << ", render, " << totalTime << endl;
        }

        //--
    }

    // ray trace an image
    void runRender(vtkh::DataSet iso_output, int rank, int numRanks, int step, bool cloverleaf, std::string fieldToOperateOn)
    {
        int height = 1024;
        int width = 1024;

        char imgNm[128];
        if (cloverleaf)
            sprintf(imgNm, "clover_out.step-%05d", step);
        else
            sprintf(imgNm, "wind_out.step-%05d", step);

        //---- Time the render
        std::chrono::steady_clock::time_point start = startTimer();
        vtkm::Bounds bounds = iso_output.GetGlobalBounds();

        vtkh::EVENT_BEGIN("createRendering");
        auto cam6 = createCamera(bounds, 0.0);

        float bg_color[4] = {0.f, 0.f, 0.f, 1.f};
        vtkh::Render render =
            vtkh::MakeRender(height, width, cam6, iso_output, imgNm, bg_color);
        render.DoRenderAnnotations(false);
        vtkh::Scene scene;
        vtkh::RayTracer tracer;
        tracer.SetInput(&iso_output);
        tracer.SetField(fieldToOperateOn);
        scene.AddRenderer(&tracer);
        scene.AddRender(render);
        vtkh::EVENT_END("createRendering");

        vtkh::EVENT_BEGIN("render");
        scene.Render();
        vtkh::EVENT_END("render");

        float totalTime = endTimer(start);
        fprintf(stderr, "\n%i, VISapp_%i_%i, render, %f", step, rank, numRanks,
                totalTime);
        (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                      << ", render, " << totalTime << endl;
        //--
    }

    vtkh::DataSet *runContour(vtkh::DataSet *data_set, int rank, int numRanks,
                              int step, bool cloverleaf, std::string fieldToOperateOn)
    {
        //---- Time the contour operation
        vtkh::EVENT_BEGIN("contour_setup");
        std::chrono::steady_clock::time_point start = startTimer();
        vtkh::MarchingCubes marcher;
        marcher.SetInput(data_set);
        marcher.SetField(fieldToOperateOn);

        const int num_levels = GLOBAL_NUM_LEVELS;
        marcher.SetLevels(num_levels);
        vtkh::EVENT_END("contour_setup");

        vtkh::EVENT_BEGIN("contour");
        marcher.Update();
        vtkh::EVENT_END("contour");

        vtkh::EVENT_BEGIN("contour_getOutput");
        vtkh::DataSet *iso_output = marcher.GetOutput();
        float totalTime = endTimer(start);
        fprintf(stderr, "\n%i, VISapp_%i_%i, contour, %f", step, rank, numRanks,
                totalTime);
        (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                      << ", contour, " << totalTime << endl;
        //--
        vtkh::EVENT_END("contour_getOutput");
        return iso_output;
    }

    vtkh::DataSet *runGhostStripper(vtkh::DataSet *data_set, int rank, int numRanks,
                                    int step, string ghostFieldName)
    {
        //---- Time the ghosts
        vtkh::EVENT_BEGIN("ghost_stripper");
        std::chrono::steady_clock::time_point start = startTimer();
        vtkh::GhostStripper ghost;
        ghost.SetInput(data_set);
        ghost.SetField(ghostFieldName);

        ghost.SetMaxValue(0);
        ghost.SetMinValue(0);
        ghost.Update();

        vtkh::DataSet *ghosts = ghost.GetOutput();
        float totalTime = endTimer(start);
        fprintf(stderr, "\n%i, VISapp_%i_%i, ghost_stripper, %f", step, rank,
                numRanks, totalTime);
        (*timingInfo) << step << ", VISapp_" << rank << "_" << numRanks
                      << ", ghost_stripper, " << totalTime << endl;
        vtkh::EVENT_END("ghost_stripper");
        //--
        return ghosts;
    }

}

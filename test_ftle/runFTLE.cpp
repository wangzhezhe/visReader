
// test code is in
// https://gitlab.kitware.com/vtk/vtk-m/-/blob/master/vtkm/filter/flow/testing/UnitTestLagrangianStructuresFilter.cxx?ref_type=heads

#include <iostream>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/filter/flow/LagrangianStructures.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/filter/flow/Tracer.h>
int main(int argc, char *argv[])
{
    vtkm::cont::Initialize(argc, argv);
    if (argc != 4)
    {
        std::cout << "<binary> <input file name> <field name> <output file name>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    std::string outputFileName = argv[3];

    vtkm::filter::flow::GetTracer().Get()->Init(0);
    vtkm::filter::flow::GetTracer().Get()->StartTimer();

    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inputData = reader.ReadDataSet();

    std::cout << "input data " << std::endl;
    inputData.PrintSummary(std::cout);

    vtkm::filter::flow::LagrangianStructures lagrangianStructures;
    lagrangianStructures.SetStepSize(0.025f);
    lagrangianStructures.SetNumberOfSteps(100);
    lagrangianStructures.SetAdvectionTime(0.025f * 100);
    lagrangianStructures.SetActiveField(fieldName);
    vtkm::cont::DataSet outputData = lagrangianStructures.Execute(inputData);
    vtkm::filter::flow::GetTracer().Get()->Finalize();

    outputData.PrintSummary(std::cout);

    //std::string outputFileName = "output_ftle.vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(outputData);
}
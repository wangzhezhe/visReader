
// test code is in
// https://gitlab.kitware.com/vtk/vtk-m/-/blob/master/vtkm/filter/flow/testing/UnitTestLagrangianStructuresFilter.cxx?ref_type=heads

#include <iostream>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/filter/flow/LagrangianStructures.h>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cout << "<binary> <file name> <field name>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];

    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inputData = reader.ReadDataSet();

    vtkm::filter::flow::LagrangianStructures lagrangianStructures;
    lagrangianStructures.SetStepSize(0.025f);
    lagrangianStructures.SetNumberOfSteps(500);
    lagrangianStructures.SetAdvectionTime(0.025f * 500);
    lagrangianStructures.SetActiveField(fieldName);
    vtkm::cont::DataSet outputData = lagrangianStructures.Execute(inputData);

    outputData.PrintSummary(std::cout);

    std::string outputFileName = "output_ftle.vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(outputData);
}
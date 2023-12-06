cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io vtkm::filter_flow)
message(STATUS "VTKm Libraries ${VTKm_LIBRARIES}")

add_executable(runFTLE runFTLE.cpp)
target_compile_definitions(runFTLE PRIVATE "MPI_ENABLED")
target_link_libraries(runFTLE ${VTKm_LIBRARIES})

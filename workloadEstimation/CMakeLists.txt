##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io vtkm::filter_flow)
message(STATUS "VTKm Libraries ${VTKm_LIBRARIES}")

add_executable(StreamlineMPI StreamlineMPI.cxx Block.cxx Block.h ../LoadData.cpp ../filter.cpp)
target_compile_definitions(StreamlineMPI PRIVATE "MPI_ENABLED")
target_link_libraries(StreamlineMPI  ${VTKm_LIBRARIES} MPI::MPI_CXX)

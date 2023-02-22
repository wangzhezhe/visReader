cmake_minimum_required(VERSION 3.3)

PROJECT(reader)
add_compile_options(-std=c++14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

##########################
find_package(MPI REQUIRED)

##########################
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
message(STATUS "VTKm Libraries ${VTKm_LIBRARIES}")

###########################
OPTION (USE_GPU OFF)
message ("USE_GPU is setting as: " ${USE_GPU})

if(USE_GPU)

  set_source_files_properties(visitReaderAdev.cpp filter.cpp PROPERTIES LANGUAGE "CUDA")
  add_executable(visitReaderAdev visitReaderAdev.cpp filter.cpp)
  set_target_properties(visitReaderAdev PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(visitReaderAdev ${VTKm_LIBRARIES} MPI::MPI_CXX)

else()
  
  add_executable(visitReaderAdev visitReaderAdev.cpp)
  #add_executable(visitReaderAdev visitReaderAdev.cpp filter.cpp)
  target_link_libraries(visitReaderAdev ${VTKm_LIBRARIES} MPI::MPI_CXX)

endif()

add_subdirectory(workloadEstimation)

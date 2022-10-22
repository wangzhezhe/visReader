cmake_minimum_required(VERSION 3.3)

PROJECT(reader)
add_compile_options(-std=c++11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

option(ENABLE_TBB         "Build TBB Support"         OFF)
option(ENABLE_OMP         "Build OMP Support"         OFF)

##########################
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

##########################
find_package(Fides REQUIRED)
message(STATUS "Fides found")
set(Fides_LIBRARIES fides)
message(STATUS "Fides Libraries ${Fides_LIBRARIES}")

##########################
find_package(VTKm REQUIRED)
set(VTKm_LIBRARIES vtkm_cont vtkm_filter vtkm_rendering)
message(STATUS "VTKm Libraries ${VTKm_LIBRARIES}")
message(STATUS "VTKm Directories ${VTKm_INCLUDE_DIRS}")
link_directories(${VTKm_LIBRARIES})
include_directories(${VTKm_INCLUDE_DIRS})

##########################
find_package(VTKh REQUIRED)
set(VTKH_LIBRARIES "vtkh_mpi")
include_directories(${VTKh_INCLUDE_DIRS})
link_directories(${VTKH_LIBRARIES})
message(STATUS "VTKh libraries ${VTKH_LIBRARIES}")
message(STATUS "VTKh directories ${VTKH_INCLUDE_DIRS}")

##########################
if(NOT ADIOS2_DIR)
    MESSAGE(FATAL_ERROR "ADIOS2 support needs explicit ADIOS2_DIR")
endif()

MESSAGE(STATUS "Looking for ADIOS2 using ADIOS2_DIR = ${ADIOS2_DIR}")

find_package(ADIOS2 REQUIRED)

if(NOT ADIOS2_FOUND)
    message(FATAL_ERROR "ADIOS2_DIR is not a path to a valid ADIOS2 install")
endif()

message(STATUS "ADIOS2 found")
message(STATUS "${ADIOS2_INCLUDE_DIRS}")
message(STATUS "${ADIOS2_LIBRARIES}")
include_directories(${ADIOS2_INCLUDE_DIRS})
link_directories(${ADIOS2_DIR}/lib)
###########################

OPTION (USE_GPU OFF)
message ("USE_GPU is setting as: " ${USE_GPU})
if(USE_GPU)

  set_source_files_properties(reader.cpp PROPERTIES LANGUAGE "CUDA")
  add_executable(reader reader.cpp)
  set_target_properties(reader PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(reader ${ADIOS2_LIBRARIES})
  target_link_libraries(reader ${MPI_LIBRARIES})
  target_link_libraries(reader ${VTKm_LIBRARIES})
  target_link_libraries(reader ${VTKH_LIBRARIES})
  target_link_libraries(reader ${Fides_LIBRARIES})
  target_compile_options(reader PRIVATE "-DVTKH_ENABLE_STATISTICS")
  target_compile_options(reader PRIVATE "-DVTKH_PARALLEL")

  if(MPI_COMPILE_FLAGS)
    set_target_properties(reader PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(reader PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()

else()

  add_executable(reader reader.cpp)
  target_link_libraries(reader ${ADIOS2_LIBRARIES})
  target_link_libraries(reader ${MPI_LIBRARIES})
  target_link_libraries(reader ${VTKm_LIBRARIES})
  target_link_libraries(reader ${VTKH_LIBRARIES})
  target_link_libraries(reader ${Fides_LIBRARIES})
  #target_link_libraries(reader dmalloc)
  target_compile_options(reader PRIVATE "-DVTKH_ENABLE_STATISTICS")
  target_compile_options(reader PRIVATE "-DVTKH_PARALLEL")

  if(MPI_COMPILE_FLAGS)
    set_target_properties(reader PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(reader PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()

endif()


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -MMD -g -O0 -Wall -std=c++11 -Wno-unused-variable -Wno-reorder")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -MMD -O3 -Wall -std=c++11 -Wno-unused-variable -Wno-reorder")

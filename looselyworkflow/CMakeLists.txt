cmake_minimum_required (VERSION 3.14)
project (lcwf C CXX)
set (CMAKE_CXX_STANDARD 14)

enable_testing()
add_definitions(-DENABLE_TIMERS)

find_package(MPI REQUIRED) 

# add our cmake module directory to the path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
     "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# link shared lib with full rpath
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# setup cache variables for ccmake
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release
         CACHE STRING "Choose the type of build." FORCE)
    set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS 
                  "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif ()
set (CMAKE_PREFIX_PATH "" CACHE STRING "External dependencies path")
set (BUILD_SHARED_LIBS "OFF" CACHE BOOL "Build a shared library")

# packages we depend on (load all the necessary depedency by spack)
find_package (cereal CONFIG REQUIRED)
get_target_property (CEREAL_INC cereal::cereal INTERFACE_INCLUDE_DIRECTORIES)

find_package (PkgConfig REQUIRED)
pkg_check_modules (margo REQUIRED IMPORTED_TARGET margo)
# for the latest thallium we do not need to link libraries
find_package (thallium REQUIRED)

include_directories(
     ${lcwf_SOURCE_DIR}/dep
     ${lcwf_SOURCE_DIR}
)

# find vtk data
find_package(VTK COMPONENTS
IOLegacy
IOXML
IOImage
CommonDataModel
FiltersSources
FiltersCore
ParallelMPI
)



find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)

add_executable(tightlyinsitu_rrb tightlyinsitu_rrb.cpp ${lcwf_SOURCE_DIR}/../LoadData.cpp ${lcwf_SOURCE_DIR}/../filter.cpp)
target_link_libraries(tightlyinsitu_rrb thallium ${VTKm_LIBRARIES} ${VTK_LIBRARIES} MPI::MPI_CXX)

add_executable(tightlyinsitu_pa tightlyinsitu_pa.cpp ${lcwf_SOURCE_DIR}/../LoadData.cpp ${lcwf_SOURCE_DIR}/../filter.cpp)
target_link_libraries(tightlyinsitu_pa thallium ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(looselyinsitu looselyinsitu.cpp ${lcwf_SOURCE_DIR}/../LoadData.cpp ${lcwf_SOURCE_DIR}/../filter.cpp)
target_link_libraries(looselyinsitu thallium ${VTKm_LIBRARIES} ${VTK_LIBRARIES} MPI::MPI_CXX)

add_executable(tightlyinsitu_webylog tightlyinsitu_webylog.cpp ${lcwf_SOURCE_DIR}/../LoadData.cpp ${lcwf_SOURCE_DIR}/../filter.cpp)
target_link_libraries(tightlyinsitu_webylog thallium ${VTKm_LIBRARIES} ${VTK_LIBRARIES} MPI::MPI_CXX)

add_executable(tightlyinsitu_webyest tightlyinsitu_webyest.cpp
${lcwf_SOURCE_DIR}/../LoadData.cpp 
${lcwf_SOURCE_DIR}/../filter.cpp 
${lcwf_SOURCE_DIR}/../workloadEstimation/Block.cxx 
${lcwf_SOURCE_DIR}/../workloadEstimation/Block.h
${lcwf_SOURCE_DIR}/../workloadEstimation/SLTool.cxx 
${lcwf_SOURCE_DIR}/../workloadEstimation/SLTool.h)
target_link_libraries(tightlyinsitu_webyest thallium ${VTKm_LIBRARIES} ${VTK_LIBRARIES} MPI::MPI_CXX)

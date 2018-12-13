cmake_minimum_required(VERSION 3.0)

project(using_with_cmake)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include("FindAscent.cmake")
include("FindConduit.cmake")
if(VTKM_DIR)
    include("FindVTKm.cmake")
    include("FindVTKh.cmake")
endif()

# setup the ascent & conduit include paths
include_directories(${ASCENT_INCLUDE_DIRS})
include_directories(${CONDUIT_INCLUDE_DIRS})

# create our example 
add_executable(sim_nyx_3D nyx_sim.cxx)

# link to ascent
target_link_libraries(sim_nyx_3D ascent ${VTK_LIBRARIES})

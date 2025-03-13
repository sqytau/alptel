#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(alptelsim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Define option for the fast simulation in SteppingAction.cc
#
option(USE_LX_FAST_SIM "Fast simulation option, when particles are discarded when they cross the surface of the beam dumps and shielding." OFF)
if (USE_LX_FAST_SIM)
  add_compile_definitions(-DLX_FAST_SIM)
endif()

# Define option for the simulation timing cut in SteppingAction.cc
#
option(USE_LX_TRACK_TIMING_CUT "Simulation option to stop particle propagation after a given time (ns)." NO)
if (USE_LX_TRACK_TIMING_CUT)
  add_compile_definitions(LX_TRACK_TIMING_CUT=${USE_LX_TRACK_TIMING_CUT})
endif()


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(alptelsim alptelsim.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(alptelsim ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES} ${HDF5_LIBRARIES})

#target_compile_options(alptelsim PUBLIC -ggdb -g)
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -g")


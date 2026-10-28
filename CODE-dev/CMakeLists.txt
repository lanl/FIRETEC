cmake_minimum_required(VERSION 3.10)

# Project Information
project(CASPM) # Coupled Atmosphere-Surface Process Modeling
set(CASPM_VERSION_MAJOR 1)

# Compiler information
enable_language(Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  message(STATUS "Using GNU Compiler")
  add_compile_options(-ffixed-line-length-none -fallow-argument-mismatch)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "INTEL")
  message(STATUS "Using INTEL Compiler")
  set(FF_FLAGS "-extend-source")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "CRAY")
  message(STATUS "Using CRAY Compiler")
  set(FF_FLAGS "-extend_source")
endif()

# Mac Silicon Speedup:
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
  message(STATUS "Using Apple Silicon compile flag: -mcpu=apple-m1")
  add_compile_options(-mcpu=apple-m1)
endif()

add_compile_options(-O2 -g)

# Option for debug flags
option(DEBUG "Enable debug flags" OFF) # default is OFF unless: cmake .. -DDEBUG=ON
if(DEBUG)
  add_compile_options(-Wall -Wextra -Werror -fbounds-check -pedantic -fcheck=all)
  message(STATUS "Building with DEBUG flags enabled")
else()
  message(STATUS "Building with DEBUG flags disabled")
endif()

# Option for precision
option(DBL_PREC "Use double precision" OFF) # default is single precision unless: cmake .. -DDBL_PREC=ON
if(DBL_PREC)
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "DBL_PREC") # COMMENT OUT IF SINGLE PRECISION
  add_compile_definitions(DBL_PREC)
  message(STATUS "Building with DOUBLE PRECISION")
else()
  message(STATUS "Building with SINGLE PRECISION")
endif()

# Find dependencies
find_package(MPI REQUIRED COMPONENTS Fortran)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Module dependecies
add_subdirectory(subGridTurb)
add_subdirectory(Vegetation)
add_subdirectory(Radiation)
add_subdirectory(Emissions)
add_subdirectory(Firebrands)
add_subdirectory(HeatSource)
add_subdirectory(Sensors)
add_subdirectory(xeVariation)

# Check if the optional LANL auxillary directory exists (different repository)
set(LANL_MOD_SUBDIR_PATH "${CMAKE_CURRENT_SOURCE_DIR}/firetec_lanl_module")
if(EXISTS "${LANL_MOD_SUBDIR_PATH}")
  message(STATUS "Building optional library in ${MY_SUBDIR_PATH}")
  # Add the subdirectory if it exists
  add_subdirectory(firetec_lanl_module)

  # You might also want to add a define to indicate its presence
  # This can be useful for conditional compilation in your source code
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "LANL_MOD")
else()
  message(STATUS "Optional LANL library directory ${MY_SUBDIR_PATH} not found. Skipping build.")
endif()

# Build
set(Higrad_Firetec_SRCS
  variables.F90 msga.F90 io.F90 advectionScheme.F90 
  higrad.F90 boundary.F90 compress.F90 definearray.F90 
  domainUpdates.F90 force.F90 input.F90 largeScalePGF.F90 
  metric.F90 MOA.F90 readgridlist.F90 rinit.F90
  rmaxmin.F90 thermodynamics.F90 tinit.F90 topo.F90 update.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/opt)
add_executable(firerad ${Higrad_Firetec_SRCS})
add_dependencies(firerad Radiation xeVariation)
set_target_properties(firerad PROPERTIES COMPILE_FLAGS "${FF_FLAGS}")
set_target_properties(firerad PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
target_link_libraries(firerad ${MPI_Fortran_LIBRARIES})
target_link_libraries(firerad Embers Emissions HeatSource Radiation
  Sensors subGridTurb Vegetation xeVariation)
  
if(EXISTS "${LANL_MOD_SUBDIR_PATH}")
  # Link the library to your main executable or another target
  target_link_libraries(firerad LANL_MOD)
endif()

# Make Clean Instructions
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
  "opt/*;firerad;Emissions/libEmissions.a;xeVariation/libxeVariation.a;
   HeatSource/libHeatSource.a;Sensors/libSensors.a;Vegetation/libVegetation.a;
   Radiation/libRadiation.a")

cmake_minimum_required(VERSION 3.10)

# Project Information
project(CASPM) # Coupled Atmosphere-Surface Process Modeling
set(CASPM_VERSION_MAJOR 1)

# Compiler information
enable_language(Fortran)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(FF_FLAGS "-O2 -g -ffixed-line-length-none -Wall -Wextra -Werror")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set(FF_FLAGS "-O2 -g -extend-source")
elseif (Fortran_COMPILER_NAME MATCHES "ftn.*")
  set(FF_FLAGS "-O2 -g -extend_source")
endif()
# add_compile_options(-Wall -Wextra)

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Module dependecies
#set_directory_properties(PROPERTIES -Wall -Wextra COMPILE_DEFINITIONS "DBL_PREC") # COMMENT OUT IF SINGLE PRECISION
add_subdirectory(Emissions)
add_subdirectory(xeVariation)
add_subdirectory(subGridTurb)
add_subdirectory(HeatSource)
add_subdirectory(Firebrands)
add_subdirectory(Sensors)

# Build
set(Higrad_Firetec_SRCS
  msga.F90 radiation.F90
  io.F90 advectionScheme.F90 higrad.F90
  boundary.F90 compress.F90 convection.F90 definearray.F90 
  diffusesies.F90 domainUpdates.F90 drag.F90 
  fuel.F90 fuel_functions.F90 filstr.F90 force.F90 
  ignite.F90 input.F90 largeScalePGF.F90 metric.F90 
  MOA.F90 readgridlist.F90 rinit.F90 rinitfire.F90 
  rmaxmin.F90 thermodynamics.F90 tinit.F90 topo.F90 update.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/opt)
add_executable(firerad ${Higrad_Firetec_SRCS})
set_target_properties(firerad PROPERTIES COMPILE_FLAGS "${FF_FLAGS}")
set_target_properties(firerad PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
target_link_libraries(firerad ${MPI_Fortran_LIBRARIES})
target_link_libraries(firerad subGridTurb Emissions xeVariation HeatSource Sensors Embers)

# Make Clean Instructions
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
  "opt/*;firerad;Emissions/libEmissions.a;xeVariation/libxeVariation.a;
     HeatSource/libHeatSource.a;Sensors/libSensors.a")

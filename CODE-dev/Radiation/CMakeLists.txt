cmake_minimum_required(VERSION 3.10)

# Project Information
project(RAD) # Radiation mechanisms
set(RAD_VERSION_MAJOR 1)
set(RAD_VERSION_MINOR 1.0)
set(RAD_VERSION "${RAD_VERSION_MAJOR}.${RAD_VERSION_MINOR}")

# Compiler Information

# Build
set(RAD_SRCS ../variables.F90                         # Other General Dependencies
  ../readgridlist.F90
  rinit_rad.F90 variables_rad.F90 input_rad.F90       # Standard Files
  definearray_rad.F90                             
  monteCarlo_rad.F90 sink_rad.F90                     # Source Calculation Files
  solarHeat_rad.F90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(Radiation STATIC ${RAD_SRCS})
add_dependencies(Radiation Vegetation)

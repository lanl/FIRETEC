cmake_minimum_required(VERSION 3.10)

# Project Information
project(HEAT) # Fire-Dependent Emissions Modeling
set(HEAT_VERSION_MAJOR 1)
set(HEAT_VERSION_MINOR 1.0)
set(HEAT_VERSION "${HEAT_VERSION_MAJOR}.${HEAT_VERSION_MINOR}")

# Build
set(HEATSOURCE_SRCS ../variables.F90 ../metric.F90       # Other Dependencies
  rinit_heat.F90 variables_heat.F90 input_heat.F90       # Standard files
  heatsource_heat.F90) # Source calculation files

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(HeatSource STATIC ${HEATSOURCE_SRCS})

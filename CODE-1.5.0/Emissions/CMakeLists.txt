cmake_minimum_required(VERSION 3.10)

# Project Information
project(EMIT) # Fire-Dependent Emissions Modeling
set(EMIT_VERSION_MAJOR 1)
set(EMIT_VERSION_MINOR 1.0)
set(EMIT_VERSION "${EMIT_VERSION_MAJOR}.${EMIT_VERSION_MINOR}")

# Build
set(EMISSIONS_SRCS ../variables.F90 ../metric.F90     
  ../fuel_functions.F90                               # Other Dependencies
  rinit_ems.F90 variables_ems.F90 input_ems.F90       # Standard files
  pointSource_ems.F90 efSource_ems.F90 zBEST_ems.F90  # Source calculation files
  stokesDrag_ems.F90 deposit_ems.F90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(Emissions STATIC ${EMISSIONS_SRCS})

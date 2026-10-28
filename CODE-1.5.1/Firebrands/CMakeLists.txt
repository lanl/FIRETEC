cmake_minimum_required(VERSION 3.10)

# Project Information
project(EMBER) # Fire-Dependent Emissions Modeling
set(EMBER_VERSION_MAJOR 1)
set(EMBER_VERSION_MINOR 1.0)
set(EMBER_VERSION "${EMBER_VERSION_MAJOR}.${EMBER_VERSION_MINOR}")

# Compiler Information
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_options(-fallow-argument-mismatch)
endif()

# Build
set(EMBER_SRCS ../variables.F90 ../metric.F90         # Other Dependencies
  ../readgridlist.F90 ../input.F90
  ../HeatSource/variables_heat.F90
  rinit_fb.F90 variables_fb.F90 input_fb.F90          # Standard files
  functions_fb.F90
  firebrand_fb.F90 generation_fb.F90 evolve_fb.F90    # Source calculation files
  spotting_fb.F90 io_fb.F90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(Embers STATIC ${EMBER_SRCS})
add_dependencies(Embers Emissions)

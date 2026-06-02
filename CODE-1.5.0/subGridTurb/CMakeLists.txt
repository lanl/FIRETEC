cmake_minimum_required(VERSION 3.10)

# Project Information
project(SUB_TURB) # Alterations to environmental variables
set(SUB_TURB_VERSION_MAJOR 1)
set(SUB_TURB_VERSION_MINOR 1.0)
set(SUB_TURB_VERSION "${SUB_TURB_VERSION_MAJOR}.${SUB_TURB_VERSION_MINOR}")

# Build
set(SUB_TURB_SRCS ../variables.F90 ../metric.F90 ../domainUpdates.F90 # Other Dependencies
  variables_turb.F90 defineArray_turb.F90 domainUpdates_turb.F90      # Standard files
  force_turb.F90 input_turb.F90 rinit_turb.F90
  drag_turb.F90 diffuse_turb.F90 rijgradu_turb.F90                    # Source calculation files
  stressrij_turb.F90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(subGridTurb STATIC ${SUB_TURB_SRCS})

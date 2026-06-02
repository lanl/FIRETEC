cmake_minimum_required(VERSION 3.10)

# Project Information
project(XEVAR) # Alterations to environmental variables
set(XEVAR_VERSION_MAJOR 1)
set(XEVAR_VERSION_MINOR 1.0)
set(XEVAR_VERSION "${XEVAR_VERSION_MAJOR}.${XEVAR_VERSION_MINOR}")

# Build
set(XEVAR_SRCS ../variables.F90 ../io.F90 ../metric.F90    # Other Dependencies
  rinit_xeVar.F90 variables_xeVar.F90 input_xeVar.F90      # Standard files
  windfield.F90 xevariation.F90) # Source calculation files

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(xeVariation STATIC ${XEVAR_SRCS})

cmake_minimum_required(VERSION 3.10)

# Project Information
project(SENSE) # Fire-Dependent Emissions Modeling
set(SENSE_VERSION_MAJOR 1)
set(SENSE_VERSION_MINOR 1.0)
set(SENSE_VERSION "${SENSE_VERSION_MAJOR}.${SENSE_VERSION_MINOR}")

# Build
set(SENSOR_SRCS ../variables.F90 ../metric.F90              # Other Dependencies
  rinit_sense.F90 variables_sense.F90 input_sense.F90       # Standard files
  sensors.F90) # Source calculation files

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(Sensors STATIC ${SENSOR_SRCS})

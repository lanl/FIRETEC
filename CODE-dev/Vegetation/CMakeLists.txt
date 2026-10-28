cmake_minimum_required(VERSION 3.10)

# Project Information
project(VEG) # Vegetation mechanisms
set(VEG_VERSION_MAJOR 1)
set(VEG_VERSION_MINOR 1.0)
set(VEG_VERSION "${VEG_VERSION_MAJOR}.${VEG_VERSION_MINOR}")

# Compiler Information

# Build
set(VEG_SRCS ../variables.F90 ../readgridlist.F90     # Other General Dependencies
  ../rmaxmin.F90
  ../subGridTurb/variables_turb.F90                   # Subgrid Turbulence Dependency
  ../Radiation/variables_rad.F90                      # Radiation Dependency
  rinit_veg.F90 variables_veg.F90 input_veg.F90       # Standard Files
  definearray_veg.F90                             
  burn_veg.F90 convection_veg.F90 diffusesies_veg.F90 # Source Calculation Files
  drag_veg.F90 evaporation_veg.F90 ignite_veg.F90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(Vegetation STATIC ${VEG_SRCS})
add_dependencies(Vegetation subGridTurb)

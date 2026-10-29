#-------------------------------------------------------------------------

© 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 
for Los Alamos National Laboratory (LANL), which is operated by Triad 
National Security, LLC for the U.S. Department of Energy/National Nuclear 
Security Administration. All rights in the program are reserved by 
Triad National Security, LLC, and the U.S. Department of Energy/National 
Nuclear Security Administration. The Government is granted for itself and 
others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, 
and to permit others to do so.

#-------------------------------------------------------------------------

FIRETEC OSS LICENSE is O5086

#-------------------------------------------------------------------------

This branch hosts the latest Check-point version of FIRETEC-1.5.1. Use the 'develop' branch for latest features and current developments.

In this repository, we have three directories:

    1) CODE-xxxx - contains the source code for the software, see the README.md in that directory to see how to configure and build the executable of the software.

    2) TOOLS - contains a number of developed example input files and post-analysis script used to help set up and interpret FIRETEC simulations.

    3) DOCS - the future destination of documentation as it's developed and posted.

#-------------------------------------------------------------------------

FIRETEC-1.5.1 builds on FIRETEC-1.5.0 with continued code modulization and restructuring to improve maintainability and extensibility.

***Major Changes from FIRETEC-1.5.0***

**Further Modulization**

Continued separation of capabilities into independent modules:
- **Radiation Module**: Radiation physics (previously a single 1000+ line file) has been restructured into a dedicated module with separate components for Monte Carlo radiation transport, solar heating, radiative sinks, initialization, and input handling.
- **Vegetation Module**: All vegetation-related physics (fuel dynamics, combustion, drag, convection, diffusion, evaporation, and ignition) consolidated into a single module, replacing the previously scattered implementation across multiple root-level files.

**Enhanced Build System**

- Improved CMake configuration with better compiler detection and platform-specific optimizations
- Added Apple Silicon optimization support (M1/M2 processors)
- Separated debug and optimization flags with cleaner option handling
- Optional double precision build flag (DBL_PREC) with improved implementation
- Support for optional LANL auxiliary modules via separate repository integration

**Code Architecture**

- Eliminated monolithic files (con.F90, convection.F90, diffusesies.F90, drag.F90, fuel.F90, fuel_functions.F90, filstr.F90, ignite.F90, rinitfire.F90, radiation.F90) in favor of modular structure
- Vegetation and Radiation modules now follow consistent naming conventions (_veg, _rad suffixes)
- Better separation of concerns with dedicated subdirectories for each major capability

**Variables Module Enhancement**

- Expanded variables.F90 to better manage global state and dependencies across modules (~19KB vs ~15KB in 1.5.0)

This version maintains backward compatibility in terms of physics implementation while providing a more maintainable and extensible codebase for future development.

***New Capabilities***

**Mirrored Boundary Capabilities**

In addition to the relaxed and cyclic boundary conditions, a mirriored boundary condition was added which will mirror internal domain conditions in domain boundary ghost cells and reflect vector conditions orthaganol to the boundary itself.

**Ambient Heating**

Solar heating of fuels with solar angles and heat fluxes calculated. Additional bulk moisture evaporation rates calculated due to this heating.

**Aerosol Physics**

Additional forcing terms introduced for particule coagulation, settling, and deposition affecting aerosols. Explicit aerosol tracking with the method of moments for mass distributions introduced.

***Bug Fixes***

Many, see the git logs

***Code Clean-up***

Whole code re-worked in general clean-up and to be compatible with non-specific compilers. Detail is too extensive to list (every file) but the end result uses the same physics in the same operations but with different architecture and orders.

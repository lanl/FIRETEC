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

This branch hosts the latest Check-point version of FIRETEC-1.5.0. Use the 'develop' branch for latest features and current developments.

In this repository, we have three directories:
    1) CODE-xxxx - contains the source code for the software, see the README.md in that directory to see how to configure and build the executable of the software.
    2) TOOLS - contains a number of developed example input files and post-analysis script used to help set up and interpret FIRETEC simulations.
    3) DOCS - the future destination of documentation as it's developed and posted.

FIRETEC-1.5.0 builds on FIRETEC-1.4.1 but consists of several major code restructurings, new capabilities, and extensive general code clean-up.

***Code Restructuring***

**HIGRAD/FIRETEC Integration**

Rather than running HIGRAD and FIRETEC as two separate codes that have partially but not fully merged together, we’ve restructured the code as a numerical solver that then can call up functionalities (not distinguishing between HIGRAD & FIRETEC at this point) to accomplish certain tasks. This numerical solver: 1) Gathers user input, 2) Sets up computer architecture, and 3) Initializes each simulation (including if a restart occurs) prior to the iterating timesteps of the simulation. After entering the iterating timesteps the solver will: 1) Update domains with any global dependencies (pressure, tempg, etc.), 2) Calculate explicit forcing functions, both for large timesteps (to be used once per timestep), and for small timesteps (to be used at each step of the implicit solver, 3) Apply an implicit solver, MOA is exclusively implemented and used for now but there is infrastructure for Runga-Kutta to be implemented in the future, 4) Domain updates and forcing calculations done as needed within the implicit solver, and 5) IO.

**Modulization**

Individual capabilities of the code are being isolated and implemented as independent modules that can be called up specifically to accomplish tasks (e.g. compute forcings), almost like libraries. The overall code uses the modules but isn’t necessarily dependent on any given one to run (compiler predirectives still to be implemented but it’s set-up to do this). This allows for easy transferability of these capabilities to other applications and codes, and also allows for the modules to be placed at various levels of export-control if needed as we pursue an open-source license.

**CMake**

We’ve replaced the building makefiles with cmake which allows for better integration of the modulization while also easing buildability across various platforms.

**Thermodynamics**

Thermodynamic variables (pressure, heat capacity, etc) are now consolidated and updated regularly rather than ad hoc when needed. These variables also now account for changing chemistry as fluid composition changes.

***New Capabilities***

**Emissions**

New module that allows a user to generate various emissions (gas or aerosols) with three distinct methodologies: 1) fixed-source emissions such as from a power-plant, gas-leak, etc, 2) fixed-emission factor emissions where a user-specifies specific fire emission factors, or 3) physics-based emission factors using Z-BEST model as before, currently limited in species.

**Artificial Heat Source**

User may specify a wide range of artificial heat sources throughout the domain.

**Sensors**

User may specify the code to output specifics variables and specific locations and time frequencies.

**xeVariation**

User may now bring in dynamic environmental conditions in either sparse or dense data formats.

**Double Precision**

Changed in the CMakeLists, one can now switch everything to double precision if needed

***Bug Fixes***

Reworked fuel bed velocities with multiple fuels, previous working was inconsistent with quality checks.

***Code Clean-up***

Whole code re-worked in general clean-up and to be compatible with non-specific compilers. Detail is too extensive to list (every file) but the end result uses the same physics in the same operations but with different architecture and orders.

All files now within a 72 character-column limit in case compiler limited.

Firebrands reworked back into the code.


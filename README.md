#-------------------------------------------------------------------------

© 2026. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

#-------------------------------------------------------------------------

FIRETEC OSS LICENSE is O5086

#-------------------------------------------------------------------------

This branch will hosts a "stable" version of FIRETEC updated from Version 1.3.2.

FIRETEC-1.4.0 builds on FIRETEC-1.3.2 but adds two major capabilities, fixes a few small bugs, and contains some general code clean-up.

#------New Capabilities
**Multiple Fuel Types**

We’ve extended the fuel arrays to account for multiple fuel types resolved simultaneously on the same grid.

To use multiple fuel types, specify the number of fuel types within the gridlist file with nfuel = x. The default number of fuel types is 1 and if nfuel is not specified in the gridlist then FIRETEC will default to 1 fuel type so old scripts and gridlists should run fine. If you use nfuel > 1, then you’ll need to use ivegread=1 as well. Multiple fuels are read into the simulation in series from the .dat files meaning each .dat file needs to contain a binary array for each fuel type stacked one on top of another. For example, if we have nfuel = 2, then the treesrhof.dat file should contain the saved array for fuel type 1 and then immediately another saved array for fuel type 2 within the same file. Likewise for treesmoist.dat, treesfueldepth.dat, and treesss.dat.

The output of each simulation will now include multiple arrays for each fuel variable (i.e. rhof_0, rhof_1, rhow_0, rho_1, etc.) and required a re-ordering of the order which those variables are outputed. To convert these output to vts use the updated convert_firetec_binary_to_VTK.py provided in REPO/TOOLS/Make_VTS_Files/Fire_Run. For now, wind_run variables have not been altered and are not as thoroughly tested.

For now, the variables that can be varied for the different fuel types are density of fuel (rhof), moisture content (fmoist), base fuel depth (afd), fuel size scale (ss). In the future, we’ll look into expanding the chemistry variables as well (heat of combustion, heat capacity, char-balance, etc.)

Multiple fuel types evaluation is accomplished by looping over the combustion, drag, and radiation of each fuel type individually and then summing there joint effect into the energy, momentum, or mass equations. This led to major changes in drag*.f files, fuel.f, and radiation.f files with minor changes in compress.f, con.f, convection.f, definearray.f, firebrand.f, firetec.f, ignite.f, io.f, rinit.f, rinitfire.f, startFromFile.f, stressrij.f, tinit.f, turb.f, update*.f, variables.f, and windfield.f.

The order of variables.f outputted changed to accommodate multiple fuels. The postprocessing script (binary to vts) was updated in the repository in the TOOLS/Make_VTS_Files/Fire_Run.

**Water Vapor**

We’ve added the ability to generate and transport water vapor (rhovapor) in the higrad array.

To use rhovapor, specify irhovapor = 1 in the gridlist. This will add another element to nv, at index 8, which represents the bulk density of water in the gas phase. To visualize this variable in vts files make sure to add rhovapor to the convert_firetec_binary_to_VTK.py to the variables list after O2 and before the density in the gas_field_names list. Unless specified, irhovapor will default to zero. 
Water vapor is generated from two sources. First, the evaporation of moisture from the fuel which simply moves water from rhow (in the solid phase) to rhovapor (in the gas phase) for transport processes. Second, as a byproduct of combustion, where 1 kg of the consumed fuel mass generates 0.56 kg of water vapor. Both of these generations are accomplished in firetec.f.

Water vapor alters the gas thermal properties (heat capacity and gas constant) of the air up to the point where the air is saturated (humidity=100%) after that then the water vapor is technically in two phases (gas and aerosol) and the aerosol no longer directly affects the air’s thermal properties, more development happening here.

#-------Bug Fixes

Removed lines from firetec.f (lines 35-44 in FIRETEC-1.3.2) which forced off the combustion of fuels until the heating time of the ignite.dat file was finished. This especially caused atv ignition and other ignition types to not work.

Changed the ambient and initial concentration of O2 to be 0.23 of air instead of 0.21 as the code is looking at mass fraction not mole fraction.

gas_dtp is only being applied on the cells where there is fuel instead of all cells. Where there is no fuel, gas_ltp (same operation but on the large timestep instead) is applied.

In radiation.f checked rhof(i,j,kreal) against min_rhof instead of rhof(i,j,k)

Fixed multiple references to wrong temperature (temg instead of temps) in radiation.f

A few other small things which will be updated as they’re found/remembered

When starting from a fire run I noticed you’d get really low oxygen concentrations at the start of your run. This was caused because we are not reading in oxygen data from the windrun and thus xvbdataold had zero for all oxygen. Thus when you did interpolations along the boundaries you were doing a linear interpolation between ambient oxygen (0.23) and zero thus giving you really low concentrations to ripple like a wave through a simulation. To fix this, changes done to rinit.f. Essentially, I’m moving the call rinitfire() which will initiate our oxygen and rhovapor before the call of startFromFile() which initializes the xvbdataold array. Now xvbdataold will be initialized with ambient oxygen (and water vapor) values.

The second problem also comes when reading a wind run. The problem was occurring if I ran with water vapor turned on, then my density was going to zero at the boundaries. The problem was that when reading in the xvbdata files we were assigning density values by the 8th index of the xvb array, but when we added water vapor (or emissions, or non-local burning) we expanded that array to more than 8 variables. To fix a change in windfield.f was made changing the 8th element call with an nv element call.

#------Code Clean-up

Changed spacing in files to be more pythonish in their indentation (pythonish not special, just wanted to be consistent)

Removed numerous redundant or archaic variables and computations that were left over from legacy operations (removed all reference to idirt and its computations)

Removed a number of old debugging codes that were commented out.

Consolidated rmaxmin functions to one file.

Updated makefile to only build relevant files and make clean to remove all produced executables.


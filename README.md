#-------------------------------------------------------------------------

© 2026. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

#-------------------------------------------------------------------------

FIRETEC OSS LICENSE is O5086

#-------------------------------------------------------------------------

This branch will hosts a "stable" version of FIRETEC updated from Version 1.4.0 with a code rewritten to be compiler non-specific.\

FIRETEC-1.4.1 builds on FIRETEC-1.4.0 but adds three major capabilities, one new pre-processing tool, several code restructurings, and contains extensive general code clean-up.

#------New Capabilities

**Diffuse Sies**

Diffusion of sies, the fuel’s internal energy, can now be turned on using idiffsies=1. This function effectively diffuses elevated fuel temperatures to surrounding fuels enabling the spread through that mechanism and not only through convective forces (advective forces for you non-engineers). If your fire requires creeping flow, use this capability.

**Multiple Fuels in Radiation**

Radiation now breaks up the bottom cell into multiple layers to properly resolve multiple fuels.

#-----Code Restructuring

**Dynamic Pointer System**

Added a dynamic pointer system to the x arrays to point to specific variables in the xv, xvb, and xe arrays. The routine is called xvpointers in the definearray.f file. Essentially, this function constructs the format of the xv, and xe arrays at simulation time and fills out a system of pointers for reference later in the code:

iuvel = u velocity

ivvel = v velocity

iwvel = w velocity

itemp = potential temperature

ika = total kinetic energy a

ikb = total kinetic energy b

iO2 = O2 mass fraction

ivapor = water vapor mass fraction

iM0 = particulate bulk number density

iM1 = particulate bulk mass density

imixfrac = gas mixture fraction (for non-local burning)

irho = gas density

Now when you write code and need to access a variable in a cell you should write xv(i,j,k,itemp) instead of xv(i,j,k,4). This was needed as we add more capabilities that are overlapping/competing for that coveted #8 spot. In the past xvb(i,j,k,8) could mean water vapor, or it could mean M0, or it could mean mixture fraction depending on what flags you were running with, and different flags were often incompatible with each other.

**XVB Gone**

Eliminated the xvb versus xv copy and executions so now there is only the transported array (xv) represented throughout the code

**Gridlist defaults**

All variables read-in from the gridlist are defaulted in the code prior to gridlist read-in meaning that if values are not included in the gridlist they will still have a value. Values found at the top of variables.f90.

**Ignition**

All ignitions are now done as aerial ignitions (no more igntypes) with a namelist at the top of the file reading target temperature, ramp rate, and number of ignition points to follow. Each point contains cartesian coordinates (x, y, and z cell) followed by time of ignition. Ignition maker included in tools along with an example ignition file in Basic Run.

#------Bug Fixes

Many

#-------Code Clean-up

Whole code re-worked in general clean-up and to be compatible with non-specific compilers. Detail is too extensive to list (every file) but the end result uses the same physics in the same operations but with different architecture and orders.


All files now updated to Fortran90 and now have .f90 extension.
Initially, firebrands not included but will be later down the road.


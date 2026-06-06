
## Commands below for pink - must have loaded modules with 
#              module load lampi intel

## Linux cluster (ie, pink) options - select only one FFLAGS definition
#----------------------------------------
  FC = mpifort   ## compiler
  #MPILIB =  -L${MPI_ROOT}/lib64 -lmpi -I${MPI_ROOT}/include/
  #MPILIB =  -I${MPI_ROOT}/include/
  FFLAGS = -O2 -ffixed-line-length-none#-extend-source
  #FFLAGS = -g -Wall -ffree-line-length-none#-extend-source
  #FFLAGS = -g -extend-source -warn all,noexternal,nointerface -check bounds -check overflow

## COMPAQ options - select only one FFLAGS definition
#----------------------------------------
# FC = f90    ## compiler
# MPILIB =  -lmpi -lfmpi -lelan
# FMACH = ${MPI_COMPILE_FLAGS} ${MPI_LD_FLAGS} 
# DEBUG_OPTIONS = -check bounds -check overflow
# FFLAGS =   -g  ${FMACH} ${DEBUG_OPTIONS} #-convert big_endian
# FFLAGS = -O3 ${FMACH} # ${DEBUG_OPTIONS} #-convert little_endian


##SGI Origin 2000 options - select only one FFLAGS definition
#----------------------------------------
# MPILIB =  -lmpi
# FMACH =  -64
# DEBUG_OPTIONS = -C -DEBUG:trap_uninitialized=ON -DEBUG:div_check
# FFLAGS =  -g -r8
# FFLAGS = ${FMACH} -g -r8 -extend_source ${DEBUG_OPTIONS}
# FFLAGS = ${FMACH} -O2  -OPT:Olimit=0 -extend_source
# FFLAGS = ${FMACH} -O3 -r8 -OPT:Olimit=0 -extend_source

OBJS = variables.o msga.o firebrand.o radiation.o  \
  windfield.o io.o advectionScheme.o higrad.o \
  compress.o boundary.o \
  definearray.o diffuse.o drag.o \
  fuel.o convection.o fieldUpdate.o filstr.o firetec.o  force.o  \
  largeScalePGF.o metric.o \
  diffusesies.o \
  startFromFile.o rijgradu.o rinit.o rinitfire.o \
  rinitturb.o rmaxmin.o \
  stressrij.o thermodynamics.o tinit.o \
  topo.o turb.o \
  update.o xevariation.o\
  soot.o chemistryTable.o \

firerad: ${OBJS}
	${FC} -o $@  ${FFLAGS} ${OBJS} ${MPILIB}

clean: 
	rm *.o *.mod firerad

# optimized suffix rules
.SUFFIXES: .f .f90 

.f.o:
	${FC} ${FFLAGS} -c $<
.f90.o:
	${FC} ${FFLAGS} -c $<

# Commands below for pink
#.f.o:
#	${FC} ${FFLAGS} -c -I${MPI_ROOT}/include $<
#.f90.o:
#	${FC} ${FFLAGS} -c -I${MPI_ROOT}/include $<

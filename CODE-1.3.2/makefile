
## Commands below for pink - must have loaded modules with 
#              module load lampi intel

## Linux cluster (ie, pink) options - select only one FFLAGS definition
#----------------------------------------
  FC = mpif90   ## compiler
  #MPILIB =  -L${MPI_ROOT}/lib64 -lmpi -I${MPI_ROOT}/include/
  #MPILIB =  -I${MPI_ROOT}/include/
  FFLAGS = -O2 -extend_source 
#  FFLAGS = -g -extend_source 
#   FFLAGS = -g -extend_source -warn all -check bounds -check overflow 

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

OBJS = variables.o msga.o firebrand.o ignite.o radiation.o  \
  nonlocalfunction.o windfield.o io.o advectionScheme.o higrad.o \
  compress.o boundary.o con.o \
  definearray.o diffuse.o draga.o dragb.o dragm.o dragm2.o\
  fuel.o convection.o fieldUpdate.o filstr.o firetec.o  force.o  \
  largeScalePGF.o metryc.o icfmewind.o \
  potflow.o readicfmew.o \
  startFromFile.o rijgradu.o rinit.o rinitfire.o \
  rinitturb.o rmaxmin.o \
  rmaxmin1.o setup.o  stressrij.o tinit.o \
  topo.o turb.o \
  update.o updated.o updated2.o  xevariation.o\

firerad: ${OBJS}
	${FC} -o $@  ${FFLAGS} ${OBJS} ${MPILIB}

clean: 
	rm *.o *.mod *__genmod.f90

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

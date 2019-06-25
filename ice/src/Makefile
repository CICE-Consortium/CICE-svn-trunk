###  CICE Makefile

### Determine the operating system being used, for compilation below
### and for the code's timers.
SYSTEM=$(shell uname -s)

### Set the first definition to "coupled" to compile the code for coupling.
### For stand-alone runs, it can be anything else (eg., "uncoupled").
### The second definition is the type of operating system determined above.
GDEFS = -Duncoupled -D$(SYSTEM)

### Set BINTYPE to be MPI, OPENMP, or HYBRID.  Set it to anything else 
### (eg., "single") for single processor, nonparallel runs.
### HYBRID uses both MPI and OPENMP for parallelization, and currently 
### works only on the IBM (system AIX).
BINTYPE = MPI

### Set the number of processors desired for each grid direction (x, y).
### This is only used for MPI grid decomposition.  NX and NY are ignored
### for OPENMP and single processor runs.
NX=4
NY=2

### SGI Origin 2000
ifeq ($(SYSTEM),IRIX64)
  DEFS = $(GDEFS) -DNPROC_X=1 -DNPROC_Y=1
  FC= f90
  FFLAGS= -r10000 -64 -r8 -O2
  LFLAGS= -r10000 -64 -r8 -O2
  CPP = /lib/cpp
  INCLUDES = -I/usr/local/include
  LIBS     = -lnetcdf -lmpi
#  LIBS     = -lnetcdf -L /usr/local/lib64/r4i4 -lmpi
  ifeq ($(BINTYPE),MPI)
    DEFS = $(GDEFS) -DNPROC_X=$(NX) -DNPROC_Y=$(NY) -D_MPI
  endif
  ifeq ($(BINTYPE),HYBRID)
    FC = **irix64/hybrid not supported**
  endif
  ifeq ($(BINTYPE),OPENMP)
    FFLAGS= -mp -r10000 -64 -r8 -O2
    LFLAGS= -mp -r10000 -64 -r8 -O2
  endif
endif

### IBM RS/6000 SP
ifeq ($(SYSTEM),AIX)
  DEFS = $(GDEFS) -DNPROC_X=1 -DNPROC_Y=1
  FC= xlf 
  FFLAGS= -O3 -qstrict -qrealsize=8 -qfixed -qarch=pwr3
  LFLAGS= -O3 -qstrict -qrealsize=8 -qfixed -qarch=pwr3
  CPP = /usr/ccs/lib/cpp
  INCLUDES = -I/usr/local/include
  LIBS = -L/usr/lpp/ppe.poe/lib -L/usr/local/lib32/r4i4 -lnetcdf -lmass
  ifeq ($(BINTYPE),MPI)
    DEFS = $(GDEFS) -DNPROC_X=$(NX) -DNPROC_Y=$(NY) -D_MPI
    FC= mpxlf
  endif
  ifeq ($(BINTYPE),HYBRID)
    DEFS = $(GDEFS) -DNPROC_X=$(NX) -DNPROC_Y=$(NY) -D_MPI
    FC= mpxlf_r -qsmp -qsmp=noauto -qnosave
    LIBS = -L/usr/lpp/ppe.poe/lib -L/usr/local/lib32/r4i4 -lnetcdf
  endif
  ifeq ($(BINTYPE),OPENMP)
    FC= xlf_r -qsmp -qsmp=noauto -qnosave
    LIBS = -L/usr/lpp/ppe.poe/lib -L/usr/local/lib32/r4i4 -lnetcdf
  endif
endif

###  various sgi debugging options that can be used above
#FFLAGS= -g -r10000 -64 -r8 
#FFLAGS= -g3 -r10000 -64 -r8 -O2 
#FFLAGS= -g -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON -r10000 -64 -r8 
#FFLAGS= -g -DEBUG:div_check=3 -r10000 -64 -r8 
#FFLAGS= -g -DEBUG:subscript_check=ON -r10000 -64 -r8
#FFLAGS= -g -DEBUG:trap_uninitialized=ON -r10000 -64 -r8 

OBJECTS= icegrid.o iceinit.o thermw.o transp.o dyn.o iceout.o ice_misc.o timers.o 

default:
	@echo Please supply a target for make

clean:	
	rm -f ice.o flux.o \
	ice_coupling.o coupling.o \
      $(OBJECTS) \
      *.err *.trace *.h *.f

.SUFFIXES:	.F .H

.F.f:
	$(CPP) $(INCLUDES) $(DEFS) -P $< > $*.f 

.f.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c $*.f

.H.h:	
	$(CPP) $(INCLUDES) $(DEFS) -P $< > $*.h

### cice:  stand-alone
### cicecpl: coupled

cice:	ice.o flux.o $(OBJECTS)
	$(FC) ice.o flux.o $(OBJECTS) \
    $(LFLAGS) $(INCLUDES) $(LIBS) \
    -o cice

cicecpl:	ice_coupling.o coupling.o $(OBJECTS)
	$(FC) ice_coupling.o coupling.o $(OBJECTS) \
    $(LFLAGS) $(INCLUDES) $(LIBS) \
    -o cicecpl

ice.h:  ice.H Makefile

ice.f:  ice.F ice.h timers.h
ice.o:  ice.f

ice_coupling.f:  ice_coupling.F ice.h timers.h ice_mpi.h
ice_coupling.o:  ice_coupling.f

coupling.f:  coupling.F ice.h coupling.h thermw.h
coupling.o:  coupling.f

flux.f:  flux.F ice.h thermw.h
flux.o:  flux.f

icegrid.f:  icegrid.F ice.h ice_mpi.h
icegrid.o:  icegrid.f

iceinit.f:  iceinit.F ice.h thermw.h ice_mpi.h history.h
iceinit.o:  iceinit.f

iceout.f:  iceout.F ice.h history.h
iceout.o:  iceout.f

ice_misc.f:  ice_misc.F ice.h ice_mpi.h
ice_misc.o:  ice_misc.f

thermw.f:  thermw.F ice.h thermw.h
thermw.o:  thermw.f

transp.f:  transp.F ice.h 
transp.o:  transp.f

dyn.f:  dyn.F ice.h 
dyn.o:  dyn.f

timers.f:  timers.F ice.h timers.h
timers.o:  timers.f


















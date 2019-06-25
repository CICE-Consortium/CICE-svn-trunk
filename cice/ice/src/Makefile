#---------------------------------------------------------------
# Makefile:  	for SGI
#-----------------------------------------------------------------
INCLUDES = -I /home/ocean/pwjones/netcdf/include \
           -I /home/ocean/pwjones/mcl/SGI64
LIBS     = -L/home/ocean/pwjones/netcdf/lib -lnetcdf \
           -L/home/ocean/pwjones/mcl/SGI64 -lMCL \
           -L/home/ocean/pwjones/pvm3/lib/SGI64 -lfpvm3 -lgpvm3 -lpvm3

FFLAGS=  -mp -r10000 -64 -r8 -O3 $(INCLUDES)
FC= f90

MPOBJECTS= icegrid.o iceinit.o flux_mp.o therm.o transp.o dyn.o iceout.o

default:
	@echo Please supply a target for make

clean:	
	rm -f ice.o ice_mp.o $(OBJECTS) $(MPOBJECTS) *.err *.trace 

.SUFFIXES:	.F

.F.o:	
	$(FC) $(FFLAGS) -c $<

cice:	ice_mp.o $(MPOBJECTS) 
	$(FC) ice_mp.o $(MPOBJECTS) $(FFLAGS) $(LIBS) -o cice


ice_mp.o:  ice_mp.F ice.h ice_mp.h

icegrid.o:  icegrid.F ice.h ice_mp.h

iceinit.o:  iceinit.F ice.h ice_mp.h

flux_mp.o:  flux_mp.F ice.h ice_mp.h

therm.o:  therm.F ice.h 

transp.o:  transp.F ice.h 

dyn.o:  dyn.F ice.h 

iceout.o:  iceout.F ice.h 












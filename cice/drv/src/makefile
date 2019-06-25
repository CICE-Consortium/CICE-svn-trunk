#=======================================================================
# SCCS module : Makefile.MCL1.0   r3.2+
# SCCS get    : 23 May 1996 / 14:13:25
# SCCS edit   : 23 May 1996 / 14:13:16
# SCCS path   : /crestone/u1/nieman/csm/drv-3.n/SCCS/s.Makefile.MCL1.0
#=======================================================================

SRCDIR = /home/ocean/eclare/release/cice/drv/src

INCLUDES = -I /home/ocean/pwjones/mcl/SGI64

MSGLIB = MCL1.0

EXEC=drv.$(MSGLIB)

FC      = f90 -r8 -r10000 -64 -O3 -mp $(INCLUDES)
CPP     = /lib/cpp 
FLINT   = flint  -gfsw

CPPOPTS = -P -DSPINUP -DRESTORE

SLIBS   = -L /usr/local/lib  

MLIBS   = -L /usr/lib64 \
          -L /home/ocean/pwjones/mcl/SGI64  -lMCL \
	  -L /home/ocean/pwjones/pvm3/lib/SGI64 -lfpvm3 -lgpvm3 -lpvm3

FILES.o = calendar.o          \
          driver.o            \
          main.o              \
          mrgf_ice.o          \
          msg_drv.$(MSGLIB).o \
          msg_ice.$(MSGLIB).o \
          restart.o           

FILES.h = buffer.h            \
	  dims.h              \
          driver.h            \
          grids.h             \
          index.h             \
          msg.$(MSGLIB).h     

#----------------------------------------------------------------------

$(EXEC): $(FILES.o) $(FILES.h)
	$(FC) -o $(EXEC) $(FILES.o) $(MLIBS) $(SLIBS) 

#----------------------------------------------------------------------
# include files
#----------------------------------------------------------------------

buffer.h: $(SRCDIR)/buffer.H
	$(CPP) $(CPPOPTS) $(SRCDIR)/buffer.H >! buffer.h

dims.h: $(SRCDIR)/dims.H
	$(CPP) $(CPPOPTS) $(SRCDIR)/dims.H >! dims.h

driver.h: $(SRCDIR)/driver.H
	$(CPP) $(CPPOPTS) $(SRCDIR)/driver.H >! driver.h

grids.h: $(SRCDIR)/grids.H
	$(CPP) $(CPPOPTS) $(SRCDIR)/grids.H >! grids.h

index.h: $(SRCDIR)/index.H
	$(CPP) $(CPPOPTS) $(SRCDIR)/index.H >! index.h

msg.$(MSGLIB).h: $(SRCDIR)/msg.$(MSGLIB).H
	$(CPP) $(CPPOPTS) $(SRCDIR)/msg.$(MSGLIB).H >! msg.$(MSGLIB).h

#----------------------------------------------------------------------
# program files
#----------------------------------------------------------------------

calendar.o: $(FILES.h) $(SRCDIR)/calendar.F
	$(CPP) $(CPPOPTS) $(SRCDIR)/calendar.F >! calendar.f 
	$(FC) -c calendar.f 

driver.o: $(FILES.h) $(SRCDIR)/driver.F
	$(CPP) $(CPPOPTS) $(SRCDIR)/driver.F >! driver.f 
	$(FC) -c driver.f 

main.o: $(FILES.h) $(SRCDIR)/main.F
	$(CPP) $(CPPOPTS) $(SRCDIR)/main.F >! main.f 
	$(FC) -c main.f 

mrgf_ice.o: $(FILES.h) $(SRCDIR)/mrgf_ice.F
	$(CPP) $(CPPOPTS) $(SRCDIR)/mrgf_ice.F >! mrgf_ice.f 
	$(FC) -c mrgf_ice.f 

msg_drv.$(MSGLIB).o: $(FILES.h) $(SRCDIR)/msg_drv.$(MSGLIB).F
	$(CPP) $(CPPOPTS) $(SRCDIR)/msg_drv.$(MSGLIB).F >! msg_drv.$(MSGLIB).f 
	$(FC) -c msg_drv.$(MSGLIB).f 

msg_ice.$(MSGLIB).o: $(FILES.h) $(SRCDIR)/msg_ice.$(MSGLIB).F
	$(CPP) $(CPPOPTS) $(SRCDIR)/msg_ice.$(MSGLIB).F >! msg_ice.$(MSGLIB).f 
	$(FC) -c msg_ice.$(MSGLIB).f 

restart.o: $(FILES.h) $(SRCDIR)/restart.F
	$(CPP) $(CPPOPTS) $(SRCDIR)/restart.F >! restart.f 
	$(FC) -c restart.f 
clean:
	/bin/rm *.o *.f *.h

#----------------------------------------------------------------------

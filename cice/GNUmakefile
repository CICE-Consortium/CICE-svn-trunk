#
# File:  GNUmakefile
#
#  The bourne shell is safest for gmake
SHELL = /bin/sh

#  Build optimized code by default. To override, add "OPTIMIZE=no" 
OPTIMIZE   = yes

### Determine the operating system being used, for compilation below
### and for the code's timers.
SYSTEM=$(shell uname -s)

#----------------------------------------------------------------------------
#  Define a symbol (TARGETX) for the executable name (csim)
#----------------------------------------------------------------------------

TARGETX = cice
ifeq ($(OPTIMIZE),no)
  #  If building a debug version, append "_db" to the executable name
  TARGETX := $(TARGETX)_db
endif

#----------------------------------------------------------------------------
#  Set valid suffixes.
#----------------------------------------------------------------------------

#  First clean out current list of suffixes, then define them
.SUFFIXES: 
.SUFFIXES: .o .f90 .f .F90 .F

#----------------------------------------------------------------------------
#  Include architecture definitions.
#----------------------------------------------------------------------------

include $(SYSTEM).defs

#  Add directories to gmake search path
VPATH = $(COMPDIR) $(SRCDIR) 

#----------------------------------------------------------------------------
#  Include source and object lists (created by the dependency script)
#  if they exist.  Otherwise force user to create them.
#----------------------------------------------------------------------------

ifeq ($(shell ls Depends),Depends)
  include Objects
  include Sources
else
  exit1:
	@echo "  Please type make Depends to create dependencies"
endif

#----------------------------------------------------------------------------
#  Make target
#----------------------------------------------------------------------------

$(TARGETX): $(OBJS)
	@echo "  Making target '$(TARGETX)'"
	$(LD) -o $(TARGETX) $(LFLAGS) $(OBJS) $(LDLIBS) 

#----------------------------------------------------------------------------
# Include the Depends files which contain all the dependencies
#----------------------------------------------------------------------------

include Depends

Depends: 
	$(PERL) Makedepends.pl $(COMPDIR) $(SRCDIR)

#----------------------------------------------------------------------------
#  Generic compilation rules
#----------------------------------------------------------------------------
 
# Cancel any implicit gmake rules for preprocessing and compiling

%.f90 : %.F90
%.o : %.F90
%.f : %.F
%.o : %.F
%.o : %.f

# Preprocessing Fortran files

$(COMPDIR)/%.f : %.F
	@echo '$(SYSTEM) preprocessing' $<
	$(CPP) $(CPPOPTS) $< >$(COMPDIR)/$*.f

$(COMPDIR)/%.f90 : %.F90
	@echo '$(SYSTEM) preprocessing' $<
	$(CPP) $(CPPOPTS) $< >$(COMPDIR)/$*.f90

# Rule for compiling fortran
# Need to move the object and module files into the compile
# directory after compilation (this could be solved by 
# changing directories into COMPDIR when compiling)

$(COMPDIR)/%.o : $(COMPDIR)/%.f
	@echo $(SYSTEM) Compiling with implicit rule $<
	$(FC) $(FFLAGS) $(FIXEDFLAGS) -c $<
	$(MV) $*.o $(COMPDIR)
	-$(MV) *.$(MOD_SUFFIX) $(COMPDIR)

$(COMPDIR)/%.o : $(COMPDIR)/%.f90
	@echo $(SYSTEM) Compiling with implicit rule $<
	$(FC) $(FFLAGS) $(FREEFLAGS) -c $<
	$(MV) $*.o $(COMPDIR)
	-$(MV) *.$(MOD_SUFFIX) $(COMPDIR)

#----------------------------------------------------------------------------
#  Make clean
#  remove just about everything related to the build process
#----------------------------------------------------------------------------

# This prevents gmake from automatically removing these files.
.PRECIOUS: %.f %.f90 %.o $(COMPDIR)/%.f90 $(COMPDIR)/%.o $(COMPDIR)/%.f

# This forces gmake to do a clean even though the dependency list is blank
.PHONY: clean

clean:
	-$(RM) $(COMPDIR)/*.f
	-$(RM) $(COMPDIR)/*.f90
	-$(RM) $(COMPDIR)/*.o
	-$(RM) $(COMPDIR)/*.$(MOD_SUFFIX)
	-$(RM) $(TARGETX)
	-$(RM) Depends Objects Sources 


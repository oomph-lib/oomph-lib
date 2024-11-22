#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

all: RESLT_no_fix/soln0.dat RESLT_fix/soln0.dat RESLT_sprittles/soln0.dat

RESLT_no_fix/soln0.dat: structured_no_correction parameters.dat
	mkdir -p RESLT_no_fix
	./structured_no_correction

RESLT_fix/soln0.dat: structured_with_correction parameters.dat
	mkdir -p RESLT_fix
	./structured_with_correction

RESLT_sprittles/soln0.dat: sprittles parameters.dat
	mkdir -p RESLT_sprittles
	./sprittles

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_with_correction \
	sprittles

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_with_correction_SOURCES = structured_with_correction.cc
sprittles_SOURCES = sprittles.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
		-lnavier_stokes \
		-laxisym_navier_stokes \
		-lfluid_interface \
		-lsolid \
		-lconstitutive \
		-lgeneric \
		$(EXTERNAL_LIBS) $(FLIBS) \
		-lgcov 

# Add the paranoid and range checking flags
AM_CPPFLAGS += --coverage#-ggdb -O0 #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 -Wno-misleading-indentation 

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

all: RESLT_no_fix/soln0.dat RESLT_fix/soln0.dat RESLT_sprittles/soln0.dat \
	RESLT_internal/soln0.dat RESLT_boundary/soln0.dat \
	RESLT_zero_boundary/soln0.dat RESLT_zero_internal/soln0.dat

RESLT_no_fix/soln0.dat: structured_no_correction parameters.dat
	mkdir -p RESLT_no_fix
	./structured_no_correction

RESLT_fix/soln0.dat: structured_with_correction parameters.dat
	mkdir -p RESLT_fix
	./structured_with_correction

RESLT_sprittles/soln0.dat: sprittles parameters.dat
	mkdir -p RESLT_sprittles
	./sprittles

RESLT_internal/soln0.dat: internal_forcing parameters.dat
	mkdir -p RESLT_internal
	./internal_forcing

RESLT_boundary/soln0.dat: boundary_forcing parameters.dat
	mkdir -p RESLT_boundary
	./boundary_forcing

RESLT_zero_internal/soln0.dat: internal_forcing_zero parameters.dat
	mkdir -p RESLT_zero
	mkdir -p RESLT_internal_zero
	./internal_forcing_zero

RESLT_zero_boundary/soln0.dat: boundary_forcing_zero parameters.dat
	mkdir -p RESLT_zero
	mkdir -p RESLT_boundary_zero
	./boundary_forcing_zero

boundaries: structured_with_correction_regions
	./structured_with_correction_regions

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_with_correction \
	structured_with_correction_regions \
	sprittles \
	internal_forcing \
	boundary_forcing \
	internal_forcing_zero \
	boundary_forcing_zero

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_with_correction_SOURCES = structured_with_correction.cc
structured_with_correction_regions_SOURCES = structured_with_correction_regions.cc
sprittles_SOURCES = sprittles.cc
internal_forcing_SOURCES = internal_forcing.cc
boundary_forcing_SOURCES = boundary_forcing.cc
internal_forcing_zero_SOURCES = internal_forcing_zero.cc
boundary_forcing_zero_SOURCES = boundary_forcing_zero.cc

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
AM_CPPFLAGS += -ggdb -O0 #--coverage #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 -Wno-misleading-indentation 

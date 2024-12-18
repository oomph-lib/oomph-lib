#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

paraview: RESLT_fix/soln0.vtu structured_with_correction parameters.dat
	paraview RESLT_fix/soln0.vtu 

paraview_region: RESLT_fix_region/soln0.vtu structured_with_correction_region parameters.dat
	paraview RESLT_fix_region/soln0.vtu 

%.vtu: %.dat
	../../../bin/oomph-convert.py $<

parameters.dat: 
	cp default_parameters.dat parameters.dat

RESLT_no_fix/soln0.dat: structured_no_correction parameters.dat
	mkdir -p RESLT_no_fix
	./structured_no_correction

RESLT_no_fix_region/soln0.dat: structured_no_correction_region parameters.dat
	mkdir -p RESLT_no_fix_region
	./structured_no_correction_region

RESLT_fix/soln0.dat: structured_with_correction parameters.dat
	mkdir -p RESLT_fix
	./structured_with_correction

RESLT_fix_region/soln0.dat: structured_with_correction_region parameters.dat
	mkdir -p RESLT_no_fix_region
	mkdir -p RESLT_fix_region
	./structured_with_correction_region

RESLT_sprittles/soln0.dat: sprittles parameters.dat
	mkdir -p RESLT_sprittles
	./sprittles

RESLT_sprittles_region/soln0.dat: sprittles_region parameters.dat
	mkdir -p RESLT_sprittles_region
	./sprittles_region

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_no_correction_region \
	structured_with_correction \
	structured_with_correction_region \
	unstructured_no_correction \
	unstructured_with_correction \
	sprittles \
	sprittles_region

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_no_correction_region_SOURCES = structured_no_correction_region.cc
structured_with_correction_SOURCES = structured_with_correction.cc
structured_with_correction_region_SOURCES = structured_with_correction_region.cc
unstructured_no_correction_SOURCES = unstructured_no_correction.cc
unstructured_with_correction_SOURCES = unstructured_with_correction.cc
sprittles_SOURCES = sprittles.cc
sprittles_region_SOURCES = sprittles_region.cc

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
AM_CPPFLAGS += #-ggdb -O0 -Wno-misleading-indentation #--coverage #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 

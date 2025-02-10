#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

FOLDER ?= ./RESLT_axi_sprittles_region

.PHONY: run

run: $(FOLDER)/soln0.vtu
	paraview $(FOLDER)/soln0.vtu

%.vtu: %.dat
	../../../bin/oomph-convert.py $<

RESLT_axi_no_fix/soln0.dat: structured_no_correction parameters.dat
	mkdir -p RESLT_axi_no_fix
	./structured_no_correction

RESLT_axi_no_fix_region/soln0.dat: structured_no_correction_region parameters.dat
	mkdir -p RESLT_axi_no_fix_region
	./structured_no_correction_region

RESLT_axi_fix/soln0.dat: structured_with_correction parameters.dat
	mkdir -p RESLT_axi_fix
	./structured_with_correction

RESLT_axi_fix_region/soln0.dat: structured_with_correction_region parameters.dat
	mkdir -p RESLT_axi_fix_region
	./structured_with_correction_region

RESLT_axi_sprittles_region/soln0.dat: sprittles_region parameters.dat
	mkdir -p RESLT_axi_sprittles_region
	./sprittles_region

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_no_correction_region \
	structured_with_correction \
	structured_with_correction_region \
	unstructured_no_correction \
	unstructured_with_correction \
	sprittles_region \
	test_parameters

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc parameter_struct.cc
structured_no_correction_region_SOURCES = structured_no_correction_region.cc parameter_struct.cc
structured_with_correction_SOURCES = structured_with_correction.cc parameter_struct.cc
structured_with_correction_region_SOURCES = structured_with_correction_region.cc parameter_struct.cc
unstructured_no_correction_SOURCES = unstructured_no_correction.cc parameter_struct.cc
unstructured_with_correction_SOURCES = unstructured_with_correction.cc parameter_struct.cc
sprittles_region_SOURCES = sprittles_region.cc parameter_struct.cc
test_parameters_SOURCES = test_parameters.cc parameter_struct.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
		-lnavier_stokes \
		-laxisym_navier_stokes \
		-lfluid_interface \
		-lsolid \
		-lconstitutive \
		-lgeneric \
		$(EXTERNAL_LIBS) $(FLIBS)
# -lgcov 

# Add the paranoid and range checking flags
AM_CPPFLAGS += -O0 -Wno-misleading-indentation #--coverage -ggdb -DPARANOID  -DRANGE_CHECKING 

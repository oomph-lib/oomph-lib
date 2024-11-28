#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

all: RESLT_no_fix/soln0.dat RESLT_fix/soln0.dat \
	RESLT_fix_region/soln0.dat

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

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_no_correction_region \
	structured_with_correction \
	structured_with_correction_region

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_no_correction_region_SOURCES = structured_no_correction_region.cc
structured_with_correction_SOURCES = structured_with_correction.cc
structured_with_correction_region_SOURCES = structured_with_correction_region.cc

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

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
	structured_no_correction \
	structured_with_correction 

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_with_correction_SOURCES = structured_with_correction.cc

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
#-lgcov 

# Add the paranoid and range checking flags
AM_CPPFLAGS += -ggdb -O0 #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 -Wno-misleading-indentation --coverage

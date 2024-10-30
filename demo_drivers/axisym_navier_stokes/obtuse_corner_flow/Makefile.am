#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

all: obtuse_corner_flow run_tests 

# Name of executables
check_PROGRAMS= \
	obtuse_corner_flow \
	run_tests

# Sources for second executable
obtuse_corner_flow_SOURCES = obtuse_corner_flow.cc
run_tests_SOURCES = run_tests.cc

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
AM_CPPFLAGS += #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 -Wno-misleading-indentation --coverage

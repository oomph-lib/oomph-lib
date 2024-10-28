#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
	linear_stability_run \
	neutral_stability_run \
	run_tests

# Sources for second executable
linear_sources = linearised_axisymmetric_fluid_interface_elements.cc \
				 full_linearised_axisym_navier_stokes_elements.cc \
				 Tfull_linearised_axisym_navier_stokes_elements.cc \
				 decomposed_pvd_elements.cc \
				 decomposed_linear_elasticity_elements.cc

linear_stability_run_SOURCES = linear_stability_run.cc $(linear_sources)
neutral_stability_run_SOURCES = neutral_stability_run.cc $(linear_sources)
run_tests_SOURCES = run_tests.cc $(linear_sources)

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ -llinear_elasticity -llinearised_axisym_navier_stokes -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric \
									$(EXTERNAL_LIBS) $(FLIBS)# -lboost_unit_test_framework -lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS +=  #-DPARANOID  -DRANGE_CHECKING #-Wfatal-errors # --coverage -O0

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
	axi_dynamic_cap \
	create_parameter_files \
	continuation_run \
	custom_continuation_run \
	custom_ca_continuation_run \
	arc_continuation_run \
	bifurcation_tracking_run \
	continue_solution_limit_run \
	height_continuation_run \
	ca_height_continuation_run \
	obtuse_bo_height_continuation_run \
	obtuse_ca_height_continuation_run \
	run_tests

axi_dynamic_cap_SOURCES = axi_dynamic_cap.cc 
create_parameter_files_SOURCES = create_parameter_files.cc
continuation_run_SOURCES = continuation_run.cc
custom_continuation_run_SOURCES = custom_continuation_run.cc
custom_ca_continuation_run_SOURCES = custom_ca_continuation_run.cc
arc_continuation_run_SOURCES = arc_continuation_run.cc
bifurcation_tracking_run_SOURCES = bifurcation_tracking_run.cc
continue_solution_limit_run_SOURCES = continue_solution_limit_run.cc
height_continuation_run_SOURCES = height_continuation_run.cc
ca_height_continuation_run_SOURCES = ca_height_continuation_run.cc
obtuse_bo_height_continuation_run_SOURCES = obtuse_bo_height_continuation_run.cc
obtuse_ca_height_continuation_run_SOURCES = obtuse_ca_height_continuation_run.cc
run_tests_SOURCES = run_tests.cc $(linear_sources)


LDADD = -L@libdir@ -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric \
									$(EXTERNAL_LIBS) $(FLIBS) -lboost_unit_test_framework -lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS += -DPARANOID -DRANGE_CHECKING --coverage -O0 # -Wfatal-errors 

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

all: unsteady_run \
	steady_run \
	obtuse_steady_run \
	axisym_stability_run \
	linear_steady_run \
	linear_unsteady_run \
	linear_stability_run \
	obtuse_linear_stability_run \
	create_parameter_files \
	neutral_stability_run \
	obtuse_neutral_stability_run \
	neutral_stability_continuation_run \
	axisym_neutral_stability_run \
	perturb_z_run \
	continuation_run \
	adapt_run \
	custom_continuation_run \
	custom_ca_continuation_run \
	arc_continuation_run \
	bifurcation_tracking_run \
	debug_jacobian_run \
	debug_linear_jacobian_run \
	perturb_by_eigensolution_run \
	continue_neutral_stability_run \
	continue_solution_limit_run \
	arc_continuation_upto_limit_point_run \
	height_continuation_run \
	ca_height_continuation_run \
	obtuse_bo_height_continuation_run \
	obtuse_ca_height_continuation_run \
	test_volume_constraint \
	run_tests \
	run_singular_tests

# Name of executables
check_PROGRAMS= \
	axi_dynamic_cap \
	create_parameter_files \
	continuation_run \
	custom_continuation_run \
	custom_ca_continuation_run \
	bifurcation_tracking_run \
	continue_solution_limit_run \
	obtuse_bo_height_continuation_run \
	obtuse_ca_height_continuation_run \
	run_tests

axi_dynamic_cap_SOURCES = axi_dynamic_cap.cc 
create_parameter_files_SOURCES = create_parameter_files.cc
continuation_run_SOURCES = continuation_run.cc
custom_continuation_run_SOURCES = custom_continuation_run.cc
custom_ca_continuation_run_SOURCES = custom_ca_continuation_run.cc
bifurcation_tracking_run_SOURCES = bifurcation_tracking_run.cc
continue_solution_limit_run_SOURCES = continue_solution_limit_run.cc
obtuse_bo_height_continuation_run_SOURCES = obtuse_bo_height_continuation_run.cc
obtuse_ca_height_continuation_run_SOURCES = obtuse_ca_height_continuation_run.cc
run_tests_SOURCES = run_tests.cc $(linear_sources)


LDADD = -L@libdir@ -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric \
									$(EXTERNAL_LIBS) $(FLIBS) #-lboost_unit_test_framework -lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS += #-DPARANOID -DRANGE_CHECKING --coverage -O0 # -Wfatal-errors 

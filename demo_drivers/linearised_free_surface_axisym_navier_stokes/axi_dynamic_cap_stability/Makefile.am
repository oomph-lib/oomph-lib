#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
	unsteady_run \
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

# Sources for second executable
linear_sources = linearised_axisymmetric_fluid_interface_elements.cc \
				 full_linearised_axisym_navier_stokes_elements.cc \
				 Tfull_linearised_axisym_navier_stokes_elements.cc \
				 decomposed_pvd_elements.cc \
				 decomposed_linear_elasticity_elements.cc

unsteady_run_SOURCES = unsteady_run.cc 
steady_run_SOURCES = steady_run.cc 
obtuse_steady_run_SOURCES = obtuse_steady_run.cc 
axisym_stability_run_SOURCES = axisym_stability_run.cc 
linear_steady_run_SOURCES = linear_steady_run.cc $(linear_sources)
linear_unsteady_run_SOURCES = linear_unsteady_run.cc $(linear_sources)
linear_stability_run_SOURCES = linear_stability_run.cc $(linear_sources)
obtuse_linear_stability_run_SOURCES = obtuse_linear_stability_run.cc $(linear_sources)
create_parameter_files_SOURCES = create_parameter_files.cc
neutral_stability_run_SOURCES = neutral_stability_run.cc $(linear_sources)
obtuse_neutral_stability_run_SOURCES = obtuse_neutral_stability_run.cc $(linear_sources)
neutral_stability_continuation_run_SOURCES = neutral_stability_continuation_run.cc $(linear_sources)
axisym_neutral_stability_run_SOURCES = axisym_neutral_stability_run.cc $(linear_sources)
perturb_z_run_SOURCES = perturb_z_run.cc 
continuation_run_SOURCES = continuation_run.cc
adapt_run_SOURCES = adapt_run.cc 
custom_continuation_run_SOURCES = custom_continuation_run.cc
custom_ca_continuation_run_SOURCES = custom_ca_continuation_run.cc
arc_continuation_run_SOURCES = arc_continuation_run.cc
bifurcation_tracking_run_SOURCES = bifurcation_tracking_run.cc
debug_jacobian_run_SOURCES = debug_jacobian_run.cc
debug_linear_jacobian_run_SOURCES = debug_linear_jacobian_run.cc $(linear_sources)
perturb_by_eigensolution_run_SOURCES = perturb_by_eigensolution_run.cc
continue_neutral_stability_run_SOURCES = continue_neutral_stability_run.cc $(linear_sources)
continue_solution_limit_run_SOURCES = continue_solution_limit_run.cc
arc_continuation_upto_limit_point_run_SOURCES = arc_continuation_upto_limit_point_run.cc
height_continuation_run_SOURCES = height_continuation_run.cc
ca_height_continuation_run_SOURCES = ca_height_continuation_run.cc
obtuse_bo_height_continuation_run_SOURCES = obtuse_bo_height_continuation_run.cc
obtuse_ca_height_continuation_run_SOURCES = obtuse_ca_height_continuation_run.cc
test_volume_constraint_SOURCES = test_volume_constraint.cc $(linear_sources)
run_tests_SOURCES = run_tests.cc $(linear_sources)
run_singular_tests_SOURCES = run_singular_tests.cc $(linear_sources)

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
linear_libraries = -L@libdir@ -llinear_elasticity -laxisym_navier_stokes \
				   					-llinearised_axisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric \
									$(EXTERNAL_LIBS) $(FLIBS) -lgcov

LDADD = -L@libdir@ -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric \
									$(EXTERNAL_LIBS) $(FLIBS) -lboost_unit_test_framework -lgcov
linear_steady_run_LDADD = $(linear_libraries)
linear_unsteady_run_LDADD = $(linear_libraries)
linear_stability_run_LDADD = $(linear_libraries)
obtuse_linear_stability_run_LDADD = $(linear_libraries)
neutral_stability_run_LDADD = $(linear_libraries)
obtuse_neutral_stability_run_LDADD = $(linear_libraries)
neutral_stability_continuation_run_LDADD = $(linear_libraries)
axisym_neutral_stability_run_LDADD = $(linear_libraries)
perturb_z_run_LDADD = $(linear_libraries)
debug_linear_jacobian_run_LDADD = $(linear_libraries)
continue_neutral_stability_run_LDADD = $(linear_libraries)
test_volume_constraint_LDADD = $(linear_libraries)
create_parameter_files_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS) -lgcov
run_tests_LDADD = $(linear_libraries)
run_singular_tests_LDADD = $(linear_libraries)

# Add the paranoid and range checking flags
AM_CPPFLAGS +=  #-DPARANOID  -DRANGE_CHECKING #-Wfatal-errors # --coverage -O0

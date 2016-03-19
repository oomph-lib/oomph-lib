#Include instructions common to all directories with validata
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
check_PROGRAMS= \
unstructured_adaptive_mesh_two_outer_boundaries \
unstructured_adaptive_mesh_two_outer_boundaries_crossed \
two_d_parallel_unstructured_adaptive_poisson

# THE EXECUTABLE:
#-----------------------

#----------------------------------------------------------------------
# Sources for  executable
unstructured_adaptive_mesh_two_outer_boundaries_SOURCES = unstructured_adaptive_mesh_two_outer_boundaries.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_adaptive_mesh_two_outer_boundaries_LDADD = \
-L@libdir@ -lpoisson -lrigid_body \
	-lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------
# Sources for  executable
unstructured_adaptive_mesh_two_outer_boundaries_crossed_SOURCES = \
	unstructured_adaptive_mesh_two_outer_boundaries_crossed.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_adaptive_mesh_two_outer_boundaries_crossed_LDADD = \
	-L@libdir@ -lpoisson -lrigid_body \
	-lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for  executable
two_d_parallel_unstructured_adaptive_poisson_SOURCES = two_d_parallel_unstructured_adaptive_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_parallel_unstructured_adaptive_poisson_LDADD = \
	-L@libdir@ -lpoisson -lrigid_body \
	-lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS = -I@includedir@

EXTRA_DIST += validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_SQUARE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_SQUARE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_SQUARE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_SQUARE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_SQUARE/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_SQUARE/input_distribution_3.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE/input_distribution_3.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_HALF_CIRCLE_INTERNAL_BOUNDARIES/input_distribution_3.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_COMPLEX_HOLES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP2_COMPLEX_HOLES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_COMPLEX_HOLES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_COMPLEX_HOLES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP3_COMPLEX_HOLES/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES/input_distribution_0.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES/input_distribution_1.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES/input_distribution_2.dat\
validata/TWO_D_PARALLEL_UNSTRUCTURED_ADAPTIVE_POISSON_NP4_COMPLEX_HOLES/input_distribution_3.dat


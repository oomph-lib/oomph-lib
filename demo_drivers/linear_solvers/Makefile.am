#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
adv_diff_iterative_linear_solver_tester \
two_d_poisson_dense_lu \
two_d_poisson_cg \
driven_cavity \
direct_solver_test \
two_d_multi_poisson \
two_d_linear_elasticity_with_simple_block_diagonal_preconditioner



# Do we have the HYPRE library?
#==============================
if OOMPH_HAS_HYPRE
#Don't add this if we have MPI 
#(HYPRE can't handle serial execution of MPI-enabled code)
#=========================================================
if !WANT_MPI
check_PROGRAMS += HypreSolver_test 
endif
endif


# Do we have the TRILINOS library?
#=================================
if OOMPH_HAS_TRILINOS
#Don't add this if we have MPI 
#(TRILINOS can't handle serial execution of MPI-enabled code)
#=========================================================
if !WANT_MPI
check_PROGRAMS += TrilinosSolver_test
endif
endif


#-------------------------------------------------------------------

# Sources for executable
#driven_cavity_with_simple_lsc_preconditioner_SOURCES = \
#driven_cavity_with_simple_lsc_preconditioner.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
#driven_cavity_with_simple_lsc_preconditioner_LDADD = \
#                -L@libdir@ -lnavier_stokes  \
#                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------------------------

# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Sources for executable
two_d_linear_elasticity_with_simple_block_diagonal_preconditioner_SOURCES = \
two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc simple_block_preconditioners.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_linear_elasticity_with_simple_block_diagonal_preconditioner_LDADD = \
                -L@libdir@ -llinear_elasticity  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


# Sources for executable
adv_diff_iterative_linear_solver_tester_SOURCES =  adv_diff_iterative_linear_solver_tester.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adv_diff_iterative_linear_solver_tester_LDADD = \
                -L@libdir@ -ladvection_diffusion  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


# Sources for executable
two_d_poisson_dense_lu_SOURCES = two_d_poisson_dense_lu.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_dense_lu_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Sources for executable
two_d_poisson_cg_SOURCES = two_d_poisson_cg.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_cg_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



#----------------------------------------------------------------------


# Sources for executable
direct_solver_test_SOURCES = direct_solver_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
direct_solver_test_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Do we have the TRILINOS library?
#=================================
if OOMPH_HAS_TRILINOS

# Sources the executable depends on:
TrilinosSolver_test_SOURCES = TrilinosSolver_test.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" library
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
TrilinosSolver_test_LDADD = -L@libdir@ -lpoisson -lgeneric \
                       $(EXTERNAL_LIBS) \
                       $(FLIBS)

endif


#----------------------------------------------------------------------


# Do we have the HYPRE library?
#===============================
if OOMPH_HAS_HYPRE

# Sources for executable
HypreSolver_test_SOURCES = HypreSolver_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
HypreSolver_test_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric  \
                 $(EXTERNAL_LIBS) \
                 $(FLIBS)

endif

#----------------------------------------------------------------------



# Sources the executable depends on:
two_d_multi_poisson_SOURCES = two_d_multi_poisson.cc \
	multi_poisson_block_preconditioners.h \
	multi_poisson_elements.h	

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" library
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
two_d_multi_poisson_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I. The slightly convoluted
# reference to generic is required to allow the multi-poisson driver code
# to include the multi-poisson header file in the serial demo driver
# directory to find the generic header files.
AM_CPPFLAGS +=  -I@includedir@   -I@includedir@/generic



EXTRA_DIST+=time_driven_cavity.bash driven_cavity_with_simple_lsc_preconditioner.cc

# Additional cleanup commands -- just in case there are serial versions
# of the Trilinos and Hypre test codes hanging around...
CLEANFILES=HypreSolver_test TrilinosSolver_test 

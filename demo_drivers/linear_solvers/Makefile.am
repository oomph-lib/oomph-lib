#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
adv_diff_iterative_linear_solver_tester \
two_d_poisson_dense_lu \
two_d_poisson_cg \
driven_cavity \
direct_solver_test



# Do we have the HYPRE library?
#==============================
if HAVE_HYPRE
#Don't add this if we have MPI 
#(HYPRE can't handle serial execution of MPI-enabled code)
#=========================================================
if !WANT_MPI
check_PROGRAMS += HypreSolver_test 
endif
endif


# Do we have the TRILINOS library?
#=================================
if HAVE_TRILINOS
#Don't add this if we have MPI 
#(TRILINOS can't handle serial execution of MPI-enabled code)
#=========================================================
if !WANT_MPI
check_PROGRAMS += TrilinosSolver_test
endif
endif

#----------------------------------------------------------------------


# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
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
if HAVE_TRILINOS

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
if HAVE_HYPRE

#----------------------------------------------------------------------

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

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  



EXTRA_DIST+=time_driven_cavity.bash

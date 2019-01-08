#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Name of executables: 
check_PROGRAMS=driven_cavity driven_cavity_without_singularity_subtracted_off \
circular_couette circular_couette_without_pseudo_singularity_subtracted_off \
adaptive_driven_cavity

#---------------------------------------------------------------------------

# Sources for executable
adaptive_driven_cavity_SOURCES = adaptive_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------------

# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc navier_stokes_elements_with_singularity.h finite_re_perturbation.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
driven_cavity_without_singularity_subtracted_off_SOURCES = driven_cavity.cc navier_stokes_elements_with_singularity.h finite_re_perturbation.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_without_singularity_subtracted_off_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

driven_cavity_without_singularity_subtracted_off_CXXFLAGS = -DDONT_USE_SINGULARITY


#---------------------------------------------------------------------------

# Sources for executable
circular_couette_SOURCES = driven_cavity.cc navier_stokes_elements_with_singularity.h 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_couette_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_couette_CXXFLAGS = -DCOUETTE


#---------------------------------------------------------------------------

# Sources for executable
circular_couette_without_pseudo_singularity_subtracted_off_SOURCES = driven_cavity.cc navier_stokes_elements_with_singularity.h finite_re_perturbation.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_couette_without_pseudo_singularity_subtracted_off_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_couette_without_pseudo_singularity_subtracted_off_CXXFLAGS = -DDONT_USE_SINGULARITY -DCOUETTE


#---------------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
#AM_CPPFLAGS += -I@includedir@  

EXTRA_DIST+=deriv_blend.map  run_adaptive_driven_cavity.bash  run_convergence_test_circular_couette.bash  run_convergence_test_driven_cavity.bash convergence_and_cpu_couette.lay compare_driven_cavity_different_re.pvsm  compare_driven_cavity.pvsm  validate_grad_u_sing.pvsm run_re_dependence.bash subtract_and_scale.bash finite_re_vs_asymptotics_along_lines.lay 



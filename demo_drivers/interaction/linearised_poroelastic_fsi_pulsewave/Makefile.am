#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=linearised_poroelastic_fsi_pulsewave \
chris

#-----------------------------------------------

# Sources for executable
linearised_poroelastic_fsi_pulsewave_SOURCES = linearised_poroelastic_fsi_pulsewave.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
linearised_poroelastic_fsi_pulsewave_LDADD = -L@libdir@ -laxisym_navier_stokes \
                                  -laxisym_poroelasticity -lgeneric \
                                  $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------

# Sources for executable
chris_SOURCES = chris.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
chris_LDADD = -L@libdir@ -laxisym_navier_stokes \
                                  -laxisym_poroelasticity -lgeneric \
                                  $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------

EXTRA_DIST += \
run_it.bash \
parameters.map \
compare_two_fluid_runs.pvsm  compare_two_poro_runs.pvsm  \
animate.py plot_it.bash animate_poro.pvsm  animate.pvsm \
animate_chris.py chris_parameters.map chris_poro_and_nst.pvsm \
chris.pvsm plot_chris.bash run_chris.bash



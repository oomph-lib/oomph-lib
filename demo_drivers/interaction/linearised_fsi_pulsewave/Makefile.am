#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=linearised_fsi_pulsewave

#-----------------------------------------------

# Sources for zeroth executable
linearised_fsi_pulsewave_SOURCES = linearised_fsi_pulsewave.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
linearised_fsi_pulsewave_LDADD = -L@libdir@ -laxisym_navier_stokes \
                                 -laxisym_linear_elasticity -lgeneric \
                                 $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------

EXTRA_DIST += plot_it.bash pulse_wave_comparison_with_wall.pvsm	\
              pulse_wave_comparison.pvsm 

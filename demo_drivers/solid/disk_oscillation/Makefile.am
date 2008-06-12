#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=disk_oscillation

# Sources for executable
disk_oscillation_SOURCES = disk_oscillation.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
disk_oscillation_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
-loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=disk_oscillation 
#disk_oscillation_demo

# Sources for executable
disk_oscillation_SOURCES = disk_oscillation.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
disk_oscillation_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
#disk_oscillation_demo_SOURCES = disk_oscillation_demo.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
#disk_oscillation_demo_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
# $(EXTERNAL_LIBS) $(FLIBS)

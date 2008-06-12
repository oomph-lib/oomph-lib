#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
thermo 

# Sources for first executable
thermo_SOURCES = thermo.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
thermo_LDADD = -L@libdir@ -lunsteady_heat \
 -lsolid -lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

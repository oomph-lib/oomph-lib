#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=adaptive_interface

# Sources for executable
adaptive_interface_SOURCES = adaptive_interface.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_interface_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lsolid -lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

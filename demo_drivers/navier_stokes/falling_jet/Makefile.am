#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=falling_jet

# Sources for executable
falling_jet_SOURCES = falling_jet.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
falling_jet_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lconstitutive -lsolid -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

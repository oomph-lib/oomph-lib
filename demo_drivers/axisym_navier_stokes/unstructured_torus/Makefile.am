#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_torus pressure_driven_torus

# Sources for executable
unstructured_torus_SOURCES = unstructured_torus.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_torus_LDADD = -L@libdir@ -laxisym_navier_stokes \
                           -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
pressure_driven_torus_SOURCES = pressure_driven_torus.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pressure_driven_torus_LDADD = -L@libdir@ -lsolid -lconstitutive \
                            -laxisym_navier_stokes \
                           -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

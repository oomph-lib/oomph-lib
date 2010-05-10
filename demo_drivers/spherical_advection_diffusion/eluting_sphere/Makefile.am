#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=eluting_sphere eluting_sphere_adapt

# Sources for executable
eluting_sphere_SOURCES = eluting_sphere.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
eluting_sphere_LDADD = -L@libdir@ -lspherical_advection_diffusion -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
eluting_sphere_adapt_SOURCES = eluting_sphere_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
eluting_sphere_adapt_LDADD = -L@libdir@ -lspherical_advection_diffusion \
                            -lgeneric \
                            $(EXTERNAL_LIBS) $(FLIBS)

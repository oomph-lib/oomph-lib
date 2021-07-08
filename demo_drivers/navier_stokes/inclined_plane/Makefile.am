#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=inclined_plane inclined_plane2

# Sources for executable
inclined_plane_SOURCES = inclined_plane.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
inclined_plane_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lsolid \
                       -lconstitutive  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
inclined_plane2_SOURCES = inclined_plane.cc

#Set compiler flag
inclined_plane2_CXXFLAGS = -DCR_ELEMENT

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
inclined_plane2_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lsolid \
                       -lconstitutive  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=vmtk_solid  

#-------------------------------------------------------------------


# Sources for executable
vmtk_solid_SOURCES = vmtk_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
vmtk_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

EXTRA_DIST += solid_iliac.1.ele  solid_iliac.1.face  solid_iliac.1.node  solid_iliac.poly



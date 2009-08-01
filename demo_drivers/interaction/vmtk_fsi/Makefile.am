#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=vmtk_fsi


#---------------------------------------------------------------------

# Sources for first executable
vmtk_fsi_SOURCES = vmtk_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
vmtk_fsi_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#vmtk_fsi_CXXFLAGS= -g -DPARANOID -DRANGE_CHECKING

#---------------------------------------------------------------------


EXTRA_DIST += solid_iliac.1.ele  solid_iliac.1.face  solid_iliac.1.node \
 solid_iliac.poly fluid_iliac.1.ele  fluid_iliac.1.face  fluid_iliac.1.node \
 fluid_iliac.poly iliac.xda





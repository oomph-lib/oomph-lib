#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=mesh_from_gmsh_solid singleElement

#-----------------------------------------------------------------------

# Sources for executable
mesh_from_gmsh_solid_SOURCES = mesh_from_gmsh_solid.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
mesh_from_gmsh_solid_LDADD = -L@libdir@ -lsolid -lconstitutive \
                             -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
# AM_CPPFLAGS += -I@includedir@  


#-----------------------------------------------------------------------

# Sources for executable
singleElement_SOURCES = singleElement.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
singleElement_LDADD = -L@libdir@ -lsolid -lconstitutive \
                      -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
# AM_CPPFLAGS += -I@includedir@  

#-----------------------------------------------------------------------

EXTRA_DIST = cantiliver.msh cube.msh


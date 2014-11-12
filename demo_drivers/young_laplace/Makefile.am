include $(top_srcdir)/config/makefile_templates/demo_drivers


# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
check_PROGRAMS= young_laplace refineable_young_laplace spherical_cap_in_cylinder barrel refineable_t_junction


#---------------------------------------------------------------------------

# Local sources that the code depends on:
refineable_t_junction_SOURCES = refineable_t_junction.cc

# Required libraries:
refineable_t_junction_LDADD = -L@libdir@ -lyoung_laplace -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

refineable_t_junction_CXXFLAGS=-g

#---------------------------------------------------------------------------

# Local sources that the code depends on:
barrel_SOURCES = barrel.cc

# Required libraries:
barrel_LDADD = -L@libdir@ -lyoung_laplace -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------------

# Local sources that the code depends on:
spherical_cap_in_cylinder_SOURCES = spherical_cap_in_cylinder.cc common_young_laplace_stuff.h

# Required libraries:
spherical_cap_in_cylinder_LDADD = -L@libdir@ -lyoung_laplace -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------------

# Local sources that the code depends on:
young_laplace_SOURCES = young_laplace.cc common_young_laplace_stuff.h

# Required libraries: 
young_laplace_LDADD = -L@libdir@ -lyoung_laplace -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Local sources that the code depends on:
refineable_young_laplace_SOURCES = refineable_young_laplace.cc common_young_laplace_stuff.h

# Required libraries:
refineable_young_laplace_LDADD = -L@libdir@ -lyoung_laplace -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  



#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=simple_shear refineable_simple_shear

# Sources for executable
simple_shear_SOURCES = simple_shear.cc precompiled_mesh.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
simple_shear_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
					$(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
refineable_simple_shear_SOURCES = refineable_simple_shear.cc precompiled_mesh.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_simple_shear_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
					$(EXTERNAL_LIBS) $(FLIBS)

#Automake doesn't have an easy way to make without optimisation
#This is the automake generated option without CXX and CPP flags
precompiled_mesh.o : precompiled_mesh.cc
	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

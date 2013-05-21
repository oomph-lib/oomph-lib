#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=pretend_melt spring_contact melt

#-------------------------------------------------

# Sources for executable
spring_contact_SOURCES = spring_contact.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spring_contact_LDADD = -L@libdir@ -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
pretend_melt_SOURCES = pretend_melt.cc unsteady_heat_flux_pseudo_melt_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pretend_melt_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
melt_SOURCES = melt.cc heat_transfer_and_melt_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
melt_LDADD = -L@libdir@ -lconstitutive -lsolid -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------




EXTRA_DIST += pseudo_melt.pvsm  plot_pseudo_melt.bash \
crop.bash plot_spring_contact_landscape.bash spring_contact_landscape.pvsm \
run_spring_contact.bash exact.map melt.pvsm plot_melt.bash re-enumerate.bash \
proper_melt.pvsm






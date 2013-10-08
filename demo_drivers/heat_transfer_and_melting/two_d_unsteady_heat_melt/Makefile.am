#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=pretend_melt spring_contact melt stefan_boltzmann stefan_boltzmann_melt

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

# Sources for executable
stefan_boltzmann_SOURCES = stefan_boltzmann.cc temporary_stefan_boltzmann_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
stefan_boltzmann_LDADD = -L@libdir@ -lconstitutive -lsolid -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
stefan_boltzmann_melt_SOURCES = stefan_boltzmann_melt.cc heat_transfer_and_melt_elements.h temporary_stefan_boltzmann_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
stefan_boltzmann_melt_LDADD = -L@libdir@ -lconstitutive -lsolid -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------




EXTRA_DIST += \
crop.bash  plot_melt.bash  plot_pseudo_melt.bash  plot_sb.bash	plot_spring_contact_landscape.bash  re-enumerate.bash  run_spring_contact.bash

EXTRA_DIST+= \
melt.pvsm  proper_melt.pvsm  pseudo_melt.pvsm  spring_contact_landscape.pvsm  stefan_boltzmann.pvsm  stefan_boltzmann_rays.pvsm stefan_boltzmann_validation.map bins.pvsm

EXTRA_DIST+= \
exact.map animate.py






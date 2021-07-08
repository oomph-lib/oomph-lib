#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=pretend_melt spring_contact melt two_d_unsteady_heat \
stefan_boltzmann stefan_boltzmann_melt \
solid_contact \
circular_boulder_melt \
circular_boulder_solar_radiation \
solid_contact_with_gravity_structured \
solid_contact_with_gravity_unstructured \
linear_solid_contact_with_gravity_structured \
linear_solid_contact_with_gravity_unstructured \
heated_linear_solid_contact_with_gravity \
spring_contact2

#-------------------------------------------------

# Sources for executable
solid_contact_SOURCES = solid_contact.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
solid_contact_LDADD = -L@libdir@  -lsolid -lconstitutive -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------

# Sources for executable
heated_linear_solid_contact_with_gravity_SOURCES = heated_linear_solid_contact_with_gravity.cc contact_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
heated_linear_solid_contact_with_gravity_LDADD = -L@libdir@ -lunsteady_heat -llinear_elasticity \
                                -lconstitutive \
                                -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
linear_solid_contact_with_gravity_unstructured_SOURCES = linear_solid_contact_with_gravity.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
linear_solid_contact_with_gravity_unstructured_LDADD = -L@libdir@ -llinear_elasticity \
                                -lconstitutive \
                                -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------

# Sources for executable
linear_solid_contact_with_gravity_structured_SOURCES = linear_solid_contact_with_gravity.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
linear_solid_contact_with_gravity_structured_LDADD = -L@libdir@ -llinear_elasticity \
                                -lconstitutive \
                                -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

linear_solid_contact_with_gravity_structured_CXXFLAGS=-DSTRUCTURED_MESH

#-------------------------------------------------

# Sources for executable
solid_contact_with_gravity_structured_SOURCES = solid_contact_with_gravity.cc \
contact_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
solid_contact_with_gravity_structured_LDADD = -L@libdir@ -lconstitutive \
                                -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

solid_contact_with_gravity_structured_CXXFLAGS=-DSTRUCTURED_MESH

#-------------------------------------------------

# Sources for executable
solid_contact_with_gravity_unstructured_SOURCES = solid_contact_with_gravity.cc \
contact_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
solid_contact_with_gravity_unstructured_LDADD = -L@libdir@ -lconstitutive \
                                -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
circular_boulder_melt_SOURCES = circular_boulder_melt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_boulder_melt_LDADD = -L@libdir@ -lunsteady_heat -lconstitutive -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
circular_boulder_solar_radiation_SOURCES = circular_boulder_solar_radiation.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_boulder_solar_radiation_LDADD = -L@libdir@ -lunsteady_heat -lconstitutive -lsolid -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
two_d_unsteady_heat_SOURCES = two_d_unsteady_heat.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
spring_contact_SOURCES = spring_contact.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spring_contact_LDADD = -L@libdir@ -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------

# Sources for executable
spring_contact2_SOURCES = spring_contact2.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spring_contact2_LDADD = -L@libdir@ -lgeneric \
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
backup.bash			     plot_sb_melt.bash \
crop.bash			     plot_solar_radiation.bash \
debug.bash			     plot_solid_contact.bash \
plot_circular_boulder_melt.bash      plot_spring_contact_landscape.bash \
plot_hertz.bash			     re-enumerate.bash \
plot_melt.bash			     run_sb_melt.bash \
plot_melt_with_free_penetrator.bash  run_solid_contact.bash \
plot_pseudo_melt.bash		     run_spring_contact.bash \
plot_sb.bash                         plot_linear_hertz.bash


EXTRA_DIST+= \
before_and_after_contact_pressure.pvsm	solar_radiation_cones.pvsm \
bins.pvsm				solar_radiation.pvsm \
circular_boulder_melt.pvsm		solid_contact.pvsm \
compare_displ.pvsm			solid_contact_with_newton_iter.pvsm \
displ.pvsm				spring_contact_landscape.pvsm \
hertz_validation.pvsm			stefan_boltzmann_melt.pvsm \
melt.pvsm				stefan_boltzmann.pvsm \
proper_melt.pvsm			stefan_boltzmann_rays.pvsm \
pseudo_melt.pvsm                        linear_hertz_validation.pvsm



EXTRA_DIST+= \
contact.map  stefan_boltzmann_melt_validation.map  view_factor.map \
exact.map    stefan_boltzmann_validation.map

EXTRA_DIST+= \
heated_linear_hertz_validation.pvsm \
plot_heated_linear_hertz.bash \
animate.py








#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=b_convection_sphere multi_domain_b_convect_sph

# Sources for executable
b_convection_sphere_SOURCES = half_rectangle_with_hole_mesh.h \
                             b_convection_sphere.cc \
                             axisym_buoyant_navier_stokes.h
# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
b_convection_sphere_LDADD = -L@libdir@ \
 -laxisym_advection_diffusion \
 -laxisym_navier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
multi_domain_b_convect_sph_SOURCES = \
 multi_domain_b_convect_sph.cc\
 half_rectangle_with_hole_mesh.h \
 multi_domain_axisym_boussinesq_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_b_convect_sph_LDADD = -L@libdir@ \
 -laxisym_advection_diffusion \
 -laxisym_navier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

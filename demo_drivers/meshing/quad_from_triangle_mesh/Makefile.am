#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
adaptive_moving_block_navier_stokes \
adaptive_unstructured_scattering_quad \
unstructured_two_d_solid \
adaptive_unstructured_two_d_solid

#---------------------------------------------------------------------------

# Sources for executable
adaptive_unstructured_two_d_solid_SOURCES = adaptive_unstructured_two_d_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_two_d_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
unstructured_two_d_solid_SOURCES = unstructured_two_d_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
adaptive_unstructured_scattering_quad_SOURCES = adaptive_unstructured_scattering_quad.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_scattering_quad_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for first executable
adaptive_moving_block_navier_stokes_SOURCES = adaptive_moving_block_navier_stokes.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_moving_block_navier_stokes_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

EXTRA_DIST += box_hole.1.ele	box_hole.1.node  box_hole.1.poly  \
              box_hole.poly  solid.1.ele  solid.1.node  solid.1.poly  solid.poly

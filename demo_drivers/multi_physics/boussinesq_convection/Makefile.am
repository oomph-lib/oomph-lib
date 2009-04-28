#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
boussinesq_convection boussinesq_convection_fd boussinesq_convection_analytic \
refineable_b_convection \
multimesh_boussinesq_convection multimesh_ref_b_convection

#---------------------------------------------------------------------

# Sources for first executable
boussinesq_convection_fd_SOURCES = boussinesq_convection.cc 
boussinesq_convection_fd_CXXFLAGS = -DUSE_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
boussinesq_convection_fd_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------


# Sources for executable
boussinesq_convection_SOURCES = boussinesq_convection.cc 
boussinesq_convection_CXXFLAGS = -DUSE_OFF_DIAGONAL_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
boussinesq_convection_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

# Sources for third executable
boussinesq_convection_analytic_SOURCES = boussinesq_convection.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
boussinesq_convection_analytic_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------


# Sources for executable
refineable_b_convection_SOURCES = refineable_b_convection.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_b_convection_LDADD = -L@libdir@ -ladvection_diffusion \
  -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------

# Sources for executable
multimesh_boussinesq_convection_SOURCES = multimesh_boussinesq_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_boussinesq_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
multimesh_ref_b_convection_SOURCES = \
 multimesh_ref_b_convection.cc \
 my_boussinesq_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_ref_b_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


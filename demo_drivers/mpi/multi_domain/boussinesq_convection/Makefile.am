#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
refineable_b_convection \
multimesh_boussinesq_convection \
multimesh_boussinesq_convection_fd \
multimesh_ind_ref_b_convection

#----------------------------------------------------------------------

# Sources for executable
refineable_b_convection_SOURCES = refineable_b_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_b_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
multimesh_boussinesq_convection_fd_SOURCES = multimesh_boussinesq_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_boussinesq_convection_fd_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

multimesh_boussinesq_convection_fd_CXXFLAGS = \
 -DUSE_FD_JACOBIAN_FOR_NAVIER_STOKES_ELEMENT \
 -DUSE_FD_JACOBIAN_FOR_ADVECTION_DIFFUSION_ELEMENT

#----------------------------------------------------------------------

# Sources for executable
multimesh_ind_ref_b_convection_SOURCES = \
 multimesh_ind_ref_b_convection.cc \
 my_boussinesq_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_ind_ref_b_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

EXTRA_DIST += multimesh_boussinesq_partition.dat \
              multimesh_ind_ref_b_partition.dat \
              refineable_b_convection_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  


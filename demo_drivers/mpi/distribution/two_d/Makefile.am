#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
two_d_mesh_dist \
fish_poisson \
adaptive_driven_cavity \
circular_driven_cavity \
airy_cantilever2_adapt \
two_d_poisson_flux_bc_adapt \
prescribed_displ_lagr_mult

#----------------------------------------------------------------------

# Sources for executable
two_d_mesh_dist_SOURCES = two_d_mesh_dist.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_mesh_dist_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
fish_poisson_SOURCES = fish_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
adaptive_driven_cavity_SOURCES = adaptive_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
circular_driven_cavity_SOURCES = circular_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
airy_cantilever2_adapt_SOURCES = airy_cantilever2.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever2_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
airy_cantilever2_adapt_CXXFLAGS = -DREFINE

#-----------------------------------------------------------------------

# Sources for executable
two_d_poisson_flux_bc_adapt_SOURCES = two_d_poisson_flux_bc_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
two_d_poisson_flux_bc_adapt_LDADD = -L@libdir@ -lpoisson \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult_SOURCES = prescribed_displ_lagr_mult.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
prescribed_displ_lagr_mult_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------




# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  


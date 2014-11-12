#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
refineable_b_convection \
multi_domain_boussinesq_convection \
multi_domain_boussinesq_convection_fd \
multi_domain_ref_b_convection

#----------------------------------------------------------------------

# Sources for executable
multi_domain_ref_b_convection_SOURCES = multi_domain_ref_b_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_ref_b_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
multi_domain_boussinesq_convection_SOURCES = multi_domain_boussinesq_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_boussinesq_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
multi_domain_boussinesq_convection_fd_SOURCES = multi_domain_boussinesq_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_boussinesq_convection_fd_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

multi_domain_boussinesq_convection_fd_CXXFLAGS = -DUSE_FD_JACOBIAN_NST_IN_MULTI_DOMAIN_BOUSSINESQ

#----------------------------------------------------------------------

EXTRA_DIST += multi_domain_boussinesq_partition.dat \
              multi_domain_ref_b_partition.dat \
              refineable_b_convection_partition.dat \
              multi_domain_boussinesq_partition_2.dat \
              prune_and_balance.mcr

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  


#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
boussinesq_convection boussinesq_convection_fd boussinesq_convection_analytic \
refineable_b_convection \
multi_domain_boussinesq_convection \
multi_domain_ref_b_convection \
multi_domain_ref_b_convection_fd_for_external_data 

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
multi_domain_boussinesq_convection_SOURCES = multi_domain_boussinesq_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_boussinesq_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
multi_domain_ref_b_convection_SOURCES = \
 multi_domain_ref_b_convection.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_ref_b_convection_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------

# Sources for executable
multi_domain_ref_b_convection_fd_for_external_data_SOURCES = \
 multi_domain_ref_b_convection.cc 

multi_domain_ref_b_convection_fd_for_external_data_CXXFLAGS = \
-DUSE_FD_FOR_DERIVATIVES_WRT_EXTERNAL_DATA_IN_MULTI_DOMAIN_BOUSSINESQ

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multi_domain_ref_b_convection_fd_for_external_data_LDADD = \
                -L@libdir@ -ladvection_diffusion -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


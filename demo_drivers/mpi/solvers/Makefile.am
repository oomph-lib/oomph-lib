#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
driven_cavity \
direct_solver_test \
airy_cantilever \
fsi_channel_with_leaflet


# Do we have the TRILINOS library?
#=================================
if OOMPH_HAS_TRILINOS
   check_PROGRAMS += TrilinosSolver_test two_d_multi_poisson
endif

#----------------------------------------------------------------------


# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Sources for executable
direct_solver_test_SOURCES = direct_solver_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
direct_solver_test_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Sources for executable
airy_cantilever_SOURCES = airy_cantilever.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
fsi_channel_with_leaflet_SOURCES = fsi_channel_with_leaflet.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_channel_with_leaflet_LDADD = \
 -L@libdir@ -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)



#-----------------------------------------------------------------------

# Do we have the TRILINOS library?
#=================================
if OOMPH_HAS_TRILINOS

# Sources the executable depends on:
TrilinosSolver_test_SOURCES = TrilinosSolver_test.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" library
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
TrilinosSolver_test_LDADD = -L@libdir@ -lpoisson -lgeneric \
                       $(EXTERNAL_LIBS) \
                       $(FLIBS)

#-----------------------------------------------------------------------

# Sources for executable
two_d_multi_poisson_SOURCES = two_d_multi_poisson.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_multi_poisson_LDADD = -L@libdir@ -lgeneric \
	$(EXTERNAL_LIBS) $(FLIBS)

endif



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I. The slightly convoluted
# reference to generic is required to allow the multi-poisson driver code
# to include the multi-poisson header file in the serial demo driver
# directory to find the generic header files.
AM_CPPFLAGS +=  -I@includedir@  -I@includedir@/generic

EXTRA_DIST += multi_poisson_partition.dat 

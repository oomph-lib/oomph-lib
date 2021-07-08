#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS = circular_driven_cavity circular_driven_cavity2 


# SINGLE DRIVER FILE
#-------------------------------------


# Sources the executable depends on:
circular_driven_cavity_SOURCES = circular_driven_cavity.cc validate.sh

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
circular_driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# SINGLE DRIVER FILE
#-------------------------------------


# Sources the executable depends on:
circular_driven_cavity2_SOURCES = circular_driven_cavity2.cc validate.sh

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
circular_driven_cavity2_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


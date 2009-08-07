#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS = rayleigh_traction_channel


# FIRST EXECUTABLE: SINGLE DRIVER FILE
#-------------------------------------


# Sources the executable depends on:
rayleigh_traction_channel_SOURCES = rayleigh_traction_channel.cc validate.sh

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
rayleigh_traction_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

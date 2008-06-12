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

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES = -I@includedir@

# Only run self test if we have python for fpdiff.py
#---------------------------------------------------
if HAVE_PYTHON
  # Test script
  TESTS =  validate.sh
else
  # Just run the executables
  TESTS = rayleigh_traction_channel
endif


# Distribute validation script and validation data
EXTRA_DIST = validate.sh validata

#Add in additional cleanup commands
clean-local: 
	rm -rf Validation *~

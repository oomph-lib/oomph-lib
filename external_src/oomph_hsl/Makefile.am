###################################################
# Frontal solver from the HSL library 
# We turn it into the library oomph_hsl to avoid
# clashes with other installations of the hsl 
# library, likely to be accessible via -lhsl
###################################################

junk:
	echo "$(AM_FFLAGS) " $(AM_FFLAGS) 
	echo "$(FFLAGS) " $(FFLAGS) 
	echo "$(F77LD) " $(F77LD)


# Define the sources 
sources = dummy.cc 

if HAVE_HSL_SOURCES
 sources += frontal.f
 EXTRA_DIST = dummy_frontal.f
else
 sources += dummy_frontal.f
endif

# Headers (none)
headers =

# Define name of library
libname = oomph_hsl

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_hsl.la

# Sources that the library depends on:
#-------------------------------------
liboomph_hsl_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_hsl

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers)

# Dummy header file - there's no real header for the fortran stuff
# we provide a cfortran-based header file for the relevant
# subroutines in src/generic
$(libname).h: 
	echo "//Empty" > $(libname).h

# Dummy for the frontal.f file -- needed during make dist even if
# should really be excluded because of the if HAVE_HSL_SOURCES block.
# Touch turns it into a zero length file which can get wiped again
# by bin/customise_dist.sh 
$(srcdir)/frontal.f:
	if [ ! -f $(srcdir)/frontal.f ]; then \
	  if [ -f $(srcdir)/../../private/external_src/hsl/frontal.f ]; then \
	      echo "Creating link to private/external_src/hsl/frontal.f"; \
	      ln -sf $(srcdir)/../../private/external_src/hsl/frontal.f $(srcdir); \
	  fi; \
	  touch $(srcdir)/frontal.f; \
	fi; 


# Tidy up: Remove symbolic link to private_external_src
clean-local:
	if [ -L $(srcdir)/frontal.f ] || ( [ -e $(srcdir)/frontal.f ] && [ ! -s $(srcdir)/frontal.f ] ); then \
	    echo "Removing empty file or symbolic link to private/external_src/hsl/frontal.f";\
	    rm $(srcdir)/frontal.f;\
	fi;

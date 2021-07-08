#---------------------------------------------------------------------------

# These are the subdirectories we are going to build:
SUBDIRS=

# Extra files for distribution
EXTRA_DIST=\
archived_discretisations

#---------------------------------------------------------------------------
# Define the headers
headers=\
mixed_order_petrov_galerkin_space_time_fluid_traction_elements.h\
mixed_order_petrov_galerkin_space_time_navier_stokes_elements.h\
refineable_mixed_order_petrov_galerkin_space_time_navier_stokes_elements.h\
mixed_order_petrov_galerkin_space_time_block_preconditionable_elements.h

# Define the sources
sources=\
mixed_order_petrov_galerkin_space_time_navier_stokes_elements.cc\
refineable_mixed_order_petrov_galerkin_space_time_navier_stokes_elements.cc\
mixed_order_petrov_galerkin_space_time_block_preconditionable_elements.cc

# Include files which shouldn't be compiled
incl_cc_files=

# Define name of library
libname = space_time_navier_stokes

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libspace_time_navier_stokes.la

# Sources that the library depends on:
#-------------------------------------
libspace_time_navier_stokes_la_SOURCES = $(headers_and_sources)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/$(libname)

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files)

# Required libraries
#-------------------
# Include path for library headers -- need to refer to include files
# in src as they will not have been installed yet!
AM_CPPFLAGS += -I$(top_srcdir)/src

# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers) $(incl_cc_files) )
	(echo "$(headers) $(incl_cc_files)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list )
else
install-data-hook:
endif

# Cleanup: Kill automatically generated combined header file
clean-local:
	rm -f $(libname).h  *~
distclean-local:
	rm -f $(libname).h  *~
#---------------------------------------------------------------------------



#---------------------------------------------------------------------------
# These are the subdirectories we are going to build:
SUBDIRS=

# Extra files for distribution
EXTRA_DIST=
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Define the headers
headers=\
space_time_unsteady_heat_elements.h\
refineable_space_time_unsteady_heat_elements.h\
space_time_unsteady_heat_block_preconditionable_elements.h

# Define the sources
sources=\
space_time_unsteady_heat_elements.cc\
refineable_space_time_unsteady_heat_elements.cc\
space_time_unsteady_heat_block_preconditionable_elements.cc

# Include files which shouldn't be compiled
incl_cc_files=

# Define name of library
libname = space_time_unsteady_heat_equal_order_galerkin

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libspace_time_unsteady_heat_equal_order_galerkin.la

# Sources that the library depends on:
#-------------------------------------
libspace_time_unsteady_heat_equal_order_galerkin_la_SOURCES = \
$(headers_and_sources)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

# The subdirectory in build/include where the headers (or symlinked headers)
# are going to go.
# NOTE: We don't normally have to create a folder two directories deep in the
# Makefiles but because the header:
#      ../../space_time_block_preconditioner/general_purpose_space_\
#                                  time_block_preconditionable_elements.h
# is referenced here, we have to. If we don't, there will be a linking problem
# when trying to make a driver code that includes this header. This is because
# the directory in build/include (where headers get placed after being built)
# won't respect the same directory structure. We have to be careful here to
# make sure this won't be a problem.
library_include_subdir = space_time_unsteady_heat/equal_order_galerkin

#Here's the subdirectory where the actual header files are placed
library_includedir = $(includedir)/$(library_include_subdir)

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
	echo $(library_include_subdir) $(headers) > all_$(libname).aux
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
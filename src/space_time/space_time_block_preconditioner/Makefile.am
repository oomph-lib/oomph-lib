#---------------------------------------------------------------------------
# These are the subdirectories we are going to build:
SUBDIRS=

# Distribute the list of block_prec names and the build script for the include files
EXTRA_DIST=
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Define the headers
headers=\
general_purpose_space_time_block_preconditioner.h\
general_purpose_space_time_subsidiary_block_preconditioner.h\
general_purpose_space_time_block_preconditionable_elements.h

# Define the non-block_prec sources
sources=\
general_purpose_space_time_block_preconditioner.cc\
general_purpose_space_time_subsidiary_block_preconditioner.cc

# Include files which shouldn't be compiled
incl_cc_files=

# Define name of library
libname = space_time_block_preconditioner

# Combine them
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libspace_time_block_preconditioner.la

# Sources that the library depends on
libspace_time_block_preconditioner_la_SOURCES=$(headers_and_sources)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

# Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/$(libname)

# These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers_and_sources)

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

dummy_$(libname).h:
	echo  $(libname) $(headers) $(non_block_prec_headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < all_$(libname).aux > $(libname).h
	rm all_$(libname).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers_and_sources) )
	(echo "$(headers_and_sources)" > include_files.list )
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
	rm -f $(libname).h *~
#---------------------------------------------------------------------------

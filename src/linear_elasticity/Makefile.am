
# A few file definitions
#-----------------------

# Define the sources
sources =  \
elasticity_tensor.cc linear_elasticity_elements.cc \
refineable_linear_elasticity_elements.cc Tlinear_elasticity_elements.cc \
linear_elasticity_traction_elements.cc 

# Include files which shouldn't be compiled
incl_cc_files =

# Define the headers
headers =  \
elasticity_tensor.h linear_elasticity_elements.h \
linear_elasticity_traction_elements.h \
Tlinear_elasticity_elements.h \
refineable_linear_elasticity_elements.h


# Template only files. These should be included in include directory
# but should not be listed in combined header file
templates = 

# Define name of library
libname = linear_elasticity

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liblinear_elasticity.la

# Sources that the library depends on: 
#-------------------------------------
liblinear_elasticity_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h 

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/linear_elasticity

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files) $(templates)

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
	(cd $(library_includedir) && rm -f $(headers) $(incl_cc_files) $(templates) )
	(echo "$(headers) $(incl_cc_files) $(templates)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Cleanup: Kill automatically generated combined header file
clean-local: 
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~


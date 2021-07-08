
# A few file definitions
#-----------------------

# Define the sources
sources =  \
linear_wave_elements.cc refineable_linear_wave_elements.cc 

# linear_wave_flux_elements.cc 

# Include files which shouldn't be compiled
incl_cc_files =

# Define the headers
headers =  \
linear_wave_elements.h \
linear_wave_flux_elements.h \
refineable_linear_wave_elements.h 

# Define name of library
libname = linear_wave

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liblinear_wave.la

# Sources that the library depends on: 
#-------------------------------------
liblinear_wave_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h 

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/linear_wave

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files)

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
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~



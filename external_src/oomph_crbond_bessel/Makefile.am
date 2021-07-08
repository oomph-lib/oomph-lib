###################################################
# Collection of Bessel functions and the like, implemented
# in C++ by C. Bond (http://www.crbond.com/), using the
# algorithms in  Zhang and Jin's book "Computation of Special Functions".
# John Wiley and Sons, 1996.
# This is based on the version dated 06/04 which appears to have
# corrected errors in complex Bessel functions.
###################################################



# Define the sources
sources = crbond_bessel.cc

# Define the headers: 
headers = crbond_bessel.h 


# Define name of library
libname = oomph_crbond_bessel

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_crbond_bessel.la


# Sources that the library depends on:
#-------------------------------------
liboomph_crbond_bessel_la_SOURCES = $(headers_and_sources) $(extra_headers)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_crbond_bessel

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(extra_headers) 


# Combined header file 
#---------------------
# Combined header file was hand-written. It is the same as 
# metis.h (provided with METIS) but the includes are redirected
# to the metis_4.0 include directory where they will be
# found after installation.
$(libname).h: dummy_$(libname).h



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
	(cd $(library_includedir) && rm -f $(headers) $(extra_headers) )
	(echo "$(headers) $(extra_headers)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Cleanup: Kill automatically generated combined header file
clean-local: 
	rm -f $(libname).h
distclean-local:
	rm -f $(libname).h
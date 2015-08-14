###################################################
# Triangle installed as library
###################################################



# Define the sources -- Shewchuk's file with compiler
# macros for  ANSI_DECLARATORS, LINUX and TRILIBRARY
sources = triangle.c

# Define the headers. Unlikely to be used directly as it's plain
# C...
headers = triangle.h

#Specific compilation for triangle with no optimisation
triangle.lo : triangle.c
	$(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(AM_CFLAGS) -O0 -Wall \
	-MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo


# Define name of library
libname = oomph_triangle

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_triangle.la


# Sources that the library depends on:
#-------------------------------------
liboomph_triangle_la_SOURCES = $(headers_and_sources) $(extra_headers)
liboomph_triangle_la_LDFLAGS = -static


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_triangle

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

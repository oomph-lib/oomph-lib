###################################################
# Tetgen installed as library
###################################################

# There is some hackery going on here because on my main test machine
# I kept getting seg faults when tetgen.cxx was compiled with full
# optimisation and predicates.cxx is supposed to be compiled with no
# optimisation. Hopefully, one day automake will sort this problem out
# with program specific flags, but it may take some time.

# The problem may be to do with some strange initialisations within loops
# and loop unrolling, but that could be a complete red herring.

# Define the sources -- Hang Si's file version 1.4 with compiler
# macros TETLIBRARY enabled in the header
sources = tetgen.cxx predicates.cxx

# Define the headers
headers = tetgen.h

#Specific compilation for tetgen at optimisation -O2 (higher fails
#on my machine)
tetgen.lo : tetgen.cxx
	$(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(AM_CXXFLAGS) -O2 -Wall \
	-MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo


#Specific compilation for predicates without any optimisation
predicates.lo : predicates.cxx
	$(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(AM_CXXFLAGS) -O0 -Wall \
	-MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo


# Define name of library
libname = oomph_tetgen

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_tetgen.la


# Sources that the library depends on:
#-------------------------------------
liboomph_tetgen_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_tetgen

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers)

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
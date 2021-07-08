#####################################################
# DOUBLE PRECISION SOURCES FROM ARPACK96 + PATCHES
# We turn these into the libraray oomph_arpack
# to avoid any clashes with existing installations
# of arpack (probably installed as -larpack)
####################################################

# Define the sources
sources = dummy.cc

if HAVE_ARPACK_SOURCES
 sources += all_arpack_sources.f
 EXTRA_DIST = dummy_arpack.f
else
 sources += dummy_arpack.f
endif

# Define the headers
headers = debug.h stat.h

#Define the name of the library
libname = oomph_arpack

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

#Define the library that goes into the lib directory
lib_LTLIBRARIES = liboomph_arpack.la

#Source on which the library depends
liboomph_arpack_la_SOURCES = $(headers_and_sources)

#Include headers
#Headers that should be included in the  $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_arpack

# Dummy header file - there's no real header for the fortran stuff
# we provide a cfortran-based header file for the relevant
# subroutines in src/generic
$(libname).h: 
	echo "//Empty" > $(libname).h

#These are the header files that are to be placed in subdirectory
#library_include_HEADERS=$(headers)


#NONE HERE
# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
#$(libname).h: dummy_$(libname).h

#NONE HERE
# When building the combined header, we also create a symbolic
# link to it from the directory above this one, so that
# it can be used by other libaries (before "make install" is run)
# in exactly the same way as after its installation in the include
# directory. (After the installation, the combined header lives
# in the directory above the actual individual header files)
#dummy_$(libname).h:  $(headers)
#	echo  $(libname) $(headers) > all_$(libname).aux
#	$(AWK) -f $(top_builddir)/bin/headers.awk < \
#	          all_$(libname).aux > $(libname).h
#	rm all_$(libname).aux
#	rm -f ../$(libname).h 
#	(cd .. && $(LN_S) $(libname)/$(libname).h $(libname).h)

#NONE HERE
# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
#if SYMBOLIC_LINKS_FOR_HEADERS
#install-data-hook: 
#	(cd $(library_includedir) && rm -f $(headers) )
#	(echo "$(headers)" > include_files.list )
#	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
#	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
#	(rm -r include_files.list.aux include_files.list ) 
#else
install-data-hook: 
#endif

#Dummy for all_arpack_sources.f file -- needed during make dist even 
#if it should not be included because the HAVE_ARPACK_SOURCES flag is false.
#Touch creates a zero-length file
$(srcdir)/all_arpack_sources.f: $(srcdir)/debug.h $(srcdir)/stat.h
	if [ ! -f $(srcdir)/all_arpack_sources.f ]; then \
	 if [ -f $(srcdir)/../../private/external_src/arpack/all_arpack_sources.f ]; then \
	    echo "Creating link to private/external_src/arpack/all_arpack_sources.f"; \
	    ln -s $(srcdir)/../../private/external_src/arpack/all_arpack_sources.f $(srcdir); \
	 fi; \
	 touch $(srcdir)/all_arpack_sources.f; \
	fi;  

#Dummy for debug.h file -- needed during make dist even 
#if it should not be included because the HAVE_ARPACK_SOURCES flag is false.
#Touch creates a zero-length file
$(srcdir)/debug.h:
	if [ ! -f $(srcdir)/debug.h ]; then \
	 if [ -f $(srcdir)/../../private/external_src/arpack/debug.h ]; then \
	    echo "Creating link to private/external_src/arpack/debug.h"; \
	    ln -s $(srcdir)/../../private/external_src/arpack/debug.h $(srcdir); \
	 fi; \
	 touch $(srcdir)/debug.h; \
	fi;  

#Dummy for stat.h file -- needed during make dist even 
#if it should not be included because the HAVE_ARPACK_SOURCES flag is false.
#Touch creates a zero-length file
$(srcdir)/stat.h:
	if [ ! -f $(srcdir)/stat.h ]; then \
	 if [ -f $(srcdir)/../../private/external_src/arpack/stat.h ]; then \
	    echo "Creating link to private/external_src/arpack/stat.h"; \
	    ln -s $(srcdir)/../../private/external_src/arpack/stat.h $(srcdir); \
	 fi; \
	 touch $(srcdir)/stat.h; \
	fi;  


# Tidy up: Remove symbolic link to private_external_src
clean-local:                                           
	if [ -L $(srcdir)/all_arpack_sources.f ] || ( [ -e $(srcdir)/all_arpack_sources.f ] && [ ! -s $(srcdir)/all_arpack_sources.f ] ); then \
	     echo "Removing empty file or symbolic link to private_external_src/arpack/all_arpack_sources.f";\
	     rm -f $(srcdir)/all_arpack_sources.f;\
	fi;         	
	if [ -L $(srcdir)/debug.h ] || ( [ -e $(srcdir)/debug.h ] && [ ! -s $(srcdir)/debug.h ] ); then \
	     echo "Removing empty file or symbolic link to private_external_src/arpack/debug.h";\
	     rm -f $(srcdir)/debug.h;\
	fi;         	
	if [ -L $(srcdir)/stat.h ] || ( [ -e $(srcdir)/stat.h ] && [ ! -s $(srcdir)/stat.h ] ); then \
	     echo "Removing empty file or symbolic link to private_external_src/arpack/stat.h";\
	     rm -f $(srcdir)/stat.h;\
	fi;         	
	rm -f $(libname).h

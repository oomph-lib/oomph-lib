###################################################
# Lapack sources required by superlu
# We turn it into the library oomph_lapack
# to avoid clashes with any existing installations
# of superlu (likely to be accessible via -llapack)
###################################################



# Define the sources
sources =  \
lsame.c  xerbla.c


# Define the headers: This is my hand-written combined header
headers =  oomph_lapack.h

# Define name of library
libname = oomph_lapack

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_lapack.la


# Sources that the library depends on:
#-------------------------------------
liboomph_lapack_la_SOURCES = $(headers_and_sources)
liboomph_lapack_la_LDFLAGS = -static


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_lapack

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers)


# When building the combined header, we also create a symbolic
# link to it from the directory above this one, so that
# it can be used by other libaries (before "make install" is run)
# in exactly the same way as after its installation in the include
# directory. (After the installation, the combined header lives
# in the directory above the actual individual header files)
dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	          all_$(libname).aux > $(libname).h
	rm all_$(libname).aux
	rm -f ../$(libname).h 
	(cd .. && $(LN_S) $(libname)/$(libname).h $(libname).h)

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_includedir) && rm -f $(headers) )
	(echo "$(headers)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


###################################################
# Serial metis as distributed with parmetis 3.1.1
#
# Turn into oomph_metis_from_parmetis_3.1.1 to avoid 
# clashes with any existing installations of metis 
# which are likely to be accessible via -lmetis
###################################################



# Define the sources
sources =  \
balance.c     graph.c	       mcoarsen.c    mkwayrefine.c  pqueue.c\
bucketsort.c  initpart.c       memory.c      mmatch.c	    refine.c\
ccgraph.c     kmetis.c	       mesh.c	     mmd.c	    separator.c\
checkgraph.c  kvmetis.c        meshpart.c    mpmetis.c	    sfm.c\
coarsen.c     kwayfm.c	       mfm2.c	     mrefine2.c     srefine.c\
compress.c    kwayrefine.c     mfm.c	     mrefine.c	    stat.c\
debug.c       kwayvolfm.c      mincover.c    mutil.c        \
estmem.c      kwayvolrefine.c  minitpart2.c  myqsort.c	    subdomains.c\
fm.c	      match.c	       minitpart.c   ometis.c	    timing.c\
fortran.c     mbalance2.c      mkmetis.c     parmetis.c     util.c\
frename.c     mbalance.c       mkwayfmh.c    pmetis.c


# Define the headers: This is the combined header provided with
# and used by metis itself
headers = metis.h 

# Additional headers that need to be distributed
extra_headers = defs.h	macros.h  proto.h  rename.h  stdheaders.h  struct.h


# Additional files that come with parmetis but aren't actually used by
# us. Included for completeness.
EXTRA_DIST = stats.c


# Define name of library
libname = oomph_metis_from_parmetis_3.1.1

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_metis_from_parmetis_3.1.1.la


# Sources that the library depends on:
#-------------------------------------
liboomph_metis_from_parmetis_3_1_1_la_SOURCES = $(headers_and_sources) $(extra_headers)
liboomph_metis_from_parmetis_3_1_1_la_LDFLAGS = -static

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_metis_from_parmetis_3.1.1

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(extra_headers) 


# Combined header file 
#---------------------
# Combined header file was hand-written. It is the same as 
# metis.h (provided with METIS) but the includes are redirected
# to the include directory where they will be
# found after installation.
$(libname).h: dummy_$(libname).h

# When building the combined header, we also create a symbolic
# link to it from the directory above this one, so that
# it can be used by other libaries (before "make install" is run)
# in exactly the same way as after its installation in the include
# directory. (After the installation, the combined header lives
# in the directory above the actual individual header files)
dummy_$(libname).h: 
	rm -f ../$(libname).h 
	(cd .. && $(LN_S) $(libname)/$(libname).h $(libname).h)


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


###################################################
#   PARMETIS 3.1.1 
# Turn into oomph_parmetis_3.1.1 to avoid clashes with any
# existing installations of metis which are 
# likely to be accessible via -lmetis
###################################################




# If we don't have malloc.h (most likely on a mac)
#-------------------------------------------------
# make an empty dummy version
#----------------------------
if !HAVE_MALLOCH

# Dummy dependency to force build of malloc.h
defs.h: malloc.h

# Make it!
malloc.h:
	rm -f malloc.h
	touch malloc.h

# Already have malloc.h  -- this is just for testing
#---------------------------------------------------
else

# Dummy dependency to force build of malloc_shite.h
defs.h: malloc_shite.h

# Make it!
malloc_shite.h:
	rm -f malloc_shite.h
	touch malloc_shite.h

# Endif for malloc
#-----------------
endif



# Define the sources
sources =  \
frename.c	 kmetis.c	mmetis.c       selectq.c \
akwayfm.c	 gkmetis.c	 kwaybalance.c	move.c	       serial.c \
ametis.c	 grsetup.c	 kwayfm.c	msetup.c       setup.c \
backcompat.c	 iidxsort.c	 kwayrefine.c	node_refine.c  stat.c \
balancemylink.c  iintsort.c	 lmatch.c	ometis.c       timer.c \
comm.c		 ikeysort.c	 match.c	order.c        util.c \
csrmatch.c	 ikeyvalsort.c	 mdiffusion.c	pspases.c      wave.c \
debug.c		 initbalance.c 	redomylink.c   weird.c \
diffutil.c	 initmsection.c  memory.c	remap.c        xyzpart.c \
fpqueue.c	 initpart.c	 mesh.c		rmetis.c      

# Define the headers: This is the combined header provided with
# and used by metis itself
headers = parmetislib.h

# Additional headers that need to be distributed
extra_headers = defs.h	macros.h  proto.h  rename.h  stdheaders.h  struct.h 

# Define name of library
libname = oomph_parmetis_3.1.1

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_parmetis_3.1.1.la


# Sources that the library depends on:
#-------------------------------------
liboomph_parmetis_3_1_1_la_SOURCES = $(headers_and_sources) $(extra_headers)
liboomph_parmetis_3_1_1_la_LDFLAGS = -static


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_parmetis_3.1.1

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(extra_headers) 


# Combined header file 
#---------------------

# Combined header file was hand-written. It is the same as 
# parmetis.h (provided with METIS) but the includes are redirected
# to the parmetis_3.1.1 include directory where they will be
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
	(cd ..; $(LN_S) parmetis.h $(libname).h) 


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

# Tidy
clean-local:
	rm -f malloc.h malloc_shite.h


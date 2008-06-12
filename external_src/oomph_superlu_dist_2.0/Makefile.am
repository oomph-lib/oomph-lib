###################################################
# DOUBLE PRECISION SOURCES FROM SuperLU_dist_2.0
# We turn it into the library oomph_superlu_dist_2.0
# to avoid clashes with any existing installations
# of superlu (likely to be accessible via -lsuperlu)
###################################################

# Define the sources
sources =  \
 symbfact.c sp_ienv.c etree.c sp_colorder.c get_perm_c.c \
	  mmd.c colamd.c comm.c memory.c util.c superlu_grid.c \
	  pxerbla.c superlu_timer.c GetDiagU.c mc64ad.f mc21.f \
          slamch.c dlangs.c dgsequ.c dlaqgs.c dutil.c \
	  dmemory.c dsp_blas2.c dsp_blas3.c dmyblas2.c \
          dldperm.c ddistribute.c pdgstrf_irecv.c pdgstrs_Bglobal.c \
	  pdgstrs1.c pdgssvx_ABglobal.c pdgsrfs_ABXglobal.c pdgsmv_AXglobal.c\
	  pdgssvx.c pdgstrs.c pddistribute.c pdlangs.c pdutil.c \
	  pdgsequ.c pdlaqgs.c pdgsrfs.c pdgsmv.c pdgstrs_lsum.c

# These files are included in superlu but are actually part of
# lapack which can cause clashes with other packages that
# also include them.
#  xerbla.c  lsame.c dlamch.c  
EXTRA_DIST =  xerbla.c.back  lsame.c.back dlamch.c.back  

# Define the headers
headers =  \
Cnames.h  superlu_ddefs.h  supermatrix.h\
colamd.h  machines.h  superlu_defs.h  util_dist.h


# Define name of library
libname = oomph_superlu_dist_2.0

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_superlu_dist_2.0.la


# Sources that the library depends on:
#-------------------------------------
liboomph_superlu_dist_2_0_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_superlu_dist_2.0

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers)

# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
$(libname).h: dummy_$(libname).h

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

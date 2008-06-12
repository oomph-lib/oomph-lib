###################################################
# DOUBLE PRECISION SOURCES FROM SuperLU_3.0
# We turn it into the library oomph_superlu_3.0
# to avoid clashes with any existing installations
# of superlu (likely to be accessible via -lsuperlu)
###################################################



# Define the sources
sources =  \
colamd.c         dgssvx.c   dmyblas2.c      dsnode_dfs.c        mmd.c \
dcolumn_bmod.c   dgstrf.c   dpanel_bmod.c   dsp_blas2.c         relax_snode.c\
dcolumn_dfs.c    dgstrs.c   dpanel_dfs.c    dsp_blas3.c         sp_coletree.c\
dcopy_to_ucol.c  dlacon.c   dpivotgrowth.c  dutil.c             sp_ienv.c\
dgscon.c         dpivotL.c       get_perm_c.c        sp_preorder.c\
dgsequ.c         dlangs.c   dpruneL.c       heap_relax_snode.c superlu_timer.c\
dgsrfs.c         dlaqgs.c   dreadhb.c       util.c\
dgssv.c          dmemory.c  dsnode_bmod.c   memory.c \
dcomplex.c       zgssv.c    zgssvx.c        zsp_blas2.c \
zsp_blas3.c      zgscon.c \
zlangs.c         zgsequ.c   zlaqgs.c zpivotgrowth.c  \
zgsrfs.c         zgstrf.c   zgstrs.c zcopy_to_ucol.c \
zsnode_dfs.c     zsnode_bmod.c \
zpanel_dfs.c     zpanel_bmod.c  zreadhb.c \
zcolumn_dfs.c    zcolumn_bmod.c zpivotL.c zpruneL.c \
zmemory.c        zutil.c        zmyblas2.c

# These files are included in superlu but are actually part of
# lapack which can cause clashes with other packages that
# also include them.
#  xerbla.c  lsame.c dlamch.c  
EXTRA_DIST =  xerbla.c.back  lsame.c.back dzsum1.c.back izmax1.c.back \
              zlacon.c.back dlamch.c.back 

# Define the headers
headers =  slu_Cnames.h colamd.h supermatrix.h slu_util.h slu_ddefs.h \
           slu_dcomplex.h slu_zdefs.h

# Define name of library
libname = oomph_superlu_3.0

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_superlu_3.0.la


# Sources that the library depends on:
#-------------------------------------
liboomph_superlu_3_0_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
#include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_superlu_3.0

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers)


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
$(libname).h: dummy_$(libname).h

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


##################################################
# DOUBLE PRECISION SOURCES FROM SuperLU_4.3
# We turn it into the library oomph_superlu_4.3
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
zmemory.c        zutil.c        zmyblas2.c \
ilu_relax_snode.c ilu_heap_relax_snode.c mark_relax.c mc64ad.c qselect.c\
dgsisx.c dgsitrf.c dldperm.c ilu_ddrop_row.c ilu_dsnode_dfs.c \
ilu_dcolumn_dfs.c ilu_dpanel_dfs.c ilu_dcopy_to_ucol.c \
ilu_dpivotL.c ddiagonal.c

# These files are included in superlu but are actually part of
# lapack which can cause clashes with other packages that
# also include them. Plus various other files that we're not compiling
# but are keeping alive for completeness so that what we have in here
# is the complete set of sources for superlu
#  xerbla.c  lsame.c dlamch.c  
EXTRA_DIST =  xerbla.c.back  lsame.c.back dzsum1.c.back izmax1.c.back \
              zlacon.c.back dlamch.c.back \
ilu_zcopy_to_ucol.c \
scsum1.c \
cmemory.c \
scolumn_dfs.c \
cdiagonal.c \
slu_scomplex.h \
smemory.c \
ilu_zdrop_row.c \
cpivotgrowth.c \
ccopy_to_ucol.c \
sdiagonal.c \
cgstrs.c \
csp_blas2.c \
cgscon.c \
spivotgrowth.c \
cutil.c \
scopy_to_ucol.c \
sgstrs.c \
cgssv.c \
cldperm.c \
ssp_blas2.c \
sgscon.c \
ilu_zpanel_dfs.c \
sutil.c \
cgsisx.c \
ilu_ccolumn_dfs.c \
sgssv.c \
sldperm.c \
ilu_cdrop_row.c \
html_mainpage.h \
zldperm.c \
sgsisx.c \
ilu_scolumn_dfs.c \
ilu_sdrop_row.c \
cpanel_dfs.c \
dreadtriple.c \
cpanel_bmod.c \
clangs.c \
csp_blas3.c \
cgssvx.c \
spanel_dfs.c \
spanel_bmod.c \
zreadtriple.c \
slangs.c \
zgsisx.c \
ssp_blas3.c \
ilu_cpivotL.c \
sgssvx.c \
cgsitrf.c \
creadrb.c \
dreadrb.c \
ilu_cpanel_dfs.c \
ilu_zsnode_dfs.c \
ilu_spivotL.c \
ilu_zpivotL.c \
sgsitrf.c \
icmax1.c \
sreadrb.c \
ilu_spanel_dfs.c \
cgsrfs.c \
zgsitrf.c \
zreadrb.c \
clacon.c \
sgsrfs.c \
slamch.c \
csnode_dfs.c \
cpruneL.c \
ilu_ccopy_to_ucol.c \
slu_cdefs.h \
slacon.c \
cpivotL.c \
cgsequ.c \
ssnode_dfs.c \
creadhb.c \
spruneL.c \
ilu_scopy_to_ucol.c \
slu_sdefs.h \
spivotL.c \
ccolumn_bmod.c \
sgsequ.c \
sreadhb.c \
csnode_bmod.c \
cmyblas2.c \
ilu_csnode_dfs.c \
scolumn_bmod.c \
cgstrf.c \
scomplex.c \
ssnode_bmod.c \
ilu_zcolumn_dfs.c \
smyblas2.c \
creadtriple.c \
ilu_ssnode_dfs.c \
sgstrf.c \
zdiagonal.c \
claqgs.c \
sreadtriple.c \
slaqgs.c \
ccolumn_dfs.c


# Define the headers
headers =  slu_Cnames.h colamd.h supermatrix.h slu_util.h slu_ddefs.h \
           slu_dcomplex.h slu_zdefs.h superlu_enum_consts.h

# Define name of library
libname = oomph_superlu_4.3

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_superlu_4.3.la


# Sources that the library depends on:
#-------------------------------------
liboomph_superlu_4_3_la_SOURCES = $(headers_and_sources)
liboomph_superlu_4_3_la_LDFLAGS = -static 


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
#include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_superlu_4.3

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


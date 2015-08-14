###################################################
# DOUBLE PRECISION SOURCES FROM SuperLU_dist_3.0
# We turn it into the library oomph_superlu_dist_3.0
# to avoid clashes with any existing installations
# of superlu (likely to be accessible via -lsuperlu)
###################################################


# Define the sources
sources =  \
sp_ienv.c etree.c sp_colorder.c get_perm_c.c \
mmd.c comm.c memory.c util.c superlu_grid.c \
pxerbla.c superlu_timer.c GetDiagU.c mc64ad.c symbfact.c \
psymbfact.c psymbfact_util.c \
get_perm_c_parmetis.c \
dlangs.c dgsequ.c dlaqgs.c dutil.c \
dmemory.c dmyblas2.c dsp_blas2.c dsp_blas3.c \
dldperm.c ddistribute.c pdgstrf.c pdgstrs_Bglobal.c \
pdgstrs1.c pdgssvx_ABglobal.c pdgsrfs_ABXglobal.c pdgsmv_AXglobal.c\
pdgssvx.c pdgstrs.c pddistribute.c pdlangs.c pdutil.c \
pdgsequ.c pdlaqgs.c pdgsrfs.c pdgsmv.c pdgstrs_lsum.c \
pdsymbfact_distdata.c dreadhb.c dreadrb.c dreadtriple.c

# These files are included in superlu but are actually part of
# lapack which can cause clashes with other packages that
# also include them. Plus various other files that we're not compiling
# but are keeping alive for completeness so that what we have in here
# is the complete set of sources for superlu
EXTRA_DIST =  xerbla.c.back  lsame.c.back dlamch.c.back slamch.c.back \
zsp_blas2.c \
pzgstrf.c \
zutil.c \
pzlangs.c \
pzlaqgs.c \
pdgstrf_X1.c \
pdgstrs_Bglobal_Bsend.c \
pzgssvx.c \
pzdistribute.c \
pzgstrs_Bglobal.c \
zlangs.c \
zlaqgs.c \
pdgstrsL.c \
zmemory.c \
pdgstrf_sherry.c \
zdistribute.c \
pzgsrfs_ABXglobal.c \
zdistribute_mark.c \
pzgsmv.c \
pzgstrf_irecv.c \
old_colamd.c \
pzgsrfs.c \
pzgsmv_AXglobal.c \
superlu_zdefs.h \
zsp_blas3.c \
zmyblas2.c \
pzgssvx_ABglobal.c \
html_mainpage.h \
pzgstrs.c \
zldperm.c \
zreadhb.c \
pzgstrs1.c \
pzgsequ.c \
pzsymbfact_distdata.c \
pdgstrf_irecv.c \
dcomplex.c \
pzgstrs_lsum.c \
zgsequ.c \
dcomplex.h \
zreadtriple.c \
pzutil.c 

# Define the headers
headers =  \
superlu_defs.h  superlu_ddefs.h Cnames.h    machines.h	  psymbfact.h \
supermatrix.h old_colamd.h   util_dist.h superlu_enum_consts.h 

# superlu_zdefs.h dcomplex.h  

# Define name of library
libname = oomph_superlu_dist_3.0

# Combine headers and sources
headers_and_sources = $(headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_superlu_dist_3.0.la


# Sources that the library depends on:
#-------------------------------------
liboomph_superlu_dist_3_0_la_SOURCES = $(headers_and_sources)
liboomph_superlu_dist_3_0_la_LDFLAGS = -static


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_superlu_dist_3.0

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

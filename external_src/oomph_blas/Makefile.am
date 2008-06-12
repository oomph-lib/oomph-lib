###################################################
#   Blas sources
###################################################

# Define the sources
sources = dummy.cc \
caxpy.f   crotg.f   dcopy.f   dsymv.f   lsame.f   sspmv.f   zaxpy.f   zhpr.f\
ccopy.f   cscal.f   ddot.f    dsyr.f    sasum.f   sspr.f    zcopy.f   zhpr2.f\
cdotc.f   csrot.f   dgbmv.f   dsyr2.f   saxpy.f   sspr2.f   zdotc.f   zrotg.f\
cdotu.f   csscal.f  dgemm.f   dsyr2k.f  scasum.f  sswap.f   zdotu.f   zscal.f\
cgbmv.f   cswap.f   dgemv.f   dsyrk.f   scnrm2.f  ssymm.f   zdrot.f   zswap.f\
cgemm.f   csymm.f   dger.f    dtbmv.f   scopy.f   ssymv.f   zdscal.f  zsymm.f\
cgemv.f   csyr2k.f  dnrm2.f   dtbsv.f   sdot.f    ssyr.f    zgbmv.f   zsyr2k.f\
cgerc.f   csyrk.f   drot.f    dtpmv.f   sdsdot.f  ssyr2.f   zgemm.f   zsyrk.f\
cgeru.f   ctbmv.f   drotg.f   dtpsv.f   sgbmv.f   ssyr2k.f  zgemv.f   ztbmv.f\
chbmv.f   ctbsv.f   drotm.f   dtrmm.f   sgemm.f   ssyrk.f   zgerc.f   ztbsv.f\
chemm.f   ctpmv.f   drotmg.f  dtrmv.f   sgemv.f   stbmv.f   zgeru.f   ztpmv.f\
chemv.f   ctpsv.f   dsbmv.f   dtrsm.f   sger.f    stbsv.f   zhbmv.f   ztpsv.f\
cher.f    ctrmm.f   dscal.f   dtrsv.f   snrm2.f   stpmv.f   zhemm.f   ztrmm.f\
cher2.f   ctrmv.f   dsdot.f   dzasum.f  srot.f    stpsv.f   zhemv.f   ztrmv.f\
cher2k.f  ctrsm.f   dspmv.f   dznrm2.f  srotg.f   strmm.f   zher.f    ztrsm.f\
cherk.f   ctrsv.f   dspr.f    icamax.f  srotm.f   strmv.f   zher2.f   ztrsv.f\
chpmv.f   dasum.f   dspr2.f   idamax.f  srotmg.f  strsm.f   zher2k.f\
chpr.f    daxpy.f   dswap.f   isamax.f  ssbmv.f   strsv.f   zherk.f\
chpr2.f   dcabs1.f  dsymm.f   izamax.f  sscal.f   xerbla.f  zhpmv.f

# Define name of library
libname = oomph_blas

# Combine headers and sources
headers_and_sources =  $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_blas.la


# Sources that the library depends on:
#-------------------------------------
liboomph_blas_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_blas

#These are the header files that are to be placed in subdirectory: **NONE**
#library_include_HEADERS=$(headers)


# Combined header file **NONE**
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
#$(libname).h: dummy_$(libname).h


# **NONE**
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


#**NONE**
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

# ==============================================================================
# Blas sources
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_blas subdirectory")

# Define the headers
set(HEADERS)

# Define the sources
set(SOURCES
    dummy.cc
    caxpy.f
    crotg.f
    dcopy.f
    dsymv.f
    lsame.f
    sspmv.f
    zaxpy.f
    zhpr.f
    ccopy.f
    cscal.f
    ddot.f
    dsyr.f
    sasum.f
    sspr.f
    zcopy.f
    zhpr2.f
    cdotc.f
    csrot.f
    dgbmv.f
    dsyr2.f
    saxpy.f
    sspr2.f
    zdotc.f
    zrotg.f
    cdotu.f
    csscal.f
    dgemm.f
    dsyr2k.f
    scasum.f
    sswap.f
    zdotu.f
    zscal.f
    cgbmv.f
    cswap.f
    dgemv.f
    dsyrk.f
    scnrm2.f
    ssymm.f
    zdrot.f
    zswap.f
    cgemm.f
    csymm.f
    dger.f
    dtbmv.f
    scopy.f
    ssymv.f
    zdscal.f
    zsymm.f
    cgemv.f
    csyr2k.f
    dnrm2.f
    dtbsv.f
    sdot.f
    ssyr.f
    zgbmv.f
    zsyr2k.f
    cgerc.f
    csyrk.f
    drot.f
    dtpmv.f
    sdsdot.f
    ssyr2.f
    zgemm.f
    zsyrk.f
    cgeru.f
    ctbmv.f
    drotg.f
    dtpsv.f
    sgbmv.f
    ssyr2k.f
    zgemv.f
    ztbmv.f
    chbmv.f
    ctbsv.f
    drotm.f
    dtrmm.f
    sgemm.f
    ssyrk.f
    zgerc.f
    ztbsv.f
    chemm.f
    ctpmv.f
    drotmg.f
    dtrmv.f
    sgemv.f
    stbmv.f
    zgeru.f
    ztpmv.f
    chemv.f
    ctpsv.f
    dsbmv.f
    dtrsm.f
    sger.f
    stbsv.f
    zhbmv.f
    ztpsv.f
    cher.f
    ctrmm.f
    dscal.f
    dtrsv.f
    snrm2.f
    stpmv.f
    zhemm.f
    ztrmm.f
    cher2.f
    ctrmv.f
    dsdot.f
    dzasum.f
    srot.f
    stpsv.f
    zhemv.f
    ztrmv.f
    cher2k.f
    ctrsm.f
    dspmv.f
    dznrm2.f
    srotg.f
    strmm.f
    zher.f
    ztrsm.f
    cherk.f
    ctrsv.f
    dspr.f
    icamax.f
    srotm.f
    strmv.f
    zher2.f
    ztrsv.f
    chpmv.f
    dasum.f
    dspr2.f
    idamax.f
    srotmg.f
    strsm.f
    zher2k.f
    chpr.f
    daxpy.f
    dswap.f
    isamax.f
    ssbmv.f
    strsv.f
    zherk.f
    chpr2.f
    dcabs1.f
    dsymm.f
    izamax.f
    sscal.f
    xerbla.f
    zhpmv.f)

# Define the name of the library
set(LIBNAME oomph_blas)

# Define the library build type
set(LIBTYPE STATIC)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LIBTYPE ${LIBTYPE}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(TARGET ${LIBNAME}
                              FORTRAN_COMPILE_FLAGS -std=legacy)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_blas subdirectory")
# ------------------------------------------------------------------------------

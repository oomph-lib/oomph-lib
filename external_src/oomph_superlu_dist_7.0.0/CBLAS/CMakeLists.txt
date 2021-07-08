set(HEADERS f2c.h)
set(SOURCES input_error_dist.c)

if(ENABLE_DOUBLE)
  list(
    APPEND
    SOURCES
    idamax.c
    dasum.c
    daxpy.c
    dcopy.c
    ddot.c
    dnrm2.c
    drot.c
    dscal.c
    dgemv.c
    dsymv.c
    dtrsv.c
    dger.c
    dsyr2.c
    dgemm.c
    dtrsm.c)
endif()

if(ENABLE_COMPLEX16)
  list(
    APPEND
    SOURCES
    izamax.c
    dzasum.c
    zaxpy.c
    zcopy.c
    dznrm2.c
    zscal.c
    dcabs1.c
    z_internal.c
    zgemv.c
    zhemv.c
    ztrsv.c
    zgerc.c
    zgeru.c
    zher2.c
    zgemm.c
    ztrsm.c)
endif()

# Define the name of the library
set(LIBNAME superlu_dist_cblas)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

/* superlu_dist_config.h.in */

/* Enable CUDA */
/* #undef HAVE_CUDA */

/* Enable HIP */
/* #undef HAVE_HIP */

/* Enable parmetis */
#define HAVE_PARMETIS ON

/* Enable colamd */
#define HAVE_COLAMD ON

/* Enable LAPACK */
#define SLU_HAVE_LAPACK ON

/* Enable CombBLAS */
/* #undef HAVE_COMBBLAS */

/* enable 64bit index mode */
#define XSDK_INDEX_SIZE 32

#if defined(XSDK_INDEX_SIZE) && (XSDK_INDEX_SIZE == 64)
#define _LONGINT 1
#endif

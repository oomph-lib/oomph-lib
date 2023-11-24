
#ifndef SUPERLU_CONFIG_H
#define SUPERLU_CONFIG_H

/* Enable metis */
#define HAVE_METIS ON

/* Enable colamd */
#define HAVE_COLAMD ON

/* enable 64bit index mode */
#define XSDK_INDEX_SIZE 32

/* Integer type for indexing sparse matrix meta structure */
#if (XSDK_INDEX_SIZE == 64)
#include <stdint.h>
#define _LONGINT 1
typedef int64_t int_t;
#else
typedef int int_t; /* default */
#endif

#endif /* SUPERLU_CONFIG_H */


/* Hand-edited version of metis.h -- exludes proto.h */
/* to avoid clashes between C and C++ interfaces */
/* (MH) */


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1 1998/11/27 17:59:21 karypis Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include <metis_4.0/defs.h>
#include <metis_4.0/struct.h>
#include <metis_4.0/macros.h>
#include <metis_4.0/rename.h>
#include <metis_4.0/proto.h>


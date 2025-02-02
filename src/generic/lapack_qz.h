//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//QZ-ehader
// The BLAS DOUBLE dot-product 
// Fortran interface : FUNCTION DDOT(N,X,INCX,Y,INCY) 
//PROTOCCALLSFFUN5( DOUBLE, DDOT, ddot, INT, DOUBLEV, INT, DOUBLEV, INT ) 
//#define BLAS_DDOT(N,X,INCX,Y,INCY) CCALLSFFUN5(DDOT,ddot,INT,DOUBLEV,INT,DOUBLEV,INT,N,X,INCX,Y,INCY) 

// The BLAS DOUBLE matrix multiplier 
// Fortran interface : SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
//PROTOCCALLSFSUB13( DGEMM, dgemm, STRING, STRING, INT, INT, INT, DOUBLE, DOUBLEV, INT, DOUBLEV, INT, DOUBLE, DOUBLEV, INT ) 
//#define BLAS_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) CCALLSFSUB13(DGEMM,dgemm,STRING,STRING,INT,INT,INT,DOUBLE,DOUBLEV,INT,DOUBLEV,INT,DOUBLE,DOUBLEV,INT,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

// The LAPACK LU solver for a DOUBLE DENSE matrix 
// Fortran interface : SUBROUTINE DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO ) 
//PROTOCCALLSFSUB8( DGESV, dgesv, INT, INT, DOUBLEV, INT, INTV, DOUBLEV, INT, PINT ) 
//#define LAPACK_DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO) CCALLSFSUB8(DGESV,dgesv,INT,INT,DOUBLEV,INT,INTV,DOUBLEV,INT,PINT,N,NRHS,A,LDA,IPIV,B,LDB,INFO) 

// The LAPACK LU solver for a DOUBLE BANDED matrix 
//  Fortran interface : SUBROUTINE DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO) 
//PROTOCCALLSFSUB10( DGBSV, dgbsv, INT, INT, INT, INT, DOUBLEV, INT, INTV, DOUBLEV, INT, PINT ) 
//#define LAPACK_DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO) CCALLSFSUB10(DGBSV,dgbsv,INT,INT,INT,INT,DOUBLEV,INT,INTV,DOUBLEV,INT,PINT,N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO) 

// The LAPACK DOUBLE GENERALISED eigenvalue solver 
// Fortran interface : SUBROUTINE DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO) 
PROTOCCALLSFSUB17( DGGEV, dggev, STRING, STRING, INT, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, INT, PINT ) 
#define LAPACK_DGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO) CCALLSFSUB17(DGGEV,dggev,STRING,STRING,INT,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,PINT,JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO) 

// The LAPACK COMPLEX GENERALISE eigenvalue solver 
// Fortran interface : SUBROUTINE ZGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO) 
PROTOCCALLSFSUB17( ZGGEV, zggev, STRING, STRING, INT, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, PINT ) 
// jobvl jobvr n a lda b ldb alpha beta vl ldvl vr ldvr work lwork rwork info 
#define LAPACK_ZGGEV(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO) CCALLSFSUB17(ZGGEV,zggev,STRING,STRING,INT,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,PINT,JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO) 

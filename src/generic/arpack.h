// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
/*cfortran.,h headers for the ARPACK routines*/

PROTOCCALLSFSUB16(DNAUPD,
                  dnaupd,
                  PINT,
                  PSTRING,
                  INT,
                  PSTRING,
                  INT,
                  PDOUBLE,
                  DOUBLEV,
                  INT,
                  DOUBLEVV,
                  INT,
                  INTV,
                  INTV,
                  DOUBLEV,
                  DOUBLEV,
                  INT,
                  PINT)
#define DNAUPD(IDO, \
               BMAT, \
               N, \
               WHICH, \
               NEV, \
               TOL, \
               RESID, \
               NCV, \
               V, \
               LDV, \
               IPARAM, \
               IPNTR, \
               WORKD, \
               WORKL, \
               LWORKL, \
               INFO) \
  CCALLSFSUB16(DNAUPD, \
               dnaupd, \
               PINT, \
               PSTRING, \
               INT, \
               PSTRING, \
               INT, \
               PDOUBLE, \
               DOUBLEV, \
               INT, \
               DOUBLEVV, \
               INT, \
               INTV, \
               INTV, \
               DOUBLEV, \
               DOUBLEV, \
               INT, \
               PINT, \
               IDO, \
               BMAT, \
               N, \
               WHICH, \
               NEV, \
               TOL, \
               RESID, \
               NCV, \
               V, \
               LDV, \
               IPARAM, \
               IPNTR, \
               WORKD, \
               WORKL, \
               LWORKL, \
               INFO)

PROTOCCALLSFSUB25(DNEUPD,
                  dneupd,
                  LOGICAL,
                  PSTRING,
                  LOGICALV,
                  DOUBLEV,
                  DOUBLEV,
                  DOUBLEVV,
                  INT,
                  DOUBLE,
                  DOUBLE,
                  DOUBLEV,
                  PSTRING,
                  INT,
                  PSTRING,
                  INT,
                  DOUBLE,
                  DOUBLEV,
                  INT,
                  DOUBLEVV,
                  INT,
                  INTV,
                  INTV,
                  DOUBLEV,
                  DOUBLEV,
                  INT,
                  PINT)

#define DNEUPD(RVEC, \
               HOWMNY, \
               SELECT, \
               DR, \
               DI, \
               Z, \
               LDZ, \
               SIGMAR, \
               SIGMAI, \
               WORKEV, \
               BMAT, \
               N, \
               WHICH, \
               NEV, \
               TOL, \
               RESID, \
               NCV, \
               V, \
               LDV, \
               IPARAM, \
               IPNTR, \
               WORKD, \
               WORKL, \
               LWORKL, \
               INFO) \
  CCALLSFSUB25(DNEUPD, \
               dneupd, \
               LOGICAL, \
               PSTRING, \
               LOGICALV, \
               DOUBLEV, \
               DOUBLEV, \
               DOUBLEVV, \
               INT, \
               DOUBLE, \
               DOUBLE, \
               DOUBLEV, \
               PSTRING, \
               INT, \
               PSTRING, \
               INT, \
               DOUBLE, \
               DOUBLEV, \
               INT, \
               DOUBLEVV, \
               INT, \
               INTV, \
               INTV, \
               DOUBLEV, \
               DOUBLEV, \
               INT, \
               PINT, \
               RVEC, \
               HOWMNY, \
               SELECT, \
               DR, \
               DI, \
               Z, \
               LDZ, \
               SIGMAR, \
               SIGMAI, \
               WORKEV, \
               BMAT, \
               N, \
               WHICH, \
               NEV, \
               TOL, \
               RESID, \
               NCV, \
               V, \
               LDV, \
               IPARAM, \
               IPNTR, \
               WORKD, \
               WORKL, \
               LWORKL, \
               INFO)

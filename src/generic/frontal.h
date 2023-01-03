// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
/*cfortran.h headers for the ma42 routines*/
PROTOCCALLSFSUB3(MA42ID, ma42id, INTV, DOUBLEV, INTV)
#define MA42ID(ICNTL, CNTL, ISAVE) \
  CCALLSFSUB3(MA42ID, ma42id, INTV, DOUBLEV, INTV, ICNTL, CNTL, ISAVE)

PROTOCCALLSFSUB8(MA42AD, ma42ad, INT, INTV, PINT, INTV, INT, INTV, INTV, INTV)
#define MA42AD(NVAR, IVAR, NDF, LAST, LENLST, ICNTL, ISAVE, INFO) \
  CCALLSFSUB8(MA42AD, \
              ma42ad, \
              INT, \
              INTV, \
              PINT, \
              INTV, \
              INT, \
              INTV, \
              INTV, \
              INTV, \
              NVAR, \
              IVAR, \
              NDF, \
              LAST, \
              LENLST, \
              ICNTL, \
              ISAVE, \
              INFO)

PROTOCCALLSFSUB9(
  MA42JD, ma42jd, INT, INTV, INT, INTV, INT, INTV, INTV, INTV, INTV)
#define MA42JD(NVAR, IVAR, NDF, LAST, NMAXE, IFSIZE, ICNTL, ISAVE, INFO) \
  CCALLSFSUB9(MA42JD, \
              ma42jd, \
              INT, \
              INTV, \
              INT, \
              INTV, \
              INT, \
              INTV, \
              INTV, \
              INTV, \
              INTV, \
              NVAR, \
              IVAR, \
              NDF, \
              LAST, \
              NMAXE, \
              IFSIZE, \
              ICNTL, \
              ISAVE, \
              INFO)

PROTOCCALLSFSUB6(MA42PD, ma42pd, INTV, INTV, INTV, INTV, INTV, INTV)
#define MA42PD(ISTRM, LENBUF, LENFLE, ICNTL, ISAVE, INFO) \
  CCALLSFSUB6(MA42PD, \
              ma42pd, \
              INTV, \
              INTV, \
              INTV, \
              INTV, \
              INTV, \
              INTV, \
              ISTRM, \
              LENBUF, \
              LENFLE, \
              ICNTL, \
              ISAVE, \
              INFO)

PROTOCCALLSFSUB22(MA42BD,
                  ma42bd,
                  INT,
                  INTV,
                  INT,
                  INTV,
                  INT,
                  DOUBLEVV,
                  INT,
                  DOUBLEVV,
                  INT,
                  INT,
                  DOUBLEVV,
                  INTV,
                  INTV,
                  INT,
                  DOUBLEV,
                  INT,
                  INTV,
                  INTV,
                  DOUBLEV,
                  INTV,
                  INTV,
                  DOUBLEV)
#define MA42BD(NVAR, \
               IVAR, \
               NDF, \
               LAST, \
               NMAXE, \
               AVAR, \
               NRHS, \
               RHS, \
               LRHS, \
               LX, \
               X, \
               NFRONT, \
               LENBUF, \
               LW, \
               W, \
               LIW, \
               IW, \
               ICNTL, \
               CNTL, \
               ISAVE, \
               INFO, \
               RINFO) \
  CCALLSFSUB22(MA42BD, \
               ma42bd, \
               INT, \
               INTV, \
               INT, \
               INTV, \
               INT, \
               DOUBLEVV, \
               INT, \
               DOUBLEVV, \
               INT, \
               INT, \
               DOUBLEVV, \
               INTV, \
               INTV, \
               INT, \
               DOUBLEV, \
               INT, \
               INTV, \
               INTV, \
               DOUBLEV, \
               INTV, \
               INTV, \
               DOUBLEV, \
               NVAR, \
               IVAR, \
               NDF, \
               LAST, \
               NMAXE, \
               AVAR, \
               NRHS, \
               RHS, \
               LRHS, \
               LX, \
               X, \
               NFRONT, \
               LENBUF, \
               LW, \
               W, \
               LIW, \
               IW, \
               ICNTL, \
               CNTL, \
               ISAVE, \
               INFO, \
               RINFO)

PROTOCCALLSFSUB12(MA42CD,
                  ma42cd,
                  LOGICAL,
                  INT,
                  INT,
                  DOUBLEVV,
                  DOUBLEVV,
                  INT,
                  DOUBLEV,
                  INT,
                  INTV,
                  INTV,
                  INTV,
                  INTV)
#define MA42CD(TRANS, NRHS, LX, B, X, LW, W, LIW, IW, ICNTL, ISAVE, INFO) \
  CCALLSFSUB12(MA42CD, \
               ma42cd, \
               LOGICAL, \
               INT, \
               INT, \
               DOUBLEVV, \
               DOUBLEVV, \
               INT, \
               DOUBLEV, \
               INT, \
               INTV, \
               INTV, \
               INTV, \
               INTV, \
               TRANS, \
               NRHS, \
               LX, \
               B, \
               X, \
               LW, \
               W, \
               LIW, \
               IW, \
               ICNTL, \
               ISAVE, \
               INFO)

PROTOCCALLSFSUB1(MC63ID, mc63id, INTV)
#define MC63ID(ICNTL) CCALLSFSUB1(MC63ID, mc63id, INTV, ICNTL)

PROTOCCALLSFSUB19(MC63AD,
                  mc63ad,
                  LOGICAL,
                  INT,
                  INT,
                  INT,
                  INTV,
                  INTV,
                  INTV,
                  INTV,
                  PINT,
                  INTV,
                  INTV,
                  DOUBLEV,
                  INT,
                  INTV,
                  INT,
                  DOUBLEV,
                  INTV,
                  INTV,
                  DOUBLEV)
#define MC63AD(DIRECT, \
               N, \
               NELT, \
               NE, \
               ELTVAR, \
               ELTPTR, \
               ORDER, \
               PERM, \
               NSUP, \
               VARS, \
               SVAR, \
               WT, \
               LIW, \
               IW, \
               LW, \
               W, \
               ICNTL, \
               INFO, \
               RINFO) \
  CCALLSFSUB19(MC63AD, \
               mc63ad, \
               LOGICAL, \
               INT, \
               INT, \
               INT, \
               INTV, \
               INTV, \
               INTV, \
               INTV, \
               PINT, \
               INTV, \
               INTV, \
               DOUBLEV, \
               INT, \
               INTV, \
               INT, \
               DOUBLEV, \
               INTV, \
               INTV, \
               DOUBLEV, \
               DIRECT, \
               N, \
               NELT, \
               NE, \
               ELTVAR, \
               ELTPTR, \
               ORDER, \
               PERM, \
               NSUP, \
               VARS, \
               SVAR, \
               WT, \
               LIW, \
               IW, \
               LW, \
               W, \
               ICNTL, \
               INFO, \
               RINFO)

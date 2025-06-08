// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
/*cfortran.h headers for the mumps routines*/


PROTOCCALLSFSUB1(MUMPS_SETUP_SOLVER_POOL, mumps_setup_solver_pool, INT)
#define mumps_setup_solver_pool(N_POOL) \
  CCALLSFSUB1(MUMPS_SETUP_SOLVER_POOL, mumps_setup_solver_pool, INT, N_POOL)

PROTOCCALLSFSUB2(MUMPS_SETUP, mumps_setup, INT, INT)
#define mumps_setup(I_POOL, S) \
  CCALLSFSUB2(MUMPS_SETUP, mumps_setup, INT, INT, I_POOL, S)


PROTOCCALLSFSUB2(MUMPS_SET_WORKSPACE_SCALING_FACTOR,
                 mumps_set_workspace_scaling_factor,
                 INT,
                 INT)
#define mumps_set_workspace_scaling_factor(I_POOL, S) \
  CCALLSFSUB2(MUMPS_SET_WORKSPACE_SCALING_FACTOR, \
              mumps_set_workspace_scaling_factor, \
              INT, \
              INT, \
              I_POOL, \
              S)


PROTOCCALLSFSUB1(MUMPS_SWITCH_ON_DOC, mumps_switch_on_doc, INT)
#define mumps_switch_on_doc(I_POOL) \
  CCALLSFSUB1(MUMPS_SWITCH_ON_DOC, mumps_switch_on_doc, INT, I_POOL)


PROTOCCALLSFSUB1(MUMPS_SWITCH_OFF_DOC, mumps_switch_off_doc, INT)
#define mumps_switch_off_doc(I_POOL) \
  CCALLSFSUB1(MUMPS_SWITCH_OFF_DOC, mumps_switch_off_doc, INT, I_POOL)


PROTOCCALLSFSUB6(
  MUMPS_FACTORISE, mumps_factorise, INT, INT, INT, INTV, INTV, DOUBLEV)
#define mumps_factorise(I_POOL, N, NZ_LOC, IRN_LOC, JCN_LOC, A_LOC) \
  CCALLSFSUB6(MUMPS_FACTORISE, \
              mumps_factorise, \
              INT, \
              INT, \
              INT, \
              INTV, \
              INTV, \
              DOUBLEV, \
              I_POOL, \
              N, \
              NZ_LOC, \
              IRN_LOC, \
              JCN_LOC, \
              A_LOC)

PROTOCCALLSFSUB7(
  MUMPS_SOLVE, mumps_solve, INT, INT, INT, INTV, INTV, DOUBLEV, DOUBLEV)
#define mumps_solve(I_POOL, N, NZ_LOC, IRN_LOC, JCN_LOC, A_LOC, RHS) \
  CCALLSFSUB7(MUMPS_SOLVE, \
              mumps_solve, \
              INT, \
              INT, \
              INT, \
              INTV, \
              INTV, \
              DOUBLEV, \
              DOUBLEV, \
              I_POOL, \
              N, \
              NZ_LOC, \
              IRN_LOC, \
              JCN_LOC, \
              A_LOC, \
              RHS)


PROTOCCALLSFSUB3(MUMPS_BACKSUB, mumps_backsub, INT, INT, DOUBLEV)
#define mumps_backsub(I_POOL, N, RHS) \
  CCALLSFSUB3(MUMPS_BACKSUB, mumps_backsub, INT, INT, DOUBLEV, I_POOL, N, RHS)


PROTOCCALLSFSUB1(MUMPS_CLEANUP_MEMORY, mumps_cleanup_memory, INT)
#define mumps_cleanup_memory(I_POOL) \
  CCALLSFSUB1(MUMPS_CLEANUP_MEMORY, mumps_cleanup_memory, INT, I_POOL)


PROTOCCALLSFSUB1(MUMPS_SHUTDOWN, mumps_shutdown, INT)
#define mumps_shutdown(I_POOL) \
  CCALLSFSUB1(MUMPS_SHUTDOWN, mumps_shutdown, INT, I_POOL)

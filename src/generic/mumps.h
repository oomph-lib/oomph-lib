//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
/*cfortran.h headers for the mumps routines*/


PROTOCCALLSFSUB0(MUMPS_SETUP,mumps_setup)
#define mumps_setup() CCALLSFSUB0(MUMPS_SETUP,mumps_setup)


PROTOCCALLSFSUB1(MUMPS_SET_WORKSPACE_SCALING_FACTOR,mumps_set_workspace_scaling_factor,INT)
#define mumps_set_workspace_scaling_factor(S) CCALLSFSUB1(MUMPS_SET_WORKSPACE_SCALING_FACTOR,mumps_set_workspace_scaling_factor,INT,S)


PROTOCCALLSFSUB0(MUMPS_SWITCH_ON_DOC,mumps_switch_on_doc)
#define mumps_switch_on_doc() CCALLSFSUB0(MUMPS_SWITCH_ON_DOC,mumps_switch_on_doc)


PROTOCCALLSFSUB0(MUMPS_SWITCH_OFF_DOC,mumps_switch_off_doc)
#define mumps_switch_off_doc() CCALLSFSUB0(MUMPS_SWITCH_OFF_DOC,mumps_switch_off_doc)


PROTOCCALLSFSUB5(MUMPS_FACTORISE,mumps_factorise,INT,INT,INTV,INTV,DOUBLEV)
#define mumps_factorise(N,NZ_LOC,IRN_LOC,JCN_LOC,A_LOC) CCALLSFSUB5(MUMPS_FACTORISE,mumps_factorise,INT,INT,INTV,INTV,DOUBLEV,N,NZ_LOC,IRN_LOC,JCN_LOC,A_LOC)

PROTOCCALLSFSUB6(MUMPS_SOLVE,mumps_solve,INT,INT,INTV,INTV,DOUBLEV,DOUBLEV)
#define mumps_solve(N,NZ_LOC,IRN_LOC,JCN_LOC,A_LOC,RHS) CCALLSFSUB6(MUMPS_SOLVE,mumps_solve,INT,INT,INTV,INTV,DOUBLEV,DOUBLEV,N,NZ_LOC,IRN_LOC,JCN_LOC,A_LOC,RHS)


PROTOCCALLSFSUB2(MUMPS_BACKSUB,mumps_backsub,INT,DOUBLEV)
#define mumps_backsub(N,RHS) CCALLSFSUB2(MUMPS_BACKSUB,mumps_backsub,INT,DOUBLEV,N,RHS)


PROTOCCALLSFSUB0(MUMPS_CLEANUP_MEMORY,mumps_cleanup_memory)
#define mumps_cleanup_memory() CCALLSFSUB0(MUMPS_CLEANUP_MEMORY,mumps_cleanup_memory)


PROTOCCALLSFSUB0(MUMPS_SHUTDOWN,mumps_shutdown)
#define mumps_shutdown() CCALLSFSUB0(MUMPS_SHUTDOWN,mumps_shutdown)


PROTOCCALLSFSUB0(MUMPS_TEST_SOLVE,mumps_test_solve)
#define mumps_test_solve() CCALLSFSUB0(MUMPS_TEST_SOLVE,mumps_test_solve)

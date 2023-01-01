//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#include <iostream>

// Test problem to generate sparse matrices for testing
// Driver for 2D rectangular driven cavity

// Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;
using namespace oomph;


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//==start_of_namespace==============================
/// Namespace for parameters
//==================================================
namespace Global_Variables
{

 /// Reynolds number
 double Re=50;

 // Number of elements in x direction
 unsigned Nx=4;

 // Number of elements in x direction
 unsigned Ny=4;

} // end_of_namespace


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class RectangularDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 RectangularDrivenCavityProblem();

 /// Destructor to clean up memory
 ~RectangularDrivenCavityProblem();

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof,
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end of fix_pressure

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve.
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  {
   // Setup tangential flow along boundary 0:
   unsigned ibound=0;
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Tangential flow
     unsigned i=0;
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
     // No penetration
     i=1;
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
    }

   // Overwrite with no flow along the other boundaries
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=1;ibound<num_bound;ibound++)
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
        }
      }
    }
  } // end_of_actions_before_newton_solve

 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt()
  {
   // Upcast from pointer to the Mesh base class to the specific
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RectangularDrivenCavity problem
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::RectangularDrivenCavityProblem()
{
 // Setup mesh

 // # of elements in x-direction
 unsigned Nx = Global_Variables::Nx;

 // # of elements in y-direction
 unsigned Ny = Global_Variables::Ny;

 // Domain length in x-direction
 double Lx=1.0;

 // Domain length in y-direction
 double Ly=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
      }
    }
  } // end loop over boundaries

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Variables::Re;
  } // end loop over elements

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;
} // end_of_constructor


//==start_of_destructor===================================================
/// Destructor for RectangularDrivenCavity problem
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::~RectangularDrivenCavityProblem()
{
 // Mesh gets killed in general problem destructor
} // end_of_destructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solution
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_main======================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element.
//=====================================================================
int main()
{
 // build the problems and
 typedef QCrouzeixRaviartElement<2> ELEMENT;
 
 Global_Variables::Nx=5;
 Global_Variables::Ny=5;
 RectangularDrivenCavityProblem<ELEMENT>* problem1_pt
  = new RectangularDrivenCavityProblem<ELEMENT>;
 problem1_pt->actions_before_newton_solve();

 Global_Variables::Nx=5;
 Global_Variables::Ny=5;
 RectangularDrivenCavityProblem<ELEMENT>* problem2_pt
   	 = new RectangularDrivenCavityProblem<ELEMENT>;
 problem2_pt->actions_before_newton_solve();
 
 // get Jacobians
 DoubleVector residual;
 CRDoubleMatrix* CR_matrix1_pt = new CRDoubleMatrix;
 CRDoubleMatrix* CR_matrix2_pt = new CRDoubleMatrix;
 CCDoubleMatrix* CC_matrix1_pt = new CCDoubleMatrix;
 CCDoubleMatrix* CC_matrix2_pt = new CCDoubleMatrix;
 problem1_pt->get_jacobian(residual,*CR_matrix1_pt);
 residual.clear();
 problem1_pt->get_jacobian(residual,*CR_matrix2_pt);
 residual.clear();
 problem1_pt->get_jacobian(residual,*CC_matrix1_pt);
 residual.clear();
 problem1_pt->get_jacobian(residual,*CC_matrix2_pt);
 residual.clear();

 DenseDoubleMatrix* D_matrix1_pt = new DenseDoubleMatrix;
 DenseDoubleMatrix* D_matrix2_pt = new DenseDoubleMatrix;
 problem2_pt->get_jacobian(residual,*D_matrix1_pt);
 problem2_pt->get_jacobian(residual,*D_matrix2_pt);

	// test matrix-matrix multiplication
	CRDoubleMatrix* CR_result1_pt = new CRDoubleMatrix;
 CRDoubleMatrix* CR_result2_pt = new CRDoubleMatrix;
 CRDoubleMatrix* CR_result3_pt = new CRDoubleMatrix;
 CCDoubleMatrix* CC_result1_pt = new CCDoubleMatrix;
 CCDoubleMatrix* CC_result2_pt = new CCDoubleMatrix;
 CCDoubleMatrix* CC_result3_pt = new CCDoubleMatrix;
 DenseDoubleMatrix* D_result_pt = new DenseDoubleMatrix;
 
 cout << "Running through matrix-matrix multiplication methods" << endl;
 cout << "----------------------------------------------------" << endl;
 
 cout << "Multiplying " << CR_matrix1_pt->ncol()
 					<< " x " << CR_matrix1_pt->nrow()
      << " matrix by a " << CR_matrix2_pt->ncol()
      << " x " << CR_matrix1_pt->nrow()
      << " matrix" << endl;
 
 cout << "CRDoubleMatrix method 1 " << endl;
 CR_matrix1_pt->serial_matrix_matrix_multiply_method()=1;
 CR_matrix1_pt->multiply(*CR_matrix2_pt, *CR_result1_pt);

 cout << "CRDoubleMatrix method 2 " << endl;
 CR_matrix1_pt->serial_matrix_matrix_multiply_method()=2;
 CR_matrix1_pt->multiply(*CR_matrix2_pt, *CR_result2_pt);

 cout << "CRDoubleMatrix method 3 " << endl;
 CR_matrix1_pt->serial_matrix_matrix_multiply_method()=3;
 CR_matrix1_pt->multiply(*CR_matrix2_pt, *CR_result3_pt);

 cout << "CCDoubleMatrix method 1 " << endl;
 CC_matrix1_pt->matrix_matrix_multiply_method()=1;
 CC_matrix1_pt->multiply(*CC_matrix2_pt, *CC_result1_pt);

 cout << "CCDoubleMatrix method 2 " << endl;
 CC_matrix1_pt->matrix_matrix_multiply_method()=2;
 CC_matrix1_pt->multiply(*CC_matrix2_pt, *CC_result2_pt);
 
 cout << "CCDoubleMatrix method 3 " << endl;
 CC_matrix1_pt->matrix_matrix_multiply_method()=3;
 CC_matrix1_pt->multiply(*CC_matrix2_pt, *CC_result3_pt);

 cout << "Multiplying " << D_matrix1_pt->ncol()
 					<< " x " << D_matrix1_pt->nrow()
      << " matrix by a " << D_matrix2_pt->ncol()
      << " x " << D_matrix1_pt->nrow()
      << " matrix" << endl;

 cout << "DenseDoubleMatrix " << endl;
 D_matrix1_pt->multiply(*D_matrix2_pt, *D_result_pt);

 // output matrices
 bool output = true;
 if (output)
 {
  ofstream file;
  file.open("CR_matrix1.dat", ios_base::out);
  CR_matrix1_pt->sparse_indexed_output(file);
  file.close();
  file.open("CR_matrix2.dat", ios_base::out);
  CR_matrix2_pt->sparse_indexed_output(file);
  file.close();

  file.open("CR_result1.dat", ios_base::out);
  CR_result1_pt->sparse_indexed_output(file);
  file.close();
  file.open("CR_result2.dat", ios_base::out);
  CR_result2_pt->sparse_indexed_output(file);
  file.close();
  file.open("CR_result3.dat", ios_base::out);
  CR_result3_pt->sparse_indexed_output(file);
  file.close();

  file.open("CC_matrix1.dat", ios_base::out);
  CC_matrix1_pt->sparse_indexed_output(file);
  file.close();
  file.open("CC_matrix2.dat", ios_base::out);
  CR_matrix2_pt->sparse_indexed_output(file);
  file.close();

  file.open("CC_result1.dat", ios_base::out);
  CC_result1_pt->sparse_indexed_output(file);
  file.close();
  file.open("CC_result2.dat", ios_base::out);
  CC_result2_pt->sparse_indexed_output(file);
  file.close();
  file.open("CC_result3.dat", ios_base::out);
  CC_result3_pt->sparse_indexed_output(file);
  file.close();

  file.open("D_matrix1.dat", ios_base::out);
  D_matrix1_pt->sparse_indexed_output(file);
  file.close();
  file.open("D_matrix2.dat", ios_base::out);
  D_matrix2_pt->sparse_indexed_output(file);
  file.close();
  file.open("D_result.dat", ios_base::out);
  D_result_pt->sparse_indexed_output(file);
  file.close();
 }

 // clear up memory
 delete problem1_pt;
 delete problem2_pt;

 delete CR_matrix1_pt;
 delete CR_matrix2_pt;
 delete CR_result1_pt;
 delete CR_result2_pt;
 delete CR_result3_pt;

 delete CC_matrix1_pt;
 delete CC_matrix2_pt;
 delete CC_result1_pt;
 delete CC_result2_pt;
 delete CC_result3_pt;

 delete D_matrix1_pt;
 delete D_matrix2_pt;
 delete D_result_pt;

} // end_of_main




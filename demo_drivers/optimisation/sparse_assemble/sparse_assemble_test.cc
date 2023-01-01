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
// Driver to compare different sparse matrix assembly strategies

#include<time.h>

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
 double Re=100;

 /// Number of elements in x direction
 unsigned Nx=4;
 
 /// Number of elements in x direction
 unsigned Ny=4;

 /// By default we're dumping the matrices for comparison
 bool Dump_matrices=true;

 /// Halt code execution at the end of the most memory-intensive 
 /// phase?
 bool Halt_code=false;

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

 /// Destructor (empty)
 ~RectangularDrivenCavityProblem()
  {
   delete Problem::mesh_pt();
  }

 /// To be killed
 void sparse_assemble_row_or_column_compressed_test(
  Vector<Vector<int> > &column_or_row_index,
  Vector<Vector<int> > &row_or_column_start,
  Vector<Vector<double> > &value,
  Vector<Vector<double>*> &residuals,
  bool compressed_row_flag,
  unsigned method=1);


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
 

 /// Do what it says...
 void compare_assembly_strategies(const unsigned& method);

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





//========================================================================
/// Compare assembly strategies
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::compare_assembly_strategies(
 const unsigned& method)
{

 // Halt code execution?
 Problem::Pause_at_end_of_sparse_assembly=Global_Variables::Halt_code;

 DoubleVector residual;
 CRDoubleMatrix matrix;
  
 Vector<int*> col_index(1);
 Vector<int* > row_start(1);
 Vector<double* > value(1);
 Vector<unsigned> nnz(1);
 Vector<double* > residuals(1);
 bool compressed_row = true;
 
 clock_t clock1 = clock();
// time_t time1 = std::time(0);
 
 switch(method)
  {
  case 1:
     
   Problem::Sparse_assembly_method=
    Problem::Perform_assembly_using_vectors_of_pairs;
     
   std::cout << std::endl << "Oomph vector of pairs" 
             << std::endl << std::endl;
   break;
     
     
  case 2:
     
   Problem::Sparse_assembly_method=
    Problem::Perform_assembly_using_two_vectors;
     

   std::cout << std::endl << "Oomph two vectors" 
             << std::endl << std::endl;
   break;
     
  case 3:
     
   Problem::Sparse_assembly_method=
    Problem::Perform_assembly_using_maps;
     

   std::cout << std::endl << "Oomph maps" 
             << std::endl << std::endl;
   break;
     
     
  case 4:
     
   Problem::Sparse_assembly_method=
    Problem::Perform_assembly_using_lists;


   std::cout << std::endl << "Oomph lists" 
             << std::endl << std::endl;
     
   break;


  case 5:
     
   Problem::Sparse_assembly_method=
    Problem::Perform_assembly_using_two_arrays;


   std::cout << std::endl << "Oomph arrays" 
             << std::endl << std::endl;
     
   break;
     
  default:
   
   std::cout << "Never get here" << std::endl;
   exit(1);
   
  }

 // Call assembly routine
 sparse_assemble_row_or_column_compressed(col_index,
                                          row_start,
                                          value,
                                          nnz,
                                          residuals,
                                          compressed_row);
  
 clock_t clock2 = clock();
 //time_t time2 = time(0);
 
 cout <<"CPU time for " 
      << Global_Variables::Nx << " x " 
      << Global_Variables::Ny << " mesh; assembly with method " 
      << method << ": "
      << double(clock2-clock1)/CLOCKS_PER_SEC << "s" << endl;
 
// cout << "Wall time for method " << method << ": "
//       << difftime(time2, time1) << "s" << endl;

//Read out the dimension of the matrix 
 unsigned long n = ndof();
 LinearAlgebraDistribution dist(this->communicator_pt(),n,false);
 matrix.build(&dist);
 matrix.build_without_copy(n,nnz[0],value[0],col_index[0],row_start[0]);
 
 /// Dump matrix?
 if (Global_Variables::Dump_matrices)
  {
   ofstream matrix_file;
   char filename[100];
   sprintf(filename,"matrix%i.dat",method);
   matrix_file.open(filename);
   matrix.sparse_indexed_output(matrix_file);
   matrix_file.close();
  }

 // clean up
 delete[] residuals[0];
}


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_main======================================================
/// Driver to compare different assembly strategies
//=====================================================================
int main(int argc, char* argv[])
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 cout << "Testing sparse assembly" << endl;
 cout << "-----------------------" << endl;

 // Loop over all assembly strategies or just do one?
 unsigned lo_method=1;
 unsigned hi_method=8;

 if (CommandLineArgs::Argc==1)
  {
   std::cout 
    << "No command line args specified -- running all assembly strategies"
    << std::endl;
  }
 else if (CommandLineArgs::Argc==6)
  {
   Global_Variables::Nx=atoi(argv[1]);
   Global_Variables::Ny=atoi(argv[1]);
   lo_method=atoi(argv[3]);
   hi_method=atoi(argv[3]);
   std::cout << "Assembly method: " << hi_method << std::endl;
   Global_Variables::Dump_matrices=atoi(argv[4]);
   Global_Variables::Halt_code=atoi(argv[5]);
  }
 else
  {
   std::cout 
    << "Wrong number of command line args specified!\n" 
    << "Specify none or five: \n"
    << "- Number of elements in x-direction\n"
    << "- Number of elements in y-direction\n"
    << "- assembly method [1-4]\n" 
    << "- output of matrices [0/1]\n" 
    << "- halt code to analyse memory usage [0/1]\n" 
    << std::endl;
   exit(1);
  }


 cout << "Mesh with " 
      << Global_Variables::Nx << " x " 
      << Global_Variables::Ny << " elements."
      << std::endl;
 
 for (unsigned method = lo_method; method<=hi_method; method++)
  {
   // Build the problem with QTaylorHoodElements
   RectangularDrivenCavityProblem<QTaylorHoodElement<2> > problem;

   // Perform actions before solve to get the problem ready
   problem.actions_before_newton_solve();

   // Run with assembly with specified method
   problem.compare_assembly_strategies(method);

  }

} // end_of_main




//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.
// This one is used to check the parallel behaviour of the line visualiser.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;
} // end_of_namespace



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class RefineableDrivenCavityProblem : public Problem
{
private:
 /// Line visualiser pointer
 LineVisualiser* LV_pt;

public:

 /// Constructor
 RefineableDrivenCavityProblem();

 /// Destructor: Empty
 ~RefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
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
  
  // Overwrite with no flow along all other boundaries
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

 /// Before adaptation
 void actions_before_adapt()
  {
   // Kill line visualiser.
   delete LV_pt;
   LV_pt=0;
  }


 /// After adaptation: Unpin pressure and pin redundant pressure dofs.
 void actions_after_adapt()
  {
   // Reincarnate line visualiser.
   
   // How many points to you want?
   unsigned int  npt=100;
   Vector<Vector<double> > coord_vec(npt);
   coord_vec[0].resize(2);
   coord_vec[0][0]=0.1;
   coord_vec[0][1]=0.9;
   for (unsigned j=1;j<npt;j++)
    {
     coord_vec[j].resize(2);
     coord_vec[j][0]=double(j)/100.0;
     coord_vec[j][1]=0.9/double(j);
    }

   // Setup line visualiser
   LV_pt=new LineVisualiser(Problem::mesh_pt(),
                            coord_vec);


   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0

   // Loop over all elements
   const unsigned n_element=mesh_pt()->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     // If the lower left node of this element is (0,0), then fix the 
     // pressure dof in this element to zero
     if (mesh_pt()->finite_element_pt(e)->node_pt(0)->x(0)==0.0 && 
         mesh_pt()->finite_element_pt(e)->node_pt(0)->x(1)==0.0) // 2d problem
      {
       // Fixing the pressure.
       // Fix the pressure in element e at pdof=0 to 0.0
       unsigned pdof=0;
       fix_pressure(e,pdof,0.0);
      }
    }



  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 
private:

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
RefineableDrivenCavityProblem<ELEMENT>::RefineableDrivenCavityProblem()
{ 
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
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

 //Find number of elements in mesh
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the first pressure dof in the first element to 0.0
 fix_pressure(0,0,0.0);
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
 // Line visualiser data allocation

 // How many points to you want?
 unsigned int  npt=100;
 Vector<Vector<double> > coord_vec(npt);
  coord_vec[0].resize(2);
  coord_vec[0][0]=0.1;
  coord_vec[0][1]=0.9;
 for (unsigned j=1;j<npt;j++)
  {
   coord_vec[j].resize(2);
   coord_vec[j][0]=double(j)/100.0;
   coord_vec[j][1]=0.9/double(j);
  }

 // Setup line visualiser
 LV_pt=
  new LineVisualiser(Problem::mesh_pt(),
                     coord_vec);


} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Output solution, only output Line Visualsier as it's all we're 
 // interested
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 LV_pt->output(some_file);
 some_file.close();
} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main(int argc, char **argv)
{

#ifdef OOMPH_HAS_MPI

 // Initialise MPI
 MPI_Helpers::init(argc,argv);

#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Set max. number of black-box adaptation
 unsigned max_adapt=3;

 // Solve problem with Taylor Hood elements
 //---------------------------------------
 {
  //Build problem
  RefineableDrivenCavityProblem<RefineableQTaylorHoodElement<2> > problem;

  //Are there command-line arguments?
  if (CommandLineArgs::Argc==1)
   {

#ifdef OOMPH_HAS_MPI

    // Provide storage for each element's partition number
    const unsigned n_element=problem.mesh_pt()->nelement();
    Vector<unsigned> out_element_partition(n_element);

    // Distribute the problem
    bool report_stats=false;
    out_element_partition=problem.distribute(report_stats);

    // Write partition to disk
    std::ofstream output_file;
    char filename[100];
    sprintf(filename,"out_adaptive_cavity_1_partition.dat");
    output_file.open(filename);
    for (unsigned e=0;e<n_element;e++)
     {
      output_file << out_element_partition[e] << std::endl;
     }

    // Check halo schemes (optional)
    problem.check_halo_schemes(doc_info);

#endif

    // Solve the problem with automatic adaptation
    problem.newton_solve(max_adapt);
  
    //Output solution
    problem.doc_solution(doc_info);

   } 
  // Validation run - read in partition from file
  else 
   {

#ifdef OOMPH_HAS_MPI

    // DocInfo object specifies directory in which we document
    // the distribution
    DocInfo mesh_doc_info;
    mesh_doc_info.set_directory("RESLT");

    // Create storage for pre-determined partitioning
    const unsigned n_element=problem.mesh_pt()->nelement();
    Vector<unsigned> element_partition(n_element);

    // Read in partitioning from disk
    std::ifstream input_file;
    char filename[100];
    sprintf(filename,"adaptive_cavity_1_partition.dat");
    input_file.open(filename);
    std::string input_string;
    for (unsigned e=0;e<n_element;e++)
     {
      getline(input_file,input_string,'\n');
      element_partition[e]=atoi(input_string.c_str());
     }

    // Now perform the distribution
    problem.distribute(element_partition);

#endif

    // solve with adaptation
    problem.newton_solve(max_adapt);

    //Output solution
    problem.doc_solution(doc_info);
   }

 } // end of Taylor Hood elements
 

 // Solve problem with Crouzeix Raviart elements
 //--------------------------------------------
 {
  // Build problem
  RefineableDrivenCavityProblem<RefineableQCrouzeixRaviartElement<2> > problem;

  //Are there command-line arguments?
  if (CommandLineArgs::Argc==1)
   {
#ifdef OOMPH_HAS_MPI
    // Distribute the problem
    problem.distribute();

    // Check halo schemes (optional)
    problem.check_halo_schemes(doc_info);
#endif

    // Solve the problem with automatic adaptation
    problem.newton_solve(max_adapt);
  
    // Step number
    doc_info.number()=1;
   
    //Output solution
    problem.doc_solution(doc_info);
   } // end of no command-line arguments
  else // Validation run - read in partition from file
   {

#ifdef OOMPH_HAS_MPI

    DocInfo mesh_doc_info;
    mesh_doc_info.set_directory("RESLT");
    mesh_doc_info.number()=0;
    std::ifstream input_file;
    char filename[100];

    // Get the partition to be used from file
    const unsigned n_element=problem.mesh_pt()->nelement();
    Vector<unsigned> element_partition(n_element);
    sprintf(filename,"adaptive_cavity_2_partition.dat");
    input_file.open(filename);
    std::string input_string;
    for (unsigned e=0;e<n_element;e++)
     {
      getline(input_file,input_string,'\n');
      element_partition[e]=atoi(input_string.c_str());
     }

    // Now perform the distribution
    problem.distribute(element_partition); 

#endif

    // Re-solve with adaptation
    problem.newton_solve(max_adapt);

    // change doc_info number
    doc_info.number()=1;

    //Output solution
    problem.doc_solution(doc_info);
   }

 } // end of Crouzeix Raviart elements


// Finalise MPI
#ifdef OOMPH_HAS_MPI

 MPI_Helpers::finalize();

#endif

} // end_of_main


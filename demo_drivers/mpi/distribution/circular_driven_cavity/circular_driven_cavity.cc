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
// Driver for adaptive 2D quarter circle driven cavity. Solved with black
// box adaptation, using Taylor Hood elements (2 problems).

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;


//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;

 /// Reynolds/Froude number
 double Re_invFr=100;

 /// Gravity vector
 Vector<double> Gravity(2);

 /// Functional body force
 void body_force(const double& time, const Vector<double>& x, 
                 Vector<double>& result)
 {
  result[0]=0.0;
  result[1]=-Re_invFr;
 }

 /// Zero functional body force
 void zero_body_force(const double& time, const Vector<double>& x, 
                      Vector<double>& result)
 {
  result[0]=0.0;
  result[1]=0.0;
 }

} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//==start_of_problem_class============================================
/// Driven cavity problem in quarter circle domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class QuarterCircleDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 QuarterCircleDrivenCavityProblem(
  NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt);

 /// Destructor: Empty -- all memory gets cleaned up in base destructor
 ~QuarterCircleDrivenCavityProblem() {}

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


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());

   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now pin the first pressure dof in the first element and set it to 0.0

   // Loop over all elements
   unsigned n_element=mesh_pt()->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     // If the lower left node of this element is (0,0), then fix the 
     // pressure dof in this element to zero
     if (mesh_pt()->finite_element_pt(e)->node_pt(0)->x(0)==0.0 && 
         mesh_pt()->finite_element_pt(e)->node_pt(0)->x(1)==0.0) // 2d problem
      {
       oomph_info << "I'm fixing the pressure " << std::endl;
       // Fix the pressure in element e at pdof=0 to 0.0
       unsigned pdof=0;
       fix_pressure(e,pdof,0.0);
      }
    }

  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Pointer to body force function
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt Body_force_fct_pt;

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for driven cavity problem in quarter circle domain
//========================================================================
template<class ELEMENT>
QuarterCircleDrivenCavityProblem<ELEMENT>::QuarterCircleDrivenCavityProblem(
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt) :
 Body_force_fct_pt(body_force_fct_pt)
{ 

 // Build geometric object that parametrises the curved boundary
 // of the domain

 // Half axes for ellipse
 double a_ellipse=1.0;
 double b_ellipse=1.0;

 // Setup elliptical ring 
 GeomObject* Wall_pt=new Ellipse(a_ellipse,b_ellipse); 

 // End points for wall
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 //Now create the mesh
 double fract_mid=0.5;
 Problem::mesh_pt() = new 
  RefineableQuarterCircleSectorMesh<ELEMENT>(
   Wall_pt,xi_lo,fract_mid,xi_hi);
 
 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
  mesh_pt())->spatial_error_estimator_pt()=error_estimator_pt;
 
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
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   //Set the Re/Fr
   el_pt->re_invfr_pt() = &Global_Physical_Variables::Re_invFr;
   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Gravity;
   //set body force function
   el_pt->body_force_fct_pt() = Body_force_fct_pt;

  } // end loop over elements
 
 // Initial refinement level
 refine_uniformly();
 refine_uniformly();

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now pin the first pressure dof in the first element and set it to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution



//==start_of_main======================================================
/// Driver for QuarterCircleDrivenCavityProblem test problem 
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Set output directory and initialise count
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 DocInfo mesh_doc_info;
 mesh_doc_info.set_directory("RESLT_MESH");
 bool report_stats=false;
 DocInfo error_doc_info;
 error_doc_info.set_directory("RESLT_ERROR");

 // Set max. number of black-box adaptation
 unsigned max_adapt=3; 

 // Solve problem 1 with Taylor-Hood elements
 //--------------------------------------------
 {
  oomph_info << " " << std::endl
             << "-----------------------" << std::endl
             << "-- Running problem 1 --" << std::endl
             << "-----------------------" << std::endl;

  // Set up downwards-Gravity vector
  Global_Physical_Variables::Gravity[0] = 0.0;
  Global_Physical_Variables::Gravity[1] = -1.0;
  
  // Set up Gamma vector for stress-divergence form
  NavierStokesEquations<2>::Gamma[0]=1;
  NavierStokesEquations<2>::Gamma[1]=1;

  // Build problem with Gravity vector in stress divergence form, 
  // using zero body force function
  QuarterCircleDrivenCavityProblem<RefineableQTaylorHoodElement<2> > 
   problem(&Global_Physical_Variables::zero_body_force);

  // output error estimates during solver
  error_doc_info.number()=0;

  RefineableMeshBase* mmesh_pt = 
   dynamic_cast<RefineableMeshBase*>(problem.mesh_pt(0));

  mmesh_pt->doc_info_pt()=&error_doc_info;

  // Initial solve
  problem.newton_solve();

#ifdef OOMPH_HAS_MPI
  mesh_doc_info.number()=0;
  std::ifstream input_file;
  std::ofstream output_file;
  char filename[100];

  // Get the partition to be used from file
  unsigned n_partition=problem.mesh_pt()->nelement();
  Vector<unsigned> element_partition(n_partition);
  sprintf(filename,"circular_cavity_1_partition.dat");
  input_file.open(filename);
  std::string input_string;
  for (unsigned e=0;e<n_partition;e++)
   {
    getline(input_file,input_string,'\n');
    element_partition[e]=atoi(input_string.c_str());
   }

  // Distribute the problem and check halo schemes
  problem.distribute(element_partition,mesh_doc_info,report_stats);
  problem.check_halo_schemes(mesh_doc_info);

  oomph_info << "---- Now solve after distribute ----" << std::endl;
#endif // (OOMPH_HAS_MPI)

  // Solve
  problem.newton_solve(max_adapt);

  // Step number
  doc_info.number()=1;
   
  //Output solution
  problem.doc_solution(doc_info);

 } // end of problem 1
 


 // Solve problem 2 with Crouzeix Raviart elements
 //--------------------------------------------
 {
  oomph_info << " " << std::endl
             << "-----------------------" << std::endl
             << "-- Running problem 2 --" << std::endl
             << "-----------------------" << std::endl;

  // Set up zero-Gravity vector
  Global_Physical_Variables::Gravity[0] = 0.0;
  Global_Physical_Variables::Gravity[1] = 0.0;

  // Set up Gamma vector for simplified form
  NavierStokesEquations<2>::Gamma[0]=0;
  NavierStokesEquations<2>::Gamma[1]=0;

  // Build problem with body force function and simplified form,
  // using body force function
  QuarterCircleDrivenCavityProblem<RefineableQCrouzeixRaviartElement<2> >
   problem(&Global_Physical_Variables::body_force);
 
  // Initial solve
  problem.newton_solve();

#ifdef OOMPH_HAS_MPI
  std::ifstream input_file;
  std::ofstream output_file;
  char filename[100];

  // Get the partition to be used from file
  unsigned n_partition=problem.mesh_pt()->nelement();
  Vector<unsigned> element_partition(n_partition);
  sprintf(filename,"circular_cavity_2_partition.dat");
  input_file.open(filename);
  std::string input_string;
  for (unsigned e=0;e<n_partition;e++)
   {
    getline(input_file,input_string,'\n');
    element_partition[e]=atoi(input_string.c_str());
   }

  // Distribute problem and check halo schemes
  problem.distribute(element_partition,mesh_doc_info,report_stats);
  problem.check_halo_schemes(mesh_doc_info);

  oomph_info << "---- Now solve after distribute ----" << std::endl;
#endif // (OOMPH_HAS_MPI)

  // Solve the problem with automatic adaptation
  problem.newton_solve(max_adapt);
  
  // Step number
  doc_info.number()=2;

  // Output solution
  problem.doc_solution(doc_info);

 } // end of problem 2

#ifdef OOMPH_HAS_MPI
// Finalise MPI
 MPI_Helpers::finalize();
#endif

} // end_of_main



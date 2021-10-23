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
///Driver for a 3D navier stokes entry flow problem in quarter tube domain

//Generic routines
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;

//=start_of_namespace================================================
/// Namespace for physical parameters
//===================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;
} // end_of_namespace


//=start_of_problem_class=============================================
/// Entry flow problem in quarter tube domain
//====================================================================
template<class ELEMENT>
class EntryFlowProblem : public Problem
{

public:

 /// Constructor: Pass DocInfo object and target errors
 EntryFlowProblem(DocInfo& doc_info, const double& min_error_target,
                  const double& max_error_target);

 /// Destructor (empty)
 ~EntryFlowProblem() {}

 /// Doc the solution after solve
 void actions_after_newton_solve() 
  {
   // Doc solution after solving
   doc_solution();

   // Increment label for output files
   Doc_info.number()++;
  }

 /// Update the problem specs before solve 
 void actions_before_newton_solve();

 /// After adaptation: Pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Pin redudant pressure dofs
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
  } 

 /// Doc the solution
 void doc_solution();

 /// Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 RefineableQuarterTubeMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableQuarterTubeMesh<ELEMENT>*>(Problem::mesh_pt());
  }

private:

 /// Exponent for bluntness of velocity profile
 int Alpha;
 
 /// Doc info object
 DocInfo Doc_info;

}; // end_of_problem_class




//=start_of_constructor===================================================
/// Constructor: Pass DocInfo object and error targets
//========================================================================
template<class ELEMENT>
EntryFlowProblem<ELEMENT>::EntryFlowProblem(DocInfo& doc_info,
                                            const double& min_error_target,
                                            const double& max_error_target) 
 : Doc_info(doc_info)
{ 

 // Setup mesh:
 //------------

 // Create geometric objects: Elliptical tube with half axes = radius = 1.0
 double radius=1.0;
 GeomObject* Wall_pt=new EllipticalTube(radius,radius);

 // Boundaries on object
 Vector<double> xi_lo(2);
 // height of inflow
 xi_lo[0]=0.0;
 // start of Wall_pt
 xi_lo[1]=0.0;

 Vector<double> xi_hi(2);
 // height of outflow
 xi_hi[0]=7.0;
 // end of Wall_pt
 xi_hi[1]=2.0*atan(1.0);

 // # of layers
 unsigned nlayer=6;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;

 // Build and assign mesh
 Problem::mesh_pt()=
  new RefineableQuarterTubeMesh<ELEMENT>(Wall_pt,xi_lo,frac_mid,xi_hi,nlayer);
 

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error()=max_error_target; 
 mesh_pt()->min_permitted_error()=min_error_target; 



 //Doc the boundaries
 ofstream some_file;
 char filename[100];
 sprintf(filename,"boundaries.dat");
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();
 
 
 // Set the boundary conditions for this problem: All nodal values are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Boundary 1 is the vertical symmetry boundary: We allow flow in 
     //the y-direction. Elsewhere, pin the y velocity
     if(ibound!=1) mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);

     // Boundary 2 is the horizontal symmetry boundary: We allow flow in 
     //the x-direction. Elsewhere, pin the x velocity
     if(ibound!=2) mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);

     // Boundaries 0 and 3 are the inflow and the wall respectively.
     // Pin the axial velocity because of the prescribed inflow
     // profile and the no slip on the stationary wall, respectively
     if((ibound==0) || (ibound==3)) 
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(2);
      }
    }
  } // end loop over boundaries

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  }

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());

 // Set Poiseuille flow as initial for solution
 // Inflow will be overwritten in actions_before_newton_solve()
 unsigned n_nod=mesh_pt()->nnode();
 for (unsigned j=0;j<n_nod;j++)
  {
   Node* node_pt=mesh_pt()->node_pt(j);
   // Recover coordinates
   double x=node_pt->x(0);
   double y=node_pt->x(1);
   double r=sqrt(x*x+y*y );  
   
   // Poiseuille flow
   node_pt->set_value(0,0.0);
   node_pt->set_value(1,0.0);
   node_pt->set_value(2,(1.0-r*r));
  }

 // Set the exponent for bluntness: Alpha=2 --> Poisseuille; anything
 // larger makes the inflow blunter
 Alpha=20;

 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=start_of_actions_before_newton_solve==========================================
/// Set the inflow boundary conditions
//========================================================================
template<class ELEMENT>
void EntryFlowProblem<ELEMENT>::actions_before_newton_solve()
{

 // (Re-)assign velocity profile at inflow values
 //--------------------------------------------

 // Setup bluntish parallel inflow on boundary 0:
 unsigned ibound=0; 
 unsigned num_nod= mesh_pt()->nboundary_node(ibound); 
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Recover coordinates
   double x=mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
   double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
   double r=sqrt(x*x+y*y);  
   
   // Bluntish profile for axial velocity (component 2)
   mesh_pt()->boundary_node_pt(ibound,inod)->
    set_value(2,(1.0-pow(r,Alpha)));
  }

} // end_of_actions_before_newton_solve


//=start_of_doc_solution==================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void EntryFlowProblem<ELEMENT>::doc_solution()
{ 
 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution

 


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=start_of_main=======================================================
/// Driver for 3D entry flow into a quarter tube. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only a single adaptation
//=====================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Allow (up to) five rounds of fully automatic adapation in response to 
 //-----------------------------------------------------------------------
 // error estimate
 //---------------
 unsigned max_adapt;
 double max_error_target,min_error_target;

 // Set max number of adaptations in black-box Newton solver and
 // error targets for adaptation
 if (CommandLineArgs::Argc==1)
  {
   // Up to five adaptations
   max_adapt=5;

   // Error targets for adaptive refinement
   max_error_target=0.005;
   min_error_target=0.0005;
  } 
 // Validation run: Only one adaptation. Relax error targets
 // to ensure that not all elements are refined so we're getting
 // some hanging nodes.
 else
  {
   // Validation run: Just one round of adaptation
   max_adapt=1;
   
   // Error targets for adaptive refinement
   max_error_target=0.02;
   min_error_target=0.002;
  }
 // end max_adapt setup
 

 // Set up doc info
 DocInfo doc_info;
 
 // Do Taylor-Hood elements
 //------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_TH");
  
  // Step number
  doc_info.number()=0;
  
  // Build problem
  EntryFlowProblem<RefineableQTaylorHoodElement<3> > 
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Taylor-Hood elements " << std::endl;


  // Doc solution after solving
  problem.doc_solution();
  
  // Solve the problem 
  problem.newton_solve(max_adapt);
 }


 // Do Crouzeix-Raviart elements
 //-----------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_CR");
  
  // Step number
  doc_info.number()=0;
  
  // Build problem
  EntryFlowProblem<RefineableQCrouzeixRaviartElement<3> >
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Crouzeix-Raviart elements " << std::endl;
  
  // Solve the problem 
  problem.newton_solve(max_adapt);
 }

} // end_of_main



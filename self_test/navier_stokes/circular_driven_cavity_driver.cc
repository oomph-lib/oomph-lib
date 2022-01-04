//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
//Driver for 2D adaptive driven cavity

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

// Mesh
#include "meshes/quarter_circle_sector_mesh.h"



using namespace std;

using namespace oomph;

//==================================================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100;

 /// Womersley number
 double ReSt=100.0;

}
//==================================================


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//====================================================================
/// Driven cavity problem in deformable ring domain.
///
/// Tecplot documentation: 
/// 
/// - 2DNavierStokesSolution.mcr documents the solution (no surprise!)
/// .
/// If refinement process was doced:
/// - 2DNavierStokesZ2Errors.mcr documents the error reduction during 
///   mesh refinement
/// - quadtree mesh refinement can be traced with various
///   other generic tecplot macros. 
/// 
//====================================================================
template<class ELEMENT, class MESH>
class DrivenCavityProblem : public Problem
{

public:

 /// Constructor
 DrivenCavityProblem();

 /// Destructor to clean up memory
 ~DrivenCavityProblem();


 //-------------------------------------------------------------------------
 /// Get pointer to wall
 //-------------------------------------------------------------------------
 GeomObject* wall_pt()
  {
   return Wall_pt;
  }
 //-------------------------------------------------------------------------


 //-------------------------------------------------------------------------
 /// Fix pressure in element ielem at pressure dof jpdof and set to pvalue
 //-------------------------------------------------------------------------
 void fix_pressure(const unsigned &ielem, const unsigned &jpdof, 
                   const double &pvalue)
  {
   //Fix the pressure at that element
   dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(ielem))->
                          fix_pressure(jpdof,pvalue);
  }
 //-------------------------------------------------------------------------

 /// Update the after solve (empty)
 void actions_after_newton_solve()
  {}

 //-------------------------------------------------------------------------
 /// Update the problem specs before solve. 
 /// Re-set velocity boundary conditions. 
 //-------------------------------------------------------------------------
 void actions_before_newton_solve()
 {

  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
   }
  
  // Overwrite with no flow along the other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned long ibound=1;ibound<num_bound;ibound++)
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
 
 }
 //-------------------------------------------------------------------------


 /// Finish problem setup: Setup element-specific things 
 /// (source fct pointers etc.)
 void actions_after_adapt();



 //-------------------------------------------------------------------------
 //Access function for the mesh
 //-------------------------------------------------------------------------
 MESH* mesh_pt() 
  {
   return dynamic_cast<MESH*>(Problem::mesh_pt());
  }
 //-------------------------------------------------------------------------

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,ofstream& trace_file);

 /// run parameter study
 void run(unsigned& next_istep);

 
private:

 /// Pointer to wall
 GeomObject* Wall_pt;

};
//========================================================================


//========================================================================
/// Constructor for DrivenCavity problem on deformable ring mesh
///
//========================================================================
template<class ELEMENT, class MESH>
DrivenCavityProblem<ELEMENT,MESH>::DrivenCavityProblem()
{ 

 // Build a linear solver: Use HSL's MA42 frontal solver
 //linear_solver_pt() = new HSL_MA42;
 
 // Switch off full doc for frontal solver
 //static_cast<HSL_MA42*>(linear_solver_pt())->disable_doc_stats();

 // Half axes for ellipse
 double a_ellipse=1.0;
 double b_ellipse=1.0;

 // Setup elliptical ring 
 Wall_pt=new Ellipse(a_ellipse,b_ellipse); 

 // End points for wall
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 //Now create the mesh
 double fract_mid=0.4;
 Problem::mesh_pt() = new MESH(Wall_pt,xi_lo,fract_mid,xi_hi);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned long ibound=0;ibound<num_bound;ibound++)
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
  }

 // Complete the build of all elements so they are fully functional
 //Find number of elements in mesh
 unsigned Nelement = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<Nelement;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() =  &Global_Physical_Variables::ReSt;

   // Free all internal pressure dofs (to make sure we don't 
   // accidentally pin two after refinement)
   unsigned num_internals=el_pt->ninternal_data();
   for (unsigned iint=0;iint<num_internals;iint++)
    {
     unsigned nvals=el_pt->internal_data_pt(iint)->nvalue();
     for (unsigned ival=0;ival<nvals;ival++)
      {
       el_pt->internal_data_pt(iint)->unpin(ival);
      }
    }
  }

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the pressure in final element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);
 
 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


}
//========================================================================


//========================================================================
/// Destructor for DrivenCavity problem on deformable ring mesh
///
//========================================================================
template<class ELEMENT, class MESH>
DrivenCavityProblem<ELEMENT,MESH>::~DrivenCavityProblem()
{ 

 // Timestepper gets killed in general problem destructor

 // Mesh gets killed in general problem destructor

}
//========================================================================




//========================================================================
/// Finish build of DrivenCavity problem on deformable ring mesh
///
/// Loop over elements and 
/// - setup pointers to physical parameters (Re, St, etc.)
/// - free internal pressure variables
/// 
/// Then choose one element and fix an internal pressure degree.
///
//========================================================================
template<class ELEMENT, class MESH>
void DrivenCavityProblem<ELEMENT,MESH>::actions_after_adapt()
{ 

 //Find number of elements in mesh
 unsigned Nelement = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<Nelement;i++)
  {

   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   // Free all internal pressure dofs (to make sure we don't 
   // accidentally pin two after refinement)
   unsigned num_internals=el_pt->ninternal_data();
   for (unsigned iint=0;iint<num_internals;iint++)
    {
     unsigned nvals=el_pt->internal_data_pt(iint)->nvalue();
     for (unsigned ival=0;ival<nvals;ival++)
      {
       el_pt->internal_data_pt(iint)->unpin(ival);
      }
    }
  }

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());

 // Now set the pressure in final element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);

}
//========================================================================





//========================================================================
/// Doc the solution
///
/// Process with
/// - 2DNavierStokesSolution.mcr
///
//========================================================================
template<class ELEMENT, class MESH>
void DrivenCavityProblem<ELEMENT,MESH>::doc_solution(DocInfo& doc_info,
                                                ofstream& trace_file)
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


 // Find first element in unrefined mesh
 ELEMENT* el_pt=dynamic_cast<ELEMENT*>(dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->
  quadtree_pt()->root_pt()->object_pt());

 // Number of nodes
 unsigned Nnode=el_pt->nnode();

 // Write trace file
 trace_file << mesh_pt()->nelement() << " " 
            << el_pt->nodal_position(Nnode-1,0) << " " 
	    << el_pt->nodal_position(Nnode-1,1) << " " 
            << el_pt->u_nst(Nnode-1,0) << " " 
            << el_pt->u_nst(Nnode-1,1) << " " 
            << std::endl;

 }
//========================================================================
 



//========================================================================
/// Run parameter study. Pass number of next step for consecutive
/// numbering of result files.
///
//========================================================================
template<class ELEMENT, class MESH>
void DrivenCavityProblem<ELEMENT,MESH>::run(unsigned& next_istep)
{ 
 // Label for output
 //-----------------
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=next_istep;
 
 
 // Doc refinement targets
 mesh_pt()->doc_adaptivity_targets(cout);

 // Initial refinement/mesh deformation
 //------------------------------------
 // Refine once uniformly
 refine_uniformly();
 
 // Refine once uniformly with doc
 refine_uniformly(doc_info);
   
 // Deform wall
 static_cast<Ellipse*>(wall_pt())->set_A_ellips(1.5);
 static_cast<Ellipse*>(wall_pt())->set_B_ellips(1.0);
   
 // Update mesh
 mesh_pt()->node_update();
 
 // Check mesh update
 cout << "Self test of algebraic-node-based mesh update: " 
      << mesh_pt()->self_test() << std::endl;
 
 
 // Setup parameters for problem adaptation
 //----------------------------------------
 
 // Don't allow refinement to drop under given level
 mesh_pt()->min_refinement_level()=2;
 
 // Don't allow refinement beyond given level 
 mesh_pt()->max_refinement_level()=5;
 
 // Get max/min refinement levels in mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 mesh_pt()->get_refinement_levels(min_refinement_level,
                                  max_refinement_level);
 
 cout << "\n Initial mesh: min/max refinement levels: " 
      << min_refinement_level << " " << max_refinement_level << std::endl << std::endl;
 
 

 
 // Open trace file
 //----------------
 ofstream trace_file;
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace%i.dat",doc_info.directory().c_str(),next_istep);
 trace_file.open(filename);
 trace_file << "# Driven cavity validation " << std::endl;
 trace_file << "# # of nodes along element edge " 
            << dynamic_cast<ELEMENT*>(mesh_pt()->finite_element_pt(0))
                ->nnode_1d() << std::endl;
 trace_file << "# err_max " << mesh_pt()->max_permitted_error() << std::endl;
 trace_file << "# err_min " << mesh_pt()->min_permitted_error() << std::endl;
 trace_file << "# min_refinement_level " 
            << mesh_pt()->min_refinement_level() << std::endl;
 trace_file << "# max_refinement_level " 
            << mesh_pt()->max_refinement_level() << std::endl;
 trace_file << "# A_ellips " <<
  static_cast<Ellipse*>(wall_pt())->a_ellips()<< std::endl; 
 trace_file << "# B_ellips " <<
  static_cast<Ellipse*>(wall_pt())->b_ellips()<< std::endl;
 trace_file << "# initial # elements " 
            << mesh_pt()->nelement() << std::endl;
 
   
   
   
 // Solve the problem on the initial mesh
 //---------------------------------------
 newton_solve();
 

 //Output solution
 doc_solution(doc_info,trace_file);
 
 
 // Loop over two different wall shapes
 //------------------------------------
 for (unsigned istep=0;istep<2;istep++)
  {
   
   // Adapt a few times
   unsigned n_adapt=2;
   
   for (unsigned iadapt=0;iadapt<n_adapt;iadapt++)
    {
     
     unsigned n_refined;
     unsigned n_unrefined;
     
     // Adapt problem
     adapt(n_refined,n_unrefined);
     
     
     // Doc max/min refinement levels
     mesh_pt()->get_refinement_levels(min_refinement_level,
                                      max_refinement_level);
     cout << "\n Adapted mesh: min/max refinement levels: " 
          << min_refinement_level << " " 
          << max_refinement_level << std::endl << std::endl;
     
     
     // Check convergence of adaptation cycle
     if ((n_refined==0)&&(n_unrefined==0))
      {
       cout << "\n\n----------------------------------------\n";
       cout <<"Solution is fully converged.\n";
       cout << "----------------------------------------\n \n";
       break;
      }    
     
     //Increment counter for solutions 
     next_istep++;
     doc_info.number()++;
     
     // Solve the problem again
     newton_solve();

     //Output solution
     doc_solution(doc_info,trace_file);
     
    }
   
   
   
   // Change wall and update mesh
   //----------------------------
   
   cout << "\n-----------------------" << std::endl;
   cout << "Changing wall shape " << std::endl;
   cout << "-----------------------" << std::endl;
   static_cast<Ellipse*>(wall_pt())->set_A_ellips(0.9); 
   static_cast<Ellipse*>(wall_pt())->set_B_ellips(1.1); 
   
   // Update mesh in response to change in wall 
   mesh_pt()->node_update();
   
   // Doc in trace file
   trace_file << "# A_ellips " <<
    static_cast<Ellipse*>(wall_pt())->a_ellips()<< std::endl; 
   trace_file << "# B_ellips " <<
    static_cast<Ellipse*>(wall_pt())->b_ellips()<< std::endl;
   
   
   // Check mesh update
   cout << "Max. error in algebraic-node-based mesh update: " 
        << mesh_pt()->self_test() << std::endl;
   
   
  } // end of loop over different orientations
 
 next_istep++;

}
//========================================================================
 





/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////









//=====================================================================
/// Driver for DrivenCavity test problem with mesh adaptation
//=====================================================================
int main()
{


 typedef AlgebraicElement<RefineableQCrouzeixRaviartElement<2> > ELEMENT0;
 typedef AlgebraicElement<RefineableQTaylorHoodElement<2> > ELEMENT1;
 typedef RefineableQCrouzeixRaviartElement<2> ELEMENT2;
 typedef RefineableQTaylorHoodElement<2> ELEMENT3;
 
 typedef AlgebraicRefineableQuarterCircleSectorMesh<ELEMENT0> MESH0;
 typedef AlgebraicRefineableQuarterCircleSectorMesh<ELEMENT1> MESH1;
 typedef RefineableQuarterCircleSectorMesh<ELEMENT2> MESH2;
 typedef RefineableQuarterCircleSectorMesh<ELEMENT3> MESH3;
 
// Counter for steps
 unsigned next_istep=0; 
  
 {
  
  cout << " " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << "Algebraic Taylor Hood " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << " " << std::endl;
  DrivenCavityProblem<ELEMENT1,MESH1>* problem_pt= new
   DrivenCavityProblem<ELEMENT1,MESH1>;
  
  problem_pt->run(next_istep);
 }
 
 
 {
  cout << " " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << "Algebraic Crouzeix Raviart " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << " " << std::endl;
  DrivenCavityProblem<ELEMENT0,MESH0>* problem_pt= new
   DrivenCavityProblem<ELEMENT0,MESH0>;
  
  problem_pt->run(next_istep);
 }
 
 {
  cout << " " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << "Taylor Hood " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << " " << std::endl;
  DrivenCavityProblem<ELEMENT3,MESH3>* problem_pt= new
   DrivenCavityProblem<ELEMENT3,MESH3>;
  
  problem_pt->run(next_istep);
 }
 
 {
  cout << " " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << "Crouzeix Raviart " << std::endl;
  cout << "==========================================================" << std::endl;
  cout << " " << std::endl;
  DrivenCavityProblem<ELEMENT2,MESH2>* problem_pt= new
   DrivenCavityProblem<ELEMENT2,MESH2>;
  
  problem_pt->run(next_istep);
 }

 
}
//=====================================================================









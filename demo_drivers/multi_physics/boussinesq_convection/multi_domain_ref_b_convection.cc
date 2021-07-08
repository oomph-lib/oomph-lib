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
//Driver for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion,
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"
#include "multi_physics.h"


// Both meshes are the standard rectangular quadmesh
#include "meshes/rectangular_quadmesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;



//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{

 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// \short Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh = 1800.0;
 
 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//=======start_of_problem_class=======================================
/// 2D Convection problem on two rectangular domains, discretised 
/// with refineable Navier-Stokes and Advection-Diffusion elements. 
/// The specific type of element is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
class RefineableConvectionProblem : public Problem
{

public:

 ///Constructor
 RefineableConvectionProblem();

 /// Destructor. Empty
 ~RefineableConvectionProblem() {}

 /// \short Update the problem specs before solve:
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Access function to the NST mesh. 
 /// Casts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<NST_ELEMENT>*>
    (Nst_mesh_pt);
  } // end_of_nst_mesh

 /// \short Access function to the AD mesh. 
 /// Casts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<AD_ELEMENT>* adv_diff_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<AD_ELEMENT>*>
    (Adv_diff_mesh_pt);
  } // end_of_ad_mesh

 /// Actions before adapt:(empty)
 void actions_before_adapt() {}

 /// \short Actions after adaptation, reset all sources, then
 /// re-pin a single pressure degree of freedom
 void actions_after_adapt()
  {
   //Unpin all the pressures in NST mesh to avoid pinning two pressures
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(nst_mesh_pt()->element_pt());

   //Pin the zero-th pressure dof in the zero-th element and set
   // its value to zero
   fix_pressure(0,0,0.0);

   // Set external elements for the multi-domain solution.
   Multi_domain_functions::
    setup_multi_domain_interactions<NST_ELEMENT,AD_ELEMENT>
    (this,nst_mesh_pt(),adv_diff_mesh_pt());
 
  } //end_of_actions_after_adapt


 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure in NST element
   if (nst_mesh_pt()->nelement()>0)
    {
     dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e))->
      fix_pressure(pdof,pvalue);
    }
  } // end_of_fix_pressure

 
 /// \short Set the
 /// boundary condition on the upper wall to be perturbed slightly
 /// to force the solution into the symmetry broken state.
 void enable_imperfection() {Imperfect = true;}

 /// \short Set th
 /// boundary condition on the upper wall to be unperturbed.
 void disable_imperfection() {Imperfect = false;}
 
 /// \short Doc the solution.
 void doc_solution();

private:
 
 /// DocInfo object
 DocInfo Doc_info;
 
 /// Is the boundary condition imperfect or not
 bool Imperfect;

protected:

 /// Navier Stokes mesh
 RefineableRectangularQuadMesh<NST_ELEMENT>* Nst_mesh_pt;

 /// Advection diffusion mesh
 RefineableRectangularQuadMesh<AD_ELEMENT>* Adv_diff_mesh_pt;

}; // end of problem class


//=======start_of_constructor=============================================
/// Constructor for adaptive thermal convection problem
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
RefineableConvectionProblem<NST_ELEMENT,AD_ELEMENT>::
RefineableConvectionProblem()
{ 
 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=9;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=3.0;

 // Domain length in y-direction
 double l_y=1.0;
 
 // Build the meshes
 Nst_mesh_pt =
  new RefineableRectangularQuadMesh<NST_ELEMENT>(n_x,n_y,l_x,l_y);
 Adv_diff_mesh_pt =
  new RefineableRectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 Nst_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Adv_diff_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set error targets for adaptive refinement
 Nst_mesh_pt->max_permitted_error()=0.5e-3; 
 Nst_mesh_pt->min_permitted_error()=0.5e-4; 
 Adv_diff_mesh_pt->max_permitted_error()=0.5e-3; 
 Adv_diff_mesh_pt->min_permitted_error()=0.5e-4; 


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries of the NST mesh
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max;

   //Loop over the number of nodes on the boundry
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {

     //If we are on the side-walls, the v-velocity and temperature
     //satisfy natural boundary conditions, so we only pin the
     //first value
     if((ibound==1) || (ibound==3)) 
      {
       val_max=1;
      }
     else
      {
       val_max=nst_mesh_pt()->boundary_node_pt(ibound,inod)->nvalue();
      }

     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }
 
 // Pin the zero-th pressure value in the zero-th element and
 // set its value to zero.
 fix_pressure(0,0,0.0);

 //Loop over the boundaries of the AD mesh
 num_bound = adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max;

   //Loop over the number of nodes on the boundry
   unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity and temperature
     //satisfy natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {
       val_max=0;
      }
     else // pin all values
      {
       val_max=adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->nvalue();
       //Loop over the desired values stored at the nodes and pin
       for(unsigned j=0;j<val_max;j++)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
        }
      }
    }
  } // end of loop over AD mesh boundaries
 
 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_nst_element = nst_mesh_pt()->nelement();
 for(unsigned i=0;i<n_nst_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   NST_ELEMENT *el_pt = dynamic_cast<NST_ELEMENT*>
    (nst_mesh_pt()->element_pt(i));

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   // We can ignore the external geometric data in the "external"
   // advection diffusion element when computing the Jacobian matrix
   // because the interaction does not involve spatial gradients of 
   // the temperature (and also because the mesh isn't moving!)
   el_pt->ignore_external_geometric_data();
  }

 unsigned n_ad_element = adv_diff_mesh_pt()->nelement();
 for(unsigned i=0;i<n_ad_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (adv_diff_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   // We can ignore the external geometric data in the "external"
   // Navier Stokes element when computing the Jacobian matrix
   // because the interaction does not involve spatial gradients of 
   // the velocities (and also because the mesh isn't moving!)
   el_pt->ignore_external_geometric_data();

  } // end of setup for all AD elements

 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Adv_diff_mesh_pt);
 build_global_mesh();
 
 // Set external elements for the multi-domain solution.
 Multi_domain_functions::
  setup_multi_domain_interactions<NST_ELEMENT,AD_ELEMENT>
  (this,nst_mesh_pt(),adv_diff_mesh_pt());
 
 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor


//===================start_actions_before_newton_solve====================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to include an imperfection (or not) depending on the control flag.
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableConvectionProblem<NST_ELEMENT,AD_ELEMENT>::actions_before_newton_solve()
{
 // Loop over the boundaries on the NST mesh
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=nst_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=nst_mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;
     //If we are on the side walls we only pin the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}
     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
      {
       //Add small velocity imperfection if desired
       if(Imperfect)
        {
         //Read out the x position
         double x = nod_pt->x(0);
         //Set a sinusoidal perturbation in the vertical velocity
         //This perturbation is mass conserving
         double value = sin(2.0*3.141592654*x/3.0);
         nod_pt->set_value(1,value);
        }
      }

    }
  }

 // Loop over all the boundaries on the AD mesh
 num_bound=adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=adv_diff_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=adv_diff_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the top boundary, set the temperature 
     //to -0.5 (cooled)
     if(ibound==2) {nod_pt->set_value(0,-0.5);}

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }


}  // end of actions before solve



//====================start_of_doc_solution===============================
/// Doc the solution
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableConvectionProblem<NST_ELEMENT,AD_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output Navier-Stokes solution
 sprintf(filename,"%s/fluid_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nst_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output advection diffusion solution
 sprintf(filename,"%s/temperature_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 adv_diff_mesh_pt()->output(some_file,npts);
 some_file.close();


 Doc_info.number()++;
} // end of doc


//===============start_of_main========================================
/// Driver code for 2D Boussinesq convection problem with 
/// adaptivity.
//====================================================================
int main()
{


 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 // Create the problem with 2D nine-node refineable elements.
 RefineableConvectionProblem<
 RefineableNavierStokesBoussinesqElement<
 RefineableQCrouzeixRaviartElement<2>,
  RefineableQAdvectionDiffusionElement<2,3> > ,
  RefineableAdvectionDiffusionBoussinesqElement<
  RefineableQAdvectionDiffusionElement<2,3>,
  RefineableQCrouzeixRaviartElement<2> >
  > problem;
 
 // Apply a perturbation to the upper boundary condition to
 // force the solution into the symmetry-broken state.
 problem.enable_imperfection();
 
 //Solve the problem with (up to) two levels of adaptation
 problem.newton_solve(2);
 
 //Document the solution
 problem.doc_solution();
 
 // Make the boundary conditions perfect and solve again. 
 // Now the slightly perturbed symmetry broken state computed
 // above is used as the initial condition and the Newton solver
 // converges to the symmetry broken solution, even without
 // the perturbation
 problem.disable_imperfection();
 problem.newton_solve(2);
 problem.doc_solution();

} // end of main










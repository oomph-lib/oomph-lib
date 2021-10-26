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

//Oomph-lib headers and derived elements are in a separate header file
#include "my_boussinesq_elements.h"
#include "multi_physics.h"

// Both meshes are the standard rectangular quadmesh
#include "meshes/simple_cubic_mesh.h"


// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

//================================================================
/// Function-type-object to perform comparison of elements
//================================================================
class ElementCmp
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x, GeneralisedElement* const &y) 
  const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);

   if((cast_x ==0) || (cast_y==0)) {return 0;}
   else
    {return cast_x->node_pt(0)->x(2) < cast_y->node_pt(0)->x(2);}
  }
};


//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{

 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=10.0;

 /// Reynolds number
 double Reynolds=10.0;

 /// Rayleigh number
 double Rayleigh = 100.0;
 
 /// Gravity vector
 Vector<double> Direction_of_gravity(3);
  
} // end_of_namespace


/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////



//=======start_of_problem_class=======================================
/// 3D Convection problem on two rectangular domains, discretised 
/// with refineable Navier-Stokes and Advection-Diffusion elements. 
/// The specific type of element is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
class RefineableConvectionProblem : public Problem
{

public:

 /// Constructor
 RefineableConvectionProblem();

 /// Destructor. clean up after the allocated memory
 ~RefineableConvectionProblem() 
  {
   //Fluid mesh
   //Delete the mesh's error estimator
   delete Nst_mesh_pt->spatial_error_estimator_pt();
   //Delete the mesh
   delete Nst_mesh_pt;
   
   //Temperature mesh
   //Delete the mesh's error estimator
   delete Adv_diff_mesh_pt->spatial_error_estimator_pt();
   //Delete the mesh
   delete Adv_diff_mesh_pt;
  }

 /// Update the problem specs before solve:
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Overloaded version of the problem's access function to 
 /// the NST mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableSimpleCubicMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RefineableSimpleCubicMesh<NST_ELEMENT>*>
    (Nst_mesh_pt);
  }

 /// Overloaded version of the problem's access function to 
 /// the AD mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableSimpleCubicMesh<AD_ELEMENT>* adv_diff_mesh_pt() 
  {
   return dynamic_cast<RefineableSimpleCubicMesh<AD_ELEMENT>*>
    (Adv_diff_mesh_pt);
  }

 /// Actions before adapt:(empty)
 void actions_before_adapt() {}

 /// Actions after adaptation, reset all sources, then
 /// re-pin a single pressure degree of freedom
 void actions_after_adapt()
  {
   // Setup all interactions
   Multi_domain_functions::setup_multi_domain_interactions
    <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),adv_diff_mesh_pt());
  }

 /// Doc the solution.
 void doc_solution();

 /// Switch to the iterative linear solve
 void switch_to_iterative_linear_solver();

private:
 
 /// DocInfo object
 DocInfo Doc_info;
 
 //Length of domain
 double Lz;

protected:

 RefineableSimpleCubicMesh<NST_ELEMENT>* Nst_mesh_pt;
 RefineableSimpleCubicMesh<AD_ELEMENT>* Adv_diff_mesh_pt;

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
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // # of elements in the z-direction
 unsigned n_z=12;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;
 
 // Domain length in the z-direction
 double l_z = 3.0;

 //Set the internal value
 Lz = l_z;
 
 // Build the meshes
 Nst_mesh_pt =
  new RefineableSimpleCubicMesh<NST_ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z);
 Adv_diff_mesh_pt =
  new RefineableSimpleCubicMesh<AD_ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z);

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

 //Loop over the boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum values of the velocity index
   unsigned vel_max=3;
   
   //If we are on the inlet or side walls pin everything
   //If we are on the outlet only pin transverse velocities
   if(ibound==5) {vel_max = 2;}

   //if on the side walls only pin the x-velocity
   if((ibound==2) || (ibound==4)) {vel_max=1;}

   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<vel_max;j++)
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }
 
 //Loop over the boundaries of the AD mesh
 num_bound = adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //If we are on the inlet or top and bottom walls pin everything
   if((ibound==0) || (ibound==1) || (ibound==3)) 
    {
     //Loop over the number of nodes on the boundry
     unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }
 
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
   el_pt->re_pt() = &Global_Physical_Variables::Reynolds;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Reynolds;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
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
  }

 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Adv_diff_mesh_pt);
 build_global_mesh();

 // Setup all interactions
 Multi_domain_functions::setup_multi_domain_interactions
  <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),adv_diff_mesh_pt());

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << endl; 
 
 
 //Set the linear solver
 /*linear_solver_pt() = new HSL_MA42;
 //Now sort the elements into z-order
 std::sort(mesh_pt()->element_pt().begin(),
 mesh_pt()->element_pt().end(),
 ElementCmp());*/

 switch_to_iterative_linear_solver();
} // end of constructor

/// Switch to an iterative linear solver
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableConvectionProblem<NST_ELEMENT,AD_ELEMENT>::
switch_to_iterative_linear_solver()
{
 // Set this to higher than default (10)
 Problem::Max_newton_iterations=20;
 //Build iterative linear solver
 GMRES<CRDoubleMatrix>* iterative_linear_solver_pt = new
  GMRES<CRDoubleMatrix>;

 // Set maximum number of iterations
 iterative_linear_solver_pt->max_iter() = 100;
 
 // Set tolerance
 iterative_linear_solver_pt->tolerance() = 1.0e-8;   

 FSIPreconditioner* prec_pt=new FSIPreconditioner(this); 
 
 // Tell preconditioner about meshes
 prec_pt->set_navier_stokes_mesh(Nst_mesh_pt);
 prec_pt->set_wall_mesh(Adv_diff_mesh_pt);

 // This is probably the most important because the temperature will
 // affect the fluid motion more than the other way round.
 prec_pt->use_block_triangular_version_with_solid_on_fluid();

 //Set the preconditioner
 iterative_linear_solver_pt->preconditioner_pt() = prec_pt;

 // Set linear solver
 linear_solver_pt() = iterative_linear_solver_pt;

#ifdef OOMPH_HAS_HYPRE
//If compiled with MPI, only use HYPRE if MPI has been initialised
#ifdef OOMPH_HAS_MPI
 if(MPI_Helpers::mpi_has_been_initialised())
#endif
  {
   // Need to create a SuperLUDistPreconditioner for the temperature ("solid")
   HyprePreconditioner* Temperature_prec_pt = new HyprePreconditioner;
   Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
    static_cast<HyprePreconditioner*>(Temperature_prec_pt)); 
   //SuperLUPreconditioner* Temperature_prec_pt =new SuperLUPreconditioner;
   prec_pt->set_solid_preconditioner_pt(Temperature_prec_pt);
   
   //Set up the internal preconditioners
   Preconditioner* P_matrix_preconditioner_pt = new HyprePreconditioner;
   
   // Set parameters for use as preconditioner on Poisson-type problem
   Hypre_default_settings::set_defaults_for_3D_poisson_problem(
    static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
   
   // Use Hypre for the Schur complement block
   prec_pt->navier_stokes_preconditioner_pt()->
    set_p_preconditioner(P_matrix_preconditioner_pt);
   
   Preconditioner* F_matrix_preconditioner_pt = new HyprePreconditioner;
   
   // Set parameters for use as preconditioner in for momentum 
   // block in Navier-Stokes problem
   Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
    static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt));
   
   // Use Hypre for momentum block 
   prec_pt->navier_stokes_preconditioner_pt()->
    set_f_preconditioner(F_matrix_preconditioner_pt);
  }
#endif
}

//===================start_actions_before_newton_solve====================
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to include an imperfection (or not) depending on the control flag.
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableConvectionProblem<NST_ELEMENT,AD_ELEMENT>::
actions_before_newton_solve()
{
 //Set the halfwith of the jet
 double x_halfwidth = 0.5;
 double y_halfwidth = 0.1875;
 //Set the peak flow of the jet
 double u_jet = 1.0;
 //Set the coordinates of the centre of the jet
 double x_jet = 0.5;  double y_jet = 0.375;
 //Bluntness of jet
 int alpha = 10;

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
     
     //If we are on the inlet
     if(ibound==0)
      {
       //No x-velocity
       nod_pt->set_value(0,0.0);
       //No y-velocity
       nod_pt->set_value(1,0.0);
       
       //Default z-velocity
       double zvel=0.0;
       
       //get x and y values
       double x = nod_pt->x(0);
       double y = nod_pt->x(1);

       //Jetlike z-velocity in square/rectangle
       if((std::abs(x-x_jet) <= x_halfwidth) &&
          (std::abs(y-y_jet) <= y_halfwidth))
        {
         zvel =   u_jet*(1.0 - pow(std::abs(x-x_jet)/x_halfwidth,alpha))*
          (1.0 - pow(std::abs(y-y_jet)/y_halfwidth,alpha));
        }
        
       
       //Set the velocity
       nod_pt->set_value(2,zvel);
      }
     
     //If we are on the upper and lower  walls
     if((ibound==1) || (ibound==3))
      {
       //No slip on all velocity components
       for(unsigned j=0;j<3;j++) {nod_pt->set_value(j,0.0);}
      }
     
     //If we are on the side walls
     if((ibound==2) || (ibound==4))
      {
       //No penetration on x-component only
       nod_pt->set_value(0,0.0);
      }

     //If we are on the outlet
     if(ibound==5)
      {
       //No slip on transverse  velocity components
       for(unsigned j=0;j<2;j++) {nod_pt->set_value(j,0.0);}
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

     //If we are on the inlet
     if(ibound==0)
      {
       //Set the temperature (hot)
       nod_pt->set_value(0,0.5);
      }

     //If we are on the upper and lower  walls, set the temperature
     if((ibound==1) || (ibound==3))
      {
       //The temperature has a discontinuity halfway along
       if(nod_pt->x(2) < 0.5*Lz)
        {
         nod_pt->set_value(0,0.5);
        }
       else
        {
         nod_pt->set_value(0,-0.5);
        }
      }
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

 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 sprintf(filename,"%s/fluid_%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nst_mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 sprintf(filename,"%s/temp_%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 adv_diff_mesh_pt()->output(some_file,npts);
 some_file.close();


 Doc_info.number()++;
} // end of doc


//===============start_of_main========================================
/// Driver code for 3D Boussinesq convection problem with 
/// adaptivity.
//====================================================================
int main(int argc, char **argv)
{


//Uncomment this to turn this into a parallel code
//#ifdef OOMPH_HAS_MPI
// // Set up MPI_Helpers
// MPI_Helpers::init(argc,argv);
// #endif

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;
 Global_Physical_Variables::Direction_of_gravity[2] = 0.0;

 // Create the problem with 3D twenty-seven-node refineable elements.
 RefineableConvectionProblem<
  RefineableQCrouzeixRaviartElementWithExternalElement<3>,
  RefineableQAdvectionDiffusionElementWithExternalElement<3> > problem;

 // Distribute the problem (including set sources)
// #ifdef OOMPH_HAS_MPI
//  DocInfo mesh_doc_info;
//  bool report_stats=true;
//  problem.distribute(mesh_doc_info,report_stats);
// #endif
  
 //Solve the problem with (up to) two levels of adaptation
 //problem.switch_to_iterative_linear_solver();
 problem.newton_solve();
 //Document the solution
 problem.doc_solution();


//#ifdef OOMPH_HAS_MPI 
// // finalize MPI
// MPI_Helpers::finalize();
//#endif

} // end of main










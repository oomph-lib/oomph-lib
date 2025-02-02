//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations and adds
//a temperature-dependent surface tension to give
//Bernard-Marangoni convection

//Oomph-lib headers, 
//We require the generic header
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"
//Include our coupling of advection-diffusion and Navier--Stokes
//the Boussinesq elements
#include "multi_physics.h"
//The fluid interface elements
#include "fluid_interface.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/single_layer_spine_mesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

namespace oomph
{

//==================================================================
/// Spine-based Marangoni surface tension elements that add
/// a linear dependence on temperature to the surface tension, 
/// which decreases with increasing temperature. In addition, 
/// this element adds a flux contribution to the advection-diffusion
/// equation to represent heat loss at the free surface. This
/// introduces the Biot number.
//=================================================================
template<class ELEMENT>
class SpineLineMarangoniFluidInterfaceElement :
public SpineLineFluidInterfaceElement<ELEMENT>
{
private:

 /// Pointer to a Biot number
 double *Bi_pt;

 /// Pointer to a Marangoni number
 double *Ma_pt;

 /// Index at which the temperature is stored at the nodes
 unsigned T_index;

 /// Default value of the physical constants
 static double Default_Physical_Constant_Value;

protected:

 /// The surface tension function is linear in the
 /// temperature with constant of proportionality equal
 /// to the Marangoni number.
 double sigma(const Vector<double> &s)
  {
   //Find the number of shape functions
   const unsigned n_node = this->nnode();
   //Now get the shape fuctions at the local coordinate
   Shape psi(n_node);
   this->shape(s,psi);

   //Now interpolate the temperature
   double T = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,T_index)*psi(l);
    }

   //Get the Marangoni and Capillary numbers
   double Ma = this->ma();
   double Ca = this->ca();
   //Return the variable surface tension
   //The additional multiplication by Ca will cancel with the 1/Ca
   //in the underlying equations
   return (1.0 - Ca*Ma*T);
  }

 /// Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
   {
    //Use finite differences to handle temperature variations
    const unsigned n_node = this->nnode();
    //Find the number of dofs in the element
    const unsigned n_dof = this->ndof();
    //Create newres vector
    Vector<double> newres(n_dof);
    
    //Integer storage for local unknown
    int local_unknown=0;
    
    //Use the default finite difference step
    const double fd_step = this->Default_fd_jacobian_step;
    
    //Loop over the nodes
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned t_index = this->T_index;

      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,t_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(t_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;
       
        //Calculate the new residuals
        this->get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
   }

   //Call the generic routine to handle the spine variables
   SpineElement<FaceGeometry<ELEMENT> >::
    fill_in_jacobian_from_geometric_data(jacobian);
  }

 
 /// Overload the Helper function to calculate the residuals and 
 /// jacobian entries. This particular function ensures that the
 /// additional entries are calculated inside the integration loop
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif, 
  const DShape &dpsifds,
  const DShape &dpsifdS,
  const DShape &dpsifdS_div,
  const Vector<double> &s,
  const Vector<double>& interpolated_x,
  const Vector<double> &interpolated_n, 
  const double &W,
  const double &J)
  {
   //Find the index at which the temperature is stored
   unsigned t_index = this->T_index;
   
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Read out the Bi number
   const double Bi = this->bi();

   //Now calculate the temperature at this point
   //Assuming the same shape functions are used (which they are)
   double T = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,t_index)*psif(l);
    }

   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;

   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,t_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       residuals[local_eqn] -= Bi*T*psif(l)*W*J;

       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,t_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= Bi*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes
  }  

 
 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Add the contribution to the jacobian
   this->fill_in_contribution_to_jacobian(residuals,jacobian);
   //No mass matrix terms, but should probably do kinematic bit here
  }

public:
 /// Constructor that passes the bulk element and face index down 
 /// to the underlying element
 SpineLineMarangoniFluidInterfaceElement(
 FiniteElement* const &element_pt, const int &face_index) : 
  SpineLineFluidInterfaceElement<ELEMENT>
  (element_pt,face_index)
  {
   //Initialise the values
   Bi_pt = &Default_Physical_Constant_Value;
   Ma_pt = &Default_Physical_Constant_Value;

   //Cast the bulk element 
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   //Now find the index at which the temperature is stored from the 
   //advection-diffusion part of the bulk element
   T_index = cast_element_pt->u_index_adv_diff();
  }

 /// Return the Biot number
 double bi() {return *Bi_pt;}
 
 /// Return the Marangoni number
 double ma() {return *Ma_pt;}

 /// Access function for pointer to the Marangoni number
 double* &ma_pt() {return Ma_pt;}

 /// Access function for pointer to the Biot number
 double* &bi_pt() {return Bi_pt;}

};


//Define the default physical value to be one
template<class ELEMENT>
double SpineLineMarangoniFluidInterfaceElement<ELEMENT>::
Default_Physical_Constant_Value = 1.0;

}

//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// Rayleigh number, set to be zero so that
 /// there are no gravitational effects
 double Rayleigh = 0.0;

 /// Scaled Bond number (Bo/Ca), set to be zero
 /// so that there are no gravitational effects
 double Scaled_Bond = 0.0;
 
 /// Biot number
 double Biot = 1.0;

 /// Marangoni number (just above the threshold for 
 /// linear instability)
 double Marangoni = 125.0;

 /// Capillary number (of which the results are independent
 /// for a pinned surface)
 double Capillary = 0.0045;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT,class INTERFACE_ELEMENT> 
class ConvectionProblem : public Problem
{

public:

 /// Constructor. The boolean indicates whether the free surface
 //should be pinned or not in the first instance
 ConvectionProblem(const bool &pin=true);

 /// Destructor. Empty
 ~ConvectionProblem() {}

 /// Release the free surface so that it can move
 void unpin_surface()
  {
   //Only bother if the surface is pinned
   if(Surface_pinned)
    {
     Surface_pinned = false;
     
     //Unpin the heights of all the spines in the middle
     unsigned n_spine = Bulk_mesh_pt->nspine();
     for(unsigned n=0;n<n_spine;n++)
      {
       Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->unpin(0);
      }
     
     //If we on the top wall, v velocity is no longer pinned
     unsigned ibound=2;
     //Loop over the number of nodes on the boundary
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //Loop over the desired values stored at the nodes and unpin
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->unpin(1);
      }
     
     //Unfix the pressure
     unfix_pressure(0,0);

     // Loop over the elements to set up element-specific 
     // and re-enable ALE
     unsigned n_element = Bulk_mesh_pt->nelement();
     for(unsigned i=0;i<n_element;i++)
      {
       // Upcast from GeneralsedElement to the present element
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

       el_pt->enable_ALE();
      }

     //Reassign the equation number
     std::cout << "Surface unpinned to give " 
               << assign_eqn_numbers() << " equation numbers\n";
    }
  }


 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Remember to update the nodes if the surface is not pinned
 void actions_before_newton_convergence_check()
  {
   if(!Surface_pinned) {Bulk_mesh_pt->node_update();}

   // This driver code cannot be allowed to use the analytical form of
   // get_dresidual_dnodal_coordinates(...) that is implemented in the
   // NavierStokesEquations class, since the elemental residuals have
   // contributions from external data which is not taken into account
   // by that routine. We therefore force the bulk elements to use the
   // fully-finite differenced version.
   const unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     ElementWithMovingNodes* el_pt =
      dynamic_cast<ElementWithMovingNodes*>(Bulk_mesh_pt->element_pt(e));
     el_pt->evaluate_shape_derivs_by_direct_fd();
    }
  }

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time_pt()->time());}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// UnFix pressure in element e at pressure dof pdof and set to pvalue
 void unfix_pressure(const unsigned &e, const unsigned &pdof)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    unfix_pressure(pdof);
  } // end_of_unfix_pressure


 /// Doc the solution.
 void doc_solution();

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 SingleLayerSpineMesh<ELEMENT>* Bulk_mesh_pt;
 
 Mesh* Surface_mesh_pt;

private:
 
 /// DocInfo object
 DocInfo Doc_info;

 /// Boolean to indicate whether the surface is pinned
 bool Surface_pinned;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::
ConvectionProblem(const bool &pin) : Surface_pinned(pin)
{
 //Allocate a timestepper
 add_time_stepper_pt(new BDF<2>);

 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=3.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build a standard rectangular quadmesh
 Bulk_mesh_pt = 
  new SingleLayerSpineMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 for(unsigned i=0;i<n_x;i++)
  {
   //Construct a new 1D line element on the face on which the local
   //coordinate 1 is fixed at its max. value (1) --- This is face 2
   FiniteElement *interface_element_pt =
    new INTERFACE_ELEMENT(
     Bulk_mesh_pt->finite_element_pt(n_x*(n_y-1)+i),2);
   
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 
  }
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 
 //Pin the heights of all the spines if the surface is pinned
 if(Surface_pinned)
  {
   unsigned n_spine = Bulk_mesh_pt->nspine();
   for(unsigned n=0;n<n_spine;n++)
    {
     Bulk_mesh_pt->spine_pt(n)->spine_height_pt()->pin(0);
    }
  }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the minimum index to be pinned (all values by default)
   unsigned val_min=0;
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=3;
   //If we are on the side-walls, the v-velocity and temperature
   //satisfy natural boundary conditions, so we only pin the
   //first value
   if((ibound==1) || (ibound==3)) {val_max=1;}

   //If we on the top wall, v velocity is pinned
   if(ibound==2) 
    {
     //If the surface is pinned, pin the v velocity
     if(Surface_pinned) {val_min=1; val_max=2;}
     //Otherwise pin nothing
     else {val_min=0; val_max=0;}
    }

   //Loop over the number of nodes on the boundary
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=val_min;j<val_max;j++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 if(Surface_pinned) {fix_pressure(0,0,0.0);}

 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Re/Fr equal to Bo/Ca, the scaled Bond number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::Scaled_Bond;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   //If the mesh is fixed, we can disable ALE
   if(Surface_pinned) {el_pt->disable_ALE();}
  }


  // Loop over the interface elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   // Upcast from GeneralsedElement to the present element
   INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(
    Surface_mesh_pt->element_pt(i));
   
   // Set the Biot number
   el_pt->bi_pt() = &Global_Physical_Variables::Biot;

   // Set the Marangoni number
   el_pt->ma_pt() =&Global_Physical_Variables::Marangoni;

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Capillary;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //Set initial temperature profile
 if(time <= 0.0)
  {
   unsigned n_node = Bulk_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = Bulk_mesh_pt->node_pt(n);
     //Set linear variation
     nod_pt->set_value(2,2.0 - nod_pt->x(1));
    }
  }

 // Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //If we are on the top boundary, do not set the velocities
     //(yet)
     if(ibound==2) {vel_max = 0;}

     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the bottom boundary peturb the velocity
     if(ibound==0) //2 
      {
       //Add small velocity imperfection if desired
       double epsilon = 0.01;

       //Read out the x position
       double x = nod_pt->x(0);

       //Set a sinusoidal perturbation in the vertical velocity
       //This perturbation is mass conserving
       double value = sin(2.0*MathematicalConstants::Pi*x/3.0)*
        epsilon*time*exp(-time);
       nod_pt->set_value(1,value);
      }

            
       //If we are on the bottom boundary, set the temperature
       //to 2 (heated)
     if(ibound==0)
      {
       nod_pt->set_value(2,2.0);
      }
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Construct our problem
 ConvectionProblem<SpineElement<BuoyantQCrouzeixRaviartElement<2> >,
  SpineLineMarangoniFluidInterfaceElement<BuoyantQCrouzeixRaviartElement<2> > > 
problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);
 
 //Perform a single steady Newton solve
 problem.steady_newton_solve();
 
 //Document the solution
 problem.doc_solution();

 //Now release the interface for real fun
 //problem.unpin_surface();

 //Set the timestep
 double dt = 0.1;

 //Initialise the value of the timestep and set initial values 
 //of previous time levels assuming an impulsive start.
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 200;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }

} // end of main










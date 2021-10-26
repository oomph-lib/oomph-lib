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
//
//Driver for a multi-physics problem that couples a Navier--Stokes
//mesh to two advection diffusion meshes, giving double diffusive convection
//in the Boussinesq approximation

//Oomph-lib headers and derived elements are in a separate header file
#include "db_nst_external_elements.h"

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
 /// The Lewis number
 double Lewis = 10.0;

 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// Thermal Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh_T = 1800.0;
 
 /// Solutal Rayleigh number
 double Rayleigh_S = -1000;

 /// Length of domain
 double Lambda = 1.414;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace


/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////



//=======start_of_problem_class=======================================
/// 2D Convection problem on two rectangular domains, discretised 
/// with refineable Navier-Stokes and Advection-Diffusion elements. 
/// The specific type of element is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
class RefineableDDConvectionProblem : public Problem
{

public:

 /// Constructor
 RefineableDDConvectionProblem();

 /// Destructor. Empty
 ~RefineableDDConvectionProblem() 
  {
   //Delete the meshes
   delete Conc_mesh_pt;
   delete Temp_mesh_pt;
   delete Nst_mesh_pt;
   //Delete the timestepper
   delete this->time_stepper_pt();
  }


 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
   set_boundary_conditions(time_pt()->time());
  }


 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);


 /// Access function to the NST mesh. 
 /// Casts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<NST_ELEMENT>*>
    (Nst_mesh_pt);
  } // end_of_nst_mesh

 /// Access function to the AD mesh. 
 /// Casts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<AD_ELEMENT>* temp_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<AD_ELEMENT>*>
    (Temp_mesh_pt);
  } // end_of_ad_mesh

 /// Access function to the AD mesh. 
 /// Casts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<AD_ELEMENT>* conc_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<AD_ELEMENT>*>
    (Conc_mesh_pt);
  } // end_of_ad_mesh


/// Get kinetic energy and kinetic energy flux
 void get_kinetic_energy(double &E, double &Edot)
  {
   //Reset values to zero
   E = 0.0; Edot=0.0;
   
   //Loop over the elements
   unsigned n_element = nst_mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     NST_ELEMENT* elem_pt = 
      dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e));
     
     E += elem_pt->kin_energy();
     Edot += elem_pt->d_kin_energy_dt();
    }
  }


 /// Actions before adapt:(empty)
 void actions_before_adapt() {}

 /// Actions after adaptation, reset all sources, then
 /// re-pin a single pressure degree of freedom
 void actions_after_adapt()
  {
   //Unpin all the pressures in NST mesh to avoid pinning two pressures
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(nst_mesh_pt()->element_pt());

   //Pin the zero-th pressure dof in the zero-th element and set
   // its value to zero
   fix_pressure(0,0,0.0);

   // Set sources for temperature
   Multi_domain_functions::setup_multi_domain_interactions
    <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),temp_mesh_pt());

   // Set sources for concentration
   Multi_domain_functions::setup_multi_domain_interactions
    <NST_ELEMENT,AD_ELEMENT>(this,nst_mesh_pt(),conc_mesh_pt(),1,0);

  } //end_of_actions_after_adapt


 /// Fix pressure in element e at pressure dof pdof and set to pvalue
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

  
 /// Doc the solution.
 void doc_solution();

private:
 
 /// DocInfo object
 DocInfo Doc_info;
 
protected:

 RefineableRectangularQuadMesh<NST_ELEMENT>* Nst_mesh_pt;
 RefineableRectangularQuadMesh<AD_ELEMENT>* Temp_mesh_pt;
 RefineableRectangularQuadMesh<AD_ELEMENT>* Conc_mesh_pt;

}; // end of problem class


//=======start_of_constructor=============================================
/// Constructor for adaptive thermal convection problem
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
RefineableDDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::
RefineableDDConvectionProblem()
{ 
 
 //Allocate a timestepper
 add_time_stepper_pt(new BDF<2>);

 // Set output directory
 Doc_info.set_directory("RESLT_ref_multimesh");
 
 // # of elements in x-direction
 unsigned n_x=9;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=Global_Physical_Variables::Lambda;

 // Domain length in y-direction
 double l_y=1.0;
 
 // Build the meshes
 Nst_mesh_pt =
  new RefineableRectangularQuadMesh<NST_ELEMENT>(n_x,n_y,l_x,l_y,
                                                 time_stepper_pt());
 Temp_mesh_pt =
  new RefineableRectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,
                                                time_stepper_pt());
 Conc_mesh_pt =
  new RefineableRectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,
                                                time_stepper_pt());

 // Create/set error estimator
 Nst_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Temp_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Conc_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set error targets for adaptive refinement
 Nst_mesh_pt->max_permitted_error()=0.5e-3; 
 Nst_mesh_pt->min_permitted_error()=0.5e-4; 
 Temp_mesh_pt->max_permitted_error()=0.5e-3; 
 Temp_mesh_pt->min_permitted_error()=0.5e-4; 
 Conc_mesh_pt->max_permitted_error()=0.5e-3; 
 Conc_mesh_pt->min_permitted_error()=0.5e-4; 


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here




 //Loop over the boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity 
     //satisfies natural boundary conditions, so we only pin the
     //first value
     if ((ibound==1) || (ibound==3)) 
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
     else // On the top and bottom walls, we have "stress-free" conditions
      //which actually corresponds to transverse stress free and normal 
      //zero velocity (symmetry)
      //Thus we pin the second value
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 //Loop over the boundaries of the AD mesh
 num_bound = temp_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= temp_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the temperature
     //satisfies natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {
      
      }
     //Otherwise pin the temperature
     else // pin all values
      {
       temp_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }


 //Loop over the boundaries of the AD mesh
 num_bound = conc_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Loop over the number of nodes on the boundry
   unsigned num_nod= conc_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the concentration
     //satisfies natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {

      }
     //Otherwiwse pin the concentration
     else // pin all values
      {
       conc_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
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
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_t_pt() = &Global_Physical_Variables::Rayleigh_T;

   // Set the Solutal Rayleigh number
   el_pt->ra_s_pt() = &Global_Physical_Variables::Rayleigh_S;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
  }


 unsigned n_temp_element = temp_mesh_pt()->nelement();
 for(unsigned i=0;i<n_temp_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (temp_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the timescale to be the same as the Navier--Stokes
   // equations (1.0)
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }

 unsigned n_conc_element = conc_mesh_pt()->nelement();
 for(unsigned i=0;i<n_conc_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (conc_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Lewis;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Lewis;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }


 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Temp_mesh_pt);
 add_sub_mesh(Conc_mesh_pt);
 build_global_mesh();

 // Setup the interaction
 actions_after_adapt();

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << endl; 

 // Set this to higher than default (10)
 Problem::Max_newton_iterations=20;

} // end of constructor


//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableDDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::
set_boundary_conditions(const double &time)
{
 // Loop over all the boundaries on the NST mesh
 unsigned num_bound=nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=nst_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=nst_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {nod_pt->set_value(0,0.0);}
     //If we are on the top and bottom walls we only set the y-velocity
     else {nod_pt->set_value(1,0.0);}
    }
  }

 // Loop over all the boundaries on the AD mesh
 num_bound=temp_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=temp_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=temp_mesh_pt()->boundary_node_pt(ibound,inod);
     
     //If we are on the top boundary, set the temperature 
     //to -0.5 (cooled)
     if(ibound==2) {nod_pt->set_value(0,-0.5);}

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }



 // Loop over all the boundaries on the AD mesh
 num_bound=conc_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=conc_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=conc_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the top boundary, set the concentration to be low
     //to -0.5 (cooled)
     if(ibound==2) 
      {
       nod_pt->set_value(0,-0.5);
       
       //Add small concentration imperfection if desired
       double epsilon = 0.01;
       
       //Read out the x position
       double x = nod_pt->x(0);
       
       //Set a sinusoidal perturbation in the concentration
       double value = sin(2.0*MathematicalConstants::Pi*x/1.5)*
        epsilon*time*exp(-10.0*time);
       nod_pt->set_value(0, -0.5 + value);
      }

     //If we are on the bottom boundary, set the concentration to be high
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }


} // end_of_set_boundary_conditions




//====================start_of_doc_solution===============================
/// Doc the solution
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void RefineableDDConvectionProblem<NST_ELEMENT,AD_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output whole solution (this will output elements from one mesh
 //----------------------  followed by the other mesh at the moment...?)
 /*sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();*/


 sprintf(filename,"%s/nst_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nst_mesh_pt()->output(some_file,npts);
 some_file.close();

 sprintf(filename,"%s/temp_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 temp_mesh_pt()->output(some_file,npts);
 some_file.close();

 sprintf(filename,"%s/conc_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 conc_mesh_pt()->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//===============start_of_main========================================
/// Driver code for 2D Boussinesq convection problem with 
/// adaptivity.
//====================================================================
int main(int argc, char **argv)
{


 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 // Create the problem with 2D nine-node refineable elements.
 RefineableDDConvectionProblem<
  RefineableQCrouzeixRaviartElementWithTwoExternalElement<2>,
  RefineableQAdvectionDiffusionElementWithExternalElement<2> > problem;
 


 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);

 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution();

  //Start a trace file
 ofstream trace("RESLT_ref_multimesh/trace.dat");
 //Local variables for the kinetic energy and its rate of change
 double E=0.0, Edot = 0.0;
 
 //Output to the trace file
 problem.get_kinetic_energy(E,Edot);
 trace << problem.time_pt()->time() << " " 
       << problem.nst_mesh_pt()->boundary_node_pt(1,8)->value(1) << " " 
       << E << " " << Edot << std::endl;
 
 //Set the timestep
 double dt = 0.01;

 //Initialise the value of the timestep and set an impulsive start
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 2000;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 bool first=true;

 unsigned n_refine = 0;

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt,n_refine,first);
   problem.doc_solution();

   first = false;

   //Refine the last case for the demo
   if((argc > 1) && (i==3)) {n_refine=1;}
   //Only start to refine after 700 steps once the
   //limit cycle has essentially settled down
   if(i > 700) {n_refine=1;}

 //Output to the trace file
   problem.get_kinetic_energy(E,Edot);
   trace << problem.time_pt()->time() << " " 
         << problem.nst_mesh_pt()->boundary_node_pt(1,8)->value(1) << " "
         << E << " " << Edot << std::endl;
  }

 trace.close();


} // end of main










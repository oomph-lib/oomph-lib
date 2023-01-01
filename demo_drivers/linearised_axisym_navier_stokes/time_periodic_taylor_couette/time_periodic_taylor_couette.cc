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
// Driver for validation problem for linearised 3D N-S equations in
// cylindrical polar coordinates which replicates the investigation of
// B.T. Murray, G.B. McFadden and S.R. Coriell in the 1990 paper
// "Stabilization of Taylor-Couette flow due to time-periodic outer
// cylinder oscillation".

// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Linearised axisymmetric Navier-Stokes headers
#include "linearised_axisym_navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

// Self-starting BDF2 timestepper
#include "../self_starting_BDF2_timestepper.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



//==start_of_namespace_for_physical_parameters===========================
/// Namespace for physical parameters
//=======================================================================
namespace GlobalPhysicalVariables
{

 /// Reynolds number
 double Re = 1.0;

 /// Product of Reynolds number and Strouhal number
 double ReSt = 1.0;

 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 0.0; // (Fr = inf)

 /// Murray, McFadden and Coriell's defn of the Reynolds number
 /// Note: Do not specify this as it is determined from the information
 /// below
 double MMC_Re_current;

 /// Loop information
 double MMC_Re_lower_limit = 20.93;
 double MMC_Re_upper_limit = 20.94;
 unsigned Nparam_steps = 5;

 /// Angular frequency of outer cylinder rotation (omega = 2*gamma^2)
 double AngularFrequency = 10.58; // corresponds to gamma=2.3

 /// Dimensionless modulation amplitude (epsilon)
 double Epsilon = 0.9;

 /// Radius ratio (radius of inner cylinder / radius of outer cylinder)
 double Eta = 0.88;

 /// Wavenumber of azimuthal waves e^iaz
 double a = 3.2;

 /// Azimuthal mode number k in e^ik(theta) decomposition
 int k = 0;

 /// Direction of gravity
 Vector<double> G(3);

} // End of namespace



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=====start_of_base_state_problem_class=================================
/// Base state problem class for the time-periodic Taylor--Couette
/// problem
//=======================================================================
template<class BASE_ELEMENT>
class BaseStateProblem : public Problem
{

public:

 /// Constructor
 BaseStateProblem(const unsigned& n_r,
                  const unsigned& n_z,
                  const double& radius_inner_cylinder,
                  const double& radius_outer_cylinder,
                  const double& l_z);

 /// Destructor (empty)
 ~BaseStateProblem() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep()
  {
   set_boundary_conditions(time_pt()->time());
  }

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// Set initial condition (incl previous timesteps)
 void set_initial_condition();

 /// Set the boundary conditions
 void set_boundary_conditions(const double& time);

 /// Access function for the specific mesh
 RectangularQuadMesh<BASE_ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RectangularQuadMesh<BASE_ELEMENT>*>(Problem::mesh_pt());
  }

 /// Access function for the specific timestepper
 SelfStartingBDF2* time_stepper_pt()
  {
   return dynamic_cast<SelfStartingBDF2*>(Problem::time_stepper_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,
                   const bool& output_soln);

 /// Create a trace file
 void create_trace_file(DocInfo& doc_info)
  {
   // Open trace file
   char filename1[256];
   sprintf(filename1,"%s/base_trace_epsilon%2.1f_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current);
   Trace_file.open(filename1);
  }

 /// Initialise trace file (print column headings)
 void initialise_trace_file()
  {
   Trace_file << "time, norm_of_dof_vector" << std::endl;
  }

 /// Clear and close trace file
 void close_trace_file()
  {
   Trace_file.clear();
   Trace_file.close();
  }

 /// Access function for trace file
 ofstream& trace_file() { return Trace_file; }

private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   dynamic_cast<BASE_ELEMENT*>(mesh_pt()->element_pt(e))
    ->fix_pressure(pdof,pvalue);
  }
 
 /// Trace file
 ofstream Trace_file;

}; // End of base_state_problem class



//========start_of_base_state_constructor================================
/// Constructor for base state time-periodic Taylor--Couette problem
//=======================================================================
template<class BASE_ELEMENT>
BaseStateProblem<BASE_ELEMENT>::
BaseStateProblem(const unsigned& n_r,
                 const unsigned& n_z,
                 const double& radius_inner_cylinder,
                 const double& radius_outer_cylinder,
                 const double& l_z)
{
 // Be less verbose during newton solve
 Problem::disable_info_in_newton_solve();

 // Be less verbose about linear solve timings
 linear_solver_pt()->disable_doc_time();

 // Overwrite maximum allowed residual to accomodate possibly
 // poor initial guess for solution
 Problem::Max_residuals=1000;

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new SelfStartingBDF2);

 // Build and assign mesh
 Problem::mesh_pt() =
  new RectangularQuadMesh<BASE_ELEMENT>(n_r,n_z,radius_inner_cylinder,
                                        radius_outer_cylinder,0.0,l_z,
                                        time_stepper_pt());
 
 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------
 
 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here
 
 // Determine number of boundaries in the mesh
 const unsigned n_boundary = mesh_pt()->nboundary();

 // Loop over boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);
   
   // Loop over nodes on boundary b
   for (unsigned n=0;n<n_node;n++)
    {
     // Pin axial velocity component on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(1);
     
     // Pin radial and azimuthal velocity components on all solid boundaries
     if(b==1 || b==3)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); // Radial
       mesh_pt()->boundary_node_pt(b,n)->pin(2); // Azimuthal
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries

 // Determine total number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();

 // Pin all radial and axial velocities throughout the bulk of the domain
 for(unsigned n=0;n<n_node;n++)
  {
//   mesh_pt()->node_pt(n)->pin(0); // Radial
//   mesh_pt()->node_pt(n)->pin(1); // Axial
  }

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of elements in the mesh
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   BASE_ELEMENT *el_pt = dynamic_cast<BASE_ELEMENT*>(mesh_pt()->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

   // The mesh remains fixed
   el_pt->disable_ALE();
   
  } // End of loop over elements

 // Set the pressure in first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);

 // Set up equation numbering scheme
 std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of base state constructor



//======start_of_base_state_set_initial_condition========================
/// Set the initial conditions to be zero everywhere for base state
//=======================================================================
template<class BASE_ELEMENT>
void BaseStateProblem<BASE_ELEMENT>::set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the three velocity components
   for(unsigned i=0;i<3;i++)
    {
     // Set velocity component i of node n to zero
     mesh_pt()->node_pt(n)->set_value(i,0.0);
    }
  }
 
 // Initialise the previous velocity values for timestepping
 // corresponding to an impulsive start
 assign_initial_values_impulsive();
 
} // End of set_initial_condition for base state



//======start_of_base_state_set_boundary_conditions======================
/// Reset the boundary conditions for the current time
//=======================================================================
template<class BASE_ELEMENT>
void BaseStateProblem<BASE_ELEMENT>::
set_boundary_conditions(const double& time)
{
 // Cache radius ratio
 const double eta = GlobalPhysicalVariables::Eta;

 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();
   
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);
   
   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       // For the inner solid boundary (boundary 3)
       if(b==3)
        {
         // Evaluate boundary condition
         double azi_vel=eta*GlobalPhysicalVariables::MMC_Re_current/(1.0-eta);
         
         // Set all velocity components to no flow along boundary
         switch(i)
          {
          case 2: // Azimuthal velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,azi_vel);
           break;
          case 1: // Axial velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
           break;
          case 0: // Radial velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
           break;
          }
        }
       
       // For the outer solid boundary (boundary 1)
       if(b==1)
        {
         // Evaluate time-dependent boundary condition
         double azi_vel
          = (GlobalPhysicalVariables::Epsilon
             *GlobalPhysicalVariables::MMC_Re_current/(1.0-eta))
          *cos(GlobalPhysicalVariables::AngularFrequency*time);
         
         // Set all velocity components to no flow along boundary
         switch(i)
          {
          case 2: // Azimuthal velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,azi_vel);
           break;
          case 1: // Axial velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
           break;
          case 0: // Radial velocity
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
           break;
          }
        }
       
       // For the symmetry boundaries (boundaries 0 and 2)
       if(b==0 || b==2)
        {
         // Set only the axial (i=1) velocity component to zero
         // (no penetration of symmetry boundary)
         if(i==1) { mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0); }
        }
      } // End of loop over velocity components
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
} // End of set_boundary_conditions for base state



//==start_of_base_state_doc_solution=====================================
/// Document the base state solution
//=======================================================================
template<class BASE_ELEMENT>
void BaseStateProblem<BASE_ELEMENT>::doc_solution(DocInfo &doc_info,
                                                  const bool& output_soln)
{
 // Create vector of dofs
 DoubleVector dofs;
 
 // Get dofs
 this->get_dofs(dofs);
 
 // Get L2 norm of dof vector
 const double dof_norm = dofs.norm();
 
 // Document in trace file
 Trace_file << time_pt()->time() << " " << dof_norm << std::endl;
 
 // If desired, output solution to file
 if(output_soln)
  {
   ofstream some_file;
   char filename[256];
   
   // Set number of plot points (in each coordinate direction)
   const unsigned npts = 2;//5;
   
   // Open solution output file
   sprintf(filename,"%s/base_soln_epsilon%2.1f_Re%4.2f_soln%i.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current,
           doc_info.number());
   some_file.open(filename);
   
   // Output solution to file
   mesh_pt()->output(some_file,npts);
   
   // Close solution output file
   some_file.close();
  }

} // End of doc_solution for base state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_perturbed_state_problem_class===============================
/// Perturbed state problem class for the time-periodic Taylor--Couette
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
class PerturbedStateProblem : public Problem
{
 
public:
 
 /// Constructor
 PerturbedStateProblem(const unsigned& n_r,
                       const unsigned& n_z,
                       const double& radius_inner_cylinder,
                       const double& radius_outer_cylinder,
                       const double& l_z,
                       Mesh* external_mesh_pt);
 
 /// Destructor (empty)
 ~PerturbedStateProblem() {}
 
 /// Update before solve (empty)
 void actions_before_newton_solve() {}
 
 /// Update after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Actions before timestep (empty)
 void actions_before_implicit_timestep() {}
 
 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// Set initial conditions to a "Poiseuille-style" profile
 void set_initial_condition();

 /// Set up the (homogeneous) boundary conditions
 void set_boundary_conditions();
 
 /// Access function for the specific mesh
 RectangularQuadMesh<PERTURBED_ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RectangularQuadMesh<PERTURBED_ELEMENT>*>
    (Problem::mesh_pt());
  }
 
 /// Access function for the base state mesh
 RectangularQuadMesh<BASE_ELEMENT>* base_state_mesh_pt()
  {
   return dynamic_cast<RectangularQuadMesh<BASE_ELEMENT>*>
    (Base_state_mesh_pt);
  }
 
 /// Access function for the specific timestepper
 SelfStartingBDF2* time_stepper_pt()
  {
   return dynamic_cast<SelfStartingBDF2*>(Problem::time_stepper_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,const bool& output_soln=true);

 /// Create a trace file
 void create_trace_file(DocInfo& doc_info)
  {
   // Open trace file
   char filename[256];
   sprintf(filename,"%s/perturbed_trace_epsilon%2.1f_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current);
   Trace_file.open(filename);
  }

 /// Initialise trace file (print column headings)
 void initialise_trace_file()
  {
   Trace_file << "time, norm_of_dof_vector" << std::endl;
  }

 /// Clear and close trace file
 void close_trace_file()
  {
   Trace_file.clear();
   Trace_file.close();
  }

 /// Access function for trace file
 ofstream& trace_file() { return Trace_file; }

private:

 /// Pointer to the base state mesh
 Mesh* Base_state_mesh_pt;

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned& e, const unsigned& pdof, 
                   const double& pvalue)
  {
   dynamic_cast<PERTURBED_ELEMENT*>(mesh_pt()->element_pt(e))
    ->fix_pressure(pdof,pvalue);
  }

 /// Trace file
 ofstream Trace_file;

 /// Inner and outer cylinder radii
 double Inner_cylinder_radius;
 double Outer_cylinder_radius;

 /// Height of domain
 double Domain_height;

}; // End of perturbed_state_problem class



//==start_of_perturbed_state_constructor=================================
/// Constructor for perturbed state time-periodic Taylor--Couette problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
PerturbedStateProblem(const unsigned& n_r,
                      const unsigned& n_z,
                      const double& radius_inner_cylinder,
                      const double& radius_outer_cylinder,
                      const double& l_z,
                      Mesh* external_mesh_pt)
 : Base_state_mesh_pt(external_mesh_pt),
   Inner_cylinder_radius(radius_inner_cylinder),
   Outer_cylinder_radius(radius_outer_cylinder),
   Domain_height(l_z)
{
 // Be less verbose during newton solve
 Problem::disable_info_in_newton_solve();

 // Be less verbose about linear solve timings
 linear_solver_pt()->disable_doc_time();
 
 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new SelfStartingBDF2);
 
 // Build and assign mesh
 Problem::mesh_pt() = 
  new RectangularQuadMesh<PERTURBED_ELEMENT>(n_r,n_z,radius_inner_cylinder,
                                             radius_outer_cylinder,0.0,l_z,
                                             time_stepper_pt());

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------
 
 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here
 
 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();

 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Pin values for axial velocity components on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(2); // Axial (real)
     mesh_pt()->boundary_node_pt(b,n)->pin(3); // Axial (imaginary)

     // Pin values for radial and azimuthal velocity components on
     // all solid boundaries
     if(b==1 || b==3)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); // Radial (real)
       mesh_pt()->boundary_node_pt(b,n)->pin(1); // Radial (imaginary)
       mesh_pt()->boundary_node_pt(b,n)->pin(4); // Azimuthal (real)
       mesh_pt()->boundary_node_pt(b,n)->pin(5); // Azimuthal (imaginary)
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of elements in mesh
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   PERTURBED_ELEMENT* el_pt
    = dynamic_cast<PERTURBED_ELEMENT*>(mesh_pt()->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;

   // The mesh remains fixed
   el_pt->disable_ALE();

   // Set the azimuthal wavenumber
   el_pt->azimuthal_mode_number_pt() = &GlobalPhysicalVariables::k;

  } // End of loop over elements
 
 // Set the pressure in first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);
 
 // ------------------------------------------------------------
 // Set up interaction between base and perturbed state problems
 // ------------------------------------------------------------

 // First interaction (base state velocities)
 Multi_domain_functions::
  setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),Base_state_mesh_pt,0);

 // Second interaction (base state velocity derivatives w.r.t. space)
 Multi_domain_functions::
   setup_multi_domain_interaction<BASE_ELEMENT>
  (this,mesh_pt(),Base_state_mesh_pt,1);

 // Set up equation numbering scheme
 std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of perturbed state constructor



//======start_of_perturbed_state_set_initial_condition===================
/// Set the initial conditions to a "Poiseuille-style" profile
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the six velocity components
   for(unsigned i=0;i<6;i++)
    {
     // Get x and y position of node
     const double x = mesh_pt()->node_pt(n)->x(0);
     const double y = mesh_pt()->node_pt(n)->x(1);
     
     // Assign "Poiseuille-style" initial conditions
     const double value = 4*(x-Inner_cylinder_radius)
      *(x-Outer_cylinder_radius)*y*(y-Domain_height);
     
     mesh_pt()->node_pt(n)->set_value(i,value);
    }
  } // End of loop over nodes in mesh
 
 // Note that we do not need to set any history values as we are
 // using a self-starting timestepper
 
 // Set the boundary conditions
 set_boundary_conditions();

} // End of set_initial_condition for perturbed state



//======start_of_perturbed_state_set_boundary_conditions=================
/// Set the (homogeneous) boundary conditions for the perturbed state
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
set_boundary_conditions()
{
 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();
 
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);
   
   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the six velocity components
     for(unsigned i=0;i<6;i++)
      {
       // For the solid boundaries set all components to zero
       if(b==1 || b==3)
        {
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
        }
       // Otherwise set only the axial components to zero
       // (no penetration of symmetry boundary)
       else
        { 
         if(i==2 || i==3)
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0); 
          }
        }
      } // End of loop over velocity components
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
} // End of set_boundary_conditions for perturbed state



//==start_of_perturbed_state_doc_solution================================
/// Document the perturbed state solution
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
doc_solution(DocInfo& doc_info,const bool& output_soln)
{
 // Create vector of dofs
 DoubleVector dofs;
 
 // Get dofs
 this->get_dofs(dofs);
 
 // Get L2 norm of dof vector
 const double dof_norm = dofs.norm();
 
 // Document in trace file
 Trace_file << time_pt()->time() << " " << dof_norm << std::endl;

 // If desired, output solution to file
 if(output_soln)
  {
   ofstream some_file;
   char filename[256];
   
   // Set number of plot points (in each coordinate direction)
   const unsigned npts = 2;//5;
   
   // Open solution output file
   sprintf(filename,"%s/perturbed_soln_epsilon%2.1f_Re%4.2f_soln%i.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current,
           doc_info.number());
   some_file.open(filename);
   
   // Output solution to file
   mesh_pt()->output(some_file,npts);
   
   // Close solution output file
   some_file.close();
  }

} // End of doc_solution for perturbed state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=====start_of_stability_problem_class==================================
/// Container class for the time-periodic Taylor--Couette stability
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
class StabilityProblem
{

public:

 /// Constructor: Build base and perturbed state problems
 StabilityProblem(const unsigned& base_n_r,
                  const unsigned& base_n_z,
                  const unsigned& perturbed_n_r,
                  const unsigned& perturbed_n_z,
                  const double& radius_inner_cylinder,
                  const double& radius_outer_cylinder,
                  const double& l_z)
  {
   // Build base state problem
   Base_state_problem_pt = new BaseStateProblem<BASE_ELEMENT>
    (base_n_r,base_n_z,radius_inner_cylinder,radius_outer_cylinder,l_z);
   
   // Build perturbed state problem
   Perturbed_state_problem_pt = 
    new PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>
    (perturbed_n_r,perturbed_n_z,radius_inner_cylinder,
     radius_outer_cylinder,l_z,Base_state_problem_pt->mesh_pt());
  }
 
 /// Destructor (empty)
 ~StabilityProblem() {}

 /// Create trace files
 void create_trace_files(DocInfo& doc_info)
  {
   // Set up base and perturbed state trace files
   Base_state_problem_pt->create_trace_file(doc_info);
   Perturbed_state_problem_pt->create_trace_file(doc_info);

   // Set up stability problem trace files
   char filename[256];
   sprintf(filename,"%s/periodic_base_flow_trace_epsilon%2.1f_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current);
   Trace_file_periodic_base_flow.open(filename);

   sprintf(filename,"%s/power_method_trace_epsilon%2.1f_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::Epsilon,
           GlobalPhysicalVariables::MMC_Re_current);
   Trace_file_power_method.open(filename);
  }

 /// Initialise trace files (creates column headings)
 void initialise_trace_files()
  {
   // Initialise base and perturbed state trace files
   Base_state_problem_pt->initialise_trace_file();
   Perturbed_state_problem_pt->initialise_trace_file();

   // Initialise stability problem trace files
   Trace_file_periodic_base_flow << "time, norm, tolerance" << std::endl;
   Trace_file_power_method << "time, convergence_criterion, "
                           << "tolerance*abs(eigenvalue)" << std::endl;
  }

 /// Clear and close trace files
 void close_trace_files()
  {
   // Initialise base and perturbed state trace files
   Base_state_problem_pt->close_trace_file();
   Perturbed_state_problem_pt->close_trace_file();

   // Initialise stability problem trace files
   Trace_file_periodic_base_flow.clear();
   Trace_file_periodic_base_flow.close();
   Trace_file_power_method.clear();
   Trace_file_power_method.close();
  }

 /// Initialise timestep (single dt version)
 void initialise_dt(const double& dt)
  {
   Base_state_problem_pt->initialise_dt(dt);
   Perturbed_state_problem_pt->initialise_dt(dt);
  }

 /// Initialise timestep (vector of dts version)
 void initialise_dt(const Vector<double>& dt)
  {
   Base_state_problem_pt->initialise_dt(dt);
   Perturbed_state_problem_pt->initialise_dt(dt);
  }

 /// Set initial conditions
 void set_initial_condition()
  {
   Base_state_problem_pt->set_initial_condition();
   Perturbed_state_problem_pt->set_initial_condition();
  }

 /// Document base and perturbed state solutions
 void doc_solution(DocInfo& doc_info,
                   const bool& output_base_soln,
                   const bool& output_perturbed_soln)
  {
   Base_state_problem_pt->
    doc_solution(doc_info,output_base_soln);
   Perturbed_state_problem_pt->doc_solution(doc_info,output_perturbed_soln);
  }

 /// Integrate forwards in time with timestep dt for n_timesteps
 /// Record an entire period's worth of the perturbed state problem
 /// solution in perturbed_state_dofs_over_period
 void unsteady_run(const double& dt,
                   const unsigned& n_timesteps,
                   DocInfo& doc_info,
                   const bool& output_base_soln,
                   Vector<DoubleVector>& perturbed_state_dofs_over_period,
                   Vector<DoubleVector>& base_state_dofs_over_period);
 
 /// Integrate forwards in time with timestep dt for n_timesteps.
 /// Re-set the initial conditions for the perturbed state problem to
 /// be input_perturbed_state_dofs and, at the end of the timestepping
 /// loop, return the perturbed state problem's dofs
 void unsteady_run_for_power_method(
  const double& dt,const unsigned& n_timesteps,DocInfo& doc_info,
  const DoubleVector& input_perturbed_state_dofs,
  DoubleVector& output_perturbed_state_dofs,
  const bool& output_base_soln,
  Vector<DoubleVector>& perturbed_state_dofs_over_period,
  Vector<DoubleVector>& base_state_dofs_over_period)
  {
   // Re-set the starting dof vector for the perturbed state problem
   Perturbed_state_problem_pt->set_dofs(input_perturbed_state_dofs);
   
   // Perform standard unsteady run
   this->unsteady_run(dt,n_timesteps,doc_info,
                      output_base_soln,perturbed_state_dofs_over_period,
                      base_state_dofs_over_period);
   
   // Get dofs from perturbed state problem
   Perturbed_state_problem_pt->get_dofs(output_perturbed_state_dofs);
  }
 
 /// Integrate just the base state problem forwards in time with
 /// timestep dt for n_timesteps. Once a periodic state has been reached,
 /// set perturbed state problem's global time to be equal to the base
 /// state problem's global time and begin solving the perturbed state
 /// problem also. This function returns the number of periods taken to
 /// satisfy the "periodic state" exit criteria.
 unsigned set_up_periodic_base_flow(const double& dt,
                                    const unsigned& n_timesteps,
                                    DocInfo& doc_info,
                                    const double& tolerance);
                                    
 /// Perform a power method with a maximum of max_iter iterations. This
 /// function requires an initial guess for the dominant eigenvector
 /// "input"and a tolerance which controls at which point the dominant
 /// eigenvalue will be deemed converged. At this point the calculated
 /// value of this eigenvalue "calc_eigenvalue" and its corresponding
 /// eigenvector "input" will be returned. The return value of this
 /// function is the number of power method iterations which were
 /// performed.
 unsigned perform_power_method(
  const double& dt,const unsigned& n_timesteps,DocInfo& doc_info,
  const double& tolerance,const unsigned& max_iter,double& calc_eigenvalue,
  DoubleVector& input,Vector<DoubleVector>& base_state_dofs_over_period);

 /// Get perturbed state dofs
 void get_perturbed_state_problem_dofs(DoubleVector& d)
  {
   Perturbed_state_problem_pt->get_dofs(d);
  }

 /// (Re-)set the global time in both problems to zero
 void reset_global_time_to_zero()
  {
   Base_state_problem_pt->time_pt()->time() = 0.0;
   Perturbed_state_problem_pt->time_pt()->time() = 0.0;

   cout << "\nReset global time to zero in both problems.\n" << endl;
  }

 /// Access function for base state problem
 const BaseStateProblem<BASE_ELEMENT>* base_state_problem_pt()
  { return Base_state_problem_pt; }
 
 /// Access function for perturbed state problem
 const PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>*
 perturbed_state_problem_pt() { return Perturbed_state_problem_pt; }
 
 /// Switch on flag
 void set_first_power_method_iteration_flag_to_true()
  { First_power_method_iteration = true; }

 /// Switch off flag
 void set_first_power_method_iteration_flag_to_false()
  { First_power_method_iteration = false; }
 
private:

 /// (Re-)set the global time in the perturbed state problem to be
 /// equal to the global time in the base state problem
 void set_perturbed_state_time_equal_to_base_state_time()
  {
   // Get current global time according to base state problem
   const double base_state_time = Base_state_problem_pt->time_pt()->time();

   // Set perturbed state problem's global time
   Perturbed_state_problem_pt->time_pt()->time() = base_state_time;

   cout << "\nReset perturbed state problem's global time to "
        << base_state_time << endl;
  }

 /// Trace file which documents the convergence of the base flow to a
 /// periodic state
 ofstream Trace_file_periodic_base_flow;

 /// Trace file which documents the convergence of the power method
 ofstream Trace_file_power_method;

 /// Pointer to base state problem class
 BaseStateProblem<BASE_ELEMENT>* Base_state_problem_pt;
 
 /// Pointer to perturbed state problem class
 PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>*
 Perturbed_state_problem_pt;

 /// Flag to indicate whether or not we are doing the first iteration
 /// of the power method
 bool First_power_method_iteration;

}; // End of stability_problem class



//==start_of_stability_problem_unsteady_run==============================
/// Integrate forwards in time with timestep dt for n_timesteps
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
void StabilityProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::unsteady_run(
 const double& dt,const unsigned& n_timesteps,
 DocInfo& doc_info,const bool& output_base_soln,
 Vector<DoubleVector>& perturbed_state_dofs_over_period,
 Vector<DoubleVector>& base_state_dofs_over_period)
{
 // "Outer" timestepping loop
 for(unsigned i=0;i<n_timesteps;i++)
  {
   // Output timestep and global time
   std::cout << "  Timestep " << i+1 << " of " << n_timesteps
             << "; Time is " << Base_state_problem_pt->time_pt()->time()
             << " (base), " << Perturbed_state_problem_pt->time_pt()->time()
             << " (perturbed)" << std::endl;
   
   // On the first "timestep", perform n_timesteps BDF1 timesteps
   if(i==0)
    {
     // Set the perturbed state problem's timestepper to work in BDF1 mode
     Base_state_problem_pt->time_stepper_pt()->turn_on_bdf1_mode();
     Perturbed_state_problem_pt->time_stepper_pt()->turn_on_bdf1_mode();
     
     // Determine the BDF1 timestep
     const double dt_bdf1 = dt/n_timesteps;
     
     // Initialise timestep (also sets timestepper weights)
     this->initialise_dt(dt_bdf1);
     
     // BDF1 timestepping loop
     for(unsigned j=0;j<n_timesteps;j++)
      {
       
       cout << "    - BDF1 loop: Timestep " << j+1 << " of " << n_timesteps
            << "; Time (before advance) is " 
            << Base_state_problem_pt->time_pt()->time()
            << " (base), "
            << Perturbed_state_problem_pt->time_pt()->time()
            << " (perturbed)" << std::endl;
       
       // Take timestep
       if(First_power_method_iteration)
        {
         Base_state_problem_pt->unsteady_newton_solve(dt_bdf1);
         Base_state_problem_pt->get_dofs(base_state_dofs_over_period[j]);
        }
       else
        {
         Base_state_problem_pt
          ->assign_eigenvector_to_dofs(base_state_dofs_over_period[j]);
        }
       Perturbed_state_problem_pt->unsteady_newton_solve(dt_bdf1);

       // Document in trace file on all steps other than the last
       // (when documentation is performed in outer loop)
       if(j!=(n_timesteps-1))
        {
         this->doc_solution(doc_info,false,false);
        }

      } // End of BDF1 timestepping loop
     
     // Now set the timestepper back to BDF2 mode
     Base_state_problem_pt->time_stepper_pt()->turn_off_bdf1_mode();
     Perturbed_state_problem_pt->time_stepper_pt()->turn_off_bdf1_mode();
    }
   
   // On the second "timestep", use a timestep that's twice as big as
   // used during the BDF1 timestepping loop
   else if(i==1)
    {
     // Determine the timestep and number of timesteps for this loop
     const double dt_thisloop = 2*dt/n_timesteps;
     const double n_timesteps_thisloop = 0.5*n_timesteps;
     
     // Determine the timestep for the previous loop
     const double dt_prevloop = dt/n_timesteps;
     
     Vector<double> dt_vector;
     dt_vector.push_back(dt_prevloop);
     dt_vector.push_back(dt_thisloop);
     
     // Initialise timestep for first time around this timestepping loop
     // (also sets timestepper weights)
     this->initialise_dt(dt_vector);
     
     // Timestepping loop
     for(unsigned j=0;j<n_timesteps_thisloop;j++)
      {
       std::cout << "    - nstep/2 loop: Timestep " << j+1
                 << " of " << n_timesteps_thisloop
                 << "; Time (before advance) is " 
                 << Base_state_problem_pt->time_pt()->time()
                 << " (base), "
                 << Perturbed_state_problem_pt->time_pt()->time()
                 << " (perturbed)" << std::endl;
       
       // Take timestep
       if(First_power_method_iteration)
        {
         Base_state_problem_pt->unsteady_newton_solve(dt_thisloop);
         Base_state_problem_pt
          ->get_dofs(base_state_dofs_over_period[(n_timesteps+j)]);
        }
       else
        {
         Base_state_problem_pt->assign_eigenvector_to_dofs
          (base_state_dofs_over_period[(n_timesteps+j)]);
        }
       Perturbed_state_problem_pt->unsteady_newton_solve(dt_thisloop);
       
       // Document in trace file on all steps other than the last
       // (when documentation is performed in outer loop)
       if(j!=(n_timesteps_thisloop-1))
        {
         this->doc_solution(doc_info,false,false);
        }

       // On the first time through this loop, initialise the timestep
       // again with equal timesteps
       if(j==0) { this->initialise_dt(dt_thisloop); }
      }
    }
   
   // On the third "timestep", use a timestep that's five times as big as
   // used during the BDF1 timestepping loop
   else if(i==2)
    {
     // Determine the timestep and number of timesteps for this loop
     const double dt_thisloop = 5*dt/n_timesteps;
     const double n_timesteps_thisloop = 0.2*n_timesteps;
     
     // Determine the timestep for the previous loop
     const double dt_prevloop = 2*dt/n_timesteps;
     
     Vector<double> dt_vector;
     dt_vector.push_back(dt_prevloop);
     dt_vector.push_back(dt_thisloop);
     
     // Initialise timestep for first time around this timestepping loop
     // (also sets timestepper weights)
     this->initialise_dt(dt_vector);
     
     // Timestepping loop
     for(unsigned j=0;j<n_timesteps_thisloop;j++)
      {
       std::cout << "    - nstep/2 loop: Timestep " << j+1
                 << " of " << n_timesteps_thisloop
                 << "; Time (before advance) is " 
                 << Base_state_problem_pt->time_pt()->time()
                 << " (base), "
                 << Perturbed_state_problem_pt->time_pt()->time()
                 << " (perturbed)" << std::endl;
       
       // Take timestep
       if(First_power_method_iteration)
        {
         Base_state_problem_pt->unsteady_newton_solve(dt_thisloop);
         Base_state_problem_pt
          ->get_dofs(base_state_dofs_over_period[((n_timesteps*1.5)+j)]);
        }
       else
        {
         Base_state_problem_pt->assign_eigenvector_to_dofs
          (base_state_dofs_over_period[((n_timesteps*1.5)+j)]);
        }
       Perturbed_state_problem_pt->unsteady_newton_solve(dt_thisloop);
       
       // Document in trace file on all steps other than the last
       // (when documentation is performed in outer loop)
       if(j!=(n_timesteps_thisloop-1))
        {
         this->doc_solution(doc_info,false,false);
        }

       // On the first time through this loop, initialise the timestep
       // again with equal timesteps
       if(j==0) { this->initialise_dt(dt_thisloop); }
      }
    }
   
   // On the fourth "timestep", use the actual timestep dt but need to
   // specify past dts as changing step size from last time
   else if(i==3)
    {
     // Determine the timestep for the previous loop
     const double dt_prevloop = 5*dt/n_timesteps;
     
     Vector<double> dt_vector;
     dt_vector.push_back(dt_prevloop);
     dt_vector.push_back(dt);
     
     // Initialise timestep for first time around this timestepping loop
     // (also sets timestepper weights)
     this->initialise_dt(dt_vector);
     
     // Take timestep
     if(First_power_method_iteration)
      {
       Base_state_problem_pt->unsteady_newton_solve(dt);
       Base_state_problem_pt
        ->get_dofs(base_state_dofs_over_period
                   [((n_timesteps*1.5)+(n_timesteps*0.2))]);
      }
     else
      {
       Base_state_problem_pt->assign_eigenvector_to_dofs
        (base_state_dofs_over_period[((n_timesteps*1.5)+(n_timesteps*0.2))]);
      }
     Perturbed_state_problem_pt->unsteady_newton_solve(dt);
     
     // Initialise timestep back to equal timesteps of length dt
     // (also sets timestepper weights)
     this->initialise_dt(dt);
    }
   
   // On all subsequent timesteps perform normal BDF2 timestepping
   else
    {
     // Take timestep
     if(First_power_method_iteration)
      {
       Base_state_problem_pt->unsteady_newton_solve(dt);
       Base_state_problem_pt
        ->get_dofs(base_state_dofs_over_period
                   [((n_timesteps*1.5)+(n_timesteps*0.2)+(i-3))]);
      }
     else
      {
       Base_state_problem_pt->assign_eigenvector_to_dofs
        (base_state_dofs_over_period
         [((n_timesteps*1.5)+(n_timesteps*0.2)+(i-3))]);
      }
     Perturbed_state_problem_pt->unsteady_newton_solve(dt);
    }

   // At end of each "outer" timestep, document info in trace file
   // and also output base state solution (if flag is set appropriately)
   this->doc_solution(doc_info,output_base_soln,false);

   // Store perturbed state problems dofs so that eigenvector over one
   // period can be outputted at the end of the power method
   Perturbed_state_problem_pt->get_dofs(perturbed_state_dofs_over_period[i]);

   // Increment counter for solutions
   doc_info.number()++;

  } // End of loop over "outer timesteps"
} // End of unsteady run



//==start_of_set_up_periodic_base_flow===================================
/// Integrate just the base state problem forwards in time with
/// timestep dt for n_timesteps. Once a periodic state has been reached,
/// set perturbed state problem's global time to be equal to the base
/// state problem's global time and begin solving the perturbed state
/// problem also. This function returns the number of periods taken to
/// satisfy the "periodic state" exit criteria.
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
unsigned StabilityProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
set_up_periodic_base_flow(const double& dt,const unsigned& n_timesteps,
                          DocInfo& doc_info,const double& tolerance)
{
 // Check that base state timestepper is working in BDF2 mode
 if(Base_state_problem_pt->time_stepper_pt()->bdf1_mode())
  {
   std::cout << "Base state problem's timestepper is not in BDF2 mode."
             << "\n Exiting...\n" << std::endl;
   exit(1);
  }
 
 // Initialise timestep (also sets timestepper weights)
 this->initialise_dt(dt);

 // Set up storage for base state problem's dof vectors at the end of the
 // previous period and at the end of the current period
 Vector<double> previous_dofs;
 Vector<double> current_dofs;

 // Determine number of nodes in base state problem
 const unsigned n_node = Base_state_problem_pt->mesh_pt()->nnode();

 // Set up storage for norm of vector of differences between current
 // dofs and dofs from one period ago (initialise to 1.0)
 double norm_over_ndof = 1.0;
 
 // Set up storage for loop counter
 unsigned loop_counter = 0;

 // While this norm is greater than the base state convergence criterion
 while(norm_over_ndof > tolerance)
  {
   // Record the dofs corresponding to the azimuthal velocity at the end
   // of the previous period (equivalent to recording them at the start
   // of the current period, as we do here)
   for(unsigned n=0;n<n_node;n++)
    {
     previous_dofs.push_back
      (Base_state_problem_pt->mesh_pt()->node_pt(n)->value(2));
    }
   
   // Timestepping loop
   for(unsigned i=0;i<n_timesteps;i++)
    {
     // Output timestep and global time
     std::cout << "  Timestep " << i+1 << " of " << n_timesteps
               << "; Time is " << Base_state_problem_pt->time_pt()->time()
               << " (base)" << std::endl;
     
     // Take timestep
     Base_state_problem_pt->unsteady_newton_solve(dt);
     
     // Output solution
     Base_state_problem_pt->doc_solution(doc_info,false);
     
     // Increment counter for solutions 
     doc_info.number()++;
    }
   
   // Get dofs from base state problem at the end of this period
   for(unsigned n=0;n<n_node;n++)
    {
     current_dofs.push_back
      (Base_state_problem_pt->mesh_pt()->node_pt(n)->value(2));
    }
   
   // Find difference between each current dof and the corresponding 
   // previous dof from one period ago
   for(unsigned i=0;i<n_node;i++) { current_dofs[i] -= previous_dofs[i]; }
   
   // Evaluate L2 norm of this vector over number of dofs
   norm_over_ndof = 0.0;
   for(unsigned i=0;i<n_node;i++)
    {
     norm_over_ndof += (current_dofs[i]*current_dofs[i]);
    }
   norm_over_ndof = sqrt(norm_over_ndof);
   norm_over_ndof = norm_over_ndof/n_node;
   
   ios_base::fmtflags flags =  cout.flags(); // Save old flags
   cout << "    ==>  norm_over_ndof = " << scientific << norm_over_ndof
        << ", tolerance = " << tolerance << endl;
   cout.flags(flags);  // Set the flags to the way they were

   // Document norm and tolerance in trace file
   Trace_file_periodic_base_flow << Base_state_problem_pt->time_pt()->time()
                                 << " " << norm_over_ndof
                                 << " " << tolerance << std::endl;

   // Clear the Vector<double>s
   previous_dofs.clear();
   current_dofs.clear();

   // Increment loop counter
   loop_counter++;
   
  } // End of while loop
 
 cout << "\nI have broken out of the base state loop after "
      << loop_counter << " periods.\n" << endl;

 // (Re-)set the perturbed state problem's global time to correspond
 // to the base state problem's global time
 this->set_perturbed_state_time_equal_to_base_state_time();
 
 // Return the number of periods required to reach a periodic flow
 return loop_counter;

} // End of set_up_periodic_base_flow



//==start_of_perform_power_method========================================
/// Perform a power method with a maximum of max_iter iterations. This
/// function requires an initial guess for the dominant eigenvector
/// "input"and a tolerance which controls at which point the dominant
/// eigenvalue will be deemed converged. At this point the calculated
/// value of this eigenvalue "calc_eigenvalue" and its corresponding
/// eigenvector "input" will be returned. The return value of this
/// function is the number of power method iterations which were
/// performed.
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
unsigned StabilityProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
perform_power_method(const double& dt,
                     const unsigned& n_timesteps,
                     DocInfo& doc_info,
                     const double& tolerance,
                     const unsigned& max_iter,
                     double& calc_eigenvalue,
                     DoubleVector& input,
                     Vector<DoubleVector>& base_state_dofs_over_period)
{
 // Initialise output_base_state_solution flag to true
 bool output_base_state_solution = true;

 // Reset the solution counter to zero
 doc_info.number()=0;

 // Determine number of degrees of freedom in perturbed state problem
 const unsigned n_dof_perturbed = Perturbed_state_problem_pt->ndof();

 // Set up external storage for the dofs over one period
 Vector<DoubleVector> dofs_over_period(n_timesteps);

 // Set flag to true
 this->set_first_power_method_iteration_flag_to_true();

 // Begin power method loop
 for(unsigned p=0;p<max_iter;p++)
  {
   cout << "\nThis is power method loop number " << p+1 
        << " of a maximum " << max_iter << "." << endl;
   
   // Normalise input vector
   const double inv_norm = 1.0/input.norm(); 
   for(unsigned i=0;i<n_dof_perturbed;i++) { input[i] *= inv_norm; }
   
   // Set up storage for output vector from Arnoldi method
   DoubleVector output;

   // Perform unsteady run over this period
   unsteady_run_for_power_method(dt,n_timesteps,doc_info,input,output,
                                 output_base_state_solution,dofs_over_period,
                                 base_state_dofs_over_period);

   // Calculate eigenvalue (scalar product of input and output vectors)
   calc_eigenvalue = input.dot(output);

   // Calculate L2 norm of (output - eigenvalue*input)
   double convergence_criterion = 0.0;
   for(unsigned i=0;i<n_dof_perturbed;i++)
    {
     convergence_criterion +=
      (output[i] - calc_eigenvalue*input[i])
      *(output[i] - calc_eigenvalue*input[i]);
    }
   convergence_criterion = sqrt(convergence_criterion);
   
   cout << "\nConvergence criterion     = " << convergence_criterion << endl;
   cout << "tolerance*abs(eigenvalue) = " << tolerance*abs(calc_eigenvalue)
        << endl;
   
   // Document convergence criterion and tolerance*eigenvalue in trace file
   Trace_file_power_method << Perturbed_state_problem_pt->time_pt()->time()
                           << " " << convergence_criterion << " "
                           << tolerance*abs(calc_eigenvalue) << std::endl;

   // We only wish to output the base state solution during the first
   // iteration of the power method (so that we can see what it looks like
   // when converged)
   output_base_state_solution = false;

   // Check for convergence
   if(convergence_criterion <= tolerance*abs(calc_eigenvalue))
    {
     cout << "\nPower method has converged to within a tolerance of " 
          << tolerance << "\nin " << p+1
          << " iterations of the power method.\n" << endl;

     // Reset the solution counter to zero
     doc_info.number()=0;

     // Output the dominant eigenvector over an entire period
     for(unsigned i=0;i<n_timesteps;i++)
      {
       // Populate perturbed state problem with dofs which have been stored
       // over the last period
       Perturbed_state_problem_pt->
        assign_eigenvector_to_dofs(dofs_over_period[i]);

       // Output the perturbed state solution
       this->doc_solution(doc_info,false,true);

       // Increment counter for solutions
       doc_info.number()++;
      }

     // Return the number of iterations of the power method
     return p+1;
    }
   else
    {
     input = output;
     output.clear();
    }

   // Set flag to false
   this->set_first_power_method_iteration_flag_to_false();
  }

 // If we reach here, the power method has failed to converge in max_iter
 // iterations.
 cout << "\nPower method has failed to converge to within a tolerance of " 
      << tolerance << "\nin " << max_iter
      << " iterations of the power method" << endl;
 
 // Reset the solution counter to zero
 doc_info.number()=0;
 
 // Output the dominant eigenvector over an entire period
 for(unsigned i=0;i<n_timesteps;i++)
  {
   // Populate perturbed state problem with dofs which have been stored
   // over the last period
   Perturbed_state_problem_pt->
    assign_eigenvector_to_dofs(dofs_over_period[i]);
   
   // Output the perturbed state solution
   this->doc_solution(doc_info,false,true);
   
   // Increment counter for solutions
   doc_info.number()++;
  }
 
 return max_iter;
 
} // End of perform_power_method



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_main========================================================
/// Driver for time-periodic Taylor--Couette problem
//=======================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set number of elements in radial (r) direction
 const unsigned base_n_r = 16;
 const unsigned perturbed_n_r = 16;

 // Set number of elements in axial (z) direction
 const unsigned base_n_z = 16;
 const unsigned perturbed_n_z = 16;

 // Cache radius ratio
 const double eta = GlobalPhysicalVariables::Eta;

 // Determine tadius of inner cylinder
 const double radius_inner_cylinder = eta/(1.0-eta);

 // Determine radius of outer cylinder
 const double radius_outer_cylinder = 1.0/(1.0-eta);

 // Length in axial (z) direction (set to fit in half a wavelength)
 const double l_z = MathematicalConstants::Pi/GlobalPhysicalVariables::a;

 // Set direction of gravity (vertically downwards)
 GlobalPhysicalVariables::G[0] = 0.0;
 GlobalPhysicalVariables::G[1] = -1.0;
 GlobalPhysicalVariables::G[2] = 0.0;
 
 // Determine period of the base state problem
 const double T = 
  (2*MathematicalConstants::Pi)/GlobalPhysicalVariables::AngularFrequency;

 // Set number of timesteps per period
 unsigned nstep = 40;

 // Set maximum number of iterations of the power method
 unsigned max_iter = 1000;

 // Set tolerance for power method convergence (this is usually machine
 // precision according to Bai, Demmel, Dongarra, Ruhe & Van der Vorst)
 const double power_method_tolerance = 1e-8;

 // Set tolerance for L2 norm of difference in dof vectors from one
 // period to the next for defining the point at which the base flow
 // has settled down to a periodic state.
 double base_flow_norm_tolerance = 1e-6;

 // Cache upper and lower limits and number of steps for loop over MMC_Re
 double lower = GlobalPhysicalVariables::MMC_Re_lower_limit;
 double upper = GlobalPhysicalVariables::MMC_Re_upper_limit;
 double n_param_steps = GlobalPhysicalVariables::Nparam_steps;

 // If we are doing a validation run, only compute one value of the
 // Reynolds number, use a lower number of timesteps per period,
 // accept a lower base flow convergence tolerance and insist that
 // only one iteration of the power method is performed
 if(CommandLineArgs::Argc>1)
  {
   n_param_steps = 1;
   nstep = 20;
   base_flow_norm_tolerance = 1e-2;
   max_iter = 1;
  }

 // Check that nstep is exactly divisible by both 2 and 5
 // (This is necessary for the self_starting_BDF2 timestepper)
 if(nstep%2 != 0 || nstep%5 !=0)
  {
   cout << "\nThe self-starting BDF2 timestepper requires that the number"
        << "\nof timesteps be exactly divisible by both 2 and 5. This is"
        << "\nnot the case for nstep = " << nstep << ". Exiting...\n"
        << endl;
   exit(2);
  }

 // Calculate length of timestep
 const double dt = T/nstep;

 // -----------------------------------------
 // LinearisedAxisymmetricQTaylorHoodElements
 // -----------------------------------------
 {
  cout << "\nDoing LinearisedAxisymmetricQTaylorHoodElements\n" << endl;

  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<AxisymmetricQTaylorHoodElement,
   LinearisedAxisymmetricQTaylorHoodMultiDomainElement>
   problem(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,
           radius_inner_cylinder,radius_outer_cylinder,l_z);
  
  // Set up labels for output
  DocInfo doc_info;
 
  // Output directory
  doc_info.set_directory("RESLT_TH");
  
  // Initialise timestep (also sets weights for all timesteppers)
  problem.initialise_dt(dt);
  
  // Set initial conditions
  problem.set_initial_condition();
  
  // Set up storage for dominant eigenvector and eigenvalue
  DoubleVector eigenvector;
  double eigenvalue = 0.0;
  
  // Get initial guess for eigenvector (these are the initial conditions
  // set up in PerturbedStateProblem::set_initial_conditions())
  problem.get_perturbed_state_problem_dofs(eigenvector);
  
  // Create and initialise global trace file
  ofstream global_trace;
  char filename[256];
  sprintf(filename,"%s/global_trace_epsilon%2.1f.dat",
          doc_info.directory().c_str(),
          GlobalPhysicalVariables::Epsilon);
  global_trace.open(filename);
  global_trace << "MMC_Re, dominating eigenvalue, "
               << "n_periods_to_set_up_periodic_base_flow, "
               << "n_power_method_iterations" << std::endl;
  
  // Begin loop over Reynolds number
  for(unsigned param=0;param<n_param_steps;param++)
   {
    // Set MMC_Re_current to the correct value
    if(lower==upper || n_param_steps<=1)
     {
      GlobalPhysicalVariables::MMC_Re_current = lower;
     }
    else
     {
      GlobalPhysicalVariables::MMC_Re_current =
       lower + (((upper-lower)/(n_param_steps-1))*param);
     }
    
    cout << "\n====================================================" << endl;
    cout << "Beginning parameter run " << param+1 << " of "
         << n_param_steps << ": MMC_Re = "
         << GlobalPhysicalVariables::MMC_Re_current << endl;
    cout << "====================================================" << endl;
    
    // Reset the global time (of both problems!) to zero
    problem.reset_global_time_to_zero();
    
    // Initialise counter for solutions
    doc_info.number()=0;
    
    // Create trace files
    problem.create_trace_files(doc_info);
    
    // Initialise trace files
    problem.initialise_trace_files();
    
    // Output initial conditions
    problem.doc_solution(doc_info,false,false);
    
    // Increment counter for solutions 
    doc_info.number()++;
    
    // Set up periodic base flow and record how many period it takes
    const unsigned n_periods_to_set_up_periodic_base_flow
     = problem.set_up_periodic_base_flow(dt,nstep,doc_info,
                                         base_flow_norm_tolerance);
    
    // Set up storage for base state dofs over one period
    const unsigned n_entries_required = (27*(nstep/10))-3;
    Vector<DoubleVector> base_state_dofs_over_period(n_entries_required);
    
    // Perform power method and store number of iterations
    const unsigned n_power_method_iterations =
     problem.perform_power_method(dt,nstep,doc_info,power_method_tolerance,
                                  max_iter,eigenvalue,eigenvector,
                                  base_state_dofs_over_period);
    
    cout << "\nDominating eigenvalue is " << eigenvalue << endl;
    
    // The corresponding dominating eigenvector is "input"
    
    // Document in the global trace file
    global_trace << GlobalPhysicalVariables::MMC_Re_current << " "
                 << eigenvalue << " "
                 << n_periods_to_set_up_periodic_base_flow << " "
                 << n_power_method_iterations << std::endl;
    
    // The initial guess for the next power method will be the eigenvector
    // calculated by this power method
    
    // Clear and close trace files
    problem.close_trace_files();
    
   } // End loop over Reynolds number
 }

 // ----------------------------------------------
 // LinearisedAxisymmetricQCrouzeixRaviartElements
 // ----------------------------------------------
 {
  cout << "\nDoing LinearisedAxisymmetricQCrouzeixRaviartElements\n" << endl;

  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<AxisymmetricQCrouzeixRaviartElement,
   LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement>
   problem(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,
           radius_inner_cylinder,radius_outer_cylinder,l_z);
  
  // Set up labels for output
  DocInfo doc_info;
 
  // Output directory
  doc_info.set_directory("RESLT_CR");
  
  // Initialise timestep (also sets weights for all timesteppers)
  problem.initialise_dt(dt);
  
  // Set initial conditions
  problem.set_initial_condition();
  
  // Set up storage for dominant eigenvector and eigenvalue
  DoubleVector eigenvector;
  double eigenvalue = 0.0;
  
  // Get initial guess for eigenvector (these are the initial conditions
  // set up in PerturbedStateProblem::set_initial_conditions())
  problem.get_perturbed_state_problem_dofs(eigenvector);
  
  // Create and initialise global trace file
  ofstream global_trace;
  char filename[256];
  sprintf(filename,"%s/global_trace_epsilon%2.1f.dat",
          doc_info.directory().c_str(),
          GlobalPhysicalVariables::Epsilon);
  global_trace.open(filename);
  global_trace << "MMC_Re, dominating eigenvalue, "
               << "n_periods_to_set_up_periodic_base_flow, "
               << "n_power_method_iterations" << std::endl;
  
  // Begin loop over Reynolds number
  for(unsigned param=0;param<n_param_steps;param++)
   {
    // Set MMC_Re_current to the correct value
    if(lower==upper || n_param_steps<=1)
     {
      GlobalPhysicalVariables::MMC_Re_current = lower;
     }
    else
     {
      GlobalPhysicalVariables::MMC_Re_current =
       lower + (((upper-lower)/(n_param_steps-1))*param);
     }
    
    cout << "\n====================================================" << endl;
    cout << "Beginning parameter run " << param+1 << " of "
         << n_param_steps << ": MMC_Re = "
         << GlobalPhysicalVariables::MMC_Re_current << endl;
    cout << "====================================================" << endl;
    
    // Reset the global time (of both problems!) to zero
    problem.reset_global_time_to_zero();
    
    // Initialise counter for solutions
    doc_info.number()=0;
    
    // Create trace files
    problem.create_trace_files(doc_info);
    
    // Initialise trace files
    problem.initialise_trace_files();
    
    // Output initial conditions
    problem.doc_solution(doc_info,false,false);
    
    // Increment counter for solutions 
    doc_info.number()++;
    
    // Set up periodic base flow and record how many period it takes
    const unsigned n_periods_to_set_up_periodic_base_flow
     = problem.set_up_periodic_base_flow(dt,nstep,doc_info,
                                         base_flow_norm_tolerance);
    
    // Set up storage for base state dofs over one period
    const unsigned n_entries_required = (27*(nstep/10))-3;
    Vector<DoubleVector> base_state_dofs_over_period(n_entries_required);
    
    // Perform power method and store number of iterations
    const unsigned n_power_method_iterations =
     problem.perform_power_method(dt,nstep,doc_info,power_method_tolerance,
                                  max_iter,eigenvalue,eigenvector,
                                  base_state_dofs_over_period);
    
    cout << "\nDominating eigenvalue is " << eigenvalue << endl;
    
    // The corresponding dominating eigenvector is "input"
    
    // Document in the global trace file
    global_trace << GlobalPhysicalVariables::MMC_Re_current << " "
                 << eigenvalue << " "
                 << n_periods_to_set_up_periodic_base_flow << " "
                 << n_power_method_iterations << std::endl;
    
    // The initial guess for the next power method will be the eigenvector
    // calculated by this power method
    
    // Clear and close trace files
    problem.close_trace_files();
    
   } // End loop over Reynolds number
 }

} // End of main







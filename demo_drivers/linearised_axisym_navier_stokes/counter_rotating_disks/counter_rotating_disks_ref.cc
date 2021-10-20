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
// Driver for validation problem for refineable linearised axisymmetric
// Navier-Stokes equations in cylindrical polar coordinates which
// replicates the investigation of C. Nore, M. Tartar, O. Daube and L.S.
// Tuckerman's paper "Survey of instability thresholds of flow between
// exactly counter-rotating disks" (2003).

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

 /// Strouhal number
 double St = 1.0;

 /// Inverse of Froude number
 double InvFr = 0.0; // (Fr = inf)

 /// Current value of the Reynolds number etc.
 /// Note: Do not specify these as they are determined from the
 /// information below
 double Re_current;
 double ReSt_current;
 double ReInvFr_current;

 /// Loop information
 double Re_lower_limit = 300;
 double Re_upper_limit = 302;
 unsigned Nparam_steps = 5;

 /// Aspect ratio (cylinder height / cylinder radius)
 double Gamma = 1.0;

 /// Azimuthal mode number k in e^ik(theta) decomposition
 int k = 2;

 /// Direction of gravity
 Vector<double> G(3);

} // End of namespace



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//=====start_of_base_state_problem_class=================================
/// Base state problem class for Tuckerman counter-rotating lids problem
//=======================================================================
template<class BASE_ELEMENT>
class BaseStateProblem : public Problem
{
 
public:
 
 /// Constructor
 BaseStateProblem(const unsigned& n_r,
                  const unsigned& n_z,
                  const double& domain_height);

 /// Destructor (empty)
 ~BaseStateProblem() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 ///  Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep()
  {
   set_boundary_conditions();
  }

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);

  } // End of actions_after_adapt

 ///  Set initial condition (incl previous timesteps)
 void set_initial_condition();

 ///  Set the boundary conditions
 void set_boundary_conditions();

 /// Access function for the specific mesh
 RefineableRectangularQuadMesh<BASE_ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableRectangularQuadMesh<BASE_ELEMENT>*>
    (Problem::mesh_pt());
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
   sprintf(filename,"%s/base_trace_k%i_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::k,
           GlobalPhysicalVariables::Re_current);
   Trace_file.open(filename);
  }

 /// Initialise trace file (print column headings)
 void initialise_trace_file()
  {
   Trace_file << "time, norm_of_dof_vector" << std::endl;
  }

 /// Clear and close trace file
 void close_trace_file() { Trace_file.clear(); Trace_file.close(); }

 /// Access function for trace file
 ofstream& trace_file() { return Trace_file; }

 void pass_updated_nondim_parameters_to_elements()
  {
   // Determine number of elements in the mesh
   const unsigned n_element = this->mesh_pt()->nelement();

   // Loop over the elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     BASE_ELEMENT *el_pt=dynamic_cast<BASE_ELEMENT*>(mesh_pt()->element_pt(e));

     // Set the Reynolds number
     el_pt->re_pt() = &GlobalPhysicalVariables::Re_current;
     
     // Set the Womersley number
     el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt_current;
     
     // Set the product of the Reynolds number and the inverse of the
     // Froude number
     el_pt->re_invfr_pt() = &GlobalPhysicalVariables::ReInvFr_current;
    }
  }

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
/// Constructor for base state Tuckerman counter-rotating lids problem
//=======================================================================
template<class BASE_ELEMENT>
BaseStateProblem<BASE_ELEMENT>::BaseStateProblem(const unsigned& n_r,
                                                 const unsigned& n_z,
                                                 const double& domain_height)
{
 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new SelfStartingBDF2);

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<BASE_ELEMENT>(n_r,n_z,0.0,1.0,0.0,
                                                  domain_height,
                                                  time_stepper_pt());
 
 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

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
     // Pin radial velocity component on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(0);

     // Pin azimuthal velocity component on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(2);

     // Pin axial velocity component on all SOLID boundaries
     if(b!=3) { mesh_pt()->boundary_node_pt(b,n)->pin(1); }
    }
  } // End of loop over mesh boundaries

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

   // Note that at this point we do not know the value of the Reynolds
   // number. Therefore we update Re, ReSt, ReInvFr using the function
   // "pass_updated_nondim_parameters_to_elements()"

   // Set the direction of gravity
   el_pt->g_pt() = &GlobalPhysicalVariables::G;

   // The mesh remains fixed
   el_pt->disable_ALE();
   
  } // End of loop over elements

 // Pin redundant pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());

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
void BaseStateProblem<BASE_ELEMENT>::set_boundary_conditions()
{
 // Set fraction along lid, "a",  after which to apply "smoothing"
 // of boundary conditions (0 < a < 1)
 const double a = 0.96;

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
     // Get the radial value
     double r = mesh_pt()->boundary_node_pt(b,n)->x(0);

     // Provide storage for azimuthal velocity value
     double azi_vel_value = 0.0;
     
     // Apply "smoothing" if r > a
     if(r>a) { azi_vel_value = a*(1.0-r)/(1.0-a); }
     else { azi_vel_value = r; }

     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       switch(b)
        {
         // Outer wall (all components = 0)
        case 1:
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
         break;

         // Top lid (azimuthal component = r, others = 0)
        case 2:
         if(i==2)
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,azi_vel_value);
          }
         else
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
          }
         break;

         // Bottom lid (azimuthal component = -r, others = 0)
        case 0:
         if(i==2)
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,-azi_vel_value);
          }
         else
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
          }
         break;

         // Symmetry boundary (axial component not set, others = 0)
        case 3:
         if(i!=1)
          {
           mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0); 
          }
         break;
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
   const unsigned npts = 5;
   
   // Open solution output file
   sprintf(filename,"%s/base_soln_k%i_Re%4.2f_soln%i.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::k,
           GlobalPhysicalVariables::Re_current,
           doc_info.number());
   some_file.open(filename);
   
   // Output solution to file
   mesh_pt()->output(some_file,npts);
   
   // Close solution output file
   some_file.close();
  }
} // End of doc_solution for base state



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//==start_of_perturbed_state_problem_class===============================
/// Perturbed state problem class for Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
class PerturbedStateProblem : public Problem
{
 
public:
 
 /// Constructor
 PerturbedStateProblem
 (const unsigned& n_r,
  const unsigned& n_z,
  const double& domain_height,
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

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableLinearisedAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableLinearisedAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Fix pressure ONLY when looking for axisymmetric modes
   if(GlobalPhysicalVariables::k==0) { fix_pressure(0,0,0.0); }

   // Set external elements for the multi-domain solution

   // First interaction (base state velocities)
   Multi_domain_functions::
    setup_multi_domain_interaction<BASE_ELEMENT>
    (this,mesh_pt(),Base_state_mesh_pt,0);

   // Second interaction (base state velocity derivatives w.r.t. space)
   Multi_domain_functions::
    setup_multi_domain_interaction<BASE_ELEMENT>
    (this,mesh_pt(),Base_state_mesh_pt,1);

  } // End of actions_after_adapt

 /// Set initial conditions to a "Poiseuille-style" profile
 void set_initial_condition();

 /// Set up the (homogeneous) boundary conditions
 void set_boundary_conditions();
 
 /// Access function for the specific mesh
 RefineableRectangularQuadMesh<PERTURBED_ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableRectangularQuadMesh<PERTURBED_ELEMENT>*>
    (Problem::mesh_pt());
  }
 
 /// Access function for the base state mesh
 RefineableRectangularQuadMesh<BASE_ELEMENT>* base_state_mesh_pt()
  {
   return dynamic_cast<RefineableRectangularQuadMesh<BASE_ELEMENT>*>
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
   sprintf(filename,"%s/perturbed_trace_k%i_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::k,
           GlobalPhysicalVariables::Re_current);
   Trace_file.open(filename);
  }

 /// Initialise trace file (print column headings)
 void initialise_trace_file()
  {
   Trace_file << "time, norm_of_dof_vector" << std::endl;
  }

 /// Clear and close trace file
 void close_trace_file() { Trace_file.clear(); Trace_file.close(); }

 /// Access function for trace file
 ofstream& trace_file() { return Trace_file; }

 void pass_updated_nondim_parameters_to_elements()
  {
   // Determine number of elements in mesh
   const unsigned n_element = this->mesh_pt()->nelement();

   // Loop over the elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PERTURBED_ELEMENT* el_pt
      = dynamic_cast<PERTURBED_ELEMENT*>(mesh_pt()->element_pt(e));

     // Set the Reynolds number
     el_pt->re_pt() = &GlobalPhysicalVariables::Re_current;
     
     // Set the Womersley number
     el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt_current;
    }
  }

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
/// Constructor for perturbed state Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
PerturbedStateProblem(const unsigned& n_r,
                      const unsigned& n_z,
                      const double& domain_height,
                      Mesh* external_mesh_pt)
 : Base_state_mesh_pt(external_mesh_pt),
   Domain_height(domain_height)
{
 // Be less verbose during newton solve
 Problem::disable_info_in_newton_solve();

 // Be less verbose about linear solve timings
 linear_solver_pt()->disable_doc_time();
 
 // Tell problem that it is linear (avoid doing unnecessary checks)
 Problem::Problem_is_nonlinear = false;

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new SelfStartingBDF2);
 
 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<PERTURBED_ELEMENT>(n_r,n_z,0.0,1.0,0.0,
                                                       domain_height,
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
     // On all solid boundaries, pin all velocity components
     if(b!=3)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); // Radial (real)
       mesh_pt()->boundary_node_pt(b,n)->pin(1); // Radial (imaginary)
       mesh_pt()->boundary_node_pt(b,n)->pin(2); // Axial (real)
       mesh_pt()->boundary_node_pt(b,n)->pin(3); // Axial (imaginary)
       mesh_pt()->boundary_node_pt(b,n)->pin(4); // Azimuthal (real)
       mesh_pt()->boundary_node_pt(b,n)->pin(5); // Azimuthal (imaginary)
      }
     // On symmetry boundary, pin only radial and azimuthal velocity
     // components, unless we are looking at the k=1 azimuthal mode.
     // In this case, do not pin any velocity components on this boundary
     else
      {
       if(GlobalPhysicalVariables::k!=1)
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(0); // Radial (real)
         mesh_pt()->boundary_node_pt(b,n)->pin(1); // Radial (imaginary)
         mesh_pt()->boundary_node_pt(b,n)->pin(4); // Azimuthal (real)
         mesh_pt()->boundary_node_pt(b,n)->pin(5); // Azimuthal (imaginary)
        }
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

   // Note that at this point we do not know the value of the Reynolds
   // number. Therefore we update Re and ReSt using the function
   // "pass_updated_nondim_parameters_to_elements()"

   // The mesh remains fixed
   el_pt->disable_ALE();

   // Set the azimuthal wavenumber
   el_pt->azimuthal_mode_number_pt() = &GlobalPhysicalVariables::k;

  } // End of loop over elements
 
 // Pin redundant pressure dofs
 RefineableLinearisedAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());

 // Fix pressure ONLY when looking for axisymmetric modes
 if(GlobalPhysicalVariables::k==0) { fix_pressure(0,0,0.0); }
 
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
     const double value = 4.0*x*(x-1.0)*y*(y-Domain_height);
     
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
     // For the solid boundaries set all components to zero
     if(b!=3)
      {
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,0,0.0); // Radial
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,1,0.0); // Radial
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,2,0.0); // Axial
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,3,0.0); // Axial
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,4,0.0); // Azimuthal
       mesh_pt()->boundary_node_pt(b,n)->set_value(0,5,0.0); // Azimuthal
      }
     // On symmetry boundary, set only radial and azimuthal components to
     // zero, unless we are looking at the k=1 azimuthal mode. In this
     // case, do not set the value of any of the velocity components
     else
      {
       if(GlobalPhysicalVariables::k!=1)
        {
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,0,0.0); // Radial
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,1,0.0); // Radial
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,4,0.0); // Azimuthal
         mesh_pt()->boundary_node_pt(b,n)->set_value(0,5,0.0); // Azimuthal
        }
      }
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
   const unsigned npts = 5;
   
   // Open solution output file
   sprintf(filename,"%s/perturbed_soln_k%i_Re%4.2f_soln%i.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::k,
           GlobalPhysicalVariables::Re_current,
           doc_info.number());
   some_file.open(filename);
   
   // Output solution to file
   mesh_pt()->output(some_file,npts);
   
   // Close solution output file
   some_file.close();
  }
} // End of doc_solution for perturbed state



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//=====start_of_stability_problem_class==================================
/// Container class for the Tuckerman counter-rotating lids stability
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
                  const double& domain_height)
  {
   // Build base state problem
   Base_state_problem_pt = new BaseStateProblem<BASE_ELEMENT>
    (base_n_r,base_n_z,domain_height);
   
   // Build perturbed state problem
   Perturbed_state_problem_pt = new PerturbedStateProblem
    <BASE_ELEMENT,PERTURBED_ELEMENT>
    (perturbed_n_r,perturbed_n_z,domain_height,
     Base_state_problem_pt->mesh_pt());
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
   sprintf(filename,"%s/power_method_trace_k%i_Re%4.2f.dat",
           doc_info.directory().c_str(),
           GlobalPhysicalVariables::k,
           GlobalPhysicalVariables::Re_current);
   Trace_file_power_method.open(filename);
  }

 /// Initialise trace files (creates column headings)
 void initialise_trace_files()
  {
   // Initialise base and perturbed state trace files
   Base_state_problem_pt->initialise_trace_file();
   Perturbed_state_problem_pt->initialise_trace_file();

   // Initialise stability problem trace files
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
 void doc_solution(DocInfo& doc_info,const bool& output_base_soln,
                   const bool& output_perturbed_soln)
  {
   Base_state_problem_pt->doc_solution(doc_info,output_base_soln);
   Perturbed_state_problem_pt->doc_solution(doc_info,output_perturbed_soln);
  }

 /// Pass updated nondimensional parameters to the elements of both
 /// problems
 void pass_updated_nondim_parameters_to_elements()
  {
   Base_state_problem_pt->pass_updated_nondim_parameters_to_elements();
   Perturbed_state_problem_pt->pass_updated_nondim_parameters_to_elements();
  }



 /// Perform a power method with a maximum of max_iter iterations. This
 /// function requires an initial guess for the dominant eigenvector
 /// "input"and a tolerance which controls at which point the dominant
 /// eigenvalue will be deemed converged. At this point the calculated
 /// value of this eigenvalue "calc_eigenvalue" and its corresponding
 /// eigenvector "input" will be returned. The return value of this
 /// function is the number of power method iterations which were
 /// performed.
 unsigned perform_power_method(const double& dt,
                               DocInfo& doc_info,
                               const double& tolerance,
                               const unsigned& max_iter,
                               double& calc_eigenvalue,
                               DoubleVector& input);

 // Solve the steady base flow problem
 void base_flow_steady_newton_solve()
  {
   Base_state_problem_pt->set_boundary_conditions();
   Base_state_problem_pt->steady_newton_solve();
  }

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
 BaseStateProblem<BASE_ELEMENT>*
 base_state_problem_pt() const { return Base_state_problem_pt; }
 
 /// Access function for perturbed state problem
 PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>*
 perturbed_state_problem_pt() const { return Perturbed_state_problem_pt; }
 
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

}; // End of stability_problem class



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
                     DocInfo& doc_info,
                     const double& tolerance,
                     const unsigned& max_iter,
                     double& calc_eigenvalue,
                     DoubleVector& input)
{
 // Initialise output_base_state_solution flag to true
 bool output_base_state_solution = true;

 // Reset the solution counter to zero
 doc_info.number()=0;

 // Determine number of degrees of freedom in perturbed state problem
 const unsigned n_dof_perturbed = Perturbed_state_problem_pt->ndof();

 // Begin power method loop
 for(unsigned p=0;p<max_iter;p++)
  {
   cout << "\nThis is power method loop number " << p+1 
        << " of a maximum " << max_iter << "." << endl;

   // Normalise input vector
   const double inv_norm = 1.0/input.norm();
   for(unsigned i=0;i<n_dof_perturbed;i++) { input[i] *= inv_norm; }

   // Set up storage for output vector
   DoubleVector output;

   // Re-set the starting dof vector for the perturbed state problem
   Perturbed_state_problem_pt->set_dofs(input);

   // Set the perturbed state problem's timestepper to work in BDF1 mode
   Perturbed_state_problem_pt->time_stepper_pt()->turn_on_bdf1_mode();
 
   // Initialise timestep (also sets timestepper weights)
   this->initialise_dt(dt);

   // Take timestep
   Perturbed_state_problem_pt->unsteady_newton_solve(dt);

   // Doc info in trace file and output base state solution
   // (if flag is set appropriately)
   this->doc_solution(doc_info,output_base_state_solution,false);

   // Increment counter for solutions
   doc_info.number()++;

   // Get dofs from perturbed state problem
   Perturbed_state_problem_pt->get_dofs(output);

   // Calculate eigenvalue (scalar product of input and output vectors)
   calc_eigenvalue = input.dot(output);

   // Calculate L2 norm of (output - eigenvalue*input)
   // divided by the square root of the number of degrees of freedom
   double convergence_criterion = 0.0;
   for(unsigned i=0;i<n_dof_perturbed;i++)
    {
     convergence_criterion +=
      (output[i] - calc_eigenvalue*input[i])
      *(output[i] - calc_eigenvalue*input[i]);
    }
   convergence_criterion = sqrt(convergence_criterion/n_dof_perturbed);

   // Print convergence criterion and threshold
   cout << "\nConvergence criterion     = " << convergence_criterion << endl;
   cout << "tolerance*abs(eigenvalue) = " << tolerance*abs(calc_eigenvalue)
        << endl;
   
   // Document convergence criterion and threshold in trace file
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

     // Output the dominant eigenvector (perturbed state solution)
     this->doc_solution(doc_info,false,true);

     // Return the number of iterations of the power method
     return p+1;
    }
   else
    {
     input = output;
     output.clear();
    }
  }

 // If we reach here, the power method has failed to converge in max_iter
 // iterations.
 cout << "\nPower method has failed to converge to within a tolerance of " 
      << tolerance << "\nin " << max_iter
      << " iterations of the power method" << endl;
 
 // Reset the solution counter to zero
 doc_info.number()=0;
 
 // Output the dominant eigenvector (perturbed state solution)
 this->doc_solution(doc_info,false,true);
 
 return max_iter;

} // End of perform_power_method



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//==start_of_main========================================================
/// Driver for Tuckerman counter-rotating lids problem (no power method)
//=======================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set duration of timestep
 const double dt = 0.005;

 // Set number of elements in radial (r) direction
 const unsigned base_n_r = 12;
 const unsigned perturbed_n_r = 12;
 
 // Set number of elements in axial (z) direction
 const unsigned base_n_z = 12;
 const unsigned perturbed_n_z = 12;
 
 // Determine height of computational domain (this is just the aspect
 // ratio since we take the cylinder radius to always be one)
 const double domain_height = GlobalPhysicalVariables::Gamma;
 
 // Set direction of gravity (vertically downwards)
 GlobalPhysicalVariables::G[0] = 0.0;
 GlobalPhysicalVariables::G[1] = -1.0;
 GlobalPhysicalVariables::G[2] = 0.0;
 
 // Set maximum number of iterations of the power method
 const unsigned max_iter = 10000;

 // Set tolerance for power method convergence (this is usually machine
 // precision according to Bai, Demmel, Dongarra, Ruhe & Van der Vorst)
 double power_method_tolerance = 1e-8;

 // Cache upper and lower limits and number of steps for loop over Re
 double lower = GlobalPhysicalVariables::Re_lower_limit;
 double upper = GlobalPhysicalVariables::Re_upper_limit;
 double n_param_steps = GlobalPhysicalVariables::Nparam_steps;

 // If we are doing a validation run, only compute one value of the
 // Reynolds number and use a lower tolerance
 if(CommandLineArgs::Argc>1)
  {
   n_param_steps = 1;
   power_method_tolerance = 1e-5;
  }

 // ---------------------------------------------------
 // RefineableLinearisedAxisymmetricQTaylorHoodElements
 // ---------------------------------------------------
 {
  cout << "\nDoing RefineableLinearisedAxisymmetricQTaylorHoodElements\n" 
       << endl;

  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<RefineableAxisymmetricQTaylorHoodElement,
   RefineableLinearisedAxisymmetricQTaylorHoodMultiDomainElement>
   problem(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,domain_height);
  
  // Set up labels for output
  DocInfo doc_info;
  
  // Output directory
  doc_info.set_directory("RESLT_TH");
  
  // Refine both base and perturbed state problems uniformly
  problem.base_state_problem_pt()->refine_uniformly();
  problem.perturbed_state_problem_pt()->refine_uniformly();

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
  sprintf(filename,"%s/global_trace_k%i.dat",
          doc_info.directory().c_str(),
          GlobalPhysicalVariables::k);
  global_trace.open(filename);
  global_trace << "Re, dominating eigenvalue, "
               << "n_base_flow_continuation_steps, "
               << "n_power_method_iterations" << std::endl;
  
  // Begin loop over Reynolds number
  for(unsigned param=0;param<n_param_steps;param++)
   {
    // Set Re_current to the correct value
    if(lower==upper || n_param_steps<=1)
     {
      GlobalPhysicalVariables::Re_current = lower;
     }
    else
     {
      GlobalPhysicalVariables::Re_current =
       lower + (((upper-lower)/(n_param_steps-1))*param);
     }
    
    // Set ReSt and ReInvFr to the correct values
    GlobalPhysicalVariables::ReSt_current = 
     GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::St;
    GlobalPhysicalVariables::ReInvFr_current = 
     GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::InvFr;
    
    cout << "\n====================================================" << endl;
    cout << "Beginning parameter run " << param+1 << " of "
         << n_param_steps << ": Re = "
         << GlobalPhysicalVariables::Re_current << endl;
    cout << "====================================================" << endl;
    
    // Pass updated Reynolds number (and ReSt and ReInvFr) to elements
    problem.pass_updated_nondim_parameters_to_elements();
    
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
    
    // Find steady base flow
    // ---------------------

    // Determine number of steps to take when ramping up to actual
    // Reynolds number
    const unsigned n_base_flow_continuation_steps = 1;
    
    // If we're on the first parameter, need to get to steady state by
    // continuation
    if(param==0)
     {
      cout << "\nEntering steady base flow loop..." << endl;
      
      // Store the current parameter value
      const double Re_current_backup = GlobalPhysicalVariables::Re_current;
      
      // Loop over continuation steps
      for(unsigned i=0;i<=n_base_flow_continuation_steps;i++)
       {
        // Update current parameters
        GlobalPhysicalVariables::Re_current =
         i*(Re_current_backup/n_base_flow_continuation_steps);
        GlobalPhysicalVariables::ReSt_current = 
         GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::St;
        GlobalPhysicalVariables::ReInvFr_current = 
         GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::InvFr;
        
        cout << "   - Setting GlobalPhysicalVariables::Re_current = "
             << GlobalPhysicalVariables::Re_current << endl;
        
        // Pass updated Reynolds number (and ReSt and ReInvFr) to elements
        problem.pass_updated_nondim_parameters_to_elements();
        
        // Perform a steady base flow solve
        problem.base_flow_steady_newton_solve();
       }
     }
    // Otherwise, just require one steady base flow solve
    else { problem.base_flow_steady_newton_solve(); }
    
    // Enable Jacobian reuse in the perturbed state problem. Note that this
    // can only be done since the problem is linear and the base state is
    // steady, and so the Jacobian will be the same at each timestep.
    // Calling this function also resets the
    // Problem::Jacobian_has_been_computed flag to "false", so that during
    // the first solve it will be recomputed. Note that this is the reason
    // that we call this function here rather than in the problem constructor.
    // Were we to call it in the constructor then the Jacobian would not be
    // recomputed when we changed the problem parameters (Reynolds number etc.)
    problem.perturbed_state_problem_pt()->enable_jacobian_reuse();
    
    // Perform power method and store number of iterations
    const unsigned n_power_method_iterations =
     problem.perform_power_method(dt,doc_info,power_method_tolerance,
                                  max_iter,eigenvalue,eigenvector);
    
    cout << "\nDominating eigenvalue is " << eigenvalue << endl;
    
    // The corresponding dominating eigenvector is "input"
    
    // Document in the global trace file
    global_trace << GlobalPhysicalVariables::Re_current << " ";
    ios_base::fmtflags flags =  cout.flags(); // Save old flags
    global_trace.precision(14); // ten decimal places
    global_trace << eigenvalue << " ";
    global_trace.flags(flags);  // Set the flags to the way they were
    global_trace << n_base_flow_continuation_steps << " "
                 << n_power_method_iterations << std::endl;
    
    // The initial guess for the next power method will be the eigenvector
    // calculated by this power method
    
    // Clear and close trace files
    problem.close_trace_files();

   } // End loop over Reynolds number
 }

 // --------------------------------------------------------
 // RefineableLinearisedAxisymmetricQCrouzeixRaviartElements
 // --------------------------------------------------------
{
  cout << "\nDoing RefineableLinearisedAxisymmetricQCrouzeixRaviartElements\n" 
       << endl;

  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<RefineableAxisymmetricQCrouzeixRaviartElement,
   RefineableLinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement>
   problem(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,domain_height);
  
  // Set up labels for output
  DocInfo doc_info;
  
  // Output directory
  doc_info.set_directory("RESLT_CR");
  
  // Refine both base and perturbed state problems uniformly
  problem.base_state_problem_pt()->refine_uniformly();
  problem.perturbed_state_problem_pt()->refine_uniformly();

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
  sprintf(filename,"%s/global_trace_k%i.dat",
          doc_info.directory().c_str(),
          GlobalPhysicalVariables::k);
  global_trace.open(filename);
  global_trace << "Re, dominating eigenvalue, "
               << "n_base_flow_continuation_steps, "
               << "n_power_method_iterations" << std::endl;
  
  // Begin loop over Reynolds number
  for(unsigned param=0;param<n_param_steps;param++)
   {
    // Set Re_current to the correct value
    if(lower==upper || n_param_steps<=1)
     {
      GlobalPhysicalVariables::Re_current = lower;
     }
    else
     {
      GlobalPhysicalVariables::Re_current =
       lower + (((upper-lower)/(n_param_steps-1))*param);
     }
    
    // Set ReSt and ReInvFr to the correct values
    GlobalPhysicalVariables::ReSt_current = 
     GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::St;
    GlobalPhysicalVariables::ReInvFr_current = 
     GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::InvFr;
    
    cout << "\n====================================================" << endl;
    cout << "Beginning parameter run " << param+1 << " of "
         << n_param_steps << ": Re = "
         << GlobalPhysicalVariables::Re_current << endl;
    cout << "====================================================" << endl;
    
    // Pass updated Reynolds number (and ReSt and ReInvFr) to elements
    problem.pass_updated_nondim_parameters_to_elements();
    
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
    
    // Find steady base flow
    // ---------------------

    // Determine number of steps to take when ramping up to actual
    // Reynolds number
    const unsigned n_base_flow_continuation_steps = 1;
    
    // If we're on the first parameter, need to get to steady state by
    // continuation
    if(param==0)
     {
      cout << "\nEntering steady base flow loop..." << endl;
      
      // Store the current parameter value
      const double Re_current_backup = GlobalPhysicalVariables::Re_current;
      
      // Loop over continuation steps
      for(unsigned i=0;i<=n_base_flow_continuation_steps;i++)
       {
        // Update current parameters
        GlobalPhysicalVariables::Re_current =
         i*(Re_current_backup/n_base_flow_continuation_steps);
        GlobalPhysicalVariables::ReSt_current = 
         GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::St;
        GlobalPhysicalVariables::ReInvFr_current = 
         GlobalPhysicalVariables::Re_current*GlobalPhysicalVariables::InvFr;
        
        cout << "   - Setting GlobalPhysicalVariables::Re_current = "
             << GlobalPhysicalVariables::Re_current << endl;
        
        // Pass updated Reynolds number (and ReSt and ReInvFr) to elements
        problem.pass_updated_nondim_parameters_to_elements();
        
        // Perform a steady base flow solve
        problem.base_flow_steady_newton_solve();
       }
     }
    // Otherwise, just require one steady base flow solve
    else { problem.base_flow_steady_newton_solve(); }
    
    // Enable Jacobian reuse in the perturbed state problem. Note that this
    // can only be done since the problem is linear and the base state is
    // steady, and so the Jacobian will be the same at each timestep.
    // Calling this function also resets the
    // Problem::Jacobian_has_been_computed flag to "false", so that during
    // the first solve it will be recomputed. Note that this is the reason
    // that we call this function here rather than in the problem constructor.
    // Were we to call it in the constructor then the Jacobian would not be
    // recomputed when we changed the problem parameters (Reynolds number etc.)
    problem.perturbed_state_problem_pt()->enable_jacobian_reuse();
    
    // Perform power method and store number of iterations
    const unsigned n_power_method_iterations =
     problem.perform_power_method(dt,doc_info,power_method_tolerance,
                                  max_iter,eigenvalue,eigenvector);
    
    cout << "\nDominating eigenvalue is " << eigenvalue << endl;
    
    // The corresponding dominating eigenvector is "input"
    
    // Document in the global trace file
    global_trace << GlobalPhysicalVariables::Re_current << " ";
    ios_base::fmtflags flags =  cout.flags(); // Save old flags
    global_trace.precision(14); // ten decimal places
    global_trace << eigenvalue << " ";
    global_trace.flags(flags);  // Set the flags to the way they were
    global_trace << n_base_flow_continuation_steps << " "
                 << n_power_method_iterations << std::endl;
    
    // The initial guess for the next power method will be the eigenvector
    // calculated by this power method
    
    // Clear and close trace files
    problem.close_trace_files();

   } // End loop over Reynolds number
 }

} // End of main

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
// Driver for validation problem for linearised axisymmetric Navier-Stokes
// equations in cylindrical polar coordinates which replicates the
// investigation of C. Nore, M. Tartar, O. Daube and L.S. Tuckerman's
// paper "Survey of instability thresholds of flow between exactly
// counter-rotating disks" (2003).

// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Linearised axisymmetric Navier-Stokes elements
#include "linearised_axisym_navier_stokes_elements.h"
#include "linearised_axisym_navier_stokes_elements.cc"
#include "multi_domain_linearised_axisym_navier_stokes_elements.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"


using namespace std;

using namespace oomph;

//using namespace MathematicalConstants;



//==start_of_namespace_for_physical_parameters===========================
/// Namespace for physical parameters
//=======================================================================
namespace GlobalPhysicalVariables
{
 /// Reynolds number
 double Re = 301.0;

 /// Product of Reynolds and Strouhal numbers
 double ReSt = 301.0; // (St = 1.0)

 /// Product of Rynolds number and inverse of Froude number
 double ReInvFr = 0.0; // (Fr = inf)

 /// Aspect ratio (cylinder height / cylinder radius)
 double Gamma = 1.0;

 /// Azimuthal mode number k in e^ik(theta) decomposition
 int k = 2;

 /// Direction of gravity
 Vector<double> G(3);

} // End of namespace



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



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

 /// Actions before the timestep (empty)
 void actions_before_implicit_timestep() {}

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// Set initial condition (incl previous timesteps)
 void set_initial_condition();

 /// Set the boundary conditions
 void set_boundary_conditions();

 /// Access function for the specific mesh
 RectangularQuadMesh<BASE_ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RectangularQuadMesh<BASE_ELEMENT>*>
    (Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   dynamic_cast<BASE_ELEMENT*>(mesh_pt()->element_pt(e))
    ->fix_pressure(pdof,pvalue);
  }

}; // End of base_state_problem class



//========start_of_base_state_constructor================================
/// Constructor for base state Tuckerman counter-rotating lids problem
//=======================================================================
template<class BASE_ELEMENT>
BaseStateProblem<BASE_ELEMENT>::
BaseStateProblem(const unsigned& n_r,
                 const unsigned& n_z,
                 const double& domain_height)
{
 // Build and assign mesh
 Problem::mesh_pt() = 
  new RectangularQuadMesh<BASE_ELEMENT>(n_r,n_z,0.0,1.0,0.0,domain_height);
 
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
void BaseStateProblem<BASE_ELEMENT>::
set_initial_condition()
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
set_boundary_conditions()
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
     // Set up pointer to node
     Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);

     // Get the radial value
     double r = mesh_pt()->boundary_node_pt(b,n)->x(0);

     // Provide storage for azimuthal velocity value
     double azi_vel = 0.0;
     
     // Apply "smoothing" if r > a
     if(r>a) { azi_vel = a*(1.0-r)/(1.0-a); }
     else { azi_vel = r; }

     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       switch(b)
        {
         // Outer wall (all components = 0)
        case 1:
         nod_pt->set_value(0,i,0.0);
         break;

         // Top lid (azimuthal component = r, others = 0)
        case 2:
         if(i==2) { nod_pt->set_value(0,i,azi_vel); }
         else { nod_pt->set_value(0,i,0.0); }
         break;

         // Bottom lid (azimuthal component = -r, others = 0)
        case 0:
         if(i==2) { nod_pt->set_value(0,i,-azi_vel); }
         else { nod_pt->set_value(0,i,0.0); }
         break;

         // Symmetry boundary (axial component not set, others = 0)
        case 3:
         if(i!=1) { nod_pt->set_value(0,i,0.0); }
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
void BaseStateProblem<BASE_ELEMENT>::
doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[256];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;
 
 // Open solution output file
 sprintf(filename,"%s/base_soln%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 // Output solution to file
 mesh_pt()->output(some_file,npts);
 
 // Close solution output file
 some_file.close();

} // End of doc_solution for base state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_perturbed_state_problem_class===============================
/// Perturbed state problem class for Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
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

 /// Set the three velocity components to a smoothly varying,
 /// non-zero profile
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
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

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

 /// Height of domain
 double Domain_height;

}; // End of perturbed_state_problem class



//==start_of_perturbed_state_constructor=================================
/// Constructor for perturbed state Tuckerman counter-rotating lids
/// problem
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
PerturbedStateProblem
<BASE_ELEMENT,PERTURBED_ELEMENT>::
PerturbedStateProblem(const unsigned& n_r,
                      const unsigned& n_z,
                      const double& domain_height,
                      Mesh* external_mesh_pt)
 : Base_state_mesh_pt(external_mesh_pt),
   Domain_height(domain_height)
{
 // Be less verbose during newton solve
 this->disable_info_in_newton_solve();

 // Be less verbose about linear solve timings
 linear_solver_pt()->disable_doc_time();
 
 // Tell problem that it is linear (avoid doing unnecessary checks)
 Problem::Problem_is_nonlinear = false;

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new BDF<1>);
 
 // Build and assign mesh
 Problem::mesh_pt() = 
  new RectangularQuadMesh<PERTURBED_ELEMENT>(n_r,n_z,0.0,1.0,0.0,domain_height,
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
     // Set up pointer to node
     Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);

     // On all solid boundaries, pin all velocity components
     if(b!=3)
      {
       nod_pt->pin(0); // Radial (real)
       nod_pt->pin(1); // Radial (imaginary)
       nod_pt->pin(2); // Axial (real)
       nod_pt->pin(3); // Axial (imaginary)
       nod_pt->pin(4); // Azimuthal (real)
       nod_pt->pin(5); // Azimuthal (imaginary)
      }
     // On symmetry boundary, pin only radial and azimuthal velocity
     // components, unless we are looking at the k=1 azimuthal mode.
     // In this case, do not pin any velocity components on this boundary
     else
      {
       if(GlobalPhysicalVariables::k!=1)
        {
         nod_pt->pin(0); // Radial (real)
         nod_pt->pin(1); // Radial (imaginary)
         nod_pt->pin(4); // Azimuthal (real)
         nod_pt->pin(5); // Azimuthal (imaginary)
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

   // Set the Reynolds number
   el_pt->re_pt() = &GlobalPhysicalVariables::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &GlobalPhysicalVariables::ReSt;
   
   // The mesh remains fixed
   el_pt->disable_ALE();

   // Set the azimuthal wavenumber
   el_pt->azimuthal_mode_number_pt() = &GlobalPhysicalVariables::k;

  } // End of loop over elements
 
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
/// Set the three velocity components to a smoothly-varying,
/// non-zero profile
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
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
 // using a self-starting timestepper (BDF<1>)
 
 // Set the boundary conditions
 set_boundary_conditions();

} // End of set_initial_condition for perturbed state



//======start_of_perturbed_state_set_boundary_conditions=================
/// Set the (homogeneous) boundary conditions for the perturbed state
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
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
     // Set up pointer to node
     Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);

     // For the solid boundaries set all components to zero
     if(b!=3)
      {
       nod_pt->set_value(0,0,0.0); // Radial
       nod_pt->set_value(0,1,0.0); // Radial
       nod_pt->set_value(0,2,0.0); // Axial
       nod_pt->set_value(0,3,0.0); // Axial
       nod_pt->set_value(0,4,0.0); // Azimuthal
       nod_pt->set_value(0,5,0.0); // Azimuthal
      }
     // On symmetry boundary, set only radial and azimuthal components to
     // zero, unless we are looking at the k=1 azimuthal mode. In this
     // case, do not set the value of any of the velocity components
     else
      {
       if(GlobalPhysicalVariables::k!=1)
        {
         nod_pt->set_value(0,0,0.0); // Radial
         nod_pt->set_value(0,1,0.0); // Radial
         nod_pt->set_value(0,4,0.0); // Azimuthal
         nod_pt->set_value(0,5,0.0); // Azimuthal
        }
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
} // End of set_boundary_conditions for perturbed state



//==start_of_perturbed_state_doc_solution================================
/// Document the perturbed state solution
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
void PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[256];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;
 
 // Open solution output file
 sprintf(filename,"%s/perturbed_soln%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 // Output solution to file
 mesh_pt()->output(some_file,npts);
 
 // Close solution output file
 some_file.close();
 
} // End of doc_solution for perturbed state



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=====start_of_stability_problem_class==================================
/// Container class for the Tuckerman counter-rotating lids
/// stability problem
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
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

 /// Set initial conditions
 void set_initial_condition()
  {
   Base_state_problem_pt->set_initial_condition();
   Perturbed_state_problem_pt->set_initial_condition();
  }

 /// Integrate forwards in time with timestep dt for n_timesteps
 void unsteady_run(const double& dt,
                   const unsigned& n_timesteps,
                   DocInfo& doc_info);
 
 /// Perform a power method with a maximum of max_iter iterations. This
 /// function requires an initial guess for the dominant eigenvector
 /// "input"and a tolerance which controls at which point the dominant
 /// eigenvalue will be deemed converged. At this point the calculated
 /// value of this eigenvalue "calc_eigenvalue" and its corresponding
 /// eigenvector "input" will be returned. The return value of this
 /// function is the number of power method iterations which were
 /// performed.
 unsigned perform_power_method(const double& dt,
                               const unsigned& n_timesteps,
                               DocInfo& doc_info,
                               const double& tolerance,
                               const unsigned& max_iter,
                               double& calc_eigenvalue,
                               DoubleVector& input);

 /// Access function for base state problem
 BaseStateProblem<BASE_ELEMENT>*
 base_state_problem_pt() const { return Base_state_problem_pt; }
 
 /// Access function for perturbed state problem
 PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>*
 perturbed_state_problem_pt() const { return Perturbed_state_problem_pt; }
 
private:

 /// Pointer to base state problem class
 BaseStateProblem<BASE_ELEMENT>* Base_state_problem_pt;
 
 /// Pointer to perturbed state problem class
 PerturbedStateProblem<BASE_ELEMENT,PERTURBED_ELEMENT>*
 Perturbed_state_problem_pt;

}; // End of stability_problem class



//==start_of_unsteady_run================================================
/// Integrate forwards in time with timestep dt for n_timesteps
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
void StabilityProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
unsteady_run(const double& dt,
             const unsigned& n_timesteps,
             DocInfo& doc_info)
{
 // Timestepping loop
 for(unsigned i=0;i<n_timesteps;i++)
  {
   // Output timestep and global time
   std::cout << "  Timestep " << i+1 << " of " << n_timesteps
             << "; Time is " << Perturbed_state_problem_pt->time_pt()->time()
             << std::endl;

   // Take timestep
   Perturbed_state_problem_pt->unsteady_newton_solve(dt);

  } // End of loop over timesteps
} // End of unsteady_run



//==start_of_perform_power_method========================================
/// Perform a power method with a maximum of max_iter iterations. This
/// function requires an initial guess for the dominant eigenvector
/// "input" and a tolerance which controls at which point the dominant
/// eigenvalue will be deemed converged. At this point the calculated
/// value of this eigenvalue "calc_eigenvalue" and its corresponding
/// eigenvector "input" will be returned. The return value of this
/// function is the number of power method iterations which were
/// performed.
//=======================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
unsigned StabilityProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
perform_power_method(const double& dt,
                     const unsigned& n_timesteps,
                     DocInfo& doc_info,
                     const double& tolerance,
                     const unsigned& max_iter,
                     double& calc_eigenvalue,
                     DoubleVector& input)
{
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

   // Perform unsteady run for n_timesteps
   this->unsteady_run(dt,n_timesteps,doc_info);

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
   cout << "\nConvergence criterion     = " << convergence_criterion
        << "\ntolerance*abs(eigenvalue) = " << tolerance*abs(calc_eigenvalue)
        << endl;

   // Check for convergence
   if(convergence_criterion <= tolerance*abs(calc_eigenvalue))
    {
     cout << "\nPower method has converged to within a tolerance of " 
          << tolerance << "\nin " << p+1
          << " iterations of the power method.\n" << endl;

     // Output the dominant eigenvector
     Perturbed_state_problem_pt->doc_solution(doc_info);

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
 
 // Output the dominant eigenvector
 Perturbed_state_problem_pt->doc_solution(doc_info);
   
 return max_iter;

} // End of perform_power_method



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//==start_of_main========================================================
/// Driver for Tuckerman counter-rotating lids problem
//=======================================================================
int main(int argc, char* argv[])
{
 // Number of elements in radial (r) direction
 const unsigned base_n_r = 36;
 const unsigned perturbed_n_r = 16;

 // Number of elements in axial (z) direction
 const unsigned base_n_z = 36;
 const unsigned perturbed_n_z = 16;
 
 // Determine height of computational domain (this is just the aspect
 // ratio since we take the cylinder radius to always be one)
 const double domain_height = GlobalPhysicalVariables::Gamma;

 // Set direction of gravity (vertically downwards)
 GlobalPhysicalVariables::G[0] = 0.0;
 GlobalPhysicalVariables::G[1] = -1.0;
 GlobalPhysicalVariables::G[2] = 0.0;
  
 // Set up labels for output
 DocInfo doc_info;

 // Determine length of timestep
 const double dt = 0.01;

 // Set number of timesteps to perform during a single iteration of the
 // power method
 const unsigned nstep = 1;

 // Cache tolerance for power method convergence
 // Note that we are using a very large tolerance for the purposes of
 // this demo driver so that it converges in a very small number of steps
 const double power_method_tolerance = 1e-4;

 // Set maximum number of iterations of the power method
 const unsigned max_iter = 1000;

 // For the purposes of code validation we will solve this problem twice.
 // The first time we shall use Taylor-Hood elements and the second time
 // Crouzeix-Raviart elements.

 // --------------------------------------------------------------------
 // Taylor-Hood elements
 // --------------------------------------------------------------------
 {
  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<AxisymmetricQTaylorHoodElement,
   LinearisedAxisymmetricQTaylorHoodMultiDomainElement>
   problem_TH(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,domain_height);
  
  // Output directory
  doc_info.set_directory("RESLT_TH");

  // Create and initialise trace file
  ofstream trace;
  char filename[256];
  sprintf(filename,"%s/trace.dat",
          doc_info.directory().c_str());
  trace.open(filename);
  trace << "Re, dominating eigenvalue, n_power_method_iterations" << std::endl;

  // Initialise timestep for perturbed state problem
  // (also sets weights for all timesteppers)
  problem_TH.perturbed_state_problem_pt()->initialise_dt(dt);
  
  // Set initial conditions
  problem_TH.set_initial_condition();
  
  // Set up storage for dominant eigenvector and eigenvalue
  DoubleVector eigenvector;
  double eigenvalue = 0.0;
  
  // Get initial guess for eigenvector (these are the initial conditions
  // set up in PerturbedStateProblem::set_initial_conditions())
  problem_TH.perturbed_state_problem_pt()->get_dofs(eigenvector);
  
  // Initialise counter for solutions
  doc_info.number()=0;
  
  // Set base state problem's boundary conditions
  problem_TH.base_state_problem_pt()->set_boundary_conditions();
  
  // Perform a steady base flow solve
  problem_TH.base_state_problem_pt()->steady_newton_solve();
  
  // Doc the base flow solution
  problem_TH.base_state_problem_pt()->doc_solution(doc_info);
  
  // Enable Jacobian reuse in the perturbed state problem. Note that this
  // can only be done since the problem is linear and the base state is
  // steady, and so the Jacobian will be the same at each timestep.
  problem_TH.perturbed_state_problem_pt()->enable_jacobian_reuse();
  
  // Perform power method and store number of iterations
  const unsigned n_power_method_iterations =
   problem_TH.perform_power_method(dt,nstep,doc_info,power_method_tolerance,
                                   max_iter,eigenvalue,eigenvector);
  
  cout << "\nDominating eigenvalue is " << eigenvalue << endl;
  
  // Document in trace file
  trace << GlobalPhysicalVariables::Re << " "
        << eigenvalue << " "
        << n_power_method_iterations << std::endl;
 }
 
 
 // --------------------------------------------------------------------
 // Crouzeix-Raviart elements
 // --------------------------------------------------------------------
 {
  // Build stability problem (this creates base and perturbed state problems)
  StabilityProblem<AxisymmetricQCrouzeixRaviartElement,
   LinearisedAxisymmetricQCrouzeixRaviartMultiDomainElement>
   problem_CR(base_n_r,base_n_z,perturbed_n_r,perturbed_n_z,domain_height);
  
  // Output directory
  doc_info.set_directory("RESLT_CR");
  
  // Create and initialise trace file
  ofstream trace;
  char filename[256];
  sprintf(filename,"%s/trace.dat",
          doc_info.directory().c_str());
  trace.open(filename);
  trace << "Re, dominating eigenvalue, n_power_method_iterations" << std::endl;

  // Initialise timestep for perturbed state problem
  // (also sets weights for all timesteppers)
  problem_CR.perturbed_state_problem_pt()->initialise_dt(dt);
  
  // Set initial conditions
  problem_CR.set_initial_condition();
  
  // Set up storage for dominant eigenvector and eigenvalue
  DoubleVector eigenvector;
  double eigenvalue = 0.0;
  
  // Get initial guess for eigenvector (these are the initial conditions
  // set up in PerturbedStateProblem::set_initial_conditions())
  problem_CR.perturbed_state_problem_pt()->get_dofs(eigenvector);
  
  // Initialise counter for solutions
  doc_info.number()=0;
  
  // Set base state problem's boundary conditions
  problem_CR.base_state_problem_pt()->set_boundary_conditions();
  
  // Perform a steady base flow solve
  problem_CR.base_state_problem_pt()->steady_newton_solve();
  
  // Doc the base flow solution
  problem_CR.base_state_problem_pt()->doc_solution(doc_info);
  
  // Enable Jacobian reuse in the perturbed state problem. Note that this
  // can only be done since the problem is linear and the base state is
  // steady, and so the Jacobian will be the same at each timestep.
  problem_CR.perturbed_state_problem_pt()->enable_jacobian_reuse();
  
  // Perform power method and store number of iterations
  const unsigned n_power_method_iterations =
   problem_CR.perform_power_method(dt,nstep,doc_info,power_method_tolerance,
                                   max_iter,eigenvalue,eigenvector);
  
  cout << "\nDominating eigenvalue is " << eigenvalue << endl;
  
  // Document in trace file
  trace << GlobalPhysicalVariables::Re << " "
        << eigenvalue << " "
        << n_power_method_iterations << std::endl;
 }

} // End of main

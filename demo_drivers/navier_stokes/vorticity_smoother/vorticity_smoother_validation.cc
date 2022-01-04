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
// Driver for Anne's MSc project

// The oomphlib headers
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;

//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace GlobalParameters
{
 /// Reynolds number
 double Re=10.0;

 /// Length of computational domain
 double Length=2.0*MathematicalConstants::Pi;

 /// Height of computational domain
 double Height=2.0*MathematicalConstants::Pi;

 //----------------------
 // Parameters for vortex
 //----------------------
 /// The value of the kinematic viscosity (assumed to be the same everywhere)
 double Kinematic_viscosity=1.0;
 
 /// Returns the velocity field associated with the Taylor-Green vortex
 /// solution (for validation of vorticity projection)
 void sin_cos_velocity_field(const Vector<double>& x,
			     Vector<double>& u)
 {
  // Time remains fixed so set to zero
  double t=0.0;
  
  // Calculate the temporal component
  double temporal_component=exp(-2.0*Kinematic_viscosity*t);
  
  // Frequency in the x-direction
  double omega_x=2.0*MathematicalConstants::Pi/Length;
  
  // Frequency in the y-direction
  double omega_y=2.0*MathematicalConstants::Pi/Height;

  // The phase shift
  double phi=MathematicalConstants::Pi/2.0;
  
  // First velocity component
  u[0]=cos(omega_x*x[0]+phi)*sin(omega_y*x[1]+phi)*temporal_component;
  
  // Second velocity component
  u[1]=-sin(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*temporal_component;
 } // End of sin_cos_velocity_field


 /// Returns the vorticity field associated with the Taylor-Green vortex
 /// solution (for validation of vorticity projection)
 void sin_cos_vorticity(const Vector<double>& x, 
                        Vector<Vector<double> >& vort_and_derivs)
 {
  // Time remains fixed so set to zero
  double t=0.0;
  
  // Calculate the temporal component
  double temporal_component=exp(-2.0*Kinematic_viscosity*t);
  
  // Frequency in the x-direction
  double omega_x=2.0*MathematicalConstants::Pi/Length;
  
  // Frequency in the y-direction
  double omega_y=2.0*MathematicalConstants::Pi/Height;

  // The phase shift
  double phi=MathematicalConstants::Pi/2.0;

  // The vorticity itself: w
  vort_and_derivs[0][0]=(-1.0*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The first derivative of the vorticity: dw/dx
  vort_and_derivs[1][0]=(omega_x*(omega_x+omega_y)*
			 sin(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);

  // The first derivative of the vorticity: dw/dy
  vort_and_derivs[1][1]=(omega_y*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*sin(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The second derivative of the vorticity: dw^2/dx^2
  vort_and_derivs[2][0]=(omega_x*omega_x*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);

  // The second derivative of the vorticity: dw^2/dxdy
  vort_and_derivs[2][1]=(-1.0*omega_x*omega_y*(omega_x+omega_y)*
			 sin(omega_x*x[0]+phi)*sin(omega_y*x[1]+phi)*
			 temporal_component);

  // The second derivative of the vorticity: dw^2/dy^2
  vort_and_derivs[2][2]=(omega_y*omega_y*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);

  // The third derivative of the vorticity: dw^3/dx^3
  vort_and_derivs[3][0]=(-1.0*omega_x*omega_x*omega_x*(omega_x+omega_y)*
			 sin(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The third derivative of the vorticity: dw^3/dx^2dy
  vort_and_derivs[3][1]=(-1.0*omega_x*omega_x*omega_y*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*sin(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The third derivative of the vorticity: dw^3/dxdy^2
  vort_and_derivs[3][2]=(-1.0*omega_x*omega_y*omega_y*(omega_x+omega_y)*
			 sin(omega_x*x[0]+phi)*cos(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The third derivative of the vorticity: dw^3/dy^3
  vort_and_derivs[3][3]=(-1.0*omega_y*omega_y*omega_y*(omega_x+omega_y)*
			 cos(omega_x*x[0]+phi)*sin(omega_y*x[1]+phi)*
			 temporal_component);
  
  // The derivatives of the velocity: du/dx
  vort_and_derivs[4][0]=(-1.0*omega_x*sin(omega_x*x[0]+phi)*
			 sin(omega_y*x[1]+phi)*temporal_component);
  
  // The derivatives of the velocity: du/dy
  vort_and_derivs[4][1]=(omega_y*cos(omega_x*x[0]+phi)*
			 cos(omega_y*x[1]+phi)*temporal_component);

  // The derivatives of the velocity: dv/dx
  vort_and_derivs[4][2]=(-1.0*omega_x*cos(omega_x*x[0]+phi)*
			 cos(omega_y*x[1]+phi)*temporal_component);
  
  // The derivatives of the velocity: dv/dy
  vort_and_derivs[4][3]=(omega_y*sin(omega_x*x[0]+phi)*
			 sin(omega_y*x[1]+phi)*temporal_component);
 } // End of sin_cos_vorticity

 
 /// Synthetic velocity field for validation
 void synthetic_velocity_field(const Vector<double>& x, 
                               Vector<double>& veloc)
 {
  // Get the sin/cos velocity field
  sin_cos_velocity_field(x,veloc);
 } // End of synthetic_velocity_field

 
 /// Synthetic vorticity field and derivs for validation
 void synthetic_vorticity(const Vector<double>& x,
                          Vector<Vector<double> >& vort_and_derivs)
 {
  // Get the sin/cos vorticity field
  sin_cos_vorticity(x,vort_and_derivs);
 } // End of synthetic_vorticity

 
 /// Initial condition for velocity
 void initial_condition(const Vector<double>& x,
			Vector<double>& u)
 {
  // Call a helper function to calculate the initial conditions
  synthetic_velocity_field(x,u);
 } // End of initial_condition
} // End of GlobalParameters


//===start_of_problem_class=============================================
/// Problem class for Anne's MSc problem
//======================================================================
template<class ELEMENT>
class VorticityRecoveryProblem : public Problem
{
public:

 /// Constructor:
 VorticityRecoveryProblem();

 // Update before solve is empty
 void actions_before_newton_solve(){}

 /// Update after solve is empty
 void actions_after_newton_solve(){}
 
 /// Actions before adapt: empty
 void actions_before_adapt(){} 

 /// After adaptation
 void actions_after_adapt()
 {
  // Apply the appropriate boundary conditions
  apply_boundary_conditions();
  
  // Do everything that needs to be done to complete the problem setup
  complete_problem_setup();
 } // End of actions_after_adapt
   
 /// Apply Dirichlet conditions on all of the boundaries
 void apply_boundary_conditions();

 /// Complete problem setup
 void complete_problem_setup();
 
 /// Assign the synthetic flow field
 void assign_synthetic_veloc_field();

 /// Check the vorticity smoothing
 void check_smoothed_vorticity(DocInfo& doc_info);

 /// Document the solution
 void doc_solution(DocInfo& doc_info,const bool& vorticity_recovered=false);

private:

 /// Vorticity recoverer
 VorticitySmoother<ELEMENT>* Vorticity_recoverer_pt;
}; // end of problem_class


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT>
VorticityRecoveryProblem<ELEMENT>::VorticityRecoveryProblem()
{
 // Set the order of vorticity recovery
 unsigned nrecovery_order=2;
 
 // Make an instance of the vorticity recoverer
 Vorticity_recoverer_pt=new VorticitySmoother<ELEMENT>(nrecovery_order);

 // Allocate the timestepper
 add_time_stepper_pt(new BDF<2>); 

 // Number of elements in x direction.
 // purely nominal initial value
 unsigned n_x=4;

 // Number of elements in y direction
 // purely nominal initial value
 unsigned n_y=4; 

 // Left end of computational domain
 double x_min=0.0;
 
 // Right end of computational domain
 double x_max=GlobalParameters::Length;
 
 // Botton of computational domain
 double y_min=0.0;

 // Height of computational domain
 double y_max=GlobalParameters::Height;
 
 // Set the maximum order of vorticity derivatives to calculate to
 // zero -- only calculate the vorticity itself!
 VorticityRecoveryHelpers::Recovery_helper.
  set_maximum_order_of_vorticity_derivative(3);
 
 // Set the maximum order of velocity derivatives to calculate to
 // zero -- don't enable the recovery of the velocity derivatives!
 VorticityRecoveryHelpers::Recovery_helper.
  set_maximum_order_of_velocity_derivative(1);
  
 // Now create the mesh 
 mesh_pt()=new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,
						      x_min,x_max,
						      y_min,y_max,
						      time_stepper_pt());
 
 // Apply the appropriate boundary conditions
 apply_boundary_conditions();
 
 // Complete the problem setup to make the elements fully functional
 complete_problem_setup();
 
 // Assign equation numbers
 assign_eqn_numbers();
 
 // Assign the initial condition
 assign_synthetic_veloc_field();
} // End of VorticityRecoveryProblem


//======start_of_assign_synthetic_veloc_field=============================
/// Assign the synthetic flow field
//========================================================================
template<class ELEMENT>
void VorticityRecoveryProblem<ELEMENT>::assign_synthetic_veloc_field()
{ 
 // The number of dimensions in the mesh
 unsigned n_dim=2;
   
 // Storage for the (Eulerian) coordinates of the node
 Vector<double> x(n_dim,0.0);
   
 // Storage for the velocity components
 Vector<double> u(n_dim,0.0);
 
 // Get the number of nodes in the mesh
 unsigned num_nod=mesh_pt()->nnode();
 
 // Loop over nodes
 for (unsigned n=0;n<num_nod;n++)
 {
  // Get a pointer to the n-th node in the mesh
  Node* node_pt=mesh_pt()->node_pt(n);
  
  // Get the nodal coordinates of the n-th node in the mesh
  node_pt->position(x);
   
  // Get initial velocity field
  GlobalParameters::synthetic_velocity_field(x,u);
   
  // Loop over the velocity components
  for (unsigned i=0;i<n_dim;i++)
  {    
   // Set the i-th velocity component value
   node_pt->set_value(i,u[i]);
  }
 } // for (unsigned n=0;n<num_nod;n++)
} // End of assign_synthetic_veloc_field

//======start_of_apply_boundary_conditions================================
/// Apply the appropriate boundary conditions
//========================================================================
template<class ELEMENT>
void VorticityRecoveryProblem<ELEMENT>::apply_boundary_conditions()
{ 
 // The number of dimensions in the mesh
 unsigned n_dim=2;
   
 // Storage for the (Eulerian) coordinates of the node
 Vector<double> x(n_dim,0.0);
   
 // Storage for the velocity components
 Vector<double> u(n_dim,0.0);
   
 // Get the number of boundaries in the mesh
 unsigned num_bound=mesh_pt()->nboundary();

 // Loop over the boundaries
 for (unsigned ibound=0;ibound<num_bound;ibound++)
 {
  // Get the number of nodes on the ibound-th boundary
  unsigned num_nod=mesh_pt()->nboundary_node(ibound);

  // Loop over the nodes on the ibound-th boundary
  for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get a pointer to the inod-th node on the ibound-th boundary
   Node* node_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   
   // Get the global coordinates of this node
   node_pt->position(x);
   
   // Get the exact velocity field values at this node
   GlobalParameters::synthetic_velocity_field(x,u);

   // Loop over the velocity components
   for (unsigned i=0;i<n_dim;i++)
   {    
    // Pin the i-th velocity component
    node_pt->pin(i);

    // Set the i-th velocity component value
    node_pt->set_value(i,u[i]);
   } 
  } // for (unsigned inod=0;inod<num_nod;inod++)
 } // for (unsigned ibound=0;ibound<num_bound;ibound++)
} // End of apply_boundary_conditions


//========================================================================
/// Complete problem setup
//========================================================================
template<class ELEMENT>
void VorticityRecoveryProblem<ELEMENT>::complete_problem_setup()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(mesh_pt()->element_pt());

 // Get the number of elements in the mesh
 unsigned n_el=mesh_pt()->nelement();
 
 // Loop over the elements
 for (unsigned e=0;e<n_el;e++)
 {
  // Upcast to a fluid element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

  // Set the Reynolds number
  el_pt->re_pt()=&GlobalParameters::Re;

  // Set the Womersley number, i.e. Reynolds*Strouhal
  el_pt->re_st_pt()=&GlobalParameters::Re;

  // Set exact solution for vorticity and derivatives (for validation)
  el_pt->exact_vorticity_fct_pt()=&GlobalParameters::synthetic_vorticity;

  // Pin the smoothed vorticity
  el_pt->pin_smoothed_vorticity();
 } 

 // Upcast the first element to a fluid element
 ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));

 // Calculate the kinematic viscosity (constant everywhere)
 double kinematic_viscosity=(el_pt->viscosity_ratio())/(el_pt->density_ratio());

 // Assign it!
 GlobalParameters::Kinematic_viscosity=kinematic_viscosity;
 
 // Pin the redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
} // End of complete_problem_setup


//==start_check_smoothed_vort=============================================
/// Check the smoothed vorticity
//========================================================================
template<class ELEMENT>
void VorticityRecoveryProblem<ELEMENT>::
check_smoothed_vorticity(DocInfo& doc_info)
{ 
 // Get any element in the mesh
 ELEMENT* const el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));

 ofstream some_file;
 char filename[10000];
 sprintf(filename,"%s/vorticity_convergence.dat",
         doc_info.directory().c_str());
 some_file.open(filename);
 some_file << "VARIABLES=\"nel\",\"sqrt(1/nel)\",";
 if (el_pt->nvorticity_derivatives_to_recover()>=1)
 {
  some_file << "\"Error(vort)\",";
 }
 if (el_pt->nvorticity_derivatives_to_recover()>=3)
 {
  some_file << "\"Error(dvort/dx)\","
	    << "\"Error(dvort/dy)\",";
 }
 if (el_pt->nvorticity_derivatives_to_recover()>=6)
 {
  some_file << "\"Error(d^2vort/dx^2)\","
	    << "\"Error(d^2vort/dxdy)\","
	    << "\"Error(d^2vort/dy^2)\",";
 }
 if (el_pt->nvorticity_derivatives_to_recover()>=10)
 {
  some_file << "\"Error(d^3vort/dx^3)\","
	    << "\"Error(d^3vort/dx^2dy)\","
	    << "\"Error(d^3vort/dxdy^2)\","
	    << "\"Error(d^3vort/dy^3)\",";
 }
 if (el_pt->nvelocity_derivatives_to_recover()==4)
 {
  some_file << "\"Error(du/dx)\","
	    << "\"Error(du/dy)\","
	    << "\"Error(dv/dx)\","
	    << "\"Error(dv/dy)\",";
 }
 some_file << "\"Area\"\n";
 
 // Number of derivatives to be recovered
 unsigned n_recovered_derivs=(el_pt->nvorticity_derivatives_to_recover()+
			      el_pt->nvelocity_derivatives_to_recover());
 
 // Uniform mesh refinements
 unsigned n=4;
 
 // Loop over mesh refinements
 for (unsigned ii=0;ii<n;ii++)
 {  
  // Refine the mesh if we've already recovered the vorticity with this mesh
  if (ii>0)
  {
   // Refine the mesh
   refine_uniformly();
  }
  
  // Assign the synthetic velocity field
  assign_synthetic_veloc_field();

  // Create some space in the output
  oomph_info << std::endl;
  
  // Recover the vorticity field from the velocity field
  Vorticity_recoverer_pt->recover_vorticity(mesh_pt());

  // Create some space in the output
  oomph_info << std::endl;
  
  // Storage for the area of the mesh
  double full_area=0.0;
  
  // Storage for the error in the projection
  Vector<double> full_error(n_recovered_derivs,0.0);

  // Get the number of elements in the mesh
  unsigned n_el=mesh_pt()->nelement();

  // Output the number of elements and the square root of its inverse
  some_file << n_el << " " << sqrt(1.0/double(n_el)) << " ";

  // Now loop over the elements in the mesh
  for (unsigned e=0;e<n_el;e++)
  {
   // Upcast the e-th element in the mesh
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   // Get the size of the element
   double size=el_pt->size();

   // Update the mesh area variable
   full_area+=size;

   // Now loop over the recovered fields
   for (unsigned i=0;i<n_recovered_derivs;i++)
   {
    // Calculate the error in the i-th recovered field
    double el_error=el_pt->vorticity_error_squared(i);

    // Update the total error value
    full_error[i]+=el_error;
   }
  } // for (unsigned e=0;e<nel;e++)

  // Now loop over the recovered derivatives
  for (unsigned i=0;i<n_recovered_derivs;i++)
  {
   // Output the error in the i-th recovered quantity
   some_file << sqrt(full_error[i]) << " ";
  }

  // Output the total area of the mesh
  some_file << full_area << " " << std::endl;
   
  // Document the solution and indicate that we've recovered the vorticity
  doc_solution(doc_info,true);
 } // for (unsigned ii=0;ii<n;ii++)

 // Close the output file
 some_file.close();
} // End of check_smoothed_vorticity


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void VorticityRecoveryProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
					const bool& vorticity_recovered)
{ 
 // Check if we've recovered the vorticity
 if (!vorticity_recovered)
 {
  // Reconstruct smooth vorticity
  Vorticity_recoverer_pt->recover_vorticity(mesh_pt());
 }

 // Create an output stream
 ofstream some_file;

 // Storage for the file name
 char filename[100];

 // Number of plot points
 unsigned npts=2; 

 // Create the file name
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());

 // Open a file with the chosen file name
 some_file.open(filename);
 
 // Output solution 
 mesh_pt()->output(some_file,npts);

 // Now close the file; we're done with it
 some_file.close();

 // Output analytical vorticity and derivatives -- uses fake
 // (zero) data for velocities and pressure
 sprintf(filename,"%s/analytical_vorticity_and_indicator%i.dat",
	 doc_info.directory().c_str(),
	 doc_info.number());

 // Open a file with the chosen file name
 some_file.open(filename);

 // Get the number of elements in the mesh
 unsigned n_el=mesh_pt()->nelement();

 // Loop over the elements in the mesh
 for (unsigned e=0;e<n_el;e++)
 {
  // Upcast the e-th element in the mesh
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

  // Output the analytical vorticity and derivatives
  el_pt->output_analytical_veloc_and_vorticity(some_file,npts);
 }

 // Now close the file
 some_file.close();

 // Increment the Doc_info object counter
 doc_info.number()++;
} // End of doc_solution   


//===start_of_main======================================================
/// Driver code for Anne channel problem
//======================================================================
int main(int argc, char* argv[]) 
{

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
 
 // Switch off output modifier
 oomph_info.output_modifier_pt()=&default_output_modifier;

 // Switch off oomph_info output for all processors but rank 0
 if (MPI_Helpers::communicator_pt()->my_rank()!=0)
 {
  oomph_info.stream_pt()=&oomph_nullstream;
  OomphLibWarning::set_stream_pt(&oomph_nullstream);
  OomphLibError::set_stream_pt(&oomph_nullstream);
 }
 else
 {
  oomph_info << "\n\n====================================================="
	     << "\nNumber of processors: "
	     << MPI_Helpers::communicator_pt()->nproc()
	     << "\n=====================================================\n\n";
 }
#endif

 // Create doc info object
 DocInfo doc_info;

 // Set the output directory
 doc_info.set_directory("RESLT");

 // Typedef element
 typedef ProjectableTaylorHoodElement<RefineableQTaylorHoodElement<2> > ELEMENT;
 
 // Set up problem
 VorticityRecoveryProblem<VorticitySmootherElement<ELEMENT> > problem;

 // Check/document the smoothed vorticity field
 problem.check_smoothed_vorticity(doc_info);
 
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
} // end of main

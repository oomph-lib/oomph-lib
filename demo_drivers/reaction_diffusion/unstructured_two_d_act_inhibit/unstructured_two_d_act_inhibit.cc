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
//Driver for a simple 2D reaction diffusioni problem with
//adaptive mesh refinement on an unstructured mesh
//The problem is the "hot spot" expansion given in Muratov and Osipov
//(see two_d_act_inhibit.cc for full reference)
//Starting with a radially symmetric spot gives slightly slower growth
//that with an eccentric one, but the general dynamics are the same.

//Generic routines
#include "generic.h"
#include "advection_diffusion_reaction.h"
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph;


namespace GlobalVariables
{
 /// The vector of timescales
 Vector<double> Tau(2,1.0);
 
 /// The vector of diffusion coefficients
 Vector<double> D(2,1.0);

 /// Bifurcation (loss) parameter
 double A = -0.1;

 //Simple reaction kinetics
 void activator_inhibitor_reaction(const Vector<double> &C, Vector<double> &R)
 {
  //Inhibitor loss is linearly proportional to concentrations of activator
  //and inhibitor
  R[0] = C[0] + C[1] - A;
  //Activator growth is linearly proportional to activator and inhibitor
  //and loss is self-catalysed at a cubic rate
  R[1] = C[1]*C[1]*C[1] - C[1] - C[0];
 }

 /// Derivative of simple reaction kinetics
 void activator_inhibitor_reaction_derivative(const Vector<double> &C,
                                 DenseMatrix<double> &dRdC)
 {
  dRdC(0,0) = 1.0; dRdC(0,1) = 1.0;
  dRdC(1,0) = -1.0; dRdC(1,1) = 3.0*C[1]*C[1] - 1.0;
 }

}


//====== start_of_problem_class=======================================
/// 2D ActivatorInhibitor problem on rectangular domain, discretised 
/// with refineable 2D QActivatorInhibitor elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableActivatorInhibitorProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source and wind functions
 RefineableActivatorInhibitorProblem();

 /// Destructor. Empty
 ~RefineableActivatorInhibitorProblem()
  {
   //Delete the mesh
   delete My_mesh_pt; My_mesh_pt=0;
   //Delete the timestepper
   delete this->time_stepper_pt(); this->time_stepper_pt() = 0;
  }

 /// Empty
 void actions_before_newton_solve() {} 

 // Empty
 void actions_after_newton_solve(){}

 //Nothing to do before adaptation
 void actions_before_adapt() {}

 //Need to reset all the parameters after adaptation
 void actions_after_adapt()
  {
   complete_problem_setup();
  }
 
 //Set the initial condition
 void set_initial_condition();

 //Assign physical properties etc
 void complete_problem_setup();
 
 //Set the timestep
 void timestep(const double &dt, const unsigned &nstep);


private:


 /// Pointers to specific mesh
 RefineableTriangleMesh<ELEMENT>* My_mesh_pt;
 
 /// Storage for the timestep
 double Dt;

}; // end of problem class



//=====start_of_constructor===============================================
/// Constructor for ActivatorInhibitor problem
//========================================================================
template<class ELEMENT>
RefineableActivatorInhibitorProblem<ELEMENT>::
RefineableActivatorInhibitorProblem() : Dt(0.1)
{ 
 //Allocate the timestepper 
 add_time_stepper_pt(new BDF<2>);

 // Setup initial mesh

 // Domain length in x-direction
 double l_x=10.0;

 // Domain length in y-direction
 double l_y=10.0;

 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polylines
 //------------------------
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 
 // Each polyline only has three vertices -- provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord(2);
 for(unsigned i=0;i<2;i++) {vertex_coord[i].resize(2);}
 
 // First polyline: Base
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=l_x;
 vertex_coord[1][1]=0.0;
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,0);
 
 // Second boundary polyline: Right-hand side
 vertex_coord[0][0]=l_x;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=l_x;
 vertex_coord[1][1]=l_y;

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,1);

 // Third boundary polyline: Upper wall
 vertex_coord[0][0]=l_x;
 vertex_coord[0][1]=l_y;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=l_y;

 // Build the 3rd boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,2);

 // Fourth boundary polyline: Left-hand side
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=l_y;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=0.0;

 // Build the 4th boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,3);
 
 // Create the triangle mesh polygon for outer boundary
 TriangleMeshPolygon* Outer_boundary_polyline_pt =
  new TriangleMeshPolygon(boundary_polyline_pt);
  
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------

 // Convert to "closed curve" objects
 TriangleMeshClosedCurve* outer_closed_curve_pt=Outer_boundary_polyline_pt;

 unsigned n_internal_closed_boundaries = 0;
 Vector<TriangleMeshClosedCurve *>
  inner_boundaries_pt(n_internal_closed_boundaries);
 
 // Target area for initial mesh
 double uniform_element_area=1.0;

 // Use the TriangleMeshParameter object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

 //Define the inner boundaries
 triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundaries_pt;
 
 // Define the maximum element area
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 My_mesh_pt = new RefineableTriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());

 // Store as the problem's one and only mesh
 Problem::mesh_pt()=My_mesh_pt;
 
 // Create/set error estimator
 My_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 //Totally free (zero-flux) boundary conditions, so we do nothing
 
 // Complete the build of all elements so they are fully functional 
 complete_problem_setup();
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//=====================================================================
/// Set the initial conditions to be a single "hot" spot
//=====================================================================
template<class ELEMENT>
void RefineableActivatorInhibitorProblem<ELEMENT>::set_initial_condition()
{
 //Set the number of hot spots
 const unsigned n_spot = 1;
 //Set the centre of the hot spot
 const double centre_x[n_spot] = {6.0};
 const double centre_y[n_spot] = {7.0};
 
 //Set the initial concentrations of the reagents --- the homogeneous state
 //Plus an exponential hot spot
 unsigned n_node = mesh_pt()->nnode();
 //Loop over the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Local pointer to the node
   Node* nod_pt = mesh_pt()->node_pt(n);
   //Get the absolute value of A
   double a13 = pow(std::abs(GlobalVariables::A),(1.0/3.0));

   //Set the value of the inhibitor from the homogeneous solution
   nod_pt->set_value(0,a13*(1.0-a13*a13));
   
   //Set a localised hot spot by making an exponential
   //hump in the activator
   double x = nod_pt->x(0);    double y = nod_pt->x(1);
   double spot = 0.0;

   for(unsigned s=0; s<n_spot; s++)
    {
     //Find the square distance from the centre of the 
     //hot spot
     double r2 = (x - centre_x[s])*(x - centre_x[s]) + 
                 (y - centre_y[s])*(y - centre_y[s]);
     
     //Add an exponential hump to the value of the spot variable
     spot += 2.0*exp(-2.0*r2);
    }
   
   //Set the value of the activator
   nod_pt->set_value(1,-a13 + spot);
  }

 //Document the initial solution
 ofstream filename("input.dat");
 mesh_pt()->output(filename,5);
 filename.close();

 //Set the initial values impulsive
 assign_initial_values_impulsive(Dt);
} 





//=====================================================================
/// Set the elemental properties
//=====================================================================
template<class ELEMENT>
void RefineableActivatorInhibitorProblem<ELEMENT>::complete_problem_setup()
{

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the timescales
   elem_pt->tau_pt() = &GlobalVariables::Tau;

   //Set the diffusion coefficients
   elem_pt->diff_pt() = &GlobalVariables::D;

   //Set the reaction terms
   elem_pt->reaction_fct_pt() = &GlobalVariables::activator_inhibitor_reaction;

    //And their derivatives
   elem_pt->reaction_deriv_fct_pt() = 
    &GlobalVariables::activator_inhibitor_reaction_derivative;
 
  }

 
}


//====================================================================
/// Timestep the problem
//===================================================================
template<class ELEMENT>
void RefineableActivatorInhibitorProblem<ELEMENT>::timestep(
 const double &dt, const unsigned &nstep)
{
//Open a trace file
 std::ofstream tracefile("trace.dat");

//Set the problem's Dt for the inital condition bit
 Dt = dt;

 //One uniform refinement so that we stand a chance of getting 
 //something
 this->refine_uniformly();
 //Maximum adaptation for the first timestep
 unsigned max_adapt = 2;

 //Take the first timestep
 bool first = true;
 
 //Set the initial condition
 set_initial_condition();

 //Compute the norm
 double norm=0.0;
 
 //Solve the first step
 unsteady_newton_solve(dt,max_adapt,first);
 {
  unsigned i=0;
  char file1[100];
  sprintf(file1,"step%i.dat",i+1);
  ofstream out1(file1);
  mesh_pt()->output(out1,5);
  out1.close();
  //Output the norm
  mesh_pt()->compute_norm(norm);
  tracefile << this->time_pt()->time() << " " << sqrt(norm) << std::endl;
 }

 //Now set so that only one round of adaptation is performed each timestep
 max_adapt = 1;
 first = false;

 //Loop over possibilities
 for(unsigned i=1;i<nstep;i++)
  {
   unsteady_newton_solve(dt,max_adapt,first);
   char file1[100];
   sprintf(file1,"step%i.dat",i+1);
   ofstream out1(file1);
   mesh_pt()->output(out1,5);
   out1.close();
   //Output the norm
   mesh_pt()->compute_norm(norm);
   tracefile << this->time_pt()->time() << " " << sqrt(norm) << std::endl;
  }
 tracefile.close();
}

//===== start_of_main=====================================================
/// Driver code for 2D ActivatorInhibitor problem
//========================================================================
int main()
{
 //Diffusive length-scale
 double epsilon = 0.05;
 //Set the control parameters
 GlobalVariables::Tau[1] = 0.2;
 GlobalVariables::D[1] = epsilon*epsilon;
 GlobalVariables::A = -0.4;
 
 //Set the (largeish) timestep
 double dt = 0.25;
 
 //Construct the problem
 RefineableActivatorInhibitorProblem 
   <ProjectableAdvectionDiffusionReactionElement<TAdvectionDiffusionReactionElement<2,2,3> > > problem;

 //Timestep it
 problem.timestep(dt,2);
}










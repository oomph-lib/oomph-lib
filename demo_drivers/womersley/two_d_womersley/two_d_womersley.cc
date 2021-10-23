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
//Driver for 2D Womersley problem

//Generic routines
#include "generic.h"

// The Womersley equations
#include "womersley.h"
#include "navier_stokes.h"

// Meshes
#include "meshes/rectangular_quadmesh.h"
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



//======start_of_GlobalPhysicalParameters=============================
/// Namespace for Womersley problem
//====================================================================
namespace GlobalPhysicalParameters
{

 /// Prescribed volume flux -- must be assigned to Prescribed_volume_flux
 double prescribed_volume_flux(const double& time)
 {
  return 2.0*cos(2.0*MathematicalConstants::Pi*time);
 }

 /// Prescribed pressure gradient
 double prescribed_pressure_gradient(const double& time)
 {
  return 2.0*cos(2.0*MathematicalConstants::Pi*time);
 }

 /// Prescribed volume flux
 double Prescribed_volume_flux=0.0;

 /// Womersley number
 double Alpha_sq=10.0;

 /// Length of impedance tube
 double L_impedance=7.0;

} // end of GlobalPhysicalParameters




/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



//====================================================================
/// Specific Womersley impedance tube for a rectangular cross-section.
//====================================================================
template<class ELEMENT>
class RectangularWomersleyImpedanceTube  : 
 public WomersleyImpedanceTubeBase<ELEMENT,2>
{

public:

 /// Constructor: Pass length and function that prescribes the volume
 /// flux to constructor of underlying base class
 RectangularWomersleyImpedanceTube(
  const double& length,  
  typename WomersleyImpedanceTubeBase<ELEMENT,2>::PrescribedVolumeFluxFctPt 
  prescribed_volume_flux_fct_pt) : 
  WomersleyImpedanceTubeBase<ELEMENT,2>(length, prescribed_volume_flux_fct_pt)
  {}


 /// Implement pure virtual fct (defined in the base class 
 /// WomersleyImpedanceTubeBase) that builds the mesh of Womersley elements
 /// (of the type specified by the template argument), using the
 /// specified timestepper. Also applies the boundary condition.
 Mesh* build_mesh_and_apply_boundary_conditions(TimeStepper* time_stepper_pt)
  {

   // Setup mesh
   
   // Number of elements in x and y directions
   unsigned nx=5;
   unsigned ny=5;
   
   // Lengths in x and y directions
   double lx=1.0;
   double ly=1.0;
   
   // Build mesh
   Mesh* my_mesh_pt = 
    new RectangularQuadMesh<ELEMENT >(nx,ny,lx,ly,time_stepper_pt);
   
 
   // Set the boundary conditions for this problem: 
   
   // All nodes are free by default -- just pin the ones that have 
   // Dirichlet conditions here. 
   unsigned n_bound = my_mesh_pt->nboundary();
   for(unsigned b=0;b<n_bound;b++)
    {
     unsigned n_node = my_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       my_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    } // end of set boundary conditions

   return my_mesh_pt;

  }

};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




 //=======start_of_main====================================================
 /// Driver code for Womersley problem driven by outflow
 /// from Navier-Stokes mesh
 //========================================================================
 void run_navier_stokes_outflow()
 {

  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT_navier_stokes");

  // Output number
  doc_info.number()=0;

  // Open a trace file
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);


  // Create Navier-Stokes mesh
  
  // Create geometric objects: Elliptical tube with half axes = radius = 1.0
  double radius=1.0;
  GeomObject* Wall_pt=new EllipticalTube(radius,radius);
  
  // Boundaries on object
  Vector<double> xi_lo(2);
  xi_lo[0]=0.0;
  xi_lo[1]=0.0;
  Vector<double> xi_hi(2);
  xi_hi[0]=7.0;
  xi_hi[1]=2.0*atan(1.0);
  
  // # of layers
  unsigned nlayer=6;
  
  //Radial divider is located half-way along the circumference
  double frac_mid=0.5;
  
  // Build volume mesh
  RefineableQuarterTubeMesh<RefineableQTaylorHoodElement<3> >* n_st_mesh_pt=
   new RefineableQuarterTubeMesh<RefineableQTaylorHoodElement<3> >(
    Wall_pt,xi_lo,frac_mid,xi_hi,nlayer);
  
  n_st_mesh_pt->refine_uniformly();
  n_st_mesh_pt->refine_uniformly();

  
  // Apply boundary conditons: Pin axial velocity on wall (boundary 3)
  unsigned bound=3;
  unsigned num_nod=n_st_mesh_pt->nboundary_node(bound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    n_st_mesh_pt->boundary_node_pt(bound,inod)->pin(2);
   }
  

  // Create outflow mesh
  Mesh* Outflow_traction_mesh_pt=new Mesh;

  // Populate it...
  unsigned b=n_st_mesh_pt->nboundary()-1;
  {
   // Loop over all elements on these boundaries
   unsigned n_bound_el=n_st_mesh_pt->nboundary_element(b);
   for (unsigned e=0;e<n_bound_el;e++)
    { 
     // Get pointer to bulk element
     RefineableQTaylorHoodElement<3>* bulk_elem_pt = 
      dynamic_cast<RefineableQTaylorHoodElement<3>*>(
       n_st_mesh_pt->boundary_element_pt(b,e));
     
     //What is the index of the face of element e along boundary b
     int face_index = n_st_mesh_pt->face_index_at_boundary(b,e);
          
     // Build the corresponding prescribed traction element
     NavierStokesImpedanceTractionElement<RefineableQTaylorHoodElement<3>,
      QWomersleyElement<2,3>,2>* 
      traction_element_pt 
      = new NavierStokesImpedanceTractionElement<
      RefineableQTaylorHoodElement<3>,
      QWomersleyElement<2,3>,2>(bulk_elem_pt,face_index);
     
     //Add the prescribed traction element to the mesh
     Outflow_traction_mesh_pt->add_element_pt(traction_element_pt);
       
    }
  }
 

  // Choose simulation interval and timestep
  double t_max=5.0;
  double dt=0.05;

  // Build Womersley impedance tube
  unsigned fixed_coordinate=2;
  unsigned w_index=2;
  WomersleyOutflowImpedanceTube<QWomersleyElement<2,3>,2>* 
   womersley_impedance_tube_pt=
   new WomersleyOutflowImpedanceTube<QWomersleyElement<2,3>,2>
   (GlobalPhysicalParameters::L_impedance,
    Outflow_traction_mesh_pt,
    fixed_coordinate,
    w_index);
  
  // Set up steady flow in impedance tube, imposing the same 
  // flow rate that "comes out of" the 3D tube.
  double q_initial=womersley_impedance_tube_pt->
   total_volume_flux_into_impedance_tube();
  womersley_impedance_tube_pt->
   setup(&GlobalPhysicalParameters::Alpha_sq,dt,q_initial);
  
  
  //Output Womersley solution
  womersley_impedance_tube_pt->womersley_problem_pt()->
   doc_solution(doc_info,trace_file, xi_hi[0]);
 

//   // Output Navier Stokes "solution"
//   ofstream some_file;
  
//   // Number of plot points
//   unsigned npts;
//   npts=5;
  
//   sprintf(filename,"%s/navier_stokes_soln%i.dat",
//           doc_info.directory().c_str(),
//           doc_info.number());
//   some_file.open(filename);
//   n_st_mesh_pt->output(some_file,npts);
//   some_file.close();
  
  //Increment counter for solutions 
  doc_info.number()++;
  
  // Find number of steps
  unsigned nstep = static_cast<unsigned>(round(t_max/dt));

  // Timestepping loop
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Shift timesteps. NOTE: In NSt problem this needs to go into
    // the 3D Problem's actions_before_implicit_time_step()
    womersley_impedance_tube_pt->shift_time_values(dt);
 

    // Crank up flow rate in 3D tube
    unsigned nnod= n_st_mesh_pt->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      Node* nod_pt=n_st_mesh_pt->node_pt(j);
      double x=nod_pt->x(0);
      double y=nod_pt->x(1);
      double time=womersley_impedance_tube_pt->womersley_problem_pt()->
       time_pt()->time();
      nod_pt->set_value(2,(1.0-x*x-y*y)*
                        sin(2.0*MathematicalConstants::Pi*time));            
     }


    // Get response
    double p_in=0.0;
    double dp_in_dq=0.0;
    womersley_impedance_tube_pt->get_response(p_in,dp_in_dq);

    //Output Womersley solution
    womersley_impedance_tube_pt->womersley_problem_pt()->
     doc_solution(doc_info,trace_file,xi_hi[0]);

//     // Output Navier Stokes "solution"
//     sprintf(filename,"%s/navier_stokes_soln%i.dat",
//             doc_info.directory().c_str(),
//             doc_info.number());
//     some_file.open(filename);
//     n_st_mesh_pt->output(some_file,npts);
//     some_file.close();
    
    //Increment counter for solutions 
    doc_info.number()++;

   }

  // Close trace file
  trace_file.close();


 }


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




 //=======start_of_main====================================================
 /// Driver code for Womersley problem based on impedance tube
 //========================================================================
 void run_impedance_tube()
 {

  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT_impedance_tube");

  // Output number
  doc_info.number()=0;

  // Open a trace file
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);

  // Choose simulation interval and timestep
  double t_max=5.0;
  double dt=0.05;

  // Build Womersley impedance tube
  RectangularWomersleyImpedanceTube<QWomersleyElement<2,4> >* 
   womersley_impedance_tube_pt=
   new RectangularWomersleyImpedanceTube<QWomersleyElement<2,4> >
   (GlobalPhysicalParameters::L_impedance,
    &GlobalPhysicalParameters::prescribed_volume_flux);

  // Set initial flux
  GlobalPhysicalParameters::Prescribed_volume_flux=
   GlobalPhysicalParameters::prescribed_volume_flux(0.0);

  // Set up steady flow in impedance tube. NOTE: In NSt problem this 
  // needs to be called when the outflow from the 3D NSt problem
  // has been assigned (implicitly by assigning its velocity dofs)
  womersley_impedance_tube_pt->
   setup(&GlobalPhysicalParameters::Alpha_sq,
         dt,GlobalPhysicalParameters::Prescribed_volume_flux);
  
  
  //Output solution
  womersley_impedance_tube_pt->womersley_problem_pt()->
   doc_solution(doc_info,trace_file);
 
  //Increment counter for solutions 
  doc_info.number()++;
  
  // Find number of steps
  unsigned nstep = static_cast<unsigned>(round(t_max/dt));

  // Timestepping loop
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Shift timesteps. NOTE: In NSt problem this needs to go into
    // the 3D Problem's actions_before_implicit_time_step().
    // Note: The volume flux is prescribed by a function pointer
    // so it updates itself.
    womersley_impedance_tube_pt->shift_time_values(dt);
 
    // Get response
    double p_in=0.0;
    double dp_in_dq=0.0;
    womersley_impedance_tube_pt->get_response(p_in,dp_in_dq);

    //Output solution
    womersley_impedance_tube_pt->womersley_problem_pt()->
     doc_solution(doc_info,trace_file);
 
    //Increment counter for solutions 
    doc_info.number()++;

   }

  // Close trace file
  trace_file.close();


 }


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////





//========================================================================
/// Normal driver code for Womersley problem with prescribed flux
//========================================================================
void run_prescribed_flux()
{

 // Create timestepper
 TimeStepper* time_stepper_pt=new BDF<2>;

 // Setup mesh
 
 // Number of elements in x and y directions
 unsigned nx=5;
 unsigned ny=5;
 
 // Lengths in x and y directions
 double lx=1.0;
 double ly=1.0;
 
 // Build mesh
 Mesh* my_mesh_pt = 
  new RectangularQuadMesh<QWomersleyElement<2,4> >(nx,ny,lx,ly,
                                                   time_stepper_pt);
 

 // Set the boundary conditions for this problem: 
 
 // All nodes are free by default -- just pin the ones that have 
 // Dirichlet conditions here. 
 unsigned n_bound = my_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = my_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     my_mesh_pt->boundary_node_pt(b,n)->pin(0); 
    }
  } // end of set boundary conditions
 


 // Build problem with specifically created mesh (BCs applied and
 // timestepper)
 WomersleyProblem<QWomersleyElement<2,4>,2> 
  problem(&GlobalPhysicalParameters::Alpha_sq,
          &GlobalPhysicalParameters::Prescribed_volume_flux,
          time_stepper_pt,my_mesh_pt);
 
 // Setup labels for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT_prescribed_volume_flux");
 
 // Output number
 doc_info.number()=0;
 
 // Open a trace file
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);

 // Choose simulation interval and timestep
 double t_max=5.0;
 double dt=0.05;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);

 // Assign current flux
 GlobalPhysicalParameters::Prescribed_volume_flux=
  GlobalPhysicalParameters::prescribed_volume_flux(0.0);
 
 // Do steady solve and impulsive start from it
 problem.steady_newton_solve();

 //Output initial condition
 problem.doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Reuse Jacobian during time-dependent run
 problem.enable_jacobian_reuse();

 // Find number of steps
 unsigned nstep = static_cast<unsigned>(round(t_max/dt));

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Assign flux at next timestep
   GlobalPhysicalParameters::Prescribed_volume_flux=
    GlobalPhysicalParameters::prescribed_volume_flux(
     problem.time_pt()->time()+dt);
   
   // Take timestep
   problem.unsteady_newton_solve(dt);
   
   //Output solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;

  }
 
 // Close trace file
 trace_file.close();

}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




//========================================================================
/// Normal driver code for Womersley problem with prescribed 
/// pressure gradient
//========================================================================
void run_prescribed_pressure_gradient()
{

 // Create timestepper
 TimeStepper* time_stepper_pt=new BDF<2>;

 // Setup mesh
 
 // Number of elements in x and y directions
 unsigned nx=5;
 unsigned ny=5;
 
 // Lengths in x and y directions
 double lx=1.0;
 double ly=1.0;
 
 // Build mesh
 Mesh* my_mesh_pt = 
  new RectangularQuadMesh<QWomersleyElement<2,4> >(nx,ny,lx,ly,
                                                   time_stepper_pt);
 

 // Set the boundary conditions for this problem: 
 
 // All nodes are free by default -- just pin the ones that have 
 // Dirichlet conditions here. 
 unsigned n_bound = my_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = my_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     my_mesh_pt->boundary_node_pt(b,n)->pin(0); 
    }
  } // end of set boundary conditions
 


 // Build problem with specifically created mesh (BCs applied and
 // timestepper)
 WomersleyProblem<QWomersleyElement<2,4>,2> 
  problem(&GlobalPhysicalParameters::Alpha_sq,
          &GlobalPhysicalParameters::prescribed_pressure_gradient,
          time_stepper_pt,my_mesh_pt);
 
 // Setup labels for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT_prescribed_pressure_gradient");
 
 // Output number
 doc_info.number()=0;
 
 // Open a trace file
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);

 // Choose simulation interval and timestep
 double t_max=5.0;
 double dt=0.05;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);

 // Assign initial pressure gradient
 problem.pressure_gradient_data_pt()->set_value(
  0,
  GlobalPhysicalParameters::prescribed_pressure_gradient(0.0));

 
 // Do steady solve and impulsive start from it
 problem.steady_newton_solve();

 //Output initial condition
 problem.doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Reuse Jacobian during time-dependent run
 problem.enable_jacobian_reuse();

 // Find number of steps
 unsigned nstep = static_cast<unsigned>(round(t_max/dt));

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {   
   // Take timestep
   problem.unsteady_newton_solve(dt);
   
   //Output solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;

  }
 
 // Close trace file
 trace_file.close();

}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////





//=======start_of_main====================================================
/// Driver code for Womersley problem
//========================================================================
int main()
{

 run_navier_stokes_outflow();
 
 run_prescribed_pressure_gradient();

 run_prescribed_flux();

 run_impedance_tube();


};
 

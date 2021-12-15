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
// Driver for axisymmetric single-layer fluid problem 
 
// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// Bessel functions
#include "oomph_crbond_bessel.h"

// The mesh
#include "meshes/single_layer_spine_mesh.h"

using namespace std;

using namespace oomph;



//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re = 5.0;

 /// Womersley number
 double ReSt = 5.0; // (St = 1)
 
 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 5.0; // (Fr = 1)

 /// Capillary number
 double Ca = 0.01;

 /// Direction of gravity
 Vector<double> G(3);

} // End of namespace



//==start_of_problem_class===============================================
/// Single axisymmetric fluid interface problem in rectangular domain
//=======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:
 
 /// Constructor: Pass the number of elements and the lengths of the
 /// domain in the r and z directions (h is the height of the fluid layer
 /// i.e. the length of the domain in the z direction)
 InterfaceProblem(const unsigned &n_r, const unsigned &n_z, 
                  const double &l_r, const double &h);
 
 /// Destructor (empty)
 ~InterfaceProblem() {}

 /// Spine heights/lengths are unknowns in the problem so their values get
 /// corrected during each Newton step. However, changing their value does
 /// not automatically change the nodal positions, so we need to update all
 /// of them here.
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 // Update before solve (empty)
 void actions_before_newton_solve() {}

 /// Update after solve can remain empty, because the update 
 /// is performed automatically after every Newton step.
 void actions_after_newton_solve() {}

 /// Set initial conditions: Set all nodal velocities to zero and
 /// initialise the previous velocities to correspond to an impulsive
 /// start
 void set_initial_condition()
  {
   // Determine number of nodes in mesh
   const unsigned n_node = Bulk_mesh_pt->nnode();

   // Loop over all nodes in mesh
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       // Set velocity component i of node n to zero
       Bulk_mesh_pt->node_pt(n)->set_value(i,0.0);
      }
    }

   // Initialise the previous velocity values for timestepping
   // corresponding to an impulsive start
   assign_initial_values_impulsive();

  } // End of set_initial_condition

 /// Access function for the specific mesh
 SingleLayerSpineMesh<ELEMENT> *Bulk_mesh_pt;

 /// Mesh for the interface elements
 Mesh* Surface_mesh_pt;

 /// Doc the solution
 void doc_solution(DocInfo &doc_info);

 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double &t_max, const double &dt); 
 

private:

 /// Deform the mesh/free surface to a prescribed function
 void deform_free_surface(const double &epsilon, const double &k)
  {
   // Initialise Bessel functions (only need the first!)
   double j0, j1, y0, y1, j0p, j1p, y0p, y1p;

   // Determine number of spines in mesh
   const unsigned n_spine = Bulk_mesh_pt->nspine();

   // Loop over spines in mesh
   for(unsigned i=0;i<n_spine;i++)
    {
     
     // Determine r coordinate of spine
     double r_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(0);

     // Get Bessel functions J_0(kr), J_1(kr), Y_0(kr), Y_1(kr)
     // and their derivatives
     CRBond_Bessel::bessjy01a(k*r_value,j0,j1,y0,y1,j0p,j1p,y0p,y1p);

     // Set spine height
     Bulk_mesh_pt->spine_pt(i)->height() = 1.0 + epsilon*j0;

    } // End of loop over spines
   
   // Update nodes in bulk mesh
   Bulk_mesh_pt->node_update();

  } // End of deform_free_surface

 /// Trace file
 ofstream Trace_file;

 /// Width of domain
 double Lr;

}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for single fluid interface problem
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::
InterfaceProblem(const unsigned &n_r, const unsigned &n_z,
                 const double &l_r, const double& h) : Lr(l_r)
{

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER);

 // Build and assign mesh (the "false" boolean flag tells the mesh
 // constructor that the domain is not periodic in r)
 Bulk_mesh_pt = new SingleLayerSpineMesh<ELEMENT>
  (n_r,n_z,l_r,h,false,time_stepper_pt());


 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 for(unsigned i=0;i<n_r;i++)
  {
   //Construct a new 1D line element on the face on which the local
   //coordinate 1 is fixed at its max. value (1) --- This is face 2
   FiniteElement *interface_element_pt =
    new SpineAxisymmetricFluidInterfaceElement<ELEMENT>(
     Bulk_mesh_pt->finite_element_pt(n_r*(n_z-1)+i),2);
   
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 
  }
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here

 // Determine number of mesh boundaries
 const unsigned n_boundary = Bulk_mesh_pt->nboundary();
 
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = Bulk_mesh_pt->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Pin radial and azimuthal velocities on all boundaries other
     // than top (no slip/penetration)
     if(b!=2)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
      }
     // Pin axial velocity on bottom wall only (no penetration)
     if(b==0)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
 // Determine total number of nodes in mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();

 // Pin all azimuthal velocities throughout the bulk of the domain
 for(unsigned n=0;n<n_node;n++)
  {
   Bulk_mesh_pt->node_pt(n)->pin(2);
  }

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------
 
 // Determine number of bulk elements in mesh
 const unsigned n_bulk = Bulk_mesh_pt->nelement();

 // Loop over the bulk elements
 for(unsigned e=0;e<n_bulk;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::G;

  } // End of loop over bulk elements

 // Create a Data object whose single value stores the external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Pin and set the external pressure to some arbitrary value
 external_pressure_data_pt->pin(0);
 external_pressure_data_pt->set_value(0,1.31);

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = Surface_mesh_pt->nelement();

 // Loop over the interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   SpineAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineAxisymmetricFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

  } // End of loop over interface elements

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor


   
//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo &doc_info)
{ 

 // Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 // Document in trace file
 Trace_file << time_pt()->time() << " "
            << Bulk_mesh_pt->spine_pt(0)->height() <<  std::endl;

 ofstream some_file;
 char filename[100];

 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;

 // Open solution output file
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Output solution to file
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);

 // Close solution output file
 some_file.close();
 
} // End of doc_solution

 

//==start_of_unsteady_run=================================================
/// Perform run up to specified time t_max with given timestep dt
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
unsteady_run(const double &t_max, const double &dt)
{

 // Set value of epsilon
 const double epsilon = 0.2;

 // Set value of k in Bessel function J_0(kr)
 const double k = 7.0156; // Value of the second zero of J_1(k)

 // Deform the mesh/free surface
 deform_free_surface(epsilon, k);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Initialise counter for solutions
 doc_info.number()=0;
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise trace file
 Trace_file << "time" << ", "
            << "edge spine height" << ", "
            << "contact angle left" << ", "
            << "contact angle right" << ", " << std::endl;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial conditions
 set_initial_condition();

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);

 // Doc initial solution
 doc_solution(doc_info);

 // Increment counter for solutions 
 doc_info.number()++;

 // Timestepping loop
 for(unsigned t=1;t<=n_timestep;t++)
  {
   // Output current timestep to screen
   cout << "\nTimestep " << t << " of " << n_timestep << std::endl;
   
   // Take one fixed timestep
   unsteady_newton_solve(dt);

   // Doc solution
   doc_solution(doc_info);

   // Increment counter for solutions 
   doc_info.number()++;

  } // End of timestepping loop

} // End of unsteady_run


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver code for single fluid axisymmetric interface problem
//=====================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 /// Maximum time
 double t_max = 0.8;

 /// Duration of timestep
 const double dt = 0.005;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::Argc>1) { t_max = 0.01; }

 // Number of elements in radial (r) direction
 const unsigned n_r = 16;
   
 // Number of elements in axial (z) direction
 const unsigned n_z = 12;

 // Width of domain
 const double l_r = 1.0;

 // Height of fluid layer
 const double h = 1.0;
 
 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = -1.0;
 Global_Physical_Variables::G[2] = 0.0; 

 // Set up the spine test problem with AxisymmetricQCrouzeixRaviartElements,
 // using the BDF<2> timestepper
 InterfaceProblem<SpineElement<AxisymmetricQCrouzeixRaviartElement >,BDF<2> >
  problem(n_r,n_z,l_r,h);
 
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);
 
} // End of main

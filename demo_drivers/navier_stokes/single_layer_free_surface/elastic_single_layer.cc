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
// Driver for two-dimensional single fluid free surface problem, where the
// mesh is deformed using a pseudo-solid node-update strategy
 
// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// Constitutive law headers
#include "constitutive.h"

// Solid headers
#include "solid.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//==start_of_namespace====================================================
/// Namespace for physical parameters
//========================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re = 5.0;

 /// Strouhal number
 double St = 1.0;

 /// Womersley number (Reynolds x Strouhal, computed automatically)
 double ReSt;
 
 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 5.0; // (Fr = 1)

 /// Capillary number
 double Ca = 0.01;

 /// Direction of gravity
 Vector<double> G(2);

 /// Pseudo-solid Poisson ratio
 double Nu = 0.1;

} // End of namespace


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//==start_of_problem_class================================================
/// Single fluid free surface problem in a rectangular domain which is
/// periodic in the x direction
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:
 
 /// Constructor: Pass the number of elements and the lengths of the
 /// domain in the x and y directions (h is the height of the fluid layer
 /// i.e. the length of the domain in the y direction)
 InterfaceProblem(const unsigned &n_x, const unsigned &n_y, 
                  const double &l_x, const double &h);
 
 /// Destructor (empty)
 ~InterfaceProblem() {}

 /// Set initial conditions
 void set_initial_condition();

 /// Set boundary conditions
 void set_boundary_conditions();

 /// Document the solution
 void doc_solution(DocInfo &doc_info);

 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double &t_max, const double &dt); 

private:

 /// No actions required before solve step
 void actions_before_newton_solve() {}

 /// No actions required after solve step
 void actions_after_newton_solve() {}

 /// Actions before the timestep: For maximum stability, reset
 /// the current nodal positions to be the "stress-free" ones.
 void actions_before_implicit_timestep()
  {
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

 /// Deform the mesh/free surface to a prescribed function
 void deform_free_surface(const double &epsilon, const unsigned &n_periods);
 
 /// Pointer to the (specific) "bulk" mesh
 ElasticRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 // Pointer to the constitutive law used to determine the mesh deformation
 ConstitutiveLaw* Constitutive_law_pt;

 /// Width of domain
 double Lx;

 /// Trace file
 ofstream Trace_file;

}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for single fluid free surface problem
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::
InterfaceProblem(const unsigned &n_x, const unsigned &n_y,
                 const double &l_x, const double& h) : Lx(l_x)
{

 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER);

 // Build and assign the "bulk" mesh (the "true" boolean flag tells
 // the mesh constructor that the domain is periodic in x)
 Bulk_mesh_pt = new ElasticRectangularQuadMesh<ELEMENT>
  (n_x,n_y,l_x,h,true,time_stepper_pt());

 // Create the "surface mesh" that will contain only the interface
 // elements. The constructor just creates the mesh without giving
 // it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // -----------------------------
 // Create the interface elements
 // -----------------------------

 // Loop over those elements adjacent to the free surface
 for(unsigned e=0;e<n_x;e++)
  {
   // Set a pointer to the bulk element we wish to our interface
   // element to
   FiniteElement* bulk_element_pt =
    Bulk_mesh_pt->finite_element_pt(n_x*(n_y-1)+e);

   // Create the interface element (on face 2 of the bulk element)
   FiniteElement* interface_element_pt =
    new ElasticLineFluidInterfaceElement<ELEMENT>(bulk_element_pt,2);

   // Add the interface element to the surface mesh
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
     // Fluid boundary conditions:
     // --------------------------

     // On lower boundary (solid wall), pin x and y components of
     // the velocity (no slip/penetration)
     if(b==0)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }

     // On left and right boundaries, pin x-component of the velocity
     // (no penetration of the periodic boundaries)
     if(b==1 || b==3)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
      }

     // Solid boundary conditions:
     // --------------------------

     // On lower boundary (solid wall), pin vertical displacement
     // (no penetration)
     if(b==0)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries

 // Pin horizontal displacement of all nodes
 const unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++) { Bulk_mesh_pt->node_pt(n)->pin_position(0); }
 
 // Define a constitutive law for the solid equations: generalised Hookean
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);

 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------
 
 // Determine number of bulk elements in mesh
 const unsigned n_element_bulk = Bulk_mesh_pt->nelement();

 // Loop over the bulk elements
 for(unsigned e=0;e<n_element_bulk;e++)
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

   // Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

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
   ElasticLineFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<ElasticLineFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   // Set the Strouhal number
   el_pt->st_pt() = &Global_Physical_Variables::St;

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

  } // End of loop over interface elements

 // Apply the boundary conditions
 set_boundary_conditions();

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor



//==start_of_set_initial_condition========================================
/// Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = mesh_pt()->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the two velocity components
   for(unsigned i=0;i<2;i++)
    {
     // Set velocity component i of node n to zero
     mesh_pt()->node_pt(n)->set_value(i,0.0);
    }
  }
 
 // Initialise the previous velocity values and nodal positions
 // for timestepping corresponding to an impulsive start
 assign_initial_values_impulsive();
 
} // End of set_initial_condition



//==start_of_set_boundary_conditions======================================
/// Set boundary conditions: Set both velocity components to zero
/// on the bottom (solid) wall and the horizontal component only to zero
/// on the side (periodic) boundaries
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::set_boundary_conditions()
{
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
     // Set x-component of the velocity to zero on all boundaries
     // other than the free surface
     if(b!=2)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(0,0.0);
      }
     
     // Set y-component of the velocity to zero on the bottom wall
     if(b==0)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(1,0.0);
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
} // End of set_boundary_conditions



//==start_of_deform_free_surface==========================================
/// Deform the mesh/free surface to a prescribed function
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
deform_free_surface(const double &epsilon,const unsigned &n_periods)
{
 // Determine number of nodes in the "bulk" mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Determine eulerian position of node
   const double current_x_pos = Bulk_mesh_pt->node_pt(n)->x(0);
   const double current_y_pos = Bulk_mesh_pt->node_pt(n)->x(1);
   
   // Determine new vertical position of node
   const double new_y_pos = current_y_pos
    + (1.0-fabs(1.0-current_y_pos))*epsilon
    *(cos(2.0*n_periods*MathematicalConstants::Pi*current_x_pos/Lx));
   
   // Set new position
   Bulk_mesh_pt->node_pt(n)->x(1) = new_y_pos;
  }
} // End of deform_free_surface



//==start_of_doc_solution=================================================
/// Document the solution
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo &doc_info)
{ 

 // Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 // Upcast from GeneralisedElement to the present element
 ElasticLineFluidInterfaceElement<ELEMENT>* el_pt = 
  dynamic_cast<ElasticLineFluidInterfaceElement<ELEMENT>*>
  (Surface_mesh_pt->element_pt(0));

 // Document time and vertical position of left hand side of interface
 // in trace file
 Trace_file << time_pt()->time() << " "
            << el_pt->node_pt(0)->x(1) << std::endl;
 
 ofstream some_file;
 char filename[100];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;
 
 // Open solution output file
 snprintf(filename, sizeof(filename), "%s/soln%i.dat",
         doc_info.directory().c_str(),doc_info.number());
 some_file.open(filename);

 // Output solution to file
 Bulk_mesh_pt->output(some_file,npts);

 // Close solution output file
 some_file.close();

 // Open interface solution output file
 snprintf(filename, sizeof(filename), "%s/interface_soln%i.dat",
         doc_info.directory().c_str(),doc_info.number());
 some_file.open(filename);
 
 // Output solution to file
 Surface_mesh_pt->output(some_file,npts);
 
 // Close solution output file
 some_file.close();
 
} // End of doc_solution

 

//==start_of_unsteady_run=================================================
/// Perform run up to specified time t_max with given timestep dt
//========================================================================
template <class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
unsteady_run(const double &t_max, const double &dt)
{

 // Set value of epsilon
 const double epsilon = 0.1;

 // Set number of periods for cosine term
 const unsigned n_periods = 1;

 // Deform the mesh/free surface
 deform_free_surface(epsilon,n_periods);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Initialise counter for solutions
 doc_info.number()=0;
 
 // Open trace file
 char filename[100];   
 snprintf(filename, sizeof(filename), "%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise trace file
 Trace_file << "time, free surface height" << std::endl;

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


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//==start_of_main=========================================================
/// Driver code for two-dimensional single fluid free surface problem
//========================================================================
int main(int argc, char* argv[]) 
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Compute the Womersley number
 Global_Physical_Variables::ReSt =
  Global_Physical_Variables::Re*Global_Physical_Variables::St;

 /// Maximum time
 double t_max = 0.6;

 /// Duration of timestep
 const double dt = 0.0025;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::Argc>1) { t_max = 0.005; }

 // Number of elements in x direction
 const unsigned n_x = 12;
   
 // Number of elements in y direction
 const unsigned n_y = 12;

 // Width of domain
 const double l_x = 1.0;

 // Height of fluid layer
 const double h = 1.0;
 
 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = -1.0;

 // Set up the elastic test problem with QCrouzeixRaviartElements,
 // using the BDF<2> timestepper
 InterfaceProblem<PseudoSolidNodeUpdateElement< QCrouzeixRaviartElement<2>,
  QPVDElement<2,3> > , BDF<2> >
  problem(n_x,n_y,l_x,h);
 
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);
 
} // End of main

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
//A demo driver that solves the classic fluid flow problem of flow
//of a fluid film along an inclined plane. Stability analysis performed by 
//Yih (1963), Benjamin (1957), ... and Blyth & Pozrikidis (2004).

//This is an example of the subtleties involved in even a seemingly simple
//free surface problem.

//Standard C++ library includes
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

//Finite-Element library routines
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "fluid_interface.h"
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables
{
 /// Reynolds number, based on the average velocity within the fluid film
 double Re=0.0;

 /// The product of Reynolds number and inverse Froude number 
 /// is set to two in this problem, which gives the free surface velocity
 /// to be sin(alpha). [Set to three in order to get the same scale as 
 /// used by Yih, Benjamin, etc]
 double ReInvFr=2.0;

 /// Angle of incline of the slope (45 degrees)
 double Alpha = 1.0*atan(1.0);

 /// The Vector direction of gravity, set in main()
 Vector<double> G(2,0.0);

 /// The Capillary number
 double Ca= 1.0;

 /// Set the wavenumber
 double K = 0.1;

 /// Set the number of waves desired in the domain
 double N_wave = 3;

 /// The length of the domain to fit the desired number of waves
 double Length = 2*N_wave*4.0*atan(1.0)/K;

 /// Direction of the wall normal vector (at the inlet)
 Vector<double> Wall_normal;

 /// Function that specifies the wall unit normal at the inlet
 void wall_unit_normal_inlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  normal=Wall_normal;
 }

 /// Function that specified the wall unit normal at the outlet
 void wall_unit_normal_outlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  //Set the normal
  normal = Wall_normal;
  //and flip the sign
  unsigned n_dim = normal.size();
  for(unsigned i=0;i<n_dim;++i) {normal[i] *= -1.0;}
 }

 /// The contact angle that is imposed at the inlet (pi)
 double Inlet_Angle = 2.0*atan(1.0);


 /// Function that prescribes the hydrostatic pressure field at the outlet
 void hydrostatic_pressure_outlet(const double& time, const Vector<double> &x, 
                                  const Vector<double> &n, 
                                  Vector<double> &traction)
 {
  traction[0] = ReInvFr*G[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet(const double& time, const Vector<double> &x, 
                                 const Vector<double> &n,
                                 Vector<double> &traction)
 {
  traction[0] = -ReInvFr*G[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }
 //end of traction functions

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

}

//=====================================================================
/// Generic problem class that will form the base class for both 
/// spine and elastic mesh-updates of the problem.
/// Templated by the bulk element and interface element types
//====================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
class InclinedPlaneProblem : public Problem
{

protected:

 /// Bulk fluid mesh
 Mesh* Bulk_mesh_pt;

 /// Mesh for the traction elements that are added at inlet and outlet
 Mesh* Traction_mesh_pt;

 /// Mesh for the free surface elements
 Mesh* Surface_mesh_pt;

 /// Mesh for the point elements at each end of the free surface
 Mesh* Point_mesh_pt;

 /// Prefix for output files
 std::string Output_prefix;

public:

 /// Generic Constructor (empty)
 InclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length) :
  Output_prefix("Unset") { }
 
 /// Solve the steady problem
 void solve_steady();
 
 /// Take n_tsteps timesteps of size dt
 void timestep(const double &dt, const unsigned &n_tsteps);

 /// Actions before the timestep 
 /// (update the time-dependent boundary conditions)
 void actions_before_implicit_timestep()
  {
   //Read out the current time
   double time = this->time_pt()->time();
   //Now add a temporary sinusoidal suction and blowing to the base
   //Amplitude of the perturbation
   double epsilon = 0.01;
   //Loop over the nodes on the base
   unsigned n_node = this->Bulk_mesh_pt->nboundary_node(0);
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = this->Bulk_mesh_pt->boundary_node_pt(0,n);
     double arg = Global_Physical_Variables::K*nod_pt->x(0);
     double value = sin(arg)*epsilon*time*exp(-time);
     nod_pt->set_value(1,value);
    }
  } //end_of_actions_before_implicit_timestep

 /// Function to add the traction boundary elements to boundaries
 /// 3(inlet) and 1(outlet) of the mesh
 void make_traction_elements()
  {
   //Create a new (empty mesh)
   Traction_mesh_pt = new Mesh;
   //Inlet boundary conditions (boundary 3)
   {
    unsigned b = 3;
    //Find the number of elements adjacent to mesh boundary
    unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
    //Loop over these elements and create the traction elements
    for(unsigned e=0;e<n_boundary_element;e++)
     {
      NavierStokesTractionElement<ELEMENT> *surface_element_pt =
       new NavierStokesTractionElement<ELEMENT>
       (Bulk_mesh_pt->boundary_element_pt(b,e),
        Bulk_mesh_pt->face_index_at_boundary(b,e));
      //Add the elements to the mesh
      Traction_mesh_pt->add_element_pt(surface_element_pt);
      //Set the traction function
      surface_element_pt->traction_fct_pt() = 
       &Global_Physical_Variables::hydrostatic_pressure_inlet;
     }
   }
   
   //Outlet boundary conditions (boundary 1)
   {
    unsigned b=1;
    //Find the number of elements adjacent to mesh boundary
    unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
    //Loop over these elements and create the traction elements
    for(unsigned e=0;e<n_boundary_element;e++)
     {
      NavierStokesTractionElement<ELEMENT> *surface_element_pt =
       new NavierStokesTractionElement<ELEMENT>
       (Bulk_mesh_pt->boundary_element_pt(b,e),
        Bulk_mesh_pt->face_index_at_boundary(b,e));
      //Add the elements to the mesh
      Traction_mesh_pt->add_element_pt(surface_element_pt);
      //Set the traction function
      surface_element_pt->traction_fct_pt() = 
       &Global_Physical_Variables::hydrostatic_pressure_outlet;
     }
   }
  } //end of make_traction_elements

 //Make the free surface elements on the top surface
 void make_free_surface_elements()
  {
   //Create the (empty) meshes
   Surface_mesh_pt = new Mesh;
   Point_mesh_pt = new Mesh;

   //The free surface is on the boundary 2
   unsigned b = 2;
   unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(b);
   //Loop over the elements and create the appropriate interface elements
   for(unsigned e=0;e<n_boundary_element;e++)
    {
     INTERFACE_ELEMENT *surface_element_pt =
      new INTERFACE_ELEMENT
      (Bulk_mesh_pt->boundary_element_pt(b,e),
       Bulk_mesh_pt->face_index_at_boundary(b,e));
     //Add elements to the mesh
     Surface_mesh_pt->add_element_pt(surface_element_pt);
     //Assign the capillary number to the free surface
     surface_element_pt->ca_pt() = 
      &Global_Physical_Variables::Ca;
     
     //Make a point element from left-hand side of the 
     //first surface element (note that this relies on knowledge of 
     //the element order within the mesh)
     if(e==0)
      {
       FluidInterfaceBoundingElement* point_element_pt =
        surface_element_pt->make_bounding_element(-1);
       //Add element to the point mesh
       Point_mesh_pt->add_element_pt(point_element_pt);
       //Set the capillary number
       point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
       //Set the wall normal
       point_element_pt->wall_unit_normal_fct_pt() = 
        &Global_Physical_Variables::wall_unit_normal_inlet_fct;
       //Set the contact angle (using the strong version of the constraint)
       point_element_pt->set_contact_angle(
        &Global_Physical_Variables::Inlet_Angle);
      }
     
     //Make another point element from the right-hand side of the 
     //last surface element (note that this relies on knowledge of 
     //the element order within the mesh)
     if(e==n_boundary_element-1)
      {
       FluidInterfaceBoundingElement* point_element_pt =
        surface_element_pt->make_bounding_element(1);
       //Add element to the mesh
       Point_mesh_pt->add_element_pt(point_element_pt);
       //Set the capillary number
       point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
       // Set the function that specifies the wall normal
       point_element_pt->wall_unit_normal_fct_pt() = 
        &Global_Physical_Variables::wall_unit_normal_outlet_fct;
      }
    }
  } //end of make_free_surface_elements

 /// Complete the build of the problem setting all standard
 /// parameters and boundary conditions
 void complete_build()
  {
   using namespace Global_Physical_Variables;
   
   //Complete the build of the fluid elements by passing physical parameters
   //Find the number of bulk elements
   unsigned n_element = Bulk_mesh_pt->nelement();
   //Loop over all the fluid elements 
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     
     //Set the Reynolds number
     temp_pt->re_pt() = &Re;
     //The Strouhal number is 1, so ReSt = Re
     temp_pt->re_st_pt() = &Re;
     //Set the Reynolds number / Froude number
     temp_pt->re_invfr_pt() = &ReInvFr;
     //Set the direction of gravity
     temp_pt->g_pt() = &G;
    }
   
   //------------Set the boundary conditions for this problem----------

   {
    //Determine whether we are solving an elastic problem or not
    bool elastic = false;
    if(dynamic_cast<SolidNode*>(Bulk_mesh_pt->node_pt(0))) {elastic=true;}

    //Loop over the bottom of the mesh (the wall of the channel)
    unsigned n_node = Bulk_mesh_pt->nboundary_node(0);
    for(unsigned j=0;j<n_node;j++)
     {
      //Pin the u- and v- velocities
      Bulk_mesh_pt->boundary_node_pt(0,j)->pin(0);
      Bulk_mesh_pt->boundary_node_pt(0,j)->pin(1);

      //If we are formulating the elastic problem pin both positions
      //of nodes
      if(elastic)
       {
        static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(0,j))
         ->pin_position(0);
        static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(0,j))
         ->pin_position(1);
       }
     }
    
    //Loop over the inlet and set the Dirichlet condition
    //of no vertical velocity
    n_node = Bulk_mesh_pt->nboundary_node(3);
    for(unsigned j=0;j<n_node;j++)
     {
      Bulk_mesh_pt->boundary_node_pt(3,j)->pin(1);

      //If elastic pin horizontal position of nodes
      if(elastic)
       { 
        static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(3,j))
         ->pin_position(0);
       }
     }
    
    //Loop over the outlet and set the Dirichlet condition
    //of no vertical velocity
    n_node = Bulk_mesh_pt->nboundary_node(1);
    for(unsigned j=0;j<n_node;j++)
     {
      Bulk_mesh_pt->boundary_node_pt(1,j)->pin(1);

      //If elastic pin horizontal position
      if(elastic)
       { 
        static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(1,j))
         ->pin_position(0);
       }
     }
   }
 
   //Attach the boundary conditions to the mesh
   std::cout << assign_eqn_numbers() << " in the main problem" << std::endl; 
  } //end of complete_build

 /// Generic desructor to clean up the memory allocated in the problem
 ~InclinedPlaneProblem()
  {
   //Clear node storage before the mesh can be deleted, 
   //to avoid double deletion
   this->Point_mesh_pt->flush_node_storage();
   //Then delete the mesh
   delete this->Point_mesh_pt;
   //Clear node storage and then delete mesh
   this->Surface_mesh_pt->flush_node_storage();
   delete this->Surface_mesh_pt;
   //Clear node storage and then delete mesh
   this->Traction_mesh_pt->flush_node_storage();
   delete this->Traction_mesh_pt;
   //Delete the bulk mesh (no need to clear node storage)
   delete this->Bulk_mesh_pt;
   //Delete the time stepper
   delete this->time_stepper_pt();
  }

};


//-------------------------------------------------------------------------
/// Solve the steady problem
//-------------------------------------------------------------------------
template<class ELEMENT,class INTERFACE_ELEMENT>
void InclinedPlaneProblem<ELEMENT,INTERFACE_ELEMENT>::solve_steady()
{
 //Load the namespace
 using namespace Global_Physical_Variables;
 
 //Initially set all nodes to the Nusselt flat-film solution
 {
  unsigned n_node = Bulk_mesh_pt->nnode();
  for(unsigned  n=0;n<n_node;n++)
   {
    double y = Bulk_mesh_pt->node_pt(n)->x(1);
    //Top row
    Bulk_mesh_pt->node_pt(n)->set_value(0,0.5*ReInvFr*sin(Alpha)*(2.0*y - y*y));
   }
 }
 
 //Do one steady solve
 steady_newton_solve();

 //Output the full flow field
 std::string filename = Output_prefix;;
 filename.append("_output.dat");
 ofstream file(filename.c_str());
 Bulk_mesh_pt->output(file,5);
 file.close();
} //end of solve_steady


//----------------------------------------------------------------------
/// Perform n_tsteps timesteps of size dt
//----------------------------------------------------------------------
template<class ELEMENT, class INTERFACE_ELEMENT>
void InclinedPlaneProblem<ELEMENT, INTERFACE_ELEMENT>::
timestep(const double &dt, const unsigned &n_tsteps)
{
 //Need to use the Global variables here
 using namespace Global_Physical_Variables;
 
 //Open an output file
 std::string filename = Output_prefix;
 filename.append("_time_trace.dat");
 ofstream trace(filename.c_str()); 
 //Counter that will be used to output the full flowfield
 //at certain timesteps
 int counter=0; 
 
 //Initial output of the time and the value of the vertical position at the
 //left and right-hand end of the free surface
 trace << time_pt()->time() << " " 
       << Bulk_mesh_pt->boundary_node_pt(2,0)->value(1) 
       << " "
       <<  Bulk_mesh_pt->
  boundary_node_pt(2, Bulk_mesh_pt->nboundary_node(2)-1)->x(1) 
       << " "
       << std::endl;
 
 //Loop over the desired number of timesteps
 for(unsigned t=1;t<=n_tsteps;t++)
  {
   //Increase the counter
   counter++;
   cout << std::endl;
   cout << "--------------TIMESTEP " << t<< " ------------------" << std::endl;
   
   //Take a timestep of size dt
   unsteady_newton_solve(dt);
   
   //Uncomment to get full solution output
   if(counter==2) //Change this number to get output every n steps
    {
     std::ofstream file;
     std::ostringstream filename;
     filename << Output_prefix << "_step" << Re << "_" << t << ".dat";
     file.open(filename.str().c_str());
     Bulk_mesh_pt->output(file,5);
     file.close();
     
     counter=0;
    }

   //Always output the interface
   {
    std::ofstream file;
     std::ostringstream filename;
     filename << Output_prefix << "_interface_" << Re << "_" << t << ".dat";
     file.open(filename.str().c_str());
     Surface_mesh_pt->output(file,5);
     file.close();
   }
   
   //Output the time and value of the vertical position of the free surface
   //at the left- and right-hand ends
   trace << time_pt()->time() << " "
         << Bulk_mesh_pt->boundary_node_pt(2,0)->x(1) << " "
         << 
    Bulk_mesh_pt->
    boundary_node_pt(2,Bulk_mesh_pt->nboundary_node(2)-1)->x(1) << " "
         << std::endl;
  }
} //end of timestep


//====================================================================
// Spine formulation of the problem
//===================================================================


//======================================================================
/// Create a spine mesh for the problem
//======================================================================
template <class ELEMENT>
class SpineInclinedPlaneMesh : 
 public SimpleRectangularQuadMesh<ELEMENT>,
 public SpineMesh
{
public:
 SpineInclinedPlaneMesh(const unsigned &nx, const unsigned &ny,
                        const double &lx, const double &ly,
                        TimeStepper* time_stepper_pt) :
  SimpleRectangularQuadMesh<ELEMENT>
 (nx,ny,lx,ly,time_stepper_pt), SpineMesh()
  {
   //Find the number of linear points in the element
   unsigned n_p =  dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();
   //Reserve storage for the number of spines
   Spine_pt.reserve((n_p-1)*nx + 1);
   
   //Create single pointer to a spine
   Spine* new_spine_pt=0;

   //Now loop over the elements horizontally
   for(unsigned long j=0;j<nx;j++)
    {
     //In most elements, we don't assign a spine to the last column,
     //beacuse that will be done by the next element
     unsigned n_pmax = n_p-1;
     //In the last element, however, we must assign the final spine
     if(j==nx-1) {n_pmax = n_p;}

     //Loop over all nodes horizontally
     for(unsigned l2=0;l2<n_pmax;l2++)
      {
       //Create a new spine with unit height and add to the mesh
       new_spine_pt=new Spine(1.0);
       Spine_pt.push_back(new_spine_pt);
       
       // Get the node
       SpineNode* nod_pt=element_node_pt(j,l2);
       //Set the pointer to spine
       nod_pt->spine_pt() = new_spine_pt;
       //Set the fraction
       nod_pt->fraction() = 0.0;
       // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = this; 
       
       //Loop vertically along the spine
       //Loop over the elements 
       for(unsigned long i=0;i<ny;i++)
        {
         //Loop over the vertical nodes, apart from the first
         for(unsigned l1=1;l1<n_p;l1++)
          {
           // Get the node
           SpineNode* nod_pt=element_node_pt(i*nx+j,l1*n_p+l2);
           //Set the pointer to the spine
           nod_pt->spine_pt() = new_spine_pt;
           //Set the fraction
           nod_pt->fraction()=(double(i)+double(l1)/double(n_p-1))/double(ny);
           // Pointer to the mesh that implements the update fct
           nod_pt->spine_mesh_pt() = this; 
          }  
        }
      }
    } //End of horizontal loop over elements  
  } //end of constructor

 /// General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific node update
 /// actions:  along vertical spines
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of y
   spine_node_pt->x(1) = W*H;
  }
};



//============================================================================
//Specific class for inclined plane problem using spines
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class SpineInclinedPlaneProblem : 
 public InclinedPlaneProblem<ELEMENT,SpineLineFluidInterfaceElement<ELEMENT> >
{
public:

 //Constructor
 SpineInclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length): 
  InclinedPlaneProblem<ELEMENT,SpineLineFluidInterfaceElement<ELEMENT> >
  (nx,ny,length) 
  {
   //Set the name
   this->Output_prefix = "spine";

   //Create our one and only timestepper, with adaptive timestepping
   this->add_time_stepper_pt(new TIMESTEPPER);

   //Create the bulk mesh
   this->Bulk_mesh_pt = new  SpineInclinedPlaneMesh<ELEMENT>(
    nx,ny,length,1.0,this->time_stepper_pt());

   //Create the traction elements
   this->make_traction_elements();
   //Create the free surface elements
   this->make_free_surface_elements();

   //Add all sub meshes to the problem
   this->add_sub_mesh(this->Bulk_mesh_pt);
   this->add_sub_mesh(this->Traction_mesh_pt);
   this->add_sub_mesh(this->Surface_mesh_pt);
   this->add_sub_mesh(this->Point_mesh_pt);
   //Create the global mesh
   this->build_global_mesh();

   //Complete the build of the problem
   this->complete_build();
  }

 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check()
  {this->Bulk_mesh_pt->node_update();}

};


//====================================================================
// Elastic formulation of the problem
//===================================================================


//======================================================================
/// Create an Elastic mesh for the problem
//======================================================================
template <class ELEMENT>
class ElasticInclinedPlaneMesh : 
 public SimpleRectangularQuadMesh<ELEMENT>,
 public SolidMesh
{
 //Public functions
 public:
 ElasticInclinedPlaneMesh(const unsigned &nx, const unsigned &ny,
                          const double &lx, const double &ly,
                          TimeStepper* time_stepper_pt) :
  SimpleRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt), SolidMesh()
  {
   //Make the current configuration the undeformed one
   set_lagrangian_nodal_coordinates();
  }
};




//============================================================================
//Specific class for inclined plane problem using pseudo-elastic formulation
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class ElasticInclinedPlaneProblem : 
 public InclinedPlaneProblem<ELEMENT,ElasticLineFluidInterfaceElement<ELEMENT> >
{
public:
 //Constructor
 ElasticInclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length) :
  InclinedPlaneProblem<ELEMENT,ElasticLineFluidInterfaceElement<ELEMENT> >
  (nx,ny,length) 
  {
   //Set the name
   this->Output_prefix = "elastic";
   
   //Create our one and only timestepper, with adaptive timestepping
   this->add_time_stepper_pt(new TIMESTEPPER);

   //Create the bulk mesh
   this->Bulk_mesh_pt = new  ElasticInclinedPlaneMesh<ELEMENT>(
    nx,ny,length,1.0,this->time_stepper_pt());

   //Set the consititutive law for the elements
   unsigned n_element = this->Bulk_mesh_pt->nelement();
   //Loop over all the fluid elements 
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->element_pt(e));
     //Set the constitutive law 
     temp_pt->constitutive_law_pt() = 
      Global_Physical_Variables::Constitutive_law_pt;
    }
   
   //Create the traction elements
   this->make_traction_elements();
   //Create the free surface element
   this->make_free_surface_elements();

   //Add all sub meshes to the problem
   this->add_sub_mesh(this->Bulk_mesh_pt);
   this->add_sub_mesh(this->Traction_mesh_pt);
   this->add_sub_mesh(this->Surface_mesh_pt);
   this->add_sub_mesh(this->Point_mesh_pt);
   //Create the global mesh
   this->build_global_mesh();

   //Complete the rest of the build
   this->complete_build();
  } //end of constructor

 /// Update Lagrangian positions after each timestep 
 /// (updated-lagrangian approach)
 void actions_after_implicit_timestep()
  {
   //Now loop over all the nodes and reset their Lagrangian coordinates
   unsigned n_node = this->Bulk_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     //Cast node to an elastic node
     SolidNode* temp_pt = 
      static_cast<SolidNode*>(this->Bulk_mesh_pt->node_pt(n));
     for(unsigned j=0;j<2;j++) {temp_pt->xi(j) = temp_pt->x(j);}
    }
  } //end of actions_after_implicit_timestep

};


//start of main
int main(int argc, char **argv)
{
 using namespace Global_Physical_Variables;
 
 //Set the constitutive law for the mesh deformation
 Constitutive_law_pt  = new GeneralisedHookean(&Global_Physical_Variables::Nu);

#ifdef CR_ELEMENT
#define FLUID_ELEMENT QCrouzeixRaviartElement<2>
#else
#define FLUID_ELEMENT QTaylorHoodElement<2>
#endif 

 //Initialise physical parameters
 //Scale Reynolds number to be independent of alpha.
 Re = 4.0/sin(Alpha); 

 //Set the direction of gravity
 G[0] = sin(Alpha);
 G[1] = -cos(Alpha);

 //The wall normal to the inlet is in the negative x direction
 Wall_normal.resize(2);
 Wall_normal[0] = -1.0;
 Wall_normal[1] = 0.0;

 //Spine problem
 {
  //Create the problem
  SpineInclinedPlaneProblem<SpineElement<FLUID_ELEMENT >, BDF<2> > 
   problem(30,4,Global_Physical_Variables::Length);
  
  //Solve the steady problem
  problem.solve_steady();
  
  //Prepare the problem for timestepping 
  //(assume that it's been at the flat-film solution for all previous time)
  double dt = 0.1;
  problem.assign_initial_values_impulsive(dt);
  
  //Timestep it 
  problem.timestep(dt,2);
 } //End of spine problem
 

 //Elastic problem
 {
  //Create the problem
  ElasticInclinedPlaneProblem<
   PseudoSolidNodeUpdateElement<FLUID_ELEMENT,QPVDElement<2,3> >, BDF<2> > 
     problem(30,4,Global_Physical_Variables::Length);
  
  //Solve the steady problem
  problem.solve_steady();
  
  //Prepare the problem for timestepping 
  //(assume that it's been at the flat-film solution for all previous time)
  double dt = 0.1;
  problem.assign_initial_values_impulsive(dt);

  //Timestep it
  problem.timestep(dt,2);
 } //End of elastic problem
}

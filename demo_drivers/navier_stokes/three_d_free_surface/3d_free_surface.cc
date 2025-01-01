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
//Three-dimensional free-surface test case

// C++ includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio> 
#include <complex>
 
// The oomphlib headers   
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"


// The mesh
#include "meshes/single_layer_cubic_spine_mesh.h"

using namespace std;

using namespace oomph;


//======================================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Physical_Variables
{
 
 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Reynolds number
 double Re;

 /// Womersley = Reynolds times Strouhal
 double ReSt; 
 
 /// Product of Reynolds and Froude number
 double ReInvFr;

 /// Capillary number
 double Ca;  

 /// Direction of gravity
 Vector<double> G(3);

}


//======================================================================
/// Single fluid interface problem
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:

 /// Constructor: Pass number of elements in x and y directions. Also lengths
 /// of the domain in x- and y-directions and the height of the layer

 InterfaceProblem(const unsigned &Nx, const unsigned &Ny, const unsigned &Nz,
                  const double &Lx, const double &Ly, const double &h);
 
 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check(){Bulk_mesh_pt->node_update();}

 /// Run an unsteady simulation with specified number of steps
 void unsteady_run(const unsigned& nstep); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 

private:

 /// Trace file
 ofstream Trace_file;

 /// Axial lengths of domain
 double Lx;

 double Ly;

 /// Pointer to bulk mesh
 SingleLayerCubicSpineMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the surface mes
 Mesh* Surface_mesh_pt;

 /// Pointer to a node for documentation purposes
 Node* Document_node_pt;

};


//====================================================================
/// Problem constructor
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::InterfaceProblem
(const unsigned &nx, const unsigned &ny,const unsigned &nz,
 const double &lx, const double &ly, const double& h)
 : Lx(lx), Ly(ly)
{  
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the bulk mesh
 Bulk_mesh_pt = new SingleLayerCubicSpineMesh<ELEMENT>
  (nx,ny,nz,lx,ly,h,time_stepper_pt());

 //Set the documented node
 {
 //Find the number of nodes per row
  unsigned n_node_row = 
   1 + nx*(Bulk_mesh_pt->finite_element_pt(0)->nnode_1d()-1);
  //Find the near central node
  unsigned central_node = (n_node_row + 1)*n_node_row/2;
  Document_node_pt = Bulk_mesh_pt->boundary_node_pt(5,central_node);
 }



 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;

 // Loop over those elements adjacent to the free surface,
 // which we shall choose to be the upper surface
 for(unsigned e1=0;e1<ny;e1++)
  {
   for(unsigned e2=0;e2<nx;e2++)
    {
     // Set a pointer to the bulk element we wish to our interface
     // element to
     FiniteElement* bulk_element_pt =
      Bulk_mesh_pt->finite_element_pt(nx*ny*(nz-1) + e2 + e1*nx);

     // Create the interface element (on face 3 of the bulk element)
     FiniteElement* interface_element_pt =
      new SpineSurfaceFluidInterfaceElement<ELEMENT>(bulk_element_pt,3);

   // Add the interface element to the surface mesh
   this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    }
  }
 
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();
 
 //Pin all nodes on the bottom
 unsigned long n_boundary_node = Bulk_mesh_pt->nboundary_node(0);
 for(unsigned long n=0;n<n_boundary_node;n++)
  {
   for(unsigned i=0;i<3;i++)
    {
     Bulk_mesh_pt->boundary_node_pt(0,n)->pin(i);
    }
  }

 //On the front and back (y=const) pin in y-direction
 for(unsigned b=1;b<4;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
    }
  }

 //On sides pin in x-direction
 for(unsigned b=2;b<5;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
    }
  }

 //Create a Data object whose single value stores the
 //external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Set and pin the external pressure to some random value
 external_pressure_data_pt->set_value(0,1.31);
 external_pressure_data_pt->pin(0);

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements in the layer
 unsigned n_bulk=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_bulk;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Loop over 2D interface elements and set capillary number and 
 //the external pressure
 unsigned long interface_element_pt_range = 
  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<interface_element_pt_range;e++)
  {
   //Cast to a interface element
   SpineSurfaceFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineSurfaceFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   //Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);
  }

 //Do equation numbering
 cout << assign_eqn_numbers() << std::endl; 
 
}

   

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 //Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 // Doc
 Trace_file << time_pt()->time();
 Trace_file << " "  << Document_node_pt->x(2);
 Trace_file << std::endl;


 // Output solution, bulk elements followed by surface elements
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();
 
}

 


//=============================================================================
/// Unsteady run with specified number of steps
//=============================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::unsteady_run(const unsigned& nstep)
{

 // Increase maximum residual
 Problem::Max_residuals=100.0;

 //Distort the mesh
 double epsilon=0.1;
 unsigned Nperiods = 1;
 unsigned Nspine = Bulk_mesh_pt->nspine();
 for(unsigned i=0;i<Nspine;i++)
  {
   double x_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(0);
   double y_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(1);
   
   Bulk_mesh_pt->spine_pt(i)->height() = 
    1.0 - epsilon*(cos(2.0*Nperiods*MathematicalConstants::Pi*x_value/Lx)
                   + cos(2.0*Nperiods*MathematicalConstants::Pi*y_value/Ly));

  }

 //Make sure to update 
 Bulk_mesh_pt->node_update();


 // Doc info object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 doc_info.number()=0;
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 Trace_file << "VARIABLES=\"time\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\"";
 


 //Set value of dt
 double  dt = 0.01;

 //Initialise all time values
 assign_initial_values_impulsive(dt);
  
 //Doc initial solution
 doc_solution(doc_info);

//Loop over the timesteps
 for(unsigned t=1;t<=nstep;t++)
  {
   cout << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   unsteady_newton_solve(dt);

   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info);
   }

}


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//=====================================================================
/// Elastic version of the simple cubic mesh
//=====================================================================


//======================================================================
/// Create an Elastic mesh for the problem
//======================================================================
template <class ELEMENT>
class ElasticSimpleCubicMesh : 
 public SimpleCubicMesh<ELEMENT>,
 public SolidMesh
{
 //Public functions
 public:
 ElasticSimpleCubicMesh(const unsigned &nx, const unsigned &ny, 
                        const unsigned &nz,
                        const double &lx, const double &ly,
                        const double &lz,
                        TimeStepper* time_stepper_pt) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,lx,ly,lz,time_stepper_pt), SolidMesh()
  {
   //Make the current configuration the undeformed one
   set_lagrangian_nodal_coordinates();
  }
};


//======================================================================
/// Elastic version of Single fluid interface problem
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class ElasticInterfaceProblem : public Problem
{
 
public:

 /// Constructor: Pass number of elements in x and y directions. Also lengths
 /// of the domain in x- and y-directions and the height of the layer
 ElasticInterfaceProblem(const unsigned &Nx, const unsigned &Ny, 
                         const unsigned &Nz,
                         const double &Lx, const double &Ly, const double &h);
 
 /// Run an unsteady simulation with specified number of steps
 void unsteady_run(const unsigned& nstep); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Trace file
 ofstream Trace_file;

 /// Axial lengths of domain
 double Lx;

 double Ly;

 /// Pointer to bulk mesh
 ElasticSimpleCubicMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the surface mes
 Mesh* Surface_mesh_pt;

 /// Pointer to a node for documentation purposes
 Node* Document_node_pt;
 
};


//====================================================================
/// Problem constructor
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
ElasticInterfaceProblem<ELEMENT,TIMESTEPPER>::ElasticInterfaceProblem
(const unsigned &nx, const unsigned &ny,const unsigned &nz,
 const double &lx, const double &ly, const double& h)
 : Lx(lx), Ly(ly)
{  
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the bulk mesh
 Bulk_mesh_pt = new ElasticSimpleCubicMesh<ELEMENT>
  (nx,ny,nz,lx,ly,h,time_stepper_pt());

 //Set the documented node
 {
 //Find the number of nodes per row
  unsigned n_node_row = 
   1 + nx*(Bulk_mesh_pt->finite_element_pt(0)->nnode_1d()-1);
  //Find the near central node
  unsigned central_node = (n_node_row + 1)*n_node_row/2;
  Document_node_pt = Bulk_mesh_pt->boundary_node_pt(5,central_node);
 }

 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
 
 // Loop over those elements adjacent to the free surface,
 // which we shall choose to be the upper surface
 for(unsigned e1=0;e1<ny;e1++)
  {
   for(unsigned e2=0;e2<nx;e2++)
    {
     // Set a pointer to the bulk element we wish to our interface
     // element to
     FiniteElement* bulk_element_pt =
      Bulk_mesh_pt->finite_element_pt(nx*ny*(nz-1) + e2 + e1*nx);

     // Create the interface element (on face 3 of the bulk element)
     FiniteElement* interface_element_pt =
      new ElasticSurfaceFluidInterfaceElement<ELEMENT>(bulk_element_pt,3);
     
     // Add the interface element to the surface mesh
     this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    }
  }
 
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 
 // Combine all sub-meshes into a single mesh
 build_global_mesh();
 
 //Pin all nodes on the bottom
 unsigned long n_boundary_node = Bulk_mesh_pt->nboundary_node(0);
 for(unsigned long n=0;n<n_boundary_node;n++)
  {
   SolidNode* nod_pt = 
    static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(0,n));
   for(unsigned i=0;i<3;i++)
    {
     nod_pt->pin(i);
     nod_pt->pin_position(i);
    }
  }
 
 //On the front and back (y=const) pin in y-direction
 for(unsigned b=1;b<4;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     SolidNode* nod_pt = 
      static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n));
     nod_pt->pin(1);
     nod_pt->pin_position(1);
    }
  }

 //On sides pin in x-direction
 for(unsigned b=2;b<5;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     SolidNode* nod_pt = 
      static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n));
     nod_pt->pin(0);
     nod_pt->pin_position(0);
    }
  }


 //Only allow the nodes to move in the vertical direction
 unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   SolidNode* nod_pt = static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n));
   nod_pt->pin_position(0);
   nod_pt->pin_position(1);
  }

 //Create a Data object whose single value stores the
 //external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Set and pin the external pressure to some random value
 external_pressure_data_pt->set_value(0,1.31);
 external_pressure_data_pt->pin(0);
 
 //Complete the problem setup to make the elements fully functional

 //Loop over the elements in the layer
 unsigned n_bulk=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_bulk;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;

   el_pt->constitutive_law_pt() = 
    Global_Physical_Variables::Constitutive_law_pt;
  }

 //Loop over 2D interface elements and set capillary number and 
 //the external pressure
 unsigned long interface_element_pt_range = 
  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<interface_element_pt_range;e++)
  {
   //Cast to a interface element
   ElasticSurfaceFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<ElasticSurfaceFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   //Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);
  }

 //Do equation numbering
 cout << assign_eqn_numbers() << std::endl; 
 
}

   

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void ElasticInterfaceProblem<ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo& doc_info)
{ 
 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 //Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 //Estimate the node in the middle of the domain

 // Doc
 Trace_file << time_pt()->time();
 Trace_file << " "  << Document_node_pt->x(2);
 Trace_file << std::endl;
 

 // Output solution, bulk elements followed by surface elements
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();
 
}

 


//=============================================================================
/// Unsteady run with specified number of steps
//=============================================================================
template<class ELEMENT, class TIMESTEPPER>
void ElasticInterfaceProblem<ELEMENT,TIMESTEPPER>::unsteady_run(const unsigned& nstep)
{

 // Increase maximum residual
 Problem::Max_residuals=100.0;

 //Distort the mesh
 double epsilon=0.1;
 unsigned Nperiods = 1;
 
 //Loop over all nodes
 unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   Node* nod_pt = Bulk_mesh_pt->node_pt(n);
   //Get the x y and z values
   double x_value = nod_pt->x(0);
   double y_value = nod_pt->x(1);
   double z_value = nod_pt->x(2);

   //Calculate the new scaled z- position
   double scaled_z = z_value*
    (1.0 - epsilon*(cos(2.0*Nperiods*MathematicalConstants::Pi*x_value/Lx)
                    + cos(2.0*Nperiods*MathematicalConstants::Pi*y_value/Ly)));
   
   //Set the new position
   nod_pt->x(2) = scaled_z;
  }

 //Reset the lagrangian positions
 Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

 // Doc info object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT_elastic");
 doc_info.number()=0;
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 Trace_file << "VARIABLES=\"time\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\"";

 //Set value of dt
 double  dt = 0.01;

 //Initialise all time values
 assign_initial_values_impulsive(dt);
  
 //Doc initial solution
 doc_solution(doc_info);

//Loop over the timesteps
 for(unsigned t=1;t<=nstep;t++)
  {
   cout << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   unsteady_newton_solve(dt);

   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info);
   }

}



//======================================================================
/// Driver code for unsteady two-layer fluid problem. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only two steps.

// In my version we will change nsteps in the programs
//======================================================================
int main(int argc, char *argv[]) 
{

 // Set physical parameters:

 //Set direction of gravity: Vertically downwards
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = 0.0;
 Global_Physical_Variables::G[2] = -1.0;
 
 // Womersley number = Reynolds number (St = 1)
 Global_Physical_Variables::ReSt = 5.0;
 Global_Physical_Variables::Re = Global_Physical_Variables::ReSt;
 
 // The Capillary number
 Global_Physical_Variables::Ca = 0.01;
 
 // Re/Fr 
 Global_Physical_Variables::ReInvFr = 0.0;
 
 //Axial lngth of domain
 double lx = 2.0;
 
 double ly = 2.0;
   
 // Number of elements in x- and y-direction
 unsigned nx=5;
 
 unsigned ny = 5;
 
 // Number of elements in layer
 unsigned nz=7;
 
 // Height of layer
 double h=1.0;

 {
  //Set up the spine test problem
  InterfaceProblem<SpineElement<QCrouzeixRaviartElement<3> >,BDF<2> >
   problem(nx,ny,nz,lx,ly,h);
  
  // Number of steps: 
  unsigned nstep;
  if(argc > 1)
   {
    // Validation run: Just two steps
    nstep=2;
   }
  else
   {
    // Full run otherwise
    nstep=160;
   }
  
  // Run the unsteady simulation
  problem.unsteady_run(nstep);
 }

 //Set the constituive law
   Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);

 {
  //Set up the spine test problem
  ElasticInterfaceProblem<PseudoSolidNodeUpdateElement
                          <QCrouzeixRaviartElement<3>,
                           QPVDElementWithPressure<3> >,BDF<2> >
   problem(nx,ny,nz,lx,ly,h);
  
  // Number of steps: 
  unsigned nstep;
  if(argc > 1)
   {
    // Validation run: Just two steps
    nstep=2;
   }
  else
   {
    // Full run otherwise
    nstep=160;
   }
  
  // Run the unsteady simulation
  problem.unsteady_run(nstep);
 }


}


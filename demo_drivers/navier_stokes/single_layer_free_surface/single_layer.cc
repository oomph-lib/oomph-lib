//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
// Two-dimensional single-layer fluid problem 
 
// The oomphlib headers   
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/single_layer_spine_mesh.h"

using namespace std;

using namespace oomph;


//======================================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re;

 /// Womersley = Reynolds times Strouhal
 double ReSt; 
 
 /// Product of Reynolds and Froude number
 double ReInvFr;

 /// Capillary number
 double Ca;  

 /// Direction of gravity
 Vector<double> G(2);

}



//======================================================================
/// Single fluid interface problem
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{


public:

 /// Constructor: Pass number of elements in x and y directions. Also length
 /// of the domain in x-directions and the height of the layer
 /// Bool flag to indicate symmetry in x.
 InterfaceProblem(const unsigned &Nx, const unsigned &Ny, 
                  const double &Lx, const double &h,
                  const bool& symmetric_in_x);
 
 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check(){Bulk_mesh_pt->node_update();}

 //Update before solve is empty
 void actions_before_newton_solve() {}

 /// \short Update after solve can remain empty, because the update 
 /// is performed automatically after every Newton step.
 void actions_after_newton_solve() {}

 ///Fix pressure value l in element e to value p_value
 void fix_pressure(const unsigned &e, const unsigned &l, 
                   const double &pvalue)
  {
   //Fix the pressure at that element
   dynamic_cast<ELEMENT *>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(l,pvalue);
  }
 

 /// Run an unsteady simulation with specified number of steps
 void unsteady_run(const unsigned& nstep); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 

private:

 /// Trace file
 ofstream Trace_file;

 /// Axial length of domain
 double Lx;

 /// Pointer to bulk mesh
 SingleLayerSpineMesh<ELEMENT,SpineLineFluidInterfaceElement<ELEMENT> >*
 Bulk_mesh_pt;

 /// Is the domain symmetric in the x-direction? 
 bool Symmetric_in_x;

};


//====================================================================
/// Problem constructor
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::InterfaceProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& h, const bool& symmetric_in_x) 
 : Lx(lx), Symmetric_in_x(symmetric_in_x)
{  

 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the mesh
 Bulk_mesh_pt = new SingleLayerSpineMesh<ELEMENT,
  SpineLineFluidInterfaceElement<ELEMENT> >(nx,ny,lx,h,symmetric_in_x,
                                                   time_stepper_pt());
 // Make bulk mesh the global mesh
 mesh_pt()=Bulk_mesh_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 unsigned num_bound=mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Boundary conditions depend on case considered:
     // Symmetric in x:
     if (Symmetric_in_x)
      {
       // Pin x-velocity: No slip/no penetration on all boundaries
       // apart from top
       if (ibound!=2)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
        }
       // Pin y-velocity only on bottom wall
       if (ibound==0)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
        }
      }
     // Not symmetric in x: Parallel inflow at bottom; symmetry
     // at x=0; no slip at x=1.
     else
      {
       // Pin x-velocity: Parallel inflow; no penetration through the
       // side walls
       if (ibound!=2)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
        }
       //Pin y-velocity only on outer wall: boundary 1
       if (ibound==1)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
        }
      }
    }
  }

 //Create a Data object whose single value stores the
 //external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Set and pin the external pressure to some random value
 external_pressure_data_pt->set_value(0,1.31);
 external_pressure_data_pt->pin(0);

 
 // Left boundary node pinned?
 if (!symmetric_in_x)
  {
   // Fix right-most spine height
   unsigned n_last_spine=Bulk_mesh_pt->nspine()-1;
   Bulk_mesh_pt->spine_pt(n_last_spine)->spine_height_pt()->pin(0);
  }


 //Complete the problem setup to make the elements fully functional

 //Loop over the elements in the layer
 unsigned n_bulk=Bulk_mesh_pt->nbulk();
 for(unsigned i=0;i<n_bulk;i++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           bulk_element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;

   //Assign the time pointer
   el_pt->time_pt() = time_pt();
  }

 //Loop over 1D interface elements and set capillary number and 
 //the external pressure
 unsigned interface_element_pt_range = Bulk_mesh_pt->ninterface_element();
 for(unsigned i=0;i<interface_element_pt_range;i++)
  {
   //Cast to a interface element
   SpineLineFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>
    (Bulk_mesh_pt->interface_element_pt(i));

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

 // Number of interface elements
 unsigned ninterface=Bulk_mesh_pt->ninterface_element();

 // Number of spines
 unsigned nspine=Bulk_mesh_pt->nspine();

 // Doc
 Trace_file << time_pt()->time();
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(0)->height();
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(nspine-1)->height();
 Trace_file << " "  
            << dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>(
               Bulk_mesh_pt->interface_element_pt(0))->
               actual_contact_angle_left()*
               180.0/MathematicalConstants::Pi << " " ;
 Trace_file << " "  
            << dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>(
               Bulk_mesh_pt->interface_element_pt(ninterface-1))->
               actual_contact_angle_right()*180.0/MathematicalConstants::Pi 
            << " ";
 Trace_file << std::endl;


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
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
   if (Symmetric_in_x)
    {
     Bulk_mesh_pt->spine_pt(i)->height() = 
      1.2 + epsilon*(cos(2.0*Nperiods*MathematicalConstants::Pi*x_value/Lx));
    }
   else
    {
     Bulk_mesh_pt->spine_pt(i)->height() = 
      1.0 + epsilon*(4.0*pow(x_value/Lx,2)+
       cos(2.0*Nperiods*MathematicalConstants::Pi*x_value/Lx));
    }
  }

 //Make sure to update 
 Bulk_mesh_pt->node_update();


 // Doc info object
 DocInfo doc_info;
 if (Symmetric_in_x)
  {
   // Set output directory
   doc_info.set_directory("RESLT_sym");
   doc_info.number()=0;
  }
 else
  {
   // Set output directory
   doc_info.set_directory("RESLT_pinned");
   doc_info.number()=0;
  }
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 Trace_file << "VARIABLES=\"time\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>left</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>right</sub>\"";


 //Set value of dt
 double dt;
 if (Symmetric_in_x)
  {
   dt = 0.01;
  }
 else
  {
   dt = 0.01;
  }

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
//======================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set physical parameters:

 //Set direction of gravity: Vertically downwards
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = -1.0;

 // Womersley number = Reynolds number (St = 1)
 Global_Physical_Variables::ReSt = 5.0;
 Global_Physical_Variables::Re = Global_Physical_Variables::ReSt;

 // The Capillary number
 Global_Physical_Variables::Ca = 0.01;

 // Re/Fr hierher check what Fr is 
 Global_Physical_Variables::ReInvFr = 0.0;

 // Apply symmetry boundary conditions in x-direction?
 bool symmetric_in_x=true; //Overwritten during validation run
 
 // During validation: Run both cases
 unsigned ncase=1;
 if (CommandLineArgs::Argc>1) ncase=2;

 // Loop over problems with and without symmetry boundary conditions
 for (unsigned i=0;i<ncase;i++)
  {
   
   // During validation: Run both cases
   if (CommandLineArgs::Argc>1)
    {
     if (i==0)
      {
       symmetric_in_x=true;
      }
     else
      {
       symmetric_in_x=false;
      }
    }
   
   //Axial lngth of domain
   double lx = 2.0;
   
   // Number of elements in x-direction
   unsigned nx=10;
   
   // Number of elements in layer
   unsigned ny=7;
   
   // Height of layer
   double h=2.0;
   
   //Set up the spine test problem
   InterfaceProblem<SpineElement<QCrouzeixRaviartElement<2> >,BDF<2> >
    problem(nx,ny,lx,h,symmetric_in_x);
   
   // Number of steps: 
   unsigned nstep;
   if (CommandLineArgs::Argc>1)
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


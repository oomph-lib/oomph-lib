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
//Three-dimensional free-surface test case including insoluble
//surfactant transport. This is a 3D implementation of the
//axisymmetric Rayleigh--Plateau problem solved in 
//rayleigh_instability_insoluble_surfactant.cc, which 
//should be a hard test
//because it's implemented in Cartesian coordinates. 
//The main aim is to test the implementation of the surface transport
//equations in 3D (2D surface).
 
// The oomphlib headers   
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

// The basic mesh
#include "meshes/single_layer_cubic_spine_mesh.h"

using namespace std;
using namespace oomph;


//==start_of_namespace==================================================
/// Namepspace for global parameters, chosen from Campana et al. as
/// in the axisymmetric problem.
//======================================================================
namespace Global_Physical_Variables
{
 //Film thickness parameter
 double Film_Thickness = 0.2;

 /// Reynolds number
 double Re = 40.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt = Re; 
 
 /// Product of Reynolds and Froude number
 double ReInvFr = 0.0;

 /// Capillary number
 double Ca =  pow(Film_Thickness,3.0);

 /// Direction of gravity
 Vector<double> G(3);

// Wavelength of the domain
 double Alpha = 1.047;
 
 /// Free surface cosine deformation parameter
 double Epsilon = 1.0e-3;

 /// Surface Elasticity number (weak case)
 double Beta = 3.6e-3;

 /// Surface Peclet number
 double Peclet_S = 4032.0;

 /// Sufrace Peclet number multiplied by Strouhal number
 double Peclet_St_S = 1.0; 
 
} // End of namespace



//=======================================================================
/// Function-type-object to perform comparison of elements in y-direction
//=======================================================================
class ElementCmp
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x, GeneralisedElement* const &y) 
  const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);

   if((cast_x ==0) || (cast_y==0)) {return 0;}
   else
    {return (cast_x->node_pt(0)->x(1)< cast_y->node_pt(0)->x(1));}
             
  }
};


namespace oomph
{
//===================================================================
/// Deform the existing cubic spine mesh into a annular section
/// with spines directed radially inwards from the wall
//===================================================================

template<class ELEMENT>
class AnnularSpineMesh : public SingleLayerCubicSpineMesh<ELEMENT>
{

public:
  AnnularSpineMesh(const unsigned &n_r, const unsigned &n_y, 
                   const unsigned &n_theta, 
		   const double &r_min, const double &r_max, 
                   const double &l_y,
		   const double &theta_min, const double &theta_max, 
		   TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper):
   //This will make a cubic mesh with n_r in the x-direction
   //n_y in the y-direction and n_theta in the z-direction
   //The coordinates will run from 0 to r_max, 0 to l_y and 0 to theta_max
   SingleLayerCubicSpineMesh<ELEMENT>(n_theta,n_y,n_r,r_max,l_y,theta_max,
                                      time_stepper_pt)
  {  
   
    //Find out how many nodes there are
    unsigned n_node = this->nnode();

    //Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
      {
	//pointer to node
	Node* nod_pt=this->node_pt(n);
	SpineNode* spine_node_pt=dynamic_cast<SpineNode*>(nod_pt);
	//Get x/y/z coordinates
	double x_old=nod_pt->x(0);
	double y_old=nod_pt->x(1);
	double z_old=nod_pt->x(2);

	//Mapping
	
	double r = r_min + (r_max-r_min)*z_old;
        double theta = (theta_min + (theta_max-theta_min)*x_old);
	double y = y_old;
	
	
	if(spine_node_pt->spine_pt()->ngeom_parameter()==0)
	  {
           spine_node_pt->h() =
            Global_Physical_Variables::Film_Thickness; 
	    spine_node_pt->spine_pt()->add_geom_parameter(theta);
	  }
	
	//cout << spine_node_pt->spine_pt()->ngeom_parameter()  << std::endl;
	//Set new nodal coordinates
	nod_pt->x(0)=r*cos(theta);
	nod_pt->x(2)=r*sin(theta);
	nod_pt->x(1)=y;
	//cout << "one" << theta << std::endl;	
      }

  }


    virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get theta along the spine
   double theta = spine_node_pt->spine_pt()->geom_parameter(0);
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of y
   spine_node_pt->x(0) = (1.0-W*H)*cos(theta);
   spine_node_pt->x(2) = (1.0-W*H)*sin(theta);
  }
};
}

//======================================================================
/// Single fluid interface problem including transport of an
/// insoluble surfactant.
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:

 /// Constructor: Pass number of elements in x and y directions. Also lengths
 /// of the domain in x- and y-directions and the height of the layer

 InterfaceProblem(const unsigned &n_r, const unsigned &n_y, 
                  const unsigned &n_theta, const double &r_min,
                  const double &r_max, const double &l_y,  
                  const double &theta_max);
 
 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check(){Bulk_mesh_pt->node_update();}

 /// Run an unsteady simulation with specified number of steps
 void unsteady_run(const unsigned& nstep); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Compute the total mass
 double compute_total_mass()
  {
   double mass = 0.0;
   
   // Determine number of 1D interface elements in mesh
   const unsigned n_interface_element = Surface_mesh_pt->nelement();
   
   // Loop over the interface elements
   for(unsigned e=0;e<n_interface_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     SpineSurfaceSurfactantTransportInterfaceElement<ELEMENT>* el_pt = 
      dynamic_cast<SpineSurfaceSurfactantTransportInterfaceElement<ELEMENT>*>
      (Surface_mesh_pt->element_pt(e));

     mass += el_pt->integrate_c();
    }
   return mass;
  }


 
 private:

 /// Trace file
 ofstream Trace_file;

 /// Axial lengths of domain
 double R_max;

 double L_y;

 /// Pointer to bulk mesh
 AnnularSpineMesh<ELEMENT>* Bulk_mesh_pt;

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
(const unsigned &n_r, const unsigned &n_y,const unsigned &n_theta,
 const double &r_min,
 const double &r_max, const double &l_y, const double &theta_max)
 : R_max(r_max), L_y(l_y)
{  
 //this->linear_solver_pt() = new HSL_MA42;

 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor0() = 3.0;
 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor1() = 3.0;
 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor2() = 3.0;
 
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the bulk mesh -- this should be your new annular mesh
 Bulk_mesh_pt = new AnnularSpineMesh<ELEMENT>
   (n_r,n_y,n_theta,r_min,r_max,l_y,0.0,theta_max,time_stepper_pt());

 //Set the documented node
 Document_node_pt = Bulk_mesh_pt->boundary_node_pt(5,0);


 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
   
 // Loop over those elements adjacent to the free surface,
 // which we shall choose to be the upper surface
 for(unsigned e1=0;e1<n_y;e1++)
  {
   for(unsigned e2=0;e2<n_theta;e2++)
    {
     // Set a pointer to the bulk element we wish to our interface
     // element to
     FiniteElement* bulk_element_pt =
      Bulk_mesh_pt->finite_element_pt(n_theta*n_y*(n_r-1) + e2 + e1*n_theta);

     // Create the interface element (on face 3 of the bulk element)
     FiniteElement* interface_element_pt =
      new SpineSurfaceSurfactantTransportInterfaceElement<ELEMENT>(bulk_element_pt,3);

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

 //On sides pin in z-direction
 for(unsigned b=4;b<5;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
    }
  }

 // pinning the top surface
 for(unsigned b=2;b<3;b++)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
    }
  }

 //Loop over the upper surface and set the surface concentration
 {
  unsigned b=5;
  n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
  for(unsigned long n=0;n<n_boundary_node;n++) 
   {
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(3,1.0);;
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
   SpineSurfaceSurfactantTransportInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineSurfaceSurfactantTransportInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   //Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;
   
   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Peclet_S;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Peclet_St_S;

  }

 //Do equation numbering
 cout << assign_eqn_numbers() << std::endl; 

 //Now sort the elements to minimise frontwidth
 std::sort(mesh_pt()->element_pt().begin(),
           mesh_pt()->element_pt().end(),
           ElementCmp());

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
 Trace_file << " "  << Document_node_pt->x(0) << " "
            << this->compute_total_mass()
            << std::endl;


 // Output solution, bulk elements followed by surface elements
 /*sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();*/
 //Change the file name and output surface in separate file
 sprintf(filename,"%s/surface%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
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
 Problem::Max_residuals=500.0;

 //Distort the mesh
 const double epsilon=Global_Physical_Variables::Epsilon;
 unsigned Nspine = Bulk_mesh_pt->nspine();
 for(unsigned i=0;i<Nspine;i++)
  {
    double y_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(1);
    
    Bulk_mesh_pt->spine_pt(i)->height() =
     Global_Physical_Variables::Film_Thickness*(1.0 +  
                                                + epsilon*cos(Global_Physical_Variables::Alpha*y_value));
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
 Trace_file << "\n";

 //Set value of dt
 double  dt = 0.1;

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


//==start_of_main========================================================
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
 Global_Physical_Variables::G[2] = 0.0;

 //Set the film thickness
 const double h = Global_Physical_Variables::Film_Thickness;
 
 //Axial lngth of domain
 double r_max = 1.0;
 double r_min = r_max - h;

 const double pi = MathematicalConstants::Pi;

 double alpha = Global_Physical_Variables::Alpha;

 double l_y = pi/alpha;
 
 double theta_max = 0.5*pi;
  
 // Number of elements in r and theta direction
 unsigned n_r = 4;
 
 unsigned n_theta = 4;
 
 // Number of elements in axial direction
 unsigned n_y = 20;
 

 {
  //Set up the spine test problem
  //The minimum values are all assumed to be zero.
  InterfaceProblem<SpineElement<QCrouzeixRaviartElement<3> >,BDF<2> >
   problem(n_r,n_y,n_theta,r_min,r_max,l_y,theta_max);
  
  // Number of steps: 
  unsigned nstep;
  if(argc > 1)
   {
    // Validation run: Just five steps
    nstep=5;
   }
  else
   {
    // Full run otherwise
    nstep=1000;
   }
  
  // Run the unsteady simulation
  problem.unsteady_run(nstep);
 }
} // End of main



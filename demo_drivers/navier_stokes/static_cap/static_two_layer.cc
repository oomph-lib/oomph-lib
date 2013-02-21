//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Driver code for a 2D interface hydrostatics problem.
// The system consists of equal volumes of two fluids with 
// the same physical properties, but a non-zero surface tension at 
// their interface, in a domain of height 2 and half-width 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// circular arcs and the pressure jump across the interface should be
// cos(alpha)/0.5 = 2 cos(alpha)/Ca.

//OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/two_layer_spine_mesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables
{

 ///Direction of the wall normal vector
 Vector<double> Wall_normal;

 /// \short Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }

}

//============================================================================
///A Problem class that solves the Navier--Stokes equations
///in an 2D geometry
//============================================================================
template<class ELEMENT>
class CapProblem : public Problem
{

public:

 //Constructor:
 //Nx: Number of elements in the x (horizontal) direction
 //Nh1: Number of elements in the y (vertical) direction in the lower fluid
 //Nh2: Number of elements in the y (vertical) direction in the upper fluid
 CapProblem(const unsigned &Nx, const unsigned &Nh1, 
            const unsigned &Nh2);

 /// Peform a parameter study: Solve problem for a range of contact angles
 void parameter_study();

 /// Finish full specification of the elements
 void finish_problem_setup();

 /// Pointer to the mesh
 TwoLayerSpineMesh<SpineElement<ELEMENT> >* Bulk_mesh_pt;

 /// Pointer to the surface mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to the point mesh
 Mesh* Point_mesh_pt;

 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {Bulk_mesh_pt->node_update();}

 /// Create the volume constraint elements
 void create_volume_constraint_elements();

 /// No other actions required after solve step
 void actions_after_newton_solve() {}

 /// No other actions required before after solve step
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// The Capillary number 
 double Ca;

 /// The volume of the fluid
 double Volume;

 /// The contact angle
 double Angle;

  /// The volume constraint mesh 
 Mesh* Volume_constraint_mesh_pt;

 /// Trace file
 ofstream Trace_file;

 /// Data that is traded for the volume constraint
 Data* Traded_pressure_data_pt;

};



//======================================================================
/// Constructor 
//======================================================================
template<class ELEMENT>
CapProblem<ELEMENT>::CapProblem
(const unsigned &Nx, const unsigned &Nh1, const unsigned &Nh2) :
 Ca(2.1),  //Initialise value of Ca to some random value
 Volume(0.5),  //Initialise the value of the volume 
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = 1.0;
 Global_Physical_Variables::Wall_normal[1] = 0.0;

 //Construct our mesh
 Bulk_mesh_pt = new TwoLayerSpineMesh<SpineElement<ELEMENT> >
 (Nx,Nh1,Nh2,0.5,1.0,1.0);

 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 unsigned n_interface_lower = Bulk_mesh_pt->ninterface_lower();
 for(unsigned i=0;i<n_interface_lower;i++)
  {
   //Construct a new 1D line element adjacent to the interface
   FiniteElement *interface_element_pt =  new 
    SpineLineFluidInterfaceElement<SpineElement<ELEMENT> >(
     Bulk_mesh_pt->interface_lower_boundary_element_pt(i),
     Bulk_mesh_pt->interface_lower_face_index_at_boundary(i));
 
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
   }

//Create the point mesh
Point_mesh_pt = new Mesh;
 {
  //Make the point (contact) element from the last surface element
  FiniteElement* point_element_pt = 
  dynamic_cast<SpineLineFluidInterfaceElement<
   SpineElement<ELEMENT> >*>(Surface_mesh_pt->element_pt(n_interface_lower-1))
   ->make_bounding_element(1);
 
 //Add it to the mesh
 this->Point_mesh_pt->add_element_pt(point_element_pt);
}

//std::ofstream mesh_file("boundaries.dat");
//this->Bulk_mesh_pt->output_boundaries(mesh_file);
// mesh_file.close();

 //Overwrite the linear solver (frontal solver) 
 //Problem::linear_solver_pt() = new HSL_MA42;

 //Hijack one of the pressure values in the upper fluid. Its value
 //will affect the residual of that element but it will not
 //be determined by it!
 Traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
  Bulk_mesh_pt->upper_layer_element_pt(0))->hijack_internal_value(0,0);

 // Loop over the elements on the interface to pass pointer to Ca and
 // the pointer to the Data item that contains the single
 // (pressure) value that is "traded" for the volume constraint 
 unsigned n_interface = Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_interface;e++)
  {
   //Cast to a 1D element
   SpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*el_pt=
   dynamic_cast<SpineLineFluidInterfaceElement<
    SpineElement<ELEMENT> >*>
    (Surface_mesh_pt->element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;
  }
 

 //Set the boundary conditions

 //Pin the velocities on all external boundaries apart from the symmetry
 //line (boundaries 4 and 5) where only the horizontal velocity is pinned
 for (unsigned b=0;b<6;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
     if ((b!=4) && (b!=5))
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }
    }
  }


 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
// dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
 // (Bulk_mesh_pt->interface_element_pt(n_interface-1))->
 // set_contact_angle_right(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->set_contact_angle(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->ca_pt() = &Ca;

 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->wall_unit_normal_fct_pt() = 
  &Global_Physical_Variables::wall_unit_normal_fct;
 
 // Pin a single pressure value: Set the pressure dof 0 in element 0
 // to zero.
 dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->fix_pressure(0,0.0);

 this->create_volume_constraint_elements();

 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Surface_mesh_pt);
 this->add_sub_mesh(Point_mesh_pt);
 this->add_sub_mesh(Volume_constraint_mesh_pt);
 
 this->build_global_mesh();

 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
}



//======================================================================
/// Create the volume constraint elements
//========================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::create_volume_constraint_elements()
{
 //The single volume constraint element
 Volume_constraint_mesh_pt = new Mesh;
 VolumeConstraintElement* vol_constraint_element = 
  new VolumeConstraintElement(&Volume,Traded_pressure_data_pt,0);
 Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);

 //Loop over lower boundaries (or just 1, why?)
 unsigned lower_boundaries[4]={0,1,5};
 for(unsigned i=0;i<3;i++)
  {
   //Read out the actual boundaries
   unsigned b = lower_boundaries[i];

   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));

     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     SpineLineVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new SpineLineVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }

 //Loop over one side of the interface
  {
   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->ninterface_lower();
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->interface_lower_boundary_element_pt(e));

     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->interface_lower_face_index_at_boundary(e);
     
     // Create new element
     SpineLineVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new SpineLineVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }

}



//======================================================================
/// Perform a parameter study
//======================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::parameter_study()
{

 // Create DocInfo object (allows checking if output directory exists)
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 doc_info.number()=0;


 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 Trace_file << "VARIABLES=\"<greek>a</greek><sub>prescribed</sub>\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>left</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>right</sub>\",";
 Trace_file << "\"p<sub>upper</sub>\",";
 Trace_file << "\"p<sub>lower</sub>\",";
 Trace_file << "\"p<sub>exact</sub>\"";
 Trace_file << std::endl;


 // Doc initial state
 doc_solution(doc_info);

 // Bump up counter
 doc_info.number()++;

 //Gradually increase the contact angle
 for(unsigned i=0;i<6;i++)
  {
   //Solve the problem
   steady_newton_solve();

   //Output result
   doc_solution(doc_info);

   // Bump up counter
   doc_info.number()++;

   //Decrease the contact angle
   Angle -= 5.0*MathematicalConstants::Pi/180.0;
  }

} 




//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 //Output domain
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();


 // Number of interface elements
// unsigned ninterface=Bulk_mesh_pt->ninterface_element();

 // Number of spines
 unsigned nspine=Bulk_mesh_pt->nspine();

 // Doc
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(0)->height();
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(nspine-1)->height();
 Trace_file << " "  
 //            << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                Bulk_mesh_pt->interface_element_pt(0))->
//                actual_contact_angle_left()*
//                180.0/MathematicalConstants::Pi << " " ;
//  Trace_file << " "  
//             << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                Bulk_mesh_pt->interface_element_pt(ninterface-1))->
//                actual_contact_angle_right()*180.0/MathematicalConstants::Pi 
            << " ";
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->upper_layer_element_pt(0))->p_nst(0);
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->lower_layer_element_pt(0))->p_nst(0);
 Trace_file << " " << 2.0*cos(Angle)/Ca;
 Trace_file << std::endl;

}

 


//======================================================================
///Main driver: Build problem and initiate parameter study
//======================================================================
int main()
{
 //Construct the problem with 4 x (4+4) elements
 CapProblem<Hijacked<QCrouzeixRaviartElement<2> >  > problem(4,4,4);

 //Solve the problem :)
 problem.parameter_study();
}









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

// The mesh
#include "meshes/two_layer_spine_mesh.h"
//Include our special fixed volume interface elements
#include "fix_vol_int_elements.h"

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
/// Specific mesh for the static cap problem: Essentially a 
/// standard two-layer mesh with an additional element that applies
/// the volume constraint.
//============================================================================
template <class ELEMENT>
class StaticCapMesh : public TwoLayerSpineMesh<SpineElement<ELEMENT>, 
 FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> > >
{

public:

 /// Constructor: Pass number of elements in axial direction and number
 /// of elements in the two layers.
 StaticCapMesh(const unsigned &nx, const unsigned &nh1,
               const unsigned &nh2);
 
 /// Return pointer to the volumetric constraint element
 FiniteElement* &point_element_pt() {return Point_element_pt;}


private:

 /// \short Pointer to the point element that is used to enforce
 /// conservation of mass
 FiniteElement* Point_element_pt;
 

};


//======================================================================
/// Constructor: Pass number of elements in horizontal direction, number
/// of elements in the two layers.
//======================================================================
template<class ELEMENT>
StaticCapMesh<ELEMENT>::StaticCapMesh(const unsigned &nx,
                                            const unsigned &nh1,
                                            const unsigned &nh2) :
 TwoLayerSpineMesh<SpineElement<ELEMENT>, 
  FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> > >
  (nx,nh1,nh2,0.5,1.0,1.0)
                                               
{
 //Reorder the elements
 this-> element_reorder();

 // Last interface element:
 FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >* el_pt=
 dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
  SpineElement<ELEMENT> >*>(this->Interface_element_pt[this->Nx-1]);

 //Make an edge element
 Point_element_pt = el_pt->make_bounding_element(1);

 //Push it back onto the stack
 this->Element_pt.push_back(Point_element_pt);
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

 /// Overload access function for the mesh
 StaticCapMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<StaticCapMesh<ELEMENT>*>(Problem::mesh_pt());}

 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {mesh_pt()->node_update();}

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
 
 /// The normal to the wall
 //Vector<double> Wall_normal;

 /// Trace file
 ofstream Trace_file;

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
 Problem::mesh_pt() = new StaticCapMesh<ELEMENT>(Nx,Nh1,Nh2);

 //Overwrite the linear solver (frontal solver) 
 //Problem::linear_solver_pt() = new HSL_MA42;

 //Hijack one of the pressure values in the upper fluid. Its value
 //will affect the residual of that element but it will not
 //be determined by it!
 Data* hijacked_pressure_data_pt = dynamic_cast<ELEMENT*>(
  mesh_pt()->upper_layer_element_pt(0))->hijack_internal_value(0,0);

 // Loop over the elements on the interface to pass pointer to Ca and
 // the pointer to the Data item that contains the single
 // (pressure) value that is "traded" for the volume constraint 
 unsigned n_interface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<n_interface;e++)
  {
   //Cast to a 1D element
   FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*el_pt=
   dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
    (mesh_pt()->interface_element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;

   //Pass the Data item that contains the single (pressure) value
   // that has been "traded" for the volume constraint to the
   // surface elements -- hacky! 
   el_pt->set_traded_pressure_data(hijacked_pressure_data_pt);
  }
 

 
 //Finally, pass the Data item that contains the single (pressure)
 // value that has been "traded" for the volume constraint
 // to the volume constraint element.
 {
  SpineVolumeConstraintPointElement<SpineElement<ELEMENT> >* el_pt =
   dynamic_cast<SpineVolumeConstraintPointElement<SpineElement<ELEMENT> >*>
   (mesh_pt()->point_element_pt());
    
  el_pt->volume_pt() = &Volume;
  el_pt->set_traded_pressure_data(hijacked_pressure_data_pt);
 }    
 

 //Set the boundary conditions

 //Pin the velocities on all boundaries apart from the symmetry
 //line (boundary 3) where only the horizontal velocity is pinned
 unsigned n_bound=mesh_pt()->nboundary();
 for (unsigned b=0;b<n_bound;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     mesh_pt()->boundary_node_pt(b,n)->pin(0);
     if (b!=3)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(1);
      }
    }
  }


 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
// dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
 // (mesh_pt()->interface_element_pt(n_interface-1))->
 // set_contact_angle_right(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->set_contact_angle(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->ca_pt() = &Ca;

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->wall_unit_normal_fct_pt() = 
  &Global_Physical_Variables::wall_unit_normal_fct;

 
 // Pin a single pressure value: Set the pressure dof 0 in element 0
 // to zero.
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);


 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
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
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Number of interface elements
// unsigned ninterface=mesh_pt()->ninterface_element();

 // Number of spines
 unsigned nspine=mesh_pt()->nspine();

 // Doc
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << mesh_pt()->spine_pt(0)->height();
 Trace_file << " "  << mesh_pt()->spine_pt(nspine-1)->height();
 Trace_file << " "  
 //            << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                mesh_pt()->interface_element_pt(0))->
//                actual_contact_angle_left()*
//                180.0/MathematicalConstants::Pi << " " ;
//  Trace_file << " "  
//             << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                mesh_pt()->interface_element_pt(ninterface-1))->
//                actual_contact_angle_right()*180.0/MathematicalConstants::Pi 
            << " ";
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             mesh_pt()->upper_layer_element_pt(0))->p_nst(0);
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             mesh_pt()->lower_layer_element_pt(0))->p_nst(0);
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









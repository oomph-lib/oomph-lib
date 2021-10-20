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
#include "constitutive.h"
#include "solid.h"

// The mesh
#include "meshes/two_layer_spine_mesh.h"
#include "meshes/rectangular_quadmesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables
{
 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 ///Direction of the wall normal vector
 Vector<double> Wall_normal;

 ///  Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }

}

//============================================================================
///A Problem class that solves the Navier--Stokes equations + free surface
///in a 2D geometry using a spine-based node update
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
 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->set_contact_angle(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->ca_pt() = &Ca;

 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(0))->wall_unit_normal_fct_pt() = 
  &Global_Physical_Variables::wall_unit_normal_fct;
 
 // Pin a single pressure value: Set the pressure dof 0 in element 0
 // of the lower layer to zero.
 dynamic_cast<ELEMENT*>(Bulk_mesh_pt->lower_layer_element_pt(0))
  ->fix_pressure(0,0.0);

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
 unsigned lower_boundaries[3]={0,1,5};
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
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->upper_layer_element_pt(0))->p_nst(0);
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->lower_layer_element_pt(0))->p_nst(0);
 Trace_file << " " << 2.0*cos(Angle)/Ca;
 Trace_file << std::endl;

}



//==start_of_specific_mesh_class==========================================
/// Two layer mesh which employs a pseudo-solid node-update strategy.
/// This class is essentially a wrapper to an 
/// ElasticRectangularQuadMesh, with an additional boundary
/// to represent the interface between the two fluid layers. In addition,
/// the mesh paritions the elements into those above and below
/// the interface and relabels boundaries so that we can impose
/// a volume constraint on the lower or upper fluid.
///
///                                 3
///               ---------------------------------------
///               |                                     |
///             4 |                                     | 2
///               |                 6                   |
///               ---------------------------------------
///               |                                     |
///             5 |                                     | 1
///               |                                     |
///               ---------------------------------------
///                                 0
//========================================================================
template <class ELEMENT>
class ElasticTwoLayerMesh :
 public ElasticRectangularQuadMesh<ELEMENT>
{

public:

 ///  Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction in bottom and top layer, respectively,
 /// axial length and height of top and bottom layers, a boolean
 /// flag to make the mesh periodic in the x-direction, and pointer 
 /// to timestepper (defaults to Steady timestepper)
 ElasticTwoLayerMesh(const unsigned &nx, 
                     const unsigned &ny1,
                     const unsigned &ny2, 
                     const double &lx,
                     const double &h1,
                     const double &h2,
                     const bool& periodic_in_x=false,
                     TimeStepper* time_stepper_pt=
                     &Mesh::Default_TimeStepper)
  : 
  RectangularQuadMesh<ELEMENT>(nx,ny1+ny2,lx,h1+h2,
                               periodic_in_x,time_stepper_pt),
  ElasticRectangularQuadMesh<ELEMENT>(nx,ny1+ny2,lx,h1+h2,
                                      periodic_in_x,time_stepper_pt)
  {
   //Set up the pointers to elements in the upper and lower fluid
   //Calculate numbers of nodes in upper and lower regions
   unsigned long n_lower = nx*ny1;
   unsigned long n_upper = nx*ny2;
   //Loop over lower elements and push back onto the array
   Lower_layer_element_pt.resize(n_lower);
   for(unsigned e=0;e<n_lower;e++)
    {
     Lower_layer_element_pt[e] = this->finite_element_pt(e);
    }
   //Loop over upper elements and push back onto the array
   Upper_layer_element_pt.resize(n_upper);
   for(unsigned e=0;e<n_upper;e++)
    {
     Upper_layer_element_pt[e] =  this->finite_element_pt(n_lower + e);
    }
   //end of upper and lower fluid element assignment

   //Set the elements adjacent to the interface on both sides
   Interface_lower_boundary_element_pt.resize(nx);
   Interface_upper_boundary_element_pt.resize(nx);
   {
    unsigned count_lower=nx*(ny1-1);
    unsigned count_upper= count_lower + nx; 
    for(unsigned e=0;e<nx;e++)
     {
      Interface_lower_boundary_element_pt[e] =
     this->finite_element_pt(count_lower); ++count_lower;
      Interface_upper_boundary_element_pt[e] = 
       this->finite_element_pt(count_upper); ++count_upper;
     }
   } //end of bulk elements next to interface setup


   // Reset the number of boundaries
   this->set_nboundary(7);

   //Relabel the boundary nodes
   //Storage for the boundary coordinates that will be transferred directly
   Vector<double> b_coord;
   {
    //Store the interface level
    const double y_interface = h1 - this->Ymin;
    
    //Nodes on boundary 3 are now on boundaries 4 and 5
    unsigned n_boundary_node = this->nboundary_node(3);
    //Loop over the nodes remove them from boundary 3
    //and add them to boundary 4 or 5 depending on their y coordinate
    for(unsigned n=0;n<n_boundary_node;n++)
     {
      //Cache pointer to the node
      Node* const nod_pt = this->boundary_node_pt(3,n);
      //Get the boundary coordinates if set
      if(this->boundary_coordinate_exists(3))
       {
        b_coord.resize(nod_pt->ncoordinates_on_boundary(3));
        nod_pt->get_coordinates_on_boundary(3,b_coord);
        //Indicate that the boundary coordinates are to be set on the 
        //new boundaries
        this->Boundary_coordinate_exists[4]=true;
        this->Boundary_coordinate_exists[5]=true;
       }
      
      //Now remove the node from the old boundary
      nod_pt->remove_from_boundary(3);
      
      //Find the height of the node
      double y = nod_pt->x(1);
      //If it's above or equal to the interface, it's on boundary 4
      if(y >= y_interface) 
       {
        this->add_boundary_node(4,nod_pt);
        //Add the boundary coordinate if it has been set up
        if(this->Boundary_coordinate_exists[4])
         {
          nod_pt->set_coordinates_on_boundary(4,b_coord);
         }
       }
      //otherwise it's on boundary 5
      if(y <=  y_interface) 
       {
        this->add_boundary_node(5,nod_pt);
        //Add the boundary coordinate if it has been set up
        if(this->Boundary_coordinate_exists[5])
         {
          nod_pt->set_coordinates_on_boundary(5,b_coord);
         }
       }
     }
    
    //Now clear the boundary node information from boundary 3
    this->Boundary_node_pt[3].clear();
    
    //Relabel the nodes on boundary 2 to be on boundary 3
    n_boundary_node = this->nboundary_node(2);
    //Loop over the nodes remove them from boundary 2
    //and add them to boundary 3 
    for(unsigned n=0;n<n_boundary_node;n++)
     {
      //Cache pointer to the node
      Node* const nod_pt = this->boundary_node_pt(2,n);
      //Get the boundary coordinates if set
      if(this->boundary_coordinate_exists(2))
       {
        b_coord.resize(nod_pt->ncoordinates_on_boundary(2));
        nod_pt->get_coordinates_on_boundary(2,b_coord);
        this->Boundary_coordinate_exists[3]=true;
       }
      
      //Now remove the node from the boundary 2
      nod_pt->remove_from_boundary(2);
      //and add to boundary 3
      this->add_boundary_node(3,nod_pt);
      if(this->Boundary_coordinate_exists[3])
       {
        nod_pt->set_coordinates_on_boundary(3,b_coord);
       }
     }
    
    //Clear the information from boundary 2
    this->Boundary_node_pt[2].clear();
    
    //Storage for nodes that should be removed from boundary 1
    std::list<Node*> nodes_to_be_removed;
    
    //Nodes on boundary 1 are now on boundaries 1 and 2
    n_boundary_node = this->nboundary_node(1);
    //Loop over the nodes remove some of them from boundary 1
    for(unsigned n=0;n<n_boundary_node;n++)
     {
      //Cache pointer to the node
      Node* const nod_pt = this->boundary_node_pt(1,n);
      
      //Find the height of the node
      double y = nod_pt->x(1);
      //If it's above or equal to the interface it's on boundary 2
      if(y >= y_interface) 
       {
        //Get the boundary coordinates if set
        if(this->boundary_coordinate_exists(1))
         {
          b_coord.resize(nod_pt->ncoordinates_on_boundary(1));
          nod_pt->get_coordinates_on_boundary(1,b_coord);
          this->Boundary_coordinate_exists[2]=true;
         }
        
        //Now remove the node from the boundary 1 if above interace
        if(y > y_interface)
         {
          nodes_to_be_removed.push_back(nod_pt);
         }
        //Always add to boundary 2
        this->add_boundary_node(2,nod_pt);
        //Add the boundary coordinate if it has been set up
        if(this->Boundary_coordinate_exists[2])
         {
          nod_pt->set_coordinates_on_boundary(2,b_coord);
         }
       }
     }
    
    //Loop over the nodes that are to be removed from boundary 1 and remove
    //them
    for(std::list<Node*>::iterator it=nodes_to_be_removed.begin();
        it!=nodes_to_be_removed.end();++it)
     {
      this->remove_boundary_node(1,*it);
     }
    nodes_to_be_removed.clear();
   } //end of boundary relabelling

   //Add the nodes to the interface

   //Read out number of linear points in the element
   unsigned n_p = dynamic_cast<ELEMENT*>
    (this->finite_element_pt(0))->nnode_1d();
   
   //Add the nodes on the interface to the boundary 6
   //Storage for boundary coordinates (x-coordinate)
   b_coord.resize(1);
   this->Boundary_coordinate_exists[6];
   //Starting index of the nodes
   unsigned n_start=0;
   for(unsigned e=0;e<nx;e++)
    {
     //If we are past the 
     if(e>0) {n_start=1;}
     //Get the layer of elements just above the interface
     FiniteElement* el_pt = this->finite_element_pt(nx*ny1+e);
     //The first n_p nodes lie on the boundary
     for(unsigned n=n_start;n<n_p;n++)
      {
       Node* nod_pt = el_pt->node_pt(n);
       this->convert_to_boundary_node(nod_pt);
       this->add_boundary_node(6,nod_pt);
       b_coord[0] = nod_pt->x(0);
       nod_pt->set_coordinates_on_boundary(6,b_coord);
      }
    }

   // Set up the boundary element information
   this->setup_boundary_element_info();
  }

 ///Access functions for pointers to elements in upper layer
 FiniteElement* &upper_layer_element_pt(const unsigned long &i) 
  {return Upper_layer_element_pt[i];}

 ///Access functions for pointers to elements in bottom layer
 FiniteElement* &lower_layer_element_pt(const unsigned long &i) 
  {return Lower_layer_element_pt[i];}

 ///Number of elements in upper layer
 unsigned long nupper() const {return Upper_layer_element_pt.size();}

 ///Number of elements in top layer
 unsigned long nlower() const {return Lower_layer_element_pt.size();}

 ///Access functions for pointers to elements in upper layer
 FiniteElement* &interface_upper_boundary_element_pt(const unsigned long &i) 
  {return Interface_upper_boundary_element_pt[i];}

 ///Access functions for pointers to elements in bottom layer
 FiniteElement* &interface_lower_boundary_element_pt(const unsigned long &i) 
  {return Interface_lower_boundary_element_pt[i];}

 ///Number of elements in upper layer
 unsigned long ninterface_upper() const 
 {return Interface_upper_boundary_element_pt.size();}

 ///Number of elements in top layer
 unsigned long ninterface_lower() const 
 {return Interface_lower_boundary_element_pt.size();}

 /// Index of the face of the elements next to the interface
 ///in the upper region (always -2)
 int interface_upper_face_index_at_boundary(const unsigned &e)
 {return -2;}

 /// Index of the face of the elements next to the interface in
 /// the lower region (always 2)
 int interface_lower_face_index_at_boundary(const unsigned &e)
 {return 2;}

private:

 /// Vector of pointers to element in the upper layer
 Vector <FiniteElement *> Lower_layer_element_pt;

 /// Vector of pointers to element in the lower layer
 Vector <FiniteElement *> Upper_layer_element_pt;

 ///  Vector of pointers to the elements adjacent to the interface
 /// on the lower layer
 Vector <FiniteElement*> Interface_lower_boundary_element_pt;

 ///  Vector of pointers to the element adjacent to the interface
 /// on the upper layer
 Vector<FiniteElement *> Interface_upper_boundary_element_pt;


}; // End of specific mesh class



//===========start_of_pseudo_elastic_class====================================
/// A class that solves the Navier--Stokes equations
///to compute the shape of a static interface between two fluids in a 
///rectangular container with an imposed contact angle at the boundary.
//============================================================================
template<class ELEMENT>
class PseudoSolidCapProblem : public Problem
{
public:

 //Constructor: Arguements are the number of elements in the x
 //direction and in each of the layers
 PseudoSolidCapProblem(
  const unsigned &Nx, const unsigned &Nh1, const unsigned &Nh2);

 /// Destructor: clean up memory allocated by the object
 ~PseudoSolidCapProblem();

 /// Peform a parameter study: Solve problem for a range of contact angles
 /// Pass name of output directory as a string
 void parameter_study(const string& dir_name);

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


private:
 
 /// Create the free surface elements
 void create_free_surface_elements();

 /// Create the volume constraint elements
 void create_volume_constraint_elements();

 /// Create the contact angle element
 void create_contact_angle_element();

 /// The Capillary number 
 double Ca;

 /// The prescribed volume of the fluid
 double Volume;

 /// The external pressure
 double Pext;

 /// The contact angle
 double Angle;

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 // Pointer to the (single valued) Data item that
 // will contain the pressure value that we're
 // trading for the volume constraint
 Data* Traded_pressure_data_pt;

 /// Trace file
 ofstream Trace_file;

 ///Storage for the bulk mesh
 ElasticTwoLayerMesh<ELEMENT>* Bulk_mesh_pt;

 /// Storage for the free surface mesh
 Mesh* Free_surface_mesh_pt;

 /// Storage for the element bounding the free surface
 Mesh* Free_surface_bounding_mesh_pt;

 /// Storage for the elements that compute the enclosed volume
 Mesh* Volume_computation_mesh_pt;

 /// Storage for the volume constraint
 Mesh* Volume_constraint_mesh_pt;

}; //end_of_pseudo_solid_problem_class



//============start_of_constructor=====================================
/// Constructor: Pass boolean flag to indicate if the volume
/// constraint is applied by hijacking an internal pressure
/// or the external pressure
//======================================================================
template<class ELEMENT>
PseudoSolidCapProblem<ELEMENT>::PseudoSolidCapProblem(
 const unsigned &Nx, const unsigned &Nh1, const unsigned &Nh2) :
 Ca(2.1),       //Initialise value of Ca to some random value
 Volume(0.5),   //Initialise the value of the volume 
 Pext(1.23),    //Initialise the external pressure to some random value
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 //Set the wall normal
 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = 1.0; 
 Global_Physical_Variables::Wall_normal[1] = 0.0;

 //Construct mesh
 Bulk_mesh_pt = new ElasticTwoLayerMesh<ELEMENT>(Nx,Nh1,Nh2,0.5,1.0,1.0);
 
 //Hijack one of the pressure values in the fluid and use it 
 //as the pressure whose value is determined by the volume constraint.
 //(Its value will affect the residual of that element but it will not
 //be determined by it, i.e. it's hijacked).
 Traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
  Bulk_mesh_pt->upper_layer_element_pt(0))->hijack_internal_value(0,0);

 //Set the constituive law
 Constitutive_law_pt =  
  new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 //Loop over the elements to set the consitutive law
 unsigned n_bulk = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_bulk;e++)
  {
   ELEMENT* el_pt = 
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   el_pt->constitutive_law_pt() = Constitutive_law_pt;
  }
 
 //Set the boundary conditions

 //Fluid velocity conditions
 //Pin the velocities on all external boundaries apart from the symmetry
 //line boundaries 4 and 5) where only the horizontal velocity is pinned
 for (unsigned b=0;b<6;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
     if((b!=4) && (b!=5))
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }
    }
  } //end_of_fluid_boundary_conditions

 //PesudoSolid boundary conditions,
 for(unsigned b=0;b<6;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     //Pin vertical displacement on the bottom and top
     if((b==0) || (b==3))
      {
       static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n))
        ->pin_position(1);
      }
     //Pin horizontal displacement on the sides
     if((b==1) || (b==2) || (b==4) || (b==5))
      {
       static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n))
        ->pin_position(0);
      }
    }
  } //end_of_solid_boundary_conditions
 
 //Constrain all nodes only to move vertically (not horizontally)
 {
  unsigned n_node = Bulk_mesh_pt->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(0);
   }
 } //end_of_constraint

 // Pin a single pressure value in the lower fluid: 
 // Set the pressure dof 0 in element 0 to zero.
 dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->fix_pressure(0,0.0);

 //Create the free surface elements
 create_free_surface_elements();

 //Create the volume constraint elements
 create_volume_constraint_elements();

 //Need to make the bounding element
 create_contact_angle_element();

 //Now need to add all the meshes
 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Free_surface_mesh_pt);
 this->add_sub_mesh(Volume_computation_mesh_pt);
 this->add_sub_mesh(Volume_constraint_mesh_pt);
 this->add_sub_mesh(Free_surface_bounding_mesh_pt);

 //and build the global mesh
 this->build_global_mesh();
 
 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
} //end_of_constructor


//==========================================================================
/// Destructor. Make sure to clean up all allocated memory, so that multiple
/// instances of the problem don't lead to excessive memory usage.
//==========================================================================
template<class ELEMENT>
PseudoSolidCapProblem<ELEMENT>::~PseudoSolidCapProblem() 
{
 //Delete the contact angle element
 delete Free_surface_bounding_mesh_pt->element_pt(0);
 Free_surface_bounding_mesh_pt->flush_element_and_node_storage();
 delete Free_surface_bounding_mesh_pt;
 //Delete the volume constraint mesh
 delete Volume_constraint_mesh_pt;
 //Delete the surface volume computation elements
 unsigned n_element = Volume_computation_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {delete Volume_computation_mesh_pt->element_pt(e);}
 //Now flush the storage
 Volume_computation_mesh_pt->flush_element_and_node_storage();
 //Now delete the mesh
 delete Volume_computation_mesh_pt;
 //Delete the free surface elements
 n_element = Free_surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {delete Free_surface_mesh_pt->element_pt(e);}
 //Now flush the storage
 Free_surface_mesh_pt->flush_element_and_node_storage();
 //Now delete the mesh
 delete Free_surface_mesh_pt;

 //Delete the constitutive law
 delete Constitutive_law_pt;

 //Then delete the bulk mesh
 delete Bulk_mesh_pt;
}

//============create_free_surface_element================================
/// Create the free surface elements
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_free_surface_elements()
{
 //Allocate storage for the free surface mesh
 Free_surface_mesh_pt = new Mesh;

 // How many bulk fluid elements are adjacent to the interface
 unsigned n_element = Bulk_mesh_pt->ninterface_lower();
 //Loop over them
 for(unsigned e=0;e<n_element;e++)
  {
   // Create new element
   ElasticLineFluidInterfaceElement<ELEMENT>* el_pt =
    new ElasticLineFluidInterfaceElement<ELEMENT>(
     Bulk_mesh_pt->interface_lower_boundary_element_pt(e),
     Bulk_mesh_pt->interface_lower_face_index_at_boundary(e));
   
   // Add it to the mesh
   Free_surface_mesh_pt->add_element_pt(el_pt);
   
   //Add the capillary number
   el_pt->ca_pt() = &Ca;
  } 

}


//============start_of_create_volume_constraint_elements==================
/// Create the volume constraint elements
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_volume_constraint_elements()
{
 //Build the single volume constraint element
 Volume_constraint_mesh_pt = new Mesh;
 VolumeConstraintElement* vol_constraint_element = 
  new VolumeConstraintElement(&Volume,Traded_pressure_data_pt,0);
 Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);
  
  //Now create the volume computation elements
  Volume_computation_mesh_pt = new Mesh;

  //Loop over lower boundaries (or just 1, why?)
  unsigned lower_boundaries[3]={0,1,5};
  //Loop over all boundaries
  for(unsigned i=0;i<3;i++)
   {
    //Read out the actual boundary
    unsigned b= lower_boundaries[i];
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
      ElasticLineVolumeConstraintBoundingElement<ELEMENT>* el_pt =
       new ElasticLineVolumeConstraintBoundingElement<ELEMENT>(
        bulk_elem_pt,face_index);   

      //Set the "master" volume control element
      el_pt->set_volume_constraint_element(vol_constraint_element);
     
      // Add it to the mesh
      Volume_computation_mesh_pt->add_element_pt(el_pt);     
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
     ElasticLineVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new ElasticLineVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }


} //end_of_create_volume_constraint_elements

//==========start_of_create_contact_angle_elements========================
/// Create the contact angle element
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_contact_angle_element()
{
 Free_surface_bounding_mesh_pt = new Mesh;

 //Find the element at the end of the free surface
 //The elements are assigned in order of increasing x coordinate
 unsigned n_free_surface = Free_surface_mesh_pt->nelement();
  
 //Make the bounding element for the contact angle constraint
 //which works because the order of elements in the mesh is known
 FluidInterfaceBoundingElement* el_pt = 
  dynamic_cast<ElasticLineFluidInterfaceElement<ELEMENT>*> 
  (Free_surface_mesh_pt->element_pt(n_free_surface-1))->
  make_bounding_element(1);
 
 //Set the contact angle (strong imposition)
 el_pt->set_contact_angle(&Angle);
 
 //Set the capillary number
 el_pt->ca_pt() = &Ca;
 
 //Set the wall normal of the external boundary
 el_pt->wall_unit_normal_fct_pt() 
  =  &Global_Physical_Variables::wall_unit_normal_fct;

 //Add the element to the mesh
 Free_surface_bounding_mesh_pt->add_element_pt(el_pt);

} //end_of_create_contact_angle_element



//================start_of_parameter_study===========================
/// Perform a parameter study. Pass name of output directory as 
/// a string
//======================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::parameter_study(const string& dir_name)
{
 // Create DocInfo object (allows checking if output directory exists)
 DocInfo doc_info;
 doc_info.set_directory(dir_name);
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 Trace_file << "VARIABLES=\"<greek>a</greek><sub>prescribed</sub>\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\",";
 Trace_file << "\"p<sub>fluid</sub>-p<sub>ext</sub>\",";
 Trace_file << "\"<greek>D</greek>p<sub>exact</sub>\"";
 Trace_file << std::endl;

 // Doc initial state
 doc_solution(doc_info);

 // Bump up counter
 doc_info.number()++;

 //Solve the problem for six different contact angles
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




//==============start_of_doc_solution=====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 //Output stream
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 //Output domain
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();


 // Number of interface elements
 unsigned ninterface=Free_surface_mesh_pt->nelement();
 //Find number of nodes in the last interface element
 unsigned np = Free_surface_mesh_pt->finite_element_pt(ninterface-1)->nnode();
 // Document the contact angle (in degrees), the height of the interface at 
 // the centre of the container, the height at the wall, the computed 
 // pressure drop across the interface and 
 // the analytic prediction of the pressure drop.
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << Free_surface_mesh_pt->finite_element_pt(0)
  ->node_pt(0)->x(1)
            << " " 
            << Free_surface_mesh_pt->finite_element_pt(ninterface-1)
  ->node_pt(np-1)->x(1);
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->upper_layer_element_pt(0))->p_nst(0);
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->lower_layer_element_pt(0))->p_nst(0);
 Trace_file << " " << 2.0*cos(Angle)/Ca;
 Trace_file << std::endl;

} //end_of_doc_solution


//======================================================================
///Main driver: Build problem and initiate parameter study
//======================================================================
int main()
{
 {
 //Construct the problem with 4 x (4+4) elements
 CapProblem<Hijacked<QCrouzeixRaviartElement<2> >  > problem(4,4,4);

 //Solve the problem :)
 problem.parameter_study();
 }

 {
 //Construct the problem with 4 x (4+4) elements
 PseudoSolidCapProblem<Hijacked<
     PseudoSolidNodeUpdateElement<QCrouzeixRaviartElement<2>,
                                  QPVDElementWithPressure<2> > > > 
  problem(4,4,4);
 
 //Solve the problem :)
 problem.parameter_study("RESLT_elastic");
 }

}

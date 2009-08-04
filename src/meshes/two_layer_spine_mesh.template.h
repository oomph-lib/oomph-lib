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
#ifndef OOMPH_TWO_LAYER_SPINE_MESH_HEADER
#define OOMPH_TWO_LAYER_SPINE_MESH_HEADER


// The mesh
#include "../generic/spines.h"
#include "rectangular_quadmesh.template.h"

namespace oomph
{


//======================================================================
/// Two-layer spine mesh class derived from standard 2D mesh.
/// The mesh contains two layers of spinified fluid elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
/// and an intermediate interface layer of corresponding Spine 
/// interface elements, of type INTERFACE_ELEMENT, e.g. 
/// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar problems.
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class TwoLayerSpineMesh : public RectangularQuadMesh<ELEMENT >, 
                          public SpineMesh
{

public:

 /// \short Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction in bottom and top layer, respectively,
 /// axial length and height of top and bottom layers and pointer 
 /// to timestepper (defaults to Steady timestepper)
 TwoLayerSpineMesh(const unsigned &nx, 
                   const unsigned &ny1,
                   const unsigned &ny2, 
                   const double &lx,
                   const double &h1,
                   const double &h2,
                   TimeStepper* time_stepper_pt=
                   &Mesh::Default_TimeStepper);


 /// \short Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction in bottom and top layer, respectively,
 /// axial length and height of top and bottom layers, a boolean
 /// flag to make the mesh periodic in the x-direction, and pointer 
 /// to timestepper (defaults to Steady timestepper)
 TwoLayerSpineMesh(const unsigned &nx, 
                   const unsigned &ny1,
                   const unsigned &ny2, 
                   const double &lx,
                   const double &h1,
                   const double &h2,
                   const bool& periodic_in_x,
                   TimeStepper* time_stepper_pt=
                   &Mesh::Default_TimeStepper);
 
 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Number of elements on interface
 unsigned long ninterface_element() const {return Interface_element_pt.size();}
 
 /// \short Reorder the elements so we loop over them vertically first
 /// (advantageous in "wide" domains if a frontal solver is used).
 void element_reorder();

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

 /// \short General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific update
 /// actions, depending on the node update fct id stored for each node.
 void spine_node_update(SpineNode* spine_node_pt)
  {
   unsigned id=spine_node_pt->node_update_fct_id();
   switch(id)
    {
    case 0:
     spine_node_update_lower(spine_node_pt);
     break;
     
    case 1:
     spine_node_update_upper(spine_node_pt);
     break;

    default:
     std::ostringstream error_message;
     error_message << "Unknown id passed to spine_node_update " << id 
                   << std::endl;
     throw OomphLibError(error_message.str(),
                         "TwoLayerSpineMesh::spine_node_update()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }



protected:

 /// Number of elements in lower layer
 unsigned Ny1;

 /// Number of elements in upper layer
 unsigned Ny2;

 /// \short  Height of the lower layer
 double H1;

 /// \short  Height of the upper layer
 double H2;

 /// Vector of pointers to element in the upper layer
 Vector <FiniteElement *> Lower_layer_element_pt;

 /// Vector of pointers to element in the lower layer
 Vector <FiniteElement *> Upper_layer_element_pt;
 
 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

 /// \short The spacing function for the x co-ordinates with two 
 /// regions.
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// \short The spacing function for the y co-ordinates with three
 /// regions in each fluid.
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// \short Update function for the lower part of the domain
 void spine_node_update_lower(SpineNode *spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of y
   spine_node_pt->x(1) = this->Ymin + W*H;
  }
 

 /// \short Update function for the upper part of the domain
 void spine_node_update_upper(SpineNode *spine_node_pt)
  {
   //Get fraction alon the spine
   double W = spine_node_pt->fraction();

   //Get spine height
   double H = spine_node_pt->h();  
 
   //Set the value of y
   spine_node_pt->x(1) = (this->Ymin+H) + W *(this->Ymax - (this->Ymin+H) );
  }

 /// \short Helper function to actually build the two-layer spine mesh 
 /// (called from various constructors)
 virtual void build_two_layer_mesh(TimeStepper* time_stepper_pt);

};



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//======================================================================
/// \short Axisymmetric-Two-layer spine mesh class derived from standard 
/// Two-Layer spine mesh.
/// This creates a finer mesh around boundaries 0, 1 and 2 (bottom, 
/// right and top) and also around the interface. This mesh was designed
/// for use with the spin-up and spin-over problems in an axisymmetric 
/// domain. Using the same spacing fractions in the both the upper and 
/// lower fluids.
//
//                     Xfraction
//      _______________________ 
//     |                 |     |
//     |                 |     |nya2
//     |-----------------+-----|      1.0-Yfraction1
//     |                 |     |
//     |                 |     |
//     |                 |     |nyb2
//     |                 |     |
//     |                 |     |
//     |-----------------+-----|      1.0-Yfraction2
//     |                 |     |nyc2
//     |_________________|_____|   interface
//     |                 |     |
//     |                 |     |nyc1
//     |-----------------+-----|      Yfraction2
//     |                 |     |
//     |                 |     |
//     |                 |     |
//     |                 |     |nyb1
//     |                 |     |
//     |                 |     |
//     |                 |     |
//     |-----------------+-----|      Yfraction1
//     |                 |     |nya1
//     |_________________|_____|
//        nxa              nxb
//
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
class Axisym2x6TwoLayerSpineMesh : 
public virtual TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>
{
private:
 //Data to hold the number of elements in each section of the mesh
 unsigned Nxa, Nxb, Nya1, Nyb1, Nyc1, Nya2, Nyb2, Nyc2;
 //Data that holds the fractions along the channel in which to introduce
 //the second zones
 double Xfraction, Yfraction1, Yfraction2;

public:

 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in right 
 /// region; nya1: number of elements in y direction in fluid 1 on the 
 /// bottom boundary 0; nyb1: number of elements in y direction in the 
 /// core region of fluid 1; nyc1: number of elements in y direction in 
 /// the region next to the interface for fluid 1; nya2: number of
 /// elements in y direction in fluid 2 on the top boundary 2; nyb2: 
 /// number of elements in y direction in the core region of fluid 2; 
 /// nyc2: number of elements in y direction in the region next to the 
 /// interface for fluid 2; lx: length of mesh in the x direction; h1: 
 /// height of fluid 1; h2: height of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym2x6TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const double &lx,
                            const double &h1, 
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);


 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in right 
 /// region; Xfrac: the spacing fraction for nxl; nya1: number of 
 /// elements in y direction in fluid 1 on the bottom boundary 0; nyb1:
 /// number of elements in y direction in the core region of fluid 1; 
 /// nyc1: number of elements in y direction in the region next to the
 /// interface for fluid 1; nya2: number of elements in y direction in 
 /// fluid 2 on the top boundary 2; nyb2: number of elements in y 
 /// direction in the core region of fluid 2; nyc2: number of elements 
 /// in y direction in the region next to the interface for fluid 2; 
 /// Yfrac1: the spacing fraction for nya; Yfrac2: the spacing fraction
 /// for nyb; lx: length of mesh in the x direction; h1: height of 
 /// fluid 1; h2: height of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym2x6TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const double &Xfrac,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const double &Yfrac1,
                            const double &Yfrac2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);


 /// \short The spacing function for the x co-ordinates with two 
 /// regions.
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// \short The spacing function for the y co-ordinates with three
 /// regions in each fluid.
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);
};


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//======================================================================
/// \short Axisymmetric-Two-layer spine mesh class derived from standard 
/// Two-Layer spine mesh.
/// This creates a finer mesh around boundaries 0, 1, 2 and 3 and also 
/// around the interface. This mesh was designed
/// for use with the spin-up and spin-over problems in an axisymmetric 
/// domain. Using the same spacing fractions in the both the upper and 
/// lower fluids.
//
//      Xfraction1     Xfraction2
//      _______________________ 
//     |    |            |     |
//     |    |            |     |nya2
//     |----+------------+-----|      1.0-Yfraction1
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |nyb2
//     |    |            |     |
//     |    |            |     |
//     |----+------------+-----|      1.0-Yfraction2
//     |    |            |     |nyc2
//     |____|____________|_____|   interface
//     |    |            |     |
//     |    |            |     |nyc1
//     |----+------------+-----|      Yfraction2
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |nyb1
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |
//     |----+------------+-----|      Yfraction1
//     |    |            |     |nya1
//     |____|____________|_____|
//      nxa     nxb        nxc
//
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
class Axisym3x6TwoLayerSpineMesh : 
public virtual TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>
{
private:
 //Data to hold the number of elements in each section of the mesh
 unsigned Nxa, Nxb, Nxc, Nya1, Nyb1, Nyc1, Nya2, Nyb2, Nyc2;
 //Data that holds the fractions along the channel in which to introduce
 //the second zones
 double Xfraction1, Xfraction2, Yfraction1, Yfraction2;

public:

 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in centre 
 /// region; nxc: number of elements in the x direction in the right 
 /// region; nya1: number of elements in y direction in fluid 1 on 
 /// the bottom boundary 0; nyb1: number of elements in y direction in
 /// the core region of fluid 1; nyc1: number of elements in y direction
 /// in the region next to the interface for fluid 1; nya2: number of 
 /// elements in the y direction in fluid 2 on the top boundary 2; nyb2:
 /// number of elements in the y direction in the core region of fluid 2;
 /// nyc2: number of elements in the y direction in the region next to 
 /// the interface for fluid 2; lx: length of mesh in the x direction; 
 /// h1: the height of fluid 1; h2: the height of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym3x6TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const unsigned &nxc,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);
 

 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in centre
 /// region; nxc: number of elements in the x direction in the right region
 /// Xfrac1: the spacing fraction for nxa; Xfrac2: the spacing fraction for
 /// nxb; nya1: number of elements in y direction in fluid 1 on the bottom
 /// boundary 0; nyb1: number of elements in y direction in the core region 
 /// of fluid 1; nyc1: number of elements in y direction in the region next
 /// to the interface for fluid 1; nya2: number of elements in the y
 /// direction in fluid 2 on the top boundary 2; nyb2: number of elements
 /// in the y direction in the core region of fluid 2; nyc2: number of
 /// elements in the y direction in the region next to the interface for 
 /// fluid 2; Yfrac1: the spacing fraction for nya; Yfrac2: the spacing 
 /// fraction for nyb; lx: length of mesh in the x direction; h1: the 
 /// height of fluid 1; h2: the height of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym3x6TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const unsigned &nxc,
                            const double &Xfrac1,
                            const double &Xfrac2,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const double &Yfrac1,
                            const double &Yfrac2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);


 /// \short The spacing function for the x co-ordinates with two 
 /// regions.
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// \short The spacing function for the y co-ordinates with three
 /// regions in each fluid.
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);
};



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


 //======================================================================
/// Axisymmetric-Two-layer spine mesh class derived from standard 
/// Two-Layer spine mesh.
/// This creates a finer mesh around boundaries 0, 1 and 2 (bottom, 
/// right and top) and also around the interface. This mesh was designed
/// for use with the spin-up and spin-over problems in an axisymmetric 
/// domain. Using the same spacing fractions in the both the upper and 
/// lower fluids.
//
//      Xfraction1     Xfraction2
//      _______________________ 
//     |    |            |     |
//     |    |            |     |nya2
//     |----+------------+-----|      1.0-Yfraction1
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |nyb2
//     |    |            |     |
//     |    |            |     |
//     |----+------------+-----|      1.0-Yfraction2
//     |    |            |     |nyc2
//     |----+------------+-----|      1.0-Yfraction3
//     |____|____________|_____|nyd2   interface
//     |    |            |     |nyd1
//     |----+------------+-----|      Yfraction3
//     |    |            |     |nyc1
//     |----+------------+-----|      Yfraction2
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |nyb1
//     |    |            |     |
//     |    |            |     |
//     |    |            |     |
//     |----+------------+-----|      Yfraction1
//     |    |            |     |nya1
//     |____|____________|_____|
//      nxa     nxb        nxc
//
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT >
class Axisym3x8TwoLayerSpineMesh : 
public virtual TwoLayerSpineMesh<ELEMENT,INTERFACE_ELEMENT>
{
private:
 //Data to hold the number of elements in each section of the mesh
 unsigned Nxa, Nxb, Nxc, Nya1, Nyb1, Nyc1, Nyd1, Nya2, Nyb2, Nyc2, Nyd2;
 //Data that holds the fractions along the channel in which to introduce
 //the second zones
 double Xfraction1, Xfraction2, Yfraction1, Yfraction2, Yfraction3;

public:

 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in centre 
 /// region; nxc: number of elements in the x direction in right region; 
 /// nya1: number of elements in y direction in fluid 1 on 
 /// the bottom boundary 0; nyb1: number of elements in y direction in
 /// the core region of fluid 1; nyc1: number of elements in y direction
 /// in the region near to the interface for fluid 1; nyd1: number of 
 /// elements in the y direction in the region next to the interface; 
 /// nya2: number of elements in the y direction in fluid 2 on the top 
 /// boundary 2; nyb2: number of elements in the y direction in the core
 /// region of fluid 2; nyc2: number of elements in the y direction in 
 /// the region near to the interface for fluid 2; nyd2: number of 
 /// elements in the y direction next to the interface; lx: length of 
 /// mesh in the x direction; h1: the height of fluid 1; h2: the height 
 /// of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym3x8TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const unsigned &nxc,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nyd1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const unsigned &nyd2,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);
 
 /// \short Simple constructor: nxa: number of elements in x direction
 /// in left region; nxb: number of elements in x direction in centre 
 /// region; nxc: number of elements in the x direction in right region; 
 /// Xfrac1: the spacing fraction for nxa; Xfrac2: the spacing fraction
 /// for nxb; nya1: number of elements in y direction in fluid 1 on 
 /// the bottom boundary 0; nyb1: number of elements in y direction in
 /// the core region of fluid 1; nyc1: number of elements in y direction
 /// in the region near to the interface for fluid 1; nyd1: number of 
 /// elements in the y direction in the region next to the interface; 
 /// nya2: number of elements in the y direction in fluid 2 on the top 
 /// boundary 2; nyb2: number of elements in the y direction in the core
 /// region of fluid 2; nyc2: number of elements in the y direction in 
 /// the region near to the interface for fluid 2; nyd2: number of 
 /// elements in the y direction next to the interface; Yfrac1: the 
 /// spacing fraction for nya1 and nya2; Yfrac2: the spacing fraction 
 /// for nyb1 and nyb2; Yfrac3: the spacing fraction for nyc1 and nyc2;
 /// lx: length of mesh in the x direction; h1: the height of fluid 1; 
 /// h2: the height of fluid 2.
 /// Also pass pointer to timestepper (defaults to Static)
 Axisym3x8TwoLayerSpineMesh(const unsigned &nxa,
                            const unsigned &nxb,
                            const unsigned &nxc,
                            const double &Xfrac1,
                            const double &Xfrac2,
                            const unsigned &nya1,
                            const unsigned &nyb1,
                            const unsigned &nyc1,
                            const unsigned &nyd1,
                            const unsigned &nya2,
                            const unsigned &nyb2,
                            const unsigned &nyc2,
                            const unsigned &nyd2,
                            const double &Yfrac1,
                            const double &Yfrac2,
                            const double &Yfrac3,
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);


 /// \short The spacing function for the x co-ordinates with two 
 /// regions.
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// \short The spacing function for the y co-ordinates with three
 /// regions in each fluid.
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);


};


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

} 

#endif


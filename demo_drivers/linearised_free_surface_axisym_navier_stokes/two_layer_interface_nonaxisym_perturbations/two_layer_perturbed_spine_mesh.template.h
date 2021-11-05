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
#ifndef OOMPH_TWO_LAYER_PERTURBED_SPINE_MESH_HEADER
#define OOMPH_TWO_LAYER_PERTURBED_SPINE_MESH_HEADER

// The mesh
#include "perturbed_spines.h"
#include "meshes.h"
#include "generic.h"

namespace oomph
{


//======================================================================
/// Two-layer spine mesh class derived from standard 2D mesh.
/// The mesh contains two layers of spinified fluid elements (of type ELEMENT;
/// e.g  PerturbedSpineElement<QCrouzeixRaviartElement<2>)
/// and an intermediate interface layer of corresponding PerturbedSpine 
/// interface elements, of type INTERFACE_ELEMENT, e.g. 
/// PerturbedSpineLineFluidInterfaceElement<ELEMENT> for 2D planar problems.
//======================================================================
template <class ELEMENT>
class TwoLayerPerturbedSpineMesh
 : public RectangularQuadMesh<ELEMENT>, public PerturbedSpineMesh
{

public:

 /// Constructor: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 TwoLayerPerturbedSpineMesh(const unsigned &nx, 
                            const unsigned &ny1,
                            const unsigned &ny2, 
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);
 
 /// Constructor: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - boolean flag to make the mesh periodic in the x-direction
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 TwoLayerPerturbedSpineMesh(const unsigned &nx, 
                            const unsigned &ny1,
                            const unsigned &ny2, 
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            const bool& periodic_in_x,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);
 
 /// Constructor: Pass in:
 ///   - number of elements in x-direction
 ///   - number of elements in y-direction in bottom layer
 ///   - number of elements in y-direction in top layer
 ///   - axial length
 ///   - height of bottom layer
 ///   - height of top layer
 ///   - boolean flag to make the mesh periodic in the x-direction
 ///   - boolean flag to specify whether or not to call the
 ///     "build_two_layer_mesh" function
 ///   - pointer to the base mesh (where we compute the base flow)
 ///   - pointer to timestepper (defaults to Steady timestepper)
 TwoLayerPerturbedSpineMesh(const unsigned &nx, 
                            const unsigned &ny1,
                            const unsigned &ny2, 
                            const double &lx,
                            const double &h1,
                            const double &h2,
                            const bool& periodic_in_x,
                            const bool& build_mesh,
                            SpineMesh* &base_mesh_pt,
                            TimeStepper* time_stepper_pt=
                            &Mesh::Default_TimeStepper);
 
 
 /// Reorder the elements so we loop over them vertically first
 /// (advantageous in "wide" domains if a frontal solver is used).
 void element_reorder();

 /// Access functions for pointers to elements in upper layer
 FiniteElement* &upper_layer_element_pt(const unsigned long &i) 
  { return Upper_layer_element_pt[i]; }

 /// Access functions for pointers to elements in lower layer
 FiniteElement* &lower_layer_element_pt(const unsigned long &i) 
  { return Lower_layer_element_pt[i]; }

 /// Number of elements in upper layer
 unsigned long nupper() const { return Upper_layer_element_pt.size(); }

 /// Number of elements in top layer
 unsigned long nlower() const { return Lower_layer_element_pt.size(); }
 
 /// General node update function implements pure virtual function 
 /// defined in PerturbedSpineMesh base class and performs specific update
 /// actions, depending on the node update fct id stored for each node.
 /// This function sets the nodal positions to be the same as in the base
 /// state problem through the perturbed spine's pointer to its
 /// corresponding base spine. Additionally this function updates the
 /// perturbations to the nodal positions, which are stored as (pinned)
 /// values at the nodes.
 void perturbed_spine_node_update(PerturbedSpineNode* spine_node_pt)
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
    error_message << "Unknown id passed to perturbed_spine_node_update " << id 
                  << std::endl;
    throw OomphLibError(error_message.str(),
                        "TwoLayerPerturbedSpineMesh::perturbed_spine_node_update()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }
 
 /// Access function for number of elements in lower layer
 const unsigned& ny1() const { return Ny1; }
 
 /// Access function for number of elements in upper layer
 const unsigned& ny2() const { return Ny2; }
 
 /// Set up the internal data YC_index and YS_index, which
 /// store the indices at which the perturbations to the nodal
 /// y-position are stored in the bulk element
 const void set_perturbation_to_nodal_positions_indices(
  const unsigned& cosine_index,
  const unsigned& sine_index)
 {
  YC_index = cosine_index;
  YS_index = sine_index;
 }

  protected:
 
 /// Number of elements in lower layer
 unsigned Ny1;

 /// Number of elements in upper layer
 unsigned Ny2;

 ///  Height of the lower layer
 double H1;

 ///  Height of the upper layer
 double H2;

 /// Pointer to corresponding mesh of base state problem
 //  Note: this has got to be a mesh that inherits from the SpineMesh
 //  base class. We store this so that the build_two_layer_mesh function
 //  can pass pointers to the corresponding base spines to each of the
 //  perturbed spines as it builds them
 SpineMesh* Base_mesh_pt;

 /// Index at which the cosine part of the perturbation to the
 /// nodal y-position is stored in the bulk element. The mesh needs to
 /// know this so that the value of the nodal position perturbation can
 /// be updated during the "perturbed spine node update" procedure
 int YC_index;

 /// Index at which the sine part of the perturbation to the
 /// nodal y-position is stored in the bulk element. The mesh needs to
 /// know this so that the value of the nodal position perturbation can
 /// be updated during the "perturbed spine node update" procedure
 int YS_index;

 /// Vector of pointers to element in the lower layer
 Vector <FiniteElement *> Lower_layer_element_pt;

 /// Vector of pointers to element in the upper layer
 Vector <FiniteElement *> Upper_layer_element_pt;
 
 /// The spacing function for the x co-ordinates with two 
 /// regions.
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// The spacing function for the y co-ordinates with three
 /// regions in each fluid.
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode);

 /// Update function for the lower part of the domain
 void spine_node_update_lower(PerturbedSpineNode *spine_node_pt)
  {
   // Get fraction along the spine
   const double w = spine_node_pt->fraction();

   // Get base spine height
   const double h =
    spine_node_pt->perturbed_spine_pt()->base_spine_pt()->height();

   // Set the nodal y-position
   spine_node_pt->x(1) = this->Ymin + w*h;

   // Get perturbed spine height (cosine and sine parts)
   const double HC = spine_node_pt->perturbed_spine_pt()->height(0);
   const double HS = spine_node_pt->perturbed_spine_pt()->height(1);

   // Set the perturbation (cosine and sine parts) to the nodal
   // y-position. Note that YC_index and YS_index are the indices at
   // which the perturbations to the nodal y-position are stored in
   // the bulk elements
   if(YC_index>=0)
    {
     spine_node_pt->set_value(YC_index,w*HC);
    }
   else
    {
     std::ostringstream error_message;
     error_message << "The cosine part of the perturbation to the nodal\n"
                   << "y-position has not been set up." << std::endl;
     throw OomphLibError(
      error_message.str(),
      "TwoLayerPerturbedSpineMesh::spine_node_update_lower(...)",
      OOMPH_EXCEPTION_LOCATION);
    }
   if(YS_index>=0)
    {
     spine_node_pt->set_value(YS_index,w*HS);
    }
   else
    {
     std::ostringstream error_message;
     error_message << "The sine part of the perturbation to the nodal\n"
                   << "y-position has not been set up." << std::endl;
     throw OomphLibError(
      error_message.str(),
      "TwoLayerPerturbedSpineMesh::spine_node_update_lower(...)",
      OOMPH_EXCEPTION_LOCATION);
    }
  }
 
 /// Update function for the upper part of the domain
 void spine_node_update_upper(PerturbedSpineNode *spine_node_pt)
  {
   // Get fraction along the spine
   const double w = spine_node_pt->fraction();

   // Get base spine height
   const double h =
    spine_node_pt->perturbed_spine_pt()->base_spine_pt()->height();
 
   // Set the nodal y-position
   spine_node_pt->x(1) = (this->Ymin + h) + w*(this->Ymax - (this->Ymin+h));

   // Get perturbed spine height (cosine and sine parts)
   const double HC = spine_node_pt->perturbed_spine_pt()->height(0);
   const double HS = spine_node_pt->perturbed_spine_pt()->height(1);

   // Set the perturbation (cosine and sine parts) to the nodal
   // y-position. Note that YC_index and YS_index are the indices at
   // which the perturbations to the nodal y-position are stored in
   // the bulk elements
   if(YC_index>=0)
    {
     // To see analysis of why this is okay, see the pdf file
     // "perturbed_node_update_procedure.pdf" (emailed to myself)
     spine_node_pt->set_value(YC_index,HC*(1.0-w));
    }
   else
    {
     std::ostringstream error_message;
     error_message << "The cosine part of the perturbation to the nodal\n"
                   << "y-position has not been set up." << std::endl;
     throw OomphLibError(
      error_message.str(),
      "TwoLayerPerturbedSpineMesh::spine_node_update_upper(...)",
      OOMPH_EXCEPTION_LOCATION);
    }
   if(YS_index>=0)
    {
     // To see analysis of why this is okay, see the pdf file
     // "perturbed_node_update_procedure.pdf" (emailed to myself)
     spine_node_pt->set_value(YS_index,HS*(1.0-w));
    }
   else
    {
     std::ostringstream error_message;
     error_message << "The sine part of the perturbation to the nodal\n"
                   << "y-position has not been set up." << std::endl;
     throw OomphLibError(
      error_message.str(),
      "TwoLayerPerturbedSpineMesh::spine_node_update_upper(...)",
      OOMPH_EXCEPTION_LOCATION);
    }
  }

 /// Helper function to actually build the two-layer spine mesh 
 /// (called from various constructors)
 virtual void build_two_layer_mesh(TimeStepper* time_stepper_pt);

  private:

 /// Static "magic" number that indicates that the indices at
 /// which the perturbations to the nodal y-positions are stored have
 /// not been set up
 static int Perturbation_to_nodal_position_indices_not_set_up;
};


} 

#endif

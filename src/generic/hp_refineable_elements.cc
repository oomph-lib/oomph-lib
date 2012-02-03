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
//Non-inline member functions for hp-refineable elements

//oomph-lib includes
#include "hp_refineable_elements.h"
#include "shape.h"

namespace oomph
{

////////////////////////////////////////////////////////////////
//       1D PRefineableQElements
////////////////////////////////////////////////////////////////

/// Get local coordinates of node j in the element; vector sets its own size
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
local_coordinate_of_node(const unsigned& n, Vector<double>& s)
 {
  s.resize(1);
  
  switch(this->nnode_1d())
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<2>::nodal_position(n);
    break;
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<3>::nodal_position(n);
    break;
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<4>::nodal_position(n);
    break;
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<5>::nodal_position(n);
    break;
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<6>::nodal_position(n);
    break;
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<7>::nodal_position(n);
    break;
  default:
    oomph_info << "\n ERROR: Exceeded maximum polynomial order for";
    oomph_info << "\n        shape functions." << std::endl;
    break;
   }
 }

/// Get the local fractino of node j in the element
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
local_fraction_of_node(const unsigned &n, Vector<double> &s_fraction)
 {
  this->local_coordinate_of_node(n,s_fraction);
  s_fraction[0] = 0.5*(s_fraction[0] + 1.0);
 }

template<unsigned INITIAL_NNODE_1D>
double PRefineableQElement<1,INITIAL_NNODE_1D>::
local_one_d_fraction_of_node(const unsigned &n1d, const unsigned &i)
 {
  switch(this->nnode_1d())
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
  default:
    std::ostringstream error_message;
    error_message <<"\nERROR: Exceeded maximum polynomial order for";
    error_message <<"\n       shape functions.\n";
    throw OomphLibError(error_message.str(),
                  "PRefineableQElement<1,INITIAL_NNODE_1D>::local_one_d_fraction_of_node()",
                  OOMPH_EXCEPTION_LOCATION);
    return 0.0;
   }
 }


//==================================================================
/// Return the node at the specified local coordinate
//==================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<1,INITIAL_NNODE_1D>::
get_node_at_local_coordinate(const Vector<double> &s)
{
 //Load the tolerance into a local variable
 double tol = FiniteElement::Node_location_tolerance;
 //There is one possible index.
 Vector<int> index(1);

 // Determine the index
 // -------------------
 
 // If we are at the lower limit, the index is zero
 if(std::abs(s[0] + 1.0) < tol)
  {
   index[0] = 0;
  }
 // If we are at the upper limit, the index is the number of nodes minus 1
 else if(std::abs(s[0] - 1.0) < tol)
  {
   index[0] = this->nnode_1d()-1;
  }
 // Otherwise, we have to calculate the index in general
 else
  {
   // Compute Gauss-Lobatto-Legendre node positions
   Vector<double> z;
   Orthpoly::gll_nodes(this->nnode_1d(), z);
   // Loop over possible internal nodal positions
   for (unsigned n=1; n<this->nnode_1d()-1; n++)
    {
     if (std::abs(z[n] - s[0]) < tol)
      {
       index[0] = n;
       break;
      }
    }
   // No matching nodes
   return 0;
  }
   // If we've got here we have a node, so let's return a pointer to it
   return this->node_pt(index[0]);
}

//===================================================================
/// If a neighbouring element has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0). If the node is
/// periodic the flag is_periodic will be true
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<1,INITIAL_NNODE_1D>::
node_created_by_neighbour(const Vector<double> &s_fraction, 
                          bool &is_periodic) 
{
 return RefineableQElement<1>::node_created_by_neighbour(s_fraction, is_periodic);
}

//===================================================================
/// If a neighbouring element's son has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0). If the node is
/// periodic the flag is_periodic will be true
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<1,INITIAL_NNODE_1D>::
node_created_by_son_of_neighbour(const Vector<double> &s_fraction, 
                                 bool &is_periodic) 
{
 // Not possible in 1D case, so return null pointer
 return 0;
}

//==================================================================
/// Set the correct p-order of the element based on that of its
/// father. Then construct an integration scheme of the correct order.
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::initial_setup()
{
 // Check if element is in a tree
 if (Tree_pt!=0)
  {
   //Pointer to my father (in quadtree impersonation)
   BinaryTree* father_pt = dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());
   
   // Check if element has father
   if (father_pt!=0)
    {
     PRefineableQElement<1,INITIAL_NNODE_1D>* father_el_pt =
            dynamic_cast<PRefineableQElement<1,INITIAL_NNODE_1D>*>
                (this->tree_pt()->father_pt()->object_pt());
     if (father_el_pt!=0)
      {
       unsigned father_p_order = father_el_pt->p_order();
       // Set p-order to that of father
       P_order = father_p_order;
      }
     
     // Now sort out the element...
     // (has p nodes)
     unsigned new_n_node = P_order;
     
     // Allocate new space for Nodes (at the element level)
     this->set_n_node(new_n_node);
     
     // Set integration scheme
     delete this->integral_pt();
     switch(P_order)
     {
     case 2:
      this->set_integration_scheme(new GaussLobattoLegendre<1,2>);
      break;
     case 3:
      this->set_integration_scheme(new GaussLobattoLegendre<1,3>);
      break;
     case 4:
      this->set_integration_scheme(new GaussLobattoLegendre<1,4>);
      break;
     case 5:
      this->set_integration_scheme(new GaussLobattoLegendre<1,5>);
      break;
     case 6:
      this->set_integration_scheme(new GaussLobattoLegendre<1,6>);
      break;
     case 7:
      this->set_integration_scheme(new GaussLobattoLegendre<1,7>);
      break;
     default:
      std::ostringstream error_message;
      error_message <<"\nERROR: Exceeded maximum polynomial order for";
      error_message <<"\n       integration scheme.\n";
      throw OomphLibError(error_message.str(),
                          "PRefineableQElement<1,INITIAL_NNODE_1D>::initial_setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    }
  }
 else
  {
   throw OomphLibError("Element not in a tree!",
                       "PRefineableQElement<1,INITIAL_NNODE_1D>::initial_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//==================================================================
/// Set the correct p-order of the element before it is built
/// and then check the father element for any required nodes which
/// already exist
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::pre_build()
{
 /*
 //Pointer to my father (in binarytree impersonation)
 BinaryTree* father_pt = dynamic_cast<BinaryTree*>(binary_tree_pt()->father_pt());
 
 // Check if element has father
 if (father_pt!=0)
  {
   PRefineableQElement<1>* father_el_pt =
          dynamic_cast<PRefineableQElement<1>*>
              (this->tree_pt()->father_pt()->object_pt());
   if (father_el_pt!=0)
    {
     // Pre-build actions
     //??
    }
   else
    {
     std::ostringstream error_message;
     error_message <<"\nERROR: Dynamic cast failed!\n";
     throw OomphLibError(error_message.str(),
                         "PRefineableQElement<1,INITIAL_NNODE_1D>::pre_build()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 */
}

//==================================================================
/// p-refine the element inc times. (If inc<0 then p-unrefinement
/// is performed.)
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::p_refine(const int &inc,
                                                       Mesh* &mesh_pt)
{
 // Timestepper should be the same for all nodes -- use it
 // to create timesteppers for new nodes
 TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
 
 // Number of history values (incl. present)
 unsigned ntstorage = time_stepper_pt->ntstorage();
 
 // Back up old vertex nodes
 unsigned n_vertex_node = this->nvertex_node();
 Vector<Node*> old_vertex_node_pt(n_vertex_node);
 for (unsigned n=0; n<n_vertex_node; n++)
  {
   old_vertex_node_pt[n] = this->vertex_node_pt(n);
  }
 
 // Compute new coordinates and projected values
 Vector<double> new_local_x(P_order + inc);
 Vector<Vector<double> > new_global_x(P_order + inc);
 for (unsigned i=0; i<P_order + inc; i++)
  {new_global_x[i].resize(1);}
 unsigned ncont = this->ncont_interpolated_values();
 Vector<Vector<Vector<double> > > projected_value(ntstorage);
 for(unsigned t=0;t<ntstorage;t++)
  {
   projected_value[t].resize(P_order + inc);
   for (unsigned i=0; i<P_order + inc; i++)
    {
     projected_value[t][i].resize(ncont);
    }
  }
 
 // Compute 1D Gauss-Lobatto-Legendre node spacing
 Orthpoly::gll_nodes(P_order + inc, new_local_x);
 
 for (unsigned n=0; n<P_order+inc; n++)
  {
   // Create coordinate vector
   Vector<double> s(1);
   s[0] = new_local_x[n];
   
   new_global_x[n][0] = this->interpolated_x(s,0);
   
   // Loop over all history values
   for(unsigned t=0;t<ntstorage;t++)
    {
     // Interpolate new nodal values
     // (while still using old integration scheme and shape functions)
     Vector<double> prev_values;
     this->get_interpolated_values(t,s,prev_values);
     for(unsigned i=0; i<ncont; i++)
      {
       projected_value[t][n][i] = prev_values[i];
      }
    }
 }
 
 // Increment p-order of the element
 p_order() += inc;
 
 // Change integration scheme
 delete this->integral_pt();
 switch(p_order())
 {
 case 2:
  this->set_integration_scheme(new GaussLobattoLegendre<1,2>);
  break;
 case 3:
  this->set_integration_scheme(new GaussLobattoLegendre<1,3>);
  break;
 case 4:
  this->set_integration_scheme(new GaussLobattoLegendre<1,4>);
  break;
 case 5:
  this->set_integration_scheme(new GaussLobattoLegendre<1,5>);
  break;
 case 6:
  this->set_integration_scheme(new GaussLobattoLegendre<1,6>);
  break;
 case 7:
  this->set_integration_scheme(new GaussLobattoLegendre<1,7>);
  break;
 default:
  oomph_info << "\n ERROR: Exceeded maximum polynomial order for";
  oomph_info << "\n        integration scheme." << std::endl;
  break;
 }
 
 // Allocate new space for Nodes (at the element level)
 this->set_n_node(P_order);
 
 // Copy back vertex nodes
 this->node_pt(0        ) = old_vertex_node_pt[0];
 this->node_pt(P_order-1) = old_vertex_node_pt[1];

 //Create new internal nodes
 for (unsigned n=1; n<P_order-1; n++)
  {
   // Build node
   Node* created_node_pt = this->construct_node(n,time_stepper_pt);
   this->node_pt(n) = created_node_pt;
   // Add node to mesh
   if (mesh_pt!=0)
    {
     mesh_pt->add_node_pt(created_node_pt);
    }
  }
 
 // Set coordinates and project data
 for(unsigned t=0;t<ntstorage;t++)
  {
   for (unsigned n=0; n<P_order; n++)
    {
     this->node_pt(n)->x(0) = new_global_x[n][0];
     for(unsigned i=0; i<ncont; i++)
      {
       this->node_pt(n)->set_value(t,i,projected_value[t][n][i]);
      }
    }
  }
}

//=======================================================================
///Shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
shape(const Vector<double> &s, Shape &psi) const
{
 switch(p_order())
 {
 case 2:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<2> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 case 3:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<3> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 case 4:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<4> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 case 5:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<5> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 case 6:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<6> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 case 7:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  // Create one-dim shape functions
  OneDimensionalLegendreShape<7> psi1(s[0]);
  //Now let's loop over the nodal points in the element
  //and copy the values back in
  for(unsigned i=0;i<p_order();i++) {psi(i) = psi1[i];}
  break;
 }
 default:
  oomph_info << "\n ERROR: PRefineableQElement::shape() exceeded maximum";
  oomph_info << "\n        polynomial order for shape functions." << std::endl;
 }
 
}

//=======================================================================
///Derivatives of shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsi) const
{
 switch(p_order())
 {
 case 2:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<2> psi1(s[0]);
  OneDimensionalLegendreDShape<2> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  
  break;
 }
 case 3:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<3> psi1(s[0]);
  OneDimensionalLegendreDShape<3> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  break;
 }
 case 4:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<4> psi1(s[0]);
  OneDimensionalLegendreDShape<4> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  break;
 }
 case 5:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<5> psi1(s[0]);
  OneDimensionalLegendreDShape<5> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  break;
 }
 case 6:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<6> psi1(s[0]);
  OneDimensionalLegendreDShape<6> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  break;
 }
 case 7:
 {
  // Calculate nodal positions
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<7> psi1(s[0]);
  OneDimensionalLegendreDShape<7> dpsi1ds(s[0]);
  // Loop over shapes and copy across
  for(unsigned i=0;i<p_order();i++) 
   {
    psi(i) = psi1[i];
    dpsi(i,0) = dpsi1ds[i];
   }
  break;
 }
 default:
  oomph_info << "\n ERROR: PRefineableQElement::dshape_local() exceeded maximum";
  oomph_info << "\n        polynomial order for shape functions." << std::endl;
 }
 
}

//=======================================================================
/// Second derivatives of shape functions for PRefineableQElement<DIM>
/// \n d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
d2shape_local(const Vector<double> &s, Shape &psi, DShape &dpsids,
              DShape &d2psids) const
{
 std::ostringstream error_message;
 error_message <<"\nd2shape_local currently not implemented for this element\n";
 throw OomphLibError(error_message.str(),
                     "PRefineableQElement<1,INITIAL_NNODE_1D>::d2shape_local()",
                     OOMPH_EXCEPTION_LOCATION);
}

//=================================================================
/// Internal function to set up the hanging nodes on a particular
/// edge of the element. (Not required in 1D.)
//=================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
binary_hang_helper(const int &value_id, const int &my_edge,
                   std::ofstream& output_hangfile)
{}

//=================================================================
/// Check inter-element continuity of 
/// - nodal positions
/// - (nodally) interpolated function values
//==================================================================== 
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<1,INITIAL_NNODE_1D>::
check_integrity(double& max_error)
{
 RefineableQElement<1>::check_integrity(max_error);
}

////////////////////////////////////////////////////////////////
//       2D PRefineableQElements
////////////////////////////////////////////////////////////////

/// Get local coordinates of node j in the element; vector sets its own size
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
local_coordinate_of_node(const unsigned& n, Vector<double>& s)
 {
  s.resize(2);
  unsigned Nnode_1d = this->nnode_1d();
  unsigned n0 = n%Nnode_1d;
  unsigned n1 = unsigned(double(n)/double(Nnode_1d));
  
  switch(Nnode_1d)
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<2>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<2>::nodal_position(n1);
    break;
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<3>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<3>::nodal_position(n1);
    break;
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<4>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<4>::nodal_position(n1);
    break;
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<5>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<5>::nodal_position(n1);
    break;
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<6>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<6>::nodal_position(n1);
    break;
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<7>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<7>::nodal_position(n1);
    break;
  default:
    std::ostringstream error_message;
    error_message <<"\nERROR: Exceeded maximum polynomial order for";
    error_message <<"\n       shape functions.\n";
    throw OomphLibError(error_message.str(),
                        "PRefineableQElement<2,INITIAL_NNODE_1D>::local_coordinate_of_node()",
                        OOMPH_EXCEPTION_LOCATION);
    break;
   }
 }
 
/// Get the local fractino of node j in the element
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
local_fraction_of_node(const unsigned &n, Vector<double> &s_fraction)
 {
  this->local_coordinate_of_node(n,s_fraction);
  s_fraction[0] = 0.5*(s_fraction[0] + 1.0);
  s_fraction[1] = 0.5*(s_fraction[1] + 1.0);
 }

template<unsigned INITIAL_NNODE_1D>
double PRefineableQElement<2,INITIAL_NNODE_1D>::
local_one_d_fraction_of_node(const unsigned &n1d, const unsigned &i)
 {
  switch(this->nnode_1d())
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
  default:
    std::ostringstream error_message;
    error_message <<"\nERROR: Exceeded maximum polynomial order for";
    error_message <<"\n       shape functions.\n";
    throw OomphLibError(error_message.str(),
                  "PRefineableQElement<2,INITIAL_NNODE_1D>::local_one_d_fraction_of_node()",
                  OOMPH_EXCEPTION_LOCATION);
    return 0.0;
   }
 }


//==================================================================
/// Return the node at the specified local coordinate
//==================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<2,INITIAL_NNODE_1D>::
get_node_at_local_coordinate(const Vector<double> &s)
{
 unsigned Nnode_1d = this->nnode_1d();
 //Load the tolerance into a local variable
 double tol = FiniteElement::Node_location_tolerance;
 //There are two possible indices.
 Vector<int> index(2);
 
 // Loop over indices
 for (unsigned i=0; i<2; i++)
  {
   // Determine the index
   // -------------------
   
   bool is_found=false;
   
   // If we are at the lower limit, the index is zero
   if(std::abs(s[i] + 1.0) < tol)
    {
     index[i] = 0;
     is_found=true;
    }
   // If we are at the upper limit, the index is the number of nodes minus 1
   else if(std::abs(s[i] - 1.0) < tol)
    {
     index[i] = Nnode_1d-1;
     is_found=true;
    }
   // Otherwise, we have to calculate the index in general
   else
    {
     // Compute Gauss-Lobatto-Legendre node positions
     Vector<double> z;
     Orthpoly::gll_nodes(Nnode_1d, z);
     // Loop over possible internal nodal positions
     for (unsigned n=1; n<Nnode_1d-1; n++)
      {
       if (std::abs(z[n] - s[i]) < tol)
        {
         index[i] = n;
         is_found=true;
         break;
        }
      }
    }
   
   if (!is_found)
    {
     // No matching nodes
     return 0;
    }
  }
 // If we've got here we have a node, so let's return a pointer to it
 return this->node_pt(index[0] + Nnode_1d*index[1]);
}

//===================================================================
/// If a neighbouring element has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0). If the node is
/// periodic the flag is_periodic will be true
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<2,INITIAL_NNODE_1D>::
node_created_by_neighbour(const Vector<double> &s_fraction, bool &is_periodic) 
{  
 using namespace QuadTreeNames;

 //Calculate the edges on which the node lies
 Vector<int> edges;
 if(s_fraction[0]==0.0) {edges.push_back(W);}
 if(s_fraction[0]==1.0) {edges.push_back(E);}
 if(s_fraction[1]==0.0) {edges.push_back(S);}
 if(s_fraction[1]==1.0) {edges.push_back(N);}

 //Find the number of edges
 unsigned n_size = edges.size();
 //If there are no edges, then there is no neighbour, return 0
 if(n_size==0) {return 0;}

 Vector<unsigned> translate_s(2);
 Vector<double> s_lo_neigh(2);
 Vector<double> s_hi_neigh(2);
 Vector<double> s(2);

 int neigh_edge, diff_level;
 bool in_neighbouring_tree;

 //Loop over the edges
 for(unsigned j=0;j<n_size;j++)
  {
   // Find pointer to neighbouring element along edge 
   QuadTree* neigh_pt;
   neigh_pt = quadtree_pt()->
    gteq_edge_neighbour(edges[j],translate_s,s_lo_neigh,s_hi_neigh,
                        neigh_edge,diff_level,in_neighbouring_tree);
   
   // Neighbour exists
   if(neigh_pt!=0)
    {
     // Have any of its vertex nodes been created yet?
     // (BS: Must look in incomplete neighbours because after the pre-build
     // they may contain pointers to the required nodes. e.g. h-refinement of
     // neighbouring linear and quadratic elements)
     bool a_vertex_node_is_built = false;
     QElement<2,INITIAL_NNODE_1D>* neigh_obj_pt =
      dynamic_cast<QElement<2,INITIAL_NNODE_1D>*>(neigh_pt->object_pt());
     if(neigh_obj_pt==0)
      {
       throw
        OomphLibError("Not a quad element!",
         "PRefineableQElement<2,INITIAL_NNODE_1D>::node_created_by_neighbour()",
         OOMPH_EXCEPTION_LOCATION);
      }
     for(unsigned vnode=0; vnode<neigh_obj_pt->nvertex_node(); vnode++)
      {
       if(neigh_obj_pt->vertex_node_pt(vnode)!=0)
        a_vertex_node_is_built = true;
       break;
      }
     if(a_vertex_node_is_built)
      {
       //We now need to translate the nodal location
       //as defined in terms of the fractional coordinates of the present
       //element into those of its neighbour
       
       //Calculate the local coordinate in the neighbour
       //Note that we need to use the translation scheme in case
       //the local coordinates are swapped in the neighbour.
       for(unsigned i=0;i<2;i++)
        {
         s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]]*
          (s_hi_neigh[i] - s_lo_neigh[i]);
        }
       
       //Find the node in the neighbour
       Node* neighbour_node_pt = 
        neigh_pt->object_pt()->get_node_at_local_coordinate(s);
       
       //If there is a node, return it
       if(neighbour_node_pt!=0) 
        {
         //Now work out whether it's a periodic boundary
         //only possible if we have moved into a neighbouring tree
         if(in_neighbouring_tree)
          {
           //Return whether the neighbour is periodic 
           is_periodic = 
            quadtree_pt()->root_pt()->is_neighbour_periodic(edges[j]);
          }
         //Return the pointer to the neighbouring node
         return neighbour_node_pt;
        }
      }
    }
  }
 //Node not found, return null
 return 0;
}

//===================================================================
/// If a neighbouring element's son has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0). If the node is
/// periodic the flag is_periodic will be true
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<2,INITIAL_NNODE_1D>::
node_created_by_son_of_neighbour(const Vector<double> &s_fraction,
                                 bool &is_periodic) 
{
 using namespace QuadTreeNames;

 //Calculate the edges on which the node lies
 Vector<int> edges;
 if(s_fraction[0]==0.0) {edges.push_back(W);}
 if(s_fraction[0]==1.0) {edges.push_back(E);}
 if(s_fraction[1]==0.0) {edges.push_back(S);}
 if(s_fraction[1]==1.0) {edges.push_back(N);}

 //Find the number of edges
 unsigned n_size = edges.size();
 //If there are no edges, then there is no neighbour, return 0
 if(n_size==0) {return 0;}

 Vector<unsigned> translate_s(2);
 Vector<double> s_lo_neigh(2);
 Vector<double> s_hi_neigh(2);
 Vector<double> s(2);

 int neigh_edge, diff_level;
 bool in_neighbouring_tree;

 //Loop over the edges
 for(unsigned j=0;j<n_size;j++)
  {
   // Find pointer to neighbouring element along edge 
   QuadTree* neigh_pt;
   neigh_pt = quadtree_pt()->
    gteq_edge_neighbour(edges[j],translate_s,s_lo_neigh,s_hi_neigh,
                        neigh_edge,diff_level,in_neighbouring_tree);
   
   // Neighbour exists
   if(neigh_pt!=0)
    {
     // Have its nodes been created yet?
     // (Must look in sons of unfinished neighbours too!!!)
     if(1)
      {
       //We now need to translate the nodal location
       //as defined in terms of the fractional coordinates of the present
       //element into those of its neighbour
       
       //Calculate the local coordinate in the neighbour
       //Note that we need to use the translation scheme in case
       //the local coordinates are swapped in the neighbour.
       for(unsigned i=0;i<2;i++)
        {
         s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]]*
          (s_hi_neigh[i] - s_lo_neigh[i]);
        }
       
       // Check if the element has sons
       if(neigh_pt->nsons()!=0)
        {
         //First, find the son element in which the node should live
         
         //Find coordinates in the sons
         Vector<double> s_in_son(2);
         int son=-10;
         //If negative on the west side
         if(s[0] < 0.0)
          {
           //On the south side
           if(s[1] < 0.0)
            {
             //It's the southwest son
             son = SW;
             s_in_son[0] =  1.0 + 2.0*s[0];
             s_in_son[1] =  1.0 + 2.0*s[1];
            }
           //On the north side
           else
            {
             //It's the northwest son
             son = NW;
             s_in_son[0] =  1.0 + 2.0*s[0];
             s_in_son[1] = -1.0 + 2.0*s[1];
            }
          }
         else
          {
           //On the south side
           if(s[1] < 0.0)
            {
             //It's the southeast son
             son = SE;
             s_in_son[0] =  -1.0 + 2.0*s[0];
             s_in_son[1] =   1.0 + 2.0*s[1];
            }
           //On the north side
           else
            {
             //It's the northeast son
             son = NE;
             s_in_son[0] = -1.0 + 2.0*s[0];
             s_in_son[1] = -1.0 + 2.0*s[1];
            }
          }
         
         //Find the node in the neighbour's son
         Node* neighbour_son_node_pt = 
          neigh_pt->son_pt(son)->object_pt()->
           get_node_at_local_coordinate(s_in_son);
         
         //If there is a node, return it
         if(neighbour_son_node_pt!=0) 
          {
           //Now work out whether it's a periodic boundary
           //only possible if we have moved into a neighbouring tree
           if(in_neighbouring_tree)
            {
             //Return whether the neighbour is periodic 
             is_periodic = 
              quadtree_pt()->root_pt()->is_neighbour_periodic(edges[j]);
            }
           //Return the pointer to the neighbouring node
           return neighbour_son_node_pt;
          }
        }
       else
        {
         // No sons to search in, so no node can be found
         return 0;
        }
      }
    }
  }
 //Node not found, return null
 return 0;
}

//==================================================================
/// Set the correct p-order of the element based on that of its
/// father. Then construct an integration scheme of the correct order.
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::initial_setup()
{
 // Check if element is in a tree
 if (Tree_pt!=0)
  {
   //Pointer to my father (in quadtree impersonation)
   QuadTree* father_pt = dynamic_cast<QuadTree*>(quadtree_pt()->father_pt());
   
   // Check if element has father
   if (father_pt!=0)
    {
     PRefineableQElement<2,INITIAL_NNODE_1D>* father_el_pt =
            dynamic_cast<PRefineableQElement<2,INITIAL_NNODE_1D>*>
                (this->tree_pt()->father_pt()->object_pt());
     if (father_el_pt!=0)
      {
       unsigned father_p_order = father_el_pt->p_order();
       // Set p-order to that of father
       P_order = father_p_order;
      }
     
     // Now sort out the element...
     // (has p^2 nodes)
     unsigned new_n_node = P_order*P_order;
     
     // Allocate new space for Nodes (at the element level)
     this->set_n_node(new_n_node);
     
     // Set integration scheme
     delete this->integral_pt();
     switch(P_order)
     {
     case 2:
      this->set_integration_scheme(new GaussLobattoLegendre<2,2>);
      break;
     case 3:
      this->set_integration_scheme(new GaussLobattoLegendre<2,3>);
      break;
     case 4:
      this->set_integration_scheme(new GaussLobattoLegendre<2,4>);
      break;
     case 5:
      this->set_integration_scheme(new GaussLobattoLegendre<2,5>);
      break;
     case 6:
      this->set_integration_scheme(new GaussLobattoLegendre<2,6>);
      break;
     case 7:
      this->set_integration_scheme(new GaussLobattoLegendre<2,7>);
      break;
     default:
      std::ostringstream error_message;
      error_message <<"\nERROR: Exceeded maximum polynomial order for";
      error_message <<"\n       integration scheme.\n";
      throw OomphLibError(error_message.str(),
                          "PRefineableQElement<2,INITIAL_NNODE_1D>::initial_setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    }
  }
 else
  {
   throw OomphLibError("Element not in a tree!",
                       "PRefineableQElement<2,INITIAL_NNODE_1D>::initial_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//==================================================================
/// Set the correct p-order of the element before it is built
/// and then check the father element for any required nodes which
/// already exist
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::pre_build()
{
 //Pointer to my father (in quadtree impersonation)
 QuadTree* father_pt = dynamic_cast<QuadTree*>(quadtree_pt()->father_pt());
 
 // Check if element has father
 if (father_pt!=0)
  {
   PRefineableQElement<2,INITIAL_NNODE_1D>* father_el_pt =
          dynamic_cast<PRefineableQElement<2,INITIAL_NNODE_1D>*>
              (this->tree_pt()->father_pt()->object_pt());
   if (father_el_pt!=0)
    {
     // Check nodes from father
     //------------------------
     unsigned n_p = this->nnode_1d();
    
     // What type of son am I? Ask my quadtree representation...
     int son_type = Tree_pt->son_type();
     
     Vector<double> s_lo(2);
     Vector<double> s_hi(2);
     Vector<double> s(2);
     Vector<double> x(2);
     
     // Setup vertex coordinates in father element:
     switch(son_type)
      {
      case QuadTreeNames::SW:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       break;
       
      case QuadTreeNames::SE:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       break;
      
      case QuadTreeNames::NE:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       break;
      
      case QuadTreeNames::NW:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       break;
      }
     unsigned jnod=0;
   
     Vector<double> s_fraction(2);
     // Loop over nodes in element
     for(unsigned i0=0;i0<n_p;i0++)
      {
       //Get the fractional position of the node in the direction of s[0]
       s_fraction[0] = local_one_d_fraction_of_node(i0,0);
       // Local coordinate in father element
       s[0] = s_lo[0] + (s_hi[0]-s_lo[0])*s_fraction[0];
   
       for(unsigned i1=0;i1<n_p;i1++)
        {
         //Get the fractional position of the node in the direction of s[1]
         s_fraction[1] = local_one_d_fraction_of_node(i1,1);
         // Local coordinate in father element
         s[1] = s_lo[1] + (s_hi[1]-s_lo[1])*s_fraction[1];
         
         // Local node number
         jnod= i0 + n_p*i1;
   
         //Get the pointer to the node in the father, returns NULL
         //if there is not node
         Node* created_node_pt = father_el_pt->get_node_at_local_coordinate(s);
   
         // Does this node already exist in father element?
         //------------------------------------------------
         if(created_node_pt!=0) 
          {
           // Copy node across
           this->node_pt(jnod) = created_node_pt;
          }
        }
      }
    }
   else
    {
     std::ostringstream error_message;
     error_message <<"\nERROR: Dynamic cast failed!\n";
     throw OomphLibError(error_message.str(),
                         "PRefineableQElement<2,INITIAL_NNODE_1D>::pre_build()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
}

//==================================================================
/// p-refine the element inc times. (If inc<0 then p-unrefinement
/// is performed.)
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::p_refine(const int &inc,
                                                       Mesh* &mesh_pt)
{
 // Timestepper should be the same for all nodes -- use it
 // to create timesteppers for new nodes
 TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
 
 // Number of history values (incl. present)
 unsigned ntstorage = time_stepper_pt->ntstorage();
 
 // Back up pointers to old vertex nodes before they are lost
 unsigned n_vertex_node = this->nvertex_node();
 Vector<Node*> old_vertex_node_pt(n_vertex_node);
 for (unsigned n=0; n<n_vertex_node; n++)
  {
   old_vertex_node_pt[n] = this->vertex_node_pt(n);
  }
 
 // Compute new coordinates and projected values
 Vector<double> new_s((P_order + inc)*(P_order + inc));
 Vector<Vector<Vector<double> > > new_global_x(ntstorage);
 for(unsigned t=0;t<ntstorage;t++)
  {
   new_global_x[t].resize((P_order + inc)*(P_order + inc));
   for (unsigned n=0; n<(P_order + inc)*(P_order + inc); n++)
    {
     new_global_x[t][n].resize(2);
    }
  }
 unsigned ncont = this->ncont_interpolated_values();
 Vector<Vector<Vector<double> > > projected_value(ntstorage);
 for(unsigned t=0;t<ntstorage;t++)
  {
   projected_value[t].resize((P_order + inc)*(P_order + inc));
   for (unsigned n=0; n<(P_order + inc)*(P_order + inc); n++)
    {
     projected_value[t][n].resize(ncont);
    }
  }
 
 // Compute Gauss-Lobatto-Legendre node spacing
 Orthpoly::gll_nodes(P_order + inc, new_s);
 for (unsigned n=0; n<(P_order + inc)*(P_order + inc); n++)
  {
   // Create coordinate vector
   Vector<double> s(2);
   unsigned n0 = n%(P_order+inc);
   unsigned n1 = unsigned(double(n)/double(P_order+inc));
   s[0] = new_s[n0];
   s[1] = new_s[n1];
   
   // Loop over all history values
   for(unsigned t=0;t<ntstorage;t++)
    {
     // Interpolate new nodal positions
     // (while still using old integration scheme and shape functions)
     new_global_x[t][n][0] = this->interpolated_x(t,s,0);
     new_global_x[t][n][1] = this->interpolated_x(t,s,1);
     
     // Interpolate new nodal values
     // (while still using old integration scheme and shape functions)
     Vector<double> prev_values;
     this->get_interpolated_values(t,s,prev_values);
     for(unsigned i=0; i<ncont; i++)
      {
       projected_value[t][n][i] = prev_values[i];
      }
    }
  }
 
 // Increment p-order of the element
 P_order += inc;
 
 // Change integration scheme
 delete this->integral_pt();
 switch(P_order)
 {
 case 2:
  this->set_integration_scheme(new GaussLobattoLegendre<2,2>);
  break;
 case 3:
  this->set_integration_scheme(new GaussLobattoLegendre<2,3>);
  break;
 case 4:
  this->set_integration_scheme(new GaussLobattoLegendre<2,4>);
  break;
 case 5:
  this->set_integration_scheme(new GaussLobattoLegendre<2,5>);
  break;
 case 6:
  this->set_integration_scheme(new GaussLobattoLegendre<2,6>);
  break;
 case 7:
  this->set_integration_scheme(new GaussLobattoLegendre<2,7>);
  break;
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       integration scheme.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<2,INITIAL_NNODE_1D>::p_refine()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
 // Allocate new space for Nodes (at the element level)
 this->set_n_node(P_order*P_order);
 unsigned new_Nnode = this->nnode();
 
 // Copy vertex nodes and create new edge and internal nodes
 //---------------------------------------------------------
 
 // Copy vertex nodes
 this->node_pt(0                  ) = old_vertex_node_pt[0];
 this->node_pt(P_order-1          ) = old_vertex_node_pt[1];
 this->node_pt(P_order*(P_order-1)) = old_vertex_node_pt[2];
 this->node_pt(P_order*P_order-1  ) = old_vertex_node_pt[3];


 //Max and min local coordinates in element
 Vector<double> s_lo(2);
 Vector<double> s_hi(2);
 s_lo[0]=-1.0;
 s_hi[0]= 1.0;
 s_lo[1]=-1.0;
 s_hi[1]= 1.0;

 //Local coordinate in element
 Vector<double> s(2);
 
 unsigned jnod=0;
 
 Vector<double> s_fraction(2);
 // Loop over nodes in element
 for(unsigned i0=0;i0<P_order;i0++)
  {
   //Get the fractional position of the node in the direction of s[0]
   s_fraction[0] = local_one_d_fraction_of_node(i0,0);
   // Local coordinate
   s[0] = s_lo[0] + (s_hi[0]-s_lo[0])*s_fraction[0];
   
   for(unsigned i1=0;i1<P_order;i1++)
    {
     //Get the fractional position of the node in the direction of s[1]
     s_fraction[1] = local_one_d_fraction_of_node(i1,1);
     // Local coordinate
     s[1] = s_lo[1] + (s_hi[1]-s_lo[1])*s_fraction[1];
     
     // Local node number
     jnod= i0 + P_order*i1;
     
     // Initialise flag: So far, this node hasn't been built
     // or copied yet
     bool node_done=false; 
     
     //Get the pointer to the node in this element, returns NULL
     //if there is not node
     Node* created_node_pt = this->get_node_at_local_coordinate(s);

     // Does this node already exist in this element?
     //----------------------------------------------
     if (created_node_pt!=0)
      {
       // Copy node across
       this->node_pt(jnod) = created_node_pt;
       
       // Node has been created by copy
       node_done = true;
      }
     // Node does not exist in this element but might already
     //------------------------------------------------------
     // have been created by neighbouring elements
     //-------------------------------------------
     else
      {
       //Was the node created by one of its neighbours
       //Whether or not the node lies on an edge can be calculated
       //by from the fractional position
       bool is_periodic = false;
       created_node_pt =
        node_created_by_neighbour(s_fraction,is_periodic);

       //If the node was so created, assign the pointers
       if (created_node_pt!=0)
        {
         //If the node is periodic
         if(is_periodic)
          {
           //Now the node must be on a boundary, but we don't know which
           //one
           //The returned created_node_pt is actually the neighbouring
           //periodic node
           Node* neighbour_node_pt = created_node_pt;
           
           // Determine the edge on which the new node will live
           //(cannot be a vertex node)
           int my_bound = Tree::OMEGA;
           if(s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
           else if(s_fraction[0] == 1.0) my_bound = QuadTreeNames::E;
           else if(s_fraction[1] == 0.0) my_bound = QuadTreeNames::S;
           else if(s_fraction[1] == 1.0) my_bound = QuadTreeNames::N;
           
           // Storage for the set of Mesh boundaries on which the 
           // appropriate edge lives.
           // [New nodes should always be mid-edge nodes and therefore
           //only live on one boundary but just to play it safe...]
           std::set<unsigned> boundaries;
           //Only get the boundaries if we are at the edge of
           //an element. Nodes in the centre of an element cannot be
           //on Mesh boundaries
           if(my_bound!=Tree::OMEGA)
            {this->get_boundaries(my_bound,boundaries);}
           
#ifdef PARANOID
           //Case where a new node lives on more than one boundary
           // seems fishy enough to flag
           if (boundaries.size()>1)
            {
             throw OomphLibError(
              "boundaries.size()!=1 seems a bit strange..\n",
              "RefineableQElement<2>::build()",
              OOMPH_EXCEPTION_LOCATION);
            }
         
           //Case when there are no boundaries, we are in big trouble
           if(boundaries.size() == 0)
            {
             std::ostringstream error_stream;
             error_stream    
              << "Periodic node is not on a boundary...\n"
              << "Coordinates: " 
              << created_node_pt->x(0) << " "
              << created_node_pt->x(1) << "\n";
             throw OomphLibError(
              error_stream.str(),
              "RefineableQElement<2>::build()",
              OOMPH_EXCEPTION_LOCATION);
            }
#endif
           
           // Create node and set the pointer to it from the element 
           created_node_pt = 
            construct_boundary_node(jnod,time_stepper_pt);
           //Make the node periodic from the neighbour
           created_node_pt->
            make_periodic(neighbour_node_pt);
           
           // Next, we Update the boundary lookup schemes
           //Loop over the boundaries stored in the set
           for(std::set<unsigned>::iterator it = boundaries.begin();
               it != boundaries.end(); ++it)
            {
             //Add the node to the boundary
             mesh_pt->add_boundary_node(*it,created_node_pt);
             
             //If we have set an intrinsic coordinate on this
             //mesh boundary then it must also be interpolated on
             //the new node
             //Now interpolate the intrinsic boundary coordinate
             if(mesh_pt->boundary_coordinate_exists(*it)==true)
              {
               Vector<double> zeta(1);
               this->interpolated_zeta_on_edge(*it,
                                               my_bound,
                                               s,zeta);
               
               created_node_pt->set_coordinates_on_boundary(*it,zeta);
              }
            }
           
           //Make sure that we add the node to the mesh
           mesh_pt->add_node_pt(created_node_pt);        
          } //End of periodic case
         //Otherwise the node is not periodic, so just set the 
         //pointer to the neighbours node
         else
          {
           node_pt(jnod) = created_node_pt;
          }
         //Node has been created
         node_done = true;
        }
       else
        {
         //Was the node created by one of its neighbours' sons
         //Whether or not the node lies on an edge can be calculated
         //by from the fractional position
         bool is_periodic = false;
         created_node_pt =
          node_created_by_son_of_neighbour(s_fraction,is_periodic);
         
         //If the node was so created, assign the pointers
         if (created_node_pt!=0)
          {
           //If the node is periodic
           if(is_periodic)
            {
             //Now the node must be on a boundary, but we don't know which
             //one
             //The returned created_node_pt is actually the neighbouring
             //periodic node
             Node* neighbour_node_pt = created_node_pt;
             
             // Determine the edge on which the new node will live
             //(cannot be a vertex node)
             int my_bound = Tree::OMEGA;
             if(s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
             else if(s_fraction[0] == 1.0) my_bound = QuadTreeNames::E;
             else if(s_fraction[1] == 0.0) my_bound = QuadTreeNames::S;
             else if(s_fraction[1] == 1.0) my_bound = QuadTreeNames::N;
             
             // Storage for the set of Mesh boundaries on which the 
             // appropriate edge lives.
             // [New nodes should always be mid-edge nodes and therefore
             //only live on one boundary but just to play it safe...]
             std::set<unsigned> boundaries;
             //Only get the boundaries if we are at the edge of
             //an element. Nodes in the centre of an element cannot be
             //on Mesh boundaries
             if(my_bound!=Tree::OMEGA)
              {this->get_boundaries(my_bound,boundaries);}
             
#ifdef PARANOID
             //Case where a new node lives on more than one boundary
             // seems fishy enough to flag
             if (boundaries.size()>1)
              {
               throw OomphLibError(
                "boundaries.size()!=1 seems a bit strange..\n",
                "RefineableQElement<2>::build()",
                OOMPH_EXCEPTION_LOCATION);
              }
             
             //Case when there are no boundaries, we are in big trouble
             if(boundaries.size() == 0)
              {
               std::ostringstream error_stream;
               error_stream    
                << "Periodic node is not on a boundary...\n"
                << "Coordinates: " 
                << created_node_pt->x(0) << " "
                << created_node_pt->x(1) << "\n";
               throw OomphLibError(
                error_stream.str(),
                "RefineableQElement<2>::build()",
                OOMPH_EXCEPTION_LOCATION);
              }
#endif
             
             // Create node and set the pointer to it from the element 
             created_node_pt = 
              construct_boundary_node(jnod,time_stepper_pt);
             //Make the node periodic from the neighbour
             created_node_pt->
              make_periodic(neighbour_node_pt);
             
             // Next, we Update the boundary lookup schemes
             //Loop over the boundaries stored in the set
             for(std::set<unsigned>::iterator it = boundaries.begin();
                 it != boundaries.end(); ++it)
              {
               //Add the node to the boundary
               mesh_pt->add_boundary_node(*it,created_node_pt);
               
               //If we have set an intrinsic coordinate on this
               //mesh boundary then it must also be interpolated on
               //the new node
               //Now interpolate the intrinsic boundary coordinate
               if(mesh_pt->boundary_coordinate_exists(*it)==true)
                {
                 Vector<double> zeta(1);
                 this->interpolated_zeta_on_edge(*it,
                                                 my_bound,
                                                 s,zeta);
                 
                 created_node_pt->set_coordinates_on_boundary(*it,zeta);
                }
              }
             
             //Make sure that we add the node to the mesh
             mesh_pt->add_node_pt(created_node_pt);        
            } //End of periodic case
           //Otherwise the node is not periodic, so just set the 
           //pointer to the neighbours node
           else
            {
             node_pt(jnod) = created_node_pt;
            }
           //Node has been created
           node_done = true;

#ifdef PARANOID
           oomph_info << "% Node exists in son of neighbouring element: "
                      << created_node_pt << std::endl;
           oomph_info << "% This functionality has not yet been fully tested."
                      << std::endl;
#endif
          } //Node does not exist in son of neighbouring element
        } //Node does not exist in neighbouring element
      } // Node does not exist in this element
     
     // Node has not been built anywhere ---> build it here
     if (!node_done)
      {
       //Firstly, we need to determine whether or not a node lies
       //on the boundary before building it, because 
       //we actually assign a different type of node on boundaries.

       // Determine the edge on which the new node will live
       //(cannot be a vertex node)
       int my_bound = Tree::OMEGA;
       if(s_fraction[0] == 0.0) my_bound = QuadTreeNames::W;
       else if(s_fraction[0] == 1.0) my_bound = QuadTreeNames::E;
       else if(s_fraction[1] == 0.0) my_bound = QuadTreeNames::S;
       else if(s_fraction[1] == 1.0) my_bound = QuadTreeNames::N;
           
       // Storage for the set of Mesh boundaries on which the 
       // appropriate edge lives.
       // [New nodes should always be mid-edge nodes and therefore
       //only live on one boundary but just to play it safe...]
       std::set<unsigned> boundaries;
       //Only get the boundaries if we are at the edge of
       //an element. Nodes in the centre of an element cannot be
       //on Mesh boundaries
       if(my_bound!=Tree::OMEGA)
        {this->get_boundaries(my_bound,boundaries);}
           
#ifdef PARANOID
       //Case where a new node lives on more than one boundary
       // seems fishy enough to flag
       if (boundaries.size()>1)
        {
         throw OomphLibError(
                             "boundaries.size()!=1 seems a bit strange..\n",
                             "RefineableQElement<2>::build()",
                             OOMPH_EXCEPTION_LOCATION);
        }
#endif
           
       //If the node lives on a mesh boundary, 
       //then we need to create a boundary node
       if(boundaries.size()> 0)
        {
         // Create node and set the pointer to it from the element 
         created_node_pt = construct_boundary_node(jnod,time_stepper_pt);

         //Now we need to work out whether to pin the values at
         //the new node based on the boundary conditions applied at
         //its Mesh boundary 

         //Get the boundary conditions from the father
         Vector<int> bound_cons(ncont_interpolated_values());
         this->get_bcs(my_bound,bound_cons);
             
         //Loop over the values and pin, if necessary
         unsigned n_value = created_node_pt->nvalue();
         for(unsigned k=0;k<n_value;k++)
          {
           if (bound_cons[k]) {created_node_pt->pin(k);}
          }
             
         // Solid node? If so, deal with the positional boundary
         // conditions:
         SolidNode* solid_node_pt = 
          dynamic_cast<SolidNode*>(created_node_pt);
         if (solid_node_pt!=0)
          {
           //Get the positional boundary conditions from the father:
           unsigned n_dim = created_node_pt->ndim();
           Vector<int> solid_bound_cons(n_dim);
           RefineableSolidQElement<2>* solid_el_pt=
            dynamic_cast<RefineableSolidQElement<2>*>(this);
#ifdef PARANOID
           if (solid_el_pt==0)
            {
             std::string error_message =
              "We have a SolidNode outside a refineable SolidElement\n";
             error_message +=
              "during mesh refinement -- this doesn't make sense";

             throw OomphLibError(error_message,
                                 "RefineableQElement<2>::build()",
                                 OOMPH_EXCEPTION_LOCATION);
            }
#endif
           solid_el_pt->
            get_solid_bcs(my_bound,solid_bound_cons);
               
           //Loop over the positions and pin, if necessary
           for(unsigned k=0;k<n_dim;k++)
            {
             if (solid_bound_cons[k]) {solid_node_pt->pin_position(k);}
            }
          } //End of if solid_node_pt

             

         // Next, we Update the boundary lookup schemes
         //Loop over the boundaries stored in the set
         for(std::set<unsigned>::iterator it = boundaries.begin();
             it != boundaries.end(); ++it)
          {
           //Add the node to the boundary
           mesh_pt->add_boundary_node(*it,created_node_pt);
               
           //If we have set an intrinsic coordinate on this
           //mesh boundary then it must also be interpolated on
           //the new node
           //Now interpolate the intrinsic boundary coordinate
           if(mesh_pt->boundary_coordinate_exists(*it)==true)
            {
             Vector<double> zeta(1);
             this->interpolated_zeta_on_edge(*it,
                                             my_bound,
                                             s,zeta);

             created_node_pt->set_coordinates_on_boundary(*it,zeta);
            }
          }
        }
       //Otherwise the node is not on a Mesh boundary and
       //we create a normal "bulk" node
       else
        {
         // Create node and set the pointer to it from the element 
         created_node_pt = construct_node(jnod,time_stepper_pt);
        }
           
       // Add new node to mesh
       mesh_pt->add_node_pt(created_node_pt);
           
      } //End of case when we build the node ourselves
    }
  }
 
 // Set coordinates and project data
 for(unsigned t=0;t<ntstorage;t++)
  {
   for (unsigned n=0; n<new_Nnode; n++)
    {
     for(unsigned i=0; i<2; i++)
      {
       this->node_pt(n)->x(t,i) = new_global_x[t][n][i];
      }
     for(unsigned i=0; i<ncont; i++)
      {
       this->node_pt(n)->set_value(t,i,projected_value[t][n][i]);
      }
    }
  }
 
 // Not necessary to delete the old nodes since all original nodes are in the
 // current mesh and so will be pruned as part of the mesh adaption process.
 
 // Do any further-build required
 this->further_build();
 
}

//=======================================================================
///Shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
shape(const Vector<double> &s, Shape &psi) const
{
 switch(p_order())
 {
 case 2:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<2;i++) 
   {
    for(unsigned j=0;j<2;j++)
     {
      psi(2*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 case 3:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<3;i++) 
   {
    for(unsigned j=0;j<3;j++)
     {
      psi(3*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 case 4:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<4;i++) 
   {
    for(unsigned j=0;j<4;j++)
     {
      psi(4*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 case 5:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<5;i++) 
   {
    for(unsigned j=0;j<5;j++)
     {
      psi(5*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 case 6:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<6;i++) 
   {
    for(unsigned j=0;j<6;j++)
     {
      psi(6*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 case 7:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<7;i++) 
   {
    for(unsigned j=0;j<7;j++)
     {
      psi(7*i + j) = psi2[i]*psi1[j];
     }
   }
  break;
 }
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       polynomial order for shape functions.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<2,INITIAL_NNODE_1D>::shape()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
}

//=======================================================================
///Derivatives of shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsids) const
{
 switch(p_order())
 {
 case 2:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<2> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<2;i++)
   {
    for(unsigned j=0;j<2;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 case 3:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<3> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<3;i++)
   {
    for(unsigned j=0;j<3;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 case 4:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<4> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<4;i++)
   {
    for(unsigned j=0;j<4;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 case 5:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<5> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<5;i++)
   {
    for(unsigned j=0;j<5;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 case 6:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<6> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<6;i++)
   {
    for(unsigned j=0;j<6;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 case 7:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]);
  OneDimensionalLegendreDShape<7> dpsi1ds(s[0]), dpsi2ds(s[1]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<7;i++)
   {
    for(unsigned j=0;j<7;j++)
     {
      //Assign the values
      dpsids(index,0) = psi2[i]*dpsi1ds[j];
      dpsids(index,1) = dpsi2ds[i]*psi1[j];
      psi[index]      = psi2[i]*psi1[j];
      //Increment the index
      ++index;
     }
   }
  break;
 }
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       polynomial order for shape functions.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<2,INITIAL_NNODE_1D>::dshape_local()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
}

//=======================================================================
/// Second derivatives of shape functions for PRefineableQElement<DIM>
/// \n d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
d2shape_local(const Vector<double> &s, Shape &psi, DShape &dpsids,
              DShape &d2psids) const
{
 std::ostringstream error_message;
 error_message <<"\nd2shape_local currently not implemented for this element\n";
 throw OomphLibError(error_message.str(),
                     "PRefineableQElement<2,INITIAL_NNODE_1D>::d2shape_local()",
                     OOMPH_EXCEPTION_LOCATION);
}

//=================================================================
/// Internal function to set up the hanging nodes on a particular
/// edge of the element.
/// Implements the mortarting method to enforce continuity weakly
/// across non-conforming element boundaries \f$\Gamma\f$ using an
/// integral matching condition
/// \f[ int_\Gamma (u_\mbox{S} - u_\mbox{M}) \psi \mbox{d} s = 0 \f]
/// for all polynomials \f$\psi\f$ on \f$\Gamma\f$ of degree at most
/// p-2 (where p is the spectral-order of the slave element) and a
/// vertex matching condition
/// \f[ (u_\mbox{S} - u_\mbox{M})\big\vert_{\partial\Gamma} = 0.\f]
/// 
/// The algorithm works as follows:
///  - First the element determines if its edge my_edge is on the
///    master or slave side of the non-conformity. At h-type non-conformities
///    short edges must be masters, and at p-type nonconformities the edge with
///    lower p-order is the master.
///  - Mortaring is performed by the slave side.
///  - The slave element constructs a list of all its neighbouring element
///    which are leaves of the tree. These elements' nodes are the masters of
///    its slave nodes.
///  - The slave element then constructs a list of all master, slave and shared
///    nodes.
///  - The integral matching condition is discretised and the mortar test
///    functions \f$ \psi \f$ are chosen to be derivatives of Legendre
///    polynomials of degree p-1.
///  - The mortar mass matrix \f$ M \f$ is constructed. Its entries are the
///    mortar test functions evaluated at the slave nodal positions, so it is
///    diagonal.
///  - Local projection matrices are constructed for each master element by
///    applying the discretised integral matching condition along its non-
///    conforming edge using the appropriate quadrature order.
///  - These are then assembled by collecting together the unknowns
///    corresponding to shared nodes and applying the vertex matching
///    condition, giving the global projection matrix \f$ P \f$.
///  - The mortar system \f$ M\xi^s = P\hat{\xi^m} \f$ is constructed,
///    where \f$ xi^m \f$ and \f$ \xi^s \f$ are the nodal values at the master
///    and slave nodes respectively.
///  - The conformity matrix \f$ C = M^{-1}P \f$ is computed. This is
///    straightforward since the mass matrix is diagonal.
///  - Finally, the master nodes and weights for each slave node are read from
///    the conformity matrix and stored in the slave's hanging scheme.
//=================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
quad_hang_helper(const int &value_id, const int &my_edge,
                 std::ofstream& output_hangfile)
{
 //RefineableQElement<2>::quad_hang_helper(value_id, my_edge, output_hangfile);
 //return;
 using namespace QuadTreeNames;
 
 Vector<unsigned> translate_s(2);
 Vector<double> s_lo_neigh(2);
 Vector<double> s_hi_neigh(2);
 int neigh_edge,diff_level;
 bool in_neighbouring_tree;
 
 // Find pointer to neighbour in this direction
 QuadTree* neigh_pt;
 neigh_pt=this->quadtree_pt()->
  gteq_edge_neighbour(my_edge, translate_s, s_lo_neigh, 
                      s_hi_neigh,neigh_edge,diff_level,in_neighbouring_tree);
 
 // Work out master/slave edges
 //----------------------------
 
 // Set up booleans
 bool h_type_master = false;
 bool h_type_slave  = false;
 bool p_type_master = false;
 bool p_type_slave  = false;
 
 // Neighbour exists and all nodes have been created
 if(neigh_pt!=0)
  {
   // Different sized element?
   if(diff_level!=0)
    {
     // Master at h-type non-conformity
     h_type_master = true;
    }
   else if(neigh_pt->nsons()==0)
    {
     unsigned my_p_order = dynamic_cast<PRefineableQElement<2,INITIAL_NNODE_1D>*>
      (this)->p_order();
     unsigned neigh_p_order = dynamic_cast<PRefineableQElement<2,INITIAL_NNODE_1D>*>
      (neigh_pt->object_pt())->p_order();
     if(neigh_p_order==my_p_order)
      {
       // At a conforming interface
      }
     else if(neigh_p_order<my_p_order)
      {
       // Slave at p-type non-conformity
       p_type_slave = true;
      }
     else
      {
       // Master at p-type non-conformity
       p_type_master = true;
      }
    }
   else
    {
     // Slave at h-type non-conformity
     h_type_slave = true;
    }
  }
 else
  {
   // Edge is on a boundary
  }
   
 // Now do the hanging nodes
 // (Amazing and efficient version)
 //--------------------------------
 if (h_type_slave || p_type_slave)
  {
   // Get pointer to neighbour objects
   Vector<PRefineableQElement<2,INITIAL_NNODE_1D>*> neigh_obj_pt;
   Vector<const QuadTree*> quadtree_neighouring_nodes;
   Vector<Vector<double> > quadtree_neighouring_s_lo,
                           quadtree_neighouring_s_hi;
   Vector<int> quadtree_neighouring_diff_level;
   neigh_pt->
    stick_neighbouring_leaves_into_vector(quadtree_neighouring_nodes,
                                          quadtree_neighouring_s_lo,
                                          quadtree_neighouring_s_hi,
                                          quadtree_neighouring_diff_level,
                                          this->quadtree_pt(),
                                          neigh_edge);
   for (unsigned e=0; e<quadtree_neighouring_nodes.size(); e++)
    {
     neigh_obj_pt.push_back(
         dynamic_cast<PRefineableQElement<2,INITIAL_NNODE_1D>*>
                (quadtree_neighouring_nodes[e]->object_pt()));
    }
   
   Vector<double> s(2);
   
   // Create vector of master, slave and shared nodes
   //------------------------------------------------
   Vector<Node*> master_node_pt, slave_node_pt, shared_node_pt,
                 slave_element_node_pt;
   std::set<Node*> master_node_pt_set;
   Vector<unsigned> slave_pos, shared_pos;
   
   // Storage for node indices in element
   unsigned first=0, space=0;
   
   // Loop over neighbouring master elements to find master nodes
   for (unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     switch(neigh_edge)
      {
      case S:
       first=0;
       space=1;
       break;
      case N:
       first=neigh_obj_pt[e]->nnode_1d()*(neigh_obj_pt[e]->nnode_1d()-1);
       space=1;
       break;
      case W:
       first=0;
       space=neigh_obj_pt[e]->nnode_1d();
       break;
      case E:
       first=neigh_obj_pt[e]->nnode_1d()-1;
       space=neigh_obj_pt[e]->nnode_1d();
       break;
      default:
       // Should never get here
       break;
      }
     for(unsigned n=0; n<neigh_obj_pt[e]->nnode_1d(); n++)
      {
       Node* curent_node_pt = neigh_obj_pt[e]->node_pt(first+n*space);
       // Add to set
       // (vertex nodes of master elements will be inserted more than once)
       if(master_node_pt_set.insert(curent_node_pt).second)
        {
         // If successfully added to the set, node is a new master node, so
         // add to end of master node vector
         // (preserves node order along the interface, rather than sorting
         // by pointer)
         master_node_pt.push_back(curent_node_pt);
        }
      }
    }
   
   // Now do slaves
   switch(my_edge)
    {
    case S:
     first=0;
     space=1;
     break;
    case N:
     first=this->nnode_1d()*(this->nnode_1d()-1);
     space=1;
     break;
    case W:
     first=0;
     space=this->nnode_1d();
     break;
    case E:
     first=this->nnode_1d()-1;
     space=this->nnode_1d();
     break;
    default:
     // Should never get here
     break;
    }
   for(unsigned n=0; n<this->nnode_1d(); n++)
    {
     // Add node to vector of slave element nodes
     slave_element_node_pt.push_back(this->node_pt(first+n*space));
     
     // Check if master node
     bool repeat=false;
     for(unsigned i=0; i<master_node_pt.size(); i++)
      {
       if(this->node_pt(first+n*space)==master_node_pt[i])
        {
         repeat=true;
         break;
        }
      }
     
     // This is a shared node
     if(repeat==true)
      {
       shared_node_pt.push_back(this->node_pt(first+n*space));
       shared_pos.push_back(n);
      }
     // Otherwise, must be a slave node
     else
      {
       slave_node_pt.push_back(this->node_pt(first+n*space));
       slave_pos.push_back(n);
      }
    }
     
#ifdef PARANOID
   // Check that there are at least two shared nodes
   // (Otherwise we cannot impose the vertex matching condition)
   if(shared_node_pt.size()<2)
    {
     oomph_info << std::endl;
     if(h_type_slave)
      {
       oomph_info << "h-type nonconformity:" << std::endl;
       oomph_info << "    diff_level = " << diff_level << std::endl;
      }
     if(p_type_slave)
      {
       oomph_info << "p-type nonconformity:" << std::endl;
       oomph_info << "    master nnode_1d = "
                  << neigh_pt->object_pt()->nnode_1d()
                  << std::endl;
       oomph_info << "     slave nnode_1d = " << this->nnode_1d() << std::endl;
      }
     oomph_info << std::endl;
         
     // Print nodal info
     for(unsigned n=0; n<master_node_pt.size(); n++)
      {
       oomph_info << "master_node_pt["<<n<<"] = " << master_node_pt[n]
                  << " at (" << master_node_pt[n]->x(0)
                  << ", " << master_node_pt[n]->x(1) << ")" << std::endl;
      }
     oomph_info << std::endl;
     for(unsigned n=0; n<slave_node_pt.size(); n++)
      {
       oomph_info << "slave_node_pt["<<n<<"] = " << slave_node_pt[n]
                  << " at (" << slave_node_pt[n]->x(0)
                  << ", " << slave_node_pt[n]->x(1)
                  << ")    slave_pos["<<n<<"] = " << slave_pos[n] << std::endl;
      }
     oomph_info << std::endl;
     for(unsigned n=0; n<shared_node_pt.size(); n++)
      {
       oomph_info << "shared_node_pt["<<n<<"] = " << shared_node_pt[n]
                  << " at (" << shared_node_pt[n]->x(0)
                  << ", " << shared_node_pt[n]->x(1) << ")" << std::endl;
      }
       
     std::string error_message =
      "There are not enough shared nodes at the interface.\n";
     error_message +=
      "The mortar method requires at least two shared nodes.\n";
       
     throw OomphLibError(
            error_message,
            "PRefineableQElement<2>::quad_hang_helper()",
            OOMPH_EXCEPTION_LOCATION);
    }
#endif
     
   // Store the number of slave and master nodes for use later
   unsigned n_slave_nodes = slave_node_pt.size();
   unsigned n_master_nodes = master_node_pt.size();
   unsigned slave_element_nnode_1d = this->nnode_1d();
   Vector<unsigned> master_element_nnode_1d(neigh_obj_pt.size());
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     master_element_nnode_1d[e] = neigh_obj_pt[e]->nnode_1d();
    }
     
   // Get master and slave nodal positions
   //-------------------------------------
   Vector<double> slave_nodal_position, slave_weight(this->nnode_1d());
   Orthpoly::gll_nodes(slave_element_nnode_1d,
                       slave_nodal_position,slave_weight);
   Vector<Vector<double> > master_nodal_position(neigh_obj_pt.size());
   Vector<Vector<double> > master_weight(neigh_obj_pt.size());
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     master_weight[e].resize(neigh_obj_pt[e]->nnode_1d());
     Orthpoly::gll_nodes(master_element_nnode_1d[e],
                         master_nodal_position[e],master_weight[e]);
    }
     
   // Assemble mass matrix for mortar
   //--------------------------------
   Vector<double> psi(n_slave_nodes);
   Shape shapes(this->nnode());
   Vector<double> diag_M(n_slave_nodes);
   Vector<Vector<double> > shared_node_M(shared_node_pt.size());
   for (unsigned i=0; i<shared_node_M.size(); i++)
    {shared_node_M[i].resize(n_slave_nodes);}
     
   for (unsigned i=0; i<n_slave_nodes; i++)
    {
     // Use L'Hosptal's rule:
     psi[i] = pow(-1.0,int(slave_element_nnode_1d-1-(slave_pos[i]-1)-1))
               *-Orthpoly::ddlegendre(slave_element_nnode_1d-1,
                             slave_nodal_position[(slave_pos[i]-1)+1]);
    }
   for (unsigned i=0; i<n_slave_nodes; i++)
    {
     diag_M[i] = psi[i]*slave_weight[(slave_pos[i]-1)+1];
    }
     
   for(unsigned v=0; v<shared_node_M.size(); v++)
    {
     for (unsigned i=0; i<n_slave_nodes; i++)
      {
       // Check if denominator is zero
       if ( fabs(slave_nodal_position[(slave_pos[i]-1)+1]
                 - slave_nodal_position[shared_pos[v]]) >= 1.0e-8 )
        {
         // We're ok
         psi[i] = pow(-1.0,int(slave_element_nnode_1d-1-(slave_pos[i]-1)-1))
                   * Orthpoly::dlegendre(slave_element_nnode_1d-1,
                                slave_nodal_position[shared_pos[v]])
                   / (slave_nodal_position[(slave_pos[i]-1)+1]
                    - slave_nodal_position[shared_pos[v]]);
        }
       // Check if numerator is zero
       else if (fabs(Orthpoly::dlegendre(slave_element_nnode_1d-1,
                                         slave_nodal_position[shared_pos[v]]))
                < 1.0e-8)
        {
         // We can use l'hopital's rule
         psi[i] = pow(-1.0,int(slave_element_nnode_1d-1-(slave_pos[i]-1)-1))
                   *-Orthpoly::ddlegendre(slave_element_nnode_1d-1,
                                 slave_nodal_position[(slave_pos[i]-1)+1]);
        }
       else
        {
         // We can't use l'hopital's rule
         throw
          OomphLibError(
           "Cannot use l'Hopital's rule. Dividing by zero is not allowed!",
           "PRefineableQElement<2>::quad_hang_helper()",
           OOMPH_EXCEPTION_LOCATION);
        }
      }
     for (unsigned i=0; i<shared_node_M[v].size(); i++)
      {
       shared_node_M[v][i] = psi[i]*slave_weight[shared_pos[v]];
      }
    }
   
   // Assemble local projection matrices for mortar
   //----------------------------------------------
   Vector<DenseDoubleMatrix*> P;
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     // Create and initialise local projection matrix P^e
     P.push_back(new DenseDoubleMatrix(n_slave_nodes,
                                       master_element_nnode_1d[e],0.0));
    }
     
     
   // Take local coordinates along bottom edge for shapes
   // (So that values are stored in the first nnode_1d() entries of shapes)
   s[1] = -1.0;
     
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     // Resize shapes
     Shape shapes(neigh_obj_pt[e]->nnode());
       
     // Sum contributions from master element shapes (quadrature)
     // (Use pointers for the quadrature knots and weights so that
     // data is not unnecessarily copied)
     unsigned quadrature_order =
      std::max(slave_element_nnode_1d,master_element_nnode_1d[e]);
     Vector<double> *quadrature_knot, *quadrature_weight;
     if (slave_element_nnode_1d > master_element_nnode_1d[e])
      {
       quadrature_knot = &slave_nodal_position;
       quadrature_weight = &slave_weight;
      }
     else
      {
       quadrature_knot = &master_nodal_position[e];
       quadrature_weight = &master_weight[e];
      }
     // Quadrature loop
     for (unsigned q=0; q<quadrature_order; q++)
      {
       // Project local master coordinates into the slave element
       unsigned active_coord_index;
       if(quadtree_neighouring_s_lo[e][0]==quadtree_neighouring_s_hi[e][0])
        {active_coord_index = 1;}
       else if(quadtree_neighouring_s_lo[e][1]==quadtree_neighouring_s_hi[e][1])
        {active_coord_index = 0;}
       else
        {
         throw OomphLibError(
                "Cannot transform coordinates",
                "PRefineableQElement<2>::quad_hang_helper()",
                OOMPH_EXCEPTION_LOCATION);
        }
       s[0] = pow(2.0,quadtree_neighouring_diff_level[e])
               * ((*quadrature_knot)[q]+1.0)
               + quadtree_neighouring_s_lo[e][active_coord_index];
       
       // Get psi
       for(unsigned k=0; k<n_slave_nodes; k++)
        {
         // Check if denominator is zero
         if (fabs(slave_nodal_position[(slave_pos[k]-1)+1]-s[0]) >= 1.0e-08)
          {
           // We're ok
           psi[k] = pow(-1.0,
                     int(slave_element_nnode_1d-1-(slave_pos[k]-1)-1))
                      * Orthpoly::dlegendre(slave_element_nnode_1d-1,s[0])
                      / (slave_nodal_position[(slave_pos[k]-1)+1]-s[0]);
          }
         // Check if numerator is zero
         else if (fabs(Orthpoly::dlegendre(slave_element_nnode_1d-1,s[0]))
                  < 1.0e-8)
          {
           // We can use l'Hopital's rule
           psi[k] = pow(-1.0,
                     int(slave_element_nnode_1d-1-(slave_pos[k]-1)-1))
                      * -Orthpoly::ddlegendre(slave_element_nnode_1d-1,s[0]);
          }
         else
          {
           // We can't use l'hopital's rule
           throw
            OomphLibError(
             "Cannot use l'Hopital's rule. Dividing by zero is not allowed!",
             "PRefineableQElement<2>::quad_hang_helper()",
             OOMPH_EXCEPTION_LOCATION);
          }
        }
       
       // Get shapes
       s[0] = (*quadrature_knot)[q];
       neigh_obj_pt[e]->shape(s,shapes);
       
       for(unsigned i=0; i<n_slave_nodes; i++)
        {
         for(unsigned j=0; j<master_element_nnode_1d[e]; j++)
          {
           P[e]->entry(i,j) += pow(2.0,
               quadtree_neighouring_diff_level[e])*shapes[j]
                 * psi[i]*(*quadrature_weight)[q];
          }
        }
      }
    }
   
   // Assemble global projection matrices for mortar
   //-----------------------------------------------
   DenseDoubleMatrix P_global(n_slave_nodes,n_master_nodes,0.0);
   
   // Assemble vertex nodes
   for(unsigned v=0; v<shared_node_pt.size(); v++)
    {
     // Find master node corresponding to shared node v
     for(unsigned k=0; k<master_node_pt.size(); k++)
      {
       if(master_node_pt[k]==shared_node_pt[v])
        {
         for(unsigned i=0; i<n_slave_nodes; i++)
          {
           P_global(i,k) -= shared_node_M[v][i];
          }
        }
      }
    }
   
   // Assemble local matrices
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     // Loop over slave element non-vertex nodes along edge
     for(unsigned i=0; i<n_slave_nodes; i++)
      {
       // Loop over local master element nodes
       for(unsigned j=0; j<master_element_nnode_1d[e]; j++)
        {
         // Get the curent local master node
         Node* cur_local_master_node;
         unsigned m_first=0, m_space=0;
         switch(neigh_edge)
          {
          case S:
           m_first=0;
           m_space=1;
           break;
          case N:
           m_first=neigh_obj_pt[e]->nnode_1d()*(neigh_obj_pt[e]->nnode_1d()-1);
           m_space=1;
           break;
          case W:
           m_first=0;
           m_space=neigh_obj_pt[e]->nnode_1d();
           break;
          case E:
           m_first=neigh_obj_pt[e]->nnode_1d()-1;
           m_space=neigh_obj_pt[e]->nnode_1d();
           break;
          default:
           // Should never get here
           break;
          }
         cur_local_master_node = neigh_obj_pt[e]->node_pt(m_first+j*m_space);
           
         // Loop over global master nodes to find match
         for(unsigned k=0; k<n_master_nodes; k++)
          {
           if (master_node_pt[k]==cur_local_master_node)
            {
             P_global(i,k) += P[e]->get_entry(i,j);
             break;
            }
          }
        }
      }
    }
     
   // Solve mortar system
   //--------------------
   for(unsigned i=0; i<n_slave_nodes; i++)
    {
     for(unsigned j=0; j<n_master_nodes; j++)
      {
       P_global(i,j)/=diag_M[i];
      }
    }
     
   // Create structures to hold the hanging info
   //-------------------------------------------
   Vector<HangInfo*> hang_info_pt(n_slave_nodes);
   for (unsigned i=0; i<n_slave_nodes; i++)
    {
     hang_info_pt[i] = new HangInfo(n_master_nodes);
    }
     
   // Set pointers to hanging info
   //-----------------------------
   for (unsigned i=0; i<n_slave_nodes; i++)
    {
     slave_node_pt[i]->set_hanging_pt(hang_info_pt[i],-1);
    }
     
   // Copy information to hanging nodes
   //----------------------------------
   for(unsigned i=0; i<n_slave_nodes; i++)
    {
     for(unsigned j=0; j<n_master_nodes; j++)
      {
       hang_info_pt[i]->set_master_node_pt(j,master_node_pt[j],P_global(i,j));
      }
    }
     
   // Free memory for local projection matrices
   for(unsigned e=0; e<neigh_obj_pt.size(); e++)
    {
     delete P[e];
    }
  }
}

//=================================================================
/// Check inter-element continuity of 
/// - nodal positions
/// - (nodally) interpolated function values
/// Overloaded to not check differences in the value. Mortaring
/// doesn't enforce strong continuity between elements.
//==================================================================== 
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<2,INITIAL_NNODE_1D>::
check_integrity(double& max_error)
{

 // BENFLAG: Overloaded to *not* check for continuity in value of interpolated
 // variables. This is necessary because mortaring does not ensure continuity
 // across element boundaries. It therefore makes no sense to test for this.

 using namespace QuadTreeNames;

 // Number of nodes along edge
 unsigned n_p=nnode_1d();

 // Number of timesteps (incl. present) for which continuity is 
 // to be checked.
 unsigned n_time=1;

 // Initialise errors
 max_error=0.0;
 Vector<double> max_error_x(2,0.0);
 double max_error_val=0.0;

 Vector<int> edges(4);
 edges[0] = S; edges[1] = N; edges[2] = W; edges[3] = E;

 //Loop over the edges
 for(unsigned edge_counter=0;edge_counter<4;edge_counter++)
  {
   Vector<unsigned> translate_s(2);
   Vector<double> s(2), s_lo_neigh(2), s_hi_neigh(2), s_fraction(2);
   int neigh_edge,diff_level;
   bool in_neighbouring_tree;
   
   // Find pointer to neighbour in this direction
   QuadTree* neigh_pt;
   neigh_pt=quadtree_pt()->gteq_edge_neighbour(edges[edge_counter], 
                                               translate_s,
                                               s_lo_neigh,s_hi_neigh,
                                               neigh_edge,diff_level,
                                               in_neighbouring_tree);
   
   // Neighbour exists and has existing nodes
   if((neigh_pt!=0) && (neigh_pt->object_pt()->nodes_built()))
    {
     //Need to exclude periodic nodes from this check
     //There are only periodic nodes if we are in a neighbouring tree
     bool is_periodic=false;
     if(in_neighbouring_tree)
      {
       //Is it periodic
       is_periodic = 
        this->tree_pt()->root_pt()->is_neighbour_periodic(edges[edge_counter]);
      }

     // Loop over nodes along the edge
     for(unsigned i0=0;i0<n_p;i0++)
      {
       //Storage for pointer to the local node
       Node* local_node_pt=0;
       
       switch(edge_counter)
        {
        case 0:
         // Local fraction of node
         s_fraction[0] = local_one_d_fraction_of_node(i0,0);
         s_fraction[1] = 0.0;
         // Get pointer to local node
         local_node_pt = this->node_pt(i0);
         break;
         
        case 1:
         // Local fraction of node
         s_fraction[0] = local_one_d_fraction_of_node(i0,0);
         s_fraction[1] = 1.0;
         // Get pointer to local node
         local_node_pt =  this->node_pt(i0 + n_p*(n_p-1));
         break;
         
         case 2:
          // Local fraction of node
          s_fraction[0] = 0.0; 
          s_fraction[1] = local_one_d_fraction_of_node(i0,1);
          // Get pointer to local node
          local_node_pt = this->node_pt(n_p*i0);
          break;
          
        case 3:
         // Local fraction of node
         s_fraction[0] = 1.0; 
         s_fraction[1] = local_one_d_fraction_of_node(i0,1);          
         // Get pointer to local node
         local_node_pt = this->node_pt(n_p-1 + n_p*i0);
         break;
        }

       //Calculate the local coordinate and the local coordinate as viewed
       //from the neighbour
       Vector<double> s_in_neighb(2);
       for(unsigned i=0;i<2;i++)
        {
         //Local coordinate in this element
         s[i] = -1.0 + 2.0*s_fraction[i];
         //Local coordinate in the neighbour
         s_in_neighb[i] = s_lo_neigh[i] + s_fraction[translate_s[i]]*
          (s_hi_neigh[i] - s_lo_neigh[i]);
        }
       
       //Loop over timesteps
       for(unsigned t=0;t<n_time;t++)
        {
         // Get the nodal position from neighbour element
         Vector<double> x_in_neighb(2);
         neigh_pt->object_pt()->interpolated_x(t,s_in_neighb,x_in_neighb);
         
         // Check error only if the node is NOT periodic
         if(is_periodic==false)
          {
           for(int i=0;i<2;i++)
            {
             //Find the spatial error
             double err = std::abs(local_node_pt->x(t,i) - x_in_neighb[i]);
             
             //If it's bigger than our tolerance, say so
             if (err>1e-9)
              {
               oomph_info << "errx " << err << " " << t << " " 
                          << local_node_pt->x(t,i) 
                          << " " <<  x_in_neighb[i]<< std::endl;
               
               oomph_info << "at " <<  local_node_pt->x(0) << " "
                          <<  local_node_pt->x(1) << std::endl;
              }
             
             //If it's bigger than the previous max error, it is the
             //new max error!
             if (err>max_error_x[i]) {max_error_x[i]=err;}
            }
          }

         //BENFLAG: This test is not required
         /*
         // Get the values from neighbour element. Note: # of values
         // gets set by routine (because in general we don't know
         // how many interpolated values a certain element has
         Vector<double> values_in_neighb;
         neigh_pt->object_pt()->
          get_interpolated_values(t,s_in_neighb,values_in_neighb);
         
         // Get the values in current element.
         Vector<double> values;
         get_interpolated_values(t,s,values);
         
         // Now figure out how many continuously interpolated values there are
         unsigned num_val=neigh_pt->object_pt()->ncont_interpolated_values();
         
         // Check error
         for(unsigned ival=0;ival<num_val;ival++)
          {
           double err=std::abs(values[ival] - values_in_neighb[ival]);
           
           if (err>1.0e-10)
             {
              oomph_info <<  local_node_pt->x(0) << " " 
                        <<  local_node_pt->x(1) << " \n# "
                        << "erru (S)" << err << " " << ival << " " 
                        << get_node_number(local_node_pt) << " "
                        << values[ival]
                        << " " << values_in_neighb[ival] << std::endl;
             }
           
           if (err>max_error_val) {max_error_val=err;}
           
          }
         */
        }
       
      }
    }
  }
 
 max_error=max_error_x[0];
 if (max_error_x[1]>max_error) max_error=max_error_x[1];
 if (max_error_val>max_error) max_error=max_error_val;
 
 if (max_error>1e-9)
  {
   oomph_info << "\n#------------------------------------ \n#Max error " ;
   oomph_info << max_error_x[0] 
        << " " << max_error_x[1] 
        << " " << max_error_val << std::endl;
   oomph_info << "#------------------------------------ \n " << std::endl;
   
  }

}

////////////////////////////////////////////////////////////////
//       3D PRefineableQElements
////////////////////////////////////////////////////////////////

/// Get local coordinates of node j in the element; vector sets its own size
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
local_coordinate_of_node(const unsigned& n, Vector<double>& s)
 {
  OomphLibWarning(
   "3D PRefineableQElements have not been fully implemented.",
   "PRefineableQElement<3,INITIAL_NNODE_1D>::local_coordinate_of_node()",
   OOMPH_EXCEPTION_LOCATION);
  
  s.resize(3);
  unsigned Nnode_1d = this->nnode_1d();
  unsigned n0 = n%Nnode_1d;
  unsigned n1 = unsigned(double(n)/double(Nnode_1d))%Nnode_1d;
  unsigned n2 = unsigned(double(n)/double(Nnode_1d*Nnode_1d));
  
  switch(Nnode_1d)
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<2>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<2>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<2>::nodal_position(n2);
    break;
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<3>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<3>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<3>::nodal_position(n2);
    break;
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<4>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<4>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<4>::nodal_position(n2);
    break;
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<5>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<5>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<5>::nodal_position(n2);
    break;
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<6>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<6>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<6>::nodal_position(n2);
    break;
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    s[0] = OneDimensionalLegendreShape<7>::nodal_position(n0);
    s[1] = OneDimensionalLegendreShape<7>::nodal_position(n1);
    s[2] = OneDimensionalLegendreShape<7>::nodal_position(n2);
    break;
  default:
    std::ostringstream error_message;
    error_message <<"\nERROR: Exceeded maximum polynomial order for";
    error_message <<"\n       shape functions.\n";
    throw OomphLibError(error_message.str(),
                        "PRefineableQElement<3,INITIAL_NNODE_1D>::local_coordinate_of_node()",
                        OOMPH_EXCEPTION_LOCATION);
    break;
   }
 }
 
/// Get the local fractino of node j in the element
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
local_fraction_of_node(const unsigned &n, Vector<double> &s_fraction)
 {
  OomphLibWarning(
   "3D PRefineableQElements have not been fully implemented.",
   "PRefineableQElement<3,INITIAL_NNODE_1D>::local_fraction_of_node()",
   OOMPH_EXCEPTION_LOCATION);
  
  this->local_coordinate_of_node(n,s_fraction);
  s_fraction[0] = 0.5*(s_fraction[0] + 1.0);
  s_fraction[1] = 0.5*(s_fraction[1] + 1.0);
  s_fraction[2] = 0.5*(s_fraction[2] + 1.0);
 }

template<unsigned INITIAL_NNODE_1D>
double PRefineableQElement<3,INITIAL_NNODE_1D>::
local_one_d_fraction_of_node(const unsigned &n1d, const unsigned &i)
 {
  OomphLibWarning(
   "3D PRefineableQElements have not been fully implemented.",
   "PRefineableQElement<3,INITIAL_NNODE_1D>::one_d_fraction_of_node()",
   OOMPH_EXCEPTION_LOCATION);
  
  switch(this->nnode_1d())
   {
  case 2:
    OneDimensionalLegendreShape<2>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<2>::nodal_position(n1d) + 1.0);
  case 3:
    OneDimensionalLegendreShape<3>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<3>::nodal_position(n1d) + 1.0);
  case 4:
    OneDimensionalLegendreShape<4>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<4>::nodal_position(n1d) + 1.0);
  case 5:
    OneDimensionalLegendreShape<5>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<5>::nodal_position(n1d) + 1.0);
  case 6:
    OneDimensionalLegendreShape<6>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<6>::nodal_position(n1d) + 1.0);
  case 7:
    OneDimensionalLegendreShape<7>::calculate_nodal_positions();
    return 0.5*(OneDimensionalLegendreShape<7>::nodal_position(n1d) + 1.0);
  default:
    std::ostringstream error_message;
    error_message <<"\nERROR: Exceeded maximum polynomial order for";
    error_message <<"\n       shape functions.\n";
    throw OomphLibError(error_message.str(),
                  "PRefineableQElement<3,INITIAL_NNODE_1D>::local_one_d_fraction_of_node()",
                  OOMPH_EXCEPTION_LOCATION);
    return 0.0;
   }
 }
  
//==================================================================
/// Return the node at the specified local coordinate
//==================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<3,INITIAL_NNODE_1D>::
get_node_at_local_coordinate(const Vector<double> &s)
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::get_node_at_local_coordinate()",
  OOMPH_EXCEPTION_LOCATION);
  
 unsigned Nnode_1d = this->nnode_1d();
 //Load the tolerance into a local variable
 double tol = FiniteElement::Node_location_tolerance;
 //There are two possible indices.
 Vector<int> index(3);
 
 // Loop over indices
 for (unsigned i=0; i<3; i++)
  {
   // Determine the index
   // -------------------
   
   bool is_found=false;
   
   // If we are at the lower limit, the index is zero
   if(std::abs(s[i] + 1.0) < tol)
    {
     index[i] = 0;
     is_found=true;
    }
   // If we are at the upper limit, the index is the number of nodes minus 1
   else if(std::abs(s[i] - 1.0) < tol)
    {
     index[i] = Nnode_1d-1;
     is_found=true;
    }
   // Otherwise, we have to calculate the index in general
   else
    {
     // Compute Gauss-Lobatto-Legendre node positions
     Vector<double> z;
     Orthpoly::gll_nodes(Nnode_1d, z);
     // Loop over possible internal nodal positions
     for (unsigned n=1; n<Nnode_1d-1; n++)
      {
       if (std::abs(z[n] - s[i]) < tol)
        {
         index[i] = n;
         is_found=true;
         break;
        }
      }
    }
   
   if (!is_found)
    {
     // No matching nodes
     return 0;
    }
  }
 // If we've got here we have a node, so let's return a pointer to it
 return this->node_pt(index[0] + Nnode_1d*index[1] + Nnode_1d*Nnode_1d*index[2]);
}

//===================================================================
/// If a neighbouring element has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0).
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<3,INITIAL_NNODE_1D>::
node_created_by_neighbour(const Vector<double> &s_fraction, bool &is_periodic) 
{  
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::node_created_by_neighbour()",
  OOMPH_EXCEPTION_LOCATION);
  
 using namespace OcTreeNames;
 
 //Calculate the faces/edges on which the node lies
 Vector<int> faces;
 Vector<int> edges;

 if(s_fraction[0]==0.0)
  {
   faces.push_back(L);
   if (s_fraction[1]==0.0) {edges.push_back(LD);}
   if (s_fraction[2]==0.0) {edges.push_back(LB);}
   if (s_fraction[1]==1.0) {edges.push_back(LU);}
   if (s_fraction[2]==1.0) {edges.push_back(LF);}
  }

 if(s_fraction[0]==1.0) 
  {
   faces.push_back(R);
   if (s_fraction[1]==0.0) {edges.push_back(RD);}
   if (s_fraction[2]==0.0) {edges.push_back(RB);}
   if (s_fraction[1]==1.0) {edges.push_back(RU);}
   if (s_fraction[2]==1.0) {edges.push_back(RF);}
  }

 if(s_fraction[1]==0.0)
  {
   faces.push_back(D);
   if (s_fraction[2]==0.0) {edges.push_back(DB);}
   if (s_fraction[2]==1.0) {edges.push_back(DF);}
  }

 if(s_fraction[1]==1.0) 
  {
   faces.push_back(U);
   if (s_fraction[2]==0.0) {edges.push_back(UB);}
   if (s_fraction[2]==1.0) {edges.push_back(UF);}
  }

 if(s_fraction[2]==0.0)
  {
   faces.push_back(B);
  }

 if(s_fraction[2]==1.0)
  {
   faces.push_back(F);
  }
 
 //Find the number of faces
 unsigned n_face = faces.size();
 
 //Find the number of edges
 unsigned n_edge = edges.size();
 
 Vector<unsigned> translate_s(3);
 Vector<double> s_lo_neigh(3);
 Vector<double> s_hi_neigh(3);
 Vector<double> s(3);

 int neigh_face, diff_level;
 
 //Loop over the faces on which the node lies
 //------------------------------------------
 for(unsigned j=0;j<n_face;j++)
  {
   // Find pointer to neighbouring element along face 
   OcTree* neigh_pt;
   neigh_pt = octree_pt()->
    gteq_face_neighbour(faces[j],translate_s,s_lo_neigh,s_hi_neigh,neigh_face,
                        diff_level);
   
   // Neighbour exists
   if(neigh_pt!=0)
    {
     // Have any of its vertex nodes been created yet?
     // (BS: Must look in incomplete neighbours because after the pre-build
     // they may contain pointers to the required nodes. e.g. h-refinement of
     // neighbouring linear and quadratic elements)
     bool a_vertex_node_is_built = false;
     QElement<3,INITIAL_NNODE_1D>* neigh_obj_pt =
      dynamic_cast<QElement<3,INITIAL_NNODE_1D>*>(neigh_pt->object_pt());
     if(neigh_obj_pt==0)
      {
       throw
        OomphLibError("Not a quad element!",
         "PRefineableQElement<3,INITIAL_NNODE_1D>::node_created_by_neighbour()",
         OOMPH_EXCEPTION_LOCATION);
      }
     for(unsigned vnode=0; vnode<neigh_obj_pt->nvertex_node(); vnode++)
      {
       if(neigh_obj_pt->vertex_node_pt(vnode)!=0)
        a_vertex_node_is_built = true;
       break;
      }
     if(a_vertex_node_is_built)
      {
       //We now need to translate the nodal location, defined in terms
       //of the fractional coordinates of the present element into
       //those of its neighbour. For this we use the information returned
       //to use from the octree function.
       
       //Calculate the local coordinate in the neighbour
       //Note that we need to use the translation scheme in case
       //the local coordinates are swapped in the neighbour.
       for(unsigned i=0;i<3;i++)
        {
         s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]]*
          (s_hi_neigh[i] - s_lo_neigh[i]);
        }
       
       //Find the node in the neighbour
       Node* neighbour_node_pt =
        neigh_pt->object_pt()->get_node_at_local_coordinate(s);
       
       //If there is a node, return it
       if(neighbour_node_pt!=0)
        {
         return neighbour_node_pt;
        }
      }
    }
  } //End of loop over faces


 //Loop over the edges on which the node lies
 //------------------------------------------
 for(unsigned j=0;j<n_edge;j++)
  {
   // Find pointer to neighbouring element along edge
   OcTree* neigh_pt;
   // Warning, need to search additional trees than these:
   unsigned i_root_edge_neighbour=0, nroot_edge_neighbour=0;
   neigh_pt = octree_pt()->
    gteq_true_edge_neighbour(edges[j],
                             i_root_edge_neighbour, nroot_edge_neighbour,
                             translate_s,s_lo_neigh,s_hi_neigh,neigh_face,
                             diff_level);
   
   // Neighbour exists
   if(neigh_pt!=0)
    {
     // Have its nodes been created yet?
     if(true || neigh_pt->object_pt()->nodes_built())
      {
       //We now need to translate the nodal location, defined in terms
       //of the fractional coordinates of the present element into
       //those of its neighbour. For this we use the information returned
       //to use from the octree function.
       
       //Calculate the local coordinate in the neighbour
       //Note that we need to use the translation scheme in case
       //the local coordinates are swapped in the neighbour.
       for(unsigned i=0;i<3;i++)
        {
         s[i] = s_lo_neigh[i] + s_fraction[translate_s[i]]*
          (s_hi_neigh[i] - s_lo_neigh[i]);
        }
       
       //Find the node in the neighbour
       Node* neighbour_node_pt =
        neigh_pt->object_pt()->get_node_at_local_coordinate(s);
       
       //If there is a node, return it
       if(neighbour_node_pt!=0)
        {
         return neighbour_node_pt;
        }
      }
    }
  } //End of loop over faces

 //Node not found, return null
 return 0;
}

//===================================================================
/// If a neighbouring element's son has already created a node at
/// a position corresponding to the local fractional position within the
/// present element, s_fraction, return
/// a pointer to that node. If not, return NULL (0). If the node is
/// periodic the flag is_periodic will be true
//===================================================================
template<unsigned INITIAL_NNODE_1D>
Node* PRefineableQElement<3,INITIAL_NNODE_1D>::
node_created_by_son_of_neighbour(const Vector<double> &s_fraction, 
                                 bool &is_periodic) 
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::node_created_by_son_of_neighbour()",
  OOMPH_EXCEPTION_LOCATION);

 throw
  OomphLibError(
   "This function has not yet been implemented.",
   "PRefineableQElement<3,INITIAL_NNODE_1D>::node_created_by_son_of_neighbour()",
   OOMPH_EXCEPTION_LOCATION);
  
 return 0;
}

//==================================================================
/// Set the correct p-order of the element based on that of its
/// father. Then construct an integration scheme of the correct order.
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::initial_setup()
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::initial_setup()",
  OOMPH_EXCEPTION_LOCATION);
  
 // Check if element is in a tree
 if (Tree_pt!=0)
  {
   //Pointer to my father (in quadtree impersonation)
   OcTree* father_pt = dynamic_cast<OcTree*>(octree_pt()->father_pt());
   
   // Check if element has father
   if (father_pt!=0)
    {
     PRefineableQElement<3,INITIAL_NNODE_1D>* father_el_pt =
            dynamic_cast<PRefineableQElement<3,INITIAL_NNODE_1D>*>
                (this->tree_pt()->father_pt()->object_pt());
     if (father_el_pt!=0)
      {
       unsigned father_p_order = father_el_pt->p_order();
       // Set p-order to that of father
       P_order = father_p_order;
      }
     
     // Now sort out the element...
     // (has p^3 nodes)
     unsigned new_n_node = P_order*P_order*P_order;
     
     // Allocate new space for Nodes (at the element level)
     this->set_n_node(new_n_node);
     
     // Set integration scheme
     delete this->integral_pt();
     switch(P_order)
     {
     case 2:
      this->set_integration_scheme(new GaussLobattoLegendre<3,2>);
      break;
     case 3:
      this->set_integration_scheme(new GaussLobattoLegendre<3,3>);
      break;
     case 4:
      this->set_integration_scheme(new GaussLobattoLegendre<3,4>);
      break;
     case 5:
      this->set_integration_scheme(new GaussLobattoLegendre<3,5>);
      break;
     case 6:
      this->set_integration_scheme(new GaussLobattoLegendre<3,6>);
      break;
     case 7:
      this->set_integration_scheme(new GaussLobattoLegendre<3,7>);
      break;
     default:
      std::ostringstream error_message;
      error_message <<"\nERROR: Exceeded maximum polynomial order for";
      error_message <<"\n       integration scheme.\n";
      throw OomphLibError(error_message.str(),
                          "PRefineableQElement<3,INITIAL_NNODE_1D>::initial_setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    }
  }
 else
  {
   throw OomphLibError("Element not in a tree!",
                       "PRefineableQElement<3,INITIAL_NNODE_1D>::initial_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//==================================================================
/// Set the correct p-order of the element before it is built
/// and then check the father element for any required nodes which
/// already exist
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::pre_build()
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::pre_build()",
  OOMPH_EXCEPTION_LOCATION);
  
 //Pointer to my father (in quadtree impersonation)
 OcTree* father_pt = dynamic_cast<OcTree*>(octree_pt()->father_pt());
 
 // Check if element has father
 if (father_pt!=0)
  {
   PRefineableQElement<3,INITIAL_NNODE_1D>* father_el_pt =
          dynamic_cast<PRefineableQElement<3,INITIAL_NNODE_1D>*>
              (this->tree_pt()->father_pt()->object_pt());
   if (father_el_pt!=0)
    {
     // Check nodes from father
     //------------------------------------
     unsigned n_p = this->nnode_1d();
    
     // What type of son am I? Ask my octree representation...
     int son_type = Tree_pt->son_type();
     
     Vector<double> s_lo(3);
     Vector<double> s_hi(3);
     Vector<double> s(3);
     Vector<double> x(3);
     
     using namespace OcTreeNames;
     
     // Setup vertex coordinates in father element:
     switch(son_type)
      {
      case LDB:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       s_lo[2]=-1.0;
       s_hi[2]= 0.0;
       break;
      
      case LDF:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       s_lo[2]= 0.0;
       s_hi[2]= 1.0;
       break;
      
      case LUB:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       s_lo[2]=-1.0;
       s_hi[2]= 0.0;
       break;
      
      case LUF:
       s_lo[0]=-1.0;
       s_hi[0]= 0.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       s_lo[2]= 0.0;
       s_hi[2]= 1.0;
       break;
      
      case RDB:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       s_lo[2]=-1.0;
       s_hi[2]= 0.0;
       break;
      
      case RDF:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]=-1.0;
       s_hi[1]= 0.0;
       s_lo[2]= 0.0;
       s_hi[2]= 1.0;
       break;
      
      case RUB:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       s_lo[2]=-1.0;
       s_hi[2]= 0.0;
       break;
      
      case RUF:
       s_lo[0]= 0.0;
       s_hi[0]= 1.0;
       s_lo[1]= 0.0;
       s_hi[1]= 1.0;
       s_lo[2]= 0.0;
       s_hi[2]= 1.0;
       break;
      }
     unsigned jnod=0;
   
     Vector<double> s_fraction(3);
     // Loop over nodes in element
     for(unsigned i0=0;i0<n_p;i0++)
      {
       //Get the fractional position of the node in the direction of s[0]
       s_fraction[0] = local_one_d_fraction_of_node(i0,0);
       // Local coordinate in father element
       s[0] = s_lo[0] + (s_hi[0]-s_lo[0])*s_fraction[0];
   
       for(unsigned i1=0;i1<n_p;i1++)
        {
         //Get the fractional position of the node in the direction of s[1]
         s_fraction[1] = local_one_d_fraction_of_node(i1,1);
         // Local coordinate in father element
         s[1] = s_lo[1] + (s_hi[1]-s_lo[1])*s_fraction[1];
   
         for(unsigned i2=0;i2<n_p;i2++)
          {
           //Get the fractional position of the node in the direction of s[1]
           s_fraction[2] = local_one_d_fraction_of_node(i2,2);
           // Local coordinate in father element
           s[2] = s_lo[2] + (s_hi[2]-s_lo[2])*s_fraction[2];
         
           // Local node number
           jnod= i0 + n_p*i1 + n_p*n_p*i2;
           
           //Get the pointer to the node in the father, returns NULL
           //if there is not node
           Node* created_node_pt = father_el_pt->get_node_at_local_coordinate(s);
           
           // Does this node already exist in father element?
           //------------------------------------------------
           if(created_node_pt!=0) 
            {
             // Copy node across
             node_pt(jnod) = created_node_pt;
            }
          }
        }
      }
    }
   else
    {
     std::ostringstream error_message;
     error_message <<"\nERROR: Dynamic cast failed!\n";
     throw OomphLibError(error_message.str(),
                         "PRefineableQElement<3,INITIAL_NNODE_1D>::pre_build()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
}

//==================================================================
/// p-refine the element inc times. (If inc<0 then p-unrefinement
/// is performed.)
//==================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine(const int &inc,
                                                       Mesh* &mesh_pt)
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
  OOMPH_EXCEPTION_LOCATION);

 // Timestepper should be the same for all nodes -- use it
 // to create timesteppers for new nodes
 TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
 
 // Number of history values (incl. present)
 unsigned ntstorage = time_stepper_pt->ntstorage();
 
 // Back up pointers to old vertex nodes before they are lost
 unsigned n_vertex_node = this->nvertex_node();
 Vector<Node*> old_vertex_node_pt(n_vertex_node);
 for (unsigned n=0; n<n_vertex_node; n++)
  {
   old_vertex_node_pt[n] = this->vertex_node_pt(n);
  }
 
 // Compute new coordinates and projected values
 Vector<double> new_local_x((P_order + inc)*(P_order + inc)*(P_order + inc));
 Vector<Vector<double> > new_global_x((P_order + inc)*(P_order + inc)*(P_order + inc));
 for (unsigned i=0; i<(P_order + inc)*(P_order + inc)*(P_order + inc); i++)
  {new_global_x[i].resize(3);}
 unsigned ncont = this->ncont_interpolated_values();
 Vector<Vector<Vector<double> > > projected_value(ntstorage);
 for(unsigned t=0;t<ntstorage;t++)
  {
   projected_value[t].resize((P_order + inc)*(P_order + inc)*(P_order + inc));
   for (unsigned i=0; i<(P_order + inc)*(P_order + inc)*(P_order + inc); i++)
    {
     projected_value[t][i].resize(ncont);
    }
  }
 
 // Compute Gauss-Lobatto-Legendre node spacing
 Orthpoly::gll_nodes(P_order + inc, new_local_x);
 
 for (unsigned n=0; n<(P_order + inc)*(P_order + inc)*(P_order + inc); n++)
  {
   // Create coordinate vector
   Vector<double> s(3);
   unsigned n0 = n%(P_order+inc);
   unsigned n1 = unsigned(double(n)/double(P_order+inc))%(P_order+inc);
   unsigned n2 = unsigned(double(n)/double((P_order+inc)*(P_order+inc)));
   s[0] = new_local_x[n0];
   s[1] = new_local_x[n1];
   s[2] = new_local_x[n2];
   
   new_global_x[n][0] = this->interpolated_x(s,0);
   new_global_x[n][1] = this->interpolated_x(s,1);
   new_global_x[n][2] = this->interpolated_x(s,2);
   
   // Loop over all history values
   for(unsigned t=0;t<ntstorage;t++)
    {
     // Interpolate new nodal values
     // (while still using old integration scheme and shape functions)
     Vector<double> values(ncont);
     this->get_interpolated_values(t,s,values);
     for(unsigned i=0; i<ncont; i++)
      {
       projected_value[t][n][i] = values[i];
      }
    }
  }
 
 // Increment p-order of the element
 P_order += inc;
 
 // Change integration scheme
 delete this->integral_pt();
 switch(P_order)
 {
 case 2:
  this->set_integration_scheme(new GaussLobattoLegendre<3,2>);
  break;
 case 3:
  this->set_integration_scheme(new GaussLobattoLegendre<3,3>);
  break;
 case 4:
  this->set_integration_scheme(new GaussLobattoLegendre<3,4>);
  break;
 case 5:
  this->set_integration_scheme(new GaussLobattoLegendre<3,5>);
  break;
 case 6:
  this->set_integration_scheme(new GaussLobattoLegendre<3,6>);
  break;
 case 7:
  this->set_integration_scheme(new GaussLobattoLegendre<3,7>);
  break;
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       integration scheme.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<3,INITIAL_NNODE_1D>::p_refine()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
 // Allocate new space for Nodes (at the element level)
 this->set_n_node(P_order*P_order*P_order);
 unsigned new_Nnode = this->nnode();
 
 // Copy vertex nodes and create new edge and internal nodes
 //---------------------------------------------------------
 
 // Copy vertex nodes
 //for(unsigned n=0; n<n_vertex_node; n++)
 // {
 //  this->vertex_node_pt(n) = old_vertex_node_pt[n];
 // }
 this->node_pt(0                              ) = old_vertex_node_pt[0];
 this->node_pt(P_order-1                      ) = old_vertex_node_pt[1];
 this->node_pt(P_order*(P_order-1)            ) = old_vertex_node_pt[2];
 this->node_pt(P_order*P_order-1              ) = old_vertex_node_pt[3];
 this->node_pt(P_order*P_order*(P_order-1)    ) = old_vertex_node_pt[4];
 this->node_pt((P_order*P_order+1)*(P_order-1)) = old_vertex_node_pt[5];
 this->node_pt(P_order*(P_order+1)*(P_order-1)) = old_vertex_node_pt[6];
 this->node_pt(P_order*P_order*P_order-1      ) = old_vertex_node_pt[7];


 //=======================================
 //Still to fill in:
 //Find/create all nodes in the element...
 //=======================================

 
 // Set coordinates and project data
 for(unsigned t=0;t<ntstorage;t++)
  {
   for (unsigned n=0; n<new_Nnode; n++)
    {
     for(unsigned i=0; i<3; i++)
      {
       this->node_pt(n)->x(i) = new_global_x[n][i];
      }
     for(unsigned i=0; i<ncont; i++)
      {
       this->node_pt(n)->set_value(t,i,projected_value[t][n][i]);
      }
    }
  }
 
 // Not necessary to delete the old nodes since all original nodes are in the
 // current mesh and so will be pruned as part of the mesh adaption process.
 
 // Do any further-build required
 further_build();
 
}

//=======================================================================
///Shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
shape(const Vector<double> &s, Shape &psi) const
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::shape()",
  OOMPH_EXCEPTION_LOCATION);
  
 switch(P_order)
 {
 case 2:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 case 3:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 case 4:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 case 5:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 case 6:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 case 7:
 {
  //Call the OneDimensional Shape functions
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]), psi3(s[2]);
 
  //Now let's loop over the nodal points in the element
  //and copy the values back in  
  for(unsigned i=0;i<P_order;i++) 
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        psi(P_order*P_order*i + P_order*j + k) = psi3[i]*psi2[j]*psi1[k];
       }
     }
   }
  break;
 }
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       polynomial order for shape functions.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<3,INITIAL_NNODE_1D>::shape()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
}

//=======================================================================
///Derivatives of shape functions for PRefineableQElement<DIM>
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsids) const
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::dshape_local()",
  OOMPH_EXCEPTION_LOCATION);
  
 switch(P_order)
 {
 case 2:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<2>::calculate_nodal_positions();
  OneDimensionalLegendreShape<2> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<2> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 case 3:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<3>::calculate_nodal_positions();
  OneDimensionalLegendreShape<3> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<3> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 case 4:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<4>::calculate_nodal_positions();
  OneDimensionalLegendreShape<4> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<4> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 case 5:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<5>::calculate_nodal_positions();
  OneDimensionalLegendreShape<5> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<5> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 case 6:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<6>::calculate_nodal_positions();
  OneDimensionalLegendreShape<6> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<6> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 case 7:
 {
  //Call the shape functions and derivatives
  OneDimensionalLegendreShape<7>::calculate_nodal_positions();
  OneDimensionalLegendreShape<7> psi1(s[0]), psi2(s[1]), psi3(s[2]);
  OneDimensionalLegendreDShape<7> dpsi1ds(s[0]), dpsi2ds(s[1]), dpsi3ds(s[2]);
 
  //Index for the shape functions
  unsigned index=0;
  //Loop over shape functions in element
  for(unsigned i=0;i<P_order;i++)
   {
    for(unsigned j=0;j<P_order;j++)
     {
      for(unsigned k=0;k<P_order;k++)
       {
        //Assign the values
        dpsids(index,0) = psi3[i]*psi2[j]*dpsi1ds[k];
        dpsids(index,1) = psi3[i]*dpsi2ds[j]*psi1[k];
        dpsids(index,2) = dpsi3ds[i]*psi2[j]*psi1[k];
        psi[index]      = psi3[i]*psi2[j]*psi1[k];
        //Increment the index
        ++index;
       }
     }
   }
  break;
 }
 default:
  std::ostringstream error_message;
  error_message <<"\nERROR: Exceeded maximum polynomial order for";
  error_message <<"\n       polynomial order for shape functions.\n";
  throw OomphLibError(error_message.str(),
                      "PRefineableQElement<3,INITIAL_NNODE_1D>::dshape_local()",
                      OOMPH_EXCEPTION_LOCATION);
 }
 
}

//=======================================================================
/// Second derivatives of shape functions for PRefineableQElement<DIM>
/// \n d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
//=======================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
d2shape_local(const Vector<double> &s, Shape &psi, DShape &dpsids,
              DShape &d2psids) const
{
 std::ostringstream error_message;
 error_message <<"\nd2shape_local currently not implemented for this element\n";
 throw OomphLibError(error_message.str(),
                     "PRefineableQElement<3,INITIAL_NNODE_1D>::d2shape_local()",
                     OOMPH_EXCEPTION_LOCATION);
}

//=================================================================
/// Internal function to set up the hanging nodes on a particular
/// edge of the element. Implements the mortar method.
//=================================================================
template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
oc_hang_helper(const int &value_id, const int &my_edge,
               std::ofstream& output_hangfile)
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::oc_hang_helper()",
  OOMPH_EXCEPTION_LOCATION);
  
 //Not implemented yet
 RefineableQElement<3>::oc_hang_helper(value_id, my_edge, output_hangfile);
}

template<unsigned INITIAL_NNODE_1D>
void PRefineableQElement<3,INITIAL_NNODE_1D>::
check_integrity(double& max_error)
{
 OomphLibWarning(
  "3D PRefineableQElements have not been fully implemented.",
  "PRefineableQElement<3,INITIAL_NNODE_1D>::node_created_by_son_of_neighbour()",
  OOMPH_EXCEPTION_LOCATION);
 
 throw
  OomphLibError("This function has not yet been implemented.",
                "PRefineableQElement<3,INITIAL_NNODE_1D>::check_integrity()",
                OOMPH_EXCEPTION_LOCATION);
 
 //Not yet implemented
}

//===================================================================
// Build required templates
//===================================================================
template class PRefineableQElement<1,2>;
template class PRefineableQElement<1,3>;
template class PRefineableQElement<1,4>;

template class PRefineableQElement<2,2>;
template class PRefineableQElement<2,3>;
template class PRefineableQElement<2,4>;

template class PRefineableQElement<3,2>;
template class PRefineableQElement<3,3>;
template class PRefineableQElement<3,4>;

}

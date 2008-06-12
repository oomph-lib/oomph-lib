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
//Header file for base class that specified common interfaces
//used in elements with moving Nodes (Spine,Algebraic,MacroElementNodeUpdate)

//Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_ELEMENT_WITH_MOVING_NODES
#define OOMPH_ELEMENT_WITH_MOVING_NODES


//oomph-lib headers
#include "elements.h"

namespace oomph
{
 
 //=======================================================================
 /// \short A policy class that serves to establish the 
 /// common interfaces for elements that contain moving nodes. This
 /// class provides storage for the geometric data that affect the 
 /// update of all the nodes of the element, i.e. USUALLY 
 /// all data that are using during a call
 /// to the Element's node_update() function. In some cases 
 /// (e.g. FluidInterfaceEdge elements), node_update() is overloaded to
 /// perform an update of the bulk element, in which case the additional
 /// bulk geometric data become external data of the element and the
 /// function GeneralisedElement::update_in_external_fd(i) is overloaded
 /// to also perform the bulk node update.
 /// The storage is populated 
 /// during the assignment of the equation numbers via the 
 /// complete_setup_of_dependencies() function and then local equations
 /// numbers are assigned to these data, accessible via 
 /// geometric_data_local_eqn(n,i). Finally, a function is provided 
 /// that calculates the terms in the jacobian matrix by due to these
 /// geometric data by finite differences.
 //=======================================================================
 class ElementWithMovingNodes : public virtual FiniteElement
{
  public:

 ///Empty constructor 
 ElementWithMovingNodes() : Geometric_data_local_eqn(0) {}
 
 /// Broken copy constructor
 ElementWithMovingNodes(const ElementWithMovingNodes&) 
  { 
   BrokenCopy::broken_copy("ElementWithMovingNodes");
  } 
 
 /// Broken assignment operator
 void operator=(const ElementWithMovingNodes&) 
  {
   BrokenCopy::broken_assign("ElementWithMovingNodes");
  }

 /// Virtual destructor (clean up and allocated memory)
 virtual ~ElementWithMovingNodes()
  {
   if(Geometric_data_local_eqn)
    {
     delete[] Geometric_data_local_eqn[0];
     delete[] Geometric_data_local_eqn;
    }
  }

 /// \short Return pointer to the i-th Data item that the element's 
 /// shape depends on 
 //Data* geom_data_pt(const unsigned& i) {return Geom_data_pt[i];}


 /// \short Return the local equation number corresponding to the i-th
 /// value at the n-th geometric data object.
 inline int geometric_data_local_eqn(const unsigned &n, const unsigned &i)
   {
#ifdef RANGE_CHECKING
   unsigned n_data = Geom_data_pt.size();
   if(n >= n_data)
    {
     std::ostringstream error_message;
     error_message << "Range Error:  Data number " << s
                   << " is not in the range (0,"
                   << n_data-1 << ")";
     throw OomphLibError(error_message.str(),
                         "SpineElement::geometric_data_local_eqn()",
                         OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     unsigned n_value = Geom_data_pt[n]->nvalue();
     if(i >= n_value)
      {
       std::ostringstream error_message;
       error_message << "Range Error: value " << i << " at data " << s 
                     << " is not in the range (0,"
                     << n_value -1 << ")";
       throw OomphLibError(error_message.str(),
                           "SpineElement::geometric_data_local_eqn()",
                           OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
#ifdef PARANOID
   //Check that the equations have been allocated
   if(Geometric_data_local_eqn==0)
    {
     throw OomphLibError(
      "Spine local equation numbers have not been allocated",
      "SpineElement::geometric_data_local_eqn()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   return Geometric_data_local_eqn[n][i];
  }


 ///Return a set of all geometric data associated with the element
 void assemble_set_of_all_geometric_data(std::set<Data*> &unique_geom_data_pt);
 
/// \short Specify Data that affects the geometry of the element 
/// by filling in the entris of the set  to the set.
/// (This functionality is required in FSI problems). 
 void identify_geometric_data(std::set<Data*> &geometric_data_pt) 
  {
   //Loop over the node update data and add to the set
   const unsigned n_geom_data = Geom_data_pt.size();
   for(unsigned n=0;n<n_geom_data;n++)
    {
     geometric_data_pt.insert(Geom_data_pt[n]);
    }
  }


  protected:

 /// Construct the vector of (unique) geometric data
 void complete_setup_of_dependencies();

/// Assign local equation numbers for the spines in the element
 virtual void assign_all_generic_local_eqn_numbers();

 ///Calculate the contributions to the jacobian matrix from the
 ///geometric data using finite differences. This version
 ///assumes that the (full) residuals vector has already been calculated
 void fill_in_jacobian_from_geometric_data_by_fd(
  Vector<double> &residuals, DenseMatrix<double> &jacobian);

 ///Calculate the contributions to the jacobian matrix from the
 ///geometric data using finite differences. This version computes
 ///the residuals vector before calculating the jacobian terms
 void fill_in_jacobian_from_geometric_data_by_fd(DenseMatrix<double> &jacobian)
  {
   //Allocate storage for the residuals
   const unsigned n_dof = ndof();
   Vector<double> residuals(n_dof);
   //Get the residuals for the entire element
   get_residuals(residuals);
   //Call the jacobian calculation
   fill_in_jacobian_from_geometric_data_by_fd(residuals,jacobian);
  }

 
 /// \short Vector that stores pointers to all Data that affect the 
 /// node update operations, i.e. the variables that can affect
 /// the position of the node. 
 Vector<Data*> Geom_data_pt;

private: 
 
 /// \short Return the number of geometric data upon which the shape
 /// of the element depends
 unsigned ngeom_data() const {return Geom_data_pt.size();}


 /// \short Array to hold local eqn number information for the
 /// Spine variables
 int **Geometric_data_local_eqn;



};


//===============================================================
/// Specific implementation of the class
//==============================================================
template<class ELEMENT,class NODE_TYPE>
 class ElementWithSpecificMovingNodes : public ELEMENT,
 public ElementWithMovingNodes
{
  public:

 /// Constructor, call the constructor of the base element
 ElementWithSpecificMovingNodes() : ELEMENT(), ElementWithMovingNodes() {}
 
 /// Constructor used for spine face elements
 ElementWithSpecificMovingNodes(FiniteElement* const &element_pt, 
                                const int &face_index) : 
  ELEMENT(element_pt, face_index), ElementWithMovingNodes() {}

 /// Empty Destructor, 
 ~ElementWithSpecificMovingNodes() {} 


 /// Overload the node assignment routine to assign spine nodes
 Node* construct_node(const unsigned &n)
  {
   //Assign a spine node to the local node pointer
   //The dimension and number of values are taken from internal element data
   //The number of timesteps to be stored comes from the problem!
   this->node_pt(n) = 
    new NODE_TYPE(this->nodal_dimension(),this->nnodal_position_type(),
                  this->required_nvalue(n));
   //Now return a pointer to the node, so that the mesh can find it
   return this->node_pt(n);
  }

 /// Overloaded node allocation for unsteady problems
 Node* construct_node(const unsigned &n, 
                      TimeStepper* const &time_stepper_pt)
  {
   //Assign a node to the local node pointer
   //The dimension and number of values are taken from internal element data
   //The number of timesteps to be stored comes from the problem!
   this->node_pt(n) = 
    new NODE_TYPE(time_stepper_pt,this->nodal_dimension(),
                  this->nnodal_position_type(),
                  this->required_nvalue(n));
   //Now return a pointer to the node, so that the mesh can find it
   return this->node_pt(n);
  }

 /// Overload the node assignment routine to assign spine boundary nodes
 Node* construct_boundary_node(const unsigned &n)
  {
   //Assign a spine node to the local node pointer
   //The dimension and number of values are taken from internal element data
   //The number of timesteps to be stored comes from the problem!
   this->node_pt(n) = 
    new BoundaryNode<NODE_TYPE>(this->nodal_dimension(),
                                this->nnodal_position_type(),
                                this->required_nvalue(n));
   //Now return a pointer to the node, so that the mesh can find it
   return this->node_pt(n);
  }

 /// Overloaded boundary node allocation for unsteady problems
 Node* construct_boundary_node(const unsigned &n, 
                               TimeStepper* const &time_stepper_pt)
  {
   //Assign a node to the local node pointer
   //The dimension and number of values are taken from internal element data
   //The number of timesteps to be stored comes from the problem!
   this->node_pt(n) = 
    new BoundaryNode<NODE_TYPE>(time_stepper_pt,this->nodal_dimension(),
                                this->nnodal_position_type(),
                                this->required_nvalue(n));
   //Now return a pointer to the node, so that the mesh can find it
   return this->node_pt(n);
  }


 /// \short Complete the setup of additional dependencies. Overloads
 /// empty virtual function in GeneralisedElement to determine the "geometric 
 /// Data", i.e. the Data that affects the element's shape.
 /// This function is called (for all elements) at the very beginning of the
 /// equation numbering procedure to ensure that all dependencies
 /// are accounted for.
 void complete_setup_of_dependencies()
  {

   // Call function of underlying element
   ELEMENT::complete_setup_of_dependencies();
   //Call function of the element with moving nodes
   ElementWithMovingNodes::complete_setup_of_dependencies();
  }


 /// Assign local equation numbers for the spines in the element
 void assign_all_generic_local_eqn_numbers()
 {
  // Call the generic local equation numbering scheme of the ELEMENT
  ELEMENT::assign_all_generic_local_eqn_numbers();
  ElementWithMovingNodes::assign_all_generic_local_eqn_numbers();
 }

 /// Compute the element's residuals vector and jacobian matrix
 void get_jacobian(Vector<double> &residuals,
                   DenseMatrix<double> &jacobian)
  {
   ///Call the element's get jacobian function
   ELEMENT::get_jacobian(residuals,jacobian);
   //Now call the additional spine jacobian terms (finite differences)
   this->fill_in_jacobian_from_geometric_data_by_fd(jacobian);
  }

 /// Add the contribution to the jacobian and residuals vector
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                   DenseMatrix<double> &jacobian)
  {
   //Call the element's add to Jacobian function 
   ELEMENT::fill_in_contribution_to_jacobian(residuals,jacobian);
   //Now call the additional spine jacobian terms (finite differences)
   this->fill_in_jacobian_from_geometric_data_by_fd(jacobian);
  }

};


}
#endif

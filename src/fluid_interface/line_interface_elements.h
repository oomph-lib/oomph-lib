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
//Header file for specific (one-dimensional) free surface elements
//Include guards, to prevent multiple includes
#ifndef OOMPH_LINE_INTERFACE_ELEMENTS_HEADER
#define OOMPH_LINE_INTERFACE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

//OOMPH-LIB headers
#include "../generic/Qelements.h"
#include "../generic/spines.h"
#include "../generic/hijacked_elements.h"
#include "interface_elements.h"

namespace oomph
{

//======================================================================
/// One-dimensional fluid interface elements that are used with a spine mesh,
/// i.e. the mesh deformation is handled by Kistler & Scriven's "method
/// of spines". These elements are FaceElements attached to 2D bulk
/// fluid elements and the particular type of fluid element is passed 
/// as a template parameter to the element. It 
/// shouldn't matter whether the passed 
/// element is the underlying (fixed) element or the templated 
/// SpineElement<Element>.
/// Optionally, an external pressure may be specified, which must be
/// passed to the element as external data. If there is no such object,
/// the external pressure is assumed to be zero. 
//======================================================================
template<class ELEMENT>
class SpineLineFluidInterfaceElement : 
public virtual Hijacked<SpineElement<FaceGeometry<ELEMENT> > >, 
public virtual LineFluidInterfaceElement 
{
  private:

 /// \short In spine elements, the kinematic condition is the equation 
 /// used to determine the unknown spine heights. Overload the
 /// function accordingly
 int kinematic_local_eqn(const unsigned &n) 
  {return this->spine_local_eqn(n);}


 /// \short Hijacking the kinematic condition corresponds to hijacking the
 /// variables associated with the spine heights --- only used when
 /// strong imposition of the contact angle condition is required.
 void hijack_kinematic_conditions(const Vector<unsigned> &bulk_node_number)
 {
  //Loop over all the node numbers that are passed in
  for(Vector<unsigned>::const_iterator it=bulk_node_number.begin();
      it!=bulk_node_number.end();++it)
   {
    //Hijack the spine heights. (and delete the returned data object)
    delete this->hijack_nodal_spine_value(*it,0);
   }
 }
 
  public:
 
 /// \short Constructor, the arguments are a pointer to the  "bulk" element 
 /// and the index of the face to be created.
 SpineLineFluidInterfaceElement(FiniteElement* const &element_pt, 
                                const int &face_index) : 
  Hijacked<SpineElement<FaceGeometry<ELEMENT> > >(), 
  LineFluidInterfaceElement()
   {
    //Attach the geometrical information to the element, by
    //making the face element from the bulk element
    element_pt->build_face_element(face_index,this);
    
    //Find the index at which the velocity unknowns are stored 
    //from the bulk element
    ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
    this->U_index_interface.resize(2);
    for(unsigned i=0;i<2;i++)
     {
      this->U_index_interface[i] = cast_element_pt->u_index_nst(i);
     }
   } //End of constructor



 /// Calculate the contribution to the residuals and the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
   
   //Call the generic routine to evaluate shape derivatives
   this->fill_in_jacobian_from_geometric_data(jacobian);
  } //End of jacobian contribution
 

 /// \short
 /// Helper function to calculate the additional contributions
 /// to be added at each integration point. Empty for this
 /// implemenetation
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif,
  const DShape &dpsifds,
  const Vector<double> &interpolated_x, 
  const Vector<double> &interpolated_n, 
  const double &W, 
  const double &J)
 {
 }
 
 /// Overload the output function
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

 /// Output the element
 void output(std::ostream &outfile, const unsigned &n_plot)
  {LineFluidInterfaceElement::output(outfile,n_plot);}

 ///Overload the C-style output function
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}

 ///C-style Output function
 void output(FILE* file_pt, const unsigned &n_plot)
  {LineFluidInterfaceElement::output(file_pt,n_plot);}


 /// \short Create an "bounding" element (here actually a 1D point element
 /// of type SpinePointFluidInterfaceBoundingElement<ELEMENT> that allows
 /// the application of a contact angle boundary condition on the
 /// the specified face.
 virtual FluidInterfaceBoundingElement* make_bounding_element(
  const int &face_index)
 {
  //Create a temporary pointer to the appropriate FaceElement
  SpinePointFluidInterfaceBoundingElement<ELEMENT> *face_el_pt = 
   new SpinePointFluidInterfaceBoundingElement<ELEMENT>;
  
  //Attach the geometrical information to the new element
  this->build_face_element(face_index,face_el_pt);
  
  //Set the index at which the unknowns are stored from the element
  face_el_pt->u_index_interface_boundary() = this->U_index_interface;
  
  //Set the value of the nbulk_value, the node is not resized
  //in this problem, so it will just be the actual nvalue
  face_el_pt->nbulk_value(0) = face_el_pt->node_pt(0)->nvalue();
  
  //Set of unique geometric data that is used to update the bulk,
  //but is not used to update the face
  std::set<Data*> unique_additional_geom_data;

  //Get all the geometric data for this (bulk) element
  this->assemble_set_of_all_geometric_data(unique_additional_geom_data);
  
  //Now assemble the set of geometric data for the face element
  std::set<Data*> unique_face_geom_data_pt;
  face_el_pt->assemble_set_of_all_geometric_data(unique_face_geom_data_pt);

  //Erase the face geometric data from the additional data
  for(std::set<Data*>::iterator it=unique_face_geom_data_pt.begin();
      it!=unique_face_geom_data_pt.end();++it)
   {unique_additional_geom_data.erase(*it);}
  
  //Finally add all unique additional data as external data
  for(std::set<Data*>::iterator it = unique_additional_geom_data.begin();
      it!= unique_additional_geom_data.end();++it)
   {
    face_el_pt->add_external_data(*it);
   }
  
  //Return the value of the pointer
  return face_el_pt;
 }
 
};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//=======================================================================
/// One-dimensional interface elements that are used when the mesh
/// deformation is handled by a set of equations that modify the nodal
/// positions. These elements are FaceElements attached to 2D bulk fluid 
/// elements and the fluid element is passed as a template parameter to 
/// the element.
/// Optionally an external pressure may be specified, which must be
/// passed to the element as external data. The default value of the external
/// pressure is zero.
//=======================================================================
template<class ELEMENT>
class ElasticLineFluidInterfaceElement : 
public virtual Hijacked<FaceGeometry<ELEMENT> >, 
public LineFluidInterfaceElement
{
  private:

 /// \short ID of the Lagrange Lagrange multiplier (in the collection of nodal
 /// values accomodated by resizing)
 unsigned Id;
 
 /// \short Equation number of the kinematic BC associated with node j.
 /// (This is the equation for the Lagrange multiplier)
 int kinematic_local_eqn(const unsigned &j)
 {
  // Get the index of the nodal value associated with Lagrange multiplier
  const unsigned lagr_index=
   dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
   index_of_first_value_assigned_by_face_element(Id);
  
  // Return nodal value
  return this->nodal_local_eqn(j,lagr_index);
 } 


 /// \short Hijacking the kinematic condition corresponds to hijacking the
 /// variables associated with the Lagrange multipliers that are assigned
 /// on construction of this element.
 void hijack_kinematic_conditions(const Vector<unsigned> &bulk_node_number)
  {
   //Loop over all the nodes that are passed in
   for(Vector<unsigned>::const_iterator it=bulk_node_number.begin();
       it!=bulk_node_number.end();++it)
    {
     //Get the index associated with the Id for each node
     //(the Lagrange multiplier)
     unsigned n_lagr = dynamic_cast<BoundaryNodeBase*>(node_pt(*it))->
      index_of_first_value_assigned_by_face_element(Id);
  
     //Hijack the appropriate value and delete the returned Node
     delete this->hijack_nodal_value(*it,n_lagr);  
    }
  }
 
  public:

 /// \short The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                          const unsigned &i) const 
  {return FaceElement::zeta_nodal(n,k,i);}     
 

 /// \short Constructor, pass a pointer to the bulk element and the face index
 /// of the bulk element to which the element is to be attached to.
 /// The optional identifier can be used
 /// to distinguish the additional nodal value (Lagr mult) created by 
 /// this element from those created by other FaceElements.
 ElasticLineFluidInterfaceElement(FiniteElement* const &element_pt, 
                                  const int &face_index,
                                  const unsigned &id=0) : 
  // Final optional argument to specify id for lagrange multiplier
  FaceGeometry<ELEMENT>(), LineFluidInterfaceElement(), Id(id)
  {
   //Attach the geometrical information to the element
   //This function also assigned nbulk_value from required_nvalue of the
   //bulk element
   element_pt->build_face_element(face_index,this);
   
   //Find the index at which the velocity unknowns are stored 
   //from the bulk element
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   this->U_index_interface.resize(2);
   for(unsigned i=0;i<2;i++)
    {
     this->U_index_interface[i] = cast_element_pt->u_index_nst(i);
    }

   //Read out the number of nodes on the face
   unsigned n_node_face = this->nnode();

   //Set the additional data values in the face
   //There is one additional values at each node --- the lagrange multiplier
   Vector<unsigned> additional_data_values(n_node_face);
   for(unsigned i=0;i<n_node_face;i++)
    {
     additional_data_values[i] = 1;
    }
   
   // Now add storage for Lagrange multipliers and set the map containing 
   // the position of the first entry of this face element's 
   // additional values.
   add_additional_values(additional_data_values,id);
   
  } //End of constructor
  
 /// Return the Lagrange multiplier at local node j
 double &lagrange(const unsigned &j)
  {
   // Get the index of the nodal value associated with Lagrange multiplier
   unsigned lagr_index=dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
    index_of_first_value_assigned_by_face_element(Id);
   return *node_pt(j)->value_pt(lagr_index); 
  }
 
 /// Fill in contribution to residuals and Jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
 {
  //Call the generic routine with the flag set to 1
  fill_in_generic_residual_contribution_interface(residuals,jacobian,1);

  //Call the generic finite difference routine for the solid variables
  this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
 }
 
 /// Overload the output function
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}
 
 /// Output the element
 void output(std::ostream &outfile, const unsigned &n_plot)
 {LineFluidInterfaceElement::output(outfile,n_plot);}
 
 ///Overload the C-style output function
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}
 
 ///C-style Output function
 void output(FILE* file_pt, const unsigned &n_plot)
 {LineFluidInterfaceElement::output(file_pt,n_plot);}
 
 
 /// \short Helper function to calculate the additional contributions
 /// to be added at each integration point. This deals with 
 /// Lagrange multiplier contribution
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif,
  const DShape &dpsifds,
  const Vector<double> &interpolated_x, 
  const Vector<double> &interpolated_n, 
  const double &W, 
  const double &J)
 {
  //Loop over the shape functions
  unsigned n_node = this->nnode();

  // Assemble Lagrange multiplier in loop over the shape functions
  double interpolated_lagrange = 0.0;
  for(unsigned l=0;l<n_node;l++)
   {
    //Note same shape functions used for lagrange multiplier field
    interpolated_lagrange += lagrange(l)*psif[l];
   }
  
  int local_eqn=0, local_unknown = 0;
  
  //Loop over the shape functions to assemble contributions
  for(unsigned l=0;l<n_node;l++)
   {
    //Loop over the directions
    for(unsigned i=0;i<2;i++)
     {
      //Now using the same shape functions for the elastic equations,
      //so we can stay in the loop
      local_eqn = this->position_local_eqn(l,0,i);
      if(local_eqn >= 0)
       {
        //Add in a "Lagrange multiplier"
        residuals[local_eqn] -= 
         interpolated_lagrange*interpolated_n[i]*psif[l]*W*J;
        
        //Do the Jacobian calculation
        if(flag)
         {
          //Loop over the nodes 
          for(unsigned l2=0;l2<n_node;l2++)
           {
            //Derivatives w.r.t. solid positions will be handled by FDing
            //That leaves the "lagrange multipliers" only
            local_unknown = kinematic_local_eqn(l2);
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) -=
               psif[l2]*interpolated_n[i]*psif[l]*W*J;
             }
           }
         } //End of Jacobian calculation
       } 
     } 
   } //End of loop over shape functions
 }
 
 
 /// \short Create an "bounding" element (here actually a 1D point element
 /// of type ElasticPointFluidInterfaceBoundingElement<ELEMENT> that allows
 /// the application of a contact angle boundary condition on the
 /// the specified face. 
 virtual FluidInterfaceBoundingElement* make_bounding_element(
  const int &face_index)
 {
  //Create a temporary pointer to the appropriate FaceElement
  ElasticPointFluidInterfaceBoundingElement<ELEMENT> *face_el_pt = 
   new ElasticPointFluidInterfaceBoundingElement<ELEMENT>;

  //Attach the geometrical information to the new element
  this->build_face_element(face_index,face_el_pt);

  //Set the index at which the unknowns are stored from the element
  face_el_pt->u_index_interface_boundary() = this->U_index_interface;

  //Set the value of the nbulk_value, the node is not resized
  //in this problem, so it will just be the actual nvalue 
  face_el_pt->nbulk_value(0) = this->node_pt(0)->nvalue();

  //Pass the ID down
  face_el_pt->set_id(Id);

  //Find the nodes
  std::set<SolidNode*> set_of_solid_nodes;
  unsigned n_node = this->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    set_of_solid_nodes.insert(static_cast<SolidNode*>(this->node_pt(n)));
   }

  //Delete the nodes from the face
  n_node = face_el_pt->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    set_of_solid_nodes.erase(static_cast<SolidNode*>(face_el_pt->node_pt(n)));
   }

  //Now add these as external data
  for(std::set<SolidNode*>::iterator it=set_of_solid_nodes.begin();
      it!=set_of_solid_nodes.end();++it)
   {
    face_el_pt->add_external_data((*it)->variable_position_pt());
   }
  
  //Return the value of the pointer
  return face_el_pt;
 }
 
};

}

#endif







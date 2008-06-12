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
///This header file contains elements that combine two element types in
//a generic way. 

#ifndef OOMPH_PSEUDO_SOLID_REMESH_ELEMENTS_HEADER
#define OOMPH_PSEUDO_SOLID_REMESH_ELEMENTS_HEADER

#include "elements.h"

namespace oomph
{

//==========================================================================
/// A templated class that permits combination two different element types,
/// for the solution of problems in deforming domains. The first template
/// paremter BASIC is the standard element and the second SOLID solves
/// the equations that are used to control the mesh deformation.
//==========================================================================
template<class BASIC, class SOLID>
 class PseudoSolidNodeUpdateElement 
 : public virtual BASIC, public virtual SOLID

{
  public: 

 ///Constructor, call the BASIC and SOLID elements' constructors
 PseudoSolidNodeUpdateElement() : BASIC(), SOLID() {}

 /// \short The required number of values is the sum of the two
 unsigned required_nvalue(const unsigned &n) const
  {return BASIC::required_nvalue(n) + SOLID::required_nvalue(n);}
 
 /// \short We assume that the solid stuff is stored at the end of
 /// the nodes, i.e. its index is the number of continuously interplated
 /// values in the BASIC equations.
 int solid_p_nodal_index() const
  {
   //At the moment, we can't handle this case in generality so throw an
   //error if the solid pressure is stored at the nodes
   if(SOLID::solid_p_nodal_index() >= 0)
    {
     throw OomphLibError(
      "Cannot handle (non-refineable) continuous solid pressure interpolation",
      "PseudoSolidNodeUpdateElement::solid_p_nodal_index()",
      OOMPH_EXCEPTION_LOCATION);
    }

   return SOLID::solid_p_nodal_index();
  }

 /// \short Final override for the residuals function. Contributions are
 /// added from both underlying element types
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the basic equations first
   BASIC::fill_in_contribution_to_residuals(residuals);
   //Add the solid equations contribution
   SOLID::fill_in_contribution_to_residuals(residuals);
  }
 
 /// \short Final override for jacobian function: Contributions are
 /// included from both the underlying element types
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the basic equations first
   BASIC::fill_in_contribution_to_jacobian(residuals,jacobian);
   //Call the solid equations
   SOLID::fill_in_contribution_to_jacobian(residuals,jacobian);
  }
 
 ///Overload the output function: Call that of the basic element
 void output(std::ostream &outfile) {BASIC::output(outfile);}

 /// \short Output function: Plot at n_p plot points using the basic element's
 /// output function
 void output(std::ostream &outfile, const unsigned &n_p)
  {BASIC::output(outfile,n_p);}

 ///Overload the output function: Call that of the basic element
 void output(FILE* file_pt) {BASIC::output(file_pt);}

 ///Output function is just the same as the basic equations
 void output(FILE* file_pt, const unsigned &n_p)
  {BASIC::output(file_pt,n_p);}


 /// \short The number of "blocks" that degrees of freedom in this element
 /// are sub-divided into.
 /// This is needed as a final overload in cases where both fluid and
 /// solid elements are block preconditionable. However, we break
 /// it here because it isn't obvious which classification we
 /// should use. This forces the user to re-implement this function
 /// if it's used
 unsigned nblock_types()
  {
   throw OomphLibError(
    "nblock_types() is deliberately broken. Provide your own final overload!",
    "PseudoSolidNodeUpdateElement::block_types()",
    OOMPH_EXCEPTION_LOCATION);
   
   // dummy return
   return 1;
  }

 /// \short Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "block" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)
 /// This is needed as a final overload in cases where both fluid and
 /// solid elements are block preconditionable. However, we break
 /// it here because it isn't obvious which classification we
 /// should use. This forces the user to re-implement this function
 /// if it's used
 void get_block_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& block_lookup_list)
  {
   throw OomphLibError(
    "get_block_numbers_for_unknowns() is deliberately broken. Provide your own final overload!",
    "PseudoSolidNodeUpdateElement::get_block_numbers_for_unknowns()",
    OOMPH_EXCEPTION_LOCATION);
  }


};

///Explicit definition of the face geometry of these elements
template<class BASIC, class SOLID>
class FaceGeometry<PseudoSolidNodeUpdateElement<BASIC,SOLID> >:
public virtual FaceGeometry<SOLID> 
{
  public:

 /// \short Constuctor calls the constructor of the SolidQElement
 /// (Only the Intel compiler seems to need this!)
 FaceGeometry() : FaceGeometry<SOLID>() {}

  protected:

};

///Explicit definition of the face geometry of these elements
template<class BASIC, class SOLID>
class FaceGeometry<FaceGeometry<PseudoSolidNodeUpdateElement<BASIC,SOLID> > >:
public virtual FaceGeometry<FaceGeometry<SOLID> >
{
  public:

 /// \short Constuctor calls the constructor of the SolidQElement
 /// (Only the Intel compiler seems to need this!)
 FaceGeometry() : FaceGeometry<FaceGeometry<SOLID> >() {}

  protected:

};


//===================================================================
/// Refineable version of the PseudoSolidNodeUpdateELement
//===================================================================
template<class BASIC, class SOLID>
class RefineablePseudoSolidNodeUpdateElement : public virtual BASIC, 
 public virtual SOLID
{

  public: 
 ///Constructor, call the BASIC and SOLID elements' constructors
 RefineablePseudoSolidNodeUpdateElement() : 
  //Need to sort this one out (at the moment both element must be quads)
  RefineableElement(),
  //RefineableQElement<2>(), 
  BASIC(), SOLID() {} 

 /// \short The required number of values is the sum of the two
 unsigned required_nvalue(const unsigned &n) const
  {return BASIC::required_nvalue(n) + SOLID::required_nvalue(n);}

 /// \short The number of continuously interpolated values is the 
 /// sum of the SOLID and BASIC values
 unsigned ncont_interpolated_values() const
  {return BASIC::ncont_interpolated_values() + 
    SOLID::ncont_interpolated_values();}

 /// \short We assume that the solid stuff is stored at the end of
 /// the nodes, i.e. its index is the number of continuously interplated
 /// values in the BASIC equations.
 int solid_p_nodal_index() const
  {
   //Find the index in the solid
   int solid_p_index = SOLID::solid_p_nodal_index();
   //If there is a solid pressure at the nodes, return the
   //index after all the BASIC stuff
   if(solid_p_index >= 0)
    {return BASIC::ncont_interpolated_values() + 
      SOLID::solid_p_nodal_index();}
   else {return solid_p_index;}
  }

 ///\short Final override for residuals function: adds contributions
 ///from both underlying element types
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the basic equations first
   BASIC::fill_in_contribution_to_residuals(residuals);
   //Call the solid equations
   SOLID::fill_in_contribution_to_residuals(residuals);
  }

 ///\short Final override for jacobian function: Calls get_jacobian() for 
 /// both of the underlying element types
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the basic equations first
   BASIC::fill_in_contribution_to_jacobian(residuals,jacobian);
   //Call the solid equations
   SOLID::fill_in_contribution_to_jacobian(residuals,jacobian);
  }
 
 /// \short Final override for the assign__additional_local_eqn_numbers():
 ///  Call the version for the BASIC element
 void assign_additional_local_eqn_numbers()
  {
   BASIC::assign_additional_local_eqn_numbers();
   SOLID::assign_additional_local_eqn_numbers();
  }

 /// Call rebuild_from_sons() for both of the underlying element types
 void rebuild_from_sons(Mesh* &mesh_pt) 
  {
   BASIC::rebuild_from_sons(mesh_pt);
   SOLID::rebuild_from_sons(mesh_pt);
  }

 /// Call get_interpolated_values(...) for both of the underlying element types
 void get_interpolated_values(const unsigned& t, const Vector<double>&s,
                              Vector<double>& values)
  {
   Vector<double> basic_values;
   BASIC::get_interpolated_values(t,s,basic_values);
   Vector<double> solid_values;
   SOLID::get_interpolated_values(t,s,solid_values);

    //Now add the basic value first
   for(Vector<double>::iterator it=basic_values.begin();
       it!=basic_values.end();++it)
    {
     values.push_back(*it);
    }
   //Then the solid
   for(Vector<double>::iterator it=solid_values.begin();
       it!=solid_values.end();++it)
    {
     values.push_back(*it);
    }
  }

 
 /// Call get_interpolated_values(...) for both of the underlying element types
 void get_interpolated_values(const Vector<double> &s, Vector<double> & values)
  { 
   Vector<double> basic_values;
   BASIC::get_interpolated_values(s,basic_values);
   Vector<double> solid_values;
   SOLID::get_interpolated_values(s,solid_values);

   //Now add the basic value first
   for(Vector<double>::iterator it=basic_values.begin();
       it!=basic_values.end();++it)
    {
     values.push_back(*it);
    }
   //Then the solid
   for(Vector<double>::iterator it=solid_values.begin();
       it!=solid_values.end();++it)
    {
     values.push_back(*it);
    }
  }

 /// \short We must compose the underlying interpolating nodes from
 /// the BASIC and SOLID equations, the BASIC ones are first
 Node* interpolating_node_pt(const unsigned &n,
                             const int &value_id) 
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {return BASIC::interpolating_node_pt(n,value_id);}
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {return SOLID::interpolating_node_pt(n,(value_id-n_basic_values));}
  }

 /// \short The pressure nodes are the corner nodes, so when value_id==0,
 /// the fraction is the same as the 1d node number, 0 or 1.
 double local_one_d_fraction_of_interpolating_node(const unsigned &n1d,
                                                   const unsigned &i, 
                                                   const int &value_id)
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {
     return BASIC::local_one_d_fraction_of_interpolating_node(n1d,i,value_id);
    }
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {
     return 
      SOLID::local_one_d_fraction_of_interpolating_node(
       n1d,i,(value_id-n_basic_values));
    }
  }


 /// \short The velocity nodes are the same as the geometric nodes. The
 /// pressure nodes must be calculated by using the same methods as
 /// the geometric nodes, but by recalling that there are only two pressure
 /// nodes per edge.
 Node* get_interpolating_node_at_local_coordinate(const Vector<double> &s,   
                                                  const int &value_id)
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {
     return BASIC::get_interpolating_node_at_local_coordinate(s,value_id);
    }
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {
     return 
      SOLID::get_interpolating_node_at_local_coordinate(
       s,(value_id-n_basic_values));
    }
  }

 /// \short The number of 1d pressure nodes is 2, otherwise we have
 /// the positional nodes
 unsigned ninterpolating_node_1d(const int &value_id)
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {
     return BASIC::ninterpolating_node_1d(value_id);
    }
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {
     return SOLID::ninterpolating_node_1d((value_id-n_basic_values));
    }
  }

 /// \short The number of pressure nodes is 2^DIM. The number of 
 /// velocity nodes is the same as the number of geometric nodes.
 unsigned ninterpolating_node(const int &value_id)
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {
     return BASIC::ninterpolating_node(value_id);
    }
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {
     return SOLID::ninterpolating_node((value_id-n_basic_values));
    }
  }
 
 /// \short The basis interpolating the pressure is given by pshape().
 //// The basis interpolating the velocity is shape().
 void interpolating_basis(const Vector<double> &s,
                          Shape &psi,
                          const int &value_id) const
  {
   //Find the number of interpolated values in the BASIC equations
   int n_basic_values = BASIC::ncont_interpolated_values();
   //If the id is below this number, we call the BASIC functon
   if(value_id < n_basic_values)
    {
     return BASIC::interpolating_basis(s,psi,value_id);
    }
   //Otherwise it's the solid and its value_id is the the current
   //it minus n_basic_values
   else
    {
     return SOLID::interpolating_basis(s,psi,(value_id-n_basic_values));
    }
  }


 /// \short Number of 'flux' terms for Z2 error estimation: Error estimation
 /// is based on error in BASIC element
 unsigned num_Z2_flux_terms() {return BASIC::num_Z2_flux_terms();}
 
 /// 'Flux' vector for Z2 error estimation: Error estimation
 /// is based on error in BASIC element
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {
   BASIC::get_Z2_flux(s,flux);
  }
 
 /// \short Perform additional hanging node procedures for variables
 /// that are not interpolated by all nodes. Done for both of the
 /// underlying element types.
 void further_setup_hanging_nodes()
  {
   BASIC::further_setup_hanging_nodes();
   SOLID::further_setup_hanging_nodes();
  }


 /// \short Build function: Call the one for the SOLID element since it
 /// calls the one basic build function automatically.
 void build(Mesh*& mesh_pt, Vector<Node*>& new_node_pt,
            bool& was_already_built,
            std::ofstream& new_nodes_file)
  {
   SOLID::build(mesh_pt,new_node_pt,
                was_already_built,new_nodes_file);
  }


 /// \short Build function: Call the one for the SOLID element since it
 /// calls the one basic build function automatically.
 void build(Mesh*& mesh_pt, Vector<Node*>& new_node_pt,
            bool& was_already_built)
  {
   SOLID::build(mesh_pt,new_node_pt,was_already_built);
  }

 ///  Further build: Done for both of the
 /// underlying element types.
 void further_build()
  {
   BASIC::further_build();
   SOLID::further_build();
  }
 

 /// \short Number of vertex nodes in the element hierher 
 //// \todo TO BE MOVED TO BASE CLASS
 unsigned nvertex_node() const
  {return BASIC::nvertex_node();}

 /// \short Pointer to the j-th vertex node in the element hierher \todo 
 /// TO BE MOVED TO BASE CLASS
 Node* vertex_node_pt(const unsigned& j) const
  {return BASIC::vertex_node_pt(j);}

 /// \short Order of recovery shape functions for Z2 error estimation: Done
 /// for BASIC element since it determines the refinement
 unsigned nrecovery_order() {return BASIC::nrecovery_order();}

 ///Overload the output function: Use that of the BASIC element
 void output(std::ostream &outfile) {BASIC::output(outfile);}

 ///Output function, plotting at n_p points: Use that of the BASIC element
 void output(std::ostream &outfile, const unsigned &n_p)
  {BASIC::output(outfile,n_p);}

 ///Overload the output function: Use that of the BASIC element
 void output(FILE* file_pt) {BASIC::output(file_pt);}

 ///Output function: Use that of the BASIC element
  void output(FILE* file_pt, const unsigned &n_p)
  {BASIC::output(file_pt,n_p);}


 /// \short The number of "blocks" that degrees of freedom in this element
 /// are sub-divided into.
 /// This is needed as a final overload in cases where both fluid and
 /// solid elements are block preconditionable. However, we break
 /// it here because it isn't obvious which classification we
 /// should use. This forces the user to re-implement this function
 /// if it's used
 unsigned nblock_types()
  {
   throw OomphLibError(
    "nblock_types() is deliberately broken. Provide your own final overload!",
    "RefineablePseudoSolidNodeUpdateElement::block_types()",
    OOMPH_EXCEPTION_LOCATION);
   
   // dummy return
   return 1;
  }

 /// \short Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "block" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)
 /// This is needed as a final overload in cases where both fluid and
 /// solid elements are block preconditionable. However, we break
 /// it here because it isn't obvious which classification we
 /// should use. This forces the user to re-implement this function
 /// if it's used
 void get_block_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& block_lookup_list)
  {
   throw OomphLibError(
    "get_block_numbers_for_unknowns() is deliberately broken. Provide your own final overload!",
    "RefineablePseudoSolidNodeUpdateElement::get_block_numbers_for_unknowns()",
    OOMPH_EXCEPTION_LOCATION);
  }



};

///Explicit definition of the face geometry of these elements
template<class BASIC, class SOLID>
class FaceGeometry<RefineablePseudoSolidNodeUpdateElement<BASIC,SOLID> >:
public virtual FaceGeometry<SOLID> 
{
  public:
 /// \short Constructor calls the constructor of the SolidQElement
 /// (Only the Intel compiler seems to need this!)
 FaceGeometry() : FaceGeometry<SOLID>() {}

  protected:
};

///Explicit definition of the face geometry of these elements
template<class BASIC, class SOLID>
class FaceGeometry<FaceGeometry<
 RefineablePseudoSolidNodeUpdateElement<BASIC,SOLID> > >:
public virtual FaceGeometry<FaceGeometry<SOLID> >
{
  public:

 /// \short Constuctor calls the constructor of the SolidQElement
 /// (Only the Intel compiler seems to need this!)
 FaceGeometry() : FaceGeometry<FaceGeometry<SOLID> >() {}

  protected:

};


}

#endif

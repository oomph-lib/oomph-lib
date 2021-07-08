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
#ifndef OOMPH_REFINEABLE_QUAD_ELEMENT_HEADER
#define OOMPH_REFINEABLE_QUAD_ELEMENT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//oomph-lib headers
#include "quadtree.h"
#include "macro_element.h"
#include "refineable_elements.h"
#include "Qelements.h"

namespace oomph
{

 //Forward definition for mesh.
 class Mesh;

//=======================================================================
/// Refineable version of QElement<2,NNODE_1D>.
///
/// Refinement is performed by quadtree procedures. When the element is
/// subdivided, the geometry of its sons is established by calls
/// to their father's
/// \code  get_x(...) \endcode
/// function which refers to
/// - the father element's geometric FE mapping  (this
///   is the default)
/// . 
///  or
/// - to a MacroElement 's MacroElement::macro_map (if the pointer 
/// to the macro element is non-NULL)
/// 
/// The class provides a generic RefineableQElement<2>::build() function 
/// which deals with generic
/// isoparametric QElements in which all values are associated with
/// nodes. The RefineableQElement<2>::further_build() function provides
/// an interface for any element-specific non-generic build operations.
/// 
//=======================================================================
template<>
class RefineableQElement<2> : public virtual RefineableElement, 
                              public virtual QuadElementBase
{

public:

 /// \short Shorthand for pointer to an argument-free void member 
 /// function of the refineable element
 typedef void (RefineableQElement<2>::*VoidMemberFctPt)();
 
 /// Constructor: Pass refinement level (default 0 = root)
 RefineableQElement() : RefineableElement()
  {
#ifdef LEAK_CHECK
   LeakCheckNames::RefineableQElement<2>_build+=1;
#endif
  } 
 

 /// Broken copy constructor
 RefineableQElement(const RefineableQElement<2>& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableQElement<2>");
  } 
 
 /// Broken assignment operator
//Commented out broken assignment operator because this can lead to a conflict warning
//when used in the virtual inheritence hierarchy. Essentially the compiler doesn't
//realise that two separate implementations of the broken function are the same and so,
//quite rightly, it shouts.
 /*void operator=(const RefineableQElement<2>&) 
  {
   BrokenCopy::broken_assign("RefineableQElement<2>");
   }*/

 /// Destructor
 virtual ~RefineableQElement()
  {
#ifdef LEAK_CHECK
   LeakCheckNames::RefineableQElement<2>_build-=1;
#endif
  } 
    
 /// A refineable quad element has four sons
 unsigned required_nsons() const {return 4;}
 
 /// \short If a neighbouring element has already created a node at
 /// a position corresponding to the local fractional position within the
 /// present element, s_fraction, return
 /// a pointer to that node. If not, return NULL (0).
 /// If the node is on a periodic boundary the flag is_periodic is true,
 /// otherwise it will be false.
 virtual Node* node_created_by_neighbour(const Vector<double> &s_fraction,
                                         bool &is_periodic);
 
 /// \short If a son of a neighbouring element has already created a node at
 /// a position corresponding to the local fractional position within the
 /// present element, s_fraction, return
 /// a pointer to that node. If not, return NULL (0).
 /// If the node is on a periodic boundary the flag is_periodic is true,
 /// otherwise it will be false.
 virtual Node* node_created_by_son_of_neighbour(const Vector<double> &s_fraction,
                                                bool &is_periodic)
  {
   // It is impossible for this situation to arise in meshes
   // containing elements of uniform p-order. This is here so
   // that it can be overloaded for p-refineable elements.
   return 0;
  }

 /// \short Build the element, i.e. give it nodal positions, apply BCs, etc. 
 /// Pointers to any new nodes will be returned in new_node_pt. If 
 /// it is open, the positions of the new 
 /// nodes will be written to the file stream new_nodes_file
 virtual void build(Mesh*& mesh_pt, Vector<Node*>& new_node_pt,
                    bool& was_already_built,
                    std::ofstream &new_nodes_file);
    
 /// \short Check the integrity of the element: ensure that the position and
 /// values are continuous across the element edges
 void check_integrity(double& max_error);
 
 ///  Print corner nodes, use colour
 void output_corners(std::ostream& outfile, const std::string& colour) const;

 /// Pointer to quadtree representation of this element
 QuadTree* quadtree_pt() 
  {return dynamic_cast<QuadTree*>(Tree_pt);}

 /// Pointer to quadtree representation of this element
 QuadTree* quadtree_pt() const  {return dynamic_cast<QuadTree*>(Tree_pt);}

 /// \short Markup all hanging nodes & document the results in
 /// the output streams contained in the vector output_stream, if they
 /// are open.
 void setup_hanging_nodes(Vector<std::ofstream*> &output_stream);
 
 /// \short Perform additional hanging node procedures for variables
 /// that are not interpolated by all nodes (e.g. lower order interpolations
 /// as for the pressure in Taylor Hood). 
 virtual void further_setup_hanging_nodes()=0;
 
  protected:
 
 /// \short Coincidence between son nodal points and father boundaries:  
 /// Father_bound[node_1d](jnod_son,son_type)={SW/SE/NW/NE/S/E/N/W/OMEGA}
 static std::map<unsigned, DenseMatrix<int> > Father_bound;

 /// \short Setup static matrix for coincidence between son 
 /// nodal points and father boundaries
 void setup_father_bounds();
  
 /// Determine Vector of boundary conditions along edge (N/S/W/E)
 void get_edge_bcs(const int& edge, Vector<int>& bound_cons) const;
 
  public:
 /// \short Determine set of (mesh) boundaries that the 
 /// element edge/vertex lives on
 void get_boundaries(const int& edge, std::set<unsigned>& boundaries) const;
  
 /// \short Determine Vector of boundary conditions along edge 
 /// (or on vertex) bound (S/W/N/E/SW/SE/NW/NE): For value ival 
 /// on this boundary, bound_cons[ival]=1 if pinned and 0 if free.
 void get_bcs(int bound, Vector<int>& bound_cons) const;
 
 /// \short Return the value of the intrinsic boundary coordinate
 /// interpolated along the edge (S/W/N/E)
 void interpolated_zeta_on_edge(const unsigned &boundary,
                                const int &edge,
                                const Vector<double> &s,
                                Vector<double> &zeta);

  protected:

 /// \short Internal helper function that is used to construct the
 /// hanging node schemes for the value_id-th interpolated value
 void setup_hang_for_value(const int &value_id);
 
 /// \short Internal helper function that is used to construct the
 /// hanging node schemes for the positions.
 virtual void quad_hang_helper(const int &value_id, const int &my_edge,
                               std::ofstream &output_hangfile);

};



//========================================================================
/// Refineable version of Solid quad elements
//========================================================================
template<>
class RefineableSolidQElement<2> : public virtual RefineableQElement<2>, 
                                   public virtual RefineableSolidElement,
                                   public virtual QSolidElementBase
{
  public:
 
 /// Constructor, just call the constructor of the RefineableQElement<2>
 RefineableSolidQElement() :
  RefineableQElement<2>(), RefineableSolidElement()
  {}


 /// Broken copy constructor
 RefineableSolidQElement(const RefineableSolidQElement<2>& dummy) 
  { 
   BrokenCopy::broken_copy("RefineableSolidQElement<2>");
  } 
 
 /// Broken assignment operator
 /*void operator=(const RefineableSolidQElement<2>&) 
  {
   BrokenCopy::broken_assign("RefineableSolidQElement<2>");
   }*/

 /// Virtual Destructor
 virtual ~RefineableSolidQElement() {}


 /// \short Final over-ride: Use version in QSolidElementBase
 void set_macro_elem_pt(MacroElement* macro_elem_pt)
  {
   QSolidElementBase::set_macro_elem_pt(macro_elem_pt);
  }
 
 /// \short Final over-ride: Use version in QSolidElementBase
 void set_macro_elem_pt(MacroElement* macro_elem_pt,
                        MacroElement* undeformed_macro_elem_pt)
  {
   QSolidElementBase::set_macro_elem_pt(macro_elem_pt,
                                        undeformed_macro_elem_pt);
  }
 
 /// \short Use the generic finite difference routine defined in 
 /// RefineableSolidElement to calculate the Jacobian matrix
 void get_jacobian(Vector<double> &residuals, 
                   DenseMatrix<double> &jacobian) 
  {RefineableSolidElement::get_jacobian(residuals,jacobian);}

 /// \short Determine vector of solid (positional) boundary conditions 
 /// along edge (N/S/W/E) [Pressure does not have to be included
 /// since it can't be subjected to bc at more than one node anyway]
 void get_edge_solid_bcs(const int& edge, Vector<int>& solid_bound_cons) const;
 
 /// \short Determine vector of solid (positional) boundary conditions 
 /// along edge (or on vertex) bound (S/W/N/E/SW/SE/NW/NE): For direction i,
 /// solid_bound_cons[i]=1 if displacement in this coordinate direction
 /// is pinned and 0 if it's free.
 void get_solid_bcs(int bound, Vector<int>& solid_bound_cons) const;


 /// \short Build the element, i.e. give it nodal positions, apply BCs, etc. 
 /// Incl. documention into new_nodes_file
 // NOTE: FOR SOME REASON THIS NEEDS TO LIVE IN *.H TO WORK ON INTEL
 void build(Mesh*& mesh_pt, Vector<Node*> &new_node_pt, 
            bool& was_already_built,
            std::ofstream &new_nodes_file)
  {
   using namespace QuadTreeNames;

   //Call the standard (non-elastic) build function
   RefineableQElement<2>::build(mesh_pt,new_node_pt,was_already_built,
                                new_nodes_file);

   // Are we done?
   if (was_already_built) return;

   //Now need to loop over the nodes again and set solid variables
 
   // What type of son am I? Ask my quadtree representation...
   int son_type=Tree_pt->son_type();
 
   // Which element (!) is my father? (We must have a father
   // since was_already_built is false...)
   RefineableSolidQElement<2>* father_el_pt=
    dynamic_cast<RefineableSolidQElement<2>*>
    (Tree_pt->father_pt()->object_pt());
 

#ifdef PARANOID
   // Currently we can't handle the case of generalised coordinates
   // since we haven't established how they should be interpolated
   // Buffer this case:
   if (static_cast<SolidNode*>(father_el_pt->node_pt(0))
       ->nlagrangian_type()!=1)
    {
     throw OomphLibError(
      "We can't handle generalised nodal positions (yet).\n",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
 
   Vector<double> s_lo(2);
   Vector<double> s_hi(2);
   Vector<double> s(2);
   Vector<double> xi(2);
   Vector<double> xi_fe(2);
   Vector<double> x(2);
   Vector<double> x_fe(2);
 
   // Setup vertex coordinates in father element:
   //--------------------------------------------
   switch(son_type)
    {
    case SW:
     s_lo[0]=-1.0;
     s_hi[0]= 0.0;
     s_lo[1]=-1.0;
     s_hi[1]= 0.0;
     break;
   
    case SE:
     s_lo[0]= 0.0;
     s_hi[0]= 1.0;
     s_lo[1]=-1.0;
     s_hi[1]= 0.0;
     break;
     
    case NE:
     s_lo[0]= 0.0;
     s_hi[0]= 1.0;
     s_lo[1]= 0.0;
     s_hi[1]= 1.0;
     break;

    case NW:
     s_lo[0]=-1.0;
     s_hi[0]= 0.0;
     s_lo[1]= 0.0;
     s_hi[1]= 1.0;
     break;
    }

   //Pass the undeformed macro element onto the son
   // hierher why can I read this? 
   if(father_el_pt->Undeformed_macro_elem_pt!=0)
    {
     Undeformed_macro_elem_pt = father_el_pt->Undeformed_macro_elem_pt;
     for(unsigned i=0;i<2;i++)
      {
       s_macro_ll(i)=      father_el_pt->s_macro_ll(i)+
        0.5*(s_lo[i]+1.0)*(father_el_pt->s_macro_ur(i)-
                           father_el_pt->s_macro_ll(i));
       s_macro_ur(i)=      father_el_pt->s_macro_ll(i)+
        0.5*(s_hi[i]+1.0)*(father_el_pt->s_macro_ur(i)-
                           father_el_pt->s_macro_ll(i));
      }
    }

   //Local node number
   unsigned n=0;

   //Find number of 1D nodes in element
   unsigned n_p = nnode_1d();
   
   // Loop over nodes in element
   for (unsigned i0=0;i0<n_p;i0++)
    {
     // Local coordinate in father element
     s[0]=s_lo[0] + (s_hi[0]-s_lo[0])*double(i0)/double(n_p-1);
     
     for (unsigned i1=0;i1<n_p;i1++)
      {
       // Local coordinate in father element
       s[1]=s_lo[1] + (s_hi[1]-s_lo[1])*double(i1)/double(n_p-1);
       
       // Local node number
       n = i0 + n_p*i1;
       
       // Get position from father element -- this uses the macro
       // element representation(s) if appropriate. If the node
       // turns out to be a hanging node later on, then
       // its position gets adjusted in line with its
       // hanging node interpolation.
       father_el_pt->get_x_and_xi(s,x_fe,x,xi_fe,xi);
       
       //Cast the node to an Solid node
       SolidNode* elastic_node_pt = 
        static_cast<SolidNode*>(node_pt(n));

       for (unsigned i=0;i<2;i++)
        {
         // x_fe is the FE representation -- this is all we can
         // work with in a solid mechanics problem. If you wish
         // to reposition nodes on curvilinear boundaries of 
         // a domain to their exact positions on those boundaries 
         // you'll have to do this yourself! [Note: We used to 
         // use the macro-element-based representation 
         // to assign the position of pinned nodes but this is not always
         // correct since pinning doesn't mean "pin in place" or
         // "pin to the curvilinear boundary". For instance, we could impose
         // the boundary displacement manually.
         // x_fe is the FE representation
         elastic_node_pt->x(i) = x_fe[i];

         // Lagrangian coordinates can come from undeformed macro element
         
         if (Use_undeformed_macro_element_for_new_lagrangian_coords)
          {
           elastic_node_pt->xi(i) = xi[i];
          }
         else
          {
           elastic_node_pt->xi(i) = xi_fe[i];
          }
        }
       
       
       // Are there any history values to be dealt with?
       TimeStepper* time_stepper_pt=father_el_pt->
        node_pt(0)->time_stepper_pt();
       
       // Number of history values (incl. present)
       unsigned ntstorage=time_stepper_pt->ntstorage();
       if (ntstorage!=1)
        {
         // Loop over # of history values (excluding present which has been
         // done above)
         for (unsigned t=1;t<ntstorage;t++)
          {
           // History values can (and in the case of Newmark timestepping,
           // the scheme most likely to be used for Solid computations, do)
           // include non-positional values, e.g. velocities and accelerations.
           
           // Set previous positions of the new node
           for(unsigned i=0;i<2;i++)
            {
             elastic_node_pt->x(t,i) = father_el_pt->interpolated_x(t,s,i);
            }
          }
        }       
       
      } // End of vertical loop over nodes in element
     
    } // End of horizontal loop over nodes in element
   
  }
 

};

}

#endif

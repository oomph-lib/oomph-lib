// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#ifndef OOMPH_REFINEABLE_BRICK_ELEMENT_HEADER
#define OOMPH_REFINEABLE_BRICK_ELEMENT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// ooomph-lib includes
#include "octree.h"
#include "refineable_elements.h"
#include "Qelements.h"

namespace oomph
{
  class Mesh;

  //=======================================================================
  /// Refineable version of QElement<3,NNODE_1D>.
  ///
  /// Refinement is performed by octree procedures. When the element is
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
  /// The class provides a generic RefineableQElement<3>::build() function
  /// which deals with generic
  /// isoparametric QElements in which all values are associated with
  /// nodes. The RefineableQElement<3>::further_build() function provides
  /// an interface for any element-specific non-generic build operations.
  ///
  //=======================================================================
  template<>
  class RefineableQElement<3> :
    public virtual RefineableElement,
    public virtual BrickElementBase
  {
  public:
    /// \short Shorthand for pointer to an argument-free void member
    /// function of the refineable element
    typedef void (RefineableQElement<3>::*VoidMemFctPt)();

    /// Constructor: Pass refinement level (default 0 = root)
    RefineableQElement<3>() : RefineableElement()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQElement<3> _build += 1;
#endif
    }

    /// Broken copy constructor
    RefineableQElement<3>(const RefineableQElement<3> &dummy)
    {
      BrokenCopy::broken_copy("RefineableQElement<3>");
    }

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const RefineableQElement<3>&)
     {
      BrokenCopy::broken_assign("RefineableQElement<3>");
      }*/

    /// Destructor
    virtual ~RefineableQElement<3>()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQElement<3> _build -= 1;
#endif
    }

    /// A refineable brick element has eight sons
    unsigned required_nsons() const
    {
      return 8;
    }

    /// \short If a neighbouring element has already created a node at
    /// a position corresponding to the local fractional position within the
    /// present element, s_fraction, return
    /// a pointer to that node. If not, return NULL (0).
    virtual Node *node_created_by_neighbour(const Vector<double> &s_fraction);

    /// \short If a neighbouring element has already created a node at
    /// a position corresponding to the local fractional position within the
    /// present element, s_fraction, return
    /// a pointer to that node. If not, return NULL (0).
    virtual Node *node_created_by_son_of_neighbour(
      const Vector<double> &s_fraction)
    {
      // It is impossible for this situation to arise in meshes
      // containing elements of uniform p-order. This is here so
      // that it can be overloaded for p-refineable elements.
      return 0;
    }

    /// \short Build the element: i.e. give it nodal positions, apply BCs, etc.
    /// Pointers to any new nodes will be returned in new_node_pt. If
    /// it is open, the positions of the new
    /// nodes will be written to the file stream new_nodes_file
    virtual void build(Mesh *&mesh_pt,
                       Vector<Node *> &new_node_pt,
                       bool &was_already_built,
                       std::ofstream &new_nodes_file);

    /// \short Check the integrity of the element: ensure that the position and
    /// values are continuous across the element faces
    void check_integrity(double &max_error);

    ///  Print corner nodes, use colour
    void output_corners(std::ostream &outfile, const std::string &colour) const;

    /// Pointer to octree representation of this element
    OcTree *octree_pt()
    {
      return dynamic_cast<OcTree *>(Tree_pt);
    }

    /// Pointer to octree representation of this element
    OcTree *octree_pt() const
    {
      return dynamic_cast<OcTree *>(Tree_pt);
    }

    /// \short Markup all hanging nodes & document the results in
    /// the output streams contained in the vector output_stream, if they
    /// are open.
    void setup_hanging_nodes(Vector<std::ofstream *> &output_stream);

    /// \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes (e.g. lower order interpolations
    /// as for the pressure in Taylor Hood).
    virtual void further_setup_hanging_nodes() = 0;

  protected:
    /// \short Coincidence between son nodal points and father boundaries:
    /// Father_bound[nnode_1d](nnode_son,son_type)={RU/RF/RD/RB/.../ OMEGA}
    /// so that node nnode_son in element of type son_type lies
    /// on boundary/vertex Father_bound[nnode_1d](nnode_son,son_type) in its
    /// father element. If the node doesn't lie on a boundary
    /// the value is OMEGA.
    static std::map<unsigned, DenseMatrix<int>> Father_bound;

    /// \short Setup static matrix for coincidence between son
    /// nodal points and father boundaries
    void setup_father_bounds();

    /// \short Determine Vector of boundary conditions along the element's
    /// face (R/L/U/D/B/F) -- BC is the least restrictive combination
    /// of all the nodes on this face.
    ///
    /// Usual convention:
    ///   - bound_cons[ival]=0 if value ival on this boundary is free
    ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
    void get_face_bcs(const int &edge, Vector<int> &bound_cons) const;

    /// Given an element edge/vertex, return a Vector which contains
    /// all the (mesh-)boundary numbers that this element edge/vertex
    /// lives on.
    ///
    /// For proper edges, the boundary is the one (if any) that is shared by
    /// both vertex nodes). For vertex nodes, we just return their
    /// boundaries.
    void get_boundaries(const int &edge, std::set<unsigned> &boundaries) const;

    /// \short Determine Vector of boundary conditions along the element's
    /// boundary (or vertex) bound (S/W/N/E/SW/SE/NW/NE).
    ///
    /// This function assumes that the same boundary condition is applied
    /// along the entire area of an element's face (of course, the
    /// vertices combine the boundary conditions of their two adjacent edges
    /// in the most restrictive combination. Hence, if we're at a vertex,
    /// we apply the most restrictive boundary condition of the
    /// two adjacent edges. If we're on an edge (in its proper interior),
    /// we apply the least restrictive boundary condition of all nodes
    /// along the edge.
    ///
    /// Usual convention:
    ///   - bound_cons[ival]=0 if value ival on this boundary is free
    ///   - bound_cons[ival]=1 if value ival on this boundary is pinned
    void get_bcs(int bound, Vector<int> &bound_cons) const;

    /// \short Return the value of the intrinsic boundary coordinate
    /// interpolated along the face
    void interpolated_zeta_on_face(const unsigned &boundary,
                                   const int &face,
                                   const Vector<double> &s,
                                   Vector<double> &zeta);

    /// \short Internal helper function that is used to construct the
    /// hanging node schemes for the value_id-th interpolated value
    void setup_hang_for_value(const int &value_id);

    /// \short Internal helper function that is used to construct the
    /// hanging node schemes for the positions.
    virtual void oc_hang_helper(const int &value_id,
                                const int &my_edge,
                                std::ofstream &output_hangfile);
  };

  //========================================================================
  /// Refineable version of Solid brick elements
  //========================================================================
  template<>
  class RefineableSolidQElement<3> :
    public virtual RefineableQElement<3>,
    public virtual RefineableSolidElement,
    public virtual QSolidElementBase
  {
  public:
    /// Constructor, just call the constructor of the RefineableQElement<2>
    RefineableSolidQElement() :
      RefineableQElement<3>(), RefineableSolidElement()
    {
    }

    /// Broken copy constructor
    RefineableSolidQElement(const RefineableSolidQElement<3> &dummy)
    {
      BrokenCopy::broken_copy("RefineableSolidQElement<3>");
    }

    /// Broken assignment operator
    /*void operator=(const RefineableSolidQElement<3>&)
     {
      BrokenCopy::broken_assign("RefineableSolidQElement<3>");
      }*/

    /// Virtual Destructor
    virtual ~RefineableSolidQElement() {}

    /// \short Final over-ride: Use version in QSolidElementBase
    void set_macro_elem_pt(MacroElement *macro_elem_pt)
    {
      QSolidElementBase::set_macro_elem_pt(macro_elem_pt);
    }

    /// \short Final over-ride: Use version in QSolidElementBase
    void set_macro_elem_pt(MacroElement *macro_elem_pt,
                           MacroElement *undeformed_macro_elem_pt)
    {
      QSolidElementBase::set_macro_elem_pt(macro_elem_pt,
                                           undeformed_macro_elem_pt);
    }

    /// \short Use the generic finite difference routine defined in
    /// RefineableSolidElement to calculate the Jacobian matrix
    void get_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian)
    {
      RefineableSolidElement::get_jacobian(residuals, jacobian);
    }

    /// \short Determine vector of solid (positional) boundary conditions
    /// along face (R/L/U/D/B/F) [Pressure does not have to be included
    /// since it can't be subjected to bc at more than one node anyway]
    void get_face_solid_bcs(const int &edge,
                            Vector<int> &solid_bound_cons) const;

    /// \short Determine vector of solid (positional) boundary conditions
    /// along edge (or on vertex) bound (S/W/N/E/SW/SE/NW/NE): For direction i,
    /// solid_bound_cons[i]=1 if displacement in this coordinate direction
    /// is pinned and 0 if it's free.
    void get_solid_bcs(int bound, Vector<int> &solid_bound_cons) const;

    /// \short Build the element, i.e. give it nodal positions, apply BCs, etc.
    /// Incl. documention into new_nodes_file
    // NOTE: FOR SOME REASON THIS NEEDS TO LIVE IN *.H TO WORK ON INTEL
    void build(Mesh *&mesh_pt,
               Vector<Node *> &new_node_pt,
               bool &was_already_built,
               std::ofstream &new_nodes_file)
    {
      using namespace OcTreeNames;

      // Call the standard (non-elastic) build function
      RefineableQElement<3>::build(
        mesh_pt, new_node_pt, was_already_built, new_nodes_file);

      // Are we done?
      if (was_already_built) return;

      // Now need to loop over the nodes again and set solid variables

      // What type of son am I? Ask my quadtree representation...
      int son_type = octree_pt()->son_type();

      // Which element (!) is my father? (We must have a father
      // since was_already_built is false...)
      RefineableSolidQElement<3> *father_el_pt =
        dynamic_cast<RefineableSolidQElement<3> *>(
          Tree_pt->father_pt()->object_pt());

#ifdef PARANOID
      // Currently we can't handle the case of generalised coordinates
      // since we haven't established how they should be interpolated
      // Buffer this case:
      if (static_cast<SolidNode *>(father_el_pt->node_pt(0))
            ->nlagrangian_type() != 1)
      {
        throw OomphLibError(
          "We can't handle generalised nodal positions (yet).\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Now get coordinates and stuff
      Vector<int> s_lo(3);
      Vector<int> s_hi(3);
      Vector<double> s(3);
      Vector<double> xi(3);
      Vector<double> xi_fe(3);
      Vector<double> x(3);
      Vector<double> x_fe(3);

      // Get the number of 1d nodes
      unsigned n_p = nnode_1d();

      // Setup vertex coordinates in father element:
      //--------------------------------------------

      // find the s_lo coordinates
      s_lo = octree_pt()->Direction_to_vector[son_type];

      // just scale them, because the Direction_to_vector
      // doesn't really gives s_lo;
      for (int i = 0; i < 3; i++)
      {
        s_lo[i] = (s_lo[i] + 1) / 2 - 1;
      }

      // setup s_hi (Actually s_hi[i]=s_lo[i]+1)
      for (int i = 0; i < 3; i++)
      {
        s_hi[i] = s_lo[i] + 1;
      }

      // Pass Undeformed macro element pointer on to sons and
      // set coordinates in macro element
      if (father_el_pt->Undeformed_macro_elem_pt != 0)
      {
        Undeformed_macro_elem_pt = father_el_pt->Undeformed_macro_elem_pt;
        for (unsigned i = 0; i < 3; i++)
        {
          s_macro_ll(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_lo[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
          s_macro_ur(i) =
            father_el_pt->s_macro_ll(i) +
            0.5 * (s_hi[i] + 1.0) *
              (father_el_pt->s_macro_ur(i) - father_el_pt->s_macro_ll(i));
        }
      }

      unsigned jnod = 0;
      Vector<double> x_small(3);
      Vector<double> x_large(3);

      Vector<double> s_fraction(3);

      // Loop over nodes in element
      for (unsigned i0 = 0; i0 < n_p; i0++)
      {
        // Get the fractional position of the node in the direction of s[0]
        s_fraction[0] = local_one_d_fraction_of_node(i0, 0);
        // Local coordinate in father element
        s[0] = s_lo[0] + (s_hi[0] - s_lo[0]) * s_fraction[0];

        for (unsigned i1 = 0; i1 < n_p; i1++)
        {
          // Get the fractional position of the node in the direction of s[1]
          s_fraction[1] = local_one_d_fraction_of_node(i1, 1);
          // Local coordinate in father element
          s[1] = s_lo[1] + (s_hi[1] - s_lo[1]) * s_fraction[1];

          for (unsigned i2 = 0; i2 < n_p; i2++)
          {
            // Get the fractional position of the node in the direction of s[2]
            s_fraction[2] = local_one_d_fraction_of_node(i2, 2);
            // Local coordinate in father element
            s[2] = s_lo[2] + (s_hi[2] - s_lo[2]) * s_fraction[2];

            // Local node number
            jnod = i0 + n_p * i1 + n_p * n_p * i2;

            // Get position from father element -- this uses the macro
            // element representation(s) if appropriate. If the node
            // turns out to be a hanging node later on, then
            // its position gets adjusted in line with its
            // hanging node interpolation.
            father_el_pt->get_x_and_xi(s, x_fe, x, xi_fe, xi);

            // Cast the node to an Solid node
            SolidNode *elastic_node_pt =
              static_cast<SolidNode *>(node_pt(jnod));

            for (unsigned i = 0; i < 3; i++)
            {
              // x_fe is the FE representation -- this is all we can
              // work with in a solid mechanics problem. If you wish
              // to reposition nodes on curvilinear boundaries of
              // a domain to their exact positions on those boundaries
              // you'll have to do this yourself! [Note: We used to
              // use the macro-element-based representation
              // to assign the position of pinned nodes but this is not always
              // correct since pinning doesn't mean "pin in place" or
              // "pin to the curvilinear boundary". For instance, we could
              // impose the boundary displacement manually.
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
            TimeStepper *time_stepper_pt =
              father_el_pt->node_pt(0)->time_stepper_pt();

            // Number of history values (incl. present)
            unsigned ntstorage = time_stepper_pt->ntstorage();
            if (ntstorage != 1)
            {
              // Loop over # of history values (excluding present which has been
              // done above)
              for (unsigned t = 1; t < ntstorage; t++)
              {
                // History values can (and in the case of Newmark timestepping,
                // the scheme most likely to be used for Solid computations, do)
                // include non-positional values, e.g. velocities and accels.

                // Set previous positions of the new node
                for (unsigned i = 0; i < 3; i++)
                {
                  elastic_node_pt->x(t, i) =
                    father_el_pt->interpolated_x(t, s, i);
                }
              }
            }
          } // End of s2 loop
        } // End of vertical loop over nodes in element

      } // End of horizontal loop over nodes in element
    }
  };

} // namespace oomph

#endif

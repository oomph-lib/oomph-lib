// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_REFINEABLE_QUAD_SPECTRAL_ELEMENT_HEADER
#define OOMPH_REFINEABLE_QUAD_SPECTRAL_ELEMENT_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib headers
#include "refineable_quad_element.h"

namespace oomph
{
  //=======================================================================
  /// Refineable version of QuadElements that add functionality for spectral
  /// Elements.
  //=======================================================================
  template<>
  class RefineableQSpectralElement<2> : public virtual RefineableQElement<2>
  {
  public:
    /// Constructor
    RefineableQSpectralElement() : RefineableElement()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQSpectralElement<2> _build += 1;
#endif
    }


    /// Broken copy constructor
    RefineableQSpectralElement(const RefineableQSpectralElement<2>& dummy) =
      delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const RefineableQSpectralElement<2>&) = delete;*/

    /// Destructor
    virtual ~RefineableQSpectralElement()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQSpectralElement<2> _build -= 1;
#endif
    }

    /// The only thing to add is rebuild from sons
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      // The timestepper should be the same for all nodes and node 0 should
      // never be deleted.
      if (this->node_pt(0) == 0)
      {
        throw OomphLibError("The Corner node (0) does not exist",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      TimeStepper* time_stepper_pt = this->node_pt(0)->time_stepper_pt();
      unsigned ntstorage = time_stepper_pt->ntstorage();

      unsigned jnod = 0;
      Vector<double> s_fraction(2), s(2);
      // Loop over the nodes in the element
      unsigned n_p = this->nnode_1d();
      for (unsigned i0 = 0; i0 < n_p; i0++)
      {
        // Get the fractional position of the node
        s_fraction[0] = this->local_one_d_fraction_of_node(i0, 0);
        // Local coordinate
        s[0] = -1.0 + 2.0 * s_fraction[0];

        for (unsigned i1 = 0; i1 < n_p; i1++)
        {
          // Get the fractional position of the node in the direction of s[1]
          s_fraction[1] = this->local_one_d_fraction_of_node(i1, 1);
          // Local coordinate in father element
          s[1] = -1.0 + 2.0 * s_fraction[1];

          // Set the local node number
          jnod = i0 + n_p * i1;

          // If the node has not been built
          if (this->node_pt(jnod) == 0)
          {
            // Has the node been created by one of its neighbours
            bool is_periodic = false;
            Node* created_node_pt =
              this->node_created_by_neighbour(s_fraction, is_periodic);

            // If it has set the pointer
            if (created_node_pt != 0)
            {
              // If the node is periodic
              if (is_periodic)
              {
                throw OomphLibError("Cannot handle periodic nodes in "
                                    "refineable spectral elements yet",
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
              // Non-periodic case, just set the pointer
              else
              {
                this->node_pt(jnod) = created_node_pt;
              }
            }
            // Otherwise, we need to build it
            else
            {
              // First, find the son element in which the node should live

              // Find coordinates in the sons
              Vector<double> s_in_son(2);
              using namespace QuadTreeNames;
              int son = -10;
              // If negative on the west side
              if (s_fraction[0] < 0.5)
              {
                // On the south side
                if (s_fraction[1] < 0.5)
                {
                  // It's the southwest son
                  son = SW;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                }
                // On the north side
                else
                {
                  // It's the northwest son
                  son = NW;
                  s_in_son[0] = -1.0 + 4.0 * s_fraction[0];
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                }
              }
              else
              {
                // On the south side
                if (s_fraction[1] < 0.5)
                {
                  // It's the southeast son
                  son = SE;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * s_fraction[1];
                }
                // On the north side
                else
                {
                  // It's the northeast son
                  son = NE;
                  s_in_son[0] = -1.0 + 4.0 * (s_fraction[0] - 0.5);
                  s_in_son[1] = -1.0 + 4.0 * (s_fraction[1] - 0.5);
                }
              }

              // Get the pointer to the son element in which the new node
              // would live
              RefineableQSpectralElement<2>* son_el_pt =
                dynamic_cast<RefineableQSpectralElement<2>*>(
                  this->tree_pt()->son_pt(son)->object_pt());

              // If we are rebuilding, then worry about the boundary conditions
              // Find the boundary of the node
              // Initially none
              int boundary = Tree::OMEGA;
              // If we are on the western boundary
              if (i0 == 0)
              {
                boundary = W;
              }
              // If we are on the eastern boundary
              else if (i0 == n_p - 1)
              {
                boundary = E;
              }

              // If we are on the southern boundary
              if (i1 == 0)
              {
                // If we already have already set the boundary, we're on a
                // corner
                switch (boundary)
                {
                  case W:
                    boundary = SW;
                    break;
                  case E:
                    boundary = SE;
                    break;
                    // Boundary not set
                  default:
                    boundary = S;
                    break;
                }
              }
              // If we are the northern bounadry
              else if (i1 == n_p - 1)
              {
                // If we already have a boundary
                switch (boundary)
                {
                  case W:
                    boundary = NW;
                    break;
                  case E:
                    boundary = NE;
                    break;
                  default:
                    boundary = N;
                    break;
                }
              }

              // set of boundaries that this edge in the son lives on
              std::set<unsigned> boundaries;

              // Now get the boundary conditions from the son
              // The boundaries will be common to the son because there can be
              // no rotations here
              if (boundary != Tree::OMEGA)
              {
                son_el_pt->get_boundaries(boundary, boundaries);
              }

              // If the node lives on a boundary:
              // Construct a boundary node,
              // Get boundary conditions and
              // update all lookup schemes
              if (boundaries.size() > 0)
              {
                // Construct the new node
                this->node_pt(jnod) =
                  this->construct_boundary_node(jnod, time_stepper_pt);

                // Get the boundary conditions from the son
                Vector<int> bound_cons(ncont_interpolated_values());
                son_el_pt->get_bcs(boundary, bound_cons);

                // Loop over the values and pin if necessary
                unsigned nval = this->node_pt(jnod)->nvalue();
                for (unsigned k = 0; k < nval; k++)
                {
                  if (bound_cons[k])
                  {
                    this->node_pt(jnod)->pin(k);
                  }
                }

                // Solid node? If so, deal with the positional boundary
                // conditions:
                SolidNode* solid_node_pt =
                  dynamic_cast<SolidNode*>(this->node_pt(jnod));
                if (solid_node_pt != 0)
                {
                  // Get the positional boundary conditions from the father:
                  unsigned n_dim = this->node_pt(jnod)->ndim();
                  Vector<int> solid_bound_cons(n_dim);
                  RefineableSolidQElement<2>* son_solid_el_pt =
                    dynamic_cast<RefineableSolidQElement<2>*>(son_el_pt);
#ifdef PARANOID
                  if (son_solid_el_pt == 0)
                  {
                    std::string error_message =
                      "We have a SolidNode outside a refineable SolidElement\n";
                    error_message +=
                      "during mesh refinement -- this doesn't make sense\n";

                    throw OomphLibError(error_message,
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
#endif
                  son_solid_el_pt->get_solid_bcs(boundary, solid_bound_cons);

                  // Loop over the positions and pin, if necessary
                  for (unsigned k = 0; k < n_dim; k++)
                  {
                    if (solid_bound_cons[k])
                    {
                      solid_node_pt->pin_position(k);
                    }
                  }
                }

                // Next we update the boundary look-up schemes
                // Loop over the boundaries stored in the set
                for (std::set<unsigned>::iterator it = boundaries.begin();
                     it != boundaries.end();
                     ++it)
                {
                  // Add the node to the boundary
                  mesh_pt->add_boundary_node(*it, this->node_pt(jnod));

                  // If we have set an intrinsic coordinate on this
                  // mesh boundary then it must also be interpolated on
                  // the new node
                  // Now interpolate the intrinsic boundary coordinate
                  if (mesh_pt->boundary_coordinate_exists(*it) == true)
                  {
                    Vector<double> zeta(1);
                    son_el_pt->interpolated_zeta_on_edge(
                      *it, boundary, s_in_son, zeta);

                    this->node_pt(jnod)->set_coordinates_on_boundary(*it, zeta);
                  }
                }
              }
              // Otherwise the node is not on a Mesh boundary
              // and we create a normal "bulk" node
              else
              {
                // Construct the new node
                this->node_pt(jnod) =
                  this->construct_node(jnod, time_stepper_pt);
              }

              // Now we set the position and values at the newly created node

              // In the first instance use macro element or FE representation
              // to create past and present nodal positions.
              // (THIS STEP SHOULD NOT BE SKIPPED FOR ALGEBRAIC
              // ELEMENTS AS NOT ALL OF THEM NECESSARILY IMPLEMENT
              // NONTRIVIAL NODE UPDATE FUNCTIONS. CALLING
              // THE NODE UPDATE FOR SUCH ELEMENTS/NODES WILL LEAVE
              // THEIR NODAL POSITIONS WHERE THEY WERE (THIS IS APPROPRIATE
              // ONCE THEY HAVE BEEN GIVEN POSITIONS) BUT WILL
              // NOT ASSIGN SENSIBLE INITIAL POSITONS!

              // Loop over # of history values
              // Loop over # of history values
              for (unsigned t = 0; t < ntstorage; t++)
              {
                using namespace QuadTreeNames;
                // Get the position from the son
                Vector<double> x_prev(2);

                // Now let's fill in the value
                son_el_pt->get_x(t, s_in_son, x_prev);
                for (unsigned i = 0; i < 2; i++)
                {
                  this->node_pt(jnod)->x(t, i) = x_prev[i];
                }
              }

              // Now set up the values
              // Loop over all history values
              for (unsigned t = 0; t < ntstorage; t++)
              {
                // Get values from father element
                // Note: get_interpolated_values() sets Vector size itself.
                Vector<double> prev_values;
                son_el_pt->get_interpolated_values(t, s_in_son, prev_values);

                // Initialise the values at the new node
                for (unsigned k = 0; k < this->node_pt(jnod)->nvalue(); k++)
                {
                  this->node_pt(jnod)->set_value(t, k, prev_values[k]);
                }
              }

              // Add the node to the mesh
              mesh_pt->add_node_pt(this->node_pt(jnod));

            } // End of the case when we build the node ourselvesx

            // Algebraic stuff here
            // Check whether the element is an algebraic element
            AlgebraicElementBase* alg_el_pt =
              dynamic_cast<AlgebraicElementBase*>(this);

            // If we do have an algebraic element
            if (alg_el_pt != 0)
            {
              std::string error_message =
                "Have not implemented rebuilding from sons for";
              error_message += "Algebraic Spectral elements yet\n";

              throw OomphLibError(
                error_message,
                "RefineableQSpectralElement::rebuild_from_sons()",
                OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
    }


    /// Overload the nodes built function
    virtual bool nodes_built()
    {
      unsigned n_node = this->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        if (node_pt(n) == 0)
        {
          return false;
        }
      }
      // If we get to here, OK
      return true;
    }
  };

} // namespace oomph

#endif

// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_BACKWARD_STEP_MESH_TEMPLATE_CC
#define OOMPH_BACKWARD_STEP_MESH_TEMPLATE_CC

#include <set>
#include <map>
#include "backward_step_mesh.template.h"

namespace oomph
{
  //=================================================================
  /// Actual build function. Pass overall number of elements in the
  /// horizontal and vertical directions, nx and ny, and the corresponding
  /// dimensions, lx and ly. nx_cut_out and ny_cut_out elements
  /// are cut out from the lower right corner to create the
  /// (reversed) backward step geometry. Timestepper defaults
  /// to Steady.
  //=================================================================
  template<class ELEMENT>
  void BackwardStepQuadMesh<ELEMENT>::build_mesh(const unsigned& nx,
                                                 const unsigned& ny,
                                                 const unsigned& nx_cut_out,
                                                 const unsigned& ny_cut_out,
                                                 const double& lx,
                                                 const double& ly)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // By default nobody's claiming any nodes
    std::map<Node*, bool> keep;

    // Get elements outside lower right block
    Vector<FiniteElement*> el_retained_pt;
    Vector<FiniteElement*> el_killed_pt;
    for (unsigned i = 0; i < nx; i++)
    {
      for (unsigned j = 0; j < ny; j++)
      {
        FiniteElement* el_pt = this->finite_element_pt(i + nx * j);
        if ((i > (nx_cut_out - 1)) && (j < ny_cut_out))
        {
          el_killed_pt.push_back(el_pt);
        }
        else
        {
          el_retained_pt.push_back(el_pt);
          unsigned nnod_el = el_pt->nnode();
          for (unsigned jj = 0; jj < nnod_el; jj++)
          {
            keep[el_pt->node_pt(jj)] = true;
          }
        }
      }
    }


    // By default nobody's claiming the nodes; also store old
    // boundary ids
    std::map<Node*, std::set<unsigned>> boundaries;
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      std::set<unsigned>* aux_pt = 0;
      this->node_pt(j)->get_boundaries_pt(aux_pt);
      if (aux_pt != 0)
      {
        boundaries[this->node_pt(j)] = (*aux_pt);
      }
    }

    // Remove information about boundary nodes
    this->remove_boundary_nodes();

    // Reset number of boundaries
    this->set_nboundary(6);

    // Kill superfluous nodes
    Vector<Node*> node_backup_pt(this->Node_pt);
    this->Node_pt.clear();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = node_backup_pt[j];
      if (keep[nod_pt])
      {
        this->Node_pt.push_back(nod_pt);
        std::set<unsigned> aux = boundaries[nod_pt];
        for (std::set<unsigned>::iterator it = boundaries[nod_pt].begin();
             it != boundaries[nod_pt].end();
             it++)
        {
          unsigned b = (*it);
          if (b > 0) b += 2;
          this->add_boundary_node(b, nod_pt);
        }
      }
      else
      {
        delete node_backup_pt[j];
      }
    }

    // Add nodes to new boundary 1
    unsigned i = nx_cut_out - 1;
    for (unsigned j = 0; j < ny_cut_out; j++)
    {
      FiniteElement* el_pt = this->finite_element_pt(i + nx * j);
      unsigned nnod_1d = el_pt->nnode_1d();
      for (unsigned jj = 0; jj < nnod_1d; jj++)
      {
        unsigned jnod = (nnod_1d - 1) + jj * nnod_1d;
        if (!(el_pt->node_pt(jnod)->is_on_boundary()))
        {
          Node* nod_pt = el_pt->node_pt(jnod);
          this->convert_to_boundary_node(nod_pt);
        }
        this->add_boundary_node(1, el_pt->node_pt(jnod));
      }
    }

    // Add nodes to new boundary 2
    unsigned j = ny_cut_out;
    for (unsigned i = nx_cut_out; i < nx; i++)
    {
      FiniteElement* el_pt = this->finite_element_pt(i + nx * j);
      unsigned nnod_1d = el_pt->nnode_1d();
      for (unsigned jj = 0; jj < nnod_1d; jj++)
      {
        unsigned jnod = jj;
        if (!(el_pt->node_pt(jnod)->is_on_boundary()))
        {
          Node* nod_pt = el_pt->node_pt(jnod);
          this->convert_to_boundary_node(nod_pt);
        }
        this->add_boundary_node(2, el_pt->node_pt(jnod));
      }
    }


    // Kill redundant elements
    this->Element_pt.clear();
    unsigned n_retained = el_retained_pt.size();
    for (unsigned e = 0; e < n_retained; e++)
    {
      this->Element_pt.push_back(el_retained_pt[e]);
    }
    unsigned n_killed = el_killed_pt.size();
    for (unsigned e = 0; e < n_killed; e++)
    {
      delete el_killed_pt[e];
    }

    // Re-setup lookup scheme that establishes which elements are located
    // on the mesh boundaries
    this->setup_boundary_element_info();
  }

} // namespace oomph

#endif

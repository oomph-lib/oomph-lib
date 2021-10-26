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
#include "geom_objects.h"
#include "algebraic_elements.h"
#include "refineable_quad_element.h"


namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Base class for Algebraic Elements
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Set up node update info for (newly created) algebraic node: Work out its
  /// node update information by interpolation from its father element,
  /// based on pointer to father element and its local coordinate in
  /// the father element. We're creating the node update info
  /// for update functions that are shared by all nodes in the
  /// father element.
  //========================================================================
  void AlgebraicElementBase::setup_algebraic_node_update(
    Node*& node_pt,
    const Vector<double>& s_father,
    FiniteElement* father_el_pt) const
  {
    // Turn node into algebraic node
    AlgebraicNode* alg_node_pt = dynamic_cast<AlgebraicNode*>(node_pt);

    // Get number of nodes in father element
    unsigned nnod_father = father_el_pt->nnode();

    // Create map that stores the number of times a node update fct id
    // has been encountered
    std::map<int, unsigned> id_count;

    // Loop over all nodes in father element to extract node update ids
    Vector<int> id;
    for (unsigned j = 0; j < nnod_father; j++)
    {
      // Get vector of ids and add them to map
      dynamic_cast<AlgebraicNode*>(father_el_pt->node_pt(j))
        ->node_update_fct_id(id);
      unsigned n_id = id.size();
      for (unsigned i = 0; i < n_id; i++)
      {
        id_count[id[i]]++;
      }
    }

    /// Now loop over the ids and check if they appear in all nodes
    Vector<int> shared_ids;
    typedef std::map<int, unsigned>::iterator IT;
    for (IT it = id_count.begin(); it != id_count.end(); it++)
    {
      if (it->second == nnod_father)
      {
        shared_ids.push_back(it->first);
      }
    }


    // How many update functions we have?
    unsigned n_update_id = shared_ids.size();

    // Loop over all node udate functions -- since it's shared by
    // nodes we may as well read the required data  from the first
    // node in the father.
    AlgebraicNode* father_node_pt =
      dynamic_cast<AlgebraicNode*>(father_el_pt->node_pt(0));
    for (unsigned i = 0; i < n_update_id; i++)
    {
      // Actual id:
      int id = shared_ids[i];

      // Is this a real update fct or the default dummy one?
      if (id >= 0)
      {
        // Get vector of geometric objects involved in specified node update
        // function (create vector by copy operation)
        Vector<GeomObject*> geom_obj_pt(
          father_node_pt->vector_geom_object_pt(id));

        // Loop over reference values and obtain the
        // ones for the current (son) element by interpolation from
        // the father

        // Number of reference values for this udpate function
        unsigned nvalue = father_node_pt->nref_value(id);
        Vector<double> ref_value(nvalue);

        // Set up the shape functions in father element
        Shape psi(nnod_father);

        // Get shape functions in father element
        father_el_pt->shape(s_father, psi);

        // Initialise reference values
        for (unsigned ivalue = 0; ivalue < nvalue; ivalue++)
        {
          ref_value[ivalue] = 0.0;
        }

        // Loop over all nodes in father element for
        // interpolation of nodes -- don't need to
        // worry about hanging nodes here as reference values
        // are only assigned once and then in a consistent way
        for (unsigned j_father = 0; j_father < nnod_father; j_father++)
        {
          // Get reference values at node in father by copy operation
          Vector<double> father_ref_value(
            dynamic_cast<AlgebraicNode*>(father_el_pt->node_pt(j_father))
              ->vector_ref_value(id));

          // Loop over reference values
          for (unsigned ivalue = 0; ivalue < nvalue; ivalue++)
          {
            ref_value[ivalue] += father_ref_value[ivalue] * psi(j_father);
          }
        }

        // Get mesh that implements the update operation
        AlgebraicMesh* mesh_pt = father_node_pt->mesh_pt(id);


        // Setup node update info for node
        alg_node_pt->add_node_update_info(id, // id
                                          mesh_pt, // mesh
                                          geom_obj_pt, // vector of geom objects
                                          ref_value); // vector of ref. values

        // Update the geometric references (e.g. in FSI) if required
        mesh_pt->update_node_update(alg_node_pt);

      } // endif for real vs. default dummy node update fct

    } // end of loop over different update fcts

    // Update the node at its current and previous positions
    bool update_all_time_levels_for_new_node = true;
    alg_node_pt->node_update(update_all_time_levels_for_new_node);
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic nodes
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Assign default value for test if different
  /// node update functions produce the same result.
  //========================================================================
  double AlgebraicNode::Max_allowed_difference_between_node_update_fcts =
    1.0e-10;


  //========================================================================
  /// Default (negative!) remesh fct id for nodes for which no remesh
  /// fct is defined
  //========================================================================
  int AlgebraicNode::Dummy_node_update_fct_id = -100;


  //======================================================================
  /// Set the dummy mesh
  //====================================================================
  DummyAlgebraicMesh AlgebraicNode::Dummy_mesh;

  //========================================================================
  /// Default dummy mesh to point to for nodes for which no remesh
  /// fct is defined
  //========================================================================
  AlgebraicMesh* AlgebraicNode::Dummy_mesh_pt = &AlgebraicNode::Dummy_mesh;


  //========================================================================
  /// Zero-sized default dummy vector of geom objects to point to for nodes
  /// for which no remesh fct is defined
  //========================================================================
  Vector<GeomObject*> AlgebraicNode::Dummy_geom_object_pt;

  //========================================================================
  /// Zero-sized default dummy vector of reference values
  /// to point to for nodes  for which no remesh fct is defined
  //========================================================================
  Vector<double> AlgebraicNode::Dummy_ref_value;


  //========================================================================
  /// Excute the node update function: Update the current (and if
  /// update_all_time_levels_for_new_node==true also the previous)
  /// nodal position. Also update the current nodal values if
  /// an auxiliary update function is defined.
  /// Note: updating of previous positions is only required (and should
  /// only be performed for) newly created AlgebraicNodes
  /// i.e. when this function is called from
  /// AlgebraicElementBase::setup_algebraic_node_update(...). We create
  /// the history of its nodal positions from the time-dependent
  /// version of the specific AlgebraicMesh's algebraic_node_update(...)
  /// function.
  //========================================================================
  void AlgebraicNode::node_update(
    const bool& update_all_time_levels_for_new_node)
  {
    // Number of time levels that need to be updated
    unsigned ntime;
    if (update_all_time_levels_for_new_node)
    {
      // Present value plus previous values
      ntime = 1 + Position_time_stepper_pt->nprev_values();
    }
    else
    {
      // Present value only
      ntime = 1;
    }

    // Is it a hanging node?
    if (is_hanging())
    {
      // Loop over all master nodes and update their position
      // That's all we need to update the position of hanging nodes!
      // (Recall that for hanging nodes Node::x(...) is not
      // guaranteed to be kept up-to-date; the (constrained!) nodal
      // position of hanging nodes must be determined via
      // Node::position() which determines the position
      // via the hanging node constraints from the position of
      // the master nodes)
      unsigned nmaster = hanging_pt()->nmaster();
      for (unsigned imaster = 0; imaster < nmaster; imaster++)
      {
        dynamic_cast<AlgebraicNode*>(hanging_pt()->master_node_pt(imaster))
          ->node_update();
      }
    }
    // Node isn't hanging --> update it directly
    else
    {
      // If no update function is defined, keep the nodal positions where
      // they were (i.e. don't do anything), else update
      if (nnode_update_fcts() != 0)
      {
        // Loop over time levels
        for (unsigned t = 0; t < ntime; t++)
        {
          // Update nodal position
          AlgebraicMesh* mesh_pt = Mesh_pt.begin()->second;
          // Not sure why I need this intermediate variable....
          AlgebraicNode* node_pt = this;
          mesh_pt->algebraic_node_update(t, node_pt);
        }
      }
    }

    /// Perform auxiliary update of function values?
    if (Aux_node_update_fct_pt != 0)
    {
      Aux_node_update_fct_pt(this);
    }
  }


  //========================================================================
  /// Perform self test: If the node has multiple update functions,
  /// check that all update functions give the same result (with a tolerance of
  /// AlgebraicNode::Max_allowed_difference_between_node_update_fcts
  //========================================================================
  unsigned AlgebraicNode::self_test()
  {
    // Initialise
    bool passed = true;

    unsigned test = Node::self_test();
    if (test != 0)
    {
      passed = false;
    }

    // Loop over all update functions
    unsigned nnode_update = nnode_update_fcts();

    // If there is just one (or no) update function, then no conflict
    // can arise and the error is zero.
    if (nnode_update <= 1)
    {
      // return 0;
    }
    // Multiple update funcions: check consistency
    else
    {
      // Initialise error
      double err_max = 0.0;

      // Spatial (Eulerian) position of the node
      unsigned ndim_node = ndim();
      Vector<double> x_0(ndim_node);
      Vector<double> x_new(ndim_node);

      // Get vector of update fct ids
      Vector<int> id;
      node_update_fct_id(id);


      // Quick consistency check
#ifdef PARANOID
      if (id.size() != nnode_update)
      {
        std::ostringstream error_stream;
        error_stream << "Inconsistency between number of node update ids:"
                     << nnode_update << " and " << id.size() << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // Update with first update function

      // Set default update function
      set_default_node_update(id[0]);

      // Update the node:
      node_update();

      // Store coordinates
      for (unsigned i = 0; i < ndim_node; i++)
      {
        x_0[i] = x(i);
      }


      // Loop over all other update functions
      for (unsigned iupdate = 1; iupdate < nnode_update; iupdate++)
      {
        // Set default update function
        set_default_node_update(id[iupdate]);

        // Update the node:
        node_update();

        // Store coordinates
        for (unsigned i = 0; i < ndim_node; i++)
        {
          x_new[i] = x(i);
        }

        // Check error
        double err = 0.0;
        for (unsigned i = 0; i < ndim_node; i++)
        {
          err += (x_new[i] - x_0[i]) * (x_new[i] - x_0[i]);
        }
        err = sqrt(err);
        if (err > err_max)
        {
          err_max = err;
        }

        // Dump out
        if (err > Max_allowed_difference_between_node_update_fcts)
        {
          oomph_info << "Discrepancy in algebraic update function " << iupdate
                     << ": " << x_0[0] << " " << x_0[1] << " " << x_new[0]
                     << " " << x_new[1] << std::endl;

          passed = false;
        }
      }


      // Update again with first update function to reset

      // Set default update function
      set_default_node_update(id[0]);

      // Update the node:
      node_update();

      // Return verdict
      if (passed)
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }
    // Catch all to remove compiler warning
    return 0;
  }

} // namespace oomph

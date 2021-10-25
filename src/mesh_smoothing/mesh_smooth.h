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

#ifndef OOMPH_MESH_SMOOTH_HEADER
#define OOMPH_MESH_SMOOTH_HEADER

#include <fstream>
#include <iostream>

#include "../linear_elasticity/elasticity_tensor.h"
#include "../constitutive/constitutive_laws.h"
#include "../solid/solid_traction_elements.h"


namespace oomph
{
  //======================================================================
  /// Helper namespace
  //======================================================================
  namespace Helper_namespace_for_mesh_smoothing
  {
    /// Poisson's ratio  (for smoothing by linear or nonlinear elasticity)
    double Nu = 0.3;

    /// Young's modulus (for smoothing by linear or nonlinear elasticity)
    double E = 1.0;

    /// The elasticity tensor  (for smoothing by linear elasticity)
    IsotropicElasticityTensor Isotropic_elasticity_tensor(Nu);

    /// Create constitutive law (for smoothing by nonlinear elasticity)
    ConstitutiveLaw* Constitutive_law_pt = new GeneralisedHookean(&Nu, &E);

    /// Scale for displacement of quadratic boundary (0.0: simplex; 1.0:
    /// quadratic)
    double Scale = 0.1;

    /// Increment for scale factor for displacement of quadratic boundary
    double Scale_increment = 0.1;

  } // namespace Helper_namespace_for_mesh_smoothing

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  //====================================================================
  /// Auxiliary Problem to smooth a SolidMesh by adjusting the internal
  /// nodal positions via the solution of a nonlinear solid mechanics problem.
  /// The mesh will typically have been created with an unstructured
  /// mesh generator that uses a low-order (simplex) representation of the
  /// element geometry; some of the nodes, typically non-vertex nodes on
  /// the domain's curvilinear boundaries, were then moved to their new
  /// position to provide a more accurate representation of the geometry.
  /// This class should be used to deal with elements that may have
  /// become inverted during the node motion.
  /// \n
  /// \b Important \b assumption:
  /// - It is important that the Lagrangian coordinates of all nodes still
  ///   indicate their original position, i.e. their position before (some of)
  ///   them were moved to their new position. This is because
  ///   we apply the boundary displacements in small increments.
  /// .
  /// Template argument specifies type of element. It must be a pure
  /// solid mechanics element! This shouldn't cause any problems
  /// since mesh smoothing operations tend to be performed off-line
  /// so the mesh may as well be built with pure solid elements
  /// even if it is ultimately to be used with other element types.
  /// (This restriction could easily be avoided but would
  /// require double templating and would generally be messy...)
  //====================================================================
  template<class ELEMENT>
  class NonLinearElasticitySmoothMesh : public Problem
  {
  public:
    /// Functor to update the nodal positions in SolidMesh pointed to by
    /// orig_mesh_pt in response to the displacement of some of its
    /// nodes relative to their original position which must still be indicated
    /// by the nodes' Lagrangian position. copy_of_mesh_pt must be a deep copy
    /// of orig_mesh_pt, with the same boundary coordinates etc. This mesh
    /// is used as workspace and can be deleted afterwards.
    /// The vector controlled_boundary_id contains the ids of
    /// the mesh boundaries in orig_mesh_pt whose position is supposed
    /// to remain fixed (while the other nodes are re-positioned to avoid
    /// the inversion of elements). The final optional argument
    /// specifies the max. number of increments in which the mesh
    /// boundary is deformed.
    void operator()(SolidMesh* orig_mesh_pt,
                    SolidMesh* copy_of_mesh_pt,
                    const Vector<unsigned>& controlled_boundary_id,
                    const unsigned& max_steps = 100000000)
    {
      // Dummy doc_info
      DocInfo doc_info;
      doc_info.disable_doc();
      NonLinearElasticitySmoothMesh<ELEMENT>()(orig_mesh_pt,
                                               copy_of_mesh_pt,
                                               controlled_boundary_id,
                                               doc_info,
                                               max_steps);
    }


    /// Functor to update the nodal positions in SolidMesh pointed to by
    /// orig_mesh_pt in response to the displacement of some of its
    /// nodes relative to their original position which must still be indicated
    /// by the nodes' Lagrangian position. copy_of_mesh_pt must be a deep copy
    /// of orig_mesh_pt, with the same boundary coordinates etc. This mesh
    /// is used as workspace and can be deleted afterwards.
    /// The vector controlled_boundary_id contains the ids of
    /// the mesh boundaries in orig_mesh_pt whose position is supposed
    /// to remain fixed (while the other nodes are re-positioned to avoid
    /// the inversion of elements). The DocInfo allows allows the output
    /// of the intermediate meshes. The final optional argument
    /// specifies the max. number of increments in which the mesh
    /// boundary is deformed.
    void operator()(SolidMesh* orig_mesh_pt,
                    SolidMesh* copy_of_mesh_pt,
                    const Vector<unsigned>& controlled_boundary_id,
                    DocInfo doc_info,
                    const unsigned& max_steps = 100000000)
    {
      // Make original mesh available to everyone...
      Orig_mesh_pt = orig_mesh_pt;
      Dummy_mesh_pt = copy_of_mesh_pt;

      unsigned nnode = orig_mesh_pt->nnode();
      unsigned nbound = orig_mesh_pt->nboundary();
      unsigned dim = orig_mesh_pt->node_pt(0)->ndim();

      // Add to problem's collection of sub-meshes
      add_sub_mesh(Dummy_mesh_pt);

      // Backup original nodal positions with boundary nodes snapped
      // into quadratic position; will soon move these back to
      // undeformed positon and gently move them back towards
      // their original position
      unsigned nnod = Orig_mesh_pt->nnode();
      Orig_node_pos.resize(nnod);
      for (unsigned j = 0; j < nnod; j++)
      {
        Orig_node_pos[j].resize(dim);
        SolidNode* nod_pt = dynamic_cast<SolidNode*>(Orig_mesh_pt->node_pt(j));
        for (unsigned i = 0; i < dim; i++)
        {
          Orig_node_pos[j][i] = nod_pt->x(i);
        }
      }

      // Meshes containing the face elements that represent the
      // quadratic surface
      Vector<SolidMesh*> quadratic_surface_mesh_pt(nbound);

      // GeomObject incarnations
      Vector<MeshAsGeomObject*> quadratic_surface_geom_obj_pt(nbound);


      // Create FaceElements on original mesh to define
      //-------------------------------------------------
      // the desired boundary shape
      //---------------------------

      unsigned n = controlled_boundary_id.size();
      for (unsigned i = 0; i < n; i++)
      {
        // Get boundary ID
        unsigned b = controlled_boundary_id[i];

        // Create mesh for surface elements
        quadratic_surface_mesh_pt[b] = new SolidMesh;

        // How many bulk elements are adjacent to boundary b?
        unsigned n_element = Orig_mesh_pt->nboundary_element(b);

        // Loop over the bulk elements adjacent to boundary b
        for (unsigned e = 0; e < n_element; e++)
        {
          // Get pointer to the bulk element that is adjacent to boundary b
          ELEMENT* bulk_elem_pt =
            dynamic_cast<ELEMENT*>(Orig_mesh_pt->boundary_element_pt(b, e));

          // What is the index of the face of the element e along boundary b
          int face_index = Orig_mesh_pt->face_index_at_boundary(b, e);

          // Create new element
          SolidTractionElement<ELEMENT>* el_pt =
            new SolidTractionElement<ELEMENT>(bulk_elem_pt, face_index);

          // Add it to the mesh
          quadratic_surface_mesh_pt[b]->add_element_pt(el_pt);

          // Specify boundary number
          el_pt->set_boundary_number_in_bulk_mesh(b);
        }

        // Create GeomObject incarnation
        quadratic_surface_geom_obj_pt[b] =
          new MeshAsGeomObject(quadratic_surface_mesh_pt[b]);
      }


      // Now create Lagrange multiplier elements on dummy mesh
      //-------------------------------------------------------
      Vector<SolidMesh*> dummy_lagrange_multiplier_mesh_pt(n);
      for (unsigned i = 0; i < n; i++)
      {
        // Get boundary ID
        unsigned b = controlled_boundary_id[i];

        // Make new mesh
        dummy_lagrange_multiplier_mesh_pt[i] = new SolidMesh;

        // How many bulk elements are adjacent to boundary b?
        unsigned n_element = Dummy_mesh_pt->nboundary_element(b);

        // Loop over the bulk fluid elements adjacent to boundary b?
        for (unsigned e = 0; e < n_element; e++)
        {
          // Get pointer to the bulk fluid element that is adjacent to boundary
          // b
          ELEMENT* bulk_elem_pt =
            dynamic_cast<ELEMENT*>(Dummy_mesh_pt->boundary_element_pt(b, e));

          // Find the index of the face of element e along boundary b
          int face_index = Dummy_mesh_pt->face_index_at_boundary(b, e);

          // Create new element
          ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt =
            new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
              bulk_elem_pt, face_index);

          // Add it to the mesh
          dummy_lagrange_multiplier_mesh_pt[i]->add_element_pt(el_pt);

          // Set the GeomObject that defines the boundary shape and set
          // which bulk boundary we are attached to (needed to extract
          // the boundary coordinate from the bulk nodes)
          el_pt->set_boundary_shape_geom_object_pt(
            quadratic_surface_geom_obj_pt[b], b);
        }

        // Add sub mesh
        add_sub_mesh(dummy_lagrange_multiplier_mesh_pt[i]);
      }


      // Combine the lot
      build_global_mesh();

      oomph_info << "Number of equations for nonlinear smoothing problem: "
                 << assign_eqn_numbers() << std::endl;


      // Complete the build of the elements so they are fully functional
      //----------------------------------------------------------------
      unsigned n_element = Dummy_mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        // Upcast from GeneralisedElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Dummy_mesh_pt->element_pt(e));

        // Set the constitutive law for pseudo-elastic mesh deformation
        el_pt->constitutive_law_pt() =
          Helper_namespace_for_mesh_smoothing::Constitutive_law_pt;

      } // end loop over elements


      // Output initial configuration
      doc_solution(doc_info);
      doc_info.number()++;

      // Initial scale
      Helper_namespace_for_mesh_smoothing::Scale = 0.0;
      Helper_namespace_for_mesh_smoothing::Scale_increment = 0.1;


      // Increase scale of deformation until full range is reached
      //----------------------------------------------------------
      bool done = false;
      unsigned count = 0;
      while (!done)
      {
        // Increase scale
        Helper_namespace_for_mesh_smoothing::Scale +=
          Helper_namespace_for_mesh_smoothing::Scale_increment;

        // Backup current nodal positions in dummy mesh
        backup();

        // Try it...
        bool success = true;
        try
        {
          // Avoid overshoot
          if (Helper_namespace_for_mesh_smoothing::Scale > 1.0)
          {
            Helper_namespace_for_mesh_smoothing::Scale = 1.0;
          }

          // Solve
          newton_solve();
        }
        catch (oomph::NewtonSolverError&)
        {
          success = false;
          Helper_namespace_for_mesh_smoothing::Scale -=
            Helper_namespace_for_mesh_smoothing::Scale_increment;
          Helper_namespace_for_mesh_smoothing::Scale_increment /= 2.0;

          // Reset current nodal positions in dummy mesh
          reset();
        }

        // Output solution
        if (success)
        {
          count++;
          doc_solution(doc_info);
          doc_info.number()++;
          if (Helper_namespace_for_mesh_smoothing::Scale >= 1.0) done = true;
          if (count == max_steps)
          {
            oomph_info << "Bailing out after " << count << " steps.\n";
            done = true;
          }
        }
      }

      oomph_info << "Done with Helper_namespace_for_mesh_smoothing::Scale="
                 << Helper_namespace_for_mesh_smoothing::Scale << std::endl;

      // Loop over nodes in actual mesh and assign new position
      for (unsigned j = 0; j < nnode; j++)
      {
        // Get nodes
        Node* orig_node_pt = orig_mesh_pt->node_pt(j);
        Node* new_node_pt = Dummy_mesh_pt->node_pt(j);

        // Assign new position
        for (unsigned i = 0; i < dim; i++)
        {
          orig_node_pt->x(i) = new_node_pt->x(i);
        }
      }

      // Now re-assign undeformed position
      orig_mesh_pt->set_lagrangian_nodal_coordinates();

      // Cleanup
      //--------
      n = controlled_boundary_id.size();
      for (unsigned i = 0; i < n; i++)
      {
        // Get boundary ID
        unsigned b = controlled_boundary_id[i];

        // Kill meshes and GeomObject representations
        delete quadratic_surface_mesh_pt[b];
        delete quadratic_surface_geom_obj_pt[b];
        delete dummy_lagrange_multiplier_mesh_pt[i];
      }
    }

    /// Destructor (empty)
    ~NonLinearElasticitySmoothMesh() {}


    /// Update nodal positions in main mesh -- also moves the
    /// nodes of the FaceElements that impose the new position
    void actions_before_newton_solve()
    {
      oomph_info << "Solving nonlinear smoothing problem for scale "
                 << Helper_namespace_for_mesh_smoothing::Scale << std::endl;
      unsigned nnod = Orig_mesh_pt->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        SolidNode* nod_pt = dynamic_cast<SolidNode*>(Orig_mesh_pt->node_pt(j));
        unsigned dim = nod_pt->ndim();
        for (unsigned i = 0; i < dim; i++)
        {
          nod_pt->x(i) =
            nod_pt->xi(i) + Helper_namespace_for_mesh_smoothing::Scale *
                              (Orig_node_pos[j][i] - nod_pt->xi(i));
        }
      }
    }


    /// Backup nodal positions in dummy mesh to allow for reset
    /// after non-convergence of Newton method
    void backup()
    {
      unsigned nnod = Dummy_mesh_pt->nnode();
      Backup_node_pos.resize(nnod);
      for (unsigned j = 0; j < nnod; j++)
      {
        SolidNode* nod_pt = dynamic_cast<SolidNode*>(Dummy_mesh_pt->node_pt(j));
        unsigned dim = nod_pt->ndim();
        Backup_node_pos[j].resize(dim);
        for (unsigned i = 0; i < dim; i++)
        {
          Backup_node_pos[j][i] = nod_pt->x(i);
        }
      }
    }


    /// Reset nodal positions in dummy mesh to allow for restart of
    /// Newton method with reduced increment in Scale
    void reset()
    {
      unsigned nnod = Dummy_mesh_pt->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        SolidNode* nod_pt = dynamic_cast<SolidNode*>(Dummy_mesh_pt->node_pt(j));
        unsigned dim = nod_pt->ndim();
        for (unsigned i = 0; i < dim; i++)
        {
          nod_pt->x(i) = Backup_node_pos[j][i];
        }
      }
    }

    /// Doc the solution
    void doc_solution(DocInfo& doc_info)
    {
      // Bail out
      if (!doc_info.is_doc_enabled()) return;

      std::ofstream some_file;
      std::ostringstream filename;

      // Number of plot points
      unsigned npts;
      npts = 5;

      filename << doc_info.directory() << "/smoothing_soln" << doc_info.number()
               << ".dat";

      some_file.open(filename.str().c_str());
      Dummy_mesh_pt->output(some_file, npts);
      some_file.close();

      // Check for inverted elements
      bool mesh_has_inverted_elements;
      std::ofstream inverted_fluid_elements;
      filename.str("");
      filename << doc_info.directory() << "/inverted_elements_during_smoothing"
               << doc_info.number() << ".dat";
      some_file.open(filename.str().c_str());
      Dummy_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                             some_file);
      some_file.close();
      oomph_info << "Dummy mesh does ";
      if (!mesh_has_inverted_elements) oomph_info << "not ";
      oomph_info << "have inverted elements. \n";
    }

  private:
    /// Original nodal positions
    Vector<Vector<double>> Orig_node_pos;

    /// Backup nodal positions
    Vector<Vector<double>> Backup_node_pos;

    /// Bulk original mesh
    SolidMesh* Orig_mesh_pt;

    /// Copy of mesh to work on
    SolidMesh* Dummy_mesh_pt;
  };


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  //====================================================================
  /// Auxiliary Problem to smooth a SolidMesh by adjusting the internal
  /// nodal positions by solving a LINEAR solid mechanics problem for the
  /// nodal displacements between the specified displacements of certain
  /// pinned nodes (usually located on boundaries). The template
  /// parameter specifies the linear elasticity element that must have
  /// the same shape (geometric element type) as the elements contained
  /// in the mesh that's to be smoothed. So, e.g. for the ten-noded
  /// three-dimensional tetrahedral TTaylorHoodElement<3>, it would be
  /// a TLinearElasticityElement<3,3>, etc.
  /// \b Important \b assumptions:
  /// - It is important that the Lagrangian coordinates of all nodes
  ///   still indicate their original position, i.e. their position
  ///   before (some of) them were moved to their new position.
  /// - It is assumed that in its original state, the mesh does not contain
  ///   any inverted elements.
  /// .
  //====================================================================
  template<class LINEAR_ELASTICITY_ELEMENT>
  class LinearElasticitySmoothMesh : public Problem
  {
  public:
    /// Constructor: Specify SolidMesh whose nodal positions are to
    /// be adjusted, and set of nodes in that mesh whose position
    /// are to remain fixed.
    void operator()(SolidMesh* orig_mesh_pt, std::set<Node*> pinned_nodes)
    {
      // Create new mesh and read out node/element numbers from old one
      mesh_pt() = new Mesh;
      unsigned nelem = orig_mesh_pt->nelement();
      unsigned nnode = orig_mesh_pt->nnode();

      // Have we already created that node?
      std::map<Node*, Node*> new_node;

      // Create new elements
      for (unsigned e = 0; e < nelem; e++)
      {
        // Make/add new element
        LINEAR_ELASTICITY_ELEMENT* el_pt = new LINEAR_ELASTICITY_ELEMENT;
        mesh_pt()->add_element_pt(el_pt);

        // Set elasticity tensor
        el_pt->elasticity_tensor_pt() =
          &Helper_namespace_for_mesh_smoothing::Isotropic_elasticity_tensor;

        // Find corresponding original element
        SolidFiniteElement* orig_elem_pt =
          dynamic_cast<SolidFiniteElement*>(orig_mesh_pt->finite_element_pt(e));
        unsigned nnod = orig_elem_pt->nnode();

        // Create nodes
        for (unsigned j = 0; j < nnod; j++)
        {
          // Does it not exist yet?
          if (new_node[orig_elem_pt->node_pt(j)] == 0)
          {
            Node* new_nod_pt =
              mesh_pt()->finite_element_pt(e)->construct_node(j);
            new_node[orig_elem_pt->node_pt(j)] = new_nod_pt;
            mesh_pt()->add_node_pt(new_nod_pt);
            unsigned dim = new_nod_pt->ndim();
            for (unsigned i = 0; i < dim; i++)
            {
              // Set new nodal position to be the old one in the
              // SolidMesh (assumed to contain no inverted elements)
              new_nod_pt->x(i) =
                dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))->xi(i);
            }
          }
          // It already exists -- copy across
          else
          {
            mesh_pt()->finite_element_pt(e)->node_pt(j) =
              new_node[orig_elem_pt->node_pt(j)];
          }
        }
      }

      // Loop over pinned nodes -- pin their positions and assign updated nodal
      // positions
      double scale = 1.0;
      for (std::set<Node*>::iterator it = pinned_nodes.begin();
           it != pinned_nodes.end();
           it++)
      {
        unsigned dim = (*it)->ndim();
        for (unsigned i = 0; i < dim; i++)
        {
          new_node[*it]->pin(i);
          new_node[*it]->set_value(i,
                                   scale *
                                     (dynamic_cast<SolidNode*>(*it)->x(i) -
                                      dynamic_cast<SolidNode*>(*it)->xi(i)));
        }
      }

      oomph_info << "Number of equations for smoothing problem: "
                 << assign_eqn_numbers() << std::endl;


      // Solve
      newton_solve();

      // Loop over nodes and assign displacement difference
      for (unsigned j = 0; j < nnode; j++)
      {
        // Get nodes
        SolidNode* orig_node_pt = orig_mesh_pt->node_pt(j);
        Node* new_node_pt = new_node[orig_node_pt];

        // Assign displacement difference
        unsigned dim = new_node_pt->ndim();
        for (unsigned i = 0; i < dim; i++)
        {
          orig_node_pt->x(i) = orig_node_pt->xi(i) + new_node_pt->value(i);
        }
      }

      // Now re-assign undeformed position
      orig_mesh_pt->set_lagrangian_nodal_coordinates();

      // Clean up -- mesh deletes nodes and elements
      delete mesh_pt();
    }

    /// Destructor (empty)
    ~LinearElasticitySmoothMesh() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //====================================================================
  /// Functor to smooth a SolidMesh by adjusting the internal
  /// nodal positions by solving a Poisson problem for the
  /// nodal displacements in the interior. The displacements of the specified
  /// pinned nodes (usually located on boundaries) remain fixed (their
  /// displacements are computed from the difference between their
  /// Lagrangian and Eulerian coordinates). The assumptions is
  /// that the Lagrangian coordinates in the SolidMesh still reflect
  /// the original nodal positions before the boundary nodes were
  /// moved.
  /// \n
  /// The template parameter specifies the Poisson element that must have
  /// the same shape (geometric element type) as the elements contained
  /// in the mesh that's to be smoothed. So, e.g. for the ten-noded
  /// three-dimensional tetrahedral TTaylorHoodElement<3>, it would be
  /// a TPoissonElement<3,3>, etc.
  //====================================================================
  template<class POISSON_ELEMENT>
  class PoissonSmoothMesh : public Problem
  {
  public:
    /// Functor: Specify SolidMesh whose nodal positions are to
    /// be adjusted, and set of nodes in that mesh whose position
    /// are to remain fixed.
    void operator()(SolidMesh* orig_mesh_pt, std::set<Node*> pinned_nodes)
    {
      // Create new mesh and read out node/element numbers from old one
      mesh_pt() = new Mesh;
      unsigned nelem = orig_mesh_pt->nelement();
      unsigned nnode = orig_mesh_pt->nnode();

      // Have we already created that node?
      std::map<Node*, Node*> new_node;

      // Create new elements
      for (unsigned e = 0; e < nelem; e++)
      {
        mesh_pt()->add_element_pt(new POISSON_ELEMENT);

        // Find corresponding original element
        SolidFiniteElement* orig_elem_pt =
          dynamic_cast<SolidFiniteElement*>(orig_mesh_pt->finite_element_pt(e));
        unsigned nnod = orig_elem_pt->nnode();

        // Create nodes
        for (unsigned j = 0; j < nnod; j++)
        {
          // Does it not exist yet?
          if (new_node[orig_elem_pt->node_pt(j)] == 0)
          {
            Node* new_nod_pt =
              mesh_pt()->finite_element_pt(e)->construct_node(j);
            new_node[orig_elem_pt->node_pt(j)] = new_nod_pt;
            mesh_pt()->add_node_pt(new_nod_pt);
            unsigned dim = new_nod_pt->ndim();
            for (unsigned i = 0; i < dim; i++)
            {
              // Set new nodal position to be the old one in the
              // SolidMesh (assumed to contain no inverted elements)
              new_nod_pt->x(i) =
                dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))->xi(i);
            }
          }
          // It already exists -- copy across
          else
          {
            mesh_pt()->finite_element_pt(e)->node_pt(j) =
              new_node[orig_elem_pt->node_pt(j)];
          }
        }
      }


      // Loop over pinned nodes
      for (std::set<Node*>::iterator it = pinned_nodes.begin();
           it != pinned_nodes.end();
           it++)
      {
        new_node[*it]->pin(0);
      }

      oomph_info << "Number of equations for Poisson displacement smoothing: "
                 << assign_eqn_numbers() << std::endl;

      // Solve separate displacement problems
      unsigned dim = orig_mesh_pt->node_pt(0)->ndim();
      for (unsigned i = 0; i < dim; i++)
      {
        // Loop over nodes and assign displacement difference
        for (unsigned j = 0; j < nnode; j++)
        {
          // Get nodes
          SolidNode* orig_node_pt = orig_mesh_pt->node_pt(j);
          Node* new_node_pt = new_node[orig_node_pt];

          // Assign displacement difference
          new_node_pt->set_value(0, orig_node_pt->x(i) - orig_node_pt->xi(i));
        }

        // Solve
        newton_solve();

        // Loop over nodes and assign displacement difference
        for (unsigned j = 0; j < nnode; j++)
        {
          // Get nodes
          SolidNode* orig_node_pt = orig_mesh_pt->node_pt(j);
          Node* new_node_pt = new_node[orig_node_pt];

          // Assign displacement difference
          orig_node_pt->x(i) = orig_node_pt->xi(i) + new_node_pt->value(0);
        }
      }

      // Now re-assign undeformed position
      orig_mesh_pt->set_lagrangian_nodal_coordinates();

      // Clean up -- mesh deletes nodes and elements
      delete mesh_pt();
    }
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

} // namespace oomph

#endif

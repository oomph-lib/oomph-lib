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
#ifndef OOMPH_REFINEABLE_MESH_HEADER
#define OOMPH_REFINEABLE_MESH_HEADER

#include "mesh.h"
#include "refineable_elements.h"
// Include the tree template to fill in the C++ split function
// Must be called after refineable_element.h
#include "tree.template.cc"
#include "error_estimator.h"

namespace oomph
{
  //=======================================================================
  /// Base class for refineable meshes. Provides standardised interfaces
  /// for the following standard mesh adaptation routines:
  /// - Uniform mesh refinement.
  /// - Adaptation, based on the elemental error estimates obtained
  ///   from an error estimator function which is accessed via a function
  ///   pointer.
  ///
  //=======================================================================
  class RefineableMeshBase : public virtual Mesh
  {
  public:
    bool adapt_flag()
    {
      return Adapt_flag;
    }

    /// Constructor sets default values for refinement targets etc.
    /// and initialises pointer to spatial error estimator to NULL.
    RefineableMeshBase()
    {
      // Initialise pointer to spatial error estimator
      Spatial_error_estimator_pt = 0;

      // Targets
      Min_permitted_error = 1.0e-5;
      Max_permitted_error = 1.0e-3;

      // Actual errors
      Min_error = 0.0;
      Max_error = 0.0;

      // Do I adapt or not?
      Adapt_flag = true;

      // Do I p-adapt or not?
      P_adapt_flag = true;

      // Do I disable additional synchronisation of hanging nodes?
      Additional_synchronisation_of_hanging_nodes_not_required = false;

      // Where do I write the documatation of the refinement process?
      Doc_info_pt = 0;

      // If only a few elements are scheduled for unrefinement, don't bother
      // By default unrefine all
      Max_keep_unrefined = 0;

      // Number of elements where the refinement is over-ruled
      Nrefinement_overruled = 0;
    };


    /// Broken copy constructor
    RefineableMeshBase(const RefineableMeshBase& dummy) = delete;

    /// Broken assignment operator
    void operator=(const RefineableMeshBase&) = delete;

    /// Empty Destructor:
    virtual ~RefineableMeshBase() {}

    ///  Access fct for number of elements that were refined
    unsigned nrefined()
    {
      return Nrefined;
    }

    /// Access fct for  number of elements that were unrefined
    unsigned nunrefined()
    {
      return Nunrefined;
    }

    /// Number of elements that would have liked to be refined further
    /// but can't because they've reached the max. refinement level
    unsigned& nrefinement_overruled()
    {
      return Nrefinement_overruled;
    }

    /// Max. number of elements that we allow to remain unrefined
    /// if no other mesh adaptation is required (to avoid
    /// mesh-adaptations that would only unrefine a few elements
    /// and then force a new solve -- this can't be worth our while!)
    unsigned& max_keep_unrefined()
    {
      return Max_keep_unrefined;
    }

    /// Doc the targets for mesh adaptation
    virtual void doc_adaptivity_targets(std::ostream& outfile)
    {
      outfile << std::endl;
      outfile << "Targets for mesh adaptation: " << std::endl;
      outfile << "---------------------------- " << std::endl;
      outfile << "Target for max. error: " << Max_permitted_error << std::endl;
      outfile << "Target for min. error: " << Min_permitted_error << std::endl;
      outfile << "Don't unrefine if less than " << Max_keep_unrefined
              << " elements need unrefinement." << std::endl;
      outfile << std::endl;
    }


    /// Access to spatial error estimator
    ErrorEstimator*& spatial_error_estimator_pt()
    {
      return Spatial_error_estimator_pt;
    }

    /// Access to spatial error estimator (const version
    ErrorEstimator* spatial_error_estimator_pt() const
    {
      return Spatial_error_estimator_pt;
    }

    /// Access fct for min. error (i.e. (try to) merge elements if
    /// their error is smaller)
    double& min_permitted_error()
    {
      return Min_permitted_error;
    }

    /// Access fct for max. error (i.e. split elements if their
    /// error is larger)
    double& max_permitted_error()
    {
      return Max_permitted_error;
    }

    /// Access fct for min. actual error in present solution (i.e. before
    /// re-solve on adapted mesh)
    double& min_error()
    {
      return Min_error;
    }

    /// Access fct for max. actual error in present solution (i.e. before
    /// re-solve on adapted mesh)
    double& max_error()
    {
      return Max_error;
    }

    ///  Access fct for pointer to DocInfo
    DocInfo*& doc_info_pt()
    {
      return Doc_info_pt;
    }

    /// Enable adaptation
    void enable_adaptation()
    {
      Adapt_flag = true;
    }

    /// Disable adaptation
    void disable_adaptation()
    {
      Adapt_flag = false;
    }

    /// Enable adaptation
    void enable_p_adaptation()
    {
      P_adapt_flag = true;
    }

    /// Disable adaptation
    void disable_p_adaptation()
    {
      P_adapt_flag = false;
    }

    /// Enable additional synchronisation of hanging nodes
    void enable_additional_synchronisation_of_hanging_nodes()
    {
      Additional_synchronisation_of_hanging_nodes_not_required = false;
    }

    /// Disable additional synchronisation of hanging nodes
    void disable_additional_synchronisation_of_hanging_nodes()
    {
      Additional_synchronisation_of_hanging_nodes_not_required = true;
    }

    /// Return whether the mesh is to be adapted
    bool is_adaptation_enabled() const
    {
      return Adapt_flag;
    }

    /// Return whether the mesh is to be adapted
    bool is_p_adaptation_enabled() const
    {
      return P_adapt_flag;
    }

    /// Return whether additional synchronisation is enabled
    bool is_additional_synchronisation_of_hanging_nodes_disabled() const
    {
      return Additional_synchronisation_of_hanging_nodes_not_required;
    }

    ///  Access fct for DocInfo
    DocInfo doc_info()
    {
      return *Doc_info_pt;
    }

    /// Adapt mesh: Refine elements whose error is lager than err_max
    /// and (try to) unrefine those whose error is smaller than err_min
    virtual void adapt(const Vector<double>& elemental_error) = 0;

    /// p-adapt mesh: Refine elements whose error is lager than err_max
    /// and (try to) unrefine those whose error is smaller than err_min
    virtual void p_adapt(const Vector<double>& elemental_error)
    {
      // Derived classes must implement this as required. Default throws an
      // error.
      std::ostringstream err_stream;
      err_stream << "p_adapt() called in base class RefineableMeshBase."
                 << std::endl
                 << "This needs to be implemented in the derived class."
                 << std::endl;
      throw OomphLibError(
        err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Refine mesh uniformly and doc process
    virtual void refine_uniformly(DocInfo& doc_info) = 0;

    /// Refine mesh uniformly
    virtual void refine_uniformly()
    {
      DocInfo doc_info;
      doc_info.directory() = "";
      doc_info.disable_doc();
      refine_uniformly(doc_info);
    }

    /// p-refine mesh uniformly and doc process
    virtual void p_refine_uniformly(DocInfo& doc_info)
    {
      // Derived classes must implement this as required. Default throws an
      // error.
      std::ostringstream err_stream;
      err_stream
        << "p_refine_uniformly() called in base class RefineableMeshBase."
        << std::endl
        << "This needs to be implemented in the derived class." << std::endl;
      throw OomphLibError(
        err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// p-refine mesh uniformly
    virtual void p_refine_uniformly()
    {
      DocInfo doc_info;
      doc_info.directory() = "";
      doc_info.disable_doc();
      p_refine_uniformly(doc_info);
    }

    /// Unrefine mesh uniformly: Return 0 for success,
    /// 1 for failure (if unrefinement has reached the coarsest permitted
    /// level)
    virtual unsigned unrefine_uniformly() = 0;

    /// p-unrefine mesh uniformly
    void p_unrefine_uniformly(DocInfo& doc_info)
    {
      // Derived classes must implement this as required. Default throws an
      // error.
      std::ostringstream err_stream;
      err_stream
        << "p_unrefine_uniformly() called in base class RefineableMeshBase."
        << std::endl
        << "This needs to be implemented in the derived class." << std::endl;
      throw OomphLibError(
        err_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

  protected:
    /// Pointer to spatial error estimator
    ErrorEstimator* Spatial_error_estimator_pt;

    /// Max. error (i.e. split elements if their error is larger)
    double Max_permitted_error;

    /// Min. error (i.e. (try to) merge elements if their error is smaller)
    double Min_permitted_error;

    /// Min.actual error
    double Min_error;

    /// Max. actual error
    double Max_error;

    /// Stats: Number of elements that were refined
    unsigned Nrefined;

    /// Stats: Number of elements that were unrefined
    unsigned Nunrefined;

    /// Flag that requests adaptation
    bool Adapt_flag;

    /// Flag that requests p-adaptation
    bool P_adapt_flag;

    /// Flag that disables additional synchronisation of hanging nodes
    bool Additional_synchronisation_of_hanging_nodes_not_required;

    /// Pointer to DocInfo
    DocInfo* Doc_info_pt;

    /// Max. number of elements that can remain unrefined
    /// if no other mesh adaptation is required (to avoid
    /// mesh-adaptations that would only unrefine a few elements
    /// and then force a new solve -- this can't be worth our while!)
    unsigned Max_keep_unrefined;

    /// Number of elements that would like to be refined further but
    /// can't because they've reached the max. refinement level
    unsigned Nrefinement_overruled;
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Base class for tree-based refineable meshes.
  //=======================================================================
  class TreeBasedRefineableMeshBase : public virtual RefineableMeshBase
  {
  public:
    /// Constructor
    TreeBasedRefineableMeshBase()
    {
      // Max/min refinement levels
      Max_refinement_level = 5;
      Min_refinement_level = 0;

      // Max/min p-refinement levels
      Max_p_refinement_level = 7;
      Min_p_refinement_level = 2;

      // Stats
      Nrefined = 0;
      Nunrefined = 0;

      // Where do I write the documatation of the refinement process?
      Doc_info_pt = 0;

      // Initialise the forest pointer to NULL
      Forest_pt = 0;

      // Mesh hasn't been pruned yet
      Uniform_refinement_level_when_pruned = 0;
    }


    /// Broken copy constructor
    TreeBasedRefineableMeshBase(const TreeBasedRefineableMeshBase& dummy) =
      delete;

    /// Broken assignment operator
    void operator=(const TreeBasedRefineableMeshBase&) = delete;

    /// Empty Destructor:
    virtual ~TreeBasedRefineableMeshBase()
    {
      // Kill the forest if there is one
      if (Forest_pt != 0)
      {
        delete Forest_pt;
        Forest_pt = 0;
      }
    }

    /// Adapt mesh: Refine elements whose error is lager than err_max
    /// and (try to) unrefine those whose error is smaller than err_min
    void adapt(const Vector<double>& elemental_error);

    /// p-adapt mesh: Refine elements whose error is lager than err_max
    /// and (try to) unrefine those whose error is smaller than err_min
    void p_adapt(const Vector<double>& elemental_error);

    /// Refine mesh uniformly and doc process
    void refine_uniformly(DocInfo& doc_info);

    /// Refine mesh uniformly
    void refine_uniformly()
    {
      RefineableMeshBase::refine_uniformly();
    }

    /// p-refine mesh uniformly and doc process
    void p_refine_uniformly(DocInfo& doc_info);

    /// p-refine mesh uniformly
    void p_refine_uniformly()
    {
      RefineableMeshBase::p_refine_uniformly();
    }

    /// Unrefine mesh uniformly: Return 0 for success,
    /// 1 for failure (if unrefinement has reached the coarsest permitted
    /// level)
    unsigned unrefine_uniformly();

    /// p-unrefine mesh uniformly
    void p_unrefine_uniformly(DocInfo& doc_info);

    /// Set up the tree forest associated with the Mesh (if any)
    virtual void setup_tree_forest() = 0;

    /// Return pointer to the Forest represenation of the mesh
    TreeForest* forest_pt()
    {
      return Forest_pt;
    }


    /// Doc the targets for mesh adaptation
    void doc_adaptivity_targets(std::ostream& outfile)
    {
      outfile << std::endl;
      outfile << "Targets for mesh adaptation: " << std::endl;
      outfile << "---------------------------- " << std::endl;
      outfile << "Target for max. error: " << Max_permitted_error << std::endl;
      outfile << "Target for min. error: " << Min_permitted_error << std::endl;
      outfile << "Min. refinement level: " << Min_refinement_level << std::endl;
      outfile << "Max. refinement level: " << Max_refinement_level << std::endl;
      outfile << "Min. p-refinement level: " << Min_p_refinement_level
              << std::endl;
      outfile << "Max. p-refinement level: " << Max_p_refinement_level
              << std::endl;
      outfile << "Don't unrefine if less than " << Max_keep_unrefined
              << " elements need unrefinement." << std::endl;
      outfile << std::endl;
    }

    /// Access fct for max. permissible refinement level (relative to base mesh)
    unsigned& max_refinement_level()
    {
      return Max_refinement_level;
    }

    /// Access fct for min. permissible refinement level (relative to base mesh)
    unsigned& min_refinement_level()
    {
      return Min_refinement_level;
    }

    /// Access fct for max. permissible p-refinement level (relative to base
    /// mesh)
    unsigned& max_p_refinement_level()
    {
      return Max_p_refinement_level;
    }

    /// Access fct for min. permissible p-refinement level (relative to base
    /// mesh)
    unsigned& min_p_refinement_level()
    {
      return Min_p_refinement_level;
    }

    /// Perform the actual tree-based mesh adaptation,
    /// documenting the progress in the directory specified in DocInfo object.
    virtual void adapt_mesh(DocInfo& doc_info);

    /// Perform the actual tree-based mesh adaptation. A simple wrapper
    /// to call the function without documentation.
    virtual void adapt_mesh()
    {
      // Create a dummy doc_info object
      DocInfo doc_info;
      doc_info.directory() = "";
      doc_info.disable_doc();
      // Call the other adapt mesh
      adapt_mesh(doc_info);
    }

    /// Perform the actual tree-based mesh p-adaptation,
    /// documenting the progress in the directory specified in DocInfo object.
    void p_adapt_mesh(DocInfo& doc_info);

    /// Perform the actual tree-based mesh p-adaptation. A simple wrapper
    /// to call the function without documentation.
    void p_adapt_mesh()
    {
      // Create a dummy doc_info object
      DocInfo doc_info;
      doc_info.directory() = "";
      doc_info.disable_doc();
      // Call the other adapt mesh
      p_adapt_mesh(doc_info);
    }

    /// Refine mesh by splitting the elements identified
    /// by their numbers.
    virtual void refine_selected_elements(
      const Vector<unsigned>& elements_to_be_refined);

    /// Refine mesh by splitting the elements identified
    /// by their pointers.
    virtual void refine_selected_elements(
      const Vector<RefineableElement*>& elements_to_be_refined);

    /// p-refine mesh by refining the elements identified
    /// by their numbers.
    void p_refine_selected_elements(
      const Vector<unsigned>& elements_to_be_refined);

    /// p-refine mesh by refining the elements identified
    /// by their pointers.
    void p_refine_selected_elements(
      const Vector<PRefineableElement*>& elements_to_be_refined_pt);


    /// Refine base mesh to same degree as reference mesh (relative
    /// to original unrefined mesh).
    virtual void refine_base_mesh_as_in_reference_mesh(
      TreeBasedRefineableMeshBase* const& ref_mesh_pt);

    /// Refine base mesh to same degree as reference mesh minus one
    /// level of refinement (relative to original unrefined mesh). Useful
    /// function for multigrid solvers; allows the easy copy of a mesh
    /// to the level of refinement just below the current one
    virtual bool refine_base_mesh_as_in_reference_mesh_minus_one(
      TreeBasedRefineableMeshBase* const& ref_mesh_pt);

    /// Refine mesh once so that its topology etc becomes that of the
    /// (finer!) reference mesh -- if possible! Useful for meshes in multigrid
    /// hierarchies. If the meshes are too different and the conversion
    /// cannot be performed, the code dies (provided PARANOID is enabled).
    virtual void refine_as_in_reference_mesh(
      TreeBasedRefineableMeshBase* const& ref_mesh_pt);

    /// Get max/min refinement levels in mesh
    virtual void get_refinement_levels(unsigned& min_refinement_level,
                                       unsigned& max_refinement_level);

    /// Extract the elements at a particular refinement level in
    /// the refinement pattern - used in Mesh::redistribute or whatever it's
    /// going to be called (RefineableMeshBase::reduce_halo_layers or something)
    virtual void get_elements_at_refinement_level(
      unsigned& refinement_level, Vector<RefineableElement*>& level_elements);

    /// Extract refinement pattern: Consider the hypothetical mesh
    /// obtained by truncating the refinement of the current mesh to a given
    /// level (where \c level=0 is the un-refined base mesh). To advance
    /// to the next refinement level, we need to refine (split) the
    /// \c to_be_refined[level].size() elements identified by the
    /// element numbers contained in \c vector to_be_refined[level][...]
    virtual void get_refinement_pattern(
      Vector<Vector<unsigned>>& to_be_refined);

    /// Refine base mesh according to specified refinement pattern
    void refine_base_mesh(Vector<Vector<unsigned>>& to_be_refined);

    /// Refine mesh according to refinement pattern in restart file
    virtual void refine(std::ifstream& restart_file);

    /// Dump refinement pattern to allow for rebuild
    virtual void dump_refinement(std::ostream& outfile);

    /// Read refinement pattern to allow for rebuild
    virtual void read_refinement(std::ifstream& restart_file,
                                 Vector<Vector<unsigned>>& to_be_refined);


    /// Level to which the mesh was uniformly refined when it was pruned
    /// (const version)
    unsigned uniform_refinement_level_when_pruned() const
    {
      return Uniform_refinement_level_when_pruned;
    }


    /// Level to which the mesh was uniformly refined when it was pruned
    unsigned& uniform_refinement_level_when_pruned()
    {
      return Uniform_refinement_level_when_pruned;
    }

#ifdef OOMPH_HAS_MPI

    /// Classify all halo and haloed information in the mesh (overloaded
    /// version from Mesh base class. Calls that one first, then synchronises
    /// hanging nodes)
    void classify_halo_and_haloed_nodes(DocInfo& doc_info,
                                        const bool& report_stats)
    {
      // Call version in base class but don't bother to call
      // resize_halo_nodes() -- we'll do it ourselves below
      bool backup = Resize_halo_nodes_not_required;
      Resize_halo_nodes_not_required = false;
      Mesh::classify_halo_and_haloed_nodes(doc_info, report_stats);
      Resize_halo_nodes_not_required = backup;

      double t_start = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_start = TimingHelpers::timer();
      }

      // Communicate positions of non-hanging nodes that depend on non-existent
      // neighbours (e.g. h-refinement of neighbouring elements with different
      // p-orders where the shared edge is on the outer edge of the halo layer)
      synchronise_nonhanging_nodes();

      // Get number of continously interpolated variables
      unsigned local_ncont_interpolated_values = 0;

      // Bypass if we have no elements and then get overall
      // value by reduction
      if (nelement() > 0)
      {
        local_ncont_interpolated_values =
          dynamic_cast<RefineableElement*>(element_pt(0))
            ->ncont_interpolated_values();
      }
      unsigned ncont_interpolated_values = 0;
      MPI_Allreduce(&local_ncont_interpolated_values,
                    &ncont_interpolated_values,
                    1,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    Comm_pt->mpi_comm());

      // Synchronise the hanging nodes
      synchronise_hanging_nodes(ncont_interpolated_values);

      if (Global_timings::Doc_comprehensive_timings)
      {
        double t_end = TimingHelpers::timer();
        oomph_info
          << "Time for "
          << "TreeBasedRefineableMeshBase::synchronise_hanging_nodes() "
          << "incl. time for initial allreduce in "
          << "TreeBasedRefineableMeshBase::classify_halo_and_haloed_nodes(): "
          << t_end - t_start << std::endl;
      }

      // Now do resize halo nodes after all (unless it's been
      // declared to be unnecessary by somebody...)
      if (!Resize_halo_nodes_not_required)
      {
        resize_halo_nodes();
      }
    }


    /// Classify the halo and haloed nodes in the mesh.
    void classify_halo_and_haloed_nodes(const bool& report_stats = false)
    {
      DocInfo doc_info;
      doc_info.disable_doc();
      classify_halo_and_haloed_nodes(doc_info, report_stats);
    }

#endif

  protected:
#ifdef OOMPH_HAS_MPI

    /// Synchronise the hanging nodes if the mesh is distributed
    void synchronise_hanging_nodes(const unsigned& ncont_interpolated_values);

    /// Additional synchronisation of hanging nodes
    /// Required for reconcilliation of hanging nodes on the outer edge of the
    /// halo layer when using elements with nonuniformly spaced nodes
    /// Must be implemented in templated derived class because ELEMENT template
    /// parameter is required for the node-creation functions in the
    /// Missing_masters_functions namespace
    virtual void additional_synchronise_hanging_nodes(
      const unsigned& ncont_interpolated_values) = 0;

    /// Synchronise the positions of non-hanging nodes that depend on
    /// non-existent neighbours (e.g. h-refinement of neighbouring elements
    /// with different p-orders where the shared edge is on the outer edge of
    /// the halo layer)
    void synchronise_nonhanging_nodes();

#endif

    /// Split all the elements in the mesh if required. This template
    /// free interface will be overloaded in RefineableMesh<ELEMENT> so that any
    /// new elements that are created will be of the correct type.
    virtual void split_elements_if_required() = 0;

    /// p-refine all the elements in the mesh if required. This template
    /// free interface will be overloaded in RefineableMesh<ELEMENT> so that
    /// any temporary copies of the element that are created will be of the
    /// correct type.
    virtual void p_refine_elements_if_required() = 0;

    /// Complete the hanging node scheme recursively
    void complete_hanging_nodes(const int& ncont_interpolated_values);


    /// Auxiliary routine for recursive hanging node completion
    void complete_hanging_nodes_recursively(Node*& nod_pt,
                                            Vector<Node*>& master_nodes,
                                            Vector<double>& hang_weights,
                                            const int& ival);

    /// Level to which the mesh was uniformly refined when it was pruned
    unsigned Uniform_refinement_level_when_pruned;

    /// Max. permissible refinement level (relative to base mesh)
    unsigned Max_refinement_level;

    /// Min. permissible refinement level (relative to base mesh)
    unsigned Min_refinement_level;

    /// Max. permissible p-refinement level (relative to base mesh)
    unsigned Max_p_refinement_level;

    /// Min. permissible p-refinement level (relative to base mesh)
    unsigned Min_p_refinement_level;

    /// Forest representation of the mesh
    TreeForest* Forest_pt;

  private:
#ifdef OOMPH_HAS_MPI

    /// Helper struct to collate data required during
    /// TreeBasedRefineableMeshBase::synchronise_hanging_nodes
    struct HangHelperStruct
    {
      unsigned Sending_processor;
      unsigned Shared_node_id_on_sending_processor;
      unsigned Shared_node_proc;
      double Weight;
      HangInfo* Hang_pt;
      unsigned Master_node_index;
      // Store pointer to node and index of value in case translation attempt
      // fails (see synchronise_hanging_nodes())
      Node* Node_pt;
      int icont;
    };

#endif
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Templated base class for refineable meshes. The use of the template
  /// parameter is required only for creating new elements during mesh
  /// adaptation. This class overloaded the template-free inteface to
  /// the function split_elements_if_required() to make use of the template
  /// parameter.
  /// All refineable meshes should inherit directly from
  /// TreeBasedRefineableMesh<ELEMENT>
  //=======================================================================
  template<class ELEMENT>
  class TreeBasedRefineableMesh : public virtual TreeBasedRefineableMeshBase
  {
  private:
    /// Split all the elements if required. Overload the template-free
    /// interface so that any new elements that are created
    /// will be of the correct type.
    void split_elements_if_required()
    {
      // Find the number of trees in the forest
      unsigned n_tree = this->Forest_pt->ntree();
      // Loop over all "active" elements in the forest and split them
      // if required
      for (unsigned long e = 0; e < n_tree; e++)
      {
        this->Forest_pt->tree_pt(e)->traverse_leaves(
          &Tree::split_if_required<ELEMENT>);
      }
    }

    /// p-refine all the elements if required. Overload the template-free
    /// interface so that any temporary copies of the element that are created
    /// will be of the correct type.
    void p_refine_elements_if_required()
    {
      // BENFLAG: Make a non const pointer to the mesh so it can be passed
      // (HACK)
      Mesh* mesh_pt = this;
      // Find the number of trees in the forest
      unsigned n_tree = this->Forest_pt->ntree();
      // Loop over all "active" elements in the forest and p-refine them
      // if required
      for (unsigned long e = 0; e < n_tree; e++)
      {
        this->Forest_pt->tree_pt(e)->traverse_leaves(
          &Tree::p_refine_if_required<ELEMENT>, mesh_pt);
      }
    }

  protected:
#ifdef OOMPH_HAS_MPI

    /// Additional setup of shared node scheme
    /// This is Required for reconcilliation of hanging nodes acrross processor
    /// boundaries when using elements with nonuniformly spaced nodes.
    /// ELEMENT template parameter is required so that
    /// MacroElementNodeUpdateNodes which are added as external halo master
    /// nodes can be made fully functional
    void additional_synchronise_hanging_nodes(
      const unsigned& ncont_interpolated_values);

#endif

  public:
    /// Constructor, call the constructor of the base class
    TreeBasedRefineableMesh() : TreeBasedRefineableMeshBase() {}

    /// Broken copy constructor
    TreeBasedRefineableMesh(const TreeBasedRefineableMesh& dummy) = delete;

    /// Empty virtual destructor
    virtual ~TreeBasedRefineableMesh() {}
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Base class for refineable tet meshes
  //=======================================================================
  class RefineableTetMeshBase : public virtual RefineableMeshBase
  {
  public:
    /// Max element size allowed during adaptation
    double& max_element_size()
    {
      return Max_element_size;
    }

    /// Min element size allowed during adaptation
    double& min_element_size()
    {
      return Min_element_size;
    }

    /// Min edge ratio before remesh gets triggered
    double& max_permitted_edge_ratio()
    {
      return Max_permitted_edge_ratio;
    }


    /// Doc the targets for mesh adaptation
    void doc_adaptivity_targets(std::ostream& outfile)
    {
      outfile << std::endl;
      outfile << "Targets for mesh adaptation: " << std::endl;
      outfile << "---------------------------- " << std::endl;
      outfile << "Target for max. error: " << Max_permitted_error << std::endl;
      outfile << "Target for min. error: " << Min_permitted_error << std::endl;
      outfile << "Target max edge ratio: " << Max_permitted_edge_ratio
              << std::endl;
      outfile << "Min. allowed element size: " << Min_element_size << std::endl;
      outfile << "Max. allowed element size: " << Max_element_size << std::endl;
      outfile << "Don't unrefine if less than " << Max_keep_unrefined
              << " elements need unrefinement." << std::endl;
      outfile << std::endl;
    }


    /// Compute target volume based on the elements' error and the
    /// error target; return max edge ratio
    double compute_volume_target(const Vector<double>& elem_error,
                                 Vector<double>& target_volume)
    {
      double max_edge_ratio = 0.0;
      unsigned count_unrefined = 0;
      unsigned count_refined = 0;
      this->Nrefinement_overruled = 0;

      unsigned nel = this->nelement();
      for (unsigned e = 0; e < nel; e++)
      {
        // Get element
        FiniteElement* el_pt = this->finite_element_pt(e);

        // Calculate the volume of the element
        double volume = el_pt->size();

        // Find the vertex coordinates
        // (vertices are enumerated first)
        double vertex[4][3];
        for (unsigned n = 0; n < 4; ++n)
        {
          for (unsigned i = 0; i < 3; ++i)
          {
            vertex[n][i] = el_pt->node_pt(n)->x(i);
          }
        }

        // Compute the radius of the circumsphere of the tetrahedron
        // Algorithm stolen from tetgen for consistency
        DenseDoubleMatrix A(3);
        for (unsigned i = 0; i < 3; ++i)
        {
          A(0, i) = vertex[1][i] - vertex[0][i];
          A(1, i) = vertex[2][i] - vertex[0][i];
          A(2, i) = vertex[3][i] - vertex[0][i];
        }

        Vector<double> rhs(3);
        // Compute the right hand side vector b (3x1).
        for (unsigned i = 0; i < 3; ++i)
        {
          rhs[i] = 0.0;
          for (unsigned k = 0; k < 3; ++k)
          {
            rhs[i] += A(i, k) * A(i, k);
          }
          rhs[i] *= 0.5;
        }

        // Solve the linear system, in which the rhs is over-written with
        // the solution
        A.solve(rhs);

        // Calculate the circum-radius
        double circum_radius =
          sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);

        // Now find the shortest edge length
        Vector<double> edge(3);
        double min_length = DBL_MAX;
        for (unsigned start = 0; start < 4; ++start)
        {
          for (unsigned end = start + 1; end < 4; ++end)
          {
            for (unsigned i = 0; i < 3; ++i)
            {
              edge[i] = vertex[start][i] - vertex[end][i];
            }
            double length =
              sqrt(edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2]);
            if (length < min_length)
            {
              min_length = length;
            }
          }
        }

        // Now calculate the minimum edge ratio for this element
        double local_max_edge_ratio = circum_radius / min_length;
        if (local_max_edge_ratio > max_edge_ratio)
        {
          max_edge_ratio = local_max_edge_ratio;
        }

        // Mimick refinement in tree-based procedure: Target volumes
        // for elements that exceed permitted error is 1/4 of their
        // current volume, corresponding to a uniform sub-division.
        if (elem_error[e] > this->max_permitted_error())
        {
          target_volume[e] = std::max(volume / 4.0, Min_element_size);
          if (target_volume[e] != Min_element_size)
          {
            count_refined++;
          }
          else
          {
            this->Nrefinement_overruled++;
          }
        }
        else if (elem_error[e] < this->min_permitted_error())
        {
          target_volume[e] = std::min(4.0 * volume, Max_element_size);
          if (target_volume[e] != Max_element_size)
          {
            count_unrefined++;
          }
        }
        else
        {
          // Leave it alone
          target_volume[e] = std::max(volume, Min_element_size);
        }

      } // End of loop over elements

      // Tell everybody
      this->Nrefined = count_refined;
      this->Nunrefined = count_unrefined;

      return max_edge_ratio;
    }

    /// Max permitted element size
    double Max_element_size;

    /// Min permitted element size
    double Min_element_size;

    /// Max edge ratio before remesh gets triggered
    double Max_permitted_edge_ratio;
  };


} // namespace oomph

#endif

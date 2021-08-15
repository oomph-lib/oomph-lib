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


#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <cstdlib>
#include <stdlib.h>
#include <limits>

#include "refineable_mesh.h"
// Include to fill in additional_synchronise_hanging_nodes() function
#include "refineable_mesh.template.cc"

namespace oomph
{
  //========================================================================
  /// Get refinement pattern of mesh: Consider the hypothetical mesh
  /// obtained by truncating the refinement of the current mesh to a given level
  /// (where \c level=0 is the un-refined base mesh). To advance
  /// to the next refinement level, we need to refine (split) the
  /// \c to_be_refined[level].size() elements identified by the
  /// element numbers contained in \c vector to_be_refined[level][...]
  //========================================================================
  void TreeBasedRefineableMeshBase::get_refinement_pattern(
    Vector<Vector<unsigned>>& to_be_refined)
  {
    // Extract *all* elements from current (fully refined) mesh.
    Vector<Tree*> all_tree_nodes_pt;
    forest_pt()->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

    // Find out maximum refinement level
    unsigned max_level = 0;
    unsigned nnodes = all_tree_nodes_pt.size();
    for (unsigned e = 0; e < nnodes; e++)
    {
      unsigned level = all_tree_nodes_pt[e]->level();
      if (level > max_level) max_level = level;
    }

    // Assign storage for refinement pattern
    to_be_refined.clear();
    to_be_refined.resize(max_level);
    Vector<unsigned> el_count(max_level);

    // Initialise count of elements that exist in mesh when refinement
    // has proceeded to this level
    for (unsigned l = 0; l < max_level; l++)
    {
      el_count[l] = 0;
    }

    // Loop over all levels and extract all elements that exist
    // in reference mesh when refinement has proceeded to this level
    for (unsigned l = 0; l < max_level; l++)
    {
      // Loop over all elements (tree nodes)
      for (unsigned e = 0; e < nnodes; e++)
      {
        // What level does this element exist on?
        unsigned level = all_tree_nodes_pt[e]->level();

        // Element is part of the mesh at this refinement level
        // if it exists at this level OR if it exists at a lower level
        // and is a leaf
        if ((level == l) || ((level < l) && (all_tree_nodes_pt[e]->is_leaf())))
        {
          // If element exsts at this level and is not a leaf it will
          // be refined when we move to the next level:
          if ((level == l) && (!all_tree_nodes_pt[e]->is_leaf()))
          {
            // Add element number (in mesh at current refinement level)
            // to the list of elements that need to be refined
            to_be_refined[l].push_back(el_count[l]);
          }
          // Element exists in this mesh: Add to counter
          el_count[l]++;
        }
      }
    }
  }

  //========================================================================
  /// \short Extract the elements at a particular refinement level in
  /// the refinement pattern (used in Mesh::redistribute or whatever it's
  /// going to be called (RefineableMeshBase::reduce_halo_layers or something)
  //========================================================================
  void TreeBasedRefineableMeshBase::get_elements_at_refinement_level(
    unsigned& refinement_level, Vector<RefineableElement*>& level_elements)
  {
    // Extract *all* elements from current (fully refined) mesh.
    Vector<Tree*> all_tree_nodes_pt;
    forest_pt()->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

    // Add the element to the vector if its level matches refinement_level
    unsigned nnodes = all_tree_nodes_pt.size();
    for (unsigned e = 0; e < nnodes; e++)
    {
      unsigned level = all_tree_nodes_pt[e]->level();
      if (level == refinement_level)
      {
        level_elements.push_back(
          dynamic_cast<RefineableElement*>(all_tree_nodes_pt[e]->object_pt()));
      }
    }
  }

  //========================================================================
  /// Refine original, unrefined mesh according to specified refinement
  /// pattern (relative to original, unrefined mesh).
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_base_mesh(
    Vector<Vector<unsigned>>& to_be_refined)
  {
    // Get mesh back to unrefined state
    unsigned my_max, my_min;
    get_refinement_levels(my_min, my_max);

    // Max refinement level:
    unsigned my_max_level = to_be_refined.size();

    unsigned global_max = 0;
    unsigned global_max_level = 0;
    Vector<unsigned> data(2, 0);
    data[0] = my_max;
    data[1] = my_max_level;
    Vector<unsigned> global_data(2, 0);
#ifdef OOMPH_HAS_MPI
    if (this->is_mesh_distributed())
    {
      MPI_Allreduce(&data[0],
                    &global_data[0],
                    2,
                    MPI_UNSIGNED,
                    MPI_MAX,
                    Comm_pt->mpi_comm());
      global_max = global_data[0];
      global_max_level = global_data[1];
    }
    else
#endif
    {
      global_max = my_max;
      global_max_level = my_max_level;
    }


    for (unsigned i = 0; i < global_max; i++)
    {
      unrefine_uniformly();
    }

    // Do refinement steps in current mesh
    for (unsigned l = 0; l < global_max_level; l++)
    {
      // Loop over elements that need to be refined at this level
      unsigned n_to_be_refined = 0;
      if (l < my_max_level) n_to_be_refined = to_be_refined[l].size();

      // Select relevant elements to be refined
      for (unsigned i = 0; i < n_to_be_refined; i++)
      {
        dynamic_cast<RefineableElement*>(this->element_pt(to_be_refined[l][i]))
          ->select_for_refinement();
      }

      // Now do the actual mesh refinement
      adapt_mesh();
    }
  }


  //========================================================================
  /// Refine base mesh according to refinement pattern in restart file
  //========================================================================
  void TreeBasedRefineableMeshBase::refine(std::ifstream& restart_file)
  {
    // Assign storage for refinement pattern
    Vector<Vector<unsigned>> to_be_refined;

    // Read refinement pattern
    read_refinement(restart_file, to_be_refined);

    // Refine
    refine_base_mesh(to_be_refined);
  }


  //========================================================================
  /// Dump refinement pattern to allow for rebuild
  ///
  //========================================================================
  void TreeBasedRefineableMeshBase::dump_refinement(std::ostream& outfile)
  {
    // Assign storage for refinement pattern
    Vector<Vector<unsigned>> to_be_refined;

    // Get refinement pattern of reference mesh:
    get_refinement_pattern(to_be_refined);

    // Dump max refinement level:
    unsigned max_level = to_be_refined.size();
    outfile << max_level << " # max. refinement level " << std::endl;

    // Doc the numbers of the elements that need to be refined at this level
    for (unsigned l = 0; l < max_level; l++)
    {
      // Loop over elements that need to be refined at this level
      unsigned n_to_be_refined = to_be_refined[l].size();
      outfile << n_to_be_refined << " # number of elements to be refined. "
              << "What follows are the numbers of the elements. " << std::endl;

      // Select relevant elements to be refined
      for (unsigned i = 0; i < n_to_be_refined; i++)
      {
        outfile << to_be_refined[l][i] << std::endl;
      }
    }
  }


  //========================================================================
  /// Read refinement pattern to allow for rebuild
  ///
  //========================================================================
  void TreeBasedRefineableMeshBase::read_refinement(
    std::ifstream& restart_file, Vector<Vector<unsigned>>& to_be_refined)
  {
    std::string input_string;

    // Read max refinement level:

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Convert
    unsigned max_level = std::atoi(input_string.c_str());

    // Assign storage for refinement pattern
    to_be_refined.resize(max_level);

    // Read the number of the elements that need to be refined at different
    // levels
    for (unsigned l = 0; l < max_level; l++)
    {
      // Read line up to termination sign
      getline(restart_file, input_string, '#');

      // Ignore rest of line
      restart_file.ignore(80, '\n');

      // Convert
      unsigned n_to_be_refined = atoi(input_string.c_str());
      ;

      // Assign storage
      to_be_refined[l].resize(n_to_be_refined);

      // Read numbers of the elements that need to be refined
      for (unsigned i = 0; i < n_to_be_refined; i++)
      {
        restart_file >> to_be_refined[l][i];
      }
    }
  }


  //========================================================================
  /// Do adaptive refinement for mesh.
  /// - Pass Vector of error estimates for all elements.
  /// - Refine those whose errors exceeds the threshold
  /// - (Try to) unrefine those whose errors is less than
  ///   threshold (only possible if the three brothers also want to be
  ///   unrefined, of course.)
  /// - Update the nodal positions in the whole lot
  /// - Store # of refined/unrefined elements.
  /// - Doc refinement process (if required)
  //========================================================================
  void TreeBasedRefineableMeshBase::adapt(const Vector<double>& elemental_error)
  {
    // Set the refinement tolerance to be the max permissible error
    double refine_tol = max_permitted_error();

    // Set the unrefinement tolerance to be the min permissible error
    double unrefine_tol = min_permitted_error();

    // Setup doc info
    DocInfo local_doc_info;
    if (doc_info_pt() == 0)
    {
      local_doc_info.disable_doc();
    }
    else
    {
      local_doc_info = doc_info();
    }


    // Check that the errors make sense
    if (refine_tol <= unrefine_tol)
    {
      std::ostringstream error_stream;
      error_stream << "Refinement tolerance <= Unrefinement tolerance"
                   << refine_tol << " " << unrefine_tol << std::endl
                   << "doesn't make sense and will almost certainly crash"
                   << std::endl
                   << "this beautiful code!" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Select elements for refinement and unrefinement
    //================================================
    // Reset counter for number of elements that would like to be
    // refined further but can't
    nrefinement_overruled() = 0;

    // Note: Yes, this needs to be a map because we'll have to check
    // the refinement wishes of brothers (who we only access via pointers)
    std::map<RefineableElement*, bool> wants_to_be_unrefined;

    // Initialise a variable to store the number of elements for refinment
    unsigned n_refine = 0;

    // Loop over all elements and mark them according to the error criterion
    unsigned long Nelement = this->nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // (Cast) pointer to the element
      RefineableElement* el_pt =
        dynamic_cast<RefineableElement*>(this->element_pt(e));

      // Initially element is not to be refined
      el_pt->deselect_for_refinement();

      // If the element error exceeds the threshold ...
      if (elemental_error[e] > refine_tol)
      {
        // ... and its refinement level is less than the maximum desired level
        // mark is to be refined
        if ((el_pt->refinement_is_enabled()) &&
            (el_pt->refinement_level() < max_refinement_level()))
        {
          el_pt->select_for_refinement();
          n_refine++;
        }
        // ... otherwise mark it as having been over-ruled
        else
        {
          nrefinement_overruled() += 1;
        }
      }

      // Now worry about unrefinement (first pass):

      // Is my error too small AND do I have a father?
      if ((elemental_error[e] < unrefine_tol) &&
          (el_pt->tree_pt()->father_pt() != 0))
      {
        // Flag to indicate whether to unrefine
        bool unrefine = true;
        unsigned n_sons = el_pt->tree_pt()->father_pt()->nsons();

        // Are all brothers leaf nodes?
        for (unsigned ison = 0; ison < n_sons; ison++)
        {
          // (At least) one brother is not a leaf: end of story; we're not doing
          // it (= the unrefinement)
          if (!(el_pt->tree_pt()->father_pt()->son_pt(ison)->is_leaf()))
          {
            unrefine = false;
          }
        }

        // Don't allow unrefinement of elements that would become larger
        // than the minimum legal refinement level
        if (el_pt->refinement_level() - 1 < min_refinement_level())
        {
          unrefine = false;
        }

        // So, all things considered, is the element eligbible for refinement?
        if (unrefine)
        {
          wants_to_be_unrefined[el_pt] = true;
        }
        else
        {
          wants_to_be_unrefined[el_pt] = false;
        }
      }
    }

    oomph_info << " \n Number of elements to be refined: " << n_refine
               << std::endl;
    oomph_info << " \n Number of elements whose refinement was overruled: "
               << nrefinement_overruled() << std::endl;

    // Second pass for unrefinement --- an element cannot be unrefined unless
    // all brothers want to be unrefined.
    // Loop over all elements again and let the first set of sons check if their
    // brothers also want to be unrefined
    unsigned n_unrefine = 0;
    for (unsigned long e = 0; e < Nelement; e++)
    {
      //(Cast) pointer to the element
      RefineableElement* el_pt =
        dynamic_cast<RefineableElement*>(this->element_pt(e));

      // hierher: This is a bit naughty... We want to put the
      // first son in charge -- the statement below assumes (correctly) that the
      // enumeration of all (!) trees starts with son types.
      // This is correct for oc and quadtrees but will bite us if we
      // ever introduce other trees if/when we accidentally break this
      // tacit assumption. Not sure what to do about it for
      // now other than leaving it hierher-ed...
      if (el_pt->tree_pt()->son_type() == OcTreeNames::LDB)
      {
        // Do all sons want to be unrefined?
        bool unrefine = true;
        unsigned n_sons = el_pt->tree_pt()->father_pt()->nsons();
        for (unsigned ison = 0; ison < n_sons; ison++)
        {
          if (!(wants_to_be_unrefined[dynamic_cast<RefineableElement*>(
                el_pt->tree_pt()->father_pt()->son_pt(ison)->object_pt())]))
          {
            // One guy isn't cooperating and spoils the party.
            unrefine = false;
          }
        }

        // Tell father that his sons need to be merged
        if (unrefine)
        {
          el_pt->tree_pt()
            ->father_pt()
            ->object_pt()
            ->select_sons_for_unrefinement();
          n_unrefine += n_sons;
        }
        // Otherwise mark the sons as not to be touched
        else
        {
          el_pt->tree_pt()
            ->father_pt()
            ->object_pt()
            ->deselect_sons_for_unrefinement();
        }
      }
    }
    oomph_info << " \n Number of elements to be merged : " << n_unrefine
               << std::endl
               << std::endl;


    // Now do the actual mesh adaptation
    //---------------------------------

    // Check whether its worth our while
    // Either some elements want to be refined,
    // or the number that want to be unrefined are greater than the
    // specified tolerance

    // In a parallel job, it is possible that one process may not have
    // any elements to refine, BUT a neighbouring process may refine an
    // element which changes the hanging status of a node that is on
    // both processes (i.e. a halo(ed) node).  To get around this issue,
    // ALL processes need to call adapt_mesh if ANY refinement is to
    // take place anywhere.

    unsigned total_n_refine = 0;
#ifdef OOMPH_HAS_MPI
    // Sum n_refine across all processors
    if (this->is_mesh_distributed())
    {
      MPI_Allreduce(&n_refine,
                    &total_n_refine,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    Comm_pt->mpi_comm());
    }
    else
    {
      total_n_refine = n_refine;
    }
#else
    total_n_refine = n_refine;
#endif

    // There may be some issues with unrefinement too, but I have not
    // been able to come up with an example (either in my head or in a
    // particular problem) where anything has arisen. I can see that
    // there may be an issue if n_unrefine differs across processes so
    // that (total_n_unrefine > max_keep_unrefined()) on some but not
    // all processes. I haven't seen any examples of this yet so the
    // following code may or may not work!  (Andy, 06/03/08)

    unsigned total_n_unrefine = 0;
#ifdef OOMPH_HAS_MPI
    // Sum n_unrefine across all processors
    if (this->is_mesh_distributed())
    {
      MPI_Allreduce(&n_unrefine,
                    &total_n_unrefine,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    Comm_pt->mpi_comm());
    }
    else
    {
      total_n_unrefine = n_unrefine;
    }
#else
    total_n_unrefine = n_unrefine;
#endif

    oomph_info << "---> " << total_n_refine << " elements to be refined, and "
               << total_n_unrefine << " to be unrefined, in total.\n"
               << std::endl;

    if ((total_n_refine > 0) || (total_n_unrefine > max_keep_unrefined()))
    {
#ifdef PARANOID
#ifdef OOMPH_HAS_MPI


      // Sanity check: Each processor checks if the enforced unrefinement of
      // its haloed element is matched by enforced unrefinement of the
      // corresponding halo elements on the other processors.
      if (this->is_mesh_distributed())
      {
        // Store number of processors and current process
        MPI_Status status;
        int n_proc = Comm_pt->nproc();
        int my_rank = Comm_pt->my_rank();

        // Loop over processes: Each processor sends unrefinement pattern
        // for halo elements with processor d to processor d where it's
        // compared against the unrefinement pattern for the corresponding
        // haloed elements
        for (int d = 0; d < n_proc; d++)
        {
          // No halo with self: Send unrefinement info to proc d
          if (d != my_rank)
          {
            // Get the vector of halo elements whose non-halo counterpart
            // are on processor d
            Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(d));

            // Create vector containing (0)1 to indicate that
            // halo element is (not) to be unrefined
            unsigned nhalo = halo_elem_pt.size();
            Vector<int> halo_to_be_unrefined(nhalo, 0);
            for (unsigned e = 0; e < nhalo; e++)
            {
              if (dynamic_cast<RefineableElement*>(halo_elem_pt[e])
                    ->sons_to_be_unrefined())
              {
                halo_to_be_unrefined[e] = 1;
              }
            }

            // Trap the case when there are no halo elements
            // so that we don't get a segfault in the MPI send
            if (nhalo > 0)
            {
              // Send it across to proc d
              MPI_Send(&halo_to_be_unrefined[0],
                       nhalo,
                       MPI_INT,
                       d,
                       0,
                       Comm_pt->mpi_comm());
            }
          }
          // else (d=my_rank): Receive unrefinement pattern from all
          // other processors (dd)
          else
          {
            // Loop over other processors
            for (int dd = 0; dd < n_proc; dd++)
            {
              // No halo with yourself
              if (dd != d)
              {
                // Get the vector of haloed elements on current processor
                // with processor dd
                Vector<GeneralisedElement*> haloed_elem_pt(
                  this->haloed_element_pt(dd));

                // Ask processor dd to send vector containing (0)1 for
                // halo element with current processor to be (not)unrefined
                unsigned nhaloed = haloed_elem_pt.size();
                Vector<int> halo_to_be_unrefined(nhaloed);
                // Trap to catch the case that there are no haloed elements
                if (nhaloed > 0)
                {
                  // Receive unrefinement pattern of haloes from proc dd
                  MPI_Recv(&halo_to_be_unrefined[0],
                           nhaloed,
                           MPI_INT,
                           dd,
                           0,
                           Comm_pt->mpi_comm(),
                           &status);
                }

                // Check it
                for (unsigned e = 0; e < nhaloed; e++)
                {
                  if (((halo_to_be_unrefined[e] == 0) &&
                       (dynamic_cast<RefineableElement*>(haloed_elem_pt[e])
                          ->sons_to_be_unrefined())) ||
                      ((halo_to_be_unrefined[e] == 1) &&
                       (!dynamic_cast<RefineableElement*>(haloed_elem_pt[e])
                           ->sons_to_be_unrefined())))
                  {
                    std::ostringstream error_message;
                    error_message
                      << "Error in unrefinement: \n"
                      << "Haloed element: " << e << " on proc " << my_rank
                      << " \n"
                      << "wants to be unrefined whereas its halo counterpart "
                         "on\n"
                      << "proc " << dd << " doesn't (or vice versa)...\n"
                      << "This is most likely because the error estimator\n"
                      << "has not assigned the same errors to halo and haloed\n"
                      << "elements -- it ought to!\n";
                    throw OomphLibError(error_message.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
                }
              }
            }
          }
        }


        // Loop over processes: Each processor sends refinement pattern
        // for halo elements with processor d to processor d where it's
        // compared against the refinement pattern for the corresponding
        // haloed elements
        for (int d = 0; d < n_proc; d++)
        {
          // No halo with self: Send refinement info to proc d
          if (d != my_rank)
          {
            // Get the vector of halo elements whose non-halo counterpart
            // are on processor d
            Vector<GeneralisedElement*> halo_elem_pt(this->halo_element_pt(d));

            // Create vector containing (0)1 to indicate that
            // halo element is (not) to be refined
            unsigned nhalo = halo_elem_pt.size();
            Vector<int> halo_to_be_refined(nhalo, 0);
            for (unsigned e = 0; e < nhalo; e++)
            {
              if (dynamic_cast<RefineableElement*>(halo_elem_pt[e])
                    ->to_be_refined())
              {
                halo_to_be_refined[e] = 1;
              }
            }

            // Trap the case when there are no halo elements
            // so that we don't get a segfault in the MPI send
            if (nhalo > 0)
            {
              // Send it across to proc d
              MPI_Send(&halo_to_be_refined[0],
                       nhalo,
                       MPI_INT,
                       d,
                       0,
                       Comm_pt->mpi_comm());
            }
          }
          // else (d=my_rank): Receive refinement pattern from all
          // other processors (dd)
          else
          {
            // Loop over other processors
            for (int dd = 0; dd < n_proc; dd++)
            {
              // No halo with yourself
              if (dd != d)
              {
                // Get the vector of haloed elements on current processor
                // with processor dd
                Vector<GeneralisedElement*> haloed_elem_pt(
                  this->haloed_element_pt(dd));

                // Ask processor dd to send vector containing (0)1 for
                // halo element with current processor to be (not)refined
                unsigned nhaloed = haloed_elem_pt.size();
                Vector<int> halo_to_be_refined(nhaloed);
                // Trap to catch the case that there are no haloed elements
                if (nhaloed > 0)
                {
                  // Receive unrefinement pattern of haloes from proc dd
                  MPI_Recv(&halo_to_be_refined[0],
                           nhaloed,
                           MPI_INT,
                           dd,
                           0,
                           Comm_pt->mpi_comm(),
                           &status);
                }

                // Check it
                for (unsigned e = 0; e < nhaloed; e++)
                {
                  if (((halo_to_be_refined[e] == 0) &&
                       (dynamic_cast<RefineableElement*>(haloed_elem_pt[e])
                          ->to_be_refined())) ||
                      ((halo_to_be_refined[e] == 1) &&
                       (!dynamic_cast<RefineableElement*>(haloed_elem_pt[e])
                           ->to_be_refined())))
                  {
                    std::ostringstream error_message;
                    error_message
                      << "Error in refinement: \n"
                      << "Haloed element: " << e << " on proc " << my_rank
                      << " \n"
                      << "wants to be refined whereas its halo counterpart on\n"
                      << "proc " << dd << " doesn't (or vice versa)...\n"
                      << "This is most likely because the error estimator\n"
                      << "has not assigned the same errors to halo and haloed\n"
                      << "elements -- it ought to!\n";
                    throw OomphLibError(error_message.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
                }
              }
            }
          }
        }
      }

#endif
#endif

      // Perform the actual adaptation
      adapt_mesh(local_doc_info);

      // The number of refineable elements is still local to each process
      Nunrefined = n_unrefine;
      Nrefined = n_refine;
    }
    // If not worthwhile, say so but still reorder nodes and kill
    // external storage for consistency in parallel computations
    else
    {
#ifdef OOMPH_HAS_MPI
      // Delete any external element storage - any interaction will still
      // be set up on the fly again, so we need to get rid of old information.
      // This particularly causes problems in multi-domain examples where
      // we decide not to refine one of the meshes
      this->delete_all_external_storage();
#endif

      // Reorder the nodes within the mesh's node vector
      // to establish a standard ordering regardless of the sequence
      // of mesh refinements -- this is required to allow dump/restart
      // on refined meshes
      this->reorder_nodes();

#ifdef OOMPH_HAS_MPI

      // Now (re-)classify halo and haloed nodes and synchronise hanging
      // nodes
      // This is required in cases where delete_all_external_storage()
      // made dependent nodes to external halo nodes nonhanging.
      if (this->is_mesh_distributed())
      {
        DocInfo doc_info;
        doc_info.disable_doc();
        classify_halo_and_haloed_nodes(doc_info, doc_info.is_doc_enabled());
      }

#endif

      if (n_refine == 0)
      {
        oomph_info << " Not enough benefit in adapting mesh." << std::endl
                   << std::endl;
      }
      Nunrefined = 0;
      Nrefined = 0;
    }
  }

  //========================================================================
  /// Get max/min refinement level
  //========================================================================
  void TreeBasedRefineableMeshBase::get_refinement_levels(
    unsigned& min_refinement_level, unsigned& max_refinement_level)
  {
    // Initialise
    min_refinement_level = UINT_MAX;
    max_refinement_level = 0;

    // Loop over all elements
    unsigned long n_element = this->nelement();
    if (n_element == 0)
    {
      min_refinement_level = 0;
      max_refinement_level = 0;
    }
    else
    {
      for (unsigned long e = 0; e < n_element; e++)
      {
        // Get the refinement level of the element
        unsigned level = dynamic_cast<RefineableElement*>(this->element_pt(e))
                           ->refinement_level();

        if (level > max_refinement_level) max_refinement_level = level;
        if (level < min_refinement_level) min_refinement_level = level;
      }
    }
  }


  //================================================================
  /// Adapt mesh, which exists in two representations,
  /// namely as:
  ///  - a FE mesh
  ///  - a forest of Oc or QuadTrees
  ///
  /// Refinement/derefinement process is documented (in tecplot-able form)
  /// if requested.
  ///
  /// Procedure:
  /// - Loop over all elements and do the refinement for those who want to
  ///   be refined. Note: Refinement/splitting only allocates elements but
  ///   doesn't build them.
  /// - Build the new elements (i.e. give them nodes (create new ones where
  ///   necessary), assign boundary conditions, and add nodes to mesh
  ///   and mesh boundaries.
  /// - For all nodes that were hanging on the previous mesh (and are still
  ///   marked as such), fill in their nodal values (consistent
  ///   with the current hanging node scheme) to make sure they are fully
  ///   functional, should they have become non-hanging during the
  ///   mesh-adaptation. Then mark the nodes as non-hanging.
  /// - Unrefine selected elements (which may cause nodes to be re-built).
  /// - Add the new elements to the mesh (by completely overwriting
  ///   the old Vector of elements).
  /// - Delete any nodes that have become obsolete.
  /// - Mark up hanging nodes and setup hanging node scheme (incl.
  ///   recursive cleanup for hanging nodes that depend on other
  ///   hanging nodes).
  /// - Adjust position of hanging nodes to make sure their position
  ///   is consistent with the FE-based represenetation of their larger
  ///   neighbours.
  /// - run a quick self-test on the neighbour finding scheme and
  ///   check the integrity of the elements (if PARANOID)
  /// - doc hanging node status, boundary conditions, neighbour
  ///   scheme if requested.
  ///
  ///
  /// After adaptation, all nodes (whether new or old) have up-to-date
  /// current and previous values.
  ///
  /// If refinement process is being documented, the following information
  /// is documented:
  /// - The files
  ///   - "neighbours.dat"
  ///   - "all_nodes.dat"
  ///   - "new_nodes.dat"
  ///   - "hang_nodes_*.dat"
  ///     where the * denotes a direction (n,s,e,w) in 2D
  ///     or (r,l,u,d,f,b) in 3D
  ///.
  ///   can be viewed with
  ///   - QHangingNodes.mcr
  ///   .
  /// - The file
  ///    - "hangnodes_withmasters.dat"
  ///    .
  ///    can be viewed with
  ///    - QHangingNodesWithMasters.mcr
  ///    .
  ///    to check the hanging node status.
  /// - The neighbour status of the elements is documented in
  ///   - "neighbours.dat"
  ///   .
  ///   and can be viewed with
  ///   - QuadTreeNeighbours.mcr
  ///   .
  //=================================================================
  void TreeBasedRefineableMeshBase::adapt_mesh(DocInfo& doc_info)
  {
#ifdef OOMPH_HAS_MPI
    // Delete any external element storage before performing the adaptation
    // (in particular, external halo nodes that are on mesh boundaries)
    this->delete_all_external_storage();
#endif

    // Only perform the adapt step if the mesh has any elements.  This is
    // relevant in a distributed problem with multiple meshes, where a
    // particular process may not have any elements on a particular submesh.
    if (this->nelement() > 0)
    {
      // Pointer to mesh needs to be passed to some functions
      Mesh* mesh_pt = this;

      double t_start = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_start = TimingHelpers::timer();
      }

      // Do refinement(=splitting) of elements that have been selected
      // This function encapsulates the template parameter
      this->split_elements_if_required();


      if (Global_timings::Doc_comprehensive_timings)
      {
        double t_end = TimingHelpers::timer();
        oomph_info << "Time for split_elements_if_required: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Now elements have been created -- build all the leaves
      //-------------------------------------------------------
      // Firstly put all the elements into a vector
      Vector<Tree*> leaf_nodes_pt;
      Forest_pt->stick_leaves_into_vector(leaf_nodes_pt);


      if (Global_timings::Doc_comprehensive_timings)
      {
        double t_end = TimingHelpers::timer();
        oomph_info << "Time for stick_leaves_into_vector: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // If we are documenting the output, create the filename
      std::ostringstream fullname;
      std::ofstream new_nodes_file;
      if (doc_info.is_doc_enabled())
      {
        fullname << doc_info.directory() << "/new_nodes" << doc_info.number()
                 << ".dat";
        new_nodes_file.open(fullname.str().c_str());
      }


      // Build all elements and store vector of pointers to new nodes
      // (Note: build() checks if the element has been built
      // already, i.e. if it's not a new element).
      Vector<Node*> new_node_pt;
      bool was_already_built;
      unsigned long num_tree_nodes = leaf_nodes_pt.size();
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        // Pre-build must be performed before any elements are built
        leaf_nodes_pt[e]->object_pt()->pre_build(mesh_pt, new_node_pt);
      }
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        // Now do the actual build of the new elements
        leaf_nodes_pt[e]->object_pt()->build(
          mesh_pt, new_node_pt, was_already_built, new_nodes_file);
      }


      double t_end = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for building " << num_tree_nodes
                   << " new elements: " << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Close the new nodes files, if it was opened
      if (doc_info.is_doc_enabled())
      {
        new_nodes_file.close();
      }

      // Loop over all nodes in mesh and free the dofs of those that were
      //-----------------------------------------------------------------
      // pinned only because they were hanging nodes. Also update their
      //-----------------------------------------------------------------
      // nodal values so that they contain data that is consistent
      //----------------------------------------------------------
      // with the hanging node representation
      //-------------------------------------
      // (Even if the nodal data isn't actually accessed because the node
      // is still hanging -- we don't know this yet, and this step makes
      // sure that all nodes are fully functional and up-to-date, should
      // they become non-hanging below).
      //
      //
      // However, if we have a fixed mesh and hanging nodes on the boundary
      // become non-hanging they will not necessarily respect the curvilinear
      // boundaries. This can only happen in 3D of course because it is not
      // possible to have hanging nodes on boundaries in 2D.
      // The solution is to store those nodes on the boundaries that are
      // currently hanging and then check to see whether they have changed
      // status at the end of the refinement procedure.
      // If it has changed, then we need to adjust their positions.
      const unsigned n_boundary = this->nboundary();
      const unsigned mesh_dim = this->finite_element_pt(0)->dim();
      Vector<std::set<Node*>> hanging_nodes_on_boundary_pt(n_boundary);

      unsigned long n_node = this->nnode();
      for (unsigned long n = 0; n < n_node; n++)
      {
        // Get the pointer to the node
        Node* nod_pt = this->node_pt(n);

        // Get the number of values in the node
        unsigned n_value = nod_pt->nvalue();

        // We need to find if any of the values are hanging
        bool is_hanging = nod_pt->is_hanging();
        // Loop over the values and find out whether any are hanging
        for (unsigned n = 0; n < n_value; n++)
        {
          is_hanging |= nod_pt->is_hanging(n);
        }

        // If the node is hanging then ...
        if (is_hanging)
        {
          // Unless they are turned into hanging nodes again below
          // (this might or might not happen), fill in all the necessary
          // data to make them 'proper' nodes again.

          // Reconstruct the nodal values/position from the node's
          // hanging node representation
          unsigned nt = nod_pt->ntstorage();
          Vector<double> values(n_value);
          unsigned n_dim = nod_pt->ndim();
          Vector<double> position(n_dim);
          // Loop over all history values
          for (unsigned t = 0; t < nt; t++)
          {
            nod_pt->value(t, values);
            for (unsigned i = 0; i < n_value; i++)
            {
              nod_pt->set_value(t, i, values[i]);
            }
            nod_pt->position(t, position);
            for (unsigned i = 0; i < n_dim; i++)
            {
              nod_pt->x(t, i) = position[i];
            }
          }

          // If it's an algebraic node: Update its previous nodal positions too
          AlgebraicNode* alg_node_pt = dynamic_cast<AlgebraicNode*>(nod_pt);
          if (alg_node_pt != 0)
          {
            bool update_all_time_levels = true;
            alg_node_pt->node_update(update_all_time_levels);
          }


          // If it's a Solid node, update Lagrangian coordinates
          // from its hanging node representation
          SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
          if (solid_node_pt != 0)
          {
            unsigned n_lagrangian = solid_node_pt->nlagrangian();
            for (unsigned i = 0; i < n_lagrangian; i++)
            {
              solid_node_pt->xi(i) = solid_node_pt->lagrangian_position(i);
            }
          }

          // Now store geometrically hanging nodes on boundaries that
          // may need updating after refinement.
          // There will only be a problem if we have 3 spatial dimensions
          if ((mesh_dim > 2) && (nod_pt->is_hanging()))
          {
            // If the node is on a boundary then add a pointer to the node
            // to our lookup scheme
            if (nod_pt->is_on_boundary())
            {
              // Storage for the boundaries on which the Node is located
              std::set<unsigned>* boundaries_pt;
              nod_pt->get_boundaries_pt(boundaries_pt);
              if (boundaries_pt != 0)
              {
                // Loop over the boundaries and add a pointer to the node
                // to the appropriate storage scheme
                for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                     it != boundaries_pt->end();
                     ++it)
                {
                  hanging_nodes_on_boundary_pt[*it].insert(nod_pt);
                }
              }
            }
          }

        } // End of is_hanging

        // Initially mark all nodes as 'non-hanging' and `obsolete'
        nod_pt->set_nonhanging();
        nod_pt->set_obsolete();
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for sorting out initial hanging status: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Unrefine all the selected elements: This needs to be
      //-----------------------------------------------------
      // all elements, because the father elements are not actually leaves.
      //-------------------------------------------------------------------

      // Unrefine
      for (unsigned long e = 0; e < Forest_pt->ntree(); e++)
      {
        Forest_pt->tree_pt(e)->traverse_all(&Tree::merge_sons_if_required,
                                            mesh_pt);
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for unrefinement: " << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Add the newly created elements to mesh
      //---------------------------------------

      // Stick all elements into a new vector
      //(note the leaves may have changed, so this is not duplicated work)
      Vector<Tree*> tree_nodes_pt;
      Forest_pt->stick_leaves_into_vector(tree_nodes_pt);

      // Copy the elements into the mesh Vector
      num_tree_nodes = tree_nodes_pt.size();
      Element_pt.resize(num_tree_nodes);
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        Element_pt[e] = tree_nodes_pt[e]->object_pt();

        // Now loop over all nodes in element and mark them as non-obsolete
        // Logic: Initially all nodes in the unrefined mesh were labeled
        // as deleteable. Then we create new elements (whose newly created
        // nodes are obviously non-obsolete), and killed some other elements (by
        // by deleting them and marking the nodes that were not shared by
        // their father as obsolete. Now we loop over all the remaining
        // elements and (re-)label all their nodes as non-obsolete. This
        // saves some nodes that were regarded as obsolete by deleted
        // elements but are still required in some surviving ones
        // from a tragic early death...
        FiniteElement* this_el_pt = this->finite_element_pt(e);
        unsigned n_node = this_el_pt->nnode(); // caching pre-loop
        for (unsigned n = 0; n < n_node; n++)
        {
          this_el_pt->node_pt(n)->set_non_obsolete();
        }
      }


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for adding elements to mesh: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Cannot delete nodes that are still marked as obsolete
      // because they may still be required to assemble the hanging schemes
      //-------------------------------------------------------------------

      // Mark up hanging nodes
      //----------------------

      // Output streams for the hanging nodes
      Vector<std::ofstream*> hanging_output_files;
      // Setup the output files for hanging nodes, this must be called
      // precisely once for the forest. Note that the files will only
      // actually be opened if doc_info.is_doc_enabled() is true
      Forest_pt->open_hanging_node_files(doc_info, hanging_output_files);

      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        // Generic setup
        tree_nodes_pt[e]->object_pt()->setup_hanging_nodes(
          hanging_output_files);
        // Element specific setup
        tree_nodes_pt[e]->object_pt()->further_setup_hanging_nodes();
      }

      // Close the hanging node files and delete the memory allocated
      // for the streams
      Forest_pt->close_hanging_node_files(doc_info, hanging_output_files);


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for setup_hanging_nodes() and "
                      "further_setup_hanging_nodes() for "
                   << num_tree_nodes << " elements: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Read out the number of continously interpolated values
      // from one of the elements (assuming it's the same in all elements)
      unsigned ncont_interpolated_values =
        tree_nodes_pt[0]->object_pt()->ncont_interpolated_values();

      // Complete the hanging nodes schemes by dealing with the
      // recursively hanging nodes
      complete_hanging_nodes(ncont_interpolated_values);

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for complete_hanging_nodes: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      /// Update the boundary element info -- this can be a costly procedure
      /// and for this reason the mesh writer might have decided not to
      /// set up this scheme. If so, we won't change this and suppress
      /// its creation...
      if (Lookup_for_elements_next_boundary_is_setup)
      {
        this->setup_boundary_element_info();
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for boundary element info: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

#ifdef PARANOID

      // Doc/check the neighbours
      //-------------------------
      Vector<Tree*> all_tree_nodes_pt;
      Forest_pt->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

      // Check the neighbours
      Forest_pt->check_all_neighbours(doc_info);

      // Check the integrity of the elements
      // -----------------------------------

      // Loop over elements and get the elemental integrity
      double max_error = 0.0;
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        double max_el_error;
        tree_nodes_pt[e]->object_pt()->check_integrity(max_el_error);
        // If the elemental error is greater than our maximum error
        // reset the maximum
        if (max_el_error > max_error)
        {
          max_error = max_el_error;
        }
      }

      if (max_error > RefineableElement::max_integrity_tolerance())
      {
        std::ostringstream error_stream;
        error_stream << "Mesh refined: Max. error in integrity check: "
                     << max_error << " is too big\n";
        error_stream
          << "i.e. bigger than RefineableElement::max_integrity_tolerance()="
          << RefineableElement::max_integrity_tolerance() << "\n"
          << std::endl;

        std::ofstream some_file;
        some_file.open("ProblemMesh.dat");
        for (unsigned long n = 0; n < n_node; n++)
        {
          // Get the pointer to the node
          Node* nod_pt = this->node_pt(n);
          // Get the dimension
          unsigned n_dim = nod_pt->ndim();
          // Output the coordinates
          for (unsigned i = 0; i < n_dim; i++)
          {
            some_file << this->node_pt(n)->x(i) << " ";
          }
          some_file << std::endl;
        }
        some_file.close();

        error_stream << "Doced problem mesh in ProblemMesh.dat" << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "Mesh refined: Max. error in integrity check: "
                   << max_error << " is OK" << std::endl;
        oomph_info
          << "i.e. less than RefineableElement::max_integrity_tolerance()="
          << RefineableElement::max_integrity_tolerance() << "\n"
          << std::endl;
      }


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for (paranoid only) checking of integrity: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

#endif

      // Loop over all elements other than the final level and deactivate the
      // objects, essentially set the pointer that point to nodes that are
      // about to be deleted to NULL. This must take place here because nodes
      // addressed by elements that are dead but still living in the tree might
      // have been made obsolete in the last round of refinement
      for (unsigned long e = 0; e < Forest_pt->ntree(); e++)
      {
        Forest_pt->tree_pt(e)->traverse_all_but_leaves(
          &Tree::deactivate_object);
      }

      // Now we can prune the dead nodes from the mesh.
      Vector<Node*> deleted_node_pt = this->prune_dead_nodes();

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for deactivating objects and pruning nodes: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Finally: Reorder the nodes within the mesh's node vector
      // to establish a standard ordering regardless of the sequence
      // of mesh refinements -- this is required to allow dump/restart
      // on refined meshes
      this->reorder_nodes();

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for reordering " << nnode()
                   << " nodes: " << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Now we can correct the nodes on boundaries that were hanging that
      // are no longer hanging
      // Only bother if we have more than two dimensions
      if (mesh_dim > 2)
      {
        // Loop over the boundaries
        for (unsigned b = 0; b < n_boundary; b++)
        {
          // Remove deleted nodes from the set
          unsigned n_del = deleted_node_pt.size();
          for (unsigned j = 0; j < n_del; j++)
          {
            hanging_nodes_on_boundary_pt[b].erase(deleted_node_pt[j]);
          }

          // If the nodes that were hanging are still hanging then remove them
          // from the set (note increment is not in for command for efficiencty)
          for (std::set<Node*>::iterator it =
                 hanging_nodes_on_boundary_pt[b].begin();
               it != hanging_nodes_on_boundary_pt[b].end();)
          {
            if ((*it)->is_hanging())
            {
              hanging_nodes_on_boundary_pt[b].erase(it++);
            }
            else
            {
              ++it;
            }
          }

          // Are there any nodes that have changed geometric hanging status
          // on the boundary
          // The slightly painful part is that we must adjust the position
          // via the macro-elements which are only available through the
          // elements and not the nodes.
          if (hanging_nodes_on_boundary_pt[b].size() > 0)
          {
            // If so we loop over all elements adjacent to the boundary
            unsigned n_boundary_element = this->nboundary_element(b);
            for (unsigned e = 0; e < n_boundary_element; ++e)
            {
              // Get a pointer to the element
              FiniteElement* el_pt = this->boundary_element_pt(b, e);

              // Do we have a solid element
              SolidFiniteElement* solid_el_pt =
                dynamic_cast<SolidFiniteElement*>(el_pt);

              // Determine whether there is a macro element
              bool macro_present = (el_pt->macro_elem_pt() != 0);
              // Or a solid macro element
              if (solid_el_pt != 0)
              {
                macro_present |= (solid_el_pt->undeformed_macro_elem_pt() != 0);
              }

              // Only bother to do anything if there is a macro element
              // or undeformed macro element in a SolidElement
              if (macro_present)
              {
                // Loop over the nodes
                // ALH: (could optimise to only loop over
                // node associated with the boundary with more effort)
                unsigned n_el_node = el_pt->nnode();
                for (unsigned n = 0; n < n_el_node; n++)
                {
                  // Cache pointer to the node
                  Node* nod_pt = el_pt->node_pt(n);
                  if (nod_pt->is_on_boundary(b))
                  {
                    // Is the Node in our set
                    std::set<Node*>::iterator it =
                      hanging_nodes_on_boundary_pt[b].find(nod_pt);

                    // If we have found the Node then update the position
                    // to be consistent with the macro-element representation
                    if (it != hanging_nodes_on_boundary_pt[b].end())
                    {
                      // Specialise local and global coordinates to 3D
                      // because there is only a problem in 3D.
                      Vector<double> s(3), x(3);

                      // Find the local coordinate of the ndoe
                      el_pt->local_coordinate_of_node(n, s);

                      // Find the number of time history values
                      const unsigned ntstorage = nod_pt->ntstorage();

                      // Do we have a solid node
                      SolidNode* solid_node_pt =
                        dynamic_cast<SolidNode*>(nod_pt);
                      if (solid_node_pt)
                      {
                        // Assign Lagrangian coordinates from undeformed
                        // macro element (if it has one -- get_x_and_xi()
                        // does "the right thing" anyway. Leave actual
                        // nodal positions alone -- we're doing a solid
                        // mechanics problem and once we're going
                        // the nodal positions are always computed, never
                        // (automatically) reset to macro-element based
                        // positions; not even on pinned boundaries
                        // because the user may have other ideas about where
                        // these should go -- pinning means "don't touch the
                        // value", not "leave where the macro-element thinks
                        // it should be"
                        Vector<double> x_fe(3), xi(3), xi_fe(3);
                        solid_el_pt->get_x_and_xi(s, x_fe, x, xi_fe, xi);
                        for (unsigned i = 0; i < 3; i++)
                        {
                          solid_node_pt->xi(i) = xi[i];
                        }
                      }
                      else
                      {
                        // Set position and history values from the
                        // macro-element representation
                        for (unsigned t = 0; t < ntstorage; t++)
                        {
                          // Get the history value from the macro element
                          el_pt->get_x(t, s, x);

                          // Set the coordinate to that of the macroelement
                          // representation
                          for (unsigned i = 0; i < 3; i++)
                          {
                            nod_pt->x(t, i) = x[i];
                          }
                        }
                      } // End of non-solid node case

                      // Now remove the node from the list
                      hanging_nodes_on_boundary_pt[b].erase(it);

                      // If there are no Nodes left then exit the loops
                      if (hanging_nodes_on_boundary_pt[b].size() == 0)
                      {
                        e = n_boundary_element;
                        break;
                      }
                    }
                  }
                }
              } // End of macro element case
            }
          }
        }
      } // End of case when we have fixed nodal positions


      // Final doc
      //-----------
      if (doc_info.is_doc_enabled())
      {
        // Doc the boundary conditions ('0' for non-existent, '1' for free,
        //----------------------------------------------------------------
        // '2' for pinned -- ideal for tecplot scatter sizing.
        //----------------------------------------------------
        // num_tree_nodes=tree_nodes_pt.size();

        // Determine maximum number of values at any node in this type of
        // element
        RefineableElement* el_pt = tree_nodes_pt[0]->object_pt();
        // Initalise max_nval
        unsigned max_nval = 0;
        for (unsigned n = 0; n < el_pt->nnode(); n++)
        {
          if (el_pt->node_pt(n)->nvalue() > max_nval)
          {
            max_nval = el_pt->node_pt(n)->nvalue();
          }
        }

        // Open the output file
        std::ofstream bcs_file;
        fullname.str("");
        fullname << doc_info.directory() << "/bcs" << doc_info.number()
                 << ".dat";
        bcs_file.open(fullname.str().c_str());

        // Loop over elements
        for (unsigned long e = 0; e < num_tree_nodes; e++)
        {
          el_pt = tree_nodes_pt[e]->object_pt();
          // Loop over nodes in element
          unsigned n_nod = el_pt->nnode();
          for (unsigned n = 0; n < n_nod; n++)
          {
            // Get pointer to the node
            Node* nod_pt = el_pt->node_pt(n);
            // Find the dimension of the node
            unsigned n_dim = nod_pt->ndim();
            // Write the nodal coordinates to the file
            for (unsigned i = 0; i < n_dim; i++)
            {
              bcs_file << nod_pt->x(i) << " ";
            }

            // Loop over all values in this element
            for (unsigned i = 0; i < max_nval; i++)
            {
              // Value exists at this node:
              if (i < nod_pt->nvalue())
              {
                bcs_file << " " << 1 + nod_pt->is_pinned(i);
              }
              // ...if not just dump out a zero
              else
              {
                bcs_file << " 0 ";
              }
            }
            bcs_file << std::endl;
          }
        }
        bcs_file.close();

        // Doc all nodes
        //---------------
        std::ofstream all_nodes_file;
        fullname.str("");
        fullname << doc_info.directory() << "/all_nodes" << doc_info.number()
                 << ".dat";
        all_nodes_file.open(fullname.str().c_str());

        all_nodes_file << "ZONE \n";

        // Need to recompute the number of nodes since it may have
        // changed during mesh refinement/unrefinement
        n_node = this->nnode();
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);
          unsigned n_dim = nod_pt->ndim();
          for (unsigned i = 0; i < n_dim; i++)
          {
            all_nodes_file << this->node_pt(n)->x(i) << " ";
          }
          all_nodes_file << std::endl;
        }

        all_nodes_file.close();


        // Doc all hanging nodes:
        //-----------------------
        std::ofstream some_file;
        fullname.str("");
        fullname << doc_info.directory() << "/all_hangnodes"
                 << doc_info.number() << ".dat";
        some_file.open(fullname.str().c_str());
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);

          if (nod_pt->is_hanging())
          {
            unsigned n_dim = nod_pt->ndim();
            for (unsigned i = 0; i < n_dim; i++)
            {
              some_file << nod_pt->x(i) << " ";
            }

            // ALH: Added this to stop Solid problems seg-faulting
            if (this->node_pt(n)->nvalue() > 0)
            {
              some_file << " " << nod_pt->raw_value(0);
            }
            some_file << std::endl;
          }
        }
        some_file.close();

        // Doc all hanging nodes and their masters
        // View with QHangingNodesWithMasters.mcr
        fullname.str("");
        fullname << doc_info.directory() << "/geometric_hangnodes_withmasters"
                 << doc_info.number() << ".dat";
        some_file.open(fullname.str().c_str());
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);
          if (nod_pt->is_hanging())
          {
            unsigned n_dim = nod_pt->ndim();
            unsigned nmaster = nod_pt->hanging_pt()->nmaster();
            some_file << "ZONE I=" << nmaster + 1 << std::endl;
            for (unsigned i = 0; i < n_dim; i++)
            {
              some_file << nod_pt->x(i) << " ";
            }
            some_file << " 2 " << std::endl;

            for (unsigned imaster = 0; imaster < nmaster; imaster++)
            {
              Node* master_nod_pt =
                nod_pt->hanging_pt()->master_node_pt(imaster);
              unsigned n_dim = master_nod_pt->ndim();
              for (unsigned i = 0; i < n_dim; i++)
              {
                some_file << master_nod_pt->x(i) << " ";
              }
              some_file << " 1 " << std::endl;
            }
          }
        }
        some_file.close();

        // Doc all hanging nodes and their masters
        // View with QHangingNodesWithMasters.mcr
        for (unsigned i = 0; i < ncont_interpolated_values; i++)
        {
          fullname.str("");
          fullname << doc_info.directory()
                   << "/nonstandard_hangnodes_withmasters" << i << "_"
                   << doc_info.number() << ".dat";
          some_file.open(fullname.str().c_str());
          unsigned n_nod = this->nnode();
          for (unsigned long n = 0; n < n_nod; n++)
          {
            Node* nod_pt = this->node_pt(n);
            if (nod_pt->is_hanging(i))
            {
              if (nod_pt->hanging_pt(i) != nod_pt->hanging_pt())
              {
                unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                some_file << "ZONE I=" << nmaster + 1 << std::endl;
                unsigned n_dim = nod_pt->ndim();
                for (unsigned j = 0; j < n_dim; j++)
                {
                  some_file << nod_pt->x(j) << " ";
                }
                some_file << " 2 " << std::endl;
                for (unsigned imaster = 0; imaster < nmaster; imaster++)
                {
                  Node* master_nod_pt =
                    nod_pt->hanging_pt(i)->master_node_pt(imaster);
                  unsigned n_dim = master_nod_pt->ndim();
                  for (unsigned j = 0; j < n_dim; j++)
                  {
                    //               some_file << master_nod_pt->x(i) << " ";
                  }
                  some_file << " 1 " << std::endl;
                }
              }
            }
          }
          some_file.close();
        }
      } // End of documentation
    } // End if (this->nelement()>0)


#ifdef OOMPH_HAS_MPI

    // Now (re-)classify halo and haloed nodes and synchronise hanging
    // nodes
    if (this->is_mesh_distributed())
    {
      classify_halo_and_haloed_nodes(doc_info, doc_info.is_doc_enabled());
    }

#endif
  }


  //========================================================================
  /// Refine mesh uniformly
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_uniformly(DocInfo& doc_info)
  {
    // Select all elements for refinement
    unsigned long Nelement = this->nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      dynamic_cast<RefineableElement*>(this->element_pt(e))
        ->select_for_refinement();
    }

    // Do the actual mesh adaptation
    adapt_mesh(doc_info);
  }


  //========================================================================
  /// p-refine mesh uniformly
  //========================================================================
  void TreeBasedRefineableMeshBase::p_refine_uniformly(DocInfo& doc_info)
  {
    // Select all elements for refinement
    unsigned long Nelement = this->nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // Get pointer to p-refineable element
      PRefineableElement* el_pt =
        dynamic_cast<PRefineableElement*>(this->element_pt(e));
      // Mark for p-refinement if possible. If not then p_adapt_mesh() will
      // report the error.
      if (el_pt != 0)
      {
        el_pt->select_for_p_refinement();
      }
    }

    // Do the actual mesh adaptation
    p_adapt_mesh(doc_info);
  }


  //========================================================================
  /// Refine mesh by splitting the elements identified
  /// by their numbers.
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_selected_elements(
    const Vector<unsigned>& elements_to_be_refined)
  {
#ifdef OOMPH_HAS_MPI
    if (this->is_mesh_distributed())
    {
      std::ostringstream warn_stream;
      warn_stream << "You are attempting to refine selected elements of a "
                  << std::endl
                  << "distributed mesh. This may have undesired effects."
                  << std::endl;

      OomphLibWarning(warn_stream.str(),
                      "TreeBasedRefineableMeshBase::refine_selected_elements()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Select elements for refinement
    unsigned long nref = elements_to_be_refined.size();
    for (unsigned long e = 0; e < nref; e++)
    {
      dynamic_cast<RefineableElement*>(
        this->element_pt(elements_to_be_refined[e]))
        ->select_for_refinement();
    }

    // Do the actual mesh adaptation
    adapt_mesh();
  }


  //========================================================================
  /// Refine mesh by splitting the elements identified
  /// by their pointers
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_selected_elements(
    const Vector<RefineableElement*>& elements_to_be_refined_pt)
  {
#ifdef OOMPH_HAS_MPI
    if (this->is_mesh_distributed())
    {
      std::ostringstream warn_stream;
      warn_stream << "You are attempting to refine selected elements of a "
                  << std::endl
                  << "distributed mesh. This may have undesired effects."
                  << std::endl;

      OomphLibWarning(warn_stream.str(),
                      "TreeBasedRefineableMeshBase::refine_selected_elements()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Select elements for refinement
    unsigned long nref = elements_to_be_refined_pt.size();
    for (unsigned long e = 0; e < nref; e++)
    {
      elements_to_be_refined_pt[e]->select_for_refinement();
    }

    // Do the actual mesh adaptation
    adapt_mesh();
  }


  //========================================================================
  /// \short Refine to same degree as the reference mesh.
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_base_mesh_as_in_reference_mesh(
    TreeBasedRefineableMeshBase* const& ref_mesh_pt)
  {
    // Assign storage for refinement pattern
    Vector<Vector<unsigned>> to_be_refined;

    // Get refinement pattern of reference mesh:
    ref_mesh_pt->get_refinement_pattern(to_be_refined);

    // Refine mesh according to given refinement pattern
    refine_base_mesh(to_be_refined);
  }


  //========================================================================
  /// \short Refine to same degree as the reference mesh minus one. Useful
  /// function for multigrid solvers; allows the easy copy of a mesh
  /// to the level of refinement just below the current one. Returns
  /// a boolean variable which indicates if the reference mesh has not
  /// been refined at all
  //========================================================================
  bool TreeBasedRefineableMeshBase::
    refine_base_mesh_as_in_reference_mesh_minus_one(
      TreeBasedRefineableMeshBase* const& ref_mesh_pt)
  {
    // Assign storage for refinement pattern
    Vector<Vector<unsigned>> to_be_refined;

    // Get refinement pattern of reference mesh:
    ref_mesh_pt->get_refinement_pattern(to_be_refined);

    // Find the length of the vector
    unsigned nrefinement_levels = to_be_refined.size();

    // If the reference mesh has not been refined a single time then
    // we cannot create an unrefined copy so stop here
    if (nrefinement_levels == 0)
    {
      return false;
    }
    // If the reference mesh has been refined at least once
    else
    {
      // Remove the last entry of the vector to make sure we refine to
      // the same level minus one
      to_be_refined.resize(nrefinement_levels - 1);

      // Refine mesh according to given refinement pattern
      refine_base_mesh(to_be_refined);

      // Indicate that it was possible to create an unrefined copy
      return true;
    }
  }


  //========================================================================
  /// Refine mesh once so that its topology etc becomes that of the
  /// (finer!) reference mesh -- if possible! Useful for meshes in multigrid
  /// hierarchies. If the meshes are too different and the conversion
  /// cannot be performed, the code dies (provided PARANOID is enabled).
  //========================================================================
  void TreeBasedRefineableMeshBase::refine_as_in_reference_mesh(
    TreeBasedRefineableMeshBase* const& ref_mesh_pt)
  {
    oomph_info << "WARNING : This has not been checked comprehensively yet"
               << std::endl
               << "Check it and remove this break " << std::endl;
    pause("Yes really pause");

#ifdef PARANOID
    // The max. refinement levels of the two meshes need to differ
    // by one, otherwise what we're doing here doesn't make sense.
    unsigned my_min, my_max;
    get_refinement_levels(my_min, my_max);

    unsigned ref_min, ref_max;
    ref_mesh_pt->get_refinement_levels(ref_min, ref_max);

    if (ref_max != my_max + 1)
    {
      std::ostringstream error_stream;
      error_stream
        << "Meshes definitely don't differ by one refinement level \n"
        << "max. refinement levels: " << ref_max << " " << my_max << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Vector storing the elements of the uniformly unrefined mesh
    Vector<Tree*> coarse_elements_pt;

    // Map storing which father elements have already been added to coarse mesh
    // (Default return is 0).
    std::map<Tree*, unsigned> father_element_included;

    // Extract active elements (=leaf nodes in the quadtree) from reference
    // mesh.
    Vector<Tree*> leaf_nodes_pt;
    ref_mesh_pt->forest_pt()->stick_leaves_into_vector(leaf_nodes_pt);

    // Loop over all elements (in their quadtree impersonation) and
    // check if their fathers's sons are all leaves too:
    unsigned nelem = leaf_nodes_pt.size();
    for (unsigned e = 0; e < nelem; e++)
    {
      // Pointer to leaf node
      Tree* leaf_pt = leaf_nodes_pt[e];

      // Get pointer to father:
      Tree* father_pt = leaf_pt->father_pt();

      // If we don't have a father we're at the root level in which
      // case this element can't be unrefined.
      if (0 == father_pt)
      {
        coarse_elements_pt.push_back(leaf_pt);
      }
      else
      {
        // Loop over the father's sons to check if they're
        // all non-leafs, i.e. if they can be unrefined
        bool can_unrefine = true;
        unsigned n_sons = father_pt->nsons();
        for (unsigned i = 0; i < n_sons; i++)
        {
          // If (at least) one of the sons is not a leaf, we can't unrefine
          if (!father_pt->son_pt(i)->is_leaf()) can_unrefine = false;
        }

        // If we can unrefine, the father element will be
        // an element in the coarse mesh, the sons won't
        if (can_unrefine)
        {
          if (father_element_included[father_pt] == 0)
          {
            coarse_elements_pt.push_back(father_pt);
            father_element_included[father_pt] = 1;
          }
        }
        // Son will still be there on the coarse mesh
        else
        {
          coarse_elements_pt.push_back(leaf_pt);
        }
      }
    }

    // Number of elements in ref mesh if it was unrefined uniformly:
    unsigned nel_coarse = coarse_elements_pt.size();


#ifdef PARANOID
    bool stop_it = false;
    // The numbers had better match otherwise we might as well stop now...
    if (nel_coarse != this->nelement())
    {
      oomph_info << "Number of elements in uniformly unrefined reference mesh: "
                 << nel_coarse << std::endl;
      oomph_info << "Number of elements in 'this' mesh: " << nel_coarse
                 << std::endl;
      oomph_info << "don't match" << std::endl;
      stop_it = true;
    }
#endif

    // Now loop over all elements in uniformly coarsened reference mesh
    // and check if add the number of any element that was created
    // by having had its sons merged to the vector of elements that
    // need to get refined if we go the other way
    Vector<unsigned> elements_to_be_refined;
    for (unsigned i = 0; i < nel_coarse; i++)
    {
      if (father_element_included[coarse_elements_pt[i]] == 1)
      {
        elements_to_be_refined.push_back(i);
      }
    }


#ifdef PARANOID
    // Doc troublesome meshes:
    if (stop_it)
    {
      std::ofstream some_file;
      some_file.open("orig_mesh.dat");
      this->output(some_file);
      some_file.close();
      oomph_info << "Documented original ('this')mesh in orig_mesh.dat"
                 << std::endl;
    }
#endif


    // Now refine precisely these elements in "this" mesh.
    refine_selected_elements(elements_to_be_refined);


#ifdef PARANOID

    // Check if the nodal positions of all element's nodes agree
    // in the two fine meshes:
    double tol = 1.0e-5;
    for (unsigned e = 0; e < nelem; e++)
    {
      // Get elements
      FiniteElement* ref_el_pt = ref_mesh_pt->finite_element_pt(e);
      FiniteElement* el_pt = this->finite_element_pt(e);

      // Loop over nodes
      unsigned nnod = ref_el_pt->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        // Get nodes
        Node* ref_node_pt = ref_el_pt->node_pt(j);
        Node* node_pt = el_pt->node_pt(j);

        // Check error in position
        double error = 0.0;
        unsigned ndim = node_pt->ndim();
        for (unsigned i = 0; i < ndim; i++)
        {
          error += pow(node_pt->x(i) - ref_node_pt->x(i), 2);
        }
        error = sqrt(error);

        if (error > tol)
        {
          oomph_info << "Error in nodal position of node " << j << ": " << error
                     << "           [tol=" << tol << "]" << std::endl;
          stop_it = true;
        }
      }
    }

    // Do we have a death wish?
    if (stop_it)
    {
      // Doc troublesome meshes:
      std::ofstream some_file;
      some_file.open("refined_mesh.dat");
      this->output(some_file);
      some_file.close();

      some_file.open("finer_mesh.dat");
      ref_mesh_pt->output(some_file);
      some_file.close();

      throw OomphLibError(
        "Bailing out. Doced refined_mesh.dat finer_mesh.dat\n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

#endif
  }


  //========================================================================
  /// Unrefine mesh uniformly. Return 0 for success,
  /// 1 for failure (if unrefinement has reached the coarsest permitted
  /// level)
  //========================================================================
  unsigned TreeBasedRefineableMeshBase::unrefine_uniformly()
  {
    // We can't just select all elements for unrefinement
    // because they need to merge with their brothers.
    // --> Rather than repeating the convoluted logic of
    // RefineableQuadMesh<ELEMENT>::adapt(Vector<double>& elemental_error)
    // here (code duplication!) hack it by filling the error
    // vector with values that ensure unrefinement for all
    // elements where this is possible

    // Create dummy vector for elemental errors
    unsigned long Nelement = this->nelement();
    Vector<double> elemental_error(Nelement);

    // In order to force unrefinement, set the min permitted error to
    // be the default and then set the actual error to be below this.
    // This avoids problems when the actual min error is zero (or small)
    // For sanity's sake, also set the max permitted error back to the default
    // so that we have a max error bigger than a min error
    const double current_min_error = this->min_permitted_error();
    const double current_max_error = this->max_permitted_error();

    this->min_permitted_error() = 1.0e-5;
    this->max_permitted_error() = 1.0e-3;

    double error = min_permitted_error() / 100.0;
    for (unsigned long e = 0; e < Nelement; e++)
    {
      elemental_error[e] = error;
    }

    // Temporarily lift any restrictions on the minimum number of
    // elements that need to be unrefined to make it worthwhile
    unsigned backup = max_keep_unrefined();
    max_keep_unrefined() = 0;

    // Do the actual mesh adaptation with fake error vector
    adapt(elemental_error);

    // Reset the minimum number of elements that need to be unrefined
    // to make it worthwhile
    max_keep_unrefined() = backup;

    // Now restore the error tolerances
    this->min_permitted_error() = current_min_error;
    this->max_permitted_error() = current_max_error;

    // Has the unrefinement actually changed anything?
    if (Nelement == this->nelement())
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }

  //==================================================================
  /// Given a node, return a vector of pointers to master nodes and a
  /// vector of the associated weights.
  /// This is done recursively, so if a node is not hanging,
  /// the node is regarded as its own master node which has weight 1.0.
  //==================================================================
  void TreeBasedRefineableMeshBase::complete_hanging_nodes_recursively(
    Node*& nod_pt,
    Vector<Node*>& master_nodes,
    Vector<double>& hang_weights,
    const int& i)
  {
    // Is the node hanging in the variable i
    if (nod_pt->is_hanging(i))
    {
      // Loop over all master nodes
      HangInfo* const hang_pt = nod_pt->hanging_pt(i);
      unsigned nmaster = hang_pt->nmaster();

      for (unsigned m = 0; m < nmaster; m++)
      {
        // Get the master node
        Node* master_nod_pt = hang_pt->master_node_pt(m);

        // Keep in memory the size of the list before adding the nodes this
        // master node depends on. This is required so that the recursion is
        // only performed on these particular master nodes. A master node
        // could contain contributions from two separate pseudo-masters.
        // These contributions must be summed, not multiplied.
        int first_new_node = master_nodes.size();

        // Now check which master nodes this master node depends on
        complete_hanging_nodes_recursively(
          master_nod_pt, master_nodes, hang_weights, i);

        // Multiply old weight by new weight for all the nodes this master
        // node depends on
        unsigned n_new_master_node = master_nodes.size();

        double mtr_weight = hang_pt->master_weight(m);

        for (unsigned k = first_new_node; k < n_new_master_node; k++)
        {
          hang_weights[k] = mtr_weight * hang_weights[k];
        }
      }
    }
    else
    // Node isn't hanging so it enters itself with the full weight
    {
      master_nodes.push_back(nod_pt);
      hang_weights.push_back(1.0);
    }
  }


  //==================================================================
  /// Complete the hanging node scheme recursively.
  /// After the initial markup scheme, hanging nodes
  /// can depend on other hanging nodes ---> AAAAAAAAARGH!
  /// Need to translate this into a scheme where all
  /// hanging  nodes only depend on non-hanging nodes...
  //==================================================================
  void TreeBasedRefineableMeshBase::complete_hanging_nodes(
    const int& ncont_interpolated_values)
  {
    // Number of nodes in mesh
    unsigned long n_node = this->nnode();
    double min_weight = 1.0e-8; // RefineableBrickElement::min_weight_value();

    // Loop over the nodes in the mesh
    for (unsigned long n = 0; n < n_node; n++)
    {
      // Assign a local pointer to the node
      Node* nod_pt = this->node_pt(n);

      // Loop over the values,
      // N.B. geometric hanging data is stored at the index -1
      for (int i = -1; i < ncont_interpolated_values; i++)
      {
        // Is the node hanging?
        if (nod_pt->is_hanging(i))
        {
          // If it is geometric OR has hanging node data that differs
          // from the geometric data, we must do some work
          if ((i == -1) || (nod_pt->hanging_pt(i) != nod_pt->hanging_pt()))
          {
            // Find out the ultimate map of dependencies: Master nodes
            // and associated weights
            Vector<Node*> master_nodes;
            Vector<double> hanging_weights;
            complete_hanging_nodes_recursively(
              nod_pt, master_nodes, hanging_weights, i);

            // put them into a map to merge all the occurences of the same node
            // (add the weights)
            std::map<Node*, double> hang_weights;
            unsigned n_master = master_nodes.size();
            for (unsigned k = 0; k < n_master; k++)
            {
              if (std::fabs(hanging_weights[k]) > min_weight)
                hang_weights[master_nodes[k]] += hanging_weights[k];
            }

            // Create new hanging data (we know how many data there are)
            HangInfo* hang_pt = new HangInfo(hang_weights.size());

            unsigned hang_weights_index = 0;
            // Copy the map into the HangInfo object
            typedef std::map<Node*, double>::iterator IT;
            for (IT it = hang_weights.begin(); it != hang_weights.end(); ++it)
            {
              hang_pt->set_master_node_pt(
                hang_weights_index, it->first, it->second);
              ++hang_weights_index;
            }

            // Assign the new hanging pointer to the appropriate value
            nod_pt->set_hanging_pt(hang_pt, i);
          }
        }
      }
    }

#ifdef PARANOID

    // Check hanging node scheme: The weights need to add up to one
    //-------------------------------------------------------------
    // Loop over all values indices
    for (int i = -1; i < ncont_interpolated_values; i++)
    {
      // Loop over all nodes in mesh
      for (unsigned long n = 0; n < n_node; n++)
      {
        // Set a local pointer to the node
        Node* nod_pt = this->node_pt(n);

        // Is it hanging?
        if (nod_pt->is_hanging(i))
        {
          unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
          double sum = 0.0;
          for (unsigned imaster = 0; imaster < nmaster; imaster++)
          {
            sum += nod_pt->hanging_pt(i)->master_weight(imaster);
          }
          if (std::fabs(sum - 1.0) > 1.0e-7)
          {
            oomph_info << "WARNING: Sum of master node weights fabs(sum-1.0) "
                       << std::fabs(sum - 1.0) << " for node number " << n
                       << " at value " << i << std::endl;
          }
        }
      }
    }
#endif
  }

  // Sorting out the cases of nodes that are hanging on at least one but not
  // all of the processors for which they are part of the local mesh.

#ifdef OOMPH_HAS_MPI

  //========================================================================
  /// Deal with nodes that are hanging on one process but not another
  /// (i.e. the hanging status of the haloed and halo layers disagrees)
  //========================================================================
  void TreeBasedRefineableMeshBase::synchronise_hanging_nodes(
    const unsigned& ncont_interpolated_values)
  {
    // Store number of processors and current process
    MPI_Status status;
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    double t_start = 0.0;
    double t_end = 0.0;

    // Storage for the hanging status of halo/haloed nodes on elements
    Vector<Vector<int>> haloed_hanging(n_proc);
    Vector<Vector<int>> halo_hanging(n_proc);

    // Counter for the number of nodes which require additional synchronisation
    unsigned nnode_still_requiring_synchronisation = 0;

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Store number of continuosly interpolated values as int
    int ncont_inter_values = ncont_interpolated_values;

    // Loop over processes: Each processor checks that is haloed nodes
    // with proc d have consistent hanging stats with halo counterparts.
    for (int d = 0; d < n_proc; d++)
    {
      // No halo with self: Setup hang info for my haloed nodes with proc d
      // then get ready to receive halo info from processor d.
      if (d != my_rank)
      {
        // Loop over haloed nodes
        unsigned nh = nhaloed_node(d);
        for (unsigned j = 0; j < nh; j++)
        {
          // Get node
          Node* nod_pt = haloed_node_pt(d, j);

          // Loop over the hanging status for each interpolated variable
          // (and the geometry)
          for (int icont = -1; icont < ncont_inter_values; icont++)
          {
            // Store the hanging status of this haloed node
            if (nod_pt->is_hanging(icont))
            {
              unsigned n_master = nod_pt->hanging_pt(icont)->nmaster();
              haloed_hanging[d].push_back(n_master);
            }
            else
            {
              haloed_hanging[d].push_back(0);
            }
          }
        }

        // Receive the hanging status information from the corresponding process
        unsigned count_haloed = haloed_hanging[d].size();

#ifdef PARANOID
        // Check that number of halo and haloed data match
        unsigned tmp = 0;
        MPI_Recv(&tmp, 1, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm(), &status);
        if (tmp != count_haloed)
        {
          std::ostringstream error_stream;
          error_stream << "Number of halo data, " << tmp
                       << ", does not match number of haloed data, "
                       << count_haloed << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Get the data (if any)
        if (count_haloed != 0)
        {
          halo_hanging[d].resize(count_haloed);
          MPI_Recv(&halo_hanging[d][0],
                   count_haloed,
                   MPI_INT,
                   d,
                   0,
                   Comm_pt->mpi_comm(),
                   &status);
        }
      }
      else // d==my_rank, i.e. current process: Send halo hanging status
           // to process dd where it's received (see above) and compared
           // and compared against the hang status of the haloed nodes
      {
        for (int dd = 0; dd < n_proc; dd++)
        {
          // No halo with yourself
          if (dd != d)
          {
            // Storage for halo hanging status and counter
            Vector<int> local_halo_hanging;

            // Loop over halo nodes
            unsigned nh = nhalo_node(dd);
            for (unsigned j = 0; j < nh; j++)
            {
              // Get node
              Node* nod_pt = halo_node_pt(dd, j);

              // Loop over the hanging status for each interpolated variable
              // (and the geometry)
              for (int icont = -1; icont < ncont_inter_values; icont++)
              {
                // Store hanging status of halo node
                if (nod_pt->is_hanging(icont))
                {
                  unsigned n_master = nod_pt->hanging_pt(icont)->nmaster();
                  local_halo_hanging.push_back(n_master);
                }
                else
                {
                  local_halo_hanging.push_back(0);
                }
              }
            }


            // Send the information to the relevant process
            unsigned count_halo = local_halo_hanging.size();

#ifdef PARANOID
            // Check that number of halo and haloed data match
            MPI_Send(&count_halo, 1, MPI_UNSIGNED, dd, 0, Comm_pt->mpi_comm());
#endif

            // Send data (if any)
            if (count_halo != 0)
            {
              MPI_Send(&local_halo_hanging[0],
                       count_halo,
                       MPI_INT,
                       dd,
                       0,
                       Comm_pt->mpi_comm());
            }
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for first all-to-all in synchronise_hanging_nodes(): "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }


    // Now compare equivalent halo and haloed vectors to find discrepancies.
    // It is possible that a master node may not be on either process involved
    // in the halo-haloed scheme; to work round this, we use the shared_node
    // storage scheme, which stores all nodes that are on each pair of
    // processors in the same order on each of the two processors

    // Vector to store info that will help us locate master nodes on "other"
    // processor
    Vector<HangHelperStruct> hang_info;

    // Copy vector-based representation of shared nodes into
    // map for faster search
    Vector<std::map<Node*, unsigned>> shared_node_map(n_proc);
    for (int d = 0; d < n_proc; d++)
    {
      unsigned n = Shared_node_pt[d].size();
      for (unsigned jj = 0; jj < n; jj++)
      {
        shared_node_map[d][Shared_node_pt[d][jj]] = jj;
      }
    }


    // Loop over domains: Each processor checks consistency of hang status
    // of its haloed nodes with proc d against the halo counterpart. Haloed
    // wins if there are any discrepancies.
    for (int d = 0; d < n_proc; d++)
    {
      // No halo with yourself
      if (d != my_rank)
      {
        unsigned discrepancy_count = 0;
        unsigned discrepancy_count_buff = 0;

        // Storage for hanging information that needs to be sent to the
        // relevant process if there is a discrepancy in the hanging status
        Vector<int> send_data;
        Vector<double> send_double_data;
        // Buffer storage for data to be sent
        //(We need this because we cannot tell until the end of the loop over
        // a node's master nodes whether we can reconcile its hanging status)
        Vector<int> send_data_buff(0);
        Vector<double> send_double_data_buff(0);

        // Counter for traversing haloed data
        unsigned count = 0;

        // Indicate presence of discrepancy. Default: there's none
        unsigned discrepancy = 0;
        unsigned discrepancy_buff = 0;

        // Loop over haloed nodes
        unsigned nh = nhaloed_node(d);
        for (unsigned j = 0; j < nh; j++)
        {
          // Get node
          Node* nod_pt = haloed_node_pt(d, j);

          // Loop over the hanging status for each interpolated variable
          // (and the geometry)
          for (int icont = -1; icont < ncont_inter_values; icont++)
          {
            // Compare hanging status of halo/haloed counterpart structure

            // Haloed is is hanging and haloed has different number
            // of master nodes (which includes none in which case it isn't
            // hanging)
            if ((haloed_hanging[d][count] > 0) &&
                (haloed_hanging[d][count] != halo_hanging[d][count]))
            {
              discrepancy_buff = 1;
              discrepancy_count_buff++;

              // Flag to check if all masters of this node have been found
              bool found_all_masters = true;

              // Find master nodes of haloed node
              HangInfo* hang_pt = nod_pt->hanging_pt(icont);
              unsigned nhd_master = hang_pt->nmaster();

              // Add the number of master nodes to the hanging_nodes vector
              send_data_buff.push_back(nhd_master);

              // Loop over master nodes required for HangInfo
              for (unsigned m = 0; m < nhd_master; m++)
              {
                // Get mth master node
                Node* master_nod_pt = hang_pt->master_node_pt(m);


                //              //------------------------------------------
                //              // Old direct search is much slower!
                //              // Keep this code alive to allow comparisons
                //
                //              double t_start_search=TimingHelpers::timer();
                //
                //              // This node will be shared: find it!
                //              bool found=false;
                //
                //              // Look in shared node storage with domain d
                //              first unsigned nnod_shared=nshared_node(d); for
                //              (unsigned k=0; k<nnod_shared; k++)
                //               {
                //                if (master_nod_pt==shared_node_pt(d,k))
                //                 {
                //                  // Found a master: Put its number in the
                //                  shared
                //                  // scheme into send data
                //                  send_data.push_back(k);
                //
                //                  // Add the weight
                //                  send_double_data.push_back(hang_pt->master_weight(m));
                //
                //                  // Done
                //                  found=true;
                //                  break;
                //                 }
                //               }
                //
                //              double t_end_search=TimingHelpers::timer();
                //              t_start_search_total+=(t_end_search-t_start_search);
                //
                //              // end direct search demo
                //              //-----------------------


                // This node will be shared: find it!
                bool found = false;

                // Which processor holds the non-halo counterpart of this
                // node?
                int non_halo_proc_id = master_nod_pt->non_halo_proc_ID();

                // Try to find node in map with proc d and get iterator to entry
                std::map<Node*, unsigned>::iterator it =
                  shared_node_map[d].find(master_nod_pt);

                // If it's not in there iterator points to end of
                // set
                if (it != shared_node_map[d].end())
                {
                  // Found a master: When looking up the node in the shared
                  // node scheme on processor d, which processor do I work
                  // with? The current one
                  send_data_buff.push_back(my_rank);

                  // Found a master: Send its index in the shared
                  // node scheme
                  send_data_buff.push_back((*it).second);

                  // Add the weight
                  send_double_data_buff.push_back(hang_pt->master_weight(m));

                  // Done
                  found = true;
                }

                // If we haven't found it in the shared node scheme with proc d
                // find it in the shared node scheme with the processor that
                // holds the non-halo version
                if (!found)
                {
                  // This is odd -- can't currently handle the case where
                  // node is owned by current processor (indicated by
                  // non_halo_proc_id being negative
                  if (non_halo_proc_id < 0)
                  {
                    // This case is now handled by the function
                    // additional_synchronise_hanging_nodes()
                    // called (if necessary) at the end
                    // OomphLibWarning(
                    // "Odd: missing master node is owned by current proc. Will
                    // crash below.",
                    // "TreeBasedRefineableMeshBase::synchronise_hanging_nodes(...)",
                    // OOMPH_EXCEPTION_LOCATION);
                  }
                  else // i.e. (non_halo_proc_id>=0)
                  {
                    if (shared_node_map[non_halo_proc_id].size() > 0)
                    {
                      std::map<Node*, unsigned>::iterator it =
                        shared_node_map[non_halo_proc_id].find(master_nod_pt);

                      // If it's not in there iterator points to end of
                      // set
                      if (it != shared_node_map[non_halo_proc_id].end())
                      {
                        // Found a master: Send ID of processor that holds
                        // non-halo (the fact that this is different from
                        // my_rank (the processor that sends this) will alert
                        // the other processor to the fact that it needs to
                        send_data_buff.push_back(non_halo_proc_id);

                        // Found a master: Send its index in the shared
                        // node scheme
                        send_data_buff.push_back((*it).second);

                        // Add the weight
                        send_double_data_buff.push_back(
                          hang_pt->master_weight(m));

                        // Done
                        found = true;

                        // It is possible that the master node found in the
                        // shared storage with processor <non_halo_proc_id> does
                        // not actually exist on processor d. If it does turn
                        // out to exist, we are ok, but if not then we must
                        // remember to create it later (in the
                        // additional_synchronise_hanging_nodes() function)
                      }
                    }
                  }
                }

                /*
                // Don't throw error, we will construct missing master nodes in
                // additional_synchronise_hanging_nodes() below

                // Paranoid check: if we haven't found the master node
                // then throw an error
                if (!found)
                 {
                  char filename[100];
                  std::ofstream some_file;
                  sprintf(filename,"sync_hanging_node_crash_mesh_proc%i.dat",
                          my_rank);
                  some_file.open(filename);
                  this->output(some_file);
                  some_file.close();

                  sprintf(filename,
                          "sync_hanging_node_crash_mesh_with_haloes_proc%i.dat",
                          my_rank);
                  some_file.open(filename);
                  this->enable_output_of_halo_elements();
                  this->output(some_file);
                  this->disable_output_of_halo_elements();
                  some_file.close();


                  std::set<unsigned> other_proc_id;
                  other_proc_id.insert(d);
                  other_proc_id.insert(non_halo_proc_id);
                  for (std::set<unsigned>::iterator it=other_proc_id.begin();
                       it!=other_proc_id.end();it++)
                   {
                    unsigned d_doc=(*it);

                    sprintf(
                     filename,
                     "sync_hanging_node_crash_halo_elements_with_proc%i_proc%i.dat",
                     d_doc,my_rank);
                    some_file.open(filename);
                    Vector<GeneralisedElement*>
                     halo_elem_pt(this->halo_element_pt(d_doc));
                    unsigned nelem=halo_elem_pt.size();
                    for (unsigned e=0;e<nelem;e++)
                     {
                      FiniteElement* el_pt=
                       dynamic_cast<FiniteElement*>(halo_elem_pt[e]);
                      el_pt->output(some_file);
                     }
                    some_file.close();

                    sprintf(
                     filename,
                     "sync_hanging_node_crash_haloed_elements_with_proc%i_proc%i.dat",
                     d_doc,my_rank);
                    some_file.open(filename);
                    Vector<GeneralisedElement*>
                     haloed_elem_pt(this->haloed_element_pt(d_doc));
                    nelem=haloed_elem_pt.size();
                    for (unsigned e=0;e<nelem;e++)
                     {
                      FiniteElement* el_pt=
                       dynamic_cast<FiniteElement*>(haloed_elem_pt[e]);
                      el_pt->output(some_file);
                     }
                    some_file.close();


                    sprintf(
                     filename,
                     "sync_hanging_node_crash_shared_nodes_with_proc%i_proc%i.dat",
                     d_doc,my_rank);
                    some_file.open(filename);
                    unsigned n=nshared_node(d_doc);
                    for (unsigned j=0;j<n;j++)
                     {
                      Node* nod_pt=shared_node_pt(d_doc,j);
                      unsigned nd=nod_pt->ndim();
                      for (unsigned i=0;i<nd;i++)
                       {
                        some_file << nod_pt->x(i) << " ";
                       }
                      some_file << "\n";
                     }
                    some_file.close();


                    sprintf(
                     filename,
                     "sync_hanging_node_crash_halo_nodes_with_proc%i_proc%i.dat",
                     d_doc,my_rank);
                    some_file.open(filename);
                    n=nhalo_node(d_doc);
                    for (unsigned j=0;j<n;j++)
                     {
                      Node* nod_pt=halo_node_pt(d_doc,j);
                      unsigned nd=nod_pt->ndim();
                      for (unsigned i=0;i<nd;i++)
                       {
                        some_file << nod_pt->x(i) << " ";
                       }
                      some_file << "\n";
                     }
                    some_file.close();


                    sprintf(
                     filename,
                     "sync_hanging_node_crash_haloed_nodes_with_proc%i_proc%i.dat",
                     d_doc,my_rank);
                    some_file.open(filename);
                    n=nhaloed_node(d_doc);
                    for (unsigned j=0;j<n;j++)
                     {
                      Node* nod_pt=haloed_node_pt(d_doc,j);
                      unsigned nd=nod_pt->ndim();
                      for (unsigned i=0;i<nd;i++)
                       {
                        some_file << nod_pt->x(i) << " ";
                       }
                      some_file << "\n";
                     }
                    some_file.close();

                   } // end of loop over all inter-processor lookup schemes

                  std::ostringstream error_stream;
                  unsigned n=master_nod_pt->ndim();
                  error_stream  << "Error: Master node at:\n\n";
                  for (unsigned i=0;i<n;i++)
                   {
                    error_stream <<  master_nod_pt->x(i) << " ";
                   }
                  error_stream   << "\n\nnot found for icont="
                                 << icont << "in  "
                                 << "shared node storage with proc " << d <<
                "\n"
                                 << "or in shared node storage with proc "
                                 << non_halo_proc_id
                                 << " which is where its non-halo counterpart
                lives.\n"
                                 << "Relevant files:
                sync_hanging_node_crash*.dat\n\n"
                                 << "Hanging node itself: \n\n";
                  n=nod_pt->ndim();
                  for (unsigned i=0;i<n;i++)
                   {
                    error_stream << nod_pt->x(i) << " ";
                   }
                  error_stream << nod_pt->non_halo_proc_ID();
                  error_stream << "\n\nMaster nodes:\n\n";
                  for (unsigned m=0; m<nhd_master; m++)
                   {
                    Node* master_nod_pt=hang_pt->master_node_pt(m);
                    n=master_nod_pt->ndim();
                    for (unsigned i=0;i<n;i++)
                     {
                      error_stream << master_nod_pt->x(i) << " ";
                     }
                    error_stream << master_nod_pt->non_halo_proc_ID();
                    error_stream << "\n";
                   }

                  // try to find it somewhere else -- sub-optimal search but
                  // we're about to die on our arses (yes, plural -- this is
                  // a parallel run!) anyway...
                  for (int dddd=0;dddd<n_proc;dddd++)
                   {
                    bool loc_found=false;
                    unsigned nnnod_shared=nshared_node(dddd);
                    for (unsigned k=0; k<nnnod_shared; k++)
                     {
                      if (master_nod_pt==shared_node_pt(dddd,k))
                       {
                        loc_found=true;
                        error_stream
                         << "Found that master node as " << k
                         << "-th entry in shared node storage with proc "
                         << dddd << "\n";
                       }
                     }
                    if (!loc_found)
                     {
                      error_stream
                       << "Did not find that master node in shared node storage
                with proc "
                       << dddd << "\n";
                     }
                   }
                  error_stream << "\n\n";

                  throw OomphLibError(
                   error_stream.str(),
                   OOMPH_CURRENT_FUNCTION,
                   OOMPH_EXCEPTION_LOCATION);
                 }
                */

                // Check if the master has been found
                if (!found)
                {
                  // If this master hasn't been found then set the flag
                  found_all_masters = false;
                  // No need to continue searching for masters
                  break;
                }


              } // loop over master nodes


              // Check if we need to send the data
              if (found_all_masters)
              {
                // All masters were found, so populate send data from buffer
                discrepancy = discrepancy_buff;
                discrepancy_count += discrepancy_count_buff;
                for (unsigned i = 0; i < send_data_buff.size(); i++)
                {
                  send_data.push_back(send_data_buff[i]);
                }
                for (unsigned i = 0; i < send_double_data_buff.size(); i++)
                {
                  send_double_data.push_back(send_double_data_buff[i]);
                }

                // Clear buffers and reset
                discrepancy_buff = 0;
                discrepancy_count_buff = 0;
                send_data_buff.clear();
                send_double_data_buff.clear();
              }
              else
              {
                // At least one master node was not found, so we can't
                // reconcile the hanging status of this node yet. We tell
                // the other processor to do nothing for now.
                send_data.push_back(0);

                // Clear buffers and reset
                discrepancy_buff = 0;
                discrepancy_count_buff = 0;
                send_data_buff.clear();
                send_double_data_buff.clear();

                // Set flag to trigger another round of synchronisation
                nnode_still_requiring_synchronisation++;
              }
            }
            // Haloed node isn't hanging but halo is: the latter
            // shouldn't so send a -1 to indicate that it's to be made
            // non-hanging
            else if ((haloed_hanging[d][count] == 0) &&
                     (halo_hanging[d][count] > 0))
            {
              discrepancy = 1;
              discrepancy_count++;
              send_data.push_back(-1);
            }
            // Both halo and haloed node have the same number of masters
            // we're happy!
            else if (haloed_hanging[d][count] == halo_hanging[d][count])
            {
              send_data.push_back(0);
            }
            else
            {
              std::ostringstream error_stream;
              error_stream << "Never get here!\n "
                           << "haloed_hanging[d][count]="
                           << haloed_hanging[d][count]
                           << "; halo_hanging[d][count]="
                           << halo_hanging[d][count] << std::endl;
              throw OomphLibError(error_stream.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
            // Increment counter for number of haloed data
            count++;
          } // end of loop over icont
        } // end of loop over haloed nodes

        // Now send all the required info to the equivalent halo layer -
        // If there are no discrepancies, no need to send anything
        Vector<unsigned> n_all_send(2, 0);
        if (discrepancy == 0)
        {
          MPI_Send(&n_all_send[0], 2, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());
        }
        else
        {
          // How much data is there to be sent?
          n_all_send[0] = send_data.size();
          n_all_send[1] = send_double_data.size();
          MPI_Send(&n_all_send[0], 2, MPI_UNSIGNED, d, 0, Comm_pt->mpi_comm());

          // Send flat-packed ints
          if (n_all_send[0] != 0)
          {
            MPI_Send(
              &send_data[0], n_all_send[0], MPI_INT, d, 1, Comm_pt->mpi_comm());
          }

          // Send flat-packed double data
          if (n_all_send[1] != 0)
          {
            MPI_Send(&send_double_data[0],
                     n_all_send[1],
                     MPI_DOUBLE,
                     d,
                     1,
                     Comm_pt->mpi_comm());
          }
        }
      }
      else // (d==my_rank), current process
      {
        // Receive the master nodes and weights in order to modify the
        // hanging status of nodes in the halo layer
        for (int dd = 0; dd < n_proc; dd++)
        {
          if (dd != d) // don't talk to yourself
          {
            // Are we going to receive anything? This is zero
            // either if there are no discrepancies or there's zero
            // data to be sent.
            Vector<unsigned> n_all_recv(2, 0);
            MPI_Recv(&n_all_recv[0],
                     2,
                     MPI_UNSIGNED,
                     dd,
                     0,
                     Comm_pt->mpi_comm(),
                     &status);

            // Storage for received information
            Vector<int> receive_data;
            Vector<double> receive_double_data;

            // Receive unsigneds (if any)
            if (n_all_recv[0] != 0)
            {
              // Receive the data
              receive_data.resize(n_all_recv[0]);
              MPI_Recv(&receive_data[0],
                       n_all_recv[0],
                       MPI_INT,
                       dd,
                       1,
                       Comm_pt->mpi_comm(),
                       &status);
            }

            // Receive doubles (if any)
            if (n_all_recv[1] != 0)
            {
              // Receive the data
              receive_double_data.resize(n_all_recv[1]);
              MPI_Recv(&receive_double_data[0],
                       n_all_recv[1],
                       MPI_DOUBLE,
                       dd,
                       1,
                       Comm_pt->mpi_comm(),
                       &status);
            }

            // If no information, no need to do anything else
            if (n_all_recv[0] != 0)
            {
              // Counters for traversing received data
              unsigned count = 0;
              unsigned count_double = 0;

              // Loop over halo nodes
              unsigned nh = nhalo_node(dd);
              for (unsigned j = 0; j < nh; j++)
              {
                // Get node
                Node* nod_pt = halo_node_pt(dd, j);

                // Loop over the hanging status for each interpolated variable
                // (and the geometry)
                for (int icont = -1; icont < ncont_inter_values; icont++)
                {
                  // Read next entry
                  int next_entry = receive_data[count++];

                  // If it's positive, then the number tells us how
                  // many master nodes we have
                  if (next_entry > 0)
                  {
                    unsigned nhd_master = unsigned(next_entry);

                    // Set up a new HangInfo for this node
                    HangInfo* hang_pt = new HangInfo(nhd_master);

                    // Now set up the master nodes and weights
                    for (unsigned m = 0; m < nhd_master; m++)
                    {
                      // Get the sent master node (a shared node) and
                      // the weight

                      // ID of proc in whose shared node lookup scheme
                      // the sending processor found the node
                      unsigned shared_node_proc =
                        unsigned(receive_data[count++]);

                      // Index of node in the shared node lookup scheme
                      unsigned shared_node_id = unsigned(receive_data[count++]);

                      // Get weight
                      double mtr_weight = receive_double_data[count_double++];

                      // If the shared node processor is the same as the
                      // the sending processor we can processor everything here
                      if (shared_node_proc == unsigned(dd))
                      {
                        // Get node
                        Node* master_nod_pt =
                          shared_node_pt(dd, shared_node_id);

                        // Set as a master node (with corresponding weight)
                        hang_pt->set_master_node_pt(
                          m, master_nod_pt, mtr_weight);
                      }
                      // ...otherwise we have do another communication with
                      // intermediate processor that holds the non-halo
                      // version of the master node -- only that processor can
                      // translate the index of the node the share node
                      // lookup scheme with the sending processor to the
                      // index in the shared node lookup scheme with this
                      // processor
                      else
                      {
                        // Store
                        HangHelperStruct tmp;
                        tmp.Sending_processor = dd;
                        tmp.Shared_node_id_on_sending_processor =
                          shared_node_id;
                        tmp.Shared_node_proc = shared_node_proc;
                        tmp.Weight = mtr_weight;
                        tmp.Hang_pt = hang_pt;
                        tmp.Master_node_index = m;
                        tmp.Node_pt = nod_pt;
                        tmp.icont = icont;
                        hang_info.push_back(tmp);
                      }
                    }

                    // Set the hanging pointer for the current halo node
                    // (does delete any previous hang data)
                    nod_pt->set_hanging_pt(hang_pt, icont);
                  }
                  // Negative entry: the hanging node already exists,
                  // but it shouldn't, so set it to nonhanging
                  else if (next_entry < 0)
                  {
                    nod_pt->set_hanging_pt(0, icont);
                  }

                } // end of loop over icont
              } // end of loop over nodes
            } // end of anything to receive
          }
        }
      }
    } // end loop over all processors


    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for second all-to-all in synchronise_hanging_nodes() "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }


    // Now identify master nodes by translating index in shared
    // node lookup scheme from the lookup scheme with the sending
    // processor to that with the current processor
    unsigned n = hang_info.size();

    // Is there anything to do be done?
    unsigned global_n = 0;
    MPI_Allreduce(&n, &global_n, 1, MPI_UNSIGNED, MPI_MAX, Comm_pt->mpi_comm());
    if (global_n == 0)
    {
      oomph_info
        << "No need for reconciliation of wrongly synchronised hang nodes\n";
    }
    // Reconcilation required
    else
    {
      oomph_info << "Need to reconcile of wrongly syncronised hang nodes\n";

      // Storage for the translated entries (in order) for/from
      // other processors
      // Must be ints so that an entry of -1 tells the other processor
      // that the node could not be found. See comment above for why
      // this may be necessary.
      Vector<Vector<int>> translated_entry(n_proc);

      // Storage for how-many-th entry in this processor's
      // hang_info vector will be completed by processor rank.
      Vector<Vector<unsigned>> hang_info_index_for_proc(n_proc);

      // Send information to intermediate processor that holds
      // non-halo version of missing master node
      {
        // Storage for number of data to be sent to each processor
        Vector<int> send_n(n_proc, 0);

        // Storage for all values to be sent to all processors
        Vector<int> send_data;

        // Start location within send_data for data to be sent to each processor
        Vector<int> send_displacement(n_proc, 0);

        // Loop over all processors
        for (int rank = 0; rank < n_proc; rank++)
        {
          // Set the offset for the current processor
          send_displacement[rank] = send_data.size();

          // Don't bother to do anything if the processor in the loop is the
          // current processor
          if (rank != my_rank)
          {
            // Search through the (typically few) entries
            // in hang info vector to find the ones that
            // must be sent to proc rank (the proc that holds
            // the non-halo version of the "missing" master node
            for (unsigned i = 0; i < n; i++)
            {
              HangHelperStruct tmp = hang_info[i];
              if (tmp.Shared_node_proc == unsigned(rank))
              {
                // Add the sending processor
                send_data.push_back(tmp.Sending_processor);

                // Add the index of the missing master node
                // in the shared node lookup scheme between
                // sending processor and processor rank
                send_data.push_back(tmp.Shared_node_id_on_sending_processor);

                // Record the how-many-th entry in this processor's
                // hang_info vector will be completed by processor rank.
                hang_info_index_for_proc[rank].push_back(i);
              }
            }
          }

          // Find the number of data added to the vector
          send_n[rank] = send_data.size() - send_displacement[rank];
        }


        // Storage for the number of data to be received from each processor
        Vector<int> receive_n(n_proc, 0);

        // Now send numbers of data to be sent between all processors
        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     Comm_pt->mpi_comm());

        // We now prepare the data to be received
        // by working out the displacements from the received data
        Vector<int> receive_displacement(n_proc, 0);
        int receive_data_count = 0;
        for (int rank = 0; rank < n_proc; ++rank)
        {
          // Displacement is number of data received so far
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors
        // Make sure that it has a size of at least one
        if (receive_data_count == 0)
        {
          ++receive_data_count;
        }
        Vector<unsigned> receive_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault
        if (send_data.size() == 0)
        {
          send_data.resize(1);
        }

        // Now send the data between all the processors
        MPI_Alltoallv(&send_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_INT,
                      &receive_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_INT,
                      Comm_pt->mpi_comm());

        // Now use the received data to update the halo nodes
        for (int send_rank = 0; send_rank < n_proc; send_rank++)
        {
          // Don't bother to do anything for the processor corresponding to the
          // current processor or if no data were received from this processor
          if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
          {
            // Counter for the data within the large array
            unsigned count = receive_displacement[send_rank];

            // We're reading two numbers per missing halo node
            unsigned n_rec = unsigned(receive_n[send_rank]);
            for (unsigned i = 0; i < n_rec / 2; i++)
            {
              // Receive orig sending proc
              unsigned orig_sending_proc = receive_data[count];
              count++;

              // Receive the index of the missing master node
              // in the shared node lookup scheme between
              // orig sending processor and current processor
              unsigned shared_node_id_on_orig_sending_proc =
                receive_data[count];
              count++;

              // Extract node from shared node lookup scheme
              Node* master_nod_pt = shared_node_pt(
                orig_sending_proc, shared_node_id_on_orig_sending_proc);

              // Now find it in shared halo scheme with the processor
              // that's sent the request

              std::map<Node*, unsigned>::iterator it =
                shared_node_map[send_rank].find(master_nod_pt);

              // If it's not in there iterator points to end of
              // set
              if (it != shared_node_map[send_rank].end())
              {
                // Store it so we can send it back
                translated_entry[send_rank].push_back((*it).second);
              }
              else
              {
                // This node has not been found in the shared scheme, so
                // the translation query has failed. We send a -1 to tell
                // the other processor the bad news.
                translated_entry[send_rank].push_back(-1);

                /*
                // We don't need to crash anymore because the function
                // additional_synchronise_hanging_nodes() will magically
                // sort out all the problems!
                std::ostringstream error_stream;
                error_stream
                 << "Received translation query for shared node"
                 << " entry " << shared_node_id_on_orig_sending_proc
                 << " with processor " << orig_sending_proc
                 << " from proc " << send_rank << std::endl
                 << "but did not find node in shared node scheme with proc "
                 << send_rank << std::endl;
                throw OomphLibError(
                 error_stream.str(),
                 OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);
                */
              }
            }
          }
        } // End of data is received

      } // end of sending stuff to intermediate processor that holds
        // non halo version of missing master node


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for third all-to-all in synchronise_hanging_nodes() "
          << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }


      // Send information back to processor that needs to identify
      // missing master node via shared node lookup scheme with
      // this processor
      {
        // Storage for number of data to be sent to each processor
        Vector<int> send_n(n_proc, 0);

        // Storage for all values to be sent to all processors
        Vector<int> send_data;

        // Start location within send_data for data to be sent to each processor
        Vector<int> send_displacement(n_proc, 0);

        // Loop over all processors
        for (int rank = 0; rank < n_proc; rank++)
        {
          // Set the offset for the current processor
          send_displacement[rank] = send_data.size();

          // Don't bother to do anything if the processor in the loop is the
          // current processor
          if (rank != my_rank)
          {
            // Put the translated entries for processor rank into
            // send data
            unsigned n = translated_entry[rank].size();
            for (unsigned j = 0; j < n; j++)
            {
              send_data.push_back(translated_entry[rank][j]);
            }
          }

          // Find the number of data added to the vector
          send_n[rank] = send_data.size() - send_displacement[rank];
        }


        // Storage for the number of data to be received from each processor
        Vector<int> receive_n(n_proc, 0);

        // Now send numbers of data to be sent between all processors
        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     Comm_pt->mpi_comm());

        // We now prepare the data to be received
        // by working out the displacements from the received data
        Vector<int> receive_displacement(n_proc, 0);
        int receive_data_count = 0;
        for (int rank = 0; rank < n_proc; ++rank)
        {
          // Displacement is number of data received so far
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors
        // Make sure that it has a size of at least one
        if (receive_data_count == 0)
        {
          ++receive_data_count;
        }
        Vector<unsigned> receive_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault
        if (send_data.size() == 0)
        {
          send_data.resize(1);
        }

        // Now send the data between all the processors
        MPI_Alltoallv(&send_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_INT,
                      &receive_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_INT,
                      Comm_pt->mpi_comm());

        // Now use the received data to update the halo nodes
        for (int send_rank = 0; send_rank < n_proc; send_rank++)
        {
          // Don't bother to do anything for the processor corresponding to the
          // current processor or if no data were received from this processor
          if ((send_rank != my_rank) && (receive_n[send_rank] != 0))
          {
            // Counter for the data within the large array
            unsigned count = receive_displacement[send_rank];

            // We're reading one number per missing halo node
            unsigned n_rec = unsigned(receive_n[send_rank]);
            for (unsigned i = 0; i < n_rec; i++)
            {
              // Index of missing master node in shared node lookup scheme
              // with processor send_rank:
              // Must be an int because failure returns -1
              int index = receive_data[count];
              count++;

              // Translation query has been successful if index >= 0
              if (index >= 0)
              {
                // Recall information associated with that missing master
                unsigned hang_info_index =
                  hang_info_index_for_proc[send_rank][i];
                HangHelperStruct tmp = hang_info[hang_info_index];

                // Extract node from shared node lookup scheme
                Node* master_nod_pt = shared_node_pt(send_rank, index);

                // Set as a master node (with corresponding weight)
                tmp.Hang_pt->set_master_node_pt(
                  tmp.Master_node_index, master_nod_pt, tmp.Weight);
              }
              else
              {
                // Translation query has failed. This is the processor
                // on which the node was a halo, so we must delete the
                // partial hang info.

                // Recall information associated with that missing master
                unsigned hang_info_index =
                  hang_info_index_for_proc[send_rank][i];
                HangHelperStruct tmp = hang_info[hang_info_index];

                // Delete partial hanging information
                tmp.Node_pt->set_hanging_pt(0, tmp.icont);

                // Set flag to trigger another round of synchronisation
                // This works even though we don't own the node that
                // still requires synchrionisation because this variable
                // is reduced over all processors at the end
                nnode_still_requiring_synchronisation++;
              }
            }
          }
        } // End of data is received

      } // end of completing hang info for missing master nodes


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info
          << "Time for fourth all-to-all in synchronise_hanging_nodes() "
          << t_end - t_start << std::endl;
      }
    } // end  of reconciliation required


    // Get global number of nodes still requiring synchronisation due to
    // missing master nodes
    // This will only be necessary for meshes involving elements
    // with nonuniformly spaced nodes. All other cases will continue to
    // work as before because all nodes will have been successfully
    // synchronised by now
    unsigned global_nnode_still_requiring_synchronisation = 0;
    MPI_Allreduce(&nnode_still_requiring_synchronisation,
                  &global_nnode_still_requiring_synchronisation,
                  1,
                  MPI_UNSIGNED,
                  MPI_MAX,
                  Comm_pt->mpi_comm());
    if (global_nnode_still_requiring_synchronisation > 0)
    {
      double tt_start = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_start = TimingHelpers::timer();
      }

      oomph_info << "Need to do additional synchronisation of hanging nodes"
                 << std::endl;

      // Do additional synchronisation
      additional_synchronise_hanging_nodes(ncont_interpolated_values);

      double tt_end = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        tt_end = TimingHelpers::timer();
        oomph_info
          << "Time for RefineableMesh::additional_synchronise_hanging_nodes() "
          << "in TreeBasedRefineableMeshBase::synchronise_hanging_nodes(): "
          << tt_end - tt_start << std::endl;
        tt_start = TimingHelpers::timer();
      }
    }
    else
    {
      oomph_info << "No need to do additional synchronisation of hanging nodes"
                 << std::endl;
    }
  }

  //========================================================================
  /// Synchronise the positions of non-hanging nodes that depend on
  /// non-existent neighbours (e.g. h-refinement of neighbouring elements
  /// with different p-orders where the shared edge is on the outer edge of
  /// the halo layer)
  //========================================================================
  void TreeBasedRefineableMeshBase::synchronise_nonhanging_nodes()
  {
    // Store number of processors and current process
    MPI_Status status;
    int n_proc = Comm_pt->nproc();
    int my_rank = Comm_pt->my_rank();

    double t_start = 0.0;
    double t_end = 0.0;

    // Storage for the hanging status of halo/haloed nodes on elements
    Vector<Vector<unsigned>> recv_unsigneds(n_proc);
    Vector<Vector<double>> recv_doubles(n_proc);

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_start = TimingHelpers::timer();
    }

    // Loop over processes: Each processor checks if its nonhaning nodes in
    // haloed elements with proc d require additional information to determine
    // their positions.
    for (int d = 0; d < n_proc; d++)
    {
      // No halo with self: Setup hang info for my haloed nodes with proc d
      // then get ready to receive halo info from processor d.
      if (d != my_rank)
      {
        // Receive the position information from the corresponding process
        unsigned recv_unsigneds_count = 0;
        MPI_Recv(&recv_unsigneds_count,
                 1,
                 MPI_UNSIGNED,
                 d,
                 0,
                 Comm_pt->mpi_comm(),
                 &status);
        unsigned recv_doubles_count = 0;
        MPI_Recv(&recv_doubles_count,
                 1,
                 MPI_UNSIGNED,
                 d,
                 1,
                 Comm_pt->mpi_comm(),
                 &status);

        // Get the data (if any)
        if (recv_unsigneds_count != 0)
        {
          recv_unsigneds[d].resize(recv_unsigneds_count);
          MPI_Recv(&recv_unsigneds[d][0],
                   recv_unsigneds_count,
                   MPI_UNSIGNED,
                   d,
                   0,
                   Comm_pt->mpi_comm(),
                   &status);
        }
        if (recv_doubles_count != 0)
        {
          recv_doubles[d].resize(recv_doubles_count);
          MPI_Recv(&recv_doubles[d][0],
                   recv_doubles_count,
                   MPI_DOUBLE,
                   d,
                   1,
                   Comm_pt->mpi_comm(),
                   &status);
        }

        // Counters for received data
        unsigned recv_unsigneds_index = 0;
        double recv_doubles_index = 0;

        // Get halo elements with processor d
        Vector<GeneralisedElement*> halo_element_pt(this->halo_element_pt(d));

        // Loop over recieved indices
        while (recv_unsigneds_index < recv_unsigneds_count)
        {
          // Get (finite) element
          FiniteElement* el_pt = dynamic_cast<FiniteElement*>(
            halo_element_pt[recv_unsigneds[d][recv_unsigneds_index++]]);

          // If we have a finite element...
          if (el_pt != 0)
          {
            // Get dimension
            unsigned n_dim = el_pt->dim();

            // Get node
            Node* nod_pt =
              el_pt->node_pt(recv_unsigneds[d][recv_unsigneds_index++]);

            // Get current position
            Vector<double> x_cur(n_dim);
            for (unsigned dir = 0; dir < n_dim; dir++)
            {
              x_cur[dir] = nod_pt->x(dir);
            }

            // Get recieved position
            Vector<double> x_rec(n_dim);
            for (unsigned dir = 0; dir < n_dim; dir++)
            {
              x_rec[dir] = recv_doubles[d][recv_doubles_index + dir];
            }

            // Compare actual and expected positions
            bool node_pos_differs = false;
            for (unsigned dir = 0; dir < n_dim; dir++)
            {
              node_pos_differs = node_pos_differs ||
                                 (std::fabs(x_cur[dir] - x_rec[dir]) > 1.0e-14);
            }

            // Set the actual position
            Vector<double> x_act(n_dim);
            for (unsigned dir = 0; dir < n_dim; dir++)
            {
              nod_pt->x(dir) = recv_doubles[d][recv_doubles_index++];
            }
          }
        }

        if (recv_unsigneds_count != recv_unsigneds_index)
        {
          std::ostringstream error_stream;
          error_stream << "recv_unsigneds_count != recv_unsigneds_index ( "
                       << recv_unsigneds_count << " != " << recv_unsigneds_index
                       << ")" << std::endl;
          throw OomphLibError(
            error_stream.str(),
            "TreeBasedRefineableMeshBase::synchronise_nonhanging_nodes()",
            OOMPH_EXCEPTION_LOCATION);
        }
      }
      else // d==my_rank, i.e. current process: Send halo hanging status
           // to process dd where it's received (see above) and compared
           // and compared against the hang status of the haloed nodes
      {
        for (int dd = 0; dd < n_proc; dd++)
        {
          // No halo with yourself
          if (dd != d)
          {
            // Storage for halo hanging status and counter
            Vector<int> send_unsigneds;
            Vector<double> send_doubles;

            // Set to store nodes whose position requires adjustment
            std::set<Node*> nodes_requiring_adjustment;

            // Get haloed elements with processor dd
            Vector<GeneralisedElement*> haloed_element_pt(
              this->haloed_element_pt(dd));

            // Loop over haloed elements with processor dd
            unsigned nh = haloed_element_pt.size();
            for (unsigned e = 0; e < nh; e++)
            {
              // Get (finite) element
              FiniteElement* el_pt =
                dynamic_cast<FiniteElement*>(haloed_element_pt[e]);

              // If we have a finite element...
              if (el_pt != 0)
              {
                // Get dimension
                unsigned n_dim = el_pt->dim();

                // Loop over element nodes
                unsigned n_node = el_pt->nnode();
                for (unsigned j = 0; j < n_node; j++)
                {
                  // Get node
                  Node* nod_pt = el_pt->node_pt(j);

                  // Only do non-hanging nodes
                  if (!nod_pt->is_hanging())
                  {
                    // Check if node's position is the same as that interpolated
                    // using its local coordinate in the haloed element

                    // Loop over all history values
                    unsigned nt = nod_pt->ntstorage();
                    for (unsigned t = 0; t < nt; t++)
                    {
                      // Get expected position
                      Vector<double> s(n_dim), x_exp(n_dim);
                      el_pt->local_coordinate_of_node(j, s);
                      el_pt->get_x(t, s, x_exp);

                      // Get actual position
                      Vector<double> x_act(n_dim);
                      for (unsigned dir = 0; dir < n_dim; dir++)
                      {
                        x_act[dir] = nod_pt->x(dir);
                      }

                      // Compare actual and expected positions
                      bool node_pos_differs = false;
                      for (unsigned dir = 0; dir < n_dim; dir++)
                      {
                        node_pos_differs =
                          node_pos_differs ||
                          (std::fabs(x_act[dir] - x_exp[dir]) > 1.0e-14);
                      }

                      // If the node's actual position differs from its
                      // expected position we need to communicate this
                      // information to processors on which this is a halo node
                      if (node_pos_differs)
                      {
                        // Check that node has not been done already
                        if (nodes_requiring_adjustment.insert(nod_pt).second)
                        {
                          // Send index of haloed element
                          send_unsigneds.push_back(e);
                          // Send index of node in the element
                          send_unsigneds.push_back(j);
                          // Send actual position of node
                          for (unsigned dir = 0; dir < n_dim; dir++)
                          {
                            send_doubles.push_back(x_act[dir]);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

            // Send the information to the relevant process
            unsigned send_unsigneds_count = send_unsigneds.size();
            unsigned send_doubles_count = send_doubles.size();

            if (send_unsigneds_count > 0)
            {
              // exit(1);
            }

            // Tell processor dd how much data to receive
            MPI_Send(&send_unsigneds_count,
                     1,
                     MPI_UNSIGNED,
                     dd,
                     0,
                     Comm_pt->mpi_comm());
            MPI_Send(
              &send_doubles_count, 1, MPI_UNSIGNED, dd, 1, Comm_pt->mpi_comm());

            // Send data (if any)
            if (send_unsigneds_count != 0)
            {
              MPI_Send(&send_unsigneds[0],
                       send_unsigneds_count,
                       MPI_UNSIGNED,
                       dd,
                       0,
                       Comm_pt->mpi_comm());
            }
            if (send_doubles_count != 0)
            {
              MPI_Send(&send_doubles[0],
                       send_doubles_count,
                       MPI_DOUBLE,
                       dd,
                       1,
                       Comm_pt->mpi_comm());
            }
          }
        }
      }
    }

    if (Global_timings::Doc_comprehensive_timings)
    {
      t_end = TimingHelpers::timer();
      oomph_info << "Time for synchronise_nonhanging_nodes(): "
                 << t_end - t_start << std::endl;
      t_start = TimingHelpers::timer();
    }
  }

#endif


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Do adaptive p-refinement for mesh.
  /// - Pass Vector of error estimates for all elements.
  /// - p-refine those whose errors exceeds the threshold
  /// - p-unrefine those whose errors is less than
  ///   threshold.
  //========================================================================
  void TreeBasedRefineableMeshBase::p_adapt(
    const Vector<double>& elemental_error)
  {
    // Set the refinement tolerance to be the max permissible error
    double refine_tol = this->max_permitted_error();

    // Set the unrefinement tolerance to be the min permissible error
    double unrefine_tol = this->min_permitted_error();

    // Setup doc info
    DocInfo local_doc_info;
    if (doc_info_pt() == 0)
    {
      local_doc_info.disable_doc();
    }
    else
    {
      local_doc_info = this->doc_info();
    }


    // Check that the errors make sense
    if (refine_tol <= unrefine_tol)
    {
      std::ostringstream error_stream;
      error_stream << "Refinement tolerance <= Unrefinement tolerance"
                   << refine_tol << " " << unrefine_tol << std::endl
                   << "doesn't make sense and will almost certainly crash"
                   << std::endl
                   << "this beautiful code!" << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Select elements for refinement and unrefinement
    //==============================================
    // Reset counter for number of elements that would like to be
    // refined further but can't
    this->nrefinement_overruled() = 0;

    unsigned n_refine = 0;
    unsigned n_unrefine = 0;
    // Loop over all elements and mark them according to the error criterion
    unsigned long Nelement = this->nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      //(Cast) pointer to the element
      PRefineableElement* el_pt =
        dynamic_cast<PRefineableElement*>(this->element_pt(e));

      // Check that we can p-refine the element
      if (el_pt != 0)
      {
        // Initially element is not to be refined
        el_pt->deselect_for_p_refinement();
        el_pt->deselect_for_p_unrefinement();

        // If the element error exceeds the threshold ...
        if (elemental_error[e] > refine_tol)
        {
          // ... and its refinement level is less than the maximum desired level
          // mark is to be refined
          if ((el_pt->p_refinement_is_enabled()) &&
              (el_pt->p_order() < this->max_p_refinement_level()))
          {
            el_pt->select_for_p_refinement();
            n_refine++;
          }
          // ... otherwise mark it as having been over-ruled
          else
          {
            this->nrefinement_overruled() += 1;
          }
        }
        if (elemental_error[e] < unrefine_tol)
        {
          // ... and its refinement level is more than the minimum desired level
          // mark is to be refined
          //(Also check that we don't unrefine past the initial refinement
          // level)
          if ((el_pt->p_refinement_is_enabled()) &&
              (el_pt->p_order() > this->min_p_refinement_level()) &&
              (el_pt->p_order() > el_pt->initial_p_order()))
          {
            el_pt->select_for_p_unrefinement();
            n_unrefine++;
          }
          // Don't mark as overruled - it's misleading
          //// ... otherwise mark it as having been over-ruled
          // else
          // {
          //  this->nrefinement_overruled()+=1;
          // }
        }
      } // End of check for p-refineability of element
      else
      {
        oomph_info << "p-refinement is not possible for these elements"
                   << std::endl;
        // Don't try to p-refine any more elements
        break;
      }
    }

    oomph_info << " \n Number of elements to be refined: " << n_refine
               << std::endl;
    oomph_info << " \n Number of elements whose refinement was overruled: "
               << this->nrefinement_overruled() << std::endl;

    oomph_info << " \n Number of elements to be unrefined : " << n_unrefine
               << std::endl
               << std::endl;


    // Now do the actual mesh adaptation
    //---------------------------------

    // Check whether its worth our while
    // Either some elements want to be refined,
    // or the number that want to be unrefined are greater than the
    // specified tolerance

    // In a parallel job, it is possible that one process may not have
    // any elements to refine, BUT a neighbouring process may refine an
    // element which changes the hanging status of a node that is on
    // both processes (i.e. a halo(ed) node).  To get around this issue,
    // ALL processes need to call adapt_mesh if ANY refinement is to
    // take place anywhere.

    unsigned total_n_refine = 0;
#ifdef OOMPH_HAS_MPI
    // Sum n_refine across all processors
    if (this->is_mesh_distributed())
    {
      MPI_Allreduce(
        &n_refine, &total_n_refine, 1, MPI_INT, MPI_SUM, Comm_pt->mpi_comm());
    }
    else
    {
      total_n_refine = n_refine;
    }
#else
    total_n_refine = n_refine;
#endif

    // There may be some issues with unrefinement too, but I have not
    // been able to come up with an example (either in my head or in a
    // particular problem) where anything has arisen.  I can see that
    // there may be an issue if n_unrefine differs across processes so
    // that (total_n_unrefine > max_keep_unrefined()) on some but not
    // all processes. I haven't seen any examples of this yet so the
    // following code may or may not work!  (Andy, 06/03/08)

    unsigned total_n_unrefine = 0;
#ifdef OOMPH_HAS_MPI
    // Sum n_unrefine across all processors
    if (this->is_mesh_distributed())
    {
      MPI_Allreduce(&n_unrefine,
                    &total_n_unrefine,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    Comm_pt->mpi_comm());
    }
    else
    {
      total_n_unrefine = n_unrefine;
    }
#else
    total_n_unrefine = n_unrefine;
#endif

    oomph_info << "---> " << total_n_refine << " elements to be refined, and "
               << total_n_unrefine << " to be unrefined, in total."
               << std::endl;

    if ((total_n_refine > 0) || (total_n_unrefine > this->max_keep_unrefined()))
    {
#ifdef PARANOID
#ifdef OOMPH_HAS_MPI

      // Sanity check: Each processor checks if the enforced unrefinement of
      // its haloed element is matched by enforced unrefinement of the
      // corresponding halo elements on the other processors.
      if (this->is_mesh_distributed())
      {
        // Store number of processors and current process
        MPI_Status status;
        int n_proc = Comm_pt->nproc();
        int my_rank = Comm_pt->my_rank();

        // Loop over all other domains/processors
        for (int d = 0; d < n_proc; d++)
        {
          // Don't talk to yourself
          if (d != my_rank)
          {
            {
              // Get the vector of halo elements whose non-halo counterpart
              // are on processor d
              Vector<GeneralisedElement*> halo_elem_pt(
                this->halo_element_pt(d));

              // Create vector containing (0)1 to indicate that
              // halo element is (not) to be unrefined
              unsigned nhalo = halo_elem_pt.size();
              Vector<int> halo_to_be_unrefined(nhalo, 0);
              for (unsigned e = 0; e < nhalo; e++)
              {
                if (dynamic_cast<PRefineableElement*>(halo_elem_pt[e])
                      ->to_be_p_unrefined())
                {
                  halo_to_be_unrefined[e] = 1;
                }
              }

              // Trap the case when there are no halo elements
              // so that we don't get a segfault in the MPI send
              if (nhalo > 0)
              {
                // Send it across
                MPI_Send(&halo_to_be_unrefined[0],
                         nhalo,
                         MPI_INT,
                         d,
                         0,
                         Comm_pt->mpi_comm());
              }
            }

            {
              // Get the vector of haloed elements on current processor
              Vector<GeneralisedElement*> haloed_elem_pt(
                this->haloed_element_pt(d));

              // Ask processor d to send vector containing (0)1 for
              // halo element with current processor to be (not)unrefined
              unsigned nhaloed = haloed_elem_pt.size();
              Vector<int> halo_to_be_unrefined(nhaloed);
              // Trap to catch the case that there are no haloed elements
              if (nhaloed > 0)
              {
                MPI_Recv(&halo_to_be_unrefined[0],
                         nhaloed,
                         MPI_INT,
                         d,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);
              }

              // Check it
              for (unsigned e = 0; e < nhaloed; e++)
              {
                if (((halo_to_be_unrefined[e] == 0) &&
                     (dynamic_cast<PRefineableElement*>(haloed_elem_pt[e])
                        ->to_be_p_unrefined())) ||
                    ((halo_to_be_unrefined[e] == 1) &&
                     (!dynamic_cast<PRefineableElement*>(haloed_elem_pt[e])
                         ->to_be_p_unrefined())))
                {
                  std::ostringstream error_message;
                  error_message
                    << "Error in refinement: \n"
                    << "Haloed element: " << e << " on proc " << my_rank
                    << " \n"
                    << "wants to be unrefined whereas its halo counterpart on\n"
                    << "proc " << d << " doesn't (or vice versa)...\n"
                    << "This is most likely because the error estimator\n"
                    << "has not assigned the same errors to halo and haloed\n"
                    << "elements -- it ought to!\n";
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
              }
            }
          }
        }


        // Loop over all other domains/processors
        for (int d = 0; d < n_proc; d++)
        {
          // Don't talk to yourself
          if (d != my_rank)
          {
            {
              // Get the vector of halo elements whose non-halo counterpart
              // are on processor d
              Vector<GeneralisedElement*> halo_elem_pt(
                this->halo_element_pt(d));

              // Create vector containing (0)1 to indicate that
              // halo element is (not) to be refined
              unsigned nhalo = halo_elem_pt.size();
              Vector<int> halo_to_be_refined(nhalo, 0);
              for (unsigned e = 0; e < nhalo; e++)
              {
                if (dynamic_cast<PRefineableElement*>(halo_elem_pt[e])
                      ->to_be_p_refined())
                {
                  halo_to_be_refined[e] = 1;
                }
              }

              // Send it across
              if (nhalo > 0)
              {
                MPI_Send(&halo_to_be_refined[0],
                         nhalo,
                         MPI_INT,
                         d,
                         0,
                         Comm_pt->mpi_comm());
              }
            }

            {
              // Get the vector of haloed elements on current processor
              Vector<GeneralisedElement*> haloed_elem_pt(
                this->haloed_element_pt(d));

              // Ask processor d to send vector containing (0)1 for
              // halo element with current processor to be (not)refined
              unsigned nhaloed = haloed_elem_pt.size();
              Vector<int> halo_to_be_refined(nhaloed);
              if (nhaloed > 0)
              {
                MPI_Recv(&halo_to_be_refined[0],
                         nhaloed,
                         MPI_INT,
                         d,
                         0,
                         Comm_pt->mpi_comm(),
                         &status);
              }

              // Check it
              for (unsigned e = 0; e < nhaloed; e++)
              {
                if (((halo_to_be_refined[e] == 0) &&
                     (dynamic_cast<PRefineableElement*>(haloed_elem_pt[e])
                        ->to_be_p_refined())) ||
                    ((halo_to_be_refined[e] == 1) &&
                     (!dynamic_cast<PRefineableElement*>(haloed_elem_pt[e])
                         ->to_be_p_refined())))
                {
                  std::ostringstream error_message;
                  error_message
                    << "Error in refinement: \n"
                    << "Haloed element: " << e << " on proc " << my_rank
                    << " \n"
                    << "wants to be refined whereas its halo counterpart on\n"
                    << "proc " << d << " doesn't (or vice versa)...\n"
                    << "This is most likely because the error estimator\n"
                    << "has not assigned the same errors to halo and haloed\n"
                    << "elements -- it ought to!\n";
                  throw OomphLibError(error_message.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
              }
            }
          }
        }
      }
#endif
#endif

      // Perform the actual adaptation
      p_adapt_mesh(local_doc_info);

      // The number of refineable elements is still local to each process
      this->Nunrefined = n_unrefine;
      this->Nrefined = n_refine;
    }
    // If not worthwhile, say so but still reorder nodes and kill external
    // storage for consistency in parallel computations
    else
    {
#ifdef OOMPH_HAS_MPI
      // Delete any external element storage - any interaction will still
      // be set up on the fly again, so we need to get rid of old information.
      // This particularly causes problems in multi-domain examples where
      // we decide not to refine one of the meshes
      this->delete_all_external_storage();
#endif

      // Reorder the nodes within the mesh's node vector
      // to establish a standard ordering regardless of the sequence
      // of mesh refinements -- this is required to allow dump/restart
      // on refined meshes
      this->reorder_nodes();

#ifdef OOMPH_HAS_MPI

      // Now (re-)classify halo and haloed nodes and synchronise hanging
      // nodes
      // This is required in cases where delete_all_external_storage()
      // made dependent nodes to external halo nodes nonhanging.
      if (this->is_mesh_distributed())
      {
        DocInfo doc_info;
        doc_info.disable_doc();
        classify_halo_and_haloed_nodes(doc_info, doc_info.is_doc_enabled());
      }

#endif

      if (n_refine == 0)
      {
        oomph_info << "\n Not enough benefit in adapting mesh. " << std::endl
                   << std::endl;
      }
      this->Nunrefined = 0;
      this->Nrefined = 0;
    }
  }


  //================================================================
  /// p-adapt mesh, which exists in two representations,
  /// namely as:
  ///  - a FE mesh
  ///  - a forest of Oc or QuadTrees
  ///
  /// p-refinement/derefinement process is documented (in tecplot-able form)
  /// if requested.
  ///
  /// Procedure:
  /// - Loop over all elements and do the p-refinement/unrefinement for
  ///   those who want to be refined. Note: p-refinement builds fully-
  ///   functional elements.
  /// - For all nodes that were hanging on the previous mesh (and are still
  ///   marked as such), fill in their nodal values (consistent
  ///   with the current hanging node scheme) to make sure they are fully
  ///   functional, should they have become non-hanging during the
  ///   mesh-adaptation. Then mark the nodes as non-hanging.
  /// - Delete any nodes that have become obsolete.
  /// - Mark up hanging nodes and setup hanging node scheme (incl.
  ///   recursive cleanup for hanging nodes that depend on other
  ///   hanging nodes).
  /// - run a quick self-test on the neighbour finding scheme and
  ///   check the integrity of the elements (if PARANOID)
  /// - doc hanging node status, boundary conditions, neighbour
  ///   scheme if requested.
  ///
  ///
  /// After adaptation, all nodes (whether new or old) have up-to-date
  /// current and previous values.
  ///
  /// If refinement process is being documented, the following information
  /// is documented:
  /// - The files
  ///   - "neighbours.dat"
  ///   - "all_nodes.dat"
  ///   - "new_nodes.dat"
  ///   - "hang_nodes_*.dat"
  ///     where the * denotes a direction (n,s,e,w) in 2D
  ///     or (r,l,u,d,f,b) in 3D
  ///.
  ///   can be viewed with
  ///   - QHangingNodes.mcr
  ///   .
  /// - The file
  ///    - "hangnodes_withmasters.dat"
  ///    .
  ///    can be viewed with
  ///    - QHangingNodesWithMasters.mcr
  ///    .
  ///    to check the hanging node status.
  /// - The neighbour status of the elements is documented in
  ///   - "neighbours.dat"
  ///   .
  ///   and can be viewed with
  ///   - QuadTreeNeighbours.mcr
  ///   .
  //=================================================================
  void TreeBasedRefineableMeshBase::p_adapt_mesh(DocInfo& doc_info)
  {
#ifdef OOMPH_HAS_MPI
    // Delete any external element storage before performing the adaptation
    // (in particular, external halo nodes that are on mesh boundaries)
    this->delete_all_external_storage();
#endif

    // Only perform the adapt step if the mesh has any elements.  This is
    // relevant in a distributed problem with multiple meshes, where a
    // particular process may not have any elements on a particular submesh.
    if (this->nelement() > 0)
    {
      double t_start = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_start = TimingHelpers::timer();
      }

      // Do refinement/unrefinement if required
      this->p_refine_elements_if_required();

      if (Global_timings::Doc_comprehensive_timings)
      {
        double t_end = TimingHelpers::timer();
        oomph_info << "Time for p-refinement/unrefinement: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Loop over all nodes in mesh and free the dofs of those that were
      //-----------------------------------------------------------------
      // pinned only because they were hanging nodes. Also update their
      //-----------------------------------------------------------------
      // nodal values so that they contain data that is consistent
      //----------------------------------------------------------
      // with the hanging node representation
      //-------------------------------------
      // (Even if the nodal data isn't actually accessed because the node
      // is still hanging -- we don't know this yet, and this step makes
      // sure that all nodes are fully functional and up-to-date, should
      // they become non-hanging below).
      //
      //
      // However, if we have a fixed mesh and hanging nodes on the boundary
      // become non-hanging they will not necessarily respect the curvilinear
      // boundaries. This can only happen in 3D of course because it is not
      // possible to have hanging nodes on boundaries in 2D.
      // The solution is to store those nodes on the boundaries that are
      // currently hanging and then check to see whether they have changed
      // status at the end of the refinement procedure.
      // If it has changed, then we need to adjust their positions.
      const unsigned n_boundary = this->nboundary();
      const unsigned mesh_dim = this->finite_element_pt(0)->dim();
      Vector<std::set<Node*>> hanging_nodes_on_boundary_pt(n_boundary);

      unsigned long n_node = this->nnode();
      for (unsigned long n = 0; n < n_node; n++)
      {
        // Get the pointer to the node
        Node* nod_pt = this->node_pt(n);

        // Get the number of values in the node
        unsigned n_value = nod_pt->nvalue();

        // We need to find if any of the values are hanging
        bool is_hanging = nod_pt->is_hanging();
        // Loop over the values and find out whether any are hanging
        for (unsigned n = 0; n < n_value; n++)
        {
          is_hanging |= nod_pt->is_hanging(n);
        }

        // If the node is hanging then ...
        if (is_hanging)
        {
          // Unless they are turned into hanging nodes again below
          // (this might or might not happen), fill in all the necessary
          // data to make them 'proper' nodes again.

          // Reconstruct the nodal values/position from the node's
          // hanging node representation
          unsigned nt = nod_pt->ntstorage();
          Vector<double> values(n_value);
          unsigned n_dim = nod_pt->ndim();
          Vector<double> position(n_dim);
          // Loop over all history values
          for (unsigned t = 0; t < nt; t++)
          {
            nod_pt->value(t, values);
            for (unsigned i = 0; i < n_value; i++)
            {
              nod_pt->set_value(t, i, values[i]);
            }
            nod_pt->position(t, position);
            for (unsigned i = 0; i < n_dim; i++)
            {
              nod_pt->x(t, i) = position[i];
            }
          }

          // If it's an algebraic node: Update its previous nodal positions too
          AlgebraicNode* alg_node_pt = dynamic_cast<AlgebraicNode*>(nod_pt);
          if (alg_node_pt != 0)
          {
            bool update_all_time_levels = true;
            alg_node_pt->node_update(update_all_time_levels);
          }


          // If it's a Solid node, update Lagrangian coordinates
          // from its hanging node representation
          SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
          if (solid_node_pt != 0)
          {
            unsigned n_lagrangian = solid_node_pt->nlagrangian();
            for (unsigned i = 0; i < n_lagrangian; i++)
            {
              solid_node_pt->xi(i) = solid_node_pt->lagrangian_position(i);
            }
          }

          // Now store geometrically hanging nodes on boundaries that
          // may need updating after refinement.
          // There will only be a problem if we have 3 spatial dimensions
          if ((mesh_dim > 2) && (nod_pt->is_hanging()))
          {
            // If the node is on a boundary then add a pointer to the node
            // to our lookup scheme
            if (nod_pt->is_on_boundary())
            {
              // Storage for the boundaries on which the Node is located
              std::set<unsigned>* boundaries_pt;
              nod_pt->get_boundaries_pt(boundaries_pt);
              if (boundaries_pt != 0)
              {
                // Loop over the boundaries and add a pointer to the node
                // to the appropriate storage scheme
                for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                     it != boundaries_pt->end();
                     ++it)
                {
                  hanging_nodes_on_boundary_pt[*it].insert(nod_pt);
                }
              }
            }
          }

        } // End of is_hanging

        // Initially mark all nodes as 'non-hanging' and `obsolete'
        nod_pt->set_nonhanging();
        nod_pt->set_obsolete();
      }

      double t_end = 0.0;
      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for sorting out initial hanging status: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Stick all elements into a vector
      Vector<Tree*> tree_nodes_pt;
      this->forest_pt()->stick_leaves_into_vector(tree_nodes_pt);

      // Copy the elements into the mesh Vector
      unsigned long num_tree_nodes = tree_nodes_pt.size();
      this->element_pt().resize(num_tree_nodes);
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        this->element_pt(e) = tree_nodes_pt[e]->object_pt();

        // Now loop over all nodes in element and mark them as non-obsolete
        FiniteElement* this_el_pt = this->finite_element_pt(e);
        unsigned n_node = this_el_pt->nnode(); // caching pre-loop
        for (unsigned n = 0; n < n_node; n++)
        {
          this_el_pt->node_pt(n)->set_non_obsolete();
        }

        // Mark up so that repeated refinements do not occur
        // (Required because refined element is the same element as the
        // original)
        PRefineableElement* cast_el_pt =
          dynamic_cast<PRefineableElement*>(this->element_pt(e));
        cast_el_pt->deselect_for_p_refinement();
        cast_el_pt->deselect_for_p_unrefinement();
      }

      // Cannot delete nodes that are still marked as obsolete
      // because they may still be required to assemble the hanging schemes
      //-------------------------------------------------------------------

      // Mark up hanging nodes
      //----------------------

      // Output streams for the hanging nodes
      Vector<std::ofstream*> hanging_output_files;
      // Setup the output files for hanging nodes, this must be called
      // precisely once for the forest. Note that the files will only
      // actually be opened if doc_info.Doc_flag is true
      this->forest_pt()->open_hanging_node_files(doc_info,
                                                 hanging_output_files);

      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        // Generic setup
        tree_nodes_pt[e]->object_pt()->setup_hanging_nodes(
          hanging_output_files);
        // Element specific setup
        tree_nodes_pt[e]->object_pt()->further_setup_hanging_nodes();
      }

      // Close the hanging node files and delete the memory allocated
      // for the streams
      this->forest_pt()->close_hanging_node_files(doc_info,
                                                  hanging_output_files);


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for setup_hanging_nodes() and "
                      "further_setup_hanging_nodes() for "
                   << num_tree_nodes << " elements: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Read out the number of continously interpolated values
      // from one of the elements (assuming it's the same in all elements)
      unsigned ncont_interpolated_values =
        tree_nodes_pt[0]->object_pt()->ncont_interpolated_values();

      // Complete the hanging nodes schemes by dealing with the
      // recursively hanging nodes
      this->complete_hanging_nodes(ncont_interpolated_values);


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for complete_hanging_nodes: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      /// Update the boundary element info -- this can be a costly procedure
      /// and for this reason the mesh writer might have decided not to set up
      /// this scheme. If so, we won't change this and suppress its creation...
      if (Lookup_for_elements_next_boundary_is_setup)
      {
        this->setup_boundary_element_info();
      }

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for boundary element info: " << t_end - t_start
                   << std::endl;
        t_start = TimingHelpers::timer();
      }

      // BENFLAG: Reset all the node update elements.
      //         This is necessary to prevent the following case: A node N is
      //         shared between two elements, A and B. The update element for
      //         the node is set to A, say. Element A is p-refined and now
      //         nolonger has N as a node. However the node update element for N
      //         is still A but the node doesn't exist in A.
      MacroElementNodeUpdateElementBase* first_macro_el_pt =
        dynamic_cast<MacroElementNodeUpdateElementBase*>(this->element_pt(0));
      if (first_macro_el_pt != 0)
      {
        // Now set the node update info elementwise
        for (unsigned e = 0; e < this->nelement(); e++)
        {
          // Cast to macro element
          MacroElementNodeUpdateElementBase* macro_el_pt =
            dynamic_cast<MacroElementNodeUpdateElementBase*>(
              this->element_pt(e));
          if (macro_el_pt != 0)
          {
            // Get vector of geometric objects from element (construct vector
            // via copy operation)
            Vector<GeomObject*> geom_object_pt(macro_el_pt->geom_object_pt());

            // (Re)set node update info for all the nodes in the element
            macro_el_pt->set_node_update_info(geom_object_pt);
          }
        }
      }

#ifdef PARANOID

      // Doc/check the neighbours
      //-------------------------
      Vector<Tree*> all_tree_nodes_pt;
      this->forest_pt()->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

      // Check the neighbours
      this->forest_pt()->check_all_neighbours(doc_info);

      // Check the integrity of the elements
      // -----------------------------------

      // Loop over elements and get the elemental integrity
      double max_error = 0.0;
      for (unsigned long e = 0; e < num_tree_nodes; e++)
      {
        double max_el_error;
        tree_nodes_pt[e]->object_pt()->check_integrity(max_el_error);
        // If the elemental error is greater than our maximum error
        // reset the maximum
        if (max_el_error > max_error)
        {
          max_error = max_el_error;
        }
      }

      if (max_error > RefineableElement::max_integrity_tolerance())
      {
        std::ostringstream error_stream;
        error_stream << "Mesh refined: Max. error in integrity check: "
                     << max_error << " is too big"
                     << "\ni.e. bigger than RefineableElement::"
                     << "max_integrity_tolerance()="
                     << RefineableElement::max_integrity_tolerance()
                     << std::endl;

        std::ofstream some_file;
        some_file.open("ProblemMesh.dat");
        for (unsigned long n = 0; n < n_node; n++)
        {
          // Get the pointer to the node
          Node* nod_pt = this->node_pt(n);
          // Get the dimension
          unsigned n_dim = nod_pt->ndim();
          // Output the coordinates
          for (unsigned i = 0; i < n_dim; i++)
          {
            some_file << this->node_pt(n)->x(i) << " ";
          }
          some_file << std::endl;
        }
        some_file.close();

        error_stream << "Documented problem mesh in ProblemMesh.dat"
                     << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "Mesh refined: Max. error in integrity check: "
                   << max_error << " is OK" << std::endl;
        oomph_info
          << "i.e. less than RefineableElement::max_integrity_tolerance()="
          << RefineableElement::max_integrity_tolerance() << "\n"
          << std::endl;
      }


      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for (paranoid only) checking of integrity: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

#endif

      // Loop over all elements other than the final level and deactivate the
      // objects, essentially set the pointer that point to nodes that are
      // about to be deleted to NULL. This must take place here because nodes
      // addressed by elements that are dead but still living in the tree might
      // have been made obsolete in the last round of refinement
      //(Not strictly required, as tree structure has not changed, but does no
      // harm)
      for (unsigned long e = 0; e < this->forest_pt()->ntree(); e++)
      {
        this->forest_pt()->tree_pt(e)->traverse_all(&Tree::deactivate_object);
      }

      // Now we can prune the dead nodes from the mesh.
      Vector<Node*> deleted_node_pt = this->prune_dead_nodes();

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for deactivating objects and pruning nodes: "
                   << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Finally: Reorder the nodes within the mesh's node vector
      // to establish a standard ordering regardless of the sequence
      // of mesh refinements -- this is required to allow dump/restart
      // on refined meshes
      this->reorder_nodes();

      if (Global_timings::Doc_comprehensive_timings)
      {
        t_end = TimingHelpers::timer();
        oomph_info << "Time for reordering " << nnode()
                   << " nodes: " << t_end - t_start << std::endl;
        t_start = TimingHelpers::timer();
      }

      // Now we can correct the nodes on boundaries that were hanging that
      // are no longer hanging
      // Only bother if we have more than two dimensions
      if (mesh_dim > 2)
      {
        // Loop over the boundaries
        for (unsigned b = 0; b < n_boundary; b++)
        {
          // Remove deleted nodes from the set
          unsigned n_del = deleted_node_pt.size();
          for (unsigned j = 0; j < n_del; j++)
          {
            hanging_nodes_on_boundary_pt[b].erase(deleted_node_pt[j]);
          }

          // If the nodes that were hanging are still hanging then remove them
          // from the set (note increment is not in for command for efficiencty)
          for (std::set<Node*>::iterator it =
                 hanging_nodes_on_boundary_pt[b].begin();
               it != hanging_nodes_on_boundary_pt[b].end();)
          {
            if ((*it)->is_hanging())
            {
              hanging_nodes_on_boundary_pt[b].erase(it++);
            }
            else
            {
              ++it;
            }
          }

          // Are there any nodes that have changed geometric hanging status
          // on the boundary
          // The slightly painful part is that we must adjust the position
          // via the macro-elements which are only available through the
          // elements and not the nodes.
          if (hanging_nodes_on_boundary_pt[b].size() > 0)
          {
            // If so we loop over all elements adjacent to the boundary
            unsigned n_boundary_element = this->nboundary_element(b);
            for (unsigned e = 0; e < n_boundary_element; ++e)
            {
              // Get a pointer to the element
              FiniteElement* el_pt = this->boundary_element_pt(b, e);

              // Do we have a solid element
              SolidFiniteElement* solid_el_pt =
                dynamic_cast<SolidFiniteElement*>(el_pt);

              // Determine whether there is a macro element
              bool macro_present = (el_pt->macro_elem_pt() != 0);
              // Or a solid macro element
              if (solid_el_pt != 0)
              {
                macro_present |= (solid_el_pt->undeformed_macro_elem_pt() != 0);
              }

              // Only bother to do anything if there is a macro element
              // or undeformed macro element in a SolidElement
              if (macro_present)
              {
                // Loop over the nodes
                // ALH: (could optimise to only loop over
                // node associated with the boundary with more effort)
                unsigned n_el_node = el_pt->nnode();
                for (unsigned n = 0; n < n_el_node; n++)
                {
                  // Cache pointer to the node
                  Node* nod_pt = el_pt->node_pt(n);
                  if (nod_pt->is_on_boundary(b))
                  {
                    // Is the Node in our set
                    std::set<Node*>::iterator it =
                      hanging_nodes_on_boundary_pt[b].find(nod_pt);
                    // If we have found the Node then update the position
                    // to be consistent with the macro-element representation
                    if (it != hanging_nodes_on_boundary_pt[b].end())
                    {
                      // Specialise local and global coordinates to 3D
                      // because there is only a problem in 3D.
                      Vector<double> s(3), x(3);
                      // Find the local coordinate of the ndoe
                      el_pt->local_coordinate_of_node(n, s);
                      // Find the number of time history values
                      const unsigned ntstorage = nod_pt->ntstorage();

                      // Do we have a solid node
                      SolidNode* solid_node_pt =
                        dynamic_cast<SolidNode*>(nod_pt);
                      if (solid_node_pt)
                      {
                        // Assign Lagrangian coordinates from undeformed
                        // macro element (if it has one -- get_x_and_xi()
                        // does "the right thing" anyway. Leave actual
                        // nodal positions alone -- we're doing a solid
                        // mechanics problem and once we're going
                        // the nodal positions are always computed, never
                        // (automatically) reset to macro-element based
                        // positions; not even on pinned boundaries
                        // because the user may have other ideas about where
                        // these should go -- pinning means "don't touch the
                        // value", not "leave where the macro-element thinks
                        // it should be"
                        Vector<double> x_fe(3), xi(3), xi_fe(3);
                        solid_el_pt->get_x_and_xi(s, x_fe, x, xi_fe, xi);
                        for (unsigned i = 0; i < 3; i++)
                        {
                          solid_node_pt->xi(i) = xi[i];
                        }
                      }
                      else
                      {
                        // Set position and history values from the
                        // macro-element representation
                        for (unsigned t = 0; t < ntstorage; t++)
                        {
                          // Get the history value from the macro element
                          el_pt->get_x(t, s, x);

                          // Set the coordinate to that of the macroelement
                          // representation
                          for (unsigned i = 0; i < 3; i++)
                          {
                            nod_pt->x(t, i) = x[i];
                          }
                        }
                      } // End of non-solid node case

                      // Now remove the node from the list
                      hanging_nodes_on_boundary_pt[b].erase(it);
                      // If there are no Nodes left then exit the loops
                      if (hanging_nodes_on_boundary_pt[b].size() == 0)
                      {
                        e = n_boundary_element;
                        break;
                      }
                    }
                  }
                }
              } // End of macro element case
            }
          }
        }
      } // End of case when we have fixed nodal positions

      // Final doc
      //-----------
      if (doc_info.is_doc_enabled())
      {
        // Doc the boundary conditions ('0' for non-existent, '1' for free,
        //----------------------------------------------------------------
        // '2' for pinned -- ideal for tecplot scatter sizing.
        //----------------------------------------------------
        // num_tree_nodes=tree_nodes_pt.size();

        // Determine maximum number of values at any node in this type of
        // element
        RefineableElement* el_pt = tree_nodes_pt[0]->object_pt();
        // Initalise max_nval
        unsigned max_nval = 0;
        for (unsigned n = 0; n < el_pt->nnode(); n++)
        {
          if (el_pt->node_pt(n)->nvalue() > max_nval)
          {
            max_nval = el_pt->node_pt(n)->nvalue();
          }
        }

        // Open the output file
        std::ostringstream fullname;
        std::ofstream bcs_file;
        fullname.str("");
        fullname << doc_info.directory() << "/bcs" << doc_info.number()
                 << ".dat";
        bcs_file.open(fullname.str().c_str());

        // Loop over elements
        for (unsigned long e = 0; e < num_tree_nodes; e++)
        {
          el_pt = tree_nodes_pt[e]->object_pt();
          // Loop over nodes in element
          unsigned n_nod = el_pt->nnode();
          for (unsigned n = 0; n < n_nod; n++)
          {
            // Get pointer to the node
            Node* nod_pt = el_pt->node_pt(n);
            // Find the dimension of the node
            unsigned n_dim = nod_pt->ndim();
            // Write the nodal coordinates to the file
            for (unsigned i = 0; i < n_dim; i++)
            {
              bcs_file << nod_pt->x(i) << " ";
            }

            // Loop over all values in this element
            for (unsigned i = 0; i < max_nval; i++)
            {
              // Value exists at this node:
              if (i < nod_pt->nvalue())
              {
                bcs_file << " " << 1 + nod_pt->is_pinned(i);
              }
              // ...if not just dump out a zero
              else
              {
                bcs_file << " 0 ";
              }
            }
            bcs_file << std::endl;
          }
        }
        bcs_file.close();

        // Doc all nodes
        //---------------
        std::ofstream all_nodes_file;
        fullname.str("");
        fullname << doc_info.directory() << "/all_nodes" << doc_info.number()
                 << ".dat";
        all_nodes_file.open(fullname.str().c_str());

        all_nodes_file << "ZONE \n";

        // Need to recompute the number of nodes since it may have
        // changed during mesh refinement/unrefinement
        n_node = this->nnode();
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);
          unsigned n_dim = nod_pt->ndim();
          for (unsigned i = 0; i < n_dim; i++)
          {
            all_nodes_file << this->node_pt(n)->x(i) << " ";
          }
          all_nodes_file << std::endl;
        }

        all_nodes_file.close();


        // Doc all hanging nodes:
        //-----------------------
        std::ofstream some_file;
        fullname.str("");
        fullname << doc_info.directory() << "/all_hangnodes"
                 << doc_info.number() << ".dat";
        some_file.open(fullname.str().c_str());
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);

          if (nod_pt->is_hanging())
          {
            unsigned n_dim = nod_pt->ndim();
            for (unsigned i = 0; i < n_dim; i++)
            {
              some_file << nod_pt->x(i) << " ";
            }

            // ALH: Added this to stop Solid problems seg-faulting
            if (this->node_pt(n)->nvalue() > 0)
            {
              some_file << " " << nod_pt->raw_value(0);
            }
            some_file << std::endl;
          }
        }
        some_file.close();

        // Doc all hanging nodes and their masters
        // View with QHangingNodesWithMasters.mcr
        fullname.str("");
        fullname << doc_info.directory() << "/geometric_hangnodes_withmasters"
                 << doc_info.number() << ".dat";
        some_file.open(fullname.str().c_str());
        for (unsigned long n = 0; n < n_node; n++)
        {
          Node* nod_pt = this->node_pt(n);
          if (nod_pt->is_hanging())
          {
            unsigned n_dim = nod_pt->ndim();
            unsigned nmaster = nod_pt->hanging_pt()->nmaster();
            some_file << "ZONE I=" << nmaster + 1 << std::endl;
            for (unsigned i = 0; i < n_dim; i++)
            {
              some_file << nod_pt->x(i) << " ";
            }
            some_file << " 2 " << std::endl;

            for (unsigned imaster = 0; imaster < nmaster; imaster++)
            {
              Node* master_nod_pt =
                nod_pt->hanging_pt()->master_node_pt(imaster);
              unsigned n_dim = master_nod_pt->ndim();
              for (unsigned i = 0; i < n_dim; i++)
              {
                some_file << master_nod_pt->x(i) << " ";
              }
              some_file << " 1 " << std::endl;
            }
          }
        }
        some_file.close();

        // Doc all hanging nodes and their masters
        // View with QHangingNodesWithMasters.mcr
        for (unsigned i = 0; i < ncont_interpolated_values; i++)
        {
          fullname.str("");
          fullname << doc_info.directory()
                   << "/nonstandard_hangnodes_withmasters" << i << "_"
                   << doc_info.number() << ".dat";
          some_file.open(fullname.str().c_str());
          unsigned n_nod = this->nnode();
          for (unsigned long n = 0; n < n_nod; n++)
          {
            Node* nod_pt = this->node_pt(n);
            if (nod_pt->is_hanging(i))
            {
              if (nod_pt->hanging_pt(i) != nod_pt->hanging_pt())
              {
                unsigned nmaster = nod_pt->hanging_pt(i)->nmaster();
                some_file << "ZONE I=" << nmaster + 1 << std::endl;
                unsigned n_dim = nod_pt->ndim();
                for (unsigned j = 0; j < n_dim; j++)
                {
                  some_file << nod_pt->x(j) << " ";
                }
                some_file << " 2 " << std::endl;
                for (unsigned imaster = 0; imaster < nmaster; imaster++)
                {
                  Node* master_nod_pt =
                    nod_pt->hanging_pt(i)->master_node_pt(imaster);
                  unsigned n_dim = master_nod_pt->ndim();
                  for (unsigned j = 0; j < n_dim; j++)
                  {
                    //               some_file << master_nod_pt->x(i) << " ";
                  }
                  some_file << " 1 " << std::endl;
                }
              }
            }
          }
          some_file.close();
        }

      } // End of documentation
    } // End if (this->nelement()>0)

    ////BENFLAG: Check that all the nodes belong to their update elements
    // std::cout << "p_adapt_mesh(): Checking stuff works!" << std::endl;
    // for(unsigned j=0; j<this->nnode(); j++)
    // {
    //  MacroElementNodeUpdateNode* macro_nod_pt =
    //  dynamic_cast<MacroElementNodeUpdateNode*>(this->node_pt(j));
    //  if(macro_nod_pt!=0)
    //   {
    //    bool big_problem = true;
    //    std::cout << "Node " << macro_nod_pt << " at [ " << macro_nod_pt->x(0)
    //    << ", " << macro_nod_pt->x(1) << " ]" << std::endl; FiniteElement*
    //    up_el_pt =
    //    dynamic_cast<FiniteElement*>(macro_nod_pt->node_update_element_pt());
    //    for(unsigned l=0; l<up_el_pt->nnode(); l++)
    //     {
    //      if(up_el_pt->node_pt(l)==macro_nod_pt)
    //       {
    //        big_problem = false;
    //        break;
    //       }
    //     }
    //    if(big_problem)
    //     {
    //      std::cout << "  This node doesn't exist in it's update element!" <<
    //      std::endl;
    //     }
    //   }
    // }

#ifdef OOMPH_HAS_MPI

    // Now (re-)classify halo and haloed nodes and synchronise hanging
    // nodes
    if (this->is_mesh_distributed())
    {
      classify_halo_and_haloed_nodes(doc_info, doc_info.is_doc_enabled());
    }

#endif
  }

  //========================================================================
  /// p-unrefine mesh uniformly
  /// Unlike in h-refinement, we can simply p-unrefine each element in the mesh
  //========================================================================
  void TreeBasedRefineableMeshBase::p_unrefine_uniformly(DocInfo& doc_info)
  {
    // Select all elements for unrefinement
    unsigned long Nelement = this->nelement();
    for (unsigned long e = 0; e < Nelement; e++)
    {
      // Get pointer to p-refineable element
      PRefineableElement* el_pt =
        dynamic_cast<PRefineableElement*>(this->element_pt(e));
      // Mark for p-refinement if possible. If not then p_adapt_mesh() will
      // report the error.
      if (el_pt != 0)
      {
        el_pt->select_for_p_unrefinement();
      }
    }

    // Do the actual mesh adaptation
    p_adapt_mesh(doc_info);
  }

  //========================================================================
  /// p-refine mesh by refining the elements identified by their numbers.
  //========================================================================
  void TreeBasedRefineableMeshBase::p_refine_selected_elements(
    const Vector<unsigned>& elements_to_be_refined)
  {
#ifdef OOMPH_HAS_MPI
    if (this->is_mesh_distributed())
    {
      std::ostringstream warn_stream;
      warn_stream << "You are attempting to refine selected elements of a "
                  << std::endl
                  << "distributed mesh. This may have undesired effects."
                  << std::endl;

      OomphLibWarning(warn_stream.str(),
                      "TreeBasedRefineableMeshBase::refine_selected_elements()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Select elements for refinement
    unsigned long nref = elements_to_be_refined.size();
    for (unsigned long e = 0; e < nref; e++)
    {
      // Get pointer to p-refineable element
      PRefineableElement* el_pt = dynamic_cast<PRefineableElement*>(
        this->element_pt(elements_to_be_refined[e]));
      // Mark for p-refinement if possible. If not then p_adapt_mesh() will
      // report the error.
      if (el_pt != 0)
      {
        el_pt->select_for_p_refinement();
      }
    }

    // Do the actual mesh adaptation
    p_adapt_mesh();
  }

  //========================================================================
  /// p-refine mesh by refining the elements identified by their pointers.
  //========================================================================
  void TreeBasedRefineableMeshBase::p_refine_selected_elements(
    const Vector<PRefineableElement*>& elements_to_be_refined_pt)
  {
#ifdef OOMPH_HAS_MPI
    if (this->is_mesh_distributed())
    {
      std::ostringstream warn_stream;
      warn_stream << "You are attempting to refine selected elements of a "
                  << std::endl
                  << "distributed mesh. This may have undesired effects."
                  << std::endl;

      OomphLibWarning(warn_stream.str(),
                      "TreeBasedRefineableMeshBase::refine_selected_elements()",
                      OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Select elements for refinement
    unsigned long nref = elements_to_be_refined_pt.size();
    for (unsigned long e = 0; e < nref; e++)
    {
      elements_to_be_refined_pt[e]->select_for_p_refinement();
    }

    // Do the actual mesh adaptation
    p_adapt_mesh();
  }

} // namespace oomph

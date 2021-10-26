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

#include "problem.h"
#include "fsi.h"
#include "elastic_problems.h"


namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Setup IC problem by:
  /// - Pinning all nodal values in the mesh
  /// - Pinning the internal data of all elements.
  /// - Freeing/unpinnning all positional data.
  /// - Flushing the pointers to the elements' external data.
  /// - Setting the pointer to the IC object for all elements to
  ///   ensure that the correct residuals/Jacobians are computed.
  //======================================================================
  void SolidICProblem::setup_problem()
  {
    // Find out how many nodes there are
    unsigned long n_node = mesh_pt()->nnode();

    // Loop over all the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cast to an elastic node
      SolidNode* node_pt = dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n));

#ifdef PARANOID
      if (node_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidNode\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get spatial dimension of node
      unsigned ndim = node_pt->ndim();

      // Find out how many positional dofs there are
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(0));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned ntype = elem_pt->nnodal_position_type();

      // Loop over coordinate directions
      for (unsigned i = 0; i < ndim; i++)
      {
        // Loop over type of dof
        for (unsigned k = 0; k < ntype; k++)
        {
          // Unpin them
          node_pt->unpin_position(k, i);
        }
      }

      // Loop over nodal values
      unsigned nval = node_pt->nvalue();
      for (unsigned ival = 0; ival < nval; ival++)
      {
        // Pin them
        node_pt->pin(ival);
      }
    }


    // Loop over the elements
    unsigned Nelement = mesh_pt()->nelement();
    for (unsigned i = 0; i < Nelement; i++)
    {
      // Cast to proper element type
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(i));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // Set flag for setting initial condition
      elem_pt->solid_ic_pt() = IC_pt;

      // We've backed up the element's external data: Flush it
      elem_pt->flush_external_data();

      // IF it's an FSI wall element then kill external stuff
      if (FSIWallElement* fsi_elem_pt = dynamic_cast<FSIWallElement*>(elem_pt))
      {
        fsi_elem_pt->exclude_external_load_data();
      }

      // Find out number of internal data
      unsigned nint = elem_pt->ninternal_data();

      // Loop over internal data
      for (unsigned j = 0; j < nint; j++)
      {
        Data* data_pt = elem_pt->internal_data_pt(j);

        // Loop over internal values
        unsigned nval = data_pt->nvalue();
        for (unsigned ival = 0; ival < nval; ival++)
        {
          // Pin internal values
          data_pt->pin(ival);
        }
      }

#ifdef PARANOID
      // Is there internal solid data
      if (elem_pt->has_internal_solid_data())
      {
        std::string error_message =
          "Automatic assignment of initial conditions doesn't work yet\n";
        error_message +=
          "for elasticity elements with internal solid dofs (pressures)\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      //    for(unsigned i=0;i<nint_solid;i++)
      //     {
      //      Data* data_pt=elem_pt->internal_solid_data_pt(i);

      //      // Loop over values
      //      unsigned nval=data_pt->nvalue();
      //      for (unsigned ival=0;ival<nval;ival++)
      //       {
      //        // Pin internal values
      //        data_pt->pin(ival);
      //       }
      //     }
    }

    // Setup equation numbers for IC problem
    oomph_info << "# of dofs for wall initial guess" << assign_eqn_numbers()
               << std::endl;
  }


  //======================================================================
  /// Backup pinned status of all data associated with the mesh.
  /// Also backup the (pointers to the) elements' external data.
  //======================================================================
  void SolidICProblem::backup_original_state()
  {
    // Find out how many nodes there are
    unsigned long n_node = mesh_pt()->nnode();

    // Flush vector which holds backup of pinned status
    Backup_pinned.clear();

    // Flush vector which holds vectors with backup of (pointers to) external
    // data
    Backup_ext_data.clear();

    // Loop over all the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cast to an elastic node
      SolidNode* node_pt = dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n));

#ifdef PARANOID
      if (node_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidNode\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get spatial dimension of node
      unsigned ndim = node_pt->ndim();

      // Find out how many positional dofs there are
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(0));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned ntype = elem_pt->nnodal_position_type();


      // Loop over coordinate directions
      for (unsigned i = 0; i < ndim; i++)
      {
        // Loop over type of dof
        for (unsigned k = 0; k < ntype; k++)
        {
          // Backup pinned status
          Backup_pinned.push_back(node_pt->position_is_pinned(k, i));
        }
      }

      // Loop over nodal values
      unsigned nval = node_pt->nvalue();
      for (unsigned ival = 0; ival < nval; ival++)
      {
        // Backup pinned status
        Backup_pinned.push_back(node_pt->is_pinned(ival));
      }
    }

    // Loop over the elements
    unsigned Nelement = mesh_pt()->nelement();
    Backup_ext_data.resize(Nelement);
    for (unsigned i = 0; i < Nelement; i++)
    {
      // Cast to proper element type
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(i));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Find out number of external data
      unsigned next = elem_pt->nexternal_data();
      Backup_ext_data[i].resize(next);

      // Loop over external data
      for (unsigned j = 0; j < next; j++)
      {
        Data* data_pt = elem_pt->external_data_pt(j);

        // Backup the pointer to external data
        Backup_ext_data[i][j] = data_pt;
      }

      // Find out number of internal data
      unsigned nint = elem_pt->ninternal_data();

      // Loop over internal data
      for (unsigned j = 0; j < nint; j++)
      {
        Data* data_pt = elem_pt->internal_data_pt(j);

        // Loop over internal values
        unsigned nval = data_pt->nvalue();
        for (unsigned ival = 0; ival < nval; ival++)
        {
          // Backup pinned status
          Backup_pinned.push_back(data_pt->is_pinned(ival));
        }
      }


#ifdef PARANOID
      // If there is internal solid data, complain
      if (elem_pt->has_internal_solid_data())
      {
        std::string error_message =
          "Automatic assignment of initial conditions doesn't work yet\n";
        error_message +=
          "for elasticity elements with internal solid dofs (pressures)\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      //    for(unsigned i=0;i<nint_solid;i++)
      //     {
      //      Data* data_pt=elem_pt->internal_solid_data_pt(i);

      //      // Loop over values
      //      unsigned nval=data_pt->nvalue();
      //      for (unsigned ival=0;ival<nval;ival++)
      //       {
      //        // Backup pinned status
      //        Backup_pinned.push_back(data_pt->is_pinned(ival));
      //       }
      //     }
    }

    // Record number of dofs whose status was pinned
    // oomph_info << "Number of backed up values " << Backup_pinned.size() <<
    // std::endl;
  }


  //======================================================================
  /// Reset pinned status of all data and re-instate the pointers
  /// to the elements' external data.
  //======================================================================
  void SolidICProblem::reset_original_state()
  {
    // Find out how many nodes there are
    unsigned long n_node = mesh_pt()->nnode();

    // Initialise counter for backed up dofs
    unsigned count = 0;

    // Loop over all the nodes
    for (unsigned n = 0; n < n_node; n++)
    {
      // Cast to an elastic node
      SolidNode* node_pt = dynamic_cast<SolidNode*>(mesh_pt()->node_pt(n));

#ifdef PARANOID
      if (node_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidNode\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get spatial dimension of node
      unsigned ndim = node_pt->ndim();

      // Find out how many positional dofs there are
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(0));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned ntype = elem_pt->nnodal_position_type();

      // Loop over coordinate directions
      for (unsigned i = 0; i < ndim; i++)
      {
        // Loop over type of dof
        for (unsigned k = 0; k < ntype; k++)
        {
          // Reset pinned status (positional dofs were all unpinned)
          if (Backup_pinned[count])
          {
            node_pt->pin_position(k, i);
          }
          count++;
        }
      }

      // Loop over nodal values
      unsigned nval = node_pt->nvalue();
      for (unsigned ival = 0; ival < nval; ival++)
      {
        // Reset pinned status (nodal values were all pinned)
        if (Backup_pinned[count])
        {
          node_pt->unpin(ival);
        }
      }
    }


    // Loop over the elements
    unsigned Nelement = mesh_pt()->nelement();
    for (unsigned i = 0; i < Nelement; i++)
    {
      // Cast to proper element type
      SolidFiniteElement* elem_pt =
        dynamic_cast<SolidFiniteElement*>(mesh_pt()->element_pt(i));

#ifdef PARANOID
      if (elem_pt == 0)
      {
        throw OomphLibError("Wasn't able to cast to SolidFiniteElement\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Switch back to normal Jacobian
      dynamic_cast<SolidFiniteElement*>(elem_pt)
        ->disable_solve_for_consistent_newmark_accel();

      // Switch off flag for setting initial condition
      dynamic_cast<SolidFiniteElement*>(elem_pt)->solid_ic_pt() = 0;

      // IF it's an FSI wall element then turn on external stuff again
      if (FSIWallElement* fsi_elem_pt = dynamic_cast<FSIWallElement*>(elem_pt))
      {
        fsi_elem_pt->include_external_load_data();
      }

      // Find out number of external data
      unsigned next = Backup_ext_data[i].size();

      // Loop over external data
      for (unsigned j = 0; j < next; j++)
      {
        // Backed up external data
        Data* data_pt = Backup_ext_data[i][j];

        // Add external data
        elem_pt->add_external_data(data_pt);
      }

      // Find out number of internal data
      unsigned nint = elem_pt->ninternal_data();

      // Loop over internal data
      for (unsigned j = 0; j < nint; j++)
      {
        Data* data_pt = elem_pt->internal_data_pt(j);

        // Loop over internal values
        unsigned nval = data_pt->nvalue();
        for (unsigned ival = 0; ival < nval; ival++)
        {
          // Restore pinned status (values were all pinned)
          if (!Backup_pinned[count])
          {
            data_pt->unpin(ival);
          }
        }
      }


#ifdef PARANOID
      // If there is internal solid data, complain
      if (elem_pt->has_internal_solid_data())
      {
        std::string error_message =
          "Automatic assignment of initial conditions doesn't work yet\n";
        error_message +=
          "for elasticity elements with internal solid dofs (pressures)\n";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      //    for(unsigned i=0;i<nint_solid;i++)
      //     {
      //      Data* data_pt=elem_pt->internal_solid_data_pt(i);

      //      // Loop over values
      //      unsigned nval=data_pt->nvalue();
      //      for (unsigned ival=0;ival<nval;ival++)
      //       {
      //        // Restore pinned status (values were all pinned)
      //        if (!Backup_pinned[count])
      //         {
      //          data_pt->unpin(ival);
      //         }
      //       }
      //     }
    }

    // Check number of recovered pinned values
    // oomph_info << "Recovered pin values " << count << std::endl;

    // Flush vector which holds backup of pinned status
    Backup_pinned.clear();

    // Flush vector which holds vectors with backup of (pointers to) external
    // data
    Backup_ext_data.clear();
  }


  //======================================================================
  /// IC problem for wall: Deform wall into the static initial shape
  /// described by the IC object at given time.
  //======================================================================
  void SolidICProblem::set_static_initial_condition(
    Problem* problem_pt,
    Mesh* wall_mesh_pt,
    SolidInitialCondition* ic_pt,
    const double& time)
  {
    // Tell this sub-problem it is distributed if the main problem is
    // distributed
#ifdef OOMPH_HAS_MPI
    if (problem_pt->problem_has_been_distributed())
    {
      Problem_has_been_distributed = true;
    }
    // This (sub-)problem needs to know about the oomph communicator
    delete Communicator_pt;
    Communicator_pt = new OomphCommunicator(problem_pt->communicator_pt());
#endif


    // Backup value of time
    double backup_time = 0.0;

    // Set value of time for IC object (needs to be backed up and
    // restored since it might be pointed to by other objects)
    TimeStepper* timestepper_pt = ic_pt->geom_object_pt()->time_stepper_pt();
    if (timestepper_pt != 0)
    {
      backup_time = timestepper_pt->time_pt()->time();
      timestepper_pt->time_pt()->time() = time;
    }

    // Delete dummy mesh
    delete mesh_pt();

    // Set pointer to mesh
    mesh_pt() = wall_mesh_pt;

    // Set pointer to initial condition object
    IC_pt = ic_pt;

    // Backup the pinned status of all dofs and remove external data
    // of all elements
    backup_original_state();

    // Now alter the pinned status so that the IC problem for the
    // positional variables can be solved; setup equation numbering
    // scheme
    setup_problem();

    // Assign displacements
    IC_pt->ic_time_deriv() = 0;

    // Solve the problem for initial shape
    newton_solve();

    // Impulsive start
    assign_initial_values_impulsive();

    // Reset the pinned status and re-attach the external data to the elements
    reset_original_state();

    // Reset time
    if (timestepper_pt != 0)
    {
      timestepper_pt->time_pt()->time() = backup_time;
    }

    // Set pointer to dummy mesh so there's something that can be deleted
    // when static problem finally goes out of scope.
    mesh_pt() = new DummyMesh;

    // We have temporarily over-written equation numbers -- need
    // to reset them now
    oomph_info << "Number of equations in big problem: "
               << problem_pt->assign_eqn_numbers() << std::endl;
  }

} // namespace oomph

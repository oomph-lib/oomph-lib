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
#include "sample_point_container.h"


namespace oomph
{
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                RefineableBin class
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //==============================================================================
  /// Destructor
  //==============================================================================
  RefineableBin::~RefineableBin()
  {
    if (Sub_bin_array_pt != 0)
    {
      delete Sub_bin_array_pt;
    }
    Sub_bin_array_pt = 0;

    if (Sample_point_pt != 0)
    {
      unsigned n = Sample_point_pt->size();
      for (unsigned i = 0; i < n; i++)
      {
        delete (*Sample_point_pt)[i];
      }
      delete Sample_point_pt;
    }
  }


  //==============================================================================
  /// Compute total number of sample points recursively
  //==============================================================================
  unsigned RefineableBin::total_number_of_sample_points_computed_recursively()
    const
  {
    unsigned count = 0;

    // Recurse?
    if (Sub_bin_array_pt != 0)
    {
      count =
        Sub_bin_array_pt->total_number_of_sample_points_computed_recursively();
    }
    else
    {
      if (Sample_point_pt != 0)
      {
        count = Sample_point_pt->size();
      }
    }
    return count;
  }


  //==============================================================================
  /// Function called for making a sub bin array in a given RefineableBin.
  /// Pass the vector of min and max coordinates of the NEW bin array.
  //==============================================================================
  void RefineableBin::make_sub_bin_array(
    const Vector<std::pair<double, double>>& min_and_max_coordinates)
  {
    // Setup parameters for sub-bin
    RefineableBinArrayParameters* ref_bin_array_parameters_pt =
      new RefineableBinArrayParameters(Bin_array_pt->mesh_pt());

    // Pass coordinates and dimensions
    ref_bin_array_parameters_pt->min_and_max_coordinates() =
      min_and_max_coordinates;

    ref_bin_array_parameters_pt->dimensions_of_bin_array() =
      Bin_array_pt->dimensions_of_bin_array();

    // Eulerian coordinates or zeta?
    if (Bin_array_pt->use_eulerian_coordinates_during_setup())
    {
      ref_bin_array_parameters_pt
        ->enable_use_eulerian_coordinates_during_setup();
    }
    else
    {
      ref_bin_array_parameters_pt
        ->disable_use_eulerian_coordinates_during_setup();
    }


#ifdef OOMPH_HAS_MPI

    // How do we handle halo elements?
    if (Bin_array_pt->ignore_halo_elements_during_locate_zeta_search())
    {
      ref_bin_array_parameters_pt
        ->enable_ignore_halo_elements_during_locate_zeta_search();
    }
    else
    {
      ref_bin_array_parameters_pt
        ->disable_ignore_halo_elements_during_locate_zeta_search();
    }

#endif

    // "Measure of" number of sample points per element
    ref_bin_array_parameters_pt->nsample_points_generated_per_element() =
      Bin_array_pt->nsample_points_generated_per_element();


    // Is it recursive?
    if (Bin_array_pt->bin_array_is_recursive())
    {
      ref_bin_array_parameters_pt->enable_bin_array_is_recursive();
    }
    else
    {
      ref_bin_array_parameters_pt->disable_bin_array_is_recursive();
    }

    // Depth
    ref_bin_array_parameters_pt->depth() = Bin_array_pt->depth() + 1;

    // Max. depth
    ref_bin_array_parameters_pt->max_depth() = Bin_array_pt->max_depth();

    // Max. number of sample points before it's subdivided
    ref_bin_array_parameters_pt->max_number_of_sample_point_per_bin() =
      Bin_array_pt->max_number_of_sample_point_per_bin();

    // Root bin array
    ref_bin_array_parameters_pt->root_bin_array_pt() =
      Bin_array_pt->root_bin_array_pt();

    // We first construct a new bin array, providing the right parameters
    BinArrayParameters* bin_array_parameters_pt = ref_bin_array_parameters_pt;
    Sub_bin_array_pt = new RefineableBinArray(bin_array_parameters_pt);
    delete ref_bin_array_parameters_pt;

    // Fill it
    Sub_bin_array_pt->fill_bin_array(*Sample_point_pt);

    // Now deleting Sample_point_pt, we no longer need it; note that
    // sample points themselves stay alive!
    delete Sample_point_pt;
    Sample_point_pt = 0;
  }


  //============================================================================
  /// Output bin; x,[y,[z]],n_sample_points.
  //============================================================================
  void RefineableBin::output(std::ofstream& outfile, const bool& don_t_recurse)
  {
    // Recurse?
    if ((Sub_bin_array_pt != 0) && (!don_t_recurse))
    {
      Sub_bin_array_pt->output_bins(outfile);
    }
    else
    {
      unsigned n_lagr = Bin_array_pt->ndim_zeta();
      Vector<std::pair<double, double>> min_and_max_coordinates(n_lagr);
      get_bin_boundaries(min_and_max_coordinates);

      // How many sample points do we have in this bin?
      unsigned n_sample_points = 0;
      if (Sample_point_pt != 0)
      {
        n_sample_points = Sample_point_pt->size();
      }

      switch (n_lagr)
      {
        case 1:
          outfile << "ZONE I=2\n"
                  << min_and_max_coordinates[0].first << " " << n_sample_points
                  << std::endl
                  << min_and_max_coordinates[0].second << " " << n_sample_points
                  << std::endl;
          break;

        case 2:

          outfile << "ZONE I=2, J=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " " << n_sample_points
                  << "\n";
          break;

        case 3:


          outfile << "ZONE I=2, J=2, K=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " " << n_sample_points
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " " << n_sample_points
                  << "\n";

          break;

        default:

          oomph_info << "n_lagr=" << n_lagr << std::endl;
          throw OomphLibError("Wrong dimension",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  //============================================================================
  /// Output bin; x,[y,[z]]
  //============================================================================
  void RefineableBin::output_bin_vertices(std::ofstream& outfile)
  {
    // Recurse?
    if (Sub_bin_array_pt != 0)
    {
      Sub_bin_array_pt->output_bin_vertices(outfile);
    }
    else
    {
      unsigned n_lagr = Bin_array_pt->ndim_zeta();
      Vector<std::pair<double, double>> min_and_max_coordinates(n_lagr);
      get_bin_boundaries(min_and_max_coordinates);

      switch (n_lagr)
      {
        case 1:
          outfile << "ZONE I=2\n"
                  << min_and_max_coordinates[0].first << std::endl
                  << min_and_max_coordinates[0].second << std::endl;
          break;

        case 2:

          outfile << "ZONE I=2, J=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << "\n";
          break;

        case 3:


          outfile << "ZONE I=2, J=2, K=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n";

          break;

        default:

          oomph_info << "n_lagr=" << n_lagr << std::endl;
          throw OomphLibError("Wrong dimension",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  //==============================================================================
  /// Add a SamplePoint* to a RefineableBin object.
  //==============================================================================
  void RefineableBin::add_sample_point(SamplePoint* new_sample_point_pt,
                                       const Vector<double>& zeta_coordinates)
  {
    // If the bin is a "leaf" (ie no sub bin array)
    if (Sub_bin_array_pt == 0)
    {
      // if there is no Sample_point_pt create it
      if (Sample_point_pt == 0)
      {
        Sample_point_pt = new Vector<SamplePoint*>;
      }
      this->Sample_point_pt->push_back(new_sample_point_pt);

      // If we are recursive (ie not at the maximum depth or not
      // in fill bin by diffusion configuration) and if the number
      // of elements there are in the RefineableBin is bigger than the maximum
      // one
      if ((Bin_array_pt->bin_array_is_recursive()) &&
          (Sample_point_pt->size() >
           Bin_array_pt->max_number_of_sample_point_per_bin()) &&
          (Bin_array_pt->depth() < Bin_array_pt->max_depth()))
      {
        // Get min and max coordinates of current bin...
        Vector<std::pair<double, double>> min_and_max_coordinates(
          Bin_array_pt->ndim_zeta());
        get_bin_boundaries(min_and_max_coordinates);

        // ...and use them as the boundaries for new sub-bin-array
        // (this transfers all the new points into the new sub-bin-array
        this->make_sub_bin_array(min_and_max_coordinates);
      }
    }
    else // if the bin has a sub bin array
    {
      // we call the corresponding method of the sub bin array
      this->Sub_bin_array_pt->add_sample_point(new_sample_point_pt,
                                               zeta_coordinates);
    }
  }


  //==============================================================================
  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  //==============================================================================
  void RefineableBin::locate_zeta(const Vector<double>& zeta,
                                  GeomObject*& sub_geom_object_pt,
                                  Vector<double>& s)
  {
    // Haven't found zeta yet!
    sub_geom_object_pt = 0;

    // Descend?
    if (Sub_bin_array_pt != 0)
    {
      Sub_bin_array_pt->locate_zeta(zeta, sub_geom_object_pt, s);
    }
    else
    {
      // Do we have to look into any sample points in this bin?
      // NOTE: There's some slight potential for overlap/duplication
      // because we always search through all the sample points in a bin
      // (unless we find the required point in which case we stop).
      bool do_it = true;
      if (
        Bin_array_pt->root_bin_array_pt()
          ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level() <
        Bin_array_pt->root_bin_array_pt()
          ->first_sample_point_to_actually_lookup_during_locate_zeta())
      {
        // oomph_info << "Not doing it (ref-bin) because counter less than
        // first" << std::endl;
        do_it = false;
      }
      if (
        Bin_array_pt->root_bin_array_pt()
          ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level() >
        Bin_array_pt->root_bin_array_pt()
          ->last_sample_point_to_actually_lookup_during_locate_zeta())
      {
        // oomph_info << "Not doing it (ref-bin) because counter more than last"
        // << std::endl;
        do_it = false;
      }
      double max_search_radius =
        Bin_array_pt->root_bin_array_pt()->max_search_radius();
      bool dont_do_it_because_of_radius = false;
      if (max_search_radius < DBL_MAX)
      {
        // Base "radius" of bin on centre of gravity
        unsigned n = zeta.size();
        double dist_squared = 0.0;
        double cog = 0.0;
        double aux = 0.0;
        Vector<std::pair<double, double>> min_and_max_coordinates(n);
        get_bin_boundaries(min_and_max_coordinates);
        for (unsigned i = 0; i < n; i++)
        {
          cog = 0.5 * (min_and_max_coordinates[i].first +
                       min_and_max_coordinates[i].second);
          aux = (cog - zeta[i]);
          dist_squared += aux * aux;
        }
        if (dist_squared > max_search_radius * max_search_radius)
        {
          do_it = false;
          dont_do_it_because_of_radius = true;
        }
      }

      if (!do_it)
      {
        if (!dont_do_it_because_of_radius)
        {
          // Skip all the sample points in this bin
          Bin_array_pt->root_bin_array_pt()
            ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level() +=
            Sample_point_pt->size();
        }
        return;
      }


      // Now search through (at most) all the sample points in this bin
      unsigned n_sample_point = Sample_point_pt->size();
      unsigned i = 0;
      while ((i < n_sample_point) && (sub_geom_object_pt == 0))
      {
        // Get the corresponding finite element
        FiniteElement* el_pt = Bin_array_pt->mesh_pt()->finite_element_pt(
          (*Sample_point_pt)[i]->element_index_in_mesh());

#ifdef OOMPH_HAS_MPI
        // We only look at the sample point if it isn't halo
        // if we are set up to ignore the halo elements
        if ((Bin_array_pt->ignore_halo_elements_during_locate_zeta_search()) &&
            (el_pt->is_halo()))
        {
          // Halo
        }
        else
        {
#endif
          // Provide initial guess for Newton search using local coordinate
          // of sample point
          bool use_equally_spaced_interior_sample_points =
            SamplePointContainer::Use_equally_spaced_interior_sample_points;
          unsigned j = (*Sample_point_pt)[i]->sample_point_index_in_element();
          el_pt->get_s_plot(
            j,
            Bin_array_pt->nsample_points_generated_per_element(),
            s,
            use_equally_spaced_interior_sample_points);


          // History of sample points visited
          if (BinArray::Visited_sample_points_file.is_open())
          {
            unsigned cached_dim_zeta = Bin_array_pt->ndim_zeta();
            Vector<double> zeta_sample(cached_dim_zeta);
            if (Bin_array_pt->use_eulerian_coordinates_during_setup())
            {
              el_pt->interpolated_x(s, zeta_sample);
            }
            else
            {
              el_pt->interpolated_zeta(s, zeta_sample);
            }
            double dist = 0.0;
            for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
            {
              BinArray::Visited_sample_points_file << zeta_sample[ii] << " ";
              dist +=
                (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
            }
            BinArray::Visited_sample_points_file
              << Bin_array_pt->root_bin_array_pt()
                   ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
              << " " << sqrt(dist) << std::endl;
          }


          // Bump counter
          Bin_array_pt->root_bin_array_pt()
            ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level()++;

          bool use_coordinate_as_initial_guess = true;
          el_pt->locate_zeta(
            zeta, sub_geom_object_pt, s, use_coordinate_as_initial_guess);

          // Always fail? (Used for debugging, e.g. to trace out
          // spiral path)
          if (BinArray::Always_fail_elemental_locate_zeta)
          {
            sub_geom_object_pt = 0;
          }

#ifdef OOMPH_HAS_MPI
        }
#endif
        // Next one please
        i++;
      }
    }
  }

  //==============================================================================
  /// Boundaries of bin in each coordinate direction. *.first = min;
  /// *.second = max.
  //==============================================================================
  void RefineableBin::get_bin_boundaries(
    Vector<std::pair<double, double>>& min_and_max_coordinates)
  {
    unsigned n_bin = Bin_index_in_bin_array;

    // temporary storage for the eulerian dim
    unsigned current_dim = Bin_array_pt->ndim_zeta();
    min_and_max_coordinates.resize(current_dim);
    for (unsigned u = 0; u < current_dim; u++)
    {
      // The number of bins there are according to the u-th dimension
      double nbin_in_dir = n_bin % Bin_array_pt->dimension_of_bin_array(u);
      n_bin /= Bin_array_pt->dimension_of_bin_array(u);

      // The range between the maximum and minimum u-th coordinates of a bin
      double range = (Bin_array_pt->min_and_max_coordinates(u).second -
                      Bin_array_pt->min_and_max_coordinates(u).first) /
                     double(Bin_array_pt->dimension_of_bin_array(u));

      // Now updating the minimum and maximum u-th coordinates for this bin.
      min_and_max_coordinates[u].first =
        Bin_array_pt->min_and_max_coordinates(u).first + nbin_in_dir * range;
      min_and_max_coordinates[u].second =
        min_and_max_coordinates[u].first + range;
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///                        SamplePointContainer base class
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /// File to record sequence of visited sample points in
  std::ofstream SamplePointContainer::Visited_sample_points_file;

  /// Boolean flag to make to make locate zeta fail
  bool SamplePointContainer::Always_fail_elemental_locate_zeta = false;

  ///  Use equally spaced sample points? (otherwise vertices are sampled
  /// repeatedly
  bool SamplePointContainer::Use_equally_spaced_interior_sample_points = true;

  /// Time setup?
  bool SamplePointContainer::Enable_timing_of_setup = false;

  /// Offset of sample point container boundaries beyond max/min coords
  double SamplePointContainer::Percentage_offset = 5.0;


  //==============================================================================
  /// Max. bin dimension (number of bins in coordinate directions)
  //==============================================================================
  unsigned BinArray::max_bin_dimension() const
  {
    unsigned dim = ndim_zeta();
    unsigned n_max_level = Dimensions_of_bin_array[0];
    if (dim >= 2)
    {
      if (Dimensions_of_bin_array[1] > n_max_level)
      {
        n_max_level = Dimensions_of_bin_array[1];
      }
    }
    if (dim == 3)
    {
      if (Dimensions_of_bin_array[2] > n_max_level)
      {
        n_max_level = Dimensions_of_bin_array[2];
      }
    }
    return n_max_level;
  }


  //========================================================================
  /// Setup the min and max coordinates for the mesh, in each dimension
  //========================================================================
  void SamplePointContainer::setup_min_and_max_coordinates()
  {
    // Get the lagrangian dimension
    int n_lagrangian = ndim_zeta();

    // Storage locally (i.e. in parallel on each processor)
    // for the minimum and maximum coordinates
    double zeta_min_local[n_lagrangian];
    double zeta_max_local[n_lagrangian];
    for (int i = 0; i < n_lagrangian; i++)
    {
      zeta_min_local[i] = DBL_MAX;
      zeta_max_local[i] = -DBL_MAX;
    }

    // Loop over the elements of the mesh
    unsigned n_el = Mesh_pt->nelement();
    for (unsigned e = 0; e < n_el; e++)
    {
      FiniteElement* el_pt = Mesh_pt->finite_element_pt(e);

      // Get the number of vertices (nplot=2 does the trick)
      unsigned n_plot = 2;
      unsigned n_plot_points = el_pt->nplot_points(n_plot);

      // Loop over the number of plot points
      for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
      {
        Vector<double> s_local(n_lagrangian);
        Vector<double> zeta_global(n_lagrangian);

        // Get the local s -- need to sample over the entire range
        // of the elements!
        bool use_equally_spaced_interior_sample_points = false;
        el_pt->get_s_plot(
          iplot, n_plot, s_local, use_equally_spaced_interior_sample_points);

        // Now interpolate to global coordinates
        if (Use_eulerian_coordinates_during_setup)
        {
          el_pt->interpolated_x(s_local, zeta_global);
        }
        else
        {
          el_pt->interpolated_zeta(s_local, zeta_global);
        }

        // Check the max and min in each direction
        for (int i = 0; i < n_lagrangian; i++)
        {
          // Is the coordinate less than the minimum?
          if (zeta_global[i] < zeta_min_local[i])
          {
            zeta_min_local[i] = zeta_global[i];
          }
          // Is the coordinate bigger than the maximum?
          if (zeta_global[i] > zeta_max_local[i])
          {
            zeta_max_local[i] = zeta_global[i];
          }
        }
      }
    }

    // Global extrema - in parallel, need to get max/min across all processors
    double zeta_min[n_lagrangian];
    double zeta_max[n_lagrangian];
    for (int i = 0; i < n_lagrangian; i++)
    {
      zeta_min[i] = 0.0;
      zeta_max[i] = 0.0;
    }

#ifdef OOMPH_HAS_MPI
    // If the mesh has been distributed and we want consistent bins
    // across all processors
    if (Mesh_pt->is_mesh_distributed())
    {
      // .. we need a non-null communicator!
      if (Mesh_pt->communicator_pt() != 0)
      {
        int n_proc = Mesh_pt->communicator_pt()->nproc();
        if (n_proc > 1)
        {
          // Get the minima and maxima over all processors
          MPI_Allreduce(zeta_min_local,
                        zeta_min,
                        n_lagrangian,
                        MPI_DOUBLE,
                        MPI_MIN,
                        Mesh_pt->communicator_pt()->mpi_comm());
          MPI_Allreduce(zeta_max_local,
                        zeta_max,
                        n_lagrangian,
                        MPI_DOUBLE,
                        MPI_MAX,
                        Mesh_pt->communicator_pt()->mpi_comm());
        }
      }
      else // Null communicator - throw an error
      {
        std::ostringstream error_message_stream;
        error_message_stream << "Communicator not set for a Mesh\n"
                             << "that was created from a distributed Mesh";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else // If the mesh hasn't been distributed then the
    // max and min are the same on all processors
    {
      for (int i = 0; i < n_lagrangian; i++)
      {
        zeta_min[i] = zeta_min_local[i];
        zeta_max[i] = zeta_max_local[i];
      }
    }
#else // If we're not using MPI then the mesh can't be distributed
    for (int i = 0; i < n_lagrangian; i++)
    {
      zeta_min[i] = zeta_min_local[i];
      zeta_max[i] = zeta_max_local[i];
    }
#endif

    // Decrease/increase min and max to allow for any overshoot in
    // meshes that may move around
    // There's no point in doing this for DIM_LAGRANGIAN==1
    for (int i = 0; i < n_lagrangian; i++)
    {
      double length = zeta_max[i] - zeta_min[i];
      zeta_min[i] -= ((Percentage_offset / 100.0) * length);
      zeta_max[i] += ((Percentage_offset / 100.0) * length);
    }

    // Set the entries as the Min/Max_coords vector
    Min_and_max_coordinates.resize(n_lagrangian);
    for (int i = 0; i < n_lagrangian; i++)
    {
      Min_and_max_coordinates[i].first = zeta_min[i];
      Min_and_max_coordinates[i].second = zeta_max[i];
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///                                RefineableBin array class
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //==============================================================================
  /// Constructor
  //==============================================================================
  RefineableBinArray::RefineableBinArray(
    SamplePointContainerParameters* sample_point_container_parameters_pt)
    : SamplePointContainer(
        sample_point_container_parameters_pt->mesh_pt(),
        sample_point_container_parameters_pt->min_and_max_coordinates(),
        sample_point_container_parameters_pt
          ->use_eulerian_coordinates_during_setup(),
        sample_point_container_parameters_pt
          ->ignore_halo_elements_during_locate_zeta_search(),
        sample_point_container_parameters_pt
          ->nsample_points_generated_per_element()),
      BinArray(
        sample_point_container_parameters_pt->mesh_pt(),
        sample_point_container_parameters_pt->min_and_max_coordinates(),
        dynamic_cast<BinArrayParameters*>(sample_point_container_parameters_pt)
          ->dimensions_of_bin_array(),
        sample_point_container_parameters_pt
          ->use_eulerian_coordinates_during_setup(),
        sample_point_container_parameters_pt
          ->ignore_halo_elements_during_locate_zeta_search(),
        sample_point_container_parameters_pt
          ->nsample_points_generated_per_element())
  {
    RefineableBinArrayParameters* ref_bin_array_parameters_pt =
      dynamic_cast<RefineableBinArrayParameters*>(
        sample_point_container_parameters_pt);

#ifdef PARANOID
    if (ref_bin_array_parameters_pt == 0)
    {
      throw OomphLibError("Wrong sample_point_container_parameters_pt",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    Bin_array_is_recursive =
      ref_bin_array_parameters_pt->bin_array_is_recursive();
    Depth = ref_bin_array_parameters_pt->depth();
    Max_depth = ref_bin_array_parameters_pt->max_depth();
    Max_number_of_sample_point_per_bin =
      ref_bin_array_parameters_pt->max_number_of_sample_point_per_bin();
    Root_bin_array_pt = ref_bin_array_parameters_pt->root_bin_array_pt();

    // Set default size of bin array (and spatial dimension!)
    if (Dimensions_of_bin_array.size() == 0)
    {
      int dim = 0;
      if (Mesh_pt->nelement() != 0)
      {
        dim = Mesh_pt->finite_element_pt(0)->dim();
      }

      // Need to do an Allreduce to ensure that the dimension is consistent
      // even when no elements are assigned to a certain processor
#ifdef OOMPH_HAS_MPI
      // Only a problem if the mesh has been distributed
      if (Mesh_pt->is_mesh_distributed())
      {
        // Need a non-null communicator
        if (Mesh_pt->communicator_pt() != 0)
        {
          int n_proc = Mesh_pt->communicator_pt()->nproc();
          if (n_proc > 1)
          {
            int dim_reduce;
            MPI_Allreduce(&dim,
                          &dim_reduce,
                          1,
                          MPI_INT,
                          MPI_MAX,
                          Mesh_pt->communicator_pt()->mpi_comm());
            dim = dim_reduce;
          }
        }
      }
#endif

      Dimensions_of_bin_array.resize(dim, Default_n_bin_1d);
    }

    // Have we specified max/min coordinates?
    // If not, compute them on the fly from mesh
    if (Min_and_max_coordinates.size() == 0)
    {
      setup_min_and_max_coordinates();
    }

    // Get total number of bins and make space
    unsigned dim = Dimensions_of_bin_array.size();
    unsigned n_bin = 1;
    for (unsigned i = 0; i < dim; i++)
    {
      n_bin *= Dimensions_of_bin_array[i];
    }
    Bin_pt.resize(n_bin, 0);

    // I'm my own root bin array
    if (Depth == 0)
    {
      Root_bin_array_pt = this;
    }

#ifdef PARANOID
    if (Depth > 0)
    {
      if (Root_bin_array_pt == 0)
      {
        throw OomphLibError(
          "Must specify root_bin_array for lower-level bin arrays\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Initialise
    Total_number_of_sample_points_visited_during_locate_zeta_from_top_level = 0;
    First_sample_point_to_actually_lookup_during_locate_zeta = 0;
    Last_sample_point_to_actually_lookup_during_locate_zeta = UINT_MAX;
    Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta =
      2; // hierher tune this and create public static default
    Initial_last_sample_point_to_actually_lookup_during_locate_zeta =
      10; // hierher UINT MAX is temporary bypass! tune this  and create public
          // static default

    // Now fill the bastard...
    if (Depth == 0)
    {
      // Time it
      double t_start = 0.0;
      if (SamplePointContainer::Enable_timing_of_setup)
      {
        t_start = TimingHelpers::timer();
      }
      fill_bin_array();
      if (SamplePointContainer::Enable_timing_of_setup)
      {
        double t_end = TimingHelpers::timer();
        unsigned npts = total_number_of_sample_points_computed_recursively();
        oomph_info << "Time for setup of " << dim
                   << "-dimensional sample point container containing " << npts
                   << " sample points: " << t_end - t_start
                   << " sec  (ref_bin); third party: 0 sec ( = 0 %)"
                   << std::endl;
      }
    }
  }


  //==============================================================================
  /// Boundaries of specified bin in each coordinate direction.
  /// *.first = min; *.second = max.
  //==============================================================================
  void RefineableBinArray::get_bin_boundaries(
    const unsigned& bin_index,
    Vector<std::pair<double, double>>& min_and_max_coordinates_of_bin)
  {
    unsigned bin_index_local = bin_index;

    // temporary storage for the eulerian dim
    unsigned current_dim = ndim_zeta();
    min_and_max_coordinates_of_bin.resize(current_dim);
    for (unsigned u = 0; u < current_dim; u++)
    {
      // The number of bins there are according to the u-th dimension
      unsigned nbin_in_dir = bin_index_local % dimension_of_bin_array(u);
      bin_index_local /= dimension_of_bin_array(u);

      // The range between the maximum and minimum u-th coordinates of a bin
      double range =
        (Min_and_max_coordinates[u].second - Min_and_max_coordinates[u].first) /
        double(Dimensions_of_bin_array[u]);

      // Now updating the minimum and maximum u-th coordinates for this bin.
      min_and_max_coordinates_of_bin[u].first =
        Min_and_max_coordinates[u].first + double(nbin_in_dir) * range;
      min_and_max_coordinates_of_bin[u].second =
        min_and_max_coordinates_of_bin[u].first + range;
    }
  }


  //==============================================================================
  /// Output neighbouring bins up to given "radius" of the specified bin
  //==============================================================================
  void RefineableBinArray::output_bin_vertices(std::ofstream& outfile)
  {
    // Loop over bins
    unsigned n_bin = Bin_pt.size();
    for (unsigned i = 0; i < n_bin; i++)
    {
      if (Bin_pt[i] != 0)
      {
        Bin_pt[i]->output_bin_vertices(outfile);
      }
    }
  }


  //==============================================================================
  /// Output neighbouring bins up to given "radius" of the specified bin
  //==============================================================================
  void RefineableBinArray::output_neighbouring_bins(const unsigned& bin_index,
                                                    const unsigned& radius,
                                                    std::ofstream& outfile)
  {
    unsigned n_lagr = ndim_zeta();

    Vector<unsigned> neighbouring_bin_index;
    get_neighbouring_bins_helper(bin_index, radius, neighbouring_bin_index);
    unsigned nneigh = neighbouring_bin_index.size();

    // Outline of bin structure
    switch (n_lagr)
    {
      case 1:
        outfile << "ZONE I=2\n"
                << Min_and_max_coordinates[0].first << std::endl
                << Min_and_max_coordinates[0].second << std::endl;
        break;

      case 2:

        outfile << "ZONE I=2, J=2\n"
                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].first << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].first << " "
                << "\n"

                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].second << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].second << " "
                << "\n";
        break;

      case 3:

        outfile << "ZONE I=2, J=2, K=2 \n"
                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].first << " "
                << Min_and_max_coordinates[2].first << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].first << " "
                << Min_and_max_coordinates[2].first << " "
                << "\n"

                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].second << " "
                << Min_and_max_coordinates[2].first << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].second << " "
                << Min_and_max_coordinates[2].first << " "
                << "\n"

                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].first << " "
                << Min_and_max_coordinates[2].second << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].first << " "
                << Min_and_max_coordinates[2].second << " "
                << "\n"

                << Min_and_max_coordinates[0].first << " "
                << Min_and_max_coordinates[1].second << " "
                << Min_and_max_coordinates[2].second << " "
                << "\n"

                << Min_and_max_coordinates[0].second << " "
                << Min_and_max_coordinates[1].second << " "
                << Min_and_max_coordinates[2].second << " "
                << "\n";
        break;

      default:

        oomph_info << "n_lagr=" << n_lagr << std::endl;
        throw OomphLibError(
          "Wrong dimension!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Loop over neighbours
    for (unsigned i = 0; i < nneigh; i++)
    {
      Vector<std::pair<double, double>> min_and_max_coordinates(n_lagr);
      get_bin_boundaries(neighbouring_bin_index[i], min_and_max_coordinates);

      switch (n_lagr)
      {
        case 1:
          outfile << "ZONE I=2\n"
                  << min_and_max_coordinates[0].first << std::endl
                  << min_and_max_coordinates[0].second << std::endl;
          break;

        case 2:

          outfile << "ZONE I=2, J=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << "\n";
          break;

        case 3:

          outfile << "ZONE I=2, J=2, K=2\n"
                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].first << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].first << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].first << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n"

                  << min_and_max_coordinates[0].second << " "
                  << min_and_max_coordinates[1].second << " "
                  << min_and_max_coordinates[2].second << " "
                  << "\n";
          break;

        default:

          oomph_info << "n_lagr=" << n_lagr << std::endl;
          throw OomphLibError("Wrong dimension!",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  //========================================================================
  /// Profiling function to compare performance of two different
  /// versions of the get_neighbouring_bins_helper(...) function
  //========================================================================
  void BinArray::profile_get_neighbouring_bins_helper()
  {
    unsigned old_faster = 0;
    unsigned new_faster = 0;
    double t_total_new = 0.0;
    double t_total_old = 0.0;

    unsigned dim = ndim_zeta();

    // Choose bin in middle of the domain
    Vector<double> zeta(dim);
    for (unsigned i = 0; i < dim; i++)
    {
      zeta[i] = 0.5 * (Min_and_max_coordinates[i].first +
                       Min_and_max_coordinates[i].second);
    }

    // Finding the bin in which the point is located
    int bin_index = coords_to_bin_index(zeta);

#ifdef PARANOID
    if (bin_index < 0)
    {
      throw OomphLibError("Negative bin index...",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Start at this radius (radius = 0 is the central bin)
    unsigned radius = 0;

    // "coordinates" of the bin which is most likely to contain the
    // point
    Vector<unsigned> bin_index_v(dim);
    coords_to_vectorial_bin_index(zeta, bin_index_v);

    // We loop over all the dimensions to find the maximum radius we have to
    // do the spiraling if we want to be sure there is no bin left at the end
    unsigned max_radius = 0;
    for (unsigned k = 0; k < dim; k++)
    {
      unsigned local =
        std::max((bin_index_v[k] + 1),
                 (Dimensions_of_bin_array[k] - bin_index_v[k] - 1));
      if (local > max_radius)
      {
        max_radius = local;
      }
    }

    // Vector which will store the indices of the neighbouring bins
    // at the current radius
    Vector<unsigned> bin_index_at_current_radius_old;
    Vector<unsigned> bin_index_at_current_radius_new;
    while (radius <= max_radius)
    {
      // Get the neighbouring bins
      bin_index_at_current_radius_old.clear();
      double t_start = TimingHelpers::timer();
      get_neighbouring_bins_helper(
        bin_index, radius, bin_index_at_current_radius_old, true);
      unsigned nbin_at_current_radius_old =
        bin_index_at_current_radius_old.size();
      double t_end = TimingHelpers::timer();
      double t_old = t_end - t_start;

      bin_index_at_current_radius_new.clear();
      double t_start_new = TimingHelpers::timer();
      get_neighbouring_bins_helper(
        bin_index, radius, bin_index_at_current_radius_new, false);
      unsigned nbin_at_current_radius_new =
        bin_index_at_current_radius_new.size();
      double t_end_new = TimingHelpers::timer();

      double t_new = t_end_new - t_start_new;

      if (nbin_at_current_radius_new != nbin_at_current_radius_old)
      {
        oomph_info << "Number of bins don't match: new = "
                   << nbin_at_current_radius_new
                   << "old = " << nbin_at_current_radius_old
                   << " radius = " << radius << std::endl;
        oomph_info << "Old: " << std::endl;
        for (unsigned i = 0; i < nbin_at_current_radius_old; i++)
        {
          oomph_info << bin_index_at_current_radius_old[i] << " ";
        }
        oomph_info << std::endl;
        oomph_info << "New: " << std::endl;
        for (unsigned i = 0; i < nbin_at_current_radius_new; i++)
        {
          oomph_info << bin_index_at_current_radius_new[i] << " ";
        }
        oomph_info << std::endl;
      }

      t_total_new += t_new;
      t_total_old += t_old;

      if (t_new < t_old)
      {
        new_faster++;
      }
      else
      {
        old_faster++;
      }
      radius++;
    }


    oomph_info << "Number of times old/new version was faster: " << old_faster
               << " " << new_faster << std::endl
               << "Total old/new time: " << t_total_old << " " << t_total_new
               << " " << std::endl;
  }


  //==============================================================================
  /// Helper function for computing the bin indices of neighbouring bins
  /// at a given "radius" of the specified bin
  //==============================================================================
  void BinArray::get_neighbouring_bins_helper(
    const unsigned& bin_index,
    const unsigned& radius,
    Vector<unsigned>& neighbouring_bin_index,
    const bool& use_old_version)
  {
    // OLD VERSION
    if (use_old_version)
    {
      neighbouring_bin_index.clear();

      // Quick return (slightly hacky -- the machinery below adds the
      // "home" bin twice...)
      if (radius == 0)
      {
        neighbouring_bin_index.push_back(bin_index);
        return;
      }

      // translate old code
      unsigned level = radius;
      unsigned bin = bin_index;

      // Dimension
      const unsigned n_lagrangian = this->ndim_zeta();

      // This will be different depending on the number of Lagrangian
      // coordinates
      if (n_lagrangian == 1)
      {
        // Reserve memory for the container where we return the indices
        // of the neighbouring bins (2 bins max, left and right)
        neighbouring_bin_index.reserve(2);

        // Single "loop" in one direction - always a vector of max size 2
        unsigned nbr_bin_left = bin - level;
        if (nbr_bin_left < Dimensions_of_bin_array[0])
        {
          unsigned nbr_bin = nbr_bin_left;
          neighbouring_bin_index.push_back(nbr_bin);
        }
        unsigned nbr_bin_right = bin + level;
        if ((nbr_bin_right < Dimensions_of_bin_array[0]) &&
            (nbr_bin_right != nbr_bin_left))
        {
          unsigned nbr_bin = nbr_bin_right;
          neighbouring_bin_index.push_back(nbr_bin);
        }
      }
      else if (n_lagrangian == 2)
      {
        // Reserve memory for the container where we return the indices
        // of the neighbouring bins
        const unsigned n_max_neighbour_bins = 8 * level;
        neighbouring_bin_index.reserve(n_max_neighbour_bins);

        const unsigned n_total_bin =
          Dimensions_of_bin_array[0] * Dimensions_of_bin_array[1];

        // Which row of the bin structure is the current bin on?
        // This is just given by the integer answer of dividing bin
        // by Nbin_x (the number of bins in a single row)
        // e.g. in a 6x6 grid, bins 6 through 11 would all have bin_row=1
        const unsigned bin_row = bin / Dimensions_of_bin_array[0];

        // The neighbour_bin vector contains all bin numbers at the
        // specified "distance" (level) away from the current bin

        // Row/column length
        const unsigned n_length = (level * 2) + 1;

        {
          // Visit all the bins at the specified distance (level) away
          // from the current bin. In order to obtain the same order in
          // the visited bins as the previous algorithm we visit all the
          // bins at the specified distance (level) as follows:

          // Suppose we want the bins at distance (level=2) of the
          // specified bin, then we visit them as indicated below

          // 01 02 03 04 05   // First row
          // 06          07
          // 08    B     09
          // 10          11
          // 12 13 14 15 16   // Last row
          // ^--------------- First column
          //              ^-- Last column

          // ----------------------------------------------------------------
          // Visit all the bins in the first row at the specified
          // distance (level) away from the current bin

          // ------------------ FIRST ROW ------------------------
          // Pre-compute the distance in the j-direction
          const unsigned j_precomputed = level * Dimensions_of_bin_array[0];
          // Pre-compute the row where the bin should lie on
          const unsigned j_initial_row = bin_row - level;

          // Loop over the columns (of the first row)
          for (unsigned i = 0; i < n_length; i++)
          {
            // --------- First row ------------------------------------------
            const unsigned initial_neighbour_bin =
              bin - (level - i) - j_precomputed;
            // This number might fall on the wrong row of the bin
            // structure; this needs to be tested? Not sure why, but leave
            // the test!

            // Which row is this number on? (see above)
            const unsigned initial_neighbour_bin_row =
              initial_neighbour_bin / Dimensions_of_bin_array[0];
            // These numbers for the rows must match; if it is then add
            // initial_neighbour_bin to the neighbour scheme (The bin
            // number must also be greater than zero and less than the
            // total number of bins)
            if ((j_initial_row == initial_neighbour_bin_row) &&
                (initial_neighbour_bin < n_total_bin))
            {
              neighbouring_bin_index.push_back(initial_neighbour_bin);
            }

          } // for (unsigned i=0;i<n_length;i++)

          // Then visit all the bins in the first and last column at the
          // specified distance (level) away from the current bin

          // ------------------ FIRST AND LAST COLUMNS ---------------------
          // Loop over the rows (of the first and last column)
          for (unsigned j = 1; j < n_length - 1; j++)
          {
            // --------- First column ---------------------------------------
            const unsigned initial_neighbour_bin =
              bin - (level) - ((level - j) * Dimensions_of_bin_array[0]);
            // This number might fall on the wrong row of the bin
            // structure; this needs to be tested? Not sure why, but leave
            // the test!

            // Which row is this number on? (see above)
            const unsigned initial_neighbour_bin_row =
              initial_neighbour_bin / Dimensions_of_bin_array[0];

            // Which row should it be on?
            const unsigned initial_row = bin_row - (level - j);

            // These numbers for the rows must match; if it is then add
            // initial_neighbour_bin to the neighbour scheme (The bin
            // number must also be greater than zero and less than the
            // total number of bins)
            if ((initial_row == initial_neighbour_bin_row) &&
                (initial_neighbour_bin < n_total_bin))
            {
              neighbouring_bin_index.push_back(initial_neighbour_bin);
            }

            // --------- Last column -----------------------------------------
            const unsigned final_neighbour_bin =
              bin + (level) - ((level - j) * Dimensions_of_bin_array[0]);
            // This number might fall on the wrong row of the bin
            // structure; this needs to be tested? Not sure why, but leave
            // the test!

            // Which row is this number on? (see above)
            const unsigned final_neighbour_bin_row =
              final_neighbour_bin / Dimensions_of_bin_array[0];

            // Which row should it be on?
            const unsigned final_row = bin_row - (level - j);

            // These numbers for the rows must match; if it is then add
            // initial_neighbour_bin to the neighbour scheme (The bin
            // number must also be greater than zero and less than the
            // total number of bins)
            if ((final_row == final_neighbour_bin_row) &&
                (final_neighbour_bin < n_total_bin))
            {
              neighbouring_bin_index.push_back(final_neighbour_bin);
            }

          } // for (unsigned j=1;j<n_length-1;j++)

          // ------------------ LAST ROW ------------------------
          // Pre-compute the row where the bin should lie on
          const unsigned j_final_row = bin_row + level;

          // Loop over the columns (of the last row)
          for (unsigned i = 0; i < n_length; i++)
          {
            // --------- Last row ------------------------------------------
            const unsigned final_neighbour_bin =
              bin - (level - i) + j_precomputed;
            // This number might fall on the wrong row of the bin
            // structure; this needs to be tested? Not sure why, but leave
            // the test!

            // Which row is this number on? (see above)
            const unsigned final_neighbour_bin_row =
              final_neighbour_bin / Dimensions_of_bin_array[0];
            // These numbers for the rows must match; if it is then add
            // initial_neighbour_bin to the neighbour scheme (The bin
            // number must also be greater than zero and less than the
            // total number of bins)
            if ((j_final_row == final_neighbour_bin_row) &&
                (final_neighbour_bin < n_total_bin))
            {
              neighbouring_bin_index.push_back(final_neighbour_bin);
            }

          } // for (unsigned i=0;i<n_length;i++)
        }
      }
      else if (n_lagrangian == 3)
      {
        // Reserve memory for the container where we return the indices
        // of the neighbouring bins
        const unsigned n_max_neighbour_bins =
          8 * level * (3 + 2 * (level - 1)) +
          2 * (2 * (level - 1) + 1) * (2 * (level - 1) + 1);
        neighbouring_bin_index.reserve(n_max_neighbour_bins);

        unsigned n_total_bin = Dimensions_of_bin_array[0] *
                               Dimensions_of_bin_array[1] *
                               Dimensions_of_bin_array[2];

        // Which layer of the bin structure is the current bin on?
        // This is just given by the integer answer of dividing bin
        // by Nbin_x*Nbin_y (the number of bins in a single layer
        // e.g. in a 6x6x6 grid, bins 72 through 107 would all have bin_layer=2
        unsigned bin_layer =
          bin / (Dimensions_of_bin_array[0] * Dimensions_of_bin_array[1]);

        // Which row in this layer is the bin number on?
        unsigned bin_row = (bin / Dimensions_of_bin_array[0]) -
                           (bin_layer * Dimensions_of_bin_array[1]);

        // The neighbour_bin vector contains all bin numbers at the
        // specified "distance" (level) away from the current bin

        // Row/column/layer length
        unsigned n_length = (level * 2) + 1;

        // Loop over the layers
        for (unsigned k = 0; k < n_length; k++)
        {
          // Loop over the rows
          for (unsigned j = 0; j < n_length; j++)
          {
            // Loop over the columns
            for (unsigned i = 0; i < n_length; i++)
            {
              // Only do this for the end points of every row/layer/column
              if ((k == 0) || (k == n_length - 1) || (j == 0) ||
                  (j == n_length - 1) || (i == 0) || (i == n_length - 1))
              {
                unsigned nbr_bin = bin - level + i -
                                   ((level - j) * Dimensions_of_bin_array[0]) -
                                   ((level - k) * Dimensions_of_bin_array[0] *
                                    Dimensions_of_bin_array[1]);
                // This number might fall on the wrong
                // row or layer of the bin structure; this needs to be tested

                // Which layer is this number on?
                unsigned nbr_bin_layer = nbr_bin / (Dimensions_of_bin_array[0] *
                                                    Dimensions_of_bin_array[1]);

                // Which row is this number on? (see above)
                unsigned nbr_bin_row =
                  (nbr_bin / Dimensions_of_bin_array[0]) -
                  (nbr_bin_layer * Dimensions_of_bin_array[1]);

                // Which layer and row should it be on, given level?
                unsigned layer = bin_layer - level + k;
                unsigned row = bin_row - level + j;

                // These layers and rows must match up:
                // if so then add nbr_bin to the neighbour schemes
                // (The bin number must also be greater than zero
                //  and less than the total number of bins)
                if ((row == nbr_bin_row) && (layer == nbr_bin_layer) &&
                    (nbr_bin < n_total_bin))
                {
                  neighbouring_bin_index.push_back(nbr_bin);
                }
              }
            }
          }
        }
      }
    }
    // LOUIS'S VERSION
    else
    {
      neighbouring_bin_index.clear();

      // Just testing if the radius is equal to 0
      if (radius == 0)
      {
        // in this case we only have to push back the bin
        neighbouring_bin_index.push_back(bin_index);
      }
      // Else, if the radius is different from 0
      else
      {
        unsigned dim = ndim_zeta();

        // The vector which will store the coordinates of the bin in the bin
        // arrray
        Vector<int> vector_of_positions(3);

        // Will store locally the dimension of the bin array
        Vector<int> vector_of_dimensions(3);

        // Will store the radiuses for each dimension
        // If it is 0, that means the dimension is not vector_of_active_dim
        Vector<int> vector_of_radiuses(3);

        // Stores true if the dimension is "active" or false if the
        // dimension is inactive (for example the third dimension is inactive in
        // a 2D mesh).
        std::vector<bool> vector_of_active_dim(3);

        // Will store the coefficients you have to multiply each bin coordinate
        // to have the correct unsigned corresponding to the bin index.
        Vector<int> vector_of_coef(3);


        // First initializing this MyArrays for the active dimensions
        unsigned coef = 1;
        for (unsigned u = 0; u < dim; u++)
        {
          vector_of_positions[u] =
            (((bin_index / coef) % (Dimensions_of_bin_array[u])));
          vector_of_coef[u] = coef;
          coef *= Dimensions_of_bin_array[u];
          vector_of_dimensions[u] = Dimensions_of_bin_array[u];
          vector_of_radiuses[u] = radius;
          vector_of_active_dim[u] = true;
        }
        // Filling the rest with default values
        for (unsigned u = dim; u < 3; u++)
        {
          vector_of_positions[u] = 0;
          vector_of_coef[u] = 0;
          vector_of_dimensions[u] = 1;
          vector_of_radiuses[u] = 0;
          vector_of_active_dim[u] = false;
        }

        // First we fill the bins which corresponds to x+radius and x-radius in
        // the bin array (x being the first coordinate/dimension).
        // For j equal to radius or -radius
        for (int j = -vector_of_radiuses[0]; j <= vector_of_radiuses[0];
             j += 2 * vector_of_radiuses[0]) // j corresponds to x
        {
          int local_tempj =
            vector_of_positions[0] + j; // Corresponding bin index
          // if we are in the bin array
          if (local_tempj >= 0 && local_tempj < vector_of_dimensions[0])
          {
            // Loop over all the bins
            for (int i = -vector_of_radiuses[1]; i <= vector_of_radiuses[1];
                 i++)
            {
              // i corresponds to y (2nd dim)
              int local_tempi = vector_of_positions[1] + i;
              // if we are still in the bin array
              if (local_tempi >= 0 && local_tempi < vector_of_dimensions[1])
              {
                for (int k = -vector_of_radiuses[2]; k <= vector_of_radiuses[2];
                     k++)
                {
                  // k corresponds to z (third dimension)
                  int local_tempk = vector_of_positions[2] + k;
                  // if we are still in the bin array
                  if (local_tempk >= 0 && local_tempk < vector_of_dimensions[2])
                  {
                    neighbouring_bin_index.push_back(
                      local_tempj * vector_of_coef[0] +
                      local_tempi * vector_of_coef[1] +
                      local_tempk * vector_of_coef[2]);
                  }
                }
              }
            }
          }
        }
        // Secondly we get the bins corresponding to y+radius and y-radius
        // only if the second dimension is active
        if (vector_of_active_dim[1])
        {
          // For i equal to radius or -radius (i corresponds to y)
          for (int i = -vector_of_radiuses[1]; i <= vector_of_radiuses[1];
               i += 2 * vector_of_radiuses[1])
          {
            int local_tempi = vector_of_positions[1] + i;
            if (local_tempi >= 0 && local_tempi < vector_of_dimensions[1])
            {
              // We loop over all the surface
              for (int j = -vector_of_radiuses[0] + 1;
                   j <= vector_of_radiuses[0] - 1;
                   j++)
              {
                int local_tempj = vector_of_positions[0] + j;
                if (local_tempj >= 0 && local_tempj < vector_of_dimensions[0])
                {
                  for (int k = -vector_of_radiuses[2];
                       k <= vector_of_radiuses[2];
                       k++)
                  {
                    int local_tempk = vector_of_positions[2] + k;
                    if (local_tempk >= 0 &&
                        local_tempk < vector_of_dimensions[2])
                    {
                      neighbouring_bin_index.push_back(
                        local_tempj * vector_of_coef[0] +
                        local_tempi * vector_of_coef[1] +
                        local_tempk * vector_of_coef[2]);
                    }
                  }
                }
              }
            }
          }
        }
        // Thirdly we get the bins corresponding to z+radius and z-radius
        // if the third dimension is active.
        if (vector_of_active_dim[2])
        {
          // for k equal to radius or -radius (k corresponds to z)
          for (int k = -vector_of_radiuses[2]; k <= vector_of_radiuses[2];
               k += 2 * vector_of_radiuses[2])
          {
            int local_tempk = vector_of_positions[2] + k;
            if (local_tempk >= 0 && local_tempk < vector_of_dimensions[2])
            {
              // We loop over all the surface
              for (int j = -vector_of_radiuses[0] + 1;
                   j <= vector_of_radiuses[0] - 1;
                   j++)
              {
                int local_tempj = vector_of_positions[0] + j;
                if (local_tempj >= 0 && local_tempj < vector_of_dimensions[0])
                {
                  for (int i = -vector_of_radiuses[1] + 1;
                       i <= vector_of_radiuses[1] - 1;
                       i++)
                  {
                    int local_tempi = vector_of_positions[1] + i;
                    if (local_tempi >= 0 &&
                        local_tempi < vector_of_dimensions[1])
                    {
                      neighbouring_bin_index.push_back(
                        local_tempj * vector_of_coef[0] +
                        local_tempi * vector_of_coef[1] +
                        local_tempk * vector_of_coef[2]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } // end new version
  }


  //==============================================================================
  /// Compute total number of sample points recursively
  //==============================================================================
  unsigned RefineableBinArray::
    total_number_of_sample_points_computed_recursively() const
  {
    unsigned count = 0;
    unsigned n_bin = nbin();
    for (unsigned i = 0; i < n_bin; i++)
    {
      if (Bin_pt[i] != 0)
      {
        count +=
          Bin_pt[i]->total_number_of_sample_points_computed_recursively();
      }
    }
    return count;
  }


  //============================================================
  /// Get (linearly enumerated) bin index of bin that
  /// contains specified zeta
  //============================================================
  unsigned BinArray::coords_to_bin_index(const Vector<double>& zeta)
  {
    unsigned coef = 1;
    unsigned n_bin = 0;
    const unsigned dim = ndim_zeta();

    // Loop over all the dimensions
    for (unsigned u = 0; u < dim; u++)
    {
      // for each one get the correct bin "coordinate"
      unsigned local_bin_number = 0;

      if (zeta[u] < Min_and_max_coordinates[u].first)
      {
        local_bin_number = 0;
      }
      else if (zeta[u] > Min_and_max_coordinates[u].second)
      {
        local_bin_number = dimensions_of_bin_array(u) - 1;
      }
      else
      {
        local_bin_number = (int(std::min(
          dimensions_of_bin_array(u) - 1,
          unsigned(floor(double(zeta[u] - Min_and_max_coordinates[u].first) /
                         double(Min_and_max_coordinates[u].second -
                                Min_and_max_coordinates[u].first) *
                         double(dimensions_of_bin_array(u)))))));
      }

      /// Update the coef and the bin index
      n_bin += local_bin_number * coef;
      coef *= dimensions_of_bin_array(u);
    }

    // return the correct bin index
    return (n_bin);
  }

  //==========================================================================
  /// Fill the bin array with sample points from FiniteElements stored in mesh
  //==========================================================================
  void RefineableBinArray::fill_bin_array()
  {
    // Fill 'em in:
    unsigned nel = Mesh_pt->nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* el_pt = Mesh_pt->finite_element_pt(e);

      // Total number of sample point we will create
      unsigned nplot =
        el_pt->nplot_points(Nsample_points_generated_per_element);

      /// For all the sample points we have to create ...
      for (unsigned j = 0; j < nplot; j++)
      {
        // ... create it: Pass element index in mesh (vector
        // of elements and index of sample point within element
        SamplePoint* new_sample_point_pt = new SamplePoint(e, j);

        // Coordinates of this point
        Vector<double> zeta(ndim_zeta());
        Vector<double> s(ndim_zeta());
        bool use_equally_spaced_interior_sample_points =
          SamplePointContainer::Use_equally_spaced_interior_sample_points;
        el_pt->get_s_plot(j,
                          Nsample_points_generated_per_element,
                          s,
                          use_equally_spaced_interior_sample_points);
        if (Use_eulerian_coordinates_during_setup)
        {
          el_pt->interpolated_x(s, zeta);
        }
        else
        {
          el_pt->interpolated_zeta(s, zeta);
        }

#ifdef PARANOID

        // Check if point is inside
        bool is_inside = true;
        std::ostringstream error_message;
        unsigned dim = ndim_zeta();
        for (unsigned i = 0; i < dim; i++)
        {
          if ((zeta[i] < Min_and_max_coordinates[i].first) ||
              (zeta[i] > Min_and_max_coordinates[i].second))
          {
            is_inside = false;
            error_message << "Sample point at zeta[" << i << "]  = " << zeta[i]
                          << " is outside limits of bin array: "
                          << Min_and_max_coordinates[i].first << " and "
                          << Min_and_max_coordinates[i].second << std::endl;
          }
        }

        if (!is_inside)
        {
          error_message << "Please correct the limits passed to the "
                        << "constructor." << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

#endif


        // Finding the correct bin to put the sample point
        unsigned bin_index = coords_to_bin_index(zeta);

        // if the bin is not yet created, create it
        if (Bin_pt[bin_index] == 0)
        {
          Bin_pt[bin_index] = new RefineableBin(this, bin_index);
        }

        // ... and then fill the bin with this new sample point
        Bin_pt[bin_index]->add_sample_point(new_sample_point_pt, zeta);
      }
    }
  }


  //==================================================================
  /// Get "coordinates" of bin that contains specified zeta
  //==================================================================
  void BinArray::coords_to_vectorial_bin_index(const Vector<double>& zeta,
                                               Vector<unsigned>& bin_index)
  {
    unsigned dim = ndim_zeta();
    bin_index.resize(dim);
    for (unsigned u = 0; u < dim; u++)
    {
      if (zeta[u] < Min_and_max_coordinates[u].first)
      {
        bin_index[u] = 0;
      }
      else if (zeta[u] > Min_and_max_coordinates[u].second)
      {
        bin_index[u] = dimensions_of_bin_array(u) - 1;
      }
      else
      {
        bin_index[u] = (int(
          std::min(dimensions_of_bin_array(u) - 1,
                   unsigned(floor((zeta[u] - Min_and_max_coordinates[u].first) /
                                  double(Min_and_max_coordinates[u].second -
                                         Min_and_max_coordinates[u].first) *
                                  double(dimensions_of_bin_array(u)))))));
      }
    }
  }


  //==============================================================================
  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  //==============================================================================
  void RefineableBinArray::locate_zeta(const Vector<double>& zeta,
                                       GeomObject*& sub_geom_object_pt,
                                       Vector<double>& s)
  {
    // Default: we've failed miserably
    sub_geom_object_pt = 0;

    unsigned dim = ndim_zeta();

    // Top level book keeping and sanity checking
    if (Depth == 0)
    {
      // Reset counter for number of sample points visited.
      // If we can't find the point we should at least make sure that
      // we've visited all the sample points before giving up.
      Total_number_of_sample_points_visited_during_locate_zeta_from_top_level =
        0;

      // Does the zeta coordinate lie within the current (top level!) bin
      // structure? Skip this test if we want to always fail because that's
      //  usually done to trace out the spiral path
      if (!BinArray::Always_fail_elemental_locate_zeta)
      {
        // Loop over the lagrangian dimension
        for (unsigned i = 0; i < dim; i++)
        {
          // If the i-th coordinate is less than the minimum
          if (zeta[i] < Min_and_max_coordinates[i].first)
          {
            return;
          }
          // Otherwise coordinate may be bigger than the maximum
          else if (zeta[i] > Min_and_max_coordinates[i].second)
          {
            return;
          }
        }
      }
    }

    // Finding the bin in which the point is located
    int bin_index = coords_to_bin_index(zeta);

#ifdef PARANOID
    if (bin_index < 0)
    {
      throw OomphLibError("Negative bin index...",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Start at this radius (radius = 0 is the central bin)
    unsigned radius = 0;

    // "coordinates" of the bin which is most likely to contain the
    // point
    Vector<unsigned> bin_index_v(dim);
    coords_to_vectorial_bin_index(zeta, bin_index_v);

    // We loop over all the dimensions to find the maximum radius we have to
    // do the spiraling if we want to be sure there is no bin left at the end
    unsigned max_radius = 0;
    for (unsigned k = 0; k < dim; k++)
    {
      unsigned local =
        std::max((bin_index_v[k] + 1),
                 (Dimensions_of_bin_array[k] - bin_index_v[k] - 1));
      if (local > max_radius)
      {
        max_radius = local;
      }
    }

    // Vector which will store the indices of the neighbouring bins
    // at the current radius
    Vector<unsigned> bin_index_at_current_radius;
    while (radius <= max_radius)
    {
      // Get the neighbouring bins
      bin_index_at_current_radius.clear();
      get_neighbouring_bins_helper(
        bin_index, radius, bin_index_at_current_radius);
      // How many are there
      unsigned nbin_at_current_radius = bin_index_at_current_radius.size();
      unsigned n_bin = nbin();

      // Keep looping over entries (stop if we found geom object)
      unsigned k = 0;
      while ((k < nbin_at_current_radius) && (sub_geom_object_pt == 0))
      {
        int neigh_bin_index = bin_index_at_current_radius[k];
#ifdef PARANOID
        if (neigh_bin_index < 0)
        {
          throw OomphLibError("Negative neighbour bin index...",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        if (neigh_bin_index < int(n_bin))
        {
          // If the bin exists
          if (Bin_pt[neigh_bin_index] != 0)
          {
            // We call the correct method to locate_zeta in the bin
            Bin_pt[neigh_bin_index]->locate_zeta(zeta, sub_geom_object_pt, s);
          }
        }
        k++;
      }

      // Reached end of loop over bins at this radius (or found the point)
      // Which one?
      if (sub_geom_object_pt != 0)
      {
        return;
      }
      // Increment radius
      else
      {
        radius++;
      }
    }


#ifdef PARANOID

    // If we still haven't found the point, check that we've at least visited
    // all the sample points
    if (Depth == 0)
    {
      if (
        Total_number_of_sample_points_visited_during_locate_zeta_from_top_level !=
        total_number_of_sample_points_computed_recursively())
      {
        if (max_search_radius() == DBL_MAX)
        {
          std::ostringstream error_message;
          error_message
            << "Zeta not found after visiting "
            << Total_number_of_sample_points_visited_during_locate_zeta_from_top_level
            << " sample points out of "
            << total_number_of_sample_points_computed_recursively() << std::endl
            << "Where are the missing sample points???\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
    }

#endif
  }

  /// Default number of bins (in each coordinate direction)
  unsigned RefineableBinArray::Default_n_bin_1d = 5;


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Constructor
  //======================================================================
  NonRefineableBinArray::NonRefineableBinArray(
    SamplePointContainerParameters* sample_point_container_parameters_pt)
    : SamplePointContainer(
        sample_point_container_parameters_pt->mesh_pt(),
        sample_point_container_parameters_pt->min_and_max_coordinates(),
        sample_point_container_parameters_pt
          ->use_eulerian_coordinates_during_setup(),
        sample_point_container_parameters_pt
          ->ignore_halo_elements_during_locate_zeta_search(),
        sample_point_container_parameters_pt
          ->nsample_points_generated_per_element()),
      BinArray(
        sample_point_container_parameters_pt->mesh_pt(),
        sample_point_container_parameters_pt->min_and_max_coordinates(),
        dynamic_cast<BinArrayParameters*>(sample_point_container_parameters_pt)
          ->dimensions_of_bin_array(),
        sample_point_container_parameters_pt
          ->use_eulerian_coordinates_during_setup(),
        sample_point_container_parameters_pt
          ->ignore_halo_elements_during_locate_zeta_search(),
        sample_point_container_parameters_pt
          ->nsample_points_generated_per_element())
  {
    // Set default size of bin array (and spatial dimension!)
    if (Dimensions_of_bin_array.size() == 0)
    {
      int dim = 0;
      if (Mesh_pt->nelement() != 0)
      {
        dim = Mesh_pt->finite_element_pt(0)->dim();
      }


      // Need to do an Allreduce to ensure that the dimension is consistent
      // even when no elements are assigned to a certain processor
#ifdef OOMPH_HAS_MPI
      // Only a problem if the mesh has been distributed
      if (Mesh_pt->is_mesh_distributed())
      {
        // Need a non-null communicator
        if (Mesh_pt->communicator_pt() != 0)
        {
          int n_proc = Mesh_pt->communicator_pt()->nproc();
          if (n_proc > 1)
          {
            int dim_reduce;
            MPI_Allreduce(&dim,
                          &dim_reduce,
                          1,
                          MPI_INT,
                          MPI_MAX,
                          Mesh_pt->communicator_pt()->mpi_comm());
            dim = dim_reduce;
          }
        }
      }
#endif

      Dimensions_of_bin_array.resize(dim, Default_n_bin_1d);
    }

    // Have we specified max/min coordinates?
    // If not, compute them on the fly from mesh
    if (Min_and_max_coordinates.size() == 0)
    {
      setup_min_and_max_coordinates();
    }

    // Spiraling parameters
    Max_spiral_level = UINT_MAX;
    Current_min_spiral_level = 0;

    NonRefineableBinArrayParameters* non_ref_bin_array_parameters_pt =
      dynamic_cast<NonRefineableBinArrayParameters*>(
        sample_point_container_parameters_pt);

#ifdef PARANOID
    if (non_ref_bin_array_parameters_pt == 0)
    {
      throw OomphLibError("Wrong sample_point_container_parameters_pt",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    Nspiral_chunk = non_ref_bin_array_parameters_pt->nspiral_chunk();
    Current_max_spiral_level = max_bin_dimension();

    // Time it
    double t_start = 0.0;
    if (SamplePointContainer::Enable_timing_of_setup)
    {
      t_start = TimingHelpers::timer();
    }

    // Now fill the bastard...
    fill_bin_array();

    if (SamplePointContainer::Enable_timing_of_setup)
    {
      double t_end = TimingHelpers::timer();
      unsigned npts = total_number_of_sample_points_computed_recursively();
      oomph_info << "Time for setup of " << Dimensions_of_bin_array.size()
                 << "-dimensional sample point container containing " << npts
                 << " sample points: " << t_end - t_start
                 << " sec  (non-ref_bin); third party: 0 sec ( = 0 %)"
                 << std::endl;
    }
  }


  //==============================================================================
  /// Compute total number of sample points recursively
  //==============================================================================
  unsigned NonRefineableBinArray::
    total_number_of_sample_points_computed_recursively() const
  {
    // Get pointer to map-based representation
    const std::map<unsigned, Vector<std::pair<FiniteElement*, Vector<double>>>>*
      map_pt = Bin_object_coord_pairs.map_pt();

    // Initialise
    unsigned count = 0;

    // loop...
    typedef std::map<
      unsigned,
      Vector<std::pair<FiniteElement*, Vector<double>>>>::const_iterator IT;
    for (IT it = map_pt->begin(); it != map_pt->end(); it++)
    {
      count += (*it).second.size();
    }
    return count;
  }


  //========================================================================
  /// Output bins
  //========================================================================
  void NonRefineableBinArray::output_bins(std::ofstream& outfile)
  {
    // Spatial dimension of bin
    const unsigned n_lagrangian = this->ndim_zeta();

    unsigned nbin_x = Dimensions_of_bin_array[0];
    unsigned nbin_y = 1;
    if (n_lagrangian > 1) nbin_y = Dimensions_of_bin_array[1];
    unsigned nbin_z = 1;
    if (n_lagrangian > 2) nbin_z = Dimensions_of_bin_array[2];

    unsigned b = 0;
    for (unsigned iz = 0; iz < nbin_z; iz++)
    {
      for (unsigned iy = 0; iy < nbin_y; iy++)
      {
        for (unsigned ix = 0; ix < nbin_x; ix++)
        {
          unsigned nentry = Bin_object_coord_pairs[b].size();
          for (unsigned e = 0; e < nentry; e++)
          {
            FiniteElement* el_pt = Bin_object_coord_pairs[b][e].first;
            Vector<double> s(Bin_object_coord_pairs[b][e].second);
            Vector<double> zeta(n_lagrangian);
            if (Use_eulerian_coordinates_during_setup)
            {
              el_pt->interpolated_x(s, zeta);
            }
            else
            {
              el_pt->interpolated_zeta(s, zeta);
            }
            for (unsigned i = 0; i < n_lagrangian; i++)
            {
              outfile << zeta[i] << " ";
            }
            outfile << ix << " " << iy << " " << iz << " " << b << " "
                    << std::endl;
          }
          b++;
        }
      }
    }
  }


  //========================================================================
  /// Output bin vertices (allowing display of bins as zones).
  //========================================================================
  void NonRefineableBinArray::output_bin_vertices(std::ofstream& outfile)
  {
    // Store the lagrangian dimension
    const unsigned n_lagrangian = this->ndim_zeta();

    unsigned nbin = Dimensions_of_bin_array[0];
    if (n_lagrangian > 1) nbin *= Dimensions_of_bin_array[1];
    if (n_lagrangian > 2) nbin *= Dimensions_of_bin_array[2];

    for (unsigned i_bin = 0; i_bin < nbin; i_bin++)
    {
      // Get bin vertices
      Vector<Vector<double>> bin_vertex;
      get_bin_vertices(i_bin, bin_vertex);
      switch (n_lagrangian)
      {
        case 1:
          outfile << "ZONE I=2\n";
          break;

        case 2:
          outfile << "ZONE I=2, J=2\n";
          break;

        case 3:
          outfile << "ZONE I=2, J=2, K=2\n";
          break;
      }

      unsigned nvertex = bin_vertex.size();
      for (unsigned i = 0; i < nvertex; i++)
      {
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          outfile << bin_vertex[i][j] << " ";
        }
        outfile << std::endl;
      }
    }
  }


  //========================================================================
  /// Fill the bin array with sample points from FiniteElements stored in mesh
  //========================================================================
  void NonRefineableBinArray::fill_bin_array()
  {
    // Store the lagrangian dimension
    const unsigned n_lagrangian = this->ndim_zeta();

    // Flush all objects out of the bin structure
    flush_bins_of_objects();

    // Temporary storage in bin
    std::map<unsigned, Vector<std::pair<FiniteElement*, Vector<double>>>>
      tmp_bin_object_coord_pairs;

    // Total number of bins
    unsigned ntotalbin = nbin();

    // Initialise the structure that keep tracks of the occupided bins
    Bin_object_coord_pairs.initialise(ntotalbin);


    // Issue warning about small number of bins
    if (!Suppress_warning_about_small_number_of_bins)
    {
      // Calculate the (nearest integer) number of elements per bin
      unsigned n_mesh_element = Mesh_pt->nelement();
      unsigned elements_per_bin = n_mesh_element / ntotalbin;
      // If it is above the threshold then issue a warning
      if (elements_per_bin > Threshold_for_elements_per_bin_warning)
      {
        if (!Already_warned_about_small_number_of_bin_cells)
        {
          Already_warned_about_small_number_of_bin_cells = true;
          std::ostringstream warn_message;
          warn_message
            << "The average (integer) number of elements per bin is \n\n"
            << elements_per_bin << ", which is more than \n\n"
            << "   "
               "NonRefineableBinArray::Threshold_for_elements_per_bin_warning="
            << Threshold_for_elements_per_bin_warning
            << "\n\nIf the lookup seems slow (and you have the memory),"
            << "consider increasing\n"
            << "BinArray::Dimensions_of_bin_array from their current\n"
            << "values of { ";
          unsigned nn = Dimensions_of_bin_array.size();
          for (unsigned ii = 0; ii < nn; ii++)
          {
            warn_message << Dimensions_of_bin_array[ii] << " ";
          }
          warn_message
            << " }.\n"
            << "\nNOTE: You can suppress this warning by increasing\n\n"
            << "   NonRefineableBinArray::"
            << "Threshold_for_elements_per_bin_warning\n\n"
            << "or by setting \n\n   "
            << "NonRefineableBinArray::Suppress_warning_about_small_number_of_"
               "bins\n\n"
            << "to true (both are public static data).\n\n";
          OomphLibWarning(warn_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
      }
    }


    // Increase overall counter
    Total_nbin_cells_counter += ntotalbin;


    // Issue warning?
    if (!Suppress_warning_about_large_total_number_of_bins)
    {
      if (Total_nbin_cells_counter >
          Threshold_for_total_bin_cell_number_warning)
      {
        if (!Already_warned_about_large_number_of_bin_cells)
        {
          Already_warned_about_large_number_of_bin_cells = true;
          std::ostringstream warn_message;
          warn_message
            << "Have allocated \n\n"
            << "   NonRefineableBinArray::Total_nbin_cells_counter="
            << NonRefineableBinArray::Total_nbin_cells_counter
            << "\n\nbin cells, which is more than \n\n"
            << "   NonRefineableBinArray::"
            << "Threshold_for_total_bin_cell_number_warning="
            << NonRefineableBinArray::
                 Threshold_for_total_bin_cell_number_warning
            << "\n\nIf you run out of memory, consider reducing\n"
            << "BinArray::Dimensions_of_bin_array from their current\n"
            << "values of { ";
          unsigned nn = Dimensions_of_bin_array.size();
          for (unsigned ii = 0; ii < nn; ii++)
          {
            warn_message << Dimensions_of_bin_array[ii] << " ";
          }
          warn_message
            << " }.\n"
            << "\nNOTE: You can suppress this warning by increasing\n\n"
            << "   NonRefineableBinArray::"
            << "Threshold_for_total_bin_cell_number_warning\n\n"
            << "or by setting \n\n   NonRefineableBinArray::"
            << "Suppress_warning_about_large_total_number_of_bins\n\n"
            << "to true (both are public static data).\n\n"
            << "NOTE: I'll only issue this warning once -- total number of\n"
            << "bins may grow yet further!\n";
          OomphLibWarning(warn_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
      }
    }

    /// Loop over subobjects (elements) to decide which bin they belong in...
    unsigned n_sub = Mesh_pt->nelement();
    for (unsigned e = 0; e < n_sub; e++)
    {
      // Cast to the element (sub-object) first
      FiniteElement* el_pt =
        dynamic_cast<FiniteElement*>(Mesh_pt->finite_element_pt(e));

      // Get specified number of points within the element
      unsigned n_plot_points =
        el_pt->nplot_points(Nsample_points_generated_per_element);

      for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
      {
        // Storage for local and global coordinates
        Vector<double> local_coord(n_lagrangian, 0.0);
        Vector<double> global_coord(n_lagrangian, 0.0);

        // Get local coordinate and interpolate to global
        bool use_equally_spaced_interior_sample_points =
          SamplePointContainer::Use_equally_spaced_interior_sample_points;
        el_pt->get_s_plot(iplot,
                          Nsample_points_generated_per_element,
                          local_coord,
                          use_equally_spaced_interior_sample_points);

        // Now get appropriate global coordinate
        if (Use_eulerian_coordinates_during_setup)
        {
          el_pt->interpolated_x(local_coord, global_coord);
        }
        else
        {
          el_pt->interpolated_zeta(local_coord, global_coord);
        }

        // Which bin are the global coordinates in?
        unsigned bin_number = 0;
        unsigned multiplier = 1;
        // Loop over the dimension
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
#ifdef PARANOID
          if ((global_coord[i] < Min_and_max_coordinates[i].first) ||
              (global_coord[i] > Min_and_max_coordinates[i].second))
          {
            std::ostringstream error_message;
            error_message
              << "Bin sample point " << iplot << " in element " << e << "\n"
              << "is outside bin limits in coordinate direction " << i << ":\n"
              << "Sample point coordinate: " << global_coord[i] << "\n"
              << "Max bin coordinate     : "
              << Min_and_max_coordinates[i].second << "\n"
              << "Min bin coordinate     : " << Min_and_max_coordinates[i].first
              << "\n"
              << "You should either setup the bin boundaries manually\n"
              << "or increase the percentage offset by which the\n"
              << "automatically computed bin limits are increased \n"
              << "beyond their sampled max/mins. This is defined in\n"
              << "the (public) namespace member\n\n"
              << "SamplePointContainer::Percentage_offset \n\n which \n"
              << "currently has the value: "
              << SamplePointContainer::Percentage_offset << "\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

#endif
          unsigned bin_number_i =
            int(Dimensions_of_bin_array[i] *
                ((global_coord[i] - Min_and_max_coordinates[i].first) /
                 (Min_and_max_coordinates[i].second -
                  Min_and_max_coordinates[i].first)));

          // Buffer the case when the global coordinate is the maximum
          // value
          if (bin_number_i == Dimensions_of_bin_array[i])
          {
            bin_number_i -= 1;
          }

          // Add to the bin number
          bin_number += multiplier * bin_number_i;

          // Sort out the multiplier
          multiplier *= Dimensions_of_bin_array[i];
        }

        // Add element-sample local coord pair to the calculated bin
        tmp_bin_object_coord_pairs[bin_number].push_back(
          std::make_pair(el_pt, local_coord));
      }
    }

    // Finally copy filled vectors across -- wiping memory from temporary
    // map as we go along is good and wouldn't be possible if we
    // filled the SparseVector's internal map within that class.
    typedef std::map<
      unsigned,
      Vector<std::pair<FiniteElement*, Vector<double>>>>::iterator IT;
    for (IT it = tmp_bin_object_coord_pairs.begin();
         it != tmp_bin_object_coord_pairs.end();
         it++)
    {
      Bin_object_coord_pairs.set_value((*it).first, (*it).second);
      // Make space immediately
      (*it).second.clear();
    }
  }


  //========================================================================
  /// Provide some stats on the fill level of the associated bin
  //========================================================================
  void NonRefineableBinArray::get_fill_stats(unsigned& n_bin,
                                             unsigned& max_n_entry,
                                             unsigned& min_n_entry,
                                             unsigned& tot_n_entry,
                                             unsigned& n_empty) const
  {
    // Total number of bins
    n_bin = nbin();
    n_empty = n_bin - Bin_object_coord_pairs.nnz();

    // Get pointer to map-based representation
    const std::map<unsigned, Vector<std::pair<FiniteElement*, Vector<double>>>>*
      map_pt = Bin_object_coord_pairs.map_pt();

    // Initialise
    max_n_entry = 0;
    min_n_entry = UINT_MAX;
    tot_n_entry = 0;

    // Do stats
    typedef std::map<
      unsigned,
      Vector<std::pair<FiniteElement*, Vector<double>>>>::const_iterator IT;
    for (IT it = map_pt->begin(); it != map_pt->end(); it++)
    {
      unsigned nentry = (*it).second.size();
      if (nentry > max_n_entry) max_n_entry = nentry;
      if (nentry < min_n_entry) min_n_entry = nentry;
      tot_n_entry += nentry;
    }
  }


  //========================================================================
  /// Fill bin by diffusion, populating each empty bin with the same content
  /// as the first non-empty bin found during a spiral-based search
  /// up to the specified "radius" (default 1)
  //========================================================================
  void NonRefineableBinArray::fill_bin_by_diffusion(
    const unsigned& bin_diffusion_radius)
  {
    // oomph_info << "PROFILING GET_NEIGHBOURING_BINS_HELPER():" << std::endl;
    // profile_get_neighbouring_bins_helper();

    // Loop over all bins to check if they're empty
    std::list<unsigned> empty_bins;
    unsigned n_bin = nbin();
    std::vector<bool> was_empty_until_current_iteration(n_bin, false);
    for (unsigned i = 0; i < n_bin; i++)
    {
      if (Bin_object_coord_pairs[i].size() == 0)
      {
        empty_bins.push_front(i);
        was_empty_until_current_iteration[i] = true;
      }
    }

    // Now keep processing the empty bins until there are none left
    unsigned iter = 0;
    Vector<unsigned> newly_filled_bin;
    while (empty_bins.size() != 0)
    {
      newly_filled_bin.clear();
      newly_filled_bin.reserve(empty_bins.size());
      for (std::list<unsigned>::iterator it = empty_bins.begin();
           it != empty_bins.end();
           it++)
      {
        unsigned bin = (*it);

        // Look for immediate neighbours within the specified "bin radius"
        unsigned level = bin_diffusion_radius;
        Vector<unsigned> neighbour_bin;
        get_neighbouring_bins_helper(bin, level, neighbour_bin);
        unsigned n_neigh = neighbour_bin.size();

        // Find closest pair
        double min_dist = DBL_MAX;
        std::pair<FiniteElement*, Vector<double>> closest_pair;
        for (unsigned i = 0; i < n_neigh; i++)
        {
          unsigned neigh_bin = neighbour_bin[i];

          // Only allow to re-populate from bins that were already filled at
          // previous iteration, otherwise things can progate too fast
          if (!was_empty_until_current_iteration[neigh_bin])
          {
            unsigned nbin_content = Bin_object_coord_pairs[neigh_bin].size();
            for (unsigned j = 0; j < nbin_content; j++)
            {
              FiniteElement* el_pt = Bin_object_coord_pairs[neigh_bin][j].first;
              Vector<double> s(Bin_object_coord_pairs[neigh_bin][j].second);
              Vector<double> x(2);
              el_pt->interpolated_x(s, x);
              // Get minimum distance of sample point from any of the vertices
              // of current bin
              // hierher Louis questions if this is the right thing to do!
              // Matthias agrees
              // with him but hasn't done anything yet...
              // should use the distance to the centre of gravity of
              // the empty bin instead!
              double dist = min_distance(bin, x);
              if (dist < min_dist)
              {
                min_dist = dist;
                closest_pair = Bin_object_coord_pairs[neigh_bin][j];
              }
            }
          }
        }

        // Have we filled the bin?
        if (min_dist != DBL_MAX)
        {
          Vector<std::pair<FiniteElement*, Vector<double>>> new_entry;
          new_entry.push_back(closest_pair);
          Bin_object_coord_pairs.set_value(bin, new_entry);

          // Record that we've filled it.
          newly_filled_bin.push_back(bin);

          // Wipe entry without breaking the linked list (Andrew's trick --
          // nice!)
          std::list<unsigned>::iterator it_to_be_deleted = it;
          it--;
          empty_bins.erase(it_to_be_deleted);
        }
      }

      // Get ready for next iteration on remaining empty bins
      iter++;
      // Now update the vector that records which bins were empty up to now
      unsigned n = newly_filled_bin.size();
      for (unsigned i = 0; i < n; i++)
      {
        was_empty_until_current_iteration[newly_filled_bin[i]] = false;
      }
    }


#ifdef PARANOID
    // Loop over all bins to check if they're empty
    n_bin = nbin();
    for (unsigned i = 0; i < n_bin; i++)
    {
      if (Bin_object_coord_pairs[i].size() == 0)
      {
        std::ostringstream error_message_stream;
        error_message_stream << "Bin " << i << " is still empty\n"
                             << "after trying to fill it by diffusion\n";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

#endif
  }


  //=================================================================
  /// Get the number of the bin containing the specified coordinate.
  /// Bin number is negative if the coordinate is outside
  /// the bin structure.
  //=================================================================
  void NonRefineableBinArray::get_bin(const Vector<double>& zeta,
                                      int& bin_number)
  {
    // Default for point not found
    bin_number = -1;

    // Get the lagrangian dimension
    const unsigned n_lagrangian = this->ndim_zeta();

    // Does the zeta coordinate lie within the current bin structure?

    // Loop over the lagrangian dimension
    for (unsigned i = 0; i < n_lagrangian; i++)
    {
      // If the i-th coordinate is less than the minimum
      if (zeta[i] < Min_and_max_coordinates[i].first)
      {
        return;
      }
      // Otherwise coordinate may be bigger than the maximum
      else if (zeta[i] > Min_and_max_coordinates[i].second)
      {
        return;
      }
    }

    // Use the min and max coords of the bin structure, to find
    // the bin structure containing the current zeta cooordinate

    // Initialise for subsequent computation
    bin_number = 0;

    // Offset for rows/matrices in higher dimensions
    unsigned multiplier = 1;

    // Loop over the dimension
    for (unsigned i = 0; i < n_lagrangian; i++)
    {
      // Find the bin number of the current coordinate
      unsigned bin_number_i =
        int(Dimensions_of_bin_array[i] *
            ((zeta[i] - Min_and_max_coordinates[i].first) /
             (Min_and_max_coordinates[i].second -
              Min_and_max_coordinates[i].first)));
      // Buffer the case when we are exactly on the edge
      if (bin_number_i == Dimensions_of_bin_array[i])
      {
        bin_number_i -= 1;
      }
      // Now add to the bin number using the multiplier
      bin_number += multiplier * bin_number_i;
      // Increase the current row/matrix multiplier for the next loop
      multiplier *= Dimensions_of_bin_array[i];
    }


#ifdef PARANOID

    // Tolerance for "out of bin" test
    double tol = 1.0e-10;

    unsigned nvertex = (int)pow(2, n_lagrangian);
    Vector<Vector<double>> bin_vertex(nvertex);
    for (unsigned j = 0; j < nvertex; j++)
    {
      bin_vertex[j].resize(n_lagrangian);
    }
    get_bin_vertices(bin_number, bin_vertex);
    for (unsigned i = 0; i < n_lagrangian; i++)
    {
      double min_vertex_coord = DBL_MAX;
      double max_vertex_coord = -DBL_MAX;
      for (unsigned j = 0; j < nvertex; j++)
      {
        if (bin_vertex[j][i] < min_vertex_coord)
        {
          min_vertex_coord = bin_vertex[j][i];
        }
        if (bin_vertex[j][i] > max_vertex_coord)
        {
          max_vertex_coord = bin_vertex[j][i];
        }
      }
      if (zeta[i] < min_vertex_coord - tol)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Trouble! " << i << " -th coordinate of sample point, " << zeta[i]
          << " , isn't actually between limits, " << min_vertex_coord << " and "
          << max_vertex_coord << " [it's below by more than " << tol << " !] "
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      if (zeta[i] > max_vertex_coord + tol)
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "Trouble! " << i << " -th coordinate of sample point, " << zeta[i]
          << " , isn't actually between limits, " << min_vertex_coord << " and "
          << max_vertex_coord << " [it's above by more than " << tol << "!] "
          << std::endl;
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
  }


  //========================================================================
  /// Get vector of vectors containing the coordinates of the
  /// vertices of the i_bin-th bin: bin_vertex[j][i] contains the
  /// i-th coordinate of the j-th vertex.
  //========================================================================
  void NonRefineableBinArray::get_bin_vertices(
    const unsigned& i_bin, Vector<Vector<double>>& bin_vertex)
  {
    // Spatial dimension of bin
    const unsigned n_lagrangian = this->ndim_zeta();

    // How many vertices do we have?
    unsigned n_vertices = 1;
    for (unsigned i = 0; i < n_lagrangian; i++)
    {
      n_vertices *= 2;
    }
    bin_vertex.resize(n_vertices);

    // Vectors for min [0] and max [1] coordinates of the bin in each
    // coordinate direction
    Vector<Vector<double>> zeta_vertex_bin(2);
    zeta_vertex_bin[0].resize(n_lagrangian);
    zeta_vertex_bin[1].resize(n_lagrangian);

    Vector<double> dzeta;
    unsigned count = 0;
    Vector<unsigned> i_1d;

    // Deal with different spatial dimensions separately
    switch (n_lagrangian)
    {
      case 1:

        // Assign vertex coordinates of the bin directly
        dzeta.resize(1);
        dzeta[0] = (Min_and_max_coordinates[0].second -
                    Min_and_max_coordinates[0].first) /
                   double(Dimensions_of_bin_array[0]);
        bin_vertex[0].resize(1);
        bin_vertex[0][0] =
          Min_and_max_coordinates[0].first + double(i_bin) * dzeta[0];
        bin_vertex[1].resize(1);
        bin_vertex[1][0] =
          Min_and_max_coordinates[0].first + double(i_bin + 1) * dzeta[0];

        break;

      case 2:

        dzeta.resize(2);
        dzeta[0] = (Min_and_max_coordinates[0].second -
                    Min_and_max_coordinates[0].first) /
                   double(Dimensions_of_bin_array[0]);
        dzeta[1] = (Min_and_max_coordinates[1].second -
                    Min_and_max_coordinates[1].first) /
                   double(Dimensions_of_bin_array[1]);

        i_1d.resize(2);
        i_1d[0] = i_bin % Dimensions_of_bin_array[0];
        i_1d[1] = (i_bin - i_1d[0]) / Dimensions_of_bin_array[0];

        // Max/min coordinates of the bin
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
          zeta_vertex_bin[0][i] =
            Min_and_max_coordinates[i].first + double(i_1d[i]) * dzeta[i];
          zeta_vertex_bin[1][i] =
            Min_and_max_coordinates[i].first + double(i_1d[i] + 1) * dzeta[i];
        }

        // Build 4 vertex coordinates
        count = 0;
        for (unsigned i_min_max = 0; i_min_max < 2; i_min_max++)
        {
          for (unsigned j_min_max = 0; j_min_max < 2; j_min_max++)
          {
            bin_vertex[count].resize(2);
            bin_vertex[count][0] = zeta_vertex_bin[i_min_max][0];
            bin_vertex[count][1] = zeta_vertex_bin[j_min_max][1];
            count++;
          }
        }

        break;

      case 3:

        dzeta.resize(3);
        dzeta[0] = (Min_and_max_coordinates[0].second -
                    Min_and_max_coordinates[0].first) /
                   double(Dimensions_of_bin_array[0]);
        dzeta[1] = (Min_and_max_coordinates[1].second -
                    Min_and_max_coordinates[1].first) /
                   double(Dimensions_of_bin_array[1]);
        dzeta[2] = (Min_and_max_coordinates[2].second -
                    Min_and_max_coordinates[2].first) /
                   double(Dimensions_of_bin_array[2]);

        i_1d.resize(3);
        i_1d[0] = i_bin % Dimensions_of_bin_array[0];
        i_1d[1] = ((i_bin - i_1d[0]) / Dimensions_of_bin_array[0]) %
                  Dimensions_of_bin_array[1];
        i_1d[2] = (i_bin - (i_1d[1] * Dimensions_of_bin_array[0] + i_1d[0])) /
                  (Dimensions_of_bin_array[0] * Dimensions_of_bin_array[1]);

        // Max/min coordinates of the bin
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
          zeta_vertex_bin[0][i] =
            Min_and_max_coordinates[i].first + double(i_1d[i]) * dzeta[i];
          zeta_vertex_bin[1][i] =
            Min_and_max_coordinates[i].first + double(i_1d[i] + 1) * dzeta[i];
        }

        // Build 8 vertex coordinates
        count = 0;
        for (unsigned i_min_max = 0; i_min_max < 2; i_min_max++)
        {
          for (unsigned j_min_max = 0; j_min_max < 2; j_min_max++)
          {
            for (unsigned k_min_max = 0; k_min_max < 2; k_min_max++)
            {
              bin_vertex[count].resize(3);
              bin_vertex[count][0] = zeta_vertex_bin[i_min_max][0];
              bin_vertex[count][1] = zeta_vertex_bin[j_min_max][1];
              bin_vertex[count][2] = zeta_vertex_bin[k_min_max][2];
              count++;
            }
          }
        }

        break;

      default:

        std::ostringstream error_message;
        error_message << "Can't deal with bins in dimension " << n_lagrangian
                      << "\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Compute the minimum distance of any vertex in the specified bin
  /// from the specified Lagrangian coordinate zeta.
  //========================================================================
  double NonRefineableBinArray::min_distance(const unsigned& i_bin,
                                             const Vector<double>& zeta)
  {
    // Spatial dimension of bin
    const unsigned n_lagrangian = this->ndim_zeta();

    // Get bin vertices
    Vector<Vector<double>> bin_vertex;
    get_bin_vertices(i_bin, bin_vertex);
    double min_dist = DBL_MAX;
    unsigned nvertex = bin_vertex.size();
    for (unsigned v = 0; v < nvertex; v++)
    {
      double dist = 0.0;
      for (unsigned i = 0; i < n_lagrangian; i++)
      {
        dist += pow(bin_vertex[v][i] - zeta[i], 2);
      }
      dist = sqrt(dist);
      if (dist < min_dist) min_dist = dist;
    }
    return min_dist;
  }


  //==============================================================================
  ///  Find the sub geometric object and local coordinate therein that
  /// corresponds to the intrinsic coordinate zeta. If sub_geom_object_pt=0
  /// on return from this function, none of the constituent sub-objects
  /// contain the required coordinate.
  //==============================================================================
  void NonRefineableBinArray::locate_zeta(const Vector<double>& zeta,
                                          GeomObject*& sub_geom_object_pt,
                                          Vector<double>& s)
  {
    // Reset counter for number of sample points visited only if we're
    // starting from the beginning
    if (current_min_spiral_level() == 0)
    {
      Total_number_of_sample_points_visited_during_locate_zeta_from_top_level =
        0;
    }

    // Initialise return to null -- if it's still null when we're
    // leaving we've failed!
    sub_geom_object_pt = 0;

    // Get the lagrangian dimension
    const unsigned n_lagrangian = this->ndim_zeta();

    // Does the zeta coordinate lie within the current bin structure?
    // Skip this test if we want to always fail because that's usually
    // done to trace out the spiral path
    if (!BinArray::Always_fail_elemental_locate_zeta)
    {
      // Loop over the lagrangian dimension
      for (unsigned i = 0; i < n_lagrangian; i++)
      {
        // If the i-th coordinate is less than the minimum
        if (zeta[i] < Min_and_max_coordinates[i].first)
        {
          return;
        }
        // Otherwise coordinate may be bigger than the maximum
        else if (zeta[i] > Min_and_max_coordinates[i].second)
        {
          return;
        }
      }
    }

    // Use the min and max coords of the bin structure, to find
    // the bin structure containing the current zeta cooordinate
    unsigned bin_number = 0;

    // Offset for rows/matrices in higher dimensions
    unsigned multiplier = 1;

    // Loop over the dimension
    for (unsigned i = 0; i < n_lagrangian; i++)
    {
      // Find the bin number of the current coordinate
      unsigned bin_number_i =
        int(Dimensions_of_bin_array[i] *
            ((zeta[i] - Min_and_max_coordinates[i].first) /
             (Min_and_max_coordinates[i].second -
              Min_and_max_coordinates[i].first)));

      // Buffer the case when we are exactly on the edge
      if (bin_number_i == Dimensions_of_bin_array[i])
      {
        bin_number_i -= 1;
      }

      // Now add to the bin number using the multiplier
      bin_number += multiplier * bin_number_i;

      // Increase the current row/matrix multiplier for the next loop
      multiplier *= Dimensions_of_bin_array[i];
    }

    // Loop over spirals that are to be visited in one go
    Vector<unsigned> neighbour_bin;
    unsigned min_level = current_min_spiral_level();
    unsigned max_level = current_max_spiral_level();
    for (unsigned i_level = min_level; i_level <= max_level; i_level++)
    {
      // Call helper function to find the neighbouring bins at this
      // level
      get_neighbouring_bins_helper(bin_number, i_level, neighbour_bin);

      // Total number of bins to be visited
      unsigned n_nbr_bin = neighbour_bin.size();

      // // keep around and add "false" as last argument to previous call
      // // to get_neighbouring_bins_helper(...) for annotation below to make
      // sense
      // {
      //  Vector<unsigned> old_neighbour_bin;
      //  get_neighbouring_bins_helper(bin_number,
      //                               i_level,
      //                               old_neighbour_bin,
      //                               true); // true=use old version!
      //  unsigned nbin_new = neighbour_bin.size();
      //  unsigned nbin_old = old_neighbour_bin.size();
      //  unsigned n=std::min(nbin_new,nbin_old);
      //  std::sort(old_neighbour_bin.begin(),
      //            old_neighbour_bin.end());
      //  std::sort(neighbour_bin.begin(),
      //            neighbour_bin.end());
      //  oomph_info << "# of neighbouring bins: "
      //             << nbin_new << " " << nbin_old;
      //  if (nbin_new!=nbin_old)
      //   {
      //    oomph_info << " DIFFERENT!";
      //   }
      //  oomph_info << std::endl;
      //  for (unsigned i=0;i<n;i++)
      //   {
      //    oomph_info << neighbour_bin[i] << " "
      //               << old_neighbour_bin[i] << " ";
      //    if (neighbour_bin[i]!=
      //        old_neighbour_bin[i])
      //     {
      //      oomph_info << "DIFFFFERENT";
      //     }
      //    oomph_info << std::endl;
      //   }
      //  if (nbin_new>nbin_old)
      //   {
      //    for (unsigned i=n;i<nbin_new;i++)
      //     {
      //      oomph_info << "NNEW:" << neighbour_bin[i]
      //                 << std::endl;
      //     }
      //   }
      //  if (nbin_old>nbin_new)
      //   {
      //    for (unsigned i=n;i<nbin_old;i++)
      //     {
      //      oomph_info << "OOLD:" << old_neighbour_bin[i]
      //                 << std::endl;
      //     }
      //   }
      // }


      // Set bool for finding zeta
      bool found_zeta = false;
      for (unsigned i_nbr = 0; i_nbr < n_nbr_bin; i_nbr++)
      {
        // Only search if bin is within the max. radius but min_distance
        // is expensive...
        bool do_it = true;
        if (Max_search_radius < DBL_MAX)
        {
          if (min_distance(neighbour_bin[i_nbr], zeta) > Max_search_radius)
          {
            do_it = false;
          }
        }
        if (do_it)
        {
          // Get the number of element-sample point pairs in this bin
          unsigned n_sample =
            Bin_object_coord_pairs[neighbour_bin[i_nbr]].size();

          // Don't do anything if this bin has no sample points
          if (n_sample > 0)
          {
            for (unsigned i_sam = 0; i_sam < n_sample; i_sam++)
            {
              // Get the element
              FiniteElement* el_pt =
                Bin_object_coord_pairs[neighbour_bin[i_nbr]][i_sam].first;

              // Get the local coordinate
              s = Bin_object_coord_pairs[neighbour_bin[i_nbr]][i_sam].second;

              // History of sample points visited
              if (BinArray::Visited_sample_points_file.is_open())
              {
                unsigned cached_dim_zeta = this->ndim_zeta();
                Vector<double> zeta_sample(cached_dim_zeta);
                if (this->use_eulerian_coordinates_during_setup())
                {
                  el_pt->interpolated_x(s, zeta_sample);
                }
                else
                {
                  el_pt->interpolated_zeta(s, zeta_sample);
                }
                double dist = 0.0;
                for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
                {
                  BinArray::Visited_sample_points_file << zeta_sample[ii]
                                                       << " ";
                  dist +=
                    (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
                }
                BinArray::Visited_sample_points_file
                  << total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
                  << " " << sqrt(dist) << std::endl;
              }

              // Use this coordinate as the initial guess
              bool use_coordinate_as_initial_guess = true;

              // Attempt to find zeta within a sub-object
              el_pt->locate_zeta(
                zeta, sub_geom_object_pt, s, use_coordinate_as_initial_guess);


              Total_number_of_sample_points_visited_during_locate_zeta_from_top_level++;

              // Always fail? (Used for debugging, e.g. to trace out
              // spiral path)
              if (BinArray::Always_fail_elemental_locate_zeta)
              {
                sub_geom_object_pt = 0;
              }


#ifdef OOMPH_HAS_MPI
              // Ignore halos?
              if (Ignore_halo_elements_during_locate_zeta_search)
              {
                // Dynamic cast the result to a FiniteElement
                FiniteElement* test_el_pt =
                  dynamic_cast<FiniteElement*>(sub_geom_object_pt);
                if (test_el_pt != 0)
                {
                  // We only want to exit if this is a non-halo element
                  if (test_el_pt->is_halo())
                  {
                    sub_geom_object_pt = 0;
                  }
                }
              }
#endif

              // If the FiniteElement is non-halo and has been located, exit
              if (sub_geom_object_pt != 0)
              {
                found_zeta = true;
                break;
              }
            } // end loop over sample points
          }


          if (found_zeta)
          {
            return; // break;
          }

        } // end of don't search if outside search radius
      } // end loop over bins at this level
    } // End of loop over levels
  }


  /// Default number of bins (in each coordinate direction)
  unsigned NonRefineableBinArray::Default_n_bin_1d = 100;

  ///  Counter for overall number of bins allocated -- used to
  /// issue warning if this exceeds a threshhold. (Default assignment
  /// of 100^DIM bins per MeshAsGeomObject can be a killer if there
  /// are huge numbers of sub-meshes (e.g. in unstructured FSI).
  unsigned long NonRefineableBinArray::Total_nbin_cells_counter = 0;

  ///  Total number of bins above which warning is issued.
  /// (Default assignment of 100^DIM bins per MeshAsGeomObject can
  /// be a killer if there are huge numbers of sub-meshes (e.g. in
  /// unstructured FSI).
  unsigned long
    NonRefineableBinArray::Threshold_for_total_bin_cell_number_warning =
      50000000;

  ///  Boolean to supppress warnings about large number of bins
  bool
    NonRefineableBinArray::Suppress_warning_about_large_total_number_of_bins =
      false;

  ///  Boolean flag to make sure that warning about large number
  /// of bin cells only gets triggered once.
  bool NonRefineableBinArray::Already_warned_about_large_number_of_bin_cells =
    false;

  ///  Fraction of elements/bin that triggers warning. Too many
  /// elements per bin can lead to very slow computations
  unsigned NonRefineableBinArray::Threshold_for_elements_per_bin_warning = 100;

  ///  Boolean to supppress warnings about small number of bins
  bool NonRefineableBinArray::Suppress_warning_about_small_number_of_bins =
    false;

  ///  Boolean flag to make sure that warning about small number
  /// of bin cells only gets triggered once.
  bool NonRefineableBinArray::Already_warned_about_small_number_of_bin_cells =
    false;


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_CGAL


  //====================================================================
  /// Constructor
  //====================================================================
  CGALSamplePointContainer::CGALSamplePointContainer(
    SamplePointContainerParameters* sample_point_container_parameters_pt)
    : SamplePointContainer(
        sample_point_container_parameters_pt->mesh_pt(),
        sample_point_container_parameters_pt->min_and_max_coordinates(),
        sample_point_container_parameters_pt
          ->use_eulerian_coordinates_during_setup(),
        sample_point_container_parameters_pt
          ->ignore_halo_elements_during_locate_zeta_search(),
        sample_point_container_parameters_pt
          ->nsample_points_generated_per_element())
  {
    // Get the spatial dimension (int because of mpi below)
    int dim = 0;
    if (Mesh_pt->nelement() != 0)
    {
      dim = Mesh_pt->finite_element_pt(0)->dim();
    }

    // Need to do an Allreduce to ensure that the dimension is consistent
    // even when no elements are assigned to a certain processor
#ifdef OOMPH_HAS_MPI
    // Only a problem if the mesh has been distributed
    if (Mesh_pt->is_mesh_distributed())
    {
      // Need a non-null communicator
      if (Mesh_pt->communicator_pt() != 0)
      {
        int n_proc = Mesh_pt->communicator_pt()->nproc();
        if (n_proc > 1)
        {
          int dim_reduce;
          MPI_Allreduce(&dim,
                        &dim_reduce,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        Mesh_pt->communicator_pt()->mpi_comm());
          dim = dim_reduce;
        }
      }
    }
#endif

    Ndim_zeta = dim;

    // Have we specified max/min coordinates?
    // If not, compute them on the fly from mesh
    if (Min_and_max_coordinates.size() == 0)
    {
      setup_min_and_max_coordinates();
    }


    // Time it
    double t_start = 0.0;
    if (SamplePointContainer::Enable_timing_of_setup)
    {
      t_start = TimingHelpers::timer();
    }

    // Fill the bastard!
    double CGAL_setup_time = get_sample_points();

    if (SamplePointContainer::Enable_timing_of_setup)
    {
      double t_end = TimingHelpers::timer();
      unsigned npts = total_number_of_sample_points_computed_recursively();
      oomph_info << "Time for setup of " << dim
                 << "-dimensional sample point container containing " << npts
                 << " sample points: " << t_end - t_start
                 << " sec (cgal); third party: " << CGAL_setup_time
                 << " sec ( = " << CGAL_setup_time / (t_end - t_start) * 100.0
                 << " %)" << std::endl;
    }

    // Initialise
    Total_number_of_sample_points_visited_during_locate_zeta_from_top_level = 0;
    First_sample_point_to_actually_lookup_during_locate_zeta = 0;
    Last_sample_point_to_actually_lookup_during_locate_zeta = UINT_MAX;
    Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta =
      2; // hierher tune this and create public static default
    Initial_last_sample_point_to_actually_lookup_during_locate_zeta =
      10; // UINT MAX is temporary bypass! tune this and create public static
          // default
  }

  //==============================================================================
  /// Get the sample points; returns time taken for setup of CGAL tree
  //==============================================================================
  double CGALSamplePointContainer::get_sample_points()
  {
    // Number of elements
    unsigned nel = Mesh_pt->nelement();

    // Estimate number of sample points
    unsigned n_sample_estimate = 0;
    if (nel > 0)
    {
      FiniteElement* el_pt = Mesh_pt->finite_element_pt(0);
      if (el_pt != 0)
      {
        // Total number of sample point we will create
        n_sample_estimate =
          nel * el_pt->nplot_points(Nsample_points_generated_per_element);
      }
    }
    CGAL_sample_point_zeta_d.reserve(n_sample_estimate);
    Sample_point_pt.reserve(n_sample_estimate);

    // Fill 'em in:
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* el_pt = Mesh_pt->finite_element_pt(e);

      // Total number of sample point we will create
      unsigned nplot =
        el_pt->nplot_points(Nsample_points_generated_per_element);

      /// For all the sample points we have to create ...
      for (unsigned j = 0; j < nplot; j++)
      {
        // ... create it: Pass element index in mesh (vector
        // of elements and index of sample point within element
        SamplePoint* new_sample_point_pt = new SamplePoint(e, j);

        // Coordinates of this point
        Vector<double> zeta(ndim_zeta());
        Vector<double> s(ndim_zeta());
        bool use_equally_spaced_interior_sample_points =
          SamplePointContainer::Use_equally_spaced_interior_sample_points;
        el_pt->get_s_plot(j,
                          Nsample_points_generated_per_element,
                          s,
                          use_equally_spaced_interior_sample_points);
        if (Use_eulerian_coordinates_during_setup)
        {
          el_pt->interpolated_x(s, zeta);
        }
        else
        {
          el_pt->interpolated_zeta(s, zeta);
        }

#ifdef PARANOID

        // Check if point is inside
        bool is_inside = true;
        std::ostringstream error_message;
        unsigned dim = ndim_zeta();
        for (unsigned i = 0; i < dim; i++)
        {
          if ((zeta[i] < Min_and_max_coordinates[i].first) ||
              (zeta[i] > Min_and_max_coordinates[i].second))
          {
            is_inside = false;
            error_message << "Sample point at zeta[" << i << "]  = " << zeta[i]
                          << " is outside limits of sample point container: "
                          << Min_and_max_coordinates[i].first << " and "
                          << Min_and_max_coordinates[i].second << std::endl;
          }
        }

        if (!is_inside)
        {
          error_message << "Please correct the limits passed to the "
                        << "constructor." << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

#endif

        CGAL_sample_point_zeta_d.push_back(
          Kernel_d::Point_d(zeta.size(), zeta.begin(), zeta.end()));
        Sample_point_pt.push_back(new_sample_point_pt);
      }
    }

    // Make tree structure
    double CGAL_setup_time = 0.0;
    if (SamplePointContainer::Enable_timing_of_setup)
    {
      CGAL_setup_time = TimingHelpers::timer();
    }

    CGAL_tree_d_pt = new K_neighbor_search_d::Tree(
      boost::make_zip_iterator(boost::make_tuple(
        CGAL_sample_point_zeta_d.begin(), Sample_point_pt.begin())),
      boost::make_zip_iterator(boost::make_tuple(CGAL_sample_point_zeta_d.end(),
                                                 Sample_point_pt.end())));
    if (SamplePointContainer::Enable_timing_of_setup)
    {
      CGAL_setup_time = TimingHelpers::timer() - CGAL_setup_time;
    }
    return CGAL_setup_time;
  }


  //==============================================================================
  /// Compute total number of sample points in sample point container
  //==============================================================================
  unsigned CGALSamplePointContainer::
    total_number_of_sample_points_computed_recursively() const
  {
    return Sample_point_pt.size();
  }


  //==============================================================================
  ///  Find the sub geometric object and local coordinate therein that
  /// corresponds to the intrinsic coordinate zeta. If sub_geom_object_pt=0
  /// on return from this function, none of the constituent sub-objects
  /// contain the required coordinate.
  //==============================================================================
  void CGALSamplePointContainer::locate_zeta(const Vector<double>& zeta,
                                             GeomObject*& sub_geom_object_pt,
                                             Vector<double>& s)
  {
    // Top level book keeping and sanity checking
    if (first_sample_point_to_actually_lookup_during_locate_zeta() == 0)
    {
      // Reset counter for number of sample points visited.
      // If we can't find the point we should at least make sure that
      // we've visited all the sample points before giving up.
      Total_number_of_sample_points_visited_during_locate_zeta_from_top_level =
        0;
    }

    // Initialise return to null -- if it's still null when we're
    // leaving we've failed!
    sub_geom_object_pt = 0;

    // Get the lagrangian dimension
    const unsigned n_lagrangian = this->ndim_zeta();

    // Does the zeta coordinate lie within the current bin structure?
    // Skip this test if we want to always fail because that's usually
    // done to trace out the spiral path
    if (!SamplePointContainer::Always_fail_elemental_locate_zeta)
    {
      // Loop over the lagrangian dimension
      for (unsigned i = 0; i < n_lagrangian; i++)
      {
        // If the i-th coordinate is less than the minimum
        if (zeta[i] < Min_and_max_coordinates[i].first)
        {
          return;
        }
        // Otherwise coordinate may be bigger than the maximum
        else if (zeta[i] > Min_and_max_coordinates[i].second)
        {
          return;
        }
      }
    }

    // Number of sample points
    unsigned n_sample = Sample_point_pt.size();

    // Create CGAL query -- this is the point we want!
    Point_d query(zeta.size(), zeta.begin(), zeta.end());

    // Max. number of nearest neighbours
    const unsigned n_nearest_neighbours_max = std::min(
      n_sample, last_sample_point_to_actually_lookup_during_locate_zeta());

    // Start with just one...
    unsigned n_nearest_neighbours = 1;
    unsigned n_neighbours_visited_last_time = 0;
    bool can_increase_n_nearest_neighbours = true;
    bool keep_going = true;
    while (keep_going)
    {
      // Find specified number of nearest neighbours only
      const unsigned n_nearest_neighbours_actual =
        std::min(n_nearest_neighbours, n_nearest_neighbours_max);
      K_neighbor_search_d search(
        *CGAL_tree_d_pt, query, n_nearest_neighbours_actual);

      // Search
      unsigned count = 0;
      for (K_neighbor_search_d::iterator it = search.begin();
           it != search.end();
           it++)
      {
        count++;

        if ((count > n_neighbours_visited_last_time) &&
            (count >
             first_sample_point_to_actually_lookup_during_locate_zeta()))
        {
          // Recover the sample point
          SamplePoint* sample_point_pt = boost::get<1>(it->first);

          // Get the element
          FiniteElement* el_pt = Mesh_pt->finite_element_pt(
            sample_point_pt->element_index_in_mesh());


#ifdef OOMPH_HAS_MPI
          // We only look at the sample point if it isn't halo
          // if we are set up to ignore the halo elements
          if (ignore_halo_elements_during_locate_zeta_search() &&
              (el_pt->is_halo()))
          {
            // Halo
          }
          else
          {
#endif

            // Provide initial guess for Newton search using local coordinate
            // of sample point
            bool use_equally_spaced_interior_sample_points =
              SamplePointContainer::Use_equally_spaced_interior_sample_points;
            unsigned i = sample_point_pt->sample_point_index_in_element();
            el_pt->get_s_plot(i,
                              nsample_points_generated_per_element(),
                              s,
                              use_equally_spaced_interior_sample_points);

            bool do_it = true;
            if (Max_search_radius < DBL_MAX)
            {
              unsigned cached_dim_zeta = ndim_zeta();
              Vector<double> zeta_sample(cached_dim_zeta);
              if (use_eulerian_coordinates_during_setup())
              {
                el_pt->interpolated_x(s, zeta_sample);
              }
              else
              {
                el_pt->interpolated_zeta(s, zeta_sample);
              }
              double dist_sq = 0.0;
              for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
              {
                dist_sq +=
                  (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
              }
              if (dist_sq > Max_search_radius * Max_search_radius)
              {
                do_it = false;
              }
            }
            if (do_it)
            {
              // History of sample points visited
              if (SamplePointContainer::Visited_sample_points_file.is_open())
              {
                unsigned cached_dim_zeta = ndim_zeta();
                Vector<double> zeta_sample(cached_dim_zeta);
                if (use_eulerian_coordinates_during_setup())
                {
                  el_pt->interpolated_x(s, zeta_sample);
                }
                else
                {
                  el_pt->interpolated_zeta(s, zeta_sample);
                }
                double dist = 0.0;
                for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
                {
                  SamplePointContainer::Visited_sample_points_file
                    << zeta_sample[ii] << " ";
                  dist +=
                    (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
                }
                SamplePointContainer::Visited_sample_points_file
                  << total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
                  << " " << sqrt(dist) << std::endl;
              }

              // Bump counter
              total_number_of_sample_points_visited_during_locate_zeta_from_top_level()++;

              bool use_coordinate_as_initial_guess = true;
              el_pt->locate_zeta(
                zeta, sub_geom_object_pt, s, use_coordinate_as_initial_guess);

              // Always fail? (Used for debugging, e.g. to trace out
              // spiral path)
              if (SamplePointContainer::Always_fail_elemental_locate_zeta)
              {
                sub_geom_object_pt = 0;
              }

              if (sub_geom_object_pt != 0)
              {
                return;
              }
            }

#ifdef OOMPH_HAS_MPI
          }
#endif
        }
      }

      n_neighbours_visited_last_time = count;

      // Can we increase the number of neighbours further?
      if (can_increase_n_nearest_neighbours)
      {
        unsigned factor_for_increase_in_nearest_neighbours = 10;
        n_nearest_neighbours *= factor_for_increase_in_nearest_neighbours;

        // There's no point going any further (next time)
        if (n_nearest_neighbours > n_nearest_neighbours_max)
        {
          can_increase_n_nearest_neighbours = false;
        }
      }
      // Bailing out; not found but we can't increase number of search pts
      // further
      else
      {
        keep_going = false;
      }
    } // while loop to increase number of nearest neighbours
  }


  //==============================================================================
  ///  Find the sub geometric object and local coordinate therein that
  /// corresponds to the intrinsic coordinate zeta, using up to the specified
  /// number of sample points as initial guess for the Newton-based search.
  /// If this fails, return the nearest sample point.
  //==============================================================================
  void CGALSamplePointContainer::limited_locate_zeta(
    const Vector<double>& zeta,
    const unsigned& max_sample_points_for_newton_based_search,
    GeomObject*& sub_geom_object_pt,
    Vector<double>& s)
  {
    // Reset counter for number of sample points visited.
    Total_number_of_sample_points_visited_during_locate_zeta_from_top_level = 0;

    // Initialise return to null -- if it's still null when we're
    // leaving we've failed!
    sub_geom_object_pt = 0;

    // Number of sample points
    unsigned n_sample = Sample_point_pt.size();

    // Create CGAL query -- this is the point we want!
    Point_d query(zeta.size(), zeta.begin(), zeta.end());

    // Max. number of nearest neighbours
    const unsigned n_nearest_neighbours_max =
      std::min(n_sample, max_sample_points_for_newton_based_search);

    // Find 'em
    K_neighbor_search_d search(
      *CGAL_tree_d_pt, query, n_nearest_neighbours_max);

    // Do Newton method starting from each of the nearest sample points
    for (K_neighbor_search_d::iterator it = search.begin(); it != search.end();
         it++)
    {
      // Recover the sample point
      SamplePoint* sample_point_pt = boost::get<1>(it->first);

      // Get the element
      FiniteElement* el_pt =
        Mesh_pt->finite_element_pt(sample_point_pt->element_index_in_mesh());


#ifdef OOMPH_HAS_MPI

      // We only look at the sample point if it isn't halo
      // if we are set up to ignore the halo elements
      if (ignore_halo_elements_during_locate_zeta_search() &&
          (el_pt->is_halo()))
      {
        // Halo
      }
      else
      { // not halo

#endif

        // Provide initial guess for Newton search using local coordinate
        // of sample point
        bool use_equally_spaced_interior_sample_points =
          SamplePointContainer::Use_equally_spaced_interior_sample_points;
        unsigned i = sample_point_pt->sample_point_index_in_element();
        el_pt->get_s_plot(i,
                          nsample_points_generated_per_element(),
                          s,
                          use_equally_spaced_interior_sample_points);

        bool do_it = true;
        if (Max_search_radius < DBL_MAX)
        {
          unsigned cached_dim_zeta = ndim_zeta();
          Vector<double> zeta_sample(cached_dim_zeta);
          if (use_eulerian_coordinates_during_setup())
          {
            el_pt->interpolated_x(s, zeta_sample);
          }
          else
          {
            el_pt->interpolated_zeta(s, zeta_sample);
          }
          double dist_sq = 0.0;
          for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
          {
            dist_sq +=
              (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
          }
          if (dist_sq > Max_search_radius * Max_search_radius)
          {
            do_it = false;
          }
        }
        if (do_it)
        {
          // History of sample points visited
          if (SamplePointContainer::Visited_sample_points_file.is_open())
          {
            unsigned cached_dim_zeta = ndim_zeta();
            Vector<double> zeta_sample(cached_dim_zeta);
            if (use_eulerian_coordinates_during_setup())
            {
              el_pt->interpolated_x(s, zeta_sample);
            }
            else
            {
              el_pt->interpolated_zeta(s, zeta_sample);
            }
            double dist = 0.0;
            for (unsigned ii = 0; ii < cached_dim_zeta; ii++)
            {
              SamplePointContainer::Visited_sample_points_file
                << zeta_sample[ii] << " ";
              dist +=
                (zeta[ii] - zeta_sample[ii]) * (zeta[ii] - zeta_sample[ii]);
            }
            SamplePointContainer::Visited_sample_points_file
              << total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
              << " " << sqrt(dist) << std::endl;
          }

          // Bump counter
          total_number_of_sample_points_visited_during_locate_zeta_from_top_level()++;

          bool use_coordinate_as_initial_guess = true;
          el_pt->locate_zeta(
            zeta, sub_geom_object_pt, s, use_coordinate_as_initial_guess);

          // Always fail? (Used for debugging, e.g. to trace out
          // spiral path)
          if (SamplePointContainer::Always_fail_elemental_locate_zeta)
          {
            sub_geom_object_pt = 0;
          }

          if (sub_geom_object_pt != 0)
          {
            return;
          }
        }

#ifdef OOMPH_HAS_MPI
      }
#endif
    }

    // We've searched over all the sample points but the Newton method
    // hasn't converged from any, so just use the nearest one
    K_neighbor_search_d::iterator it = search.begin();

    // Recover the sample point
    SamplePoint* sample_point_pt = boost::get<1>(it->first);

    // Get the element
    FiniteElement* el_pt =
      Mesh_pt->finite_element_pt(sample_point_pt->element_index_in_mesh());

    // Get local coordinate of sample point in element
    bool use_equally_spaced_interior_sample_points =
      SamplePointContainer::Use_equally_spaced_interior_sample_points;
    unsigned i = sample_point_pt->sample_point_index_in_element();
    el_pt->get_s_plot(i,
                      nsample_points_generated_per_element(),
                      s,
                      use_equally_spaced_interior_sample_points);

    sub_geom_object_pt = el_pt;
  }


#endif // cgal

} // namespace oomph

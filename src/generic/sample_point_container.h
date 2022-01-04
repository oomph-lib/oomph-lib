// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_SAMPLE_POINT_CONTAINER_HEADER
#define OOMPH_SAMPLE_POINT_CONTAINER_HEADER

#ifdef OOMPH_HAS_CGAL

#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#endif

// oomph-lib includes
#include "sample_point_parameters.h"
#include "sparse_vector.h"

/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////


//=============================================================================
/// Class for containing sample points: Number of finite element in
/// its mesh and index of sample point within that element.
//=============================================================================
class SamplePoint
{
public:
  /// Construct SamplePoint object from number of finite element
  /// in its mesh, and index of sample point within that element
  SamplePoint(const unsigned& element_index_in_mesh,
              const unsigned& sample_point_index_in_element)
    : Element_index_in_mesh(element_index_in_mesh),
      Sample_point_index_in_element(sample_point_index_in_element)
  {
  }

  /// Broken copy constructor.
  SamplePoint(const SamplePoint& data) = delete;

  /// Broken assignment operator.
  void operator=(const SamplePoint&) = delete;

  /// Access function to the index of finite element in its mesh
  unsigned element_index_in_mesh() const
  {
    return Element_index_in_mesh;
  }

  /// Index of sample point within element
  unsigned sample_point_index_in_element() const
  {
    return Sample_point_index_in_element;
  }


private:
  /// Index of finite element in its mesh
  unsigned Element_index_in_mesh;

  /// Index of the sample point within element
  unsigned Sample_point_index_in_element;
};


/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////


//  Forward declaration of the RefineableBinArray class.
class RefineableBinArray;


//==============================================================================
/// RefineableBin class. Contains sample points and is embedded in a
/// RefineableBinArray. May itself be represented by a RefineableBinArray to
/// make it recursive.
//==============================================================================
class RefineableBin
{
public:
  /// Constructor. Pass pointer to bin array that
  /// contains this bin and the index of the newly created bin in that
  /// RefineableBinArray
  RefineableBin(RefineableBinArray* bin_array_pt,
                const unsigned& bin_index_in_bin_array)
    : Sample_point_pt(0),
      Sub_bin_array_pt(0),
      Bin_array_pt(bin_array_pt),
      Bin_index_in_bin_array(bin_index_in_bin_array)
  {
  }


  /// Broken copy constructor.
  RefineableBin(const RefineableBin& data) = delete;

  /// Broken assignment operator.
  void operator=(const RefineableBin&) = delete;

  /// Destructor
  ~RefineableBin();

  /// Compute total number of sample points recursively
  unsigned total_number_of_sample_points_computed_recursively() const;

  /// Add a new sample point to RefineableBin
  void add_sample_point(SamplePoint* new_sample_point_pt,
                        const Vector<double>& zeta_coordinates);

  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  void locate_zeta(const Vector<double>& zeta,
                   GeomObject*& sub_geom_object_pt,
                   Vector<double>& s);

  /// Output bin; x,[y,[z]],n_sample_points
  void output(std::ofstream& outfile, const bool& don_t_recurse = false);

  /// Output bin vertices (allowing display of bins as zones).
  void output_bins(std::ofstream& outfile);

  /// Output bin vertices (allowing display of bins as zones).
  void output_bin_vertices(std::ofstream& outfile);

  /// Number of sample points stored in bin
  unsigned nsample_points_in_bin()
  {
    if (Sample_point_pt == 0)
    {
      return 0;
    }
    else
    {
      return Sample_point_pt->size();
    }
  }

protected:
  /// Container of SamplePoints. Pointer to vector because it's shorter
  /// than an empty vector! (Not all RefineableBins have sample
  /// points -- the ones that are subdivided don't!)
  Vector<SamplePoint*>* Sample_point_pt;

  /// Pointer to a possible sub-BinArray. Null by default
  RefineableBinArray* Sub_bin_array_pt;

  /// Pointer to the bin array which "owns" this RefineableBin.
  RefineableBinArray* Bin_array_pt;

  /// Index of bin in its bin array
  unsigned Bin_index_in_bin_array;

  /// Method for building a new subbin_array (called when the Bin size
  /// is greater than the Max_number_of_sample_point_per_bin (and the Bin is
  /// recursive). Pass in the extremal coordinates of the bin which is
  /// being subdivided. Redistributes all existing sample points to
  /// newly made sub-bin-array and empties its own storage. Pass
  /// Max./min. coordinates of new bin array for efficiency.
  void make_sub_bin_array(
    const Vector<std::pair<double, double>>& min_and_max_coordinates);

  /// Boundaries of bin in each coordinate direction.
  /// *.first = min; *.second = max
  void get_bin_boundaries(
    Vector<std::pair<double, double>>& min_and_max_coordinates);
};


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//=========================================================================
/// Base class for all sample point containers
//=========================================================================
class SamplePointContainer
{
public:
  /// Constructor
  SamplePointContainer(
    Mesh* mesh_pt,
    const Vector<std::pair<double, double>>& min_and_max_coordinates,
    const bool& use_eulerian_coordinates_during_setup,
    const bool& ignore_halo_elements_during_locate_zeta_search,
    const unsigned& nsample_points_generated_per_element)
    : Mesh_pt(mesh_pt),
      Min_and_max_coordinates(min_and_max_coordinates),
      Use_eulerian_coordinates_during_setup(
        use_eulerian_coordinates_during_setup),
#ifdef OOMPH_HAS_MPI
      Ignore_halo_elements_during_locate_zeta_search(
        ignore_halo_elements_during_locate_zeta_search),
#endif
      Nsample_points_generated_per_element(
        nsample_points_generated_per_element),
      Total_number_of_sample_points_visited_during_locate_zeta_from_top_level(0)
  {
    // Don't limit max. search radius
    Max_search_radius = DBL_MAX;
  }

  /// Broken default constructor; needed for broken
  /// copy constructors. Don't call. It will die.
  SamplePointContainer()
  {
    // Throw the error
    throw OomphLibError("Broken default constructor. Don't call this!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  /// Broken copy constructor.
  SamplePointContainer(const SamplePointContainer& data) = delete;

  /// Broken assignment operator.
  void operator=(const SamplePointContainer&) = delete;

  /// Virtual destructor
  virtual ~SamplePointContainer() {}

  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  virtual void locate_zeta(const Vector<double>& zeta,
                           GeomObject*& sub_geom_object_pt,
                           Vector<double>& s) = 0;


  /// Counter to keep track of how many sample points we've
  /// visited during top level call to locate_zeta. Virtual so it can be
  /// overloaded for different versions.
  virtual unsigned& total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
  {
    return Total_number_of_sample_points_visited_during_locate_zeta_from_top_level;
  }

  /// Total number of sample points in sample point container, possibly
  /// computed recursively.
  virtual unsigned total_number_of_sample_points_computed_recursively()
    const = 0;

  /// Dimension of the zeta ( =  dim of local coordinate of elements)
  virtual unsigned ndim_zeta() const = 0;

  /// Pointer to mesh from whose FiniteElements sample points are created
  Mesh* mesh_pt() const
  {
    return Mesh_pt;
  }

  /// Pair of doubles for min and maximum coordinates in i-th direction:
  /// min (first) and max. (second) coordinates
  const std::pair<double, double>& min_and_max_coordinates(
    const unsigned& i) const
  {
    return Min_and_max_coordinates[i];
  }

  /// Vector of pair of doubles for min and maximum coordinates.
  /// min (first) and max. (second) coordinates
  const Vector<std::pair<double, double>>& min_and_max_coordinates() const
  {
    return Min_and_max_coordinates;
  }


#ifdef OOMPH_HAS_MPI

  /// Ignore halo elements?
  bool ignore_halo_elements_during_locate_zeta_search() const
  {
    return Ignore_halo_elements_during_locate_zeta_search;
  }

#endif

  /// Use Eulerian coordinates (i.e. interpolated_x) rather than
  /// zeta itself (i.e. interpolated_zeta) to identify point.
  bool use_eulerian_coordinates_during_setup() const
  {
    return Use_eulerian_coordinates_during_setup;
  }

  /// "Measure of" number of sample points generated in each element
  unsigned& nsample_points_generated_per_element()
  {
    return Nsample_points_generated_per_element;
  }

  /// Set maximum search radius for locate zeta. This is initialised
  /// do DBL_MAX so we brutally search through the entire bin structure,
  /// no matter how big it is until we've found the required point (or
  /// failed to do so. This can be VERY costly with fine meshes.
  /// Here the user takes full responsibility and states that we have
  /// no chance in hell to find the required point in
  /// a bin whose closest vertex is further than the specified
  /// max search radius.
  double& max_search_radius()
  {
    return Max_search_radius;
  }


  /// File to record sequence of visited sample points in. Used for debugging/
  /// illustration of search procedures.
  static std::ofstream Visited_sample_points_file;

  /// Boolean flag to make to make locate zeta fail. Used for debugging/
  /// illustration of search procedures.
  static bool Always_fail_elemental_locate_zeta;

  /// Use equally spaced sample points? (otherwise vertices are sampled
  /// repeatedly
  static bool Use_equally_spaced_interior_sample_points;

  /// Time setup?
  static bool Enable_timing_of_setup;

  /// Offset of sample point container boundaries beyond max/min coords
  static double Percentage_offset;

protected:
  /// Helper function to compute the min and max coordinates for the
  /// mesh, in each dimension
  void setup_min_and_max_coordinates();

  /// Pointer to mesh from whose FiniteElements sample points are created
  Mesh* Mesh_pt;

  /// Vector of pairs of doubles for min and maximum coordinates.
  /// Call: Min_and_max_coordinates[j] gives me the
  /// pair of min (first) and max. (second) coordinates in the j-th
  /// coordinate direction.
  Vector<std::pair<double, double>> Min_and_max_coordinates;

  /// Use Eulerian coordinates (i.e. interpolated_x) rather than
  /// zeta itself (i.e. interpolated_zeta) to identify point.
  bool Use_eulerian_coordinates_during_setup;

#ifdef OOMPH_HAS_MPI

  /// Ignore halo elements?
  bool Ignore_halo_elements_during_locate_zeta_search;

#endif

  /// "Measure of" number of sample points generated in each element
  unsigned Nsample_points_generated_per_element;

  /// Counter to keep track of how many sample points we've
  /// visited during top level call to locate_zeta
  unsigned
    Total_number_of_sample_points_visited_during_locate_zeta_from_top_level;

  /// Max radius beyond which we stop searching the bin. Initialised
  /// to DBL_MAX so keep going until the point is found or until
  /// we've searched every single bin. Overwriting this means we won't search
  /// in bins whose closest vertex is at a distance greater than
  /// Max_search_radius from the point to be located.
  double Max_search_radius;
};


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=========================================================================
/// Base class for all bin arrays
//=========================================================================
class BinArray : public virtual SamplePointContainer
{
public:
  /// Constructor
  BinArray(Mesh* mesh_pt,
           const Vector<std::pair<double, double>>& min_and_max_coordinates,
           const Vector<unsigned>& dimensions_of_bin_array,
           const bool& use_eulerian_coordinates_during_setup,
           const bool& ignore_halo_elements_during_locate_zeta_search,
           const unsigned& nsample_points_generated_per_element)
    : SamplePointContainer(mesh_pt,
                           min_and_max_coordinates,
                           use_eulerian_coordinates_during_setup,
                           ignore_halo_elements_during_locate_zeta_search,
                           nsample_points_generated_per_element),
      Dimensions_of_bin_array(dimensions_of_bin_array)
  {
    // Note: Resizing of Dimensions_of_bin_array if no sizes are specified
    // is delayed to derived class since refineable and nonrefineable
    // bin arrays have different defaults.
  }

  /// Broken default constructor; needed for broken
  /// copy constructors. Don't call. It will die.
  BinArray()
  {
    // Throw the error
    throw OomphLibError("Broken default constructor. Don't call this!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  /// Broken copy constructor.
  BinArray(const BinArray& data) = delete;

  /// Broken assignment operator.
  void operator=(const BinArray&) = delete;

  /// Virtual destructor
  virtual ~BinArray() {}

  /// Helper function for computing the bin indices of neighbouring bins
  /// at a given "radius" of the specified bin. Final, optional boolean
  /// (default: true) chooses to use the old version which appears to be
  /// faster than Louis' new one after all (in the few cases where this
  /// functionality is still used -- not all if we have cgal!)
  void get_neighbouring_bins_helper(const unsigned& bin_index,
                                    const unsigned& radius,
                                    Vector<unsigned>& neighbouring_bin_index,
                                    const bool& use_old_version = true);

  /// Profiling function to compare performance of two different
  /// versions of the get_neighbouring_bins_helper(...) function
  void profile_get_neighbouring_bins_helper();


  /// Get (linearly enumerated) bin index of bin that
  /// contains specified zeta
  unsigned coords_to_bin_index(const Vector<double>& zeta);


  ///  Get "coordinates" of bin that contains specified zeta
  void coords_to_vectorial_bin_index(const Vector<double>& zeta,
                                     Vector<unsigned>& bin_index);

  /// Output bins (boundaries and number of sample points in them)
  virtual void output_bins(std::ofstream& outfile) = 0;

  /// Output bin vertices (allowing display of bin boundaries as zones).
  virtual void output_bin_vertices(std::ofstream& outfile) = 0;

  /// Number of bins (not taking recursion into account for refineable
  /// versions)
  virtual unsigned nbin() const = 0;

  /// Max. bin dimension (number of bins in coordinate directions)
  unsigned max_bin_dimension() const;

  /// Dimension of the zeta ( =  dim of local coordinate of elements)
  unsigned ndim_zeta() const
  {
    return Dimensions_of_bin_array.size();
  }

  /// Number of bins in coordinate direction i
  unsigned dimension_of_bin_array(const unsigned& i) const
  {
    return Dimensions_of_bin_array[i];
  }


  /// Number of bins in coordinate directions. Const vector-based
  /// version
  Vector<unsigned> dimensions_of_bin_array() const
  {
    return Dimensions_of_bin_array;
  }


  /// Number of bins in specified coordinate direction
  unsigned dimensions_of_bin_array(const unsigned& i) const
  {
    return Dimensions_of_bin_array[i];
  }

protected:
  /// Number of bins in each coordinate direction
  Vector<unsigned> Dimensions_of_bin_array;
};


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==============================================================================
/// RefineableBinArray class.
//==============================================================================
class RefineableBinArray : public virtual BinArray
{
public:
  /// Constructor
  RefineableBinArray(SamplePointContainerParameters* bin_array_parameters_pt);

  /// Broken copy constructor.
  RefineableBinArray(const RefineableBinArray& data) = delete;

  /// Broken assignment operator.
  void operator=(const RefineableBinArray&) = delete;

  /// Destructor
  ~RefineableBinArray()
  {
    unsigned n = Bin_pt.size();
    for (unsigned i = 0; i < n; i++)
    {
      if (Bin_pt[i] != 0)
      {
        delete Bin_pt[i];
        Bin_pt[i] = 0;
      }
    }
  }

  /// Root bin array
  RefineableBinArray* root_bin_array_pt() const
  {
    return Root_bin_array_pt;
  }

  /// Pointer to i-th bin; can be null if bin is empty
  RefineableBin* bin_pt(const unsigned& i) const
  {
    return Bin_pt[i];
  }

  /// Number of bins (not taking recursion into account)
  unsigned nbin() const
  {
    return Bin_pt.size();
  }

  /// Default number of bins (in each coordinate direction)
  /// (Note: don't move this into a common base class because
  /// each derived class has its own value; we'll want far fewer
  /// in the refineable version!)
  static unsigned Default_n_bin_1d;

  /// Compute total number of sample points recursively
  unsigned total_number_of_sample_points_computed_recursively() const;

  /// Fill the bin array with specified SamplePoints
  void fill_bin_array(const Vector<SamplePoint*>& sample_point_pt)
  {
    unsigned n_dim = ndim_zeta();

    unsigned n = sample_point_pt.size();
    for (unsigned i = 0; i < n; i++)
    {
      // Coordinates of this point
      Vector<double> zeta(n_dim);

      // Which element is the point in?
      unsigned e = sample_point_pt[i]->element_index_in_mesh();
      FiniteElement* el_pt = Mesh_pt->finite_element_pt(e);

      // Which sample point is it at?
      unsigned j = sample_point_pt[i]->sample_point_index_in_element();
      Vector<double> s(n_dim);
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

      // Add it
      add_sample_point(sample_point_pt[i], zeta);
    }
  }

  /// Add specified SamplePoint to RefineableBinArray
  void add_sample_point(SamplePoint* new_sample_point_pt,
                        const Vector<double>& zeta)
  {
    // Find the correct bin
    unsigned bin_index = coords_to_bin_index(zeta);

    // if the bin is not yet created, create it...
    if (Bin_pt[bin_index] == 0)
    {
      Bin_pt[bin_index] = new RefineableBin(this, bin_index);
    }
    // Then add the SamplePoint
    Bin_pt[bin_index]->add_sample_point(new_sample_point_pt, zeta);
  }

  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  void locate_zeta(const Vector<double>& zeta,
                   GeomObject*& sub_geom_object_pt,
                   Vector<double>& s);

  /// Boundaries of specified bin in each coordinate direction.
  /// *.first = min; *.second = max.
  void get_bin_boundaries(
    const unsigned& bin_index,
    Vector<std::pair<double, double>>& min_and_max_coordinates);

  /// Depth of the hierarchical bin_array.
  unsigned depth() const
  {
    return Depth;
  }

  /// Max depth of the hierarchical bin_array; const version
  unsigned max_depth() const
  {
    return Max_depth;
  }

  /// Max depth of the hierarchical bin_array
  unsigned& max_depth()
  {
    return Max_depth;
  }

  /// Is the BinArray recursive?
  bool bin_array_is_recursive() const
  {
    return Bin_array_is_recursive;
  }

  /// Maximum number of sample points in bin (before its subdivided
  /// recursively)
  unsigned max_number_of_sample_point_per_bin() const
  {
    return Max_number_of_sample_point_per_bin;
  }

  /// Output bins
  void output_bins(std::ofstream& outfile)
  {
    /// Loop over bins
    unsigned n_bin = Bin_pt.size();
    for (unsigned i = 0; i < n_bin; i++)
    {
      if (Bin_pt[i] != 0)
      {
        Bin_pt[i]->output(outfile);
      }
    }
  }

  /// Output bin vertices (allowing display of bins as zones).
  void output_bin_vertices(std::ofstream& outfile);

  /// Output neighbouring bins at given "radius" of the specified bin
  void output_neighbouring_bins(const unsigned& bin_index,
                                const unsigned& radius,
                                std::ofstream& outfile);


  /// Counter to keep track of how many sample points we've
  /// visited during top level call to locate_zeta
  unsigned& total_number_of_sample_points_visited_during_locate_zeta_from_top_level()
  {
    if (Depth == 0)
    {
      return Total_number_of_sample_points_visited_during_locate_zeta_from_top_level;
    }
    else
    {
      return Root_bin_array_pt
        ->total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
    }
  }

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  unsigned& first_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return First_sample_point_to_actually_lookup_during_locate_zeta;
  }

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter is less than this value.
  unsigned& last_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Last_sample_point_to_actually_lookup_during_locate_zeta;
  }

  /// Every time we've completed a "spiral", visiting a finite
  /// number of sample points in a deterministic order, use this
  /// multiplier to increase the max. number of sample points to be visited.
  /// Using a multiplier rather than a constant increment increases
  /// the amount of (more and more unlikely to yield anything!) work
  /// done locally before doing another costly mpi round trip
  /// when we're already far from the point we're trying to find.
  unsigned& multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta;
  }

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  /// This is the initial value when starting the spiral based search.
  unsigned& initial_last_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Initial_last_sample_point_to_actually_lookup_during_locate_zeta;
  }


private:
  /// Fill the bin array with sample points from FiniteElements stored in mesh
  void fill_bin_array();

  /// Loop over all sample points in the element specified via the
  /// pointer and create a SamplePoint for each. Also specify the index of the
  /// element in its mesh.
  void create_sample_points_from_element(FiniteElement* const element_pt,
                                         const unsigned& n_element);


  /// Vector of pointers to constituent RefineableBins.
  Vector<RefineableBin*> Bin_pt;

  /// Variable which stores if the RefineableBinArray is
  /// recursive or not.
  bool Bin_array_is_recursive;

  /// Variable which stores the Depth value of the bin_array. Useful for
  /// debugging and for preventing "infinite" recursion in case if there is
  /// a problem.
  unsigned Depth;

  /// Max depth of the hierarchical bin_array
  unsigned Max_depth;

  /// Maximum number of sample points in bin (before it's subdivided
  /// recursively)
  unsigned Max_number_of_sample_point_per_bin;

  /// Pointer to root bin array
  RefineableBinArray* Root_bin_array_pt;

  // hierher only used in root

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  unsigned First_sample_point_to_actually_lookup_during_locate_zeta;

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter is less than this value.
  unsigned Last_sample_point_to_actually_lookup_during_locate_zeta;

  /// Every time we've completed a "spiral", visiting a finite
  /// number of sample points in a deterministic order, use this
  /// multiplier to increase the max. number of sample points to be visited.
  /// Using a multiplier rather than a constant increment increases
  /// the amount of (more and more unlikely to yield anything!) work
  /// done locally before doing another costly mpi round trip
  /// when we're already far from the point we're trying to find.
  unsigned
    Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta;

  /// When searching through sample points recursively from the top
  /// level RefineableBinArray (in deterministic order!) only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  /// This is the initial value when starting the spiral based search.
  unsigned Initial_last_sample_point_to_actually_lookup_during_locate_zeta;
};


/// /////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////


//==============================================================================
/// NonRefineableBinArray class.
//==============================================================================
class NonRefineableBinArray : public virtual BinArray
{
public:
  /// Constructor
  NonRefineableBinArray(
    SamplePointContainerParameters* bin_array_parameters_pt);

  /// Destructor:
  ~NonRefineableBinArray()
  {
    flush_bins_of_objects();
  }

  /// Broken copy constructor.
  NonRefineableBinArray(const NonRefineableBinArray& data) = delete;

  /// Broken assignment operator.
  void operator=(const NonRefineableBinArray&) = delete;

  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  void locate_zeta(const Vector<double>& zeta,
                   GeomObject*& sub_geom_object_pt,
                   Vector<double>& s);

  /// Total number of bins (empty or not)
  unsigned nbin() const
  {
    const unsigned n_lagrangian = ndim_zeta();
    unsigned ntotalbin = Dimensions_of_bin_array[0];
    for (unsigned i = 1; i < n_lagrangian; i++)
    {
      ntotalbin *= Dimensions_of_bin_array[i];
    }
    return ntotalbin;
  }

  /// Default number of bins (in each coordinate direction).
  /// (Note: don't move this into a common base class because
  /// each derived class has its own value; nonrefineable bin
  /// wants a much larger value than the refineable one!)
  static unsigned Default_n_bin_1d;

  /// Compute total number of sample points recursively
  unsigned total_number_of_sample_points_computed_recursively() const;

  /// Number of spirals to be searched in one go const version
  unsigned n_spiral_chunk() const
  {
    return Nspiral_chunk;
  }

  /// Number of spirals to be searched in one go
  unsigned& n_spiral_chunk()
  {
    return Nspiral_chunk;
  }

  /// Access function to max. spiral level during straight locate_zeta
  /// search (for efficiency; similar to max_search_radius())
  unsigned& max_spiral_level()
  {
    return Max_spiral_level;
  }

  /// Access function to current min. spiral level
  unsigned& current_min_spiral_level()
  {
    return Current_min_spiral_level;
  }

  /// Access function to current max. spiral level
  unsigned& current_max_spiral_level()
  {
    return Current_max_spiral_level;
  }

  /// Provide some stats on the fill level of the associated bin
  void get_fill_stats(unsigned& n_bin,
                      unsigned& max_n_entry,
                      unsigned& min_n_entry,
                      unsigned& tot_n_entry,
                      unsigned& n_empty) const;

  /// Compute the minimum distance of any vertex in the specified bin
  /// from the specified Lagrangian coordinate zeta
  double min_distance(const unsigned& i_bin, const Vector<double>& zeta);

  /// Output bin vertices (allowing display of bins as zones).
  void output_bin_vertices(std::ofstream& outfile);

  /// Get vector of vectors containing the coordinates of the
  /// vertices of the i_bin-th bin: bin_vertex[j][i] contains the
  /// i-th coordinate of the j-th vertex.
  void get_bin_vertices(const unsigned& i_bin,
                        Vector<Vector<double>>& bin_vertex);

  /// Get the number of the bin containing the specified coordinate.
  /// Bin number is negative if the coordinate is outside
  /// the bin structure.
  void get_bin(const Vector<double>& zeta, int& bin_number);

  /// Get the number of the bin containing the specified coordinate; also
  /// return the contents of that bin. Bin number is negative if the
  /// coordinate is outside the bin structure.
  void get_bin(
    const Vector<double>& zeta,
    int& bin_number,
    Vector<std::pair<FiniteElement*, Vector<double>>>& sample_point_pairs);

  /// Get the contents of all bins in vector
  Vector<Vector<std::pair<FiniteElement*, Vector<double>>>> bin_content() const
  {
    Vector<Vector<std::pair<FiniteElement*, Vector<double>>>> all_vals;
    Bin_object_coord_pairs.get_all_values(all_vals);
    return all_vals;
  }

  /// Get the contents of all bins in vector
  const std::map<unsigned, Vector<std::pair<FiniteElement*, Vector<double>>>>* get_all_bins_content()
    const
  {
    // Return the content of the bins
    return Bin_object_coord_pairs.map_pt();
  }

  /// Fill bin by diffusion, populating each empty bin with the
  /// same content as the first non-empty bin found during a spiral-based search
  /// up to the specified "radius" (default 1)
  void fill_bin_by_diffusion(const unsigned& bin_diffusion_radius = 1);


  /// Output bins
  void output_bins(std::ofstream& outfile);

  /// Output bins
  void output_bins(std::string& filename)
  {
    std::ofstream outfile;
    outfile.open(filename.c_str());
    output_bins(outfile);
    outfile.close();
  }

  /// Counter for overall number of bins allocated -- used to
  /// issue warning if this exceeds a threshhold. (Default assignment
  /// of 100^DIM bins per MeshAsGeomObject can be a killer if there
  /// are huge numbers of sub-meshes (e.g. in unstructured FSI).
  static unsigned long Total_nbin_cells_counter;

  /// Total number of bins above which warning is issued.
  /// (Default assignment of 100^DIM bins per MeshAsGeomObject can
  /// be a killer if there are huge numbers of sub-meshes (e.g. in
  /// unstructured FSI).
  static unsigned long Threshold_for_total_bin_cell_number_warning;

  /// Boolean to supppress warnings about large number of bins
  static bool Suppress_warning_about_large_total_number_of_bins;

  /// Boolean flag to make sure that warning about large number
  /// of bin cells only gets triggered once.
  static bool Already_warned_about_large_number_of_bin_cells;

  /// Fraction of elements/bin that triggers warning. Too many
  /// elements per bin can lead to very slow computations
  static unsigned Threshold_for_elements_per_bin_warning;

  /// Boolean to supppress warnings about small number of bins
  static bool Suppress_warning_about_small_number_of_bins;

  /// Boolean flag to make sure that warning about small number
  /// of bin cells only gets triggered once.
  static bool Already_warned_about_small_number_of_bin_cells;

private:
  /// Fill the bin array with sample points from FiniteElements stored in mesh
  void fill_bin_array();

  /// Initialise and populate the "bin" structure for locating coordinates
  /// and increment counter for total number of bins in active use by any
  /// MeshAsGeomObject)
  void create_bins_of_objects();

  /// Flush the storage for the binning method (and decrement counter
  /// for total number of bins in active use by any MeshAsGeomObject)
  void flush_bins_of_objects()
  {
    Total_nbin_cells_counter -= Bin_object_coord_pairs.nnz();
    Bin_object_coord_pairs.clear();
  }

  /// Storage for paired objects and coords in each bin
  SparseVector<Vector<std::pair<FiniteElement*, Vector<double>>>>
    Bin_object_coord_pairs;

  /// Max. spiralling level (for efficiency; effect similar to
  /// max_search_radius)
  unsigned Max_spiral_level;

  /// Current min. spiralling level
  unsigned Current_min_spiral_level;

  /// Current max. spiralling level
  unsigned Current_max_spiral_level;

  /// Number of spirals to be searched in one go
  unsigned Nspiral_chunk;
};


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_CGAL


//====================================================================
/// CGAL-based SamplePointContainer
//====================================================================
class CGALSamplePointContainer : public virtual SamplePointContainer
{
public:
  /// Constructor
  CGALSamplePointContainer(
    SamplePointContainerParameters* sample_point_container_parameters_pt);

  /// Broken copy constructor.
  CGALSamplePointContainer(const CGALSamplePointContainer& data) = delete;

  /// Broken assignment operator.
  void operator=(const CGALSamplePointContainer&) = delete;

  /// Virtual destructor
  virtual ~CGALSamplePointContainer()
  {
    unsigned n = Sample_point_pt.size();
    for (unsigned i = 0; i < n; i++)
    {
      delete Sample_point_pt[i];
      Sample_point_pt[i] = 0;
    }
    delete CGAL_tree_d_pt;
    CGAL_tree_d_pt = 0;
  }

  /// When searching through sample points only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  unsigned& first_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return First_sample_point_to_actually_lookup_during_locate_zeta;
  }

  /// When searching through sample points  only actually do the
  /// locate_zeta calls when when the counter is less than this value.
  unsigned& last_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Last_sample_point_to_actually_lookup_during_locate_zeta;
  }


  /// Every time we've completed a "spiral", visiting a finite
  /// number of sample points in a deterministic order, use this
  /// multiplier to increase the max. number of sample points to be visited.
  /// Using a multiplier rather than a constant increment increases
  /// the amount of (more and more unlikely to yield anything!) work
  /// done locally before doing another costly mpi round trip
  /// when we're already far from the point we're trying to find.
  unsigned& multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta;
  }

  /// When searching through sample points only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  /// This is the initial value when starting the spiral based search.
  unsigned& initial_last_sample_point_to_actually_lookup_during_locate_zeta()
  {
    return Initial_last_sample_point_to_actually_lookup_during_locate_zeta;
  }


  /// Find sub-GeomObject (finite element) and the local coordinate
  /// s within it that contains point with global coordinate zeta.
  /// sub_geom_object_pt=0 if point can't be found.
  void locate_zeta(const Vector<double>& zeta,
                   GeomObject*& sub_geom_object_pt,
                   Vector<double>& s);


  /// Find the sub geometric object and local coordinate therein that
  /// corresponds to the intrinsic coordinate zeta, using up to the specified
  /// number of sample points as initial guess for the Newton-based search.
  /// If this fails, return the nearest sample point.
  void limited_locate_zeta(
    const Vector<double>& zeta,
    const unsigned& max_sample_points_for_newton_based_search,
    GeomObject*& sub_geom_object_pt,
    Vector<double>& s);


  /// Dimension of the zeta ( =  dim of local coordinate of elements)
  unsigned ndim_zeta() const
  {
    return Ndim_zeta;
  }

  /// Compute total number of sample points in sample point container
  unsigned total_number_of_sample_points_computed_recursively() const;

private:
  /// Get the sample points; return time for setup of CGAL tree.
  double get_sample_points();

  /// Dimension of the zeta ( =  dim of local coordinate of elements)
  unsigned Ndim_zeta;

  /// typedefs for cgal stuff
  typedef CGAL::Cartesian_d<double> Kernel_d;
  typedef Kernel_d::Point_d Point_d;
  typedef boost::tuple<Point_d, SamplePoint*> Point_d_and_pointer;
  typedef CGAL::Search_traits_d<Kernel_d> Traits_base_d;
  typedef CGAL::Search_traits_adapter<
    Point_d_and_pointer,
    CGAL::Nth_of_tuple_property_map<0, Point_d_and_pointer>,
    Traits_base_d>
    Traits_d;
  typedef CGAL::Orthogonal_k_neighbor_search<Traits_d> K_neighbor_search_d;

  /// Vector containing sample point coordinates
  Vector<Point_d> CGAL_sample_point_zeta_d;

  /// Pointer to tree-based representation of sample points
  K_neighbor_search_d::Tree* CGAL_tree_d_pt;

  /// Vector storing pointers to sample point objects (which represent
  /// sample point in terms of number of element
  /// in its mesh and number of sample point)
  Vector<SamplePoint*> Sample_point_pt;

  /// When searching through sample points only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  unsigned First_sample_point_to_actually_lookup_during_locate_zeta;

  /// When searching through sample points  only actually do the
  /// locate_zeta calls when when the counter is less than this value.
  unsigned Last_sample_point_to_actually_lookup_during_locate_zeta;

  /// Every time we've completed a "spiral", visiting a finite
  /// number of sample points in a deterministic order, use this
  /// multiplier to increase the max. number of sample points to be visited.
  /// Using a multiplier rather than a constant increment increases
  /// the amount of (more and more unlikely to yield anything!) work
  /// done locally before doing another costly mpi round trip
  /// when we're already far from the point we're trying to find.
  unsigned
    Multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta;

  /// When searching through sample points only actually do the
  /// locate_zeta calls when when the counter exceeds this value.
  /// This is the initial value when starting the "spiral based" search.
  unsigned Initial_last_sample_point_to_actually_lookup_during_locate_zeta;
};

#endif // endif oomph has cgal

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


} // end of namespace extension

#endif

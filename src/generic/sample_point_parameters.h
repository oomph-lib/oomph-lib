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
#ifndef SAMPLE_POINT_PARAMETERS_HEADER
#define SAMPLE_POINT_PARAMETERS_HEADER


// oomph-lib includes
#include "elements.h"
#include "mesh.h"

namespace oomph
{
  //=========================================================
  /// Enumeration to identify type of sample point container
  //=========================================================
  enum Sample_Point_Container_Type
  {
    UseRefineableBinArray = 1,
    UseNonRefineableBinArray = 2
#ifdef OOMPH_HAS_CGAL
    ,
    UseCGALSamplePointContainer = 3
#endif
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  // Forward references
  class RefineableBinArray;
  class NonRefineableBinArray;


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Helper object for dealing with the parameters used for the
  /// SamplePointContainer objects
  //=========================================================================
  class SamplePointContainerParameters
  {
  public:
    /// Constructor is private and only accessible by friends
    /// to ensure that parameters are set correctly
    friend class BinArrayParameters;
    friend class RefineableBinArrayParameters;
    friend class NonRefineableBinArrayParameters;
#ifdef OOMPH_HAS_CGAL
    friend class CGALSamplePointContainerParameters;
#endif

    /// Broken copy constructor.
    SamplePointContainerParameters(const SamplePointContainerParameters& data) =
      delete;

    /// Broken assignment operator.
    void operator=(const SamplePointContainerParameters&) = delete;

    /// Empty destructor
    virtual ~SamplePointContainerParameters() {}

    /// Pointer to mesh from whose FiniteElements sample points are created
    Mesh* mesh_pt() const
    {
      return Mesh_pt;
    }

    /// Vector of pairs of doubles for min and maximum coordinates.
    /// Call: Min_and_max_coordinates[j] gives me the
    /// pair of min (first) and max. (second) coordinates in the j-th
    /// coordinate direction.
    Vector<std::pair<double, double>>& min_and_max_coordinates()
    {
      return Min_and_max_coordinates;
    }

    /// Vector of pairs of doubles for min and maximum coordinates.
    /// Call: Min_and_max_coordinates[j] gives me the
    /// pair of min (first) and max. (second) coordinates in the j-th
    /// coordinate direction. Const version
    Vector<std::pair<double, double>> min_and_max_coordinates() const
    {
      return Min_and_max_coordinates;
    }

    /// "Measure of" number of sample points generated in each element
    /// const version
    unsigned nsample_points_generated_per_element() const
    {
      return Nsample_points_generated_per_element;
    }

    /// "Measure of" number of sample points generated in each element
    unsigned& nsample_points_generated_per_element()
    {
      return Nsample_points_generated_per_element;
    }

    /// Use eulerian coordinates (via interpolated_x) during
    /// setup (otherwise use interpolated_zeta())?
    bool use_eulerian_coordinates_during_setup() const
    {
      return Use_eulerian_coordinates_during_setup;
    }

    /// Enable use of eulerian coordinates (via interpolated_x) during
    /// setup (otherwise use interpolated_zeta())
    void enable_use_eulerian_coordinates_during_setup()
    {
      Use_eulerian_coordinates_during_setup = true;
    }

    /// Disable use of eulerian coordinates (via interpolated_x) during
    /// setup (otherwise use interpolated_zeta())
    void disable_use_eulerian_coordinates_during_setup()
    {
      Use_eulerian_coordinates_during_setup = false;
    }

    /// Ignore halo elements? (MPI only)
    bool ignore_halo_elements_during_locate_zeta_search() const
    {
      return Ignore_halo_elements_during_locate_zeta_search;
    }

    /// Enable Ignore halo elements? (MPI only)
    void enable_ignore_halo_elements_during_locate_zeta_search()
    {
      Ignore_halo_elements_during_locate_zeta_search = true;
    }

    /// Disable Ignore halo elements? (MPI only)
    void disable_ignore_halo_elements_during_locate_zeta_search()
    {
      Ignore_halo_elements_during_locate_zeta_search = false;
    }

    /// Default for "measure of" number of sample points per element
    static unsigned Default_nsample_points_generated_per_element;

  protected:
    /// Pointer to mesh from whose FiniteElements sample points are created
    Mesh* Mesh_pt;

    /// Vector of pairs of doubles for min and maximum coordinates.
    /// Call: Min_and_max_coordinates[j] gives me the
    /// pair of min (first) and max. (second) coordinates in the j-th
    /// coordinate direction.
    Vector<std::pair<double, double>> Min_and_max_coordinates;

    /// "Measure of" number of sample points generated in each element
    unsigned Nsample_points_generated_per_element;

    /// Use Eulerian coordinates to setup bin (i.e. use interpolated_x()
    /// rather than interpolated_zeta() when setting up and searching sample
    /// point container)
    bool Use_eulerian_coordinates_during_setup;

    /// Ignore halo elements? Accepting halo elements can drastically
    /// reduce the number of external halo elements in multidomain
    /// problems -- currently not aware of any problems with doing this
    /// therefore set to false by default but retention
    /// of this flag allows easy return to previous implementation.
    bool Ignore_halo_elements_during_locate_zeta_search;


  private:
    /// Constructor: Pass mesh.
    /// Constructor is private and can only be called
    /// by the derived friends.
    SamplePointContainerParameters(Mesh* mesh_pt)
      : Mesh_pt(mesh_pt),
        Nsample_points_generated_per_element(
          Default_nsample_points_generated_per_element),
        Use_eulerian_coordinates_during_setup(false),
        Ignore_halo_elements_during_locate_zeta_search(false)
    {
    }

    /// Broken default constructor; needed for broken
    /// copy constructors. Don't call. It will die.
    SamplePointContainerParameters()
    {
      // Throw the error
      throw OomphLibError("Broken default constructor. Don't call this!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_CGAL

  //=========================================================================
  /// Helper object for dealing with the parameters used for the
  /// CGALSamplePointContainer objects
  //=========================================================================
  class CGALSamplePointContainerParameters
    : public virtual SamplePointContainerParameters
  {
  public:
    /// Constructor: Pass mesh.
    CGALSamplePointContainerParameters(Mesh* mesh_pt)
      : SamplePointContainerParameters(mesh_pt)
    {
    }


    /// Broken copy constructor.
    CGALSamplePointContainerParameters(
      const CGALSamplePointContainerParameters& data) = delete;

    /// Broken assignment operator.
    void operator=(const CGALSamplePointContainerParameters&) = delete;
  };

#endif

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Helper object for dealing with the parameters used for the
  /// BinArray objects
  //=========================================================================
  class BinArrayParameters : public virtual SamplePointContainerParameters
  {
  public:
    /// Constructor is private and only accessible by friends
    friend class RefineableBinArrayParameters;
    friend class NonRefineableBinArrayParameters;

    /// Broken copy constructor.
    BinArrayParameters(const BinArrayParameters& data) = delete;

    /// Broken assignment operator.
    void operator=(const BinArrayParameters&) = delete;

    /// Empty destructor
    virtual ~BinArrayParameters() {}

    /// Number of bins in each coordinate direction
    Vector<unsigned>& dimensions_of_bin_array()
    {
      return Dimensions_of_bin_array;
    }

    /// Number of bins in each coordinate direction. Const version
    Vector<unsigned> dimensions_of_bin_array() const
    {
      return Dimensions_of_bin_array;
    }

  protected:
    /// Number of bins in each coordinate direction
    Vector<unsigned> Dimensions_of_bin_array;

  private:
    /// Constructor: Pass mesh. Constructor is private and can only
    /// be called by the derived friends.
    BinArrayParameters(Mesh* mesh_pt) : SamplePointContainerParameters(mesh_pt)
    {
    }

    /// Broken default constructor; needed for broken
    /// copy constructors. Don't call. It will die.
    BinArrayParameters()
    {
      // Throw the error
      throw OomphLibError("Broken default constructor. Don't call this!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Helper object for dealing with the parameters used for the
  /// RefineableBinArray objects
  //=========================================================================
  class RefineableBinArrayParameters : public virtual BinArrayParameters
  {
  public:
    /// Constructor: Pass mesh
    RefineableBinArrayParameters(Mesh* mesh_pt)
      : SamplePointContainerParameters(mesh_pt),
        BinArrayParameters(mesh_pt),
        Bin_array_is_recursive(true),
        Depth(0),
        Max_depth(Default_max_depth),
        Max_number_of_sample_point_per_bin(
          Default_max_number_of_sample_point_per_bin),
        Root_bin_array_pt(0)
    {
    }

    /// Broken copy constructor.
    RefineableBinArrayParameters(const RefineableBinArrayParameters& data) =
      delete;

    /// Broken assignment operator.
    void operator=(const RefineableBinArrayParameters&) = delete;

    /// Empty destructor
    virtual ~RefineableBinArrayParameters() {}


    /// Is bin recursive?
    bool bin_array_is_recursive() const
    {
      return Bin_array_is_recursive;
    }

    /// Enable recursiveness
    void enable_bin_array_is_recursive()
    {
      Bin_array_is_recursive = true;
    }

    /// Disable recursiveness
    void disable_bin_array_is_recursive()
    {
      Bin_array_is_recursive = false;
    }

    /// Variable which stores the Depth value of the bin_array.
    unsigned& depth()
    {
      return Depth;
    }

    /// Variable which stores the Depth value of the bin_array.
    /// const version
    unsigned depth() const
    {
      return Depth;
    }

    /// Max. depth value of the bin_array.
    unsigned& max_depth()
    {
      return Max_depth;
    }

    /// Max. depth value of the bin_array.
    /// const version
    unsigned max_depth() const
    {
      return Max_depth;
    }

    /// Maximum number of sample points in bin (before it's subdivided
    /// recursively)
    unsigned& max_number_of_sample_point_per_bin()
    {
      return Max_number_of_sample_point_per_bin;
    }

    /// Maximum number of sample points in bin (before it's subdivided
    /// recursively; const version
    unsigned max_number_of_sample_point_per_bin() const
    {
      return Max_number_of_sample_point_per_bin;
    }

    /// Pointer to root bin array
    RefineableBinArray*& root_bin_array_pt()
    {
      return Root_bin_array_pt;
    }

    /// Pointer to root bin array; const version
    RefineableBinArray* root_bin_array_pt() const
    {
      return Root_bin_array_pt;
    }

    /// Default value for max. depth
    static unsigned Default_max_depth;

    /// Default value for max. number of sample points before refinement
    static unsigned Default_max_number_of_sample_point_per_bin;

  private:
    /// Variable which stores if the RefineableBinArray is
    /// recursive or not.
    bool Bin_array_is_recursive;

    /// Variable which stores the Depth value of the bin_array. Useful
    /// for debugging and for preventing "infinite" recursion in case if there
    /// is a problem.
    unsigned Depth;

    /// Max. depth value of the bin_array.
    unsigned Max_depth;

    /// Maximum number of sample points in bin (before its subdivided
    /// recursively
    unsigned Max_number_of_sample_point_per_bin;

    /// Pointer to root bin array
    RefineableBinArray* Root_bin_array_pt;
  };


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //=========================================================================
  /// Helper object for dealing with the parameters used for the
  /// NonRefineableBinArray objects
  //=========================================================================
  class NonRefineableBinArrayParameters : public virtual BinArrayParameters
  {
  public:
    /// Constructor: Pass mesh
    NonRefineableBinArrayParameters(Mesh* mesh_pt)
      : SamplePointContainerParameters(mesh_pt),
        BinArrayParameters(mesh_pt),
        Nspiral_chunk(Default_nspiral_chunk)
    {
    }

    /// Broken copy constructor.
    NonRefineableBinArrayParameters(
      const NonRefineableBinArrayParameters& data) = delete;

    /// Broken assignment operator.
    void operator=(const NonRefineableBinArrayParameters&) = delete;

    /// Empty destructor
    virtual ~NonRefineableBinArrayParameters() {}

    /// Number of spirals that are being
    /// visited before doing another circular mpi communication
    /// const version
    unsigned nspiral_chunk() const
    {
      return Nspiral_chunk;
    }
    /// Number of spirals that are being
    /// visited before doing another circular mpi communication
    unsigned& nspiral_chunk()
    {
      return Nspiral_chunk;
    }

    /// Default value for number of spirals that are being
    /// visited before doing another circular mpi communication
    static unsigned Default_nspiral_chunk;

  private:
    /// Number of spirals that are being
    /// visited before doing another circular mpi communication
    unsigned Nspiral_chunk;
  };

} // namespace oomph

#endif

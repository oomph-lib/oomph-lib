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

#include "sample_point_parameters.h"


namespace oomph
{
  /// Default for "measure of" number of sample points per element
  unsigned SamplePointContainerParameters::
    Default_nsample_points_generated_per_element = 5;

  /// Default value for max. depth
  unsigned RefineableBinArrayParameters::Default_max_depth =
    100; // hierher explore

  /// Default value for max. number of sample points before refinement
  unsigned
    RefineableBinArrayParameters::Default_max_number_of_sample_point_per_bin =
      15; // hierher explore. Can be 1 if points are uniformly spaced

  /// Default value for number of spirals that are being
  /// visited before doing another circular mpi communication
  unsigned NonRefineableBinArrayParameters::Default_nspiral_chunk =
    10; // hierher explore 1;

} // namespace oomph

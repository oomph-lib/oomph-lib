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
// Non-inline function for the Elasticity Tensor
#include "time_harmonic_elasticity_tensor.h"

namespace oomph
{
  /// Translation scheme that takes account of the symmetries of the
  /// tensor. The independent coefficients are related to the coefficients of
  /// the elasticity tensor as follows:
  /**\f[\begin{array}{cc}    0 & E_{1111} \\
                             1 & E_{1112} \\
                             2 & E_{1122} \\
                             3 & E_{1212} \\
                             4 & E_{1222} \\
                             5 & E_{2222} \\
                             6 & E_{1113} \\
                             7 & E_{1123} \\
                             8 & E_{1133} \\
                             9 & E_{1213} \\
                             10 & E_{1223} \\
                             11 & E_{1233} \\
                             12 & E_{1313} \\
                             13 & E_{1322} \\
                             14 & E_{1323} \\
                             15 & E_{1333} \\
                             16 & E_{2223} \\
                             17 & E_{2233} \\
                             18 & E_{2323} \\
                             19 & E_{2333} \\
                             20 & E_{3333}
                             \end{array}\f] **/
  const unsigned TimeHarmonicElasticityTensor::Index[3][3][3][3] = {
    {{{0, 1, 6}, {1, 2, 7}, {6, 7, 8}},
     {{1, 3, 9}, {3, 4, 10}, {9, 10, 11}},
     {{6, 9, 12}, {9, 13, 14}, {12, 14, 15}}},

    {{{1, 3, 9}, {3, 4, 10}, {9, 10, 11}},
     {{2, 4, 13}, {4, 5, 16}, {13, 16, 17}},
     {{7, 10, 14}, {10, 16, 18}, {14, 18, 19}}},

    {{{6, 9, 12}, {9, 13, 14}, {12, 14, 15}},
     {{7, 10, 14}, {10, 16, 18}, {14, 18, 19}},
     {{8, 11, 15}, {11, 17, 19}, {15, 19, 20}}}};


  /// Translation scheme for the isotropic elasticity tensor
  const unsigned TimeHarmonicIsotropicElasticityTensor::StaticIndex[21] = {
    1, 0, 2, 3, 0, 1, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 2, 3, 0, 1};

} // namespace oomph

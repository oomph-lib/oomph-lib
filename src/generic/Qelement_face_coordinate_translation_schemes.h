// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header containing functions that map face coordinates on
// Qelements to bulk coordinates in those elements

// Include guards
#ifndef OOMPH_QELEMENT_FACE_COORDINATE_TRANSLATION_HEADER
#define OOMPH_QELEMENT_FACE_COORDINATE_TRANSLATION_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// oomph-lib headers
#include "Vector.h"
#include "matrices.h"

namespace oomph
{
  //=============================================================
  /// Namespace for helper functions that return the local
  /// coordinates in the bulk elements
  //==============================================================
  namespace QElement1FaceToBulkCoordinates
  {
    /// The translation scheme for the face s0 = -1.0
    void face0(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the face s0 = 1.0
    void face1(const Vector<double> &s, Vector<double> &s_bulk);
  } // namespace QElement1FaceToBulkCoordinates

  //=============================================================
  /// Namespace for helper functions that calculate derivatives
  /// of the local coordinates in the bulk elements wrt the
  /// local coordinates in the face element.
  //=============================================================
  namespace QElement1BulkCoordinateDerivatives
  {
    /// Function for both faces -- the bulk coordinate is fixed on both
    void faces0(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);
  } // namespace QElement1BulkCoordinateDerivatives

  //===================================================================
  /// Namespace for the functions that translate local face coordinates
  /// to the coordinates in the bulk element
  //==================================================================
  namespace QElement2FaceToBulkCoordinates
  {
    /// The translation scheme for the west face (s0 = -1.0)
    void face0(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the south face (s1 = -1.0)
    void face1(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the east face (s0 = 1.0)
    void face2(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the north face (s1 = 1.0)
    void face3(const Vector<double> &s, Vector<double> &s_bulk);
  } // namespace QElement2FaceToBulkCoordinates

  //=============================================================
  /// Namespace for helper functions that calculate derivatives
  /// of the local coordinates in the bulk elements wrt the
  /// local coordinates in the face element.
  //=============================================================
  namespace QElement2BulkCoordinateDerivatives
  {
    /// Function for the east and west faces, along which s0 is fixed
    void faces0(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);

    /// Function for the north and south faces, along which s1 is fixed
    void faces1(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);
  } // namespace QElement2BulkCoordinateDerivatives

  //===================================================================
  /// Namespace for the functions that translate local face coordinates
  /// to the coordinates in the bulk element
  //==================================================================
  namespace QElement3FaceToBulkCoordinates
  {
    /// The translation scheme for the left face (s0 = -1.0)
    void face0(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the down face (s1 = -1.0)
    void face1(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the back face (s2 = -1.0)
    void face2(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the right face (s0 = 1.0)
    void face3(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the up face (s1 = 1.0)
    void face4(const Vector<double> &s, Vector<double> &s_bulk);

    /// The translation scheme for the front face (s2 = 1.0)
    void face5(const Vector<double> &s, Vector<double> &s_bulk);
  } // namespace QElement3FaceToBulkCoordinates

  //=============================================================
  /// Namespace for helper functions that calculate derivatives
  /// of the local coordinates in the bulk elements wrt the
  /// local coordinates in the face element.
  //=============================================================
  namespace QElement3BulkCoordinateDerivatives
  {
    /// Function for the back and front  faces, along which s0 is fixed
    void faces0(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);

    /// Function for the up and down  faces, along which s1 is fixed
    void faces1(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);

    /// Function for the left and right  faces, along which s2 is fixed
    void faces2(const Vector<double> &s,
                DenseMatrix<double> &dsbulk_dsface,
                unsigned &interior_direction);
  } // namespace QElement3BulkCoordinateDerivatives

} // namespace oomph

#endif

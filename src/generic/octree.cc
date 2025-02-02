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
#include <algorithm>

#include "octree.h"
#include "refineable_brick_element.h"

namespace oomph
{
  //========================================================================
  /// Bool indicating that static member data has been setup
  //========================================================================
  bool OcTree::Static_data_has_been_setup = false;


  // Public static member data:
  //---------------------------

  //====================================================================
  /// Translate (enumerated) directions into strings
  //====================================================================
  Vector<std::string> OcTree::Direct_string;

  //====================================================================
  /// Get opposite face, e.g. Reflect_face[L]=R
  //====================================================================
  Vector<int> OcTree::Reflect_face;

  //====================================================================
  /// Get opposite edge, e.g. Reflect_edge[DB]=UF
  //====================================================================
  Vector<int> OcTree::Reflect_edge;

  //====================================================================
  /// Get opposite vertex, e.g. Reflect_vertex[LDB]=RUF
  //====================================================================
  Vector<int> OcTree::Reflect_vertex;

  //====================================================================
  /// Map of vectors containing the two vertices for each edge
  //====================================================================
  Vector<Vector<int>> OcTree::Vertex_at_end_of_edge;

  //====================================================================
  /// Map storing the information to translate a vector of directions
  /// (in the three coordinate directions) into a direction
  //====================================================================
  std::map<Vector<int>, int> OcTree::Vector_to_direction;

  //====================================================================
  /// Vector storing the information to translate a direction into a
  /// vector of directions (in the three coordinate directions)
  //====================================================================
  Vector<Vector<int>> OcTree::Direction_to_vector;

  //====================================================================
  /// Storage for the up/right-equivalents corresponding to two
  /// pairs of vertices along an element edge:
  /// - The first pair contains
  ///   -# the vertex in the reference element
  ///   -# the corresponding vertex in the edge neighbour (i.e. the
  ///      vertex in the edge neighbour that is located at the same
  ///      position as that first vertex).
  ///   .
  /// - The second pair contains
  ///    -# the vertex at the other end of the edge in the reference element
  ///    -# the corresponding vertex in the edge neighbour.
  ///    .
  /// .
  /// These two pairs completely define the relative rotation
  /// between the reference element and its edge neighbour. The map
  /// returns a pair which contains the up_equivalent and the
  /// right_equivalent of the edge neighbour, i.e. it tells us
  /// which direction in the edge neighbour coincides with the
  /// up (or right) direction in the reference element.
  //====================================================================
  std::map<std::pair<std::pair<int, int>, std::pair<int, int>>,
           std::pair<int, int>>
    OcTree::Up_and_right_equivalent_for_pairs_of_vertices;


  // Protected static member data:
  //------------------------------


  //====================================================================
  /// Entry in rotation matrix: cos(i*90)
  //====================================================================
  Vector<int> OcTree::Cosi;


  //====================================================================
  /// Entry in rotation matrix: sin(i*90)
  //====================================================================
  Vector<int> OcTree::Sini;

  //====================================================================
  /// Array of direction/octant adjacency scheme:
  /// Is_adjacent(direction,octant): Is face/edge \c direction
  /// adjacent to octant \c octant ? Table in Samet's book.
  //====================================================================
  DenseMatrix<bool> OcTree::Is_adjacent;

  //====================================================================
  /// Reflection scheme: Reflect(direction,octant): Get mirror
  /// of octant/edge in specified direction. E.g. Reflect(LDF,L)=RDF
  //====================================================================
  DenseMatrix<int> OcTree::Reflect;

  //====================================================================
  /// Determine common face of edges or octants.
  /// Slightly bizarre lookup scheme from Samet's book.
  //====================================================================
  DenseMatrix<int> OcTree::Common_face;

  //====================================================================
  /// Colours for neighbours in various directions
  //====================================================================
  Vector<std::string> OcTree::Colour;

  //====================================================================
  /// s_base(i,direction):  Initial value for coordinate s[i] on
  /// the face indicated by direction (L/R/U/D/F/B)
  //====================================================================
  DenseMatrix<double> OcTree::S_base;

  //====================================================================
  /// Each face of the RefineableQElement<3> that is represented
  /// by the octree is parametrised by two (of the three)
  /// local coordinates that parametrise the entire 3D element. E.g.
  /// the B[ack] face is parametrised by (s[0], s[1]); the D[own] face
  /// is parametrised by (s[0],s[2]); etc. We always identify the
  /// in-face coordinate with the lower (3D) index with the subscript
  /// \c _lo and the one with the larger (3D) index with the subscript \c _hi.
  ///  Here we set up the translation scheme between the 2D in-face
  /// coordinates (s_lo,s_hi) and the corresponding 3D coordinates:
  /// If we're located on face \c face [L/R/F/B/U/D], then
  /// an increase in s_lo from -1 to +1 corresponds to a change
  /// of \c S_steplo(i,face) in the 3D coordinate \c s[i]. S_steplo(i,direction)
  //====================================================================
  DenseMatrix<double> OcTree::S_steplo;


  //====================================================================
  /// If we're located on face \c face [L/R/F/B/U/D], then
  /// an increase in s_hi from -1 to +1 corresponds to a change
  /// of \c S_stephi(i,face) in the 3D coordinate \ s[i].
  /// [Read the discussion of \c S_steplo for an explanation of
  /// the subscripts \c _hi and \c _lo.]
  //====================================================================
  DenseMatrix<double> OcTree::S_stephi;

  //====================================================================
  ///  Relative to the left/down/back vertex in any (father) octree, the
  /// corresponding vertex in the son specified by \c son_octant has an offset.
  /// If we project the son_octant's left/down/back vertex onto the
  /// father's face \c face, it is located at the in-face coordinate
  ///  \c s_lo = h/2 \c S_directlo(face,son_octant). [See discussion of
  ///  \c S_steplo for an explanation of the subscripts \c _hi and \c _lo.]
  //====================================================================
  DenseMatrix<double> OcTree::S_directlo;

  //====================================================================
  /// Relative to the left/down/back vertex in any (father) octree, the
  /// corresponding vertex in the son specified by \c son_octant has an offset.
  /// If we project the son_octant's left/down/back vertex onto the
  /// father's face \c face, it is located at the in-face coordinate
  /// \c s_hi = h/2 \c S_directlhi(face,son_octant). [See discussion of
  /// \c S_steplo for an explanation of the subscripts \c _hi and \c _lo.]
  //====================================================================
  DenseMatrix<double> OcTree::S_directhi;

  //====================================================================
  /// S_base_edge(i,edge):  Initial value for coordinate s[i] on
  /// the specified edge (LF/RF/...).
  //====================================================================
  DenseMatrix<double> OcTree::S_base_edge;

  //====================================================================
  ///  Each edge of the RefineableQElement<3> that is represented
  /// by the octree  is parametrised by one (of the three)
  /// local coordinates that parametrise the entire 3D element.
  /// If we're located on edge \c edge [DB,UB,...], then
  /// an increase in s from -1 to +1 corresponds to a change
  /// of \c s_step_edge(i,edge) in the 3D coordinates \c s[i].
  //====================================================================
  DenseMatrix<double> OcTree::S_step_edge;

  //====================================================================
  ///  Relative to the left/down/back vertex in any (father) octree, the
  /// corresponding vertex in the son specified by \c son_octant has an offset.
  /// If we project the son_octant's left/down/back vertex onto the
  /// father's edge \c edge, it is located at the in-face coordinate
  ///  \c s_lo = h/2 \c S_direct_edge(edge,son_octant).
  //====================================================================
  DenseMatrix<double> OcTree::S_direct_edge;


  // End static data
  //----------------


  //===================================================================
  /// This function is used to translate the position of a vertex node
  /// (given by his local number n into a vector giving the position of
  /// this node in the local coordinates system.
  /// It also needs the value of nnode1d to work.
  //===================================================================
  Vector<int> OcTree::vertex_node_to_vector(const unsigned& n,
                                            const unsigned& nnode1d)
  {
#ifdef PARANOID
    if ((n != 0) && (n != nnode1d - 1) && (n != (nnode1d - 1) * nnode1d) &&
        (n != nnode1d * nnode1d - 1) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + 0) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + nnode1d - 1) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + (nnode1d - 1) * nnode1d) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + nnode1d * nnode1d - 1))
    {
      std::ostringstream error_stream;
      error_stream << "Node " << n
                   << " is not a vertex node in a brick element with "
                   << nnode1d << " nodes along each edge!";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    int a, b, c;
    Vector<int> result_vect(3);
    a = n / (nnode1d * nnode1d);
    b = (n - a * nnode1d * nnode1d) / nnode1d;
    c = n - a * nnode1d * nnode1d - b * nnode1d;

    result_vect[0] = 2 * c / (nnode1d - 1) - 1;
    result_vect[1] = 2 * b / (nnode1d - 1) - 1;
    result_vect[2] = 2 * a / (nnode1d - 1) - 1;

    return result_vect;
  }

  //==================================================================
  /// Given an edge, this function returns the faces on which it lies.
  //==================================================================
  Vector<int> OcTree::faces_of_common_edge(const int& edge)
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((edge != LB) && (edge != RB) && (edge != DB) && (edge != UB) &&
        (edge != LD) && (edge != RD) && (edge != LU) && (edge != RU) &&
        (edge != LF) && (edge != RF) && (edge != DF) && (edge != UF))
    {
      std::ostringstream error_stream;
      error_stream << "Edge " << Direct_string[edge] << "is not valid!";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Allocate space for the result
    Vector<int> faces(2, 0);

    // Which faces does this edge lie on?
    switch (edge)
    {
      case LB:
        // First entry
        faces[0] = L;

        // Second entry
        faces[1] = B;

        // Break
        break;
      case RB:
        // First entry
        faces[0] = R;

        // Second entry
        faces[1] = B;

        // Break
        break;
      case DB:
        // First entry
        faces[0] = D;

        // Second entry
        faces[1] = B;

        // Break
        break;
      case UB:
        // First entry
        faces[0] = U;

        // Second entry
        faces[1] = B;

        // Break
        break;
      case LD:
        // First entry
        faces[0] = L;

        // Second entry
        faces[1] = D;

        // Break
        break;
      case RD:
        // First entry
        faces[0] = R;

        // Second entry
        faces[1] = D;

        // Break
        break;
      case LU:
        // First entry
        faces[0] = L;

        // Second entry
        faces[1] = U;

        // Break
        break;
      case RU:
        // First entry
        faces[0] = R;

        // Second entry
        faces[1] = U;

        // Break
        break;
      case LF:
        // First entry
        faces[0] = L;

        // Second entry
        faces[1] = F;

        // Break
        break;
      case RF:
        // First entry
        faces[0] = R;

        // Second entry
        faces[1] = F;

        // Break
        break;
      case DF:
        // First entry
        faces[0] = D;

        // Second entry
        faces[1] = F;

        // Break
        break;
      case UF:
        // First entry
        faces[0] = U;

        // Second entry
        faces[1] = F;

        // Break
        break;
      default:
        // Throw an error
        throw OomphLibError("Incorrect edge input!",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }

    // Return the faces
    return faces;
  } // End of faces_of_common_edge


  //==================================================================
  /// Return the local node number of given vertex
  /// in an element with nnode1d nodes in each coordinate direction.
  //==================================================================
  unsigned OcTree::vertex_to_node_number(const int& vertex,
                                         const unsigned& nnode1d)
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((vertex != LDB) && (vertex != RDB) && (vertex != LUB) &&
        (vertex != RUB) && (vertex != LDF) && (vertex != RDF) &&
        (vertex != LUF) && (vertex != RUF))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong vertex: " << Direct_string[vertex] << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    switch (vertex)
    {
      case LDB:
        return 0;
        break;

      case RDB:
        return nnode1d - 1;
        break;


      case LUB:
        return nnode1d * (nnode1d - 1);
        break;

      case RUB:
        return nnode1d * nnode1d - 1;
        break;


      case LDF:
        return nnode1d * nnode1d * (nnode1d - 1);
        break;


      case RDF:
        return (nnode1d * nnode1d + 1) * (nnode1d - 1);
        break;

      case LUF:
        return nnode1d * nnode1d * nnode1d - nnode1d;
        break;

      case RUF:
        return nnode1d * nnode1d * nnode1d - 1;
        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here. vertex: " << Direct_string[vertex]
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        break;
    }


    std::ostringstream error_stream;
    error_stream << "Never get here. vertex: " << Direct_string[vertex]
                 << std::endl;
    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //==================================================================
  /// Return the vertex of local (vertex) node n
  /// in an element with nnode1d nodes in each coordinate direction.
  //==================================================================
  int OcTree::node_number_to_vertex(const unsigned& n, const unsigned& nnode1d)
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((n != 0) && (n != nnode1d - 1) && (n != (nnode1d - 1) * nnode1d) &&
        (n != nnode1d * nnode1d - 1) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + 0) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + nnode1d - 1) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + (nnode1d - 1) * nnode1d) &&
        (n != (nnode1d * nnode1d) * (nnode1d - 1) + nnode1d * nnode1d - 1))
    {
      std::ostringstream error_stream;
      error_stream << "Node " << n
                   << " is not a vertex node in a brick element with "
                   << nnode1d << " nodes along each edge!";

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    if (n == 0)
    {
      return LDB;
    }
    else if (n == nnode1d - 1)
    {
      return RDB;
    }
    else if (n == nnode1d * (nnode1d - 1))
    {
      return LUB;
    }
    else if (n == nnode1d * nnode1d - 1)
    {
      return RUB;
    }
    else if (n == nnode1d * nnode1d * (nnode1d - 1))
    {
      return LDF;
    }
    else if (n == (nnode1d * nnode1d + 1) * (nnode1d - 1))
    {
      return RDF;
    }
    else if (n == nnode1d * nnode1d * nnode1d - nnode1d)
    {
      return LUF;
    }
    else if (n == nnode1d * nnode1d * nnode1d - 1)
    {
      return RUF;
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "Never get here. local node number: " << n
                   << " is not a vertex node in a brick element with "
                   << nnode1d << " nodes along each edge!" << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //==================================================================
  /// This function takes as argument two node numbers of two nodes
  /// delimiting an edge, and one face of this edge and returns the
  /// other face that is sharing this edge. The node numbers given to
  /// this function MUST be vertices nodes to work.
  /// it also need the value of nnode1d to work.
  /// (\c face is a direction in the set U,D,F,B,L,R).
  //==================================================================
  int OcTree::get_the_other_face(const unsigned& n1,
                                 const unsigned& n2,
                                 const unsigned& nnode1d,
                                 const int& face)
  {
    Vector<int> vect_node1(3);
    Vector<int> vect_node2(3);
    Vector<int> vect_face(3);
    Vector<int> vect_other_face(3);

    // Translate the nodes to vectors
    vect_node1 = vertex_node_to_vector(n1, nnode1d);
    vect_node2 = vertex_node_to_vector(n2, nnode1d);

    // Translate the face to a vector
    vect_face = Direction_to_vector[face];

    // Compute the vector to the other face -- magic, courtesy of Renaud
    // Schleck!
    for (unsigned i = 0; i < 3; i++)
    {
      // Calculate the i-th entry
      vect_other_face[i] = (vect_node1[i] + vect_node2[i]) / 2 - vect_face[i];

#ifdef PARANOID
      if ((vect_other_face[i] != 1) && (vect_other_face[i] != -1) &&
          (vect_other_face[i] != 0))
      {
        throw OomphLibError("The nodes given are not vertices",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // return the corresponding face
    return Vector_to_direction[vect_other_face];
  }


  //====================================================================
  /// Build the rotation matrix for a rotation around the axis \c axis of
  /// an angle \c angle*90
  //====================================================================
  void OcTree::construct_rotation_matrix(int& axis,
                                         int& angle,
                                         DenseMatrix<int>& mat)
  {
    using namespace OcTreeNames;
    int a = 0, b = 0, c = 0, i, j;

    switch (axis)
    {
      case R:
        a = 1;
        b = 2;
        c = 0;
        break;
      case U:
        a = 2;
        b = 0;
        c = 1;
        break;
      case F:
        a = 0;
        b = 1;
        c = 2;
        break;
      default:
        // Object to create error message
        std::ostringstream error_message_stream;

        // Create the error message
        error_message_stream << "Bad axis (" << axis << "). Expected "
                             << OcTreeNames::R << ", " << OcTreeNames::U
                             << " or " << OcTreeNames::F << "." << std::endl;

        // Throw the error message
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        mat(i, j) = 0;
      }
    }
    mat(a, a) = Cosi[angle];
    mat(b, b) = Cosi[angle];
    mat(a, b) = -Sini[angle];
    mat(b, a) = Sini[angle];
    mat(c, c) = 1;
  }

  //===================================================================
  /// Helper: Performs the operation Mat3=Mat1*Mat2
  //===================================================================
  void OcTree::mult_mat_mat(const DenseMatrix<int>& mat1,
                            const DenseMatrix<int>& mat2,
                            DenseMatrix<int>& mat3)
  {
    int Sum, i, j, k;
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        Sum = 0;
        for (k = 0; k < 3; k++)
        {
          Sum += mat1(i, k) * mat2(k, j);
        }
        mat3(i, j) = Sum;
      }
    }
  }


  //===================================================================
  /// Helper: Performs the operation Vect2=Mat*Vect1
  //===================================================================
  void OcTree::mult_mat_vect(const DenseMatrix<int>& mat,
                             const Vector<int>& vect1,
                             Vector<int>& vect2)
  {
    int i, k, sum;
    for (i = 0; i < 3; i++)
    {
      sum = 0;
      for (k = 0; k < 3; k++)
      {
        sum += mat(i, k) * vect1[k];
      }
      vect2[i] = sum;
    }
  }


  //===================================================================
  /// A rotation is defined by the newUp and newRight
  /// directions; so if Up becomes newUp and Right becomes newRight
  /// then dir becomes rotate(newUp,newRight,dir);
  //===================================================================
  int OcTree::rotate(const int& new_up, const int& new_right, const int& dir)
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((new_up != L) && (new_up != R) && (new_up != F) && (new_up != B) &&
        (new_up != U) && (new_up != D))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong new_up: " << Direct_string[new_up] << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if ((new_right != L) && (new_right != R) && (new_right != F) &&
        (new_right != B) && (new_right != U) && (new_right != D))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong new_right: " << Direct_string[new_right]
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    Vector<int> vect_dir(3);
    Vector<int> vect_new_dir(3);
    // translate the direction to a vector
    vect_dir = Direction_to_vector[dir];

    // rotate it
    vect_new_dir = rotate(new_up, new_right, vect_dir);

    // translate the vector back to a direction
    return Vector_to_direction[vect_new_dir];
  }


  //==================================================================
  /// This function rotates a vector according to a rotation of
  /// the axes that changes up to new_up and right to new_right.
  //==================================================================
  Vector<int> OcTree::rotate(const int& new_up,
                             const int& new_right,
                             const Vector<int>& dir)
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((new_up != L) && (new_up != R) && (new_up != F) && (new_up != B) &&
        (new_up != U) && (new_up != D))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong new_up: " << Direct_string[new_up] << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if ((new_right != L) && (new_right != R) && (new_right != F) &&
        (new_right != B) && (new_right != U) && (new_right != D))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong new_right: " << Direct_string[new_right]
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Every possible rotation of an element can be produced by the
    // compostition of two (or one) rotations about one of the main
    // axes of the element (R, U, F). So for each (newRight,newUp)
    // we define the (angle1,axis1) of the first rotation and the
    // (angle2, axis2) of the second one (if needed)

    int axis1, axis2, angle1, angle2, nrot;
    nrot = 2;
    if (new_up == U)
    {
      nrot = 1;
      axis1 = U;
      switch (new_right)
      {
        case R:
          angle1 = 0;
          break;
        case B:
          angle1 = 1;
          break;
        case L:
          angle1 = 2;
          break;
        case F:
          angle1 = 3;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be R, B, L, or F" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    if (new_up == D)
    {
      switch (new_right)
      {
        case R:
          nrot = 1;
          axis1 = R;
          angle1 = 2;
          break;
        case B:
          axis1 = R;
          angle1 = 2;
          axis2 = U;
          angle2 = 1;
          break;
        case L:
          nrot = 1;
          axis1 = F;
          angle1 = 2;
          break;
        case F:
          axis1 = R;
          angle1 = 2;
          axis2 = U;
          angle2 = 3;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be R, B, L, or F" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    if (new_up == R)
    {
      switch (new_right)
      {
        case D:
          nrot = 1;
          axis1 = F;
          angle1 = 3;
          break;
        case B:
          axis1 = F;
          angle1 = 3;
          axis2 = R;
          angle2 = 1;
          break;
        case U:
          axis1 = F;
          angle1 = 1;
          axis2 = U;
          angle2 = 2;
          break;
        case F:
          axis1 = F;
          angle1 = 3;
          axis2 = R;
          angle2 = 3;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be D, B, U, or F" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    if (new_up == L)
    {
      switch (new_right)
      {
        case D:
          axis1 = F;
          angle1 = 1;
          axis2 = R;
          angle2 = 2;
          break;
        case B:
          axis1 = F;
          angle1 = 1;
          axis2 = R;
          angle2 = 3;
          break;
        case U:
          nrot = 1;
          axis1 = F;
          angle1 = 1;
          break;
        case F:
          axis1 = F;
          angle1 = 1;
          axis2 = R;
          angle2 = 1;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be D, B, U, or F" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    if (new_up == F)
    {
      switch (new_right)
      {
        case R:
          nrot = 1;
          axis1 = R;
          angle1 = 1;
          break;
        case L:
          axis1 = R;
          angle1 = 1;
          axis2 = F;
          angle2 = 2;
          break;
        case U:
          axis1 = R;
          angle1 = 1;
          axis2 = F;
          angle2 = 1;
          break;
        case D:
          axis1 = R;
          angle1 = 1;
          axis2 = F;
          angle2 = 3;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be R, L, U, or D" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }
    if (new_up == B)
    {
      switch (new_right)
      {
        case R:
          nrot = 1;
          axis1 = R;
          angle1 = 3;
          break;
        case L:
          axis1 = R;
          angle1 = 3;
          axis2 = F;
          angle2 = 2;
          break;
        case U:
          axis1 = R;
          angle1 = 3;
          axis2 = F;
          angle2 = 1;
          break;
        case D:
          axis1 = R;
          angle1 = 3;
          axis2 = F;
          angle2 = 3;
          break;
        default:
          std::ostringstream error_stream;
          error_stream << "New_right is " << new_right << " ("
                       << Direct_string[new_right] << "). "
                       << "It should be R, L, U, or D" << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    Vector<int> vect_new_dir(3);
    DenseMatrix<int> mat1(3);
    DenseMatrix<int> mat2(3);
    DenseMatrix<int> mat3(3);

    // Then we build the first rotation matrix
    construct_rotation_matrix(axis1, angle1, mat1);

    // If needed we build the second one
    if (nrot == 2)
    {
      // And we make the composition of the two rotations
      // which is stored in Mat3
      construct_rotation_matrix(axis2, angle2, mat2);
      mult_mat_mat(mat2, mat1, mat3);
    }
    else
    {
      // Else we just copy Mat1 into Mat3
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          mat3(i, j) = mat1(i, j);
        }
      }
    }

    // Rotate
    mult_mat_vect(mat3, dir, vect_new_dir);

    // Return the corresponding vector
    return vect_new_dir;
  }


  //====================================================================
  /// Setup static data for OcTree
  //====================================================================
  void OcTree::setup_static_data()
  {
    using namespace OcTreeNames;


#ifdef PARANOID
    if (Tree::OMEGA != OcTree::OMEGA)
    {
      std::ostringstream error_stream;
      error_stream << "Inconsistent enumeration!  \n  Tree::OMEGA="
                   << Tree::OMEGA << "\nOcTree::OMEGA=" << OcTree::OMEGA
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Private static data
    //--------------------

    // Set flag to indicate that static data has been setup
    Static_data_has_been_setup = true;

    // Initialisation of cosi and sini used in rotation matrix
    Cosi.resize(4);
    Sini.resize(4);

    Cosi[0] = 1;
    Cosi[1] = 0;
    Cosi[2] = -1;
    Cosi[3] = 0;

    Sini[0] = 0;
    Sini[1] = 1;
    Sini[2] = 0;
    Sini[3] = -1;


    // Build direction/octant adjacency scheme
    // Is_adjacent(i_face,j_octant):
    // Is face adjacent to octant?
    // See the table in Samet book
    Is_adjacent.resize(27, 27);
    Is_adjacent(L, LDB) = true;
    Is_adjacent(R, LDB) = false;
    Is_adjacent(D, LDB) = true;
    Is_adjacent(U, LDB) = false;
    Is_adjacent(B, LDB) = true;
    Is_adjacent(F, LDB) = false;

    Is_adjacent(L, LDF) = true;
    Is_adjacent(R, LDF) = false;
    Is_adjacent(D, LDF) = true;
    Is_adjacent(U, LDF) = false;
    Is_adjacent(B, LDF) = false;
    Is_adjacent(F, LDF) = true;

    Is_adjacent(L, LUB) = true;
    Is_adjacent(R, LUB) = false;
    Is_adjacent(D, LUB) = false;
    Is_adjacent(U, LUB) = true;
    Is_adjacent(B, LUB) = true;
    Is_adjacent(F, LUB) = false;

    Is_adjacent(L, LUF) = true;
    Is_adjacent(R, LUF) = false;
    Is_adjacent(D, LUF) = false;
    Is_adjacent(U, LUF) = true;
    Is_adjacent(B, LUF) = false;
    Is_adjacent(F, LUF) = true;

    Is_adjacent(L, RDB) = false;
    Is_adjacent(R, RDB) = true;
    Is_adjacent(D, RDB) = true;
    Is_adjacent(U, RDB) = false;
    Is_adjacent(B, RDB) = true;
    Is_adjacent(F, RDB) = false;

    Is_adjacent(L, RDF) = false;
    Is_adjacent(R, RDF) = true;
    Is_adjacent(D, RDF) = true;
    Is_adjacent(U, RDF) = false;
    Is_adjacent(B, RDF) = false;
    Is_adjacent(F, RDF) = true;

    Is_adjacent(L, RUB) = false;
    Is_adjacent(R, RUB) = true;
    Is_adjacent(D, RUB) = false;
    Is_adjacent(U, RUB) = true;
    Is_adjacent(B, RUB) = true;
    Is_adjacent(F, RUB) = false;

    Is_adjacent(L, RUF) = false;
    Is_adjacent(R, RUF) = true;
    Is_adjacent(D, RUF) = false;
    Is_adjacent(U, RUF) = true;
    Is_adjacent(B, RUF) = false;
    Is_adjacent(F, RUF) = true;


    // Build direction/octant adjacency scheme
    // Is_adjacent(i_edge,j_octant):
    // Is edge adjacent to octant?
    // See the table in Samet book
    Is_adjacent(LD, LDB) = true;
    Is_adjacent(LD, LDF) = true;
    Is_adjacent(LD, LUB) = false;
    Is_adjacent(LD, LUF) = false;
    Is_adjacent(LD, RDB) = false;
    Is_adjacent(LD, RDF) = false;
    Is_adjacent(LD, RUB) = false;
    Is_adjacent(LD, RUF) = false;

    Is_adjacent(LU, LDB) = false;
    Is_adjacent(LU, LDF) = false;
    Is_adjacent(LU, LUB) = true;
    Is_adjacent(LU, LUF) = true;
    Is_adjacent(LU, RDB) = false;
    Is_adjacent(LU, RDF) = false;
    Is_adjacent(LU, RUB) = false;
    Is_adjacent(LU, RUF) = false;


    Is_adjacent(LB, LDB) = true;
    Is_adjacent(LB, LDF) = false;
    Is_adjacent(LB, LUB) = true;
    Is_adjacent(LB, LUF) = false;
    Is_adjacent(LB, RDB) = false;
    Is_adjacent(LB, RDF) = false;
    Is_adjacent(LB, RUB) = false;
    Is_adjacent(LB, RUF) = false;


    Is_adjacent(LF, LDB) = false;
    Is_adjacent(LF, LDF) = true;
    Is_adjacent(LF, LUB) = false;
    Is_adjacent(LF, LUF) = true;
    Is_adjacent(LF, RDB) = false;
    Is_adjacent(LF, RDF) = false;
    Is_adjacent(LF, RUB) = false;
    Is_adjacent(LF, RUF) = false;


    Is_adjacent(RD, LDB) = false;
    Is_adjacent(RD, LDF) = false;
    Is_adjacent(RD, LUB) = false;
    Is_adjacent(RD, LUF) = false;
    Is_adjacent(RD, RDB) = true;
    Is_adjacent(RD, RDF) = true;
    Is_adjacent(RD, RUB) = false;
    Is_adjacent(RD, RUF) = false;


    Is_adjacent(RU, LDB) = false;
    Is_adjacent(RU, LDF) = false;
    Is_adjacent(RU, LUB) = false;
    Is_adjacent(RU, LUF) = false;
    Is_adjacent(RU, RDB) = false;
    Is_adjacent(RU, RDF) = false;
    Is_adjacent(RU, RUB) = true;
    Is_adjacent(RU, RUF) = true;

    Is_adjacent(RB, LDB) = false;
    Is_adjacent(RB, LDF) = false;
    Is_adjacent(RB, LUB) = false;
    Is_adjacent(RB, LUF) = false;
    Is_adjacent(RB, RDB) = true;
    Is_adjacent(RB, RDF) = false;
    Is_adjacent(RB, RUB) = true;
    Is_adjacent(RB, RUF) = false;

    Is_adjacent(RF, LDB) = false;
    Is_adjacent(RF, LDF) = false;
    Is_adjacent(RF, LUB) = false;
    Is_adjacent(RF, LUF) = false;
    Is_adjacent(RF, RDB) = false;
    Is_adjacent(RF, RDF) = true;
    Is_adjacent(RF, RUB) = false;
    Is_adjacent(RF, RUF) = true;


    Is_adjacent(DB, LDB) = true;
    Is_adjacent(DB, LDF) = false;
    Is_adjacent(DB, LUB) = false;
    Is_adjacent(DB, LUF) = false;
    Is_adjacent(DB, RDB) = true;
    Is_adjacent(DB, RDF) = false;
    Is_adjacent(DB, RUB) = false;
    Is_adjacent(DB, RUF) = false;


    Is_adjacent(DF, LDB) = false;
    Is_adjacent(DF, LDF) = true;
    Is_adjacent(DF, LUB) = false;
    Is_adjacent(DF, LUF) = false;
    Is_adjacent(DF, RDB) = false;
    Is_adjacent(DF, RDF) = true;
    Is_adjacent(DF, RUB) = false;
    Is_adjacent(DF, RUF) = false;

    Is_adjacent(UB, LDB) = false;
    Is_adjacent(UB, LDF) = false;
    Is_adjacent(UB, LUB) = true;
    Is_adjacent(UB, LUF) = false;
    Is_adjacent(UB, RDB) = false;
    Is_adjacent(UB, RDF) = false;
    Is_adjacent(UB, RUB) = true;
    Is_adjacent(UB, RUF) = false;

    Is_adjacent(UF, LDB) = false;
    Is_adjacent(UF, LDF) = false;
    Is_adjacent(UF, LUB) = false;
    Is_adjacent(UF, LUF) = true;
    Is_adjacent(UF, RDB) = false;
    Is_adjacent(UF, RDF) = false;
    Is_adjacent(UF, RUB) = false;
    Is_adjacent(UF, RUF) = true;


    // Common face of various octants from Samet's book
    Common_face.resize(27, 27);
    Common_face(LDB, LDB) = OMEGA;
    Common_face(LDB, LDF) = OMEGA;
    Common_face(LDB, LUB) = OMEGA;
    Common_face(LDB, LUF) = L;
    Common_face(LDB, RDB) = OMEGA;
    Common_face(LDB, RDF) = D;
    Common_face(LDB, RUB) = B;
    Common_face(LDB, RUF) = OMEGA;

    Common_face(LDF, LDB) = OMEGA;
    Common_face(LDF, LDF) = OMEGA;
    Common_face(LDF, LUB) = L;
    Common_face(LDF, LUF) = OMEGA;
    Common_face(LDF, RDB) = D;
    Common_face(LDF, RDF) = OMEGA;
    Common_face(LDF, RUB) = OMEGA;
    Common_face(LDF, RUF) = F;

    Common_face(LUB, LDB) = OMEGA;
    Common_face(LUB, LDF) = L;
    Common_face(LUB, LUB) = OMEGA;
    Common_face(LUB, LUF) = OMEGA;
    Common_face(LUB, RDB) = B;
    Common_face(LUB, RDF) = OMEGA;
    Common_face(LUB, RUB) = OMEGA;
    Common_face(LUB, RUF) = U;

    Common_face(LUF, LDB) = L;
    Common_face(LUF, LDF) = OMEGA;
    Common_face(LUF, LUB) = OMEGA;
    Common_face(LUF, LUF) = OMEGA;
    Common_face(LUF, RDB) = OMEGA;
    Common_face(LUF, RDF) = F;
    Common_face(LUF, RUB) = U;
    Common_face(LUF, RUF) = OMEGA;

    Common_face(RDB, LDB) = OMEGA;
    Common_face(RDB, LDF) = D;
    Common_face(RDB, LUB) = B;
    Common_face(RDB, LUF) = OMEGA;
    Common_face(RDB, RDB) = OMEGA;
    Common_face(RDB, RDF) = OMEGA;
    Common_face(RDB, RUB) = OMEGA;
    Common_face(RDB, RUF) = R;

    Common_face(RDF, LDB) = D;
    Common_face(RDF, LDF) = OMEGA;
    Common_face(RDF, LUB) = OMEGA;
    Common_face(RDF, LUF) = F;
    Common_face(RDF, RDB) = OMEGA;
    Common_face(RDF, RDF) = OMEGA;
    Common_face(RDF, RUB) = R;
    Common_face(RDF, RUF) = OMEGA;

    Common_face(RUB, LDB) = B;
    Common_face(RUB, LDF) = OMEGA;
    Common_face(RUB, LUB) = OMEGA;
    Common_face(RUB, LUF) = U;
    Common_face(RUB, RDB) = OMEGA;
    Common_face(RUB, RDF) = R;
    Common_face(RUB, RUB) = OMEGA;
    Common_face(RUB, RUF) = OMEGA;

    Common_face(RUF, LDB) = OMEGA;
    Common_face(RUF, LDF) = F;
    Common_face(RUF, LUB) = U;
    Common_face(RUF, LUF) = OMEGA;
    Common_face(RUF, RDB) = R;
    Common_face(RUF, RDF) = OMEGA;
    Common_face(RUF, RUB) = OMEGA;
    Common_face(RUF, RUF) = OMEGA;


    // Common face of various edges/octants from Samet's book
    Common_face(LD, LDB) = OMEGA;
    Common_face(LD, LDF) = OMEGA;
    Common_face(LD, LUB) = L;
    Common_face(LD, LUF) = L;
    Common_face(LD, RDB) = D;
    Common_face(LD, RDF) = D;
    Common_face(LD, RUB) = OMEGA;
    Common_face(LD, RUF) = OMEGA;

    Common_face(LU, LDB) = L;
    Common_face(LU, LDF) = L;
    Common_face(LU, LUB) = OMEGA;
    Common_face(LU, LUF) = OMEGA;
    Common_face(LU, RDB) = OMEGA;
    Common_face(LU, RDF) = OMEGA;
    Common_face(LU, RUB) = U;
    Common_face(LU, RUF) = U;

    Common_face(LB, LDB) = OMEGA;
    Common_face(LB, LDF) = L;
    Common_face(LB, LUB) = OMEGA;
    Common_face(LB, LUF) = L;
    Common_face(LB, RDB) = B;
    Common_face(LB, RDF) = OMEGA;
    Common_face(LB, RUB) = B;
    Common_face(LB, RUF) = OMEGA;

    Common_face(LF, LDB) = L;
    Common_face(LF, LDF) = OMEGA;
    Common_face(LF, LUB) = L;
    Common_face(LF, LUF) = OMEGA;
    Common_face(LF, RDB) = OMEGA;
    Common_face(LF, RDF) = F;
    Common_face(LF, RUB) = OMEGA;
    Common_face(LF, RUF) = F;

    Common_face(RD, LDB) = D;
    Common_face(RD, LDF) = D;
    Common_face(RD, LUB) = OMEGA;
    Common_face(RD, LUF) = OMEGA;
    Common_face(RD, RDB) = OMEGA;
    Common_face(RD, RDF) = OMEGA;
    Common_face(RD, RUB) = R;
    Common_face(RD, RUF) = R;

    Common_face(RU, LDB) = OMEGA;
    Common_face(RU, LDF) = OMEGA;
    Common_face(RU, LUB) = U;
    Common_face(RU, LUF) = U;
    Common_face(RU, RDB) = R;
    Common_face(RU, RDF) = R;
    Common_face(RU, RUB) = OMEGA;
    Common_face(RU, RUF) = OMEGA;

    Common_face(RB, LDB) = B;
    Common_face(RB, LDF) = OMEGA;
    Common_face(RB, LUB) = B;
    Common_face(RB, LUF) = OMEGA;
    Common_face(RB, RDB) = OMEGA;
    Common_face(RB, RDF) = R;
    Common_face(RB, RUB) = OMEGA;
    Common_face(RB, RUF) = R;

    Common_face(RF, LDB) = OMEGA;
    Common_face(RF, LDF) = F;
    Common_face(RF, LUB) = OMEGA;
    Common_face(RF, LUF) = F;
    Common_face(RF, RDB) = R;
    Common_face(RF, RDF) = OMEGA;
    Common_face(RF, RUB) = R;
    Common_face(RF, RUF) = OMEGA;


    Common_face(DB, LDB) = OMEGA;
    Common_face(DB, LDF) = D;
    Common_face(DB, LUB) = B;
    Common_face(DB, LUF) = OMEGA;
    Common_face(DB, RDB) = OMEGA;
    Common_face(DB, RDF) = D;
    Common_face(DB, RUB) = B;
    Common_face(DB, RUF) = OMEGA;

    Common_face(DF, LDB) = D;
    Common_face(DF, LDF) = OMEGA;
    Common_face(DF, LUB) = OMEGA;
    Common_face(DF, LUF) = F;
    Common_face(DF, RDB) = D;
    Common_face(DF, RDF) = OMEGA;
    Common_face(DF, RUB) = OMEGA;
    Common_face(DF, RUF) = F;

    Common_face(UB, LDB) = B;
    Common_face(UB, LDF) = OMEGA;
    Common_face(UB, LUB) = OMEGA;
    Common_face(UB, LUF) = U;
    Common_face(UB, RDB) = B;
    Common_face(UB, RDF) = OMEGA;
    Common_face(UB, RUB) = OMEGA;
    Common_face(UB, RUF) = U;

    Common_face(UF, LDB) = OMEGA;
    Common_face(UF, LDF) = F;
    Common_face(UF, LUB) = U;
    Common_face(UF, LUF) = OMEGA;
    Common_face(UF, RDB) = OMEGA;
    Common_face(UF, RDF) = F;
    Common_face(UF, RUB) = U;
    Common_face(UF, RUF) = OMEGA;


    // Tecplot colours for neighbours in various directions
    Colour.resize(27);
    Colour[R] = "CYAN";
    Colour[L] = "RED";
    Colour[U] = "GREEN";
    Colour[D] = "BLUE";
    Colour[F] = "CUSTOM3";
    Colour[B] = "PURPLE";
    Colour[OMEGA] = "YELLOW";

    Colour[LB] = "BLUE";
    Colour[RB] = "BLUE";
    Colour[DB] = "BLUE";
    Colour[UB] = "BLUE";

    Colour[LD] = "GREEN";
    Colour[RD] = "GREEN";
    Colour[LU] = "GREEN";
    Colour[RU] = "GREEN";


    Colour[LF] = "RED";
    Colour[RF] = "RED";
    Colour[DF] = "RED";
    Colour[UF] = "RED";

    Colour[OMEGA] = "YELLOW";


    // Reflection scheme:
    // Reflect(direction,octant): Get mirror of octant in direction
    // Table in Samet book as well
    Reflect.resize(27, 27);

    Reflect(L, LDB) = RDB;
    Reflect(R, LDB) = RDB;
    Reflect(D, LDB) = LUB;
    Reflect(U, LDB) = LUB;
    Reflect(B, LDB) = LDF;
    Reflect(F, LDB) = LDF;

    Reflect(L, LDF) = RDF;
    Reflect(R, LDF) = RDF;
    Reflect(D, LDF) = LUF;
    Reflect(U, LDF) = LUF;
    Reflect(B, LDF) = LDB;
    Reflect(F, LDF) = LDB;

    Reflect(L, LUB) = RUB;
    Reflect(R, LUB) = RUB;
    Reflect(D, LUB) = LDB;
    Reflect(U, LUB) = LDB;
    Reflect(B, LUB) = LUF;
    Reflect(F, LUB) = LUF;

    Reflect(L, LUF) = RUF;
    Reflect(R, LUF) = RUF;
    Reflect(D, LUF) = LDF;
    Reflect(U, LUF) = LDF;
    Reflect(B, LUF) = LUB;
    Reflect(F, LUF) = LUB;

    Reflect(L, RDB) = LDB;
    Reflect(R, RDB) = LDB;
    Reflect(D, RDB) = RUB;
    Reflect(U, RDB) = RUB;
    Reflect(B, RDB) = RDF;
    Reflect(F, RDB) = RDF;

    Reflect(L, RDF) = LDF;
    Reflect(R, RDF) = LDF;
    Reflect(D, RDF) = RUF;
    Reflect(U, RDF) = RUF;
    Reflect(B, RDF) = RDB;
    Reflect(F, RDF) = RDB;

    Reflect(L, RUB) = LUB;
    Reflect(R, RUB) = LUB;
    Reflect(D, RUB) = RDB;
    Reflect(U, RUB) = RDB;
    Reflect(B, RUB) = RUF;
    Reflect(F, RUB) = RUF;

    Reflect(L, RUF) = LUF;
    Reflect(R, RUF) = LUF;
    Reflect(D, RUF) = RDF;
    Reflect(U, RUF) = RDF;
    Reflect(B, RUF) = RUB;
    Reflect(F, RUF) = RUB;


    // Reflection scheme:
    // Reflect(direction (edge) ,octant): Get mirror of edge in direction
    // Table in Samet book as well

    Reflect(LD, LDB) = RUB;
    Reflect(LU, LDB) = RUB;
    Reflect(RD, LDB) = RUB;
    Reflect(RU, LDB) = RUB;

    Reflect(LD, LDF) = RUF;
    Reflect(LU, LDF) = RUF;
    Reflect(RD, LDF) = RUF;
    Reflect(RU, LDF) = RUF;

    Reflect(LD, LUB) = RDB;
    Reflect(LU, LUB) = RDB;
    Reflect(RD, LUB) = RDB;
    Reflect(RU, LUB) = RDB;

    Reflect(LD, LUF) = RDF;
    Reflect(LU, LUF) = RDF;
    Reflect(RD, LUF) = RDF;
    Reflect(RU, LUF) = RDF;


    Reflect(LD, RDB) = LUB;
    Reflect(LU, RDB) = LUB;
    Reflect(RD, RDB) = LUB;
    Reflect(RU, RDB) = LUB;

    Reflect(LD, RDF) = LUF;
    Reflect(LU, RDF) = LUF;
    Reflect(RD, RDF) = LUF;
    Reflect(RU, RDF) = LUF;

    Reflect(LD, RUB) = LDB;
    Reflect(LU, RUB) = LDB;
    Reflect(RD, RUB) = LDB;
    Reflect(RU, RUB) = LDB;

    Reflect(LD, RUF) = LDF;
    Reflect(LU, RUF) = LDF;
    Reflect(RD, RUF) = LDF;
    Reflect(RU, RUF) = LDF;


    Reflect(LB, LDB) = RDF;
    Reflect(LF, LDB) = RDF;
    Reflect(RB, LDB) = RDF;
    Reflect(RF, LDB) = RDF;

    Reflect(LB, LDF) = RDB;
    Reflect(LF, LDF) = RDB;
    Reflect(RB, LDF) = RDB;
    Reflect(RF, LDF) = RDB;

    Reflect(LB, LUB) = RUF;
    Reflect(LF, LUB) = RUF;
    Reflect(RB, LUB) = RUF;
    Reflect(RF, LUB) = RUF;

    Reflect(LB, LUF) = RUB;
    Reflect(LF, LUF) = RUB;
    Reflect(RB, LUF) = RUB;
    Reflect(RF, LUF) = RUB;

    Reflect(LB, RDB) = LDF;
    Reflect(LF, RDB) = LDF;
    Reflect(RB, RDB) = LDF;
    Reflect(RF, RDB) = LDF;

    Reflect(LB, RDF) = LDB;
    Reflect(LF, RDF) = LDB;
    Reflect(RB, RDF) = LDB;
    Reflect(RF, RDF) = LDB;

    Reflect(LB, RUB) = LUF;
    Reflect(LF, RUB) = LUF;
    Reflect(RB, RUB) = LUF;
    Reflect(RF, RUB) = LUF;

    Reflect(LB, RUF) = LUB;
    Reflect(LF, RUF) = LUB;
    Reflect(RB, RUF) = LUB;
    Reflect(RF, RUF) = LUB;


    Reflect(DB, LDB) = LUF;
    Reflect(DF, LDB) = LUF;
    Reflect(UB, LDB) = LUF;
    Reflect(UF, LDB) = LUF;

    Reflect(DB, LDF) = LUB;
    Reflect(DF, LDF) = LUB;
    Reflect(UB, LDF) = LUB;
    Reflect(UF, LDF) = LUB;

    Reflect(DB, LUB) = LDF;
    Reflect(DF, LUB) = LDF;
    Reflect(UB, LUB) = LDF;
    Reflect(UF, LUB) = LDF;

    Reflect(DB, LUF) = LDB;
    Reflect(DF, LUF) = LDB;
    Reflect(UB, LUF) = LDB;
    Reflect(UF, LUF) = LDB;

    Reflect(DB, RDB) = RUF;
    Reflect(DF, RDB) = RUF;
    Reflect(UB, RDB) = RUF;
    Reflect(UF, RDB) = RUF;

    Reflect(DB, RDF) = RUB;
    Reflect(DF, RDF) = RUB;
    Reflect(UB, RDF) = RUB;
    Reflect(UF, RDF) = RUB;

    Reflect(DB, RUB) = RDF;
    Reflect(DF, RUB) = RDF;
    Reflect(UB, RUB) = RDF;
    Reflect(UF, RUB) = RDF;

    Reflect(DB, RUF) = RDB;
    Reflect(DF, RUF) = RDB;
    Reflect(UB, RUF) = RDB;
    Reflect(UF, RUF) = RDB;


    // S_base(i,direction) etc.:  Initial value/increment for coordinate s[i] on
    // the face indicated by direction (L/R/U/D/F/B)
    S_base.resize(3, 27);
    S_steplo.resize(3, 27);
    S_stephi.resize(3, 27);

    S_base(0, L) = -1.0;
    S_base(1, L) = -1.0;
    S_base(2, L) = -1.0;
    S_steplo(0, L) = 0.0;
    S_steplo(1, L) = 2.0;
    S_steplo(2, L) = 0.0;
    S_stephi(0, L) = 0.0;
    S_stephi(1, L) = 0.0;
    S_stephi(2, L) = 2.0;

    S_base(0, R) = 1.0;
    S_base(1, R) = -1.0;
    S_base(2, R) = -1.0;
    S_steplo(0, R) = 0.0;
    S_steplo(1, R) = 2.0;
    S_steplo(2, R) = 0.0;
    S_stephi(0, R) = 0.0;
    S_stephi(1, R) = 0.0;
    S_stephi(2, R) = 2.0;

    S_base(0, U) = -1.0;
    S_base(1, U) = 1.0;
    S_base(2, U) = -1.0;
    S_steplo(0, U) = 2.0;
    S_steplo(1, U) = 0.0;
    S_steplo(2, U) = 0.0;
    S_stephi(0, U) = 0.0;
    S_stephi(1, U) = 0.0;
    S_stephi(2, U) = 2.0;

    S_base(0, D) = -1.0;
    S_base(1, D) = -1.0;
    S_base(2, D) = -1.0;
    S_steplo(0, D) = 2.0;
    S_steplo(1, D) = 0.0;
    S_steplo(2, D) = 0.0;
    S_stephi(0, D) = 0.0;
    S_stephi(1, D) = 0.0;
    S_stephi(2, D) = 2.0;

    S_base(0, F) = -1.0;
    S_base(1, F) = -1.0;
    S_base(2, F) = 1.0;
    S_steplo(0, F) = 2.0;
    S_steplo(1, F) = 0.0;
    S_steplo(2, F) = 0.0;
    S_stephi(0, F) = 0.0;
    S_stephi(1, F) = 2.0;
    S_stephi(2, F) = 0.0;

    S_base(0, B) = -1.0;
    S_base(1, B) = -1.0;
    S_base(2, B) = -1.0;
    S_steplo(0, B) = 2.0;
    S_steplo(1, B) = 0.0;
    S_steplo(2, B) = 0.0;
    S_stephi(0, B) = 0.0;
    S_stephi(1, B) = 2.0;
    S_stephi(2, B) = 0.0;


    // Relative to the left/down/back vertex in any (father) octree, the
    // corresponding vertex in the son specified by \c son_octant has an offset.
    // If we project the son_octant's left/down/back vertex onto the
    // father's face \c face, it is located at the in-face coordinate
    //  \c s_lo = h/2 \c S_directlo(face,son_octant). [See discussion of
    //  \c S_steplo for an explanation of the subscripts \c _hi and \c _lo.]
    S_directlo.resize(27, 27);
    S_directhi.resize(27, 27);

    S_directlo(L, LDB) = 0.0;
    S_directlo(R, LDB) = 0.0;
    S_directlo(U, LDB) = 0.0;
    S_directlo(D, LDB) = 0.0;
    S_directlo(F, LDB) = 0.0;
    S_directlo(B, LDB) = 0.0;

    S_directlo(L, LDF) = 0.0;
    S_directlo(R, LDF) = 0.0;
    S_directlo(U, LDF) = 0.0;
    S_directlo(D, LDF) = 0.0;
    S_directlo(F, LDF) = 0.0;
    S_directlo(B, LDF) = 0.0;

    S_directlo(L, LUB) = 1.0;
    S_directlo(R, LUB) = 1.0;
    S_directlo(U, LUB) = 0.0;
    S_directlo(D, LUB) = 0.0;
    S_directlo(F, LUB) = 0.0;
    S_directlo(B, LUB) = 0.0;

    S_directlo(L, LUF) = 1.0;
    S_directlo(R, LUF) = 1.0;
    S_directlo(U, LUF) = 0.0;
    S_directlo(D, LUF) = 0.0;
    S_directlo(F, LUF) = 0.0;
    S_directlo(B, LUF) = 0.0;

    S_directlo(L, RDB) = 0.0;
    S_directlo(R, RDB) = 0.0;
    S_directlo(U, RDB) = 1.0;
    S_directlo(D, RDB) = 1.0;
    S_directlo(F, RDB) = 1.0;
    S_directlo(B, RDB) = 1.0;

    S_directlo(L, RDF) = 0.0;
    S_directlo(R, RDF) = 0.0;
    S_directlo(U, RDF) = 1.0;
    S_directlo(D, RDF) = 1.0;
    S_directlo(F, RDF) = 1.0;
    S_directlo(B, RDF) = 1.0;

    S_directlo(L, RUB) = 1.0;
    S_directlo(R, RUB) = 1.0;
    S_directlo(U, RUB) = 1.0;
    S_directlo(D, RUB) = 1.0;
    S_directlo(F, RUB) = 1.0;
    S_directlo(B, RUB) = 1.0;

    S_directlo(L, RUF) = 1.0;
    S_directlo(R, RUF) = 1.0;
    S_directlo(U, RUF) = 1.0;
    S_directlo(D, RUF) = 1.0;
    S_directlo(F, RUF) = 1.0;
    S_directlo(B, RUF) = 1.0;


    S_directhi(L, LDB) = 0.0;
    S_directhi(R, LDB) = 0.0;
    S_directhi(U, LDB) = 0.0;
    S_directhi(D, LDB) = 0.0;
    S_directhi(F, LDB) = 0.0;
    S_directhi(B, LDB) = 0.0;

    S_directhi(L, LDF) = 1.0;
    S_directhi(R, LDF) = 1.0;
    S_directhi(U, LDF) = 1.0;
    S_directhi(D, LDF) = 1.0;
    S_directhi(F, LDF) = 0.0;
    S_directhi(B, LDF) = 0.0;

    S_directhi(L, LUB) = 0.0;
    S_directhi(R, LUB) = 0.0;
    S_directhi(U, LUB) = 0.0;
    S_directhi(D, LUB) = 0.0;
    S_directhi(F, LUB) = 1.0;
    S_directhi(B, LUB) = 1.0;

    S_directhi(L, LUF) = 1.0;
    S_directhi(R, LUF) = 1.0;
    S_directhi(U, LUF) = 1.0;
    S_directhi(D, LUF) = 1.0;
    S_directhi(F, LUF) = 1.0;
    S_directhi(B, LUF) = 1.0;

    S_directhi(L, RDB) = 0.0;
    S_directhi(R, RDB) = 0.0;
    S_directhi(U, RDB) = 0.0;
    S_directhi(D, RDB) = 0.0;
    S_directhi(F, RDB) = 0.0;
    S_directhi(B, RDB) = 0.0;

    S_directhi(L, RDF) = 1.0;
    S_directhi(R, RDF) = 1.0;
    S_directhi(U, RDF) = 1.0;
    S_directhi(D, RDF) = 1.0;
    S_directhi(F, RDF) = 0.0;
    S_directhi(B, RDF) = 0.0;

    S_directhi(L, RUB) = 0.0;
    S_directhi(R, RUB) = 0.0;
    S_directhi(U, RUB) = 0.0;
    S_directhi(D, RUB) = 0.0;
    S_directhi(F, RUB) = 1.0;
    S_directhi(B, RUB) = 1.0;

    S_directhi(L, RUF) = 1.0;
    S_directhi(R, RUF) = 1.0;
    S_directhi(U, RUF) = 1.0;
    S_directhi(D, RUF) = 1.0;
    S_directhi(F, RUF) = 1.0;
    S_directhi(B, RUF) = 1.0;


    // S_base_edge(i,direction):  Initial value/increment for coordinate s[i] on
    // the edge indicated by direction (LB,RB,...)
    S_base_edge.resize(3, 27);
    S_step_edge.resize(3, 27);

    S_base_edge(0, LB) = -1.0;
    S_base_edge(1, LB) = -1.0;
    S_base_edge(2, LB) = -1.0;
    S_step_edge(0, LB) = 0.0;
    S_step_edge(1, LB) = 2.0;
    S_step_edge(2, LB) = 0.0;

    S_base_edge(0, RB) = 1.0;
    S_base_edge(1, RB) = -1.0;
    S_base_edge(2, RB) = -1.0;
    S_step_edge(0, RB) = 0.0;
    S_step_edge(1, RB) = 2.0;
    S_step_edge(2, RB) = 0.0;

    S_base_edge(0, DB) = -1.0;
    S_base_edge(1, DB) = -1.0;
    S_base_edge(2, DB) = -1.0;
    S_step_edge(0, DB) = 2.0;
    S_step_edge(1, DB) = 0.0;
    S_step_edge(2, DB) = 0.0;

    S_base_edge(0, UB) = -1.0;
    S_base_edge(1, UB) = 1.0;
    S_base_edge(2, UB) = -1.0;
    S_step_edge(0, UB) = 2.0;
    S_step_edge(1, UB) = 0.0;
    S_step_edge(2, UB) = 0.0;

    S_base_edge(0, LD) = -1.0;
    S_base_edge(1, LD) = -1.0;
    S_base_edge(2, LD) = -1.0;
    S_step_edge(0, LD) = 0.0;
    S_step_edge(1, LD) = 0.0;
    S_step_edge(2, LD) = 2.0;

    S_base_edge(0, RD) = 1.0;
    S_base_edge(1, RD) = -1.0;
    S_base_edge(2, RD) = -1.0;
    S_step_edge(0, RD) = 0.0;
    S_step_edge(1, RD) = 0.0;
    S_step_edge(2, RD) = 2.0;

    S_base_edge(0, LU) = -1.0;
    S_base_edge(1, LU) = 1.0;
    S_base_edge(2, LU) = -1.0;
    S_step_edge(0, LU) = 0.0;
    S_step_edge(1, LU) = 0.0;
    S_step_edge(2, LU) = 2.0;

    S_base_edge(0, RU) = 1.0;
    S_base_edge(1, RU) = 1.0;
    S_base_edge(2, RU) = -1.0;
    S_step_edge(0, RU) = 0.0;
    S_step_edge(1, RU) = 0.0;
    S_step_edge(2, RU) = 2.0;


    S_base_edge(0, LF) = -1.0;
    S_base_edge(1, LF) = -1.0;
    S_base_edge(2, LF) = 1.0;
    S_step_edge(0, LF) = 0.0;
    S_step_edge(1, LF) = 2.0;
    S_step_edge(2, LF) = 0.0;

    S_base_edge(0, RF) = 1.0;
    S_base_edge(1, RF) = -1.0;
    S_base_edge(2, RF) = 1.0;
    S_step_edge(0, RF) = 0.0;
    S_step_edge(1, RF) = 2.0;
    S_step_edge(2, RF) = 0.0;

    S_base_edge(0, DF) = -1.0;
    S_base_edge(1, DF) = -1.0;
    S_base_edge(2, DF) = 1.0;
    S_step_edge(0, DF) = 2.0;
    S_step_edge(1, DF) = 0.0;
    S_step_edge(2, DF) = 0.0;

    S_base_edge(0, UF) = -1.0;
    S_base_edge(1, UF) = 1.0;
    S_base_edge(2, UF) = 1.0;
    S_step_edge(0, UF) = 2.0;
    S_step_edge(1, UF) = 0.0;
    S_step_edge(2, UF) = 0.0;


    // Relative to the left/down/back vertex in any (father) octree, the
    // corresponding vertex in the son specified by \c son_octant has an offset.
    // If we project the son_octant's left/down/back vertex onto the
    // father's edge \c edge, it is located at the in-face coordinate
    //  \c s_lo = h/2 \c S_direct_edge(edge,son_octant).
    S_direct_edge.resize(27, 27);
    S_direct_edge(LB, LDB) = 0.0;
    S_direct_edge(RB, LDB) = 0.0;
    S_direct_edge(DB, LDB) = 0.0;
    S_direct_edge(UB, LDB) = 0.0;
    S_direct_edge(LD, LDB) = 0.0;
    S_direct_edge(RD, LDB) = 0.0;
    S_direct_edge(LU, LDB) = 0.0;
    S_direct_edge(RU, LDB) = 0.0;
    S_direct_edge(LF, LDB) = 0.0;
    S_direct_edge(RF, LDB) = 0.0;
    S_direct_edge(DF, LDB) = 0.0;
    S_direct_edge(UF, LDB) = 0.0;

    S_direct_edge(LB, RDB) = 0.0;
    S_direct_edge(RB, RDB) = 0.0;
    S_direct_edge(DB, RDB) = 1.0;
    S_direct_edge(UB, RDB) = 1.0;
    S_direct_edge(LD, RDB) = 0.0;
    S_direct_edge(RD, RDB) = 0.0;
    S_direct_edge(LU, RDB) = 0.0;
    S_direct_edge(RU, RDB) = 0.0;
    S_direct_edge(LF, RDB) = 0.0;
    S_direct_edge(RF, RDB) = 0.0;
    S_direct_edge(DF, RDB) = 1.0;
    S_direct_edge(UF, RDB) = 1.0;

    S_direct_edge(LB, LUB) = 1.0;
    S_direct_edge(RB, LUB) = 1.0;
    S_direct_edge(DB, LUB) = 0.0;
    S_direct_edge(UB, LUB) = 0.0;
    S_direct_edge(LD, LUB) = 0.0;
    S_direct_edge(RD, LUB) = 0.0;
    S_direct_edge(LU, LUB) = 0.0;
    S_direct_edge(RU, LUB) = 0.0;
    S_direct_edge(LF, LUB) = 1.0;
    S_direct_edge(RF, LUB) = 1.0;
    S_direct_edge(DF, LUB) = 0.0;
    S_direct_edge(UF, LUB) = 0.0;

    S_direct_edge(LB, RUB) = 1.0;
    S_direct_edge(RB, RUB) = 1.0;
    S_direct_edge(DB, RUB) = 1.0;
    S_direct_edge(UB, RUB) = 1.0;
    S_direct_edge(LD, RUB) = 0.0;
    S_direct_edge(RD, RUB) = 0.0;
    S_direct_edge(LU, RUB) = 0.0;
    S_direct_edge(RU, RUB) = 0.0;
    S_direct_edge(LF, RUB) = 1.0;
    S_direct_edge(RF, RUB) = 1.0;
    S_direct_edge(DF, RUB) = 1.0;
    S_direct_edge(UF, RUB) = 1.0;


    S_direct_edge(LB, LDF) = 0.0;
    S_direct_edge(RB, LDF) = 0.0;
    S_direct_edge(DB, LDF) = 0.0;
    S_direct_edge(UB, LDF) = 0.0;
    S_direct_edge(LD, LDF) = 1.0;
    S_direct_edge(RD, LDF) = 1.0;
    S_direct_edge(LU, LDF) = 1.0;
    S_direct_edge(RU, LDF) = 1.0;
    S_direct_edge(LF, LDF) = 0.0;
    S_direct_edge(RF, LDF) = 0.0;
    S_direct_edge(DF, LDF) = 0.0;
    S_direct_edge(UF, LDF) = 0.0;

    S_direct_edge(LB, RDF) = 0.0;
    S_direct_edge(RB, RDF) = 0.0;
    S_direct_edge(DB, RDF) = 1.0;
    S_direct_edge(UB, RDF) = 1.0;
    S_direct_edge(LD, RDF) = 1.0;
    S_direct_edge(RD, RDF) = 1.0;
    S_direct_edge(LU, RDF) = 1.0;
    S_direct_edge(RU, RDF) = 1.0;
    S_direct_edge(LF, RDF) = 0.0;
    S_direct_edge(RF, RDF) = 0.0;
    S_direct_edge(DF, RDF) = 1.0;
    S_direct_edge(UF, RDF) = 1.0;

    S_direct_edge(LB, LUF) = 1.0;
    S_direct_edge(RB, LUF) = 1.0;
    S_direct_edge(DB, LUF) = 0.0;
    S_direct_edge(UB, LUF) = 0.0;
    S_direct_edge(LD, LUF) = 1.0;
    S_direct_edge(RD, LUF) = 1.0;
    S_direct_edge(LU, LUF) = 1.0;
    S_direct_edge(RU, LUF) = 1.0;
    S_direct_edge(LF, LUF) = 1.0;
    S_direct_edge(RF, LUF) = 1.0;
    S_direct_edge(DF, LUF) = 0.0;
    S_direct_edge(UF, LUF) = 0.0;

    S_direct_edge(LB, RUF) = 1.0;
    S_direct_edge(RB, RUF) = 1.0;
    S_direct_edge(DB, RUF) = 1.0;
    S_direct_edge(UB, RUF) = 1.0;
    S_direct_edge(LD, RUF) = 1.0;
    S_direct_edge(RD, RUF) = 1.0;
    S_direct_edge(LU, RUF) = 1.0;
    S_direct_edge(RU, RUF) = 1.0;
    S_direct_edge(LF, RUF) = 1.0;
    S_direct_edge(RF, RUF) = 1.0;
    S_direct_edge(DF, RUF) = 1.0;
    S_direct_edge(UF, RUF) = 1.0;


    // Public static data:
    //-------------------

    // Translate (enumerated) directions into strings
    Direct_string.resize(27);
    Direct_string[LDB] = "LDB";
    Direct_string[LDF] = "LDF";
    Direct_string[LUB] = "LUB";
    Direct_string[LUF] = "LUF";
    Direct_string[RDB] = "RDB";
    Direct_string[RDF] = "RDF";
    Direct_string[RUB] = "RUB";
    Direct_string[RUF] = "RUF";


    Direct_string[L] = "L";
    Direct_string[R] = "R";
    Direct_string[U] = "U";
    Direct_string[D] = "D";
    Direct_string[F] = "F";
    Direct_string[B] = "B";

    Direct_string[LU] = "LU";
    Direct_string[LD] = "LD";
    Direct_string[LF] = "LF";
    Direct_string[LB] = "LB";
    Direct_string[RU] = "RU";
    Direct_string[RD] = "RD";
    Direct_string[RF] = "RF";
    Direct_string[RB] = "RB";
    Direct_string[UF] = "UF";
    Direct_string[UB] = "UB";
    Direct_string[DF] = "DF";
    Direct_string[DB] = "DB";

    Direct_string[OMEGA] = "OMEGA";


    // Get opposite face, e.g. Reflect_face(L)=R
    Reflect_face.resize(27);
    Reflect_face[L] = R;
    Reflect_face[R] = L;
    Reflect_face[U] = D;
    Reflect_face[D] = U;
    Reflect_face[B] = F;
    Reflect_face[F] = B;

    // Get opposite edge, e.g. Reflect_edge(DB)=UF
    Reflect_edge.resize(27);
    Reflect_edge[LB] = RF;
    Reflect_edge[RB] = LF;
    Reflect_edge[DB] = UF;
    Reflect_edge[UB] = DF;

    Reflect_edge[LD] = RU;
    Reflect_edge[RD] = LU;
    Reflect_edge[LU] = RD;
    Reflect_edge[RU] = LD;

    Reflect_edge[LF] = RB;
    Reflect_edge[RF] = LB;
    Reflect_edge[DF] = UB;
    Reflect_edge[UF] = DB;

    // Get opposite vertex, e.g. Reflect_vertex(LDB)=RUF
    Reflect_vertex.resize(27);
    Reflect_vertex[LDB] = RUF;
    Reflect_vertex[RUF] = LDB;
    Reflect_vertex[RDB] = LUF;
    Reflect_vertex[LUF] = RDB;
    Reflect_vertex[LUB] = RDF;
    Reflect_vertex[RDF] = LUB;
    Reflect_vertex[RUB] = LDF;
    Reflect_vertex[LDF] = RUB;


    // Vertices at ends of edges
    Vertex_at_end_of_edge.resize(27);

    Vertex_at_end_of_edge[DB].resize(2);
    Vertex_at_end_of_edge[DB][0] = LDB; // Pattern: both other indices
    Vertex_at_end_of_edge[DB][1] = RDB;

    Vertex_at_end_of_edge[UB].resize(2);
    Vertex_at_end_of_edge[UB][0] = LUB; // Pattern: both other indices
    Vertex_at_end_of_edge[UB][1] = RUB;


    Vertex_at_end_of_edge[LB].resize(2);
    Vertex_at_end_of_edge[LB][0] = LUB; // Pattern: both other indices
    Vertex_at_end_of_edge[LB][1] = LDB;

    Vertex_at_end_of_edge[RB].resize(2);
    Vertex_at_end_of_edge[RB][0] = RUB; // Pattern: both other indices
    Vertex_at_end_of_edge[RB][1] = RDB;


    Vertex_at_end_of_edge[LD].resize(2);
    Vertex_at_end_of_edge[LD][0] = LDF; // Pattern: both other indices
    Vertex_at_end_of_edge[LD][1] = LDB;

    Vertex_at_end_of_edge[RD].resize(2);
    Vertex_at_end_of_edge[RD][0] = RDF; // Pattern: both other indices
    Vertex_at_end_of_edge[RD][1] = RDB;

    Vertex_at_end_of_edge[LU].resize(2);
    Vertex_at_end_of_edge[LU][0] = LUF; // Pattern: both other indices
    Vertex_at_end_of_edge[LU][1] = LUB;

    Vertex_at_end_of_edge[RU].resize(2);
    Vertex_at_end_of_edge[RU][0] = RUF; // Pattern: both other indices
    Vertex_at_end_of_edge[RU][1] = RUB;


    Vertex_at_end_of_edge[DF].resize(2);
    Vertex_at_end_of_edge[DF][0] = LDF; // Pattern: both other indices
    Vertex_at_end_of_edge[DF][1] = RDF;

    Vertex_at_end_of_edge[UF].resize(2);
    Vertex_at_end_of_edge[UF][0] = LUF; // Pattern: both other indices
    Vertex_at_end_of_edge[UF][1] = RUF;

    Vertex_at_end_of_edge[LF].resize(2);
    Vertex_at_end_of_edge[LF][0] = LUF; // Pattern: both other indices
    Vertex_at_end_of_edge[LF][1] = LDF;

    Vertex_at_end_of_edge[RF].resize(2);
    Vertex_at_end_of_edge[RF][0] = RUF; // Pattern: both other indices
    Vertex_at_end_of_edge[RF][1] = RDF;


    // Initialisation of the values of Vector_to_direction
    Vector<int> vect(3);
    int elem;

    for (int i = -1; i < 2; i++)
    {
      for (int j = -1; j < 2; j++)
      {
        for (int k = -1; k < 2; k++)
        {
          vect[0] = i;
          vect[1] = j;
          vect[2] = k;
          int num_elem = 0;

          // To put a number on the vector (i,j,k), we assume that that
          // the vector (i+1,j+1,k+1) represents the decomposition
          // of the number of the corresponding direction in base 3.
          num_elem = (i + 1) * 9 + (j + 1) * 3 + (k + 1);

          // for each number we have the corresponding element
          switch (num_elem)
          {
            case 6:
              elem = LUB;
              break;
            case 24:
              elem = RUB;
              break;
            case 26:
              elem = RUF;
              break;
            case 8:
              elem = LUF;
              break;
            case 0:
              elem = LDB;
              break;
            case 18:
              elem = RDB;
              break;
            case 20:
              elem = RDF;
              break;
            case 2:
              elem = LDF;
              break;
            case 25:
              elem = RU;
              break;
            case 23:
              elem = RF;
              break;
            case 19:
              elem = RD;
              break;
            case 21:
              elem = RB;
              break;
            case 7:
              elem = LU;
              break;
            case 5:
              elem = LF;
              break;
            case 1:
              elem = LD;
              break;
            case 3:
              elem = LB;
              break;
            case 17:
              elem = UF;
              break;
            case 15:
              elem = UB;
              break;
            case 11:
              elem = DF;
              break;
            case 9:
              elem = DB;
              break;
            case 16:
              elem = U;
              break;
            case 10:
              elem = D;
              break;
            case 22:
              elem = R;
              break;
            case 4:
              elem = L;
              break;
            case 14:
              elem = F;
              break;
            case 12:
              elem = B;
              break;
            case 13:
              elem = OMEGA;
              break;
            default:
              elem = OMEGA;
              oomph_info << "there might be a problem with Vector_to_direction"
                         << std::endl;
              break;
          }
          Vector_to_direction[vect] = elem;
        }
      }
    }


    // Initialisation of Direction_to_vector:
    // Translate Octant, face, edge into direction vector using the
    // value of Direct_string; Direction_to_vector[U]={0,1,0};
    Direction_to_vector.resize(27);
    for (int i = LDB; i <= F; i++)
    {
      Direction_to_vector[i].resize(3);
      // Initialisation to 0;
      for (int j = 0; j < 3; j++)
      {
        Direction_to_vector[i][j] = 0;
      }

      // Use +1 or -1 to indicate the relevant components of
      // the vector Direction
      for (unsigned j = 0; j < 3; j++)
      {
        if (Direct_string[i].length() > j)
        {
          switch (Direct_string[i][j])
          {
            case 'R':
              Direction_to_vector[i][0] = 1;
              break;
            case 'L':
              Direction_to_vector[i][0] = -1;
              break;
            case 'U':
              Direction_to_vector[i][1] = 1;
              break;
            case 'D':
              Direction_to_vector[i][1] = -1;
              break;
            case 'F':
              Direction_to_vector[i][2] = 1;
              break;
            case 'B':
              Direction_to_vector[i][2] = -1;
              break;
            default:
              oomph_info << "Direction Error !!" << std::endl;
          }
        }
      }
    }


    // Setup map that works out required rotations based on
    //-----------------------------------------------------
    // adjacent edge vertices
    //-----------------------

    int new_up, new_right;
    int new_vertex;


    // Map that stores the set of rotations (as a pairs consisting of
    // the up_equivalent and the right_equivalent) that move the vertex
    // specified by the first entry in key pair to the position of the second
    // one:
    std::map<std::pair<int, int>, std::set<std::pair<int, int>>>
      required_rotation;

    // Loop over all vertices
    for (int vertex = LDB; vertex <= RUF; vertex++)
    {
      new_up = U;
      new_right = R;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = U;
      new_right = B;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = U;
      new_right = L;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = U;
      new_right = F;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));


      new_up = D;
      new_right = R;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = D;
      new_right = B;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = D;
      new_right = L;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = D;
      new_right = F;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));


      new_up = R;
      new_right = D;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = R;
      new_right = B;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = R;
      new_right = U;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = R;
      new_right = F;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));


      new_up = L;
      new_right = D;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = L;
      new_right = B;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = L;
      new_right = U;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = L;
      new_right = F;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));


      new_up = F;
      new_right = R;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = F;
      new_right = L;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = F;
      new_right = U;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = F;
      new_right = D;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));


      new_up = B;
      new_right = R;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = B;
      new_right = L;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = B;
      new_right = U;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));

      new_up = B;
      new_right = D;
      new_vertex = OcTree::rotate(new_up, new_right, vertex);
      required_rotation[std::make_pair(vertex, new_vertex)].insert(
        std::make_pair(new_up, new_right));
    }


    // Each vertex is part of three edges. This container stores the
    // vertices in each of the three edge neighbours that are
    // adjacent to this node if there's no relative rotation between
    // the elements.
    std::map<int, Vector<int>> vertex_in_edge_neighbour;


    // Each vertex is part of three edges. This container stores the
    // vertices at the other end of the edge
    std::map<int, Vector<int>> other_vertex_on_edge;

    // Each vertex is part of three edges. This container stores the
    // vertices in the adjacent element at the other end of the edge
    // assuming there are no rotations between the elements.
    std::map<int, Vector<int>> other_vertex_in_edge_neighbour;


    vertex_in_edge_neighbour[LDB].resize(3);
    vertex_in_edge_neighbour[LDB][0] =
      LUF; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[LDB][1] = RDF;
    vertex_in_edge_neighbour[LDB][2] = RUB;

    other_vertex_on_edge[LDB].resize(3);
    other_vertex_on_edge[LDB][0] =
      RDB; // Pattern: opposite of the matching letter
    other_vertex_on_edge[LDB][1] = LUB;
    other_vertex_on_edge[LDB][2] = LDF;

    other_vertex_in_edge_neighbour[LDB].resize(3);
    other_vertex_in_edge_neighbour[LDB][0] = RUF; // Pattern: full reflection
    other_vertex_in_edge_neighbour[LDB][1] = RUF;
    other_vertex_in_edge_neighbour[LDB][2] = RUF;


    vertex_in_edge_neighbour[RDB].resize(3);
    vertex_in_edge_neighbour[RDB][0] =
      RUF; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[RDB][1] = LDF;
    vertex_in_edge_neighbour[RDB][2] = LUB;

    other_vertex_on_edge[RDB].resize(3);
    other_vertex_on_edge[RDB][0] =
      LDB; // Pattern: opposite of the matching letter
    other_vertex_on_edge[RDB][1] = RUB;
    other_vertex_on_edge[RDB][2] = RDF;

    other_vertex_in_edge_neighbour[RDB].resize(3);
    other_vertex_in_edge_neighbour[RDB][0] = LUF; // Pattern: full reflection
    other_vertex_in_edge_neighbour[RDB][1] = LUF;
    other_vertex_in_edge_neighbour[RDB][2] = LUF;


    vertex_in_edge_neighbour[LUB].resize(3);
    vertex_in_edge_neighbour[LUB][0] =
      LDF; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[LUB][1] = RUF;
    vertex_in_edge_neighbour[LUB][2] = RDB;

    other_vertex_on_edge[LUB].resize(3);
    other_vertex_on_edge[LUB][0] =
      RUB; // Pattern: opposite of the matching letter
    other_vertex_on_edge[LUB][1] = LDB;
    other_vertex_on_edge[LUB][2] = LUF;

    other_vertex_in_edge_neighbour[LUB].resize(3);
    other_vertex_in_edge_neighbour[LUB][0] = RDF; // Pattern: full reflection
    other_vertex_in_edge_neighbour[LUB][1] = RDF;
    other_vertex_in_edge_neighbour[LUB][2] = RDF;


    vertex_in_edge_neighbour[RUB].resize(3);
    vertex_in_edge_neighbour[RUB][0] =
      RDF; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[RUB][1] = LUF;
    vertex_in_edge_neighbour[RUB][2] = LDB;

    other_vertex_on_edge[RUB].resize(3);
    other_vertex_on_edge[RUB][0] =
      LUB; // Pattern: opposite of the matching letter
    other_vertex_on_edge[RUB][1] = RDB;
    other_vertex_on_edge[RUB][2] = RUF;

    other_vertex_in_edge_neighbour[RUB].resize(3);
    other_vertex_in_edge_neighbour[RUB][0] = LDF; // Pattern: full reflection
    other_vertex_in_edge_neighbour[RUB][1] = LDF;
    other_vertex_in_edge_neighbour[RUB][2] = LDF;


    vertex_in_edge_neighbour[LDF].resize(3);
    vertex_in_edge_neighbour[LDF][0] =
      LUB; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[LDF][1] = RDB;
    vertex_in_edge_neighbour[LDF][2] = RUF;

    other_vertex_on_edge[LDF].resize(3);
    other_vertex_on_edge[LDF][0] =
      RDF; // Pattern: opposite of the matching letter
    other_vertex_on_edge[LDF][1] = LUF;
    other_vertex_on_edge[LDF][2] = LDB;

    other_vertex_in_edge_neighbour[LDF].resize(3);
    other_vertex_in_edge_neighbour[LDF][0] = RUB; // Pattern: full reflection
    other_vertex_in_edge_neighbour[LDF][1] = RUB;
    other_vertex_in_edge_neighbour[LDF][2] = RUB;


    vertex_in_edge_neighbour[RDF].resize(3);
    vertex_in_edge_neighbour[RDF][0] =
      RUB; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[RDF][1] = LDB;
    vertex_in_edge_neighbour[RDF][2] = LUF;

    other_vertex_on_edge[RDF].resize(3);
    other_vertex_on_edge[RDF][0] =
      LDF; // Pattern: opposite of the matching letter
    other_vertex_on_edge[RDF][1] = RUF;
    other_vertex_on_edge[RDF][2] = RDB;

    other_vertex_in_edge_neighbour[RDF].resize(3);
    other_vertex_in_edge_neighbour[RDF][0] = LUB; // Pattern: full reflection
    other_vertex_in_edge_neighbour[RDF][1] = LUB;
    other_vertex_in_edge_neighbour[RDF][2] = LUB;


    vertex_in_edge_neighbour[LUF].resize(3);
    vertex_in_edge_neighbour[LUF][0] =
      LDB; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[LUF][1] = RUB;
    vertex_in_edge_neighbour[LUF][2] = RDF;

    other_vertex_on_edge[LUF].resize(3);
    other_vertex_on_edge[LUF][0] =
      RUF; // Pattern: opposite of the matching letter
    other_vertex_on_edge[LUF][1] = LDF;
    other_vertex_on_edge[LUF][2] = LUB;

    other_vertex_in_edge_neighbour[LUF].resize(3);
    other_vertex_in_edge_neighbour[LUF][0] = RDB; // Pattern: full reflection
    other_vertex_in_edge_neighbour[LUF][1] = RDB;
    other_vertex_in_edge_neighbour[LUF][2] = RDB;


    vertex_in_edge_neighbour[RUF].resize(3);
    vertex_in_edge_neighbour[RUF][0] =
      RDB; // Pattern: exactly one letter matches
    vertex_in_edge_neighbour[RUF][1] = LUB;
    vertex_in_edge_neighbour[RUF][2] = LDF;

    other_vertex_on_edge[RUF].resize(3);
    other_vertex_on_edge[RUF][0] =
      LUF; // Pattern: opposite of the matching letter
    other_vertex_on_edge[RUF][1] = RDF;
    other_vertex_on_edge[RUF][2] = RUB;

    other_vertex_in_edge_neighbour[RUF].resize(3);
    other_vertex_in_edge_neighbour[RUF][0] = LDB; // Pattern: full reflection
    other_vertex_in_edge_neighbour[RUF][1] = LDB;
    other_vertex_in_edge_neighbour[RUF][2] = LDB;


    // Loop over all vertices in the reference element
    for (int vertex = LDB; vertex <= RUF; vertex++)
    {
      // Loop over the three edges that are connected to this vertex
      for (unsigned i = 0; i < 3; i++)
      {
        // This is the other vertex along this edge
        int other_vertex = other_vertex_on_edge[vertex][i];

        // This is the vertex in the edge neighbour that corresponds
        // to the present vertex in the reference element (in
        // the absence of rotations)
        int unrotated_neigh_vertex = vertex_in_edge_neighbour[vertex][i];

        // This is the vertex in the edge neighbour that corresponds
        // to the other vertex in the reference element (in
        // the absence of rotations)
        int unrotated_neigh_other_vertex =
          other_vertex_in_edge_neighbour[vertex][i];

        // Loop over all vertices in the neighbour element
        for (int neigh_vertex = LDB; neigh_vertex <= RUF; neigh_vertex++)
        {
          // What rotations would turn the neigh_vertex
          // into the unrotated_neigh_vertex?
          std::set<std::pair<int, int>> vertex_rot =
            required_rotation[std::make_pair(neigh_vertex,
                                             unrotated_neigh_vertex)];


          // Loop over all "other" vertices in the neighbour element
          for (int neigh_other_vertex = LDB; neigh_other_vertex <= RUF;
               neigh_other_vertex++)
          {
            // What rotations would turn the other_neigh_vertex
            // into the unrotated_other_neigh_vertex?
            std::set<std::pair<int, int>> other_vertex_rot =
              required_rotation[std::make_pair(neigh_other_vertex,
                                               unrotated_neigh_other_vertex)];

            // What are the common rotations?
            std::set<std::pair<int, int>> common_rotations;

            // Get the intersection of the two sets
            std::set_intersection(
              vertex_rot.begin(),
              vertex_rot.end(),
              other_vertex_rot.begin(),
              other_vertex_rot.end(),
              inserter(common_rotations, common_rotations.begin()));


            if (common_rotations.size() > 0)
            {
              for (std::set<std::pair<int, int>>::iterator it =
                     common_rotations.begin();
                   it != common_rotations.end();
                   it++)
              {
                // Copy into container

                // First: up equivalent:
                Up_and_right_equivalent_for_pairs_of_vertices
                  [std::make_pair(
                     std::make_pair(vertex, neigh_vertex),
                     std::make_pair(other_vertex, neigh_other_vertex))]
                    .first = it->first;

                // Second: Right equivalent
                Up_and_right_equivalent_for_pairs_of_vertices
                  [std::make_pair(
                     std::make_pair(vertex, neigh_vertex),
                     std::make_pair(other_vertex, neigh_other_vertex))]
                    .second = it->second;
              }
            }
          }
        }
      }
    }
  }


  //================================================================
  /// Is the edge neighbour (for edge "edge")  specified via the pointer
  /// also a face neighbour for one of the two adjacent faces?
  //================================================================
  bool OcTree::edge_neighbour_is_face_neighbour(const int& edge,
                                                OcTree* edge_neigh_pt) const
  {
#ifdef PARANOID
    // No paranoid check needed -- the default for the switch statement
    // catches illegal values for edge
#endif


    // Catch stupid case: Null doesn't have a face neighbour...
    if (edge_neigh_pt == 0)
    {
      return false;
    }

    using namespace OcTreeNames;

    // Auxiliary variables
    int face;
    Vector<unsigned> translate_s(3);
    Vector<double> s_sw(3);
    Vector<double> s_ne(3);
    int reflected_face;
    int diff_level;
    bool in_neighbouring_tree = false;

    OcTree* face_neigh_pt = 0;

    switch (edge)
    {
      case LB:

        // Get first face neighbour
        face = L;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = B;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case RB:


        // Get first face neighbour
        face = R;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = B;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);
        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case DB:

        // Get first face neighbour
        face = D;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = B;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);
        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case UB:

        // Get first face neighbour
        face = U;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = B;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;

      case LD:


        // Get first face neighbour
        face = L;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = D;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);
        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;

      case RD:


        // Get first face neighbour
        face = R;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = D;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);
        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;

      case LU:

        // Get first face neighbour
        face = L;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = U;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case RU:

        // Get first face neighbour
        face = R;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = U;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case LF:


        // Get first face neighbour
        face = L;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = F;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;

      case RF:

        // Get first face neighbour
        face = R;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        // Get second face neighbour
        face = F;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case DF:

        // Get first face neighbour
        face = D;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }


        // Get second face neighbour
        face = F;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;


      case UF:

        // Get first face neighbour
        face = U;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }


        // Get second face neighbour
        face = F;
        face_neigh_pt = gteq_face_neighbour(face,
                                            translate_s,
                                            s_sw,
                                            s_ne,
                                            reflected_face,
                                            diff_level,
                                            in_neighbouring_tree);

        // Check if they agree...
        if (face_neigh_pt != 0)
        {
          if (face_neigh_pt == edge_neigh_pt)
          {
            return true;
          }
        }

        break;

      default:

        // There is no face neighbour in this direction so they can't
        // agree:
        std::ostringstream error_stream;
        error_stream << "Never get here! Edge:" << Direct_string[edge] << " "
                     << edge << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // If we've made it to here then we've located the requested edge
    // but found that none of its two face neighbours match the specified
    // edge neighbour:
    return false;
  }


  //================================================================
  /// Find (pointer to) `greater-or-equal-sized face neighbour' in
  /// given direction (L/R/U/D/F/B).
  /// Another way of interpreting this is that we're looking for
  /// the neighbour across the present element's face 'direction'.
  /// The various arguments return additional information about the
  /// size and relative orientation of the neighbouring octree.
  /// To interpret these we use the following
  /// <B>General convention:</B>
  /// - Each face of the element that is represented by the octree
  ///   is parametrised by two (of the three)
  ///   local coordinates that parametrise the entire 3D element. E.g.
  ///   the B[ack] face is parametrised by (s[0], s[1]); the D[own] face
  ///   is parametrised by (s[0],s[2]); etc. We always identify the
  ///   in-face coordinate with the lower (3D) index with the subscript
  ///   _lo and the one with the larger (3D) index with the subscript _hi.
  /// .
  /// With this convention, the interpretation of the arguments is
  /// as follows:
  /// - The vector \c translate_s turns the index of the local coordinate
  ///   in the present octree into that of the neighbour. If there are no
  ///   rotations then \c translate_s[i] = i.
  /// - In the present octree, the "south west" vertex of the face
  ///   between the present octree and its neighbour is located at
  ///   S_lo=-1, S_hi=-1. This point is located at the (3D) local
  ///   coordinates (\c s_sw[0], \c s_sw[1], \c s_sw[2])
  ///   in the neighbouring octree.
  /// - ditto with s_ne: In the present octree, the "north east" vertex
  ///   of the face between the present octree and its neighbour is located at
  ///   S_lo=+1, S_hi=+1. This point is located
  ///   at the (3D) local coordinates (\c s_ne[0], \c s_ne[1], \c s_ne[2])
  ///   in the neighbouring octree.
  /// - We're looking for a neighbour in the specified \c direction. When
  ///   viewed from the neighbouring octree, the face that separates
  ///   the present octree from its neighbour is the neighbour's face
  ///   \c face. If there's no rotation between the two octrees, this is a
  ///   simple reflection: For instance, if we're looking
  ///   for a neighhbour in the \c R [ight] \c direction, \c face will
  ///   be \c L [eft]
  /// - \c diff_level <= 0 indicates the difference in refinement levels between
  ///   the two neighbours. If \c diff_level==0, the neighbour has the
  ///   same size as the current octree.
  //=================================================================
  OcTree* OcTree::gteq_face_neighbour(const int& direction,
                                      Vector<unsigned>& translate_s,
                                      Vector<double>& s_sw,
                                      Vector<double>& s_ne,
                                      int& face,
                                      int& diff_level,
                                      bool& in_neighbouring_tree) const
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((direction != L) && (direction != R) && (direction != F) &&
        (direction != B) && (direction != U) && (direction != D))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong direction: " << Direct_string[direction]
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise in_neighbouring tree to false. It will be set to true
    // during the recursion if we do actually hop over in to the neighbour
    in_neighbouring_tree = false;

    // Maximum level to which we're allowed to descend (we only want
    // greater-or-equal-sized neighbours)
    int max_level = Level;

    // Current element has the following root:
    OcTreeRoot* orig_root_pt = dynamic_cast<OcTreeRoot*>(Root_pt);

    // Initialise offset in local coordinate
    double s_difflo = 0;
    double s_diffhi = 0;

    // Initialise difference in level
    diff_level = 0;

    // Find neighbour
    OcTree* return_pt = gteq_face_neighbour(direction,
                                            s_difflo,
                                            s_diffhi,
                                            diff_level,
                                            in_neighbouring_tree,
                                            max_level,
                                            orig_root_pt);

    OcTree* neighb_pt = return_pt;

    // Initialise the translation scheme
    for (unsigned i = 0; i < 3; i++)
    {
      translate_s[i] = i;
    }

    // If neighbour exists: What's the direction of the interfacial
    // face when viewed from within the neighbour element?
    if (neighb_pt != 0)
    {
      // Find the reflection of the original direction, which will be the
      // direction to the face in the neighbour, if there are no rotations.
      int reflected_dir = Reflect_face[direction];

      // These coordinates are the coordinates of the ne and sw points
      // in the neighbour (provided there are no rotations -- we'll correct
      // for these below)
      s_sw[0] = S_base(0, reflected_dir) +
                S_steplo(0, reflected_dir) * s_difflo +
                S_stephi(0, reflected_dir) * s_diffhi;

      s_sw[1] = S_base(1, reflected_dir) +
                S_steplo(1, reflected_dir) * s_difflo +
                S_stephi(1, reflected_dir) * s_diffhi;

      s_sw[2] = S_base(2, reflected_dir) +
                S_steplo(2, reflected_dir) * s_difflo +
                S_stephi(2, reflected_dir) * s_diffhi;

      s_ne[0] = S_base(0, reflected_dir) +
                S_steplo(0, reflected_dir) * pow(2.0, diff_level) +
                S_steplo(0, reflected_dir) * s_difflo +
                S_stephi(0, reflected_dir) * pow(2.0, diff_level) +
                S_stephi(0, reflected_dir) * s_diffhi;

      s_ne[1] = S_base(1, reflected_dir) +
                S_steplo(1, reflected_dir) * pow(2.0, diff_level) +
                S_steplo(1, reflected_dir) * s_difflo +
                S_stephi(1, reflected_dir) * pow(2.0, diff_level) +
                S_stephi(1, reflected_dir) * s_diffhi;

      s_ne[2] = S_base(2, reflected_dir) +
                S_steplo(2, reflected_dir) * pow(2.0, diff_level) +
                S_steplo(2, reflected_dir) * s_difflo +
                S_stephi(2, reflected_dir) * pow(2.0, diff_level) +
                S_stephi(2, reflected_dir) * s_diffhi;

      // If there is no rotation then the new direction is the same as the
      // old direction
      int new_dir = direction;

      // If necessary, rotate things around (the orientation of Up in the
      // neighbour might be be different from that in the present element)
      // If the root of the neighbour is the not same as the one of the present
      // element then their orientations may not be the same and the new
      // direction is given by :
      if (neighb_pt->Root_pt != Root_pt)
      {
        new_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         direction);
      }

      // What's the direction of the interfacial face when viewed from within
      // the neighbour element?
      face = Reflect_face[new_dir];

      Vector<double> s_sw_new(3), s_ne_new(3);

      // If the root of the present element is different from the root
      // of his neighbour, we have to rotate the RUF and LDB coordinates
      // to have their equivalents in the neighbour's point of view.
      if (neighb_pt->Root_pt != Root_pt)
      {
        int tmp_dir;
        Vector<int> vect1(3);
        Vector<int> vect2(3);
        Vector<int> vect3(3);
        DenseMatrix<int> Mat_rot(3);

        // All this is just to compute the rotation matrix
        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         R);
        vect1 = Direction_to_vector[tmp_dir];

        // All this is just to compute the rotation matrix
        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         U);
        vect2 = Direction_to_vector[tmp_dir];

        // All this is just to compute the rotation matrix
        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         F);
        vect3 = Direction_to_vector[tmp_dir];

        // Setup the inverse rotation matrix
        for (int i = 0; i < 3; i++)
        {
          Mat_rot(i, 0) = vect1[i];
          Mat_rot(i, 1) = vect2[i];
          Mat_rot(i, 2) = vect3[i];
        }

        // Initialise the translation scheme
        Vector<int> translate_s_new(3);

        // Then the rotation of the coordinates
        for (unsigned i = 0; i < 3; i++)
        {
          s_ne_new[i] = 0.0;
          s_sw_new[i] = 0.0;
          translate_s_new[i] = 0;
          for (unsigned k = 0; k < 3; k++)
          {
            s_ne_new[i] += Mat_rot(i, k) * s_ne[k];
            s_sw_new[i] += Mat_rot(i, k) * s_sw[k];
            translate_s_new[i] += Mat_rot(i, k) * translate_s[k];
          }
        }

        s_ne = s_ne_new;
        s_sw = s_sw_new;

        // Set the translation scheme
        for (unsigned i = 0; i < 3; i++)
        {
          // abs is ok here; not fabs!
          translate_s[i] = std::abs(translate_s_new[i]);
        }
      }

    } // end of if(neighb_pt!=0)

    return return_pt;
  }

  //================================================================
  /// Find (pointer to) `greater-or-equal-sized true edge neighbour' in
  /// the given direction (LB,RB,DB,UB [the back edges],
  /// LD,RD,LU,RU [the side edges], LF,RF,DF,UF [the front edges]).
  /// Another way of interpreting this is that we're looking for
  /// the neighbour across the present element's edge 'direction'.
  /// The various arguments return additional information about the
  /// size and relative orientation of the neighbouring octree.
  /// Each edge of the element that is represented by the octree
  /// is parametrised by one (of the three) local coordinates that
  /// parametrise the entire 3D element. E.g. the L[eft]B[ack] edge
  /// is parametrised by s[1]; the "low" vertex of this edge
  /// (located at the low value of this coordinate, i.e. at s[1]=-1)
  /// is L[eft]D[own]B[ack]. The "high" vertex of this edge (located
  /// at the high value of this coordinate, i.e. at s[1]=1) is
  /// L[eft]U[p]B[ack]; etc
  /// The interpretation of the arguments is as follows:
  /// - In a forest, an OcTree can have multiple edge neighbours
  ///   (across an edge where multiple trees meet). \c i_root_edge_neighbour
  ///   specifies which of these is used. Use this as "reverse communication":
  ///   First call with \c i_root_edge_neighbour=0 and \c n_root_edge_neighour
  ///   initialised to anything you want (zero, ideally). On return from
  ///   the fct, \c n_root_edge_neighour contains the total number of true
  ///   edge neighbours, so additional calls to the fct with
  ///   \c i_root_edge_neighbour>0 can be made until they've all been visited.
  /// - The vector \c translate_s turns the index of the local coordinate
  ///   in the present octree into that of the neighbour. If there are no
  ///   rotations then \c translate_s[i] = i.
  /// - The "low" vertex of the edge in the present octree
  ///   coincides with a certain vertex in the edge neighbour.
  ///   In terms of the neighbour's local coordinates, this point is
  ///   located at the (3D) local coordinates (\c s_lo[0], \c s_lo[1],
  ///   \c s_lo[2])
  /// - ditto with s_hi: The "high" vertex of the edge in the present octree
  ///   coincides with a certain vertex in the edge neighbour.
  ///   In terms of the neighbour's local coordinates, this point is
  ///   located at the (3D) local coordinates (\c s_hi[0], \c s_hi[1],
  ///   \c s_hi[2])
  /// - We're looking for a neighbour in the specified \c direction. When
  ///   viewed from the neighbouring octree, the edge that separates
  ///   the present octree from its neighbour is the neighbour's edge
  ///   \c edge. If there's no rotation between the two octrees, this is a
  ///   simple reflection: For instance, if we're looking
  ///   for a neighhbour in the \c DB \c direction, \c edge will
  ///   be \c UF.
  /// - \c diff_level <= 0 indicates the difference in refinement levels between
  ///   the two neighbours. If \c diff_level==0, the neighbour has the
  ///   same size as the current octree.
  /// .
  /// \b Important: We're only looking for \b true edge neighbours
  /// i.e. edge neigbours that are not also face neighbours. This is an
  /// important difference to Samet's terminology. If the neighbour
  /// in a certain direction is not a true edge neighbour, or if there
  /// is no neighbour, then this function returns NULL.
  //=================================================================
  OcTree* OcTree::gteq_true_edge_neighbour(
    const int& direction,
    const unsigned& i_root_edge_neighbour,
    unsigned& nroot_edge_neighbour,
    Vector<unsigned>& translate_s,
    Vector<double>& s_lo,
    Vector<double>& s_hi,
    int& edge,
    int& diff_level) const
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((direction != LB) && (direction != RB) && (direction != DB) &&
        (direction != UB) && (direction != LD) && (direction != RD) &&
        (direction != LU) && (direction != RU) && (direction != LF) &&
        (direction != RF) && (direction != DF) && (direction != UF))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong direction: " << Direct_string[direction]
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Maximum level to which we're allowed to descend (we only want
    // greater-or-equal-sized neighbours)
    int max_level = Level;

    // Current element has the following root:
    OcTreeRoot* orig_root_pt = dynamic_cast<OcTreeRoot*>(Root_pt);

    // Initialise offset in local coordinate along edge
    double s_diff = 0;

    // Initialise difference in level
    diff_level = 0;

    // Find edge neighbour
    OcTree* return_pt = gteq_edge_neighbour(direction,
                                            i_root_edge_neighbour,
                                            nroot_edge_neighbour,
                                            s_diff,
                                            diff_level,
                                            max_level,
                                            orig_root_pt);

    // Only use "true" edge neighbours
    if (edge_neighbour_is_face_neighbour(direction, return_pt))
    {
      return_pt = 0;
    }

    // By default, we return what was returned as the true edge neighbour.
    OcTree* neighb_pt = return_pt;

    // Initialise the translation scheme
    for (unsigned i = 0; i < 3; i++)
    {
      translate_s[i] = i;
    }

    // If neighbour exists: What's the direction of the interfacial
    // edge when viewed from within the neighbour element?
    if (neighb_pt != 0)
    {
      // Find the reflection of the original direction, which will be the
      // direction to the edge in the neighbour, if there are no rotations.
      int reflected_dir = Reflect_edge[direction];

      // These coordinates are the coordinates of the "low" and "high" points
      // in the neighbour (provided there are no rotations -- we'll correct
      // for these below)
      s_lo[0] =
        S_base_edge(0, reflected_dir) + S_step_edge(0, reflected_dir) * s_diff;

      s_lo[1] =
        S_base_edge(1, reflected_dir) + S_step_edge(1, reflected_dir) * s_diff;

      s_lo[2] =
        S_base_edge(2, reflected_dir) + S_step_edge(2, reflected_dir) * s_diff;

      s_hi[0] = S_base_edge(0, reflected_dir) +
                S_step_edge(0, reflected_dir) * pow(2.0, diff_level) +
                S_step_edge(0, reflected_dir) * s_diff;

      s_hi[1] = S_base_edge(1, reflected_dir) +
                S_step_edge(1, reflected_dir) * pow(2.0, diff_level) +
                S_step_edge(1, reflected_dir) * s_diff;

      s_hi[2] = S_base_edge(2, reflected_dir) +
                S_step_edge(2, reflected_dir) * pow(2.0, diff_level) +
                S_step_edge(2, reflected_dir) * s_diff;

      // If there is no rotation then the new direction is the same as the
      // old direction
      int new_dir = direction;


      // If necessary, rotate things around (the orientation of Up in the
      // neighbour might be be different from that in the present element)
      // If the root of the neighbour is the not same as the one of the present
      // element then their orientations may not be the same and the new
      // direction is given by :
      if (neighb_pt->Root_pt != Root_pt)
      {
        int new_up = orig_root_pt->up_equivalent(neighb_pt->Root_pt);

        int new_right = orig_root_pt->right_equivalent(neighb_pt->Root_pt);

        new_dir = rotate(new_up, new_right, direction);
      }

      // What's the direction of the interfacial edge when viewed from within
      // the neighbour element (including rotations!)
      edge = Reflect_edge[new_dir];

      // Get ready to rotate the local coordinates
      Vector<double> s_lo_new(3), s_hi_new(3);

      // If the root of the present element is different from the root
      // of his neighbour, we have to rotate the lo and hi coordinates
      // to have their equivalents from the neighbour's point of view.
      if ((neighb_pt->Root_pt != Root_pt))
      {
        int tmp_dir;
        Vector<int> vect1(3);
        Vector<int> vect2(3);
        Vector<int> vect3(3);
        DenseMatrix<int> Mat_rot(3);

        // All this is just to compute the rotation matrix
        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         R);
        vect1 = Direction_to_vector[tmp_dir];


        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         U);
        vect2 = Direction_to_vector[tmp_dir];


        tmp_dir = rotate(orig_root_pt->up_equivalent(neighb_pt->Root_pt),
                         orig_root_pt->right_equivalent(neighb_pt->Root_pt),
                         F);
        vect3 = Direction_to_vector[tmp_dir];


        // Setup the inverse rotation matrix
        for (int i = 0; i < 3; i++)
        {
          Mat_rot(i, 0) = vect1[i];
          Mat_rot(i, 1) = vect2[i];
          Mat_rot(i, 2) = vect3[i];
        }

        // Initialise the translation scheme
        Vector<int> translate_s_new(3);

        // Then the rotation of the coordinates
        for (unsigned i = 0; i < 3; i++)
        {
          s_hi_new[i] = 0.0;
          s_lo_new[i] = 0.0;
          translate_s_new[i] = 0;
          for (unsigned k = 0; k < 3; k++)
          {
            s_hi_new[i] += Mat_rot(i, k) * s_hi[k];
            s_lo_new[i] += Mat_rot(i, k) * s_lo[k];
            translate_s_new[i] += Mat_rot(i, k) * translate_s[k];
          }
        }

        s_lo = s_lo_new;
        s_hi = s_hi_new;

        // Set the translation scheme
        for (unsigned i = 0; i < 3; i++)
        {
          // abs is ok here; not fabs!
          translate_s[i] = std::abs(translate_s_new[i]);
        }
      }

    } // end if for (neighb_pt!=0)

    return return_pt;
  }


  //==========================================================================
  /// Find `greater-or-equal-sized face neighbour' in given direction
  /// (L/R/U/D/B/F).
  ///
  /// This is an auxiliary routine which allows neighbour finding in adjacent
  /// octrees. Needs to keep track of previous son types and
  /// the maximum level to which search is performed.
  ///
  /// Parameters:
  ///
  /// - direction: (L/R/U/D/B/F) Direction in which neighbour has to be found.
  /// - s_difflo/s_diffhi: Offset of left/down/back vertex from
  ///   corresponding vertex in neighbour. Note that this is input/output
  ///   as it needs to be incremented/decremented during the recursive calls
  ///   to this function.
  /// - face: We're looking for the neighbour across our face 'direction'
  ///   (L/R/U/D/B/F). When viewed from the neighbour, this face is
  ///   `face' (L/R/U/D/B/F). [If there's no relative rotation between
  ///   neighbours then this is a mere reflection, e.g. direction=F  --> face=B
  ///   etc.]
  /// - diff_level <= 0 indicates the difference in octree levels
  ///   between the current element and its neighbour.
  /// - max_level is the maximum level to which the neighbour search is
  ///   allowed to proceed. This is necessary because in a forest,
  ///   the neighbour search isn't based on pure recursion.
  /// - orig_root_pt identifies the root node of the element whose
  ///   neighbour we're really trying to find by all these recursive calls.
  //===========================================================================
  OcTree* OcTree::gteq_face_neighbour(const int& direction,
                                      double& s_difflo,
                                      double& s_diffhi,
                                      int& diff_level,
                                      bool& in_neighbouring_tree,
                                      int max_level,
                                      OcTreeRoot* orig_root_pt) const
  {
    using namespace OcTreeNames;

#ifdef PARANOID
    if ((direction != L) && (direction != R) && (direction != F) &&
        (direction != B) && (direction != U) && (direction != D))
    {
      std::ostringstream error_stream;
      error_stream << "Direction " << Direct_string[direction]
                   << " is not L, R, B, F, D or U." << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    OcTree* next_el_pt;
    OcTree* return_el_pt;

    // Initialise in_neighbouring tree to false. It will be set to true
    // during the recursion if we do actually hop over in to the neighbour
    in_neighbouring_tree = false;

    // STEP 1: Find the neighbour's father
    //--------
    // Does the element have a father?
    if (Father_pt != 0)
    {
      // If the present octant (whose location inside its
      // father element is specified by Son_type) is adjacent to the
      // father's face in the required direction, then its neighbour has
      // a different father ---> we need to climb up the tree to
      // the father and find his neighbour in the required direction
      if (Is_adjacent(direction, Son_type))
      {
        next_el_pt = dynamic_cast<OcTree*>(Father_pt)->gteq_face_neighbour(
          direction,
          s_difflo,
          s_diffhi,
          diff_level,
          in_neighbouring_tree,
          max_level,
          orig_root_pt);
      }
      // If the present octant is not adjacent to the
      // father's face in the required direction, then the
      // neighbour has the same father and will be obtained
      // by the appropriate reflection inside the father element
      else
      {
        next_el_pt = dynamic_cast<OcTree*>(Father_pt);
      }

      // We're about to ascend one level:
      diff_level -= 1;

      // Work out position of lower-left corner of present face
      // in its father element
      s_difflo += pow(0.5, -diff_level) * S_directlo(direction, Son_type);
      s_diffhi += pow(0.5, -diff_level) * S_directhi(direction, Son_type);

      // STEP 2:  We have now located the neighbour's father and need to
      // -------  find the appropriate son.
      // Buffer cases where the neighbour (and hence its father) lie outside
      // the boundary
      if (next_el_pt != 0)
      {
        // If the father is a leaf then we can't descend to the same
        // level as the present node ---> simply return the father himself
        // as the (greater) neighbour. Same applies if we are about
        // to descend lower than the max_level (in a neighbouring tree)
        if ((next_el_pt->Son_pt.size() == 0) ||
            (next_el_pt->Level > max_level - 1))
        {
          return_el_pt = next_el_pt;
        }
        // We have located the neighbour's father: The position of the
        // neighbour is obtained by `reflecting' the position of the
        // node itself.
        else
        {
          // By default (in the absence of rotations) we obtain the
          // son octant by reflecting
          int son_octant = Reflect(direction, Son_type);

          // If there is a rotation, we rotate the son octant
          if (orig_root_pt != next_el_pt->Root_pt)
          {
            int my_up = dynamic_cast<OcTreeRoot*>(Root_pt)->up_equivalent(
              next_el_pt->Root_pt);
            int my_right = dynamic_cast<OcTreeRoot*>(Root_pt)->right_equivalent(
              next_el_pt->Root_pt);
            son_octant = rotate(my_up, my_right, son_octant);
          }

          return_el_pt = dynamic_cast<OcTree*>(next_el_pt->Son_pt[son_octant]);

          // Work out position of lower-left corner of present face
          // in next higher element
          s_difflo -= pow(0.5, -diff_level) * S_directlo(direction, Son_type);
          s_diffhi -= pow(0.5, -diff_level) * S_directhi(direction, Son_type);

          // We have just descended one level
          diff_level += 1;
        }
      }
      // The neighbour's father lies outside the boundary --> the neighbour
      // itself does too --> return NULL.
      else
      {
        return_el_pt = 0;
      }
    }
    // Element does not have a father --> check if it has a neighbouring
    // tree in the appropriate direction
    else
    {
      // Find neighbouring root
      if (Root_pt->neighbour_pt(direction) != 0)
      {
        // If we're in the neighbouring tree
        in_neighbouring_tree = true;

        // Return
        return_el_pt = dynamic_cast<OcTree*>(Root_pt->neighbour_pt(direction));
      }
      // No neighbouring tree, so there really is no neighbour --> return NULL
      else
      {
        return_el_pt = 0;
      }
    }

    // Return the appropriate OcTree pointer
    return return_el_pt;
  } // End of gteq_face_neighbour


  //==========================================================================
  /// Find `greater-or-equal-sized edge neighbour' in given direction
  ///  (LB,RB,DB,UB [the back edges],
  /// LD,RD,LU,RU [the side edges], LF,RF,DF,UF [the front edges]).
  ///
  /// This is an auxiliary routine which allows neighbour finding in adjacent
  /// octrees. Needs to keep track of previous son types and
  /// the maximum level to which search is performed.
  ///
  /// Parameters:
  ///
  /// - direction: (LB,RB/...) Direction in which neighbour has to be found.
  /// - In a forest, an OcTree can have multiple edge neighbours
  ///   (across an edge where multiple trees meet). \c i_root_edge_neighbour
  ///   specifies which of these is used. Use this as "reverse communication":
  ///   First call with \c i_root_edge_neighbour=0 and \c n_root_edge_neighour
  ///   initialised to anything you want (zero, ideally). On return from
  ///   the fct, \c n_root_edge_neighour contains the total number of true
  ///   edge neighbours, so additional calls to the fct with
  ///   \c i_root_edge_neighbour>0 can be made until they've all been visited.
  /// - s_diff: Offset of left/down/back vertex from
  ///   corresponding vertex in
  ///   neighbour. Note that this is input/output as it needs to be incremented/
  ///   decremented during the recursive calls to this function.
  /// - diff_level <= 0 indicates the difference in octree levels
  ///   between the current element and its neighbour.
  /// - max_level is the maximum level to which the neighbour search is
  ///   allowed to proceed. This is necessary because in a forest,
  ///   the neighbour search isn't based on pure recursion.
  /// - orig_root_pt identifies the root node of the element whose
  ///   neighbour we're really trying to find by all these recursive calls.
  //===========================================================================
  OcTree* OcTree::gteq_edge_neighbour(const int& direction,
                                      const unsigned& i_root_edge_neighbour,
                                      unsigned& nroot_edge_neighbour,
                                      double& s_diff,
                                      int& diff_level,
                                      int max_level,
                                      OcTreeRoot* orig_root_pt) const
  {
    using namespace OcTreeNames;


#ifdef PARANOID
    if ((direction != LB) && (direction != RB) && (direction != DB) &&
        (direction != UB) && (direction != LD) && (direction != RD) &&
        (direction != LU) && (direction != RU) && (direction != LF) &&
        (direction != RF) && (direction != DF) && (direction != UF))
    {
      std::ostringstream error_stream;
      error_stream << "Wrong direction: " << Direct_string[direction]
                   << std::endl;
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise total number of edge neighbours available across
    // edges of octree roots
    nroot_edge_neighbour = 0;

    OcTree* next_el_pt = 0;
    OcTree* return_el_pt = 0;

    // STEP 1: Find the common ancestor
    //--------
    // Does the element have a father?
    if (Father_pt != 0)
    {
      // If the present octant (whose location inside its
      // father element is specified by Son_type) is adjacent to the
      // father's edge in the required direction, then its neighbour has
      // a different father ---> we need to climb up the tree to
      // the father and find his neighbour in the required direction
      if (Is_adjacent(direction, Son_type))
      {
        next_el_pt = dynamic_cast<OcTree*>(Father_pt)->gteq_edge_neighbour(
          direction,
          i_root_edge_neighbour,
          nroot_edge_neighbour,
          s_diff,
          diff_level,
          max_level,
          orig_root_pt);
      }
      // If the present octant (whose location inside its
      // father element is specified by Son_type) is adjacent to the
      // father's face in the required direction, then its neighbour has
      // a different father ---> we need to climb up the tree to
      // the father and find his neighbour in the required direction,
      // crossing the face as we do so.
      else if (Common_face(direction, Son_type) != OMEGA)
      {
        // Initialise bool
        bool in_neighbouring_tree = false;

        // We're going across a face:
        double s_difflo = 0.0;
        double s_diffhi = 0.0;
        int diff_level_edge = 0;

        next_el_pt = dynamic_cast<OcTree*>(Father_pt)->gteq_face_neighbour(
          Common_face(direction, Son_type),
          s_difflo,
          s_diffhi,
          diff_level_edge,
          in_neighbouring_tree,
          max_level,
          orig_root_pt);
      }
      // If the present octant is not adjacent to the
      // father's face/edge in the required direction, then the
      // neighbour has the same father and will be obtained
      // by the appropriate reflection inside the father element
      else
      {
        next_el_pt = dynamic_cast<OcTree*>(Father_pt);
      }

      // We're about to ascend one level:
      diff_level -= 1;

      // Work out position of "low" vertex of present edge
      // in its father element
      s_diff += pow(0.5, -diff_level) * S_direct_edge(direction, Son_type);

      // STEP 2:  We have now located the neighbour's father and need to
      // -------  find the appropriate son.
      // Buffer cases where ...
      if (next_el_pt != 0)
      {
        // If the father is a leaf then we can't descend to the same
        // level as the present node ---> simply return the father himself
        // as the (greater) neighbour. Same applies if we are about
        // to descend lower than the max_level (in a  neighbouring tree)
        if ((next_el_pt->Son_pt.size() == 0) ||
            (next_el_pt->Level > max_level - 1))
        {
          return_el_pt = next_el_pt;
        }
        // We have located the neighbour's father: The position of the
        // neighbour is obtained by `reflecting' the position of the
        // node itself.
        else
        {
          // By default (in the absence of rotations) we obtain the
          // son octant by reflecting
          int son_octant = Reflect(direction, Son_type);

          // If there is a rotation, we rotate the son octant
          if (orig_root_pt != next_el_pt->Root_pt)
          {
            int my_up = dynamic_cast<OcTreeRoot*>(Root_pt)->up_equivalent(
              next_el_pt->Root_pt);
            int my_right = dynamic_cast<OcTreeRoot*>(Root_pt)->right_equivalent(
              next_el_pt->Root_pt);
            son_octant = rotate(my_up, my_right, son_octant);
          }

          return_el_pt = dynamic_cast<OcTree*>(next_el_pt->Son_pt[son_octant]);

          // Work out position of "low" vertex of present edge
          // in next higher element
          s_diff -= pow(0.5, -diff_level) * S_direct_edge(direction, Son_type);

          // We have just descended one level
          diff_level += 1;
        }
      }
      // The neighbour's father lies outside the boundary --> the neighbour
      // itself does too --> return NULL.
      else
      {
        return_el_pt = 0;
      }
    }
    // Element does not have a father --> check if it has a neighbouring
    // tree in the appropriate direction
    else
    {
      // Get total number of edge neighbours available across
      // edges of octree roots for return
      nroot_edge_neighbour =
        dynamic_cast<OcTreeRoot*>(Root_pt)->nedge_neighbour(direction);

      // Get vector of edge neighbours (if any) in appropriate direction
      Vector<TreeRoot*> edge_neighbour_pt =
        dynamic_cast<OcTreeRoot*>(Root_pt)->edge_neighbour_pt(direction);

      // Get the number of edge neighbours
      unsigned n_neigh = edge_neighbour_pt.size();

      // If there are any edge neighbours
      if (n_neigh > 0)
      {
        // Return the appropriate edge neighbour
        return_el_pt =
          dynamic_cast<OcTree*>(edge_neighbour_pt[i_root_edge_neighbour]);
      }
      else
      {
        return_el_pt = 0;
      }
    }

    return return_el_pt;
  } // End of gteq_edge_neighbour


  //================================================================
  /// Self-test: Check neighbour finding routine. For each element
  /// in the tree and for each vertex, determine the
  /// distance between the vertex and its position in the
  /// neigbour. . If the difference is less than
  /// Tree::Max_neighbour_finding_tolerance.
  /// return success (0), otherwise failure (1)
  //=================================================================
  unsigned OcTree::self_test()
  {
    // Stick pointers to all nodes into Vector and number elements
    // in the process
    Vector<Tree*> all_nodes_pt;
    stick_all_tree_nodes_into_vector(all_nodes_pt);

    long int count = 0;
    unsigned long num_nodes = all_nodes_pt.size();

    for (unsigned long i = 0; i < num_nodes; i++)
    {
      all_nodes_pt[i]->object_pt()->set_number(++count);
    }

    // Check neighbours vertices and their opposite points
    std::ofstream neighbours_file;
    std::ofstream no_true_edge_file;
    std::ofstream neighbours_txt_file;

    double max_error_face = 0.0;
    OcTree::doc_face_neighbours(
      all_nodes_pt, neighbours_file, neighbours_txt_file, max_error_face);

    double max_error_edge = 0.0;
    OcTree::doc_true_edge_neighbours(all_nodes_pt,
                                     neighbours_file,
                                     no_true_edge_file,
                                     neighbours_txt_file,
                                     max_error_edge);
    bool failed = false;
    if (max_error_face > Max_neighbour_finding_tolerance)
    {
      oomph_info
        << "\n \n Failed self_test() for OcTree because of faces: Max. error "
        << max_error_face << std::endl
        << std::endl;
      failed = true;
    }

    if (max_error_edge > Max_neighbour_finding_tolerance)
    {
      oomph_info
        << "\n \n Failed self_test() for OcTree because of edges: Max. error "
        << max_error_edge << std::endl
        << std::endl;
      failed = true;
    }

    double max_error = max_error_face;
    if (max_error_edge > max_error) max_error = max_error_edge;

    if (failed)
    {
      return 1;
    }
    else
    {
      oomph_info << "Passed self_test() for OcTree: Max. error " << max_error
                 << std::endl;
      return 0;
    }
  }


  //=================================================================
  /// Doc/check all face neighbours of octree (nodes) contained in the
  /// Vector forest_node_pt. Output into neighbours_file which can
  /// be viewed from tecplot with OcTreeNeighbours.mcr
  /// Neighbour info and errors are displayed on
  /// neighbours_txt_file.  Finally, compute the max. error between
  /// vertices when viewed from neighbouring element.
  /// If the two filestreams are closed, output is suppressed.
  /// (Static function.)
  //=================================================================
  void OcTree::doc_face_neighbours(Vector<Tree*> forest_nodes_pt,
                                   std::ofstream& neighbours_file,
                                   std::ofstream& neighbours_txt_file,
                                   double& max_error)
  {
    using namespace OcTreeNames;

    int diff_level;
    int face = OMEGA;
    bool in_neighbouring_tree;

    Vector<double> s(3);
    Vector<double> x(3);

    Vector<double> s_sw(3);
    Vector<double> s_ne(3);
    Vector<unsigned> translate_s(3);

    Vector<double> x_small(3);
    Vector<double> x_large(3);


    // Initialise error in vertex positions
    max_error = 0.0;

    // Loop over all elements
    // ----------------------
    unsigned long num_nodes = forest_nodes_pt.size();

    for (unsigned long i = 0; i < num_nodes; i++)
    {
      // Doc the element itself
      OcTree* el_pt = dynamic_cast<OcTree*>(forest_nodes_pt[i]);

      // If the object is incomplete omit
      if (el_pt->object_pt()->nodes_built())
      {
        // Print it
        if (neighbours_file.is_open())
        {
          neighbours_file << "#---------------------------------" << std::endl;
          neighbours_file << "#The element itself: " << i << std::endl;
          neighbours_file << "#---------------------------------" << std::endl;
          dynamic_cast<RefineableQElement<3>*>(el_pt->object_pt())
            ->output_corners(neighbours_file, "BLACK");
        }

        // Loop over directions to find neighbours
        // ----------------------------------------
        for (int direction = L; direction <= F; direction++)
        {
          // Initialise difference in levels and coordinate offset
          diff_level = 0;

          // Find greater-or-equal-sized neighbour...
          OcTree* neighb_pt = el_pt->gteq_face_neighbour(direction,
                                                         translate_s,
                                                         s_sw,
                                                         s_ne,
                                                         face,
                                                         diff_level,
                                                         in_neighbouring_tree);

          // If neighbour exists and nodes are created
          if ((neighb_pt != 0) && (neighb_pt->object_pt()->nodes_built()))
          {
            // Doc neighbour stats
            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file
                << Direct_string[direction] << " neighbour of "
                << el_pt->object_pt()->number() << " is "
                << neighb_pt->object_pt()->number() << " diff_level "
                << diff_level << ". Inside the neighbour the face is "
                << Direct_string[face] << std::endl;
            }

            // Plot neighbour in the appropriate colour
            if (neighbours_file.is_open())
            {
              neighbours_file << "#---------------------------------"
                              << std::endl;
              neighbours_file
                << "#Neighbour element: " << Direct_string[direction]
                << "\n#---------------------------------" << std::endl;
              dynamic_cast<RefineableQElement<3>*>(neighb_pt->object_pt())
                ->output_corners(neighbours_file, Colour[direction]);
            }

            // Check that local coordinates in the larger element
            // lead to the same spatial point as the node vertices
            // in the current element
            if (neighbours_file.is_open())
            {
              neighbours_file << "ZONE I=2  C=" << Colour[direction]
                              << std::endl;
            }

            // "South west" vertex in the interfacial face
            //--------------------------------------------

            // Get coordinates in large (neighbour) element
            s[0] = s_sw[0];
            s[1] = s_sw[1];
            s[2] = s_sw[2];
            neighb_pt->object_pt()->get_x(s, x_large);

            // Get coordinates in small element
            Vector<double> s(3);
            s[0] = S_base(0, direction);
            s[1] = S_base(1, direction);
            s[2] = S_base(2, direction);
            el_pt->object_pt()->get_x(s, x_small);

            // Need to exclude periodic nodes from this check
            // There can only be periodic nodes if we have moved into the
            // neighbour
            bool is_periodic = false;
            if (in_neighbouring_tree)
            {
              // is the node periodic
              is_periodic = el_pt->root_pt()->is_neighbour_periodic(direction);
            }

            double error = 0.0;
            // Only bother to calculate the error if the node is NOT periodic
            if (is_periodic == false)
            {
              error = sqrt(pow(x_small[0] - x_large[0], 2) +
                           pow(x_small[1] - x_large[1], 2) +
                           pow(x_small[2] - x_large[2], 2));
            }

            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file << "Error (1) " << error << std::endl;
            }

            if (std::fabs(error) > max_error)
            {
              max_error = std::fabs(error);
            }

            // Check error and doc mismatch if required
            bool stop = false;
            std::ofstream mismatch_file;
            if (std::fabs(error) > Max_neighbour_finding_tolerance)
            {
              stop = true;
              mismatch_file.open("mismatch.dat");
              mismatch_file << "ZONE" << std::endl;
              mismatch_file << x_large[0] << " " << x_large[1] << " "
                            << x_large[2] << " 2 \n";
              mismatch_file << x_small[0] << " " << x_small[1] << " "
                            << x_small[2] << " 3 \n";
            }

            if (neighbours_file.is_open())
            {
              neighbours_file << "#SOUTH WEST: " << std::endl;
              neighbours_file << x_large[0] << " " << x_large[1] << " "
                              << x_large[2] << "  40 \n";
            }


            // "North east" vertex in the interfacial face
            //--------------------------------------------

            // Get coordinates in large (neighbour) element
            s[0] = s_ne[0];
            s[1] = s_ne[1];
            s[2] = s_ne[2];
            neighb_pt->object_pt()->get_x(s, x_large);

            // Get coordinates in small element
            s[0] = S_base(0, direction) + S_steplo(0, direction) +
                   S_stephi(0, direction);
            s[1] = S_base(1, direction) + S_steplo(1, direction) +
                   S_stephi(1, direction);
            s[2] = S_base(2, direction) + S_steplo(2, direction) +
                   S_stephi(2, direction);
            el_pt->object_pt()->get_x(s, x_small);

            error = 0.0;
            // Only do this if we are NOT periodic
            if (is_periodic == false)
            {
              error = sqrt(pow(x_small[0] - x_large[0], 2) +
                           pow(x_small[1] - x_large[1], 2) +
                           pow(x_small[2] - x_large[2], 2));
            }

            // Output
            if (neighbours_file.is_open())
            {
              neighbours_file << "#NORTH EAST: " << std::endl;
              neighbours_file << x_large[0] << " " << x_large[1] << " "
                              << x_large[2] << "  80 \n";
            }

            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file << "Error (2) " << error << std::endl;
            }

            if (std::fabs(error) > max_error)
            {
              max_error = std::fabs(error);
            }

            // Check error and doc mismatch if required
            if (std::fabs(error) > Max_neighbour_finding_tolerance)
            {
              stop = true;
            }
            if (stop)
            {
              if (!mismatch_file.is_open())
              {
                mismatch_file.open("mismatch.dat");
              }
              mismatch_file << "ZONE" << std::endl;
              mismatch_file << x_large[0] << " " << x_large[1] << " "
                            << x_large[2] << "  2 " << std::fabs(error)
                            << " \n";
              mismatch_file << x_small[0] << " " << x_small[1] << " "
                            << x_small[2] << "  3 " << std::fabs(error)
                            << " \n";
              mismatch_file.close();
              pause("Error");
            }
          }

          // If neighbour does not exist: Insert blank zones into file
          // so that tecplot can find six neighbours for every element
          else
          {
            if (neighbours_file.is_open())
            {
              neighbours_file
                << "#---------------------------------\n"
                << "# No neighbour in direction: " << Direct_string[direction]
                << "\n"
                << "#---------------------------------" << std::endl;

              dynamic_cast<RefineableQElement<3>*>(el_pt->object_pt())
                ->output_corners(neighbours_file, "WHITE");
              neighbours_file << "ZONE I=2 \n";
              neighbours_file << "-0.05 -0.05 -0.05  0 \n";
              neighbours_file << "-0.05 -0.05 -0.05  0 \n";
            }
          }

          if (neighbours_file.is_open())
          {
            neighbours_file << std::endl << std::endl << std::endl;
          }
        }
      } // End of case when element can be documented
    }
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //=================================================================
  /// Doc/check all true edge neighbours of octree (nodes) contained
  /// in the Vector forest_node_pt. Output into neighbours_file which can
  /// be viewed from tecplot with OcTreeNeighbours.mcr
  /// Neighbour info and errors are displayed on
  /// neighbours_txt_file.  Finally, compute the max. error between
  /// vertices when viewed from neighbouring element.
  /// If the two filestreams are closed, output is suppressed.
  /// (Static function).
  //=================================================================
  void OcTree::doc_true_edge_neighbours(Vector<Tree*> forest_nodes_pt,
                                        std::ofstream& neighbours_file,
                                        std::ofstream& no_true_edge_file,
                                        std::ofstream& neighbours_txt_file,
                                        double& max_error)
  {
    using namespace OcTreeNames;

    int diff_level;
    int edge = OMEGA;

    Vector<double> s(3);
    Vector<double> x(3);

    Vector<double> s_lo(3);
    Vector<double> s_hi(3);
    Vector<unsigned> translate_s(3);

    Vector<double> x_small(3);
    Vector<double> x_large(3);


    // Initialise error in vertex positions
    max_error = 0.0;

    // Loop over all elements
    // ----------------------
    unsigned long num_nodes = forest_nodes_pt.size();

    for (unsigned long i = 0; i < num_nodes; i++)
    {
      // Doc the element itself
      OcTree* el_pt = dynamic_cast<OcTree*>(forest_nodes_pt[i]);

      // If the object is incomplete omit
      if (el_pt->object_pt()->nodes_built())
      {
        // Print it
        if (neighbours_file.is_open())
        {
          neighbours_file << "#---------------------------------" << std::endl;
          neighbours_file << "# The element itself: " << i << std::endl;
          neighbours_file << "#---------------------------------" << std::endl;
          dynamic_cast<RefineableQElement<3>*>(el_pt->object_pt())
            ->output_corners(neighbours_file, "BLACK");
        }

        // Loop over directions to find edge neighbours
        // --------------------------------------------
        for (int direction = LB; direction <= UF; direction++)
        {
          // Initialise difference in levels
          diff_level = 0;

          // For now simply doc the zero-th edge neighbour (if any)
          unsigned i_root_edge_neighbour = 0;
          unsigned nroot_edge_neighbour = 0;

          // Find greater-or-equal-sized edge neighbour...
          OcTree* neighb_pt =
            el_pt->gteq_true_edge_neighbour(direction,
                                            i_root_edge_neighbour,
                                            nroot_edge_neighbour,
                                            translate_s,
                                            s_lo,
                                            s_hi,
                                            edge,
                                            diff_level);

          // If neighbour exists and nodes are created
          if ((neighb_pt != 0) && (neighb_pt->object_pt()->nodes_built()))
          {
            // Doc neighbour stats
            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file
                << Direct_string[direction] << " neighbour of "
                << el_pt->object_pt()->number() << " is "
                << neighb_pt->object_pt()->number() << " diff_level "
                << diff_level << ". Inside the neighbour the edge is "
                << Direct_string[edge] << std::endl;
            }

            // Plot neighbour in the appropriate colour
            if (neighbours_file.is_open())
            {
              neighbours_file
                << "#---------------------------------"
                << "\n# Neighbour element: " << Direct_string[direction]
                << "\n#---------------------------------" << std::endl;
              dynamic_cast<RefineableQElement<3>*>(neighb_pt->object_pt())
                ->output_corners(neighbours_file, Colour[direction]);
            }

            // Check that local coordinates in the larger element
            // lead to the same spatial point as the node vertices
            // in the current element
            if (neighbours_file.is_open())
            {
              neighbours_file << "ZONE I=2  C=" << Colour[direction]
                              << std::endl;
            }

            // "Low" vertex in the edge
            //-------------------------
            // Get coordinates in large (neighbour) element
            s[0] = s_lo[0];
            s[1] = s_lo[1];
            s[2] = s_lo[2];
            neighb_pt->object_pt()->get_x(s, x_large);

            // Get coordinates in small element
            Vector<double> s(3);
            s[0] = S_base_edge(0, direction);
            s[1] = S_base_edge(1, direction);
            s[2] = S_base_edge(2, direction);
            el_pt->object_pt()->get_x(s, x_small);

            // Need to exclude periodic nodes from this check
            // There can only be periodic nodes if we have moved into the
            // neighbour
            bool is_periodic = false;

            // Get the faces on which the edge lies
            Vector<int> faces_attached_to_edge =
              faces_of_common_edge(direction);

            // Get the number of entries in the vector
            unsigned n_faces_attached_to_edge = faces_attached_to_edge.size();

            // Loop over the faces
            for (unsigned i_face = 0; i_face < n_faces_attached_to_edge;
                 i_face++)
            {
              // Is the node periodic in the face direction?
              is_periodic = el_pt->root_pt()->is_neighbour_periodic(
                faces_attached_to_edge[i_face]);

              // Check if the edge is periodic in the i_face-th face direction
              if (is_periodic)
              {
                // We're done!
                break;
              }
            } // for (unsigned
              // i_face=0;i_face<n_faces_attached_to_edge;i_face++)

            double error = 0.0;
            // Only bother to calculate the error if the node is NOT periodic
            if (is_periodic == false)
            {
              error = sqrt(pow(x_small[0] - x_large[0], 2) +
                           pow(x_small[1] - x_large[1], 2) +
                           pow(x_small[2] - x_large[2], 2));
            }

            if (std::fabs(error) > max_error)
            {
              max_error = std::fabs(error);
            }

            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file << "Error (1) " << error << std::endl;
            }

            // Check error and doc mismatch if required
            bool stop = false;
            std::ofstream mismatch_file;
            if (std::fabs(error) > Max_neighbour_finding_tolerance)
            {
              stop = true;
              mismatch_file.open("mismatch.dat");
              mismatch_file << "ZONE" << std::endl;
              mismatch_file << x_large[0] << " " << x_large[1] << " "
                            << x_large[2] << " 2 \n";
              mismatch_file << x_small[0] << " " << x_small[1] << " "
                            << x_small[2] << " 3 \n";
            }

            if (neighbours_file.is_open())
            {
              neighbours_file << "# LOW VERTEX: " << std::endl;
              neighbours_file << x_large[0] << " " << x_large[1] << " "
                              << x_large[2] << "  40 \n";
            }


            // "High" vertex in the edge
            //--------------------------
            // Get coordinates in large (neighbour) element
            s[0] = s_hi[0];
            s[1] = s_hi[1];
            s[2] = s_hi[2];
            neighb_pt->object_pt()->get_x(s, x_large);

            // Get coordinates in small element
            s[0] = S_base_edge(0, direction) + S_step_edge(0, direction);
            s[1] = S_base_edge(1, direction) + S_step_edge(1, direction);
            s[2] = S_base_edge(2, direction) + S_step_edge(2, direction);
            el_pt->object_pt()->get_x(s, x_small);

            // Output
            if (neighbours_file.is_open())
            {
              neighbours_file << "# HI VERTEX: " << std::endl;
              neighbours_file << x_large[0] << " " << x_large[1] << " "
                              << x_large[2] << "  80 \n";
            }

            // Reset the error value
            error = 0.0;

            // Only do this if we are NOT periodic
            if (is_periodic == false)
            {
              error = sqrt(pow(x_small[0] - x_large[0], 2) +
                           pow(x_small[1] - x_large[1], 2) +
                           pow(x_small[2] - x_large[2], 2));
            }

            if (neighbours_txt_file.is_open())
            {
              neighbours_txt_file << "Error (2) " << error << std::endl;
            }

            if (std::fabs(error) > max_error)
            {
              max_error = std::fabs(error);
            }

            // Check error and doc mismatch if required
            if (std::fabs(error) > Max_neighbour_finding_tolerance)
            {
              stop = true;
            }
            if (stop)
            {
              if (!mismatch_file.is_open())
              {
                mismatch_file.open("mismatch.dat");
              }
              mismatch_file << "ZONE" << std::endl;
              mismatch_file << x_large[0] << " " << x_large[1] << " "
                            << x_large[2] << "  2 \n";
              mismatch_file << x_small[0] << " " << x_small[1] << " "
                            << x_small[2] << "  3 \n";
              mismatch_file.close();
              pause("Error");
            }
          }

          // If neighbour does not exist: Insert blank zones into file
          // so that tecplot can find twelve neighbours for every element
          else
          {
            if (neighbours_file.is_open())
            {
              neighbours_file << "#---------------------------------"
                              << std::endl;
              neighbours_file
                << "# No neighbour in direction: " << Direct_string[direction]
                << std::endl;
              neighbours_file << "#---------------------------------"
                              << std::endl;

              dynamic_cast<RefineableQElement<3>*>(el_pt->object_pt())
                ->output_corners(neighbours_file, "WHITE");
              neighbours_file << "ZONE I=2 \n";
              neighbours_file << "-0.05 -0.05 -0.05  0 \n";
              neighbours_file << "-0.05 -0.05 -0.05  0 \n";


              // Doc edge for which no neighbour exists but only
              // for the smallest elements
              if (el_pt->is_leaf())
              {
                // Get start coordinates in small element
                Vector<double> s(3), x_start(3), x_end(3);
                s[0] = S_base_edge(0, direction);
                s[1] = S_base_edge(1, direction);
                s[2] = S_base_edge(2, direction);
                el_pt->object_pt()->get_x(s, x_start);


                // Get coordinates in small element
                s[0] = S_base_edge(0, direction) + S_step_edge(0, direction);
                s[1] = S_base_edge(1, direction) + S_step_edge(1, direction);
                s[2] = S_base_edge(2, direction) + S_step_edge(2, direction);
                el_pt->object_pt()->get_x(s, x_end);

                no_true_edge_file << "ZONE I=2" << std::endl;
                no_true_edge_file << x_start[0] << " " << x_start[1] << " "
                                  << x_start[2] << " " << std::endl;
                no_true_edge_file << x_end[0] << " " << x_end[1] << " "
                                  << x_end[2] << " " << std::endl;
              }
            }

            if (neighbours_file.is_open())
            {
              neighbours_file << std::endl << std::endl << std::endl;
            }
          }
        } // End of case when element can be documented
      }
    }
  }


  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////


  //================================================================
  /// Constructor for OcTreeForest:
  ///
  /// Pass:
  ///  - trees_pt[], the Vector of pointers to the constituent trees
  ///    (OcTreeRoot objects)
  ///
  //=================================================================
  OcTreeForest::OcTreeForest(Vector<TreeRoot*>& trees_pt) : TreeForest(trees_pt)
  {
#ifdef LEAK_CHECK
    LeakCheckNames::OcTreeForest_build += 1;
#endif

    // Don't setup neighbours etc. if forest is empty
    if (trees_pt.size() == 0)
    {
      return;
    }

    using namespace OcTreeNames;

    // Construct the neighbour and rotation scheme, note that all neighbour
    // pointers must be set before the constructor is called

    //  MemoryUsage::doc_memory_usage(
    //   "before find_neighbours in octree forest constr");

    // setup the neighbour scheme
    find_neighbours();

    //  MemoryUsage::doc_memory_usage(
    //   "after find_neighbours in octree forest constr");

    // setup the rotation scheme
    construct_up_right_equivalents();

    //  MemoryUsage::doc_memory_usage(
    //   "after construct_up_right_equivalents in octree forest constr");
  }

  //================================================================
  /// setup the neighbour scheme : tells to each element in the
  /// forest which (if any) element is its {R,L,U,D,B,F} face neighbour
  /// and which is its {LB,RB,...,UF} edge neighbour.
  //================================================================
  void OcTreeForest::find_neighbours()
  {
    using namespace OcTreeNames;
    unsigned numtrees = ntree();
    int n = 0; // to store nnode1d
    if (numtrees > 0)
    {
      n = Trees_pt[0]->object_pt()->nnode_1d();
    }
    else
    {
      throw OomphLibError(
        "Trying to setup the neighbour scheme for an empty forest",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    // Number of vertex nodes: 8
    unsigned n_vertex_node = 8;

    // Find potentially connected trees by identifying
    // those whose associated elements share a common vertex node
    std::map<Node*, std::set<unsigned>> tree_assoc_with_vertex_node;

    // Loop over all trees
    for (unsigned i = 0; i < numtrees; i++)
    {
      // Loop over the vertex nodes of the associated element
      for (unsigned j = 0; j < n_vertex_node; j++)
      {
        Node* nod_pt = dynamic_cast<BrickElementBase*>(Trees_pt[i]->object_pt())
                         ->vertex_node_pt(j);
        tree_assoc_with_vertex_node[nod_pt].insert(i);
      }
    }


    // For each tree we store a set of potentially neighbouring trees
    // i.e. trees that share at least one node
    Vector<std::set<unsigned>> potentially_neighb_tree(numtrees);

    // Loop over vertex nodes
    for (std::map<Node*, std::set<unsigned>>::iterator it =
           tree_assoc_with_vertex_node.begin();
         it != tree_assoc_with_vertex_node.end();
         it++)
    {
      // Loop over connected elements twice
      for (std::set<unsigned>::iterator it_el1 = it->second.begin();
           it_el1 != it->second.end();
           it_el1++)
      {
        unsigned i = (*it_el1);
        for (std::set<unsigned>::iterator it_el2 = it->second.begin();
             it_el2 != it->second.end();
             it_el2++)
        {
          unsigned j = (*it_el2);
          // These two elements are potentially connected
          if (i != j)
          {
            potentially_neighb_tree[i].insert(j);
          }
        }
      }
    }


    // Loop over all trees
    for (unsigned i = 0; i < numtrees; i++)
    {
      // Cast to OcTreeRoot
      OcTreeRoot* octree_root_i_pt = dynamic_cast<OcTreeRoot*>(Trees_pt[i]);

      // Loop over their potential neighbours
      for (std::set<unsigned>::iterator it = potentially_neighb_tree[i].begin();
           it != potentially_neighb_tree[i].end();
           it++)
      {
        unsigned j = (*it);

        // So far, we haven't identified the candidate element as a face
        // neighbour of element/tree i
        bool is_face_neighbour = false;

        // is it the Up neighbour ?
        bool is_Up_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * (n * n - 1))) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n - 1)) != -1));

        // is it the Down neighbour ?
        bool is_Down_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n * (n - 1))) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1));

        // is it the Right neighbour ?
        bool is_Right_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n * n - 1)) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1));

        // is it the Left neighbour ?
        bool is_Left_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * (n * n - 1))) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(0)) != -1));

        // is it the Back neighbour ?
        bool is_Back_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n - 1)) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(0)) != -1));

        // is it the Front neighbour ?
        bool is_Front_neighbour =
          ((Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n * n - 1)) != -1) &&
           (Trees_pt[j]->object_pt()->get_node_number(
              Trees_pt[i]->object_pt()->node_pt(n * n * (n - 1))) != -1));


        if (is_Down_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(D) = Trees_pt[j];
        }
        if (is_Up_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(U) = Trees_pt[j];
        }
        if (is_Right_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(R) = Trees_pt[j];
        }
        if (is_Left_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(L) = Trees_pt[j];
        }
        if (is_Back_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(B) = Trees_pt[j];
        }
        if (is_Front_neighbour)
        {
          is_face_neighbour = true;
          Trees_pt[i]->neighbour_pt(F) = Trees_pt[j];
        }


        // If it's not a face neighbour, it may still be an edge
        // neighbour. We check this by checking if the
        // vertex nodes coincide. Note: This test would also
        // evaluate to true for face neighbours but we've already
        // determined that the element is not a face neighbour!
        if (!is_face_neighbour)
        {
          // is it the left back neighbour ?
          bool is_left_back_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(0)) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * (n - 1))) != -1));

          // is it the right back neighbour ?
          bool is_right_back_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n - 1)) != -1));


          // is it the down back neighbour ?
          bool is_down_back_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(0)) != -1));

          // is it the up back neighbour ?
          bool is_up_back_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * (n - 1))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n - 1)) != -1));


          // is it the left down neighbour ?
          bool is_left_down_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n * (n - 1))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(0)) != -1));


          // is it the right down neighbour ?
          bool is_right_down_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n + 1) * (n - 1))) !=
              -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n - 1)) != -1));


          // is it the left up neighbour ?
          bool is_left_up_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - n))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * (n - 1))) != -1));


          // is it the right up neighbour ?
          bool is_right_up_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - 1))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n - 1)) != -1));


          // is it the left front neighbour ?
          bool is_left_front_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - n))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n * (n - 1))) != -1));


          // is it the right front neighbour ?
          bool is_right_front_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - 1))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n + 1) * (n - 1))) !=
              -1));


          // is it the down front neighbour ?
          bool is_down_front_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt(n * n * (n - 1))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n + 1) * (n - 1))) !=
              -1));


          // is it the up front neighbour ?
          bool is_up_front_neighbour =
            ((Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - n))) != -1) &&
             (Trees_pt[j]->object_pt()->get_node_number(
                Trees_pt[i]->object_pt()->node_pt((n * n * n - 1))) != -1));


          // Add to storage scheme for edge neighbours (only!)

          if (is_left_back_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(LB)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], LB);
          }
          if (is_right_back_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(RB)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], RB);
          }
          if (is_down_back_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(DB)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], DB);
          }
          if (is_up_back_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(UB)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], UB);
          }


          if (is_left_down_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(LD)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], LD);
          }
          if (is_right_down_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(RD)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], RD);
          }
          if (is_left_up_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(LU)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], LU);
          }
          if (is_right_up_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(RU)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], RU);
          }


          if (is_left_front_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(LF)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], LF);
          }
          if (is_right_front_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(RF)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], RF);
          }
          if (is_down_front_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(DF)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], DF);
          }
          if (is_up_front_neighbour)
          {
            // Trees_pt[i]->neighbour_pt(UF)=Trees_pt[j];
            octree_root_i_pt->add_edge_neighbour_pt(Trees_pt[j], UF);
          }
        }
      }
    }
  }


  //================================================================
  /// Construct the rotation scheme for the octree forest.
  /// Note that all pointers to neighbours must have been allocated
  /// for this to work.
  //================================================================
  void OcTreeForest::construct_up_right_equivalents()
  {
    using namespace OcTreeNames;

    // Number of trees in forest
    unsigned numtrees = ntree();

    // nnode1d from first element (if it exists!)
    int nnode1d = 0;
    if (numtrees > 0)
    {
      nnode1d = Trees_pt[0]->object_pt()->nnode_1d();
    }

    OcTreeRoot* neigh_pt = 0;

    // Loop over all the trees
    //------------------------
    for (unsigned i = 0; i < numtrees; i++)
    {
      // Find the pointer to the Up neighbour
      //------------------------------------
      neigh_pt = oc_face_neigh_pt(i, U);

      // If there is a neighbour?
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If up neighbour has a pointer to this tree
        if (direction != Tree::OMEGA)
        {
          // The direction to the element in the neighbour
          // must be equivalent to the down direction in this element
          // Hence, the up equivalent is the reflection of that direction
          octree_pt(i)->set_up_equivalent(neigh_pt,
                                          OcTree::Reflect_face[direction]);

          // The right equivalent is the direction to the equivalent of
          // the right face in the present element, which will
          // be connected to the UR edge of the present element,
          // but will not be the face adjacent to the U boundary
          // i.e. the "direction" face.
          // We find the local node numbers, in the neighbour, of
          // the UR edge of the present element
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d * nnode1d -
                                               1));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d - 1));

          // Now get the other face connected to that edge in the
          // neighbour. It is the right equivalent
          octree_pt(i)->set_right_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        // If U neighbour does not have pointer to this tree, die
        {
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Up neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Find the pointer to the Right neighbour
      //---------------------------------------
      neigh_pt = oc_face_neigh_pt(i, R);

      // If there is a neighbour?
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If the neighbour has a pointer to this tree
        if (direction != Tree::OMEGA)
        {
          // The direction to the element in the neighbour
          // must be equivalent to the left direction in this element
          // Hence, the right equivalent is the reflection of that direction
          octree_pt(i)->set_right_equivalent(neigh_pt,
                                             OcTree::Reflect_face[direction]);

          // The up equivalent will be the direction to the equivalent of
          // the up face in the neighbour, which will be connected to the
          // UR edge of the present element, but will not be the face
          // adjacent to the R boundary, i.e. the "direction" face
          // We find the local node numbers, in the neighbour, of the
          // UR edge of the present element
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d * nnode1d -
                                               1));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d - 1));

          // Now get the other face connected to that edge in the
          // neighbour. It is the up equivalent
          octree_pt(i)->set_up_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        {
          // If R neighbour does not have pointer to this tree, die
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Right neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Find the pointer to the Down neighbour
      //--------------------------------------
      neigh_pt = oc_face_neigh_pt(i, D);

      // If there is a neighbour?
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If the neighbour has a pointer to this element
        if (direction != Tree::OMEGA)
        {
          // The direction to the element in the neighbour must be
          // equivalent to the up direction
          octree_pt(i)->set_up_equivalent(neigh_pt, direction);

          // The right equivalent is the direction to the equivalent of
          // the right face in the present element, which will be
          // connected to the DR edge of the present element, but will
          // not be the face adjacent to the D boundary,
          // i.e. the "direction" face.
          // We find the local node numbers, in the neighbour, of
          // the RD edge of the present element
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d - 1));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d * nnode1d + 1) *
                                               (nnode1d - 1)));

          // Now get the other face connected to that edge in the neighbour.
          // It is the right equivalent
          octree_pt(i)->set_right_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        {
          // If D neighbour does not have pointer to this tree, die
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Down neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Find the pointer to the Left neighbour
      //--------------------------------------
      neigh_pt = oc_face_neigh_pt(i, L);

      // If there is a neighbour
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If the neighbour has a pointer to the present element
        if (direction != Tree::OMEGA)
        {
          // The direction to the element in the neighbour is
          // must be equivalent to the right direction
          octree_pt(i)->set_right_equivalent(neigh_pt, direction);

          // The up equivalent is the direction to the equivalent of the
          // up face in the present element, which will be connected to
          // the UL edge of the present element, but will not
          // be the face adjacent to the L boundary, i.e. the "direction"
          // face.
          // We find the local node numbers, in the neighbour, of the UL
          // edge
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d - 1) * nnode1d));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d * nnode1d - 1) *
                                               nnode1d));

          // Now get the other face connected to that edge in the
          // neighbour.It is the up equivalent
          octree_pt(i)->set_up_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        {
          // If L neighbour does not have pointer to this tree, die
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Left neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Find the pointer to the Front neighbour
      //---------------------------------------
      neigh_pt = oc_face_neigh_pt(i, F);

      // If there is a neighbour?
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If the neighbour has a pointer to the present element
        if (direction != Tree::OMEGA)
        {
          // The direction to the up face will be given by the
          // other face connected to the UF edge in the neighbour
          // i.e. the face that is not given by direction
          // We obtain the local node numbers, in the neighbour,
          // of the UF edge in the present element
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d * nnode1d -
                                               1));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d * nnode1d - 1) *
                                               nnode1d));

          // We now get the other face connected to that edge
          // It is the up equivalent.
          octree_pt(i)->set_up_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));

          // The direction to the right face will be given by the
          // other face connected to the RF edge in the neighbour
          // i.e. the face that is not given by the direction
          // We get the local node numbers, in the neighbour,
          // of the RF edge in the present element
          nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d * nnode1d -
                                               1));
          nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d * nnode1d + 1) *
                                               (nnode1d - 1)));

          // We now get the other face connected to that edge
          // It is the right equivalent
          octree_pt(i)->set_right_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        {
          // If F neighbour does not have pointer to this tree, die
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Front neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // Find the pointer to the Back neighbour
      //--------------------------------------
      neigh_pt = oc_face_neigh_pt(i, B);

      // If there is a neighbour?
      if (neigh_pt != 0)
      {
        // Find the direction of the present tree, as viewed from the neighbour
        int direction = neigh_pt->direction_of_neighbour(octree_pt(i));

        // If the neighbour has a pointer to the present element
        if (direction != Tree::OMEGA)
        {
          // The direction to the up face will be given by the
          // other face connected to the UB edge of the present element
          // i.e. not the "direction" face
          // We obtain the local node numbers, in the neighbour,
          // of the UB edge in the present element
          int nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt((nnode1d - 1) * nnode1d));
          int nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d - 1));

          // We now get the other face connected to that edge
          // It is the up equivalent
          octree_pt(i)->set_up_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));

          // The direction to the right face will be given by
          // the other face connected to the RB edge of the present
          // element, i.e. not the direction face.
          // We obtain local node numbers, in the neighbour,
          // of the RB edge in the present element
          nod1 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d * nnode1d - 1));
          nod2 = neigh_pt->object_pt()->get_node_number(
            octree_pt(i)->object_pt()->node_pt(nnode1d - 1));

          // We now get the other face connected to that edge
          // It is the right equivalent
          octree_pt(i)->set_right_equivalent(
            neigh_pt,
            OcTree::get_the_other_face(nod1, nod2, nnode1d, direction));
        }
        else
        {
          // If B neighbour does not have pointer to this tree, die
          std::ostringstream error_stream;
          error_stream << "Tree " << i
                       << "'s Back neighbour has no neighbour pointer to Tree "
                       << i << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      // EDGE NEIGHBOURS
      //----------------

      // Loop over all edges:
      for (int edge = LB; edge <= UF; edge++)
      {
        // Get vector to pointers of edge neighbours
        Vector<TreeRoot*> edge_neigh_pt = oc_edge_neigh_pt(i, edge);

        // Loop over all edge neighbours
        unsigned n_neigh = edge_neigh_pt.size();
        for (unsigned e = 0; e < n_neigh; e++)
        {
          // If there is a neighbour
          if (edge_neigh_pt[e] != 0)
          {
            // Here are the two vertices at the two ends of the edge
            // in the present element
            int vertex = OcTree::Vertex_at_end_of_edge[edge][0];
            int other_vertex = OcTree::Vertex_at_end_of_edge[edge][1];

            // Get node numbers of the two vertices on edge
            unsigned nod1 = OcTree::vertex_to_node_number(vertex, nnode1d);
            unsigned nod2 =
              OcTree::vertex_to_node_number(other_vertex, nnode1d);


            // Here are the local node numbers (in the neighbouring element)
            // of the start/end nodes of present element's edge:
            unsigned neighb_nod1 =
              edge_neigh_pt[e]->object_pt()->get_node_number(
                octree_pt(i)->object_pt()->node_pt(nod1));

            unsigned neighb_nod2 =
              edge_neigh_pt[e]->object_pt()->get_node_number(
                octree_pt(i)->object_pt()->node_pt(nod2));

            // Convert to vertices
            int neighb_vertex =
              OcTree::node_number_to_vertex(neighb_nod1, nnode1d);
            int neighb_other_vertex =
              OcTree::node_number_to_vertex(neighb_nod2, nnode1d);

            // Up equivalent is stored first in the pair that's returned
            // by the lookup table
            octree_pt(i)->set_up_equivalent(
              edge_neigh_pt[e],
              OcTree::Up_and_right_equivalent_for_pairs_of_vertices
                [std::make_pair(
                   std::make_pair(neighb_vertex, vertex),
                   std::make_pair(neighb_other_vertex, other_vertex))]
                  .first);

            // Right equivalent is stored second in the pair that's returned
            // by the lookup table
            octree_pt(i)->set_right_equivalent(
              edge_neigh_pt[e],
              OcTree::Up_and_right_equivalent_for_pairs_of_vertices
                [std::make_pair(
                   std::make_pair(neighb_vertex, vertex),
                   std::make_pair(neighb_other_vertex, other_vertex))]
                  .second);
          }
        }
      }
    } // end of loop over all trees.
  }


  //================================================================
  /// Self test: Check neighbour finding routine
  //=================================================================
  unsigned OcTreeForest::self_test()
  {
    // Stick pointers to all nodes into Vector and number elements
    // in the process
    Vector<Tree*> all_nodes_pt;
    stick_all_tree_nodes_into_vector(all_nodes_pt);
    long int count = 0;
    unsigned long num_nodes = all_nodes_pt.size();
    for (unsigned long i = 0; i < num_nodes; i++)
    {
      all_nodes_pt[i]->object_pt()->set_number(++count);
    }

    // Check neighbours vertices and their opposite points
    std::ofstream neighbours_file;
    std::ofstream no_true_edge_file;
    std::ofstream neighbours_txt_file;

    double max_error_face = 0.0;
    OcTree::doc_face_neighbours(
      all_nodes_pt, neighbours_file, neighbours_txt_file, max_error_face);

    double max_error_edge = 0.0;
    OcTree::doc_true_edge_neighbours(all_nodes_pt,
                                     neighbours_file,
                                     no_true_edge_file,
                                     neighbours_txt_file,
                                     max_error_edge);

    bool failed = false;
    if (max_error_face > OcTree::max_neighbour_finding_tolerance())
    {
      oomph_info << "\n\n Failed self_test() for OcTreeForest because of "
                    "faces: Max. error "
                 << max_error_face << std::endl
                 << std::endl;
      failed = true;
    }

    if (max_error_edge > OcTree::max_neighbour_finding_tolerance())
    {
      oomph_info << "\n\n Failed self_test() for OcTreeForest because of "
                    "edges: Max. error "
                 << max_error_edge << std::endl
                 << std::endl;
      failed = true;
    }

    double max_error = max_error_face;
    if (max_error_edge > max_error) max_error = max_error_edge;

    if (failed)
    {
      return 1;
    }
    else
    {
      oomph_info << "\nPassed self_test() for OcTreeForest: Max. error "
                 << max_error << std::endl;
      return 0;
    }
  }


  //================================================================
  /// Document and check all the neighbours of all the nodes
  /// in the forest. DocInfo object specifies the output directory
  /// and file numbers for the various files. If
  /// \c doc_info.is_doc_enabled()=false
  /// no output is created.
  //================================================================
  void OcTreeForest::check_all_neighbours(DocInfo& doc_info)
  {
    // Create vector of all elements in the tree
    Vector<Tree*> all_tree_nodes_pt;
    this->stick_all_tree_nodes_into_vector(all_tree_nodes_pt);

    // Face neighbours
    //----------------
    {
      // Create storage for the files
      std::ofstream neigh_file;
      std::ofstream neigh_txt_file;
      // If we are documenting, then do so
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream fullname;
        fullname << doc_info.directory() << "/neighbours" << doc_info.number()
                 << ".dat";
        oomph_info << "opened " << fullname.str() << " to doc neighbours"
                   << std::endl;
        neigh_file.open(fullname.str().c_str());
        fullname.str("");
        fullname << doc_info.directory() << "/neighbours" << doc_info.number()
                 << ".txt";
        oomph_info << "opened " << fullname.str() << " to doc neighbours"
                   << std::endl;
        neigh_txt_file.open(fullname.str().c_str());
      }

      // Call the static member of OcTree function
      double max_error = 0.0;
      OcTree::doc_face_neighbours(
        all_tree_nodes_pt, neigh_file, neigh_txt_file, max_error);
      if (max_error > Tree::max_neighbour_finding_tolerance())
      {
        std::ostringstream error_stream;
        error_stream << "\nMax. error in octree neighbour finding: "
                     << max_error << " is too big" << std::endl;
        error_stream
          << "i.e. bigger than Tree::max_neighbour_finding_tolerance()="
          << Tree::max_neighbour_finding_tolerance() << std::endl;

        // Close the files if they were opened
        if (doc_info.is_doc_enabled())
        {
          neigh_file.close();
          neigh_txt_file.close();
        }

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "\nMax. error in octree neighbour finding: " << max_error
                   << " is OK" << std::endl;
        oomph_info
          << "i.e. less than OcTree::max_neighbour_finding_tolerance()="
          << OcTree::max_neighbour_finding_tolerance() << std::endl;
      }

      // Close the documentation files if they were opened
      if (doc_info.is_doc_enabled())
      {
        neigh_file.close();
        neigh_txt_file.close();
      }
    }

    // Edge neighbours
    //----------------
    {
      // Create storage for the files
      std::ofstream neigh_file;
      std::ofstream no_true_edge_file;
      std::ofstream neigh_txt_file;
      // If we are documenting, then do so
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream fullname;
        fullname << doc_info.directory() << "/edge_neighbours"
                 << doc_info.number() << ".dat";
        neigh_file.open(fullname.str().c_str());
        fullname.str("");
        fullname << doc_info.directory() << "/no_true_edge" << doc_info.number()
                 << ".dat";
        no_true_edge_file.open(fullname.str().c_str());
        fullname.str("");
        fullname << doc_info.directory() << "/edge_neighbours"
                 << doc_info.number() << ".txt";
        neigh_txt_file.open(fullname.str().c_str());
      }

      // Call the static member of OcTree function
      double max_error = 0.0;
      // Get the maximum error between edge neighbours
      OcTree::doc_true_edge_neighbours(all_tree_nodes_pt,
                                       neigh_file,
                                       no_true_edge_file,
                                       neigh_txt_file,
                                       max_error);
      if (max_error > Tree::max_neighbour_finding_tolerance())
      {
        std::ostringstream error_stream;
        error_stream << "Max. error in octree edge neighbour finding: "
                     << max_error << " is too big" << std::endl;
        error_stream
          << "i.e. bigger than Tree::max_neighbour_finding_tolerance()="
          << Tree::max_neighbour_finding_tolerance() << std::endl;

        // Close the files if they were opened
        if (doc_info.is_doc_enabled())
        {
          neigh_file.close();
          no_true_edge_file.close();
          neigh_txt_file.close();
        }

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        oomph_info << "Max. error in octree edge neighbour finding: "
                   << max_error << " is OK" << std::endl;
        oomph_info
          << "i.e. less than OcTree::max_neighbour_finding_tolerance()="
          << OcTree::max_neighbour_finding_tolerance() << std::endl;
      }

      // Close the documentation files if they were opened
      if (doc_info.is_doc_enabled())
      {
        neigh_file.close();
        no_true_edge_file.close();
        neigh_txt_file.close();
      }
    }
  }


  //================================================================
  /// Open output files that will stored any hanging nodes that are
  /// created in the mesh refinement process.
  /// ===============================================================
  void OcTreeForest::open_hanging_node_files(
    DocInfo& doc_info, Vector<std::ofstream*>& output_stream)
  {
    // In 3D, there will be six output files
    for (unsigned i = 0; i < 6; i++)
    {
      output_stream.push_back(new std::ofstream);
    }

    // If we are documenting the output, open the files
    if (doc_info.is_doc_enabled())
    {
      std::ostringstream fullname;
      fullname << doc_info.directory() << "/hang_nodes_u" << doc_info.number()
               << ".dat";
      output_stream[0]->open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/hang_nodes_d" << doc_info.number()
               << ".dat";
      output_stream[1]->open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/hang_nodes_l" << doc_info.number()
               << ".dat";
      output_stream[2]->open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/hang_nodes_r" << doc_info.number()
               << ".dat";
      output_stream[3]->open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/hang_nodes_b" << doc_info.number()
               << ".dat";
      output_stream[4]->open(fullname.str().c_str());
      fullname.str("");
      fullname << doc_info.directory() << "/hang_nodes_f" << doc_info.number()
               << ".dat";
      output_stream[5]->open(fullname.str().c_str());
    }
  }


} // namespace oomph

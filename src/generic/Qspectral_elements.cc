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
#include "Qspectral_elements.h"
#include "Qelement_face_coordinate_translation_schemes.h"


namespace oomph
{
  std::map<unsigned, Vector<double>> OneDLegendreShapeParam::z;


  /// ///////////////////////////////////////////////////////////////////////
  ///                   1D QLegendreElements
  /// ///////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Assign the static integral
  //=======================================================================
  template<unsigned NNODE_1D>
  GaussLobattoLegendre<1, NNODE_1D> QSpectralElement<1, NNODE_1D>::integral;

  //=======================================================================
  /// The output function for general 1D QSpectralElements
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<1, NNODE_1D>::output(std::ostream& outfile)
  {
    // Tecplot header info
    outfile << "ZONE I=" << NNODE_1D << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l = 0; l < NNODE_1D; l++)
    {
      // Loop over the dimensions and output the position
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << this->node_pt(l)->x(i) << " ";
      }
      // Find out how many data values at the node
      unsigned initial_nvalue = this->node_pt(l)->nvalue();
      // Lopp over the data and output whether pinned or not
      for (unsigned i = 0; i < initial_nvalue; i++)
      {
        outfile << this->node_pt(l)->is_pinned(i) << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }

  //=======================================================================
  /// The output function for nplot points in each coordinate direction
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<1, NNODE_1D>::output(std::ostream& outfile,
                                             const unsigned& nplot)
  {
    // Local variables
    Vector<double> s(1);
    // Shape functions
    Shape psi(NNODE_1D);

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Tecplot header info
    outfile << "ZONE I=" << nplot << std::endl;
    // Loop over element nodes
    for (unsigned l = 0; l < nplot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (nplot - 1);
      shape(s, psi);
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Output the x and y positions
        outfile << this->interpolated_x(s, i) << " ";
      }
      for (unsigned i = 0; i < NNODE_1D; i++)
      {
        outfile << psi(i) << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }

  //===========================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type QSpectralElement<0,1>).
  //===========================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<1, NNODE_1D>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    /*throw OomphLibError("Untested",
      OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);*/

    // Overload the nodal dimension by reading out the value from the node
    face_element_pt->set_nodal_dimension(this->node_pt(0)->ndim());

    // Set the pointer to the "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // If the bulk element is halo then the face element must be too
    // if (this->is_halo())
    {
      face_element_pt->set_halo(Non_halo_proc_ID);
    }
#endif

    // Resize the storage for the original number of values at the (one and
    // only) node of the face element.
    face_element_pt->nbulk_value_resize(1);

    // Resize the storage for the bulk node number corresponding to the (one
    // and only) node of the face element
    face_element_pt->bulk_node_number_resize(1);

    // Set the face index in the face element
    face_element_pt->face_index() = face_index;

    // Now set up the node pointer
    // The convention is that the "normal", should always point
    // out of the element
    switch (face_index)
    {
        // Bottom, normal sign is negative (coordinate points into element)
      case (-1):
        face_element_pt->node_pt(0) = this->node_pt(0);
        face_element_pt->bulk_node_number(0) = 0;
        face_element_pt->normal_sign() = -1;

        // Set the pointer to the function that determines the bulk coordinates
        // in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement1FaceToBulkCoordinates::face0;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement1BulkCoordinateDerivatives::faces0;

        // Set the number of values stored when the node is part of the "bulk"
        // element. The required_nvalue() must be used, rather than nvalue(),
        // because otherwise nodes on boundaries will be resized multiple
        // times. If you want any other behaviour, you MUST set nbulk_value()
        // manually after construction of your specific face element.
        face_element_pt->nbulk_value(0) = this->required_nvalue(0);
        break;

        // Top, normal sign is positive (coordinate points out of element)
      case (1):
        face_element_pt->node_pt(0) = this->node_pt(NNODE_1D - 1);
        face_element_pt->bulk_node_number(0) = NNODE_1D - 1;
        face_element_pt->normal_sign() = +1;


        // Set the pointer to the function that determines the bulk coordinates
        // in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement1FaceToBulkCoordinates::face1;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement1BulkCoordinateDerivatives::faces0;


        // Set the number of values stored when the node is part of the "bulk"
        // element.
        face_element_pt->nbulk_value(0) = this->required_nvalue(NNODE_1D - 1);
        break;

        // Other cases
      default:
        std::ostringstream error_message;
        error_message << "Face_index should only take "
                      << "the values +/-1, not " << face_index << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  /// ///////////////////////////////////////////////////////////////////////
  ///                   2D QLegendreElements
  /// ///////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Assign the static integral
  //=======================================================================
  template<unsigned NNODE_1D>
  GaussLobattoLegendre<2, NNODE_1D> QSpectralElement<2, NNODE_1D>::integral;

  //=======================================================================
  /// The output function for general 1D QSpectralElements
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<2, NNODE_1D>::output(std::ostream& outfile)
  {
    // Tecplot header info
    outfile << "ZONE I=" << NNODE_1D << ", J=" << NNODE_1D << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l2 = 0; l2 < NNODE_1D; l2++)
    {
      for (unsigned l1 = 0; l1 < NNODE_1D; l1++)
      {
        unsigned l = l2 * NNODE_1D + l1;

        // Loop over the dimensions and output the position
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << this->node_pt(l)->x(i) << " ";
        }
        // Find out how many data values at the node
        unsigned initial_nvalue = this->node_pt(l)->nvalue();
        // Loop over the data and output whether pinned or not
        for (unsigned i = 0; i < initial_nvalue; i++)
        {
          outfile << this->node_pt(l)->is_pinned(i) << " ";
        }
        outfile << std::endl;
      }
    }
    outfile << std::endl;
  }


  //=======================================================================
  /// The output function for n+plot points in each coordinate direction
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<2, NNODE_1D>::output(std::ostream& outfile,
                                             const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;

    // Find the dimension of the first node
    unsigned n_dim = this->nodal_dimension();

    // Loop over plot points
    for (unsigned l2 = 0; l2 < n_plot; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
      for (unsigned l1 = 0; l1 < n_plot; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

        // Output the coordinates
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }
        outfile << std::endl;
      }
    }
    outfile << std::endl;
  }


  //===========================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type QSpectralElement<0,1>).
  //===========================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<2, NNODE_1D>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    // Set the nodal dimension from the "bulk"
    face_element_pt->set_nodal_dimension(this->node_pt(0)->ndim());

    // Set the pointer to the "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // If the bulk element is halo then the face element must be too
    // if (this->is_halo())
    {
      face_element_pt->set_halo(Non_halo_proc_ID);
    }
#endif

    // Resize the storage for the original number of values at
    // NNODE_1D nodes of the FaceElement
    face_element_pt->nbulk_value_resize(NNODE_1D);

    // Resize the storage for the bulk node numbers corresponding
    // to the NNODE_1D nodes of the FaceElement
    face_element_pt->bulk_node_number_resize(NNODE_1D);

    // Set the face index in the face element
    // The faces are
    //                       +1    East
    //                       -1    West
    //                       +2    North
    //                       -2    South

    // Set the face index in the face element
    face_element_pt->face_index() = face_index;

    // Now set up the node pointers
    // =================================
    // The convention here is that interior_tangent X tangent X tangent
    // is the OUTWARD normal
    switch (face_index)
    {
      unsigned bulk_number;
        // West face, normal sign is positive
      case (-1):
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement2FaceToBulkCoordinates::face0;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement2BulkCoordinateDerivatives::faces0;

        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          bulk_number = i * NNODE_1D;
          face_element_pt->node_pt(i) = this->node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = 1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = this->required_nvalue(bulk_number);
        }
        break;
        // South face, normal sign is positive
      case (-2):
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement2FaceToBulkCoordinates::face1;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement2BulkCoordinateDerivatives::faces1;

        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          bulk_number = i;
          face_element_pt->node_pt(i) = this->node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = 1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = this->required_nvalue(bulk_number);
        }
        break;
        // East face, normal sign is negative
      case (1):
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement2FaceToBulkCoordinates::face2;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement2BulkCoordinateDerivatives::faces0;

        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          bulk_number = NNODE_1D * i + NNODE_1D - 1;
          face_element_pt->node_pt(i) = this->node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = -1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = this->required_nvalue(bulk_number);
        }
        break;
        // North face, normal sign is negative
      case (2):
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement2FaceToBulkCoordinates::face3;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement2BulkCoordinateDerivatives::faces1;

        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          bulk_number = NNODE_1D * (NNODE_1D - 1) + i;
          face_element_pt->node_pt(i) = this->node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = -1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = this->required_nvalue(bulk_number);
        }
        break;

        // Now cover the other cases
      default:
        std::ostringstream error_message;
        error_message
          << "Face index should only take the values +/- 1 or +/- 2,"
          << " not " << face_index << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  /// ///////////////////////////////////////////////////////////////////////
  ///                   3D QLegendreElements
  /// ///////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Assign the static integral
  //=======================================================================
  template<unsigned NNODE_1D>
  GaussLobattoLegendre<3, NNODE_1D> QSpectralElement<3, NNODE_1D>::integral;

  //=======================================================================
  /// The output function for general 1D QSpectralElements
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<3, NNODE_1D>::output(std::ostream& outfile)
  {
    // Tecplot header info
    outfile << "ZONE I=" << NNODE_1D << ", J=" << NNODE_1D << ", K=" << NNODE_1D
            << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l3 = 0; l3 < NNODE_1D; l3++)
    {
      for (unsigned l2 = 0; l2 < NNODE_1D; l2++)
      {
        for (unsigned l1 = 0; l1 < NNODE_1D; l1++)
        {
          unsigned l = l3 * NNODE_1D * NNODE_1D + l2 * NNODE_1D + l1;

          // Loop over the dimensions and output the position
          for (unsigned i = 0; i < n_dim; i++)
          {
            outfile << this->node_pt(l)->x(i) << " ";
          }
          // Find out how many data values at the node
          unsigned initial_nvalue = this->node_pt(l)->nvalue();
          // Loop over the data and output whether pinned or not
          for (unsigned i = 0; i < initial_nvalue; i++)
          {
            outfile << this->node_pt(l)->is_pinned(i) << " ";
          }
          outfile << std::endl;
        }
      }
    }
    outfile << std::endl;
  }

  //=======================================================================
  /// The output function for n_plot points in each coordinate direction
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<3, NNODE_1D>::output(std::ostream& outfile,
                                             const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(3);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << ", K=" << n_plot
            << std::endl;

    // Find the dimension of the first node
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l3 = 0; l3 < n_plot; l3++)
    {
      s[2] = -1.0 + l3 * 2.0 / (n_plot - 1);
      for (unsigned l2 = 0; l2 < n_plot; l2++)
      {
        s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
        for (unsigned l1 = 0; l1 < n_plot; l1++)
        {
          s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

          // Output the coordinates
          for (unsigned i = 0; i < n_dim; i++)
          {
            outfile << this->interpolated_x(s, i) << " ";
          }
          outfile << std::endl;
        }
      }
    }
    outfile << std::endl;
  }


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type QElement<3,NNODE_1D>).
  //=======================================================================
  template<unsigned NNODE_1D>
  void QSpectralElement<3, NNODE_1D>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    oomph_info << " WARNING UNTESTED CODE" << std::endl;

    // Set the nodal dimension from the "bulk"
    face_element_pt->set_nodal_dimension(this->node_pt(0)->ndim());

    // Set the pointer to the orginal "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // If the bulk element is halo then the face element must be too
    // if (this->is_halo())
    {
      face_element_pt->set_halo(Non_halo_proc_ID);
    }
#endif

    // Resize storage for the number of values originally stored
    // at the face element's NNODE_1D*NNODE_1D nodes.
    face_element_pt->nbulk_value_resize(NNODE_1D * NNODE_1D);

    // Set the face index in the element
    // The faces are
    // -3 : BACK  (OLD: Bottom
    // -2 : DOWN  (OLD: Front
    // -1 : LEFT  (OLD: Left Side
    //  1 : RIGHT (OLD: Right Side
    //  2 : UP    (OLD: Back
    //  3 : FRONT (OLD: Top

    face_element_pt->face_index() = face_index;

    // Now set up the node pointers and the normal vectors
    switch (face_index)
    {
        // BACK
        //-----
      case -3:
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face2;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces2;

        // Copy nodes
        for (unsigned i = 0; i < (NNODE_1D * NNODE_1D); i++)
        {
          face_element_pt->node_pt(i) = this->node_pt(i);
        }
        // Outer unit normal is negative of cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;

        break;

        // FRONT
        //------
      case 3:

        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face5;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces2;

        // Copy nodes
        for (unsigned i = 0; i < (NNODE_1D * NNODE_1D); i++)
        {
          face_element_pt->node_pt(i) =
            this->node_pt(i + (NNODE_1D * NNODE_1D) * (NNODE_1D - 1));
        }
        // Outer unit normal is cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = 1;


        break;

        // DOWN:
        //------
      case -2:

      {
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face1;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces1;

        // Copy nodes
        unsigned count = 0;
        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          for (unsigned j = 0; j < NNODE_1D; j++)
          {
            face_element_pt->node_pt(count) =
              this->node_pt(j + i * (NNODE_1D * NNODE_1D));
            count++;
          }
        }

        // Outer unit normal is cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = 1;
      }
      break;


      // UP:
      //----
      case 2:

      {
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face4;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces1;

        // Copy nodes
        unsigned count = 0;
        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          for (unsigned j = 0; j < NNODE_1D; j++)
          {
            face_element_pt->node_pt(count) = this->node_pt(
              j + i * (NNODE_1D * NNODE_1D) + (NNODE_1D * (NNODE_1D - 1)));
            count++;
          }
        }

        // Outer unit normal is negative of cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;
      }
      break;

        // LEFT:
        //------
      case -1:

      {
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face0;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces0;

        // Copy nodes
        unsigned count = 0;
        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          for (unsigned j = 0; j < NNODE_1D; j++)
          {
            unsigned jj = j * NNODE_1D + i * (NNODE_1D * NNODE_1D);
            face_element_pt->node_pt(count) = this->node_pt(jj);
            count++;
          }
        }

        // Outer unit normal is negative of cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;
      }
      break;


      // RIGHT:
      //-------
      case 1:

      {
        // Set the pointer to the bulk coordinate translation scheme
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &QElement3FaceToBulkCoordinates::face3;

        // Set the pointer to the derivative mappings
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &QElement3BulkCoordinateDerivatives::faces0;

        // Copy nodes
        unsigned count = 0;
        for (unsigned i = 0; i < NNODE_1D; i++)
        {
          for (unsigned j = 0; j < NNODE_1D; j++)
          {
            unsigned jj =
              j * NNODE_1D + i * (NNODE_1D * NNODE_1D) + (NNODE_1D - 1);
            face_element_pt->node_pt(count) = this->node_pt(jj);
            count++;
          }
        }

        // Outer unit normal is cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = 1;
      }
      break;


      // Cover all other cases
      default:
        std::ostringstream error_message;
        error_message
          << "Face index should only take the values +/- 1, +/- 2 or +/- 3,"
          << " not " << face_index << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    } // end switch
  }


  template class QSpectralElement<1, 2>;
  template class QSpectralElement<1, 3>;
  template class QSpectralElement<1, 4>;
  template class QSpectralElement<1, 5>;
  template class QSpectralElement<1, 6>;
  template class QSpectralElement<1, 7>;
  template class QSpectralElement<1, 8>;
  template class QSpectralElement<1, 9>;
  template class QSpectralElement<1, 10>;
  template class QSpectralElement<1, 11>;
  template class QSpectralElement<1, 12>;
  template class QSpectralElement<1, 13>;
  template class QSpectralElement<1, 14>;

  template class QSpectralElement<2, 2>;
  template class QSpectralElement<2, 3>;
  template class QSpectralElement<2, 4>;
  template class QSpectralElement<2, 5>;
  template class QSpectralElement<2, 6>;
  template class QSpectralElement<2, 7>;
  template class QSpectralElement<2, 8>;

  template class QSpectralElement<3, 2>;
  template class QSpectralElement<3, 3>;
  template class QSpectralElement<3, 4>;
  template class QSpectralElement<3, 5>;
  template class QSpectralElement<3, 6>;
  template class QSpectralElement<3, 7>;
  template class QSpectralElement<3, 8>;

} // namespace oomph

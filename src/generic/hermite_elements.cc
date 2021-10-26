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
// Non-inline functions for the QHermiteElement classes

// oomph-lib header
#include "hermite_elements.h"
#include "Qelement_face_coordinate_translation_schemes.h"

namespace oomph
{
  /// /////////////////////////////////////////////////////////////
  //       1D Hermite elements
  /// /////////////////////////////////////////////////////////////

  //=======================================================================
  /// Assign the static Default_integration_scheme
  //=======================================================================
  template<unsigned DIM>
  Gauss<DIM, 3> QHermiteElement<DIM>::Default_integration_scheme;

  //=======================================================================
  /// Shape function for specific QHermiteElement<1>
  //=======================================================================
  template<>
  void QHermiteElement<1>::shape(const Vector<double>& s, Shape& psi) const
  {
    // Local storage
    double Psi[2][2];
    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi);

    // Loop over the number of nodes
    for (unsigned l = 0; l < 2; l++)
    {
      // Loop over the number of dofs
      for (unsigned k = 0; k < 2; k++)
      {
        psi(l, k) = Psi[l][k];
      }
    }
  }

  //=======================================================================
  /// Derivatives of shape functions for specific  QHermiteElement<1>
  //=======================================================================
  template<>
  void QHermiteElement<1>::dshape_local(const Vector<double>& s,
                                        Shape& psi,
                                        DShape& dpsids) const
  {
    // Local storage
    double Psi[2][2], DPsi[2][2];
    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi);
    OneDimHermite::dshape(s[0], DPsi);

    // Loop over number of nodes
    for (unsigned l = 0; l < 2; l++)
    {
      // Loop over number of dofs
      for (unsigned k = 0; k < 2; k++)
      {
        psi(l, k) = Psi[l][k];
        dpsids(l, k, 0) = DPsi[l][k];
      }
    }
  }

  //=======================================================================
  /// Derivatives and second derivatives of shape functions for specific
  /// QHermiteElement<1>
  /// d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
  //=======================================================================
  template<>
  void QHermiteElement<1>::d2shape_local(const Vector<double>& s,
                                         Shape& psi,
                                         DShape& dpsids,
                                         DShape& d2psids) const
  {
    // Local storage
    double Psi[2][2], DPsi[2][2], D2Psi[2][2];
    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi);
    OneDimHermite::dshape(s[0], DPsi);
    OneDimHermite::d2shape(s[0], D2Psi);

    // Loop over number of nodes
    for (unsigned l = 0; l < 2; l++)
    {
      // Loop over number of dofs
      for (unsigned k = 0; k < 2; k++)
      {
        psi(l, k) = Psi[l][k];
        dpsids(l, k, 0) = DPsi[l][k];
        d2psids(l, k, 0) = D2Psi[l][k];
      }
    }
  }


  //=======================================================================
  /// The output function for general 1D QHermiteElements
  //=======================================================================
  template<>
  void QHermiteElement<1>::output(std::ostream& outfile)
  {
    // Tecplot header info
    outfile << "ZONE I=" << 2 << std::endl;

    // Find the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l = 0; l < 2; l++)
    {
      // Loop over the dimensions and output the position
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << node_pt(l)->x(i) << " ";
      }

      // Find the number of types of dof stored at each node
      unsigned n_position_type = node_pt(l)->nposition_type();
      // Loop over the additional positional dofs
      for (unsigned k = 1; k < n_position_type; k++)
      {
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << node_pt(l)->x_gen(k, i) << " ";
        }
      }

      // Find out how many data values at the node
      unsigned initial_nvalue = node_pt(l)->nvalue();
      // Lopp over the data and output whether pinned or not
      for (unsigned i = 0; i < initial_nvalue; i++)
      {
        outfile << node_pt(l)->is_pinned(i) << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }

  //=======================================================================
  /// The output function for n_plot points in each coordinate direction
  //=======================================================================
  template<>
  void QHermiteElement<1>::output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Find the dimension of the first node
    unsigned n_dim = this->nodal_dimension();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Output the coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }


  //=======================================================================
  /// The C-style output function for general 1D QHermiteElements
  //=======================================================================
  template<>
  void QHermiteElement<1>::output(FILE* file_pt)
  {
    // Tecplot header info
    fprintf(file_pt, "ZONE I=2\n");

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l = 0; l < 2; l++)
    {
      // Loop over the dimensions and output the position
      for (unsigned i = 0; i < n_dim; i++)
      {
        fprintf(file_pt, "%g ", node_pt(l)->x(i));
      }

      // Find the number of types of dof stored at each node
      unsigned n_position_type = node_pt(l)->nposition_type();
      // Loop over the additional positional dofs
      for (unsigned k = 1; k < n_position_type; k++)
      {
        for (unsigned i = 0; i < n_dim; i++)
        {
          fprintf(file_pt, "%g ", node_pt(l)->x_gen(k, i));
        }
      }

      // Find out how many data values at the node
      unsigned initial_nvalue = node_pt(l)->nvalue();
      // Lopp over the data and output whether pinned or not
      for (unsigned i = 0; i < initial_nvalue; i++)
      {
        fprintf(file_pt, "%i ", node_pt(l)->is_pinned(i));
      }
      fprintf(file_pt, "\n");
    }
    fprintf(file_pt, "\n");
  }


  //=======================================================================
  /// The C-style output function for n_plot points in each coordinate direction
  //=======================================================================
  template<>
  void QHermiteElement<1>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    fprintf(file_pt, "ZONE I=%i \n", n_plot);

    // Find the dimension of the first node
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);
      // Output the coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      fprintf(file_pt, "\n");
    }
  }


  //===========================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (single node elements)
  //===========================================================
  template<>
  void QHermiteElement<1>::build_face_element(const int& face_index,
                                              FaceElement* face_element_pt)
  {
    // Set the nodal dimension from the "bulk"
    face_element_pt->set_nodal_dimension(node_pt(0)->ndim());

    // Set the pointer to the "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // Pass on non-halo ID
    face_element_pt->set_halo(Non_halo_proc_ID);
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
    // The convention is that the "normal", the coordinate, should always point
    // out of the element
    switch (face_index)
    {
        // Bottom, normal sign is negative (coordinate points into element)
      case (-1):
        face_element_pt->node_pt(0) = node_pt(0);
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
        // element.
        face_element_pt->nbulk_value(0) = required_nvalue(0);
        break;

        // Top, normal sign is positive (coordinate points out of element)
      case (1):
        face_element_pt->node_pt(0) = node_pt(1);
        face_element_pt->bulk_node_number(0) = 1;
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
        face_element_pt->nbulk_value(0) = required_nvalue(1);
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

  /// /////////////////////////////////////////////////////////////
  //       2D Hermite elements
  /// /////////////////////////////////////////////////////////////


  //=======================================================================
  /// Shape function for specific QHermiteElement<2>
  //=======================================================================
  template<>
  void QHermiteElement<2>::shape(const Vector<double>& s, Shape& psi) const
  {
    // Local storage
    double Psi[2][2][2];

    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi[0]);
    OneDimHermite::shape(s[1], Psi[1]);

    // Set up the two dimensional shape functions
    // Set up the functions at corner 0
    // psi_0 = 1 when s1 = 0, s2 = 0
    psi(0, 0) = Psi[0][0][0] * Psi[1][0][0];
    // dpsi_0/ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 1) = Psi[0][0][1] * Psi[1][0][0];
    // dpsi_0/ds2 = 1 when s1 = 0, s2 = 0
    psi(0, 2) = Psi[0][0][0] * Psi[1][0][1];
    // dpsi_0/ds2ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 3) = Psi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    // psi_1 = 1 when s1 = 1, s2 = 0
    psi(1, 0) = Psi[0][1][0] * Psi[1][0][0];
    // dpsi_1/ds1 = 1 when s1 = 1, s2 = 0
    psi(1, 1) = Psi[0][1][1] * Psi[1][0][0];
    // dpsi_1/ds2 = 1 when s1 = 1, s2 = 0
    psi(1, 2) = Psi[0][1][0] * Psi[1][0][1];
    // dpsi_1/ds1ds2  = 1 when s1 = 1, s2 = 0
    psi(1, 3) = Psi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    // psi_2 = 1 when s1 = 0, s2 = 1
    psi(2, 0) = Psi[0][0][0] * Psi[1][1][0];
    // dpsi_2/ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 1) = Psi[0][0][1] * Psi[1][1][0];
    // dpsi_2/ds2 = 1 when s1 = 0, s2 = 1
    psi(2, 2) = Psi[0][0][0] * Psi[1][1][1];
    // dpsi_2/ds2ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 3) = Psi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    // psi_3 = 1 when s1 = 1, s2 = 1
    psi(3, 0) = Psi[0][1][0] * Psi[1][1][0];
    // dpsi_3/ds1 = 1 when s1 = 1, s2 = 1
    psi(3, 1) = Psi[0][1][1] * Psi[1][1][0];
    // dpsi_3/ds2 = 1 when s1 = 1, s2 = 1
    psi(3, 2) = Psi[0][1][0] * Psi[1][1][1];
    // dpsi_3/ds1ds2  = 1 when s1 = 1, s2 = 1
    psi(3, 3) = Psi[0][1][1] * Psi[1][1][1];
  }

  //=======================================================================
  /// Derivatives of shape functions for specific QHermiteElement<2>
  //=======================================================================
  template<>
  void QHermiteElement<2>::dshape_local(const Vector<double>& s,
                                        Shape& psi,
                                        DShape& dpsids) const
  {
    // Local storage
    double Psi[2][2][2];
    double DPsi[2][2][2];
    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi[0]);
    OneDimHermite::shape(s[1], Psi[1]);
    OneDimHermite::dshape(s[0], DPsi[0]);
    OneDimHermite::dshape(s[1], DPsi[1]);


    // Set up the two dimensional shape functions
    // Set up the functions at corner 0
    // psi_0 = 1 when s1 = 0, s2 = 0
    psi(0, 0) = Psi[0][0][0] * Psi[1][0][0];
    // dpsi_0/ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 1) = Psi[0][0][1] * Psi[1][0][0];
    // dpsi_0/ds2 = 1 when s1 = 0, s2 = 0
    psi(0, 2) = Psi[0][0][0] * Psi[1][0][1];
    // dpsi_0/ds2ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 3) = Psi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    // psi_1 = 1 when s1 = 1, s2 = 0
    psi(1, 0) = Psi[0][1][0] * Psi[1][0][0];
    // dpsi_1/ds1 = 1 when s1 = 1, s2 = 0
    psi(1, 1) = Psi[0][1][1] * Psi[1][0][0];
    // dpsi_1/ds2 = 1 when s1 = 1, s2 = 0
    psi(1, 2) = Psi[0][1][0] * Psi[1][0][1];
    // dpsi_1/ds1ds2  = 1 when s1 = 1, s2 = 0
    psi(1, 3) = Psi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    // psi_2 = 1 when s1 = 0, s2 = 1
    psi(2, 0) = Psi[0][0][0] * Psi[1][1][0];
    // dpsi_2/ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 1) = Psi[0][0][1] * Psi[1][1][0];
    // dpsi_2/ds2 = 1 when s1 = 0, s2 = 1
    psi(2, 2) = Psi[0][0][0] * Psi[1][1][1];
    // dpsi_2/ds2ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 3) = Psi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    // psi_3 = 1 when s1 = 1, s2 = 1
    psi(3, 0) = Psi[0][1][0] * Psi[1][1][0];
    // dpsi_3/ds1 = 1 when s1 = 1, s2 = 1
    psi(3, 1) = Psi[0][1][1] * Psi[1][1][0];
    // dpsi_3/ds2 = 1 when s1 = 1, s2 = 1
    psi(3, 2) = Psi[0][1][0] * Psi[1][1][1];
    // dpsi_3/ds1ds2  = 1 when s1 = 1, s2 = 1
    psi(3, 3) = Psi[0][1][1] * Psi[1][1][1];

    // FIRST DERIVATIVES

    // D/Ds[0]

    // Set up the functions at corner 0
    dpsids(0, 0, 0) = DPsi[0][0][0] * Psi[1][0][0];
    dpsids(0, 1, 0) = DPsi[0][0][1] * Psi[1][0][0];
    dpsids(0, 2, 0) = DPsi[0][0][0] * Psi[1][0][1];
    dpsids(0, 3, 0) = DPsi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    dpsids(1, 0, 0) = DPsi[0][1][0] * Psi[1][0][0];
    dpsids(1, 1, 0) = DPsi[0][1][1] * Psi[1][0][0];
    dpsids(1, 2, 0) = DPsi[0][1][0] * Psi[1][0][1];
    dpsids(1, 3, 0) = DPsi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    dpsids(2, 0, 0) = DPsi[0][0][0] * Psi[1][1][0];
    dpsids(2, 1, 0) = DPsi[0][0][1] * Psi[1][1][0];
    dpsids(2, 2, 0) = DPsi[0][0][0] * Psi[1][1][1];
    dpsids(2, 3, 0) = DPsi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    dpsids(3, 0, 0) = DPsi[0][1][0] * Psi[1][1][0];
    dpsids(3, 1, 0) = DPsi[0][1][1] * Psi[1][1][0];
    dpsids(3, 2, 0) = DPsi[0][1][0] * Psi[1][1][1];
    dpsids(3, 3, 0) = DPsi[0][1][1] * Psi[1][1][1];

    // D/Ds[1]

    // Set up the functions at corner 0
    dpsids(0, 0, 1) = Psi[0][0][0] * DPsi[1][0][0];
    dpsids(0, 1, 1) = Psi[0][0][1] * DPsi[1][0][0];
    dpsids(0, 2, 1) = Psi[0][0][0] * DPsi[1][0][1];
    dpsids(0, 3, 1) = Psi[0][0][1] * DPsi[1][0][1];

    // Set up the functions at corner 1
    dpsids(1, 0, 1) = Psi[0][1][0] * DPsi[1][0][0];
    dpsids(1, 1, 1) = Psi[0][1][1] * DPsi[1][0][0];
    dpsids(1, 2, 1) = Psi[0][1][0] * DPsi[1][0][1];
    dpsids(1, 3, 1) = Psi[0][1][1] * DPsi[1][0][1];

    // Set up the functions at the corner 2
    dpsids(2, 0, 1) = Psi[0][0][0] * DPsi[1][1][0];
    dpsids(2, 1, 1) = Psi[0][0][1] * DPsi[1][1][0];
    dpsids(2, 2, 1) = Psi[0][0][0] * DPsi[1][1][1];
    dpsids(2, 3, 1) = Psi[0][0][1] * DPsi[1][1][1];

    // Set up the functions at corner 3
    dpsids(3, 0, 1) = Psi[0][1][0] * DPsi[1][1][0];
    dpsids(3, 1, 1) = Psi[0][1][1] * DPsi[1][1][0];
    dpsids(3, 2, 1) = Psi[0][1][0] * DPsi[1][1][1];
    dpsids(3, 3, 1) = Psi[0][1][1] * DPsi[1][1][1];
  }

  //======================================================================
  /// Second derivatives of the shape functions wrt local coordinates.
  /// d2psids(i,0) = \f$ \partial^2 \psi_j / \partial s_0^2 \f$
  /// d2psids(i,1) = \f$ \partial^2 \psi_j / \partial s_1^2 \f$
  /// d2psids(i,2) = \f$ \partial^2 \psi_j / \partial s_0 \partial s_1 \f$
  //======================================================================
  template<>
  void QHermiteElement<2>::d2shape_local(const Vector<double>& s,
                                         Shape& psi,
                                         DShape& dpsids,
                                         DShape& d2psids) const
  {
    // Local storage
    double Psi[2][2][2];
    double DPsi[2][2][2];
    double D2Psi[2][2][2];

    // Call the OneDimensional Shape functions
    OneDimHermite::shape(s[0], Psi[0]);
    OneDimHermite::shape(s[1], Psi[1]);
    OneDimHermite::dshape(s[0], DPsi[0]);
    OneDimHermite::dshape(s[1], DPsi[1]);
    OneDimHermite::d2shape(s[0], D2Psi[0]);
    OneDimHermite::d2shape(s[1], D2Psi[1]);

    // Set up the two dimensional shape functions
    // Set up the functions at corner 0
    // psi_0 = 1 when s1 = 0, s2 = 0
    psi(0, 0) = Psi[0][0][0] * Psi[1][0][0];
    // dpsi_0/ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 1) = Psi[0][0][1] * Psi[1][0][0];
    // dpsi_0/ds2 = 1 when s1 = 0, s2 = 0
    psi(0, 2) = Psi[0][0][0] * Psi[1][0][1];
    // dpsi_0/ds2ds1 = 1 when s1 = 0, s2 = 0
    psi(0, 3) = Psi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    // psi_1 = 1 when s1 = 1, s2 = 0
    psi(1, 0) = Psi[0][1][0] * Psi[1][0][0];
    // dpsi_1/ds1 = 1 when s1 = 1, s2 = 0
    psi(1, 1) = Psi[0][1][1] * Psi[1][0][0];
    // dpsi_1/ds2 = 1 when s1 = 1, s2 = 0
    psi(1, 2) = Psi[0][1][0] * Psi[1][0][1];
    // dpsi_1/ds1ds2  = 1 when s1 = 1, s2 = 0
    psi(1, 3) = Psi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    // psi_2 = 1 when s1 = 0, s2 = 1
    psi(2, 0) = Psi[0][0][0] * Psi[1][1][0];
    // dpsi_2/ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 1) = Psi[0][0][1] * Psi[1][1][0];
    // dpsi_2/ds2 = 1 when s1 = 0, s2 = 1
    psi(2, 2) = Psi[0][0][0] * Psi[1][1][1];
    // dpsi_2/ds2ds1 = 1 when s1 = 0, s2 = 1
    psi(2, 3) = Psi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    // psi_3 = 1 when s1 = 1, s2 = 1
    psi(3, 0) = Psi[0][1][0] * Psi[1][1][0];
    // dpsi_3/ds1 = 1 when s1 = 1, s2 = 1
    psi(3, 1) = Psi[0][1][1] * Psi[1][1][0];
    // dpsi_3/ds2 = 1 when s1 = 1, s2 = 1
    psi(3, 2) = Psi[0][1][0] * Psi[1][1][1];
    // dpsi_3/ds1ds2  = 1 when s1 = 1, s2 = 1
    psi(3, 3) = Psi[0][1][1] * Psi[1][1][1];

    // FIRST DERIVATIVES

    // D/Ds[0]

    // Set up the functions at corner 0
    dpsids(0, 0, 0) = DPsi[0][0][0] * Psi[1][0][0];
    dpsids(0, 1, 0) = DPsi[0][0][1] * Psi[1][0][0];
    dpsids(0, 2, 0) = DPsi[0][0][0] * Psi[1][0][1];
    dpsids(0, 3, 0) = DPsi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    dpsids(1, 0, 0) = DPsi[0][1][0] * Psi[1][0][0];
    dpsids(1, 1, 0) = DPsi[0][1][1] * Psi[1][0][0];
    dpsids(1, 2, 0) = DPsi[0][1][0] * Psi[1][0][1];
    dpsids(1, 3, 0) = DPsi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    dpsids(2, 0, 0) = DPsi[0][0][0] * Psi[1][1][0];
    dpsids(2, 1, 0) = DPsi[0][0][1] * Psi[1][1][0];
    dpsids(2, 2, 0) = DPsi[0][0][0] * Psi[1][1][1];
    dpsids(2, 3, 0) = DPsi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    dpsids(3, 0, 0) = DPsi[0][1][0] * Psi[1][1][0];
    dpsids(3, 1, 0) = DPsi[0][1][1] * Psi[1][1][0];
    dpsids(3, 2, 0) = DPsi[0][1][0] * Psi[1][1][1];
    dpsids(3, 3, 0) = DPsi[0][1][1] * Psi[1][1][1];

    // D/Ds[1]

    // Set up the functions at corner 0
    dpsids(0, 0, 1) = Psi[0][0][0] * DPsi[1][0][0];
    dpsids(0, 1, 1) = Psi[0][0][1] * DPsi[1][0][0];
    dpsids(0, 2, 1) = Psi[0][0][0] * DPsi[1][0][1];
    dpsids(0, 3, 1) = Psi[0][0][1] * DPsi[1][0][1];

    // Set up the functions at corner 1
    dpsids(1, 0, 1) = Psi[0][1][0] * DPsi[1][0][0];
    dpsids(1, 1, 1) = Psi[0][1][1] * DPsi[1][0][0];
    dpsids(1, 2, 1) = Psi[0][1][0] * DPsi[1][0][1];
    dpsids(1, 3, 1) = Psi[0][1][1] * DPsi[1][0][1];

    // Set up the functions at the corner 2
    dpsids(2, 0, 1) = Psi[0][0][0] * DPsi[1][1][0];
    dpsids(2, 1, 1) = Psi[0][0][1] * DPsi[1][1][0];
    dpsids(2, 2, 1) = Psi[0][0][0] * DPsi[1][1][1];
    dpsids(2, 3, 1) = Psi[0][0][1] * DPsi[1][1][1];

    // Set up the functions at corner 3
    dpsids(3, 0, 1) = Psi[0][1][0] * DPsi[1][1][0];
    dpsids(3, 1, 1) = Psi[0][1][1] * DPsi[1][1][0];
    dpsids(3, 2, 1) = Psi[0][1][0] * DPsi[1][1][1];
    dpsids(3, 3, 1) = Psi[0][1][1] * DPsi[1][1][1];

    // SECOND DERIVATIVES
    // Convention: index 0 is d^2/ds[0]^2,
    //            index 1 is d^2/ds[1]^2,
    //            index 2 is the mixed derivative

    // D^2/Ds[0]^2

    // Set up the functions at corner 0
    d2psids(0, 0, 0) = D2Psi[0][0][0] * Psi[1][0][0];
    d2psids(0, 1, 0) = D2Psi[0][0][1] * Psi[1][0][0];
    d2psids(0, 2, 0) = D2Psi[0][0][0] * Psi[1][0][1];
    d2psids(0, 3, 0) = D2Psi[0][0][1] * Psi[1][0][1];

    // Set up the functions at corner 1
    d2psids(1, 0, 0) = D2Psi[0][1][0] * Psi[1][0][0];
    d2psids(1, 1, 0) = D2Psi[0][1][1] * Psi[1][0][0];
    d2psids(1, 2, 0) = D2Psi[0][1][0] * Psi[1][0][1];
    d2psids(1, 3, 0) = D2Psi[0][1][1] * Psi[1][0][1];

    // Set up the functions at the corner 2
    d2psids(2, 0, 0) = D2Psi[0][0][0] * Psi[1][1][0];
    d2psids(2, 1, 0) = D2Psi[0][0][1] * Psi[1][1][0];
    d2psids(2, 2, 0) = D2Psi[0][0][0] * Psi[1][1][1];
    d2psids(2, 3, 0) = D2Psi[0][0][1] * Psi[1][1][1];

    // Set up the functions at corner 3
    d2psids(3, 0, 0) = D2Psi[0][1][0] * Psi[1][1][0];
    d2psids(3, 1, 0) = D2Psi[0][1][1] * Psi[1][1][0];
    d2psids(3, 2, 0) = D2Psi[0][1][0] * Psi[1][1][1];
    d2psids(3, 3, 0) = D2Psi[0][1][1] * Psi[1][1][1];

    // D^2/Ds[1]^2

    // Set up the functions at corner 0
    d2psids(0, 0, 1) = Psi[0][0][0] * D2Psi[1][0][0];
    d2psids(0, 1, 1) = Psi[0][0][1] * D2Psi[1][0][0];
    d2psids(0, 2, 1) = Psi[0][0][0] * D2Psi[1][0][1];
    d2psids(0, 3, 1) = Psi[0][0][1] * D2Psi[1][0][1];

    // Set up the functions at corner 1
    d2psids(1, 0, 1) = Psi[0][1][0] * D2Psi[1][0][0];
    d2psids(1, 1, 1) = Psi[0][1][1] * D2Psi[1][0][0];
    d2psids(1, 2, 1) = Psi[0][1][0] * D2Psi[1][0][1];
    d2psids(1, 3, 1) = Psi[0][1][1] * D2Psi[1][0][1];

    // Set up the functions at the corner 2
    d2psids(2, 0, 1) = Psi[0][0][0] * D2Psi[1][1][0];
    d2psids(2, 1, 1) = Psi[0][0][1] * D2Psi[1][1][0];
    d2psids(2, 2, 1) = Psi[0][0][0] * D2Psi[1][1][1];
    d2psids(2, 3, 1) = Psi[0][0][1] * D2Psi[1][1][1];

    // Set up the functions at corner 3
    d2psids(3, 0, 1) = Psi[0][1][0] * D2Psi[1][1][0];
    d2psids(3, 1, 1) = Psi[0][1][1] * D2Psi[1][1][0];
    d2psids(3, 2, 1) = Psi[0][1][0] * D2Psi[1][1][1];
    d2psids(3, 3, 1) = Psi[0][1][1] * D2Psi[1][1][1];

    // D^2/Ds[0]Ds[1]

    // Set up the functions at corner 0
    d2psids(0, 0, 2) = DPsi[0][0][0] * DPsi[1][0][0];
    d2psids(0, 1, 2) = DPsi[0][0][1] * DPsi[1][0][0];
    d2psids(0, 2, 2) = DPsi[0][0][0] * DPsi[1][0][1];
    d2psids(0, 3, 2) = DPsi[0][0][1] * DPsi[1][0][1];

    // Set up the functions at corner 1
    d2psids(1, 0, 2) = DPsi[0][1][0] * DPsi[1][0][0];
    d2psids(1, 1, 2) = DPsi[0][1][1] * DPsi[1][0][0];
    d2psids(1, 2, 2) = DPsi[0][1][0] * DPsi[1][0][1];
    d2psids(1, 3, 2) = DPsi[0][1][1] * DPsi[1][0][1];

    // Set up the functions at the corner 2
    d2psids(2, 0, 2) = DPsi[0][0][0] * DPsi[1][1][0];
    d2psids(2, 1, 2) = DPsi[0][0][1] * DPsi[1][1][0];
    d2psids(2, 2, 2) = DPsi[0][0][0] * DPsi[1][1][1];
    d2psids(2, 3, 2) = DPsi[0][0][1] * DPsi[1][1][1];

    // Set up the functions at corner 3
    d2psids(3, 0, 2) = DPsi[0][1][0] * DPsi[1][1][0];
    d2psids(3, 1, 2) = DPsi[0][1][1] * DPsi[1][1][0];
    d2psids(3, 2, 2) = DPsi[0][1][0] * DPsi[1][1][1];
    d2psids(3, 3, 2) = DPsi[0][1][1] * DPsi[1][1][1];
  }


  //===========================================================
  /// The output function for QHermiteElement<2,ORDER>
  //===========================================================
  template<>
  void QHermiteElement<2>::output(std::ostream& outfile)
  {
    // Tecplot header info
    outfile << "ZONE I=" << 2 << ", J=" << 2 << std::endl;

    // Find the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l2 = 0; l2 < 2; l2++)
    {
      for (unsigned l1 = 0; l1 < 2; l1++)
      {
        unsigned l = l2 * 2 + l1;

        // Loop over the dimensions and output the position
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << node_pt(l)->x(i) << " ";
        }

        // Find out number of types of dof stored at each node
        unsigned n_position_type = node_pt(l)->nposition_type();
        // Loop over the additional positional dofs
        for (unsigned k = 1; k < n_position_type; k++)
        {
          for (unsigned i = 0; i < n_dim; i++)
          {
            outfile << node_pt(l)->x_gen(k, i) << " ";
          }
        }

        // Find out how many data values at the node
        unsigned initial_nvalue = node_pt(l)->nvalue();
        // Loop over the data and output whether pinned or not
        for (unsigned i = 0; i < initial_nvalue; i++)
        {
          outfile << node_pt(l)->is_pinned(i) << " ";
        }
        outfile << std::endl;
      }
    }
    outfile << std::endl;
  }

  //=======================================================================
  /// The output function for n_plot points in each coordinate direction
  //=======================================================================
  template<>
  void QHermiteElement<2>::output(std::ostream& outfile, const unsigned& n_plot)
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
          outfile << interpolated_x(s, i) << " ";
        }
      }
    }
    outfile << std::endl;
  }


  //===========================================================
  /// The C-style output function for QHermiteElement<2,ORDER>
  //===========================================================
  template<>
  void QHermiteElement<2>::output(FILE* file_pt)
  {
    // Tecplot header info
    fprintf(file_pt, "ZONE I=2, J=2");

    // Find the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Loop over element nodes
    for (unsigned l2 = 0; l2 < 2; l2++)
    {
      for (unsigned l1 = 0; l1 < 2; l1++)
      {
        unsigned l = l2 * 2 + l1;

        // Loop over the dimensions and output the position
        for (unsigned i = 0; i < n_dim; i++)
        {
          fprintf(file_pt, "%g ", node_pt(l)->x(i));
        }

        // Find out number of types of dof stored at each node
        unsigned n_position_type = node_pt(l)->nposition_type();
        // Loop over the additional positional dofs
        for (unsigned k = 1; k < n_position_type; k++)
        {
          for (unsigned i = 0; i < n_dim; i++)
          {
            fprintf(file_pt, "%g ", node_pt(l)->x_gen(k, i));
          }
        }

        // Find out how many data values at the node
        unsigned initial_nvalue = node_pt(l)->nvalue();
        // Loop over the data and output whether pinned or not
        for (unsigned i = 0; i < initial_nvalue; i++)
        {
          fprintf(file_pt, "%i ", node_pt(l)->is_pinned(i));
        }
        fprintf(file_pt, "\n");
      }
    }
    fprintf(file_pt, "\n");
  }

  //=======================================================================
  /// The C-style output function for n_plot points in each coordinate direction
  //=======================================================================
  template<>
  void QHermiteElement<2>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "ZONE I=%i, J=%i \n", n_plot, n_plot);

    // Find the dimension of the nodes
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
          fprintf(file_pt, "%g ", interpolated_x(s, i));
        }
      }
    }
    fprintf(file_pt, "\n");
  }


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (of type QHermiteElement<1>).
  //=======================================================================
  template<>
  void QHermiteElement<2>::build_face_element(const int& face_index,
                                              FaceElement* face_element_pt)
  {
    // Set the nodal dimension from the "bulk"
    face_element_pt->set_nodal_dimension(node_pt(0)->ndim());

    // Set the pointer to the "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // Pass on non-halo proc ID
    face_element_pt->set_halo(Non_halo_proc_ID);
#endif

    // Resize the bulk_position_type translation scheme to the number of
    // position types in the 1D element: 2 (position and slope)
    face_element_pt->bulk_position_type_resize(2);

    // Resize the storage for the original number of values at
    // the two nodes of the FaceElement
    face_element_pt->nbulk_value_resize(2);

    // Resize the storage for the bulk node numbers coressponding to the
    // two nodes of the FaceElement
    face_element_pt->bulk_node_number_resize(2);

    // Set the face index in the face element
    // The faces are
    //                       +1    East
    //                       -1    West
    //                       +2    North
    //                       -2    South

    // Set the face index in the face element
    face_element_pt->face_index() = face_index;

    // Now set up the node pointers
    // The convention here is that interior_tangent X tangent X tangent
    // is the OUTWARD normal
    // IMPORTANT NOTE: Need to ensure that numbering is consistent here
    // i.e. node numbers increase in positive x and y directions as they should
    // If not, normals will be inward.
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

        for (unsigned i = 0; i < 2; i++)
        {
          bulk_number = i * 2;
          face_element_pt->node_pt(i) = node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = 1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
          // Set the position type for the slope, which is in the s[1]
          // direction, so is position_type 2 in the "bulk" element
          face_element_pt->bulk_position_type(1) = 2;
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

        for (unsigned i = 0; i < 2; i++)
        {
          bulk_number = i;
          face_element_pt->node_pt(i) = node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = 1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
          // Set the position type for the slope, which is in the s[0]
          // direction, so is position_type 1 in "bulk" element
          face_element_pt->bulk_position_type(1) = 1;
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

        for (unsigned i = 0; i < 2; i++)
        {
          bulk_number = 2 * i + 1;
          face_element_pt->node_pt(i) = node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = -1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
          // Set the position type for the slope, which is in the s[1]
          // direction, so is position_type 2 in the bulk element
          face_element_pt->bulk_position_type(1) = 2;
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

        for (unsigned i = 0; i < 2; i++)
        {
          bulk_number = 2 + i;
          face_element_pt->node_pt(i) = node_pt(bulk_number);
          face_element_pt->bulk_node_number(i) = bulk_number;
          face_element_pt->normal_sign() = -1;
          // Set the number of values originally stored at this node
          face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
          // Set the position type for the slope, which is in the s[0]
          // direction, so is position_type 1 in the "bulk" element
          face_element_pt->bulk_position_type(1) = 1;
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

  /// //////////////////////////////////////////////////////////////////
  //     1D SolidQHermiteElements
  /// ///////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Overload the output function
  //====================================================================
  template<unsigned DIM>
  void SolidQHermiteElement<DIM>::output(std::ostream& outfile)
  {
    QHermiteElement<DIM>::output(outfile);
  }


  //=======================================================================
  /// The output function for n_plot points in each coordinate direction
  /// for the 1D element
  //=======================================================================
  template<>
  void SolidQHermiteElement<1>::output(std::ostream& outfile,
                                       const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Find the Lagrangian dimension of the first node
    unsigned n_lagr = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Output the Eulerian coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }

      // Output the Lagrangian coordinates
      for (unsigned i = 0; i < n_lagr; i++)
      {
        outfile << interpolated_xi(s, i) << " ";
      }
      outfile << std::endl;
    }
  }


  //=====================================================================
  /// Overload the C-style output function
  //====================================================================
  template<unsigned DIM>
  void SolidQHermiteElement<DIM>::output(FILE* file_pt)
  {
    QHermiteElement<DIM>::output(file_pt);
  }


  //=======================================================================
  /// The C-style output function for n_plot points in each coordinate direction
  /// for the 1D element
  //=======================================================================
  template<>
  void SolidQHermiteElement<1>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    fprintf(file_pt, "ZONE I=%i\n", n_plot);

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Find the Lagrangian dimension of the first node
    unsigned n_lagr = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = -1.0 + l * 2.0 / (n_plot - 1);

      // Output the Eulerian coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }

      // Output the Lagrangian coordinates
      for (unsigned i = 0; i < n_lagr; i++)
      {
        fprintf(file_pt, "%g ", interpolated_xi(s, i));
      }
      fprintf(file_pt, "\n");
    }
  }


  /// ////////////////////////////////////////////////////////////////////////
  //   2D SolidQHermiteElements
  /// ///////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// The output function for any number of points per element
  //=====================================================================
  template<>
  void SolidQHermiteElement<2>::output(std::ostream& outfile,
                                       const unsigned& n_p)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    outfile << "ZONE I=" << n_p << ", J=" << n_p << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Find the Lagrangian dimension of the first node
    unsigned n_lagr = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Loop over element nodes
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_p - 1);
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_p - 1);

        // Output the Eulerian coordinates
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_x(s, i) << " ";
        }

        // Output the Lagrangian coordinates
        for (unsigned i = 0; i < n_lagr; i++)
        {
          outfile << interpolated_xi(s, i) << " ";
        }
        outfile << std::endl;
      }
    }
    outfile << std::endl;
  }


  //=====================================================================
  /// The C-style output function for any number of points per element
  //=====================================================================
  template<>
  void SolidQHermiteElement<2>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(2);

    // Tecplot header info
    fprintf(file_pt, "ZONE I=%i, J=%i\n", n_plot, n_plot);

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Find the Lagrangian dimension of the first node
    unsigned n_lagr = static_cast<SolidNode*>(node_pt(0))->nlagrangian();

    // Loop over element nodes
    for (unsigned l2 = 0; l2 < n_plot; l2++)
    {
      s[1] = -1.0 + l2 * 2.0 / (n_plot - 1);
      for (unsigned l1 = 0; l1 < n_plot; l1++)
      {
        s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

        // Output the Eulerian coordinates
        for (unsigned i = 0; i < n_dim; i++)
        {
          fprintf(file_pt, "%g ", interpolated_x(s, i));
        }

        // Output the Lagrangian coordinates
        for (unsigned i = 0; i < n_lagr; i++)
        {
          fprintf(file_pt, "%g ", interpolated_xi(s, i));
        }
        fprintf(file_pt, "\n");
      }
    }
    fprintf(file_pt, "\n");
  }


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements for the solid hermite elements. We need to
  /// construct the basic  element and then sort out the Lagrangian
  /// coordinates
  //=======================================================================
  template<unsigned DIM>
  void SolidQHermiteElement<DIM>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    // Build the standard non-solid FaceElement
    QHermiteElement<DIM>::build_face_element(face_index, face_element_pt);

    // Set the Lagrangian dimension from the first node of the present element
    dynamic_cast<SolidFiniteElement*>(face_element_pt)
      ->set_lagrangian_dimension(
        static_cast<SolidNode*>(node_pt(0))->nlagrangian());
  }

  //==================================================================
  /// Force build of templates
  //==================================================================
  template class QHermiteElement<1>;
  template class QHermiteElement<2>;
  template class DiagQHermiteElement<1>;
  template class DiagQHermiteElement<2>;


  template class SolidQHermiteElement<1>;
  template class SolidQHermiteElement<2>;
  template class SolidDiagQHermiteElement<1>;
  template class SolidDiagQHermiteElement<2>;

} // namespace oomph

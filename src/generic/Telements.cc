// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Non-inline member functions for Telements


// oomph-lib headers
#include "Telements.h"


namespace oomph
{
  //=======================================================================
  /// Assign the static integral
  //=======================================================================
  template<unsigned NNODE_1D>
  TGauss<1, NNODE_1D> TElement<1, NNODE_1D>::Default_integration_scheme;
  template<unsigned NNODE_1D>
  TGauss<2, NNODE_1D> TElement<2, NNODE_1D>::Default_integration_scheme;
  template<unsigned NNODE_1D>
  TGauss<3, NNODE_1D> TElement<3, NNODE_1D>::Default_integration_scheme;

  /// ///////////////////////////////////////////////////////////////
  // 1D Telements
  /// ///////////////////////////////////////////////////////////////

  //=======================================================================
  /// The output function for general 1D TElements
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<1, NNODE_1D>::output(std::ostream& outfile)
  {
    output(outfile, NNODE_1D);
  }

  //=======================================================================
  /// The output function for n_plot points in each coordinate direction
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<1, NNODE_1D>::output(std::ostream& outfile,
                                     const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = l * 1.0 / (n_plot - 1);
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
  /// C style output function for general 1D TElements
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<1, NNODE_1D>::output(FILE* file_pt)
  {
    output(file_pt, NNODE_1D);
  }

  //=======================================================================
  /// C style output function for n_plot points in each coordinate direction
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<1, NNODE_1D>::output(FILE* file_pt, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    // outfile << "ZONE I=" << n_plot << std::endl;
    fprintf(file_pt, "ZONE I=%i\n", n_plot);

    // Find the dimension of the first node
    unsigned n_dim = this->nodal_dimension();

    // Loop over plot points
    for (unsigned l = 0; l < n_plot; l++)
    {
      s[0] = l * 1.0 / (n_plot - 1);

      // Output the coordinates
      for (unsigned i = 0; i < n_dim; i++)
      {
        // outfile << interpolated_x(s,i) << " " ;
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      // outfile <<  std::endl;
      fprintf(file_pt, "\n");
    }
    // outfile << std::endl;
    fprintf(file_pt, "\n");
  }

  //=============================================================
  /// Namespace for helper functions that return the local
  /// coordinates in the bulk elements
  //==============================================================
  namespace TElement1FaceToBulkCoordinates
  {
    /// The translation scheme for the face s0 = 0.0
    void face0(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 0.0;
    }

    /// The translation scheme for the face s0 = 1.0
    void face1(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 1.0;
    }
  } // namespace TElement1FaceToBulkCoordinates


  //=============================================================
  /// Namespace for helper functions that calculate derivatives
  /// of the local coordinates in the bulk elements wrt the
  /// local coordinates in the face element.
  //=============================================================
  namespace TElement1BulkCoordinateDerivatives
  {
    /// Function for both faces -- the bulk coordinate is fixed on both
    void faces0(const Vector<double>& s,
                DenseMatrix<double>& dsbulk_dsface,
                unsigned& interior_direction)
    {
      // Bulk coordinate s[0] does not vary along the face
      dsbulk_dsface(0, 0) = 0.0;

      // The interior direction is given by s[0]
      interior_direction = 0;
    }
  } // namespace TElement1BulkCoordinateDerivatives

  //===========================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are Point Elements).
  //===========================================================
  template<unsigned NNODE_1D>
  void TElement<1, NNODE_1D>::build_face_element(const int& face_index,
                                                 FaceElement* face_element_pt)
  {
    // Overload the nodal dimension
    face_element_pt->set_nodal_dimension(nodal_dimension());

    // Set the pointer to the "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // Pass on non-halo proc ID
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
    // The convention is that the "normal", should always point
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
          &TElement1FaceToBulkCoordinates::face0;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &TElement1BulkCoordinateDerivatives::faces0;

        // Set the number of values stored when the node is part of the "bulk"
        // element. The required_nvalue() must be used, rather than nvalue(),
        // because otherwise nodes on boundaries will be resized multiple
        // times. If you want any other behaviour, you MUST set nbulk_value()
        // manually after construction of your specific face element.
        face_element_pt->nbulk_value(0) = required_nvalue(0);
        break;

        // Top, normal sign is positive (coordinate points out of element)
      case (1):
        face_element_pt->node_pt(0) = node_pt(NNODE_1D - 1);
        face_element_pt->bulk_node_number(0) = NNODE_1D - 1;
        face_element_pt->normal_sign() = +1;

        // Set the pointer to the function that determines the bulk coordinates
        // in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement1FaceToBulkCoordinates::face1;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &TElement1BulkCoordinateDerivatives::faces0;


        // Set the number of values stored when the node is part of the "bulk"
        // element.
        face_element_pt->nbulk_value(0) = required_nvalue(NNODE_1D - 1);
        break;

        // All other cases throw an error
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
  //       2D Telements
  /// /////////////////////////////////////////////////////////////

  /// Assign the nodal translation schemes
  template<>
  const unsigned TElement<2, 2>::Node_on_face[3][2] = {{2, 1}, {2, 0}, {0, 1}};

  template<>
  const unsigned TElement<2, 3>::Node_on_face[3][3] = {
    {2, 4, 1}, {2, 5, 0}, {0, 3, 1}};

  template<>
  const unsigned TElement<2, 4>::Node_on_face[3][4] = {
    {2, 6, 5, 1}, {2, 7, 8, 0}, {0, 3, 4, 1}};


  //===================================================================
  /// Namespace for the functions that translate local face coordinates
  /// to the coordinates in the bulk element
  //==================================================================
  namespace TElement2FaceToBulkCoordinates
  {
    /// The translation scheme for the face s0 = 0
    void face0(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 0.0;
      s_bulk[1] = s[0];
    }

    /// The translation scheme for the face s1 = 0
    void face1(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = s[0];
      s_bulk[1] = 0.0;
    }

    /// The translation scheme for the face s2 = 0
    void face2(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 1.0 - s[0];
      s_bulk[1] = s[0];
    }
  } // namespace TElement2FaceToBulkCoordinates


  //=============================================================
  /// Namespace for helper functions that calculate derivatives
  /// of the local coordinates in the bulk elements wrt the
  /// local coordinates in the face element.
  //=============================================================
  namespace TElement2BulkCoordinateDerivatives
  {
    /// Function for the "left" face along which s0 is fixed
    void face0(const Vector<double>& s,
               DenseMatrix<double>& dsbulk_dsface,
               unsigned& interior_direction)
    {
      // Bulk coordinate s[0] does not vary along the face
      dsbulk_dsface(0, 0) = 0.0;
      // Bulk coordinate s[1] is the face coordinate
      dsbulk_dsface(1, 0) = 1.0;

      // The interior direction is given by s[0]
      interior_direction = 0;
    }


    /// Function for the "bottom" face along which s1 is fixed
    void face1(const Vector<double>& s,
               DenseMatrix<double>& dsbulk_dsface,
               unsigned& interior_direction)
    {
      // Bulk coordinate s[0] is the face coordinate
      dsbulk_dsface(0, 0) = 1.0;
      // Bulk coordinate s[1] does not vary along the face
      dsbulk_dsface(1, 0) = 0.0;

      // The interior direction is given by s[1]
      interior_direction = 1;
    }

    /// Function for the sloping face
    void face2(const Vector<double>& s,
               DenseMatrix<double>& dsbulk_dsface,
               unsigned& interior_direction)
    {
      // Bulk coordinate s[0] decreases along the face
      dsbulk_dsface(0, 0) = -1.0;
      // Bulk coordinate s[1] increases along the face
      dsbulk_dsface(1, 0) = 1.0;

      // The interior direction is given by s[0] (or s[1])
      interior_direction = 0;
    }

  } // namespace TElement2BulkCoordinateDerivatives


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type TElement<2,NNODE_1D>).
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<2, NNODE_1D>::build_face_element(const int& face_index,
                                                 FaceElement* face_element_pt)
  {
    // Set the nodal dimension from the first node
    face_element_pt->set_nodal_dimension(nodal_dimension());

    // Set the pointer to the orginal "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // Pass on non-halo proc ID
    face_element_pt->set_halo(Non_halo_proc_ID);
#endif

    // Calculate the number of nodes in the face element
    const unsigned n_face_nodes = NNODE_1D;

    // Resize storage for the number of values originally stored
    // at the face element's nodes.
    face_element_pt->nbulk_value_resize(n_face_nodes);

    // Resize storage for the bulk node numbers corresponding to
    // the nodes of the face
    face_element_pt->bulk_node_number_resize(n_face_nodes);

    // Set the face index in the face element
    face_element_pt->face_index() = face_index;

    // So the faces are
    // 0 : s_0 fixed
    // 1 : s_1 fixed
    // 2 : sloping face

    // Copy nodes across to the face
    for (unsigned i = 0; i < n_face_nodes; i++)
    {
      // The number is just offset by one
      unsigned bulk_number = Node_on_face[face_index][i];
      face_element_pt->node_pt(i) = node_pt(bulk_number);
      face_element_pt->bulk_node_number(i) = bulk_number;
      // set the number of values originally stored at this node
      face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
    }

    // Now set up the node pointers and the normal vectors
    switch (face_index)
    {
        //
        //-----The face s0 is constant
      case 0:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement2FaceToBulkCoordinates::face0;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &TElement2BulkCoordinateDerivatives::face0;

        // Outer unit normal is the positive cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = +1;

        break;

        //-----The face s1 is constant
      case 1:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement2FaceToBulkCoordinates::face1;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &TElement2BulkCoordinateDerivatives::face1;

        // Outer unit normal is the positive cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = +1;

        break;

        //
        //-----The face s2 is constant
      case 2:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement2FaceToBulkCoordinates::face2;

        // Set the pointer to the function that determines the mapping of
        // derivatives
        face_element_pt->bulk_coordinate_derivatives_fct_pt() =
          &TElement2BulkCoordinateDerivatives::face2;

        // Outer unit normal is the negative cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;

        break;

      default:

        std::ostringstream error_message;
        error_message << "Face_index should only take "
                      << "the values 0, 1 or 2 not " << face_index << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    } // end switch
  }


  //=======================================================================
  /// The output function for TElement<2,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<2, NNODE_1D>::output(std::ostream& outfile)
  {
    // QUEHACERES want to perform same output, but at each node
    output(outfile, NNODE_1D);
  }


  //=======================================================================
  /// The output function for TElement<2,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<2, NNODE_1D>::output(std::ostream& outfile,
                                     const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Get the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //=======================================================================
  /// The C-style output function for TElement<2,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<2, NNODE_1D>::output(FILE* file_pt)
  {
    output(file_pt, NNODE_1D);
  }


  //=======================================================================
  /// The C-style output function for TElement<2,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<2, NNODE_1D>::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(2);

    // Find the dimensions of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Tecplot header info
    fprintf(file_pt, "%s \n", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < n_dim; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
        // outfile << interpolated_x(s,i) << " ";
      }
      fprintf(file_pt, "\n");
      // outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// The output function for TElement<3,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<3, NNODE_1D>::output(std::ostream& outfile)
  {
    output(outfile, NNODE_1D);
  }


  //=======================================================================
  /// The output function for TElement<3,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<3, NNODE_1D>::output(std::ostream& outfile,
                                     const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(3);

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << "\n";
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //=======================================================================
  /// The C-style output function for TElement<3,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<3, NNODE_1D>::output(FILE* file_pt)
  {
    output(file_pt, NNODE_1D);
  }


  //=======================================================================
  /// The C-style output function for TElement<3,NNODE_1D>
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<3, NNODE_1D>::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(3);

    // Find the dimension of the nodes
    unsigned n_dim = this->nodal_dimension();

    // Tecplot header info
    fprintf(file_pt, "%s \n", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < n_dim; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
        // outfile << interpolated_x(s,i) << " ";
      }
      fprintf(file_pt, "\n");
      // outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }

  //===================================================================
  /// Namespace for the functions that translate local face coordinates
  /// to the coordinates in the bulk element
  //==================================================================
  namespace TElement3FaceToBulkCoordinates
  {
    /// The translation scheme for the face s0 = 0
    void face0(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 0.0;
      s_bulk[1] = s[0];
      s_bulk[2] = s[1];
    }

    /// The translation scheme for the face s1 = 0
    void face1(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = s[0];
      s_bulk[1] = 0.0;
      s_bulk[2] = s[1];
    }

    /// The translation scheme for the face s2 = 0
    void face2(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = s[0];
      s_bulk[1] = s[1];
      s_bulk[2] = 0.0;
    }

    /// The translation scheme for the sloping face
    void face3(const Vector<double>& s, Vector<double>& s_bulk)
    {
      s_bulk[0] = 1 - s[0] - s[1];
      s_bulk[1] = s[0];
      s_bulk[2] = s[1];
    }

  } // namespace TElement3FaceToBulkCoordinates

  /// Assign the nodal translation scheme for linear elements
  template<>
  const unsigned TElement<3, 2>::Node_on_face[4][3] = {
    {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {1, 2, 0}};

  /// Assign the nodal translation scheme for quadratic elements
  template<>
  const unsigned TElement<3, 3>::Node_on_face[4][6] = {{1, 2, 3, 7, 8, 9},
                                                       {0, 2, 3, 5, 8, 6},
                                                       {0, 1, 3, 4, 9, 6},
                                                       {1, 2, 0, 7, 5, 4}};


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type TElement<2,NNODE_1D>).
  //=======================================================================
  template<unsigned NNODE_1D>
  void TElement<3, NNODE_1D>::build_face_element(const int& face_index,
                                                 FaceElement* face_element_pt)
  {
    // Set the nodal dimension
    face_element_pt->set_nodal_dimension(nodal_dimension());

    // Set the pointer to the orginal "bulk" element
    face_element_pt->bulk_element_pt() = this;

#ifdef OOMPH_HAS_MPI
    // Pass on non-halo proc ID
    face_element_pt->set_halo(Non_halo_proc_ID);
#endif

    // Calculate the number of nodes in the face element
    const unsigned n_face_nodes = (NNODE_1D * (NNODE_1D + 1)) / 2;

    // Resize storage for the number of values originally stored
    // at the face element's nodes.
    face_element_pt->nbulk_value_resize(n_face_nodes);

    // Resize storage for the bulk node numbers corresponding to
    // the nodes of the face
    face_element_pt->bulk_node_number_resize(n_face_nodes);

    // Set the face index in the face element
    face_element_pt->face_index() = face_index;

    // So the faces are
    // 0 : s_0 fixed
    // 1 : s_1 fixed
    // 2 : s_2 fixed
    // 3 : sloping face

    // Copy nodes across to the face
    for (unsigned i = 0; i < n_face_nodes; i++)
    {
      // The number is just offset by one
      unsigned bulk_number = Node_on_face[face_index][i];
      face_element_pt->node_pt(i) = node_pt(bulk_number);
      face_element_pt->bulk_node_number(i) = bulk_number;
      // set the number of values originally stored at this node
      face_element_pt->nbulk_value(i) = required_nvalue(bulk_number);
    }

    // Now set up the node pointers and the normal vectors
    switch (face_index)
    {
        //
        //-----The face s0 is constant
      case 0:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement3FaceToBulkCoordinates::face0;

        // Outer unit normal is the negative cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;

        break;

        //-----The face s1 is constant
      case 1:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement3FaceToBulkCoordinates::face1;

        // Outer unit normal is the positive cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = +1;

        break;

        //
        //-----The face s2 is constant
      case 2:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement3FaceToBulkCoordinates::face2;

        // Outer unit normal is the negative cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = -1;

        break;

        //-----The sloping face of the tetrahedron
      case 3:

        // Set the pointer to the function that determines the bulk
        // coordinates in the face element
        face_element_pt->face_to_bulk_coordinate_fct_pt() =
          &TElement3FaceToBulkCoordinates::face3;

        // Outer unit normal is the positive cross product of two in plane
        // tangent vectors
        face_element_pt->normal_sign() = +1;

        break;


      default:

        std::ostringstream error_message;
        error_message << "Face_index should only take "
                      << "the values 0, 1, 2 or 3, not " << face_index
                      << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    } // end switch
  }


  //==================================================================
  // Default integration scheme for the TBubbleEnrichedElement<2,3>
  //===================================================================
  template<unsigned DIM>
  TBubbleEnrichedGauss<DIM, 3>
    TBubbleEnrichedElement<DIM, 3>::Default_enriched_integration_scheme;

  //===================================================================
  // Central node on the face of the TBubbleEnrichedelement<2,3>
  //===================================================================
  template<>
  const unsigned TBubbleEnrichedElement<2, 3>::Central_node_on_face[3] = {
    4, 5, 3};


  //================================================================
  /// The face element for is the same in the two-dimesional case
  //================================================================
  template<>
  void TBubbleEnrichedElement<2, 3>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    TElement<2, 3>::build_face_element(face_index, face_element_pt);
  }


  //===================================================================
  // Central node on the face of the TBubbleEnrichedelement<3,3>
  //===================================================================
  template<>
  const unsigned TBubbleEnrichedElement<3, 3>::Central_node_on_face[4] = {
    13, 12, 10, 11};


  //=======================================================================
  /// Function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type TBubbleEnrichedElement<2,3>).
  //=======================================================================
  template<>
  void TBubbleEnrichedElement<3, 3>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    // Call the standard unenriched build function
    TElement<3, 3>::build_face_element(face_index, face_element_pt);

    // Set the enriched number of total face nodes
    const unsigned n_face_nodes = 7;

    // Resize storage for the number of values originally stored
    // at the face element's nodes.
    face_element_pt->nbulk_value_resize(n_face_nodes);

    // Resize storage for the bulk node numbers corresponding to
    // the nodes of the face
    face_element_pt->bulk_node_number_resize(n_face_nodes);

    // So the faces are
    // 0 : s_0 fixed
    // 1 : s_1 fixed
    // 2 : s_2 fixed
    // 3 : sloping face

    // Copy central node across
    unsigned bulk_number = Central_node_on_face[face_index];
    face_element_pt->node_pt(n_face_nodes - 1) = node_pt(bulk_number);
    face_element_pt->bulk_node_number(n_face_nodes - 1) = bulk_number;
    // set the number of values originally stored at this node
    face_element_pt->nbulk_value(n_face_nodes - 1) =
      required_nvalue(bulk_number);
  }

  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////

  //===========================================================
  /// Final override for
  /// function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type SolidTBubbleEnrichedElement<1,3>).
  //===========================================================
  template<>
  void SolidTBubbleEnrichedElement<2, 3>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    // Build the standard non-solid FaceElement
    TBubbleEnrichedElement<2, 3>::build_face_element(face_index,
                                                     face_element_pt);

    // Set the Lagrangian dimension from the first node of the present element
    dynamic_cast<SolidFiniteElement*>(face_element_pt)
      ->set_lagrangian_dimension(
        static_cast<SolidNode*>(node_pt(0))->nlagrangian());
  }


  //===========================================================
  /// Final override for
  /// function to setup geometrical information for lower-dimensional
  /// FaceElements (which are of type SolidTBubbleEnrichedElement<2,3>).
  //===========================================================
  template<>
  void SolidTBubbleEnrichedElement<3, 3>::build_face_element(
    const int& face_index, FaceElement* face_element_pt)
  {
    // Build the standard non-solid FaceElement
    TBubbleEnrichedElement<3, 3>::build_face_element(face_index,
                                                     face_element_pt);

    // Set the Lagrangian dimension from the first node of the present element
    dynamic_cast<SolidFiniteElement*>(face_element_pt)
      ->set_lagrangian_dimension(
        static_cast<SolidNode*>(node_pt(0))->nlagrangian());
  }


  //===================================================================
  // Build required templates
  //===================================================================
  template class TElement<1, 2>;
  template class TElement<1, 3>;
  template class TElement<1, 4>;
  template class TElement<2, 2>;
  template class TElement<2, 3>;
  template class TElement<2, 4>;
  template class TElement<3, 2>;
  template class TElement<3, 3>;
  template class TBubbleEnrichedElement<2, 3>;
  template class TBubbleEnrichedElement<3, 3>;
  template class SolidTBubbleEnrichedElement<2, 3>;
  template class SolidTBubbleEnrichedElement<3, 3>;
  template class TBubbleEnrichedGauss<2, 3>;
  template class TBubbleEnrichedGauss<3, 3>;
} // namespace oomph

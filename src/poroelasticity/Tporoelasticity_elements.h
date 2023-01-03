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
#ifndef OOMPH_TPOROELASTICITY_ELEMENTS_HEADER
#define OOMPH_TPOROELASTICITY_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "poroelasticity_elements.h"
#include "generic/Telements.h"

namespace oomph
{
  /// Element which solves the Darcy equations using TElements
  template<unsigned ORDER>
  class TPoroelasticityElement : public TElement<2, 3>,
                                 public PoroelasticityEquations<2>
  {
  private:
    /// The number of values stored at each node
    static const unsigned Initial_Nvalue[];

    /// Conversion scheme from an edge degree of freedom to the node it's stored
    /// at
    static const unsigned Q_edge_conv[];

    /// The points along each edge where the fluxes are taken to be
    static const double Gauss_point[];

    /// The internal data index where the internal q degrees of freedom are
    /// stored
    unsigned Q_internal_data_index;

    /// The internal data index where the p degrees of freedom are stored
    unsigned P_internal_data_index;

    /// Unit normal signs associated with each edge to ensure
    /// inter-element continuity of the flux
    std::vector<short> Sign_edge;

  public:
    /// Constructor
    TPoroelasticityElement();

    /// Destructor
    ~TPoroelasticityElement();

    /// Number of values required at node n
    unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue[n];
    }

    unsigned u_index(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= 2)
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,1)";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return n;
    }

    /// Return the equation number of the n-th edge (flux) degree of freedom
    int q_edge_local_eqn(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= nq_basis_edge())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nq_basis_edge() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->nodal_local_eqn(q_edge_node_number(n), q_edge_index(n));
    }

    /// Return the equation number of the n-th internal (moment) degree of
    /// freedom
    int q_internal_local_eqn(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis() - nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << (nq_basis() - nq_basis_edge()) - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return internal_local_eqn(q_internal_index(), n);
    }

    /// Return the index of the internal data where the q_internal
    /// degrees of freedom are stored
    unsigned q_internal_index() const
    {
      return Q_internal_data_index;
    }

    /// Return the nodal index at which the nth edge unknown is stored
    unsigned q_edge_index(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nq_basis_edge() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return n % (ORDER + 1) + 2;
    }

    /// Return the number of the node where the nth edge unknown is stored
    unsigned q_edge_node_number(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nq_basis_edge() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Q_edge_conv[n / (ORDER + 1)];
    }

    /// Return the values of the edge (flux) degrees of freedom
    double q_edge(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nq_basis_edge() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return nodal_value(q_edge_node_number(n), q_edge_index(n));
    }

    /// Return the values of the edge (flux) degrees of freedom at time
    /// history level t
    double q_edge(const unsigned& t, const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nq_basis_edge() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return nodal_value(t, q_edge_node_number(n), q_edge_index(n));
    }

    /// Return the values of the internal (moment) degrees of freedom
    double q_internal(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis() - nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << (nq_basis() - nq_basis_edge()) - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->internal_data_pt(q_internal_index())->value(n);
    }

    /// Return the values of the internal (moment) degrees of freedom at
    /// time history level t
    double q_internal(const unsigned& t, const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      // mjr TODO add time history level range checking
      if (n >= (nq_basis() - nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << (nq_basis() - nq_basis_edge()) - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->internal_data_pt(q_internal_index())->value(t, n);
    }

    /// Return the total number of computational basis functions for u
    unsigned nq_basis() const;

    /// Return the number of edge basis functions for u
    unsigned nq_basis_edge() const;

    /// Returns the local form of the q basis at local coordinate s
    void get_q_basis_local(const Vector<double>& s, Shape& q_basis) const;

    /// Returns the local form of the q basis and dbasis/ds at local coordinate
    /// s
    void get_div_q_basis_local(const Vector<double>& s,
                               Shape& div_q_basis_ds) const;

    /// Returns the number of gauss points along each edge of the element
    unsigned nedge_gauss_point() const;

    /// Returns the nth gauss point along an edge: if sign_edge(edge)==1,
    /// returns regular gauss point; if sign_edge(edge)==-1, returns 1-(gauss
    /// point)
    double edge_gauss_point(const unsigned& edge, const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (edge >= 3)
      {
        std::ostringstream error_message;
        error_message << "Range Error: edge " << edge
                      << " is not in the range (0,2)";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (n >= nedge_gauss_point())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nedge_gauss_point() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return (1 - sign_edge(edge)) / 2.0 + sign_edge(edge) * Gauss_point[n];
    }

    /// Returns the global coordinates of the nth gauss point along an edge
    void edge_gauss_point_global(const unsigned& edge,
                                 const unsigned& n,
                                 Vector<double>& x) const
    {
#ifdef RANGE_CHECKING
      if (edge >= 3)
      {
        std::ostringstream error_message;
        error_message << "Range Error: edge " << edge
                      << " is not in the range (0,2)";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (n >= nedge_gauss_point())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << nedge_gauss_point() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get the location of the n-th gauss point along the edge in terms of the
      // distance along the edge itself
      double gauss_point = edge_gauss_point(edge, n);

      // Convert the edge number to the number of the mid-edge node along that
      // edge
      unsigned node_number = Q_edge_conv[edge];

      // Storage for the local coords of the gauss point
      Vector<double> s_gauss(2, 0);

      // The edge basis functions are defined in a clockwise manner, so we have
      // to effectively "flip" the coordinates along edges 0 and 1 to match this
      switch (node_number)
      {
        case 3:
          s_gauss[0] = 1 - gauss_point;
          s_gauss[1] = gauss_point;
          break;
        case 4:
          s_gauss[0] = 0;
          s_gauss[1] = 1 - gauss_point;
          break;
        case 5:
          s_gauss[0] = gauss_point;
          s_gauss[1] = 0;
          break;
      }

      // Calculate the global coordinates from the local ones
      interpolated_x(s_gauss, x);
    }

    /// Pin the nth internal q value
    void pin_q_internal_value(const unsigned& n)
    {
#ifdef RANGE_CHECKING
      if (n >= (nq_basis() - nq_basis_edge()))
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << (nq_basis() - nq_basis_edge()) - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      this->internal_data_pt(q_internal_index())->pin(n);
    }

    /// Return the equation number of the n-th pressure degree of freedom
    int p_local_eqn(const unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= np_basis())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << np_basis() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->internal_local_eqn(P_internal_data_index, n);
    }

    /// Return the nth pressure value
    double p_value(unsigned& n) const
    {
#ifdef RANGE_CHECKING
      if (n >= np_basis())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << np_basis() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->internal_data_pt(P_internal_data_index)->value(n);
    }

    /// Return the total number of pressure basis functions
    unsigned np_basis() const;

    /// Return the pressure basis
    void get_p_basis(const Vector<double>& s, Shape& p_basis) const;

    /// Pin the nth pressure value
    void pin_p_value(const unsigned& n, const double& p)
    {
#ifdef RANGE_CHECKING
      if (n >= np_basis())
      {
        std::ostringstream error_message;
        error_message << "Range Error: n " << n << " is not in the range (0,"
                      << np_basis() - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      this->internal_data_pt(P_internal_data_index)->pin(n);
      this->internal_data_pt(P_internal_data_index)->set_value(n, p);
    }

    /// Scale the edge basis to allow arbitrary edge mappings
    void scale_basis(Shape& basis) const
    {
      // Storage for the lengths of the edges of the element
      Vector<double> length(3, 0.0);

      // Temporary storage for the vertex positions
      double x0, y0, x1, y1;

      // loop over the edges of the element and calculate their lengths (in x-y
      // space)
      for (unsigned i = 0; i < 3; i++)
      {
        x0 = this->node_pt(i)->x(0);
        y0 = this->node_pt(i)->x(1);
        x1 = this->node_pt((i + 1) % 3)->x(0);
        y1 = this->node_pt((i + 1) % 3)->x(1);

        length[i] = std::sqrt(std::pow(y1 - y0, 2) + std::pow(x1 - x0, 2));
      }

      // lengths of the sides of the reference element (in the same order as the
      // basis functions)
      const double ref_length[3] = {std::sqrt(2.0), 1, 1};

      // get the number of basis functions associated with the edges
      unsigned n_q_basis_edge = nq_basis_edge();

      // rescale the edge basis functions to allow arbitrary edge mappings from
      // element to ref. element
      const unsigned n_index2 = basis.nindex2();
      for (unsigned i = 0; i < n_index2; i++)
      {
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          basis(l, i) *=
            (length[l / (ORDER + 1)] / ref_length[l / (ORDER + 1)]);
        }
      }
    }

    /// Accessor for the unit normal sign of edge n (const version)
    const short& sign_edge(const unsigned& n) const
    {
      return Sign_edge[n];
    }

    /// Accessor for the unit normal sign of edge n
    short& sign_edge(const unsigned& n)
    {
      return Sign_edge[n];
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      PoroelasticityEquations<2>::output(outfile);
    }

    /// Output FE representation of soln: x,y,u1,u2,div_q,p at
    /// Nplot^DIM plot points
    void output(std::ostream& outfile, const unsigned& Nplot)
    {
      PoroelasticityEquations<2>::output(outfile, Nplot);
    }

  protected:
    /// Returns the geometric basis, and the u, p and divergence basis
    /// functions and test functions at local coordinate s
    double shape_basis_test_local(const Vector<double>& s,
                                  Shape& psi,
                                  DShape& dpsi,
                                  Shape& u_basis,
                                  Shape& u_test,
                                  DShape& du_basis_dx,
                                  DShape& du_test_dx,
                                  Shape& q_basis,
                                  Shape& q_test,
                                  Shape& p_basis,
                                  Shape& p_test,
                                  Shape& div_q_basis_ds,
                                  Shape& div_q_test_ds) const
    {
      const unsigned n_q_basis = this->nq_basis();

      Shape q_basis_local(n_q_basis, 2);
      this->get_q_basis_local(s, q_basis_local);
      this->get_p_basis(s, p_basis);
      this->get_div_q_basis_local(s, div_q_basis_ds);

      double J = this->transform_basis(s, q_basis_local, psi, dpsi, q_basis);

      // u_basis consists of the normal Lagrangian shape functions
      u_basis = psi;
      du_basis_dx = dpsi;

      u_test = psi;
      du_test_dx = dpsi;

      q_test = q_basis;
      p_test = p_basis;
      div_q_test_ds = div_q_basis_ds;

      return J;
    }

    /// Returns the geometric basis, and the u, p and divergence basis
    /// functions and test functions at integration point ipt
    double shape_basis_test_local_at_knot(const unsigned& ipt,
                                          Shape& psi,
                                          DShape& dpsi,
                                          Shape& u_basis,
                                          Shape& u_test,
                                          DShape& du_basis_dx,
                                          DShape& du_test_dx,
                                          Shape& q_basis,
                                          Shape& q_test,
                                          Shape& p_basis,
                                          Shape& p_test,
                                          Shape& div_q_basis_ds,
                                          Shape& div_q_test_ds) const
    {
      Vector<double> s(2);
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      return shape_basis_test_local(s,
                                    psi,
                                    dpsi,
                                    u_basis,
                                    u_test,
                                    du_basis_dx,
                                    du_test_dx,
                                    q_basis,
                                    q_test,
                                    p_basis,
                                    p_test,
                                    div_q_basis_ds,
                                    div_q_test_ds);
    }
  };

  /// Face geometry for TPoroelasticityElement<0>
  template<>
  class FaceGeometry<TPoroelasticityElement<0>> : public virtual TElement<1, 3>
  {
  public:
    /// Constructor: Call constructor of base
    FaceGeometry() : TElement<1, 3>() {}
  };

  /// Face geometry for TPoroelasticityElement<1>
  template<>
  class FaceGeometry<TPoroelasticityElement<1>> : public virtual TElement<1, 3>
  {
  public:
    /// Constructor: Call constructor of base class
    FaceGeometry() : TElement<1, 3>() {}
  };

} // namespace oomph

#endif

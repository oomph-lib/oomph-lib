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
#ifndef OOMPH_SUPG_ADV_DIFF_ELEMENTS_HEADER
#define OOMPH_SUPG_ADV_DIFF_ELEMENTS_HEADER

#include "../advection_diffusion/refineable_advection_diffusion_elements.h"

namespace oomph
{
  //======================================================================
  /// \short QSUPGAdvectionDiffusionElement<DIM,NNODE_1D> elements are
  /// SUPG-stabilised Advection Diffusion elements with
  /// NNODE_1D nodal points in each coordinate direction. Inherits
  /// from QAdvectionDiffusionElement and overwrites their
  /// test functions
  ///
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QSUPGAdvectionDiffusionElement
    : public virtual QAdvectionDiffusionElement<DIM, NNODE_1D>
  {
  public:
    ///\short  Constructor: Call constructors for underlying
    /// QAdvectionDiffusion element. Initialise stabilisation parameter
    /// to zero
    QSUPGAdvectionDiffusionElement()
      : QAdvectionDiffusionElement<DIM, NNODE_1D>()
    {
      Tau_SUPG = 0.0;
    }

    /// Get stabilisation parameter for the element
    double get_Tau_SUPG()
    {
      return Tau_SUPG;
    }


    /// Set stabilisation parameter for the element to zero
    void switch_off_stabilisation()
    {
      Tau_SUPG = 0.0;
    }


    /// Compute stabilisation parameter for the element
    void compute_stabilisation_parameter()
    {
      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions and their derivatives
      Shape psi(n_node), test(n_node);
      DShape dpsidx(n_node, DIM);

      // Evaluate everything at the element centroid
      Vector<double> s(DIM, 0.0);

      // Call the geometrical shape functions and derivatives
      (void)QElement<DIM, NNODE_1D>::dshape_eulerian(s, psi, dpsidx);

      // Calculate Eulerian coordinates
      Vector<double> interpolated_x(DIM, 0.0);

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions (we're in 2D)
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += this->nodal_position(l, j) * psi[l];
        }
      }

      // Element size: Choose the max. diagonal
      double h = 0;
      if (DIM == 1)
      {
        h = std::fabs(this->vertex_node_pt(1)->x(0) -
                      this->vertex_node_pt(0)->x(0));
      }
      else if (DIM == 2)
      {
        h =
          pow(this->vertex_node_pt(3)->x(0) - this->vertex_node_pt(0)->x(0),
              2) +
          pow(this->vertex_node_pt(3)->x(1) - this->vertex_node_pt(0)->x(1), 2);
        double h1 =
          pow(this->vertex_node_pt(2)->x(0) - this->vertex_node_pt(1)->x(0),
              2) +
          pow(this->vertex_node_pt(2)->x(1) - this->vertex_node_pt(1)->x(1), 2);
        if (h1 > h) h = h1;
        h = sqrt(h);
      }
      else if (DIM == 3)
      {
        // diagonals are from nodes 0-7, 1-6, 2-5, 3-4
        unsigned n1 = 0;
        unsigned n2 = 7;
        for (unsigned i = 0; i < 4; i++)
        {
          double h1 =
            pow(this->vertex_node_pt(n1)->x(0) - this->vertex_node_pt(n2)->x(0),
                2) +
            pow(this->vertex_node_pt(n1)->x(1) - this->vertex_node_pt(n2)->x(1),
                2) +
            pow(this->vertex_node_pt(n1)->x(2) - this->vertex_node_pt(n2)->x(2),
                2);
          if (h1 > h) h = h1;
          n1++;
          n2--;
        }
        h = sqrt(h);
      }

      // Get wind
      Vector<double> wind(DIM);
      // Dummy ipt argument?
      unsigned ipt = 0;
      this->get_wind_adv_diff(ipt, s, interpolated_x, wind);
      double abs_wind = 0;
      for (unsigned j = 0; j < DIM; j++)
      {
        abs_wind += wind[j] * wind[j];
      }
      abs_wind = sqrt(abs_wind);

      // Mesh Peclet number
      double Pe_mesh = 0.5 * this->pe() * h * abs_wind;

      // Stabilisation parameter
      if (Pe_mesh > 1.0)
      {
        Tau_SUPG = h / (2.0 * abs_wind) * (1.0 - 1.0 / Pe_mesh);
      }
      else
      {
        Tau_SUPG = 0.0;
      }
    }


    /// \short Output function:
    /// x,y,u,w_x,w_y,tau_supg  or    x,y,z,u,w_x,w_y,w_z,tau_supg
    /// nplot points in each coordinate direction
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        // Get Eulerian coordinate of plot point
        Vector<double> x(DIM);
        this->interpolated_x(s, x);

        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << x[i] << " ";
        }
        outfile << this->interpolated_u_adv_diff(s) << " ";

        // Get the wind
        Vector<double> wind(DIM);
        // Dummy ipt argument
        unsigned ipt = 0;
        this->get_wind_adv_diff(ipt, s, x, wind);
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << wind[i] << " ";
        }

        // Output stabilisation parameter
        outfile << Tau_SUPG << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }

    /// Output at default number of plot points
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// C-style output
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// C_style output at n_plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }


  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    double dshape_and_dtest_eulerian_adv_diff(const Vector<double>& s,
                                              Shape& psi,
                                              DShape& dpsidx,
                                              Shape& test,
                                              DShape& dtestdx) const;


    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    double dshape_and_dtest_eulerian_at_knot_adv_diff(const unsigned& ipt,
                                                      Shape& psi,
                                                      DShape& dpsidx,
                                                      Shape& test,
                                                      DShape& dtestdx) const;

    /// SUPG stabilisation parameter
    double Tau_SUPG;
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// \short Refineable version of QSUPGAdvectionDiffusionElement.
  /// Inherit from the standard QSUPGAdvectionDiffusionElement and the
  /// appropriate refineable geometric element and the refineable equations.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class RefineableQSUPGAdvectionDiffusionElement
    : public QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>,
      public virtual RefineableAdvectionDiffusionEquations<DIM>,
      public virtual RefineableQElement<DIM>
  {
  public:
    /// \short Constructor: Pass refinement level to refineable quad element
    /// (default 0 = root)
    RefineableQSUPGAdvectionDiffusionElement()
      : RefineableElement(),
        RefineableAdvectionDiffusionEquations<DIM>(),
        RefineableQElement<DIM>(),
        QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>()
    {
    }


    /// Broken copy constructor
    RefineableQSUPGAdvectionDiffusionElement(
      const RefineableQSUPGAdvectionDiffusionElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("RefineableQSUPGAdvectionDiffusionElement");
    }

    /// Broken assignment operator
    void operator=(
      const RefineableQSUPGAdvectionDiffusionElement<DIM, NNODE_1D>&)
    {
      BrokenCopy::broken_assign("RefineableQSUPGAdvectionDiffusionElement");
    }

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// \short Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>::nvertex_node();
    }

    /// \short Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QSUPGAdvectionDiffusionElement<DIM, NNODE_1D>::vertex_node_pt(j);
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// \short Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    ///  \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty.
    void further_setup_hanging_nodes() {}
  };


} // namespace oomph

#endif

// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// Header functions for functions used to generate orthogonal polynomials
// Include guards to prevent multiple inclusion of the header

#ifndef OOMPH_ORTHPOLY_HEADER
#define OOMPH_ORTHPOLY_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// Oomph-lib include
#include "Vector.h"
#include "oomph_utilities.h"

#include <cmath>

namespace oomph
{
  // Let's put these things in a namespace
  namespace Orthpoly
  {
    const double eps = 1.0e-15;

    /// Calculates Legendre polynomial of degree p at x
    /// using the three term recurrence relation
    /// \f$ (n+1) P_{n+1} = (2n+1)xP_{n} - nP_{n-1} \f$
    inline double legendre(const unsigned& p, const double& x)
    {
      // Return the constant value
      if (p == 0) return 1.0;
      // Return the linear polynomial
      else if (p == 1)
        return x;
      // Otherwise use the recurrence relation
      else
      {
        // Initialise the terms in the recurrence relation, we're going
        // to shift down before using the relation.
        double L0 = 1.0, L1 = 1.0, L2 = x;
        // Loop over the remaining polynomials
        for (unsigned n = 1; n < p; n++)
        {
          // Shift down the values
          L0 = L1;
          L1 = L2;
          // Use the three term recurrence
          L2 = ((2 * n + 1) * x * L1 - n * L0) / (n + 1);
        }
        // Once we've finished return the final value
        return L2;
      }
    }


    /// Calculates Legendre polynomial of degree p at x
    /// using three term recursive formula. Returns all polynomials up to
    /// order p in the vector
    inline void legendre_vector(const unsigned& p,
                                const double& x,
                                Vector<double>& polys)
    {
      // Set the constant term
      polys[0] = 1.0;
      // If we're only asked for the constant term return
      if (p == 0)
      {
        return;
      }
      // Set the linear polynomial
      polys[1] = x;
      // Initialise terms for the recurrence relation
      double L0 = 1.0, L1 = 1.0, L2 = x;
      // Loop over the remaining terms
      for (unsigned n = 1; n < p; n++)
      {
        // Shift down the values
        L0 = L1;
        L1 = L2;
        // Use the recurrence relation
        L2 = ((2 * n + 1) * x * L1 - n * L0) / (n + 1);
        // Set the newly calculated polynomial
        polys[n + 1] = L2;
      }
    }


    ///  Calculates first derivative of Legendre
    /// polynomial of degree p at x
    /// using three term recursive formula.
    /// \f$ nP_{n+1}^{'} = (2n+1)xP_{n}^{'} - (n+1)P_{n-1}^{'} \f$
    inline double dlegendre(const unsigned& p, const double& x)
    {
      double dL1 = 1.0, dL2 = 3 * x, dL3 = 0.0;
      if (p == 0) return 0.0;
      else if (p == 1)
        return dL1;
      else if (p == 2)
        return dL2;
      else
      {
        for (unsigned n = 2; n < p; n++)
        {
          dL3 = 1.0 / n * ((2.0 * n + 1.0) * x * dL2 - (n + 1.0) * dL1);
          dL1 = dL2;
          dL2 = dL3;
        }
        return dL3;
      }
    }

    /// Calculates second derivative of Legendre
    /// polynomial of degree p at x
    /// using three term recursive formula.
    inline double ddlegendre(const unsigned& p, const double& x)
    {
      double ddL2 = 3.0, ddL3 = 15 * x, ddL4 = 0.0;
      if (p == 0) return 0.0;
      else if (p == 1)
        return 0.0;
      else if (p == 2)
        return ddL2;
      else if (p == 3)
        return ddL3;
      else
      {
        for (unsigned i = 3; i < p; i++)
        {
          ddL4 =
            1.0 / (i - 1.0) * ((2.0 * i + 1.0) * x * ddL3 - (i + 2.0) * ddL2);
          ddL2 = ddL3;
          ddL3 = ddL4;
        }
        return ddL4;
      }
    }

    /// Calculate the Jacobi polnomials
    inline double jacobi(const int& alpha,
                         const int& beta,
                         const unsigned& p,
                         const double& x)
    {
      double P0 = 1.0;
      double P1 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
      double P2;
      if (p == 0) return P0;
      else if (p == 1)
        return P1;
      else
      {
        for (unsigned n = 1; n < p; n++)
        {
          double a1n =
            2 * (n + 1) * (n + alpha + beta + 1) * (2 * n + alpha + beta);
          double a2n =
            (2 * n + alpha + beta + 1) * (alpha * alpha - beta * beta);
          double a3n = (2 * n + alpha + beta) * (2 * n + alpha + beta + 1) *
                       (2 * n + alpha + beta + 2);
          double a4n =
            2 * (n + alpha) * (n + beta) * (2 * n + alpha + beta + 2);
          P2 = ((a2n + a3n * x) * P1 - a4n * P0) / a1n;
          P0 = P1;
          P1 = P2;
        }
        return P2;
      }
    }

    /// Calculate the Jacobi polnomials all in one goe
    inline void jacobi(const int& alpha,
                       const int& beta,
                       const unsigned& p,
                       const double& x,
                       Vector<double>& polys)
    {
      // Set the constant term
      polys[0] = 1.0;
      // If we've only been asked for the constant term, bin out
      if (p == 0)
      {
        return;
      }
      // Set the linear polynomial
      polys[1] = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
      // Initialise the terms for the recurrence relation
      double P0 = 1.0;
      double P1 = 1.0;
      double P2 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
      // Loop over the remaining terms
      for (unsigned n = 1; n < p; n++)
      {
        // Shift down the terms
        P0 = P1;
        P1 = P2;
        // Set the constants
        double a1n =
          2 * (n + 1) * (n + alpha + beta + 1) * (2 * n + alpha + beta);
        double a2n = (2 * n + alpha + beta + 1) * (alpha * alpha - beta * beta);
        double a3n = (2 * n + alpha + beta) * (2 * n + alpha + beta + 1) *
                     (2 * n + alpha + beta + 2);
        double a4n = 2 * (n + alpha) * (n + beta) * (2 * n + alpha + beta + 2);
        // Set the latest term
        P2 = ((a2n + a3n * x) * P1 - a4n * P0) / a1n;
        // Set the newly calculate polynomial
        polys[n + 1] = P2;
      }
    }


    /// Calculates the Gauss Lobatto Legendre abscissas for degree p = NNode-1
    void gll_nodes(const unsigned& Nnode, Vector<double>& x);

    // This version of gll_nodes calculates the abscissas AND weights
    void gll_nodes(const unsigned& Nnode, Vector<double>& x, Vector<double>& w);

    // Calculates the Gauss Legendre abscissas of degree p=Nnode-1
    void gl_nodes(const unsigned& Nnode, Vector<double>& x);

    // This version of gl_nodes calculates the abscissas AND weights
    void gl_nodes(const unsigned& Nnode, Vector<double>& x, Vector<double>& w);

  } // namespace Orthpoly

} // namespace oomph

#endif

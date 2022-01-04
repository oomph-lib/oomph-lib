// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#include "orthpoly.h"

namespace oomph
{
  namespace Orthpoly
  {
    // Calculates the Gauss Lobatto Legendre abscissas for degree p = NNode-1
    void gll_nodes(const unsigned& Nnode, Vector<double>& x)
    {
      double z, zold, del;
      unsigned p = Nnode - 1;
      x.resize(Nnode);
      if (Nnode < 2)
      {
        throw OomphLibError("Invalid number of nodes",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else if (Nnode == 2)
      {
        x[0] = -1.0;
        x[1] = 1.0;
      }
      else if (Nnode == 3)
      {
        x[0] = -1;
        x[1] = 0.0;
        x[2] = 1.0;
      }
      else if (Nnode == 4)
      {
        x[0] = -1.0;
        x[1] = -std::sqrt(5.0) / 5.0;
        x[2] = -x[1];
        x[3] = 1.0;
      }
      else
      {
        unsigned mid;
        if (p % 2 == 0)
        {
          mid = p / 2;
          x[mid] = 0.0;
        }
        else
          mid = p / 2 + 1;
        for (unsigned j = 1; j < mid; j++)
        {
          z = -std::cos(j * MathematicalConstants::Pi / double(p));
          do
          {
            del = dlegendre(p, z) / ddlegendre(p, z);
            zold = z;
            z = zold - del;
          } while (std::fabs(z - zold) > eps);
          x[j] = z;
          x[p - j] = -z;
        }
        x[0] = -1.0;
        x[p] = 1.0;
      }
    }

    // This version of gll_nodes calculates the abscissas AND weights
    void gll_nodes(const unsigned& Nnode, Vector<double>& x, Vector<double>& w)
    {
      gll_nodes(Nnode, x);
      // Now calculate the corresponding weights
      double l_z;
      unsigned p = Nnode - 1;
      for (unsigned i = 0; i < p + 1; i++)
      {
        l_z = legendre(p, x[i]);
        w[i] = 2.0 / (p * (p + 1) * l_z * l_z);
      }
    }


    // Calculates the Gauss Legendre abscissas of degree p=Nnode-1
    void gl_nodes(const unsigned& Nnode, Vector<double>& x)
    {
      double z, zold, del;
      unsigned p = Nnode - 1;
      x.resize(Nnode);
      if (Nnode < 2)
      {
        throw OomphLibError("Invalid number of nodes",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else if (Nnode == 2)
      {
        x[0] = -1.0 / 3.0 * std::sqrt(3.0);
        x[1] = -x[0];
      }
      else
      {
        unsigned mid;
        if (p % 2 == 0)
        {
          mid = p / 2;
          x[mid] = 0.0;
        }
        else
          mid = p / 2 + 1;
        for (unsigned j = 0; j < mid; j++)
        {
          z = -std::cos((2.0 * j + 1.0) * MathematicalConstants::Pi /
                        (2 * double(p) + 2.0));
          do
          {
            del = legendre(p + 1, z) / dlegendre(p + 1, z);
            zold = z;
            z = zold - del;
          } while (std::fabs(z - zold) > eps);
          x[j] = z;
          x[p - j] = -z;
        }
      }
    }

    // This version of gl_nodes calculates the abscissas AND weights
    void gl_nodes(const unsigned& Nnode, Vector<double>& x, Vector<double>& w)
    {
      gl_nodes(Nnode, x);
      // Now calculate the corresponding weights
      double dl_z;
      unsigned p = Nnode - 1;
      for (unsigned i = 0; i < p + 1; i++)
      {
        dl_z = dlegendre(p + 1, x[i]);
        w[i] = 2.0 / (1 - x[i] * x[i]) / (dl_z * dl_z);
      }
    }
  } // namespace Orthpoly

} // namespace oomph

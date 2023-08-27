//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef QINTEGRAL_ELEMENTS_HEADER
#define QINTEGRAL_ELEMENTS_HEADER

#include "../generic/elements.h"
#include "../generic/Qelements.h"
#include "integral_equations.h"

namespace oomph
{
  template<unsigned DIM, unsigned NNODE_1D>
  class QIntegralElement : public virtual QElement<DIM, NNODE_1D>,
                           public virtual IntegralEquations<DIM>
  {
  private:
    /// Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    /// Constructor: Call constructors for QElement and
    /// Poisson equations
    QIntegralElement() : QElement<DIM, NNODE_1D>(), IntegralEquations<DIM>() {}

    ///  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    virtual void output(std::ostream& outfile)
    {
      output(outfile, NNODE_1D);
    }

    void output(std::ostream& outfile, const unsigned& nplot)
    {
      IntegralEquations<DIM>::output(outfile, nplot);
    }
  };

  //=======================================================================
  /// Face geometry for the QIntegralElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QIntegralElement<DIM, NNODE_1D>>
    : public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QIntegralElement<DIM, NNODE_1D>::Initial_Nvalue = 0;
} // namespace oomph
#endif

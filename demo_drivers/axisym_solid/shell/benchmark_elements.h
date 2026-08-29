// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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


// OOMPH-LIB headers
#include "generic/nodes.h"
#include "generic/Qelements.h"
#include "generic/oomph_utilities.h"


namespace oomph
{
  //=============================================================
  /// A class for elements to solve the 1D radial displacement ODE
  /// in the spherical shell compression benchmark problem.
  //=============================================================
  template<unsigned N_NODE1D>
  class ShellBenchmarkElement : public virtual QElement<1, N_NODE1D>,
                                public virtual FiniteElement
  {
  public:
    /// Constructor ()
    ShellBenchmarkElement() : Nu_pt(0) {}

    /// Broken copy constructor
    ShellBenchmarkElement(const ShellBenchmarkElement& dummy) = delete;

    /// Broken assignment operator
    void operator=(const ShellBenchmarkElement&) = delete;

    ///  Required  # of `values' (pinned or dofs) at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return 1;
    }

    /// Access function: Poisson ratio
    const double& nu() const
    {
      return *Nu_pt;
    }

    /// Pointer to Poisson ratio
    double*& nu_pt()
    {
      return Nu_pt;
    }

    /// Output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_shell_benchmark(residuals);
    }

    /// Return FE representation of function value u at local coordinate s
    virtual inline double interpolated_u_shell_benchmark(
      const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the poisson unknown is stored
      const unsigned u_nodal_index = 0; // u_index_shell_benchmark();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }

  protected:
    /// Compute element residual Vector only (if flag=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_shell_benchmark(
      Vector<double>& residuals);

    /// Pointer to the Poisson ratio
    double* Nu_pt;
  };


  //======================================================================
  /// Face elements for the application of an external pressure
  //======================================================================
  class ShellBenchmarkPressureElement : public virtual PointElement,
                                        public virtual FaceElement
  {
  public:
    /// Constructor, takes the pointer to the "bulk" element and the
    /// index of the face to be created
    ShellBenchmarkPressureElement(FiniteElement* const& bulk_el_pt,
                                  const int& face_index)
    {
      // Let the bulk element build the FaceElement, i.e. setup the pointers
      // to its nodes (by referring to the appropriate nodes in the bulk
      // element), etc.
      bulk_el_pt->build_face_element(face_index, this);
    }

    /// Broken copy constructor
    ShellBenchmarkPressureElement(const ShellBenchmarkPressureElement& dummy) =
      delete;

    /// Compute the element residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_shell_benchmark_pressure(residuals);
    }

    /// Access function: Poisson ratio
    const double& pressure() const
    {
      return *Pressure_pt;
    }

    /// Pointer to Poisson ratio
    double*& pressure_pt()
    {
      return Pressure_pt;
    }

  protected:
    /// Compute the element residual vector.
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_shell_benchmark_pressure(
      Vector<double>& residuals);

    /// Pointer to the applied pressure
    double* Pressure_pt;

    /// The index at which the unknown is stored at the nodes
    unsigned U_index_shell_benchmark_pressure;
  };
} // namespace oomph

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
// Header file for specific surface elements

#ifndef OOMPH_NAVIER_STOKES_SURFACE_POWER_ELEMENTS_HEADER
#define OOMPH_NAVIER_STOKES_SURFACE_POWER_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "../generic/Qelements.h"

namespace oomph
{
  //======================================================================
  /// A class of elements that allow the determination of the power
  /// input and various other fluxes over the domain boundaries.
  /// The element operates as a FaceElement and attaches itself
  /// to a bulk element of the type specified by the template
  /// argument.
  //======================================================================
  template<class ELEMENT>
  class NavierStokesSurfacePowerElement : public virtual FaceGeometry<ELEMENT>,
                                          public virtual FaceElement
  {
  public:
    /// Constructor, which takes a "bulk" element and the value of the index
    /// and its limit
    NavierStokesSurfacePowerElement(FiniteElement* const& element_pt,
                                    const int& face_index)
      : FaceGeometry<ELEMENT>(), FaceElement()
    {
      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Set the dimension from the dimension of the first node
      Dim = node_pt(0)->ndim();
    }


    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    /// This final over-ride is required for cases where the
    /// FaceElement is a SolidFiniteElement because both SolidFiniteElements
    /// and FaceElements overload zeta_nodal.
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }


    /// Get drag force (traction acting on fluid)
    Vector<double> drag_force()
    {
      std::ofstream dummy_file;
      return drag_force(dummy_file);
    }


    /// Get  drag force (traction acting on fluid)
    /// Doc in outfile.
    Vector<double> drag_force(std::ofstream& outfile)
    {
      // Spatial dimension of element
      unsigned ndim = dim();

      // Initialise
      Vector<double> drag(ndim + 1, 0.0);

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();


      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Hacky: This is only appropriate for 3 point integration of
      // 1D line elements
      if (outfile.is_open()) outfile << "ZONE I=3" << std::endl;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);

        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        // Get traction from bulk element
        Vector<double> traction(ndim + 1);
        bulk_el_pt->get_traction(s_bulk, normal, traction);

        // Integrate
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          drag[i] += traction[i] * W;
        }

        if (outfile.is_open())
        {
          // Output x,y,...,
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << x[i] << " ";
          }

          // Output traction
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << traction[i] << " ";
          }


          // Output normal
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << normal[i] << " ";
          }

          outfile << std::endl;
        }
      }
      return drag;
    }

    /// Get integral of instantaneous rate of work done by
    /// the traction that's exerted onto the fluid.
    double get_rate_of_traction_work()
    {
      std::ofstream dummy_file;
      return get_rate_of_traction_work(dummy_file);
    }


    /// Get integral of instantaneous rate of work done by
    /// the traction that's exerted onto the fluid. Doc in outfile.
    double get_rate_of_traction_work(std::ofstream& outfile)
    {
      // Initialise
      double rate_of_work_integral = 0.0;

      // Spatial dimension of element
      unsigned ndim = dim();

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();


      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Hacky: This is only appropriate for 3x3 integration of
      // 2D quad elements
      if (outfile.is_open()) outfile << "ZONE I=3, J=3" << std::endl;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);


        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        // Get traction from bulk element
        Vector<double> traction(ndim + 1);
        bulk_el_pt->get_traction(s_bulk, normal, traction);


        // Local rate of work:
        double rate_of_work = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          rate_of_work += traction[i] * veloc[i];
        }

        // Add rate of work
        rate_of_work_integral += rate_of_work * W;

        if (outfile.is_open())
        {
          // Output x,y,...,
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << x[i] << " ";
          }

          // Output traction
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << traction[i] << " ";
          }

          // Output veloc
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << veloc[i] << " ";
          }

          // Output normal
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << normal[i] << " ";
          }

          // Output local rate of work
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << rate_of_work << " ";
          }

          outfile << std::endl;
        }
      }

      return rate_of_work_integral;
    }


    /// Get integral of instantaneous rate of work done by
    /// the traction that's exerted onto the fluid, decomposed into pressure
    /// and normal and tangential viscous components.
    void get_rate_of_traction_work_components(double& rate_of_work_integral_p,
                                              double& rate_of_work_integral_n,
                                              double& rate_of_work_integral_t)
    {
      std::ofstream dummy_file;
      get_rate_of_traction_work_components(dummy_file,
                                           rate_of_work_integral_p,
                                           rate_of_work_integral_n,
                                           rate_of_work_integral_t);
    }


    /// Get integral of instantaneous rate of work done by
    /// the traction that's exerted onto the fluid, decomposed into pressure
    /// and normal and tangential viscous components.  Doc in outfile.
    void get_rate_of_traction_work_components(std::ofstream& outfile,
                                              double& rate_of_work_integral_p,
                                              double& rate_of_work_integral_n,
                                              double& rate_of_work_integral_t)
    {
      // Initialise
      rate_of_work_integral_p = 0;
      rate_of_work_integral_n = 0;
      rate_of_work_integral_t = 0;

      // Spatial dimension of element
      unsigned ndim = dim();

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();


      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Hacky: This is only appropriate for 3x3 integration of
      // 2D quad elements
      if (outfile.is_open()) outfile << "ZONE I=3, J=3" << std::endl;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);


        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        // Get traction from bulk element
        Vector<double> traction_p(ndim + 1);
        Vector<double> traction_n(ndim + 1);
        Vector<double> traction_t(ndim + 1);
        bulk_el_pt->get_traction(
          s_bulk, normal, traction_p, traction_n, traction_t);


        // Local rate of work:
        double rate_of_work_p = 0.0;
        double rate_of_work_n = 0.0;
        double rate_of_work_t = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          rate_of_work_p += traction_p[i] * veloc[i];
          rate_of_work_n += traction_n[i] * veloc[i];
          rate_of_work_t += traction_t[i] * veloc[i];
        }

        // Add rate of work
        rate_of_work_integral_p += rate_of_work_p * W;
        rate_of_work_integral_n += rate_of_work_n * W;
        rate_of_work_integral_t += rate_of_work_t * W;

        if (outfile.is_open())
        {
          // Output x,y,...,
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << x[i] << " ";
          }

          // Output traction due to pressure
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << traction_p[i] << " ";
          }

          // Output traction due to viscous normal stress
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << traction_n[i] << " ";
          }

          // Output traction due to viscous tangential stress
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << traction_t[i] << " ";
          }

          // Output veloc
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << veloc[i] << " ";
          }

          // Output normal
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << normal[i] << " ";
          }

          // Output local rate of work due to pressure
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << rate_of_work_p << " ";
          }

          // Output local rate of work due to viscous normal stress
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << rate_of_work_n << " ";
          }

          // Output local rate of work due to viscous tangential stress
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << rate_of_work_t << " ";
          }

          outfile << std::endl;
        }
      }
    }


    /// Get integral of kinetic energy flux
    double get_kinetic_energy_flux()
    {
      std::ofstream dummy_file;
      return get_kinetic_energy_flux(dummy_file);
    }


    /// Get integral of kinetic energy flux and doc
    double get_kinetic_energy_flux(std::ofstream& outfile)
    {
      // Initialise
      double kinetic_energy_flux_integral = 0.0;

      // Spatial dimension of element
      unsigned ndim = dim();

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Hacky: This is only appropriate for 3x3 integration of
      // 2D quad elements
      if (outfile.is_open()) outfile << "ZONE I=3, J=3" << std::endl;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);

        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        double kin_energy = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          kin_energy += veloc[i] * veloc[i];
        }
        kin_energy *= 0.5;

        // Kinetic energy flux
        double kin_energy_flux = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          kin_energy_flux += kin_energy * normal[i] * veloc[i];
        }

        // Add to integral
        kinetic_energy_flux_integral += kin_energy_flux * W;

        if (outfile.is_open())
        {
          // Output x,y,...,
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << x[i] << " ";
          }

          // Output veloc
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << veloc[i] << " ";
          }

          // Output normal
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << normal[i] << " ";
          }

          // Output local kin energy flux
          outfile << kin_energy_flux << " ";

          outfile << std::endl;
        }
      }

      return kinetic_energy_flux_integral;
    }


    /// Get integral of volume flux
    double get_volume_flux()
    {
      std::ofstream dummy_file;
      return get_volume_flux(dummy_file);
    }


    /// Get integral of volume flux and doc
    double get_volume_flux(std::ofstream& outfile)
    {
      // Initialise
      double volume_flux_integral = 0.0;

      // Spatial dimension of element
      unsigned ndim = dim();

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Hacky: This is only appropriate for 3x3 integration of
      // 2D quad elements
      if (outfile.is_open()) outfile << "ZONE I=3, J=3" << std::endl;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }


        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);

        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        // Volume flux
        double volume_flux = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          volume_flux += normal[i] * veloc[i];
        }

        // Add to integral
        volume_flux_integral += volume_flux * W;

        if (outfile.is_open())
        {
          // Output x,y,...,
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << x[i] << " ";
          }

          // Output veloc
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << veloc[i] << " ";
          }

          // Output normal
          for (unsigned i = 0; i < ndim + 1; i++)
          {
            outfile << normal[i] << " ";
          }

          // Output local volume flux
          outfile << volume_flux << " ";

          outfile << std::endl;
        }
      }

      return volume_flux_integral;
    }


  private:
    /// The highest dimension of the problem
    unsigned Dim;
  };


} // namespace oomph

#endif

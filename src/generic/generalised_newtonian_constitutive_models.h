// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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


#ifndef OOMPH_GENERALISED_NEWTONIAN_CONSTITUTIVE_MODELS_HEADER
#define OOMPH_GENERALISED_NEWTONIAN_CONSTITUTIVE_MODELS_HEADER


// Oomph-lib includes
//#include "generic.h"

namespace oomph
{
  //======================================================================
  /// A Base class defining the generalise Newtonian constitutive relation
  //======================================================================
  template<unsigned DIM>
  class GeneralisedNewtonianConstitutiveEquation
  {
  public:
    /// Empty constructor
    GeneralisedNewtonianConstitutiveEquation() {}


    /// Empty virtual destructor
    virtual ~GeneralisedNewtonianConstitutiveEquation() {}

    /// Function implementing the constitutive model
    /// Input: second invariant of the rate of strain
    /// Output: the viscosity
    /// For Newtonian behaviour this returns 1
    virtual double viscosity(
      const double& second_invariant_of_rate_of_strain_tensor) = 0;

    /// Function returning the derivative of the viscosity w.r.t.
    /// the second invariant of the rate of strain tensor
    /// For Newtonian behaviour this returns 0.0
    virtual double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor) = 0;
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Newtonian fluid
  //===================================================================
  template<unsigned DIM>
  class NewtonianConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  public:
    /// Constructor: specify viscosity ratio (defaults to one)
    NewtonianConstitutiveEquation(const double& viscosity_ratio = 1.0)
      : Viscosity_ratio(viscosity_ratio)
    {
    }

    /// in the Newtonian case the viscosity is constant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      return Viscosity_ratio;
    }

    /// the derivative w.r.t. I2 is zero
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      return 0.0;
    }

  private:
    /// Viscosity ratio
    double Viscosity_ratio;
  };


  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class defining a power-law
  /// fluid regularised according to Bercovier and Engelman (1980) to allow for
  /// n < 1
  //==================================================================
  template<unsigned DIM>
  class PowerLawBerEngRegConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// power law index n
    double* Power_pt;

    /// regularisation parameter e << 1
    double* Regularisation_parameter_pt;


  public:
    PowerLawBerEngRegConstitutiveEquation(double* power_pt, double* reg_par_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Power_pt(power_pt),
        Regularisation_parameter_pt(reg_par_pt)
    {
    }

    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // Calculate the square root of the absolute value of the
      // second invariant of the rate of strain tensor
      double measure_of_rate_of_strain =
        sqrt(sign * second_invariant_of_rate_of_strain_tensor);

      return pow(
        (2.0 * measure_of_rate_of_strain + *Regularisation_parameter_pt),
        *Power_pt - 1.0);
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Herschel-Bulkley fluid
  /// using Bercovier and Engelman's (1980) regularisation
  //==================================================================
  template<unsigned DIM>
  class HerschelBulkleyBerEngRegConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// yield stress tau_y
    double* Yield_stress_pt;

    /// power law index n
    double* Flow_index_pt;

    /// regularisation parameter e << 1
    double* Regularisation_parameter_pt;


  public:
    HerschelBulkleyBerEngRegConstitutiveEquation(
      double* yield_stress_pt,
      double* flow_index_pt,
      double* regularisation_parameter_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Flow_index_pt(flow_index_pt),
        Regularisation_parameter_pt(regularisation_parameter_pt)
    {
    }

    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // Calculate the square root of the absolute value of the
      // second invariant of the rate of strain tensor
      double measure_of_rate_of_strain =
        sqrt(sign * second_invariant_of_rate_of_strain_tensor +
             (*Regularisation_parameter_pt));

      return (*Yield_stress_pt) / (2.0 * measure_of_rate_of_strain) +
             pow(2.0 * measure_of_rate_of_strain, *Flow_index_pt - 1.0);
    }

    // not implemented yet
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // std::ostringstream error_stream;
      // error_stream << "This has not been implemented yet!";
      // throw OomphLibError(
      //  error_stream.str(),
      //  OOMPH_CURRENT_FUNCTION,
      //  OOMPH_EXCEPTION_LOCATION);

      return sign * pow(2.0, (*Flow_index_pt) - 2.0) *
               ((*Flow_index_pt) - 1.0) *
               pow(sign * second_invariant_of_rate_of_strain_tensor +
                     (*Regularisation_parameter_pt),
                   ((*Flow_index_pt) - 1.0) / 2.0 - 1.0) -
             sign * (*Yield_stress_pt) /
               (4.0 * pow(sign * second_invariant_of_rate_of_strain_tensor +
                            (*Regularisation_parameter_pt),
                          3.0 / 2.0));
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Herschel-Bulkley fluid
  /// using Tanner and Milthorpe's (1983) regularisation
  //==================================================================
  template<unsigned DIM>
  class HerschelBulkleyTanMilRegConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// yield stress tau_y
    double* Yield_stress_pt;

    /// power law index n
    double* Flow_index_pt;

    /// value of the second invariant below which we have
    /// constant (Newtonian) viscosity -- assumed to be always positive
    double* Critical_second_invariant_pt;

  public:
    /// "Cutoff regularised" Herschel Bulkley constitutive equation
    HerschelBulkleyTanMilRegConstitutiveEquation(
      double* yield_stress_pt,
      double* flow_index_pt,
      double* critical_second_invariant_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Flow_index_pt(flow_index_pt),
        Critical_second_invariant_pt(critical_second_invariant_pt)
    {
      // Calculate the Newtonian cutoff viscosity from the constitutive
      // equation and the cutoff value of the second invariant
      double cut_off_viscosity = calculate_cut_off_viscosity();

      oomph_info << "HerschelBulkleyTanMilRegConstitutiveEquation: "
                 << " cutoff viscosity = " << cut_off_viscosity << std::endl;

      oomph_info << "                                              "
                 << " cutoff invariant = " << *Critical_second_invariant_pt
                 << std::endl;
    }

    /// Function that calculates the cut off viscosity
    double calculate_cut_off_viscosity()
    {
      return (*Yield_stress_pt) /
               (2.0 * sqrt((*Critical_second_invariant_pt))) +
             pow((2.0 * sqrt((*Critical_second_invariant_pt))),
                 *Flow_index_pt - 1.0);
    }

    /// Report cutoff values
    void report_cut_off_values(double& cut_off_invariant,
                               double& cut_off_viscosity)
    {
      cut_off_invariant = *Critical_second_invariant_pt;
      cut_off_viscosity = calculate_cut_off_viscosity();
    }


    /// Viscosity ratio as a fct of strain rate invariant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // std::cout<<"I2 "<<second_invariant_of_rate_of_strain_tensor<<std::endl;

      // if the second invariant is below the cutoff we have a constant,
      // Newtonian viscosity
      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return calculate_cut_off_viscosity();
      }
      else
      {
        // Calculate the square root of the absolute value of the
        // second invariant of the rate of strain tensor
        double measure_of_rate_of_strain =
          sqrt(sign * second_invariant_of_rate_of_strain_tensor);

        return (*Yield_stress_pt) / (2.0 * measure_of_rate_of_strain) +
               pow((2.0 * measure_of_rate_of_strain), *Flow_index_pt - 1.0);
      }
    }

    /// Deriv of viscosity w.r.t. strain rate invariant
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return 0.0;
      }
      else
      {
        return sign * pow(2.0, (*Flow_index_pt) - 2.0) *
                 ((*Flow_index_pt) - 1.0) *
                 pow(sign * second_invariant_of_rate_of_strain_tensor,
                     ((*Flow_index_pt) - 1.0) / 2.0 - 1.0) -
               sign * (*Yield_stress_pt) /
                 (4.0 * pow(sign * second_invariant_of_rate_of_strain_tensor,
                            3.0 / 2.0));
      }
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Herschel-Bulkley fluid
  /// using Tanner and Milthorpe's (1983) regularisation
  /// with a smooth transition using a quadratic
  //==================================================================
  template<unsigned DIM>
  class HerschelBulkleyTanMilRegWithBlendingConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// yield stress tau_y
    double* Yield_stress_pt;

    /// power law index n
    double* Flow_index_pt;

    /// value of the second invariant below which we have
    /// constant (Newtonian) viscosity -- assumed to be always positive
    double* Critical_second_invariant_pt;

    /// We use a quadratic function to smoothly blend from the Herschel Bulkley
    /// model at the cut-off to a constant viscosity as the strain rate
    /// /approaches zero; this way we avoid the
    /// discontinuity of the gradient at the cut-off, which is present in the
    /// classic Tanner Milthorpe regularisation
    double a;
    double b;
    double c;

    /// Fraction of the cut-off strain rate below which the viscosity is
    /// constant 0 <= \alpha < 1
    double alpha;


  public:
    /// "Cutoff regularised" Herschel Bulkley constitutive equation
    HerschelBulkleyTanMilRegWithBlendingConstitutiveEquation(
      double* yield_stress_pt,
      double* flow_index_pt,
      double* critical_second_invariant_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Flow_index_pt(flow_index_pt),
        Critical_second_invariant_pt(critical_second_invariant_pt),
        a(0.0),
        b(0.0),
        c(0.0)
    {
      /// Blend over one order of magnitude
      alpha = 0.1;

      /// Calculate the coefficients of the quadratic equation
      if (fabs(*Critical_second_invariant_pt) > 0.0)
      {
        a = (pow(2.0, *Flow_index_pt - 3.0) * (*Flow_index_pt - 1.0) *
               pow(fabs(*Critical_second_invariant_pt),
                   (*Flow_index_pt - 1.0) / 2.0 - 2.0) -
             (*Yield_stress_pt) /
               (8.0 * pow(fabs(*Critical_second_invariant_pt), 5.0 / 2.0))) /
            (1.0 - alpha);
        b = -2.0 * a * alpha * fabs(*Critical_second_invariant_pt);
        c = (*Yield_stress_pt) /
              (2.0 * sqrt(fabs(*Critical_second_invariant_pt))) +
            pow(2.0 * sqrt(fabs(*Critical_second_invariant_pt)),
                *Flow_index_pt - 1.0) -
            a * pow(fabs(*Critical_second_invariant_pt), 2.0) -
            b * fabs(*Critical_second_invariant_pt);
      }

      /// get the cutoff viscosity
      double cut_off_viscosity = calculate_cutoff_viscosity();

      /// get the zero shear viscosity
      double zero_shear_viscosity = calculate_zero_shear_viscosity();

      oomph_info << "HerschelBulkleyTanMilRegWithBlendingConstitutiveEquation: "
                 << " zero shear viscosity = " << zero_shear_viscosity
                 << std::endl;

      oomph_info << "HerschelBulkleyTanMilRegWithBlendingConstitutiveEquation: "
                 << " cut off viscosity = " << cut_off_viscosity << std::endl;

      oomph_info << "                                              "
                 << " cutoff invariant = " << *Critical_second_invariant_pt
                 << std::endl;
    }

    /// Function that calculates the viscosity at the cutoff invariant
    /// Note: this is NOT the viscosity at zero I2
    double calculate_cutoff_viscosity()
    {
      // Calculate the Newtonian cutoff viscosity from the constitutive
      // equation and the cutoff value of the second invariant
      return (*Yield_stress_pt) /
               (2.0 * sqrt(std::max(fabs(*Critical_second_invariant_pt),
                                    DBL_EPSILON))) +
             pow((2.0 * sqrt(std::max(fabs(*Critical_second_invariant_pt),
                                      DBL_EPSILON))),
                 *Flow_index_pt - 1.0);
    }

    /// Function that calculates the viscosity at zero I2
    double calculate_zero_shear_viscosity()
    {
      // The maximum of either the viscosity at alpha*cut-off or the cut-off
      // viscosity for a cut-off of zero
      return std::max(a * pow(alpha, 2.0) *
                          pow(fabs(*Critical_second_invariant_pt), 2.0) +
                        b * alpha * fabs(*Critical_second_invariant_pt) + c,
                      calculate_cutoff_viscosity());
    }

    /// Report cutoff values
    void report_cut_off_values(double& cut_off_invariant,
                               double& cut_off_viscosity,
                               double& zero_shear_viscosity)
    {
      cut_off_invariant = *Critical_second_invariant_pt;
      cut_off_viscosity = calculate_cutoff_viscosity();
      zero_shear_viscosity = calculate_zero_shear_viscosity();
    }

    /// Viscosity ratio as a fct of strain rate invariant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      if (fabs(second_invariant_of_rate_of_strain_tensor) <
          alpha * fabs(*Critical_second_invariant_pt))
      {
        return calculate_zero_shear_viscosity();
      }
      else if (fabs(second_invariant_of_rate_of_strain_tensor) <
               fabs(*Critical_second_invariant_pt))
      {
        return a * pow(fabs(second_invariant_of_rate_of_strain_tensor), 2.0) +
               b * fabs(second_invariant_of_rate_of_strain_tensor) + c;
      }

      return (*Yield_stress_pt) /
               (2.0 * sqrt(fabs(second_invariant_of_rate_of_strain_tensor))) +
             pow(2.0 * sqrt(fabs(second_invariant_of_rate_of_strain_tensor)),
                 *Flow_index_pt - 1.0);
    }

    /// Deriv of viscosity w.r.t. strain rate invariant
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      if (fabs(second_invariant_of_rate_of_strain_tensor) <
          alpha * fabs(*Critical_second_invariant_pt))
      {
        return 0.0;
      }
      else if (fabs(second_invariant_of_rate_of_strain_tensor) <
               fabs(*Critical_second_invariant_pt))
      {
        return 2.0 * a * fabs(second_invariant_of_rate_of_strain_tensor) + b;
      }

      return pow(2.0, *Flow_index_pt - 2.0) * (*Flow_index_pt - 1.0) *
               pow(fabs(second_invariant_of_rate_of_strain_tensor),
                   (*Flow_index_pt - 1.0) / 2.0 - 1.0) -
             (*Yield_stress_pt) /
               (4.0 * pow(fabs(second_invariant_of_rate_of_strain_tensor),
                          3.0 / 2.0));
    }
  };


  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Herschel-Bulkley fluid
  /// using Papanastasiou's (1987) regularisation
  //==================================================================
  template<unsigned DIM>
  class HerschelBulkleyPapRegConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// Yield stress tau_y
    double* Yield_stress_pt;

    /// Power law index n
    double* Flow_index_pt;

    /// Regularisation parameter m >> 1
    double* Exponential_parameter_pt;


  public:
    HerschelBulkleyPapRegConstitutiveEquation(double* yield_stress_pt,
                                              double* flow_index_pt,
                                              double* exponential_parameter_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Flow_index_pt(flow_index_pt),
        Exponential_parameter_pt(exponential_parameter_pt)
    {
    }

    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Calculate magnitude of the rate of strain tensor
      double measure_of_rate_of_strain =
        sqrt(std::fabs(second_invariant_of_rate_of_strain_tensor));

      return (1.0 - exp(-2.0 * (*Exponential_parameter_pt) *
                        measure_of_rate_of_strain)) /
               (2.0 * measure_of_rate_of_strain) * (*Yield_stress_pt) +
             (1.0 - exp(-2.0 * (*Exponential_parameter_pt) *
                        measure_of_rate_of_strain)) /
               (2.0 * measure_of_rate_of_strain) *
               pow((2.0 * measure_of_rate_of_strain), *Flow_index_pt);
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Herschel-Bulkley fluid
  /// using Mendes and Dutra's (2004) regularisation
  //==================================================================
  template<unsigned DIM>
  class HerschelBulkleyMenDutRegConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// yield stress tau_y
    double* Yield_stress_pt;

    /// power law index n
    double* Flow_index_pt;

    /// the viscosity at zero shear rate
    double* Zero_shear_viscosity_pt;

  public:
    /// "Exponentially regularised" Herschel Bulkley constitutive equation
    HerschelBulkleyMenDutRegConstitutiveEquation(
      double* yield_stress_pt,
      double* flow_index_pt,
      double* zero_shear_viscosity_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Flow_index_pt(flow_index_pt),
        Zero_shear_viscosity_pt(zero_shear_viscosity_pt)
    {
    }

    /// Viscosity ratio as a fct of strain rate invariant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      // Also, because the viscosity is exactly zero at zero invariant,
      // we have to add a small value to it
      double sign = -1.0;
      double eps = 0.0;
      if (second_invariant_of_rate_of_strain_tensor == 0.0)
      {
        eps = 1.0e-30;
        sign = 1.0;
      }
      else if (second_invariant_of_rate_of_strain_tensor > 0.0)
      {
        sign = 1.0;
      }

      // Calculate the square root of the absolute value of the
      // second invariant of the rate of strain tensor
      double measure_of_rate_of_strain =
        sqrt(sign * (second_invariant_of_rate_of_strain_tensor + eps));

      return (1.0 - exp(-(*Zero_shear_viscosity_pt) * 2.0 *
                        measure_of_rate_of_strain / (*Yield_stress_pt))) *
             ((*Yield_stress_pt) / (2.0 * measure_of_rate_of_strain) +
              pow((2.0 * measure_of_rate_of_strain), *Flow_index_pt - 1.0));
    }

    /// Deriv of viscosity w.r.t. strain rate invariant
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      // Also, because the viscosity is exactly zero at zero invariant,
      // we have to add a small value to it
      double sign = -1.0;
      double eps = 0.0;
      if (second_invariant_of_rate_of_strain_tensor == 0.0)
      {
        eps = 1.0e-30;
        sign = 1.0;
      }
      else if (second_invariant_of_rate_of_strain_tensor > 0.0)
      {
        sign = 1.0;
      }

      // Calculate the square root of the absolute value of the
      // second invariant of the rate of strain tensor
      double measure_of_rate_of_strain =
        sqrt(sign * (second_invariant_of_rate_of_strain_tensor + eps));

      return (1.0 - exp(-(*Zero_shear_viscosity_pt) * 2.0 *
                        measure_of_rate_of_strain / (*Yield_stress_pt))) *
               (sign * pow(2.0, (*Flow_index_pt) - 2.0) *
                  ((*Flow_index_pt) - 1.0) *
                  pow(sign * (second_invariant_of_rate_of_strain_tensor + eps),
                      ((*Flow_index_pt) - 1.0) / 2.0 - 1.0) -
                sign * (*Yield_stress_pt) /
                  (4.0 *
                   pow(sign * (second_invariant_of_rate_of_strain_tensor + eps),
                       3.0 / 2.0))) +
             sign *
               (((*Zero_shear_viscosity_pt) *
                 exp(-(*Zero_shear_viscosity_pt) * 2.0 *
                     measure_of_rate_of_strain / (*Yield_stress_pt))) /
                ((*Yield_stress_pt) * measure_of_rate_of_strain)) *
               (pow(2.0, (*Flow_index_pt) - 1.0) *
                  pow(sign * (second_invariant_of_rate_of_strain_tensor + eps),
                      ((*Flow_index_pt) - 1.0) / 2.0) +
                (*Yield_stress_pt) / (2.0 * measure_of_rate_of_strain));
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Sisko fluid
  /// using Tanner and Milthorpe's (1983) regularisation
  /// with a smooth transition using a cubic (for n < 1)
  //==================================================================
  template<unsigned DIM>
  class SiskoTanMilRegWithBlendingConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// pre-factor alpha
    double* Alpha_pt;

    /// power law index n
    double* Flow_index_pt;

    /// value of the second invariant below which we have
    /// constant (Newtonian) viscosity -- assumed to be always positive
    double* Critical_second_invariant_pt;


  public:
    /// "Cutoff regularised" Sisko constitutive equation
    SiskoTanMilRegWithBlendingConstitutiveEquation(
      double* alpha_pt,
      double* flow_index_pt,
      double* critical_second_invariant_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Alpha_pt(alpha_pt),
        Flow_index_pt(flow_index_pt),
        Critical_second_invariant_pt(critical_second_invariant_pt)
    {
      /// get the cutoff viscosity
      double cut_off_viscosity = calculate_cutoff_viscosity();

      /// get the zero shear viscosity
      double zero_shear_viscosity = calculate_zero_shear_viscosity();

      oomph_info << "SiskoTanMilRegWithBlendingConstitutiveEquation: "
                 << " zero shear viscosity = " << zero_shear_viscosity
                 << std::endl;

      oomph_info << "SiskoTanMilRegWithBlendingConstitutiveEquation: "
                 << " cut off viscosity = " << cut_off_viscosity << std::endl;

      oomph_info << "                                              "
                 << " cutoff invariant = " << *Critical_second_invariant_pt
                 << std::endl;
    }

    /// Function that calculates the viscosity at the cutoff invariant
    /// Note: this is NOT the viscosity at zero I2
    double calculate_cutoff_viscosity()
    {
      // Calculate the Newtonian cutoff viscosity from the constitutive
      // equation and the cutoff value of the second invariant
      return 1.0 +
             (*Alpha_pt) * pow((2.0 * sqrt((*Critical_second_invariant_pt))),
                               *Flow_index_pt - 1.0);
    }

    /// Offset by how much the zero shear rate viscosity lies
    /// above the viscosity at I2_cutoff
    /// Hard-coded to a value that ensures a smooth transition
    double calculate_viscosity_offset_at_zero_shear(double& cut_off_viscosity)
    {
      return cut_off_viscosity / 5.0;
    }

    /// Function that calculates the viscosity at zero I2
    double calculate_zero_shear_viscosity()
    {
      // get the viscosity at the cutoff point
      double cut_off_viscosity = calculate_cutoff_viscosity();

      /// Offset by how much the zero shear rate viscosity lies
      /// above the viscosity at I2_cutoff
      double epsilon =
        calculate_viscosity_offset_at_zero_shear(cut_off_viscosity);

      return cut_off_viscosity + epsilon;
    }

    /// Report cutoff values
    void report_cut_off_values(double& cut_off_invariant,
                               double& cut_off_viscosity,
                               double& zero_shear_viscosity)
    {
      cut_off_invariant = *Critical_second_invariant_pt;
      cut_off_viscosity = calculate_cutoff_viscosity();
      zero_shear_viscosity = calculate_zero_shear_viscosity();
    }

    // Calculate the fitting parameters for the cubic blending
    void calculate_fitting_parameters_of_cubic(double& a, double& b)
    {
      // get the viscosity at the cutoff invariant
      double Cut_off_viscosity = calculate_cutoff_viscosity();

      // calculate the offset at zero shear
      double epsilon =
        calculate_viscosity_offset_at_zero_shear(Cut_off_viscosity);

      a =
        1.0 /
        (4.0 * pow((*Critical_second_invariant_pt), 3.0) *
         pow((*Critical_second_invariant_pt), 5.0 / 2.0)) *
        (8.0 * (Cut_off_viscosity + epsilon - 1.0) *
           pow((*Critical_second_invariant_pt), 5.0 / 2.0) +
         pow(2.0, (*Flow_index_pt)) * ((*Flow_index_pt) - 1.0) * (*Alpha_pt) *
           pow((*Critical_second_invariant_pt), 2.0) *
           pow((*Critical_second_invariant_pt), (*Flow_index_pt) / 2.0) -
         pow(2.0, (*Flow_index_pt) + 2.0) * (*Alpha_pt) *
           pow((*Critical_second_invariant_pt), (*Flow_index_pt) / 2.0 + 2.0));

      b =
        1.0 /
        (4.0 * pow((*Critical_second_invariant_pt), 2.0) *
         pow((*Critical_second_invariant_pt), 5.0 / 2.0)) *
        (-12.0 * (Cut_off_viscosity + epsilon - 1.0) *
           pow((*Critical_second_invariant_pt), 5.0 / 2.0) -
         pow(2.0, (*Flow_index_pt)) * ((*Flow_index_pt) - 1.0) * (*Alpha_pt) *
           pow((*Critical_second_invariant_pt), 2.0) *
           pow((*Critical_second_invariant_pt), (*Flow_index_pt) / 2.0) +
         3.0 * pow(2.0, (*Flow_index_pt) + 1.0) * (*Alpha_pt) *
           pow((*Critical_second_invariant_pt), (*Flow_index_pt) / 2.0 + 2.0));
    }


    /// Viscosity ratio as a fct of strain rate invariant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Get the parameters of the cubic
      double a;
      double b;

      calculate_fitting_parameters_of_cubic(a, b);

      double zero_shear_viscosity = calculate_zero_shear_viscosity();

      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // if the second invariant is below the cutoff we have a constant,
      // Newtonian viscosity
      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return a * pow(sign * second_invariant_of_rate_of_strain_tensor, 3.0) +
               b * pow(sign * second_invariant_of_rate_of_strain_tensor, 2.0) +
               zero_shear_viscosity;
      }
      else
      {
        // Calculate the square root of the absolute value of the
        // second invariant of the rate of strain tensor
        double measure_of_rate_of_strain =
          sqrt(sign * second_invariant_of_rate_of_strain_tensor);

        return 1.0 + (*Alpha_pt) * pow((2.0 * measure_of_rate_of_strain),
                                       *Flow_index_pt - 1.0);
      }
    }

    /// Deriv of viscosity w.r.t. strain rate invariant
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Get the parameters of the cubic
      double a;
      double b;

      calculate_fitting_parameters_of_cubic(a, b);

      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return sign * 3.0 * a *
                 pow(sign * second_invariant_of_rate_of_strain_tensor, 2.0) +
               2.0 * b * second_invariant_of_rate_of_strain_tensor;
      }
      else
      {
        return (*Alpha_pt) * pow(2.0, (*Flow_index_pt) - 2.0) *
               ((*Flow_index_pt) - 1.0) *
               second_invariant_of_rate_of_strain_tensor *
               pow(sign * second_invariant_of_rate_of_strain_tensor,
                   ((*Flow_index_pt) - 1.0) / 2.0 - 2.0);
      }
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a Casson model fluid
  /// using Tanner and Milthorpe's (1983) regularisation
  /// with a smooth transition using a cubic
  //==================================================================
  template<unsigned DIM>
  class CassonTanMilRegWithBlendingConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// Yield stress
    double* Yield_stress_pt;

    /// value of the second invariant below which we have
    /// constant (Newtonian) viscosity -- assumed to be always positive
    double* Critical_second_invariant_pt;


  public:
    /// "Cutoff regularised" Casson constitutive equation
    CassonTanMilRegWithBlendingConstitutiveEquation(
      double* yield_stress_pt, double* critical_second_invariant_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Yield_stress_pt(yield_stress_pt),
        Critical_second_invariant_pt(critical_second_invariant_pt)
    {
      /// get the cutoff viscosity
      double cut_off_viscosity = calculate_cutoff_viscosity();

      /// get the zero shear viscosity
      double zero_shear_viscosity = calculate_zero_shear_viscosity();

      oomph_info << "CassonTanMilRegWithBlendingConstitutiveEquation: "
                 << " zero shear viscosity = " << zero_shear_viscosity
                 << std::endl;

      oomph_info << "CassonTanMilRegWithBlendingConstitutiveEquation: "
                 << " cut off viscosity = " << cut_off_viscosity << std::endl;

      oomph_info << "                                              "
                 << " cutoff invariant = " << *Critical_second_invariant_pt
                 << std::endl;
    }

    /// Function that calculates the viscosity at the cutoff invariant
    /// Note: this is NOT the viscosity at zero I2
    double calculate_cutoff_viscosity()
    {
      // Calculate the Newtonian cutoff viscosity from the constitutive
      // equation and the cutoff value of the second invariant
      return (*Yield_stress_pt) / (2.0 * sqrt(*Critical_second_invariant_pt)) +
             2.0 * sqrt(*Yield_stress_pt /
                        (2.0 * sqrt(*Critical_second_invariant_pt))) +
             1.0;
    }

    /// Offset by how much the zero shear rate viscosity lies
    /// above the viscosity at I2_cutoff
    /// Hard-coded to a value that ensures a smooth transition
    double calculate_viscosity_offset_at_zero_shear(double& cut_off_viscosity)
    {
      return cut_off_viscosity / 5.0;
    }

    /// Function that calculates the viscosity at zero I2
    double calculate_zero_shear_viscosity()
    {
      // get the viscosity at the cutoff point
      double cut_off_viscosity = calculate_cutoff_viscosity();

      /// Offset by how much the zero shear rate viscosity lies
      /// above the viscosity at I2_cutoff
      double epsilon =
        calculate_viscosity_offset_at_zero_shear(cut_off_viscosity);

      return cut_off_viscosity + epsilon;
    }

    /// Report cutoff values
    void report_cut_off_values(double& cut_off_invariant,
                               double& cut_off_viscosity,
                               double& zero_shear_viscosity)
    {
      cut_off_invariant = *Critical_second_invariant_pt;
      cut_off_viscosity = calculate_cutoff_viscosity();
      zero_shear_viscosity = calculate_zero_shear_viscosity();
    }

    // Calculate the fitting parameters for the cubic blending
    void calculate_fitting_parameters_of_cubic(double& a, double& b)
    {
      // get the viscosity at the cutoff invariant
      double Cut_off_viscosity = calculate_cutoff_viscosity();

      // calculate the offset at zero shear
      double epsilon =
        calculate_viscosity_offset_at_zero_shear(Cut_off_viscosity);

      a = 1.0 / pow(*Critical_second_invariant_pt, 39.0 / 4.0) *
          (pow(*Critical_second_invariant_pt, 27.0 / 4.0) *
             (2.0 * (Cut_off_viscosity + epsilon) -
              2.0 * sqrt(2.0) *
                sqrt(*Yield_stress_pt / sqrt(*Critical_second_invariant_pt)) -
              2.0) -
           5.0 / 4.0 * (*Yield_stress_pt) *
             pow(*Critical_second_invariant_pt, 25.0 / 4.0) -
           1.0 / (2.0 * sqrt(2.0)) * sqrt(*Yield_stress_pt) *
             pow(*Critical_second_invariant_pt, 13.0 / 2.0));

      b = 1.0 / pow(*Critical_second_invariant_pt, 27.0 / 4.0) *
          (pow(*Critical_second_invariant_pt, 19.0 / 4.0) *
             (-3.0 * (Cut_off_viscosity + epsilon) +
              3.0 * sqrt(2.0) *
                sqrt(*Yield_stress_pt / sqrt(*Critical_second_invariant_pt)) +
              3.0) +
           7.0 / 4.0 * (*Yield_stress_pt) *
             pow(*Critical_second_invariant_pt, 17.0 / 4.0) +
           1.0 / (2.0 * sqrt(2.0)) * sqrt(*Yield_stress_pt) *
             pow(*Critical_second_invariant_pt, 9.0 / 2.0));
    }


    /// Viscosity ratio as a fct of strain rate invariant
    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Get the parameters of the cubic
      double a;
      double b;

      calculate_fitting_parameters_of_cubic(a, b);

      double zero_shear_viscosity = calculate_zero_shear_viscosity();

      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      // if the second invariant is below the cutoff we have a constant,
      // Newtonian viscosity
      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return a * pow(sign * second_invariant_of_rate_of_strain_tensor, 3.0) +
               b * pow(sign * second_invariant_of_rate_of_strain_tensor, 2.0) +
               zero_shear_viscosity;
      }
      else
      {
        // Calculate the square root of the absolute value of the
        // second invariant of the rate of strain tensor
        double measure_of_rate_of_strain =
          sqrt(sign * second_invariant_of_rate_of_strain_tensor);

        return *Yield_stress_pt / (2.0 * measure_of_rate_of_strain) +
               2.0 *
                 sqrt(*Yield_stress_pt / (2.0 * measure_of_rate_of_strain)) +
               1.0;
      }
    }

    /// Deriv of viscosity w.r.t. strain rate invariant
    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Get the parameters of the cubic
      double a;
      double b;

      calculate_fitting_parameters_of_cubic(a, b);

      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      if (sign * second_invariant_of_rate_of_strain_tensor <
          (*Critical_second_invariant_pt))
      {
        return sign * 3.0 * a *
                 pow(sign * second_invariant_of_rate_of_strain_tensor, 2.0) +
               2.0 * b * second_invariant_of_rate_of_strain_tensor;
      }
      else
      {
        return -sqrt(*Yield_stress_pt) *
                 second_invariant_of_rate_of_strain_tensor /
                 (2.0 * sqrt(2.0) *
                  pow(sign * second_invariant_of_rate_of_strain_tensor,
                      9.0 / 4.0)) -
               (*Yield_stress_pt) * second_invariant_of_rate_of_strain_tensor /
                 (4.0 * pow(sign * second_invariant_of_rate_of_strain_tensor,
                            5.0 / 2.0));
      }
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining an arbitrary shear-thinning fluid
  //==================================================================
  template<unsigned DIM>
  class NicosConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// high shear rate viscosity
    double* Mu_inf_pt;

    /// zero shear rate viscosity
    double* Mu_0_pt;

    /// parameter that controls the steepness of the curve
    double* Alpha_pt;


  public:
    NicosConstitutiveEquation(double* mu_inf_pt,
                              double* mu_0_pt,
                              double* alpha_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Mu_inf_pt(mu_inf_pt),
        Mu_0_pt(mu_0_pt),
        Alpha_pt(alpha_pt)
    {
    }

    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      return (*Mu_inf_pt) +
             ((*Mu_0_pt) - (*Mu_inf_pt)) *
               exp(-(*Alpha_pt) * second_invariant_of_rate_of_strain_tensor);
    }

    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // std::ostringstream error_stream;
      // error_stream << "This has not been implemented yet!";
      // throw OomphLibError(
      // error_stream.str(),
      // OOMPH_CURRENT_FUNCTION,
      // OOMPH_EXCEPTION_LOCATION);

      // return 0.0;

      return (*Alpha_pt) * ((*Mu_inf_pt) - (*Mu_0_pt)) *
             exp(-(*Alpha_pt) * second_invariant_of_rate_of_strain_tensor);
    }
  };

  //===================================================================
  /// A GeneralisedNewtonianConstitutiveEquation class
  /// defining a fluid following a tanh-profile
  //==================================================================
  template<unsigned DIM>
  class TanhProfileConstitutiveEquation
    : public GeneralisedNewtonianConstitutiveEquation<DIM>
  {
  private:
    /// high shear rate viscosity
    double* Mu_inf_pt;

    /// zero shear rate viscosity
    double* Mu_0_pt;

    /// parameter controlling the steepness of the step
    /// nb -- I used 10.0/(*Critical_second_invariant_pt)
    double* Alpha_pt;

    /// parameter that controls the location of the step -- assumed
    /// to be always positive
    double* Critical_second_invariant_pt;

  public:
    TanhProfileConstitutiveEquation(double* mu_inf_pt,
                                    double* mu_0_pt,
                                    double* alpha_pt,
                                    double* critical_second_invariant_pt)
      : GeneralisedNewtonianConstitutiveEquation<DIM>(),
        Mu_inf_pt(mu_inf_pt),
        Mu_0_pt(mu_0_pt),
        Alpha_pt(alpha_pt),
        Critical_second_invariant_pt(critical_second_invariant_pt)
    {
    }

    double viscosity(const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      return ((*Mu_0_pt) - (*Mu_inf_pt)) / 2.0 *
               tanh(((*Critical_second_invariant_pt) -
                     sign * second_invariant_of_rate_of_strain_tensor) *
                    (*Alpha_pt)) -
             ((*Mu_0_pt) - (*Mu_inf_pt)) / 2.0 + (*Mu_0_pt);
    }

    double dviscosity_dinvariant(
      const double& second_invariant_of_rate_of_strain_tensor)
    {
      // Pre-multiply the second invariant with +/-1 depending on whether it's
      // positive or not
      double sign = -1.0;
      if (second_invariant_of_rate_of_strain_tensor >= 0.0)
      {
        sign = 1.0;
      }

      return -sign * ((*Mu_0_pt) - (*Mu_inf_pt)) * 10.0 /
             (2.0 * (*Critical_second_invariant_pt)) * 1.0 /
             cosh(((*Critical_second_invariant_pt) -
                   sign * second_invariant_of_rate_of_strain_tensor) *
                  (*Alpha_pt)) *
             1.0 /
             cosh(((*Critical_second_invariant_pt) -
                   sign * second_invariant_of_rate_of_strain_tensor) *
                  (*Alpha_pt));
    }
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////

} // namespace oomph


#endif

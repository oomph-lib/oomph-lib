#include "poly_interp.h"

namespace oomph
{
 using namespace PolyInterpHelpers;

 /// Get value at point x.
 void BarycentricLagrangeInterpolator::
 eval(const double& x, Vector<double>& result) const
 {
  eval_checks(x);

  // Look for this x in the list of locations
  int matching_location = PolyInterpHelpers::fp_find(x, Locations);

  // If it's not at a point in the list of known points:
  if(matching_location == -1)
   {
    // Calculate the polynomial according to (4.2) of Berrut2004.
    Vector<double> numerator(Values[0].size(), 0.0);
    double denominator = 0;
    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
     {
      double temp = Weights[i]  / (x - Locations[i]);
      vector_add_ax(temp, Values[i], numerator);
      denominator += temp;
     }

    result.reserve(numerator.size());
    for(unsigned i=0, ni=numerator.size(); i<ni; i++)
     {
      result.push_back(numerator[i] / denominator);
     }
   }
  // Otherwise we already know the answer
  else
   {
    result = Values[matching_location];
   }
 }

 /// Get derivative of the value at point x.
 void BarycentricLagrangeInterpolator::
 eval_derivative(const double& x, const unsigned &deriv_order,
                 Vector<double>& result) const
 {
  eval_checks(x);

  // Only implemented for first derivatives...
#ifdef PARANOID
  if(deriv_order != 1)
   {
    std::string error_msg = "Only implemented for first derivatives.";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // Look for this x in the list of locations
  int matching_location = fp_find(x, Locations);

  // If it's not at a point in the list of known points:
  if(matching_location == -1)
   {
    Vector<double> g(Values[0].size(), 0.0), dg(Values[0].size(), 0.0);
    double h = 0, dh = 0;

    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
     {
      double temp = Weights[i]  / (x - Locations[i]);

      vector_add_ax(temp, Values[i], g);
      vector_add_ax(-1* temp / (x - Locations[i]), Values[i], dg);
      h += temp;
      dh += -1 * temp / (x - Locations[i]);
     }

#ifdef PARANOID
    if(h == 0)
     {
      std::string error_msg = "About to divide by zero...";
      throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    result.assign(Values[0].size(), 0.0);
    for(unsigned i=0, ni=g.size(); i<ni; i++)
     {
      result[i] = (dg[i] * h - g[i]*dh) / (h*h);
     }
   }

  // Otherwise we need to do something else to avoid numerically dodgy
  // behaviour. See Berrut2004 pg 513.
  else
   {
    // The index of where we want to calculation the derivative (rename
    // to be consistent with the paper).
    unsigned i = matching_location;

    // Compute special weights
    Vector<double> deriv_weights(Locations.size(), 0.0);
    for(unsigned j=0, nj = Locations.size(); j < nj; j++)
     {
      if(i != j)
       {
        deriv_weights[j] = (Weights[j] / Weights[i]) / (Locations[i] - Locations[j]);

        // Maths trick, see the paper
        deriv_weights[i] -= deriv_weights[j];
       }
     }


    // For each value that we are interpolating
    result.assign(Values[0].size(), 0.0);
    for(unsigned k=0, nk=Values[0].size(); k<nk; k++)
     {
      // Interpolate it
      for(unsigned j=0, nj=Values.size(); j<nj; j++)
       {
        result[k] += deriv_weights[j] * Values[j][k];
       }
     }
   }
 }

 /// \short Check the everything has been set up right.
 void BarycentricLagrangeInterpolator::
 eval_checks(const double &x) const
 {
#ifdef PARANOID
  if(Locations.size() == 0)
   {
    std::string error_msg = "We need some values/locations to interpolate from!";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

#ifdef PARANOID
  if(Locations.size() != Values.size())
   {
    std::string error_msg = "Locations and values must be the same size.";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

#ifdef PARANOID
  if(Weights.size() != Locations.size())
   {
    std::string error_msg = "Weights and locations must be the same size.";
    throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

 }

}

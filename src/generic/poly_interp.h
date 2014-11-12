#ifndef OOMPH_POLY_INTERP_H
#define OOMPH_POLY_INTERP_H

#include "Vector.h"
#include <set>

namespace oomph
{

namespace PolyInterpHelpers
{
  /// Check if a vector contains any duplicate values
  template <typename T>
  inline bool contains_duplicates(const std::vector<T> &v)
  {
    // Construct a set (which has no duplicates by definition) and compare
    // sizes.
    return std::set<T>(v.begin(), v.end()).size() != v.size();
  }

 /// \short Equivalent to std::find but for floating point values. Return
 /// -1 if not found.
 inline int fp_find(double search_value, const Vector<double> &vec, double tol=1e-12)
 {
  int found_location = -1;
  for(unsigned j=0, nj=vec.size(); j<nj; j++)
   {
    if(std::abs(vec[j] - search_value) < tol)
     {
      found_location = j;
      break;
     }
   }

  return found_location;
 }


 /// \short y += ax operation for vectors x and y.
 inline void vector_add_ax(double a, const Vector<double>& x,
                            Vector<double>& y)
  {
    // my_assert(y.size() == x.size());
   #ifdef PARANOID
   if(y.size() != x.size())
    {
     std::string error_msg = "y and x must have the same size.";
     throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

    for (unsigned i = 0, ni=y.size(); i<ni; i++)
      {
        y[i] += a * x[i];
      }
  }


 /// \short Check if two doubles are with a tol of each other.
 inline bool almost_equal(double a, double b, double tol=1e-10)
 {
  return std::abs(a - b) < tol;
 }


 /// \short Check if vectors x and y are equal to within a tolerance.
 inline bool almost_equal(const Vector<double>& x,
                           const Vector<double> &y,
                           const double &tol=1e-10)
  {
    if(x.size() != y.size())
      {
        return false;
      }
    else
      {
        for(unsigned i=0, ni=x.size(); i<ni; i++)
          {
            if(!almost_equal(x[i], y[i], tol))
              {
                return false;
              }
          }
      }
    return true;
  }
}



  // ============================================================
  /// Interpolate a (list of) polynomial function(s) through the given
  /// values.
  // ============================================================
  class BarycentricLagrangeInterpolator
  {
  public:

    /// Construct from locations, values lists.
    BarycentricLagrangeInterpolator(const Vector<double>& locations,
                                    const Vector<Vector<double> >& values)
      : Locations(locations), Values(values)
    { build(); }

    /// Empty virtual destructor
    virtual ~BarycentricLagrangeInterpolator() {}

    /// Get value at point x.
    void eval(const double& x, Vector<double>& result) const;

    /// Get nth derivative at point x.
    void eval_derivative(const double& x, const unsigned &deriv_order,
                         Vector<double>& result) const;

  private:

    /// Construct the weights.
    void build()
    {
#ifdef PARANOID
     if(PolyInterpHelpers::contains_duplicates(Locations))
      {
       std::string error_msg = "";
       throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif

#ifdef PARANOID
     if(Locations.size() != Values.size())
      {
       std::string error_msg = "";
       throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // Make the right sized weights vector
      Weights.clear();
      Weights.resize(Locations.size());

      // Calculate each weight according to (3.2) in Berrut2004 (O(N^2) where
      // N is the number of interpolation points).
      for(unsigned j=0, nj=Locations.size(); j<nj; j++)
        {
          double product = 1.0;
          for(unsigned k=0, nk = Locations.size(); k<nk; k++)
            {
              if(k != j) product *= (Locations[j] - Locations[k]);
            }

          Weights[j] = 1 / product;
        }

    }

   /// Check that everything has been set up properly.
    void eval_checks(const double &x) const;

    /// \short Location of the known points that we are using to
    /// interpolate.
    Vector<double> Locations;

   /// \short Values of the known points that we are using to interpolate. First
   /// axis is "places to interpolate" second axis is a list of values to
   /// interpolate at that place. (it's more efficient to interpolate a
   /// list of functions at once).
   Vector<Vector<double> > Values;

   /// \short Data storage for the weights used in the interpolation (calculated
   /// by build()).
   Vector<double> Weights;

  };


} // End of oomph namespace

#endif

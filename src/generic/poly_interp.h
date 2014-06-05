#ifndef OOMPH_POLY_INTERP_H
#define OOMPH_POLY_INTERP_H


#include "Vector.h"


namespace oomph
{

namespace PolyInterpHelpers
{
 // Equivalent to std::find but for floating point values. Return -1 if
 // not found.
 int fp_find(double search_value, const Vector<double> &vec, double tol=1e-12);
 
 /// Check if two doubles are with a tol of each other
 bool almost_equal(double a, double b, double tol=1e-10);
 
 void vector_add_ax(double a, const Vector<double>& x,
                    Vector<double>& y);

 bool almost_equal(const Vector<double>& x,
                   const Vector<double> &y,
                   const double &tol);
}

  // ============================================================
  ///
  // ============================================================
  class PolynomialInterpolatorBase
  {
  public:
    /// Get value at point
    virtual void eval(const double& x, Vector<double>& result) const =0;

    /// Get nth derivative at point
    virtual void eval_derivative(const double& x, const unsigned &deriv_order,
                                 Vector<double>& result) const =0;

    /// Empty virtual destructor
    virtual ~PolynomialInterpolatorBase() {}

  };


  // ============================================================
  ///
  // ============================================================
  class BarycentricLagrangeInterpolator : public PolynomialInterpolatorBase
  {
  public:

    /// Construct from locations, values lists.
    BarycentricLagrangeInterpolator(const Vector<double>& locations,
                                    const Vector<Vector<double> >& values)
      : Locations(locations), Values(values)
    { build(); }

    /// Empty virtual destructor
    virtual ~BarycentricLagrangeInterpolator() {}

    /// Get value at point
    void eval(const double& x, Vector<double>& result) const;

    /// Get nth derivative at point
    void eval_derivative(const double& x, const unsigned &deriv_order,
                         Vector<double>& result) const;

  private:

    /// Construct the weights
    void build()
    {
      // assert(!VectorOps::contains_duplicates(Locations));
      // assert(Locations.size() == Values.size());

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

    void eval_checks(const double &x) const;

    /// Data storage
    Vector<double> Locations;

    /// Values for interpolation. First axis is "places to interpolate"
    /// second axis is "values to interpolate at that place".
    Vector<Vector<double> > Values;
    Vector<double> Weights;

  };


} // End of oomph namespace

#endif

#ifndef COMPLEX_LESS_HEADER
#define COMPLEX_LESS_HEADER

#include "complex"

namespace oomph
{
  //===================================================================
  /// Function-type-object to perform comparison of complex data types
  /// Needed to sort the complex eigenvalues into order based on the
  /// size of the real part.
  //==================================================================
  template<class T>
  class ComplexLess
  {
  public:
    /// Comparison in terms of magnitude of complex number
    bool operator()(const std::complex<T>& x, const std::complex<T>& y) const
    {
      //// Return the order in terms of magnitude if they are not equal
      //// Include a tolerance to avoid processor specific ordering
      // if (std::abs(std::abs(x) - std::abs(y)) > 1.0e-10)
      //{
      //  return std::abs(x) < std::abs(y);
      //}
      // Otherwise sort on real part, again with a tolerance
      // else
      if (std::abs(x.real() - y.real()) > 1.0e-10)
      {
        return x.real() < y.real();
      }
      // Otherwise sort on imaginary part
      else
      {
        return x.imag() < y.imag();
      }
    }
  };
} // namespace oomph

#endif

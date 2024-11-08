#ifndef STORE_ERROR_ELEMENT_HEADER
#define STORE_ERROR_ELEMENT_HEADER

namespace oomph
{
  /// An element to store the error generated from the error estimator, e.g.
  /// Z2ErrorEstimator
  template<class ELEMENT>
  class StoreErrorElement : public ELEMENT
  {
  private:
    /// Store the error
    double Error;

  public:
    /// Constructor
    StoreErrorElement() : ELEMENT(), Error(0.0) {}

    /// Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    /// Get error value for post-processing
    double get_error()
    {
      return Error;
    }
  };

} // namespace oomph

#endif

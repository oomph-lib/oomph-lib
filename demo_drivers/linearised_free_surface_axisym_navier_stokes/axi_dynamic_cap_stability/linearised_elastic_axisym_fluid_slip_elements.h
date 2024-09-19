#ifndef LINEARISED_ELASTIC_AXISYM_FLUID_SLIP_ELEMENTS_HEADER
#define LINEARISED_ELASTIC_AXISYM_FLUID_SLIP_ELEMENTS_HEADER

#include "linearised_axisym_fluid_slip_elements.h"

namespace oomph
{
  //======================================================================
  /// A class for elements that allow the imposition of an applied slip
  /// in the axisym Navier Stokes eqns.
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELEMENT>
  class LinearisedElasticAxisymmetricNavierStokesSlipElement
    : public LinearisedAxisymmetricNavierStokesSlipElement<ELEMENT>
  {
  public:
    LinearisedElasticAxisymmetricNavierStokesSlipElement(
      FiniteElement* const& element_pt, const int& face_index)
      : LinearisedAxisymmetricNavierStokesSlipElement<ELEMENT>(element_pt,
                                                               face_index)
    {
    }

  private:
    // Abstract function to get the solid perturbations
    virtual double Xhat(const unsigned& n, const unsigned& i, const unsigned& j)
    {
      return this->nodal_value(n, i * 2 + j);
    }
  };

} // namespace oomph
#endif

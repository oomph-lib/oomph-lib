#ifndef OVERLAYING_LINEARISED_ELASTIC_AXISYM_FLUID_SLIP_ELEMENT_HEADER
#define OVERLAYING_LINEARISED_ELASTIC_AXISYM_FLUID_SLIP_ELEMENT_HEADER

//#include "generic.h"
#include "linearised_elastic_axisym_fluid_slip_elements.h"
#include "axisym_fluid_slip_elements.h"
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"

namespace oomph
{
  template<class BASE_ELEMENT, class BULK_ELEMENT>
  class OverlayingLinearisedElasticAxisymmetricFluidSlipElement
    : public virtual LinearisedElasticAxisymmetricNavierStokesSlipElement<
        BULK_ELEMENT>
  {
  private:
    typedef AxisymmetricNavierStokesSlipElement<BASE_ELEMENT> BASE_FACE_ELEMENT;

    BASE_FACE_ELEMENT* Base_element_pt;

  public:
    OverlayingLinearisedElasticAxisymmetricFluidSlipElement(
      FiniteElement* const& element_pt, const int& face_index)
      : LinearisedElasticAxisymmetricNavierStokesSlipElement<BULK_ELEMENT>(
          element_pt, face_index)
    {
    }

    void set_base_element_pt(GeneralisedElement* const& el_pt)
    {
      Base_element_pt = dynamic_cast<BASE_FACE_ELEMENT*>(el_pt);
      if (!Base_element_pt)
      {
        throw(OomphLibError("Base element is not the correct type.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION));
      }
    }

    virtual Vector<double> get_base_u(const Vector<double>& s)
    {
      return Base_element_pt->interpolated_u(s);
    }

    virtual Vector<double> get_base_wall_velocity(const Vector<double>& s)
    {
      return Base_element_pt->interpolated_wall_velocity(s);
    }
  };
} // namespace oomph
#endif

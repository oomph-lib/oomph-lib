#ifndef OVERLAYING_LINEARISED_ELASTIC_AXISYM_FLUID_INTERFACE_ELEMENT_HEADER
#define OVERLAYING_LINEARISED_ELASTIC_AXISYM_FLUID_INTERFACE_ELEMENT_HEADER

//#include "generic.h"
#include "linearised_elastic_axisym_fluid_interface_element.h"
#include "free_surface_elements.h"

namespace oomph
{
  template<class BASE_ELEMENT, class BULK_ELEMENT>
  class OverlayingLinearisedElasticAxisymmetricFluidInterfaceElement
    : public virtual Hijacked<
        LinearisedElasticAxisymmetricFluidInterfaceElement<BULK_ELEMENT>>
  {
  private:
    typedef FreeSurfaceElement<BASE_ELEMENT> BASE_FACE_ELEMENT;

    BASE_FACE_ELEMENT* Base_element_pt;

  public:
    OverlayingLinearisedElasticAxisymmetricFluidInterfaceElement(
      FiniteElement* const& element_pt,
      const int& face_index,
      TimeStepper* const& time_stepper_pt,
      const unsigned& id = 0)
      : LinearisedElasticAxisymmetricFluidInterfaceElement<BULK_ELEMENT>(
          element_pt, face_index, time_stepper_pt, id)
    {
    }

    virtual double pext()
    {
      // return 0;
      return Base_element_pt->pext();
    }

    void set_base_element_pt(GeneralisedElement* const& el_pt)
    {
      Base_element_pt = dynamic_cast<BASE_FACE_ELEMENT*>(el_pt);
    }

    /// Overload get_base_flow_u(...) to return the external
    /// element's velocity components at the integration point
    virtual void get_base_lagrange_multiplier(const double& time,
                                              const unsigned& ipt,
                                              const Vector<double>& x,
                                              double& result)
    {
      // result = Base_element_pt->lagrange(ipt);
      result = 0;
    } // End of overloaded get_base_flow_u function

    virtual double get_base_external_pressure()
    {
      return Base_element_pt->pext();
    }

    /// Overload get_base_flow_u(...) to return the external
    /// element's velocity components at the integration point
    virtual void get_base_pressure(const double& time,
                                   const unsigned& ipt,
                                   const Vector<double>& x,
                                   double& result)
    {
      const unsigned p_index = 3;

      Shape psif(this->nnode());
      double interpolated_p = 0;
      for (unsigned n = 0; n < 2; n++)
      {
        this->shape_at_knot(n, psif);
        const double psif_ = psif(n);
        interpolated_p = Base_element_pt->node_pt(ipt)->value(p_index) * psif_;
      }

      result = interpolated_p;
    } // End of overloaded get_base_flow_u function

    /// Overload get_base_flow_u(...) to return the external
    /// element's velocity components at the integration point
    virtual double get_base_u(const Vector<double>& s, const unsigned& i)
    {
      return Base_element_pt->interpolated_u(s, i);
    } // End of overloaded get_base_flow_u function

    /// Compute the element's residual vector and the Jacobian matrix
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Get the analytical contribution from the basic patricklinearised
      // element
      // LinearisedElasticAxisymmetricFluidInterfaceElement::
      //   fill_in_contribution_to_jacobian(residuals, jacobian);

      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }
  };

  template<class BASE_ELEMENT, class BULK_ELEMENT>
  class FaceGeometry<
    OverlayingLinearisedElasticAxisymmetricFluidInterfaceElement<BASE_ELEMENT,
                                                                 BULK_ELEMENT>>
    : public PointElement
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : PointElement() {}
  };


} // namespace oomph
#endif

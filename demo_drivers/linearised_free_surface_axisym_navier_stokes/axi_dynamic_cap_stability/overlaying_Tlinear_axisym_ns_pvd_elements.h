#ifndef OVERLAYING_TLINEAR_AXISYM_NS_PVD_ELEMENTS_HEADER
#define OVERLAYING_TLINEAR_AXISYM_NS_PVD_ELEMENTS_HEADER

//#include "generic.h"
#include "projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "Tlinear_axisym_ns_pvd_elements.h"

namespace oomph
{
  class OverlayingTLinearisedAxisymNSPVDElement
    : public virtual Hijacked<TLinearisedAxisymNSPVDElement>
  {
  private:
    typedef ProjectableAxisymmetricTTaylorHoodPVDElement BASE_ELEMENT;

    BASE_ELEMENT* Base_element_pt;

  public:
    OverlayingTLinearisedAxisymNSPVDElement() : TLinearisedAxisymNSPVDElement()
    {
    }

    void set_base_element_pt(GeneralisedElement* const& el_pt)
    {
      Base_element_pt = dynamic_cast<BASE_ELEMENT*>(el_pt);
    }

    /// Overload get_base_flow_u(...) to return the external
    /// element's velocity components at the integration point
    virtual void get_base_flow_u(const double& time,
                                 const unsigned& ipt,
                                 const Vector<double>& x,
                                 Vector<double>& result) const
    {
      Vector<double> s(2);

      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      Base_element_pt->interpolated_u_axi_nst(s, result);
    } // End of overloaded get_base_flow_u function

    virtual void get_base_flow_u(const Vector<double>& s,
                                 Vector<double>& result) const
    {
      Base_element_pt->interpolated_u_axi_nst(s, result);
    } // End of overloaded get_base_flow_u function

    /// Overload get_base_flow_dudx(...) to return the derivatives of
    /// the external element's velocity components w.r.t. global coordinates
    /// at the integration point
    virtual void get_base_flow_dudx(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& result) const
    {
      Vector<double> s(2);

      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          result(i, j) = Base_element_pt->interpolated_dudx_axi_nst(s, i, j);
        }
      }

    } // End of overloaded get_base_flow_dudx function

    virtual void get_base_flow_dudx(const Vector<double>& s,
                                    DenseMatrix<double>& result) const
    {
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          result(i, j) = Base_element_pt->interpolated_dudx_axi_nst(s, i, j);
        }
      }
    } // End of overloaded get_base_flow_u function

    /// Overload get_base_flow_p(...) to return the external
    /// element's pressure at the integration point
    virtual void get_base_flow_p(const double& time,
                                 const unsigned& ipt,
                                 const Vector<double>& x,
                                 double& result) const
    {
      Vector<double> s(2);

      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      result = Base_element_pt->interpolated_p_axi_nst(s);
    } // End of overloaded get_base_flow_p function

    virtual void get_base_flow_p(const Vector<double>& s, double& result) const
    {
      result = Base_element_pt->interpolated_p_axi_nst(s);
    }

    /// Overload get_base_flow_dudt(...) to return the derivatives of
    /// the external element's velocity components w.r.t. time
    /// at the integration point
    virtual void get_base_flow_dudt(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    Vector<double>& result) const
    {
      Vector<double> s(2);

      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      for (unsigned i = 0; i < 2; i++)
      {
        result[i] = Base_element_pt->interpolated_dudt_axi_nst(s, i);
      }
    } // End of overloaded get_base_flow_dudt function

    virtual void get_base_flow_dudt(const Vector<double>& s,
                                    Vector<double>& result) const
    {
      for (unsigned i = 0; i < 2; i++)
      {
        result[i] = Base_element_pt->interpolated_dudt_axi_nst(s, i);
      }
    }

    /// Overload get_base_flow_duds(...) to return the derivatives of
    /// the external element's velocity components w.r.t. local coordinates
    /// at the integration point
    virtual void get_base_flow_duds(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& result) const
    {
      Vector<double> s(2);

      // Assign values of the local coordinates s
      for (unsigned i = 0; i < 2; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          result(i, j) = Base_element_pt->interpolated_duds_axi_nst(s, i, j);
        }
      }

    } // End of overloaded get_base_flow_duds function

    virtual void get_base_flow_duds(const Vector<double>& s,
                                    DenseMatrix<double>& result) const
    {
      for (unsigned i = 0; i < 2; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          result(i, j) = Base_element_pt->interpolated_duds_axi_nst(s, i, j);
        }
      }
    }
  };

  template<>
  class FaceGeometry<OverlayingTLinearisedAxisymNSPVDElement>
    : public TElement<1, 3>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : TElement<1, 3>() {}
  };

} // namespace oomph
#endif

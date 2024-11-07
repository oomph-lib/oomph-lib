#ifndef DEBUG_ELASTIC_AXISYMMETRIC_VOLUME_CONSTRAINT_BOUNDARY_ELEMENTS_HEADER
#define DEBUG_ELASTIC_AXISYMMETRIC_VOLUME_CONSTRAINT_BOUNDARY_ELEMENTS_HEADER

#include "solid/solid_elements.h"
#include "fluid_interface/constrained_volume_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  template<class ELEMENT>
  class DebugElasticAxisymmetricVolumeConstraintBoundingElement
    : public virtual ElasticAxisymmetricVolumeConstraintBoundingElement<
        ELEMENT>,
      public virtual DebugJacobianSolidFiniteElement,
      public virtual SolidFaceElement
  {
  public:
    DebugElasticAxisymmetricVolumeConstraintBoundingElement(
      FiniteElement* const& element_pt, const int& face_index)
      : ElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(element_pt,
                                                                    face_index),
        SolidFaceElement()
    {
      this->add_other_bulk_node_positions_as_external_data();
    }

    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    /// This final over-ride is required because both SolidFiniteElements
    /// and FaceElements overload zeta_nodal
    virtual double zeta_nodal(const unsigned& n,
                              const unsigned& k,
                              const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// Set pointer to MacroElement -- overloads generic version
    /// and uses the MacroElement
    /// also as the default for the "undeformed" configuration.
    /// This assignment must be overwritten with
    /// set_undeformed_macro_elem_pt(...) if the deformation of
    /// the solid body is driven by a deformation of the
    /// "current" Domain/MacroElement representation of it's boundary.
    /// Can be overloaded in derived classes to perform additional
    /// tasks
    virtual void set_macro_elem_pt(MacroElement* macro_elem_pt)
    {
      Macro_elem_pt = macro_elem_pt;
      Undeformed_macro_elem_pt = macro_elem_pt;
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Compute the element's residual Vector and the jacobian matrix
    /// Virtual function can be overloaded by hanging-node version
    void fill_in_contribution_to_djacobian_dparameter(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam)
    {
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }
  };
} // namespace oomph

#endif

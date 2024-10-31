#ifndef FREE_SURFACE_ELEMENT_HEADER
#define FREE_SURFACE_ELEMENT_HEADER

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"

namespace oomph
{
  // My free surface element, overloads the sigma function to remove the
  // surface tension term
  template<class ELEMENT>
  class MyFreeSurfaceElement
    : public ElasticAxisymmetricFluidInterfaceElement<ELEMENT>
  {
  public:
    MyFreeSurfaceElement(FiniteElement* const& element_pt,
                         const int& face_index,
                         const unsigned& id = 0)
      : ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
          element_pt, face_index, id)
    {
    }

    double sigma(const Vector<double>& s_local)
    {
      return 0.0;
    }
  };

  // My free surface element, overloads the sigma function to remove the
  // surface tension term
  template<class ELEMENT>
  class MyAxisymmmetricFreeSurfaceElement
    : public ElasticAxisymmetricFluidInterfaceElement<ELEMENT>
  {
  public:
    MyAxisymmmetricFreeSurfaceElement(FiniteElement* const& element_pt,
                         const int& face_index,
                         const unsigned& id = 0)
      : ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
          element_pt, face_index, id)
    {
    }
    
    double sigma(const Vector<double>& s_local)
    {
      return 0.0;
    }
  };
} // namespace oomph
#endif

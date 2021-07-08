//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Header file for elements that are used to apply contact angle
// BC for Young Laplace eqns

#ifndef OOMPH_YOUNGLAPLACE_FLUX_ELEMENTS_HEADER
#define OOMPH_YOUNGLAPLACE_FLUX_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


// oomph-lib ncludes
#include "../generic/Qelements.h"

namespace oomph
{

//======================================================================
/// \short A class for elements that allow the imposition of an 
/// contact angle bcs for Young Laplace elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class. Jacobian is evaluated by finite differencing.
//======================================================================
template <class ELEMENT>
class YoungLaplaceContactAngleElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
public:

 /// \short Constructor, takes the pointer to the "bulk" element and the 
 /// index of the face to which the element is attached.
 YoungLaplaceContactAngleElement(FiniteElement* const &bulk_el_pt, 
                         const int &face_index);

 ///\short  Broken empty constructor
 YoungLaplaceContactAngleElement()
  {
   throw OomphLibError(
    "Don't call empty constructor for YoungLaplaceContactAngleElement",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Broken copy constructor
 YoungLaplaceContactAngleElement(const YoungLaplaceContactAngleElement& dummy) 
  { 
   BrokenCopy::broken_copy("YoungLaplaceContactAngleElement");
  } 
 
 /// Broken assignment operator
 void operator=(const YoungLaplaceContactAngleElement&) 
  {
   BrokenCopy::broken_assign("YoungLaplaceContactAngleElement");
  }


 /// \short Access function for the pointer to the prescribed contact angle
 /// (const version)
 double* prescribed_cos_gamma_pt() const
  {
   return Prescribed_cos_gamma_pt;
  }


 /// Access function for the pointer to the prescribed contact angle
 double*& prescribed_cos_gamma_pt() 
  {
   return Prescribed_cos_gamma_pt;
  }


 /// Add the element's contribution to its residual vector
 void fill_in_contribution_to_residuals(Vector<double> &residuals);

 /// \short Get the local equation number of the (one and only) unknown
 /// stored at local node n (returns -1 if value is pinned).
 /// Can be overloaded in derived multiphysics elements.
 inline int u_local_eqn(const unsigned& n)
  {
   //Local equation number is the first value stored at the node
   return nodal_local_eqn(n,0);
  }

 
 /// Specify the value of nodal zeta from the face geometry
 /// \short The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                   const unsigned &i) const 
 {return FaceElement::zeta_nodal(n,k,i);}     

 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

 /// \short Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile, const unsigned &n_plot)
  {FiniteElement::output(outfile,n_plot);}

 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}

 /// \short C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}

 /// Compute cosinus of actual contact angle
 double actual_cos_contact_angle(const Vector<double>& s);

 /// Get unit tangent and normal vectors to contact line
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal)
  {
   // Dummy
   double norm_of_drds;
   Vector<double> spine(3);
   contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
  }

 /// \short Get unit tangent and normal vectors to contact line. Final
 /// argument gives the norm of dR/ds where R is the vector to the
 /// contact line and s the local coordinate in the element.
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal,
			   double& norm_of_drds)
  {
   // Dummy
   Vector<double> spine(3);
   contact_line_vectors(s,tangent,normal,spine,norm_of_drds);
  }


 /// \short Get tangent and normal vectors to contact line and spine vector 
 /// (wall normal can then be obtained by cross-product). Final
 /// argument gives the norm of dR/ds where R is the vector to the
 /// contact line and s the local coordinate in the element.
 void contact_line_vectors(const Vector<double>& s,
                           Vector<double>& tangent,
                           Vector<double>& normal,
                           Vector<double>& spine,
                           double& norm_of_drds);

protected:

 /// \short Define an access function to the first data value stored 
 /// at each node. In a more general "Equation" element, such abstraction
 /// is essential, because different Elements will store the same variables
 /// in different locations.
 double &u(const unsigned int &n) {return *this->node_pt(n)->value_pt(0);}

 /// Function to calculate the cos of the prescribed contact angle
 double prescribed_cos_gamma()
  {
   //If the function pointer is zero return zero
   if(Prescribed_cos_gamma_pt == 0)
    {
     return 0.0;
    }
   //Otherwise de-reference pointer
   else
    {
     return *Prescribed_cos_gamma_pt;
    }
  }


private:
 
 /// Pointer to prescribed cos gamma
 double* Prescribed_cos_gamma_pt;


}; 


}

#endif

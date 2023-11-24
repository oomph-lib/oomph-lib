//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#include<fenv.h>

#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The equations equations
#include "unsteady_heat.h"

#include "heat_transfer_and_melt_elements.h"
#include "temporary_stefan_boltzmann_elements.h"

//(Pseudo-)solid
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace oomph;
using namespace std;




/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////


//======================================================================
/// A class for elements that allow the imposition of (mutual)
/// Stefan Boltzmann heat flux on the boundaries of UnsteadyHeat elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT>
/// policy class.
//======================================================================
template <class ELEMENT>
class StefanBoltzmannMeltElement :
 public virtual StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>,
 public virtual SurfaceMeltElement<ELEMENT>
{
 
public:


 /// Constructor, takes the pointer to the "bulk" element and the
 /// index of the face to be created. Calls the constructors of the 
 /// underlying classes
 StefanBoltzmannMeltElement(FiniteElement* const &bulk_el_pt,
                            const int &face_index) :
  UnsteadyHeatBaseFaceElement<ELEMENT>(bulk_el_pt,face_index),
  StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>(bulk_el_pt,face_index),
  SurfaceMeltElement<ELEMENT>(bulk_el_pt,face_index)
  {
  }
 

 /// Broken copy constructor
 StefanBoltzmannMeltElement(
  const StefanBoltzmannMeltElement& dummy)
  {
   BrokenCopy::broken_copy("StefanBoltzmannMeltElement");
  }
 
  /// Specify the value of nodal zeta from the face geometry:
  /// The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default (final overload)
  double zeta_nodal(const unsigned &n, const unsigned &k,           
                    const unsigned &i) const 
  {return FaceElement::zeta_nodal(n,k,i);}     
  

};





/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//=============================================================================
/// One-dimensional integration scheme that locates integration points
/// uniformly along boundary (if elements are equal sized, that is...)
//=============================================================================
class MyIntegral : public Integral
{

public:

 /// Constructor: Specify number of uniformly spaced sample points in 
 /// unit interval
 MyIntegral(const unsigned& n_knot){N_knot=n_knot;}

 /// Broken copy constructor
 MyIntegral(const MyIntegral& dummy) 
  { 
   BrokenCopy::broken_copy("MyIntegral");
  } 
 
 /// Broken assignment operator
 void operator=(const MyIntegral&) 
  {
   BrokenCopy::broken_assign("MyIntegral");
  }

 /// Number of integration points of the scheme
 unsigned nweight() const {return N_knot;}

 /// Return coordinate s[j] (j=0) of integration point i -- 
 double knot(const unsigned &i, const unsigned &j) const 
  {
   double dx=1.0/double(N_knot);
   return ((0.5+double(i)))*dx;
  }
 
 /// Return weight of integration point i
 double weight(const unsigned &i) const 
  {
   //return 2.0/double(N_knot);
   return 1.0/double(N_knot);
  }
 
private:
 
 /// Number of integration points
 unsigned N_knot;
}; 



/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the problem parameters
//=====================================================================
namespace GlobalParameters
{
 /// Output directory
 string Directory="RESLT";

 /// Number of integration points for new integration scheme (if used)
 unsigned Nintpt=10;
 
 /// Non-dim density for pseudo-solid
 double Lambda_sq=0.0;

 /// Poisson's ratio for pseudo-solid
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Melt-temperature
 double Melt_temperature=0.8288627710; 
  
 /// Non-dim Stefan Boltzmann constant
 double Sigma= 1.0e-2;

 /// Zero degrees Celsius offset in Stefan Boltzmann law
 double Theta_0=1.0;

 /// Thermal inertia in inner region
 double Alpha0=1.0;

 /// Thermal inertia in outer annular region
 double Alpha1=1.0;

 /// Thermal conductivity in inner region
 double Beta0=0.05;

 /// Thermal conductivity in outer annular region
 double Beta1=1.5;

 /// Target element size
 double Target_area=0.05;

 /// Initial radius of inner circle
 double Radius_innermost=0.5;

 /// Temporal variation in inner radius (for exact solution)
 double R_hat=0.1; 

 /// Inner radius of annular region
 double Radius_inner=1.0;

 /// Outer radius of annular region
 double Radius_outer=1.5;

 /// Temperature on boundary of inner circle
 double U0=0.8288627710; 

 /// Coeff for (time-)constant variation of temperature in inner circle
 double V0=1.0;

 /// Coeff for time variation inner circle
 double V0_hat=0.5;

 /// Source function
 void get_source(const double& time, const Vector<double>& x, double& source)
 {

  double t=time;
  double t0=0.0;

  double MapleGenVar1 = 0.0;     
  double MapleGenVar4 = 0.0;
  double MapleGenVar6 = 0.0;
  double MapleGenVar9 = 0.0;
  double MapleGenVar8 = 0.0;
  double MapleGenVar7 = 0.0;
  double MapleGenVar5 = 0.0;
  double MapleGenVar3 = 0.0;
  double MapleGenVar11 = 0.0;
  double MapleGenVar13 = 0.0;
  double MapleGenVar14 = 0.0;
  double MapleGenVar12 = 0.0;
  double MapleGenVar10 = 0.0;
  double MapleGenVar2  = 0.0;

  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  if (r<0.5*(Radius_inner+Radius_innermost))
   {

//#--------------------------------------------------------------------
//# Source fct in inner region
//#--------------------------------------------------------------------
//> C(eval(S0));

      t0 = Beta0*(2.0*V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t)))-Alpha0*(
(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/0.3141592653589793E1
/2.0))*R_hat*(1.0-cos(2.0*0.3141592653589793E1*t))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)-(r*r-pow(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0),2.0))*V0_hat*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1/2.0);

   }
  else
   {

//#--------------------------------------------------------------------
//# Source fct in outer region
//#--------------------------------------------------------------------
//> C(eval(S1));

      MapleGenVar2 = Beta1/r;
      MapleGenVar8 = -Beta0*0.3141592653589793E1*0.3141592653589793E1*cos(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*t*Radius_inner/2.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*R_hat*V0_hat*t*t*Radius_outer/2.0+Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_inner-Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_outer+
Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*V0_hat*t
*Radius_inner/2.0-Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*
R_hat*R_hat*V0_hat*t*Radius_outer/2.0-Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0*Radius_inner/4.0+Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0*Radius_outer/4.0-Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0_hat*Radius_inner/8.0+Beta0*pow(sin(2.0*
0.3141592653589793E1*t),2.0)*R_hat*R_hat*V0_hat*Radius_outer/8.0-Beta0*
0.3141592653589793E1*0.3141592653589793E1*V0*Radius_inner*Radius_innermost*
Radius_innermost+Beta0*0.3141592653589793E1*0.3141592653589793E1*V0*
Radius_outer*Radius_innermost*Radius_innermost;
      MapleGenVar9 = MapleGenVar8-Beta0*0.3141592653589793E1*
0.3141592653589793E1*V0_hat*Radius_inner*Radius_innermost*Radius_innermost/2.0+
Beta0*0.3141592653589793E1*0.3141592653589793E1*V0_hat*Radius_outer*
Radius_innermost*Radius_innermost/2.0-0.3141592653589793E1*0.3141592653589793E1
*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*t*Radius_inner+
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*R_hat*t*Radius_outer+sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*cos(
2.0*0.3141592653589793E1*t)*R_hat*R_hat*Radius_inner/2.0-sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
R_hat*R_hat*Radius_outer/2.0;
      MapleGenVar7 = MapleGenVar9+0.3141592653589793E1*0.3141592653589793E1*cos
(2.0*0.3141592653589793E1*t)*R_hat*Radius_inner*Radius_innermost-
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*Radius_outer*Radius_innermost+0.3141592653589793E1*0.3141592653589793E1*R_hat*
R_hat*t*Radius_inner-0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*t*
Radius_outer+Beta1*0.3141592653589793E1*0.3141592653589793E1*Radius_inner*pow((
Beta0*(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/
0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0
)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/Sigma+pow(Theta_0+U0,4.0),0.25)-
Beta1*0.3141592653589793E1*0.3141592653589793E1*Radius_inner*Theta_0-sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*Radius_inner/2.0;
      MapleGenVar9 = sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*R_hat
*R_hat*Radius_outer/2.0-0.3141592653589793E1*0.3141592653589793E1*R_hat*
Radius_inner*Radius_innermost+0.3141592653589793E1*0.3141592653589793E1*R_hat*
Radius_outer*Radius_innermost+Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*
Radius_inner/2.0-Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*cos
(2.0*0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*Radius_outer/2.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*V0_hat*t*Radius_inner*Radius_innermost;
      MapleGenVar8 = MapleGenVar9-Beta0*0.3141592653589793E1*
0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat*V0_hat*t*
Radius_outer*Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat*V0_hat*Radius_inner*
Radius_innermost/2.0+Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1
*cos(2.0*0.3141592653589793E1*t)*R_hat*V0_hat*Radius_outer*Radius_innermost/2.0
+MapleGenVar7-Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t*
t*Radius_inner+Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t
*t*Radius_outer-Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*
V0_hat*t*t*Radius_inner/2.0;
      MapleGenVar9 = MapleGenVar8+Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*R_hat*V0_hat*t*t*Radius_outer/2.0-Beta0*pow(sin(2.0*
0.3141592653589793E1*t),2.0)*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*V0_hat
*Radius_inner/8.0+Beta0*pow(sin(2.0*0.3141592653589793E1*t),2.0)*cos(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*Radius_outer/8.0-Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
V0_hat*Radius_inner*Radius_innermost*Radius_innermost/2.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
V0_hat*Radius_outer*Radius_innermost*Radius_innermost/2.0+2.0*Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0*t*Radius_inner*
Radius_innermost;
      MapleGenVar6 = MapleGenVar9-2.0*Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*V0*t*Radius_outer*Radius_innermost+Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*t*Radius_inner*
Radius_innermost-Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*t
*Radius_outer*Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*V0*Radius_inner*Radius_innermost+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*V0*Radius_outer*
Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*
R_hat*V0_hat*Radius_inner*Radius_innermost/2.0+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*V0_hat*Radius_outer*
Radius_innermost/2.0;
      MapleGenVar7 = 1/(0.3141592653589793E1*0.3141592653589793E1)/Radius_inner
/Beta1;
      MapleGenVar5 = MapleGenVar6*MapleGenVar7;
      MapleGenVar6 = -pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.25)+Theta_0;
      MapleGenVar4 = MapleGenVar5+MapleGenVar6;
      MapleGenVar5 = 1/(Radius_outer-Radius_inner);
      MapleGenVar3 = MapleGenVar4*MapleGenVar5;
      MapleGenVar1 = MapleGenVar2*MapleGenVar3;
      MapleGenVar3 = -1.0;
      MapleGenVar5 = Alpha1;
      MapleGenVar7 = 1/(pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.75))*(-Beta0*R_hat*(1.0-cos(2.0*
0.3141592653589793E1*t))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0)-
Beta0*(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/
0.3141592653589793E1/2.0))*V0_hat*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1+2.0*R_hat*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1)/Sigma/4.0;
      MapleGenVar14 = pow(cos(2.0*0.3141592653589793E1*t),2.0)*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*Radius_inner-pow(sin(2.0*
0.3141592653589793E1*t),2.0)*0.3141592653589793E1*0.3141592653589793E1*R_hat*
R_hat*Radius_inner-pow(cos(2.0*0.3141592653589793E1*t),2.0)*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*Radius_outer+pow(sin(2.0*
0.3141592653589793E1*t),2.0)*0.3141592653589793E1*0.3141592653589793E1*R_hat*
R_hat*Radius_outer+0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*
Radius_inner-0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*Radius_outer
;
      MapleGenVar13 = MapleGenVar14-2.0*Beta0*0.3141592653589793E1*
0.3141592653589793E1*0.3141592653589793E1*sin(2.0*0.3141592653589793E1*t)*R_hat
*V0_hat*t*Radius_inner*Radius_innermost+2.0*Beta0*0.3141592653589793E1*
0.3141592653589793E1*0.3141592653589793E1*sin(2.0*0.3141592653589793E1*t)*R_hat
*V0_hat*t*Radius_outer*Radius_innermost+Beta1*0.3141592653589793E1*
0.3141592653589793E1*Radius_inner/pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0
*0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.75)*(-Beta0*R_hat*(1.0-cos(2.0*0.3141592653589793E1
*t))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0)-Beta0*(
Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/0.3141592653589793E1/
2.0))*V0_hat*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1+2.0*R_hat*sin
(2.0*0.3141592653589793E1*t)*0.3141592653589793E1)/Sigma/4.0-2.0*Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_inner+2.0*
Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_outer-
Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0_hat*t*
Radius_inner;
      MapleGenVar12 = MapleGenVar13+Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*R_hat*V0_hat*t*Radius_outer+Beta0*pow(sin(2.0*
0.3141592653589793E1*t),3.0)*0.3141592653589793E1*R_hat*R_hat*V0_hat*
Radius_inner/4.0-Beta0*pow(sin(2.0*0.3141592653589793E1*t),3.0)*
0.3141592653589793E1*R_hat*R_hat*V0_hat*Radius_outer/4.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*V0_hat*Radius_inner*Radius_innermost*Radius_innermost-
Beta0*0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*V0_hat*Radius_outer*Radius_innermost*Radius_innermost+
Beta0*0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*t*Radius_inner-Beta0*
0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*t*Radius_outer+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*V0*Radius_inner+2.0*
Beta0*cos(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*0.3141592653589793E1
*R_hat*R_hat*V0*t*Radius_inner-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*R_hat*V0*Radius_outer-2.0*Beta0*cos(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*
V0*t*Radius_outer+Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*
R_hat*R_hat*V0_hat*Radius_inner/2.0;
      MapleGenVar14 = MapleGenVar12-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*R_hat*V0_hat*Radius_outer/2.0-Beta0*sin(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0*Radius_inner*cos(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1+Beta0*sin(2.0*0.3141592653589793E1
*t)*R_hat*R_hat*V0*Radius_outer*cos(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1+Beta0*pow(cos(2.0*0.3141592653589793E1*t),2.0)*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0_hat*t*Radius_inner-
Beta0*pow(sin(2.0*0.3141592653589793E1*t),2.0)*0.3141592653589793E1*
0.3141592653589793E1*R_hat*R_hat*V0_hat*t*Radius_inner;
      MapleGenVar13 = MapleGenVar14-Beta0*pow(cos(2.0*0.3141592653589793E1*t),
2.0)*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0_hat*t*
Radius_outer+Beta0*pow(sin(2.0*0.3141592653589793E1*t),2.0)*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0_hat*t*Radius_outer-
Beta0*pow(cos(2.0*0.3141592653589793E1*t),2.0)*0.3141592653589793E1*
0.3141592653589793E1*R_hat*V0_hat*Radius_inner*Radius_innermost+Beta0*pow(sin(
2.0*0.3141592653589793E1*t),2.0)*0.3141592653589793E1*0.3141592653589793E1*
R_hat*V0_hat*Radius_inner*Radius_innermost+Beta0*pow(cos(2.0*
0.3141592653589793E1*t),2.0)*0.3141592653589793E1*0.3141592653589793E1*R_hat*
V0_hat*Radius_outer*Radius_innermost-Beta0*pow(sin(2.0*0.3141592653589793E1*t),
2.0)*0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*Radius_outer*
Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*pow(cos(2.0*
0.3141592653589793E1*t),2.0)*R_hat*R_hat*V0_hat*Radius_inner*
0.3141592653589793E1/2.0;
      MapleGenVar11 = MapleGenVar13+Beta0*sin(2.0*0.3141592653589793E1*t)*pow(
cos(2.0*0.3141592653589793E1*t),2.0)*R_hat*R_hat*V0_hat*Radius_outer*
0.3141592653589793E1/2.0+2.0*Beta0*0.3141592653589793E1*0.3141592653589793E1*
R_hat*V0*Radius_inner*Radius_innermost-2.0*Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*V0*Radius_outer*Radius_innermost+Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*Radius_inner*
Radius_innermost-Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*
Radius_outer*Radius_innermost-2.0*Beta0*cos(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0*Radius_inner*
Radius_innermost+2.0*Beta0*cos(2.0*0.3141592653589793E1*t)*0.3141592653589793E1
*0.3141592653589793E1*R_hat*V0*Radius_outer*Radius_innermost-2.0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*R_hat*Radius_inner+2.0*0.3141592653589793E1*0.3141592653589793E1*
0.3141592653589793E1*sin(2.0*0.3141592653589793E1*t)*R_hat*R_hat*t*Radius_inner
+2.0*0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
R_hat*R_hat*Radius_outer-2.0*0.3141592653589793E1*0.3141592653589793E1*
0.3141592653589793E1*sin(2.0*0.3141592653589793E1*t)*R_hat*R_hat*t*Radius_outer
-2.0*0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*R_hat*Radius_inner*Radius_innermost+2.0*
0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*sin(2.0*
0.3141592653589793E1*t)*R_hat*Radius_outer*Radius_innermost;
      MapleGenVar12 = 1/(0.3141592653589793E1*0.3141592653589793E1)/
Radius_inner/Beta1;
      MapleGenVar10 = MapleGenVar11*MapleGenVar12;
      MapleGenVar11 = -1/(pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.75))*(-Beta0*R_hat*(1.0-cos(2.0*
0.3141592653589793E1*t))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0)-
Beta0*(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/
0.3141592653589793E1/2.0))*V0_hat*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1+2.0*R_hat*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1)/Sigma/4.0;
      MapleGenVar9 = MapleGenVar10+MapleGenVar11;
      MapleGenVar10 = (r-Radius_inner)/(Radius_outer-Radius_inner);
      MapleGenVar8 = MapleGenVar9*MapleGenVar10;
      MapleGenVar6 = MapleGenVar7+MapleGenVar8;
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      t0 = MapleGenVar1+MapleGenVar2;



   }

  source = t0;
   
 }

 /// Exact solution as a Vector
 void get_exact_u(const double& time, const Vector<double>& x, 
                  Vector<double>& u)
 {
  double t=time;
  double t0=0.0;

  double MapleGenVar1 = 0.0;     
  double MapleGenVar4 = 0.0;
  double MapleGenVar6 = 0.0;
  double MapleGenVar9 = 0.0;
  double MapleGenVar8 = 0.0;
  double MapleGenVar7 = 0.0;
  double MapleGenVar5 = 0.0;
  double MapleGenVar3 = 0.0;
  double MapleGenVar10 = 0.0;
  double MapleGenVar2  = 0.0;

  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  if (r<0.5*(Radius_inner+Radius_innermost))
   {
//#--------------------------------------------------------------------
//# Solution in inner region
//#--------------------------------------------------------------------
//> C(eval(u0));

      t0 = U0+(r*r-pow(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t
)/0.3141592653589793E1/2.0),2.0))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*
t))/2.0)/2.0;


   }
  else
   {

//#--------------------------------------------------------------------
//# Solution in outer region
//#--------------------------------------------------------------------
//> C(eval(u1));

      MapleGenVar1 = pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.25);
      MapleGenVar3 = -Theta_0;
      MapleGenVar9 = -Beta0*0.3141592653589793E1*0.3141592653589793E1*cos(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*t*Radius_inner/2.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*R_hat*V0_hat*t*t*Radius_outer/2.0+Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_inner-Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*V0*t*Radius_outer+
Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*V0_hat*t
*Radius_inner/2.0-Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*
R_hat*R_hat*V0_hat*t*Radius_outer/2.0-Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0*Radius_inner/4.0+Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0*Radius_outer/4.0-Beta0*pow(sin(2.0*0.3141592653589793E1*t)
,2.0)*R_hat*R_hat*V0_hat*Radius_inner/8.0+Beta0*pow(sin(2.0*
0.3141592653589793E1*t),2.0)*R_hat*R_hat*V0_hat*Radius_outer/8.0-Beta0*
0.3141592653589793E1*0.3141592653589793E1*V0*Radius_inner*Radius_innermost*
Radius_innermost+Beta0*0.3141592653589793E1*0.3141592653589793E1*V0*
Radius_outer*Radius_innermost*Radius_innermost;
      MapleGenVar10 = MapleGenVar9-Beta0*0.3141592653589793E1*
0.3141592653589793E1*V0_hat*Radius_inner*Radius_innermost*Radius_innermost/2.0+
Beta0*0.3141592653589793E1*0.3141592653589793E1*V0_hat*Radius_outer*
Radius_innermost*Radius_innermost/2.0-0.3141592653589793E1*0.3141592653589793E1
*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*t*Radius_inner+
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*R_hat*t*Radius_outer+sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*cos(
2.0*0.3141592653589793E1*t)*R_hat*R_hat*Radius_inner/2.0-sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
R_hat*R_hat*Radius_outer/2.0;
      MapleGenVar8 = MapleGenVar10+0.3141592653589793E1*0.3141592653589793E1*
cos(2.0*0.3141592653589793E1*t)*R_hat*Radius_inner*Radius_innermost-
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat
*Radius_outer*Radius_innermost+0.3141592653589793E1*0.3141592653589793E1*R_hat*
R_hat*t*Radius_inner-0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*t*
Radius_outer+Beta1*0.3141592653589793E1*0.3141592653589793E1*Radius_inner*pow((
Beta0*(Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/
0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0
)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/Sigma+pow(Theta_0+U0,4.0),0.25)-
Beta1*0.3141592653589793E1*0.3141592653589793E1*Radius_inner*Theta_0-sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*R_hat*Radius_inner/2.0;
      MapleGenVar10 = MapleGenVar8+sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*R_hat*Radius_outer/2.0-0.3141592653589793E1*
0.3141592653589793E1*R_hat*Radius_inner*Radius_innermost+0.3141592653589793E1*
0.3141592653589793E1*R_hat*Radius_outer*Radius_innermost+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
R_hat*R_hat*V0_hat*t*Radius_inner/2.0-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*t*
Radius_outer/2.0;
      MapleGenVar9 = MapleGenVar10+Beta0*0.3141592653589793E1*
0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*R_hat*V0_hat*t*
Radius_inner*Radius_innermost-Beta0*0.3141592653589793E1*0.3141592653589793E1*
cos(2.0*0.3141592653589793E1*t)*R_hat*V0_hat*t*Radius_outer*Radius_innermost-
Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*cos(2.0*
0.3141592653589793E1*t)*R_hat*V0_hat*Radius_inner*Radius_innermost/2.0+Beta0*
sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*cos(2.0*
0.3141592653589793E1*t)*R_hat*V0_hat*Radius_outer*Radius_innermost/2.0-Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t*t*Radius_inner+Beta0
*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0*t*t*Radius_outer-
Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*R_hat*V0_hat*t*t*
Radius_inner/2.0;
      MapleGenVar10 = MapleGenVar9+Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*R_hat*V0_hat*t*t*Radius_outer/2.0-Beta0*pow(sin(2.0*
0.3141592653589793E1*t),2.0)*cos(2.0*0.3141592653589793E1*t)*R_hat*R_hat*V0_hat
*Radius_inner/8.0+Beta0*pow(sin(2.0*0.3141592653589793E1*t),2.0)*cos(2.0*
0.3141592653589793E1*t)*R_hat*R_hat*V0_hat*Radius_outer/8.0-Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
V0_hat*Radius_inner*Radius_innermost*Radius_innermost/2.0+Beta0*
0.3141592653589793E1*0.3141592653589793E1*cos(2.0*0.3141592653589793E1*t)*
V0_hat*Radius_outer*Radius_innermost*Radius_innermost/2.0+2.0*Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0*t*Radius_inner*
Radius_innermost;
      MapleGenVar7 = MapleGenVar10-2.0*Beta0*0.3141592653589793E1*
0.3141592653589793E1*R_hat*V0*t*Radius_outer*Radius_innermost+Beta0*
0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*t*Radius_inner*
Radius_innermost-Beta0*0.3141592653589793E1*0.3141592653589793E1*R_hat*V0_hat*t
*Radius_outer*Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*
0.3141592653589793E1*R_hat*V0*Radius_inner*Radius_innermost+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*V0*Radius_outer*
Radius_innermost-Beta0*sin(2.0*0.3141592653589793E1*t)*0.3141592653589793E1*
R_hat*V0_hat*Radius_inner*Radius_innermost/2.0+Beta0*sin(2.0*
0.3141592653589793E1*t)*0.3141592653589793E1*R_hat*V0_hat*Radius_outer*
Radius_innermost/2.0;
      MapleGenVar8 = 1/(0.3141592653589793E1*0.3141592653589793E1)/Radius_inner
/Beta1;
      MapleGenVar6 = MapleGenVar7*MapleGenVar8;
      MapleGenVar7 = -pow((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0),0.25)+Theta_0;
      MapleGenVar5 = MapleGenVar6+MapleGenVar7;
      MapleGenVar6 = (r-Radius_inner)/(Radius_outer-Radius_inner);
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3+MapleGenVar4;
      t0 = MapleGenVar1+MapleGenVar2;

   }

  u[0] = t0;

 }
 
 /// Melt flux for exact solution (to test melting without feedback
 /// with S.B.)
 double melt_flux(const double& t)
 {
  double t0=0;
  
//#--------------------------------------------------------------------
//# Melt flux
//#--------------------------------------------------------------------
//> C(eval(melt_flux));

      t0 = R_hat*(1.0-cos(2.0*0.3141592653589793E1*t));
  
      return t0;
 }


 /// Exact radius of inner circle
 double radius(const double& t)
 {
  double t0=0.0;

//#--------------------------------------------------------------------
//# Radius
//#--------------------------------------------------------------------
//> C(eval(r0(t)));   

      t0 = Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/
0.3141592653589793E1/2.0);

  
  return t0;
 }

 
 /// Incoming sb radiation on inside 
 double incoming_sb_inside(const double& t)
 {
  double t0=0.0;
  
// #--------------------------------------------------------------------
// # Incoming sb radiation inside
// #--------------------------------------------------------------------
// > C(eval(Sigma*((Theta_0+U1)^4)));

      t0 = Sigma*((Beta0*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))*(V0+V0_hat*(1.0+cos(2.0*
0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*0.3141592653589793E1*t)))/
Sigma+pow(Theta_0+U0,4.0));


  
  return t0;
 }



 /// Incoming sb radiation on outside 
 double incoming_sb_outside(const double& t)
 {
  double t0=0.0;

// #--------------------------------------------------------------------
// # Incoming sb radiation outside
// #--------------------------------------------------------------------
// > C(eval((Theta_0+U0)^4*r0(t)/Radius_inner+
// >        (Theta_0+U1)^4*(1-(r0(t)/Radius_inner))));

      t0 = Sigma*(pow(Theta_0+U0,4.0)*(Radius_innermost-R_hat*(t-sin(2.0*
0.3141592653589793E1*t)/0.3141592653589793E1/2.0))/Radius_inner+((Beta0*(
Radius_innermost-R_hat*(t-sin(2.0*0.3141592653589793E1*t)/0.3141592653589793E1/
2.0))*(V0+V0_hat*(1.0+cos(2.0*0.3141592653589793E1*t))/2.0)+R_hat*(1.0-cos(2.0*
0.3141592653589793E1*t)))/Sigma+pow(Theta_0+U0,4.0))*(1.0-(Radius_innermost-
R_hat*(t-sin(2.0*0.3141592653589793E1*t)/0.3141592653589793E1/2.0))/
Radius_inner));
  
  return t0;
 }

} // end of namespace


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class 
//=====================================================================
template<class ELEMENT> 
class StefanBoltzmannProblem : public Problem
{

public:
 
 /// Constructor
 StefanBoltzmannProblem();
 
 /// Destructor (empty)
 ~StefanBoltzmannProblem(){}

 /// Update the problem specs before solve: Update Stefan Boltzmann
 /// radiation
 void actions_before_newton_solve()
  {
   oomph_info << "Re-setting up Stefan Boltzmann radiation\n";
   setup_sb_radiation();

   // Re-setup equation numbering scheme
   cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  
  } 

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Update the problem specs before next timestep: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_implicit_timestep()
  {
   // Update pinned boundary values
   for(unsigned b=0;b<2;b++) // hierher enumerate
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       double time=time_pt()->time();
       Vector<double> u(1);
       GlobalParameters::get_exact_u(time,x,u);
       nod_pt->set_value(0,u[0]);
      }
    }
  }

 /// Setup Stefan Boltzmann radiation
 void setup_sb_radiation();

 /// Create the Stefan Boltzmann elements
 void create_sb_elements()
  {
   // Loop over inner boundaries of outer annulus
   for (unsigned b=2;b<=3;b++) // hierher enumerate
    {
     
     // How many bulk elements are adjacent to boundary b?
     unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
       
       // Create flux element
       StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT> *el_pt = 
        new StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>(bulk_elem_pt,
                                                            face_index);
       
       // Set different integration scheme
       if (CommandLineArgs::command_line_flag_has_been_set("--nintpt"))
        {
         el_pt->set_integration_scheme
          (new MyIntegral(GlobalParameters::Nintpt));
        }
       
       // Set non-dim Stefan Boltzmann constant
       el_pt->sigma_pt()= &GlobalParameters::Sigma;
       
       // Set zero-centrigrade offset in Stefan Boltzmann law
       el_pt->theta_0_pt()= &GlobalParameters::Theta_0;
       
       // Add to mesh
       Unsteady_heat_flux_mesh_pt->add_element_pt(el_pt);
      }  
    }
  }


 /// Create the melt elements and also identify nodes to be pinned
 /// on inner boundary
 void create_melt_elements()
  {
   // Storage for rightmost, leftmost and topmost node on inner boundary
   double x_min=DBL_MAX;
   double x_max=-DBL_MAX;
   double y_max=-DBL_MAX;
   SolidNode* x_min_node_pt=0;
   SolidNode* x_max_node_pt=0;
   SolidNode* y_max_node_pt=0;
   
   // Loop over boundaries of inner circle
   for (unsigned b=4;b<=5;b++) // hierher enumerate
    {
     // Loop over the bulk elements adjacent to boundary b?
     unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_pt(b,e));

       //What is the face index of element e along boundary b
       int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
       
       // // Build the corresponding melt element
       // SurfaceMeltElement<ELEMENT>* melt_element_pt = new 
       //  SurfaceMeltElement<ELEMENT>(bulk_elem_pt,face_index);
       

       // Build the corresponding melt element
       StefanBoltzmannMeltElement<ELEMENT>* melt_element_pt = new 
        StefanBoltzmannMeltElement<ELEMENT>(bulk_elem_pt,face_index);
       
       //Add the melt element to the surface mesh
       Surface_melt_mesh_pt->add_element_pt(melt_element_pt);
              
       // Set melt temperature
       melt_element_pt->melt_temperature_pt()=
        &GlobalParameters::Melt_temperature;
       
       // Set different integration scheme
       if (CommandLineArgs::command_line_flag_has_been_set("--nintpt"))
        {
         melt_element_pt->set_integration_scheme
          (new MyIntegral(GlobalParameters::Nintpt));
        }
       
       // Set non-dim Stefan Boltzmann constant
       melt_element_pt->sigma_pt()= &GlobalParameters::Sigma;
       
       // Set zero-centrigrade offset in Stefan Boltzmann law
       melt_element_pt->theta_0_pt()= &GlobalParameters::Theta_0;
       
       // Find rightmost/leftmost/topmost node on inner boundary
       if (b==4)
        {
         unsigned nnod=melt_element_pt->nnode();
         for (unsigned j=0;j<nnod;j++)
          {
           SolidNode* nod_pt=
            dynamic_cast<SolidNode*>(melt_element_pt->node_pt(j));
           double x=nod_pt->x(0);
           double y=nod_pt->x(1);
           if (x<x_min) 
            {
             x_min=x;
             x_min_node_pt=nod_pt;
            }
           if (x>x_max) 
            {
             x_max=x;
             x_max_node_pt=nod_pt;
            }
           if (y>y_max) 
            {
             y_max=y;
             y_max_node_pt=nod_pt;
            }
          }
        }
      }
    } //end of loop over bulk elements adjacent to boundary b
  
   // Pin to avoid rigid body translation/rotation of inner circle
   oomph_info << "Pinning vertical   displacement at " 
              << x_min_node_pt->x(0) << " " 
              << x_min_node_pt->x(1) << std::endl;
   x_min_node_pt->pin_position(1);
   oomph_info << "Pinning vertical   displacement at " 
              << x_max_node_pt->x(0) << " " 
              << x_max_node_pt->x(1) << std::endl;
   x_max_node_pt->pin_position(1);
   oomph_info << "Pinning horizontal displacement at " 
              << y_max_node_pt->x(0) << " " 
              << y_max_node_pt->x(1) << std::endl;
   y_max_node_pt->pin_position(0);

  }

 /// Doc the solution. 
 void doc_solution();
 
private:

 /// Pointer to the "bulk" mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to surface mesh of radiative flux elements
 Mesh* Unsteady_heat_flux_mesh_pt;

 /// Pointer to the surface melt mesh
 Mesh* Surface_melt_mesh_pt;

 /// Trace file
 ofstream Trace_file;

 /// Integration point of leftmost point on outer boundary
 unsigned Ipt_leftmost_outer;

 /// Element containing leftmost point on outer boundary
 StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>* Leftmost_outer_el_pt;

 /// Integration point of rightmost point on inner boundary
 unsigned Ipt_rightmost_inner;

 /// Element containing rightmost point on inner boundary
 StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>* Rightmost_inner_el_pt;

 /// DocInfo object stores flags/labels for where the output gets written to
 DocInfo Doc_info;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor
//========================================================================
template<class ELEMENT>
StefanBoltzmannProblem<ELEMENT>::StefanBoltzmannProblem()
{ 

 // Set output directory
 Doc_info.set_directory(GlobalParameters::Directory);

 // Open trace file
 Trace_file.open((GlobalParameters::Directory+"/trace.dat").c_str());
 
 // Number of segments in (half) the innermost boundary
 unsigned n_innermost=20;

 // Scale number of segments for convergence test
 unsigned n_inner=unsigned(GlobalParameters::Radius_inner/
                           GlobalParameters::Radius_innermost*
                           double(n_innermost));
 unsigned n_outer=20; 

 // Create circle representing outer boundary
 double a= GlobalParameters::Radius_outer;
 double x_c=0.0;
 double y_c=0.0;
 Circle* outer_circle_pt=new Circle(x_c,y_c,a);

 // Create circle representing inner boundary
 a=GlobalParameters::Radius_inner;
 x_c=0.0;
 y_c=0.0;
 Circle* inner_circle_pt=new Circle(x_c,y_c,a);

 // Create circle representing boundary of central region
 a=GlobalParameters::Radius_innermost;
 Circle* central_circle_pt=new Circle(x_c,y_c,a);


 // Outer boundary
 //---------------

 // Provide storage for pointers to the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
 
 // First bit
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned boundary_id=0;
 outer_curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  outer_circle_pt,zeta_start,zeta_end,n_outer,boundary_id);
 
 // Second bit
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
 boundary_id=1;
 outer_curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
  outer_circle_pt,zeta_start,zeta_end,n_outer,boundary_id);
 
 // Combine to curvilinear boundary and define the
 // outer boundary
 TriangleMeshClosedCurve* outer_boundary_pt=
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);


 // Inner circular boundaries
 //--------------------------
 Vector<TriangleMeshCurveSection*> inner_boundary_line_pt(2);
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start = 0.0;
 double s_end   = MathematicalConstants::Pi;
 boundary_id = 2;
 inner_boundary_line_pt[0]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_inner, 
                            boundary_id);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = MathematicalConstants::Pi;
 s_end   = 2.0*MathematicalConstants::Pi;
 boundary_id = 3;
 inner_boundary_line_pt[1]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_inner,
                            boundary_id);
 
 // Combine to hole
 Vector<TriangleMeshClosedCurve*> internal_closed_curve_pt;
 Vector<double> hole_coords(2);
 hole_coords[0]=0.5*(GlobalParameters::Radius_inner+
                     GlobalParameters::Radius_innermost);
 hole_coords[1]=0.0;
 internal_closed_curve_pt.push_back(
  new TriangleMeshClosedCurve(inner_boundary_line_pt,
                              hole_coords));
  
 
 // Boundary of innermost region
 //------------------------------
 Vector<TriangleMeshCurveSection*> central_boundary_line_pt(2);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = 0.0;
 s_end   = MathematicalConstants::Pi;
 boundary_id = 4;
 central_boundary_line_pt[0]=
  new TriangleMeshCurviLine(central_circle_pt,
                            s_start,
                            s_end,
                            n_innermost,
                            boundary_id);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = MathematicalConstants::Pi;
 s_end   = 2.0*MathematicalConstants::Pi;
 boundary_id = 5;
 central_boundary_line_pt[1]=
  new TriangleMeshCurviLine(central_circle_pt,
                            s_start,
                            s_end,
                            n_innermost,
                            boundary_id);
 
 // Define hole
 internal_closed_curve_pt.push_back(
  new TriangleMeshClosedCurve(central_boundary_line_pt));

 // Use the TriangleMeshParameters object for helping on the manage 
 // of the TriangleMesh parameters. The only parameter that needs to take 
 // is the outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() = internal_closed_curve_pt;

 // Target element size in bulk mesh
 triangle_mesh_parameters.element_area() = GlobalParameters::Target_area;
 
 /// Identify outer annulus as region 1
 Vector<double> outer_annulus_region_coord(2); 
 outer_annulus_region_coord[0]=0.0;
 outer_annulus_region_coord[1]=0.5*(GlobalParameters::Radius_inner+
                                    GlobalParameters::Radius_outer);
 triangle_mesh_parameters.add_region_coordinates(1,outer_annulus_region_coord);
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Build "bulk" mesh
 Bulk_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                       time_stepper_pt());


 // Complete setup of bulk elements

 // Central region
 unsigned r=0;
 unsigned nel=Bulk_mesh_pt->nregion_element(r);
 for (unsigned e=0;e<nel;e++)
  {
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->region_element_pt(r,e));

   // Thermal inertia
   el_pt->alpha_pt()=&GlobalParameters::Alpha0;

   // Thermal conductivity
   el_pt->beta_pt()=&GlobalParameters::Beta0;

   //Set the source function pointer
   if (!CommandLineArgs::command_line_flag_has_been_set("--no_source"))
    {
     el_pt->source_fct_pt() = &GlobalParameters::get_source;
    }

   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    GlobalParameters::Constitutive_law_pt;

   // Set density to zero
   el_pt->lambda_sq_pt()=&GlobalParameters::Lambda_sq;

   // Disable inertia
   el_pt->disable_inertia();
  }


 // Outer annular region
 r=1;
 nel=Bulk_mesh_pt->nregion_element(r);
 for (unsigned e=0;e<nel;e++)
  {
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->region_element_pt(r,e));

   // Thermal inertia
   el_pt->alpha_pt()=&GlobalParameters::Alpha1;

   // Thermal conductivity
   el_pt->beta_pt()=&GlobalParameters::Beta1;

   //Set the source function pointer
   if (!CommandLineArgs::command_line_flag_has_been_set("--no_source"))
    {
     el_pt->source_fct_pt() = &GlobalParameters::get_source;
    }

   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    GlobalParameters::Constitutive_law_pt;

   // Set density to zero
   el_pt->lambda_sq_pt()=&GlobalParameters::Lambda_sq;

   // Disable inertia
   el_pt->disable_inertia();

   // Pin the positional dofs of all nodes
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(j));
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
    }
  }

 // Doc use of integration scheme
 if (CommandLineArgs::command_line_flag_has_been_set("--nintpt"))
  {
   oomph_info << "Setting new integration scheme with nintpt="
              << GlobalParameters::Nintpt << std::endl;
  }
 else
  {
   oomph_info << "Using Gauss scheme" << std::endl;
  }

 // Create mesh with Stefan Boltmann radiation elements
 Unsteady_heat_flux_mesh_pt=new Mesh;
 create_sb_elements();
 
 // Create the surface melt mesh
 Surface_melt_mesh_pt=new Mesh;
 create_melt_elements();
 
 // Setup Stefan Boltzmann radiation
 setup_sb_radiation();
 
 // Build the entire mesh from its submeshes
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_melt_mesh_pt);
 add_sub_mesh(Unsteady_heat_flux_mesh_pt);
 build_global_mesh();
 
 // Assign exact solution as initial guess
 unsigned nnod=Bulk_mesh_pt->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=Bulk_mesh_pt->node_pt(j);
   Vector<double> x(2);
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   double time=time_pt()->time();
   Vector<double> u(1);
   GlobalParameters::get_exact_u(time,x,u);
   nod_pt->set_value(0,u[0]);
  }
 
 // Pin boundary values for temperature on very outside
 for(unsigned b=0;b<2;b++)
  {
   unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
     nod_pt->pin(0);       
    }
  }


 // Find rightmost point on inner boundary; leftmost point on outer boundary
 //-------------------------------------------------------------------------
 {
  double x_leftmost_outer = DBL_MAX;
  Vector<double> s(1);
  Vector<double> r(2);
  unsigned nel=Unsteady_heat_flux_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>* el_pt=
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>
     (Unsteady_heat_flux_mesh_pt->element_pt(e));
    unsigned nintpt=el_pt->integral_pt()->nweight();
    for (unsigned ipt=0;ipt<nintpt;ipt++)
     {
      // Local coordinate of integration point
      s[0]=el_pt->integral_pt()->knot(ipt,0);
      
      // Get position
      el_pt->interpolated_x(s,r);
      
      if (r[0]<x_leftmost_outer)
       {
        x_leftmost_outer=r[0];
        Leftmost_outer_el_pt=el_pt;
        Ipt_leftmost_outer=ipt;
       }
     }
   }
  
  
  double x_rightmost_inner=-DBL_MAX;
  nel=Surface_melt_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    StefanBoltzmannMeltElement<ELEMENT>* el_pt=
     dynamic_cast<StefanBoltzmannMeltElement<ELEMENT>*>
     (Surface_melt_mesh_pt->element_pt(e));
    unsigned nintpt=el_pt->integral_pt()->nweight();
    for (unsigned ipt=0;ipt<nintpt;ipt++)
     {
      // Local coordinate of integration point
      s[0]=el_pt->integral_pt()->knot(ipt,0);
      
      // Get position
      el_pt->interpolated_x(s,r);
      
      if (r[0]>x_rightmost_inner)
       {
        x_rightmost_inner=r[0];
        Rightmost_inner_el_pt=el_pt;
        Ipt_rightmost_inner=ipt;
       }
     }
   }
 }


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  
} // end of constructor


//=====================start_of_setup_sb==================================
/// Setup Stefan Boltzmann radiation
//========================================================================
template<class ELEMENT>
void StefanBoltzmannProblem<ELEMENT>::setup_sb_radiation()
{
 // Identify face elements on potentially sun-exposed boundaries
 Vector<FiniteElement*> shielding_face_element_pt;
 unsigned nel=Unsteady_heat_flux_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Add to shielding surface
   shielding_face_element_pt.push_back(
    Unsteady_heat_flux_mesh_pt->finite_element_pt(e));     
  }
 nel=Surface_melt_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Add to shielding surface
   shielding_face_element_pt.push_back(
    Surface_melt_mesh_pt->finite_element_pt(e));     
  }

 // Doc it?
 bool doc_it=false;
 if (doc_it)
  {
   StefanBoltzmannHelper::Nsample=2;
   StefanBoltzmannHelper::Nx_bin=20;
   StefanBoltzmannHelper::Ny_bin=20;
  }

 //Setup mutual visibility for all Stefan Boltzmann elements
 StefanBoltzmannHelper::setup_stefan_boltzmann_visibility(
  shielding_face_element_pt);
 
 // Doc the populated bins
 bool doc_bins=true;
 if (doc_bins)
  {
   ofstream some_file;
   char filename[100];

   oomph_info << "Docing bins for step " << Doc_info.number() << std::endl;

   sprintf(filename,"%s/populated_bins%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   StefanBoltzmannHelper::doc_bins(some_file);
   some_file.close();
   
   // Doc the sample points used to assess visibility 
   sprintf(filename,"%s/sample_points%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   StefanBoltzmannHelper::doc_sample_points(some_file,
                                            shielding_face_element_pt);
   some_file.close();
  }

}
 
 

//=====================start_of_doc=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StefanBoltzmannProblem<ELEMENT>::doc_solution()
{ 

 oomph_info << "Outputting step " << Doc_info.number()
            << " for time " << time_pt()->time() << std::endl;

 ofstream some_file,some_file2;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 
 
 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 
 // Output solution 
 //-----------------
 sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 unsigned npts_coarse=2;
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();
 

 // Output Stefan Boltzmann radiation rays on rightmost point
 //----------------------------------------------------------
 // on inner boundary; leftmost point on outer boundary
 //----------------------------------------------------

 // Output rays
 sprintf(filename,"%s/stefan_boltzmann_rays_left%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Leftmost_outer_el_pt->
  output_stefan_boltzmann_radiation_rays(some_file,Ipt_leftmost_outer);
 some_file.close();
 

 // Output rays
 sprintf(filename,"%s/stefan_boltzmann_rays_right%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Rightmost_inner_el_pt->
  output_stefan_boltzmann_radiation_rays(some_file,Ipt_rightmost_inner);
 some_file.close();
 

 // Output exact outer radius
 //--------------------------
 double r=GlobalParameters::radius(time_pt()->time());
 sprintf(filename,"%s/exact_melt_surface%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned np=300;
 for (unsigned i=0;i<np;i++)
  {
   double phi=2.0*MathematicalConstants::Pi*double(i)/double(np-1);
   some_file << r*cos(phi) << " " << r*sin(phi) << std::endl;
  }
 some_file.close();


// Output heat flux (x,y,u,flux,n_x,n_y)
//--------------------------------------
sprintf(filename,"%s/flux%i.dat",
        Doc_info.directory().c_str(),
        Doc_info.number());
some_file.open(filename);
Unsteady_heat_flux_mesh_pt->output(some_file,5);
some_file.close();


// Output Stefan Boltzmann radiation (x, y, net incoming, incoming, 
//-----------------------------------------------------------------
// outgoing, n_x, n_y)
//--------------------
double av_inc_outer=0.0;
double av_inc_inner=0.0;
double inc_outer_exact=
 GlobalParameters::incoming_sb_outside(time_pt()->time());
double inc_inner_exact= 
 GlobalParameters::incoming_sb_inside(time_pt()->time());
{
 double average_incoming_inner=0.0;
 double average_incoming_outer=0.0;
 unsigned count_inner=0;
 unsigned count_outer=0;
 Vector<double> s(1);
 Vector<double> x(2);
 Vector<double> unit_normal(2);
 sprintf(filename,"%s/sb_radiation%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Unsteady_heat_flux_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>* el_pt=
   dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
    Unsteady_heat_flux_mesh_pt->element_pt(e));

   // Get contribution to avergage (over integration points) of incoming rad
   unsigned nintpt=el_pt->integral_pt()->nweight();
   for (unsigned ipt=0;ipt<nintpt;ipt++)
    {
     // Local coordinate of integration point
     s[0]=el_pt->integral_pt()->knot(ipt,0);
     
     // Get Eulerian coordinates
     el_pt->interpolated_x(s,x);
          
     // Outer unit normal
     el_pt->outer_unit_normal(s,unit_normal);
        
     // Get incoming radiation
     double my_rate=el_pt->incoming_stefan_boltzmann_radiation(ipt,x,
                                                               unit_normal);

     // Add to average
     average_incoming_outer+=my_rate;
     count_outer++;
    }

   // Output
   el_pt->output_stefan_boltzmann_radiation(some_file);
  }
 nel=Surface_melt_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Cast to element
   StefanBoltzmannMeltElement<ELEMENT>* el_pt=
   dynamic_cast<StefanBoltzmannMeltElement<ELEMENT>*>(
    Surface_melt_mesh_pt->element_pt(e));

   // Get contribution to avergage (over integration points) of incoming rad
   unsigned nintpt=el_pt->integral_pt()->nweight();
   for (unsigned ipt=0;ipt<nintpt;ipt++)
    {
     // Local coordinate of integration point
     s[0]=el_pt->integral_pt()->knot(ipt,0);
     
     // Get Eulerian coordinates
     el_pt->interpolated_x(s,x);
          
     // Outer unit normal
     el_pt->outer_unit_normal(s,unit_normal);
        
     // Get incoming radiation
     double my_rate=el_pt->incoming_stefan_boltzmann_radiation(ipt,x,
                                                               unit_normal);

     // Add to average
     average_incoming_inner+=my_rate;
     count_inner++;
    }
   
   // Output
   el_pt->output_stefan_boltzmann_radiation(some_file);
  }
 some_file.close();

 av_inc_inner=average_incoming_inner/double(count_inner);
 oomph_info << "Average inc. radiation (inner): " << av_inc_inner
            << " ; exact: " << inc_inner_exact 
            << " ; abs difference: " << fabs(inc_inner_exact-av_inc_inner);
 if (std::max(inc_inner_exact,av_inc_inner)!=0)
  {
   oomph_info 
    << " ; rel difference: " 
    << fabs(inc_inner_exact-av_inc_inner)/
    double(std::max(inc_inner_exact,av_inc_inner))*100.0
    << " % ";
  }
 oomph_info << std::endl;
 
 av_inc_outer=average_incoming_outer/double(count_outer);
 oomph_info << "Average inc. radiation (outer): " << av_inc_outer
            << " ; exact: " << inc_outer_exact 
            << " ; abs difference: " << fabs(inc_outer_exact-av_inc_outer);
 if (std::max(inc_outer_exact,av_inc_outer)!=0)
  {
   oomph_info 
    << " ; rel difference: " 
    << fabs(inc_outer_exact-av_inc_outer)/
    double(std::max(inc_outer_exact,av_inc_outer))*100.0
    << " % ";
  }
 oomph_info << std::endl;

}

// Output melt (x, y, u, net incoming flux, melt rate, n_x, n_y)
//--------------------------------------------------------------
double av_melt=0.0;
double melt_exact=GlobalParameters::melt_flux(time_pt()->time());
{
 double average_melt=0.0;
 unsigned count=0;
 Vector<double> s(1);
 sprintf(filename,"%s/melt%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_melt_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   StefanBoltzmannMeltElement<ELEMENT>* el_pt=
   dynamic_cast<StefanBoltzmannMeltElement<ELEMENT>*>(
    Surface_melt_mesh_pt->element_pt(e));

   // Get contribution to avergage (over integration points) of melt rate
   unsigned nintpt=el_pt->integral_pt()->nweight();
   for (unsigned ipt=0;ipt<nintpt;ipt++)
    {
     // Local coordinate of integration point
     s[0]=el_pt->integral_pt()->knot(ipt,0);
     
     // Get melt rate
     double my_melt=0.0;
     el_pt->interpolated_melt_rate(s,my_melt);

     // Add to average
     average_melt+=my_melt;
     count++;
    }

   // Output
   el_pt->output_melt(some_file);
   
  }
 some_file.close();
 
 av_melt=average_melt/double(count);
 oomph_info << "Average melt rate: " << av_melt
            << " ; exact: " << melt_exact 
            << " ; abs difference: " << fabs(melt_exact-av_melt);
 if (std::max(melt_exact,av_melt)!=0)
  {
   oomph_info 
    << " ; rel difference: " 
    << fabs(melt_exact-av_melt)/double(std::max(melt_exact,av_melt))*100.0
    << " % ";
  }
 oomph_info << std::endl;

}

// Get average radius
double av_radius=0.0;
std::map<Node*,bool> done;
unsigned count=0;
unsigned nel=Surface_melt_mesh_pt->nelement();
for (unsigned e=0;e<nel;e++)
 {
  FiniteElement* el_pt=Surface_melt_mesh_pt->finite_element_pt(e);
  unsigned nnod=el_pt->nnode();
  for (unsigned j=0;j<nnod;j++)
   {
    Node* nod_pt=el_pt->node_pt(j);
    if (!done[nod_pt])
     {
      done[nod_pt]=true;
      count++;
      av_radius+=sqrt(nod_pt->x(0)*nod_pt->x(0)+
                      nod_pt->x(1)*nod_pt->x(1));
     }
   }
 }
av_radius/=double(count);


Trace_file << time_pt()->time() << " " 
           << av_radius << " " 
           << GlobalParameters::radius(time_pt()->time()) << " "
           << av_inc_outer << " " 
           << inc_outer_exact << " " 
           << av_inc_inner << " " 
           << inc_inner_exact << " " 
           << av_melt << " " 
           << melt_exact << " ";

// Output exact solution 
//----------------------
sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
        Doc_info.number());
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                         GlobalParameters::get_exact_u); 
some_file.close();

// Doc error
//----------
double error,norm;
sprintf(filename,"%s/error%i.dat",Doc_info.directory().c_str(),
        Doc_info.number());
some_file.open(filename);
Bulk_mesh_pt->compute_error(some_file,
                            GlobalParameters::get_exact_u,
                            time_pt()->time(),
                            error,norm); 
some_file.close();

// Doc solution and error
//-----------------------
cout << "error: " << error << std::endl; 
cout << "norm : " << norm << std::endl << std::endl;

Trace_file  << error<< " " 
<< norm << " ";

Trace_file << std::endl;

//Increment counter for solutions 
Doc_info.number()++;

} // end of doc



//==========start_of_main=================================================
/// Solve 2D problem 
//========================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Number of integration points for new integration scheme. 
 // Use normal Gauss rule if not specified.
 CommandLineArgs::specify_command_line_flag("--nintpt",
                                            &GlobalParameters::Nintpt);

 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &GlobalParameters::Directory);

 // Timestep
 double dt=0.05;
 CommandLineArgs::specify_command_line_flag("--dt",&dt);

 // Max time
 double t_max=2.0;
 CommandLineArgs::specify_command_line_flag("--t_max",&t_max);

 // No source terms
 CommandLineArgs::specify_command_line_flag("--no_source");
  
 // Target area
 CommandLineArgs::specify_command_line_flag("--el_area",
                                            &GlobalParameters::Target_area);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
 // Create generalised Hookean constitutive equations
 GlobalParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&GlobalParameters::Nu);
 
 
 oomph_info << "Building/solving for  GlobalParameters::Target_area = "
            << GlobalParameters::Target_area << std::endl;
 
 // Set up the problem with 2D six-node elements 
 StefanBoltzmannProblem
  <ProjectableUnsteadyHeatElement<
   PseudoSolidNodeUpdateElement<TUnsteadyHeatElement<2,3>,
                                TPVDElement<2,3> > > >
  problem;
 
 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);
 problem.assign_initial_values_impulsive();
  
 //Output initial condition
 problem.doc_solution();
  
 // Timestepping loop
 while (problem.time_pt()->time()<=t_max)
  {     
   // Take timestep
   problem.unsteady_newton_solve(dt);
   
   //Output solution
   problem.doc_solution();
  }
    
} //end of main


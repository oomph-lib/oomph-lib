//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Time-harmonic deformation of an elastic annulus subject to 
//  displacement and pressure

//Oomph-lib includes
#include "generic.h"
#include "time_harmonic_linear_elasticity.h"
#include "oomph_crbond_bessel.h"

//The meshes
//#include "annular_meshes.h"
#include "meshes/rectangular_quadmesh.h"
 

using namespace std;

using namespace oomph;

//#define ADAPTIVE
//#undef ADAPTIVE


namespace oomph
{

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//===================================================================
/// 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it.
//==================================================================
template<class ELEMENT>
class TwoDAnnularMesh : public virtual RectangularQuadMesh<ELEMENT>
{

public:


 ///Constructor
 TwoDAnnularMesh(const bool& periodic, const double& azimuthal_fraction,
                 const unsigned &ntheta, const unsigned &nr, 
                 const double& a, const double &h, TimeStepper* time_stepper_pt
                 = &Mesh::Default_TimeStepper) :
  RectangularQuadMesh<ELEMENT>(ntheta,nr,1.0,1.0,
                               periodic,time_stepper_pt)
  {
   // Wrap mesh into annular shape
   wrap_into_annular_shape(a,h,azimuthal_fraction);
  }


  private:
  
  /// Wrap mesh into annular shape
  void wrap_into_annular_shape(const double& a, const double& h,
                               const double& azimuthal_fraction)
  {
   //Create the hole
   Ellipse ellipse(a,a);
   
   //Set all the initial positions of the nodes
   Vector<double> xi(1);
   Vector<double> base(2);
   Vector<double> N(2);
   const unsigned n_node = this->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     // Pointer to node
     Node* nod_pt = this->node_pt(n);

     // Get the angle of the node -- rotate such that jump in angle
     // appears at periodic boundary. Shrink domain slightly
     // to keep angle unique
     xi[0] = (1.0-1.0e-10)*(-azimuthal_fraction*nod_pt->x(0))*
      2.0*MathematicalConstants::Pi+
      MathematicalConstants::Pi-
      (1.0-azimuthal_fraction)*2.0*MathematicalConstants::Pi;
   
     //Get the node's fraction in the radial direction
     double w = nod_pt->x(1);
     
     //Get the position on the ellipse base
     ellipse.position(xi,base);
     
     //Get the unit normal, if it were a circle , by normalising the base
     double norm = sqrt(base[0]*base[0] + base[1]*base[1]);
     N[0] = base[0]/norm; 
     N[1] = base[1]/norm;
     
     //Set circular film from the ellipse
     nod_pt->x(0) = base[0] + w*(h+a-norm)*N[0];
     nod_pt->x(1) = base[1] + w*(h+a-norm)*N[1];

     // Set boundary coordinates
     Vector<double> xi_bound(1);

     // Polar angle for boundary coordinate on boundary 0
     if(nod_pt->is_on_boundary(0))
      {       
       xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
       nod_pt->set_coordinates_on_boundary(0,xi_bound);
      }


     // Radius for boundary coordinate on boundary 1
     if(nod_pt->is_on_boundary(1))
      {      
       xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(1,xi_bound);
      }

     // Polar angle for boundary coordinate on boundary 2
     if(nod_pt->is_on_boundary(2))
      {
       xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
       nod_pt->set_coordinates_on_boundary(2,xi_bound);
      }


     // Radius for boundary coordinate on boundary 3
     if(nod_pt->is_on_boundary(3))
      {      
       xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(3,xi_bound);
      }
    }

   this->Boundary_coordinate_exists[0]=true;
   this->Boundary_coordinate_exists[1]=true;
   this->Boundary_coordinate_exists[2]=true;
   this->Boundary_coordinate_exists[3]=true;
  }
 
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//===================================================================
/// Refineable 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it.
//==================================================================
template<class ELEMENT>
 class RefineableTwoDAnnularMesh : public virtual TwoDAnnularMesh<ELEMENT>,
  public virtual RefineableQuadMesh<ELEMENT>
{

public:


 ///Constructor
 RefineableTwoDAnnularMesh(const bool& periodic, 
                           const double& azimuthal_fraction,
                           const unsigned &ntheta, 
                           const unsigned &nr, 
                           const double& a, 
                           const double &h, 
                           TimeStepper* time_stepper_pt
                           = &Mesh::Default_TimeStepper) :
 RectangularQuadMesh<ELEMENT>(ntheta,nr,1.0,1.0,
                              periodic,time_stepper_pt),
  TwoDAnnularMesh<ELEMENT>(periodic,azimuthal_fraction,ntheta, 
                           nr,a,h,time_stepper_pt)  
   {
    // Nodal positions etc. were created in constructor for
    // RectangularMesh<...>. Only need to setup quadtree forest
    this->setup_quadtree_forest();
    
    //Setup the periodic neighbour information of the TreeRoots
    //Cast to specific elements
    if (periodic)
     {
      Vector<TreeRoot*> left_root_pt(nr);
      Vector<TreeRoot*> right_root_pt(nr);
      for(unsigned i=0;i<nr;i++) 
       {
        left_root_pt[i] = 
         dynamic_cast<ELEMENT*>(this->element_pt(i*ntheta))->
         tree_pt()->root_pt();
        
        right_root_pt[i] = 
         dynamic_cast<ELEMENT*>(this->element_pt((i+1)*ntheta-1))->
         tree_pt()->root_pt();
       }
      
      using namespace QuadTreeNames;
      //Set the neighbour and periodicity
      for(unsigned i=0;i<nr;i++) 
       {
        left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
        left_root_pt[i]->set_neighbour_periodic(W);
        
        right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
        right_root_pt[i]->set_neighbour_periodic(E);
       }
     }
   }
 
};
 


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 // Output directory
 string Directory="RESLT";
 
 // Number of elements in azimuthal direction in solid
 unsigned Ntheta_solid=20;

 // Number of elements in radial direction in solid mesh
 unsigned Nr_solid=10;

 /// Poisson's ratio
 double Nu = 0.3;

 /// Non-dim frequency
 double Omega=10.0;
  
 /// The elasticity tensor
 TimeHarmonicIsotropicElasticityTensor E(Nu);
 
 /// Thickness of annulus
 double H_annulus=0.5;

 /// Displacement amplitude on inner radius
 double Displacement_amplitude=0.1;

 /// Displacement field on inner boundary
 void solid_boundary_displacement(const Vector<double>& x,
                                  Vector<double>& u)
 {
  Vector<double> normal(2);
  double norm=sqrt(x[0]*x[0]+x[1]*x[1]);
  normal[0]=x[0]/norm;
  normal[1]=x[1]/norm;

  u[0]=Displacement_amplitude*normal[0];
  u[1]=Displacement_amplitude*normal[1];
 }

 /// Uniform pressure
 double P = 0.0;

 /// Constant pressure load
 void constant_pressure(const Vector<double> &x,
                        const Vector<double> &n, 
                        Vector<std::complex<double> >&traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    //traction[i].real() = -P*n[i];
    //traction[i].imag() =  P*n[i];
    traction[i] = complex<double>(-P*n[i],P*n[i]);
   }
 } // end of pressure load
 

 /// Helper function to evaluate Y_n(x) from bloody maple output
 double BesselY(const double& n, const double& x)
 {
  // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
  double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
  CRBond_Bessel::bessjy01a(x,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
  
  if (n==0.0)
   {
    return y0;
   }
  else if (n==1.0)
   {
    return y1;
   }
  else
   {
    cout << "Never get here...";
    abort();
   }
   
 }


 /// Helper function to evaluate J_n(x) from bloody maple output
 double BesselJ(const double& n, const double& x)
 {
  
  // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
  double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
  CRBond_Bessel::bessjy01a(x,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
  
  if (n==0.0)
   {
    return j0;
   }
  else if (n==1.0)
   {
    return j1;
   }
  else
   {
    cout << "Never get here...";
    abort();
   }
 }

 /// Exact solution as a Vector
 void exact_u(const Vector<double>& x, Vector<double>& u)
 {
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);

  double MapleGenVar1;
  double MapleGenVar2;
  double MapleGenVar3;
  double MapleGenVar4;
  double MapleGenVar5;
  double MapleGenVar6;
  double MapleGenVar7;
  double t0;

      MapleGenVar3 = Displacement_amplitude*Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/
(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))+Displacement_amplitude
*Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(
1.0+H_annulus))*H_annulus-Displacement_amplitude*BesselY(0.0,Omega/sqrt(Nu/(1.0
+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega-
Displacement_amplitude*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0
+2.0*Nu))*(1.0+H_annulus))*Omega*H_annulus+Displacement_amplitude*sqrt(-(1.0-Nu
)/(-1.0+Nu+2.0*Nu*Nu))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0
+2.0*Nu))*(1.0+H_annulus))-2.0*Displacement_amplitude*sqrt(-(1.0-Nu)/(-1.0+Nu+
2.0*Nu*Nu))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(
1.0+H_annulus))*Nu;
      MapleGenVar2 = MapleGenVar3-BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu)))*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))-BesselY(1.0,Omega/
sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*
Nu*Nu))*H_annulus+BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*
Nu)))*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu+BesselY(1.0,Omega/sqrt(Nu/(1.0+
Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu*
H_annulus+2.0*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))
)*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu*Nu+2.0*BesselY(1.0,Omega/sqrt(Nu/(
1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*
Nu*Nu*H_annulus;
      MapleGenVar7 = Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+
2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu)))+Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)
+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*H_annulus*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu
)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)
/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar6 = MapleGenVar7-BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*H_annulus*BesselJ(1.0,Omega/sqrt(
Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))
*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+
H_annulus))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))
-2.0*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(
1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Nu*BesselJ(1.0,Omega/sqrt(Nu/(
1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar7 = MapleGenVar6+BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)
/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu
)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*H_annulus*BesselY(1.0,Omega/sqrt(Nu/
(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar5 = MapleGenVar7-sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(
1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*
BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-Nu*Omega*
BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+
H_annulus))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+
2.0*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0
-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Nu*BesselY(1.0,Omega/sqrt(Nu/(1.0+
Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-Nu*Omega*BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)
/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*H_annulus*BesselY(1.0,Omega/
sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar4 = 1/MapleGenVar5;
      MapleGenVar5 = BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+
2.0*Nu))*r);
      MapleGenVar3 = MapleGenVar4*MapleGenVar5;
      MapleGenVar1 = MapleGenVar2*MapleGenVar3;
      MapleGenVar6 = -Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+
2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu)))-Nu*Omega*BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)
+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*H_annulus*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu
)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)
/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar5 = MapleGenVar6+BesselY(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*H_annulus*BesselJ(1.0,Omega/sqrt(
Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))
*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+
H_annulus))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+
2.0*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0
-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Nu*BesselJ(1.0,Omega/sqrt(Nu/(1.0+
Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar6 = MapleGenVar5-BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)
/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu
)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*H_annulus*BesselY(1.0,Omega/sqrt(Nu/
(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar4 = MapleGenVar6+sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(
1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*
BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+Nu*Omega*
BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+
H_annulus))*BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))
-2.0*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(
1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Nu*BesselY(1.0,Omega/sqrt(Nu/(
1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+Nu*Omega*BesselJ(0.0,Omega/sqrt(Nu/(1.0
+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*H_annulus*BesselY(1.0,
Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)));
      MapleGenVar3 = 1/MapleGenVar4;
      MapleGenVar6 = -P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(1.0,Omega/
sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))-P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*
Nu*Nu))*H_annulus*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*
Nu)))+P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu*BesselJ(1.0,Omega/sqrt(Nu/(1.0+
Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu*
H_annulus*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+
2.0*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*Nu*Nu*BesselJ(1.0,Omega/sqrt(Nu/(1.0+
Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)))+2.0*P*sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*
Nu*Nu*H_annulus*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu
)));
      MapleGenVar7 = MapleGenVar6-BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*Displacement_amplitude-BesselJ(0.0
,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Omega*
H_annulus*Displacement_amplitude;
      MapleGenVar5 = MapleGenVar7+sqrt(-(1.0-Nu)/(-1.0+Nu+2.0*Nu*Nu))*BesselJ(
1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*
Displacement_amplitude+Nu*Omega*BesselJ(0.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)
+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*Displacement_amplitude-2.0*sqrt(-(1.0-Nu)/(
-1.0+Nu+2.0*Nu*Nu))*BesselJ(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+
2.0*Nu))*(1.0+H_annulus))*Nu*Displacement_amplitude+Nu*Omega*BesselJ(0.0,Omega/
sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu))*(1.0+H_annulus))*H_annulus*
Displacement_amplitude;
      MapleGenVar6 = BesselY(1.0,Omega/sqrt(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+
2.0*Nu))*r);
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      t0 = MapleGenVar1+MapleGenVar2;

  
  u[0]=t0*x[0]/r;
  u[1]=t0*x[1]/r;
  
  u[2]=-t0*x[0]/r;
  u[3]=-t0*x[1]/r;

 }

} //end namespace



//=============begin_problem============================================ 
/// Annular disk
//====================================================================== 
template<class ELASTICITY_ELEMENT>
class AnnularDiskProblem : public Problem
{

public:

 /// Constructor:
 AnnularDiskProblem();
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

#ifdef ADAPTIVE

 /// Access function for the solid mesh
 TreeBasedRefineableMeshBase*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#else

 /// Access function for the solid mesh
 Mesh*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#endif

 /// Actions before adapt: Wipe the mesh of traction elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of traction elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

private:

 /// \short Create traction elements
 void create_traction_elements();

 /// Delete traction elements
 void delete_traction_elements();

#ifdef ADAPTIVE

 /// Pointer to solid mesh
 TreeBasedRefineableMeshBase* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 Mesh* Solid_mesh_pt;

#endif

 /// Pointer to mesh of traction elements
 Vector<Mesh*> Traction_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
 /// Do we have a hole in the elastic coating?
 bool Have_hole;

};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELASTICITY_ELEMENT>
AnnularDiskProblem<ELASTICITY_ELEMENT>::AnnularDiskProblem() 
{

 // Do we have a hole in the coating?
 if (CommandLineArgs::command_line_flag_has_been_set("--have_hole"))
  {
   Have_hole=true;
  }
 else
  {
   Have_hole=false;
  }

 // The coating mesh is periodic if there's no hole
 bool periodic=true;
 if (Have_hole) periodic=false;

 // Azimuthal fraction of elastic coating
 double azimuthal_fraction_of_coating=0.9;
 if (periodic)
  {
   azimuthal_fraction_of_coating=1.0;
  }
 
 // Solid mesh
 //-----------

 // Innermost radius for solid mesh
 double a=1.0;
 
#ifdef ADAPTIVE

 // Build solid mesh
 solid_mesh_pt() = new 
  RefineableTwoDAnnularMesh<ELASTICITY_ELEMENT>(
   periodic,azimuthal_fraction_of_coating,
   Global_Parameters::Ntheta_solid,
   Global_Parameters::Nr_solid,a,
   Global_Parameters::H_annulus);
 
 // Set error estimators
 solid_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Set some refinement targets
 solid_mesh_pt()->max_refinement_level()=2;
 
#else
 
 // Build solid mesh
 solid_mesh_pt() = new 
  TwoDAnnularMesh<ELASTICITY_ELEMENT>(
   periodic,azimuthal_fraction_of_coating,
   Global_Parameters::Ntheta_solid,
   Global_Parameters::Nr_solid,a,
   Global_Parameters::H_annulus);
 
#endif

 // Let's have a look where the boundaries are
 solid_mesh_pt()->output("solid_mesh.dat");
 solid_mesh_pt()->output_boundaries("solid_mesh_boundary.dat");

 //Assign the physical properties to the elements
 //Loop over the elements in the main mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELASTICITY_ELEMENT *el_pt = 
    dynamic_cast<ELASTICITY_ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->elasticity_tensor_pt() = &Global_Parameters::E;

   // Non-dim frequency
   el_pt->omega_pt()= &Global_Parameters::Omega;
  }
 
 // Number of meshes on the FSI boundary
 unsigned n_bound_mesh=1;
 if (Have_hole)
  {
   n_bound_mesh=3;
  }
 
 // Construct the traction element mesh
 Traction_mesh_pt.resize(n_bound_mesh);
 for (unsigned i=0;i<n_bound_mesh;i++)
  {
   Traction_mesh_pt[i]=new Mesh;
  }
 create_traction_elements(); 
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Add traction sub-meshes
 unsigned n_tract=Traction_mesh_pt.size();
 for (unsigned i=0;i<n_tract;i++)
  {
   add_sub_mesh(Traction_mesh_pt[i]);
  }

 // Build combined "global" mesh
 build_global_mesh();
 

 // Solid boundary conditions:
 //---------------------------
 // Pin the bottom boundary (boundary 0) in both directions
 unsigned n_node = solid_mesh_pt()->nboundary_node(0);
 
 //Loop over the nodes to pin and assign boundary displacements on 
 //solid boundary
 Vector<double> u(2);
 Vector<double> x(2);
 for(unsigned i=0;i<n_node;i++)
  {
   Node* nod_pt=solid_mesh_pt()->boundary_node_pt(0,i);
   nod_pt->pin(0);
   nod_pt->pin(1);
   nod_pt->pin(2);
   nod_pt->pin(3);
   
   // Assign displacements
   x[0]=nod_pt->x(0);   
   x[1]=nod_pt->x(1);   
   Global_Parameters::solid_boundary_displacement(x,u);
   nod_pt->set_value(0,u[0]);
   nod_pt->set_value(1,u[1]);
   nod_pt->set_value(2,-u[0]);
   nod_pt->set_value(3,-u[1]);
  }

 //Assign equation numbers
 cout << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory(Global_Parameters::Directory);

} //end of constructor




//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of traction elements
//========================================================================
template<class ELASTICITY_ELEMENT>
void AnnularDiskProblem<ELASTICITY_ELEMENT>::actions_before_adapt()
{
 // Kill the traction elements and wipe surface mesh
 delete_traction_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the mesh of traction elements
//========================================================================
template<class ELASTICITY_ELEMENT>
void AnnularDiskProblem<ELASTICITY_ELEMENT>::actions_after_adapt()
{
 // Create traction elements from all elements that are 
 // adjacent to FSI boundaries and add them to surface meshes
 create_traction_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
}// end of actions_after_adapt


//============start_of_create_traction_elements==============================
/// Create traction elements 
//=======================================================================
template<class ELASTICITY_ELEMENT>
void AnnularDiskProblem<ELASTICITY_ELEMENT>::create_traction_elements()
{

 // Load outer surface (2) and both "ends" (1 and 3) if there's a hole
 unsigned b_lo=1;
 unsigned b_hi=3;
 // ...otherwise load only the outside (2)
 if (!Have_hole)
  {     
   b_lo=2;
   b_hi=2;
  }

 unsigned count=0;
 for (unsigned b=b_lo;b<=b_hi;b++)
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned n_element = solid_mesh_pt()->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELASTICITY_ELEMENT* bulk_elem_pt = dynamic_cast<ELASTICITY_ELEMENT*>(
      solid_mesh_pt()->boundary_element_pt(b,e));
     
     //Find the index of the face of element e along boundary b
     int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
     
     // Create element
     TimeHarmonicLinearElasticityTractionElement<ELASTICITY_ELEMENT>* el_pt=
      new TimeHarmonicLinearElasticityTractionElement<ELASTICITY_ELEMENT>
      (bulk_elem_pt,face_index);   

     // Add to mesh
     Traction_mesh_pt[count]->add_element_pt(el_pt);
     
     // Associate element with bulk boundary (to allow it to access
     // the boundary coordinates in the bulk mesh)
     el_pt->set_boundary_number_in_bulk_mesh(b); 
     
     //Set the traction function
     el_pt->traction_fct_pt() = Global_Parameters::constant_pressure;

    }
   count++;
  }
 
} // end of create_traction_elements




//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the  traction meshes
//=======================================================================
template<class ELASTICITY_ELEMENT>
void AnnularDiskProblem<ELASTICITY_ELEMENT>::delete_traction_elements()
{

 // Loop over traction meshes
 unsigned n=Traction_mesh_pt.size();
 for (unsigned m=0;m<n;m++)
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Traction_mesh_pt[m]->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Traction_mesh_pt[m]->element_pt(e);
    }
   
   // Wipe the mesh
   Traction_mesh_pt[m]->flush_element_and_node_storage();
  }

} // end of delete_traction_elements





//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELASTICITY_ELEMENT>
void AnnularDiskProblem<ELASTICITY_ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot=5; 

 // Output displacement field
 //--------------------------
 sprintf(filename,"%s/elast_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();


 // Output traction elements
 //-------------------------
 sprintf(filename,"%s/traction_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n_tract=Traction_mesh_pt.size();
 for (unsigned i=0;i<n_tract;i++)
  {
   Traction_mesh_pt[i]->output(some_file,n_plot);
  }
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output_fct(some_file,n_plot,Global_Parameters::exact_u); 
 some_file.close();

 // Increment label for output files
 Doc_info.number()++;
 
} //end doc



//=======start_of_main==================================================
/// Driver for annular disk loaded by pressure
//======================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &Global_Parameters::Directory);
  
 // Do have a hole in the annulus?
 CommandLineArgs::specify_command_line_flag("--have_hole");
  
 // Number of elements in azimuthal direction
 CommandLineArgs::specify_command_line_flag(
  "--ntheta_solid",
  &Global_Parameters::Ntheta_solid);

 // Number of elements in radial direction
 CommandLineArgs::specify_command_line_flag(
  "--nr_solid",
  &Global_Parameters::Nr_solid);
 
#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=3;
 CommandLineArgs::specify_command_line_flag("--max_adapt",&max_adapt);

 // Number of uniform refinements
 unsigned nrefine=0; 
 CommandLineArgs::specify_command_line_flag("--nrefine",&nrefine);

#endif
 
 // P increment
 double p_increment=0.1;
 CommandLineArgs::specify_command_line_flag("--p_increment",&p_increment);
 
 // Number of steps
 unsigned nstep=3; 
 CommandLineArgs::specify_command_line_flag("--nstep",&nstep);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

#ifdef ADAPTIVE

 //Set up the problem
 AnnularDiskProblem<RefineableQTimeHarmonicLinearElasticityElement<2,3> >
  problem;

 // Refine unformly
 for (unsigned i=0;i<nrefine;i++)
  {
   problem.refine_uniformly();
  }

#else

 //Set up the problem
 AnnularDiskProblem<QTimeHarmonicLinearElasticityElement<2,3> > problem;

#endif

 
 // Doc initial configuration
 problem.doc_solution();
 
 // Initial values for parameter values
 Global_Parameters::P=0.0; 
 
 //Parameter incrementation
 for(unsigned i=0;i<nstep;i++)
  {

#ifdef ADAPTIVE

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   problem.newton_solve(max_adapt);

#else

   // Solve the problem with Newton's method
   problem.newton_solve();

#endif

   
   // Doc solution
   problem.doc_solution();
   
   // Increment pressure
   Global_Parameters::P+=p_increment;
  }
 
} //end of main









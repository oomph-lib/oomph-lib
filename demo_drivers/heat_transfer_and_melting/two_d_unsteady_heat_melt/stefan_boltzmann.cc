//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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

// The mesh
#include "meshes/triangle_mesh.h"

using namespace oomph;
using namespace std;




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

 /// Radius of inner circle
 double Radius_innermost=0.5;

 /// Inner radius of annular region
 double Radius_inner=1.0;

 /// Outer radius of annular region
 double Radius_outer=1.5;

 /// Temperature on boundary of inner circle
 double U0 = 0.8288627710;

 /// Strength of source function in inner region
 double S0=0.1;

 /// Strength of source function in outer region
 double S1=1.0;

 /// Source function
 void get_source(const double& time, const Vector<double>& x, double& source)
 {
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  if (r>0.5*(Radius_inner+Radius_innermost))
   {
    source = S1;
   }
  else
   {
    source = S0;
   }
   
 }

 /// Exact solution as a Vector
 void get_exact_u(const double& time, const Vector<double>& x, 
                  Vector<double>& u)
 {
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);

  // Solutions from maple
  double t0=0.0;
  if (r<0.5*(Radius_inner+Radius_innermost))
   {
    t0 = 0.25*S0/Beta0*(r*r-1.0*Radius_innermost*Radius_innermost)+U0;
   }
  else
   {
    
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
    double MapleGenVar15 = 0.0;
    double MapleGenVar16 = 0.0;
    double MapleGenVar14 = 0.0;
    double MapleGenVar12 = 0.0;
    double MapleGenVar10 = 0.0;
    double MapleGenVar2  = 0.0;
    
      MapleGenVar1 = S1/Beta1*r*r/4.0;
      MapleGenVar4 = 1/Beta1/4.0;
      MapleGenVar6 = 1/(log(Radius_inner)-log(Radius_outer));
      MapleGenVar9 = -2.0*log(Radius_inner)*S0*Radius_innermost*
Radius_innermost+2.0*log(Radius_outer)*S0*Radius_innermost*Radius_innermost+4.0
*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*Theta_0*Radius_inner-8.0*log(
Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*Theta_0*Radius_innermost+4.0*log(
Radius_inner)*Sigma*U0*U0*U0*U0*Radius_inner-8.0*log(Radius_inner)*Sigma*U0*U0*
U0*U0*Radius_innermost-4.0*log(Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*
Theta_0*Radius_inner+8.0*log(Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*
Theta_0*Radius_innermost-4.0*log(Radius_outer)*Sigma*U0*U0*U0*U0*Radius_inner+
8.0*log(Radius_outer)*Sigma*U0*U0*U0*U0*Radius_innermost+2.0*log(Radius_inner)*
S0*Radius_inner*Radius_innermost-2.0*log(Radius_outer)*S0*Radius_inner*
Radius_innermost-2.0*S1*Radius_inner*Radius_inner*log(Radius_inner);
      MapleGenVar8 = MapleGenVar9+2.0*log(Radius_outer)*S1*Radius_inner*
Radius_inner+16.0*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*U0*
Radius_inner-32.0*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*U0*
Radius_innermost+24.0*log(Radius_inner)*Sigma*Theta_0*Theta_0*U0*U0*
Radius_inner-48.0*log(Radius_inner)*Sigma*Theta_0*Theta_0*U0*U0*
Radius_innermost+16.0*log(Radius_inner)*Sigma*Theta_0*U0*U0*U0*Radius_inner
-32.0*log(Radius_inner)*Sigma*Theta_0*U0*U0*U0*Radius_innermost-16.0*log(
Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_inner+32.0*log(
Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_innermost-24.0*log(
Radius_outer)*Sigma*Theta_0*Theta_0*U0*U0*Radius_inner+48.0*log(Radius_outer)*
Sigma*Theta_0*Theta_0*U0*U0*Radius_innermost-16.0*log(Radius_outer)*Sigma*
Theta_0*U0*U0*U0*Radius_inner+32.0*log(Radius_outer)*Sigma*Theta_0*U0*U0*U0*
Radius_innermost;
      MapleGenVar9 = log(r);
      MapleGenVar7 = MapleGenVar8*MapleGenVar9;
      MapleGenVar5 = MapleGenVar6*MapleGenVar7;
      MapleGenVar3 = MapleGenVar4*MapleGenVar5;
      MapleGenVar6 = -Beta1*pow(2.0,0.75)*pow((2.0*Sigma*Theta_0*Theta_0*
Theta_0*Theta_0+8.0*Sigma*Theta_0*Theta_0*Theta_0*U0+12.0*Sigma*Theta_0*Theta_0
*U0*U0+8.0*Sigma*Theta_0*U0*U0*U0+2.0*Sigma*U0*U0*U0*U0+S0*Radius_innermost)/
Sigma,0.25)*log(Radius_outer)/2.0-log(Radius_inner)*S1*Radius_outer*
Radius_outer/4.0;
      MapleGenVar7 = MapleGenVar6+log(Radius_outer)*S1*Radius_inner*
Radius_inner/4.0;
      MapleGenVar8 = MapleGenVar7;
      MapleGenVar11 = 1.0/4.0;
      MapleGenVar13 = log(Radius_inner);
      MapleGenVar15 = 2.0*log(Radius_inner)*S0*Radius_innermost*
Radius_innermost-2.0*log(Radius_outer)*S0*Radius_innermost*Radius_innermost-4.0
*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*Theta_0*Radius_inner+8.0*log(
Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*Theta_0*Radius_innermost-4.0*log(
Radius_inner)*Sigma*U0*U0*U0*U0*Radius_inner+8.0*log(Radius_inner)*Sigma*U0*U0*
U0*U0*Radius_innermost+4.0*log(Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*
Theta_0*Radius_inner-8.0*log(Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*
Theta_0*Radius_innermost+4.0*log(Radius_outer)*Sigma*U0*U0*U0*U0*Radius_inner
-8.0*log(Radius_outer)*Sigma*U0*U0*U0*U0*Radius_innermost-2.0*log(Radius_inner)
*S0*Radius_inner*Radius_innermost+2.0*log(Radius_outer)*S0*Radius_inner*
Radius_innermost+2.0*S1*Radius_inner*Radius_inner*log(Radius_inner)-2.0*log(
Radius_outer)*S1*Radius_inner*Radius_inner-S1*Radius_inner*Radius_inner;
      MapleGenVar16 = MapleGenVar15+2.0*pow(2.0,0.75)*pow((2.0*Sigma*Theta_0*
Theta_0*Theta_0*Theta_0+8.0*Sigma*Theta_0*Theta_0*Theta_0*U0+12.0*Sigma*Theta_0
*Theta_0*U0*U0+8.0*Sigma*Theta_0*U0*U0*U0+2.0*Sigma*U0*U0*U0*U0+S0*
Radius_innermost)/Sigma,0.25)*Beta1+S1*Radius_outer*Radius_outer-4.0*Theta_0*
Beta1-16.0*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_inner+32.0
*log(Radius_inner)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_innermost-24.0*log(
Radius_inner)*Sigma*Theta_0*Theta_0*U0*U0*Radius_inner+48.0*log(Radius_inner)*
Sigma*Theta_0*Theta_0*U0*U0*Radius_innermost;
      MapleGenVar14 = MapleGenVar16-16.0*log(Radius_inner)*Sigma*Theta_0*U0*U0*
U0*Radius_inner+32.0*log(Radius_inner)*Sigma*Theta_0*U0*U0*U0*Radius_innermost+
16.0*log(Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_inner-32.0*log(
Radius_outer)*Sigma*Theta_0*Theta_0*Theta_0*U0*Radius_innermost+24.0*log(
Radius_outer)*Sigma*Theta_0*Theta_0*U0*U0*Radius_inner-48.0*log(Radius_outer)*
Sigma*Theta_0*Theta_0*U0*U0*Radius_innermost+16.0*log(Radius_outer)*Sigma*
Theta_0*U0*U0*U0*Radius_inner-32.0*log(Radius_outer)*Sigma*Theta_0*U0*U0*U0*
Radius_innermost;
      MapleGenVar12 = MapleGenVar13*MapleGenVar14;
      MapleGenVar10 = MapleGenVar11*MapleGenVar12;
      MapleGenVar11 = Beta1*log(Radius_outer)*Theta_0;
      MapleGenVar9 = MapleGenVar10+MapleGenVar11;
      MapleGenVar5 = MapleGenVar8+MapleGenVar9;
      MapleGenVar6 = 1/Beta1/(log(Radius_inner)-log(Radius_outer));
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3+MapleGenVar4;
      t0 = MapleGenVar1+MapleGenVar2;
   }

  u[0] = t0;

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

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){} 

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to mesh of radiative flux elements
 Mesh* Unsteady_heat_flux_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor
//========================================================================
template<class ELEMENT>
StefanBoltzmannProblem<ELEMENT>::StefanBoltzmannProblem()
{ 

 // Open trace file
 Trace_file.open((GlobalParameters::Directory+std::string("/trace.dat")).c_str());
 
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
                            n_inner, // hierher
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
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(triangle_mesh_parameters,
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
   el_pt->source_fct_pt() = &GlobalParameters::get_source;
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
   el_pt->source_fct_pt() = &GlobalParameters::get_source;
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

 // Create mesh with heat flux elements
 //------------------------------------
 Unsteady_heat_flux_mesh_pt=new Mesh;

 // Face elements on potentially sun-exposed boundaries
 Vector<FiniteElement*> shielding_face_element_pt;

 // Boundaries that participate in radiative heat exchange
 std::set<unsigned> bc_set;
 bc_set.insert(2);
 bc_set.insert(3);
 bc_set.insert(4);
 bc_set.insert(5);   
 for (std::set<unsigned>::iterator it=bc_set.begin();it!=bc_set.end();it++)
  {
   unsigned b=(*it);
   
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
       el_pt->set_integration_scheme(new MyIntegral(GlobalParameters::Nintpt));
      }

     // Set non-dim Stefan Boltzmann constant
     el_pt->sigma_pt()= &GlobalParameters::Sigma;
     
     // Set zero-centrigrade offset in Stefan Boltzmann law
     el_pt->theta_0_pt()= &GlobalParameters::Theta_0;
     
     // Add to mesh
     Unsteady_heat_flux_mesh_pt->add_element_pt(el_pt);
     
     // Add to shielding surface
     shielding_face_element_pt.push_back(el_pt);     
    }  
  }


 StefanBoltzmannHelper::Nsample=2;
 StefanBoltzmannHelper::Nx_bin=20;
 StefanBoltzmannHelper::Ny_bin=20;

 //Setup mutual visibility for all Stefan Boltzmann elements
 StefanBoltzmannHelper::setup_stefan_boltzmann_visibility(
  shielding_face_element_pt);

 // Doc the populated bins
 ofstream some_file;
 some_file.open((GlobalParameters::Directory+
                 std::string("/populated_bins.dat")).c_str());
 StefanBoltzmannHelper::doc_bins(some_file);
 some_file.close();

 // Doc the sample points used to assess visibility 
 some_file.open((GlobalParameters::Directory+
                 std::string("/sample_points.dat")).c_str());
 StefanBoltzmannHelper::doc_sample_points(some_file,shielding_face_element_pt);
 some_file.close();
 
 // Build the entire mesh from its submeshes
 add_sub_mesh(Bulk_mesh_pt);
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
 
 
 // Pin boundary values
 for(unsigned b=0;b<2;b++)
  {
   unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
     nod_pt->pin(0);       
     // oomph_info << "Pinning at r="
     //            << sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2))
     //            << std::endl;
    }
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end of constructor



//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void StefanBoltzmannProblem<ELEMENT>::doc_solution(DocInfo& doc_info) 
{ 

 ofstream some_file,some_file2;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output solution 
 //-----------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 unsigned npts_coarse=2;
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Output Stefan Boltzmann radiation (best suited for tecplot)
 //------------------------------------------------------------
 sprintf(filename,"%s/stefan_boltzmann_radiation%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned nel=Unsteady_heat_flux_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>
  (Unsteady_heat_flux_mesh_pt->element_pt(e))->
   output_stefan_boltzmann_radiation_rays(some_file);
  }
 some_file.close();


// Output Stefan Boltzmann radiation rays (best suited for paraview)
//------------------------------------------------------------------
bool output_rays=true;
if (output_rays)
 {
  unsigned count=0;
  unsigned nel=Unsteady_heat_flux_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    sprintf(filename,"%s/stefan_boltzmann_rays%i_intpt%i.dat",
            doc_info.directory().c_str(),
            doc_info.number(),count);
    some_file.open(filename);
    StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>* el_pt=
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>
     (Unsteady_heat_flux_mesh_pt->element_pt(e));
    unsigned nintpt=el_pt->integral_pt()->nweight();
    for (unsigned ipt=0;ipt<nintpt;ipt++)
     {
      el_pt->output_stefan_boltzmann_radiation_rays(some_file,ipt);
     }
    some_file.close();
    count++;
   }
 }


// Output heat flux etc.
//----------------------
sprintf(filename,"%s/flux%i.dat",
        doc_info.directory().c_str(),
        doc_info.number());
some_file.open(filename);
Unsteady_heat_flux_mesh_pt->output(some_file,5);
some_file.close();




 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                          GlobalParameters::get_exact_u); 
 some_file.close();
 
 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
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


 // Get norm of solution
 Bulk_mesh_pt->compute_norm(norm); 
 Trace_file  << norm << std::endl;
 
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
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
 
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory(GlobalParameters::Directory);


 // hierher refinement loop
 for (unsigned i=0;i<1;i++)
  {

   oomph_info << "Building/solving for  GlobalParameters::Target_area = "
              << GlobalParameters::Target_area << std::endl;
   
   // Set up the problem with 2D six-node elements 
   StefanBoltzmannProblem<TUnsteadyHeatElement<2,3> >  problem;

   // Initialise timestep -- also sets the weights for all timesteppers
   // in the problem.
   double dt=0.01;
   problem.initialise_dt(dt);
   problem.assign_initial_values_impulsive();
   
   // Steady solve only
   {
    // Do steady solve
    problem.steady_newton_solve();
    
    //Output solution
    problem.doc_solution(doc_info);
    
    exit(0);
   }

   //Output initial condition
   problem.doc_solution(doc_info);

   //Increment counter for solutions 
   doc_info.number()++;
   
   // Timestepping loop
   unsigned nstep = 20;
   for (unsigned istep=0;istep<nstep;istep++)
    {
     cout << " Timestep " << istep << std::endl;
     
     // Take timestep
     problem.unsteady_newton_solve(dt);
     
     //Output solution
     problem.doc_solution(doc_info);
     
     //Increment counter for solutions 
     doc_info.number()++;
    }

   // Adjust
   GlobalParameters::Target_area/=4.0;
  }
 
    
} //end of main


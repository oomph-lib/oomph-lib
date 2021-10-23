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
// Driver for Helmholtz/TimeHarmonicTimeHarmonicLinElast coupling

//Oomph-lib includes
#include "generic.h"
#include "helmholtz.h"
#include "time_harmonic_linear_elasticity.h"
#include "multi_physics.h"

//The mesh
#include "meshes/annular_mesh.h"

using namespace std;
using namespace oomph;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Square of wavenumber for the Helmholtz equation
 double K_squared=10.0;
 
 /// Radius of outer boundary of Helmholtz domain
 double Outer_radius=4.0; 

 /// FSI parameter
 double Q=0.0;

 /// Non-dim thickness of elastic coating
 double H_coating=0.3; 
 
 /// Poisson's ratio
 double Nu = 0.3;

 /// The elasticity tensor for the solid
 TimeHarmonicIsotropicElasticityTensor E(Nu);
  
 /// Density ratio: solid to fluid
 double Density_ratio=0.0;

 /// Non-dim square of frequency for solid -- dependent variable!
 double Omega_sq=0.0;

 /// Function to update dependent parameter values
 void update_parameter_values()
 {
  Omega_sq=Density_ratio*Q;
 }

 /// Azimuthal wavenumber for imposed displacement of coating
 /// on inner boundary 
 unsigned N=0; 

 /// Displacement field on inner boundary of solid
 void solid_boundary_displacement(const Vector<double>& x,
                                  Vector<std::complex<double> >& u)
 {
  Vector<double> normal(2);
  double norm=sqrt(x[0]*x[0]+x[1]*x[1]);
  double phi=atan2(x[1],x[0]);
  normal[0]=x[0]/norm;
  normal[1]=x[1]/norm;

  u[0]=complex<double>(normal[0]*cos(double(N)*phi),0.0);
  u[1]=complex<double>(normal[1]*cos(double(N)*phi),0.0);
 }


 /// Output directory
 string Directory="RESLT";
 
 /// Multiplier for number of elements
 unsigned El_multiplier=1;

 /// Interface to Hankel function in maple style
 std::complex<double> HankelH1(const double& k, const double& x)
 {
  Vector<std::complex<double> > h(2);
  Vector<std::complex<double> > hp(2);
  Hankel_functions_for_helmholtz_problem::Hankel_first(2,x,h,hp);
  
  if (k==0.0)
   {
    return h[0];
   }
  else if (k==1.0)
   {
    return h[1];
   }
  else
   {
    cout << "Never get here. k=" << k << std::endl;
    assert(false);
    // Dummy return
    return std::complex<double>(1.0,1.0); 
   }
 }
 

 /// Coefficient in front of Hankel function for axisymmetric solution 
 /// of Helmholtz potential
 std::complex<double> axisym_coefficient()
 {
  std::complex<double> MapleGenVar1;
  std::complex<double> MapleGenVar2;
  std::complex<double> MapleGenVar3;
  std::complex<double> MapleGenVar4;
  std::complex<double> MapleGenVar5;
  std::complex<double> MapleGenVar6;
  std::complex<double> MapleGenVar7;
  std::complex<double> MapleGenVar8;
  std::complex<double> t0;
  

      MapleGenVar3 = (-2.0/(2.0+2.0*Nu)*1.0+2.0/(2.0+2.0*Nu)*1.0*
H_coating)/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*
H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0;
      MapleGenVar5 = -1/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0
*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)*Q/2.0;
      MapleGenVar6 = HankelH1(0.0,sqrt(K_squared))*(-(-2.0/(2.0+2.0*Nu)*1.0
+2.0/(2.0+2.0*Nu)*1.0*H_coating)/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)
-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0+(2.0/(2.0+
2.0*Nu)*1.0-2.0/(2.0+2.0*Nu)*1.0*H_coating+2.0*1.0*Nu/(1.0+Nu)/(1.0
-2.0*Nu)-2.0*1.0*H_coating*Nu/(1.0+Nu)/(1.0-2.0*Nu))/(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*
H_coating)/2.0)/(HankelH1(1.0,sqrt(K_squared))*sqrt(K_squared)-Q*HankelH1(0.0,
sqrt(K_squared))/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*
H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0+Q*HankelH1(0.0,sqrt(K_squared
))*(1.0-2.0*H_coating+H_coating*H_coating)/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+
2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0);
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3+MapleGenVar4;
      MapleGenVar3 = MapleGenVar2;
      MapleGenVar5 = -(2.0/(2.0+2.0*Nu)*1.0-2.0/(2.0+2.0*Nu)*1.0*
H_coating+2.0*1.0*Nu/(1.0+Nu)/(1.0-2.0*Nu)-2.0*1.0*H_coating*Nu/(1.0+Nu
)/(1.0-2.0*Nu))/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*
H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0;
      MapleGenVar7 = (1.0-2.0*H_coating+H_coating*H_coating)/(Nu/(1.0+Nu)/(1.0
-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*
H_coating)/2.0;
      MapleGenVar8 = Q*HankelH1(0.0,sqrt(K_squared))*(-(-2.0/(2.0+2.0*Nu)*
1.0+2.0/(2.0+2.0*Nu)*1.0*H_coating)/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+
2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0+(2.0
/(2.0+2.0*Nu)*1.0-2.0/(2.0+2.0*Nu)*1.0*H_coating+2.0*1.0*Nu/(1.0+Nu
)/(1.0-2.0*Nu)-2.0*1.0*H_coating*Nu/(1.0+Nu)/(1.0-2.0*Nu))/(Nu/(1.0+Nu)/(
1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*
H_coating*H_coating)/2.0)/(HankelH1(1.0,sqrt(K_squared))*sqrt(K_squared)-Q*
HankelH1(0.0,sqrt(K_squared))/(Nu/(1.0+Nu)/(1.0-2.0*Nu)+2.0/(2.0+2.0*Nu)-2.0/(
2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*H_coating)/2.0+Q*HankelH1(0.0,
sqrt(K_squared))*(1.0-2.0*H_coating+H_coating*H_coating)/(Nu/(1.0+Nu)/(1.0-2.0*
Nu)+2.0/(2.0+2.0*Nu)-2.0/(2.0+2.0*Nu)*H_coating+1/(2.0+2.0*Nu)*H_coating*
H_coating)/2.0);
      MapleGenVar6 = MapleGenVar7*MapleGenVar8;
      MapleGenVar4 = MapleGenVar5+MapleGenVar6;
      MapleGenVar1 = MapleGenVar3+MapleGenVar4;
      MapleGenVar2 = 1.0/HankelH1(1.0,sqrt(K_squared))/sqrt(K_squared);
      t0 = MapleGenVar1*MapleGenVar2;

      return t0;
 }


 /// Exact solution for Helmholtz potential for axisymmetric solution
 void exact_axisym_potential(const Vector<double>& x, 
                             Vector<double>& soln)
 {
  std::complex<double> C=axisym_coefficient();
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  soln[0]=real(HankelH1(0,sqrt(K_squared)*r)*C);
  soln[1]=imag(HankelH1(0,sqrt(K_squared)*r)*C);
 }
 
 /// Exact radiated power for axisymmetric solution
 double exact_axisym_radiated_power()
 {
  // Solution is independent of where it's evaluated (have checked!)
  double r=1.0;

  // Get the coefficient for the axisymmetric potential
  std::complex<double> C=axisym_coefficient();
  
  // Argument for Bessel/Hankel functions
  double rr=sqrt(K_squared)*r;  
  
  // Evaluate Hankel functions -- only need zero-th term
  Vector<std::complex<double> > h(2),hp(2);
  Hankel_functions_for_helmholtz_problem::Hankel_first(1,rr,h,hp);
  
  // Compute time-average radiated power
  double power=sqrt(K_squared)*0.5*r*2.0*MathematicalConstants::Pi*
   (imag(C*hp[0])*real(C*h[0])-real(C*hp[0])*imag(C*h[0]));       
 
  return power;
 }

} //end namespace



//=============begin_problem============================================ 
/// Coated disk FSI
//====================================================================== 
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
class CoatedDiskProblem : public Problem
{

public:

 /// Constructor:
 CoatedDiskProblem();
 
 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Recompute gamma integral before checking Newton residuals
 void actions_before_newton_convergence_check()
  {
   Helmholtz_outer_boundary_mesh_pt->setup_gamma();
  }
 
 /// Actions before adapt: Wipe the mesh of traction elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of traction elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

private:

 /// Create FSI traction elements
 void create_fsi_traction_elements();
 
 /// Create Helmholtz FSI flux elements
 void create_helmholtz_fsi_flux_elements();
 
 /// Delete (face) elements in specified mesh 
 void delete_face_elements(Mesh* const & boundary_mesh_pt);
 
 /// Create DtN face elements 
 void create_helmholtz_DtN_elements();

 /// Setup interaction
 void setup_interaction();

 /// Pointer to solid mesh
 TreeBasedRefineableMeshBase* Solid_mesh_pt;

 /// Pointer to mesh of FSI traction elements
 Mesh* FSI_traction_mesh_pt;

 /// Pointer to Helmholtz mesh
 TreeBasedRefineableMeshBase* Helmholtz_mesh_pt;

 /// Pointer to mesh of Helmholtz FSI flux elements
 Mesh* Helmholtz_fsi_flux_mesh_pt;
 
 /// Pointer to mesh containing the DtN elements
 HelmholtzDtNMesh<HELMHOLTZ_ELEMENT>* Helmholtz_outer_boundary_mesh_pt;
 
 /// DocInfo object for output
 DocInfo Doc_info;

 /// Trace file
 ofstream Trace_file;

};


//===========start_of_constructor======================================= 
/// Constructor
//====================================================================== 
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
CoatedDiskProblem<ELASTICITY_ELEMENT, HELMHOLTZ_ELEMENT>::CoatedDiskProblem() 
{

 // The coating mesh is periodic
 bool periodic=true;
 double azimuthal_fraction_of_coating=1.0;
 
 // Solid mesh
 //-----------
 // Number of elements in azimuthal direction
 unsigned ntheta_solid=10*Global_Parameters::El_multiplier;

 // Number of elements in radial direction 
 unsigned nr_solid=3*Global_Parameters::El_multiplier;

 // Innermost radius for solid mesh
 double a=1.0-Global_Parameters::H_coating;
 
 // Build solid mesh
 Solid_mesh_pt = new 
  RefineableTwoDAnnularMesh<ELASTICITY_ELEMENT>
  (periodic,azimuthal_fraction_of_coating,
   ntheta_solid,nr_solid,a,Global_Parameters::H_coating);
 
 
 // Helmholtz mesh
 //---------------

 // Number of elements in azimuthal direction in Helmholtz mesh
 unsigned ntheta_helmholtz=11*Global_Parameters::El_multiplier;

 // Number of elements in radial direction in Helmholtz mesh
 unsigned nr_helmholtz=3*Global_Parameters::El_multiplier;

 // Innermost radius of Helmholtz mesh
 a=1.0;
 
 // Thickness of Helmholtz mesh
 double h_thick_helmholtz=Global_Parameters::Outer_radius-a;

 // Build mesh
 Helmholtz_mesh_pt = new 
  RefineableTwoDAnnularMesh<HELMHOLTZ_ELEMENT>
  (periodic,azimuthal_fraction_of_coating,
   ntheta_helmholtz,nr_helmholtz,a,h_thick_helmholtz);


 // Set error estimators
 Solid_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Helmholtz_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;


 // Mesh containing the Helmholtz DtN
 // elements. Specify outer radius and number of Fourier terms to be
 // used in gamma integral
 unsigned nfourier=20;
 Helmholtz_outer_boundary_mesh_pt = 
  new HelmholtzDtNMesh<HELMHOLTZ_ELEMENT>(Global_Parameters::Outer_radius,
                                          nfourier);
 
 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the solid mesh
 unsigned n_element=Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELASTICITY_ELEMENT *el_pt = 
    dynamic_cast<ELASTICITY_ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law
   el_pt->elasticity_tensor_pt() = &Global_Parameters::E;

   // Square of non-dim frequency
   el_pt->omega_sq_pt()= &Global_Parameters::Omega_sq;
  }
 

 // Same for Helmholtz mesh
 n_element =Helmholtz_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   HELMHOLTZ_ELEMENT *el_pt = 
    dynamic_cast<HELMHOLTZ_ELEMENT*>(Helmholtz_mesh_pt->element_pt(i));

   //Set the pointer to square of Helmholtz wavenumber
   el_pt->k_squared_pt() = &Global_Parameters::K_squared;
  }


 // Output meshes and their boundaries so far so we can double 
 // check the boundary enumeration
 Solid_mesh_pt->output("solid_mesh.dat");
 Helmholtz_mesh_pt->output("helmholtz_mesh.dat");
 Solid_mesh_pt->output_boundaries("solid_mesh_boundary.dat");
 Helmholtz_mesh_pt->output_boundaries("helmholtz_mesh_boundary.dat");


 // Create FaceElement meshes for boundary conditions
 //---------------------------------------------------

 // Construct the fsi traction element mesh
 FSI_traction_mesh_pt=new Mesh;
 create_fsi_traction_elements(); 

 // Construct the Helmholtz fsi flux element mesh
 Helmholtz_fsi_flux_mesh_pt=new Mesh;
 create_helmholtz_fsi_flux_elements();

 // Create DtN elements on outer boundary of Helmholtz mesh
 create_helmholtz_DtN_elements();


 // Combine sub meshes
 //-------------------

 // Solid mesh is first sub-mesh
 add_sub_mesh(Solid_mesh_pt);

 // Add traction sub-mesh
 add_sub_mesh(FSI_traction_mesh_pt);

 // Add Helmholtz mesh
 add_sub_mesh(Helmholtz_mesh_pt);

 // Add Helmholtz FSI flux mesh
 add_sub_mesh(Helmholtz_fsi_flux_mesh_pt);

 // Add Helmholtz DtN mesh
 add_sub_mesh(Helmholtz_outer_boundary_mesh_pt); 

 // Build combined "global" mesh
 build_global_mesh();
 

 // Solid boundary conditions:
 //---------------------------
 // Pin displacements on innermost boundary (boundary 0) of solid mesh
 unsigned n_node = Solid_mesh_pt->nboundary_node(0);
 Vector<std::complex<double> > u(2);
 Vector<double> x(2);
 for(unsigned i=0;i<n_node;i++)
  {
   Node* nod_pt=Solid_mesh_pt->boundary_node_pt(0,i);
   nod_pt->pin(0);
   nod_pt->pin(1);
   nod_pt->pin(2);
   nod_pt->pin(3);

   // Assign displacements
   x[0]=nod_pt->x(0);   
   x[1]=nod_pt->x(1);   
   Global_Parameters::solid_boundary_displacement(x,u);

   // Real part of x-displacement
   nod_pt->set_value(0,u[0].real());

   // Imag part of x-displacement
   nod_pt->set_value(1,u[1].real());

   // Real part of y-displacement
   nod_pt->set_value(2,u[0].imag());

   //Imag part of y-displacement
   nod_pt->set_value(3,u[1].imag());
  }


 // Setup fluid-structure interaction
 //----------------------------------
 setup_interaction();

 // Assign equation  numbers
 oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory(Global_Parameters::Directory);

 // Open trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
  

} //end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the meshes face elements
//========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
actions_before_adapt()
{
 // Kill the fsi traction elements and wipe surface mesh
 delete_face_elements(FSI_traction_mesh_pt);
 
 // Kill Helmholtz FSI flux elements
 delete_face_elements(Helmholtz_fsi_flux_mesh_pt);
 
 // Kill Helmholtz BC elements 
 delete_face_elements(Helmholtz_outer_boundary_mesh_pt);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the meshes of face elements
//========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
actions_after_adapt()
{
 // Create fsi traction elements from all elements that are 
 // adjacent to FSI boundaries and add them to surface meshes
 create_fsi_traction_elements();
 
 // Create Helmholtz fsi flux elements
 create_helmholtz_fsi_flux_elements();
 
 // Create DtN elements from all elements that are 
 // adjacent to the outer boundary of Helmholtz mesh
 create_helmholtz_DtN_elements();
 
 // Setup interaction
 setup_interaction();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
}// end of actions_after_adapt


//============start_of_delete_face_elements================
/// Delete face elements and wipe the mesh
//==========================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
delete_face_elements(Mesh* const & boundary_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = boundary_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   //   Kill surface element
   delete boundary_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 boundary_mesh_pt->flush_element_and_node_storage();
 
} // end of delete_face_elements



//============start_of_create_fsi_traction_elements======================
/// Create fsi traction elements 
//=======================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
create_fsi_traction_elements()
{
 // We're on boundary 2 of the solid mesh
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Solid_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELASTICITY_ELEMENT* bulk_elem_pt = dynamic_cast<ELASTICITY_ELEMENT*>(
    Solid_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
   
   // Create element
   TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement
    <ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>* el_pt=
    new TimeHarmonicLinElastLoadedByHelmholtzPressureBCElement
    <ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>(bulk_elem_pt,
                                           face_index);   
   // Add to mesh
   FSI_traction_mesh_pt->add_element_pt(el_pt);
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   el_pt->set_boundary_number_in_bulk_mesh(b); 
   
   // Set FSI parameter
   el_pt->q_pt()=&Global_Parameters::Q;          
  }
 
} // end of create_traction_elements





//============start_of_create_helmholtz_fsi_flux_elements================
/// Create Helmholtz fsii flux elements 
//=======================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
create_helmholtz_fsi_flux_elements()
{
 
 // Attach to inner boundary of Helmholtz mesh (0)
 unsigned b=0;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Helmholtz_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   HELMHOLTZ_ELEMENT* bulk_elem_pt = dynamic_cast<HELMHOLTZ_ELEMENT*>(
    Helmholtz_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Helmholtz_mesh_pt->face_index_at_boundary(b,e);
   
   // Create element
   HelmholtzFluxFromNormalDisplacementBCElement
    <HELMHOLTZ_ELEMENT,ELASTICITY_ELEMENT>* el_pt=
    new HelmholtzFluxFromNormalDisplacementBCElement
    <HELMHOLTZ_ELEMENT,ELASTICITY_ELEMENT>(bulk_elem_pt,
                                           face_index);
   
   // Add to mesh
   Helmholtz_fsi_flux_mesh_pt->add_element_pt(el_pt);
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   el_pt->set_boundary_number_in_bulk_mesh(b); 
  }  
  
} // end of create_helmholtz_flux_elements



//============start_of_create_DtN_elements==============================
/// Create DtN elements on the outer boundary of 
/// the Helmholtz mesh
//===========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
create_helmholtz_DtN_elements()
{
 // We're on boundary 2 of the Helmholtz mesh
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Helmholtz_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   HELMHOLTZ_ELEMENT* bulk_elem_pt = dynamic_cast<HELMHOLTZ_ELEMENT*>(
    Helmholtz_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = Helmholtz_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding DtN element
   HelmholtzDtNBoundaryElement<HELMHOLTZ_ELEMENT>* flux_element_pt = new 
    HelmholtzDtNBoundaryElement<HELMHOLTZ_ELEMENT>(bulk_elem_pt,face_index);
   
   // Set pointer to the mesh that contains all the boundary condition
   // elements on this boundary
   flux_element_pt->set_outer_boundary_mesh_pt(
    Helmholtz_outer_boundary_mesh_pt);  

   //Add the DtN element to the mesh
   Helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);

  }

} // end of create_helmholtz_DtN_elements



//=====================start_of_setup_interaction======================
/// Setup interaction between two fields
//========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
setup_interaction()
{
 // Setup Helmholtz "pressure" load on traction elements
 unsigned boundary_in_helmholtz_mesh=0;
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <HELMHOLTZ_ELEMENT,2>
  (this,boundary_in_helmholtz_mesh,Helmholtz_mesh_pt,FSI_traction_mesh_pt);

 // Setup Helmholtz flux from normal displacement interaction
 unsigned boundary_in_solid_mesh=2;
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <ELASTICITY_ELEMENT,2>(
   this,boundary_in_solid_mesh,Solid_mesh_pt,Helmholtz_fsi_flux_mesh_pt);
}



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedDiskProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::doc_solution()
{

 ofstream some_file,some_file2;
 char filename[100];

 // Number of plot points
 unsigned n_plot=5; 

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 // Accumulate contribution from elements
 double power=0.0;
 unsigned nn_element=Helmholtz_outer_boundary_mesh_pt->nelement(); 
 for(unsigned e=0;e<nn_element;e++)
  {
   HelmholtzBCElementBase<HELMHOLTZ_ELEMENT> *el_pt = 
    dynamic_cast<HelmholtzBCElementBase<HELMHOLTZ_ELEMENT>*>(
     Helmholtz_outer_boundary_mesh_pt->element_pt(e)); 
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();
 oomph_info << "Step: " << Doc_info.number() 
            << ": Q=" << Global_Parameters::Q  << "\n"
            << " k_squared=" << Global_Parameters::K_squared  << "\n"
            << " density ratio=" << Global_Parameters::Density_ratio  << "\n"
            << " omega_sq=" << Global_Parameters::Omega_sq  << "\n"
            << " Total radiated power " << power  << "\n"
            << " Axisymmetric radiated power "  << "\n"
            <<  Global_Parameters::exact_axisym_radiated_power()  << "\n"
            << std::endl; 


 // Write trace file
 Trace_file << Global_Parameters::Q << " " 
            << Global_Parameters::K_squared << " "
            << Global_Parameters::Density_ratio << " "
            << Global_Parameters::Omega_sq << " "
            << power << " " 
            << Global_Parameters::exact_axisym_radiated_power() << " " 
            << std::endl;
  
 std::ostringstream case_string;
 case_string << "TEXT X=10,Y=90, T=\"Q=" 
             <<  Global_Parameters::Q 
             << ",  k<sup>2</sup>="
             <<  Global_Parameters::K_squared
             << ",  density ratio="
             <<  Global_Parameters::Density_ratio
             << ",  omega_sq="
             <<  Global_Parameters::Omega_sq
             << "\"\n";


 // Output displacement field
 //--------------------------
 sprintf(filename,"%s/elast_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,n_plot);
 some_file.close();


 // Output fsi traction elements
 //----------------------------- 
 sprintf(filename,"%s/traction_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 FSI_traction_mesh_pt->output(some_file,n_plot);
 some_file.close();


 // Output Helmholtz fsi flux elements
 //----------------------------------- 
 sprintf(filename,"%s/flux_bc_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Helmholtz_fsi_flux_mesh_pt->output(some_file,n_plot);
 some_file.close();


 // Output Helmholtz
 //-----------------
 sprintf(filename,"%s/helmholtz_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Helmholtz_mesh_pt->output(some_file,n_plot);
 some_file << case_string.str();
 some_file.close();


 // Output exact solution for Helmholtz 
 //------------------------------------
 sprintf(filename,"%s/exact_helmholtz_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Helmholtz_mesh_pt->output_fct(some_file,n_plot,
                                 Global_Parameters::exact_axisym_potential); 
 some_file.close();
 

 cout << "Doced for Q=" << Global_Parameters::Q << " (step "
      << Doc_info.number() << ")\n";

 // Increment label for output files
 Doc_info.number()++;
 
} //end doc



//=======start_of_main==================================================
/// Driver for acoustic fsi problem
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
 
 // Azimuthal wavenumber of forcing
 CommandLineArgs::specify_command_line_flag("--n",&Global_Parameters::N);

 // Minimum refinement level
 CommandLineArgs::specify_command_line_flag("--el_multiplier",
                            &Global_Parameters::El_multiplier);
 
 // Outer radius of Helmholtz domain
 CommandLineArgs::specify_command_line_flag("--outer_radius",
                            &Global_Parameters::Outer_radius);
 
 // Number of steps in parameter study
 unsigned nstep=2;
 CommandLineArgs::specify_command_line_flag("--nstep",&nstep);
 
 // Increment in FSI parameter in parameter study
 double q_increment=5.0;
 CommandLineArgs::specify_command_line_flag("--q_increment",&q_increment);
 
 // Max. number of adaptations
 unsigned max_adapt=3;
 CommandLineArgs::specify_command_line_flag("--max_adapt",&max_adapt);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Set up the problem
 CoatedDiskProblem<RefineableQTimeHarmonicLinearElasticityElement<2,3>,
                   RefineableQHelmholtzElement<2,3> > problem; 

  
 // Initial values for parameter values
 Global_Parameters::Q=0.0; 
 Global_Parameters::update_parameter_values();

 //Parameter incrementation
 for(unsigned i=0;i<nstep;i++)
  {
   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   problem.newton_solve(max_adapt);

   // Doc solution
   problem.doc_solution();

   // Increment FSI parameter
   Global_Parameters::Q+=q_increment;
   Global_Parameters::update_parameter_values();
  }
 
} //end of main









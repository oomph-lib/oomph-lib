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
//Driver for Fourier-decomposed  Helmholtz problem 

#include <complex>
#include <cmath>

//Generic routines
#include "generic.h"

// The Helmholtz equations
#include "fourier_decomposed_helmholtz.h"
 
// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//========================================================================
// AnnularQuadMesh, derived from SimpleRectangularQuadMesh.
//========================================================================
template<class ELEMENT> 
class AnnularQuadMesh : public SimpleRectangularQuadMesh<ELEMENT>
{
 
  public:

 // Constructor for angular mesh with n_r x n_phi 
 // 2D quad elements. Calls constructor for the underlying 
 // SimpleRectangularQuadMesh; then deforms the mesh so that it fits 
 // into the annular region bounded by the radii r_min and r_max
 // and angles (in degree) of phi_min and phi_max.
 AnnularQuadMesh(const unsigned& n_r, const unsigned& n_phi,
                 const double& r_min, const double& r_max,
                 const double& phi_min, const double& phi_max) :
   SimpleRectangularQuadMesh<ELEMENT>(n_r,n_phi,1.0,1.0)
  {

   // The constructor for the  SimpleRectangularQuadMesh has
   // built the mesh with n_x x n_y = n_r x n_phi elements in the unit
   // square. Let's reposition the nodal points so that the mesh
   // gets mapped into the required annular region:

   // Find out how many nodes there are
   unsigned n_node=this->nnode();
   
   // Loop over all nodes
   for (unsigned n=0;n<n_node;n++)
    {
     // Pointer to node:
     Node* nod_pt=this->node_pt(n);
     
     // Get the x/y coordinates
     double x_old=nod_pt->x(0);
     double y_old=nod_pt->x(1);

     // Map from the old x/y to the new r/phi:
     double r=r_min+(r_max-r_min)*x_old;
     double phi=(phi_min+(phi_max-phi_min)*y_old)*
      MathematicalConstants::Pi/180.0;

     // Set new nodal coordinates
     nod_pt->x(0)=r*cos(phi);
     nod_pt->x(1)=r*sin(phi);
    }
  }

};



//===== start_of_namespace_planar_wave=================================
/// Namespace to test representation of planar wave in spherical
/// polars
//=====================================================================
namespace PlanarWave
{

 /// Number of terms in series
 unsigned N_terms=100;

 /// Wave number
 double K=3.0*MathematicalConstants::Pi;

 /// Imaginary unit 
 std::complex<double> I(0.0,1.0); 

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);
  
  double theta;
  theta=atan2(x[0],x[1]);
  
  // Argument for Bessel/Hankel functions
  double kr = K*R;  
  
  // Need half-order Bessel functions
  double bessel_offset=0.5;

  // Evaluate Bessel/Hankel functions
  Vector<double> jv(N_terms);
  Vector<double> yv(N_terms);
  Vector<double> djv(N_terms);
  Vector<double> dyv(N_terms);
  double order_max_in=double(N_terms-1)+bessel_offset;
  double order_max_out=0;
  
  // This function returns vectors containing 
  // J_k(x), Y_k(x) and their derivatives
  // up to k=order_max, with k increasing in
  // integer increments starting with smallest
  // positive value. So, e.g. for order_max=3.5
  // jv[0] contains J_{1/2}(x),
  // jv[1] contains J_{3/2}(x),
  // jv[2] contains J_{5/2}(x),
  // jv[3] contains J_{7/2}(x).
  CRBond_Bessel::bessjyv(order_max_in,
                         kr,
                         order_max_out,
                         &jv[0],&yv[0],
                         &djv[0],&dyv[0]);
  
  
  // Assemble  exact solution (actually no need to add terms
  // below i=N_fourier as Legendre polynomial would be zero anyway)
  complex<double> u_ex(0.0,0.0);
  for(unsigned i=0;i<N_terms;i++)
   {
    //Associated_legendre_functions
    double p=Legendre_functions_helper::plgndr2(i,0,cos(theta));
    
    // Set exact solution
    u_ex+=(2.0*i+1.0)*pow(I,i)*
     sqrt(MathematicalConstants::Pi/(2.0*kr))*jv[i]*p;
   }
  
  // Get the real & imaginary part of the result
  u[0]=u_ex.real();
  u[1]=u_ex.imag();
  
 }//end of get_exact_u


 /// Plot 
 void plot()
 {
  unsigned nr=20;
  unsigned nz=100;
  unsigned nt=40;

  ofstream some_file("planar_wave.dat");

  for (unsigned i_t=0;i_t<nt;i_t++)
   {
    double t=2.0*MathematicalConstants::Pi*double(i_t)/double(nt-1);

    some_file << "ZONE I="<< nz << ", J="<< nr << std::endl;
    
    Vector<double> x(2);
    Vector<double> u(2);
    for (unsigned i=0;i<nr;i++)
     {
      x[0]=0.001+double(i)/double(nr-1);
      for (unsigned j=0;j<nz;j++)
       {
        x[1]=double(j)/double(nz-1);
        get_exact_u(x,u); 
        complex<double> uu=complex<double>(u[0],u[1])*exp(-I*t);
        some_file << x[0] << " " << x[1] << " " 
                  << uu.real() << " " << uu.imag() << "\n";
       }
     } 
   }
 }
 
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//===== start_of_namespace=============================================
/// Namespace for the Fourier decomposed Helmholtz problem parameters
//=====================================================================
namespace ProblemParameters
{
 /// Square of the wavenumber
 double K_squared=10.0;
 
 /// Fourier wave number
 int N_fourier=3;
 
 /// Number of terms in computation of DtN boundary condition
 unsigned Nterms_for_DtN=6;

 /// Number of terms in the exact solution
 unsigned N_terms=6; 
 
 /// Coefficients in the exact solution
 Vector<double> Coeff(N_terms,1.0);

 /// Imaginary unit 
 std::complex<double> I(0.0,1.0); 

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);
  
  double theta;
  theta=atan2(x[0],x[1]);
  
  // Argument for Bessel/Hankel functions
  double kr = sqrt(K_squared)*R;  
  
  // Need half-order Bessel functions
  double bessel_offset=0.5;

  // Evaluate Bessel/Hankel functions
  Vector<double> jv(N_terms);
  Vector<double> yv(N_terms);
  Vector<double> djv(N_terms);
  Vector<double> dyv(N_terms);
  double order_max_in=double(N_terms-1)+bessel_offset;
  double order_max_out=0;
  
  // This function returns vectors containing 
  // J_k(x), Y_k(x) and their derivatives
  // up to k=order_max, with k increasing in
  // integer increments starting with smallest
  // positive value. So, e.g. for order_max=3.5
  // jv[0] contains J_{1/2}(x),
  // jv[1] contains J_{3/2}(x),
  // jv[2] contains J_{5/2}(x),
  // jv[3] contains J_{7/2}(x).
  CRBond_Bessel::bessjyv(order_max_in,
                         kr,
                         order_max_out,
                         &jv[0],&yv[0],
                         &djv[0],&dyv[0]);
  
  
  // Assemble  exact solution (actually no need to add terms
  // below i=N_fourier as Legendre polynomial would be zero anyway)
  complex<double> u_ex(0.0,0.0);
  for(unsigned i=N_fourier;i<N_terms;i++)
   {
    //Associated_legendre_functions
    double p=Legendre_functions_helper::plgndr2(i,N_fourier,
                                                cos(theta));
    // Set exact solution
    u_ex+=Coeff[i]*sqrt(MathematicalConstants::Pi/(2.0*kr))*(jv[i]+I*yv[i])*p;
   }
  
  // Get the real & imaginary part of the result
  u[0]=u_ex.real();
  u[1]=u_ex.imag();
  
 }//end of get_exact_u

 
 /// Get -du/dr (spherical r) for exact solution. Equal to prescribed
 /// flux on inner boundary.
 void exact_minus_dudr(const Vector<double>& x, std::complex<double>& flux)
 {
  // Initialise flux
  flux=std::complex<double>(0.0,0.0);
  
  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);
  
  double theta;
  theta=atan2(x[0],x[1]);
  
  // Argument for Bessel/Hankel functions
  double kr=sqrt(K_squared)*R;  

  // Helmholtz wavenumber
  double k=sqrt(K_squared);

  // Need half-order Bessel functions
  double bessel_offset=0.5;

  // Evaluate Bessel/Hankel functions
  Vector<double> jv(N_terms);
  Vector<double> yv(N_terms);
  Vector<double> djv(N_terms);
  Vector<double> dyv(N_terms);
  double order_max_in=double(N_terms-1)+bessel_offset;
  double order_max_out=0;
  
  // This function returns vectors containing 
  // J_k(x), Y_k(x) and their derivatives
  // up to k=order_max, with k increasing in
  // integer increments starting with smallest
  // positive value. So, e.g. for order_max=3.5
  // jv[0] contains J_{1/2}(x),
  // jv[1] contains J_{3/2}(x),
  // jv[2] contains J_{5/2}(x),
  // jv[3] contains J_{7/2}(x).
  CRBond_Bessel::bessjyv(order_max_in,
                         kr,
                         order_max_out,
                         &jv[0],&yv[0],
                         &djv[0],&dyv[0]);
  
  
  // Assemble  exact solution (actually no need to add terms
  // below i=N_fourier as Legendre polynomial would be zero anyway)
  complex<double> u_ex(0.0,0.0);
  for(unsigned i=N_fourier;i<N_terms;i++)
   {
    //Associated_legendre_functions
    double p=Legendre_functions_helper::plgndr2(i,N_fourier,
                                                cos(theta));
    // Set flux of exact solution
    flux-=Coeff[i]*sqrt(MathematicalConstants::Pi/(2.0*kr))*p*
     ( k*(djv[i]+I*dyv[i]) - (0.5*(jv[i]+I*yv[i])/R) );
   }
  
 }// end of exact_normal_derivative
 
 
 /// Multiplier for number of elements
 unsigned El_multiplier=1;
 

} // end of namespace



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//========= start_of_problem_class=====================================
/// Problem class 
//=====================================================================
template<class ELEMENT> 
class FourierDecomposedHelmholtzProblem : public Problem
{
 
public:
 
 /// Constructor
 FourierDecomposedHelmholtzProblem();
 
 /// Destructor (empty)
 ~FourierDecomposedHelmholtzProblem(){}
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
 
 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
  
 /// Recompute gamma integral before checking Newton residuals
 void actions_before_newton_convergence_check()
  {
   Helmholtz_outer_boundary_mesh_pt->setup_gamma();
  }

 /// Check gamma computation
 void check_gamma(DocInfo& doc_info);
  
private:
 
 /// Create BC elements on outer boundary
 void create_outer_bc_elements();
 
 /// Create flux elements on inner boundary
 void create_flux_elements_on_inner_boundary();

 /// Pointer to bulk mesh
 AnnularQuadMesh<ELEMENT>* Bulk_mesh_pt;
  
 /// Pointer to mesh containing the DtN boundary
 /// condition elements
 FourierDecomposedHelmholtzDtNMesh<ELEMENT>* Helmholtz_outer_boundary_mesh_pt;

 /// Mesh of face elements that apply the prescribed flux
 /// on the inner boundary
 Mesh* Helmholtz_inner_boundary_mesh_pt;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Fourier-decomposed Helmholtz problem
//========================================================================
template<class ELEMENT>
FourierDecomposedHelmholtzProblem<ELEMENT>::
FourierDecomposedHelmholtzProblem()
{ 

// Build annular mesh
  
// # of elements in r-direction 
 unsigned n_r=10*ProblemParameters::El_multiplier;
 
 // # of elements in theta-direction 
 unsigned n_theta=10*ProblemParameters::El_multiplier;
 
 // Domain boundaries in theta-direction
 double theta_min=-90.0;
 double theta_max=90.0;
 
 // Domain boundaries in r-direction
 double r_min=1.0;
 double r_max=3.0;
 
 // Build and assign mesh
 Bulk_mesh_pt = 
  new AnnularQuadMesh<ELEMENT>(n_r,n_theta,r_min,r_max,theta_min,theta_max);

 // Create mesh for DtN elements on outer boundary
 Helmholtz_outer_boundary_mesh_pt=
  new FourierDecomposedHelmholtzDtNMesh<ELEMENT>(
   r_max,ProblemParameters::Nterms_for_DtN);

 // Populate it with elements
 create_outer_bc_elements();

 // Create flux elements on inner boundary
 Helmholtz_inner_boundary_mesh_pt=new Mesh;
 create_flux_elements_on_inner_boundary();
 
 // Add the several  sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt); 
 add_sub_mesh(Helmholtz_inner_boundary_mesh_pt); 
 add_sub_mesh(Helmholtz_outer_boundary_mesh_pt); 
  
 // Build the Problem's global mesh from its various sub-meshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
   
   //Set the k_squared pointer
   el_pt->k_squared_pt()=&ProblemParameters::K_squared;
   
   // Set pointer to Fourier wave number
   el_pt->fourier_wavenumber_pt()=&ProblemParameters::N_fourier;
  }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=================================start_of_check_gamma===================
/// Check gamma computation: \f$ \gamma = -du/dn \f$
//========================================================================
template<class ELEMENT>
void FourierDecomposedHelmholtzProblem<ELEMENT>::check_gamma(DocInfo& doc_info)
{
 

 // Compute gamma stuff
 Helmholtz_outer_boundary_mesh_pt->setup_gamma();
 
 ofstream some_file;
 char filename[100];

 sprintf(filename,"%s/gamma_test%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
  
 //first loop over elements e
 unsigned nel=Helmholtz_outer_boundary_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Get a pointer to element
   FourierDecomposedHelmholtzDtNBoundaryElement<ELEMENT>* el_pt=
    dynamic_cast<FourierDecomposedHelmholtzDtNBoundaryElement<ELEMENT>*>
    (Helmholtz_outer_boundary_mesh_pt->element_pt(e));
   
   //Set the value of n_intpt
   const unsigned n_intpt =el_pt->integral_pt()->nweight();
   
   // Get gamma at all gauss points in element
   Vector<std::complex<double> > gamma(
    Helmholtz_outer_boundary_mesh_pt->gamma_at_gauss_point(el_pt));
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Allocate and initialise coordiante
     Vector<double> x(el_pt->dim()+1,0.0);
     
     //Set the Vector to hold local coordinates
     unsigned n=el_pt->dim();
     Vector<double> s(n,0.0);
     for(unsigned i=0;i<n;i++)
      {
       s[i]=el_pt->integral_pt()->knot(ipt,i);
      }
     
     //Get the coordinates of the integration point
     el_pt->interpolated_x(s,x);
     
     complex<double> flux;
     ProblemParameters::exact_minus_dudr(x,flux);
     some_file << atan2(x[0],x[1]) << " " 
               << gamma[ipt].real() << " "
               << gamma[ipt].imag() << " "
               << flux.real() << " " 
               << flux.imag() << " " 
               << std::endl;
     
    }// end of loop over integration points
   
  }// end of loop over elements
 
 some_file.close();
  
}//end of output_gamma


//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void FourierDecomposedHelmholtzProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;
  
 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,ProblemParameters::get_exact_u); 
 some_file.close();
 
 
 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,ProblemParameters::get_exact_u,
                             error,norm); 
 some_file.close();
 
 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;
 
 // Check gamma computation
 check_gamma(doc_info);


 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 sprintf(filename,"%s/total_power%i.dat",doc_info.directory().c_str(),
         doc_info.number());

 // Accumulate contribution from elements
 double power=0.0;
 unsigned nn_element=Helmholtz_outer_boundary_mesh_pt->nelement(); 
 for(unsigned e=0;e<nn_element;e++)
  {
   FourierDecomposedHelmholtzBCElementBase<ELEMENT> *el_pt = 
    dynamic_cast<FourierDecomposedHelmholtzBCElementBase<ELEMENT>*>(
     Helmholtz_outer_boundary_mesh_pt->element_pt(e)); 
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();

 // Output total power
 oomph_info << "Radiated power: " 
            << ProblemParameters::N_fourier << " " 
            << ProblemParameters::K_squared << " " 
            << power << std::endl;
 some_file.open(filename); 
 some_file << ProblemParameters::N_fourier << " " 
           << ProblemParameters::K_squared << " " 
           << power << std::endl;
 some_file.close();

} // end of doc



//============start_of_create_outer_bc_elements==============================
/// Create BC elements on outer boundary
//========================================================================
template<class ELEMENT>
void FourierDecomposedHelmholtzProblem<ELEMENT>::create_outer_bc_elements()
{
 // Outer boundary is boundary 1:
 unsigned b=1;

 // Loop over the bulk elements adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding DtN element
   FourierDecomposedHelmholtzDtNBoundaryElement<ELEMENT>* flux_element_pt = new 
    FourierDecomposedHelmholtzDtNBoundaryElement<ELEMENT>(bulk_elem_pt,
                                                          face_index);
   
   //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
   Helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);

   // Set pointer to the mesh that contains all the boundary condition
   // elements on this boundary
   flux_element_pt->
    set_outer_boundary_mesh_pt(Helmholtz_outer_boundary_mesh_pt);
  }

} // end of create_outer_bc_elements



//============start_of_create_flux_elements=================
/// Create flux elements on inner boundary
//==========================================================
template<class ELEMENT>
void  FourierDecomposedHelmholtzProblem<ELEMENT>::
create_flux_elements_on_inner_boundary()
{
 // Apply flux bc on inner boundary (boundary 3)
 unsigned b=3;

// Loop over the bulk elements adjacent to boundary b
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding prescribed incoming-flux element
   FourierDecomposedHelmholtzFluxElement<ELEMENT>* flux_element_pt = new 
    FourierDecomposedHelmholtzFluxElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed incoming-flux element to the surface mesh
   Helmholtz_inner_boundary_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set the pointer to the prescribed flux function
   flux_element_pt->flux_fct_pt() = &ProblemParameters::exact_minus_dudr;

  } //end of loop over bulk elements adjacent to boundary b
 
} // end of create flux elements on inner boundary



//===== start_of_main=====================================================
/// Driver code for Fourier decomposed Helmholtz problem
//========================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Multiplier for number of elements
 CommandLineArgs::specify_command_line_flag("--el_multiplier",
                                            &ProblemParameters::El_multiplier);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();


 // Check if the claimed representation of a planar wave in
 // the tutorial is correct -- of course it is!
 //PlanarWave::plot();

 // Test Bessel/Hankel functions
 //-----------------------------
 {
  // Number of Bessel functions to be computed
  unsigned n=3;
  
  // Offset of Bessel function order (less than 1!)
  double bessel_offset=0.5;
  
  ofstream bessely_file("besselY.dat");
  ofstream bessely_deriv_file("dbesselY.dat");
  
  ofstream besselj_file("besselJ.dat");
  ofstream besselj_deriv_file("dbesselJ.dat");
  
  // Evaluate Bessel/Hankel functions
  Vector<double> jv(n+1);
  Vector<double> yv(n+1);
  Vector<double> djv(n+1);
  Vector<double> dyv(n+1);
  double x_min=0.5;
  double x_max=5.0;
  unsigned nplot=100;
  for (unsigned i=0;i<nplot;i++)
   {
    double x=x_min+(x_max-x_min)*double(i)/double(nplot-1);
    double order_max_in=double(n)+bessel_offset;
    double order_max_out=0;
    
    // This function returns vectors containing 
    // J_k(x), Y_k(x) and their derivatives
    // up to k=order_max, with k increasing in
    // integer increments starting with smallest
    // positive value. So, e.g. for order_max=3.5
    // jv[0] contains J_{1/2}(x),
    // jv[1] contains J_{3/2}(x),
    // jv[2] contains J_{5/2}(x),
    // jv[3] contains J_{7/2}(x).
    CRBond_Bessel::bessjyv(order_max_in,x,
                           order_max_out,
                           &jv[0],&yv[0],
                           &djv[0],&dyv[0]);
    bessely_file << x << " ";
    for (unsigned j=0;j<=n;j++)
     {
      bessely_file << yv[j] << " ";
     }
    bessely_file << std::endl;
    
    besselj_file << x << " ";
    for (unsigned j=0;j<=n;j++)
     {
      besselj_file << jv[j] << " ";
     }
    besselj_file << std::endl;
    
    bessely_deriv_file << x << " ";
    for (unsigned j=0;j<=n;j++)
     {
      bessely_deriv_file << dyv[j] << " ";
     }
    bessely_deriv_file << std::endl;
    
    besselj_deriv_file << x << " ";
    for (unsigned j=0;j<=n;j++)
     {
      besselj_deriv_file << djv[j] << " ";
     }
    besselj_deriv_file << std::endl;
    
   }
  bessely_file.close();
  besselj_file.close();
  bessely_deriv_file.close();
  besselj_deriv_file.close();
 }
 
 
 // Test Legendre Polynomials
 //--------------------------
 {
  // Fourier wavenumber
  unsigned n=3;
    
  ofstream some_file("legendre3.dat");
  unsigned nplot=100;
  for (unsigned i=0;i<nplot;i++)
   {
    double x=double(i)/double(nplot-1)*2.0*MathematicalConstants::Pi;

    some_file << x << " ";
    for (unsigned j=n;j<=5;j++)
     {
      some_file <<  Legendre_functions_helper::plgndr2(j,n,cos(x)) << " ";
     }
    some_file << std::endl;
   }
  some_file.close();
 }


 {
  ofstream some_file("legendre.dat");
  unsigned nplot=100;
  for (unsigned i=0;i<nplot;i++)
   {
    double x=double(i)/double(nplot-1);

    some_file << x << " ";
    for (unsigned j=0;j<=3;j++)
     {
      some_file <<  Legendre_functions_helper::plgndr2(j,0,x) << " ";
     }
    some_file << std::endl;
   }
  some_file.close();
 }


 
 // Create the problem with 2D nine-node elements from the
 // QFourierDecomposedHelmholtzElement family. 
 FourierDecomposedHelmholtzProblem<QFourierDecomposedHelmholtzElement<3> > 
  problem;
 
 // Create label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Solve for a few Fourier wavenumbers
 for (ProblemParameters::N_fourier=0;ProblemParameters::N_fourier<4;
      ProblemParameters::N_fourier++)
  {
   // Step number
   doc_info.number()=ProblemParameters::N_fourier;
   
   // Solve the problem
   problem.newton_solve();
   
   //Output the solution
   problem.doc_solution(doc_info);
  }
 
} //end of main








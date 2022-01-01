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
// Driver for Helmholtz/TimeHarmonicTimeHarmonicLinElast coupling
#include <complex>
#include <cmath>

//Oomph-lib includes
#include "generic.h"

//The Helmholtz equation
#include "fourier_decomposed_helmholtz.h"

//The Elasticity equation
#include "time_harmonic_fourier_decomposed_linear_elasticity.h"

// The interaction elements
#include "multi_physics.h"

// The mesh
#include "meshes/annular_mesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;


/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////


//=======start_of_namespace==========================================
/// Global variables
//===================================================================
namespace Global_Parameters
{

 /// Square of wavenumber for the Helmholtz equation
 double K_squared=10.0;
 
 /// Radius of outer boundary of Helmholtz domain
 double Outer_radius=2.0; 

 /// FSI parameter
 double Q=10.0;

 /// Non-dim thickness of elastic coating
 double H_coating=0.2; 

 /// Define azimuthal Fourier wavenumber
 int Fourier_wavenumber=0;
   
 /// Poisson's ratio Nu
 std::complex<double> Nu(std::complex<double>(0.3,0.0));

 /// Non-dim square of frequency for solid -- dependent variable!
 std::complex<double>  Omega_sq(std::complex<double>(100.0,0.0));

 /// Density ratio: solid to fluid
 double Density_ratio=1.0; 

 /// Function to update dependent parameter values
 void update_parameter_values()
 {
  Omega_sq=Density_ratio*Q;
 }

 /// Wavenumber "zenith"-variation of imposed displacement of coating
 /// on inner boundary 
 unsigned M=4; 

 /// Displacement field on inner boundary of solid
 void solid_boundary_displacement(const Vector<double>& x,
                                  Vector<std::complex<double> >& u)
 {
  Vector<double> normal(2);
  double norm=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta=atan2(x[1],x[0]);
  normal[0]=x[0]/norm;
  normal[1]=x[1]/norm;

  u[0]=complex<double>(normal[0]*cos(double(M)*theta),0.0);
  u[1]=complex<double>(normal[1]*cos(double(M)*theta),0.0);
 }


 /// Output directory
 string Directory="RESLT";
 
 /// Multiplier for number of elements
 unsigned El_multiplier=1;

} //end_of_namespace



//=============start_of_problem_class===================================
/// Coated sphere FSI
//====================================================================== 
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
class CoatedSphereProblem : public Problem
{

public:

 /// Constructor:
 CoatedSphereProblem();
 
 /// Update function (empty)
 void actions_before_newton_solve(){}

 /// Update function (empty)
 void actions_after_newton_solve() {}
 
 /// Recompute gamma integral before checking Newton residuals
 void actions_before_newton_convergence_check()
  {
   Helmholtz_DtN_mesh_pt->setup_gamma();
  }
  
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Create FSI traction elements
 void create_fsi_traction_elements();

 /// Create Helmholtz FSI flux elements
 void create_helmholtz_fsi_flux_elements(); 

 /// Setup interaction
 void setup_interaction();

 /// Create DtN elements on outer boundary
 void create_helmholtz_DtN_elements();

 /// Pointer to solid mesh
 TwoDAnnularMesh<ELASTICITY_ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of FSI traction elements
 Mesh* FSI_traction_mesh_pt;

 /// Pointer to Helmholtz mesh
 TwoDAnnularMesh<HELMHOLTZ_ELEMENT>* Helmholtz_mesh_pt;

 /// Pointer to mesh of Helmholtz FSI flux elements
 Mesh* Helmholtz_fsi_flux_mesh_pt;
 
 /// Pointer to mesh containing the DtN elements
 FourierDecomposedHelmholtzDtNMesh<HELMHOLTZ_ELEMENT>* Helmholtz_DtN_mesh_pt;
 
 /// Trace file
 ofstream Trace_file;

};// end_of_problem_class


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
CoatedSphereProblem<ELASTICITY_ELEMENT, HELMHOLTZ_ELEMENT>::
CoatedSphereProblem() 
{

 // Parameters for meshes
 bool periodic=false;
 double azimuthal_fraction_of_coating=0.5;
 double phi=0.5*MathematicalConstants::Pi;

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
  TwoDAnnularMesh<ELASTICITY_ELEMENT>(periodic,azimuthal_fraction_of_coating,
                                      ntheta_solid,nr_solid,a,
                                      Global_Parameters::H_coating,phi);
 
 
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
 Helmholtz_mesh_pt = new TwoDAnnularMesh<HELMHOLTZ_ELEMENT>
  (periodic,azimuthal_fraction_of_coating,
   ntheta_helmholtz,nr_helmholtz,a,h_thick_helmholtz,phi);


 // Create mesh for DtN elements on outer boundary
 unsigned nfourier=20;
 Helmholtz_DtN_mesh_pt=
  new FourierDecomposedHelmholtzDtNMesh<HELMHOLTZ_ELEMENT>(
   Global_Parameters::Outer_radius,nfourier);

 // Complete the solid problem setup to make the elements fully functional
 unsigned nel=Solid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {     
   // Cast to a bulk element
   ELASTICITY_ELEMENT* el_pt=dynamic_cast<ELASTICITY_ELEMENT*>(
    Solid_mesh_pt->element_pt(e));

   // Set the pointer to Fourier wavenumber
   el_pt->fourier_wavenumber_pt() = &Global_Parameters::Fourier_wavenumber;
   
   // Set the pointer to Poisson's ratio
   el_pt->nu_pt() = &Global_Parameters::Nu;
      
   // Set the pointer to square of the angular frequency
   el_pt->omega_sq_pt() = &Global_Parameters::Omega_sq;
  }

 // Complete the build of all Helmholtz elements so they are fully functional
 unsigned n_element = Helmholtz_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   HELMHOLTZ_ELEMENT *el_pt = dynamic_cast<HELMHOLTZ_ELEMENT*>(
    Helmholtz_mesh_pt->element_pt(i));
   
   //Set the k_squared pointer
   el_pt->k_squared_pt()=&Global_Parameters::K_squared;
   
   // Set pointer to Fourier wave number
   el_pt->fourier_wavenumber_pt()=&Global_Parameters::Fourier_wavenumber;
  }

 // Output meshes and their boundaries so far so we can double 
 // check the boundary enumeration
 Solid_mesh_pt->output("solid_mesh.dat");
 Helmholtz_mesh_pt->output("helmholtz_mesh.dat");
 Solid_mesh_pt->output_boundaries("solid_mesh_boundary.dat");
 Helmholtz_mesh_pt->output_boundaries("helmholtz_mesh_boundary.dat");

 // Create FaceElement meshes for boundary conditions
 //--------------------------------------------------
 
 // Construct the fsi traction element mesh
 FSI_traction_mesh_pt=new Mesh;
 create_fsi_traction_elements();
 
 // Construct the Helmholtz fsi flux element mesh
 Helmholtz_fsi_flux_mesh_pt=new Mesh;
 create_helmholtz_fsi_flux_elements();
 
 // Create DtN elements
 create_helmholtz_DtN_elements();


 // Combine sub meshes
 //-------------------
 add_sub_mesh(Solid_mesh_pt);
 add_sub_mesh(FSI_traction_mesh_pt);
 add_sub_mesh(Helmholtz_mesh_pt);
 add_sub_mesh(Helmholtz_fsi_flux_mesh_pt);
 add_sub_mesh(Helmholtz_DtN_mesh_pt); 
 
 // Build the Problem's global mesh from its various sub-meshes
 build_global_mesh();
 

 // Solid boundary conditions:
 //---------------------------

 // Pin the solid inner boundary (boundary 0) in all directions
 unsigned b=0;
 unsigned n_node = Solid_mesh_pt->nboundary_node(b);
 
 Vector<std::complex<double> > u(2);
 Vector<double> x(2);

 //Loop over the nodes to pin and assign boundary displacements on 
 //solid boundary
 for(unsigned i=0;i<n_node;i++)
  {
   Node* nod_pt=Solid_mesh_pt->boundary_node_pt(b,i);
   nod_pt->pin(0);
   nod_pt->pin(1);
   nod_pt->pin(2);
   nod_pt->pin(3);
   nod_pt->pin(4);
   nod_pt->pin(5);

   // Assign prescribed displacements
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   Global_Parameters::solid_boundary_displacement(x,u);

   // Real part of radial displacement
   nod_pt->set_value(0,u[0].real());
   // Real part of axial displacement
   nod_pt->set_value(1,u[1].real());
   // Real part of azimuthal displacement
   nod_pt->set_value(2,0.0);
   // Imag part of radial displacement
   nod_pt->set_value(3,u[0].imag());
   // Imag part of axial displacement
   nod_pt->set_value(4,u[1].imag());
   // Imag part of azimuthal displacement
   nod_pt->set_value(5,0.0);
  }

 // Vertical Symmetry boundary (r=0 and z<0)
 {
  unsigned ibound=1;
  {
   unsigned num_nod= Solid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin radial displacement (u_0 (real) and u_3 (imag))
     nod_pt->pin(0);
     nod_pt->set_value(0,0.0);
     nod_pt->pin(3);
     nod_pt->set_value(3,0.0);
     
     // Pin azimuthal displacement (u_2 (real) and u_5 (imag))
     nod_pt->pin(2);
     nod_pt->set_value(2,0.0);
     nod_pt->pin(5);
     nod_pt->set_value(5,0.0);
    }
  }
 }
 

 // Vertical Symmetry boundary (r=0 and z>0)
 {
  unsigned ibound=3;
  {
   unsigned num_nod= Solid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin radial displacement (u_0 (real) and u_3 (imag))
     nod_pt->pin(0);
     nod_pt->set_value(0,0.0);
     nod_pt->pin(3);
     nod_pt->set_value(3,0.0);
     
     // Pin azimuthal displacement (u_2 (real) and u_5 (imag))
     nod_pt->pin(2);
     nod_pt->set_value(2,0.0);
     nod_pt->pin(5);
     nod_pt->set_value(5,0.0);
    }
  }
 } // done sym bc

 // Setup fluid-structure interaction
 //----------------------------------
 setup_interaction();

 // Open trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",Global_Parameters::Directory.c_str());
 Trace_file.open(filename);
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}//end_of_constructor



//============start_of_create_outer_bc_elements===========================
/// Create BC elements on outer boundary
//========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedSphereProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
create_helmholtz_DtN_elements()
{
 // Outer boundary is boundary 2:
 unsigned b=2;

 // Loop over the bulk elements adjacent to boundary b?
 unsigned n_element = Helmholtz_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   HELMHOLTZ_ELEMENT* bulk_elem_pt = dynamic_cast<HELMHOLTZ_ELEMENT*>(
    Helmholtz_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = Helmholtz_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding DtN element
   FourierDecomposedHelmholtzDtNBoundaryElement<HELMHOLTZ_ELEMENT>* 
    flux_element_pt = new 
    FourierDecomposedHelmholtzDtNBoundaryElement<HELMHOLTZ_ELEMENT>
    (bulk_elem_pt,face_index);
   
   //Add the flux boundary element to the  helmholtz_DtN_mesh
   Helmholtz_DtN_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set pointer to the mesh that contains all the boundary condition
   // elements on this boundary
   flux_element_pt->set_outer_boundary_mesh_pt(Helmholtz_DtN_mesh_pt);
  }
 
} // end_of_create_outer_bc_elements





//=====================start_of_setup_interaction======================
/// Setup interaction between two fields
//========================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedSphereProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
setup_interaction()
{

 // Setup Helmholtz "pressure" load on traction elements
 unsigned boundary_in_helmholtz_mesh=0;

 // Doc boundary coordinate for Helmholtz
 ofstream the_file;
 the_file.open("boundary_coordinate_hh.dat");
 Helmholtz_mesh_pt->Mesh::template doc_boundary_coordinates<HELMHOLTZ_ELEMENT>
  (boundary_in_helmholtz_mesh, the_file);
 the_file.close();

 // Setup interaction
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <HELMHOLTZ_ELEMENT,2>
  (this,boundary_in_helmholtz_mesh,Helmholtz_mesh_pt,FSI_traction_mesh_pt);

 // Setup Helmholtz flux from normal displacement interaction
 unsigned boundary_in_solid_mesh=2;

 // Doc boundary coordinate for solid mesh
 the_file.open("boundary_coordinate_solid.dat");
 Solid_mesh_pt->Mesh::template doc_boundary_coordinates<ELASTICITY_ELEMENT>
  (boundary_in_solid_mesh, the_file);
 the_file.close();

 // Setup interaction
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <ELASTICITY_ELEMENT,2>(
   this,boundary_in_solid_mesh,Solid_mesh_pt,Helmholtz_fsi_flux_mesh_pt);

}// end_of_setup_interaction





//============start_of_create_fsi_traction_elements======================
/// Create fsi traction elements 
//=======================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedSphereProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
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
   FourierDecomposedTimeHarmonicLinElastLoadedByHelmholtzPressureBCElement
    <ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>* el_pt=
    new FourierDecomposedTimeHarmonicLinElastLoadedByHelmholtzPressureBCElement
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
 
} // end_of_create_fsi_traction_elements



//============start_of_create_helmholtz_fsi_flux_elements================
/// Create Helmholtz fsi flux elements 
//=======================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedSphereProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
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
   FourierDecomposedHelmholtzFluxFromNormalDisplacementBCElement
    <HELMHOLTZ_ELEMENT,ELASTICITY_ELEMENT>* el_pt=
    new FourierDecomposedHelmholtzFluxFromNormalDisplacementBCElement
    <HELMHOLTZ_ELEMENT,ELASTICITY_ELEMENT>(bulk_elem_pt,
                                           face_index);
   
   // Add to mesh
   Helmholtz_fsi_flux_mesh_pt->add_element_pt(el_pt);
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   el_pt->set_boundary_number_in_bulk_mesh(b); 
  }  
  
} // end_of_create_helmholtz_fsi_flux_elements



//==============start_of_doc_solution===============================
/// Doc the solution
//==================================================================
template<class ELASTICITY_ELEMENT, class HELMHOLTZ_ELEMENT>
void CoatedSphereProblem<ELASTICITY_ELEMENT,HELMHOLTZ_ELEMENT>::
doc_solution(DocInfo& doc_info)
{

 // Doc parameters
 oomph_info << "Writing result for step " << doc_info.number() 
            << ". Parameters: "<< std::endl;
 oomph_info << "Fourier mode number : N = "
            << Global_Parameters::Fourier_wavenumber << std::endl;
 oomph_info << "FSI parameter : Q = " << Global_Parameters::Q << std::endl;
 oomph_info << "Fluid outer radius : R = " << Global_Parameters::Outer_radius
            << std::endl;
 oomph_info << "Fluid wavenumber : k^2 = " << Global_Parameters::K_squared
            << std::endl;
 oomph_info << "Solid wavenumber : Omega_sq = " << Global_Parameters::Omega_sq 
            << std::endl << std::endl; 


 ofstream some_file,some_file2;
 char filename[100];

 // Number of plot points
 unsigned n_plot=5; 

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Accumulate contribution from elements
 double power=0.0;
 unsigned nn_element=Helmholtz_DtN_mesh_pt->nelement(); 
 for(unsigned e=0;e<nn_element;e++)
  {
   FourierDecomposedHelmholtzBCElementBase<HELMHOLTZ_ELEMENT> *el_pt = 
    dynamic_cast<FourierDecomposedHelmholtzBCElementBase<HELMHOLTZ_ELEMENT>*>(
     Helmholtz_DtN_mesh_pt->element_pt(e)); 
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();
 oomph_info << "Radiated power: " << power << std::endl;

 // Output displacement field
 //--------------------------
 sprintf(filename,"%s/elast_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,n_plot);
 some_file.close();

 // Output Helmholtz
 //-----------------
 sprintf(filename,"%s/helmholtz_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Helmholtz_mesh_pt->output(some_file,n_plot);
 some_file.close();


 // Output fsi traction elements
 //----------------------------- 
 sprintf(filename,"%s/fsi_traction_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 FSI_traction_mesh_pt->output(some_file,n_plot);
 some_file.close();


 // Output Helmholtz fsi flux elements
 //----------------------------------- 
 sprintf(filename,"%s/fsi_flux_bc_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Helmholtz_fsi_flux_mesh_pt->output(some_file,n_plot);
 some_file.close();

 // Write trace file
 Trace_file << Global_Parameters::Q << " " 
            << Global_Parameters::K_squared << " "
            << Global_Parameters::Density_ratio << " "
            << Global_Parameters::Omega_sq.real() << " "
            << power << " " 
            << std::endl;
   
 // Bump up counter
 doc_info.number()++;

} //end_of_doc_solution



//=======start_of_main==================================================
/// Driver for coated sphere loaded by lineared fluid loading
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

 // Parameter for the Helmholtz equation
 CommandLineArgs::specify_command_line_flag("--k_squared",
                                            &Global_Parameters::K_squared);
 
 // Initial value of Q 
 CommandLineArgs::specify_command_line_flag("--q_initial", 
                                            &Global_Parameters::Q);
 
 // Number of steps in parameter study
 unsigned nstep=2;
 CommandLineArgs::specify_command_line_flag("--nstep",&nstep);
 
 // Increment in FSI parameter in parameter study
 double q_increment=5.0;
 CommandLineArgs::specify_command_line_flag("--q_increment",&q_increment);
 
 
 // Wavenumber "zenith"-variation of imposed displacement of coating
 // on inner boundary 
 CommandLineArgs::specify_command_line_flag("--M",
                                            &Global_Parameters::M);
 
 // Multiplier for number of elements
 CommandLineArgs::specify_command_line_flag("--el_multiplier",
                                            &Global_Parameters::El_multiplier);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Update dependent parameters values
 Global_Parameters::update_parameter_values();

 // Set up doc info
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory(Global_Parameters::Directory);
 
 // Set up the problem
 CoatedSphereProblem<QTimeHarmonicFourierDecomposedLinearElasticityElement<3>,
                     QFourierDecomposedHelmholtzElement<3> > problem;

 //Parameter incrementation
 for(unsigned i=0;i<nstep;i++)
  {
   // Solve the problem with Newton's method
   problem.newton_solve();

   // Doc solution
   problem.doc_solution(doc_info);

   // Increment FSI parameter
   Global_Parameters::Q+=q_increment;
   Global_Parameters::update_parameter_values();
  }

} //end_of_main









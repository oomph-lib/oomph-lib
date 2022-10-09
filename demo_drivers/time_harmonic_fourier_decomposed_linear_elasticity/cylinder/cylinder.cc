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
// Driver 

// The oomphlib headers
#include "generic.h"
#include "time_harmonic_fourier_decomposed_linear_elasticity.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Define Poisson's ratio Nu
 std::complex<double> Nu(0.3,0.05);

 /// Define the non-dimensional Young's modulus
 std::complex<double> E(1.0,0.01);

 // Lame parameters
 std::complex<double> lambda = E*Nu/(1.0+Nu)/(1.0-2.0*Nu);
 std::complex<double> mu = E/2.0/(1.0+Nu);

 /// Define Fourier wavenumber
 int Fourier_wavenumber = 3;

 /// Define the non-dimensional square angular frequency of 
 /// time-harmonic motion
 std::complex<double> Omega_sq (10.0,5.0);

 /// Length of domain in r direction
 double Lr = 1.0;

 /// Length of domain in z-direction
 double Lz = 2.0;

 // Set up min & max (r,z) coordinates
 double rmin = 0.1;
 double zmin = 0.3;
 double rmax = rmin+Lr;
 double zmax = zmin+Lz;
 
 /// Define the imaginary unit
 const std::complex<double> I(0.0,1.0);

 /// The traction function at r=rmin: (t_r, t_z, t_theta)
 void boundary_traction(const Vector<double> &x,
                      const Vector<double> &n,
                      Vector<std::complex<double> > &result)
 {
  result[0] = -6.0*pow(x[0],2)*mu*cos(x[1])-
   lambda*(I*double(Fourier_wavenumber)*pow(x[0],2)*pow(x[1],3)+
           (4.0*pow(x[0],2)+pow(x[0],3))*cos(x[1]));
  result[1] = -mu*(3.0*pow(x[0],2)-pow(x[0],3))*sin(x[1]);
  result[2] = -mu*pow(x[0],2)*(2*pow(x[1],3)+I*double(Fourier_wavenumber)*
                               cos(x[1]));
 }
 
 
 /// The body force function; returns vector of complex doubles
 /// in the order (b_r, b_z, b_theta)
 void body_force(const Vector<double> &x,
                 Vector<std::complex<double> > &result)
 {
  result[0] = 
   x[0]*(-2.0*I*lambda*double(Fourier_wavenumber)*pow(x[1],3)-cos(x[1])*
         (lambda*(8.0+3.0*x[0])-
          mu*(pow(double(Fourier_wavenumber),2)
              -16.0+x[0]*(x[0]-3.0))+pow(x[0],2)*Omega_sq));
  result[1] = 
   x[0]*sin(x[1])*(mu*(pow(double(Fourier_wavenumber),2)-9.0)+
                   4.0*x[0]*(lambda+mu)+pow(x[0],2)*
                   (lambda+2.0*mu-Omega_sq))-
   3.0*I*double(Fourier_wavenumber)*pow(x[0],2)*pow(x[1],2)*(lambda+mu);
  result[2] = 
   -x[0]*(8.0*mu*pow(x[1],3)-pow(double(Fourier_wavenumber),2)*pow(x[1],3)*
          (lambda+2.0*mu)+pow(x[0],2)*(pow(x[1],3)*Omega_sq+6.0*mu*x[1])+
          I*cos(x[1])*double(Fourier_wavenumber)*
          (lambda*(4.0+x[0])+mu*(6.0+x[0])));
 }
 
 /// The exact solution in a flat-packed vector:
 // 0: u_r[real], 1: u_z[real],..., 5: u_theta[imag]
 void exact_solution(const Vector<double> &x,
                     Vector<double> &u)
 {
  u[0] = pow(x[0],3)*cos(x[1]);
  u[1] = pow(x[0],3)*sin(x[1]);
  u[2] = pow(x[0],3)*pow(x[1],3);
  u[3] = 0.0;
  u[4] = 0.0;
  u[5] = 0.0;
 }

} // end_of_namespace


//===start_of_problem_class=============================================
/// Class to validate time harmonic linear elasticity (Fourier 
/// decomposed)
//======================================================================
template<class ELEMENT>
class FourierDecomposedTimeHarmonicLinearElasticityProblem : public Problem
{
public:

 /// Constructor: Pass number of elements in r and z directions 
 /// and boundary locations
 FourierDecomposedTimeHarmonicLinearElasticityProblem(
 const unsigned &nr, const unsigned &nz,
 const double &rmin, const double& rmax,
 const double &zmin, const double& zmax);

 
 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:
 
 /// Allocate traction elements on the bottom surface
 void assign_traction_elements();
 
 /// Pointer to the bulk mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the mesh of traction elements
 Mesh* Surface_mesh_pt;
}; // end_of_problem_class


//===start_of_constructor=============================================
/// Problem constructor: Pass number of elements in coordinate
/// directions and size of domain.
//====================================================================
template<class ELEMENT>
FourierDecomposedTimeHarmonicLinearElasticityProblem<ELEMENT>::
FourierDecomposedTimeHarmonicLinearElasticityProblem
(const unsigned &nr, const unsigned &nz,
 const double &rmin, const double& rmax,
 const double &zmin, const double& zmax)
{
 //Now create the mesh
 Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(nr,nz,rmin,rmax,zmin,zmax);

 //Create the surface mesh of traction elements
 Surface_mesh_pt=new Mesh;
 assign_traction_elements();
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin & set the ones that have Dirichlet 
 // conditions here
 
 // storage for nodal position
 Vector<double> x(2);

 // Storage for prescribed displacements
 Vector<double> u(6);

 // Now set displacements on boundaries 0 (z=zmin),
 //------------------------------------------------
 // 1 (r=rmax) and 2 (z=zmax)
 //--------------------------
 for (unsigned ibound=0;ibound<=2;ibound++)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++) 
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     
     // get r and z coordinates
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     
     // Pinned in r, z and theta
     nod_pt->pin(0);nod_pt->pin(1);nod_pt->pin(2);
     nod_pt->pin(3);nod_pt->pin(4);nod_pt->pin(5);
     
     // Compute the value of the exact solution at the nodal point
     Vector<double> u(6);
     Global_Parameters::exact_solution(x,u);
     
     // Set the displacements
     nod_pt->set_value(0,u[0]);
     nod_pt->set_value(1,u[1]);
     nod_pt->set_value(2,u[2]);
     nod_pt->set_value(3,u[3]);
     nod_pt->set_value(4,u[4]);
     nod_pt->set_value(5,u[5]);
    }
  } // end_of_loop_over_boundary_nodes


 // Complete the problem setup to make the elements fully functional

 // Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Cast to a bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the body force
   el_pt->body_force_fct_pt() = &Global_Parameters::body_force;

   // Set the pointer to Poisson's ratio
   el_pt->nu_pt() = &Global_Parameters::Nu;

   // Set the pointer to Fourier wavenumber
   el_pt->fourier_wavenumber_pt() = &Global_Parameters::Fourier_wavenumber;

   // Set the pointer to non-dim Young's modulus
   el_pt->youngs_modulus_pt() = &Global_Parameters::E;

   // Set the pointer to square of the angular frequency
   el_pt->omega_sq_pt() = &Global_Parameters::Omega_sq;

  }// end loop over elements

 // Loop over the traction elements
 unsigned n_traction =  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_traction;e++)
  {
   // Cast to a surface element
   TimeHarmonicFourierDecomposedLinearElasticityTractionElement<ELEMENT>*
    el_pt = 
    dynamic_cast<TimeHarmonicFourierDecomposedLinearElasticityTractionElement
    <ELEMENT>* >(Surface_mesh_pt->element_pt(e));
   
   // Set the applied traction
   el_pt->traction_fct_pt() = &Global_Parameters::boundary_traction;
   
  }// end loop over traction elements
 
 // Add the submeshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Now build the global mesh
 build_global_mesh();

 // Assign equation numbers
 cout << assign_eqn_numbers() << " equations assigned" << std::endl; 

} // end of constructor


//===start_of_traction===============================================
/// Make traction elements along the boundary r=rmin
//===================================================================
template<class ELEMENT>
void FourierDecomposedTimeHarmonicLinearElasticityProblem<ELEMENT>::
assign_traction_elements()
{
 unsigned bound, n_neigh;

 // How many bulk elements are next to boundary 3
 bound=3;
 n_neigh = Bulk_mesh_pt->nboundary_element(bound); 

 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   FiniteElement *traction_element_pt 
    = new TimeHarmonicFourierDecomposedLinearElasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
 
   // Add to mesh
   Surface_mesh_pt->add_element_pt(traction_element_pt);
  }

} // end of assign_traction_elements


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FourierDecomposedTimeHarmonicLinearElasticityProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts=5; 
 
 // Output solution 
 sprintf(filename,"%s/soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 sprintf(filename,"%s/exact_soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,
                          Global_Parameters::exact_solution);
 some_file.close();

 // Doc error
 double error=0.0;
 double norm=0.0;
 sprintf(filename,"%s/error.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             Global_Parameters::exact_solution, 
                             error,norm);
 some_file.close();

 // Doc error norm:
 cout << "\nNorm of error:    " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end_of_doc_solution   


//===start_of_main======================================================
/// Driver code 
//======================================================================
int main(int argc, char* argv[]) 
{
 // Number of elements in r-direction
 unsigned nr=5;
 
 // Number of elements in z-direction (for (approximately) square elements)
 unsigned nz=unsigned(double(nr)*Global_Parameters::Lz/Global_Parameters::Lr);

 // Set up doc info
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");

 // Set up problem
 FourierDecomposedTimeHarmonicLinearElasticityProblem
  <QTimeHarmonicFourierDecomposedLinearElasticityElement<3> > 
  problem(nr,nz,Global_Parameters::rmin,Global_Parameters::rmax,
          Global_Parameters::zmin,Global_Parameters::zmax);
 
 // Solve
 problem.newton_solve();
 
 // Output the solution
 problem.doc_solution(doc_info);
  
} // end_of_main

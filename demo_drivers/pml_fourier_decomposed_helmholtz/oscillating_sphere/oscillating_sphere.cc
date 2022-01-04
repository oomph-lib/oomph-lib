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
//Driver for pml Fourier-decomposed  Helmholtz problem

#include <complex>
#include <cmath>

//Generic routines
#include "generic.h"

// The Helmholtz equations
#include "pml_fourier_decomposed_helmholtz.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//===== start_of_namespace=============================================
/// Namespace for the Fourier decomposed Helmholtz problem parameters
//=====================================================================
namespace ProblemParameters
{
 /// Output directory
 string Directory="RESLT";

 /// Frequency
 double K_squared = 10.0;

 /// Default physical PML thickness
 double PML_thickness=4.0;

 /// Default number of elements within PMLs
 unsigned Nel_pml=15;

 /// Target area for initial mesh
 double Element_area = 0.1;

 /// The default Fourier wave number
 int N_fourier=0;

 /// Number of terms in the exact solution
 unsigned N_terms=6;

 /// Coefficients in the exact solution
 Vector<double> Coeff(N_terms,1.0);

 /// Imaginary unit
 std::complex<double> I(0.0,1.0);

 /// Exact solution as a Vector of size 2, containing real and imag parts
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  double K = sqrt(K_squared);

  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);

  double theta;
  theta=atan2(x[0],x[1]);

  // Argument for Bessel/Hankel functions
  double kr= K*R;

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
  double K = sqrt(K_squared);

  // Initialise flux
  flux=std::complex<double>(0.0,0.0);

  // Switch to spherical coordinates
  double R=sqrt(x[0]*x[0]+x[1]*x[1]);

  double theta;
  theta=atan2(x[0],x[1]);

  // Argument for Bessel/Hankel functions
  double kr=K*R;

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
     ( K*(djv[i]+I*dyv[i]) - (0.5*(jv[i]+I*yv[i])/R) );
   }
 } // end of exact_normal_derivative


 /// Radial position of point source
 double R_source = 2.0;

 /// Axial position of point source
 double Z_source = 2.0;

 /// Point source magnitude (Complex)
 std::complex<double> Magnitude(100.0,100.0);

} // end of namespace


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////


namespace oomph
{

//========= start_of_point_source_wrapper==============================
/// Class to impose point source to (wrapped) Helmholtz element
//=====================================================================
template<class ELEMENT>
class PMLHelmholtzPointSourceElement : public virtual ELEMENT
{

public:

 /// Constructor
 PMLHelmholtzPointSourceElement()
  {
   // Initialise
   Point_source_magnitude=std::complex<double>(0.0,0.0);
  }

 /// Destructor (empty)
 ~PMLHelmholtzPointSourceElement(){}

 /// Set local coordinate and magnitude of point source
 void setup(const Vector<double>& s_point_source,
            const std::complex<double>& magnitude)
  {
   S_point_source=s_point_source;
   Point_source_magnitude=magnitude;
  }


 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   ELEMENT::fill_in_generic_residual_contribution_pml_fourier_decomposed_helmholtz(
    residuals,GeneralisedElement::Dummy_matrix,0);

   // Add point source contribution
   fill_in_point_source_contribution_to_residuals(residuals);
  }


 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   ELEMENT::fill_in_generic_residual_contribution_pml_fourier_decomposed_helmholtz(residuals,
                                                            jacobian,1);

   // Add point source contribution
   fill_in_point_source_contribution_to_residuals(residuals);
  }


private:


 /// Add the point source contribution to the residual vector
 void fill_in_point_source_contribution_to_residuals(Vector<double> &residuals)
  {
   // No further action
   if (S_point_source.size()==0) return;

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();

   //Set up memory for the shape/test functions
   Shape psi(n_node);

   //Integers to store the local equation and unknown numbers
   int local_eqn_real=0;
   int local_eqn_imag=0;

   // Get shape/test fcts
   this->shape(S_point_source,psi);

   // Assemble residuals
   //--------------------------------

   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     // first, compute the real part contribution
     //-------------------------------------------

     //Get the local equation
     local_eqn_real =
      this->nodal_local_eqn
      (l,this->u_index_pml_fourier_decomposed_helmholtz().real());

     /*IF it's not a boundary condition*/
     if(local_eqn_real >= 0)
      {
       residuals[local_eqn_real] -= 2.0*Point_source_magnitude.real()*psi(l);
      }

     // Second, compute the imaginary part contribution
     //------------------------------------------------

     //Get the local equation
     local_eqn_imag =
      this->nodal_local_eqn
      (l,this->u_index_pml_fourier_decomposed_helmholtz().imag());

     /*IF it's not a boundary condition*/
     if(local_eqn_imag >= 0)
      {
       // Add body force/source term and Helmholtz bit
       residuals[local_eqn_imag] -= 2.0*Point_source_magnitude.imag()*psi(l);
      }
    }
  }

 /// Local coordinates of point at which point source is applied
 Vector<double> S_point_source;

 /// Magnitude of point source (complex!)
 std::complex<double> Point_source_magnitude;

 };




//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<class ELEMENT>
 class FaceGeometry<PMLHelmholtzPointSourceElement<ELEMENT> >
  : public virtual FaceGeometry<ELEMENT>
 {
 public:
  FaceGeometry() : FaceGeometry<ELEMENT>() {}
 };


//=======================================================================
/// Face geometry of the Face Geometry for element is the same as
/// that for the underlying wrapped element
//=======================================================================
 template<class ELEMENT>
 class FaceGeometry<FaceGeometry<PMLHelmholtzPointSourceElement<ELEMENT> > >
  : public virtual FaceGeometry<FaceGeometry<ELEMENT> >
 {
 public:
  FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
 };


//=======================================================================
/// Policy class defining the elements to be used in the actual
/// PML layers.
//=======================================================================
 template<class ELEMENT>
class PMLLayerElement<PMLHelmholtzPointSourceElement<ELEMENT> > :
 public virtual PMLLayerElement<ELEMENT>
{

  public:

 /// Constructor: Call the constructor for the
 /// appropriate Element
 PMLLayerElement() : PMLLayerElement<ELEMENT>()
  {}

};



//=======================================================================
/// Policy class defining the elements to be used in the actual
/// PML layers.
//=======================================================================
 template<class ELEMENT>
class PMLLayerElement<
  ProjectablePMLFourierDecomposedHelmholtzElement<
  PMLHelmholtzPointSourceElement<ELEMENT> > >:
 public virtual PMLLayerElement<ELEMENT>
{

  public:

 /// Constructor: Call the constructor for the
 /// appropriate Element
 PMLLayerElement() : PMLLayerElement<ELEMENT>()
  {}

};

}




/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class
//=====================================================================
template<class ELEMENT>
class PMLFourierDecomposedHelmholtzProblem : public Problem
{

public:

 /// Constructor
 PMLFourierDecomposedHelmholtzProblem();

 /// Destructor (empty)
 ~PMLFourierDecomposedHelmholtzProblem(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// Create PML meshes
 void create_pml_meshes();

 /// Create mesh of face elements that monitor the radiated power
 void create_power_monitor_mesh();

 #ifdef ADAPTIVE

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();

 #endif


  // Apply boundary condtions for odd Fourier wavenumber
 void complete_problem_setup();

private:

 /// Create flux elements on inner boundary
 void create_flux_elements_on_inner_boundary();

 /// Delete boundary face elements and wipe the surface mesh
 void delete_face_elements(Mesh* const & boundary_mesh_pt)
  {
   // Loop over the surface elements
   unsigned n_element = boundary_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete  boundary_mesh_pt->element_pt(e);
    }

   // Wipe the mesh
   boundary_mesh_pt->flush_element_and_node_storage();
  }

 // Apply boundary condtions for odd Fourier wavenumber
 void apply_zero_dirichlet_boundary_conditions();

 /// Pointer to mesh that stores the power monitor elements
 Mesh* Power_monitor_mesh_pt;


#ifdef ADAPTIVE

 // Create point source element (only used in adaptive run)
 void setup_point_source();

 /// Pointer to the refineable "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Mesh as geometric object representation of bulk mesh
 MeshAsGeomObject* Mesh_as_geom_obj_pt;

#else

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif

 /// Mesh of FaceElements that apply the flux bc on the inner boundary
 Mesh* Helmholtz_inner_boundary_mesh_pt;

 /// Pointer to the right PML mesh
 Mesh* PML_right_mesh_pt;

 /// Pointer to the top PML mesh
 Mesh* PML_top_mesh_pt;

 /// Pointer to the bottom PML mesh
 Mesh* PML_bottom_mesh_pt;

 /// Pointer to the top right corner PML mesh
 Mesh* PML_top_right_mesh_pt;

 /// Pointer to the bottom right corner PML mesh
 Mesh* PML_bottom_right_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class


//===================start_of_create_power_monitor_mesh===================
/// Create BC elements on outer boundary
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
create_power_monitor_mesh()
{
 // Loop over outer boundaries
 for (unsigned b=1;b<4;b++)
  {
   // Loop over the bulk elements adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));

     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);

     // Build the corresponding element
     PMLFourierDecomposedHelmholtzPowerMonitorElement<ELEMENT>*
      flux_element_pt = new
      PMLFourierDecomposedHelmholtzPowerMonitorElement<ELEMENT>
      (bulk_elem_pt,face_index);

     //Add the flux boundary element
     Power_monitor_mesh_pt->add_element_pt(flux_element_pt);
    }
  }
} // end of create_power_monitor_mesh



#ifdef ADAPTIVE

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
actions_before_adapt()
{
 // Before adapting the added PML meshes must be removed
 // as they are not refineable and are to be rebuilt from the
 // newly refined triangular mesh
 delete PML_right_mesh_pt;
 PML_right_mesh_pt=0;
 delete PML_top_mesh_pt;
 PML_top_mesh_pt=0;
 delete PML_bottom_mesh_pt;
 PML_bottom_mesh_pt=0;
 delete PML_top_right_mesh_pt;
 PML_top_right_mesh_pt=0;
 delete PML_bottom_right_mesh_pt;
 PML_bottom_right_mesh_pt=0;

 // Wipe the power monitor elements
 delete_face_elements(Power_monitor_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 // but first flush all its submeshes
 flush_sub_meshes();

 // Then add the triangular mesh back
 add_sub_mesh(Bulk_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
actions_after_adapt()
{

 // Build PML meshes  and add them to the global mesh
 create_pml_meshes();

 // Re-attach the power monitor elements
 create_power_monitor_mesh();

 // Build the entire mesh from its submeshes
 rebuild_global_mesh();

 // Complete the build of all elements
 complete_problem_setup();

 // Setup point source
 setup_point_source();

}// end of actions_after_adapt

#endif



//=================start_of_complete_problem_setup==================
// Complete the build of all elements so that they are fully
// functional
//==================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
complete_problem_setup()
{
 // Complete the build of all elements so they are fully functional
 unsigned n_element = this->mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   PMLFourierDecomposedHelmholtzEquations *el_pt = dynamic_cast<
    PMLFourierDecomposedHelmholtzEquations*>(
     mesh_pt()->element_pt(i));

   if (!(el_pt==0))
    {
     //Set the frequency pointer
     el_pt->k_squared_pt()=&ProblemParameters::K_squared;

     // Set pointer to Fourier wave number
     el_pt->pml_fourier_wavenumber_pt()=&ProblemParameters::N_fourier;
    }
  }

 // If the Fourier wavenumber is odd, then apply zero dirichlet boundary
 // conditions on the two straight boundaries on the symmetry line.
 if (ProblemParameters::N_fourier % 2 == 1)
  {
   cout
    << "Zero Dirichlet boundary condition has been applied on symmetry line\n";
   cout << "due to an odd Fourier wavenumber\n" << std::endl;
   apply_zero_dirichlet_boundary_conditions();
  }

} // end of complete_problem_setup


#ifdef ADAPTIVE

//==================start_of_setup_point_source===========================
/// Set point source
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
setup_point_source()
{
 // Create mesh as geometric object
 delete Mesh_as_geom_obj_pt;
 Mesh_as_geom_obj_pt= new MeshAsGeomObject(Bulk_mesh_pt);

 // Position of point source
 Vector<double> x_point_source;
 x_point_source.resize(2);
 x_point_source[0]=ProblemParameters::R_source;
 x_point_source[1]=ProblemParameters::Z_source;

 GeomObject* sub_geom_object_pt=0;
 Vector<double> s_point_source(2);
 Mesh_as_geom_obj_pt->locate_zeta(x_point_source,sub_geom_object_pt,
                                  s_point_source);

 // Set point force
 if(x_point_source[0]==0.0)
  {
   if(ProblemParameters::N_fourier>0)
    {
     // if source on z axis, only contribution to residual comes
     // from Fourier wavenumber zero
     ProblemParameters::Magnitude=complex<double>(0.0,0.0);
    }
  }

 // Set point force
 dynamic_cast<ELEMENT*>(sub_geom_object_pt)->
  setup(s_point_source,ProblemParameters::Magnitude);

}


#endif

//=========start_of_apply_zero_dirichlet_boundary_conditions========
// Apply extra bounday conditions if given an odd Fourier wavenumber
//==================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
apply_zero_dirichlet_boundary_conditions()
{
 // Apply zero dirichlet conditions on the bottom straight boundary
 // and the top straight boundary located on the symmetry line.

 // Bottom straight boundary on symmetry line:
 {
  //Boundary id
  unsigned b=0;

  // How many nodes are there?
  unsigned n_node=Bulk_mesh_pt->nboundary_node(b);
  for (unsigned n=0;n<n_node;n++)
   {
    // Get the node
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

    // Pin the node
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Set the node's value
    nod_pt->set_value(0, 0.0);
    nod_pt->set_value(1, 0.0);
   }
 }

// Top straight boundary on symmetry line:
 {
  //Boundary id
  unsigned b=4;

  // How many nodes are there?
  unsigned n_node=Bulk_mesh_pt->nboundary_node(b);
  for (unsigned n=0;n<n_node;n++)
   {
    // Get the node
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

    // Pin the node
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Set the node's value
    nod_pt->set_value(0, 0.0);
    nod_pt->set_value(1, 0.0);
   }
 }


} // end of apply_zero_dirichlet_boundary_conditions


//=======start_of_constructor=============================================
/// Constructor for Pml Fourier-decomposed Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
PMLFourierDecomposedHelmholtzProblem()
{
  string trace_file_location = ProblemParameters::Directory + "/trace.dat";

 // Open trace file
 Trace_file.open(trace_file_location.c_str());

 /// Setup "bulk" mesh

 // Create the circle that represents the inner boundary
 double x_c=0.0;
 double y_c=0.0;
 double r_min=1.0;
 double r_max=3.0;
 Circle* inner_circle_pt=new Circle(x_c,y_c,r_min);


 // Edges/boundary segments making up outer boundary
 //-------------------------------------------------
 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(6);

 // All poly boundaries are defined by two vertices
 Vector<Vector<double> > boundary_vertices(2);


 // Bottom straight boundary on symmetry line
 //------------------------------------------
 boundary_vertices[0].resize(2);
 boundary_vertices[0][0]=0.0;
 boundary_vertices[0][1]=-r_min;
 boundary_vertices[1].resize(2);
 boundary_vertices[1][0]=0.0;
 boundary_vertices[1][1]=-r_max;

 unsigned boundary_id=0;
 outer_boundary_line_pt[0]=
  new TriangleMeshPolyLine(boundary_vertices,boundary_id);


 // Bottom boundary of bulk mesh
 //-----------------------------
 boundary_vertices[0][0]=0.0;
 boundary_vertices[0][1]=-r_max;
 boundary_vertices[1][0]=r_max;;
 boundary_vertices[1][1]=-r_max;

 boundary_id=1;
 outer_boundary_line_pt[1]=
  new TriangleMeshPolyLine(boundary_vertices,boundary_id);


 // Right boundary of bulk mesh
 //----------------------------
 boundary_vertices[0][0]=r_max;
 boundary_vertices[0][1]=-r_max;
 boundary_vertices[1][0]=r_max;;
 boundary_vertices[1][1]=r_max;

 boundary_id=2;
 outer_boundary_line_pt[2]=
  new TriangleMeshPolyLine(boundary_vertices,boundary_id);


 // Top boundary of bulk mesh
 //--------------------------
 boundary_vertices[0][0]=r_max;
 boundary_vertices[0][1]=r_max;
 boundary_vertices[1][0]=0.0;
 boundary_vertices[1][1]=r_max;

 boundary_id=3;
 outer_boundary_line_pt[3]=
  new TriangleMeshPolyLine(boundary_vertices,boundary_id);

// Top straight boundary on symmetry line
 //---------------------------------------
 boundary_vertices[0][0]=0.0;
 boundary_vertices[0][1]=r_max;
 boundary_vertices[1][0]=0.0;
 boundary_vertices[1][1]=r_min;

 boundary_id=4;
 outer_boundary_line_pt[4]=
  new TriangleMeshPolyLine(boundary_vertices,boundary_id);


 // Inner circular boundary:
 //-------------------------

 // Number of segments used for representing the curvilinear boundary
 unsigned n_segments = 20;

 // The intrinsic coordinates for the beginning and end of the curve
 double s_start =  0.5*MathematicalConstants::Pi;
 double s_end   =  -0.5*MathematicalConstants::Pi;

 boundary_id = 5;
 outer_boundary_line_pt[5]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);


 // Create closed curve that defines outer boundary
 //------------------------------------------------
 TriangleMeshClosedCurve *outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_boundary_line_pt);


 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters. The only parameter that needs to take is the
 // outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);


 // Specify maximum element area
 double element_area = ProblemParameters::Element_area;
 triangle_mesh_parameters.element_area() = element_area;

#ifdef ADAPTIVE

 // Build "bulk" mesh
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Add the bulk mesh to the problem
 add_sub_mesh(Bulk_mesh_pt);

 // Initialise mesh as geom object
 Mesh_as_geom_obj_pt=0;

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=0.00004;
 Bulk_mesh_pt->max_permitted_error()=0.0001;

 // Reduce wavenumber to make effect of singularity more prominent
 ProblemParameters::K_squared = 0.1*0.1;

#else

 // Create the bulk mesh
 Bulk_mesh_pt= new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Add the bulk mesh to the problem
 add_sub_mesh(Bulk_mesh_pt);

  // Create flux elements on inner boundary
 Helmholtz_inner_boundary_mesh_pt=new Mesh;
 create_flux_elements_on_inner_boundary();

 // ...and add the mesh to the problem
 add_sub_mesh(Helmholtz_inner_boundary_mesh_pt);

#endif


 // Attach the power monitor elements
 Power_monitor_mesh_pt=new Mesh;
 create_power_monitor_mesh();

 // Create the pml meshes
 create_pml_meshes();

 // Build the Problem's global mesh from its various sub-meshes
 build_global_mesh();

 // Complete the build of all elements
 complete_problem_setup();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;


 // Total radiated power
 double power=0.0;

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Accumulate contribution from elements
 unsigned nn_element=Power_monitor_mesh_pt->nelement();
 for(unsigned e=0;e<nn_element;e++)
  {
   PMLFourierDecomposedHelmholtzPowerMonitorElement<ELEMENT> *el_pt =
    dynamic_cast<PMLFourierDecomposedHelmholtzPowerMonitorElement
    <ELEMENT>*>(Power_monitor_mesh_pt->element_pt(e));
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();
 oomph_info << "Total radiated power: " << power << std::endl;


 // Output solution within the bulk mesh
 //-------------------------------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output solution within pml domains
 //-----------------------------------
 sprintf(filename,"%s/pml_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 PML_top_mesh_pt->output(some_file,npts);
 PML_right_mesh_pt->output(some_file,npts);
 PML_bottom_mesh_pt->output(some_file,npts);
 PML_top_right_mesh_pt->output(some_file,npts);
 PML_bottom_right_mesh_pt->output(some_file,npts);
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

 sprintf(filename,"%s/int_error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 // Doc L2 error and norm of solution
 some_file << "\nNorm of error   : " << sqrt(error) << std::endl;
 some_file << "Norm of solution: " << sqrt(norm) << std::endl <<std::endl;
 some_file << "Relative error: " << sqrt(error)/sqrt(norm) << std::endl;
 some_file.close();

 // Write norm and radiated power of solution to trace file
 Bulk_mesh_pt->compute_norm(norm);
 Trace_file  << norm << " "  << power << std::endl;

} // end of doc


//============start_of_create_flux_elements=================
/// Create flux elements on inner boundary
//==========================================================
template<class ELEMENT>
void  PMLFourierDecomposedHelmholtzProblem<ELEMENT>::
create_flux_elements_on_inner_boundary()
{
 // Apply flux bc on inner boundary (boundary 5)
 unsigned b=5;

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
   PMLFourierDecomposedHelmholtzFluxElement<ELEMENT>*
    flux_element_pt = new
    PMLFourierDecomposedHelmholtzFluxElement<ELEMENT>
    (bulk_elem_pt,face_index);

   //Add the prescribed incoming-flux element to the surface mesh
   Helmholtz_inner_boundary_mesh_pt->add_element_pt(flux_element_pt);

   // Set the pointer to the prescribed flux function
   flux_element_pt->flux_fct_pt() = &ProblemParameters::exact_minus_dudr;

  } //end of loop over bulk elements adjacent to boundary b

} // end of create flux elements on inner boundary



//============start_of_create_pml_meshes======================================
/// Create PML meshes and add them to the problem's sub-meshes
//============================================================================
template<class ELEMENT>
void PMLFourierDecomposedHelmholtzProblem<ELEMENT>::create_pml_meshes()
{
 // Bulk mesh bottom boundary id
 unsigned int bottom_boundary_id = 1;

 // Bulk mesh right boundary id
 unsigned int right_boundary_id = 2;

 // Bulk mesh top boundary id
 unsigned int top_boundary_id = 3;

 // PML width in elements for the right layer
 unsigned n_x_right_pml = ProblemParameters::Nel_pml;

 // PML width in elements for the top layer
 unsigned n_y_top_pml = ProblemParameters::Nel_pml;

 // PML width in elements for the bottom layer
 unsigned n_y_bottom_pml = ProblemParameters::Nel_pml;


 // Outer physical length of the PML layers
 // defaults to 4.0
 double width_x_right_pml  = ProblemParameters::PML_thickness;
 double width_y_top_pml    = ProblemParameters::PML_thickness;
 double width_y_bottom_pml = ProblemParameters::PML_thickness;

 // Build the PML meshes based on the new adapted triangular mesh
 PML_right_mesh_pt = TwoDimensionalPMLHelper::create_right_pml_mesh
  <PMLLayerElement<ELEMENT> >(Bulk_mesh_pt,
                              right_boundary_id,
                              n_x_right_pml,
                              width_x_right_pml);
 PML_top_mesh_pt   = TwoDimensionalPMLHelper::create_top_pml_mesh
  <PMLLayerElement<ELEMENT> >(Bulk_mesh_pt,
                              top_boundary_id,
                              n_y_top_pml,
                              width_y_top_pml);
 PML_bottom_mesh_pt= TwoDimensionalPMLHelper::create_bottom_pml_mesh
  <PMLLayerElement<ELEMENT> >(Bulk_mesh_pt,
                              bottom_boundary_id,
                              n_y_bottom_pml,
                              width_y_bottom_pml);

 // Add submeshes to the global mesh
 add_sub_mesh(PML_right_mesh_pt);
 add_sub_mesh(PML_top_mesh_pt);
 add_sub_mesh(PML_bottom_mesh_pt);

 // Rebuild corner PML meshes
 PML_top_right_mesh_pt    =
  TwoDimensionalPMLHelper::create_top_right_pml_mesh
  <PMLLayerElement<ELEMENT> >(PML_right_mesh_pt,
                              PML_top_mesh_pt,
                              Bulk_mesh_pt,
                              right_boundary_id);

 PML_bottom_right_mesh_pt =
  TwoDimensionalPMLHelper::create_bottom_right_pml_mesh
  <PMLLayerElement<ELEMENT> >(PML_right_mesh_pt,
                              PML_bottom_mesh_pt,
                              Bulk_mesh_pt,
                              right_boundary_id);

 // Add submeshes to the global mesh
 add_sub_mesh(PML_top_right_mesh_pt);
 add_sub_mesh(PML_bottom_right_mesh_pt);

} // end of create_pml_meshes


//===== start_of_main=====================================================
/// Driver code for Pml Fourier decomposed Helmholtz problem
//========================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Fourier wavenumber
 CommandLineArgs::specify_command_line_flag("--Fourier_wavenumber",
                                            &ProblemParameters::N_fourier);

 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &ProblemParameters::Directory);

 // PML thickness
 CommandLineArgs::specify_command_line_flag("--pml_thick",
                                            &ProblemParameters::PML_thickness);

 // Number of elements within PMLs
 CommandLineArgs::specify_command_line_flag("--npml_element",
                                            &ProblemParameters::Nel_pml);

 // Target Element size on first mesh generation
 CommandLineArgs::specify_command_line_flag("--element_area",
                                            &ProblemParameters::Element_area);
 // k squared (wavenumber squared)
 CommandLineArgs::specify_command_line_flag("--k_squared",
                                            &ProblemParameters::K_squared);

 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validate");


  // Demonstrate across a range of omega (or k) values with good mesh for a
  // visual test of accuracy (put in by Jonathan Deakin)
  CommandLineArgs::specify_command_line_flag("--demonstrate");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();


 // Create label for output
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory(ProblemParameters::Directory);

#ifdef ADAPTIVE

 // Create the problem with 2D projectable six-node elements from the
 // TPMLFourierDecomposedHelmholtzElement family,
 // allowing for the imposition of a point source (via another
 // templated wrapper)
 PMLFourierDecomposedHelmholtzProblem<
 ProjectablePMLFourierDecomposedHelmholtzElement<
 PMLHelmholtzPointSourceElement<
 TPMLFourierDecomposedHelmholtzElement<3> > > > problem;

#else

 // Create the problem with 2D six-node elements from the
 // TPMLFourierDecomposedHelmholtzElement family.
 PMLFourierDecomposedHelmholtzProblem
  <TPMLFourierDecomposedHelmholtzElement<3> >
  problem;

#endif

 // Step number
 doc_info.number()=ProblemParameters::N_fourier;

#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=4;

 // Validation run?
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   max_adapt=1;
  }

 // Solve the problem, allowing
 // up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);

#else


// Solve the problem with Newton's method
problem.newton_solve();



#endif

 //Output the solution
 problem.doc_solution(doc_info);


} //end of main

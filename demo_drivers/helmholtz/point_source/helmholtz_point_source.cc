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
// Helmholtz equation with point source
#include<fenv.h>

#include "math.h"
#include <complex>

//Generic routines
#include "generic.h"

// The Helmholtz equations
#include "helmholtz.h"

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
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{

 /// Square of the wavenumber
 double K_squared=10.0; 
  
 /// Flag to choose the Dirichlet to Neumann BC
 /// or ABC BC
 bool DtN_BC=false;

 /// Flag to choose wich order to use
 // in the ABCs BC: 1 for ABC 1st order...
 unsigned ABC_order=3;

 /// Radius of outer boundary (must be a circle!)
 double Outer_radius=1.5;

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
class HelmholtzPointSourceElement : public virtual ELEMENT
{

public:

 /// Constructor
 HelmholtzPointSourceElement()
  {
   // Initialise
   Point_source_magnitude=std::complex<double>(0.0,0.0);
  }
 
 /// Destructor (empty)
 ~HelmholtzPointSourceElement(){}
 
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
   ELEMENT::fill_in_generic_residual_contribution_helmholtz(
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
   ELEMENT::fill_in_generic_residual_contribution_helmholtz(residuals,
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
     local_eqn_real = this->nodal_local_eqn(l,this->u_index_helmholtz().real());
     
     /*IF it's not a boundary condition*/
     if(local_eqn_real >= 0)
      {
       residuals[local_eqn_real] += Point_source_magnitude.real()*psi(l);
      }     
     
     // Second, compute the imaginary part contribution 
     //------------------------------------------------
     
     //Get the local equation
     local_eqn_imag = this->nodal_local_eqn(l,this->u_index_helmholtz().imag());
     
     /*IF it's not a boundary condition*/
     if(local_eqn_imag >= 0)
      {
       // Add body force/source term and Helmholtz bit
       residuals[local_eqn_imag] += Point_source_magnitude.imag()*psi(l);
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
 class FaceGeometry<HelmholtzPointSourceElement<ELEMENT> > 
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
 class FaceGeometry<FaceGeometry<HelmholtzPointSourceElement<ELEMENT> > >
  : public virtual FaceGeometry<FaceGeometry<ELEMENT> >
 {
 public:
  FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
 };


}

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//========= start_of_problem_class=====================================
/// Problem class
//=====================================================================
template<class ELEMENT> 
class HelmholtzPointSourceProblem : public Problem
{

public:
 
 /// Constructor
 HelmholtzPointSourceProblem();
 
 /// Destructor (empty)
 ~HelmholtzPointSourceProblem(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){} 

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Recompute gamma integral before checking Newton residuals
 void actions_before_newton_convergence_check()
  {
   if (GlobalParameters::DtN_BC)
    {
     Helmholtz_outer_boundary_mesh_pt->setup_gamma();
    }
  }

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();
 
private:

 /// Create BC elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified survace Mesh 
 void create_outer_bc_elements(
  const unsigned &b, Mesh* const &bulk_mesh_pt,
  Mesh* const & helmholtz_outer_boundary_mesh_pt);
 
 /// Delete boundary face elements and wipe the surface mesh
 void delete_face_elements( Mesh* const & boundary_mesh_pt);

 /// Set up boundary condition elements on outer boundary
 void setup_outer_boundary();

 /// Set up point source
 void setup_point_source();

#ifdef ADAPTIVE

 /// Pointer to the "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif

 /// Mesh as geometric object representation of bulk mesh
 MeshAsGeomObject* Mesh_as_geom_obj_pt;

 /// Pointer to mesh containing the DtN (or ABC) boundary
 /// condition elements
 HelmholtzDtNMesh<ELEMENT>* Helmholtz_outer_boundary_mesh_pt;
 
 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
HelmholtzPointSourceProblem<ELEMENT>::
HelmholtzPointSourceProblem()
{ 

 // Open trace file
 Trace_file.open("RESLT/trace.dat");

 /// Initialise
 Mesh_as_geom_obj_pt=0;
 
 // Setup "bulk" mesh
 
 // Create circle representing outer boundary
 double x_c=0.0;
 double y_c=0.0;
 Circle* outer_circle_pt=new Circle(x_c,y_c,
                                    GlobalParameters::Outer_radius);
 
 // Outer boundary
 //---------------
 TriangleMeshClosedCurve* outer_boundary_pt=0;

 unsigned n_segments=40;
 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(2);
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start = 0.0;
 double s_end   = MathematicalConstants::Pi;
 unsigned boundary_id = 0;
 outer_boundary_line_pt[0]=
  new TriangleMeshCurviLine(outer_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = MathematicalConstants::Pi;
 s_end   = 2.0*MathematicalConstants::Pi;
 boundary_id = 1;
 outer_boundary_line_pt[1]=
  new TriangleMeshCurviLine(outer_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);
 
 // Create closed curve for outer boundary
 outer_boundary_pt=new TriangleMeshClosedCurve(outer_boundary_line_pt);
 

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters. The only parameter that needs to take is the
 // outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

#ifdef ADAPTIVE
 
 // Build "bulk" mesh
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=0.00004;
 Bulk_mesh_pt->max_permitted_error()=0.0001;

#else
 
 // Build "bulk" mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

#endif
 
 // Pointer to mesh containing the Helmholtz outer boundary condition
 // elements. Specify outer radius and number of Fourier terms to be
 // used in gamma integral
 unsigned n_fourier=10;
 Helmholtz_outer_boundary_mesh_pt = 
  new HelmholtzDtNMesh<ELEMENT>(GlobalParameters::Outer_radius,
                                n_fourier);
 
 // Create outer boundary elements from all elements that are 
 // adjacent to the outer boundary , but add them to a separate mesh.
 create_outer_bc_elements(0,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);
 create_outer_bc_elements(1,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);
 
 // Add the several  sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Helmholtz_outer_boundary_mesh_pt); 
  
 // Build the Problem's global mesh from its various sub-meshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional
 
 // Loop over the Helmholtz bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wave number squared
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the k_squared  pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Set up elements on outer boundary
 setup_outer_boundary();

 // Setup point source
 setup_point_source();

  // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::actions_before_adapt()
{ 
 // Wipe the boundary mesh
 delete_face_elements(Helmholtz_outer_boundary_mesh_pt);


 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::actions_after_adapt()
{

 // Complete the build of all elements so they are fully functional
 
 // Loop over the Helmholtz bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wave number squared
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the k_squared  pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Create BC elements 
 // from all elements that are adjacent to the boundaries and add them to 
 // Helmholtz_boundary_meshes
 create_outer_bc_elements(0,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);
 create_outer_bc_elements(1,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Set up DtN mesh
 setup_outer_boundary();
  
 // Setup point source
 setup_point_source();
 
}// end of actions_after_adapt



//==================start_of_setup_point_source===========================
/// Set point source
//========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::setup_point_source()
{ 

 // Create mesh as geometric object 
 delete Mesh_as_geom_obj_pt;
 Mesh_as_geom_obj_pt= new MeshAsGeomObject(Bulk_mesh_pt);

 // Position of point source in Eulerian coordinates
 Vector<double> x_point_source(2);
 x_point_source[0]=0.2;
 x_point_source[1]=0.3;

 GeomObject* sub_geom_object_pt=0;
 Vector<double> s_point_source(2);
 Mesh_as_geom_obj_pt->locate_zeta(x_point_source,sub_geom_object_pt,
                                  s_point_source);

 // Set point force
 std::complex<double> magnitude(2.0,3.0); 
 dynamic_cast<ELEMENT*>(sub_geom_object_pt)->setup(s_point_source,magnitude);
}




//==================start_of_setup_outer_boundary=========================
/// Set pointers for elements on outer boundary
//========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::setup_outer_boundary()
{ 
 // Loop over the flux elements to pass pointer to DtN
 // BC for the outer boundary
 unsigned n_element=Helmholtz_outer_boundary_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Dirichlet to Neumann BC
   if (GlobalParameters::DtN_BC)
    {
     // Upcast from GeneralisedElement to Helmholtz flux element
     HelmholtzDtNBoundaryElement<ELEMENT> *el_pt = 
      dynamic_cast< HelmholtzDtNBoundaryElement<ELEMENT>*>(
       Helmholtz_outer_boundary_mesh_pt->element_pt(e));
     
     // Set pointer to the mesh that contains all the boundary condition
     // elements on this boundary
     el_pt->set_outer_boundary_mesh_pt(Helmholtz_outer_boundary_mesh_pt);
    }
   // ABCs BC
   else
    {
     // Upcast from GeneralisedElement to appropriate type
     HelmholtzAbsorbingBCElement<ELEMENT> *el_pt = 
      dynamic_cast< HelmholtzAbsorbingBCElement<ELEMENT>*>(
       Helmholtz_outer_boundary_mesh_pt->element_pt(e));
     
     // Set pointer to outer radius of artificial boundary
     el_pt->outer_radius_pt()=&GlobalParameters::Outer_radius;
     
     // Set order of absorbing boundary condition
     el_pt->abc_order_pt()=&GlobalParameters::ABC_order;
    }    
  }
}



//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::doc_solution(DocInfo& 
                                              doc_info) 
{ 

 ofstream some_file,some_file2;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 

 // Total radiated power
 double power=0.0;

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 // Accumulate contribution from elements
 unsigned nn_element=Helmholtz_outer_boundary_mesh_pt->nelement(); 
 for(unsigned e=0;e<nn_element;e++)
  {
   HelmholtzBCElementBase<ELEMENT> *el_pt = 
    dynamic_cast< HelmholtzBCElementBase<ELEMENT>*>(
     Helmholtz_outer_boundary_mesh_pt->element_pt(e)); 
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();
 oomph_info << "Total radiated power: " << power << std::endl; 

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output solution in Paraview
 //---------------------------
 sprintf(filename,"%s/soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_paraview(some_file,npts);
 some_file.close();

 // Write power to trace file
 Trace_file  << power << std::endl;

} // end of doc


//============start_of_create_outer_bc_elements==============================
/// Create outer BC elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointed to by helmholtz_outer_boundary_mesh_pt.
//===========================================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::
create_outer_bc_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                         Mesh* const & helmholtz_outer_boundary_mesh_pt)
{
 // Loop over the bulk elements adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding outer flux element
   
   // Dirichlet to Neumann boundary conditon
   if (GlobalParameters::DtN_BC)
    {
     HelmholtzDtNBoundaryElement<ELEMENT>* flux_element_pt = new 
      HelmholtzDtNBoundaryElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
     helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);
    }
   //  ABCs BC
   else
    {
     HelmholtzAbsorbingBCElement<ELEMENT>* flux_element_pt = new 
      HelmholtzAbsorbingBCElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
     helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);
    }
  } //end of loop over bulk elements adjacent to boundary b
} // end of create_outer_bc_elements


//============start_of_delete_face_elements================
/// Delete face elements and wipe the boundary mesh
//==========================================================
template<class ELEMENT>
void HelmholtzPointSourceProblem<ELEMENT>::
delete_face_elements(Mesh* const & boundary_mesh_pt)
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
 
} // end of delete_outer_face_elements



//==========start_of_main=================================================
/// Solve 2D Helmholtz problem for scattering of a planar wave from a 
/// unit disk 
//========================================================================
int main(int argc, char **argv)
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define case to be run
 unsigned i_case=0;
 CommandLineArgs::specify_command_line_flag("--case",&i_case);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Now set flags accordingly
 switch(i_case)
  {
  case 0:
   GlobalParameters::DtN_BC=true;
   break;
   
  case 1:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=1;
   break;

  case 2:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=2;
   break;

  case 3:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=3;
   break;
  }
 
 
 //Set up the problem
 //------------------
 
#ifdef ADAPTIVE
 
 //Set up the problem with 2D nine-node elements from the
 //QHelmholtzElement family.
 HelmholtzPointSourceProblem<ProjectableHelmholtzElement
                             <HelmholtzPointSourceElement
                              <THelmholtzElement<2,3> > > >
  problem;
 
#else
 
 // Set up the problem with 2D six-node elements from the
 // THelmholtzElement family. 
 HelmholtzPointSourceProblem<HelmholtzPointSourceElement
                             <THelmholtzElement<2,3> > > problem;
 

#endif
 
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");


#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=5;
 
   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   problem.newton_solve(max_adapt);

#else

   // Solve the problem with Newton's method
   problem.newton_solve();

#endif

 //Output solution
 problem.doc_solution(doc_info);
    
} //end of main



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
#include <fenv.h> 

//Generic routines
#include "generic.h" 

// The equations
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

#include "locate_zeta_tester.h"


//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=10.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }
 

 ///  Zero function -- used to compute norm of the computed solution by 
 /// computing the norm of the error when compared against this.
 void zero(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=0.0;
 }

} // end of namespace



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredPoissonProblem : public virtual Problem
{

public:

 /// Constructor
 UnstructuredPoissonProblem();
    
 /// Destructor
 ~UnstructuredPoissonProblem(){};

 /// Actions before adapt. Empty
 void actions_before_adapt() {}
 
 ///  Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   complete_problem_setup();
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve: Re-apply boundary conditons
 void actions_before_newton_solve()
  {
   apply_boundary_conditions();
  }
  
 /// Doc the solution
 void doc_solution(const std::string& comment="");
 

private:

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 ///  Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Pointers to specific mesh
 RefineableTriangleMesh<ELEMENT>* My_mesh_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

}; // end_of_problem_class





//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
UnstructuredPoissonProblem<ELEMENT>::UnstructuredPoissonProblem()
{          
 // Intrinsic coordinate along GeomObject
 Vector<double> zeta(1);

 // Position vector on GeomObject
 Vector<double> posn(2);
 
 // Ellipse defining the outer boundary
 double x_center = 0.0;
 double y_center = 0.0;
 double A = 1.0;
 double B = 1.0;
 Ellipse * outer_boundary_ellipse_pt = new Ellipse(A,B);

 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon?
 //---------------------------------
 bool polygon_for_outer_boundary=false;
#ifdef OUTER_POLYGON
 polygon_for_outer_boundary=true;
#endif
 if (polygon_for_outer_boundary)
  { 
   // Number of segments that make up the boundary
   unsigned n_seg = 5; 
   double unit_zeta = 0.5*MathematicalConstants::Pi/double(n_seg);
   
   // The boundary is bounded by two distinct boundaries, each
   // represented by its own polyline
   Vector<TriangleMeshCurveSection*> boundary_polyline_pt(2);
   
   // Vertex coordinates on boundary
   Vector<Vector<double> > bound_coords(n_seg+1);
   
   // First part of the boundary 
   //---------------------------
   for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
    {
     // Resize the vector 
     bound_coords[ipoint].resize(2);
     
     // Get the coordinates
     zeta[0]=unit_zeta*double(ipoint);
     outer_boundary_ellipse_pt->position(zeta,posn);
     bound_coords[ipoint][0]=posn[0]+x_center;
     bound_coords[ipoint][1]=posn[1]+y_center;
    }
   
   // Build the 1st boundary polyline
   unsigned boundary_id=0;
   boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,boundary_id);
   
   // Second part of the boundary
   //----------------------------
   unit_zeta*=3.0;
   for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
    {
     // Resize the vector 
     bound_coords[ipoint].resize(2);
     
     // Get the coordinates
     zeta[0]=(unit_zeta*double(ipoint))+0.5*MathematicalConstants::Pi;
     outer_boundary_ellipse_pt->position(zeta,posn);
     bound_coords[ipoint][0]=posn[0]+x_center;
     bound_coords[ipoint][1]=posn[1]+y_center;
    }
   
   // Build the 2nd boundary polyline
   boundary_id=1;
   boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,boundary_id);
   

   // Create the triangle mesh polygon for outer boundary
   //----------------------------------------------------
   TriangleMeshPolygon *outer_polygon =
    new TriangleMeshPolygon(boundary_polyline_pt);

   // Enable redistribution of polylines
   outer_polygon->
    enable_redistribution_of_segments_between_polylines();

   // Set the pointer
   closed_curve_pt = outer_polygon;

  }
 // Build outer boundary as curvilinear
 //------------------------------------
 else
  {   

   // Provide storage for pointers to the two parts of the curvilinear boundary
   Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
   
   // First bit
   //----------
   double zeta_start=0.0;
   double zeta_end=MathematicalConstants::Pi;
   unsigned nsegment=5;
   unsigned boundary_id=0;
   outer_curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
    outer_boundary_ellipse_pt,zeta_start,zeta_end,nsegment,boundary_id);
   
   // Second bit
   //-----------
   zeta_start=MathematicalConstants::Pi;
   zeta_end=2.0*MathematicalConstants::Pi;
   nsegment=8;
   boundary_id=1;
   outer_curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
    outer_boundary_ellipse_pt,zeta_start,zeta_end,nsegment,boundary_id);
   
   // Combine to curvilinear boundary and define the
   //--------------------------------
   // outer boundary
   //--------------------------------
   closed_curve_pt=
     new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
   
  }
 
 
 // Now build the holes
 //====================
 Vector<TriangleMeshClosedCurve*> hole_pt(2);

 // Build polygonal hole
 //=====================
 
 // Build first hole: A circle
 x_center = 0.0;
 y_center = 0.5;
 A = 0.1;
 B = 0.1;
 Ellipse* polygon_ellipse_pt=new Ellipse(A,B);
 
 // Number of segments defining upper and lower half of the hole
 unsigned n_seg = 6; 
 double unit_zeta = MathematicalConstants::Pi/double(n_seg);
 
 // This hole is bounded by two distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> hole_polyline_pt(2);
 

 // First boundary of polygonal hole
 //---------------------------------

 // Vertex coordinates
 Vector<Vector<double> > bound_hole(n_seg+1);
 for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=unit_zeta*double(ipoint);
   polygon_ellipse_pt->position(zeta,posn);
   bound_hole[ipoint][0]=posn[0]+x_center;
   bound_hole[ipoint][1]=posn[1]+y_center;
  }
 
 // Specify the hole boundary id
 unsigned boundary_id=2;
 
 // Build the 1st hole polyline
 hole_polyline_pt[0] = new TriangleMeshPolyLine(bound_hole,boundary_id);
 

 // Second boundary of polygonal hole
 //----------------------------------
 for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=(unit_zeta*double(ipoint))+MathematicalConstants::Pi;
   polygon_ellipse_pt->position(zeta,posn);
   bound_hole[ipoint][0]=posn[0]+x_center;
   bound_hole[ipoint][1]=posn[1]+y_center;
  }
 
 // Specify the hole boundary id
 boundary_id=3;
 
 // Build the 2nd hole polyline
 hole_polyline_pt[1] = new TriangleMeshPolyLine(bound_hole,boundary_id);


 // Build the polygonal hole 
 //-------------------------
 
 // Inner hole center coordinates
 Vector<double> hole_center(2);
 hole_center[0]=x_center;
 hole_center[1]=y_center;

 hole_pt[0] = new TriangleMeshPolygon(hole_polyline_pt, hole_center);
 

 // Build curvilinear hole
 //======================
 
 // Build second hole: Another ellipse
 A = 0.2;
 B = 0.1;
 Ellipse* ellipse_pt=new Ellipse(A,B);
 
 // Build the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);
 

 // First part of curvilinear boundary
 //-----------------------------------
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned nsegment=10;
 boundary_id=4;
 curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  ellipse_pt,zeta_start,zeta_end, 
  nsegment,boundary_id);
 
 // Second part of curvilinear boundary
 //-------------------------------------
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
 nsegment=15;
 boundary_id=5;
 curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
  ellipse_pt,zeta_start,zeta_end, 
  nsegment,boundary_id);
 
 
 // Combine to hole
 //----------------
 Vector<double> hole_coords(2);
 hole_coords[0]=0.0;
 hole_coords[1]=0.0;
 Vector<TriangleMeshClosedCurve*> curvilinear_hole_pt(1);
 hole_pt[1]=
  new TriangleMeshClosedCurve(curvilinear_boundary_pt,
                                                 hole_coords);
 
 // Uncomment this as an exercise to observe how a
 // layer of fine elements get left behind near the boundary
 // once the tanh step has swept past: 

 // closed_curve_pt->disable_polyline_refinement();
 // closed_curve_pt->disable_polyline_unrefinement();
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;

 // Specify the maximum area element
 double uniform_element_area=0.2;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 My_mesh_pt=new 
  RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
 
 // Store as the problem's one and only mesh
 Problem::mesh_pt()=My_mesh_pt;

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 My_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // // Set element size limits
 // My_mesh_pt->max_element_size()=0.2;
 // My_mesh_pt->min_element_size()=0.002; 

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Open trace file
 char filename[100];
 sprintf(filename,"RESLT/trace.dat");
 Trace_file.open(filename);

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " 
            << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor




//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::complete_problem_setup()
{   

 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=My_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=My_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin one-and-only unknown value
     nod_pt->pin(0);
    }   
  } // end loop over boundaries
 
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = My_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(e));
   
   //Set the source function pointer
   el_pt->source_fct_pt() = &TanhSolnForPoisson::get_source;
  }
 
 // Re-apply Dirichlet boundary conditions (projection ignores
 // boundary conditions!)
 apply_boundary_conditions();
}




//==start_of_apply_bc=====================================================
 /// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::apply_boundary_conditions()
{
 
 // Loop over all boundary nodes
 unsigned nbound=this->My_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=this->My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=this->My_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     
     // Compute the value of the exact solution at the nodal point
     Vector<double> u(1);
     TanhSolnForPoisson::get_exact_u(x,u);
     
     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  } 

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::doc_solution(const 
                                                       std::string& comment)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 
 
 sprintf(filename,"RESLT/soln%i.dat",Doc_info.number());
 some_file.open(filename);
 this->My_mesh_pt->output(some_file,npts); 
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
           << comment << "\"\n";
 some_file.close();
 
 // Output exact solution 
 //----------------------
 sprintf(filename,"RESLT/exact_soln%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();
 
 // Output boundaries
 //------------------
 sprintf(filename,"RESLT/boundaries%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm,dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                           error,norm); 
 
 My_mesh_pt->compute_error(some_file,TanhSolnForPoisson::zero,
                           dummy_error,zero_norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 oomph_info << "\nNorm of error   : " << sqrt(error) << std::endl; 
 oomph_info << "Norm of exact solution: " << sqrt(norm) << std::endl;
 oomph_info << "Norm of computed solution: " << sqrt(dummy_error) << std::endl;
 Trace_file << sqrt(norm) << " " << sqrt(dummy_error) << std::endl;

 // Increment the doc_info number
 Doc_info.number()++;


 // Only used to doc sample points
 if (false)
  {
   unsigned nplot=npts;
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(0));
  
   el_pt->node_pt(0)->x(0)=0.0;
   el_pt->node_pt(0)->x(1)=0.0;

   el_pt->node_pt(1)->x(0)=1.0;
   el_pt->node_pt(1)->x(1)=0.0;

   el_pt->node_pt(2)->x(0)=0.0;
   el_pt->node_pt(2)->x(1)=1.0;

   el_pt->node_pt(3)->x(0)=0.5;
   el_pt->node_pt(3)->x(1)=0.0;

   el_pt->node_pt(4)->x(0)=0.5;
   el_pt->node_pt(4)->x(1)=0.5;

   el_pt->node_pt(5)->x(0)=0.0;
   el_pt->node_pt(5)->x(1)=0.5;

   sprintf(filename,"RESLT/element0.dat");
   some_file.open(filename);
   el_pt->output(some_file,nplot);
   some_file.close();

   //Vector of local coordinates
   Vector<double> s(2);
   Vector<double> s_orig(2);
  
   sprintf(filename,"RESLT/sample_points_in_element0.dat");
   some_file.open(filename);
  
   // Tecplot header info
   some_file << el_pt->tecplot_zone_string(nplot);
  
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
    
     // Get unshifted local coordinates of plot point
     el_pt->get_s_plot(iplot,nplot,s_orig);
    
     // Get local coordinates of plot point
     bool shift=true;
     el_pt->get_s_plot(iplot,nplot,s,shift);
    
     oomph_info << "Sample points: ";
     for(unsigned i=0;i<2;i++) 
      {
       oomph_info << s[i] << " [ " << s_orig[i] << " ] "  ;
       some_file << el_pt->interpolated_x(s,i) << " ";
      }
     oomph_info << std::endl;
     some_file << el_pt->interpolated_u_poisson(s) << std::endl;   
    
    }
  
   // Write tecplot footer (e.g. FE connectivity lists)
   el_pt->write_tecplot_zone_footer(some_file,nplot);

   some_file.close();
   exit(0);
  }

} // end of doc


//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Create problem
 UnstructuredPoissonProblem<ProjectablePoissonElement<TPoissonElement<2,3> > >
  problem;
 
 
 // Doc the initial mesh
 //=====================
 std::stringstream comment_stream;
 comment_stream << "Initial mesh ";
 problem.doc_solution(comment_stream.str());
   
 
 // Solve with spatial adaptation
 //==============================
 unsigned max_adapt=3;
 problem.newton_solve(max_adapt);
 
 // Doc the solution
 //=================
 problem.doc_solution();

 // Check it out!
 check_locate_zeta(problem.mesh_pt());
    
   
 
} //End of main

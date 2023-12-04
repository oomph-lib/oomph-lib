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
#include <fenv.h> 

//Generic routines
#include "generic.h" 

// The equations
#include "foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;




//========================================================
/// Global parameters
//========================================================
namespace GlobalParameters
{

 /// FvK parameter
 double Eta = 2.39e6;

 /// The "bubble" radius
 double R_b = 0.1;

 /// The pressure
 double Pressure=0.0;

 /// Pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& x, double& pressure)
  {
   pressure=Pressure;
  }

 /// Function to compute norm of solution itself (we treat this
 /// as the "exact" solution) 
 void zero(const Vector<double>& x, Vector<double>& u)
  {
   u[0]=0.0;
  }


}


/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

 /// Constructor
 UnstructuredFvKProblem(double element_area = 0.2);
    
 /// Destructor
 ~UnstructuredFvKProblem()
  {
   delete My_mesh_pt;
  }

 
 /// Update after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution();
 

private:

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Pointers to specific mesh
 TriangleMesh<ELEMENT>* My_mesh_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

 // Keep track of boundary ids
 enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Inner_boundary0 = 2,
   Inner_boundary1 = 3,
   Inner_boundary2 = 4,
   Inner_boundary3 = 5,
   Inner_boundary4 = 6,
   Inner_boundary5 = 7
  };


 /// Element area
 double Element_area;

  /// Geom object made of mesh
  MeshAsGeomObject* Mesh_as_geom_object_pt;
  
  /// Plot points along radial line are located in this element/local coord
 Vector<std::pair<GeomObject*,Vector<double> > > Radial_sample_point_pt;

}; // end_of_problem_class


//========================================================================
/// Constructor: Pass in element area
//========================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
 : Element_area(element_area)
{
 Vector<double> zeta(1);
 Vector<double> posn(2);

 //Outer boundary
 //--------------
 Ellipse* outer_boundary_ellipse_pt = new Ellipse(1.0,1.0);

 TriangleMeshClosedCurve* outer_boundary_pt = 0;

 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
 
 //First bit
 double zeta_start = 0.0;
 double zeta_end = MathematicalConstants::Pi;
 unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
 outer_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Outer_boundary0);

 //Second bit
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
 outer_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Outer_boundary1);

 outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);


 //Inner boundary
 //--------------
  Ellipse* inner_boundary_ellipse_pt = new Ellipse(GlobalParameters::R_b,
                                                  GlobalParameters::R_b);

 Vector<TriangleMeshClosedCurve*> inner_boundary_pt(1);
 Vector<TriangleMeshCurveSection*> inner_curvilinear_boundary_pt(2);

 //First part
 zeta_start = 0.0;
 zeta_end = MathematicalConstants::Pi;
 nsegment = (int)(GlobalParameters::R_b*
                  MathematicalConstants::Pi/sqrt(element_area));
 inner_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(inner_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Inner_boundary0);

 //Second part
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(GlobalParameters::R_b*
                  MathematicalConstants::Pi/sqrt(element_area));
 inner_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(inner_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Inner_boundary1);

 //Combine to internal curvilinear boundary
 Vector<double> left_region_coords(2);
 Vector<double> right_region_coords(2);
 left_region_coords[0] = -0.5*GlobalParameters::R_b;
 left_region_coords[1] = 0.0;

 right_region_coords[0] = 0.5*GlobalParameters::R_b;
 right_region_coords[1] = 0.0;

 inner_boundary_pt[0] =
  new TriangleMeshClosedCurve(inner_curvilinear_boundary_pt);

 //Diameter boundary line
 Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt(4);

 Vector<Vector<double> > vertex_coord(2, Vector<double>(2));

 //Section 1 of line
 vertex_coord[0][0] = 0.0;
 vertex_coord[0][1] = -1.0;
 vertex_coord[1][0] = 0.0;
 vertex_coord[1][1] = -GlobalParameters::R_b;

 TriangleMeshPolyLine *inner_open_polyline1_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary2);

 inner_open_polyline1_pt->connect_initial_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(outer_curvilinear_boundary_pt[1]),
   1.5*MathematicalConstants::Pi);

 inner_open_polyline1_pt->connect_final_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(inner_curvilinear_boundary_pt[1]),
   1.5*MathematicalConstants::Pi);

 Vector<TriangleMeshCurveSection *> inner_boundary_line1_pt(1);
 inner_boundary_line1_pt[0] = inner_open_polyline1_pt;

 inner_open_boundary_pt[0] =
  new TriangleMeshOpenCurve(inner_boundary_line1_pt);

 //Section 2 of line
 vertex_coord[0][0] = 0.0;
 vertex_coord[0][1] = -GlobalParameters::R_b;
 vertex_coord[1][0] = 0.0;
 vertex_coord[1][1] = 0.0;

 TriangleMeshPolyLine *inner_open_polyline2_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary3);

 inner_open_polyline2_pt->connect_initial_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(inner_curvilinear_boundary_pt[1]),
   1.5*MathematicalConstants::Pi);

 Vector<TriangleMeshCurveSection *> inner_boundary_line2_pt(1);
 inner_boundary_line2_pt[0] = inner_open_polyline2_pt;

 inner_open_boundary_pt[1] =
  new TriangleMeshOpenCurve(inner_boundary_line2_pt);

 //Section 3 of line
 vertex_coord[0][0] = 0.0;
 vertex_coord[0][1] = 0.0;
 vertex_coord[1][0] = 0.0;
 vertex_coord[1][1] = GlobalParameters::R_b;

 TriangleMeshPolyLine *inner_open_polyline3_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary4);

 inner_open_polyline3_pt->connect_initial_vertex_to_polyline(
   inner_open_polyline2_pt,1);

 inner_open_polyline3_pt->connect_final_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(inner_curvilinear_boundary_pt[0]),
   0.5*MathematicalConstants::Pi);

 Vector<TriangleMeshCurveSection *> inner_boundary_line3_pt(1);
 inner_boundary_line3_pt[0] = inner_open_polyline3_pt;

 inner_open_boundary_pt[2] =
  new TriangleMeshOpenCurve(inner_boundary_line3_pt);

 //Section 4 of line
 vertex_coord[0][0] = 0.0;
 vertex_coord[0][1] = GlobalParameters::R_b;
 vertex_coord[1][0] = 0.0;
 vertex_coord[1][1] = 1.0;

 TriangleMeshPolyLine *inner_open_polyline4_pt =
  new TriangleMeshPolyLine(vertex_coord, Inner_boundary5);

 inner_open_polyline4_pt->connect_initial_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(inner_curvilinear_boundary_pt[0]),
   0.5*MathematicalConstants::Pi);

 inner_open_polyline4_pt->connect_final_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine *>(outer_curvilinear_boundary_pt[0]),
   0.5*MathematicalConstants::Pi);

 Vector<TriangleMeshCurveSection *> inner_boundary_line4_pt(1);
 inner_boundary_line4_pt[0] = inner_open_polyline4_pt;

 inner_open_boundary_pt[3] =
  new TriangleMeshOpenCurve(inner_boundary_line4_pt);

 //Create the mesh
 //---------------

 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(outer_boundary_pt);
 mesh_parameters.internal_closed_curve_pt() = inner_boundary_pt;
 mesh_parameters.internal_open_curves_pt() = inner_open_boundary_pt;
 mesh_parameters.add_region_coordinates(1,left_region_coords);
 mesh_parameters.add_region_coordinates(2,right_region_coords);
 mesh_parameters.element_area() = element_area;

 // Build the bloody thing
 My_mesh_pt = new TriangleMesh<ELEMENT>(mesh_parameters);


 complete_problem_setup();
 add_sub_mesh(My_mesh_pt);
 build_global_mesh();

 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 oomph_info << "Number of equations: "
            << this->assign_eqn_numbers() << '\n';
}



//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{   

 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet 
 // conditions here. 
 unsigned nbound = Outer_boundary1 + 1;
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=My_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin unknown values
     nod_pt->pin(0);  // pin vertical displacement; Laplacian of w unpinned
                      // imposes clamping

     // Leave in-plane displacements unpinned; note that disk randomly
     // moves/rotates around xy plane as a result
    }   
  } // end loop over boundaries
 


 // Pin in-plane displacement at origin (still allows for rotations!)
 // -----------------------------------------------------------------
 bool done_pin=false;
 unsigned n_boundary3_node = My_mesh_pt->nboundary_node(Inner_boundary3);
 unsigned n_boundary4_node = My_mesh_pt->nboundary_node(Inner_boundary4);
 for(unsigned inode=0; inode < n_boundary3_node; inode++)
  {
   for(unsigned jnode=0; jnode < n_boundary4_node; jnode++)
    {
     if(My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)
        == My_mesh_pt->boundary_node_pt(Inner_boundary4,jnode))
      {
       My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)->pin(2);
       My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)->pin(3);
       done_pin=true;
       oomph_info 
        << "Pinned horizontal/vertical displacement at [ "
        << My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)->x(0) << " "
        << My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)->x(1) << " ] "
        << "to suppress rigid body motion.\n";
      }
    }
  }
 if (!done_pin)
  {   
   oomph_info << "Not found node to suppress rigid body displacement\n";
   abort();
  }

 // Pin horizontal displacement at point at the very top
 // ----------------------------------------------------
 bool done=false;
 unsigned n_boundary0_node = My_mesh_pt->nboundary_node(Outer_boundary0);
 unsigned n_boundary5_node = My_mesh_pt->nboundary_node(Inner_boundary5);
 for(unsigned inode=0; inode < n_boundary0_node; inode++)
  {
   for(unsigned jnode=0; jnode < n_boundary5_node; jnode++)
    {
     if(My_mesh_pt->boundary_node_pt(Outer_boundary0,inode)
        == My_mesh_pt->boundary_node_pt(Inner_boundary5,jnode))
      {
       My_mesh_pt->boundary_node_pt(Outer_boundary0,inode)->pin(2);
       done=true;
       oomph_info 
        << "Pinned horizontal displacement at [ "
        << My_mesh_pt->boundary_node_pt(Outer_boundary0,inode)->x(0) << " "
        << My_mesh_pt->boundary_node_pt(Outer_boundary0,inode)->x(1) << " ] "
        << "to suppress rotation.\n";
      }
    }
  }
 if (!done)
  {
   oomph_info << "Not found node to suppress rigid body rotation\n";
   abort();
  }

 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = My_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(e));
   
   //Set the pressure function pointers and the physical constants
   el_pt->eta_pt() = &GlobalParameters::Eta;
   el_pt->pressure_fct_pt() = &GlobalParameters::get_pressure;

  }

  // Create the mesh as Geom Object
  Mesh_as_geom_object_pt=new MeshAsGeomObject(My_mesh_pt);


  
  /// Number of sample points
  unsigned n_sample=50;
  const double dr = 1.0/double(n_sample);

  // Get 'em
  Vector<double> x(2,0.0);
  Radial_sample_point_pt.resize(n_sample);
  for (unsigned j=0;j<n_sample;j++)
   {
    Radial_sample_point_pt[j].second.resize(2);
    x[0]=double(j)*dr;
    
    // Get the element and its local coordinates
    Mesh_as_geom_object_pt->locate_zeta(x,
                                        Radial_sample_point_pt[j].first,
                                        Radial_sample_point_pt[j].second);
   }
  
}


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution()
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts = 5;
 
 sprintf(filename,"RESLT/soln%i.dat",Doc_info.number());
 some_file.open(filename);
 this->My_mesh_pt->output(some_file,npts); 
 some_file.close();
 
 // Output boundaries
 //------------------
 sprintf(filename,"RESLT/boundaries%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_boundaries(some_file);
 some_file.close();

 // Find solution at r = 0
 // ----------------------
 unsigned n_boundary3_node = My_mesh_pt->nboundary_node(Inner_boundary3);
 unsigned n_boundary4_node = My_mesh_pt->nboundary_node(Inner_boundary4);
 double w_0 = 0;
 for(unsigned inode=0; inode < n_boundary3_node; inode++)
  {
   for(unsigned jnode=0; jnode < n_boundary4_node; jnode++)
    {
     if(My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)
       == My_mesh_pt->boundary_node_pt(Inner_boundary4,jnode))
      {
       w_0 = My_mesh_pt->boundary_node_pt(Inner_boundary3,inode)->value(0);
      }
    }
  }


 // Plot solution along radial line
 sprintf(filename,"RESLT/soln_along_radial_line%i.dat",Doc_info.number());
 some_file.open(filename);
 Vector<double> x(2);
 unsigned nplot=Radial_sample_point_pt.size();
 for (unsigned j=0;j<nplot;j++)
  {
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Radial_sample_point_pt[j].first);
   Vector<double> s=Radial_sample_point_pt[j].second;
   el_pt->interpolated_x(s,x);
   some_file << x[0] << " " 
             << x[1] << " " 
             << el_pt->interpolated_w_fvk(s,0) << " " // w 
             << el_pt->interpolated_w_fvk(s,2) << " " // u_x
             << el_pt->interpolated_w_fvk(s,3) << " " // u_y
             << std::endl;
  }
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double dummy_error,zero_norm;
 sprintf(filename,"RESLT/norm%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->compute_error(some_file,GlobalParameters::zero,
                           dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Norm of computed solution: " 
            << sqrt(dummy_error) << std::endl;

 Trace_file << GlobalParameters::Pressure << " "
            << w_0 << '\n';

 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<TDisplacementBasedFoepplvonKarmanElement<3> >
  problem(0.01);

 double dp=0.1;
 unsigned n_step=10;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   n_step=2;
  }
 for (unsigned i=0;i<n_step;i++)
  {
   oomph_info << "Solving for p = " << GlobalParameters::Pressure << std::endl;

   // Solve the problem
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution();
   
   // Increment pressure
   GlobalParameters::Pressure+=dp; 
  }


} //End of main


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

namespace TestSoln
{
 double A = 1.0;
 double B = 2.0;
 double eta = 1;
 // The bubble radius
 double r_b = 0.1;

 // Storage for the pressures in the bubble and outside of the bubble
 Data *p_b_pt = 0;
 Data *p_0_pt = 0;

 // The volume of the bubble region under the membrane
 double prescribed_volume = 0.01;

 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& x, double& pressure)
  {
   double r = sqrt(x[0]*x[0] + x[1]*x[1]);

   if (r < r_b)
    {
     pressure = p_b_pt->value(0);
    }
   else
    {
     pressure = p_0_pt->value(0);
    }
  }

 void get_airy_forcing(const Vector<double>& x, double& airy_forcing)
  {
   airy_forcing = 0;
  }
 
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
  };

 /// Actions before adapt. Delete old volume constraint element and update
 /// global mesh
 void actions_before_adapt()
  {
   *Temp_pressure_pt =
    Volume_constraint_element_pt->pressure_data_pt()->value(0);
   Volume_constraint_mesh_pt->flush_element_and_node_storage();
   delete Volume_constraint_element_pt;
   rebuild_global_mesh();
  }
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 /// Also create new volume constraint element with previous pressure value
 void actions_after_adapt()
  {
   Volume_constraint_element_pt =
    new FoepplvonKarmanVolumeConstraintElement
    <ELEMENT, RefineableTriangleMesh> (My_mesh_pt,
                                       Bubble_regions,
                                       *Temp_pressure_pt);
   Volume_constraint_element_pt->set_prescribed_volume(
     &TestSoln::prescribed_volume);
   Volume_constraint_mesh_pt->add_element_pt(Volume_constraint_element_pt);
   TestSoln::p_b_pt = Volume_constraint_element_pt->pressure_data_pt();
   rebuild_global_mesh();
   complete_problem_setup();
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve()
  {
  }

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

 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Pointers to specific mesh
 //RefineableTriangleMesh<ELEMENT>* My_mesh_pt;
 RefineableTriangleMesh<ELEMENT>* My_mesh_pt;

 /// Mesh for the volume constraint element
 Mesh *Volume_constraint_mesh_pt;

 /// Single volume constraint element instance
 FoepplvonKarmanVolumeConstraintElement<ELEMENT, RefineableTriangleMesh>
  *Volume_constraint_element_pt;

 /// Temporary storage for the pressure used when switching between instances
 /// of volume constraint element
 double *Temp_pressure_pt;

 Vector<unsigned> Bubble_regions;

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

 double Element_area;
}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
 :
   Volume_constraint_mesh_pt(0),
   Volume_constraint_element_pt(0),
   Temp_pressure_pt(new double(0)),
   Element_area(element_area)
{
 Vector<double> zeta(1);
 Vector<double> posn(2);

 //Outer boundary
 //--------------

 double A = 1.0;
 double B = 1.0;
 Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);

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
 
 A = TestSoln::r_b;
 B = TestSoln::r_b;
 Ellipse* inner_boundary_ellipse_pt = new Ellipse(A, B);

 Vector<TriangleMeshClosedCurve*> inner_boundary_pt(1);
 Vector<TriangleMeshCurveSection*> inner_curvilinear_boundary_pt(2);

 //First part
 zeta_start = 0.0;
 zeta_end = MathematicalConstants::Pi;
 nsegment = (int)(TestSoln::r_b*MathematicalConstants::Pi/sqrt(element_area));
 inner_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(inner_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Inner_boundary0);

 //Second part
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(TestSoln::r_b*MathematicalConstants::Pi/sqrt(element_area));
 inner_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(inner_boundary_ellipse_pt, zeta_start,
    zeta_end, nsegment, Inner_boundary1);

 //Combine to internal curvilinear boundary
 Vector<double> left_region_coords(2);
 Vector<double> right_region_coords(2);
 left_region_coords[0] = -0.5*TestSoln::r_b;
 left_region_coords[1] = 0.0;

 right_region_coords[0] = 0.5*TestSoln::r_b;
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
 vertex_coord[1][1] = -TestSoln::r_b;

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
 vertex_coord[0][1] = -TestSoln::r_b;
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
 vertex_coord[1][1] = TestSoln::r_b;

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
 vertex_coord[0][1] = TestSoln::r_b;
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
 Bubble_regions.push_back(1);
 Bubble_regions.push_back(2);
 mesh_parameters.element_area() = element_area;

 My_mesh_pt = new RefineableTriangleMesh<ELEMENT>(mesh_parameters);

 Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
 My_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

 My_mesh_pt->max_element_size() = 0.2;
 My_mesh_pt->min_element_size() = 0.002;

 //My_mesh_pt->max_permitted_error()=0.00005;
 //My_mesh_pt->min_permitted_error()=0.00001;

 TestSoln::p_0_pt = new Data(1);
 TestSoln::p_0_pt->set_value(0,0);
 TestSoln::p_0_pt->pin(0);

 /// Create the initial volume constraint element
 Volume_constraint_element_pt
  = new FoepplvonKarmanVolumeConstraintElement
  <ELEMENT, RefineableTriangleMesh> (My_mesh_pt, Bubble_regions);
 Volume_constraint_element_pt->set_prescribed_volume(
   &TestSoln::prescribed_volume);

 /// Set the problem pressure to be that of the volume constraint element
 TestSoln::p_b_pt = Volume_constraint_element_pt->pressure_data_pt();

 /// Set initial pressure guess
 Volume_constraint_element_pt->pressure_data_pt()->set_value(0,0.1);

 complete_problem_setup();
 
 /// Create the volume constraint mesh and add the element to the mesh
 Volume_constraint_mesh_pt = new Mesh;
 Volume_constraint_mesh_pt->add_element_pt(Volume_constraint_element_pt);
 TestSoln::p_b_pt = Volume_constraint_element_pt->pressure_data_pt();

 /// Add the sub meshes and build the global mesh
 add_sub_mesh(Volume_constraint_mesh_pt);

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
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 //unsigned nbound=My_mesh_pt->nboundary();
 //Just loop over outer boundary since inner boundary doesn't have boundary
 //conditions
 unsigned nbound = Outer_boundary1 + 1;

 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=My_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin unknown values
     nod_pt->pin(0);
     nod_pt->pin(2);
    }   
  } // end loop over boundaries
 

 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = My_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(e));
   
   //Set the pressure function pointers and the physical constants
   el_pt->eta_pt() = &TestSoln::eta;
   el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
   el_pt->airy_forcing_fct_pt() = &TestSoln::get_airy_forcing;
  }

 // Add the bubble pressure as external data to the nodes in the bubble
 // region
 unsigned n_bubble_regions = Bubble_regions.size();
 
 for(unsigned r = 0; r < n_bubble_regions; r++)
  {
   n_element = My_mesh_pt->nregion_element(Bubble_regions[r]);
   for(unsigned e = 0; e < n_element; e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      My_mesh_pt->region_element_pt(Bubble_regions[r],e));
     
     el_pt->set_volume_constraint_pressure_data_as_external_data(
      Volume_constraint_element_pt->pressure_data_pt());
    }
  }
  
 // Re-apply Dirichlet boundary conditions (projection ignores
 // boundary conditions!)
 apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
 
 // Loop over all boundary nodes
 //unsigned nbound=this->My_mesh_pt->nboundary();
 //Just loop over outer boundary since inner boundary doesn't have boundary
 //conditions
 unsigned nbound = Outer_boundary1 + 1;

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
    }
  } 

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                        std::string& comment)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts = 5;
 
 sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 this->My_mesh_pt->output(some_file,npts); 
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
           << comment << "\"\n";
 some_file.close();
 
 // Output boundaries
 //------------------
 sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 My_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Output along radial boundary
 // ----------------------------
 sprintf(filename,"RESLT/radial_soln%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 for(unsigned b = Inner_boundary2; b <= Inner_boundary5; b++)
  {
   unsigned n_boundary_el = My_mesh_pt->nboundary_element(b);
   for(unsigned e = 0; e < n_boundary_el; e++)
    {
     My_mesh_pt->boundary_element_pt(b,e)->output(some_file,npts);
    }
  }
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


 // Output regions
 unsigned n_region = My_mesh_pt->nregion();
 if (n_region > 1)
  {
   for (unsigned r = 0; r < n_region; r++)
    {
     //Attempt to output elements in different regions
     sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
       Element_area);
     some_file.open(filename);
     unsigned nel = My_mesh_pt->nregion_element(r);
     for (unsigned e = 0; e < nel; e++)
      {
       My_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
      }
     some_file.close();
    }
  }

 // Calculate the integral (i.e. volume) over the inner region
 double bubble_volume = 0;
 if (n_region > 1)
  {
   sprintf(filename,"RESLT/bubble_volume%i.dat",Doc_info.number());
   some_file.open(filename);
   unsigned n_inner_el;
   unsigned n_bubble_regions = Bubble_regions.size();

   for(unsigned r = 0;r < n_bubble_regions; r++)
    {
     n_inner_el = My_mesh_pt->nregion_element(Bubble_regions[r]);
     for(unsigned e = 0; e < n_inner_el; e++)
      {
       ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
         My_mesh_pt->region_element_pt(Bubble_regions[r],e));
       if(el_pt != 0)
        {
         bubble_volume += el_pt->get_bounded_volume();
        }
      }
    }
   some_file << bubble_volume << '\n';
   some_file.close();
  }


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 //double error,norm,dummy_error,zero_norm;
 double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 My_mesh_pt->compute_error(some_file,TestSoln::zero,
                           dummy_error,zero_norm);
 some_file.close();

 // Doc L2 error and norm of solution
 oomph_info << "Norm of computed solution: " << sqrt(dummy_error) << std::endl;

 Trace_file << TestSoln::p_b_pt->value(0) << " "
            << TestSoln::p_0_pt->value(0) << " "
            << w_0 << " " << bubble_volume << '\n';

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

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<
  ProjectableFoepplvonKarmanElement<TFoepplvonKarmanElement<3> > >
  problem(0.01);

 problem.newton_solve();
 problem.doc_solution();

} //End of main


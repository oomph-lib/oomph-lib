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
//
// Linearised poro-elastic fsi problem -- flow in a cylindrical poroelastic tube

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "axisym_poroelasticity.h"
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

//=================================================================
/// A namespace for the physical parameters in the problem
//=================================================================
namespace Global_Physical_Variables
{

 // The directory in which the solution is output
 std::string Directory("RESLT");
 
 // Factor for reduction in inlet pressure
 double Chris_p_factor=1.0;

 // Fluid parameters
 // ----------------

 /// The Reynolds number
 double Re = 500.0;

 /// The Strouhal number
 double St = 1.0;

 /// Ratio of poro-elastic fluid to Navier-Stokes fluid densities -- 
 /// almost certainly one; we allow a variable value to be able to switch off
 /// wall inertia without having to set the Reynolds number to zero.
 double Poro_elastic_to_navier_stokes_fluid_density=1.0;

 /// The Womersley number -- dependent parameter; compute
 /// from Re and St
 double Wo = 0.0;

 /// Storage for Moens Korteweg wavespeed -- dependent parameter!
 double Wavespeed=0.0;
 
 // Poroelasticity parameters
 // -------------------------

 /// Define the poroelasticity inertia parameter -- dependent parameter; compute
 /// from Re, Q and density ratio
 double Lambda_sq=0.0;

 /// Non-dim permeability -- ratio of typical porous flux to fluid veloc
 /// scale
 double Permeability=5.0; 

 /// Poisson's ratio of drained poro-elastic medium
 double Nu = 0.35;

 /// Define the non-dimensional Young's modulus (ratio of actual 
 /// Young's modulus to the one used to non-dimensionalise (poro-elastic)
 /// stresses
 double E = 1.0;

 /// Biot parameter 
 double Alpha=1.0; 

 /// Porosity
 double Porosity=0.18; 

 /// Ratio of the densities of the fluid and solid phases of
 /// the poroelastic material
 double Density_ratio_poro=1.0; 

 /// Inverse slip rate coefficient
 double Inverse_slip_rate_coefficient=0.0; 

 // Derived poroelastic parameters
 // ------------------------------

 /// Lambda - first Lame parameter -- dependent parameter; compute from nu
 double Lambda_lame = 0.0;

 /// mu - second Lame parameter  -- dependent parameter; compute from nu
 double Mu_lame = 0.0;

 /// Ratio of the pore fluid density to the compound density --
 /// dependent parameter compute from density ratio and porosity.
 double Rho_f_over_rho = 0.0;

 // Mesh/problem parameters
 // -----------------------

 /// Target element area for fluid mesh
 double Element_area_fluid=0.002;

 /// Target element area for poro-elasticmesh
 double Element_area_solid=0.002;

 /// Fluid mesh boundary layer thickness
 double Fluid_mesh_bl_thickness=0.01;

 /// Poro-elastic mesh boundary layer thickness
 double Solid_mesh_bl_thickness=0.01;

 /// Inner radius of tube
 double Inner_radius=0.5;

 /// Outer radius of tube
 double Outer_radius=1.0;

 /// Length of tube
 double Length=3.0;

 /// Constant inlet pressure (for steady Poiseuille flow)
 double P_inlet_initial=500.0; // Length*4.0/(Inner_radius*Inner_radius);

 /// Actual (possibly time varying) inlet pressure -- initialised to
 /// P_inlet_initial.
 double P_inlet_const=P_inlet_initial;

 /// Increment for pressure (default: double the inlet pressure 
 /// over duration of then tanh step)
 double P_inlet_step=P_inlet_initial;

 /// Parameter for tanh origin for pressure incrementation
 double T_tanh=0.25; 

 /// Steepness parameter for tanh for pressure incrementation
 double Alpha_tanh=100.0;

 /// FSI parameter
 double Q=1.0e-9; 
  
 /// Inflow traction applied to the fluid mesh
 void fluid_inflow_boundary_traction(const double &time,
                                     const Vector<double> &x,
                                     const Vector<double> &n,
                                     Vector<double> &result)
 {
  result[0]=0.0;
  result[1]=P_inlet_const;
  result[2]=0.0;
 }



 //-----------------------------------------------------------------------------
 
 /// Helper function to update dependent parameters
 void update_dependent_parameters()
 {
#ifdef PARANOID
  if (Inner_radius!=0.5)
   {
    throw OomphLibError(
     "Inner radius must be 1/2 for non-dimensionalisation to be consistent.",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // Scale from Chris' reference values?
  if (CommandLineArgs::command_line_flag_has_been_set("--chris_p_factor"))
   {
    oomph_info 
     << "\n\nRESCALING ACCORDING TO A REDUCTION IN INLET PRESSURE OF "
     << Chris_p_factor << std::endl;

    Re*=sqrt(Chris_p_factor);
    Q*=sqrt(Chris_p_factor);
    Permeability/=sqrt(Chris_p_factor);

    // hierher
    P_inlet_initial*=Chris_p_factor;
    P_inlet_step*=Chris_p_factor;
   }

  // The Womersley number
  Wo = Re*St;
  
  // Lambda - first Lame parameter
  Lambda_lame = E*Nu/((1.0+Nu)*(1.0-2.0*Nu));
  
  // mu - second Lame parameter
  Mu_lame = E/(2.0*(1.0+Nu));
  
  // Ratio of the pore fluid density to the compound density
  Rho_f_over_rho = Density_ratio_poro/(1.0+Porosity*(Density_ratio_poro-1.0));

  // Poro-elastic timescale ratio
  Lambda_sq=Poro_elastic_to_navier_stokes_fluid_density*Re*Q/Rho_f_over_rho;

  // Non-dim wall thickness
  double h=Outer_radius-Inner_radius;


   
  //Wavespeed
  if ((Q*Re)!=0.0)
   {
    Wavespeed=sqrt(h/(2.0*Inner_radius)/(Q*Re));
   }
  else
   {
    Wavespeed=999.9e99;
   }
 }

 /// Doc dependent parameters
 void doc_dependent_parameters()
 {
  oomph_info << std::endl;
  oomph_info << "Problem parameters" << std::endl;
  oomph_info << "==================" << std::endl;
  oomph_info << "Reynolds number (Re)                           : "
             << Re << std::endl;
  oomph_info << "Strouhal number (St)                           : "
             << St << std::endl;
  oomph_info << "Womersley number (ReSt)                        : "
             << Wo << std::endl;
  oomph_info << "FSI parameter (Q)                              : "
             << Q << std::endl;
  oomph_info << "Wall thickness                                 : "
             << Outer_radius-Inner_radius << std::endl;
  oomph_info << "Poro-elastic wall inertia parameter (Lambda^2) : "
             << Lambda_sq << std::endl;
  oomph_info << "Non-dim permeability                           : "
             << Permeability << std::endl;
  oomph_info << "Biot parameter (alpha)                         : "
             << Alpha << std::endl;
  oomph_info << "Porosity                                       : "
             << Porosity << std::endl;
  oomph_info << "Inverse slip rate coefficient                  : "
             <<Inverse_slip_rate_coefficient << std::endl;
  oomph_info << "Density ratio (rho_f/rho_s)                    : "
             << Density_ratio_poro << std::endl; 
  oomph_info << "Fluid density ratio (rho_f/rho_f_nst)          : "
             << Poro_elastic_to_navier_stokes_fluid_density 
             << std::endl;
  oomph_info << "rho_fluid/rho_compound                         : "
             << Rho_f_over_rho << std::endl;
  oomph_info << "Poisson's ratio (Nu)                           : "
             << Nu << std::endl;
  oomph_info << "First Lame parameter (Lambda_lame)             : "
             << Lambda_lame << std::endl;
  oomph_info << "Second Lame parameter (Mu_lame)                : "
             << Mu_lame << std::endl;
  oomph_info << std::endl << std::endl;


  // For information:
  {
   double nt=100;
   double dt=Global_Physical_Variables::St*
    (Global_Physical_Variables::Length/double(nt))/
    Global_Physical_Variables::Wavespeed;
   oomph_info 
    << "Timestep required to sub-divide pulse-wave propagation \n"
    << "along tube length into " << nt << " steps: " << dt << std::endl;
  }

 }

 /// Global function that completes the edge sign setup
 template<class ELEMENT>
 void edge_sign_setup(Mesh* mesh_pt)
 {
  // The dictionary keeping track of edge signs
  std::map<Edge,unsigned> assignments;
  
  // Loop over all elements
  unsigned n_element = mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
   {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt->element_pt(e));
  
  // Assign edge signs: Loop over the vertex nodes (always 
    // first 3 nodes for triangles)
    for(unsigned i=0;i<3;i++)
     {
      Node *nod_pt1, *nod_pt2;
      nod_pt1 = el_pt->node_pt(i);
      nod_pt2 = el_pt->node_pt((i+1)%3);
      Edge edge(nod_pt1,nod_pt2);
      unsigned status = assignments[edge];
      
      // This relies on the default value for an int being 0 in a map
      switch(status)
       {
        // If not assigned on either side, give it a + on current side
       case 0:
        assignments[edge]=1;
        break;
        // If assigned + on other side, give it a - on current side
       case 1:
        assignments[edge]=2;
        el_pt->sign_edge(i)=-1;
        break;
        // If assigned - on other side, give it a + on current side
       case 2:
        assignments[edge]=1;
        break;
       }
     } // end of loop over vertex nodes

   } // end of loop over elements
 }
 
}



//==========================================================================
/// Problem class
//==========================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
class PressureWaveFSIProblem : public Problem
{

public:

 /// Constructor
 PressureWaveFSIProblem();

 /// Create the fluid traction elements
 void create_fluid_traction_elements();

 /// Create the poroelasticity traction/pressure elements
 void create_poro_face_elements();

 /// Update the problem specs before next timestep. Increase
 /// inlet pressure
 void actions_before_implicit_timestep()
  {
   double time=time_pt()->time();
   Global_Physical_Variables::P_inlet_const =
    Global_Physical_Variables::P_inlet_initial + 
    Global_Physical_Variables::P_inlet_step*
    0.5*(1.0+tanh(Global_Physical_Variables::Alpha_tanh*
                  (time-Global_Physical_Variables::T_tanh)));

   oomph_info << "Updating inflow pressure for time " << time << " to " 
              << Global_Physical_Variables::P_inlet_const << std::endl;
  }

 /// Doc solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Setup fsi
 void setup_fsi();

 /// Fluid mesh
 TriangleMesh<FLUID_ELEMENT> *Fluid_mesh_pt;

 /// Poroelasticity mesh
 TriangleMesh<POROELASTICITY_ELEMENT> *Poro_mesh_pt;

 /// Poroelasticity surface mesh at FSI interface
 Mesh* FSI_poro_surface_mesh_pt;

 /// Inflow fluid surface mesh
 Mesh* Inflow_fluid_surface_mesh_pt;

 /// FSI fluid surface mesh
 Mesh* FSI_fluid_surface_mesh_pt;

 /// Pointer to the poroelasticity timestepper
 TimeStepper* Poro_time_stepper_pt;

 /// Pointer to the fluid timestepper
 TimeStepper* Fluid_time_stepper_pt;

 /// Id for Lagrange multiplier that enforces (no-)slip on fluid
 unsigned Lagrange_id;

 /// Mesh as geom object representation of fluid mesh
 MeshAsGeomObject* Fluid_mesh_geom_obj_pt;

 /// Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<FLUID_ELEMENT*,Vector<double> > > 
 Fluid_regularly_spaced_plot_point;

 /// Mesh as geom object representation of solid mesh
 MeshAsGeomObject* Solid_mesh_geom_obj_pt;

 /// Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<POROELASTICITY_ELEMENT*,Vector<double> > > 
 Solid_regularly_spaced_plot_point;

 /// Enumeration of fluid boundaries
 enum {Bottom_fluid_boundary, 
       FSI_fluid_boundary,
       Top_fluid_boundary,
       Centreline_fluid_boundary,
       Internal_boundary_layer_fluid_boundary};


 /// Enumeration of poro-elastic boundaries
 enum {Bottom_poro_boundary, 
       Outer_poro_boundary,
       Top_poro_boundary,
       FSI_poro_boundary,
       Bottom_internal_poro_boundary,
       Top_internal_poro_boundary};



};

//============================================================================
/// Constructor
//============================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
PressureWaveFSIProblem<FLUID_ELEMENT, POROELASTICITY_ELEMENT>::
PressureWaveFSIProblem()
{

  // Create timesteppers
 Poro_time_stepper_pt = new Newmark<2>;
 add_time_stepper_pt(Poro_time_stepper_pt);
 Fluid_time_stepper_pt = new BDF<2>;
 add_time_stepper_pt(Fluid_time_stepper_pt);

 // Create the problem geometry - a hollow poroelastic cylinder filled 
 //-------------------------------------------------------------------
 // with fluid
 //-----------

 // Fluid mesh
 // ----------

 // Vector for the vertices
 Vector<Vector<double> > fluid_vertex_coords(7, Vector<double>(2));
 fluid_vertex_coords[0][0]=0.0;
 fluid_vertex_coords[0][1]=0.0;
 
 fluid_vertex_coords[1][0]=Global_Physical_Variables::Inner_radius-
  Global_Physical_Variables::Fluid_mesh_bl_thickness;
 fluid_vertex_coords[1][1]=0.0;
 
 fluid_vertex_coords[2][0]=Global_Physical_Variables::Inner_radius;
 fluid_vertex_coords[2][1]=0.0;
 
 fluid_vertex_coords[3][0]=Global_Physical_Variables::Inner_radius;
 fluid_vertex_coords[3][1]=Global_Physical_Variables::Length;

 fluid_vertex_coords[4][0]=Global_Physical_Variables::Inner_radius-
  Global_Physical_Variables::Fluid_mesh_bl_thickness;
 fluid_vertex_coords[4][1]=Global_Physical_Variables::Length;

 fluid_vertex_coords[5][0]=0.0;
 fluid_vertex_coords[5][1]=Global_Physical_Variables::Length;

 fluid_vertex_coords[6][0]=0.0;
 fluid_vertex_coords[6][1]=0.0;

 // Loop over the vertices and create a polyline between each consecutive pair
 Vector<TriangleMeshCurveSection*> fluid_outer_polyline_boundary_pt(4);
 unsigned vertex_count=0;
 for(unsigned i=0;i<4;i++)
  {

   // =============================================================
   // Careful -- don't change the loop structure here as it has
   // to remain consistent with the enum of the boundary ids -- 
   // awkward...
   // =============================================================

   unsigned n=2;
   if ((i==0)||(i==2))  n=3;
   Vector<Vector<double> > temp_coord(n,Vector<double>(2));
   for (unsigned j=0;j<n;j++)
    {
     temp_coord[j][0]=fluid_vertex_coords[vertex_count][0];
     temp_coord[j][1]=fluid_vertex_coords[vertex_count][1];
     vertex_count++;
    }
   vertex_count--;
   fluid_outer_polyline_boundary_pt[i]=new TriangleMeshPolyLine(temp_coord,i);
  }


 // Internal boundary forcing at least one layer of elements
 // within a nominal boundary layer near the wall
 Vector<TriangleMeshOpenCurve*> bl_boundary_pt(1);
 Vector<Vector<double> > bl_coord(2,Vector<double>(2));
 bl_coord[0][0]=fluid_vertex_coords[1][0];
 bl_coord[0][1]=fluid_vertex_coords[1][1];
 bl_coord[1][0]=fluid_vertex_coords[4][0];
 bl_coord[1][1]=fluid_vertex_coords[4][1];
 TriangleMeshPolyLine* bl_polyline_pt=
  new TriangleMeshPolyLine(bl_coord,
                           Internal_boundary_layer_fluid_boundary);

 // Connect start point to first node on lower outer boundary polyline
 bl_polyline_pt->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (fluid_outer_polyline_boundary_pt[0]),1);

 // Connect end point to first node on upper outer boundary polyline
 bl_polyline_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (fluid_outer_polyline_boundary_pt[2]),1);

 Vector<TriangleMeshCurveSection*> bl_curve_section_pt(1);
 bl_curve_section_pt[0]=bl_polyline_pt;

 Vector<TriangleMeshOpenCurve*> bl_open_boundary_pt(1);
 bl_open_boundary_pt[0]=new TriangleMeshOpenCurve(bl_curve_section_pt);

 // Create the outer boundary closed curve
 TriangleMeshClosedCurve* fluid_outer_boundary_pt =
  new TriangleMeshClosedCurve(fluid_outer_polyline_boundary_pt);

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters fluid_triangle_mesh_parameters(
   fluid_outer_boundary_pt);

 // Set element area
 fluid_triangle_mesh_parameters.element_area()=
  Global_Physical_Variables::Element_area_fluid;

 // Set the internal boundary
 fluid_triangle_mesh_parameters.internal_open_curves_pt()=bl_open_boundary_pt;

 // Create the mesh
 Fluid_mesh_pt = new TriangleMesh<FLUID_ELEMENT>(
   fluid_triangle_mesh_parameters, Fluid_time_stepper_pt);
 

 // Mesh as geom object representation of fluid mesh
 Fluid_mesh_geom_obj_pt=new MeshAsGeomObject(Fluid_mesh_pt);
 
 // Extract points
 unsigned nplot_r=5;
 unsigned nplot_z=unsigned(double(nplot_r)*Global_Physical_Variables::Length/
                           Global_Physical_Variables::Inner_radius);
 Fluid_regularly_spaced_plot_point.resize(nplot_r*nplot_z);
 Vector<double> x(2);
 Vector<double> s(2);
 unsigned count=0;
 for (unsigned ir=0;ir<nplot_r;ir++)
  {
   x[0]=double(ir)/double(nplot_r-1)*Global_Physical_Variables::Inner_radius;
   for (unsigned iz=0;iz<nplot_z;iz++)
    {
     x[1]=double(iz)/double(nplot_z-1)*Global_Physical_Variables::Length;

     // Pointer to GeomObject that contains this point
     GeomObject* geom_obj_pt=0;

     // Get it
     Fluid_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);

     // Store it
     Fluid_regularly_spaced_plot_point[count].first=
      dynamic_cast<FLUID_ELEMENT*>(geom_obj_pt);
     Fluid_regularly_spaced_plot_point[count].second=s;

     count++;
    }
  }

 Fluid_mesh_pt->output_boundaries("fluid_mesh_boundaries.dat");


 // Poroelasticity mesh
 //-------------------

 // Vector for the vertices
 Vector<Vector<double> > poro_vertex_coords(5, Vector<double>(2));
 poro_vertex_coords[0][0]=Global_Physical_Variables::Inner_radius;
 poro_vertex_coords[0][1]=0.0;

 poro_vertex_coords[1][0]=Global_Physical_Variables::Outer_radius;
 poro_vertex_coords[1][1]=0.0;

 poro_vertex_coords[2][0]=Global_Physical_Variables::Outer_radius;
 poro_vertex_coords[2][1]=Global_Physical_Variables::Length;

 poro_vertex_coords[3][0]=Global_Physical_Variables::Inner_radius;
 poro_vertex_coords[3][1]=Global_Physical_Variables::Length;

 poro_vertex_coords[4][0]=Global_Physical_Variables::Inner_radius;
 poro_vertex_coords[4][1]=0.0;

 // Loop over the vertices and create a polyline between each consecutive pair
 Vector<Vector<double> > temp_coord(2,Vector<double>(2));
 Vector<TriangleMeshCurveSection*> poro_outer_polyline_boundary_pt(4);
 for(unsigned i=0;i<4;i++)
  {

   // =============================================================
   // Careful -- don't change the loop structure here as it has
   // to remain consistent with the enum of the boundary ids -- 
   // awkward...
   // =============================================================

   temp_coord[0][0]=poro_vertex_coords[i][0];
   temp_coord[0][1]=poro_vertex_coords[i][1];
   temp_coord[1][0]=poro_vertex_coords[i+1][0];
   temp_coord[1][1]=poro_vertex_coords[i+1][1];
   poro_outer_polyline_boundary_pt[i] =
    new TriangleMeshPolyLine(temp_coord,i);
  }

 // Create the outer boundary closed curve
 TriangleMeshClosedCurve* poro_outer_boundary_pt =
  new TriangleMeshClosedCurve(poro_outer_polyline_boundary_pt);

 // Internal boundary to prevent problem with flux B.C.s
 Vector<TriangleMeshOpenCurve*> poro_inner_open_boundary_pt;

 // We need 2 inner boundaries to split the corners
 poro_inner_open_boundary_pt.resize(2);

 // Temporary storage for coords in required format
 Vector<Vector<double> > temp_inner_coord(4,Vector<double>(2));

 // First inner boundary
 temp_inner_coord[0][0]=Global_Physical_Variables::Inner_radius;
 temp_inner_coord[0][1]=0.0;
 temp_inner_coord[1][0]=Global_Physical_Variables::Inner_radius+
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[1][1]=
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[2][0]=Global_Physical_Variables::Inner_radius+
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[2][1]=Global_Physical_Variables::Length-
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[3][0]=Global_Physical_Variables::Inner_radius;
 temp_inner_coord[3][1]=Global_Physical_Variables::Length;

 TriangleMeshPolyLine *poro_inner_open_polyline1_pt=
  new TriangleMeshPolyLine(temp_inner_coord,Bottom_internal_poro_boundary);

 poro_inner_open_polyline1_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (poro_outer_polyline_boundary_pt[0]),0);

 poro_inner_open_polyline1_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (poro_outer_polyline_boundary_pt[3]),0);

 Vector<TriangleMeshCurveSection*> 
  poro_inner_open_polyline1_curve_section_pt(1);
 poro_inner_open_polyline1_curve_section_pt[0]=poro_inner_open_polyline1_pt;

 poro_inner_open_boundary_pt[0]=
  new TriangleMeshOpenCurve(poro_inner_open_polyline1_curve_section_pt);

 // Second inner boundary
 temp_inner_coord[0][0]=Global_Physical_Variables::Outer_radius;
 temp_inner_coord[0][1]=Global_Physical_Variables::Length;
 temp_inner_coord[1][0]=Global_Physical_Variables::Outer_radius-
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[1][1]=Global_Physical_Variables::Length-
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[2][0]=Global_Physical_Variables::Outer_radius-
  Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[2][1]= Global_Physical_Variables::Solid_mesh_bl_thickness;
 temp_inner_coord[3][0]=Global_Physical_Variables::Outer_radius;
 temp_inner_coord[3][1]=0.0;

 TriangleMeshPolyLine *poro_inner_open_polyline2_pt=
  new TriangleMeshPolyLine(temp_inner_coord,Top_internal_poro_boundary);
 
 poro_inner_open_polyline2_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (poro_outer_polyline_boundary_pt[2]),0);

 // // hierher JULIO: WHY DOES THIS: not work?
 // poro_inner_open_polyline2_pt->connect_final_vertex_to_polyline(
 //   dynamic_cast<TriangleMeshPolyLine*>
 //   (poro_outer_polyline_boundary_pt[2]),1);

 poro_inner_open_polyline2_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (poro_outer_polyline_boundary_pt[0]),1);

 Vector<TriangleMeshCurveSection*> 
  poro_inner_open_polyline2_curve_section_pt(1);
 poro_inner_open_polyline2_curve_section_pt[0]=poro_inner_open_polyline2_pt;

 poro_inner_open_boundary_pt[1]=
  new TriangleMeshOpenCurve(poro_inner_open_polyline2_curve_section_pt);

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters poro_triangle_mesh_parameters(
   poro_outer_boundary_pt);

 // Set the inner boundaries
 poro_triangle_mesh_parameters.internal_open_curves_pt()=
  poro_inner_open_boundary_pt;

 // Set element area
 poro_triangle_mesh_parameters.element_area()=
  Global_Physical_Variables::Element_area_solid;

 // Create the mesh
 Poro_mesh_pt = new TriangleMesh<POROELASTICITY_ELEMENT>(
   poro_triangle_mesh_parameters, Poro_time_stepper_pt);

 // Mesh as geom object representation of solid mesh
 Solid_mesh_geom_obj_pt=new MeshAsGeomObject(Poro_mesh_pt);
 
 // Extract points
 nplot_r=5;
 nplot_z=unsigned(double(nplot_r)*
                  Global_Physical_Variables::Length/
                  (Global_Physical_Variables::Outer_radius-
                   Global_Physical_Variables::Inner_radius));
 Solid_regularly_spaced_plot_point.resize(nplot_r*nplot_z);
 count=0;
 for (unsigned ir=0;ir<nplot_r;ir++)
  {
   x[0]=Global_Physical_Variables::Inner_radius+
    double(ir)/double(nplot_r-1)*(Global_Physical_Variables::Outer_radius-
                                  Global_Physical_Variables::Inner_radius);
   for (unsigned iz=0;iz<nplot_z;iz++)
    {
     x[1]=double(iz)/double(nplot_z-1)*Global_Physical_Variables::Length;

     // Pointer to GeomObject that contains this point
     GeomObject* geom_obj_pt=0;

     // Get it
     Solid_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);

     // Store it
     Solid_regularly_spaced_plot_point[count].first=
      dynamic_cast<POROELASTICITY_ELEMENT*>(geom_obj_pt);
     Solid_regularly_spaced_plot_point[count].second=s;

     count++;
    }
  }

 Poro_mesh_pt->output_boundaries("poro_mesh_boundaries.dat");


 //------------------------------------------------------------------------

 // Create surface meshes
 Inflow_fluid_surface_mesh_pt = new Mesh;
 FSI_fluid_surface_mesh_pt = new Mesh;
 FSI_poro_surface_mesh_pt = new Mesh;

 // Set Id for Lagrange multiplier that imposes (no-)slip on fluid
 Lagrange_id=99;

 // Make surface meshes
 create_fluid_traction_elements();
 create_poro_face_elements();

 //------------------------------------------------------------------------

 // Loop over the elements in fluid mesh to set physical parameters etc
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Cast to the particular element type
   FLUID_ELEMENT *el_pt = dynamic_cast<FLUID_ELEMENT*>
    (Fluid_mesh_pt->element_pt(e));
   
   // There is no need for ALE
   el_pt->disable_ALE();
   
   // Set the Reynolds number for each element
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the product of Reynolds and Strouhal numbers
   el_pt->re_st_pt() = &Global_Physical_Variables::Wo;

   // Pin u_theta at all nodes (no swirl velocity)
   unsigned n_node = el_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node *nod_pt = el_pt->node_pt(n);
     nod_pt->pin(2);
    }
  }

 // Do edge sign setup for darcy elements
 Global_Physical_Variables::edge_sign_setup<POROELASTICITY_ELEMENT>(
  Poro_mesh_pt);

 // Loop over the elements in poroelasticity mesh for setup
 n_element = Poro_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Cast to a bulk element
   POROELASTICITY_ELEMENT *el_pt = dynamic_cast<POROELASTICITY_ELEMENT*>
    (Poro_mesh_pt->element_pt(e));

   // Set the pointer to Poisson's ratio
   el_pt->nu_pt() = &Global_Physical_Variables::Nu;

   // Set the pointer to non-dim Young's modulus
   el_pt->youngs_modulus_pt() = &Global_Physical_Variables::E;

   // Set the pointer to the Lambda parameter
   el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   // Set the pointer to the density ratio
   el_pt->density_ratio_pt()=&Global_Physical_Variables::Density_ratio_poro;

   // Set the pointer to permeability
   el_pt->permeability_pt()=&Global_Physical_Variables::Permeability;

   // Set the pointer to Biot parameter
   el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;

   // Set the pointer to the porosity
   el_pt->porosity_pt()=&Global_Physical_Variables::Porosity;

   // Set the internal q dofs' timestepper to the problem timestepper
   el_pt->set_q_internal_timestepper(time_stepper_pt());

   // Switch off Darcy flow?
   if (CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
    {
     el_pt->switch_off_darcy();
    }
  } // end_loop_over_elements
 
 
 // Set the fluid boundary conditions
 //----------------------------------

 // Pin the boundary nodes
 unsigned n_boundary = Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary-1;b++) // ignore last one -- it's the 
                                      // internal boundary
  {
   // Pin radial velocity at in and outflow (parallel flow)
   // and on centreline (symm); leave FSI boundary free
   // [Doing the loop like this ensures that we pin nodes
   // that are simultaneously on FSI boundary and in/outflow boundaries]
   if (b!=FSI_fluid_boundary)
    {
     unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;++n)
      {
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,n);
       nod_pt->pin(0);


       // If node is also on FSI boundary it must be on in/outflow:
       // Pin axial velocity too
       if (nod_pt->is_on_boundary(FSI_fluid_boundary))
        {
         nod_pt->pin(1);
        }
      }
    }
  }

 // Set the poroelasticity boundary conditions
 //-------------------------------------------

 // Pin the boundary nodes
 n_boundary = Poro_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary;b++)
  {
   if((b==Bottom_poro_boundary)||(b==Top_poro_boundary)) 
    {
     unsigned n_boundary_node = Poro_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       Node* nod_pt=Poro_mesh_pt->boundary_node_pt(b,n);
       
       // pin the solid displacement
       nod_pt->pin(0);
       nod_pt->pin(1);


       // Pin the normal component of the porous flux (at those nodes 
       // that store it). Bit hacky...
       if(nod_pt->nvalue()==4)
        {
         nod_pt->pin(2);
         nod_pt->pin(3);
        }
      }
    }
  }

 // Pin theta Lagrange multipliers (no swirl) everywhere and radial and 
 // axial Lagrange multiplier for nodes on in and outflow boundaries
 // because velocity is already constrained there
 unsigned nel=FSI_fluid_surface_mesh_pt->nelement();
 for(unsigned e=0;e<nel;e++)
  {
   // Get face element
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,POROELASTICITY_ELEMENT>*
    traction_element_pt=dynamic_cast<LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,POROELASTICITY_ELEMENT>*>
    (FSI_fluid_surface_mesh_pt->element_pt(e));

   unsigned n_node=traction_element_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node *nod_pt=traction_element_pt->node_pt(n);

     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);

     // Get the index of the first nodal value associated with
     // this FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(Lagrange_id);

     // Pin swirl component of Lagrange multiplier
     nod_pt->pin(first_index+2);

     // Pin radial and axial component if on in or outflow boundary
     if (bnod_pt->is_on_boundary(Bottom_fluid_boundary))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
     if (bnod_pt->is_on_boundary(Top_fluid_boundary))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
    }
  }

 // Setup fluid structure interaction
 setup_fsi();

 // Add submeshes and build global mesh
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Poro_mesh_pt);
 add_sub_mesh(Inflow_fluid_surface_mesh_pt);
 add_sub_mesh(FSI_fluid_surface_mesh_pt);
 add_sub_mesh(FSI_poro_surface_mesh_pt);
 build_global_mesh();

 //Setup all the equation numbering and look-up schemes
 oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

} // end_of_constructor


//========================================================================
/// Create the traction elements to the appropriate boundaries of the fluid mesh
//========================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, POROELASTICITY_ELEMENT>::
create_fluid_traction_elements()
{

 // Inflow: 
 //--------
 unsigned bound=Top_fluid_boundary;

 // Now loop over bulk elements and create the face elements
 unsigned nel=Fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies pressure boundary condition
   // at inflow
   AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>* traction_element_pt
    = new AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>
    (Fluid_mesh_pt->boundary_element_pt(bound,e),
     Fluid_mesh_pt->face_index_at_boundary(bound,e));

   // Add to mesh
   Inflow_fluid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the applied traction
   traction_element_pt->traction_fct_pt() =
    &Global_Physical_Variables::fluid_inflow_boundary_traction;
  }

 // FSI boundary: 
 //--------------
 bound=FSI_fluid_boundary;

 // Now loop over bulk elements and create the face elements
 nel=Fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes [no] slip (Beavers Joseph Saffman)
   // condition at FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,POROELASTICITY_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,POROELASTICITY_ELEMENT>
    (Fluid_mesh_pt->boundary_element_pt(bound,e),
     Fluid_mesh_pt->face_index_at_boundary(bound,e),Lagrange_id);

   // Add to mesh
   FSI_fluid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);

   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;

   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }

}

//========================================================================
/// Create the traction/pressure elements to the appropriate boundaries
/// of the poroelasticity mesh
//========================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, POROELASTICITY_ELEMENT>::
create_poro_face_elements()
{

 // FSI boundary 
 //-------------
 unsigned bound=FSI_poro_boundary;

 // Now loop over bulk elements and create the face elements
 unsigned nel=Poro_mesh_pt->nboundary_element(bound);
 for(unsigned n=0;n<nel;n++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <POROELASTICITY_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <POROELASTICITY_ELEMENT,FLUID_ELEMENT>
    (Poro_mesh_pt->boundary_element_pt(bound,n),
     Poro_mesh_pt->face_index_at_boundary(bound,n));

   // Add to mesh
   FSI_poro_surface_mesh_pt->add_element_pt(face_element_pt);

   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;

   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(bound);
  }

} // end of create poro face elements


//=====================start_of_setup_fsi================================
/// Setup interaction between two fields
//========================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, POROELASTICITY_ELEMENT>::
setup_fsi()
{
 // Setup fluid traction on FSI wall elements
 // ------------------------------------------
 unsigned fsi_boundary_in_fluid_mesh= FSI_fluid_boundary;
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,fsi_boundary_in_fluid_mesh,Fluid_mesh_pt,FSI_poro_surface_mesh_pt);

 // Setup BJS bc for axisym Navier Stokes from adjacent axisym
 // ----------------------------------------------------------
 // poroelasticity elements
 // -----------------------
 unsigned fsi_boundary_in_poro_mesh=FSI_poro_boundary;
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <POROELASTICITY_ELEMENT,2>
  (this,fsi_boundary_in_poro_mesh,Poro_mesh_pt,FSI_fluid_surface_mesh_pt);

}

//==========================================================================
/// Doc solution
//==========================================================================
template<class FLUID_ELEMENT, class POROELASTICITY_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, POROELASTICITY_ELEMENT>::
doc_solution(DocInfo& doc_info)
{


 oomph_info << "Outputting step " << doc_info.number() 
            << " for time " << time_pt()->time() << std::endl;

 // Define a string that we can set to be the name of the output file
 char filename[100];

 // Define an output filestream
 ofstream file;

 // Fluid
 sprintf(filename,"%s/soln-fluid%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Fluid_mesh_pt->output(file,5);
 file.close();

 // Output solution to file in paraview format
 string pad="";
 if (doc_info.number()<10)
  {
   pad="0000";
  }
 else if (doc_info.number()<100)
  {
   pad="000";
  }
 else if (doc_info.number()<1000)
  {
   pad="00";
  }
 else if (doc_info.number()<10000)
  {
   pad="0";
  }

 sprintf(filename,"%s/soln-fluid%s%i.vtu",doc_info.directory().c_str(),
         pad.c_str(),doc_info.number());
 file.open(filename);
 Fluid_mesh_pt->output_paraview(file,5);
 file.close();
 

 // Fluid at regularly spaced points
 sprintf(filename,"%s/regular_fluid%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 unsigned npts=Fluid_regularly_spaced_plot_point.size();
 Vector<double> x(2);
 Vector<double> s(2);
 Vector<double> veloc(3);
 for (unsigned i=0;i<npts;i++)
  {
   // Pointer to element
   FLUID_ELEMENT* el_pt=Fluid_regularly_spaced_plot_point[i].first;
   
   // Coordinates in it
   s=Fluid_regularly_spaced_plot_point[i].second;
   
   // Get coords
   el_pt->interpolated_x(s,x);

   // Get velocity
   el_pt->interpolated_u_axi_nst(s,veloc);

   file << x[0] << " "
        << x[1] << " "
        << veloc[0] << " "
        << veloc[1] << " "
        << veloc[2] << " "
        << el_pt->interpolated_p_axi_nst(s) << " "
        << std::endl;
  }
 file.close();

 // Fluid inflow traction
 sprintf(filename,"%s/fluid_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Inflow_fluid_surface_mesh_pt->output(file,5);
 file.close();

 // Fluid FSI traction
 sprintf(filename,"%s/fluid_fsi_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 FSI_fluid_surface_mesh_pt->output(file,5);
 file.close();

 // Poroelasticity
 sprintf(filename,"%s/soln-poro%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Poro_mesh_pt->output(file,5);
 file.close();


 // Output solution to file in paraview format
 sprintf(filename,"%s/soln-poro%s%i.vtu",
         doc_info.directory().c_str(),
         pad.c_str(),
         doc_info.number());
 file.open(filename);
 Poro_mesh_pt->output_paraview(file,5);
 file.close();
 
 // Solid at regularly spaced points
 sprintf(filename,"%s/regular_poro%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 npts=Solid_regularly_spaced_plot_point.size();
 for (unsigned i=0;i<npts;i++)
  {
   // Pointer to element
   POROELASTICITY_ELEMENT* el_pt=Solid_regularly_spaced_plot_point[i].first;
   
   // Coordinates in it
   s=Solid_regularly_spaced_plot_point[i].second;
   
   // Get coords
   el_pt->interpolated_x(s,x);

   file << x[0] << " "
        << x[1] << " "
        << el_pt->interpolated_u(s,0) << " "
        << el_pt->interpolated_u(s,1) << " "
        << el_pt->interpolated_q(s,0) << " "
        << el_pt->interpolated_q(s,1) << " "
        << el_pt->interpolated_p(s) << " "
        << std::endl;
  }
 file.close();

 // Poroelasticity FSI traction
 sprintf(filename,"%s/poro_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 FSI_poro_surface_mesh_pt->output(file,5);
 file.close();

 doc_info.number()++;

}

//================================================================
// Driver code
//================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Reynolds number
 CommandLineArgs::specify_command_line_flag("--re",
                                            &Global_Physical_Variables::Re);

 // FSI parameter
 CommandLineArgs::specify_command_line_flag("--q",
                                            &Global_Physical_Variables::Q);

 // Drained Poisson's ratio
 CommandLineArgs::specify_command_line_flag("--nu",
                                            &Global_Physical_Variables::Nu);

 // Initial inlet pressure
 CommandLineArgs::specify_command_line_flag(
  "--p_inlet_initial",
  &Global_Physical_Variables::P_inlet_initial);

 // Step in inlet pressure (applied with tanh profile)
 CommandLineArgs::specify_command_line_flag(
  "--p_inlet_step",
  &Global_Physical_Variables::P_inlet_step);

 // Non-dim permeability
 CommandLineArgs::specify_command_line_flag(
  "--permeability",
  &Global_Physical_Variables::Permeability);

 // Ratio of poro-elastic fluid to Navier-Stokes fluid densities
 CommandLineArgs::specify_command_line_flag(
  "--fluid_density_ratio",
  &Global_Physical_Variables::Poro_elastic_to_navier_stokes_fluid_density);

 // Biot parameter
 CommandLineArgs::specify_command_line_flag(
  "--alpha",
  &Global_Physical_Variables::Alpha);

 // Porosity
 CommandLineArgs::specify_command_line_flag(
  "--porosity",
  &Global_Physical_Variables::Porosity);

 // Factor for reduction in inlet pressure
 CommandLineArgs::specify_command_line_flag(
  "--chris_p_factor", &Global_Physical_Variables::Chris_p_factor);

 // Outer radius
 CommandLineArgs::specify_command_line_flag(
  "--outer_radius",&Global_Physical_Variables::Outer_radius);

 // Outer radius
 CommandLineArgs::specify_command_line_flag(
  "--length",&Global_Physical_Variables::Length);

 // Pin darcy?
 CommandLineArgs::specify_command_line_flag("--pin_darcy");

 // Start from rest?
 CommandLineArgs::specify_command_line_flag("--start_from_steady_presolve");

 // Element area for fluid mesh
 CommandLineArgs::specify_command_line_flag(
  "--el_area_fluid",
  &Global_Physical_Variables::Element_area_fluid);

 // Element area for solid mesh
 CommandLineArgs::specify_command_line_flag(
  "--el_area_solid",
  &Global_Physical_Variables::Element_area_solid);

 // Output directory
 CommandLineArgs::specify_command_line_flag(
  "--dir",&Global_Physical_Variables::Directory);

 // Timestep
 double dt=0.00001;
 CommandLineArgs::specify_command_line_flag("--dt",&dt);

 // Number of timesteps to perform
 unsigned nstep=1000;
 CommandLineArgs::specify_command_line_flag("--nstep",&nstep);

 // Steepness parameter for tanh increment of pressure load 
 CommandLineArgs::specify_command_line_flag(
  "--alpha_tanh",
  &Global_Physical_Variables::Alpha_tanh);
 
 // "Time" for tanh increment of pressure load 
 CommandLineArgs::specify_command_line_flag(
  "--t_tanh",
  &Global_Physical_Variables::T_tanh);

 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create doc info object
 DocInfo doc_info;
 doc_info.set_directory(Global_Physical_Variables::Directory);

 // Create problem
 PressureWaveFSIProblem<AxisymmetricTTaylorHoodElement,
                         TAxisymmetricPoroelasticityElement<1> > problem;

 // Update dependent problem parameters
 Global_Physical_Variables::update_dependent_parameters();

 // Doc dependent parameters
 Global_Physical_Variables::doc_dependent_parameters();

 // Set up impulsive start from rest
 problem.assign_initial_values_impulsive(dt);

 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   // Do only 3 steps if validating
   nstep=3;
  }

 // Output current time
 oomph_info << "Output initial state.\n";

 // Output mesh and initial conditions
 problem.doc_solution(doc_info);


 if (CommandLineArgs::command_line_flag_has_been_set
     ("--start_from_steady_presolve"))
  {
   oomph_info << "Doing steady initial solve.\n";
   
   // Do a steady solve
   problem.steady_newton_solve();
   
   // Output solution
   problem.doc_solution(doc_info);
  }
 else
  {
   oomph_info << "Bypassing steady initial solve and starting from rest.\n";   
  }


 // Timestep
 oomph_info << "About to do " << nstep << " timesteps with dt = "
            << dt << std::endl;
 for(unsigned n=0;n<nstep;n++)
  {
   oomph_info << "Solving at time " << problem.time_pt()->time()+dt 
              << std::endl;

   //Solve the problem
   problem.unsteady_newton_solve(dt);

   // Output
   problem.doc_solution(doc_info);
  }

 return 0;
}

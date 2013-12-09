//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.90. August 3, 2009.
//LIC//
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#include <fenv.h>

#include "axisym_poroelasticity.h"
#include "axisym_navier_stokes.h"
#include "meshes/triangle_mesh.h"

using namespace oomph;
using namespace std;

//===Start_of_Global_Physical_Variables_namespace======================
/// Namespace for global parameters
//======================================================================
namespace Global_Physical_Variables
{
 /// Output directory
 std::string Directory = "RESLT";

 /// Global FSI parameter
 double Q = 1.0e-9;

 // Fluid parameters
 // ----------------

 /// The Reynolds number
 double Re = 10.0;

 /// The Strouhal number
 double St = 1.0;

 /// \short The Womersley number -- dependent parameter; compute
 /// from Re and St
 double Wo = 0.0;

  /// Factor for wall inertia parameter (based on dura properties)
 double Lambda_factor=1.414213562e-6; 
 
 /// Factor for FSI parameter (based on dura properties)
 double Q_factor=2.0e-12;
 
 /// Factor for non-dim permeability (based on dura properties)
 /// Note: Is multiplied by Permeability_multiplier before
 /// working out actual numbers)
 double K_factor=125.0; 

 /// Multiplier for permeability
 double Permeability_multiplier=1.0;

 // Dura and pia properties
 //------------------------

 /// \short Stiffness ratio for dura and pia: ratio of actual Young's modulus
 /// to Young's modulus used in the non-dim of the equations (unity because 
 /// we've scaled things on their properties)
 double Stiffness_ratio_dura_and_pia=1.0;

 /// \short Permeability ratio; ratio of actual permeability
 /// to permeability used in the non-dim of the equations (unity because 
 /// we've scaled things on their properties)
 double Permeability_ratio_dura_and_pia=1.0;

 /// Dura and Pia Poisson's ratio for drained poro-elastic medium
 double Nu_dura_and_pia = 0.35;

 /// Porosity for dura and pia
 double Porosity_dura_and_pia = 0.3;


 // Cord properties
 //----------------

 /// \short Stiffness ratio for cord: ratio of actual Young's modulus
 /// to Young's modulus used in the non-dim of the equations
 double Stiffness_ratio_cord=0.004;

 /// \short Permeability ratio for cord; ratio of actual permeability
 /// to permeability used in the non-dim of the equations
 double Permeability_ratio_cord=1.0;

 /// Cord Poisson's ratio for drained poro-elastic medium
 double Nu_cord = 0.35;

 /// Porosity for cord
 double Porosity_cord = 0.3;

 // Filum and block properties
 //---------------------------

 /// \short Stiffness ratio for filum and block: ratio of actual Young's modulus
 /// to Young's modulus used in the non-dim of the equations
 double Stiffness_ratio_filum_and_block=0.050;

 /// \short Permeability ratio for filum and block; ratio of actual 
 /// permeability to permeability used in the non-dim of the equations
 double Permeability_ratio_filum_and_block=1.0;

 /// Filum and block Poisson's ratio for drained poro-elastic medium
 double Nu_filum_and_block = 0.35;

 /// Porosity for filum and block
 double Porosity_filum_and_block= 0.3;


 // Generic poroelasticity parameters
 // ---------------------------------

 /// \short The poroelasticity inertia parameter -- dependent parameter;
 /// compute from Re lambda factor
 double Lambda_sq = 0;

 /// \short Non-dim permeability -- ratio of typical porous flux to fluid veloc
 /// scale -- dependent parameter; compute from Re and k factor
 double Permeability = 0.0;

 // Alpha, the Biot parameter (same for all; incompressible fluid and solid)
 double Alpha = 1.0;

 /// \short Ratio of the densities of the fluid and solid phases of
 /// the poroelastic material (same for all)
 double Density_ratio = 1.0;

 /// Inverse slip rate coefficient
 double Inverse_slip_rate_coefficient = 0.0;

 // /// \short Ratio of the pore fluid density to the compound density --
 // /// dependent parameter compute from density ratio and porosity.
 // double Rho_f_over_rho = 0.0;

 // Mesh/problem parameters
 // -----------------------

 /// Target element area for fluid meshes
 double Element_area_fluid = 0.001;

 /// Target element area for poro-elastic meshes
 double Element_area_poro = 0.01;

 /// Raw boundary layer thickness (in Chris' original coordinates)
 double BL_thick=0.1;

 /// \short Moens Korteweg wavespeed based on dura properties -- dependent
 /// parameter
 double Dura_wavespeed=0.0;

 /// \short Overall scaling factor make inner diameter of outer porous
 /// region at inlet equal to one so Re etc. mean what they're
 /// supposed to mean
 double Length_scale = 1.0/20.0;

 /// Radial scaling -- vanilla length scale
 double R_scale = Length_scale;
 
 /// Aspect ratio factor to shrink domain in z direction
 double Z_shrink_factor=1.0;

 /// Axial scaling -- dependent parameter; update
 double Z_scale = Length_scale/Z_shrink_factor;

 /// Suppress regularly spaced output
 bool Suppress_regularly_spaced_output=false;
 
 /// Regular output with approx 3 per thickness of dura
 unsigned Nplot_r_regular=33; 
 
 /// Regular output with equivalent axial spacing
 unsigned Nplot_z_regular=
  unsigned(double(Global_Physical_Variables::Nplot_r_regular)*600.0/11.0)
  /Global_Physical_Variables::Z_shrink_factor;
  

 //-----------------------------------------------------------------------------

 /// \short Helper function to update dependent parameters (based on the
 /// linearised_poroelastic_fsi_pulsewave version)
 void update_dependent_parameters()
  {
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--suppress_wall_inertia"))
    {
     Lambda_factor=0.0;
     oomph_info <<"\n\nNOTE: Set wall inertial to zero!\n\n";
    }

   // The Womersley number
   Wo = Re*St;

   // (Square of) inertia parameter
   Lambda_sq = pow(Lambda_factor*Re,2);
   
   // FSI parameter
   Q = Q_factor*Re;
   
   // Non-dim permeability
   Permeability=K_factor/Re*Permeability_multiplier;

   // Wavespeed in dura (constant thickness and radius as at inlet)
   double non_dim_thickness=1.0/20.0; // from Chris' sketch -- inner diameter
                                      // 20mm; wall thickness 1mm
   //Wavespeed
   if ((Q*Re)!=0.0)
    {
     Dura_wavespeed=sqrt(non_dim_thickness/(Q*Re));
    }
   else
    {
     Dura_wavespeed=999.9e99;
    }
  }

 /// Doc dependent parameters
 void doc_dependent_parameters()
  {
   oomph_info << std::endl;
   oomph_info << "Global problem parameters" << std::endl;
   oomph_info << "=========================" << std::endl;
   oomph_info << "Reynolds number (Re)                           : "
              << Re << std::endl;
   oomph_info << "Strouhal number (St)                           : "
              << St << std::endl;
   oomph_info << "Womersley number (ReSt)                        : "
              << Wo << std::endl;
   oomph_info << "FSI parameter (Q)                              : "
              << Q << std::endl << std::endl;

   oomph_info << "Factor for wall inertia parameter              : "
              << Lambda_factor << std::endl;
   oomph_info << "Factor for FSI parameter                       : "
              << Q_factor << std::endl;
   oomph_info << "Factor for non-dim permeability                : " 
              << K_factor << std::endl << std::endl;

   oomph_info << "Extra multiplier for non-dim permeability      : " 
              << Permeability_multiplier << std::endl << std::endl;

   oomph_info << "Dura wavespeed                                 : " 
              << Dura_wavespeed << std::endl << std::endl;

   oomph_info << "\n\nCommon poro-elastic properties:\n";
   oomph_info <<     "===============================\n\n";

   oomph_info << "(Square of) wall inertia parameter (Lambda^2)  : "
              << Lambda_sq << std::endl;
   oomph_info << "Non-dim permeability                           : "
              << Permeability << std::endl;
   oomph_info << "Biot parameter (alpha) in all porous media     : "
              << Alpha << std::endl;
   oomph_info << "Density ratio (rho_f/rho_s) in all porous media: "
              << Density_ratio << std::endl;
   oomph_info << "Inverse slip rate coefficient                  : "
              << Inverse_slip_rate_coefficient << std::endl;

   oomph_info << "\n\nDura and pia properties:\n";
   oomph_info <<     "========================\n\n";
   oomph_info << "Stiffness ratio (non-dim Young's modulus)      : "
              << Stiffness_ratio_dura_and_pia << std::endl;
   oomph_info << "Permeability ratio                             : "
              << Permeability_ratio_dura_and_pia << std::endl;
   oomph_info << "Non-dim permeability                           : "
              << Permeability*Permeability_ratio_dura_and_pia << std::endl;
   oomph_info << "Poisson's ratio for drained material           : "
              << Nu_dura_and_pia << std::endl;
   oomph_info << "Porosity                                       : "
              << Porosity_dura_and_pia <<std::endl;

   oomph_info << "\n\nCord properties:\n";
   oomph_info <<     "================\n\n";
   oomph_info << "Stiffness ratio (non-dim Young's modulus)      : "
              << Stiffness_ratio_cord << std::endl;
   oomph_info << "Permeability ratio                             : "
              << Permeability_ratio_cord << std::endl;
   oomph_info << "Non-dim permeability                           : "
              << Permeability*Permeability_ratio_cord << std::endl;
   oomph_info << "Poisson's ratio for drained material           : "
              << Nu_cord << std::endl;
   oomph_info << "Porosity                                       : "
              << Porosity_cord <<std::endl;

   oomph_info << "\n\nFilum and block properties:\n";
   oomph_info <<     "===========================\n\n";
   oomph_info << "Stiffness ratio (non-dim Young's modulus)      : "
              << Stiffness_ratio_filum_and_block << std::endl;
   oomph_info << "Permeability ratio                             : "
              << Permeability_ratio_filum_and_block << std::endl;
   oomph_info << "Non-dim permeability                           : "
              << Permeability*Permeability_ratio_filum_and_block << std::endl;
   oomph_info << "Poisson's ratio for drained material           : "
              << Nu_filum_and_block << std::endl;
   oomph_info << "Porosity                                       : "
              << Porosity_filum_and_block <<std::endl;
   oomph_info << std::endl << std::endl;
  }

 //-----------------------------------------------------------------------------

 /// Inflow traction applied to the fluid mesh
 void fluid_inflow_boundary_traction(const double &time,
                                     const Vector<double> &x,
                                     const Vector<double> &n,
                                     Vector<double> &result)
  {
   result[0]=Re*n[0];
   result[1]=Re*n[1];
   result[2]=0.0;
  }


 //-----------------------------------------------------------------------------
 
 /// Helper function to create equally spaced vertices between
 /// between final existing vertex and specified end point. 
 /// Vertices are pushed back into Vector of vertices.
 void push_back_vertices(const Vector<double>& end,
                         const double& edge_length,
                         Vector<Vector<double> >& vertex)
 {
  Vector<double> start(2);
  double n_existing=vertex.size();
  start[0]=vertex[n_existing-1][0];
  start[1]=vertex[n_existing-1][1];
  double total_length=sqrt(pow(end[0]-start[0],2)+
                           pow(end[1]-start[1],2));
  unsigned n_seg=unsigned(total_length/edge_length);
  Vector<double> x(2);
  // Skip first one!
  for (unsigned j=1;j<n_seg;j++)
   {
    x[0]=start[0]+(end[0]-start[0])*double(j)/double(n_seg-1);
    x[1]=start[1]+(end[1]-start[1])*double(j)/double(n_seg-1);
    vertex.push_back(x);
   }
 }

 /// Helper function to create equally spaced vertices between
 /// between start and end points. Vertices are pushed back
 /// into Vector of vertices.
 void push_back_vertices(const Vector<double>& start,
                         const Vector<double>& end,
                         const double& edge_length,
                         Vector<Vector<double> >& vertex)
 {
  double total_length=sqrt(pow(end[0]-start[0],2)+
                           pow(end[1]-start[1],2));
  unsigned n_seg=unsigned(total_length/edge_length);
  Vector<double> x(2);
  for (unsigned j=0;j<n_seg;j++)
   {
    x[0]=start[0]+(end[0]-start[0])*double(j)/double(n_seg-1);
    x[1]=start[1]+(end[1]-start[1])*double(j)/double(n_seg-1);
    vertex.push_back(x);
   }
 }

 //-----------------------------------------------------------------------------
 
 /// Helper function to get intersection between two lines
 Vector<double> intersection(const Vector<Vector<double> >& line1,
                             const Vector<Vector<double> >& line2)
 {
  Vector<double> intersect(2);
  double a2 =  line2[1][1] - line2[0][1];
  double b2 = -line2[1][0] + line2[0][0];
  double c2 = a2*line2[0][0]+b2*line2[0][1];
  
  double a1 =  line1[1][1] - line1[0][1];
  double b1 = -line1[1][0] + line1[0][0];
  double c1 = a1*line1[0][0]+b1*line1[0][1];
  
  double det = a1*b2 - a2*b1;
  if (abs(det)<1.0e-16)
   {
    oomph_info << "Trouble";
    abort();
   }
  intersect[0]= (b2*c1-b1*c2)/det;
  intersect[1]= (a1*c2-a2*c1)/det;
  return intersect;
  
 }

 /// Helper function to get intersection between boundary layer
 /// to the right of the two lines
 Vector<double> bl_intersection(const Vector<Vector<double> >& line1,
                                const Vector<Vector<double> >& line2,
                                const double& bl_thick)
 {
  // Get normal to right for first line
  Vector<double> normal(2);
  normal[0]=line1[1][1]-line1[0][1];
  normal[1]=-(line1[1][0]-line1[0][0]);
  double norm=normal[0]*normal[0]+normal[1]*normal[1];
  normal[0]/=sqrt(norm);
  normal[1]/=sqrt(norm);

  // Get boundary layer line
  Vector<Vector<double> >bl_line1(2,Vector<double>(2));
  bl_line1[0][0]=line1[0][0]+bl_thick*normal[0];
  bl_line1[0][1]=line1[0][1]+bl_thick*normal[1]*
   Global_Physical_Variables::Z_shrink_factor;
  bl_line1[1][0]=line1[1][0]+bl_thick*normal[0];
  bl_line1[1][1]=line1[1][1]+bl_thick*normal[1]*
   Global_Physical_Variables::Z_shrink_factor;

  // Get normal to right for second line
  normal[0]=line2[1][1]-line2[0][1];
  normal[1]=-(line2[1][0]-line2[0][0]);
  norm=normal[0]*normal[0]+normal[1]*normal[1];
  normal[0]/=sqrt(norm);
  normal[1]/=sqrt(norm);

  // Get boundary layer line
  Vector<Vector<double> >bl_line2(2,Vector<double>(2));
  bl_line2[0][0]=line2[0][0]+bl_thick*normal[0];
  bl_line2[0][1]=line2[0][1]+bl_thick*normal[1]*
   Global_Physical_Variables::Z_shrink_factor;
  bl_line2[1][0]=line2[1][0]+bl_thick*normal[0];
  bl_line2[1][1]=line2[1][1]+bl_thick*normal[1]*
   Global_Physical_Variables::Z_shrink_factor;

  Vector<double> intersect=intersection(bl_line1,bl_line2);

  return intersect;
 }

 //-----------------------------------------------------------------------------

 /// \short Global function that completes the edge sign setup
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

} // end_of_Global_Physical_Variables_namespace

//===start_of_problem_class=============================================
/// Problem class
//======================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
class AxisymmetricSpineProblem : public Problem
{
public:

 /// Constructor
 AxisymmetricSpineProblem();

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 /// Actions before implicit timestep
 void actions_before_implicit_timestep() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,unsigned npts=5);

private:

 /// Create the poroelasticity traction/pressure elements
 void create_poro_face_elements();

 /// Create the fluid traction elements
 void create_fluid_traction_elements();

 /// Setup FSI
 void setup_fsi();

 /// Pointers to the inner poro mesh
 Mesh* Inner_poro_mesh_pt;

 /// Pointers to the outer poro mesh
 Mesh* Outer_poro_mesh_pt;

 /// Pointer to the syrinx fluid mesh
 Mesh* Syrinx_fluid_mesh_pt;

 /// Pointer to the SSS fluid mesh
 Mesh* SSS_fluid_mesh_pt;

 /// Inflow fluid surface mesh
 Mesh* Inflow_fluid_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to inner poro mesh from syrinx
 /// fluid. Attaches to inner poro mesh.
 Mesh* Inner_poro_syrinx_fluid_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to syrinx mesh from inner
 /// poro. Attaches to syrinx fluid mesh.
 Mesh* Syrinx_fluid_inner_poro_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to inner poro mesh from SSS
 /// fluid. Attaches to inner poro mesh.
 Mesh* Inner_poro_SSS_fluid_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to SSS mesh from inner
 /// poro. Attaches to SSS fluid mesh.
 Mesh* SSS_fluid_inner_poro_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to outer poro mesh from SSS
 /// fluid near inlet. Attaches to outer poro mesh 
 Mesh* Outer_poro_SSS_fluid_near_inlet_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to SSS mesh from outer
 /// poro near inlet. Attaches to SSS fluid mesh.
 Mesh* SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to outer poro mesh from SSS
 /// fluid near outlet. Attaches to outer poro mesh 
 Mesh* Outer_poro_SSS_fluid_near_outlet_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to SSS mesh from outer
 /// poro near outlet. Attaches to SSS fluid mesh.
 Mesh* SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to outer poro mesh from SSS
 /// fluid. Attaches to outer poro mesh [only block!]
 Mesh* Outer_poro_SSS_fluid_block_FSI_surface_mesh_pt;

 /// \short FSI mesh to apply traction/pressure to SSS mesh from outer
 /// poro. Attaches to SSS fluid mesh. [only block!]
 Mesh* SSS_fluid_outer_poro_block_FSI_surface_mesh_pt;

 /// \short Mesh as geom object representation of SSS fluid mesh
 MeshAsGeomObject* SSS_fluid_mesh_geom_obj_pt;
 
 /// \short Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<FLUID_ELEMENT*,Vector<double> > > 
 SSS_fluid_regularly_spaced_plot_point;

 /// \short Mesh as geom object representation of syrinx fluid mesh
 MeshAsGeomObject* Syrinx_fluid_mesh_geom_obj_pt;
 
 /// \short Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<FLUID_ELEMENT*,Vector<double> > > 
 Syrinx_fluid_regularly_spaced_plot_point;

 /// \short Mesh as geom object representation of outer poro mesh
 MeshAsGeomObject* Outer_poro_mesh_geom_obj_pt;
 
 /// \short Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<PORO_ELEMENT*,Vector<double> > > 
 Outer_poro_regularly_spaced_plot_point;

 /// \short Mesh as geom object representation of inner poro mesh
 MeshAsGeomObject* Inner_poro_mesh_geom_obj_pt;
 
 /// \short Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<PORO_ELEMENT*,Vector<double> > > 
 Inner_poro_regularly_spaced_plot_point;

 /// Pointer to the poroelasticity timestepper
 TimeStepper* Poro_time_stepper_pt;

 /// Pointer to the fluid timestepper
 TimeStepper* Fluid_time_stepper_pt;

 /// \short Enumeration for Lagrange multiplier IDs for Lagrange multipliers 
 /// that enforce BJS boundary condition on fluid from poro meshes
 enum
 {
  Lagrange_id_syrinx_fluid_inner_poro,
  Lagrange_id_SSS_fluid_inner_poro,
  Lagrange_id_SSS_fluid_outer_poro
 };


 /// Enumeration of the inner poro boundaries
 enum
  {
   Inner_poro_inlet_boundary_id,
   Inner_poro_SSS_interface_boundary_id,
   Inner_poro_outlet_boundary_id,
   Inner_poro_symmetry_near_outlet_boundary_id,
   Inner_poro_syrinx_interface_boundary_id,
   Inner_poro_symmetry_near_inlet_boundary_id,
   Inner_poro_pia_internal_boundary_id,
   Inner_poro_filum_internal_boundary_id
  };


 /// Enumeration of the outer poro boundaries
 enum
  {
   Outer_poro_inlet_boundary_id,
   Outer_poro_outer_boundary_id,
   Outer_poro_outlet_boundary_id,
   Outer_poro_SSS_interface_near_outlet_boundary_id,
   Outer_poro_SSS_interface_block_boundary_id,
   Outer_poro_SSS_interface_near_inlet_boundary_id,
   Outer_poro_inner_block_boundary_id
  };

 /// Enumeration of the syrinx fluid boundaries
 enum
  {
   Syrinx_fluid_inner_poro_interface_boundary_id,
   Syrinx_fluid_symmetry_boundary_id,
   Syrinx_bl_boundary_id
  };

 /// Enumeration of the SSS fluid boundaries
 enum
  {
   SSS_fluid_outer_poro_interface_block_boundary_id,
   SSS_fluid_outer_poro_interface_near_outlet_boundary_id,
   SSS_fluid_outflow_boundary_id,
   SSS_fluid_inner_poro_interface_boundary_id,
   SSS_fluid_inflow_boundary_id,
   SSS_fluid_outer_poro_interface_near_inlet_boundary_id,
   SSS_bl_near_inlet_boundary_id,
   SSS_bl_near_outlet_boundary_id
  };


 /// Enumeration of regions
 enum
  {
   Cord_region_id, // Note: not specified therefore implied!
   Dura_region_id,
   Block_region_id,
   Pia_region_id,
   Filum_region_id
  };


}; // end_of_problem_class

//===start_of_constructor=============================================
/// Problem constructor:
//====================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
AxisymmetricSpineProblem<PORO_ELEMENT, FLUID_ELEMENT>::
AxisymmetricSpineProblem()
{
 // Create timesteppers
 Poro_time_stepper_pt = new Newmark<2>;
 add_time_stepper_pt(Poro_time_stepper_pt);
 Fluid_time_stepper_pt = new BDF<2>;
 add_time_stepper_pt(Fluid_time_stepper_pt);

 // Create the problem geometry - a polygonal model of the spinal cavity with
 // -------------------------------------------------------------------------
 // poroelastic walls filled with fluid
 // -----------------------------------

 // Update axial scale factor
 Global_Physical_Variables::Z_scale = 
  Global_Physical_Variables::Length_scale/
  Global_Physical_Variables::Z_shrink_factor;
 
 // Update regular output with equivalent axial spacing
 Global_Physical_Variables::Nplot_z_regular=
  unsigned(double(Global_Physical_Variables::Nplot_r_regular)*600.0/11.0)
  /Global_Physical_Variables::Z_shrink_factor;
 if (!Global_Physical_Variables::Suppress_regularly_spaced_output)
 {
  oomph_info << "Regular output with: Nplot_r_regular Nplot_z_regular " 
             << Global_Physical_Variables::Nplot_r_regular << " " 
             << Global_Physical_Variables::Nplot_z_regular << std::endl;
 }


 // Poro mesh (inner part)
 // ---------------------
 Vector<TriangleMeshCurveSection*> inner_poro_outer_polyline_boundary_pt(6);

 // Inlet
 {
  Vector<Vector<double> > inner_coords(3, Vector<double>(2));
  inner_coords[0][0]=0.0;
  inner_coords[0][1]=0.0;
  
  inner_coords[1][0]=5.8;
  inner_coords[1][1]=0.0;
  
  inner_coords[2][0]=6.0;
  inner_coords[2][1]=0.0;

  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }

  // Build
  inner_poro_outer_polyline_boundary_pt[Inner_poro_inlet_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,Inner_poro_inlet_boundary_id);
 }

 // FSI interface with SSS
 {
  Vector<Vector<double> > inner_coords(4, Vector<double>(2));

  inner_coords[0][0]=6.0;
  inner_coords[0][1]=0.0;
  
  inner_coords[1][0]=4.24;
  inner_coords[1][1]=-440.0;
  
  inner_coords[2][0]=1.25;
  inner_coords[2][1]=-480.0;
  
  inner_coords[3][0]=1.25; 
  inner_coords[3][1]=-600.0;
  
  
  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_poro_outer_polyline_boundary_pt[Inner_poro_SSS_interface_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,Inner_poro_SSS_interface_boundary_id);
 }


 // Outlet
 {
  Vector<Vector<double> > inner_coords(3, Vector<double>(2));

  inner_coords[0][0]=1.25;
  inner_coords[0][1]=-600.0;
  
  inner_coords[1][0]=1.05;
  inner_coords[1][1]=-600.0;

  inner_coords[2][0]=0.0;
  inner_coords[2][1]=-600.0;
  
  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_poro_outer_polyline_boundary_pt[Inner_poro_outlet_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,Inner_poro_outlet_boundary_id);
 }



 // Symm line near outlet
 {
  Vector<Vector<double> > inner_coords(3, Vector<double>(2));

  inner_coords[0][0]=0.0;
  inner_coords[0][1]=-600.0;
  
  // Need this one to define boundary with filum!
  inner_coords[1][0]=0.0;
  inner_coords[1][1]=-480.0;

  inner_coords[2][0]=0.0;
  inner_coords[2][1]=-220.0;

  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_poro_outer_polyline_boundary_pt
   [Inner_poro_symmetry_near_outlet_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,
                            Inner_poro_symmetry_near_outlet_boundary_id);
 }


 // FSI interface in syrinx
 {
  Vector<Vector<double> > inner_coords(4, Vector<double>(2));
  
  inner_coords[0][0]=0.0;
  inner_coords[0][1]=-220.0;
  
  inner_coords[1][0]=4.128;
  inner_coords[1][1]=-210.0;
  
  inner_coords[2][0]=4.512;
  inner_coords[2][1]=-90.0;
  
  inner_coords[3][0]=0.0;
  inner_coords[3][1]=-80.0;
  
  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_poro_outer_polyline_boundary_pt
   [Inner_poro_syrinx_interface_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,
                            Inner_poro_syrinx_interface_boundary_id);
 }
 
 
 // Symm line near inlet
 {
  Vector<Vector<double> > inner_coords(2, Vector<double>(2));
  
  inner_coords[0][0]=0.0;
  inner_coords[0][1]=-80.0;
  
  inner_coords[1][0]=0.0;
  inner_coords[1][1]=0.0;
  
  // Rescale
  unsigned n=inner_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    inner_coords[j][0]*=Global_Physical_Variables::R_scale;
    inner_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_poro_outer_polyline_boundary_pt
   [Inner_poro_symmetry_near_inlet_boundary_id] =
   new TriangleMeshPolyLine(inner_coords,
                            Inner_poro_symmetry_near_inlet_boundary_id);

 }
 TriangleMeshClosedCurve* inner_poro_outer_boundary_pt =
  new TriangleMeshClosedCurve(inner_poro_outer_polyline_boundary_pt);

 // We have two internal boundaries
 Vector<TriangleMeshOpenCurve*> inner_inner_open_boundary_pt(2);
 TriangleMeshPolyLine *inner_inner_polyline1_pt = 0;
 TriangleMeshPolyLine *inner_inner_polyline2_pt = 0;
   
 // First internal boundary -- boundary between "pia" and genuine interior
 {
  Vector<Vector<double> > inner_inner_coords(4, Vector<double>(2,0.0));
  
  inner_inner_coords[0][0]=5.8;
  inner_inner_coords[0][1]=0.0;
  
  inner_inner_coords[1][0]=4.04;
  inner_inner_coords[1][1]=-440.0;
  
  inner_inner_coords[2][0]=1.05;
  inner_inner_coords[2][1]=-480.0;
  
  inner_inner_coords[3][0]=1.05;
  inner_inner_coords[3][1]=-600.0;
  
  // Rescale
  for(unsigned i=0;i<4;i++)
   {
    inner_inner_coords[i][0]*=Global_Physical_Variables::R_scale;
    inner_inner_coords[i][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_inner_polyline1_pt =
   new TriangleMeshPolyLine(inner_inner_coords,
                            Inner_poro_pia_internal_boundary_id);
  
  // Connect initial vertex to middle vertex on inlet boundary
  inner_inner_polyline1_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (inner_poro_outer_polyline_boundary_pt[Inner_poro_inlet_boundary_id]),1);

  // Connect final vertex to middle vertex on outlet boundary
  inner_inner_polyline1_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (inner_poro_outer_polyline_boundary_pt[Inner_poro_outlet_boundary_id]),1);
  

  // Store in vector
  Vector<TriangleMeshCurveSection*> inner_inner_curve_section1_pt(1);
  inner_inner_curve_section1_pt[0] = inner_inner_polyline1_pt;
  
  // Create first internal open curve
  inner_inner_open_boundary_pt[0] =
   new TriangleMeshOpenCurve(inner_inner_curve_section1_pt);
 }

 // Second internal boundary -- boundary between "filum" and genuine interior
 {
  Vector<Vector<double> > inner_inner_coords(2, Vector<double>(2,0.0));

  inner_inner_coords[0][0]=0.0;
  inner_inner_coords[0][1]=-480.0;
  
  inner_inner_coords[1][0]=1.05;
  inner_inner_coords[1][1]=-480.0;
  
  // Rescale
  for(unsigned i=0;i<2;i++)
   {
    inner_inner_coords[i][0]*=Global_Physical_Variables::R_scale;
    inner_inner_coords[i][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_inner_polyline2_pt =
   new TriangleMeshPolyLine(inner_inner_coords,
                            Inner_poro_filum_internal_boundary_id);
  
  // Connect initial vertex of separating line with first vertex on 
  // symmetry line
  inner_inner_polyline2_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>(
    inner_poro_outer_polyline_boundary_pt
    [Inner_poro_symmetry_near_outlet_boundary_id]),1);

  // Conect final vertex to second vertex on pia internal bouundary
  inner_inner_polyline2_pt->connect_final_vertex_to_polyline(
   inner_inner_polyline1_pt,2);
  
  // Stick into vector
  Vector<TriangleMeshCurveSection*> inner_inner_curve_section2_pt(1);
  inner_inner_curve_section2_pt[0] = inner_inner_polyline2_pt;
  
  // Build
  inner_inner_open_boundary_pt[1] =
   new TriangleMeshOpenCurve(inner_inner_curve_section2_pt);
  
 }

 // Mesh parameters
 TriangleMeshParameters
  inner_poro_mesh_parameters(inner_poro_outer_boundary_pt);
 inner_poro_mesh_parameters.element_area() =
  Global_Physical_Variables::Element_area_poro;
 
 // Define inner boundaries
 inner_poro_mesh_parameters.internal_open_curves_pt() =
  inner_inner_open_boundary_pt;

 // Define region coordinates
 Vector<double> region_coords(2);

 // Point inside the "pia"
 region_coords[0] = Global_Physical_Variables::R_scale*5.9;
 region_coords[1] = Global_Physical_Variables::Z_scale*(-1.0);
 inner_poro_mesh_parameters.add_region_coordinates(Pia_region_id, 
                                                   region_coords);

 // Point inside the filum
 region_coords[0] = Global_Physical_Variables::R_scale*0.5;
 region_coords[1] = Global_Physical_Variables::Z_scale*(-500.0);
 inner_poro_mesh_parameters.add_region_coordinates(Filum_region_id, 
                                                   region_coords);

 // Build the mesh
 Inner_poro_mesh_pt = new TriangleMesh<PORO_ELEMENT>(
   inner_poro_mesh_parameters, Poro_time_stepper_pt);

 // Mesh as geom object representation of inner poro mesh
 Inner_poro_mesh_geom_obj_pt=new MeshAsGeomObject(Inner_poro_mesh_pt);
 
 // Extract regularly spaced points
 if (!Global_Physical_Variables::Suppress_regularly_spaced_output)
  {
   double t_start = TimingHelpers::timer();
   
   // Speed up search (somewhat hand-tuned; adjust if this takes too long
   // or misses too many points)
   Inner_poro_mesh_geom_obj_pt->max_spiral_level()=4;
   unsigned nx_back= Multi_domain_functions::Nx_bin;
   Multi_domain_functions::Nx_bin=10;
   unsigned ny_back= Multi_domain_functions::Ny_bin;
   Multi_domain_functions::Ny_bin=10;
   
   // Populate it...
   Vector<double> x(2);
   Vector<double> s(2);
   for (unsigned ir=0;ir<Global_Physical_Variables::Nplot_r_regular;ir++)
    {
     // outermost radius of dura
     x[0]=double(ir)/double(Global_Physical_Variables::Nplot_r_regular-1)*0.55;
     for (unsigned iz=0;iz<Global_Physical_Variables::Nplot_z_regular;iz++)
      {
       x[1]=double(iz)/double(Global_Physical_Variables::Nplot_z_regular-1)*
        (-600.0)*Global_Physical_Variables::Z_scale;
       
       // Pointer to GeomObject that contains this point
       GeomObject* geom_obj_pt=0;
       
       // Get it
       Inner_poro_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);
       
       // Store it
       if (geom_obj_pt!=0)
        {
         std::pair<PORO_ELEMENT*,Vector<double> > tmp;
         tmp = std::make_pair(dynamic_cast<PORO_ELEMENT*>(geom_obj_pt),s);
         Inner_poro_regularly_spaced_plot_point.push_back(tmp);
        }
      }
    }

   // Reset number of bins in binning method 
   Multi_domain_functions::Nx_bin=nx_back;
   Multi_domain_functions::Ny_bin=ny_back;
   oomph_info << "Took: " << TimingHelpers::timer()-t_start
              << " sec to setup regularly spaced points for inner poro mesh\n";
  }

 // Output mesh and boundaries
 Inner_poro_mesh_pt->output
  (Global_Physical_Variables::Directory+"/inner_poro_mesh.dat");
 Inner_poro_mesh_pt->output_boundaries
  (Global_Physical_Variables::Directory+"/inner_poro_mesh_boundaries.dat");
 
 //-----------------------------------------------------------------------------

 // Poro mesh (outer part)
 // ----------------------
 Vector<TriangleMeshCurveSection*> outer_poro_outer_polyline_boundary_pt(6);

 // Inlet
 {
  Vector<Vector<double> > outer_coords(2, Vector<double>(2));
  outer_coords[0][0]=10.0;
  outer_coords[0][1]=0.0;
  
  outer_coords[1][0]=11.0;
  outer_coords[1][1]=0.0;
  
  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }

  // Build
  outer_poro_outer_polyline_boundary_pt[Outer_poro_inlet_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,Outer_poro_inlet_boundary_id);
 }

 // Outer boundary
 {

  Vector<Vector<double> > outer_coords(3, Vector<double>(2));

  outer_coords[0][0]=11.0;
  outer_coords[0][1]=0.0;
  
  outer_coords[1][0]=8.0;
  outer_coords[1][1]=-500.0;

  outer_coords[2][0]=3.5;
  outer_coords[2][1]=-600.0;

  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }

  // Build
  outer_poro_outer_polyline_boundary_pt[Outer_poro_outer_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,Outer_poro_outer_boundary_id);
 }


 // Outlet
 {
  Vector<Vector<double> > outer_coords(2, Vector<double>(2));

  outer_coords[0][0]=3.5;
  outer_coords[0][1]=-600.0;

  outer_coords[1][0]=2.5;
  outer_coords[1][1]=-600.0;

  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  outer_poro_outer_polyline_boundary_pt[Outer_poro_outlet_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,Outer_poro_outlet_boundary_id);
 }

 // Inner/FSI near outlet
 {
  Vector<Vector<double> > outer_coords(3, Vector<double>(2));

  outer_coords[0][0]=2.5;
  outer_coords[0][1]=-600.0;
  
  outer_coords[1][0]=7.0;
  outer_coords[1][1]=-500.0;

  outer_coords[2][0]=9.01;
  outer_coords[2][1]=-165.0;
  
  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  outer_poro_outer_polyline_boundary_pt[
   Outer_poro_SSS_interface_near_outlet_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,
                            Outer_poro_SSS_interface_near_outlet_boundary_id);
 }



 // Block
 unsigned n_vertex_block=0;
 {
  Vector<Vector<double> > outer_coords(4, Vector<double>(2));

  outer_coords[0][0]=9.01;
  outer_coords[0][1]=-165.0;

  outer_coords[1][0]=5.82;
  outer_coords[1][1]=-160.0;
  
  outer_coords[2][0]=5.905;
  outer_coords[2][1]=-140.0;
  
  outer_coords[3][0]=9.19;
  outer_coords[3][1]=-135.0;

  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  outer_poro_outer_polyline_boundary_pt[
   Outer_poro_SSS_interface_block_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,
                            Outer_poro_SSS_interface_block_boundary_id);

  // Remember for connection
  n_vertex_block=n;
 }


 // Inner/FSI near inlet
 {
  Vector<Vector<double> > outer_coords(3, Vector<double>(2));
  
  outer_coords[0][0]=9.19;
  outer_coords[0][1]=-135.0;
  
  // mjr extra vertex necessary because of a possible bug which means you cannot
  // connect an internal boundary line to the second last (!?) vertex of
  // the outer boundary!!! // hierher check
  outer_coords[1][0]=9.595;
  outer_coords[1][1]=-67.5;
  
  outer_coords[2][0]=10.0;
  outer_coords[2][1]=0.0;
  
  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  outer_poro_outer_polyline_boundary_pt[
   Outer_poro_SSS_interface_near_inlet_boundary_id] =
   new TriangleMeshPolyLine(outer_coords,
                            Outer_poro_SSS_interface_near_inlet_boundary_id);
 }

 TriangleMeshClosedCurve * outer_poro_outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_poro_outer_polyline_boundary_pt);

 // Inner boundary -- separation of block from rest
 Vector<TriangleMeshOpenCurve*> outer_inner_open_boundary_pt(1);

 // Boundary of block
 {
  Vector<Vector<double> > outer_coords(2, Vector<double>(2,0.0));
  
  outer_coords[0][0]=9.01;
  outer_coords[0][1]=-165.0;
  
  outer_coords[1][0]=9.19;
  outer_coords[1][1]=-135.0;
  
  // Rescale
  unsigned n=outer_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    outer_coords[j][0]*=Global_Physical_Variables::R_scale;
    outer_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }

  // Build
  TriangleMeshPolyLine *outer_inner_polyline_pt =
   new TriangleMeshPolyLine(outer_coords,
                            Outer_poro_inner_block_boundary_id);
  
  // Connect initial vertex to second (starting from zeroth) vertex 
  // on fsi boundary near outlet
  outer_inner_polyline_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (outer_poro_outer_polyline_boundary_pt
    [Outer_poro_SSS_interface_near_outlet_boundary_id]),2);
  
  // Connect final vertex to final vertex of
  // the block
  outer_inner_polyline_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>
   (outer_poro_outer_polyline_boundary_pt[
    Outer_poro_SSS_interface_block_boundary_id]),n_vertex_block-1);

  // Store in vector
  Vector<TriangleMeshCurveSection*> outer_inner_curve_section_pt(1);
  outer_inner_curve_section_pt[0] = outer_inner_polyline_pt;
  
  // Create internal open curve
  outer_inner_open_boundary_pt[0] =
   new TriangleMeshOpenCurve(outer_inner_curve_section_pt);  
 }
 
 // Mesh parameters
 TriangleMeshParameters
  outer_poro_mesh_parameters(outer_poro_outer_boundary_pt);
 outer_poro_mesh_parameters.element_area() =
  Global_Physical_Variables::Element_area_poro;
 
 // Define inner boundaries
 outer_poro_mesh_parameters.internal_open_curves_pt() =
  outer_inner_open_boundary_pt;


 // Define region coordinates for outer region (dura)
 region_coords[0] = Global_Physical_Variables::R_scale*10.0;
 region_coords[1] = Global_Physical_Variables::Z_scale*(-1.0);
 outer_poro_mesh_parameters.add_region_coordinates(Dura_region_id
                                                   ,region_coords);

 // Define region coordinates for block
 region_coords[0] = Global_Physical_Variables::R_scale*7.3;
 region_coords[1] = Global_Physical_Variables::Z_scale*(-160.0);
 outer_poro_mesh_parameters.add_region_coordinates(Block_region_id,
                                                   region_coords);

 // Build mesh
 Outer_poro_mesh_pt = new TriangleMesh<PORO_ELEMENT>(
   outer_poro_mesh_parameters, Poro_time_stepper_pt);

 // Mesh as geom object representation of outer poro mesh
 Outer_poro_mesh_geom_obj_pt=new MeshAsGeomObject(Outer_poro_mesh_pt);
 
 // Extract regularly spaced points
 if (!Global_Physical_Variables::Suppress_regularly_spaced_output)
  {
   double t_start = TimingHelpers::timer();

   // Speed up search (somewhat hand-tuned; adjust if this takes too long
   // or misses too many points)
   Outer_poro_mesh_geom_obj_pt->max_spiral_level()=4;
   unsigned nx_back= Multi_domain_functions::Nx_bin;
   Multi_domain_functions::Nx_bin=10;
   unsigned ny_back= Multi_domain_functions::Ny_bin;
   Multi_domain_functions::Ny_bin=10;
   
   // Populate it...
   Vector<double> x(2);
   Vector<double> s(2);
   for (unsigned ir=0;ir<Global_Physical_Variables::Nplot_r_regular;ir++)
    {
     // outermost radius of dura
     x[0]=double(ir)/double(Global_Physical_Variables::Nplot_r_regular-1)*0.55;
     for (unsigned iz=0;iz<Global_Physical_Variables::Nplot_z_regular;iz++)
      {
       x[1]=double(iz)/double(Global_Physical_Variables::Nplot_z_regular-1)*
        (-600.0)*Global_Physical_Variables::Z_scale;
       
       // Pointer to GeomObject that contains this point
       GeomObject* geom_obj_pt=0;
       
       // Get it
       Outer_poro_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);
       
       // Store it
       if (geom_obj_pt!=0)
        {
         std::pair<PORO_ELEMENT*,Vector<double> > tmp;
         tmp = std::make_pair(dynamic_cast<PORO_ELEMENT*>(geom_obj_pt),s);
         Outer_poro_regularly_spaced_plot_point.push_back(tmp);
        }
      }
    }

   // Reset number of bins in binning method 
   Multi_domain_functions::Nx_bin=nx_back;
   Multi_domain_functions::Ny_bin=ny_back;
   oomph_info << "Took: " << TimingHelpers::timer()-t_start
              << " sec to setup regularly spaced points for outer poro mesh\n";
  }
 

 // Output mesh and boundaries
 Outer_poro_mesh_pt->output
  (Global_Physical_Variables::Directory+"/outer_poro_mesh.dat");
 Outer_poro_mesh_pt->output_boundaries
  (Global_Physical_Variables::Directory+"/outer_poro_mesh_boundaries.dat");
 
 // ---------------------------------------------------------------------------

 // Syrinx fluid mesh
 // -----------------
 Vector<TriangleMeshCurveSection*> syrinx_fluid_outer_polyline_boundary_pt(2);

 // FSI boundary
 {
  Vector<Vector<double> > syrinx_coords(4, Vector<double>(2));
  syrinx_coords[0][0]=0.0;
  syrinx_coords[0][1]=-80.0;
  
  syrinx_coords[1][0]=4.512;
  syrinx_coords[1][1]=-90.0;
  
  syrinx_coords[2][0]=4.128;
  syrinx_coords[2][1]=-210.0;
  
  syrinx_coords[3][0]=0.0;
  syrinx_coords[3][1]=-220.0;
  
  unsigned n=syrinx_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    syrinx_coords[j][0]*=Global_Physical_Variables::R_scale;
    syrinx_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  syrinx_fluid_outer_polyline_boundary_pt
   [Syrinx_fluid_inner_poro_interface_boundary_id] =
   new TriangleMeshPolyLine(syrinx_coords,
                            Syrinx_fluid_inner_poro_interface_boundary_id);
 }


 double z_bl_intersect_lo=0.0;
 double z_bl_intersect_hi=0.0;

 // Symm boundary
 {
  Vector<Vector<double> > syrinx_coords(4, Vector<double>(2));

  Vector<Vector<double> > line1(2,Vector<double>(2));
  Vector<Vector<double> > line2(2,Vector<double>(2));

  syrinx_coords[0][0]=0.0;
  syrinx_coords[0][1]=-220.0;

  // Symm line (offset; BL to right!)
  line1[0][0]=-Global_Physical_Variables::BL_thick;
  line1[0][1]=-100.0;
  line1[1][0]=-Global_Physical_Variables::BL_thick;
  line1[1][1]=0.0;

  // Downstream face of syrinx (BL to right)
  line2[0][0]=4.128;
  line2[0][1]=-210.0;
  line2[1][0]=0.0;
  line2[1][1]=-220.0;

  // Get intersection
  Vector<double> bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  syrinx_coords[1][0]=0.0;
  syrinx_coords[1][1]=bl_intersect[1];
  z_bl_intersect_lo=bl_intersect[1];

  // Symm line (offset; BL to right!)
  line1[0][0]=-Global_Physical_Variables::BL_thick;
  line1[0][1]=-100.0;
  line1[1][0]=-Global_Physical_Variables::BL_thick;
  line1[1][1]=0.0;

  // Upstream face of syrinx (BL to right)
  line2[0][0]=0.0;
  line2[0][1]=-80.0;
  line2[1][0]=4.512;
  line2[1][1]=-90.0;

  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  syrinx_coords[2][0]=0.0;
  syrinx_coords[2][1]=bl_intersect[1];
  z_bl_intersect_hi=bl_intersect[1];

  syrinx_coords[3][0]=0.0;
  syrinx_coords[3][1]=-80.0;
  
  unsigned n=syrinx_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    syrinx_coords[j][0]*=Global_Physical_Variables::R_scale;
    syrinx_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  syrinx_fluid_outer_polyline_boundary_pt
   [Syrinx_fluid_symmetry_boundary_id] =
   new TriangleMeshPolyLine(syrinx_coords,
                            Syrinx_fluid_symmetry_boundary_id);
 }

 TriangleMeshClosedCurve * syrinx_fluid_outer_boundary_pt =
  new TriangleMeshClosedCurve(syrinx_fluid_outer_polyline_boundary_pt);


 // One internal boundary for BL
 Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt(1);
 TriangleMeshPolyLine* inner_polyline_pt = 0;
 Vector<Vector<double> > inner_coords(4, Vector<double>(2,0.0));

  Vector<Vector<double> > line1(2,Vector<double>(2));
  Vector<Vector<double> > line2(2,Vector<double>(2));

 inner_coords[0][0]=0.0;
 inner_coords[0][1]=z_bl_intersect_lo;
 
 // Inner surface of syrinx (BL to right!)
 line1[0][0]=4.512;
 line1[0][1]=-90.0;
 line1[1][0]=4.128;
 line1[1][1]=-210.0;
 
 // Downstream face of syrinx (BL to right)
 line2[0][0]=4.128;
 line2[0][1]=-210.0;
 line2[1][0]=0.0;
 line2[1][1]=-220.0;
 
 // Get intersection
 Vector<double> bl_intersect=
  Global_Physical_Variables::bl_intersection(
   line1,line2,Global_Physical_Variables::BL_thick);
 
 inner_coords[1][0]=bl_intersect[0];
 inner_coords[1][1]=bl_intersect[1];
 

 // Inner surface of syrinx (BL to right!)
 line1[0][0]=4.512;
 line1[0][1]=-90.0;
 line1[1][0]=4.128;
 line1[1][1]=-210.0;
 
 // Upstream face of syrinx (BL to right)
 line2[0][0]=0.0;
 line2[0][1]=-80.0;
 line2[1][0]=4.512;
 line2[1][1]=-90.0;
 
 // Get intersection
 bl_intersect=
  Global_Physical_Variables::bl_intersection(
   line1,line2,Global_Physical_Variables::BL_thick);
 
 inner_coords[2][0]=bl_intersect[0];
 inner_coords[2][1]=bl_intersect[1];
 
 inner_coords[3][0]=0.0;
 inner_coords[3][1]=z_bl_intersect_hi;
 
 // Rescale
 for(unsigned i=0;i<4;i++)
  {
   inner_coords[i][0]*=Global_Physical_Variables::R_scale;
   inner_coords[i][1]*=Global_Physical_Variables::Z_scale;
  }
 
 // Build
 inner_polyline_pt =
  new TriangleMeshPolyLine(inner_coords,
                           Syrinx_bl_boundary_id);
 
 // Connect initial vertex to first (from zeroth) vertex on sym boundary
 inner_polyline_pt->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (syrinx_fluid_outer_polyline_boundary_pt[Syrinx_fluid_symmetry_boundary_id]),
  1);
 
 // Connect final vertex to second (from zeroth) vertex on outlet boundary
 inner_polyline_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (syrinx_fluid_outer_polyline_boundary_pt[Syrinx_fluid_symmetry_boundary_id]),
  2);
 
 // Store in vector
 Vector<TriangleMeshCurveSection*> inner_curve_section_pt(1);
 inner_curve_section_pt[0] = inner_polyline_pt;
 
 // Create internal open curve
 inner_open_boundary_pt[0] =
  new TriangleMeshOpenCurve(inner_curve_section_pt);


 // Mesh parameters
 TriangleMeshParameters
  syrinx_fluid_mesh_parameters(syrinx_fluid_outer_boundary_pt);

 syrinx_fluid_mesh_parameters.element_area() =
  Global_Physical_Variables::Element_area_fluid;

 // Define inner boundaries
 syrinx_fluid_mesh_parameters.internal_open_curves_pt()=inner_open_boundary_pt;
 
 // Build mesh
 Syrinx_fluid_mesh_pt = new TriangleMesh<FLUID_ELEMENT>(
   syrinx_fluid_mesh_parameters, Fluid_time_stepper_pt);

 // Mesh as geom object representation of syrinx fluid mesh
 Syrinx_fluid_mesh_geom_obj_pt=new MeshAsGeomObject(Syrinx_fluid_mesh_pt);
 
 // Extract regularly spaced points
 if (!Global_Physical_Variables::Suppress_regularly_spaced_output)
  {
   double t_start = TimingHelpers::timer();

   // Speed up search (somewhat hand-tuned; adjust if this takes too long
   // or misses too many points)
   Syrinx_fluid_mesh_geom_obj_pt->max_spiral_level()=4;
   unsigned nx_back= Multi_domain_functions::Nx_bin;
   Multi_domain_functions::Nx_bin=10;
   unsigned ny_back= Multi_domain_functions::Ny_bin;
   Multi_domain_functions::Ny_bin=10;
   
   // Populate it...
   Vector<double> x(2);
   Vector<double> s(2);
   for (unsigned ir=0;ir<Global_Physical_Variables::Nplot_r_regular;ir++)
    {
     // outermost radius of dura
     x[0]=double(ir)/double(Global_Physical_Variables::Nplot_r_regular-1)*0.55;
     for (unsigned iz=0;iz<Global_Physical_Variables::Nplot_z_regular;iz++)
      {
       x[1]=double(iz)/double(Global_Physical_Variables::Nplot_z_regular-1)*
        (-600.0)*Global_Physical_Variables::Z_scale;
       
       // Pointer to GeomObject that contains this point
       GeomObject* geom_obj_pt=0;
       
       // Get it
       Syrinx_fluid_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);
       
       // Store it
       if (geom_obj_pt!=0)
        {
         std::pair<FLUID_ELEMENT*,Vector<double> > tmp;
         tmp = std::make_pair(dynamic_cast<FLUID_ELEMENT*>(geom_obj_pt),s);
         Syrinx_fluid_regularly_spaced_plot_point.push_back(tmp);
        }
      }
    }

   // Reset number of bins in binning method 
   Multi_domain_functions::Nx_bin=nx_back;
   Multi_domain_functions::Ny_bin=ny_back;
   oomph_info << "Took: " << TimingHelpers::timer()-t_start
              << " sec to setup regularly spaced points for syrinx mesh\n";
  }

 // Output mesh and boundaries
 Syrinx_fluid_mesh_pt->output
  (Global_Physical_Variables::Directory+"/syrinx_fluid_mesh.dat");
 Syrinx_fluid_mesh_pt->output_boundaries
  (Global_Physical_Variables::Directory+"/syrinx_fluid_mesh_boundaries.dat");
 
 // ----------------------------------------------------------------------------

 // SSS fluid mesh
 // --------------
 Vector<TriangleMeshCurveSection*> sss_fluid_outer_polyline_boundary_pt(6);

 // Block
 {
  Vector<Vector<double> > sss_coords;
  Vector<double> vert(2);
  vert[0]=9.19;
  vert[1]=-135.0;
  sss_coords.push_back(vert);

  Vector<double> start(2);
  start[0]=5.905;
  start[1]=-140.0;

  Vector<double> end(2);
  end[0]=5.82;
  end[1]=-160.0;
  
  Global_Physical_Variables::push_back_vertices(
   start,end,Global_Physical_Variables::BL_thick*
   Global_Physical_Variables::Z_shrink_factor,sss_coords);

  vert[0]=9.01;
  vert[1]=-165.0;
  sss_coords.push_back(vert);

  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_outer_poro_interface_block_boundary_id]=
   new TriangleMeshPolyLine(sss_coords,
                            SSS_fluid_outer_poro_interface_block_boundary_id);
 }

 // Outer FSI interface near outlet
 {
  Vector<Vector<double> > sss_coords(3, Vector<double>(2));
  
  sss_coords[0][0]=9.01;
  sss_coords[0][1]=-165.0;

  sss_coords[1][0]=7.0;
  sss_coords[1][1]=-500.0;
  
  sss_coords[2][0]=2.5;
  sss_coords[2][1]=-600.0;

  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_outer_poro_interface_near_outlet_boundary_id]=
   new TriangleMeshPolyLine(
    sss_coords,
    SSS_fluid_outer_poro_interface_near_outlet_boundary_id);
 }
 

 // Outlet
 {
  Vector<Vector<double> > sss_coords(4, Vector<double>(2));

  sss_coords[0][0]=2.5;
  sss_coords[0][1]=-600.0;

  sss_coords[1][0]=2.5-Global_Physical_Variables::BL_thick;
  sss_coords[1][1]=-600.0;

  sss_coords[2][0]=1.25+Global_Physical_Variables::BL_thick;
  sss_coords[2][1]=-600.0;

  sss_coords[3][0]=1.25;
  sss_coords[3][1]=-600.0;

  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_outflow_boundary_id]=
   new TriangleMeshPolyLine(sss_coords,SSS_fluid_outflow_boundary_id);
 }
 
 // Inner FSI
 {


  Vector<Vector<double> > sss_coords;
  Vector<double> vert(2);
  vert[0]=1.25;
  vert[1]=-600.0;
  sss_coords.push_back(vert);

  vert[0]=1.25;
  vert[1]=-480.0;
  sss_coords.push_back(vert);

  Vector<double> start(2);
  start[0]=4.24;
  start[1]=-440.0;
  
  Vector<double> end(2);
  end[0]=6.0;
  end[1]=0.0;

  Global_Physical_Variables::push_back_vertices(
   start,end,Global_Physical_Variables::BL_thick*
   Global_Physical_Variables::Z_shrink_factor,sss_coords);

  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_inner_poro_interface_boundary_id]=
   new TriangleMeshPolyLine(sss_coords,
                            SSS_fluid_inner_poro_interface_boundary_id);
 }
 

 // Inlet
 {
  Vector<Vector<double> > sss_coords(4, Vector<double>(2));

  sss_coords[0][0]=6.0;
  sss_coords[0][1]=0.0;

  sss_coords[1][0]=6.0+Global_Physical_Variables::BL_thick;
  sss_coords[1][1]=0.0;

  sss_coords[2][0]=10.0-Global_Physical_Variables::BL_thick;
  sss_coords[2][1]=0.0;
  
  sss_coords[3][0]=10.0;
  sss_coords[3][1]=0.0;

  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_inflow_boundary_id]=
   new TriangleMeshPolyLine(sss_coords,
                            SSS_fluid_inflow_boundary_id);
 }



 // Outer FSI interface near inlet
 {
  Vector<Vector<double> > sss_coords(3, Vector<double>(2));

  sss_coords[0][0]=10.0;
  sss_coords[0][1]=0.0;

  sss_coords[1][0]=9.595;
  sss_coords[1][1]=-67.5;
  
  sss_coords[2][0]=9.19;
  sss_coords[2][1]=-135.0;
  
  // Rescale
  unsigned n=sss_coords.size();
  for (unsigned j=0;j<n;j++)
   {
    sss_coords[j][0]*=Global_Physical_Variables::R_scale;
    sss_coords[j][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  sss_fluid_outer_polyline_boundary_pt[
   SSS_fluid_outer_poro_interface_near_inlet_boundary_id]=
   new TriangleMeshPolyLine(
    sss_coords,
    SSS_fluid_outer_poro_interface_near_inlet_boundary_id);
 }
 
 TriangleMeshClosedCurve * sss_fluid_outer_boundary_pt =
  new TriangleMeshClosedCurve(sss_fluid_outer_polyline_boundary_pt);


 // Two internal boundaries for BL
 inner_open_boundary_pt.resize(2);
 TriangleMeshPolyLine* inner_polyline1_pt = 0;
 TriangleMeshPolyLine* inner_polyline2_pt = 0;

 // First one
 {
  // Boundary lines
  Vector<Vector<double> > line1(2,Vector<double>(2));
  Vector<Vector<double> > line2(2,Vector<double>(2));
  
  Vector<Vector<double> > inner_coords(4, Vector<double>(2,0.0));
  
  inner_coords[0][0]=6.0+Global_Physical_Variables::BL_thick;
  inner_coords[0][1]=0.0;
  
  // Inner wall (BL to right!)
  line1[0][0]=4.24;
  line1[0][1]=-440.0;
  line1[1][0]=6.0;
  line1[1][1]=0.0;

  // Upstream edge of block (BL to right!)
  line2[0][0]=9.19;
  line2[0][1]=-135.0;
  line2[1][0]=5.905;
  line2[1][1]=-140.0;

  // Get intersection
  Vector<double> bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[1][0]=bl_intersect[0];
  inner_coords[1][1]=bl_intersect[1];


  // Upstream edge of block (BL to right!)
  line1[0][0]=9.19;
  line1[0][1]=-135.0;
  line1[1][0]=5.905;
  line1[1][1]=-140.0;

  // Outer wall (BL to right)
  line2[0][0]=10.0;
  line2[0][1]=0.0;
  line2[1][0]=9.19;
  line2[1][1]=-135.0;

  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it
  inner_coords[2][0]=bl_intersect[0];
  inner_coords[2][1]=bl_intersect[1];
  
  // Final one
  inner_coords[3][0]=10.0-Global_Physical_Variables::BL_thick;
  inner_coords[3][1]=0.0;
  
  // Rescale
  for(unsigned i=0;i<4;i++)
   {
    inner_coords[i][0]*=Global_Physical_Variables::R_scale;
    inner_coords[i][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_polyline1_pt =
   new TriangleMeshPolyLine(inner_coords,
                            SSS_bl_near_inlet_boundary_id);
 }

 // Connect initial vertex to first (from zeroth) vertex on inlet boundary
 inner_polyline1_pt->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (sss_fluid_outer_polyline_boundary_pt[SSS_fluid_inflow_boundary_id]),1);
 
 // Connect final vertex to second (from zeroth) vertex on inlet boundary
 inner_polyline1_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (sss_fluid_outer_polyline_boundary_pt[SSS_fluid_inflow_boundary_id]),2);
 

 // Second one
 {
  // Boundary lines
  Vector<Vector<double> > line1(2,Vector<double>(2));
  Vector<Vector<double> > line2(2,Vector<double>(2));
  
  Vector<Vector<double> > inner_coords(7, Vector<double>(2,0.0));
  
  inner_coords[0][0]=2.5-Global_Physical_Variables::BL_thick;
  inner_coords[0][1]=-600.0;
  
  // Lowest part of inner wall (BL to right!)
  line1[0][0]=7.0;
  line1[0][1]=-500.0;
  line1[1][0]=2.5-Global_Physical_Variables::BL_thick;
  line1[1][1]=-600.0;
  
  // Upper part of inner wall (BL to right)
  line2[0][0]=9.01;
  line2[0][1]=-165.0;
  line2[1][0]=7.0;
  line2[1][1]=-500.0;
  
  // Get intersection
  Vector<double> bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[1][0]=bl_intersect[0];
  inner_coords[1][1]=bl_intersect[1];

  // Downstream edge of block (BL to right!)
  line1[0][0]=5.82;
  line1[0][1]=-160.0;
  line1[1][0]=9.01;
  line1[1][1]=-165.0;

  // Upper part of inner wall (BL to right)
  line2[0][0]=9.01;
  line2[0][1]=-165.0;
  line2[1][0]=7.0;
  line2[1][1]=-500.0;

  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[2][0]=bl_intersect[0];
  inner_coords[2][1]=bl_intersect[1];
  
  // Downstream edge of block (BL to right!)
  line1[0][0]=5.82;
  line1[0][1]=-160.0;
  line1[1][0]=9.01;
  line1[1][1]=-165.0;

  // Upper part of inner wall (BL to right!)
  line2[0][0]=4.24;
  line2[0][1]=-440.0;
  line2[1][0]=6.0;
  line2[1][1]=0.0;
  
  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[3][0]=bl_intersect[0];
  inner_coords[3][1]=bl_intersect[1];
  



  // Middle part of inner wall (BL to right!)
  line1[0][0]=1.25;
  line1[0][1]=-480.0;
  line1[1][0]=4.24;
  line1[1][1]=-440.0;

  // Upper part of inner wall (BL to right!)
  line2[0][0]=4.24;
  line2[0][1]=-440.0;
  line2[1][0]=6.0;
  line2[1][1]=0.0;
  
  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[4][0]=bl_intersect[0];
  inner_coords[4][1]=bl_intersect[1];
  

  // Middle part of inner wall (BL to right!)
  line1[0][0]=1.25;
  line1[0][1]=-480.0;
  line1[1][0]=4.24;
  line1[1][1]=-440.0;

  // Lower part of inner wall (BL to right!)
  line2[0][0]=1.25;
  line2[0][1]=-600.0;
  line2[1][0]=1.25;
  line2[1][1]=-480.0;
  
  // Get intersection
  bl_intersect=
   Global_Physical_Variables::bl_intersection(
    line1,line2,Global_Physical_Variables::BL_thick);

  // Set it...
  inner_coords[5][0]=bl_intersect[0];
  inner_coords[5][1]=bl_intersect[1];
  

  // Final one
  inner_coords[6][0]=1.25+Global_Physical_Variables::BL_thick;
  inner_coords[6][1]=-600.0;
  

  // Rescale
  for(unsigned i=0;i<7;i++)
   {
    inner_coords[i][0]*=Global_Physical_Variables::R_scale;
    inner_coords[i][1]*=Global_Physical_Variables::Z_scale;
   }
  
  // Build
  inner_polyline2_pt =
   new TriangleMeshPolyLine(inner_coords,
                            SSS_bl_near_outlet_boundary_id);
 }

 // Connect initial vertex to first (from zeroth) vertex on outlet boundary
 inner_polyline2_pt->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (sss_fluid_outer_polyline_boundary_pt[SSS_fluid_outflow_boundary_id]),1);
 
 // Connect final vertex to second (from zeroth) vertex on outlet boundary
 inner_polyline2_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>
  (sss_fluid_outer_polyline_boundary_pt[SSS_fluid_outflow_boundary_id]),2);
 
 // Store in vector
 inner_curve_section_pt.resize(1);
 inner_curve_section_pt[0] = inner_polyline1_pt;
 
 // Create internal open curve
 inner_open_boundary_pt[0] = new TriangleMeshOpenCurve(inner_curve_section_pt);

 // Store in vector
 inner_curve_section_pt.resize(1);
 inner_curve_section_pt[0] = inner_polyline2_pt;
 
 // Create internal open curve
 inner_open_boundary_pt[1] = new TriangleMeshOpenCurve(inner_curve_section_pt);


 
 // Mesh parameters
 TriangleMeshParameters
  sss_fluid_mesh_parameters(sss_fluid_outer_boundary_pt);

 sss_fluid_mesh_parameters.element_area() =
  Global_Physical_Variables::Element_area_fluid;

 // Define inner boundaries
 sss_fluid_mesh_parameters.internal_open_curves_pt()=inner_open_boundary_pt;
 
 // Build mesh
 SSS_fluid_mesh_pt = new TriangleMesh<FLUID_ELEMENT>(
   sss_fluid_mesh_parameters, Fluid_time_stepper_pt);

 // Mesh as geom object representation of syrinx fluid mesh
 SSS_fluid_mesh_geom_obj_pt=new MeshAsGeomObject(SSS_fluid_mesh_pt);
 
 // Extract regularly spaced points
 if (!Global_Physical_Variables::Suppress_regularly_spaced_output)
  {
   double t_start = TimingHelpers::timer();

   // Speed up search (somewhat hand-tuned; adjust if this takes too long
   // or misses too many points)
   SSS_fluid_mesh_geom_obj_pt->max_spiral_level()=4;
   unsigned nx_back= Multi_domain_functions::Nx_bin;
   Multi_domain_functions::Nx_bin=10;
   unsigned ny_back= Multi_domain_functions::Ny_bin;
   Multi_domain_functions::Ny_bin=10;

   // Populate it...
   Vector<double> x(2);
   Vector<double> s(2);
   for (unsigned ir=0;ir<Global_Physical_Variables::Nplot_r_regular;ir++)
    {
     // outermost radius of dura
     x[0]=double(ir)/double(Global_Physical_Variables::Nplot_r_regular-1)*0.55; 
     for (unsigned iz=0;iz<Global_Physical_Variables::Nplot_z_regular;iz++)
      {
       x[1]=double(iz)/double(Global_Physical_Variables::Nplot_z_regular-1)*
        (-600.0)*Global_Physical_Variables::Z_scale;
       
       // Pointer to GeomObject that contains this point
       GeomObject* geom_obj_pt=0;
       
       // Get it
       SSS_fluid_mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s);

       // Store it
       if (geom_obj_pt!=0)
        {
         std::pair<FLUID_ELEMENT*,Vector<double> > tmp;
         tmp = std::make_pair(dynamic_cast<FLUID_ELEMENT*>(geom_obj_pt),s);
         SSS_fluid_regularly_spaced_plot_point.push_back(tmp);
        }
      }
    }

   // Reset number of bins in binning method 
   Multi_domain_functions::Nx_bin=nx_back;
   Multi_domain_functions::Ny_bin=ny_back;
   oomph_info << "Took: " << TimingHelpers::timer()-t_start
              << " sec to setup regularly spaced points for sss mesh\n";
  }
 
  // Output mesh and boundaries
 SSS_fluid_mesh_pt->output
  (Global_Physical_Variables::Directory+"/sss_fluid_mesh.dat");
 SSS_fluid_mesh_pt->output_boundaries
  (Global_Physical_Variables::Directory+"/sss_fluid_mesh_boundaries.dat");
 


 // ---------------------------------------------------------------------------

 // Undo z-scaling?
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--undo_z_scaling_by_element_stretching"))
  {
   Vector<Mesh*> all_bulk_mesh_pt;
   all_bulk_mesh_pt.push_back(SSS_fluid_mesh_pt);
   all_bulk_mesh_pt.push_back(Syrinx_fluid_mesh_pt);
   all_bulk_mesh_pt.push_back(Outer_poro_mesh_pt);
   all_bulk_mesh_pt.push_back(Inner_poro_mesh_pt);
   unsigned n=all_bulk_mesh_pt.size();
   for (unsigned m=0;m<n;m++)
    {
     Mesh* my_mesh_pt=all_bulk_mesh_pt[m];
     unsigned nnod=my_mesh_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       my_mesh_pt->node_pt(j)->x(1)*=Global_Physical_Variables::Z_shrink_factor;
      }
    }
  }


 // ---------------------------------------------------------------------------

 // Create the FSI meshes between inner poro and syrinx fluid
 Inner_poro_syrinx_fluid_FSI_surface_mesh_pt=new Mesh;
 Syrinx_fluid_inner_poro_FSI_surface_mesh_pt=new Mesh;
 
 // Create the FSI meshes between inner poro and SSS fluid
 Inner_poro_SSS_fluid_FSI_surface_mesh_pt=new Mesh;
 SSS_fluid_inner_poro_FSI_surface_mesh_pt=new Mesh;

 // Create the FSI meshes between outer poro and SSS fluid
 Outer_poro_SSS_fluid_near_inlet_FSI_surface_mesh_pt=new Mesh;
 SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt=new Mesh;
 Outer_poro_SSS_fluid_near_outlet_FSI_surface_mesh_pt=new Mesh;
 SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt=new Mesh;

 // Create the FSI meshes beteween outer poro and SSS fluid for block
 Outer_poro_SSS_fluid_block_FSI_surface_mesh_pt=new Mesh;
 SSS_fluid_outer_poro_block_FSI_surface_mesh_pt=new Mesh;

 // Create the surface mesh for the inflow boundary
 Inflow_fluid_surface_mesh_pt = new Mesh;
 
 // ---------------------------------------------------------------------------

 // Make surface meshes
 create_fluid_traction_elements();
 create_poro_face_elements();

 // ----------------------------------------------------------------------------

 // Complete the problem setup to make the elements fully functional
 
 // Do edge sign setup for inner poro mesh
 Global_Physical_Variables::edge_sign_setup<PORO_ELEMENT>(Inner_poro_mesh_pt);
 
 // Loop over the elements in cord region of inner poro mesh
 //---------------------------------------------------------
 {

  ofstream some_file;
  char filename[100];
  sprintf(filename,"cord_region.dat");
  some_file.open(filename);

  unsigned r=Cord_region_id;
  unsigned nel=dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)->
   nregion_element(r);
  for(unsigned e=0;e<nel;e++)
   {
    // Cast to a bulk element
    PORO_ELEMENT *el_pt=dynamic_cast<PORO_ELEMENT*>(
     dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)
     ->region_element_pt(r,e));
    
    // Output region just using vertices
    el_pt->output(some_file,2);

    // Set the pointer to the Lambda^2 parameter (universal)
    el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    
    // Set the pointer to non-dim permeability (universal)
    el_pt->permeability_pt()= &Global_Physical_Variables::Permeability;
    
    // Set the pointer to the Biot parameter (same for all; incompressible
    // fluid and solid)
    el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;
    
    // Set the pointer to the density ratio (pore fluid to solid) (same for all)
    el_pt->density_ratio_pt()= &Global_Physical_Variables::Density_ratio;
    
    // Set the pointer to Poisson's ratio
    el_pt->nu_pt() = &Global_Physical_Variables::Nu_cord;
    
    // Set the pointer to non-dim Young's modulus (ratio of Young's modulus
    // to Young's modulus used in the non-dimensionalisation of the
    // equations
    el_pt->youngs_modulus_pt() = 
     &Global_Physical_Variables::Stiffness_ratio_cord;
    
    // Set the pointer to permeability ratio (ratio of actual permeability
    // to the one used in the non-dim of the equations)
    el_pt->permeability_ratio_pt()=
     &Global_Physical_Variables::Permeability_ratio_cord;
    
    // Set the pointer to the porosity
    el_pt->porosity_pt()=&Global_Physical_Variables::Porosity_cord;
    
    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Poro_time_stepper_pt);
    
    // Switch off Darcy flow?
    if(CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
     {
      el_pt->switch_off_darcy();
     }
   }// end loop over cord elements

  some_file.close();
 }

 // Loop over the elements in pia region of inner poro mesh
 //---------------------------------------------------------
 {

  ofstream some_file;
  char filename[100];
  sprintf(filename,"pia_region.dat");
  some_file.open(filename);

  unsigned r=Pia_region_id;
  unsigned nel=dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)->
   nregion_element(r);
  for(unsigned e=0;e<nel;e++)
   {
    // Cast to a bulk element
    PORO_ELEMENT *el_pt=dynamic_cast<PORO_ELEMENT*>(
     dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)
     ->region_element_pt(r,e));
        
    // Output region just using vertices
    el_pt->output(some_file,2);

    // Set the pointer to the Lambda^2 parameter (universal)
    el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    
    // Set the pointer to non-dim permeability (universal)
    el_pt->permeability_pt()= &Global_Physical_Variables::Permeability;
    
    // Set the pointer to the Biot parameter (same for all; incompressible
    // fluid and solid)
    el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;
    
    // Set the pointer to the density ratio (pore fluid to solid) (same for all)
    el_pt->density_ratio_pt()= &Global_Physical_Variables::Density_ratio;
    
    // Set the pointer to Poisson's ratio
    el_pt->nu_pt() = &Global_Physical_Variables::Nu_dura_and_pia;
    
    // Set the pointer to non-dim Young's modulus (ratio of Young's modulus
    // to Young's modulus used in the non-dimensionalisation of the
    // equations
    el_pt->youngs_modulus_pt() = 
     &Global_Physical_Variables::Stiffness_ratio_dura_and_pia;
    
    // Set the pointer to permeability ratio (ratio of actual permeability
    // to the one used in the non-dim of the equations)
    el_pt->permeability_ratio_pt()=
     &Global_Physical_Variables::Permeability_ratio_dura_and_pia;
    
    // Set the pointer to the porosity
    el_pt->porosity_pt()=&Global_Physical_Variables::Porosity_dura_and_pia;
    
    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Poro_time_stepper_pt);
    
    // Switch off Darcy flow?
    if(CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
     {
      el_pt->switch_off_darcy();
     }
  }// end loop over pia elements

  some_file.close();
 }

 // Loop over the elements in filum region of inner poro mesh
 //----------------------------------------------------------
 {
  ofstream some_file;
  char filename[100];
  sprintf(filename,"filum_region.dat");
  some_file.open(filename);

  unsigned r=Filum_region_id;
  unsigned nel=dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)->
   nregion_element(r);
  for(unsigned e=0;e<nel;e++)
   {
    // Cast to a bulk element
    PORO_ELEMENT *el_pt=dynamic_cast<PORO_ELEMENT*>(
     dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Inner_poro_mesh_pt)
     ->region_element_pt(r,e));
    
    // Output region just using vertices
    el_pt->output(some_file,2);

    // Set the pointer to the Lambda^2 parameter (universal)
    el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    
    // Set the pointer to non-dim permeability (universal)
    el_pt->permeability_pt()= &Global_Physical_Variables::Permeability;
    
    // Set the pointer to the Biot parameter (same for all; incompressible
    // fluid and solid)
    el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;
    
    // Set the pointer to the density ratio (pore fluid to solid) (same for all)
    el_pt->density_ratio_pt()= &Global_Physical_Variables::Density_ratio;
    
    // Set the pointer to Poisson's ratio
    el_pt->nu_pt() = &Global_Physical_Variables::Nu_filum_and_block;
    
    // Set the pointer to non-dim Young's modulus (ratio of Young's modulus
    // to Young's modulus used in the non-dimensionalisation of the
    // equations
    el_pt->youngs_modulus_pt() = 
     &Global_Physical_Variables::Stiffness_ratio_filum_and_block;
    
    // Set the pointer to permeability ratio (ratio of actual permeability
    // to the one used in the non-dim of the equations)
    el_pt->permeability_ratio_pt()=
     &Global_Physical_Variables::Permeability_ratio_filum_and_block;
    
    // Set the pointer to the porosity
    el_pt->porosity_pt()=&Global_Physical_Variables::Porosity_filum_and_block;
    
    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Poro_time_stepper_pt);
    
    // Switch off Darcy flow?
    if(CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
     {
      el_pt->switch_off_darcy();
     }
  }// end loop over cord elements

  some_file.close();
 }


 // Do edge sign setup for outer poro mesh
 Global_Physical_Variables::edge_sign_setup<PORO_ELEMENT>(Outer_poro_mesh_pt);


 // Loop over the elements in block region of outer poro mesh
 //----------------------------------------------------------
 {

  ofstream some_file;
  char filename[100];
  sprintf(filename,"block_region.dat");
  some_file.open(filename);

  unsigned r=Block_region_id;
  unsigned nel=dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Outer_poro_mesh_pt)->
   nregion_element(r);
  for(unsigned e=0;e<nel;e++)
   {
    // Cast to a bulk element
    PORO_ELEMENT *el_pt=dynamic_cast<PORO_ELEMENT*>(
     dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Outer_poro_mesh_pt)
     ->region_element_pt(r,e));
        
    // Output region just using vertices
    el_pt->output(some_file,2);

    // Set the pointer to the Lambda^2 parameter (universal)
    el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    
    // Set the pointer to non-dim permeability (universal)
    el_pt->permeability_pt()= &Global_Physical_Variables::Permeability;
    
    // Set the pointer to the Biot parameter (same for all; incompressible
    // fluid and solid)
    el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;
    
    // Set the pointer to the density ratio (pore fluid to solid) (same for all)
    el_pt->density_ratio_pt()= &Global_Physical_Variables::Density_ratio;
    
    // Set the pointer to Poisson's ratio
    el_pt->nu_pt() = &Global_Physical_Variables::Nu_filum_and_block;
    
    // Set the pointer to non-dim Young's modulus (ratio of Young's modulus
    // to Young's modulus used in the non-dimensionalisation of the
    // equations
    el_pt->youngs_modulus_pt() = 
     &Global_Physical_Variables::Stiffness_ratio_filum_and_block;
    
    // Set the pointer to permeability ratio (ratio of actual permeability
    // to the one used in the non-dim of the equations)
    el_pt->permeability_ratio_pt()=
     &Global_Physical_Variables::Permeability_ratio_filum_and_block;
    
    // Set the pointer to the porosity
    el_pt->porosity_pt()=&Global_Physical_Variables::Porosity_filum_and_block;
    
    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Poro_time_stepper_pt);
    
    // Switch off Darcy flow?
    if(CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
     {
      el_pt->switch_off_darcy();
     }
  }// end loop over block elements

  some_file.close();
 }

 // Loop over the elements in dura region of outer poro mesh
 //----------------------------------------------------------
 {
  ofstream some_file;
  char filename[100];
  sprintf(filename,"dura_region.dat");
  some_file.open(filename);

  unsigned r=Dura_region_id;
  unsigned nel=dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Outer_poro_mesh_pt)->
   nregion_element(r);
  for(unsigned e=0;e<nel;e++)
   {
    // Cast to a bulk element
    PORO_ELEMENT *el_pt=dynamic_cast<PORO_ELEMENT*>(
     dynamic_cast<TriangleMesh<PORO_ELEMENT>*>(Outer_poro_mesh_pt)
     ->region_element_pt(r,e));

    // Output region just using vertices
    el_pt->output(some_file,2);

    // Set the pointer to the Lambda^2 parameter (universal)
    el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    
    // Set the pointer to non-dim permeability (universal)
    el_pt->permeability_pt()= &Global_Physical_Variables::Permeability;
    
    // Set the pointer to the Biot parameter (same for all; incompressible
    // fluid and solid)
    el_pt->alpha_pt()=&Global_Physical_Variables::Alpha;
    
    // Set the pointer to the density ratio (pore fluid to solid) (same for all)
    el_pt->density_ratio_pt()= &Global_Physical_Variables::Density_ratio;
    
    // Set the pointer to Poisson's ratio
    el_pt->nu_pt() = &Global_Physical_Variables::Nu_dura_and_pia;
    
    // Set the pointer to non-dim Young's modulus (ratio of Young's modulus
    // to Young's modulus used in the non-dimensionalisation of the
    // equations
    el_pt->youngs_modulus_pt() = 
     &Global_Physical_Variables::Stiffness_ratio_dura_and_pia;
    
    // Set the pointer to permeability ratio (ratio of actual permeability
    // to the one used in the non-dim of the equations)
    el_pt->permeability_ratio_pt()=
     &Global_Physical_Variables::Permeability_ratio_dura_and_pia;
    
    // Set the pointer to the porosity
    el_pt->porosity_pt()=&Global_Physical_Variables::Porosity_dura_and_pia;
    
    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Poro_time_stepper_pt);
    
    // Switch off Darcy flow?
    if(CommandLineArgs::command_line_flag_has_been_set("--pin_darcy"))
     {
      el_pt->switch_off_darcy();
     }
  }// end loop over dura elements

  some_file.close();
 }


  
 // Loop over the elements in syrinx fluid mesh to set physical parameters etc
  //--------------------------------------------------------------------------
 unsigned n_el_syrinx_fluid = Syrinx_fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_el_syrinx_fluid;e++)
  {
   // Cast to the particular element type
   FLUID_ELEMENT *el_pt = dynamic_cast<FLUID_ELEMENT*>
    (Syrinx_fluid_mesh_pt->element_pt(e));

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

 // Loop over the elements in SSS fluid mesh to set physical parameters etc
 //------------------------------------------------------------------------
 unsigned n_el_sss_fluid = SSS_fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_el_sss_fluid;e++)
  {
   // Cast to the particular element type
   FLUID_ELEMENT *el_pt = dynamic_cast<FLUID_ELEMENT*>
    (SSS_fluid_mesh_pt->element_pt(e));

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


 // ------------------------------------------------------------------------

 // Boundary conditions for the inner poroelasticity mesh
 // ----------------------------------------------------

 // Axisymmetric conditions on these boundaries
 {

 // Vector of the boundaries we will pin on the inner poro mesh
  Vector<unsigned> inner_poro_pinned_boundaries;
  inner_poro_pinned_boundaries.push_back(
   Inner_poro_symmetry_near_inlet_boundary_id);
  inner_poro_pinned_boundaries.push_back(
   Inner_poro_symmetry_near_outlet_boundary_id);

  // Get the number of boundaries from the vector
  unsigned n_b=inner_poro_pinned_boundaries.size();
  
  // Set boundary conditions for the inner poro mesh
  for(unsigned ibound=0;ibound<n_b;ibound++)
   {
    unsigned num_nod=Inner_poro_mesh_pt->nboundary_node(
     inner_poro_pinned_boundaries[ibound]);
    for(unsigned inod=0;inod<num_nod;inod++)
     {
      // Get pointer to node
      Node* nod_pt=Inner_poro_mesh_pt->boundary_node_pt(
       inner_poro_pinned_boundaries[ibound],inod);
      
      // Axisymmetric conditions and u_z = 0 on r = 0
      
      // Pin u_r
      nod_pt->pin(0);
      
      // If this node stores edge flux dofs, pin these
      if(nod_pt->nvalue()==4)
       {
        nod_pt->pin(2); 
       }
      
     } // end_of_loop_over_nodes
   } // end_of_loop_over_inner_poro_boundaries
 }
 
 // No porous flux or solid displacement on these
 {
  Vector<unsigned> inner_poro_pinned_boundaries;
  inner_poro_pinned_boundaries.push_back(Inner_poro_inlet_boundary_id);
  inner_poro_pinned_boundaries.push_back(Inner_poro_outlet_boundary_id);
  
  // Get the number of boundaries from the vector
  unsigned n_b=inner_poro_pinned_boundaries.size();
  
  // Set boundary conditions for the inner poro mesh
  for(unsigned ibound=0;ibound<n_b;ibound++)
   {
    unsigned num_nod=Inner_poro_mesh_pt->nboundary_node(
     inner_poro_pinned_boundaries[ibound]);
    for(unsigned inod=0;inod<num_nod;inod++)
     {
      // Get pointer to node
      Node* nod_pt=Inner_poro_mesh_pt->boundary_node_pt(
       inner_poro_pinned_boundaries[ibound],inod);
      
      // Axisymmetric conditions and u_z = 0 on r = 0
      // and no porous flux or solid displacement on the inlet and outlet
      
      // Pin u_r, u_z
      nod_pt->pin(0);
      nod_pt->pin(1);
           
      // If this node stores edge flux dofs, pin these
      if(nod_pt->nvalue()==4)
       {
        nod_pt->pin(2);
        nod_pt->pin(3);
       }
      
     } // end_of_loop_over_nodes
   } // end_of_loop_over_inner_poro_boundaries
 }

 // Boundary conditions for the outer poroelasticity mesh
 // --------------------------------------------------------

 // Vector of the boundaries we will pin on the outer poro mesh
 std::vector<unsigned> outer_poro_pinned_boundaries;

 // No porous flux or solid displacement on these boundaries
 outer_poro_pinned_boundaries.push_back(Outer_poro_inlet_boundary_id);
 outer_poro_pinned_boundaries.push_back(Outer_poro_outlet_boundary_id);

 // Get the number of boundaries from the vector
 unsigned n_b=outer_poro_pinned_boundaries.size();

 // Set boundary conditions for the outer poro mesh
 for(unsigned ibound=0;ibound<n_b;ibound++)
  {
   unsigned num_nod=Outer_poro_mesh_pt->nboundary_node(
     outer_poro_pinned_boundaries[ibound]);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Outer_poro_mesh_pt->boundary_node_pt(
       outer_poro_pinned_boundaries[ibound],inod);

     // Pin in r and z
     nod_pt->pin(0);
     nod_pt->pin(1);

     // If this node stores edge flux dofs, pin these
     if(nod_pt->nvalue()==4)
      {
       nod_pt->pin(2);
       nod_pt->pin(3);
      }
    }
  } // end_of_loop_over_outer_poro_boundaries


 // Boundary conditions for the SSS fluid mesh
 // ------------------------------------------

 // Inlet
 {
  unsigned ibound=SSS_fluid_inflow_boundary_id;
  unsigned num_nod=SSS_fluid_mesh_pt->nboundary_node(ibound);
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    // Get pointer to node
    Node* nod_pt=SSS_fluid_mesh_pt->boundary_node_pt(ibound,inod);
    
    // Pin radial velocity (parallel in/outflow)
    nod_pt->pin(0);
        
    // If the nod is also on either the inner or outer poro-elsatic 
    // boundary then pin the axial velocity too
    if( (nod_pt->is_on_boundary(
          SSS_fluid_outer_poro_interface_near_inlet_boundary_id))||
        (nod_pt->is_on_boundary(
         SSS_fluid_inner_poro_interface_boundary_id)))
     {
      nod_pt->pin(1);
     }
   }
 }

 // Outlet -- pin everything
 {
  unsigned ibound=SSS_fluid_outflow_boundary_id;
  unsigned num_nod=SSS_fluid_mesh_pt->nboundary_node(ibound);
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    // Get pointer to node
    Node* nod_pt=SSS_fluid_mesh_pt->boundary_node_pt(ibound,inod);
    
    // Pin radial and axial velocity
    nod_pt->pin(0);
    nod_pt->pin(1);
   }
 }


 // Boundary conditions for the syrinx fluid mesh
 // ---------------------------------------------

 // Vector of the boundaries we will pin on the syrinx fluid mesh
 std::vector<unsigned> syrinx_fluid_pinned_boundaries;

 // Add the symmetry boundary to the vector
 syrinx_fluid_pinned_boundaries.push_back(Syrinx_fluid_symmetry_boundary_id);

 // Get the number of boundaries from the vector
 n_b=syrinx_fluid_pinned_boundaries.size();

 // Set the boundary conditions on the syrinx fluid mesh
 for(unsigned ibound=0;ibound<n_b;ibound++)
  {
   unsigned num_nod=Syrinx_fluid_mesh_pt->nboundary_node(
     syrinx_fluid_pinned_boundaries[ibound]);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Syrinx_fluid_mesh_pt->boundary_node_pt(
       syrinx_fluid_pinned_boundaries[ibound],inod);

     // Pin u_r
     nod_pt->pin(0);
    }
  }

 // ----------------------------------------------------------------------------

 // Pin relevant Lagrange multipliers in the
 // syrinx fluid <-> inner poro FSI surface meshe
  
 // Loop over the elements
 unsigned nel=Syrinx_fluid_inner_poro_FSI_surface_mesh_pt->nelement(); 
 for(unsigned e=0;e<nel;e++)
  {
   // Get a face element
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=dynamic_cast<LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*>
    (Syrinx_fluid_inner_poro_FSI_surface_mesh_pt->element_pt(e));
   
   // Get the number of nodes in element e
   unsigned n_node=traction_element_pt->nnode();
   
   // Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     // Get a node
     Node *nod_pt=traction_element_pt->node_pt(n);
     
     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
     
     // Get the index of the first nodal value associated with this
     // FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(
       Lagrange_id_syrinx_fluid_inner_poro);
     
     // Pin swirl component of Lagrange multiplier
     nod_pt->pin(first_index+2);
     
     // Pin radial component if on the symmetry boundary
     if(bnod_pt->is_on_boundary(Syrinx_fluid_symmetry_boundary_id))
      {       
       nod_pt->pin(first_index+0);
      }
    } //end_of_loop_over_nodes
  } //end_of_loop_over_elements


 // Pin relevant Lagrange multipliers in the
 // SSS fluid <-> inner poro FSI surface mesh
 
 // Loop over the elements
 nel=SSS_fluid_inner_poro_FSI_surface_mesh_pt->nelement();
 for(unsigned e=0;e<nel;e++)
  {
   // Get a face element
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=dynamic_cast<LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*>
    (SSS_fluid_inner_poro_FSI_surface_mesh_pt->element_pt(e));
   
   // Get the number of nodes in element e
   unsigned n_node=traction_element_pt->nnode();
   
   // Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     // Get a node
     Node *nod_pt=traction_element_pt->node_pt(n);
     
     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
     
     // Get the index of the first nodal value associated with this
     // FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(
       Lagrange_id_SSS_fluid_inner_poro);
     
     // Pin swirl component of Lagrange multiplier
     nod_pt->pin(first_index+2);
     
     // Pin radial and axial component if on inflow boundary
     if(bnod_pt->is_on_boundary(SSS_fluid_inflow_boundary_id))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
     
     // Pin radial and axial component if on outflow boundary
     if(bnod_pt->is_on_boundary(SSS_fluid_outflow_boundary_id))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
    } //end_of_loop_over_nodes
  } //end_of_loop_over_elements


 // Pin relevant Lagrange multipliers in each of the 
 // SSS fluid <-> outer poro FSI surface meshes [non-block]
 Vector<Mesh*> tmp_mesh_pt;
 tmp_mesh_pt.push_back(SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt);
 tmp_mesh_pt.push_back(SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt);
 unsigned n_mesh=tmp_mesh_pt.size();
 
 // Loop over the meshes
 for(unsigned i=0;i<n_mesh;i++)
  {
   // Get the number of elements
   unsigned nel=tmp_mesh_pt[i]->nelement();
   
   // Loop over the elements
   for(unsigned e=0;e<nel;e++)
    {
     // Get a face element
     LinearisedAxisymPoroelasticBJS_FSIElement
      <FLUID_ELEMENT,PORO_ELEMENT>*
      traction_element_pt=dynamic_cast<LinearisedAxisymPoroelasticBJS_FSIElement
      <FLUID_ELEMENT,PORO_ELEMENT>*>
      (tmp_mesh_pt[i]->element_pt(e));
     
     // Get the number of nodes in element e
     unsigned n_node=traction_element_pt->nnode();
     
     // Loop over the nodes
     for(unsigned n=0;n<n_node;n++)
      {
       // Get a node
       Node *nod_pt=traction_element_pt->node_pt(n);
       
       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
       
       // Get the index of the first nodal value associated with this
       // FaceElement
       unsigned first_index=
        bnod_pt->index_of_first_value_assigned_by_face_element(
         Lagrange_id_SSS_fluid_outer_poro);
       
       // Pin swirl component of Lagrange multiplier
       nod_pt->pin(first_index+2);
       
       // Pin radial and axial component if on inflow boundary
       if(bnod_pt->is_on_boundary(SSS_fluid_inflow_boundary_id))
        {
         nod_pt->pin(first_index+0);
         nod_pt->pin(first_index+1);
        }
       
       // Pin radial and axial component if on outflow boundary
       if(bnod_pt->is_on_boundary(SSS_fluid_outflow_boundary_id))
        {
         nod_pt->pin(first_index+0);
         nod_pt->pin(first_index+1);
        }
      } //end_of_loop_over_nodes
    } //end_of_loop_over_elements
  } //end_of_loop_over_meshes


 // Pin relevant Lagrange multipliers in the
 // SSS fluid <-> inner poro FSI surface mesh [block]
 
 // Loop over the elements
 nel=SSS_fluid_outer_poro_block_FSI_surface_mesh_pt->nelement();
 for(unsigned e=0;e<nel;e++)
  {
   // Get a face element
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=dynamic_cast<LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*>
    (SSS_fluid_outer_poro_block_FSI_surface_mesh_pt->element_pt(e));
   
   // Get the number of nodes in element e
   unsigned n_node=traction_element_pt->nnode();
   
   // Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     // Get a node
     Node *nod_pt=traction_element_pt->node_pt(n);
     
     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
     
     // Get the index of the first nodal value associated with this
     // FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(
       Lagrange_id_SSS_fluid_outer_poro);
     
     // Pin swirl component of Lagrange multiplier
     nod_pt->pin(first_index+2);
     
    } //end_of_loop_over_nodes
  } //end_of_loop_over_elements

 
 // ----------------------------------------------------------------------------

 // Setup FSI
 setup_fsi();

 // Add the poro submeshes to the problem
 add_sub_mesh(Inner_poro_mesh_pt);
 add_sub_mesh(Outer_poro_mesh_pt);

 // Add the fluid submeshes to the problem
 add_sub_mesh(Syrinx_fluid_mesh_pt);
 add_sub_mesh(SSS_fluid_mesh_pt);

 // Add the fluid inflow traction mesh to the problem
 add_sub_mesh(Inflow_fluid_surface_mesh_pt);

 // Add the FSI meshes to the problem
 add_sub_mesh(Inner_poro_syrinx_fluid_FSI_surface_mesh_pt);
 add_sub_mesh(Syrinx_fluid_inner_poro_FSI_surface_mesh_pt);

 add_sub_mesh(Inner_poro_SSS_fluid_FSI_surface_mesh_pt);
 add_sub_mesh(SSS_fluid_inner_poro_FSI_surface_mesh_pt);

 add_sub_mesh(Outer_poro_SSS_fluid_near_inlet_FSI_surface_mesh_pt);
 add_sub_mesh(SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt);

 add_sub_mesh(Outer_poro_SSS_fluid_near_outlet_FSI_surface_mesh_pt);
 add_sub_mesh(SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt);

 add_sub_mesh(Outer_poro_SSS_fluid_block_FSI_surface_mesh_pt);
 add_sub_mesh(SSS_fluid_outer_poro_block_FSI_surface_mesh_pt);
 
 // Now build the global mesh
 build_global_mesh();
 
 // Assign equation numbers
 oomph_info << assign_eqn_numbers() << " equations assigned" << std::endl;
 
} // end_of_constructor

//===start_of_create_poro_face_elements==============================
/// Make poro face elements
//===================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
void AxisymmetricSpineProblem<PORO_ELEMENT, FLUID_ELEMENT>::
create_poro_face_elements()
{

 // ----------------------------------
 // Create the poro FSI surface meshes
 // ----------------------------------

 // Inner poro mesh (from syrinx fluid mesh)
 // ---------------------------------------
 unsigned b=Inner_poro_syrinx_interface_boundary_id;
 
// Now loop over the bulk elements and create the face elements
 unsigned nel = Inner_poro_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>
    (Inner_poro_mesh_pt->boundary_element_pt(b,e),
     Inner_poro_mesh_pt->face_index_at_boundary(b,e));
   
   // Add to the appropriate surface mesh
   Inner_poro_syrinx_fluid_FSI_surface_mesh_pt->add_element_pt(face_element_pt);
   
   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(b);
  }
 
 // Inner poro mesh (from SSS fluid mesh)
 // -------------------------------------
 b = Inner_poro_SSS_interface_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = Inner_poro_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>
    (Inner_poro_mesh_pt->boundary_element_pt(b,e),
     Inner_poro_mesh_pt->face_index_at_boundary(b,e));
   
   // Add to the appropriate surface mesh
   Inner_poro_SSS_fluid_FSI_surface_mesh_pt->add_element_pt(face_element_pt);
   
   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(b);
  }
  

 // Outer poro mesh near inlet (from SSS fluid mesh)
 // ------------------------------------------------
 b = Outer_poro_SSS_interface_near_inlet_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = Outer_poro_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>
    (Outer_poro_mesh_pt->boundary_element_pt(b,e),
     Outer_poro_mesh_pt->face_index_at_boundary(b,e));
   
   // Add to mesh
   Outer_poro_SSS_fluid_near_inlet_FSI_surface_mesh_pt->add_element_pt(
    face_element_pt);
   
   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(b);
  }


 // Outer poro mesh near outlet (from SSS fluid mesh)
 // ------------------------------------------------
 b = Outer_poro_SSS_interface_near_outlet_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = Outer_poro_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>
    (Outer_poro_mesh_pt->boundary_element_pt(b,e),
     Outer_poro_mesh_pt->face_index_at_boundary(b,e));
   
   // Add to mesh
   Outer_poro_SSS_fluid_near_outlet_FSI_surface_mesh_pt->add_element_pt(
    face_element_pt);
   
   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(b);
  }


 // Outer poro mesh (block) (from SSS fluid mesh)
 // ---------------------------------------------
 b = Outer_poro_SSS_interface_block_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = Outer_poro_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies fluid traction to poroelastic
   // solid
   FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>* face_element_pt=
    new FSILinearisedAxisymPoroelasticTractionElement
    <PORO_ELEMENT,FLUID_ELEMENT>
    (Outer_poro_mesh_pt->boundary_element_pt(b,e),
     Outer_poro_mesh_pt->face_index_at_boundary(b,e));
   
   // Add to mesh
   Outer_poro_SSS_fluid_block_FSI_surface_mesh_pt->add_element_pt(
    face_element_pt);
   
   // Set the FSI parameter
   face_element_pt->q_pt()=&Global_Physical_Variables::Q;
   
   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   face_element_pt->set_boundary_number_in_bulk_mesh(b);
  }
 

} // end_of_create_poro_face_elements



//==start_of_create_fluid_traction_elements===============================
/// Create the fluid traction elements
//========================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
void AxisymmetricSpineProblem<PORO_ELEMENT, FLUID_ELEMENT>::
create_fluid_traction_elements()
{

 // Inflow:
 //--------
 unsigned bound=SSS_fluid_inflow_boundary_id;

 // Now loop over bulk elements and create the face elements
 unsigned nel=SSS_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that applies pressure boundary condition
   // at inflow
   AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>* traction_element_pt
    = new AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>
    (SSS_fluid_mesh_pt->boundary_element_pt(bound,e),
     SSS_fluid_mesh_pt->face_index_at_boundary(bound,e));

   // Add to mesh
   Inflow_fluid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the applied traction
   traction_element_pt->traction_fct_pt() =
    &Global_Physical_Variables::fluid_inflow_boundary_traction;
  }


 // -----------------------------------
 // Create the fluid FSI surface meshes
 // -----------------------------------


 // Syrinx fluid mesh (from inner poro mesh)
 // ----------------------------------------
 bound=Syrinx_fluid_inner_poro_interface_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = Syrinx_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes Beavers Joseph Saffman condition
   // at FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>
    (Syrinx_fluid_mesh_pt->boundary_element_pt(bound,e),
     Syrinx_fluid_mesh_pt->face_index_at_boundary(bound,e),
     Lagrange_id_syrinx_fluid_inner_poro);
   
   // Add to mesh
   Syrinx_fluid_inner_poro_FSI_surface_mesh_pt->add_element_pt(
    traction_element_pt);
   
   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);
   
   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;
   
   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }


 // SSS fluid mesh (from inner poro mesh)
 // ------------------------------------ 
 bound = SSS_fluid_inner_poro_interface_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = SSS_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes Beavers Joseph Saffman condition
   // at FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>
    (SSS_fluid_mesh_pt->boundary_element_pt(bound,e),
     SSS_fluid_mesh_pt->face_index_at_boundary(bound,e),
     Lagrange_id_SSS_fluid_inner_poro);
   
   // Add to mesh
   SSS_fluid_inner_poro_FSI_surface_mesh_pt->add_element_pt(
    traction_element_pt);
   
   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);
   
   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;
   
   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }


 // SSS fluid mesh near inlet (from outer poro mesh)
 // ------------------------------------------------
 bound=SSS_fluid_outer_poro_interface_near_inlet_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = SSS_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes Beavers Joseph Saffman condition at
   // FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>
    (SSS_fluid_mesh_pt->boundary_element_pt(bound,e),
     SSS_fluid_mesh_pt->face_index_at_boundary(bound,e),
     Lagrange_id_SSS_fluid_outer_poro);
   
   // Add to mesh
   SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt->add_element_pt(
    traction_element_pt);
   
   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);
   
   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;
   
   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }


 // SSS fluid mesh near outlet (from outer poro mesh)
 // --------------------------------------------------
 bound=SSS_fluid_outer_poro_interface_near_outlet_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = SSS_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes Beavers Joseph Saffman condition at
   // FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>
    (SSS_fluid_mesh_pt->boundary_element_pt(bound,e),
     SSS_fluid_mesh_pt->face_index_at_boundary(bound,e),
     Lagrange_id_SSS_fluid_outer_poro);
   
   // Add to mesh
   SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt->add_element_pt(
    traction_element_pt);
   
   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);
   
   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;
   
   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }


 // SSS fluid mesh (block) (from outer poro mesh)
 // --------------------------------------------------
 bound=SSS_fluid_outer_poro_interface_block_boundary_id;
 
 // Now loop over the bulk elements and create the face elements
 nel = SSS_fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element that imposes Beavers Joseph Saffman condition at
   // FSI interface
   LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>*
    traction_element_pt=new
    LinearisedAxisymPoroelasticBJS_FSIElement
    <FLUID_ELEMENT,PORO_ELEMENT>
    (SSS_fluid_mesh_pt->boundary_element_pt(bound,e),
     SSS_fluid_mesh_pt->face_index_at_boundary(bound,e),
     Lagrange_id_SSS_fluid_outer_poro);
   
   // Add to mesh
   SSS_fluid_outer_poro_block_FSI_surface_mesh_pt->add_element_pt(
    traction_element_pt);
   
   // Associate element with bulk boundary (to allow it to access the boundary
   // coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);
   
   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;
   
   // Set inverse slip coefficient
   traction_element_pt->inverse_slip_rate_coefficient_pt()=
    &Global_Physical_Variables::Inverse_slip_rate_coefficient;
  }
 
} // end_of_create_fluid_traction_elements

//==start_of_setup_fsi====================================================
/// Setup interaction between the poroelasticity and fluid meshes
//========================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
void AxisymmetricSpineProblem<PORO_ELEMENT, FLUID_ELEMENT>::
setup_fsi()
{

 // Inner poro <-> Syrinx fluid coupling
 // -----------------------------------
 
 // Setup syrinx fluid traction on inner poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,
   Syrinx_fluid_inner_poro_interface_boundary_id,
   Syrinx_fluid_mesh_pt,
   Inner_poro_syrinx_fluid_FSI_surface_mesh_pt);
 
 // Setup BJS bc for syrinx fluid from inner poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <PORO_ELEMENT,2>
  (this,
   Inner_poro_syrinx_interface_boundary_id,
   Inner_poro_mesh_pt,
   Syrinx_fluid_inner_poro_FSI_surface_mesh_pt);
 
 // Inner poro <-> SSS fluid coupling
 // --------------------------------

 // Setup SSS fluid traction on inner poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,
   SSS_fluid_inner_poro_interface_boundary_id,
   SSS_fluid_mesh_pt,
   Inner_poro_SSS_fluid_FSI_surface_mesh_pt);

 // Setup BJS bc for SSS fluid from inner poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <PORO_ELEMENT,2>
  (this,
   Inner_poro_SSS_interface_boundary_id,
   Inner_poro_mesh_pt,
   SSS_fluid_inner_poro_FSI_surface_mesh_pt);


 // Outer poro <-> SSS fluid coupling
 // ---------------------------------

 // Near inlet

 // Setup SSS fluid traction on outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,
   SSS_fluid_outer_poro_interface_near_inlet_boundary_id,
   SSS_fluid_mesh_pt,
   Outer_poro_SSS_fluid_near_inlet_FSI_surface_mesh_pt);


 // Setup BJS bc for SSS fluid from outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <PORO_ELEMENT,2>
  (this,
   Outer_poro_SSS_interface_near_inlet_boundary_id,
   Outer_poro_mesh_pt,
   SSS_fluid_outer_poro_near_inlet_FSI_surface_mesh_pt);


 // Near outlet

 // Setup SSS fluid traction on outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,
   SSS_fluid_outer_poro_interface_near_outlet_boundary_id,
   SSS_fluid_mesh_pt,
   Outer_poro_SSS_fluid_near_outlet_FSI_surface_mesh_pt);

 // Setup BJS bc for SSS fluid from outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <PORO_ELEMENT,2>
  (this,
   Outer_poro_SSS_interface_near_outlet_boundary_id,
   Outer_poro_mesh_pt,
   SSS_fluid_outer_poro_near_outlet_FSI_surface_mesh_pt);



 // Block

 // Setup SSS fluid traction on outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,
   SSS_fluid_outer_poro_interface_block_boundary_id,
   SSS_fluid_mesh_pt,
   Outer_poro_SSS_fluid_block_FSI_surface_mesh_pt);

 // Setup BJS bc for SSS fluid from outer poro wall elements
 Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <PORO_ELEMENT,2>
  (this,
   Outer_poro_SSS_interface_block_boundary_id,
   Outer_poro_mesh_pt,
   SSS_fluid_outer_poro_block_FSI_surface_mesh_pt);



}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class PORO_ELEMENT, class FLUID_ELEMENT>
void AxisymmetricSpineProblem<PORO_ELEMENT, FLUID_ELEMENT>::
doc_solution(DocInfo& doc_info, unsigned npts)
{
 ofstream some_file;
 char filename[100];

 // Output results in Tecplot format
 bool output_in_tecplot_format=true;
 bool output_in_paraview_format=true;

 if(output_in_tecplot_format)
  {
   // Output solution on inner poro mesh
   sprintf(filename,"%s/soln-inner-poro%i.dat",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Inner_poro_mesh_pt->output(some_file,npts);
   some_file.close();

   // Output solution on outer poro mesh
   sprintf(filename,"%s/soln-outer-poro%i.dat",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Outer_poro_mesh_pt->output(some_file,npts);
   some_file.close();

   // Output solution on syrinx fluid mesh
   sprintf(filename,"%s/soln-syrinx-fluid%i.dat",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Syrinx_fluid_mesh_pt->output(some_file,npts);
   some_file.close();

   // Output solution on SSS fluid mesh
   sprintf(filename,"%s/soln-sss-fluid%i.dat",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   SSS_fluid_mesh_pt->output(some_file,npts);
   some_file.close();
  }

 
 
 // Outer poro at regularly spaced points
 {
  sprintf(filename,"%s/regular_outer_poro%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Vector<double> x(2);
  Vector<double> s(2);
  unsigned n=Outer_poro_regularly_spaced_plot_point.size();
  for (unsigned i=0;i<n;i++)
   {
    // Pointer to element
    PORO_ELEMENT* el_pt=Outer_poro_regularly_spaced_plot_point[i].first;
    
    // Coordinates in it
    s=Outer_poro_regularly_spaced_plot_point[i].second;
    
    // Get coords
    el_pt->interpolated_x(s,x);
    
    some_file << x[0] << " "
              << x[1] << " "
              << el_pt->interpolated_u(s,0) << " "
              << el_pt->interpolated_u(s,1) << " "
              << el_pt->interpolated_q(s,0) << " "
              << el_pt->interpolated_q(s,1) << " "
              << el_pt->interpolated_p(s) << " "
              << std::endl;
   }
  some_file.close();
 }

 // Inner poro at regularly spaced points
 {
  sprintf(filename,"%s/regular_inner_poro%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Vector<double> x(2);
  Vector<double> s(2);
  unsigned n=Inner_poro_regularly_spaced_plot_point.size();
  for (unsigned i=0;i<n;i++)
   {
    // Pointer to element
    PORO_ELEMENT* el_pt=Inner_poro_regularly_spaced_plot_point[i].first;
    
    // Coordinates in it
    s=Inner_poro_regularly_spaced_plot_point[i].second;
    
    // Get coords
    el_pt->interpolated_x(s,x);
    
    some_file << x[0] << " "
              << x[1] << " "
              << el_pt->interpolated_u(s,0) << " "
              << el_pt->interpolated_u(s,1) << " "
              << el_pt->interpolated_q(s,0) << " "
              << el_pt->interpolated_q(s,1) << " "
              << el_pt->interpolated_p(s) << " "
              << std::endl;
   }
  some_file.close();
 }


 // Syrinx Fluid at regularly spaced points
 {
  sprintf(filename,"%s/regular_syrinx_fluid%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned npts=Syrinx_fluid_regularly_spaced_plot_point.size();
  Vector<double> x(2);
  Vector<double> s(2);
  Vector<double> veloc(3);
  for (unsigned i=0;i<npts;i++)
   {
    // Pointer to element
    FLUID_ELEMENT* el_pt=Syrinx_fluid_regularly_spaced_plot_point[i].first;
    
    // Coordinates in it
    s=Syrinx_fluid_regularly_spaced_plot_point[i].second;
    
    // Get coords
    el_pt->interpolated_x(s,x);
    
    // Get velocity
    el_pt->interpolated_u_axi_nst(s,veloc);
    
    some_file << x[0] << " "
              << x[1] << " "
              << veloc[0] << " "
              << veloc[1] << " "
              << veloc[2] << " "
              << el_pt->interpolated_p_axi_nst(s) << " "
              << std::endl;
   }
  some_file.close();
 }


 // SSS Fluid at regularly spaced points
 {
  sprintf(filename,"%s/regular_sss_fluid%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned npts=SSS_fluid_regularly_spaced_plot_point.size();
  Vector<double> x(2);
  Vector<double> s(2);
  Vector<double> veloc(3);
  for (unsigned i=0;i<npts;i++)
   {
    // Pointer to element
    FLUID_ELEMENT* el_pt=SSS_fluid_regularly_spaced_plot_point[i].first;
    
    // Coordinates in it
    s=SSS_fluid_regularly_spaced_plot_point[i].second;
    
    // Get coords
    el_pt->interpolated_x(s,x);
    
    // Get velocity
    el_pt->interpolated_u_axi_nst(s,veloc);
    
    some_file << x[0] << " "
              << x[1] << " "
              << veloc[0] << " "
              << veloc[1] << " "
              << veloc[2] << " "
              << el_pt->interpolated_p_axi_nst(s) << " "
              << std::endl;
   }
  some_file.close();
 }

// hierher reinstate once they've been combined
 
 // // Output inner traction elements
 // sprintf(filename,"%s/traction-inner%i.dat",doc_info.directory().c_str(),
 //   doc_info.number());
 // some_file.open(filename);
 // Inner_surface_mesh_pt->output(some_file,npts);
 // some_file.close();

 // // Output outer traction elements
 // sprintf(filename,"%s/traction-outer%i.dat",doc_info.directory().c_str(),
 //   doc_info.number());
 // some_file.open(filename);
 // Outer_surface_mesh_pt->output(some_file,npts);
 // some_file.close();

 // Output results in Paraview format
 if(output_in_paraview_format)
  {
   // Output solution on inner poro mesh in Paraview format (with padded 0s)
   sprintf(filename,"%s/soln-inner-poro%05i.vtu",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Inner_poro_mesh_pt->output_paraview(some_file,npts);
   some_file.close();

   // Output solution on outer poro mesh in Paraview format (with padded 0s)
   sprintf(filename,"%s/soln-outer-poro%05i.vtu",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Outer_poro_mesh_pt->output_paraview(some_file,npts);
   some_file.close();

   // Output solution on syrinx fluid mesh in Paraview format (with padded 0s)
   sprintf(filename,"%s/soln-syrinx-fluid%05i.vtu",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   Syrinx_fluid_mesh_pt->output_paraview(some_file,npts);
   some_file.close();

   // Output solution on SSS fluid mesh in Paraview format (with padded 0s)
   sprintf(filename,"%s/soln-sss-fluid%05i.vtu",doc_info.directory().c_str(),
     doc_info.number());
   some_file.open(filename);
   SSS_fluid_mesh_pt->output_paraview(some_file,npts);
   some_file.close();
  }

 // Bump counter
 doc_info.number()++;

} // end_of_doc_solution



//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Reynolds number
 CommandLineArgs::specify_command_line_flag(
   "--re",
   &Global_Physical_Variables::Re);

 // Multiplier for permeability
 CommandLineArgs::specify_command_line_flag(
   "--permeability_multiplier",
   &Global_Physical_Variables::Permeability_multiplier);

 // Pin darcy?
 CommandLineArgs::specify_command_line_flag(
   "--pin_darcy");

 // No wall inertia
 CommandLineArgs::specify_command_line_flag(
   "--suppress_wall_inertia");

 // Start from rest?
 CommandLineArgs::specify_command_line_flag(
   "--steady_solve");

 // Target element area for fluid meshes
 CommandLineArgs::specify_command_line_flag(
   "--el_area_fluid",
   &Global_Physical_Variables::Element_area_fluid);

 // Target element area for poro meshes
 CommandLineArgs::specify_command_line_flag(
   "--el_area_poro",
   &Global_Physical_Variables::Element_area_poro);

 // Axial scaling parameter -- shrinks aspect ratio of
 // domain by this factor
 CommandLineArgs::specify_command_line_flag(
   "--z_shrink_factor",
   &Global_Physical_Variables::Z_shrink_factor);


 // Axial scaling parameter -- shrinks aspect ratio of
 // domain by this factor
 CommandLineArgs::specify_command_line_flag(
  "--undo_z_scaling_by_element_stretching");

 // Output directory
 CommandLineArgs::specify_command_line_flag(
   "--dir",
   &Global_Physical_Variables::Directory);

 // Timestep
 double dt=0.00001;
 CommandLineArgs::specify_command_line_flag(
   "--dt",
   &dt);

 // Set timestep from number of steps per unit distance of wavetravel
 unsigned nstep_per_unit_wave_travel=0;
 CommandLineArgs::specify_command_line_flag("--nsteps_per_unit_wave_travel",
                                            &nstep_per_unit_wave_travel);

 // Number of timesteps to perform
 unsigned nstep=1000;
 CommandLineArgs::specify_command_line_flag(
   "--nstep",
   &nstep);

 // Validation?
 CommandLineArgs::specify_command_line_flag(
   "--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set up doc info
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory(Global_Physical_Variables::Directory);

 // Create problem
 AxisymmetricSpineProblem<TAxisymmetricPoroelasticityElement<1>,
                          AxisymmetricTTaylorHoodElement> problem;

 // Update dependent problem parameters
 Global_Physical_Variables::update_dependent_parameters();

 // Doc dependent parameters
 Global_Physical_Variables::doc_dependent_parameters();

 if (CommandLineArgs::command_line_flag_has_been_set
     ("--nsteps_per_unit_wave_travel"))
  {
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dt"))
    {
     throw OomphLibError(
      "Cannot set both --dt and --nsteps_per_unit_wave_travel",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     dt=1.0/Global_Physical_Variables::Dura_wavespeed/
      double(nstep_per_unit_wave_travel);
     oomph_info << "Setting dt so that we perform "
                << nstep_per_unit_wave_travel 
                << " timesteps in the in the time it takes \n"
                << "a pulswave travelling with the dura-based wavespeed\n"
                << "one diameter: dt="
                << dt << std::endl;
    }
  }


 // Set the initial time to t=0
 problem.time_pt()->time()=0.0;

 // Set up impulsive start from rest (which also initialises the timestep)
 problem.assign_initial_values_impulsive(dt);

 // Do only 3 steps if validating
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   nstep=3;
  }

 // Check if we are doing a steady solve only
 if(CommandLineArgs::command_line_flag_has_been_set(
     "--steady_solve"))
  {
   oomph_info << "Doing steady solve only\n";
   
   // Do a steady solve
   problem.steady_newton_solve();
   
   // Output solution
   problem.doc_solution(doc_info);
   exit(0);
  }


 // Doc the initial conditions
 problem.doc_solution(doc_info);

 // Timestep
 oomph_info << "About to do " << nstep << " timesteps with dt = "
            << dt << std::endl;
 // Do the timestepping
 for(unsigned istep=0;istep<nstep;istep++)
  {
   oomph_info << "Solving at time " << problem.time_pt()->time()+dt
              << std::endl;

   // Solve for this timestep
   problem.unsteady_newton_solve(dt);

   // Doc the solution
   problem.doc_solution(doc_info);
  }

 return 0;
} // end_of_main


//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Driver for an axisymmetric poro-elastic test problem

#include "axisym_poroelasticity.h"
#include "meshes/triangle_mesh.h"


using namespace oomph;
using namespace std;

//========================================================================
// NOTE: CODE FOR ADAPTIVITY IS IN HERE BUT IT'S TYPICALLY NOT WORTHWHILE
//       IT'S TOO MUCH OF A HYPERBOLIC PROBLEM...
//========================================================================


//==================================================================
/// Namespace for exact solution and problem parameters
//==================================================================
namespace ProblemParameters
{

 /// Steady flag
 bool Steady=false;

 /// Timestep
 double Dt=0.01;

 /// Timescale ratio (non-dim density)
 double Lambda_sq=0.7;
 
 /// Permeability
 double Permeability=1.0; // ONLY 1.0 MAKES SENSE IN A PURE PORO-ELASTIC 
                          // PROBLEM!

 /// Non-dim Young's modulus -- set to one implying that equations
 /// were scaled on actual young's modulus.
 double E_mod=1.0;
 
 /// Poisson's ratio
 double Nu=0.3;

 /// Alpha, the Biot parameter
 double Alpha=0.5;

 /// Porosity
 double Porosity=0.3;

 /// Ratio of pore fluid density to solid matrix density
 double Density_ratio=0.6;

 /// Scaling parameter for pressure load on 
 /// upper boundary (for non-validation case)
 double P=1.0;

 /// Parameter for tanh origin for pressure incrementation
 double T_tanh=0.25; 

 /// Steepness parameter for tanh for pressure incrementation
 double Alpha_tanh=100.0;

 // First Lame parameter -- dependent parameter!
 double Lambda = 0.0;

 // Second Lame parameter -- dependent parameter!
 double Mu = 0.0;

 /// Precalculate ratio of pore fluid density to compound density 
 /// for convenience -- dependent parameter!
 double Rho_f_over_rho = 0.0;

 /// Helper function to update dependent parameters
 void update_dependent_parameters()
 {
  Lambda = E_mod*Nu/((1.0+Nu)*(1.0-2.0*Nu));
  Mu = E_mod/(2.0*(1.0+Nu));
  Rho_f_over_rho = Density_ratio/(1.0+Porosity*(Density_ratio-1.0));
 }

 /// Radius of the smaller arcs in the curved mesh
 double Domain_radius=1.0;

 // Inner radius of annular region
 double Inner_radius=0.3;

 /// Function that returns zero (for assigment of initial conditions/
 /// history values)
 double zero_fct(const double &time,
                 const Vector<double> &x)
 {
  return 0.0;
 }

 /// Imposed boundary displacement in r-direction
 double boundary_displ_0(const double &time,
                         const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    return 1.0/3.0*pow(time,3)*(r*z*z);
   }
  else
   {
    return 0.0;
   }
 }

 /// Imposed boundary displacement in z-direction
 double boundary_displ_1(const double &time,
                         const Vector<double> &x)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     const double r=x[0],z=x[1];
     return 1.0/3.0*pow(time,3)*(-z/r);
    }
   else
    {
     return 0.0;
    }
  }
 /// Imposed boundary velocity in r-direction
 double boundary_veloc_0(const double &time,
                         const Vector<double> &x)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     const double r=x[0],z=x[1];
     return time*time*(r*z*z);
    }
   else
    {
     return 0.0;
    }
  }

 /// Imposed boundary velocity in z-direction
 double boundary_veloc_1(const double &time,
                         const Vector<double> &x)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     const double r=x[0],z=x[1];
     return time*time*(-z/r);
    }
   else
    {
     return 0.0;
    }
  }

 /// Imposed boundary accel in r-direction
 double boundary_accel_0(const double &time,
                         const Vector<double> &x)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     const double r=x[0],z=x[1];
     return 2*time*(r*z*z);
    }
   else
    {
     return 0.0;
    }
  }

 /// Imposed boundary accel in z-direction
 double boundary_accel_1(const double &time,
                         const Vector<double> &x)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     const double r=x[0],z=x[1];
     return 2*time*(-z/r);
    }
   else
    {
     return 0.0;
    }
  }

 /// Imposed boundary flux in r-direction
 double boundary_flux_0(const double &time,
                        const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    return 1.0/Permeability*Alpha*time*time*(1+r*z);
   }
  else
   {
    return 0.0;
   }
 }

 /// Imposed boundary flux in z-direction
 double boundary_flux_1(const double &time,
                        const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double z=x[1];
    return 1.0/Permeability*Alpha*time*time*(-2.0/3.0*pow(z,3)-z*z);
   }
  else
   {
    return 0.0;
   }
 }



 /// Imposed boundary d/dt flux in r-direction
 double boundary_dfluxdt_0(const double &time,
                           const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    return 1.0/Permeability*2.0*Alpha*time*(1+r*z);
   }
  else
   {
    return 0.0;
   }
 }
 
 /// Imposed boundary d/dt flux in z-direction
 double boundary_dfluxdt_1(const double &time,
                           const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double z=x[1];
    return 1.0/Permeability*2.0*Alpha*time*(-2.0/3.0*pow(z,3)-z*z);
   }
  else
   {
    return 0.0;
   }
 }
 

 /// Imposed boundary d2/dt2 flux in r-direction
 double boundary_d2fluxdt2_0(const double &time,
                             const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    return 1.0/Permeability*2.0*Alpha*(1+r*z);
   }
  else
   {
    return 0.0;
   }
 }
 
 /// Imposed boundary d2/dt2 flux in z-direction
 double boundary_d2fluxdt2_1(const double &time,
                             const Vector<double> &x)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double z=x[1];
    return 1.0/Permeability*2.0*Alpha*(-2.0/3.0*pow(z,3)-z*z);
   }
  else
   {
    return 0.0;
   }
 }
  
 // Exact solution: u1,u2,q1,q2,div_q,p
 void exact_soln(const double &time,
                 const Vector<double> &x,
                 Vector<double> &soln)
  {
   if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     // u[0] -- radial displacement
     soln[0]=boundary_displ_0(time,x);

     // u[1] -- axial displacement
     soln[1]=boundary_displ_1(time,x);

     // Convert coordinates
     double r=x[0];
     double z=x[1];

     // q[0] -- radial flux
     soln[2]=boundary_flux_0(time,x);

     // q[1] -- axial flux
     soln[3]=boundary_flux_1(time,x);
     
     // div q
     soln[4]=1.0/Permeability*Alpha*pow(time,2)*(1/r-2*z*z);

     // p -- pore pressure
     soln[5]=time*time*(r*r*z-2*z);

     // dudt[0] radial veloc
     soln[6]=boundary_veloc_0(time,x);

     // dudt[1] axial veloc
     soln[7]=boundary_veloc_1(time,x);

     // div of solid veloc
     soln[8]=boundary_veloc_0(time,x)+boundary_veloc_1(time,x);

     // d2udt2[0] radial accel
     soln[9]=boundary_accel_0(time,x);

     // d2udt2[1] axial accel
     soln[10]=boundary_accel_1(time,x);

     // dqdt[0] -- radial flux accel
     soln[11]=boundary_dfluxdt_0(time,x);

     // dqdt[1] -- axial flux accel
     soln[12]=boundary_dfluxdt_1(time,x);
    }
   else
    {
     unsigned n=soln.size();
     for (unsigned i=0;i<n;i++)
      {
       soln[i]=0.0;
      }
    }
  }

 /// Solid body force
 void Solid_body_force(const double &time,
                       const Vector<double> &x,
                       Vector<double> &b)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    b[0]=(time/3.0)*(6*r*time*z*Alpha+6*
                     (Rho_f_over_rho*Alpha+
                      r*z*(z+Rho_f_over_rho*Alpha))*Lambda_sq-
                     2*r*time*time*Mu-time*time*(Lambda+Mu)/(r*r));
    b[1]=(-time/(3.0*pow(r,3)))*
     (-3*pow(r,5)*time*Alpha+6*r*r*z*Lambda_sq-
      time*time*z*Mu+pow(r,3)*
      (6*time*Alpha+
       2*Rho_f_over_rho*z*z*(3+2*z)*Alpha*Lambda_sq+4*time*time*z*(Lambda+Mu)));
   }
  else
   {
    b[0]=0;
    b[1]=0;
   }
 }
 
 /// Fluid body force
 void Fluid_body_force(const double &time,
                       const Vector<double> &x,
                       Vector<double> &f)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    f[0]=(time*(2*Rho_f_over_rho*(1+r*z)*Alpha*Lambda_sq+
                Porosity*(1.0/Permeability*time*Alpha+r*z*
                          (2*time+1.0/Permeability*time*Alpha+
                           2*Rho_f_over_rho*z*Lambda_sq))))/
     (Porosity*Rho_f_over_rho);
    f[1]=(time*(-2*r*Rho_f_over_rho*z*z*(3+2*z)*Alpha*Lambda_sq+
                Porosity*(3*pow(r,3)*time-r*time*
                          (6+1.0/Permeability*z*z*(3+2*z)*Alpha)-
                          6*Rho_f_over_rho*z*Lambda_sq)))/
     (3*Porosity*r*Rho_f_over_rho);
   }
  else
   {
    f[0]=0;
    f[1]=0;
   }
 }
 

 /// Source term for continuity
 void Mass_source(const double &time, const Vector<double> &x, double &f)
 {
  f=0;
 }


 /// Get time-dep pressure magnitude
 double pressure_magnitude(const double& time)
 {
  return P*0.5*(1.0+tanh(Alpha_tanh*(time-T_tanh)));
 }

 /// Pressure around the boundary of the domain
 void boundary_pressure(const double &time,
                        const Vector<double> &x,
                        const Vector<double> &n,
                        double &result)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    Vector<double> soln(13);
    exact_soln(time,x,soln);
    result=soln[5];
   }
  else
   {
    result=pressure_magnitude(time);
   }
 }
 

 /// Boundary traction
 void boundary_traction(const double &time,
                        const Vector<double> &x,
                        const Vector<double> &n,
                        Vector<double> &traction)
 {
  if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
    const double r=x[0],z=x[1];
    
    traction[0]=
     n[0]*(1.0/3.0*Lambda*pow(time,3)*(2*z*z-1/r)+
           2.0/3.0*Mu*pow(time,3)*z*z-Alpha*time*time*(r*r*z-2*z))+
     n[1]*1.0/3.0*Mu*pow(time,3)*(z/(r*r)+2*r*z);
    traction[1]=
     n[0]*1.0/3.0*Mu*pow(time,3)*(z/(r*r)+2*r*z)+
     n[1]*(1.0/3.0*Lambda*pow(time,3)*(2*z*z-1/r)-
           2.0/3.0*Mu*pow(time,3)/r-Alpha*time*time*(r*r*z-2*z));
   }
  // Pure pressure load
  else
   {
    double p=0.0;
    boundary_pressure(time,x,n,p);
    traction[0]=-p*n[0];
    traction[1]=-p*n[1];
   }
 }
 
 // Target element area for Triangle
 double Element_area = 0.01;
 
 // The directory in which the solution is output
 std::string Directory("RESLT");

 /// Pointer to timestepper for internal dofs
 TimeStepper* Internal_time_stepper_pt=0;

 /// Global function that completes the edge sign setup --
 /// has to be called before projection in unstructured
 /// adaptation
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

    // Set the internal q dofs' timestepper to the problem timestepper
    el_pt->set_q_internal_timestepper(Internal_time_stepper_pt);

   } // end of loop over elements
 }
}




//===================================================================
/// Problem class
//===================================================================
template<class ELEMENT,class TIMESTEPPER>
class AxiPoroProblem : public Problem
{
public:

 /// Constructor
 AxiPoroProblem();

 /// Set the time-dependent boundary values
 void set_boundary_values();

 /// Set the initial conditions
 void set_initial_condition()
  {

   if(!CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
     assign_initial_values_impulsive();
    }
   else
    {
     Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
      initial_value_fct(2);
     Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
      initial_veloc_fct(2);
     Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
      initial_accel_fct(2);
     
     // Assign values for analytical value, veloc and accel:
     initial_value_fct[0]=&ProblemParameters::boundary_displ_0;
     initial_value_fct[1]=&ProblemParameters::boundary_displ_1;
     
     initial_veloc_fct[0]=&ProblemParameters::boundary_veloc_0;
     initial_veloc_fct[1]=&ProblemParameters::boundary_veloc_1;
     
     initial_accel_fct[0]=&ProblemParameters::boundary_accel_0;
     initial_accel_fct[1]=&ProblemParameters::boundary_accel_1;
     
     TIMESTEPPER* timestepper_pt=dynamic_cast<TIMESTEPPER*>(time_stepper_pt());
     
     unsigned n_node=Bulk_mesh_pt->nnode();
     for(unsigned n=0;n<n_node;n++)
      {
       Node *nod_pt=Bulk_mesh_pt->node_pt(n);
       timestepper_pt->assign_initial_data_values(nod_pt,
                                                  initial_value_fct,
                                                  initial_veloc_fct,
                                                  initial_accel_fct);
      }
    }
  }

 /// Doc the solution
 void doc_solution(const unsigned &label);

 /// Actions before newton solve (empty)
 void actions_before_newton_solve(){}

 /// Actions after newton solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before implicit timestep -- update boundary conditions
 void actions_before_implicit_timestep()
  {
   set_boundary_values();
  }

 /// Complete problem setup
 void complete_problem_setup();


#ifdef ADAPTIVE

 /// Actions before adapt
 void actions_before_adapt()
  {
   // Loop over the surface elements
   unsigned n_element = Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_mesh_pt->flush_element_and_node_storage();
   
   //  Rebuild the global mesh 
   rebuild_global_mesh();
  }
 
 /// Actions after adapt
 void actions_after_adapt()
  {
   // Re-create pressure bc elements
   create_pressure_elements();

   //  Rebuild the global mesh 
   rebuild_global_mesh();

   // Complete problem setup
   complete_problem_setup();
  }

#endif

 /// Access to timestepper
 TIMESTEPPER* my_time_stepper_pt()
  {
   return My_time_stepper_pt;
  }

private:

 /// Allocate traction/pressure elements
 void create_pressure_elements();
 
#ifdef ADAPTIVE

 /// Pointer to the refineable "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif
 
 /// Mesh for traction/pressure elements
 Mesh *Surface_mesh_pt;

 /// Trace file
 std::ofstream Trace_file;

 // Remember timestepper
 TIMESTEPPER* My_time_stepper_pt;

};

//==============================================================================
/// Constructor
//==============================================================================
template<class ELEMENT,class TIMESTEPPER>
AxiPoroProblem<ELEMENT,TIMESTEPPER>::AxiPoroProblem()
{
 
 // Problem is linear
 Problem::Problem_is_nonlinear=false;
 
 // Allocate the timestepper
 My_time_stepper_pt=new TIMESTEPPER;
 add_time_stepper_pt(My_time_stepper_pt);

 // disable warnings in assign initial data values
 My_time_stepper_pt->disable_warning_in_assign_initial_data_values();
 
 // Remember it
 ProblemParameters::Internal_time_stepper_pt=My_time_stepper_pt;
  
 // Solid mesh
 //-----------
 
 TriangleMeshClosedCurve *outer_boundary_pt=0;
 
 // Internal boundary to prevent problem with flux B.C.s
 Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt;
 
 // Temporary storage for coords in required format
 Vector<Vector<double> > temp_coord(2,Vector<double>(2));
 
 // Bottom small arc
 Circle* outer_boundary_lower_circle_pt=
  new Circle(ProblemParameters::Inner_radius,
             -2.0*ProblemParameters::Domain_radius,
             ProblemParameters::Domain_radius);
 
 // Large arc
 Circle* outer_boundary_right_circle_pt=
  new Circle(ProblemParameters::Inner_radius,0.0,sqrt(5)*
             ProblemParameters::Domain_radius);
 
 // Top small arc
 Circle* outer_boundary_upper_circle_pt=
  new Circle(ProblemParameters::Inner_radius, 
             2.0*ProblemParameters::Domain_radius,
             ProblemParameters::Domain_radius);
 
 // Storage for the individual boundary sections
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(4);
 
 // First bit
 double zeta_start=MathematicalConstants::Pi/2.0;
 double zeta_end=0;
 unsigned nsegment=
  (unsigned)max(
   3.0,
   MathematicalConstants::Pi*ProblemParameters::Domain_radius/
   (2*std::sqrt(ProblemParameters::Element_area)));
 outer_curvilinear_boundary_pt[0]=
  new TriangleMeshCurviLine(outer_boundary_lower_circle_pt, zeta_start,
                            zeta_end, nsegment, 0);
 
 // Second bit
 zeta_start=-atan(2.0);
 zeta_end=atan(2.0);
 nsegment=
  (unsigned)max(
   3.0,
   (2.0*MathematicalConstants::Pi*sqrt(5.0)*ProblemParameters::Domain_radius*
    (atan(2.0)/MathematicalConstants::Pi))/
   sqrt(ProblemParameters::Element_area));
 outer_curvilinear_boundary_pt[1]=
  new TriangleMeshCurviLine(outer_boundary_right_circle_pt, zeta_start,
                            zeta_end, nsegment, 1);
 
 // Third bit
 zeta_start=2.0*MathematicalConstants::Pi;
 zeta_end=3.0*MathematicalConstants::Pi/2.0;
 nsegment=
  (unsigned)max(
   3.0,
   MathematicalConstants::Pi*ProblemParameters::Domain_radius
   /(2.0*std::sqrt(ProblemParameters::Element_area)));
 outer_curvilinear_boundary_pt[2]=
  new TriangleMeshCurviLine(outer_boundary_upper_circle_pt, zeta_start,
                            zeta_end, nsegment, 2);
 
 // Fourth bit
 temp_coord[0][0]=ProblemParameters::Inner_radius;
 temp_coord[0][1]=ProblemParameters::Domain_radius;
 temp_coord[1][0]=ProblemParameters::Inner_radius;
 temp_coord[1][1]=-ProblemParameters::Domain_radius;
 
 outer_curvilinear_boundary_pt[3]=
  new TriangleMeshPolyLine(temp_coord,3);
 
 // Make the solid outer boundary
 outer_boundary_pt=
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
 
 
 // Add four internal boundaries emanating from corners to make sure
 // that no elements have faces on two boundaries.
 inner_open_boundary_pt.resize(4);
 
 // First inner boundary
 temp_coord[0][0]=ProblemParameters::Inner_radius;
 temp_coord[0][1]=-ProblemParameters::Domain_radius;
 temp_coord[1][0]=ProblemParameters::Inner_radius+
  0.25*ProblemParameters::Domain_radius;
 temp_coord[1][1]=-0.75*ProblemParameters::Domain_radius;
 
 Vector<TriangleMeshCurveSection*> inner_open_polyline1_curve_section_pt(1);
 inner_open_polyline1_curve_section_pt[0]=
  new TriangleMeshPolyLine(temp_coord,4);
 
 inner_open_polyline1_curve_section_pt[0]->
  connect_initial_vertex_to_curviline(
   dynamic_cast<TriangleMeshCurviLine*>(outer_curvilinear_boundary_pt[0]),
   MathematicalConstants::Pi/2.0);
 
 inner_open_boundary_pt[0]=
  new TriangleMeshOpenCurve(inner_open_polyline1_curve_section_pt);
 
 // Second inner boundary
 temp_coord[0][0]=ProblemParameters::Inner_radius+
  ProblemParameters::Domain_radius;
 temp_coord[0][1]=-2.0*ProblemParameters::Domain_radius;
 temp_coord[1][0]=ProblemParameters::Inner_radius+
  1.3*ProblemParameters::Domain_radius;
 temp_coord[1][1]=-1.5*ProblemParameters::Domain_radius;
 
 TriangleMeshPolyLine *inner_open_polyline2_pt=
  new TriangleMeshPolyLine(temp_coord,5);
 
 inner_open_polyline2_pt->connect_initial_vertex_to_curviline(
  dynamic_cast<TriangleMeshCurviLine*>(outer_curvilinear_boundary_pt[1]),
  -atan(2.0));
 
 Vector<TriangleMeshCurveSection*> inner_open_polyline2_curve_section_pt(1);
 inner_open_polyline2_curve_section_pt[0]=inner_open_polyline2_pt;
 
 inner_open_boundary_pt[1]=
  new TriangleMeshOpenCurve(inner_open_polyline2_curve_section_pt);
 
 // Third inner boundary
 temp_coord[0][0]=ProblemParameters::Inner_radius+
  ProblemParameters::Domain_radius;
 temp_coord[0][1]=2.0*ProblemParameters::Domain_radius;
 temp_coord[1][0]=ProblemParameters::Inner_radius+
  1.3*ProblemParameters::Domain_radius;
 temp_coord[1][1]=1.5*ProblemParameters::Domain_radius;
 
 TriangleMeshPolyLine *inner_open_polyline3_pt=
  new TriangleMeshPolyLine(temp_coord,6);
 
 inner_open_polyline3_pt->connect_initial_vertex_to_curviline(
  dynamic_cast<TriangleMeshCurviLine*>(outer_curvilinear_boundary_pt[1]),
  atan(2.0));
 
 Vector<TriangleMeshCurveSection*> inner_open_polyline3_curve_section_pt(1);
 inner_open_polyline3_curve_section_pt[0]=inner_open_polyline3_pt;
 
 inner_open_boundary_pt[2]=
  new TriangleMeshOpenCurve(inner_open_polyline3_curve_section_pt);
 
 // Fourth inner boundary
 temp_coord[0][0]=ProblemParameters::Inner_radius;
 temp_coord[0][1]=ProblemParameters::Domain_radius;
 temp_coord[1][0]=ProblemParameters::Inner_radius+
  0.25*ProblemParameters::Domain_radius;
 temp_coord[1][1]=0.75*ProblemParameters::Domain_radius;
 
 Vector<TriangleMeshCurveSection*> inner_open_polyline4_curve_section_pt(1);
 inner_open_polyline4_curve_section_pt[0]=
  new TriangleMeshPolyLine(temp_coord,7);
 
 inner_open_polyline4_curve_section_pt[0]->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>(outer_curvilinear_boundary_pt[3]),
  0);
 
 inner_open_boundary_pt[3]=
  new TriangleMeshOpenCurve(inner_open_polyline4_curve_section_pt);
 
 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
  outer_boundary_pt);
 
 // Set the inner boundaries
 triangle_mesh_parameters.internal_open_curves_pt()=
  inner_open_boundary_pt;
 
 // Target area for triangle
 triangle_mesh_parameters.element_area() = ProblemParameters::Element_area;

#ifdef ADAPTIVE
 
 // Build adaptive "bulk" mesh
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                  time_stepper_pt());

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Choose error tolerances 
 Bulk_mesh_pt->min_permitted_error()=1.0e-4;
 Bulk_mesh_pt->max_permitted_error()=1.0e-3;
 
 Bulk_mesh_pt->max_element_size()=ProblemParameters::Element_area;
 Bulk_mesh_pt->min_element_size()=1.0e-20;

 // Actions before projection
 Bulk_mesh_pt->mesh_update_fct_pt()=
  &ProblemParameters::edge_sign_setup<ELEMENT>;
 
#else
 
 // Create the buk mesh
 Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                          time_stepper_pt());
 
#endif
  
 // Setup the signs for the fluxes
 ProblemParameters::edge_sign_setup<ELEMENT>(Bulk_mesh_pt);

 // Create the surface mesh
 Surface_mesh_pt = new Mesh;

 // Assign the traction/pressure elements
 create_pressure_elements();

 // Add the submeshes to the global mesh
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Build the global mesh
 build_global_mesh();

 // Complete the problem set up
 complete_problem_setup();

 // Temporary filename string
 char filename[500];

 // Trace filename
 sprintf(filename,"%s/trace.dat",ProblemParameters::Directory.c_str());

 // Open the trace file
 Trace_file.open(filename);

 // Assign and doc number of eqns
 oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
 
} // end of problem constructor


//===================================================================
/// Complete problem setup
//===================================================================
template<class ELEMENT,class TIMESTEPPER>
void AxiPoroProblem<ELEMENT,TIMESTEPPER>::complete_problem_setup()
{

 // Loop over all elements to set parameters
 //-----------------------------------------
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   el_pt->nu_pt()=&ProblemParameters::Nu;
   el_pt->solid_body_force_fct_pt()=ProblemParameters::Solid_body_force;
   el_pt->fluid_body_force_fct_pt()=ProblemParameters::Fluid_body_force;
   el_pt->mass_source_fct_pt()=ProblemParameters::Mass_source;
   el_pt->lambda_sq_pt()=&ProblemParameters::Lambda_sq;
   el_pt->density_ratio_pt()=&ProblemParameters::Density_ratio;
   el_pt->permeability_pt()=&ProblemParameters::Permeability;
   el_pt->alpha_pt()=&ProblemParameters::Alpha;
   el_pt->porosity_pt()=&ProblemParameters::Porosity;

   // Set the internal q dofs' timestepper to the problem timestepper
   el_pt->set_q_internal_timestepper(time_stepper_pt());

  } // end of loop over elements


 // Pin flux and solid displacement along selected boundaries
 // ----------------------------------------------------------
  std::vector<unsigned> pinned_boundaries;

 // Leave outflow boundary traction free ("do nothing") for "real"
 // run. Pin and apply flux and displacement for validation case
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   pinned_boundaries.push_back(0);
  }
 pinned_boundaries.push_back(1);
 pinned_boundaries.push_back(3);
 unsigned n_boundary=pinned_boundaries.size();
 for(unsigned ibound=0;ibound<n_boundary;ibound++)
  {
   unsigned n_node = Bulk_mesh_pt->nboundary_node(pinned_boundaries[ibound]);
   for(unsigned i=0;i<n_node;i++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(pinned_boundaries[ibound],i);

     // Pin all of them: i=0: u_r; i=1: u_z. Subsequent entries
     // only exist for midside nodes where the remaining entries
     // represent the Darcy flux at the "flux interpolation points".
     unsigned n_value=nod_pt->nvalue();
     for(unsigned i=0;i<n_value;i++)
      {
       nod_pt->pin(i);
      }
    }
  }
}


//===start_of_traction===============================================
/// Make traction elements along the top boundary of the bulk mesh
//===================================================================
template<class ELEMENT,class TIMESTEPPER>
void AxiPoroProblem<ELEMENT,TIMESTEPPER>::create_pressure_elements()
{
 // How many bulk elements are next to boundary 2 (the top boundary)?
 unsigned ibound=2;
 unsigned n_neigh = Bulk_mesh_pt->nboundary_element(ibound);
 
 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   AxisymmetricPoroelasticityTractionElement<ELEMENT>* pressure_element_pt
    = new AxisymmetricPoroelasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(ibound,n),
     Bulk_mesh_pt->face_index_at_boundary(ibound,n));
   
   // Add to mesh
   Surface_mesh_pt->add_element_pt(pressure_element_pt);

   // Set function pointers
   pressure_element_pt->pressure_fct_pt()=&ProblemParameters::boundary_pressure;
   pressure_element_pt->traction_fct_pt()=&ProblemParameters::boundary_traction;
  }

} // end of assign_traction_elements


//===================================================================
/// Set the time dependent boundary values
//===================================================================
template<class ELEMENT,class TIMESTEPPER>
void AxiPoroProblem<ELEMENT,TIMESTEPPER>::set_boundary_values()
{
 Vector<double> x(2);
 
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  value_fct(2);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  veloc_fct(2);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  accel_fct(2);
 
 // Assign values for analytical value, veloc and accel:
 value_fct[0]=&ProblemParameters::boundary_displ_0;
 value_fct[1]=&ProblemParameters::boundary_displ_1;
 
 veloc_fct[0]=&ProblemParameters::boundary_veloc_0;
 veloc_fct[1]=&ProblemParameters::boundary_veloc_1;
 
 accel_fct[0]=&ProblemParameters::boundary_accel_0;
 accel_fct[1]=&ProblemParameters::boundary_accel_1;
 
 
 // Get the nodal index at which values representing edge fluxes
 // at flux interpolation points are stored
 ELEMENT *el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
 
 // How many flux interpolation points do we have?
 unsigned n=el_pt->nedge_flux_interpolation_point();
 
 // Provide storage; null out the ones we don't need
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  flux_fct(2+n,0);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  flux_ddt_fct(2+n,0);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  flux_d2dt2_fct(2+n,0);
  
 // Where are the flux values stored?
 Vector<unsigned> q_index(n);
 for (unsigned j=0;j<n;j++)
  {
   q_index[j]=el_pt->q_edge_index(j);
   flux_fct[q_index[j]]=&ProblemParameters::zero_fct;
   flux_ddt_fct[q_index[j]]=&ProblemParameters::zero_fct;
   flux_d2dt2_fct[q_index[j]]=&ProblemParameters::zero_fct;
  }
 

 TIMESTEPPER* timestepper_pt=
  dynamic_cast<TIMESTEPPER*>(time_stepper_pt());
 
 // Assign current (and history) values for solid displacements as well as
 // ----------------------------------------------------------------------
 // flux along pinned boundaries
 //-----------------------------
 std::vector<unsigned> pinned_boundaries;

 // Leave outer boundary traction free ("do nothing") for "real"
 // run. Pin and apply flux and displacement for validation case
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   pinned_boundaries.push_back(0);
  }
 pinned_boundaries.push_back(1);
 pinned_boundaries.push_back(3);
 unsigned n_boundary=pinned_boundaries.size();
 for(unsigned ibound=0;ibound<n_boundary;ibound++)
  {
   
   // Assign current (and history) values for solid displacements
   //------------------------------------------------------------
   // (ignores the flux ones)
   //------------------------
   unsigned n_node = Bulk_mesh_pt->nboundary_node(pinned_boundaries[ibound]);
   for(unsigned i=0;i<n_node;i++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(pinned_boundaries[ibound],i);
     timestepper_pt->assign_initial_data_values(nod_pt,
                                                value_fct,
                                                veloc_fct,
                                                accel_fct);
    }
   
   // Assign flux
   //------------
   
   // Coordinate vector
   Vector<double> x(2);
   
   // Get the number of elements along boundary 
   unsigned n_boundary_element=
    Bulk_mesh_pt->nboundary_element(pinned_boundaries[ibound]);
   
   // Loop over the elements along boundary ibound
   for(unsigned e=0;e<n_boundary_element;e++)
    {
     // Upcast the current element to the actual type
     ELEMENT *el_pt=
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                             boundary_element_pt(pinned_boundaries[ibound],e));
     
     // Loop over the edges
     for(unsigned edge=0;edge<3;edge++)
      {
       // Get pointer to node that stores the edge flux dofs for this edge
       Node* nod_pt=el_pt->edge_flux_node_pt(edge);
       
       // Set values for the flux degrees of freedom
       if (nod_pt->is_on_boundary(pinned_boundaries[ibound]))
        {
         // Get face index of face associated with edge
         unsigned f=el_pt->face_index_of_edge(edge);
         
         // Build a temporary face element from which we'll extract
         // the outer unit normal
         AxisymmetricPoroelasticityTractionElement<ELEMENT>* face_el_pt=
          new AxisymmetricPoroelasticityTractionElement<ELEMENT>(el_pt,f);   
         
         // Loop over the flux interpolation points
         unsigned n_flux_interpolation_points=
          el_pt->nedge_flux_interpolation_point();
         for(unsigned g=0;g<n_flux_interpolation_points;g++)
          {
           // Get the global coords of the flux_interpolation point
           el_pt->edge_flux_interpolation_point_global(edge,g,x);
           
           // Get the exact solution
           Vector<double> exact_soln(13,0.0);
           ProblemParameters::exact_soln(time(),x,exact_soln);
           
           // Get unit normal at this flux interpolation point
           Vector<double> s(1);
           el_pt->face_local_coordinate_of_flux_interpolation_point(edge,g,s);
           Vector<double> unit_normal(2);
           face_el_pt->outer_unit_normal(s,unit_normal);
           
#ifdef PARANOID          
           // Sanity check
           Vector<double> x_face(2);
           face_el_pt->interpolated_x(s,x_face);
           if ((x_face[0]-x[0])*(x_face[0]-x[0])+
               (x_face[1]-x[1])*(x_face[1]-x[1])>1.0e-3)
            {
             std::stringstream error;
             error 
              << "Discrepancy in coordinate of flux interpolation point\n"
              << "(computed by bulk and face elements) for edge " << e 
              << " and flux int pt " << g << "\n"
              << "Face thinks node is at: "
              << x_face[0] << " " << x_face[1] << "\n"
              << "Bulk thinks node is at: "
              << x[0] << " " << x[1] << "\n";
             throw OomphLibError(
              error.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
            }
#endif
           if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
            {
             // Set the boundary flux -- only does the current values;
             // should really do the history values, but it's messy...
             nod_pt->set_value(q_index[g],
                               exact_soln[2]*unit_normal[0]+
                               exact_soln[3]*unit_normal[1]);
            }
           else
            {
             // assign zero 
             timestepper_pt->assign_initial_data_values(nod_pt,
                                                        flux_fct,
                                                        flux_ddt_fct,
                                                        flux_d2dt2_fct);
            }
           
          } // End of loop over flux interpolation points
         
         // Don't need face element on that edge any more
         delete face_el_pt;
         
        } // End if for edge on required boundary
      } // End of loop over edges
     
    } // End of loop over boundary elements
  } // End of loop over boundaries
 }



//===================================================================
/// Write the solution and exact solution to file, and calculate the error
//===================================================================
template<class ELEMENT,class TIMESTEPPER>
void AxiPoroProblem<ELEMENT,TIMESTEPPER>::doc_solution(const unsigned &label)
{

 oomph_info << "Outputting solution " << label 
            << " for time " << time_pt()->time() << std::endl;

 ofstream some_file;
 char filename[100];

 unsigned npts=5;
 sprintf(filename,"%s/soln%i.dat",ProblemParameters::Directory.c_str(),label);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output coarse solution
 //-----------------------
 unsigned npts_coarse=2;
 sprintf(filename,"%s/coarse_soln%i.dat",ProblemParameters::Directory.c_str(),label);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Output boundary condition elements
 //-----------------------------------
 sprintf(filename,"%s/bc_elements%i.dat",ProblemParameters::Directory.c_str(),label);
 some_file.open(filename);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   
   // Output exact solution
   //----------------------
   sprintf(filename,"%s/exact_soln%i.dat",ProblemParameters::Directory.c_str(),label);
   some_file.open(filename);
   Bulk_mesh_pt->output_fct(some_file,
                            npts,
                            time(),
                            ProblemParameters::exact_soln);
   some_file.close();
   
   // Doc error
   //----------
   Vector<double> norm(3,0.0);
   Vector<double> error(3,0.0);
   sprintf(filename,"%s/error%i.dat",ProblemParameters::Directory.c_str(),label);
   some_file.open(filename);
   Bulk_mesh_pt->compute_error(some_file,
                               ProblemParameters::exact_soln,
                               time(),
                               error,norm);
   some_file.close();
   
   // Doc error norm:
   oomph_info << std::endl;
   oomph_info << "Norm of exact u : " << sqrt(norm[0]) << std::endl;
   oomph_info << "Norm of exact q : " << sqrt(norm[1]) << std::endl;
   oomph_info << "Norm of exact p : " << sqrt(norm[2]) << std::endl 
              << std::endl;
   
   oomph_info << "Norm of u error : " << sqrt(error[0]) << std::endl;
   oomph_info << "Norm of q error : " << sqrt(error[1]) << std::endl;
   oomph_info << "Norm of p error : " << sqrt(error[2]) << std::endl;
   oomph_info << std::endl << std::endl;
   
   Trace_file << ndof() << " "
              << ProblemParameters::Element_area << " "
              << ProblemParameters::Dt << " "
              << sqrt(error[0]) << " "
              << sqrt(error[1]) << " "
              << sqrt(error[2]) << " "
              << sqrt(norm[0]) << " "
              << sqrt(norm[1]) << " "
              << sqrt(norm[2]) << std::endl;
  }
 else
  {
   Trace_file << time_pt()->time() << " "
              << ProblemParameters::pressure_magnitude(time_pt()->time()) << " "
              << std::endl;
  }

}

//===================================================================
/// Main function
//===================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &ProblemParameters::Directory);

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Target element area for triangle
 CommandLineArgs::specify_command_line_flag("--element_area",
                                            &ProblemParameters::Element_area);

 // Steady solve?
 CommandLineArgs::specify_command_line_flag("--steady");

 // Biot parameter
 CommandLineArgs::specify_command_line_flag("--alpha",
                                            &ProblemParameters::Alpha);

 // Density ratio (fluid to solid)
 CommandLineArgs::specify_command_line_flag("--density_ratio",
                                            &ProblemParameters::Density_ratio);

 // Intertia parameter
 CommandLineArgs::specify_command_line_flag("--lambda_sq",
                                            &ProblemParameters::Lambda_sq);

 // Timestep
 CommandLineArgs::specify_command_line_flag("--dt",
                                            &ProblemParameters::Dt);

 // Steepness parameter for tanh increment of pressure load 
 CommandLineArgs::specify_command_line_flag("--alpha_tanh",
                                            &ProblemParameters::Alpha_tanh);

 // "Time" for tanh increment of pressure load 
 CommandLineArgs::specify_command_line_flag("--t_tanh",
                                            &ProblemParameters::T_tanh);

 // Precalculate n_steps based on defaults
 unsigned n_steps=unsigned(1.0/ProblemParameters::Dt);

 // Override n_steps if provided
 CommandLineArgs::specify_command_line_flag("--n_steps",
                                            &n_steps);

 // Max number of adaptations per Newton solve
 unsigned max_adapt=0;
 CommandLineArgs::specify_command_line_flag("--max_adapt",
                                            &max_adapt);

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Don't override the value of Lambda^2 set on the command line
 if(!CommandLineArgs::command_line_flag_has_been_set("--lambda_sq"))
  {
   // And if it's steady
   if(CommandLineArgs::command_line_flag_has_been_set("--steady"))
    {
     // Switch on steady flag
     ProblemParameters::Steady=true;

     // Disable inertia
     ProblemParameters::Lambda_sq=0;
    }
   else // Unsteady
    {
     // Disable steady flag
     ProblemParameters::Steady=false;

     // Set default inertia parameter value
     ProblemParameters::Lambda_sq=1e-5;
    }
  }

 /// Set dependent parameters
 ProblemParameters::update_dependent_parameters();
  
#ifdef ADAPTIVE

 // The problem instance
 AxiPoroProblem<ProjectableAxisymmetricPoroelasticityElement<
  TAxisymmetricPoroelasticityElement<1> >, Newmark<1> > problem;

#else

 // The problem instance
 AxiPoroProblem<TAxisymmetricPoroelasticityElement<1>, Newmark<1> > problem;

#endif

 // Doc initial configuration
 problem.doc_solution(0);

 // If we're doing a steady solve
 if(CommandLineArgs::command_line_flag_has_been_set("--steady"))
  {
   // Set time to 1 so that the exact solution reduces to the steady solution
   problem.time_pt()->time()=1.0;

   // Set the BCs
   problem.set_boundary_values();

   // Do the steady solve
   problem.steady_newton_solve();

   // Doc the solution
   problem.doc_solution(1);
  }
 else
  {
   // Set the initial time to t=0
   problem.time_pt()->time()=0.0;

   // Initialise the timestep
   problem.initialise_dt(ProblemParameters::Dt);

   // Set up the initial conditions
   problem.set_initial_condition();

#ifdef ADAPTIVE

   // We're doing the first timestep
   bool first=true;

#endif

   // Loop over timesteps
   for(unsigned i=1;i<=n_steps;i++)
    {
     // Output current time
     oomph_info << "Solving at time "
                << problem.time_pt()->time()+ProblemParameters::Dt << std::endl;
     
     
#ifdef ADAPTIVE
     
     // Solve the problem with the adaptive Newton's method, allowing
     // up to max_adapt mesh adaptations after every solve.
     problem.unsteady_newton_solve(ProblemParameters::Dt,max_adapt,first);

     // Now we've done the first timestep...
     first=false;
          
#else

     // Do the unsteady solve
     problem.unsteady_newton_solve(ProblemParameters::Dt);

#endif

     // Doc the solution
     problem.doc_solution(i);
    }
  }
 
 return 0;
}


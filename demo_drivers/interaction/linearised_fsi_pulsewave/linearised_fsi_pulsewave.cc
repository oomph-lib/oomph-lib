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
//
// Linearised FSI problem -- flow in a cylindrical tube

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "axisym_linear_elasticity.h"
#include "meshes/triangle_mesh.h"

#define DO_FSI

using namespace std;
using namespace oomph;

//=================================================================
//A namespace for the physical parameters in the problem
//=================================================================
namespace Global_Physical_Variables
{

 /// For validation only
 bool Pin_fluid_on_fsi=false;

 /// The Reynolds number
 double Re = 10.0;

 /// The Strouhal number
 double St = 1.0;

 /// The Womersley number -- dependent parameter; compute
 /// from Re and St
 double Wo = 0.0;

 /// Define Poisson's ratio Nu
 double Nu = 0.3;

 /// Define the non-dimensional Young's modulus
 double E = 1.0;

 /// Density ratio (zero for zero wall inertia)
 double Density_ratio=0.0;

 /// Square of the nondimensional "density" -- dependent parameter; compute
 double Lambda_sq = 0.0;

 /// Inner radius of tube
 double Inner_radius=0.5;

 /// Outer radius of tube
 double Outer_radius=1.0;

 /// Length of tube
 double Length=3.0;

 /// Constant inlet pressure (for steady Poiseuille flow)
 double P_inlet_const=Length*4.0/(Inner_radius*Inner_radius);

 /// Timescale for ramping up pressure via cos profile (if zero
 /// we have an impulsive start)
 double T_cos=0.0;

 /// Pressure growth factor for pressure acting on the outside
 /// of the solid wall
 double P_outside_scale=0.0;

 /// FSI parameter
 double Q=1.0e-6;

 /// Traction applied to the outside of the solid mesh
 void outside_solid_boundary_traction(const double &time,
                                      const Vector<double> &x,
                                      const Vector<double> &n,
                                      Vector<double> &result)
 {
  result[0]=-P_outside_scale*time*time*n[0];
  result[1]=-P_outside_scale*time*time*n[1];
  result[2]=0;
 }

 /// Inflow traction applied to the fluid mesh
 void fluid_inflow_boundary_traction(const double &time,
                                     const Vector<double> &x,
                                     const Vector<double> &n,
                                     Vector<double> &result)
 {
  double ramp_factor=1.0;
  if (T_cos!=0.0)
   {
    if (time<=T_cos)
     {
      ramp_factor=0.5*(1.0-cos(2.0*MathematicalConstants::Pi*time/T_cos));
     }
   }
  //oomph_info << "ramp factor: " << ramp_factor << std::endl;
  result[0]=0.0;
  result[1]=P_inlet_const*ramp_factor;
  result[2]=0;
 }


 /// Constant wall pressure for validation
 double P_wall=1.0;

 /// Traction applied to the solid mesh at fsi interface (for validation only)
 void validation_solid_fsi_boundary_traction(const double &time,
                                             const Vector<double> &x,
                                             const Vector<double> &n,
                                             Vector<double> &result)
 {
  result[0]=-P_wall*n[0];
  result[1]=-P_wall*n[1];
  result[2]=0;
 }


 /// Shear stress (for steady Poiseuille flow)
 double Tau=2.0/Inner_radius;

 /// "fsi" traction applied to the fluid mesh (for validation case
 /// in which fluid is driven by prescribed traction on fsi boundary
 /// rather than "lagrange multiplier traction" that enforces no slip
 void validation_fluid_fsi_boundary_traction(const double &time,
                                             const Vector<double> &x,
                                             const Vector<double> &n,
                                             Vector<double> &result)
 {
  // Pressure increases linearly
  result[0]=-P_inlet_const*x[1]/Length;

  // Shear stress is constant
  result[1]=Tau;

  // No swirl
  result[2]=0.0;
 }

 /// Storage for Moens Korteweg wavespeed -- dependent parameter!
 double Wavespeed=0.0;

 /// Storage for pressure wavespeed in solid -- dependent parameter!
 double Pressure_wavespeed=0.0;

 /// "Exact" solution for propagating pulse wave
 void pulse_wave_solution(const double& time,
                          const Vector<double>& x,
                          Vector<double>& soln)
 {
  // Three velocities and pressure
  soln.resize(4,0.0);

  // Wave position
  double z_wave=Length-Wavespeed*time;

  // Just do the pressure
  if (x[1]>z_wave)
   {
    soln[3]=P_inlet_const;
   }
  else
   {
    soln[3]=0.0;
   }
 }


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

  // Non-dim wall thickness
  double h=Outer_radius-Inner_radius;

  //Wavespeed
  Wavespeed=sqrt(h/(2.0*Inner_radius)/(Q*Re));

  // The Womersley number
  Wo = Re*St;

  // Wall inertia parameter
  Lambda_sq=St*St*Re*Q*Density_ratio;

  // Pressure wavespeed in solid
  // Commented out to avoid a division by zero
  //Pressure_wavespeed=sqrt( ( (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu)) ) /
  //                         (Density_ratio*Q*Re) );

 }


 // Doc dependent parameters
 void doc_dependent_parameters()
 {
  oomph_info << std::endl;
  oomph_info << "Problem parameters" << std::endl;
  oomph_info << "==================" << std::endl;
  oomph_info << "Womersley number (ReSt)           : "
             << Wo << std::endl;
  oomph_info << "Wall inertia parameter (Lambda^2) : "
             << Lambda_sq << std::endl;
  oomph_info << "Moens-Korteweg wavespeed          : "
             << Wavespeed << std::endl;
  oomph_info << "Pressure wavespeed in solid       : "
             << Pressure_wavespeed << std::endl;
  oomph_info << std::endl;

 }


};


//==========================================================================
/// Problem class
//==========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
class PressureWaveFSIProblem : public Problem
{

public:

 /// Constructor
 PressureWaveFSIProblem();

 /// Create the fluid traction elements
 void create_fluid_traction_elements();

 /// Create the solid traction elements
 void create_solid_traction_elements();

 /// Update the problem specs before next timestep: Emtpy
 void actions_before_implicit_timestep() {}

 /// Doc solution
 void doc_solution(DocInfo& doc_info);

private:


 /// setup fsi
 void setup_fsi();

 /// Solid mesh
 TriangleMesh<FLUID_ELEMENT> *Fluid_mesh_pt;

 /// Fluid mesh
 TriangleMesh<SOLID_ELEMENT> *Solid_mesh_pt;

 /// Solid surface mesh on outside
 Mesh* Solid_surface_mesh_pt;

 /// Solid surface mesh at FSI interface
 Mesh* FSI_solid_surface_mesh_pt;

 /// Inflow fluid surface mesh
 Mesh* Inflow_fluid_surface_mesh_pt;

 /// FSI fluid surface mesh
 Mesh* FSI_fluid_surface_mesh_pt;

 /// Pointer to wall timestepper
 TimeStepper* Solid_time_stepper_pt;

 /// Pointer to fluid timestepper
 TimeStepper* Fluid_time_stepper_pt;

 /// Mesh as geom object representation of fluid mesh
 MeshAsGeomObject* Fluid_mesh_geom_obj_pt;

 /// Vector of pairs containing pointers to elements and
 /// local coordinates within them for regularly spaced plot points
 Vector<std::pair<FLUID_ELEMENT*,Vector<double> > > Regularly_spaced_plot_point;

 /// Id for Lagrange multiplier constraint
 unsigned Lagrange_id;

};


//============================================================================
/// Constructor
//============================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
PressureWaveFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::PressureWaveFSIProblem()
{

 // Create timesteppers
 Solid_time_stepper_pt = new NewmarkBDF<2>;
 add_time_stepper_pt(Solid_time_stepper_pt);
 Fluid_time_stepper_pt = new BDF<2>;
 add_time_stepper_pt(Fluid_time_stepper_pt);

 //Create the problem geometry - a hollow solid cylinder filled with fluid
 //-----------------------------------------------------------------------

 //Temporary storage for coords in required format
 Vector<Vector<double> > temp_coord(2,Vector<double>(2));

 //Fluid mesh
 //----------

 //Vector for the vertices
 Vector<Vector<double> > fluid_vertex_coords(5, Vector<double>(2));
 fluid_vertex_coords[0][0]=0.0;
 fluid_vertex_coords[0][1]=0.0;

 fluid_vertex_coords[1][0]=Global_Physical_Variables::Inner_radius;
 fluid_vertex_coords[1][1]=0.0;

 fluid_vertex_coords[2][0]=Global_Physical_Variables::Inner_radius;
 fluid_vertex_coords[2][1]=Global_Physical_Variables::Length;

 fluid_vertex_coords[3][0]=0.0;
 fluid_vertex_coords[3][1]=Global_Physical_Variables::Length;

 fluid_vertex_coords[4][0]=0.0;
 fluid_vertex_coords[4][1]=0.0;


 //Loop over the vertices and create a polyline between each consecutive pair
 Vector<TriangleMeshCurveSection*> fluid_outer_polyline_boundary_pt(4);
 for(unsigned i=0;i<4;i++)
  {
   temp_coord[0][0]=fluid_vertex_coords[i][0];
   temp_coord[0][1]=fluid_vertex_coords[i][1];
   temp_coord[1][0]=fluid_vertex_coords[i+1][0];
   temp_coord[1][1]=fluid_vertex_coords[i+1][1];
   fluid_outer_polyline_boundary_pt[i] =
    new TriangleMeshPolyLine(temp_coord,i);
  }

 //Create the outer boundary closed curve
 TriangleMeshClosedCurve* fluid_outer_boundary_pt =
  new TriangleMeshClosedCurve(fluid_outer_polyline_boundary_pt);

 //Target element area for Triangle
 double uniform_element_area = 0.002;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters fluid_triangle_mesh_parameters(
   fluid_outer_boundary_pt);

 // Take the maximum element area
 fluid_triangle_mesh_parameters.element_area() =
   uniform_element_area;

 // Create the mesh
 Fluid_mesh_pt = new TriangleMesh<FLUID_ELEMENT>(
   fluid_triangle_mesh_parameters, Fluid_time_stepper_pt);

 // Mesh as geom object representation of fluid mesh
 Fluid_mesh_geom_obj_pt=new MeshAsGeomObject(Fluid_mesh_pt);

 // Extract points
 unsigned nplot_r=5;
 unsigned nplot_z=unsigned(double(nplot_r)*Global_Physical_Variables::Length/
                           Global_Physical_Variables::Inner_radius);
 Regularly_spaced_plot_point.resize(nplot_r*nplot_z);
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
     Regularly_spaced_plot_point[count].first=
      dynamic_cast<FLUID_ELEMENT*>(geom_obj_pt);
     Regularly_spaced_plot_point[count].second=s;

     count++;
    }
  }


 //-----------------------------------------------------------------------------

 //Solid mesh
 //----------

 //Vector for the vertices
 Vector<Vector<double> > solid_vertex_coords(5, Vector<double>(2));
 solid_vertex_coords[0][0]=Global_Physical_Variables::Inner_radius;
 solid_vertex_coords[0][1]=0.0;

 solid_vertex_coords[1][0]=Global_Physical_Variables::Outer_radius;
 solid_vertex_coords[1][1]=0.0;

 solid_vertex_coords[2][0]=Global_Physical_Variables::Outer_radius;
 solid_vertex_coords[2][1]=Global_Physical_Variables::Length;

 solid_vertex_coords[3][0]=Global_Physical_Variables::Inner_radius;
 solid_vertex_coords[3][1]=Global_Physical_Variables::Length;

 solid_vertex_coords[4][0]=Global_Physical_Variables::Inner_radius;
 solid_vertex_coords[4][1]=0.0;


 //Loop over the vertices and create a polyline between each consecutive pair
 Vector<TriangleMeshCurveSection*> solid_outer_polyline_boundary_pt(4);
 for(unsigned i=0;i<4;i++)
  {
   temp_coord[0][0]=solid_vertex_coords[i][0];
   temp_coord[0][1]=solid_vertex_coords[i][1];
   temp_coord[1][0]=solid_vertex_coords[i+1][0];
   temp_coord[1][1]=solid_vertex_coords[i+1][1];
   solid_outer_polyline_boundary_pt[i] =
    new TriangleMeshPolyLine(temp_coord,i);
  }

 //Create the outer boundary closed curve
 TriangleMeshClosedCurve* solid_outer_boundary_pt =
  new TriangleMeshClosedCurve(solid_outer_polyline_boundary_pt);

 //Target element area for Triangle
 //uniform_element_area = 0.02;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters solid_triangle_mesh_parameters(
   solid_outer_boundary_pt);

 // Take the maximum element area
 solid_triangle_mesh_parameters.element_area() =
   uniform_element_area;

 // Create the mesh
 Solid_mesh_pt = new TriangleMesh<SOLID_ELEMENT>(
   solid_triangle_mesh_parameters, Solid_time_stepper_pt);

 //-----------------------------------------------------------------------------

 // Create surface meshes
 Inflow_fluid_surface_mesh_pt = new Mesh;
 FSI_fluid_surface_mesh_pt = new Mesh;
 Solid_surface_mesh_pt = new Mesh;
 FSI_solid_surface_mesh_pt = new Mesh;

 // Set Id for Lagrange multiplier constraint
 Lagrange_id=99;

 // Make surface meshes
 create_fluid_traction_elements();
 create_solid_traction_elements();

 //-----------------------------------------------------------------------------

 //Loop over the elements in fluid mesh
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to the particular element type
   FLUID_ELEMENT *el_pt = dynamic_cast<FLUID_ELEMENT*>
    (Fluid_mesh_pt->element_pt(e));

   //There is no need for ALE
   el_pt->disable_ALE();

   //Set the Reynolds number for each element
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   //Set the product of Reynolds and Strouhal numbers
   el_pt->re_st_pt() = &Global_Physical_Variables::Wo;

   //Pin u_theta at all nodes (no swirl velocity)
   unsigned n_node = el_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node *nod_pt = el_pt->node_pt(n);
     nod_pt->pin(2);
    }
  }

 // Loop over the elements in solid mesh
 n_element = Solid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Cast to a bulk element
   SOLID_ELEMENT *el_pt = dynamic_cast<SOLID_ELEMENT*>
    (Solid_mesh_pt->element_pt(e));

   // Set the pointer to Poisson's ratio
   el_pt->nu_pt() = &Global_Physical_Variables::Nu;

   // Set the pointer to non-dim Young's modulus
   el_pt->youngs_modulus_pt() = &Global_Physical_Variables::E;

   // Set the pointer to the Lambda parameter
   el_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   // Loop over all nodes in the element and pin the theta displacement
   unsigned n_node = el_pt->nnode();
   for(unsigned i_node=0;i_node<n_node;i_node++)
    {
     Node *nod_pt = el_pt->node_pt(i_node);
     nod_pt->pin(2);
    }
  }// end_loop_over_elements


 //Set the fluid boundary conditions

 //Loop over the nodes on the fluid mesh boundary
 unsigned n_boundary = Fluid_mesh_pt->nboundary();

 //Pin the boundary nodes
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Pin radial velocity at in and outflow (parallel flow)
   // and on centreline (symm); leave FSI boundary (b=1) free
   // [Doing the loop like this ensures that we pin nodes
   // that are simultaneously on fsi boundary and in/outflow boundaries]
   if (b!=1)
    {
     unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;++n)
      {
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,n);
       nod_pt->pin(0);


       // If node is also on FSI boundary it must be on in/outflow:
       // Pin axial velocity too
       if (nod_pt->is_on_boundary(1))
        {
         nod_pt->pin(1);
        }
      }
    }
   // Completely pin the velocity on the FSI boundary
   else if (Global_Physical_Variables::Pin_fluid_on_fsi)
    {
     unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;++n)
      {
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,n);
       nod_pt->pin(0);
       nod_pt->pin(1);
      }
    }
  }


 //Set the solid boundary conditions
 n_boundary = Solid_mesh_pt->nboundary();

 //Pin the boundary nodes
 for(unsigned b=0;b<n_boundary;b++)
  {
   unsigned n_boundary_node = Solid_mesh_pt->nboundary_node(b);
   for(unsigned n=0;n<n_boundary_node;++n)
    {
     // Pin solid at in and outflow
     if ((b==0)||(b==2))
      {
       Solid_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Solid_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }
    }
  }


 // Pin theta Lagrange multipliers (no swirl) and radial and axial Lagrange
 // multiplier for nodes on in and outflow boundaries (b=0 and 2)
 // because velocity is already constrained there
 unsigned nel=FSI_fluid_surface_mesh_pt->nelement();
 for(unsigned e=0;e<nel;e++)
  {
   // Get face element
   LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_ELEMENT,SOLID_ELEMENT>*
    traction_element_pt=dynamic_cast<
    LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_ELEMENT,SOLID_ELEMENT>*>
    (FSI_fluid_surface_mesh_pt->element_pt(e));

   // Loop over nodes
   unsigned n_node = traction_element_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node *nod_pt = traction_element_pt->node_pt(n);

     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);

     // Get the index of the first nodal value associated with
     // this FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(Lagrange_id);

     // Pin swirl component
     nod_pt->pin(first_index+2);

     // Pin all Lagrange multipliers if we're pinning the fluid
     // on the fsi boundary (for validation)
     if (Global_Physical_Variables::Pin_fluid_on_fsi)
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }

     // Pin radial and axial component if on in or outflow boundary
     if (bnod_pt->is_on_boundary(0))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
     if (bnod_pt->is_on_boundary(2))
      {
       nod_pt->pin(first_index+0);
       nod_pt->pin(first_index+1);
      }
    }
  }


 // Setup fluid structure interaction
 setup_fsi();

 //Add submeshes and build global mesh
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Solid_mesh_pt);
 add_sub_mesh(Inflow_fluid_surface_mesh_pt);
 add_sub_mesh(FSI_fluid_surface_mesh_pt);
 add_sub_mesh(Solid_surface_mesh_pt);
 add_sub_mesh(FSI_solid_surface_mesh_pt);

 build_global_mesh();

 //Setup all the equation numbering and look-up schemes
 assign_eqn_numbers();

}


//========================================================================
/// Create the traction elements to the appropriate boundaries of the solid mesh
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::
create_solid_traction_elements()
{

 //Prescribed traction on fsi boundary (b=3)
 unsigned bound=3;

 // Now loop over bulk elements and create the face elements
 unsigned nel=Solid_mesh_pt->nboundary_element(bound);
 for(unsigned n=0;n<nel;n++)
  {

#ifdef DO_FSI

   // Create the face element
   FSIAxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT,FLUID_ELEMENT>*
    traction_element_pt = new
    FSIAxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT,FLUID_ELEMENT>
    (Solid_mesh_pt->boundary_element_pt(bound,n),
     Solid_mesh_pt->face_index_at_boundary(bound,n));

   // Add to mesh
   FSI_solid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the FSI parameter
   traction_element_pt->q_pt() = &Global_Physical_Variables::Q;

   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);

#else

   // Create the face element
   AxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT>*
    traction_element_pt
    = new AxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT>
    (Solid_mesh_pt->boundary_element_pt(bound,n),
     Solid_mesh_pt->face_index_at_boundary(bound,n));

   // Add to mesh
   FSI_solid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the applied traction
   traction_element_pt->traction_fct_pt() =
    &Global_Physical_Variables::validation_solid_fsi_boundary_traction;

#endif

  }


 //Prescribed traction on outside wall (b=1)
 bound=1;

 // Now loop over bulk elements and create the face elements
 nel=Solid_mesh_pt->nboundary_element(bound);
 for(unsigned n=0;n<nel;n++)
  {

   // Create the face element
   AxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT>*
    traction_element_pt
    = new AxisymmetricLinearElasticityTractionElement<SOLID_ELEMENT>
    (Solid_mesh_pt->boundary_element_pt(bound,n),
     Solid_mesh_pt->face_index_at_boundary(bound,n));

   // Add to mesh
   Solid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the applied traction
   traction_element_pt->traction_fct_pt() =
    &Global_Physical_Variables::outside_solid_boundary_traction;
  }


}


//========================================================================
/// Create the traction elements to the appropriate boundaries of the fluid mesh
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::
create_fluid_traction_elements()
{
 // Inflow:
 unsigned bound=2;

 // Now loop over bulk elements and create the face elements
 unsigned nel=Fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {
   // Create the face element
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
 bound=1;

 // Now loop over bulk elements and create the face elements
 nel=Fluid_mesh_pt->nboundary_element(bound);
 for(unsigned e=0;e<nel;e++)
  {

#ifdef DO_FSI

   // Create the face element
   LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_ELEMENT,SOLID_ELEMENT>*
    traction_element_pt = new
    LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_ELEMENT,SOLID_ELEMENT>
    (Fluid_mesh_pt->boundary_element_pt(bound,e),
     Fluid_mesh_pt->face_index_at_boundary(bound,e),Lagrange_id);

   // Add to mesh
   FSI_fluid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Associate element with bulk boundary (to allow it to access
   // the boundary coordinates in the bulk mesh)
   traction_element_pt->set_boundary_number_in_bulk_mesh(bound);

   // Set Strouhal number
   traction_element_pt->st_pt()=&Global_Physical_Variables::St;

#else

   // Create the face element
   AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>* traction_element_pt
    = new AxisymmetricNavierStokesTractionElement<FLUID_ELEMENT>
    (Fluid_mesh_pt->boundary_element_pt(bound,e),
     Fluid_mesh_pt->face_index_at_boundary(bound,e));

   // Add to mesh
   FSI_fluid_surface_mesh_pt->add_element_pt(traction_element_pt);

   // Set the applied traction
   traction_element_pt->traction_fct_pt() =
    &Global_Physical_Variables::fluid_fsi_boundary_traction;

#endif

  }

}


//=====================start_of_setup_interaction======================
/// Setup interaction between two fields
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::setup_fsi()
{
 // Setup fluid traction on FSI wall elements
 //------------------------------------------
 unsigned boundary_in_fluid_mesh=1;
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
  <FLUID_ELEMENT,2>
  (this,boundary_in_fluid_mesh,Fluid_mesh_pt,FSI_solid_surface_mesh_pt);

 // Setup no slip bc for axisym Navier Stokes from adjacent axisym
 //---------------------------------------------------------------
 // elasticity elements
 //-------------------
  unsigned boundary_in_solid_mesh=3;
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh
   <SOLID_ELEMENT,2>(
    this,boundary_in_solid_mesh,Solid_mesh_pt,FSI_fluid_surface_mesh_pt);

}


//==========================================================================
/// Doc solution
//==========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PressureWaveFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::doc_solution(
 DocInfo& doc_info)
{

 //Define a string that we can set to be the name of the output file
 char filename[100];

 //Define an output filestream
 ofstream file;

 // Fluid
 sprintf(filename,"%s/soln-fluid%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Fluid_mesh_pt->output(file,5);
 file.close();

 // Fluid traction
 sprintf(filename,"%s/fluid_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Inflow_fluid_surface_mesh_pt->output(file,5);
 file.close();


 // Fluid traction
 sprintf(filename,"%s/fluid_fsi_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 FSI_fluid_surface_mesh_pt->output(file,5);
 file.close();


 // // Fluid traction
 // sprintf(filename,"%s/tmp_fluid_fsi_traction%i.dat",
 //         doc_info.directory().c_str(),
 //         doc_info.number());
 // file.open(filename);
 // tmp_FSI_fluid_surface_mesh_pt->output(file,5);
 // file.close();

 // Solid
 sprintf(filename,"%s/soln-solid%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Solid_mesh_pt->output(file,5);
 file.close();


 // Solid traction
 sprintf(filename,"%s/solid_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 FSI_solid_surface_mesh_pt->output(file,5);
 file.close();


 // Solid traction
 sprintf(filename,"%s/outside_solid_traction%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Solid_surface_mesh_pt->output(file,5);
 file.close();


 // Fluid at regularly spaced points
 sprintf(filename,"%s/regular_fluid%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 unsigned npts=Regularly_spaced_plot_point.size();
 Vector<double> x(2);
 Vector<double> s(2);
 Vector<double> veloc(3);
 for (unsigned i=0;i<npts;i++)
  {
   // Pointer to element
   FLUID_ELEMENT* el_pt=Regularly_spaced_plot_point[i].first;

   // Coordinates in it
   s=Regularly_spaced_plot_point[i].second;

   // Get coords
   el_pt->interpolated_x(s,x);

   // Get velocity
   el_pt->interpolated_u_axi_nst(s,veloc);

   file << x[0] << " "
        << x[1] << " "
        << veloc[0] << " "
        << veloc[1] << " "
        << veloc[2] << " "
        << std::endl;
  }
 file.close();


 // Output "exact" pulse wave solution
 //-----------------------------------
 sprintf(filename,"%s/pulse_wave%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 file.open(filename);
 Fluid_mesh_pt->output_fct(file,5,time_pt()->time(),
                           Global_Physical_Variables::pulse_wave_solution);
 file.close();

 // Bump up counter
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

 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create doc info object
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Create problem
 PressureWaveFSIProblem<AxisymmetricTTaylorHoodElement,
                        TAxisymmetricLinearElasticityElement<3> > problem;


 // Update dependent problem parameters
 Global_Physical_Variables::update_dependent_parameters();

 // Doc dependent parameters
  Global_Physical_Variables::doc_dependent_parameters();

 // Setup timestep based on Moens Korteweg velocity
 //------------------------------------------------

 // Desired number of timesteps for pressure wave to propagate
 // along tube
 unsigned nt=100;

 // Choose timestep
 double dt=Global_Physical_Variables::St*
  (Global_Physical_Variables::Length/double(nt))/
  Global_Physical_Variables::Wavespeed;
 oomph_info << "Timestep: " << dt << std::endl;

 // // Choose ramp-up over 30 timesteps
 // unsigned nramp=30;
 // Global_Physical_Variables::T_cos=double(nramp)*dt;

 //Set an impulsive start from rest
 problem.assign_initial_values_impulsive(dt);

 // Output initial condition
 problem.doc_solution(doc_info);

 // Do timestepping loop
 unsigned nstep=500;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   nstep=3;
  }
 for(unsigned n=0;n<nstep;n++)
  {
   //Solve the problem
   problem.unsteady_newton_solve(dt);

   // Output
   problem.doc_solution(doc_info);

   if (Global_Physical_Variables::Pin_fluid_on_fsi)
    {
     Global_Physical_Variables::P_inlet_const*=1.1;
     oomph_info << "incrementing P_inlet_const to "
                <<  Global_Physical_Variables::P_inlet_const << std::endl;
    }

  }

}

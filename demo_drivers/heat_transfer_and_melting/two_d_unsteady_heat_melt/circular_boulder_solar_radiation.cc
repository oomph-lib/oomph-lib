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

// The unsteady heat equations
#include "unsteady_heat.h"

// The solid elements
#include "solid.h"

// Contact stuff
#include "contact_elements.h"

// The melt equations
#include "heat_transfer_and_melt_elements.h"
#include "temporary_stefan_boltzmann_elements.h"

// Mesh
#include "meshes/triangle_mesh.h"

#define ADAPTIVE

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;


/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 

//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{
 
 /// Non-dim density for solid
 double Lambda_sq=0.0;

 /// Poisson's ratio for solid
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Radius of penetrator
 double Radius=0.2;

 /// Initial y position of centre of penetrator
 double Y_c_initial=1.05; // 1.1;

 /// Position of centre of penetrator
 Vector<double> Centre;

 /// Penetrator
 Penetrator* Penetrator_pt=0;

 /// Diffuse radiation (based on 300 W/m^2)
 double Diffuse_radiation=0.4147111183;

 /// Direct radiation for peak of 1000 W/m^2
 double Direct_radiation=1.382370394; // peak of 250W/m^2 gives: 0.3455925986;

 /// One day in our non-dim units
 double One_day=0.4084910633;

 /// Function that specifies atmospheric
 /// radiation in terms of directional solar flux vector
 /// (which has direction and magnitude; well, like most vectors, actually)
 /// and total diffusive radiation (which is later weighted by diffuse limiting
 /// angles). Input argument: time.
 void atmospheric_radiation(const double& time,
                            double& solar_flux_magnitude, 
                            Vector<double>& solar_flux_unit_vector, 
                            double& total_diffuse_radiation)
 {
  // Diffuse background radiation
  total_diffuse_radiation=Diffuse_radiation;

  // Solar flux magnitude; constant (hierher: probably needs further 
  // modulation, because this obviously ain't true...)
  solar_flux_magnitude=Direct_radiation;

  // t=0: sunrise: Sun is horizontal from the left
  solar_flux_unit_vector[0]= cos(2.0*MathematicalConstants::Pi*time/One_day);
  solar_flux_unit_vector[1]=-sin(2.0*MathematicalConstants::Pi*time/One_day);
 }


 /// Nondim thermal inertia for rock
 double Alpha_rock=0.001; 

 /// Nondim thermal conductivity for rock
 double Beta_rock=10.0; 


} // end of ProblemParameters


//======start_of_ExactSolution========================================
/// Namespace for exact solution 
//====================================================================
namespace ExactSolution
{

 /// Constant/initial temperature
 double U0=0.0;

 /// Growth rate for interface
 double Growth_rate=1.0;
 
 /// Frequency of co-sinusoidal oscillation of incoming heat flux
 /// (to assess suppression of re-freezing). Set to zero for validation.
 double Omega_cos=0.0;

 /// Exact solution as a Vector
 void get_exact_u_for_unsteady_heat_validation(const double& t, 
                                               const Vector<double>& x, 
                                               Vector<double>& u)
 {
  double X=x[0];
  double Y=x[1];
  u[0]=U0+t*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1);
 }

 /// Exact solution as a scalar
 void get_exact_u_for_unsteady_heat_validation(const double& t, 
                                               const Vector<double>& x, 
                                               double& u)
 {
  Vector<double> u_vect(1);
  get_exact_u_for_unsteady_heat_validation(t,x,u_vect);
  u=u_vect[0];
 }

 /// Source function to make it an exact solution 
 void get_source_for_unsteady_heat_validation(const double& t, 
                                              const Vector<double>& x, 
                                              double& source)
 {
  double X=x[0];
  double Y=x[1];
  source = -2.0*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)-2.0*t*t*t*Y*Y*
Growth_rate*(1.0-cos(2.0*X*0.3141592653589793E1))*cos(2.0*X*
0.3141592653589793E1)+4.0*t*t*t*t*Y*Y*Growth_rate*pow(cos(2.0*X*
0.3141592653589793E1),2.0)*0.3141592653589793E1*0.3141592653589793E1-8.0*t*t*t*
t*Y*Y*Growth_rate*pow(sin(2.0*X*0.3141592653589793E1),2.0)*0.3141592653589793E1
*0.3141592653589793E1-4.0*t*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)*0.3141592653589793E1*
0.3141592653589793E1+2.0*t*t*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)+4.0*t*t*Y*cos(2.0*X*
0.3141592653589793E1);

 }


 /// Flux into bulk required by the exact solution on a 
 /// boundary with outer unit normal n. No dependence on temperature u.
 void flux_into_bulk(const double& t,
                     const Vector<double>& x, 
                     const Vector<double>& n, 
                     const double& u,
                     double& flux)
 {
  double X=x[0];
  double Y=x[1];

  //The outer unit normal 
  double Nx =  n[0];
  double Ny =  n[1];
  
  //The flux in terms of the normal is
  flux=(2.0*t*t*t*t*Y*Y*Growth_rate*sin(2.0*X*0.3141592653589793E1)*
0.3141592653589793E1*cos(2.0*X*0.3141592653589793E1)-2.0*t*t*Y*Y*(Y-1.0+
Growth_rate*t*t*(1.0-cos(2.0*X*0.3141592653589793E1)))*sin(2.0*X*
0.3141592653589793E1)*0.3141592653589793E1)*Nx+(2.0*t*t*Y*(Y-1.0+Growth_rate*t*
t*(1.0-cos(2.0*X*0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)+t*t*Y*
Y*cos(2.0*X*0.3141592653589793E1))*Ny;

 }

 /// Total flux (into bulk + melt) required by the exact 
 /// solution on a boundary with outer unit
 /// normal n. No dependence on temperature u.
 void prescribed_flux_for_unsteady_heat_validation(const double& t,
                                                   const Vector<double>& x, 
                                                   const Vector<double>& n, 
                                                   const double& u,
                                                   double& flux)
 {
  double X=x[0];

  //The flux into bulk
  flux_into_bulk(t,x,n,u,flux);

  // Add melt flux
  double melt_flux=
   2.0*Growth_rate*t*(1.0-cos(2.0*X*0.3141592653589793E1))/sqrt(1.0+4.0
*Growth_rate*Growth_rate*t*t*t*t*pow(sin(2.0*X*0.3141592653589793E1),2.0)*
0.3141592653589793E1*0.3141592653589793E1);

  flux+=melt_flux*cos(Omega_cos*t);
 }

 /// Height of melting surface
 double melting_surface_height(const double& t, const double& x)
 {
  return 1.0-Growth_rate*t*t*(1.0-cos(2.0*x*0.3141592653589793E1));
 }
 
} // end of ExactSolution


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// Problem class
//====================================================================
template<class ELEMENT>
class SolarRadiationProblem : public Problem
{

public:

 /// Constructor
 SolarRadiationProblem();
 
 /// Destructor (empty)
 ~SolarRadiationProblem(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
  
 /// Actions before next timestep
 void actions_before_implicit_timestep()
  {
   // Amplitude of oscillation
   double amplitude=0.2; //-0.2; 
   
   ProblemParameters::Centre[1]=ProblemParameters::Y_c_initial-amplitude*
    (1.0+0.05*time_pt()->time())*
    0.5*(1.0-cos(2.0*MathematicalConstants::Pi*time_pt()->time()));
   
   oomph_info << "Solving for y_c = " 
              << ProblemParameters::Centre[1] << std::endl;

   // If update of limiting angles is skipped during Newton iteration
   // update it at least before every timestep
   if (CommandLineArgs::command_line_flag_has_been_set(
        "--skip_update_limiting_angles_during_newton_iteration"))
    {
     update_limiting_angles();
    }

  }


 /// Newton convergence check
 void actions_before_newton_convergence_check() 
  {
   // Sometimes the Newton iteration 
   // stagnates, probably because of the shadow jumping
   // over integration points. Not updating limiting angles
   // is a bit naughty but the error incurred should be
   // small-ish if timesteps (i.e. change in surface shape remain
   // small between solves).
   if (!CommandLineArgs::command_line_flag_has_been_set(
        "--skip_update_limiting_angles_during_newton_iteration"))
    {
     update_limiting_angles();
    }
  }


 /// Actions before adapt: wipe contact elements
 void actions_before_adapt() 
  {
   // Kill the  elements and wipe surface mesh
   delete_flux_elements();
   delete_contact_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   // Create elements
   create_flux_elements();
   create_contact_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();

   // Rebuild elements
   complete_problem_setup();
  }
 

/// Update limiting angles for diffuse radiation
 void update_limiting_angles()
  {
   // Update limiting angles for diffuse radiation, given the
   // Vector of pointers to face elements that make up the "upper boundary"
   // that can potentially shield the integration points from diffuse
   // radiation
   unsigned nel=Surface_flux_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
      Surface_flux_mesh_pt->
      element_pt(e))->update_limiting_angles(Ordered_shielding_node_pt);
    }
  } 



 // Setup shielding nodes for solar radiation
 void setup_shielding_nodes()
  {
   // Wipe
   Ordered_shielding_node_pt.clear();
   Ordered_shielding_face_element_pt.clear();

   // Left horizon
   // hierher Ordered_shielding_node_pt.push_back(Left_horizon_node_pt);
   
   // Associate elements with vertex nodes
   std::map<Node*,std::set<FiniteElement*> > adj_el_pt;
   unsigned nel=Surface_flux_mesh_pt->nelement();

   oomph_info << "Number of elements: " << nel << std::endl;

   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=Surface_flux_mesh_pt->finite_element_pt(e);
     unsigned nnod=el_pt->nnode();
     
     // Associate element with first vertex node
     Node* nod_pt=el_pt->node_pt(0);
     adj_el_pt[nod_pt].insert(el_pt);
     
     // Associate element with last vertex node
     nod_pt=el_pt->node_pt(nnod-1);
     adj_el_pt[nod_pt].insert(el_pt);
    }
   
   // Find node that is only associated with single element
   unsigned nfound=0;
   Node* first_node_pt=0;
   FiniteElement* first_element_pt=0;
   for (std::map<Node*,std::set<FiniteElement*> >::iterator it=
         adj_el_pt.begin();
        it!=adj_el_pt.end();it++)
    {
     if ((*it).second.size()==1)
      {
       if (nfound==0)
        {
         first_node_pt=(*it).first;
         first_element_pt=*(((*it).second).begin());
         nfound++;
        }
       else
        {
         Node* other_node_pt=(*it).first;
         if (other_node_pt->x(0)<=first_node_pt->x(0))
          {
           first_node_pt=other_node_pt;
           first_element_pt=*(((*it).second).begin());
          }
         break;
        }
      }
    }
   
   
   // Current element/node
   FiniteElement* current_element_pt=first_element_pt;
   Node* current_node_pt=first_node_pt;
   
   oomph_info << "FIRST NODE: " 
              << first_node_pt->x(0) << " " 
              << first_node_pt->x(1) << " " << std::endl;

   // Add it
   Ordered_shielding_node_pt.push_back(current_node_pt);
   Ordered_shielding_face_element_pt.push_back(current_element_pt);

   // Keep going until we reach the end
   unsigned n_associated_els=2;
   while (n_associated_els==2)
    {
     // Next node
     Node* next_node_pt=current_element_pt->node_pt(0);
     if (next_node_pt==current_node_pt)
      {
       unsigned nnod=current_element_pt->nnode();
       next_node_pt=current_element_pt->node_pt(nnod-1);
      }
     
     // Get next element
     std::set<FiniteElement*>::iterator it=adj_el_pt[next_node_pt].begin();
     FiniteElement* next_element_pt=(*it);
     if (next_element_pt==current_element_pt)
      {
       it++;
       next_element_pt=(*it);
      }
     
     // Add it
     Ordered_shielding_node_pt.push_back(next_node_pt);
     Ordered_shielding_face_element_pt.push_back(next_element_pt);

     // bump up
     current_element_pt=next_element_pt;
     current_node_pt=next_node_pt;
     
     n_associated_els=adj_el_pt[current_node_pt].size();
    }
   
   // Kill last one
   Ordered_shielding_face_element_pt.pop_back();
    
   oomph_info << "Done setup shielding nodes: " 
              << Ordered_shielding_node_pt.size() << std::endl;

   // Right horizon
   // hierher Ordered_shielding_node_pt.push_back(Right_horizon_node_pt);
  }

 /// Doc the solution
 void doc_solution();
 
 /// Dummy global error norm for adaptive time-stepping
 double global_temporal_error_norm(){return 0.0;}

private:

 
 /// Create contact elements
 void create_contact_elements()
  {
   // Which boundaries are in contact with penetrator?
   Vector<unsigned> contact_boundaries;
   contact_boundaries.push_back(Upper_contact_boundary_id);
   contact_boundaries.push_back(Lower_contact_boundary_id);
   unsigned n=contact_boundaries.size();
   for (unsigned bb=0;bb<n;bb++)
    {
     // How many bulk elements are adjacent to boundary b?
     unsigned b=contact_boundaries[bb];
     unsigned n_element = Bulk_mesh_pt->nboundary_element_in_region(b,1);
     
     // Loop over the bulk elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_in_region_pt(b,1,e));
       
       //What is the face index of element e along boundary b
       int face_index = Bulk_mesh_pt->face_index_at_boundary_in_region(b,1,e);
       
       // Id of additional dofs created by this element
       unsigned id=0; // hierher make member

       // Build the corresponding contact element
       NonlinearSurfaceContactElement<ELEMENT>* contact_element_pt = new 
        NonlinearSurfaceContactElement<ELEMENT>(bulk_elem_pt,face_index,id);
       
       //Add the contact element to the surface mesh
       Surface_contact_mesh_pt->add_element_pt(contact_element_pt);
       
      } //end of loop over bulk elements adjacent to boundary b
    }
  }



 /// Delete contact elements
 void delete_contact_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_contact_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_contact_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_contact_mesh_pt->flush_element_and_node_storage();
  }

 /// Create flux elements
 void create_flux_elements()
  {
   // Which boundaries are subject to incoming radiation
   Vector<unsigned> flux_boundaries;
   flux_boundaries.push_back(Upper_contact_boundary_id);
   flux_boundaries.push_back(Left_melt_boundary_id);
   flux_boundaries.push_back(Right_melt_boundary_id);
   unsigned n=flux_boundaries.size();
   for (unsigned bb=0;bb<n;bb++)
    {
     // How many bulk elements are adjacent to boundary b?
     unsigned b=flux_boundaries[bb];
     unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_pt(b,e));
       
       //What is the face index of element e along boundary b
       int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
       
       // Create flux element
       StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT> *el_pt = 
        new StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>(bulk_elem_pt,
                                                            face_index);
                     
       // Add to mesh
       Surface_flux_mesh_pt->add_element_pt(el_pt);
       
      } //end of loop over bulk elements adjacent to boundary b
    }
  }
 


 /// Delete flux elements
 void delete_flux_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_flux_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_flux_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_flux_mesh_pt->flush_element_and_node_storage();
  }


 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup()
  {
   
   
   // Set material properties for elements in rock region
   //---------------------------------------------------
   unsigned r=1; // hierher enumerate globally
   unsigned nel=Bulk_mesh_pt->nregion_element(r);
   for (unsigned e=0;e<nel;e++)
    {
     ELEMENT *el_pt=dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(r,e));
     
     // Non-dim thermal inertia
     el_pt->alpha_pt()=&ProblemParameters::Alpha_rock;
     
     // Non-dim thermal inertia
     el_pt->beta_pt()=&ProblemParameters::Beta_rock;
    }


   // Set (pseudo-)solid mechanics properties for all elements
   //---------------------------------------------------------
   unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a solid element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     
     // Set the constitutive law
     el_pt->constitutive_law_pt() =
      ProblemParameters::Constitutive_law_pt;
     
     // Set density to zero
     el_pt->lambda_sq_pt()=&ProblemParameters::Lambda_sq;
     
     // Disable inertia
     el_pt->disable_inertia();

     // Set source function
     //el_pt->source_fct_pt() = 
     // &ExactSolution::get_source_for_unsteady_heat_validation;
    }

   // Apply boundary conditions for solid
   //------------------------------------

   // Bottom: completely pinned
   unsigned b=Bottom_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
    }

   // Sides: Symmetry bcs
   b=Left_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }
   b=Right_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }

   // Assign the Lagrangian coordinates -- sensible
   // because we've completely rebuilt the mesh 
   // and haven't copied across any Lagrange multipliers
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
   
   // Loop over the contact elements, pass pointer to penetrator and make sticky
   //---------------------------------------------------------------------------
   n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     NonlinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);

     // Make it sticky to enforce permanent contact
     el_pt->enable_stick();
    }


   // Loop over the flux elements to pass pointer to prescribed flux function
   //------------------------------------------------------------------------
   n_element=Surface_flux_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement
     StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT> *el_pt = 
      dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
       Surface_flux_mesh_pt->element_pt(e));
     
       // Set atmospheric radiation fct
       el_pt->atmospheric_radiation_fct_pt()=
        &ProblemParameters::atmospheric_radiation;
    }


   // Setup shielding nodes
   setup_shielding_nodes();

   // hierher
   update_limiting_angles();

  }
 


#ifdef ADAPTIVE

 /// Pointer to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to bulk mesh
 SolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_contact_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_flux_mesh_pt;

 /// ID of left melt boundary
 unsigned Left_melt_boundary_id;

 /// ID of right melt boundary
 unsigned Right_melt_boundary_id;

 /// ID of lower contact boundary
 unsigned Lower_contact_boundary_id;

 /// ID of upper contact boundary
 unsigned Upper_contact_boundary_id;

 /// ID of bottom boundary
 unsigned Bottom_boundary_id;

 /// ID of left boundary
 unsigned Left_boundary_id;

 /// ID of right boundary
 unsigned Right_boundary_id;

 /// Trace file
 ofstream Trace_file;

 // Setup labels for output
 DocInfo Doc_info;

 /// Storage for ordered shielding nodes
 Vector<Node*> Ordered_shielding_node_pt;

 /// Face elements on potentially sun-exposed boundaries 
 Vector<FiniteElement*> Ordered_shielding_face_element_pt;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for contact problem in square domain
//========================================================================
template<class ELEMENT>
SolarRadiationProblem<ELEMENT>::SolarRadiationProblem()
{ 

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0;

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Allow for crap initial guess
 Problem::Max_residuals=10000.0;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon
  
 // The boundary, represented by polylines
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt;

 // Store coordinates of left and right point at which penetrator
 // penetrates flat surface
 Vector<double> right_penetrator(2);
 Vector<double> left_penetrator(2);

 // Vertex coordinates on boundary
 Vector<Vector<double> > bound_coords(2);
 
 // Left boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=1.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=0.0;
 bound_coords[1][1]=0.0;

 // Build the boundary polyline
 Left_boundary_id=0;
 boundary_polyline_pt.push_back(new TriangleMeshPolyLine(bound_coords,
                                                         Left_boundary_id));

 // Bottom boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=0.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=1.0;
 bound_coords[1][1]=0.0;

 // Build the boundary polyline
 Bottom_boundary_id=1;
 boundary_polyline_pt.push_back(new TriangleMeshPolyLine(bound_coords,
                                                         Bottom_boundary_id));
 

 // Right boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=1.0;
 bound_coords[0][1]=0.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=1.0;
 bound_coords[1][1]=1.0;
 
 // Build the boundary polyline
 Right_boundary_id=2;
 boundary_polyline_pt.push_back(new TriangleMeshPolyLine(bound_coords,
                                                         Right_boundary_id));


 // Half opening angle of "boulder"
 double half_opening_angle=acos((ProblemParameters::Centre[1]-1.0)/
                                ProblemParameters::Radius);
 

 // Points at which the penetrator meets the flat surface
 right_penetrator[0]=ProblemParameters::Centre[0]+
  ProblemParameters::Radius*sin(half_opening_angle);
 right_penetrator[1]=1.0;

 left_penetrator[0]=ProblemParameters::Centre[0]-
  ProblemParameters::Radius*sin(half_opening_angle);
 left_penetrator[1]=1.0;

 // Right melt boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=1.0;
 bound_coords[0][1]=1.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=right_penetrator[0];
 bound_coords[1][1]=right_penetrator[1];
 
 
 // Build boundary poly line
 Right_melt_boundary_id=3;
 boundary_polyline_pt.push_back
  (new TriangleMeshPolyLine(bound_coords,Right_melt_boundary_id));
 
 

 // Vertex coordinates on lower contact boundary
 unsigned nvertex_contact=20;
 Vector<Vector<double> > circle_bound_coords(nvertex_contact);
 
 // First point
 circle_bound_coords[0].resize(2);
 circle_bound_coords[0][0]=right_penetrator[0];
 circle_bound_coords[0][1]=right_penetrator[1];
 
 // Internal points
 for (unsigned j=1;j<nvertex_contact-1;j++)
  {
   double phi=1.5*MathematicalConstants::Pi+half_opening_angle
    -2.0*half_opening_angle*double(j)/double(nvertex_contact-1);
   
   circle_bound_coords[j].resize(2);
   circle_bound_coords[j][0]=ProblemParameters::Centre[0]+
    ProblemParameters::Radius*cos(phi);
   
   circle_bound_coords[j][1]=ProblemParameters::Centre[1]+
    ProblemParameters::Radius*sin(phi);
  }
 
 // Last point
 circle_bound_coords[nvertex_contact-1].resize(2);
 circle_bound_coords[nvertex_contact-1][0]=left_penetrator[0];
 circle_bound_coords[nvertex_contact-1][1]=left_penetrator[1];


 // Build boundary poly line
 Lower_contact_boundary_id=6;
 TriangleMeshPolyLine* lower_contact_boundary_pt=
  new TriangleMeshPolyLine(circle_bound_coords,Lower_contact_boundary_id);
 
 // Store as internal poly line 
 Vector<TriangleMeshCurveSection*> internal_polyline_pt(1);
 internal_polyline_pt[0]=lower_contact_boundary_pt;

 // Vertex coordinates on upper contact boundary
 unsigned nvertex_contact_upper=unsigned(double(nvertex_contact)*
                                         0.5*MathematicalConstants::Pi/
                                         half_opening_angle);
 circle_bound_coords.clear();
 circle_bound_coords.resize(nvertex_contact_upper);
 
 // First point: recycle
 circle_bound_coords[0].resize(2);
 circle_bound_coords[0][0]=right_penetrator[0];
 circle_bound_coords[0][1]=right_penetrator[1];
 for (unsigned j=1;j<nvertex_contact_upper-1;j++)
  {
   double phi=-(0.5*MathematicalConstants::Pi-half_opening_angle)
    +(2.0*MathematicalConstants::Pi-2.0*half_opening_angle)*
    double(j)/double(nvertex_contact_upper-1);
      
   circle_bound_coords[j].resize(2);
   circle_bound_coords[j][0]=ProblemParameters::Centre[0]+
    ProblemParameters::Radius*cos(phi);
   
   circle_bound_coords[j][1]=ProblemParameters::Centre[1]+
    ProblemParameters::Radius*sin(phi);
  }
 
 // Last point
 circle_bound_coords[nvertex_contact_upper-1].resize(2);
 circle_bound_coords[nvertex_contact_upper-1][0]=left_penetrator[0];
 circle_bound_coords[nvertex_contact_upper-1][1]=left_penetrator[1];
 
 // Build boundary poly line
 Upper_contact_boundary_id=4; 
 TriangleMeshPolyLine* upper_contact_boundary_pt=
  new TriangleMeshPolyLine(circle_bound_coords,Upper_contact_boundary_id);
 boundary_polyline_pt.push_back(upper_contact_boundary_pt);
 
 // Left melt boundary
 Vector<Vector<double> > left_melt_bound_coords;
 Vector<double> tmp_vector(2);
 tmp_vector[0]=left_penetrator[0];
 tmp_vector[1]=left_penetrator[1];
 left_melt_bound_coords.push_back(tmp_vector);


 // Funky bit
 {
  tmp_vector[0]=0.2;
  tmp_vector[1]=1.0;
  left_melt_bound_coords.push_back(tmp_vector);
  
  tmp_vector[0]=0.2;
  tmp_vector[1]=2.2;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.5;
  tmp_vector[1]=2.2;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.5;
  tmp_vector[1]=1.7;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.4;
  tmp_vector[1]=1.7;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.4;
  tmp_vector[1]=2.0;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.35;
  tmp_vector[1]=2.0;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.35;
  tmp_vector[1]=1.65; 
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.55;
  tmp_vector[1]=1.65;
  left_melt_bound_coords.push_back(tmp_vector);


  tmp_vector[0]=0.55;
  tmp_vector[1]=2.25;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.15;
  tmp_vector[1]=2.25;
  left_melt_bound_coords.push_back(tmp_vector);

  tmp_vector[0]=0.15;
  tmp_vector[1]=1.0;
  left_melt_bound_coords.push_back(tmp_vector);
 }


 tmp_vector[0]=0.0;
 tmp_vector[1]=1.0;
 left_melt_bound_coords.push_back(tmp_vector);
  
 // Build boundary poly line
 Left_melt_boundary_id=5;
 TriangleMeshPolyLine* left_melt_boundary_pt=
  new TriangleMeshPolyLine(left_melt_bound_coords,Left_melt_boundary_id);
 boundary_polyline_pt.push_back(left_melt_boundary_pt);
 
 // Connect first vertex (vertex 0) of upper contact boundary
 // with initial vertex of lower contact boundary
 unsigned vertex_id_of_connection=0;
 lower_contact_boundary_pt->connect_initial_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>(upper_contact_boundary_pt),
  vertex_id_of_connection);


//#define USE_BROKEN
#ifdef USE_BROKEN

 // hierher this misses out final node on lower contact boundary
 // in boundary lookup scheme!

 // Connect last vertex  of upper contact boundary
 // with final vertex of lower contact boundary
 vertex_id_of_connection=nvertex_contact_upper-1;
 lower_contact_boundary_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>(upper_contact_boundary_pt),
  vertex_id_of_connection);

#else 

 // Connect first vertex left melt boundary
 // with final vertex of lower contact boundary
 vertex_id_of_connection=0;
 lower_contact_boundary_pt->connect_final_vertex_to_polyline(
  dynamic_cast<TriangleMeshPolyLine*>(left_melt_boundary_pt),
  vertex_id_of_connection);

#endif
 
 // Create open curve that defines boulder/ice interface
 Vector<TriangleMeshOpenCurve*> inner_boundary_pt;
 inner_boundary_pt.push_back(new TriangleMeshOpenCurve(internal_polyline_pt));
 

 // Create the triangle mesh polygon for outer boundary
 //----------------------------------------------------
 TriangleMeshPolygon *outer_polygon =
  new TriangleMeshPolygon(boundary_polyline_pt);
  
 // Set the pointer
 closed_curve_pt = outer_polygon;
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the maximum area element
 double uniform_element_area=0.002; 
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Specify the internal open boundary
 triangle_mesh_parameters.internal_open_curves_pt()=inner_boundary_pt;

 /// Identify penetrator as region 1
 triangle_mesh_parameters.add_region_coordinates(1,ProblemParameters::Centre);
 

#ifdef ADAPTIVE

 // Create the mesh
 Bulk_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                       time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // // Set element size limits
 // Bulk_mesh_pt->max_element_size()=0.2;
 // Bulk_mesh_pt->min_element_size()=0.002; 

#else

 // Build mesh
 Bulk_mesh_pt=new SolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                             time_stepper_pt());

#endif
 
 Bulk_mesh_pt->output("mesh.dat");
 Bulk_mesh_pt->output_boundaries("boundaries.dat");
 
 // Create the surface mesh as an empty mesh
 Surface_contact_mesh_pt=new Mesh;
 
 // Build 'em 
 create_contact_elements();
 
 // Create the surface mesh as an empty mesh
 Surface_flux_mesh_pt=new Mesh;
 
 // Build 'em 
 create_flux_elements();

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_contact_mesh_pt);
 add_sub_mesh(Surface_flux_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // Set the initial conditions
 unsigned nnod = Bulk_mesh_pt->nnode();
 for(unsigned j=0;j<nnod;j++)
  {
   Bulk_mesh_pt->node_pt(j)->set_value(0,ExactSolution::U0); 
  } 

 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void SolarRadiationProblem<ELEMENT>::doc_solution()
{ 

 oomph_info << "outputting step: " << Doc_info.number() << std::endl;

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();
 
 // Output contact elements
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
  }
 some_file.close();
 
 // Output exact solution
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(
  some_file,npts,time_pt()->time(),
  ExactSolution::get_exact_u_for_unsteady_heat_validation); 
 some_file.close();

 // Output Solar radiation
  double solar_flux_magnitude=0.0; 
 Vector<double> solar_flux_unit_vector(2); 
 double total_diffuse_radiation=0.0;
 ProblemParameters::atmospheric_radiation(time_pt()->time(),
                                          solar_flux_magnitude, 
                                          solar_flux_unit_vector, 
                                          total_diffuse_radiation);
 sprintf(filename,"%s/solar_radiation%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 some_file << "0.0 3.5 " 
           << solar_flux_magnitude  << " "
           << solar_flux_unit_vector[0] << " " 
           << solar_flux_unit_vector[1] << " "
           << "\n"; 
 some_file.close();
 

 // Output atmospheric radiation along all exposed surfaces
 sprintf(filename,"%s/atmospheric_radiation%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nel=Surface_flux_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
    Surface_flux_mesh_pt->
    element_pt(e))->output_atmospheric_radiation(some_file);
  }
 some_file.close();


 // Output illumination angles for all integration points
 sprintf(filename,"%s/illumination_angles%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nel=Surface_flux_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
    Surface_flux_mesh_pt->
    element_pt(e))->output_limiting_angles(some_file);
  }
 some_file.close();
 
 bool plot_individual=true;
 if (plot_individual)
  {
   
   // Output cone of diffuse radiation for all integration points
   // (suitable for paraview animation -- each element is in a separate
   // file)
   double radius=2.0;
   nel=Ordered_shielding_face_element_pt.size();
   unsigned count=0;
   for (unsigned e=0;e<nel;e++)
    {
     sprintf(filename,"%s/diffuse_radiation_cone%i_%i.dat",
             Doc_info.directory().c_str(),
             Doc_info.number(),count++);
     some_file.open(filename);
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
      Ordered_shielding_face_element_pt[e])
      ->output_diffuse_radiation_cone(some_file,radius);
     some_file.close();
    }
   count=0;
   for (unsigned e=0;e<nel;e++)
    {
     sprintf(filename,"%s/diffuse_radiation_cone_max_angle%i_%i.dat",
             Doc_info.directory().c_str(),
             Doc_info.number(),count++);
     some_file.open(filename);
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
      Ordered_shielding_face_element_pt[e])
      ->output_diffuse_radiation_cone_max_angle(some_file,radius);
     some_file.close();
    }
   count=0;
   for (unsigned e=0;e<nel;e++)
    {
     sprintf(filename,"%s/diffuse_radiation_cone_min_angle%i_%i.dat",
             Doc_info.directory().c_str(),
             Doc_info.number(),count++);
     some_file.open(filename);
     dynamic_cast<StefanBoltzmannUnsteadyHeatFluxElement<ELEMENT>*>(
      Ordered_shielding_face_element_pt[e])
      ->output_diffuse_radiation_cone_min_angle(some_file,radius);
     some_file.close();
    }
  }


 // Output penetrator
 sprintf(filename,"%s/penetrator%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n=100;
 ProblemParameters::Penetrator_pt->output(some_file,n);
 some_file.close();
 
 // Output Number of Newton iterations in form that can be visualised
 // as vector in paraview
 sprintf(filename,"%s/newton_iter%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 some_file << "0 0 0 " << Nnewton_iter_taken << std::endl;
 some_file.close();
 
 // Write norm of solution to trace file
 double norm=0.0;
 //Bulk_mesh_pt->compute_norm(norm); 
 Trace_file  << norm << std::endl;
 
 //Increment counter for solutions 
 Doc_info.number()++;

} // end of doc_solution



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// Driver code
//========================================================================
int main(int argc, char* argv[])
{

 // hierher
 FiniteElement::Accept_negative_jacobian=true;
 FiniteElement::Tolerance_for_singular_jacobian=0.0;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--no_adapt");
  
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 ProblemParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&ProblemParameters::Nu);
  

 // Define centre of penetrator
 ProblemParameters::Centre.resize(2);
 ProblemParameters::Centre[0]=0.5;
 ProblemParameters::Centre[1]=ProblemParameters::Y_c_initial; 

 // Create penetrator
 ProblemParameters::Penetrator_pt = new
  CircularPenetrator(&ProblemParameters::Centre,
                     ProblemParameters::Radius);

#ifdef ADAPTIVE
  
 // Build problem
 SolarRadiationProblem<ProjectableUnsteadyHeatElement<
  PseudoSolidNodeUpdateElement<TUnsteadyHeatElement<2,3>,TPVDElement<2,3> > > >
                    problem;
                    
 #else

 // Build problem
 SolarRadiationProblem<TPVDElement<2,3> > problem;


#endif

  //Output initial condition
 problem.doc_solution();
 
 unsigned max_adapt=1; 
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }

 
 // Number of parameter increments per period
 unsigned nstep_for_period=100; 

 // Parameter variation
 unsigned nperiod=3;

 // Initial timestep
 double dt=1.0/double(nstep_for_period);

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);

 // Set impulsive IC
 problem.assign_initial_values_impulsive(dt);

 while (problem.time_pt()->time()<double(nperiod))
  {
#ifdef ADAPTIVE

   // Dummy double adaptivity (timestep is always accepted because
   // tolerance is set to huge value; mainly used to automatically
   // re-solve with smaller timestep increment after non-convergence
   double epsilon_t=DBL_MAX;
   bool first=false;
   double next_dt=
    problem.doubly_adaptive_unsteady_newton_solve(dt,epsilon_t,
                                                  max_adapt,first);
   dt = next_dt; 
   
#else

   // hierher 
   abort();

   // Solve
   problem.newton_solve();
   
#endif
   
   
   //Output solution
   problem.doc_solution();

  }
 
} // end of main

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


//Generic libraries
#include "generic.h"
#include "solid.h"
#include "constitutive.h"
#include "navier_stokes.h"
#include "multi_physics.h"

// Get the mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/quarter_tube_mesh.h"

// include oomph and std
using namespace std;
using namespace oomph;



namespace oomph {

//============================================================================
/// Pseudo-Elastic Solid element class to overload the Block Preconditioner
/// methods ndof_types() and get_dof_numbers_for_unknowns() to differentiate
/// between DOFs subject to Lagrange multiplier and those that are not.
//============================================================================
template <class ELEMENT>
class PseudoElasticBulkElement : 
 public virtual ELEMENT
{

public:

 /// default constructor
 PseudoElasticBulkElement() : ELEMENT() {}
 
 /// returns the number of DOF types associated with this element. 
 unsigned ndof_types() const
  {
   return 2*ELEMENT::dim();
  }
 
 /// Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "DOF" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)\n
 /// E.g. in a 3D problem there are 6 types of DOF:\n
 /// 0 - x displacement (without lagr mult traction)\n
 /// 1 - y displacement (without lagr mult traction)\n
 /// 2 - z displacement (without lagr mult traction)\n
 /// 4 - x displacement (with lagr mult traction)\n
 /// 5 - y displacement (with lagr mult traction)\n
 /// 6 - z displacement (with lagr mult traction)\n
 void get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // temporary pair (used to store dof lookup prior to being added to list
   std::pair<unsigned,unsigned> dof_lookup;
   
   // number of nodes
   const unsigned n_node = this->nnode();
   
   //Get the number of position dofs and dimensions at the node
   const unsigned n_position_type = ELEMENT::nnodal_position_type();
   const unsigned nodal_dim = ELEMENT::nodal_dimension();
   
   //Integer storage for local unknown
   int local_unknown=0;
   
   //Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     unsigned offset = 0;
     if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
      {
       offset = ELEMENT::dim();
      }
     
     //Loop over position dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over dimension
       for(unsigned i=0;i<nodal_dim;i++)
        {
         //If the variable is free
         local_unknown = ELEMENT::position_local_eqn(n,k,i);
         
         // ignore pinned values
         if (local_unknown >= 0)
          {
           // store dof lookup in temporary pair: First entry in pair
           // is global equation number; second entry is dof type
           dof_lookup.first = this->eqn_number(local_unknown);
           dof_lookup.second = offset+i;
           
           // add to list
           dof_lookup_list.push_front(dof_lookup);
           
          }
        }
      }
    }
  }
};





//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<PseudoElasticBulkElement<ELEMENT> > :
 public virtual FaceGeometry<ELEMENT>
{
};



}

#ifdef OOMPH_HAS_HYPRE

//=========================================================================
/// Navier Stokes LSC Preconditioner Subsidiary Helper
//=========================================================================
namespace LSC_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  HyprePreconditioner* hypre_preconditioner_pt
   = new HyprePreconditioner;
  hypre_preconditioner_pt->set_amg_iterations(2);
  hypre_preconditioner_pt->amg_using_simple_smoothing();
  hypre_preconditioner_pt->amg_simple_smoother() = 0;
  hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
  hypre_preconditioner_pt->amg_strength() = 0.25;
  hypre_preconditioner_pt->amg_coarsening() = 6;
  hypre_preconditioner_pt->amg_damping() = 0.5;
  hypre_preconditioner_pt->disable_doc_time();
  return hypre_preconditioner_pt;
 }
}



//=========================================================================
/// Real Solid Preconditioner Subsidiary Helper
//=========================================================================
namespace Real_Solid_Preconditioner_Helper
{
 Preconditioner* get_preconditioner()
 {
  HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;
  hypre_preconditioner_pt->set_amg_iterations(2);
  hypre_preconditioner_pt->amg_using_simple_smoothing();
  hypre_preconditioner_pt->amg_simple_smoother() = 3;
  hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
  hypre_preconditioner_pt->amg_strength() = 0.25;
  hypre_preconditioner_pt->amg_coarsening() = 6;
  hypre_preconditioner_pt->disable_doc_time();
  return hypre_preconditioner_pt;
 }
}

#endif

/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////




//=========================================================================
// the wall mesh
//=========================================================================
template<class ELEMENT>
class WallMesh : public virtual RefineableSimpleCubicMesh<ELEMENT>, 
                 public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 WallMesh(const double& r_min, const double& r_max, const double& lz,
          const int& nr, const int& nphi,const int& nz,
          TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper)
  : SimpleCubicMesh<ELEMENT>(nr,nphi,nz,1.0,1.0,lz,time_stepper_pt),
    RefineableSimpleCubicMesh<ELEMENT>(nr,nphi,nz,1.0,1.0,lz,time_stepper_pt)
  {

   // Important: Currently wall mesh will not refine
   // towards cylindrical geometry! Need Domain!

   // move the nodes

   // Find out how many nodes there are
   unsigned n_node=nnode();

   // Calculate the value of pi
   const double pi = MathematicalConstants::Pi;

   // Loop over all nodes
   for (unsigned n=0;n<n_node;n++)
    {
     // Pointer to node
     Node* nod_pt=node_pt(n);

     // Get the x/y coordinates
     double x_old=nod_pt->x(0);
     double y_old=nod_pt->x(1);

     // Map from the old x/y to the new r/phi:
     double r=r_min+(r_max-r_min)*x_old;
     double phi=(pi/2.0)*(y_old);

     // Set new nodal coordinates
     nod_pt->x(0)=r*cos(phi);
     nod_pt->x(1)=r*sin(phi);

     // set boundary coordinates
     if (nod_pt->is_on_boundary(4))
      {
       Vector<double> zeta(2);
       zeta[0]=nod_pt->x(2);
       zeta[1]=phi;
       nod_pt->set_coordinates_on_boundary(4,zeta);
      }
    }
   this->Boundary_coordinate_exists[4]=true;
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~WallMesh() { }


 /// return the minimum x-position of a fluid-wall interface node
 /// on the x / z face 
 double min_deformation()
  {
   double xmin = 1.0;
   unsigned n_node=nnode();
   for (unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt=node_pt(n);
     if (nod_pt->is_on_boundary(1))
      {
       xmin=min(xmin,nod_pt->x(0));
      }
    }
   return xmin;
  }

 /// return the maximum y-position of a fluid-wall interface node
 /// on the y / z face 
 double max_deformation()
  {
   double ymax = 1.0;
   unsigned n_node=nnode();
   for (unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt=node_pt(n);
     if (nod_pt->is_on_boundary(3)&&nod_pt->is_on_boundary(4))
      {
       ymax=max(ymax,nod_pt->x(1));
      }
    }
   return ymax;
  }
};


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//==============start_fluid_mesh===========================================
// the fluid mesh
//=========================================================================
template<class ELEMENT>
class FluidMesh : public virtual RefineableQuarterTubeMesh<ELEMENT>,
                  public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 FluidMesh(GeomObject* wall_pt,
           const Vector<double>& xi_lo,
           const double& fract_mid,
           const Vector<double>& xi_hi,
           const unsigned& nlayer,
           TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper)
  : QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,nlayer,
                             time_stepper_pt),
    RefineableQuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,nlayer,
                                       time_stepper_pt)
  {
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~FluidMesh() { }
};
 

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////


//===================================================================
/// Global variables
//===================================================================
namespace Global_Parameters
{
 /// tube length
 double L=12.0;

 /// upstream length
 double Lup=1.5;
 
 /// downstream length
 double Ldown=3.0;

 /// wall thickness
 double H=0.15;

 /// number of axial slices in fluid mesh
 unsigned N_slice=4;

 /// reynolds number
 double Re=100.0; 

 /// Default FSI parameter
 double Q=1e-5;

 /// Pointer to constitutive law for the wall
 ConstitutiveLaw* Constitutive_law_wall_pt=0;

 /// Pointer to constitutive law for the pseudo elastic node update elements
 ConstitutiveLaw* Constitutive_law_pseudo_elastic_pt=0;

 /// Poisson's ratio for generalised Hookean constitutive equation for the
 /// wall
 double Nu_wall=0.3;

 /// Poisson's ratio for generalised Hookean constitutive equation for the
 /// pseudo elastic bulk
 double Nu_pseudo_elastic=0.1;

 /// Uniform pressure acting on wall
 double P = 2.0e-4;

 /// Perturbation pressure
 double Pcos = 2.0e-4;

//  /// external pressure increment
//  double Pincrement=5.0e-4;

//  /// maximum external pressure
//  double Pmax=1.4;

 /// Load function for wall
 void press_load(const Vector<double> &xi,
                 const Vector<double> &x,
                 const Vector<double> &N,
                 Vector<double>& load)
 {  
  // Get polar angle
  double phi=atan2(xi[1],xi[0]);

//  if (xi[0]==0.0) assert(false);

  // Pressure with pcos perturbation
  for(unsigned i=0;i<3;i++)
   {
    load[i] = (-P - Pcos*cos(2.0*phi))*N[i];
   }
 } 
 
 /// Timescale ratio (non-dim density) for solid
 double Lambda_sq=0.0;

 /// Period of periodic variation in inflow pressure
 double T=1.0;

 /// Timestep
 double Dt=0.01;

 /// Number of steps per period
 unsigned Nstep_per_period=40;

 /// Number of periods
 unsigned Nperiod=5;

} //end_of_namespace


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////


//===============start_of_problem_class===============================
// pseudo elastic collapsible channel fsi problem
//====================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
class PseudoElasticCollapsibleChannelProblem : public Problem
{

public:

 /// Constructor: 
 PseudoElasticCollapsibleChannelProblem();

 /// Destructor (empty)
 ~PseudoElasticCollapsibleChannelProblem()
  {
   delete Fluid_time_stepper_pt;
   delete Wall_pt;
   delete Wall_time_stepper_pt;
   delete Fluid_mesh_pt->spatial_error_estimator_pt();
   delete Fluid_mesh_pt;
   delete Solid_mesh_pt->spatial_error_estimator_pt();
   delete Solid_mesh_pt;
   delete Solid_fsi_traction_mesh_pt;
   delete Lagrange_multiplier_mesh_pt;
   delete Solid_traction_mesh_pt;
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Do steady run 
 void steady_run(DocInfo& doc_info);

 /// Do unsteady run 
 void unsteady_run(DocInfo& doc_info);

 /// Create FSI traction elements
 void create_fsi_traction_elements();

 /// Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 /// create solid traction elements
 /// (to apply pressure to external solid wall)
 void create_solid_traction_elements();

 /// Actions before timestep 
 void actions_before_implicit_timestep();

 /// Update no slip bc after update of unknowns
 void actions_before_newton_convergence_check()
  {
   Fluid_mesh_pt->node_update();
  }


 /// Actions before Newton solve: Reset pseudo-elastic Lagrangian coordinates
 void actions_before_newton_solve()
  {
   Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  }


 /// Before adapt: Flush face element submeshes
 void actions_before_adapt();
 
 /// After adapt: Rebuild face element submeshes and re-setup FSI
 void actions_after_adapt();

 /// Before distribute: Flush face submeshes and keep all
 /// solid elements adjacent to the fluid as halos.
 void actions_before_distribute();
 
 /// After distribute: Rebuild face element submeshes and re-setup
 /// FSI
 void actions_after_distribute();

 /// Helper function to delete elements from a mesh
 void empty_mesh(Mesh* const &surface_mesh_pt);


/// Switch to GMRES preconditioned with the pseudo-elastic 
/// FSI preconditioner
 void set_pseudo_elastic_fsi_solver();

 /// Doc parameters
 void doc_parameters()
  {
   std::cout << "\n\n=================================================\n";
   std::cout << "Q                     : " << Global_Parameters::Q 
             << std::endl;
   std::cout << "Re                    : " << Global_Parameters::Re 
             << std::endl;
   std::cout << "P                     : " << Global_Parameters::P 
             << std::endl;
   std::cout << "Pcos                   : " << Global_Parameters::Pcos 
             << std::endl;
   std::cout << "T                     : " << Global_Parameters::T 
             <<std::endl;
   std::cout << "t                     : " << time_pt()->time()
             << std::endl;
   std::cout << "=================================================\n\n";
  }
 
  /// Bulk solid mesh
 WallMesh<SOLID_ELEMENT>* wall_mesh_pt() {return Solid_mesh_pt; }

 /// Bulk fluid mesh
 FluidMesh<FLUID_ELEMENT>* fluid_mesh_pt() {return Fluid_mesh_pt; }

 /// Meshes of Lagrange multiplier elements
 SolidMesh*& lagrange_multiplier_mesh_pt() 
  {return Lagrange_multiplier_mesh_pt; }

 /// Set initial conditions 
 void set_initial_condition();

private:
 
 /// Sanity check: Doc boundary coordinates on i-th solid FSI interface
 void doc_solid_boundary_coordinates();
 
 /// Bulk solid mesh
 WallMesh<SOLID_ELEMENT>* Solid_mesh_pt;

 /// Meshes of FSI traction elements (fluid traction onto solid)
 SolidMesh* Solid_fsi_traction_mesh_pt;

 /// Bulk fluid mesh
 FluidMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

 /// Mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// Solid traction elements (prescibed external pressure on solid wall)
 SolidMesh* Solid_traction_mesh_pt;

 /// Geometric Object defining the undeformed boundary of the
 /// fluid mesh
 GeomObject* Wall_pt;
 
 /// Timestepper for the Navier-Stokes equations
 BDF<2>* Fluid_time_stepper_pt;
 
 /// Steady timestepper that stores two history values
 Steady<2>* Wall_time_stepper_pt; 
 
};



//==========start_of_constructor==========================================
/// Constructor for unstructured 3D FSI problem
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
PseudoElasticCollapsibleChannelProblem()
{
 // DEBUG: doc full stats for the multi-domain interaction
 //Multi_domain_functions::Doc_full_stats=true;
 //Multi_domain_functions::Doc_stats=true;

 // Allocate the timestepper for the Navier-Stokes equations
 Fluid_time_stepper_pt=new BDF<2>;

 // Add the fluid timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(Fluid_time_stepper_pt);
 
 // Create a Steady timestepper that stores two history values
 Wall_time_stepper_pt = new Steady<2>; 

 // Add the wall timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(Wall_time_stepper_pt);

 // Domain size and shape parameters
 //---------------------------------

 // Create geometric object that defines curvilinear boundary of
 // fluid domain: Elliptical tube with half axes = radius = 1.0
 Wall_pt=new EllipticalTube(1.0,1.0);
 
 // Bounding Lagrangian coordinates
 Vector<double> xi_lo(2);
 xi_lo[0]=0.0;
 xi_lo[1]=0.0;
 Vector<double> xi_hi(2);
 xi_hi[0]= Global_Parameters::L;
 xi_hi[1]= 0.5*MathematicalConstants::Pi;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;
 
 // Define fluid mesh
 //------------------
 
 // create the fluid mesh
 Fluid_mesh_pt =  new FluidMesh<FLUID_ELEMENT>
  (Wall_pt,xi_lo,frac_mid,xi_hi,Global_Parameters::N_slice,
   Fluid_time_stepper_pt);
  

 // Setup error estimator
 if (!CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   Fluid_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
   Fluid_mesh_pt->max_permitted_error()=0.004;
   Fluid_mesh_pt->min_permitted_error()=0.0001; 
  }
 // Validation: Fake
 else
  {
   // Setup fake error estimator:  Refine one element on the outer
   // and one on the inner surface of the wall
   Vector<unsigned> elements_to_refine(2);
   elements_to_refine[0]=3;
   elements_to_refine[1]=8;

   unsigned central_node_number=13;
   bool use_lagrangian_coordinates=false;
   Fluid_mesh_pt->spatial_error_estimator_pt()=
    new DummyErrorEstimator(Fluid_mesh_pt,elements_to_refine,
                            central_node_number,
                            use_lagrangian_coordinates);
  }


 // Define solid mesh
 //------------------
 
 //Create solid bulk mesh with a slightly different number of elements
 //in the axial direction
 double r_min=1.0;
 double r_max=1.0+Global_Parameters::H;
 double& lz=Global_Parameters::L;
 unsigned nr=2;
 unsigned nphi=2; 
 unsigned nz=Global_Parameters::N_slice;
 // Try non-matching discretisation if validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   nz=unsigned(double(nz)*1.5);
  }
 Solid_mesh_pt = new WallMesh<SOLID_ELEMENT>(r_min,r_max,lz,
                                             nr,nphi,nz,Wall_time_stepper_pt);
 

 // Setup error estimator
 if (!CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   Solid_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
   Solid_mesh_pt->max_permitted_error()=0.21;
   Solid_mesh_pt->min_permitted_error()=0.001; 
  }
 // Validation: Fake
 else
  {     
   // Setup fake error estimator: This refines two elements
   Vector<unsigned> elements_to_refine(1);
   elements_to_refine[0]=9;

   unsigned central_node_number=13;
   bool use_lagrangian_coordinates=true;
   Solid_mesh_pt->spatial_error_estimator_pt()=
    new DummyErrorEstimator(Solid_mesh_pt,elements_to_refine,
                            central_node_number,
                            use_lagrangian_coordinates);
  }


 // Create FSI Traction elements
 //-----------------------------
 
 // Build the FSI traction elements
 Solid_fsi_traction_mesh_pt=new SolidMesh;
 create_fsi_traction_elements();

 // Create Lagrange multiplier mesh for boundary motion of fluid mesh
 //------------------------------------------------------------------
 
 // Create elements
 Lagrange_multiplier_mesh_pt = new SolidMesh;
 create_lagrange_multiplier_elements();

 // Create solid traction elements for external pressure
 //-----------------------------------------------------

 Solid_traction_mesh_pt = new SolidMesh;
 create_solid_traction_elements();

 // Combine the lot
 //----------------
 
 // Add sub meshes:

 // The fluid bulk mesh
 add_sub_mesh(Fluid_mesh_pt);

 // Solid bulk mesh
 add_sub_mesh(Solid_mesh_pt);
  
 // The solid fsi traction mesh
 add_sub_mesh(Solid_fsi_traction_mesh_pt);

 // the solid traction elements
 add_sub_mesh(Solid_traction_mesh_pt);

 // The Lagrange multiplier mesh for the fluid
 add_sub_mesh(Lagrange_multiplier_mesh_pt);

 // Build global mesh
 build_global_mesh();
 
 // Apply BCs for fluid/pseudo-elastic elements
 //--------------------------------------------
  
 unsigned nbound_fluid = Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<nbound_fluid;b++)
  {
   unsigned nnode = Fluid_mesh_pt->nboundary_node(b);
   for (unsigned i = 0; i < nnode; i++)
    {

     // Fluid velocity BCs

     // Boundary 1 is the vertical symmetry boundary (y=const) where
     // we allow flow in the y-direction. Elsewhere, pin the y velocity
     if(b!=1) Fluid_mesh_pt->boundary_node_pt(b,i)->pin(1);

     // Boundary 2 is the horizontal symmetry boundary (x=const) where
     // we allow flow in the x-direction. Elsewhere, pin the x velocity
     if(b!=2) Fluid_mesh_pt->boundary_node_pt(b,i)->pin(0);

     // Boundaries 0 and 3 are the inflow and the wall respectively.
     // Pin the axial velocity because of the prescribed inflow
     // profile and the no slip on the wall (where the velocity is imposed
     // by the node update), respectively
     if((b==0)||(b==3)) 
      {
       Fluid_mesh_pt->boundary_node_pt(b,i)->pin(2);
      }

     // Pseudo solid BCs
     SolidNode* solid_node_pt = dynamic_cast<SolidNode*>
      (Fluid_mesh_pt->boundary_node_pt(b,i));
     
     // Inflow, vertical symmetry BC (y=const), outflow: Suppress x-displ
     if((b==0)||(b==1)||(b==4)) { solid_node_pt->pin_position(0); }

     // Inflow, horizontal symmetry BC (x=const), outflow: Suppress y-displ
     if((b==0)||(b==2)||(b==4)) { solid_node_pt->pin_position(1); }

     // Inflow and outflow: Suppress z-displacement
     if((b==0)||(b==4)) { solid_node_pt->pin_position(2); }
    }
  }


 // Complete the build of the fluid elements so they are fully functional
 //----------------------------------------------------------------------

 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   FLUID_ELEMENT* el_pt = 
    dynamic_cast<FLUID_ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;
   el_pt->re_st_pt() = &Global_Parameters::Re; 
   
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pseudo_elastic_pt;

   // Set the timescale ratio 
   el_pt->lambda_sq_pt() =
    &Global_Parameters::Lambda_sq;

  } // end loop over elements

 // pin the rundundant nodal pressures
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());

 // Apply BCs for solid
 //--------------------
 
 unsigned nbound_solid = Solid_mesh_pt->nboundary();
 for(unsigned b=0;b<nbound_solid;b++)
  {
   unsigned nnode = Solid_mesh_pt->nboundary_node(b);
   for (unsigned i = 0; i < nnode; i++)
    {
     SolidNode* solid_node_pt=
      dynamic_cast<SolidNode*>(Solid_mesh_pt->boundary_node_pt(b,i));

     // Inflow, vertical symmetry (y=const) and outflow: Pin x-displ
     if((b==0)||(b==3)||(b==5)){solid_node_pt->pin_position(0);}

     // Inflow, horizontal symmetry (x=const) and outflow: Pin y-displ
     if((b==0)||(b==1)||(b==5)){solid_node_pt->pin_position(1);}

     // Inflow and outflow: Pin z-displ
     if((b==0)||(b==5)){solid_node_pt->pin_position(2);}

     // Pin all positions in "rigid sections"
     double z = solid_node_pt->x(2);
     if (z<=Global_Parameters::Lup||
         z>=(Global_Parameters::L-Global_Parameters::Ldown))
      {
       solid_node_pt->pin_position(0);
       solid_node_pt->pin_position(1);
       solid_node_pt->pin_position(2);
      }
    }
  }
 
 // Complete the build of Solid elements so they are fully functional
 //----------------------------------------------------------------

 n_element = Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   SOLID_ELEMENT *el_pt = dynamic_cast<SOLID_ELEMENT*>(
    Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law   
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_wall_pt;

   // Set the timescale ratio 
   el_pt->lambda_sq_pt() =
    &Global_Parameters::Lambda_sq;
  }

 // Setup FSI
 //----------
     
 // The velocity of the fluid nodes on the fsi boundary
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:

 // Get boundary ID
 unsigned b=3;
   
 // How many nodes on boundary b?
 unsigned nnode_fsi_fluid =  Fluid_mesh_pt->nboundary_node(b);
 
 for (unsigned i=0;i<nnode_fsi_fluid;i++)
  {
   Fluid_mesh_pt->
    boundary_node_pt(b,i)->set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }

 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. 
   
 //Doc boundary coordinates in fluid
 Multi_domain_functions::Doc_boundary_coordinate_file.open(
 	"RESLT/fluid_boundary_coordinates.dat");
 
 // Setup FSI: Pass ID of fluid FSI boundary and associated
 // mesh of solid fsi traction elements.
 FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,3>
  (this,3,Fluid_mesh_pt,Solid_fsi_traction_mesh_pt);
 
 // Close the doc file
 Multi_domain_functions::Doc_boundary_coordinate_file.close();

 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
}


//====start_of_apply_initial_condition========================================
/// Apply initial conditions
//============================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
set_initial_condition()
{ 
 // Set Poiseuille flow
 unsigned nnode = Fluid_mesh_pt->nnode();
 for (unsigned i = 0; i < nnode; i++)
  {
   double x = Fluid_mesh_pt->node_pt(i)->x(0);
   double y = Fluid_mesh_pt->node_pt(i)->x(1);
   double r = sqrt(x*x+y*y);
   double v = 2.0*(1.0-r*r);
   Fluid_mesh_pt->node_pt(i)->set_value(0,0.0);
   Fluid_mesh_pt->node_pt(i)->set_value(1,0.0);
   Fluid_mesh_pt->node_pt(i)->set_value(2,v);
  }


 // Assign initial values for an impulsive start
 Fluid_mesh_pt->assign_initial_values_impulsive();
 Solid_mesh_pt->assign_initial_values_impulsive();

}

//============start_of_create_fsi_traction_elements=========================
/// Create FSI traction elements that apply fluid traction onto wall
//==========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_fsi_traction_elements()
{
 // Boundary 4 is the is the fsi boundary of solid mesh
 unsigned b=4;
   
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Solid_mesh_pt->nboundary_element(b);
   
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
    Solid_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of the element e along boundary b
   int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
     
   // Create new element 
   RefineableFSISolidTractionElement<SOLID_ELEMENT,3>* el_pt=
    new RefineableFSISolidTractionElement<SOLID_ELEMENT,3>(bulk_elem_pt,
                                                           face_index);
   
   // Add it to the mesh
   Solid_fsi_traction_mesh_pt->add_element_pt(el_pt);
     
   // Specify boundary number
   el_pt->set_boundary_number_in_bulk_mesh(b);
   
   // Specify FSI parameter
   el_pt->q_pt() = &Global_Parameters::Q; 
  }
} // end of create_fsi_traction_elements


//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Create Lagrange multiplier elements on boundary 3 of fluid mesh
 unsigned b=3;
 
 // How many bulk fluid elements are adjacent to boundary b?
 unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk fluid elements adjacent to boundary b 
 // to create elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk fluid element that is adjacent to boundary b
   FLUID_ELEMENT* bulk_elem_pt = dynamic_cast<FLUID_ELEMENT*>(
    Fluid_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
   
   // Create new element
   RefineableFSIImposeDisplacementByLagrangeMultiplierElement<FLUID_ELEMENT>* 
    el_pt =
    new RefineableFSIImposeDisplacementByLagrangeMultiplierElement
    <FLUID_ELEMENT>(bulk_elem_pt,face_index);
   
   // Add it to the mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);
   
   // Specify boundary number
   el_pt->set_boundary_number_in_bulk_mesh(b);
  }
 
 
 // Locate bulk solid elements next to boundary b_solid_fsi (the FSI boundary) 
 // in the mesh pointed to by Solid_mesh_pt that drive the deformation of the 
 // pseudo-solid fluid mesh via the Lagrange multiplier elements
 // stored in Lagrange_multiplier_mesh_pt
 unsigned b_solid=4;
 FSI_functions::setup_solid_elements_for_displacement_bc<SOLID_ELEMENT,3>
  (this,b_solid,Solid_mesh_pt,Lagrange_multiplier_mesh_pt);
 
 // Loop over the bulk fluid elements adjacent to boundary b
 // to apply boundary conditions
 for(unsigned e=0;e<n_element;e++)
  {
   // Get element
   RefineableFSIImposeDisplacementByLagrangeMultiplierElement<FLUID_ELEMENT>* 
    el_pt=
    dynamic_cast<RefineableFSIImposeDisplacementByLagrangeMultiplierElement
    <FLUID_ELEMENT>*>(
     Lagrange_multiplier_mesh_pt->finite_element_pt(e));
   
   // Loop over the nodes 
   unsigned nnod=el_pt->nnode();
   for (unsigned i=0;i<nnod;i++)
    {
     // How many nodal values were used by the "bulk" element
     // that originally created this node?
     unsigned n_bulk_value=el_pt->nbulk_value(i);
     Node* node_pt = el_pt->node_pt(i);
     
     // The remaining ones are Lagrange multipliers and we pin them.
     unsigned nval=node_pt->nvalue();
     for (unsigned j=n_bulk_value;j<nval;j++)
      {
       // Direction of Lagrange multiplier force
       int d=j-n_bulk_value;
       
       
       // Boundary 0 (inflow): Pseudo-solid displacement is fixed so pin
       // all Lagrange multipliers
       if (node_pt->is_on_boundary(0))
        {
         if (d==0||d==1||d==2) node_pt->pin(j);
        }
       // Boundary 1 (symm BC, x=const.): Pseudo-solid displacement is 
       // fixed in x-direction so pin corresponding Lagrange multiplier
       if (node_pt->is_on_boundary(1))
        {
         if (d==0) node_pt->pin(j);
        }
       // Boundary 2 (symm BC, y=const.): Pseudo-solid displacement is 
       // fixed in y-direction so pin corresponding Lagrange multiplier
       if (node_pt->is_on_boundary(2))
        {
         if (d==1) node_pt->pin(j);
        }
       // Boundary 4 (outflow): Pseudo-solid displacement is fixed so pin
       // all Lagrange multipliers
       if (node_pt->is_on_boundary(4))
        {
         if (d==0||d==1||d==2) node_pt->pin(j);
        }
      }       
    }
  }

} // end of create_lagrange_multiplier_elements



//==================start_of_create_solid_traction_elements===============
/// Create face elements that apply external pressure load onto
/// wall
//========================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_solid_traction_elements()
{
 // Set boundary
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Solid_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
    Solid_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
      
   // Create new element
   RefineableSolidTractionElement<SOLID_ELEMENT>* el_pt = 
    new RefineableSolidTractionElement<SOLID_ELEMENT>(bulk_elem_pt,
                                                      face_index);   

   //Set the traction function
   el_pt->traction_fct_pt() = &Global_Parameters::press_load;
   
   // add to mesh
   Solid_traction_mesh_pt->add_element_pt(el_pt);
  }
}

//===========================================================================
/// Flush mesh of face elements.
//===========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_before_adapt()
{
 // Flush all the submeshes out but keep the meshes of 
 // RefineableFSISolidTractionElements alive (i.e. don't delete them)
 flush_sub_meshes();
 
 // Add the fluid mesh and the solid mesh back again
 // Remember that it's important that the fluid mesh is
 // added before the solid mesh!
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Solid_mesh_pt);

 // Rebuild the global mesh
 this->rebuild_global_mesh();
} 





//===========================================================================
/// Re-create face elements and re-setup problem
//===========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_after_adapt()
{
 // delete and recreate meshes
 this->empty_mesh(Lagrange_multiplier_mesh_pt);
 this->empty_mesh(Solid_fsi_traction_mesh_pt);
 this->empty_mesh(Solid_traction_mesh_pt);

 // Build the FSI traction elements
 create_fsi_traction_elements();

 // Create Lagrange multiplier mesh for boundary motion of fluid mesh
 create_lagrange_multiplier_elements();

 // create solid traction elements for external pressure
 create_solid_traction_elements();

 // Add these meshes as submeshes
 add_sub_mesh(Solid_fsi_traction_mesh_pt);
 add_sub_mesh(Solid_traction_mesh_pt);
 add_sub_mesh(Lagrange_multiplier_mesh_pt);

 // Rebuild the global mesh
 rebuild_global_mesh();

 // Unpin all pressure dofs
 RefineableNavierStokesEquations<3>::
  unpin_all_pressure_dofs(Fluid_mesh_pt->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());


 // Re-apply BCs for solid (just to be on the safe side, because of the
 //--------------------------------------------------------------------
 // way we're keeping the upstream and downstream sections rigid)
 //--------------------------------------------------------------
 
 unsigned nbound_solid = Solid_mesh_pt->nboundary();
 for(unsigned b=0;b<nbound_solid;b++)
  {
   unsigned nnode = Solid_mesh_pt->nboundary_node(b);
   for (unsigned i = 0; i < nnode; i++)
    {
     SolidNode* solid_node_pt=
      dynamic_cast<SolidNode*>(Solid_mesh_pt->boundary_node_pt(b,i));

     // Inflow, vertical symmetry (y=const) and outflow: Pin x-displ
     if((b==0)||(b==3)||(b==5)){solid_node_pt->pin_position(0);}

     // Inflow, horizontal symmetry (x=const) and outflow: Pin y-displ
     if((b==0)||(b==1)||(b==5)){solid_node_pt->pin_position(1);}

     // Inflow and outflow: Pin z-displ
     if((b==0)||(b==5)){solid_node_pt->pin_position(2);}

     // Pin all positions in "rigid sections"
     double z = solid_node_pt->x(2);
     if (z<=Global_Parameters::Lup||
         z>=(Global_Parameters::L-Global_Parameters::Ldown))
      {
       solid_node_pt->pin_position(0);
       solid_node_pt->pin_position(1);
       solid_node_pt->pin_position(2);
      }
    }
  }

 // Setup FSI
 //----------

 unsigned b=3;

 FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,3>
  (this,3,Fluid_mesh_pt,Solid_fsi_traction_mesh_pt);

 unsigned nnode_fsi_fluid =  Fluid_mesh_pt->nboundary_node(b);
 for (unsigned i=0;i<nnode_fsi_fluid;i++)
  {
   Fluid_mesh_pt->
    boundary_node_pt(b,i)->set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }
} 


//===========================================================================
/// Retain all bulk solid elements adjacent to fluid mesh as halos
/// then remove face-elements from problem.
//===========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_before_distribute()
{

 // Flush all the submeshes out but keep the meshes of 
 // RefineableFSISolidTractionElements alive (i.e. don't delete them)
 flush_sub_meshes();
 
 // Add the fluid mesh and the solid mesh back again
 // Remember that it's important that the fluid mesh is
 // added before the solid mesh!
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Solid_mesh_pt);
 
 // Rebuild global mesh
 rebuild_global_mesh();
} 

//===========================================================================
/// Recreate face elements and rebuild problem
//===========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_after_distribute()
{
 // delete and recreate meshes
 this->empty_mesh(Lagrange_multiplier_mesh_pt);
 this->empty_mesh(Solid_fsi_traction_mesh_pt);
 this->empty_mesh(Solid_traction_mesh_pt);
 create_fsi_traction_elements();
 create_lagrange_multiplier_elements();
 create_solid_traction_elements();
 add_sub_mesh(Solid_fsi_traction_mesh_pt);
 add_sub_mesh(Solid_traction_mesh_pt);
 add_sub_mesh(Lagrange_multiplier_mesh_pt);

 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<3>::
  unpin_all_pressure_dofs(Fluid_mesh_pt->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());

 // Setup FSI
 //----------

 unsigned b=3;

 FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,3>
  (this,3,Fluid_mesh_pt,Solid_fsi_traction_mesh_pt);

 unsigned nnode_fsi_fluid =  Fluid_mesh_pt->nboundary_node(b);
 for (unsigned i=0;i<nnode_fsi_fluid;i++)
  {
   Fluid_mesh_pt->
    boundary_node_pt(b,i)->set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }
} 


//============start_of_empty_mesh========================================
/// Delete elements in specified mesh of face elements
//=======================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::empty_mesh(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of empty_mesh




//=====================start_of_actions_before_implicit_timestep==========
/// Actions before implicit timestep
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>
::actions_before_implicit_timestep()
{
 // Get current time
 double t=time_pt()->time();

 // Factor for volume flux
 double fact=1.0+0.5*sin(MathematicalConstants::Pi*t/Global_Parameters::T);
 oomph_info << "Time, volume flux factor: " << t << " " << fact << std::endl;

 // Set inflow (bc 0)
 unsigned nnode= fluid_mesh_pt()->nboundary_node(0);
 for (unsigned i=0;i<nnode;i++)
  {
   double x = Fluid_mesh_pt->boundary_node_pt(0,i)->x(0);
   double y = Fluid_mesh_pt->boundary_node_pt(0,i)->x(1);
   double r = sqrt(x*x+y*y);
   double v = 2.0*(1.0-r*r);
   v*=fact;
   fluid_mesh_pt()->boundary_node_pt(0,i)->set_value(2,v);
  }

}// end of actions_before_implicit_timestep


//============start_doc_solid_zeta=======================================
/// Doc boundary coordinates of i-th solid FSI boundary.
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
doc_solid_boundary_coordinates()
{
 
 //Doc boundary coordinates in fluid
 std::ofstream the_file("RESLT/solid_boundary_coordinates.dat");
 
 // Loop over traction elements
 unsigned n_face_element = Solid_fsi_traction_mesh_pt->nelement();
 for(unsigned e=0;e<n_face_element;e++)
  {
   //Cast the element pointer
   RefineableFSISolidTractionElement<SOLID_ELEMENT,3>* el_pt=
    dynamic_cast<RefineableFSISolidTractionElement<SOLID_ELEMENT,3>*>
    (Solid_fsi_traction_mesh_pt->element_pt(e));
   
   // Doc boundary coordinate
   Vector<double> s(2);
   Vector<double> zeta(2);
   Vector<double> x(3);
   unsigned n_plot=3;
   the_file << el_pt->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {         
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,n_plot,s);         
     el_pt->interpolated_zeta(s,zeta);
     el_pt->interpolated_x(s,x);
     for (unsigned i=0;i<3;i++)
      {
       the_file << x[i] << " ";
      }
     for (unsigned i=0;i<2;i++)
      {
       the_file << zeta[i] << " ";
      }

     the_file << std::endl;
    }
   el_pt->write_tecplot_zone_footer(the_file,n_plot);
  } 

 // Close doc file
 the_file.close();

} // end doc_solid_zeta


//========start_of_doc_solution===========================================
/// Doc the solution
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 
 std::ofstream some_file;
 std::ostringstream filename;

 // Number of plot points
 unsigned npts;
 npts=5;
 
 // Output solid boundaries
 //------------------------
 filename << doc_info.directory() << "/solid_boundaries"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Solid_mesh_pt->output_boundaries(some_file);
 some_file.close();
 
 
 // Output solid solution
 //-----------------------
 filename.str("");
 filename << doc_info.directory() << "/solid_soln"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 
 // Output fluid boundaries
 //------------------------
 filename.str("");
 filename << doc_info.directory() << "/fluid_boundaries"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Fluid_mesh_pt->output_boundaries(some_file);
 some_file.close();
 
 
 // Output fluid solution
 //-----------------------
 filename.str("");
 filename << doc_info.directory() << "/fluid_soln"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();
  
   
 // Output fsi traction
 //--------------------
 filename.str("");
 filename << doc_info.directory() << "/fsi_traction"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Solid_fsi_traction_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output fsi traction
 //--------------------
 filename.str("");
 filename << doc_info.directory() << "/solid_traction"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Solid_traction_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output Lagrange multipliers
 //----------------------------
 filename.str("");
 filename << doc_info.directory() << "/lagrange"
          << doc_info.number() << ".dat";
 some_file.open(filename.str().c_str());
 Lagrange_multiplier_mesh_pt->output(some_file,npts);
 some_file.close();

} // end_of_doc




//=========================================================================
// steady run
//=========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
steady_run(DocInfo& doc_info)
{
  
 //Output initial configuration
 this->doc_solution(doc_info);
 doc_info.number()++;   
 
 // increment P
 double p_increment=5.0e-4;
 unsigned count=0;
 while (this->wall_mesh_pt()->min_deformation() > 0.9)
  {
   oomph_info << "Min wall position: "  
              << this->wall_mesh_pt()->min_deformation() << std::endl;
   Global_Parameters::P+=p_increment;

   cout << "Doing steady solve: " << count << std::endl;
   if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
    {
     this->steady_newton_solve();
    }
   else
    {
     unsigned max_adapt=1;
     this->steady_newton_solve(max_adapt);
    }
   this->doc_parameters();
   this->doc_solution(doc_info);
   doc_info.number()++;
   count++;
   if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
    {
     oomph_info << "Just doing one steady solve during validation.\n";
     break;
    }
  }

}


//=========================================================================
// unsteady run
//=========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
unsteady_run(DocInfo& doc_info)
{

 // Set number of timesteps
 unsigned nstep=Global_Parameters::Nstep_per_period*Global_Parameters::Nperiod;
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   oomph_info << "Just doing one timestep during validation.\n";
   nstep=1;
  }

 // Do unsteady run
 for (unsigned r = 0; r < nstep; r++)
  {
   cout << "Doing unsteady solve: " << r << std::endl;
   this->unsteady_newton_solve(Global_Parameters::Dt);
   this->doc_parameters();
   this->doc_solution(doc_info);
   doc_info.number()++;
  }

}


//==============================================================================
/// Switch to GMRES preconditioned with the pseudo-elastic 
/// FSI preconditioner
//==============================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void PseudoElasticCollapsibleChannelProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
set_pseudo_elastic_fsi_solver()
{

//setup the solver
#ifdef OOMPH_HAS_TRILINOS
  
  TrilinosAztecOOSolver* solver_pt = new
   TrilinosAztecOOSolver;
  solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
  
#else
  
  GMRES<CRDoubleMatrix>* solver_pt = new GMRES<CRDoubleMatrix>;
  
#endif

  // Set tolerance
  solver_pt->tolerance()=1.0e-8;

  // Max. number of iterations
  solver_pt->max_iter()=200;
  
  // preconditioner
  PseudoElasticFSIPreconditioner* prec_pt = new
   PseudoElasticFSIPreconditioner(3,this);
  
 // meshes
  prec_pt->set_fluid_and_pseudo_elastic_mesh_pt(fluid_mesh_pt());
  prec_pt->set_solid_mesh_pt(wall_mesh_pt());
  prec_pt->set_lagrange_multiplier_mesh_pt(lagrange_multiplier_mesh_pt());
  
  // inexact pseudo-solid preconditioning
  prec_pt->pseudo_elastic_preconditioner_pt()->elastic_preconditioner_type()
   = PseudoElasticPreconditioner::
   Block_upper_triangular_preconditioner;

#ifdef OOMPH_HAS_HYPRE

  prec_pt->pseudo_elastic_preconditioner_pt()
   ->set_elastic_subsidiary_preconditioner
   (Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper
    ::get_elastic_preconditioner);

#endif


#ifdef OOMPH_HAS_TRILINOS

  prec_pt->pseudo_elastic_preconditioner_pt()
   ->set_lagrange_multiplier_subsidiary_preconditioner
   (Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper
    ::get_lagrange_multiplier_preconditioner);
  
#endif

  // inexact "real" solid preconditioning
  BlockTriangularPreconditioner<CRDoubleMatrix>*
   solid_prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;
  prec_pt->set_solid_preconditioner(solid_prec_pt);

#ifdef OOMPH_HAS_HYPRE

  solid_prec_pt->set_subsidiary_preconditioner_function
   (Real_Solid_Preconditioner_Helper::get_preconditioner);

#endif
  
  // inexact navier stokes preconditioning
  NavierStokesSchurComplementPreconditioner*
   ns_prec_pt = prec_pt->navier_stokes_schur_complement_preconditioner_pt();
  prec_pt->enable_navier_stokes_schur_complement_preconditioner();
  
  // ns momentum
  BlockDiagonalPreconditioner<CRDoubleMatrix>*
   f_prec_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;

#ifdef OOMPH_HAS_HYPRE

  f_prec_pt->set_subsidiary_preconditioner_function
   (LSC_Preconditioner_Helper::set_hypre_preconditioner);

#endif

  ns_prec_pt->set_f_preconditioner(f_prec_pt);
  
#ifdef OOMPH_HAS_HYPRE

  // ns pressure poisson
  HyprePreconditioner* p_prec_pt = new HyprePreconditioner;
  p_prec_pt->set_amg_iterations(2);
  p_prec_pt->amg_using_simple_smoothing();
  p_prec_pt->amg_simple_smoother() = 3;
  p_prec_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
  p_prec_pt->amg_strength() = 0.25;
  // Use 6 for parallel (Falgout) or 3 (Ruge Stuben) for serial
  p_prec_pt->amg_coarsening() = 6; 
  p_prec_pt->disable_doc_time();
  ns_prec_pt->set_p_preconditioner(p_prec_pt);

#endif
  
  // and pass prec to solver and solver to problem
  solver_pt->preconditioner_pt() = prec_pt;
  linear_solver_pt() = solver_pt;
 }
 
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////



//========================= start_of_main=================================
/// Demonstrate how to solve a 3D FSI problem with pseudo-solid
/// node update for the fluid mesh.
//========================================================================
int main(int argc, char* argv[])                                        
{       
                                                                
#ifdef OOMPH_HAS_MPI                                                    
   MPI_Helpers::init(argc,argv);                                          
#endif         

 // Switch off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Doc info
 DocInfo doc_info;

 // Use separate directory for output from each processor
 std::ostringstream dir_name;
 dir_name << "RESLT_proc" << MPI_Helpers::communicator_pt()->my_rank();
 doc_info.set_directory(dir_name.str());

 // Define processor-labeled output file for all on-screen stuff
 std::ofstream output_stream;
 char filename[1000];
#ifdef OOMPH_HAS_MPI
 sprintf(filename,"%s/OUTPUT.%i",doc_info.directory().c_str(),
         MPI_Helpers::communicator_pt()->my_rank());
#else
 sprintf(filename,"%s/OUTPUT.%i",doc_info.directory().c_str(),0);
#endif

 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);   

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validate");
 
 // Use iterative solver?
 CommandLineArgs::specify_command_line_flag("--use_iterative_solver");
  
 // Doc available command line flags
 CommandLineArgs::doc_available_flags();

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 Global_Parameters::Constitutive_law_wall_pt = 
  new GeneralisedHookean(&Global_Parameters::Nu_wall);
 
Global_Parameters::Constitutive_law_pseudo_elastic_pt = 
  new GeneralisedHookean(&Global_Parameters::Nu_pseudo_elastic);

//Set up the problem
PseudoElasticCollapsibleChannelProblem
 <RefineablePseudoSolidNodeUpdateElement
  <RefineableQTaylorHoodElement<3>, 
   PseudoElasticBulkElement<RefineableQPVDElement<3,3> > >,
  RefineableQPVDElement<3,3> >  problem; 


#ifdef OOMPH_HAS_MPI
 // Distribute the problem
 oomph_info << "Problem is being distributed." << std::endl;

 // Validation: manufacture distribution so that domain is split in two in a 
 // controllable way
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   unsigned nel=problem.mesh_pt()->nelement();
   Vector<unsigned> element_partition(nel);
   for (unsigned e=0;e<nel;e++)
    {
     if (problem.mesh_pt()->finite_element_pt(e)->node_pt(0)->x(2)<
         0.5*Global_Parameters::L)
      {
       element_partition[e]=0;
      }
     else
      {
       element_partition[e]=1;
      }
    }
   problem.distribute(element_partition);
  }
 // Use METIS to determine the partitioning
 else
  {
   problem.distribute(); 
  }
 oomph_info << "Problem has been distributed." << std::endl;
#endif
 // Doc initial configuration
 problem.doc_solution(doc_info);
 oomph_info << "Before manual refine: doc_info.number()=" 
            << doc_info.number() << std::endl;
 doc_info.number()++;

 // Do manual non-uniform refinement
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   problem.adapt(); 
  }

 // Doc afterwards
 problem.doc_solution(doc_info);
 oomph_info << "After manual refine: doc_info.number()=" 
            << doc_info.number() << std::endl;
 doc_info.number()++;
 
 // Setup timestepping
 double dt = Global_Parameters::T/double(Global_Parameters::Nstep_per_period);
 Global_Parameters::Dt=dt;
 problem.initialise_dt(dt);
 problem.set_initial_condition();


 // Switch solver if required
 if (CommandLineArgs::command_line_flag_has_been_set("--use_iterative_solver"))
  {
   oomph_info << "switching to iterative solver\n";
   problem.set_pseudo_elastic_fsi_solver();
  }
 else
  {
   oomph_info << "Using direct solver\n";
  }

 // Steady run
 problem.steady_run(doc_info);

 // Unteady run
 problem.unsteady_run(doc_info);

                                                         
#ifdef OOMPH_HAS_MPI                                                    
 MPI_Helpers::finalize();
#endif         
 
} // end_of_main

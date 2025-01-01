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

// Generic oomph-lib includes
#include "generic.h"
#include "beam.h"
#include "navier_stokes.h"
#include "multi_physics.h"
#include "constitutive.h"
#include "solid.h"

// The meshes
#include "meshes/channel_with_leaflet_mesh.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;
using namespace oomph;


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

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
 
 /// Return the number of DOF types associated with this element. 
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


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



#ifdef OOMPH_HAS_HYPRE

//=========================================================================
/// Namespace for Navier Stokes LSC Preconditioner
//=========================================================================
namespace LSC_Preconditioner_Helper
{

 /// Create instance of Hypre preconditioner with settings that are
 /// appropriate for serial solution of Navier-Stokes momentum block
 Preconditioner* set_hypre_preconditioner()
 {
  HyprePreconditioner* hypre_preconditioner_pt
   = new HyprePreconditioner;
  hypre_preconditioner_pt->set_amg_iterations(2);
  hypre_preconditioner_pt->amg_using_simple_smoothing();
  hypre_preconditioner_pt->amg_simple_smoother() = 0;
  hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
  hypre_preconditioner_pt->amg_strength() = 0.25;
  hypre_preconditioner_pt->amg_coarsening() = 3;
  hypre_preconditioner_pt->amg_damping() = 2.0/3.0;
  return hypre_preconditioner_pt;
 }
}

#endif

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////




//==========================================================================
/// Global parameters
//==========================================================================
namespace Global_Parameters
{
 /// x-position of root of leaflet
 double Leaflet_x0 = 1.0; 

 /// height of leaflet
 double Leaflet_height=0.5;

 /// length of fluid mesh to left of leaflet
 double Fluid_length_left=1.0; 

 /// length of fluid mesh to right of leaflet
 double Fluid_length_right=3.0; 

 /// height of fluid mesh
 double Fluid_height=1.0;

 /// Num elements to left of leaflet in coarse mesh
 unsigned Mesh_nleft=4;

 /// Num elements to right of leaflet in coarse mesh
 unsigned Mesh_nright=12;

 /// Num elements in fluid mesh in y dirn adjacent to leaflet
 unsigned Mesh_ny1=2;

 /// Num elements in fluid mesh in y dirn  above leaflet
 unsigned Mesh_ny2=2; 

 /// Reynolds number
 double Re=50.0;

 /// Womersley number: Product of Reynolds and Strouhal numbers
 double ReSt=50.0;

 /// Non-dimensional wall thickness.
 double H=0.05;
 
 /// Fluid structure interaction parameter: Ratio of stresses used 
 /// for non-dimensionalisation of fluid to solid stresses. 
 double Q=2.0e-7;

 /// Period for fluctuations in flux
 double T=1.0;

 /// Min. flux
 double U_base=1.0;
 
 /// Max. flux
 double U_perturbation=0.5;
 
 /// Flux: Pulsatile flow 
 double flux(const double& t)
 {  
  return U_base+U_perturbation*cos(2.0*MathematicalConstants::Pi*t/T);
 }

 /// Pseudo-solid mass density
 double Lambda_sq=0.0;

 /// Beam mass density
 double Lambda_sq_beam=0.0;

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 /// Timestep for simulation: 40 steps per period
 double Dt = T/40.0;
 
} // end_of_namespace


/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////


//==========================================================================
/// GeomObject: Undeformed straight, vertical leaflet
//==========================================================================
class UndeformedLeaflet : public GeomObject
{

public:

 /// Constructor: argument is the x-coordinate of the leaflet
 UndeformedLeaflet(const double& x0): GeomObject(1,2)
  {
   X0=x0;
  }
 
 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = X0;
   r[1] = zeta[0];
  }


 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Calls steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Use the steady version
   position(zeta,r);
  } // end of position


 /// Posn vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 void d2position(const Vector<double>& zeta,
                 Vector<double>& r,
                 DenseMatrix<double> &drdzeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   // Position vector
   r[0] = X0;
   r[1] = zeta[0];

   // Tangent vector
   drdzeta(0,0)=0.0;
   drdzeta(0,1)=1.0;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=0.0;
  } // end of d2position

 /// Number of geometric Data in GeomObject: None.
 unsigned ngeom_data() const {return 0;}  

private :

 /// x position of the undeformed leaflet's origin. 
 double X0;

}; //end_of_undeformed_wall


/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////


//==========================================================================
/// FSI leaflet in channel. Mesh update with pseudo-elasticity and
/// solved with pseudo-elastic fsi preconditioner.
//==========================================================================
template<class ELEMENT>
class FSIChannelWithLeafletProblem : public Problem
{

public:

 /// Constructor: Pass multiplier for uniform mesh refinement
 FSIChannelWithLeafletProblem(const unsigned& mesh_multiplier);

 /// Destructor empty
 ~FSIChannelWithLeafletProblem()
  {
   delete Bulk_mesh_pt;
   delete Lagrange_multiplier_mesh_pt;
   delete Wall_mesh_pt;
   delete Bulk_time_stepper_pt;
   delete Wall_time_stepper_pt;
   delete Wall_geom_object_pt;
   delete Undeformed_wall_pt;
   delete Constitutive_law_pt;
  }
 
 /// Actions after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before Newton solve:
 /// Reset the  pseudo-elastic undeformed configuration
 void actions_before_newton_solve()
  {
   // Reset undeformed configuration for pseudo-solid
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }
 
 /// Update no slip before Newton convergence check
 void actions_before_newton_convergence_check()
  {
   // Loop over the nodes to perform auxiliary node update (no slip) 
   unsigned nnod=Bulk_mesh_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Bulk_mesh_pt->node_pt(j)->perform_auxiliary_node_update_fct();
    }
     
  }
 
 /// Actions before implicit timestep: Update the inflow velocity
 void actions_before_implicit_timestep()
  {
   // Actual time
   double t=time_pt()->time();
   
   // Amplitude of flow
   double ampl=Global_Parameters::flux(t);
   
   // Update parabolic flow along boundary 3
   unsigned ibound=3; 
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     double ycoord = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
     double uy = ampl*4.0*ycoord/Global_Parameters::Fluid_height*
      (1.0-ycoord/Global_Parameters::Fluid_height);
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
    }
  } // end of actions_before_implicit_timestep
 
 /// Set iterative solver 
 void set_iterative_solver();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Create elements that enforce prescribed boundary motion
 /// by Lagrange multipliers
 void create_lagrange_multiplier_elements();
 
 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multipliers
 void delete_lagrange_multiplier_elements();
 
 /// Doc parameters
 void doc_parameters()
  {
   oomph_info << "\n\n=================================================\n";
   oomph_info << "Q                    : " << Global_Parameters::Q 
              << std::endl;
   oomph_info << "Re                   : " << Global_Parameters::Re 
              << std::endl;
   oomph_info << "Lambda_sq_beam       : " << Global_Parameters::Lambda_sq_beam 
              << std::endl;
   oomph_info << "T                    : " << Global_Parameters::T 
              <<std::endl;
   oomph_info << "t                    : " << time_pt()->time()
              << std::endl;
   oomph_info << "tip x                : " 
              << tip_node_pt()->x(0) << std::endl;
   oomph_info << "tip y                : " 
              << tip_node_pt()->x(1) << std::endl;
   oomph_info << "=================================================\n\n";
   }

private:
 
 /// Helper fct; returns the node at the tip of the wall mesh
 Node* tip_node_pt()
  {
   unsigned n_el_wall=Wall_mesh_pt->nelement();
   return Wall_mesh_pt->finite_element_pt(n_el_wall-1)->node_pt(1);
  }


 /// Pointer to the fluid mesh
 PseudoElasticChannelWithLeafletMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt;
 
 /// Bulk timestepper
 BDF<2>* Bulk_time_stepper_pt;
 
 /// Wall time stepper pt
 Newmark<2>* Wall_time_stepper_pt;
 
 /// Pointers to mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;
 
 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Geometric object for the leaflet (to apply lagrange mult)
 MeshAsGeomObject* Wall_geom_object_pt;
 
 /// Geom object for the leaflet
 UndeformedLeaflet* Undeformed_wall_pt;
  
};


//==========================================================================
/// Constructor
//==========================================================================
template <class ELEMENT>
FSIChannelWithLeafletProblem<ELEMENT>::FSIChannelWithLeafletProblem
(const unsigned& mesh_multiplier)
{

  // Allocate the timesteppers
 Bulk_time_stepper_pt=new BDF<2>;
 add_time_stepper_pt(Bulk_time_stepper_pt);
 Wall_time_stepper_pt=new Newmark<2>;
 add_time_stepper_pt(Wall_time_stepper_pt);
 
 // Wall mesh
 //----------

 // Geometric object that represents the undeformed leaflet
 Undeformed_wall_pt=new UndeformedLeaflet(Global_Parameters::Leaflet_x0);

 //Create the "wall" mesh with FSI Hermite beam elements
 unsigned n_wall_el=Global_Parameters::Mesh_ny1*mesh_multiplier;
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (n_wall_el,Global_Parameters::Leaflet_height,
   Undeformed_wall_pt,Wall_time_stepper_pt);


 // Fluid mesh
 // ----------

 // Build a geometric object from the wall mesh
 Wall_geom_object_pt=new MeshAsGeomObject(Wall_mesh_pt); 

 //Build the fluid mesh
 Bulk_mesh_pt =new PseudoElasticChannelWithLeafletMesh<ELEMENT>(
  Wall_geom_object_pt,
  Global_Parameters::Fluid_length_left,
  Global_Parameters::Fluid_length_right,
  Global_Parameters::Leaflet_height,
  Global_Parameters::Fluid_height,
  Global_Parameters::Mesh_nleft*mesh_multiplier,
  Global_Parameters::Mesh_nright*mesh_multiplier,
  Global_Parameters::Mesh_ny1*mesh_multiplier,
  Global_Parameters::Mesh_ny2*mesh_multiplier,
  Bulk_time_stepper_pt);


 // Construct the mesh of elements that enforce prescribed boundary motion
 //-----------------------------------------------------------------------
 // by Lagrange multipliers
 //------------------------
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();


 // Build global mesh 
 //------------------
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Lagrange_multiplier_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);
 build_global_mesh();
 


 // Fluid boundary conditions
 //--------------------------
 
 // Apply no-slip condition on all boundary nodes of the fluid mesh
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
     
     // Do not pin the x velocity of the outflow
     if( ibound != 1)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }      
    }
  }
 
 // Pin the nodal position on all boundaries apart from 
 // the moving wall
 for (unsigned ibound=0;ibound<7;ibound++)
  {
   if (ibound==0||ibound==1||ibound==2||ibound==3||ibound==6)
    {
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<2;i++)
        {
         if (!( (ibound==2)&&(i==0) ))
          {
           dynamic_cast<SolidNode*>(Bulk_mesh_pt->
                                    boundary_node_pt(ibound, inod))
            ->pin_position(i);
          }
        }
      }
    }
  }

 // Setup parabolic flow along boundary 3 (everything else that's 
 // pinned has homogeneous boundary conditions so no action is required
 // as that's the default assignment). 
 unsigned ibound=3; 
 unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
   double uy = 4.0*ycoord/Global_Parameters::Fluid_height*
    (1.0-ycoord/Global_Parameters::Fluid_height);
   Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
   Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }
 
 // Boundary conditions for wall mesh
 // ---------------------------------
 
 // Set the boundary conditions: the lower end of the beam is fixed in space
 unsigned b=0; 
 
 // Pin displacements in both x and y directions
 Wall_mesh_pt->boundary_node_pt(b,0)->pin_position(0); 
 Wall_mesh_pt->boundary_node_pt(b,0)->pin_position(1);
 
 // Infinite slope: Pin type 1 (slope) dof for displacement direction 0 
 Wall_mesh_pt->boundary_node_pt(b,0)->pin_position(1,0);


 // Complete build of fluid elements 
 // --------------------------------

 //Set the constitutive law
 Constitutive_law_pt = new GeneralisedHookean(&Global_Parameters::Nu);
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;
   
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Parameters::ReSt;
      
   //Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;
   
   // Density of pseudo-solid
   el_pt->lambda_sq_pt()=&Global_Parameters::Lambda_sq;
   
  }// end loop over elements
 
 

 // Complete build of wall elements 
 // -------------------------------
 n_element = Wall_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(Wall_mesh_pt->element_pt(e));
   
   // Set physical parameters for each element:
   elem_pt->h_pt() = &Global_Parameters::H;
   
   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Parameters::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = Undeformed_wall_pt;

   // Density of beam
   elem_pt->lambda_sq_pt()=&Global_Parameters::Lambda_sq_beam;
   
   // Leaflet is immersed and loaded by fluid on both sides
   elem_pt->enable_fluid_loading_on_both_sides();
   
   // The normal to the leaflet, as computed by the 
   // FSIHermiteElements points out of the fluid rather than 
   // into the fluid (as assumed by default) when viewed from
   // the "front" (face 0).
   elem_pt->set_normal_pointing_out_of_fluid();
   
  }
 
 // Setup FSI
 // ---------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 4,5)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 for(unsigned ibound=4;ibound<6;ibound++ )
  { 
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {   
     Bulk_mesh_pt->boundary_node_pt(ibound, inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  }// aux node update fct has been set



 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 4 and 5
 // of the 2D fluid mesh.

 // Front of leaflet (face=0, which is also the default so this argument
 // could be omitted) meets boundary 4 of the fluid mesh.
 unsigned face=0; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,
   4,
   Bulk_mesh_pt,
   Wall_mesh_pt,
   face); 
 
 // Back of leaflet: face 1, meets boundary 5 of fluid mesh
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Bulk_mesh_pt,Wall_mesh_pt,face); 
 
 // Setup equation numbering scheme
 //--------------------------------
 oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl; 
 
}//end of constructor


//===========start_iterative_solver=========================================
/// Set iterative solver 
//==========================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::set_iterative_solver()
{
 // Create the linear solver
 IterativeLinearSolver* solver_pt=0;
 
 // If we have trilinos, use it
#ifdef OOMPH_HAS_TRILINOS
 
 // Create solver
 solver_pt = new TrilinosAztecOOSolver;
 
 // Use GMRES
 dynamic_cast<TrilinosAztecOOSolver*>(solver_pt)->solver_type() 
  = TrilinosAztecOOSolver::GMRES;
 
#else
 
 // Use oomph-lib's own GMRES
 solver_pt = new GMRES<CRDoubleMatrix>;
 
#endif

 // Set solver
 linear_solver_pt() = solver_pt;
  
 // Create preconditioner for 2D problem
 unsigned dim=2;
 PseudoElasticFSIPreconditioner* prec_pt=
  new PseudoElasticFSIPreconditioner(dim, this);
 
 // Set preconditioner
 solver_pt->preconditioner_pt() = prec_pt;
 
 
 // Specify meshes that contain elements which classify the various
 // degrees of freedom:
 prec_pt->set_fluid_and_pseudo_elastic_mesh_pt(Bulk_mesh_pt);
 prec_pt->set_solid_mesh_pt(Wall_mesh_pt);
 prec_pt->set_lagrange_multiplier_mesh_pt(Lagrange_multiplier_mesh_pt);
 
 
 // Use oomph-lib's Schur complement preconditioner as Navier-Stokes
 // subsidiary preconditioner
 if (!CommandLineArgs::command_line_flag_has_been_set("--suppress_lsc"))
  {
   oomph_info << "Enabling LSC preconditioner\n";
   prec_pt->enable_navier_stokes_schur_complement_preconditioner();
  }
 else
  {
   prec_pt->disable_navier_stokes_schur_complement_preconditioner();
   oomph_info << "Not using LSC preconditioner\n";
  } // done disable lsc


 // Use approximate block solves?
 //------------------------------
 if (CommandLineArgs::command_line_flag_has_been_set("--superlu_for_blocks"))
  {
   oomph_info << "Use SuperLU for block solves\n";
  }
 else
  {
   oomph_info << "Use optimal block solves\n";

   // Get pointer to Navier-Stokes Schur complement preconditioner
   NavierStokesSchurComplementPreconditioner* ns_prec_pt = 
    prec_pt->navier_stokes_schur_complement_preconditioner_pt();

   // Navier Stokes momentum block
   //-----------------------------

   // Block triangular for momentum block in LSC precond
   BlockTriangularPreconditioner<CRDoubleMatrix>*
    f_prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;

   // Set it
   ns_prec_pt->set_f_preconditioner(f_prec_pt);  
   
#ifdef OOMPH_HAS_HYPRE

   // Use Hypre for diagonal blocks
   f_prec_pt->set_subsidiary_preconditioner_function
    (LSC_Preconditioner_Helper::set_hypre_preconditioner);


   // Navier Stokes Schur complement/pressure block
   //----------------------------------------------

   // Build/set Hypre for Schur complement (pressure) block
   HyprePreconditioner* p_prec_pt = new HyprePreconditioner;
   p_prec_pt->disable_doc_time();   
   Hypre_default_settings::set_defaults_for_2D_poisson_problem(p_prec_pt);   
   ns_prec_pt->set_p_preconditioner(p_prec_pt); 

#endif

   // Pseudo elastic block
   //---------------------

   // Use block upper triangular preconditioner for (pseudo-)elastic block
   prec_pt->pseudo_elastic_preconditioner_pt()->elastic_preconditioner_type()
    = PseudoElasticPreconditioner::Block_upper_triangular_preconditioner;

#ifdef OOMPH_HAS_HYPRE

   // Use Hypre for diagonal blocks of (pseudo-)elastic preconditioner
   prec_pt->pseudo_elastic_preconditioner_pt()->
    set_elastic_subsidiary_preconditioner(
     Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper::
     get_elastic_preconditioner_hypre);

#endif
   
#ifdef OOMPH_HAS_TRILINOS

   // Use Trilinos CG as subsidiary preconditioner (inexact solver) for
   // linear (sub-)systems to be solved in the Lagrange multiplier block
   prec_pt->pseudo_elastic_preconditioner_pt()->
    set_lagrange_multiplier_subsidiary_preconditioner
    (Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper::
     get_lagrange_multiplier_preconditioner);

#endif
  }

} //end set_iterative_solver


//==========================================================================
/// Create elements that impose the prescribed boundary displacement
/// by Lagrange multipliers
//==========================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Lagrange multiplier elements are located on boundary 4 and 5
 for (unsigned b=4;b<=5;b++)
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned n_bulk_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_bulk_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     ImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
      new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
       bulk_elem_pt,face_index);
     
     // Add to mesh 
     Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);
   
     // Set the GeomObject that defines the boundary shape and set
     // which bulk boundary we are attached to(needed to extract
     // the boundary coordinate from the bulk nodes)
     el_pt->set_boundary_shape_geom_object_pt(Wall_geom_object_pt,b);
    }  
  }

 // Pin Lagrange multiplier unknowns for fluid nodes whose position
 // is already pinned
 unsigned n_element=Lagrange_multiplier_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a Lagrange multiplier element
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(i));

   // Loop over the nodes 
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = el_pt->node_pt(j);
     
     // Is the node also on boundary 0 or 6 (i.e. on bottom wall)>
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(6)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=el_pt->nbulk_value(j);
       
       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned k=n_bulk_value;k<nval;k++)
        {
         nod_pt->pin(k);
        }
      }
    }
  }
} 



//====start_of_delete_lagrange_multiplier_elements=======================
/// Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//=======================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::
delete_lagrange_multiplier_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Lagrange_multiplier_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();

} // end of delete_lagrange_multiplier_elements



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output fluid solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output wall solution 
 sprintf(filename,"%s/wall_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Wall_mesh_pt->output(some_file,npts);
 some_file.close();

 // Help me figure out what the "front" and "back" faces of the leaflet are
 //------------------------------------------------------------------------

 // Output fluid elements on fluid mesh boundary 4 (associated with
 // the "front")
 unsigned bound=4;
 sprintf(filename,"%s/bulk_boundary_elements_front_%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned nel= Bulk_mesh_pt->nboundary_element(bound);
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(bound,e))
    ->output(some_file,npts);
  }
 some_file.close();


 // Output fluid elements on fluid mesh boundary 5 (associated with
 // the "back")
 bound=5;
 sprintf(filename,"%s/bulk_boundary_elements_back_%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 nel= Bulk_mesh_pt->nboundary_element(bound);
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(bound,e))
    ->output(some_file,npts);
  }
 some_file.close();


 // Output normal vector on wall elements
 sprintf(filename,"%s/wall_normal_%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 nel=Wall_mesh_pt->nelement();
 Vector<double> s(1);
 Vector<double> x(2);
 Vector<double> xi(1);
 Vector<double> N(2);
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   FSIHermiteBeamElement* el_pt=
    dynamic_cast<FSIHermiteBeamElement*>(Wall_mesh_pt->element_pt(e));

   // Loop over plot points
   for (unsigned i=0;i<npts;i++)
    {
     s[0]=-1.0+2.0*double(i)/double(npts-1);

     // Get Eulerian position
     el_pt->interpolated_x(s,x);

     // Get unit normal
     el_pt->get_normal(s,N);

     // Get Lagrangian coordinate
     el_pt->interpolated_xi(s,xi);
     
     some_file << x[0] << " " << x[1] << " " 
               << N[0] << " " << N[1] << " " 
               << xi[0] << std::endl;
    }
  }
 some_file.close();

} // end_of_doc_solution



//=======start_of_main=====================================================
/// Driver code 
//=========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Switch off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Multiplier for number of elements in coordinate directions.
 // Used for uniform mesh refinement studies.
 unsigned mesh_multiplier = 2; 
 CommandLineArgs::specify_command_line_flag("--mesh_multiplier",
                                            &mesh_multiplier);

 // Suppress use of LSC preconditioner for Navier Stokes block
 CommandLineArgs::specify_command_line_flag("--suppress_lsc");

 // Use direct solver (SuperLU)
 CommandLineArgs::specify_command_line_flag("--use_direct_solver");
 
 // Use SuperLU for all block solves
 CommandLineArgs::specify_command_line_flag("--superlu_for_blocks");

 // Validation only?
 CommandLineArgs::specify_command_line_flag("--validate");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Set up the problem
 FSIChannelWithLeafletProblem<PseudoSolidNodeUpdateElement
                              <QTaylorHoodElement<2>,
                               PseudoElasticBulkElement<QPVDElement<2,3> > > >* 
  problem_pt = new  
  FSIChannelWithLeafletProblem<PseudoSolidNodeUpdateElement
                               <QTaylorHoodElement<2>,
                                PseudoElasticBulkElement<QPVDElement<2,3> > > >
  (mesh_multiplier);

 // Initialise timestep
 problem_pt->initialise_dt(Global_Parameters::Dt);
 
 // Label for output
 DocInfo doc_info;
 doc_info.set_directory("RESLT");


 // Define processor-labeled output file for all on-screen stuff
 std::ofstream output_stream;
 char filename[1000];
#ifdef OOMPH_HAS_MPI
 sprintf(filename,"%s/OUTPUT_STEADY.%i",doc_info.directory().c_str(),
         MPI_Helpers::communicator_pt()->my_rank());
#else
 sprintf(filename,"%s/OUTPUT_STEADY.%i",doc_info.directory().c_str(),0);
#endif

 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);   

 
 // Output initial configuration
 problem_pt->doc_solution(doc_info);
 doc_info.number()++;   

 // Switch to iterative solver?
 if (!CommandLineArgs::command_line_flag_has_been_set("--use_direct_solver"))
  {
   problem_pt->set_iterative_solver();
  }


 // Steady solves
 //--------------

 // Increment Re and Womersley numbers in increments of 25. 
 double target_re = Global_Parameters::Re;
 Global_Parameters::Re=25.0;
 Global_Parameters::ReSt=25.0;
 while (Global_Parameters::Re<target_re)
  {
   problem_pt->steady_newton_solve();
   problem_pt->doc_parameters();
   Global_Parameters::Re+=25.0;
   Global_Parameters::ReSt+=25.0;
   problem_pt->doc_solution(doc_info);
   doc_info.number()++;
  }

 // Do final solve at desired Re
 Global_Parameters::Re=target_re;
 Global_Parameters::ReSt=target_re;
 problem_pt->steady_newton_solve();
 problem_pt->doc_parameters();
 problem_pt->doc_solution(doc_info);
 doc_info.number()++;

 // Unsteady solves
 //----------------

 // Define processor-labeled output file for all on-screen stuff
 output_stream.close();
#ifdef OOMPH_HAS_MPI
 sprintf(filename,"%s/OUTPUT_UNSTEADY.%i",doc_info.directory().c_str(),
         MPI_Helpers::communicator_pt()->my_rank());
#else
 sprintf(filename,"%s/OUTPUT_UNSTEADY.%i",doc_info.directory().c_str(),0);
#endif
 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);   

 // Loop over timesteps for specified number of periods of fluctuating
 // inflow
 unsigned n_period=1;

 unsigned nstep=unsigned(double(n_period)
                         *Global_Parameters::T/Global_Parameters::Dt);

 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   nstep=3;
  }
 for (unsigned r = 0; r < nstep; r++)
  {
   problem_pt->unsteady_newton_solve(Global_Parameters::Dt);
   problem_pt->doc_parameters();
   problem_pt->doc_solution(doc_info);
   doc_info.number()++;
  }

 // clean up
 delete problem_pt;

 // Shut down
 oomph_info << "Done\n";

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} // end_of_main



 

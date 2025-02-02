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

//Include the mesh
#include "meshes/channel_with_leaflet_mesh.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"


using namespace std;
using namespace oomph;


//==== start_of_global_parameters================================
/// Global parameters
//===============================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=50.0;

 /// Womersley number: Product of Reynolds and Strouhal numbers
 double ReSt=50.0;

 /// Non-dimensional wall thickness.
 double H=0.05;
 
 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=1.0e-6;

 /// Period for fluctuations in flux
 double Period=2.0;

 /// Min. flux
 double Min_flux=1.0;
 
 /// Max. flux
 double Max_flux=2.0;
 
 /// Flux: Pulsatile flow fluctuating between Min_flux and Max_flux 
 /// with period Period
 double flux(const double& t)
 {  
  return Min_flux+
   (Max_flux-Min_flux)*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*t/Period));
 }

} // end_of_namespace



/// ////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=====start_of_undeformed_leaflet=================================
/// GeomObject: Undeformed straight, vertical leaflet
//=================================================================
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


/// ////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=====start_of_problem_class========================================
/// FSI leaflet in channel
//===================================================================
template<class ELEMENT>
class FSIChannelWithLeafletProblem : public Problem
{

public:

 /// Constructor: Pass the lenght of the domain at the left
 /// of the leaflet lleft,the lenght of the domain at the right of the
 /// leaflet lright,the height of the leaflet hleaflet, the total height
 /// of the domain htot, the number of macro-elements at the left of the
 /// leaflet nleft, the number of macro-elements at the right of the
 /// leaflet nright, the number of macro-elements under hleaflet ny1,
 /// the number of macro-elements above hleaflet ny2, the abscissa
 /// of the origin of the leaflet x_0.
 FSIChannelWithLeafletProblem(const double& lleft,
                              const double& lright, const double& hleaflet,
                              const double& htot,
                              const unsigned& nleft, const unsigned& nright,
                              const unsigned& ny1, const unsigned&  ny2,
                              const double& x_0);  

 /// Destructor empty
 ~FSIChannelWithLeafletProblem(){}
 
 /// Actions after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before solve (empty) 
 void actions_before_newton_solve(){}

 /// Actions after adaptation
 void actions_after_adapt();

 /// Access function to the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;
  } 

 /// Access function to fluid mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* fluid_mesh_pt()
  {
   return Fluid_mesh_pt;
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace);


/// Update the inflow velocity
 void actions_before_implicit_timestep()
  {
   // Actual time
   double t=time_pt()->time();

   // Amplitude of flow
   double ampl=Global_Physical_Variables::flux(t);

   // Update parabolic flow along boundary 3
   unsigned ibound=3; 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
    double uy = ampl*6.0*ycoord/Htot*(1.0-ycoord/Htot);
    Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
    Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
   }
  } // end of actions_before_implicit_timestep

 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape
 void actions_before_newton_convergence_check()
  {
   Fluid_mesh_pt->node_update();
  }

private:
 
 /// Pointer to the fluid mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* Fluid_mesh_pt;

 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt;

 /// Pointer to the GeomObject that represents the wall
 GeomObject* Leaflet_pt;
 
 /// Total height of the domain
 double Htot;

};





//=====start_of_constructor==============================================
/// Constructor
//=======================================================================
template <class ELEMENT>
FSIChannelWithLeafletProblem<ELEMENT>::FSIChannelWithLeafletProblem(
 const double& lleft,
 const double& lright,
 const double& hleaflet,
 const double& htot,
 const unsigned& nleft,
 const unsigned& nright,
 const unsigned& ny1,
 const unsigned&  ny2,
 const double& x_0) : Htot(htot)
{
 // Timesteppers:
 //--------------

 // Allocate the timestepper
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;
 add_time_stepper_pt(fluid_time_stepper_pt);

 // Allocate the wall timestepper
 Steady<2>* wall_time_stepper_pt=new Steady<2>;
 add_time_stepper_pt(wall_time_stepper_pt);


 // Discretise leaflet
 //-------------------

 // Geometric object that represents the undeformed leaflet
 UndeformedLeaflet* undeformed_wall_pt=new UndeformedLeaflet(x_0);

 //Create the "wall" mesh with FSI Hermite beam elements
 unsigned n_wall_el=5;
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (n_wall_el,hleaflet,undeformed_wall_pt,wall_time_stepper_pt);


 // Provide GeomObject representation of leaflet mesh and build fluid mesh
 //-----------------------------------------------------------------------

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 MeshAsGeomObject* wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 

//Build the mesh
 Fluid_mesh_pt =new RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>(
  wall_geom_object_pt,
  lleft, lright,
  hleaflet,
  htot,nleft,
  nright,ny1,ny2,
  fluid_time_stepper_pt);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;



 // Build global mesh
 //------------------
 
 // Add the sub meshes to the problem
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();



 // Fluid boundary conditions
 //--------------------------

 //Pin the boundary nodes of the fluid mesh
 unsigned num_bound = Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
      Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
     
      // Do not pin the x velocity of the outflow
      if( ibound != 1)
      {
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }      
    }
  }
 // end loop over boundaries
   

 // Setup parabolic flow along boundary 3 (everything else that's 
 // pinned has homogenous boundary conditions so no action is required
 // as that's the default assignment). Inflow profile is parabolic
 // and this is interpolated correctly during mesh refinement so
 // no re-assignment necessary after adaptation.
 unsigned ibound=3; 
 unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
   double uy = 6.0*ycoord/htot*(1.0-ycoord/htot);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }// end of setup boundary condition


 
 // Boundary conditions for wall mesh
 //----------------------------------

 // Set the boundary conditions: the lower end of the beam is fixed in space
 unsigned b=0; 

 // Pin displacements in both x and y directions
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
 
 // Infinite slope: Pin type 1 (slope) dof for displacement direction 0 
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1,0);
 


 // Complete build of fluid elements
 //---------------------------------
 unsigned n_element = Fluid_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
      
  }// end loop over elements


 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());
 
 
 // Complete build of wall elements
 //--------------------------------
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));
    
   // Set physical parameters for each element:
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;

   // Leaflet is immersed and loaded by fluid on both sides
   elem_pt->enable_fluid_loading_on_both_sides();

   // The normal to the leaflet, as computed by the 
   // FSIHermiteElements points away from the fluid rather than 
   // into the fluid (as assumed by default) when viewed from
   // the "front" (face 0).
   elem_pt->set_normal_pointing_out_of_fluid();

  } // end of loop over elements

 
 // Setup FSI
 //----------
 
 // The velocity of the fluid nodes on the wall (fluid mesh boundary 4,5)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 for(unsigned ibound=4;ibound<6;ibound++ )
  { 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {   
     Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  }// aux node update fct has been set
 
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 4 and 5
 // of the 2D fluid mesh.

 // Front of leaflet: Set face=0 (which is also the default so this argument
 // could be omitted)
 unsigned face=0; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,4,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 



 // Choose iterative solvers with FSI preconditioner (only if not test)
 //====================================================================
 if (CommandLineArgs::Argc==1)
  {
   
   // Build iterative linear solver
   //------------------------------
   GMRES<CRDoubleMatrix>* iterative_linear_solver_pt =
    new GMRES<CRDoubleMatrix>;
   
   // Set maximum number of iterations
   iterative_linear_solver_pt->max_iter() = 200;
   
   // Pass solver to problem:
   linear_solver_pt()=iterative_linear_solver_pt;
   
   
   // Build preconditioner
   //---------------------
   FSIPreconditioner* prec_pt=new FSIPreconditioner(this);
   
   // Set Navier Stokes mesh:
   prec_pt->set_navier_stokes_mesh(Fluid_mesh_pt);
   
   // Set solid mesh:
   prec_pt->set_wall_mesh(Wall_mesh_pt);
   
   // By default, the Schur complement preconditioner uses SuperLU as
   // an exact preconditioner (i.e. a solver) for the
   // momentum and Schur complement blocks.
   // Can overwrite this by passing pointers to
   // other preconditioners that perform the (approximate)
   // solves of these blocks

  
#ifdef OOMPH_HAS_HYPRE

   // Create internal preconditioner used on Schur block
   //---------------------------------------------------
   HyprePreconditioner* p_preconditioner_pt = new HyprePreconditioner;
   
   // Shut up!
   p_preconditioner_pt->disable_doc_time();
   
   // Set defaults parameters for use as preconditioner on Poisson-type problem
   Hypre_default_settings::
    set_defaults_for_2D_poisson_problem(p_preconditioner_pt);
   
   // Pass to preconditioner
   //prec_pt->set_p_preconditioner(p_preconditioner_pt);
   
   
   // Create internal preconditioner used on momentum block
   //------------------------------------------------------
   HyprePreconditioner* f_preconditioner_pt = new HyprePreconditioner;
   
   // Shut up!
   f_preconditioner_pt->disable_doc_time();
   
   // Set default parameters for use as preconditioner in for momentum 
   // block in Navier-Stokes problem
   Hypre_default_settings::
    set_defaults_for_navier_stokes_momentum_block(f_preconditioner_pt);
   
   // Use Hypre for momentum block
   //prec_pt->set_f_preconditioner(f_preconditioner_pt);

#endif
   
   // Retain fluid onto solid terms in FSI preconditioner
   prec_pt->use_block_triangular_version_with_fluid_on_solid();

   // Pass preconditioner to iterative linear solver
   iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
     
  }

  
}//end of constructor





//==== start_of_actions_after_adapt=================================
/// Actions_after_adapt()
//==================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::actions_after_adapt()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(Fluid_mesh_pt->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());
 

 // (Re-)apply the no slip condition on the moving wall
 //-----------------------------------------------------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 4,5)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 for(unsigned ibound=4;ibound<6;ibound++ )
  { 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {   
     Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  } // aux node update fct has been (re-)set

 
 
 // Re-setup FSI
 //-------------
 
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 4 and 5
 // of the 2D fluid mesh.

 // Front of leaflet: Set face=0 (which is also the default so this argument
 // could be omitted)
 unsigned face=0; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,4,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 
  
} // end_of_actions_after_adapt





//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                                         ofstream& trace)
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
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output wall solution 
 sprintf(filename,"%s/wall_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Wall_mesh_pt->output(some_file,npts);
 some_file.close();

 // Get node at tip of leaflet
 unsigned n_el_wall=Wall_mesh_pt->nelement();
 Node* tip_node_pt=Wall_mesh_pt->finite_element_pt(n_el_wall-1)->node_pt(1);

 // Get time:
 double time=time_pt()->time();

 // Write trace file
 trace << time << " " 
       << Global_Physical_Variables::flux(time) << " " 
       << tip_node_pt->x(0) << " "
       << tip_node_pt->x(1) << " "
       << tip_node_pt->dposition_dt(0) << " "
       << tip_node_pt->dposition_dt(1) << " "
       << doc_info.number() << " " 
       << std::endl;


 // Help me figure out what the "front" and "back" faces of the leaflet are
 //------------------------------------------------------------------------

 // Output fluid elements on fluid mesh boundary 4 (associated with
 // the "front")
 unsigned bound=4;
 sprintf(filename,"%s/fluid_boundary_elements_front_%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned nel= Fluid_mesh_pt->nboundary_element(bound);
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Fluid_mesh_pt->boundary_element_pt(bound,e))
    ->output(some_file,npts);
  }
 some_file.close();


 // Output fluid elements on fluid mesh boundary 5 (associated with
 // the "back")
 bound=5;
 sprintf(filename,"%s/fluid_boundary_elements_back_%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 nel= Fluid_mesh_pt->nboundary_element(bound);
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Fluid_mesh_pt->boundary_element_pt(bound,e))
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
 



//======= start_of_main================================================
/// Driver code  -- pass a command line argument if you want to run
/// the code in validation mode where it only performs a few steps
//=====================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Parameters for the leaflet: x-position of root and height
 double x_0 = 1.0; 
 double hleaflet=0.5;

 // Number of elements in various regions of mesh
 unsigned nleft=6; 
 unsigned nright=18;
 unsigned ny1=3; 
 unsigned ny2=3; 

 // Dimensions of fluid mesh: length to the left and right of leaflet
 // and total height
 double lleft =1.0; 
 double lright=3.0; 
 double htot=1.0;
  
 //Build the problem
 FSIChannelWithLeafletProblem<
  AlgebraicElement<RefineableQTaylorHoodElement<2> > >
  problem(lleft,lright,hleaflet,
          htot,nleft,nright,ny1,ny2,x_0); 

 // Set up doc info
 DocInfo doc_info; 
 doc_info.set_directory("RESLT");

 // Trace file
 ofstream trace;
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace.open(filename);


 // Number of timesteps (reduced for validation)
 unsigned nstep=200;
 if (CommandLineArgs::Argc>1)
  {
   nstep=2;
  }

 //Timestep: 
 double dt=0.05;
 
 // Initialise timestep 
 problem.initialise_dt(dt);

 // Doc initial guess for steady solve
 problem.doc_solution(doc_info,trace);
 doc_info.number()++; 

 
 // Initial loop to increment the Reynolds number in sequence of steady solves
 //---------------------------------------------------------------------------
 unsigned n_increment=4;
 // Just to one step for validation run
 if (CommandLineArgs::Argc>1)
  {
   n_increment=1;
  }

 // Set max. number of adaptations
 unsigned max_adapt=3;

 Global_Physical_Variables::Re=0.0;
 for (unsigned i=0;i<n_increment;i++)
  {
   // Increase Re and ReSt (for St=1)
   Global_Physical_Variables::Re+=50.0;
   Global_Physical_Variables::ReSt=Global_Physical_Variables::Re;

   // Solve the steady problem 
   std::cout << "Computing a steady solution for Re=" 
             <<  Global_Physical_Variables::Re << std::endl;
   problem.steady_newton_solve(max_adapt);
   problem.doc_solution(doc_info,trace);
   doc_info.number()++; 
  } // reached final Reynolds number 



 // Proper time-dependent run
 //--------------------------

 // Limit the number of adaptations during unsteady run to one per timestep
 max_adapt=1;
 
 // Don't re-set the initial conditions when adapting the mesh
 bool first = false;

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  { 
   // Solve the problem
   problem.unsteady_newton_solve(dt,max_adapt,first);
   
   // Output the solution
   problem.doc_solution(doc_info,trace);
   
   // Step number
   doc_info.number()++;

  }

}//end of main



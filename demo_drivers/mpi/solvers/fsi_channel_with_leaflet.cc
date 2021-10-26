//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
 OneDLagrangianMesh<FSIHermiteBeamElement>* solid_mesh_pt() 
  {
   return Solid_mesh_pt;
  } 

 /// Access function to fluid mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* fluid_mesh_pt()
  {
   return Fluid_mesh_pt;
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


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
 OneDLagrangianMesh<FSIHermiteBeamElement>* Solid_mesh_pt;

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
 Solid_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (n_wall_el,hleaflet,undeformed_wall_pt,wall_time_stepper_pt);


 // Provide GeomObject representation of leaflet mesh and build fluid mesh
 //-----------------------------------------------------------------------

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 MeshAsGeomObject* wall_geom_object_pt=
  new MeshAsGeomObject(Solid_mesh_pt); 

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
 add_sub_mesh(Solid_mesh_pt);

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
 solid_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
 solid_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
 
 // Infinite slope: Pin type 1 (slope) dof for displacement direction 0 
 solid_mesh_pt()->boundary_node_pt(b,0)->pin_position(1,0);
 


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
 n_element = solid_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(solid_mesh_pt()->element_pt(e));
    
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
  (this,4,Fluid_mesh_pt,Solid_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Solid_mesh_pt,face); 
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
}





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
  (this,4,Fluid_mesh_pt,Solid_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Solid_mesh_pt,face); 
  
} // end_of_actions_after_adapt





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
 sprintf(filename,"%s/fsi_fluid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output wall solution 
 sprintf(filename,"%s/fsi_wall_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

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
 nel=Solid_mesh_pt->nelement();
 Vector<double> s(1);
 Vector<double> x(2);
 Vector<double> xi(1);
 Vector<double> N(2);
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   FSIHermiteBeamElement* el_pt=
    dynamic_cast<FSIHermiteBeamElement*>(Solid_mesh_pt->element_pt(e));

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
 

//=start_of_simple_fsi_preconditioner========================================== 
/// Simple FSI preconditioner. A block uppper triangular preconditioner
/// for the 2x2 FSI block system -- DOFs are decomposed into fluid DOFs and 
/// solid DOFs. The fluid subisidiary system is solved with the 
/// Navier Stokes Preconditioner and the solid subsidiary system with the
//=============================================================================
 template<typename MATRIX> 
class SimpleFSIPreconditioner 
 : public virtual BlockPreconditioner<MATRIX>
{
 
public :
 
 /// Constructor for SimpleFSIPreconditioner
 SimpleFSIPreconditioner(Problem* problem_pt)
  : BlockPreconditioner<MATRIX>(), Navier_stokes_preconditioner_pt(0),
    Solid_preconditioner_pt(0), Fluid_solid_coupling_matvec_pt(0),
    Navier_stokes_mesh_pt(0), Solid_mesh_pt(0)
  {
   // Create the Navier Stokes Schur Complement preconditioner
   Navier_stokes_preconditioner_pt = 
    new NavierStokesSchurComplementPreconditioner(problem_pt);

   // Create the Solid preconditioner
   Solid_preconditioner_pt = new SuperLUPreconditioner;

   // Create the matrix-vector product operator
   Fluid_solid_coupling_matvec_pt = new MatrixVectorProduct;

  }// end_of_constructor
 
 
 /// Destructor: Clean up.
 ~SimpleFSIPreconditioner()
  {
   //Delete the Navier-Stokes preconditioner
   delete Navier_stokes_preconditioner_pt; Navier_stokes_preconditioner_pt = 0;
   
   //Delete the solid preconditioner
   delete Solid_preconditioner_pt; Solid_preconditioner_pt = 0;
   
   // Delete the matrix vector product operator
   delete Fluid_solid_coupling_matvec_pt; Fluid_solid_coupling_matvec_pt = 0;
  }
 
 /// Broken copy constructor
 SimpleFSIPreconditioner(const SimpleFSIPreconditioner&)
  {
   BrokenCopy::broken_copy("SimpleFSIPreconditioner");
  }
 
  
 /// Access function to mesh containing the block-preconditionable
 /// Navier-Stokes elements. 
 void set_navier_stokes_mesh(Mesh* mesh_pt) 
  {
   Navier_stokes_mesh_pt = mesh_pt;
  }

 /// Access function to mesh containing the block-preconditionable
 /// FSI solid elements. 
 void set_solid_mesh(Mesh* mesh_pt) 
  {
   Solid_mesh_pt = mesh_pt;
  }

 /// Setup the preconditioner
 void setup();
 
 /// Apply preconditioner to r
 void preconditioner_solve(const DoubleVector &r,
                           DoubleVector &z);

private:

 /// Pointer the Navier Stokes preconditioner.
 NavierStokesSchurComplementPreconditioner* Navier_stokes_preconditioner_pt;

 /// Pointer to the solid preconditioner.
 Preconditioner* Solid_preconditioner_pt;

 /// Pointer to the fluid onto solid matrix vector product.
 MatrixVectorProduct* Fluid_solid_coupling_matvec_pt;

 /// Pointer to the navier stokes mesh.
 Mesh* Navier_stokes_mesh_pt;

 /// Pointer to the solid mesh.
 Mesh* Solid_mesh_pt;

};



/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////
// FSI preconditioner member functions
/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////




//=start_of_setup===============================================================
/// Setup the preconditioner.
//==============================================================================
 template<typename MATRIX> 
 void SimpleFSIPreconditioner<MATRIX>::setup()
{
 // setup the meshes for BlockPreconditioner and get the number of types of
 // DOF assoicated with each Mesh.
 // Mesh 0 is the fluid mesh, and hence DOFs 0 to n_fluid_dof_type-1 
 // are the fluid DOFs. Mesh 1 is the solid mesh and therefore DOFs 
 // n_fluid_dof_type to n_total_dof_type-1 are solid DOFs
 // set the mesh pointers
 this->set_nmesh(2);
 this->set_mesh(0,Navier_stokes_mesh_pt);
 this->set_mesh(1,Solid_mesh_pt);
 
 unsigned n_fluid_dof_type = this->ndof_types_in_mesh(0);
 unsigned n_total_dof_type = n_fluid_dof_type + this->ndof_types_in_mesh(1);
 
 // This fsi preconditioner has two types of block -- fluid and solid.
 // Create a map from DOF number to block type. The fluid block is labelled
 // 0 and the solid block 1.
 Vector<unsigned> dof_to_block_map(n_total_dof_type,0);
 for (unsigned i = n_fluid_dof_type; i < n_total_dof_type; i++)
  {
   dof_to_block_map[i] = 1;
  }
 
 // Call the BlockPreconditioner method block_setup(...) to assemble the data
 // structures required for block preconditioning.
 this->block_setup(dof_to_block_map);
 
 // First the solid preconditioner
 //===============================
 
 // get the solid block matrix (1,1)
 CRDoubleMatrix* solid_matrix_pt = new CRDoubleMatrix;
 this->get_block(1,1,*solid_matrix_pt);
 
 // setup the solid preconditioner
 // (perform the LU decomposition)
 Solid_preconditioner_pt->setup(solid_matrix_pt);
 delete solid_matrix_pt; solid_matrix_pt = 0;
 
 // Next the fluid preconditioner
 //==============================
 
 // Specify the relationship between the enumeration of DOF types in the 
 // master preconditioner and the Schur complement subsidiary preconditioner 
 // so that ns_dof_type[i_nst] contains i_master
 Vector<unsigned> ns_dof_list(n_fluid_dof_type);
 for (unsigned i = 0; i < n_fluid_dof_type; i++)
  {
   ns_dof_list[i] = i;
  }
 
 // Turn the NavierStokesSchurComplement preconditioner into a subsidiary 
 // preconditioner of this (FSI) preconditioner
 Navier_stokes_preconditioner_pt->
  turn_into_subsidiary_block_preconditioner(this,ns_dof_list);
 
 // Set up the NavierStokesSchurComplement preconditioner. 
 // (Pass it a pointer to the Navier Stokes mesh)
 Navier_stokes_preconditioner_pt->
  set_navier_stokes_mesh(Navier_stokes_mesh_pt);

 // Navier Stokes preconditioner is a subsidiary block preconditioner.
 // It therefore needs a pointer to the full matrix.
 Navier_stokes_preconditioner_pt->setup(this->matrix_pt());
 
 // Finally the fluid onto solid matrix vector product operator
 //============================================================
 
 // Similar to the solid preconditioner get the matrix
 CRDoubleMatrix* fluid_onto_solid_matrix_pt = new CRDoubleMatrix;
 this->get_block(1,0,*fluid_onto_solid_matrix_pt);
 
 // And setup the matrix vector product operator
 this->setup_matrix_vector_product(Fluid_solid_coupling_matvec_pt,
                                   fluid_onto_solid_matrix_pt,
                                   0);
 // Clean up
 delete fluid_onto_solid_matrix_pt; fluid_onto_solid_matrix_pt = 0;
}


//=start_of_preconditioner_solve================================================
/// Apply preconditioner.
//==============================================================================
 template<typename MATRIX> 
 void SimpleFSIPreconditioner<MATRIX>::preconditioner_solve(
  const DoubleVector &y, DoubleVector &z)
{
 // Fluid Subsidiary Preconditioner
 //=================================

 // Start by applying the Fluid subsidiary preconditioner
 // The fluid subsidiary preconditioner is a block preconditioner and
 // hence we pass it the global residual and solution vectors (y and z)
 Navier_stokes_preconditioner_pt->preconditioner_solve(y,z);

 // Fluid Onto Solid Matrix Vector Product Operator
 //================================================

 // The vector z_f contains the result of the action of the 
 // NavierStokesPreconditioner on a subset of  the elements of z.
 // Remember the fluid block index is 0 and the solid block index is 1.
 DoubleVector z_f;
 this->get_block_vector(0,z,z_f);

 // Apply the matrix vector product to z_f and store the results in w
 DoubleVector w;
 Fluid_solid_coupling_matvec_pt->multiply(z_f,w);

 // The vector y_s contains the solid residuals
 DoubleVector y_s;
 this->get_block_vector(1,y,y_s);

 // Subtract the action of the fluid onto solid matrix vector product from y_s
 y_s -= w;
 w = y_s; 

 // Solid Subsidiary Preconditioner
 //================================

 // Apply the solid preconditioner to s and return the result to the 
 // global solution vector z
 DoubleVector z_s;
 Solid_preconditioner_pt->preconditioner_solve(w,z_s);
 this->return_block_vector(1,z_s,z);
}



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////


//======= start_of_main================================================
/// Driver code  -- pass a command line argument if you want to run
/// the code in validation mode where it only performs a few steps
//=====================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

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

 // Create the solver. 
#ifdef OOMPH_HAS_TRILINOS
 TrilinosAztecOOSolver* solver_pt = new TrilinosAztecOOSolver;
 solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
#else
 GMRES<CRDoubleMatrix>* solver_pt = new GMRES<CRDoubleMatrix>;
#endif

 // Pass the solver to the problem.
 problem.linear_solver_pt() = solver_pt;
 
 // Create the preconditioner
 SimpleFSIPreconditioner<CRDoubleMatrix>* preconditioner_pt 
  = new SimpleFSIPreconditioner<CRDoubleMatrix>(&problem);

 // Pass the meshes to the preconditioner.
 preconditioner_pt->set_navier_stokes_mesh(problem.fluid_mesh_pt());
 preconditioner_pt->set_solid_mesh(problem.solid_mesh_pt());

 // Pass the preconditioner to the solver
 solver_pt->preconditioner_pt() = preconditioner_pt;

 // Set up doc info
 DocInfo doc_info; 
 doc_info.set_directory("RESLT");

 // Set max. number of adaptations
 unsigned max_adapt=3;

 // solve and document
 problem.steady_newton_solve(max_adapt);
 problem.doc_solution(doc_info);

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
}//end of main



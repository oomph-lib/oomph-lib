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

// Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Namespace for global parameters
//========================================================================
namespace Global_Parameters
{
 /// The number of nodes in one direction (default=2)
 unsigned Nnode_1d=2;
 
 /// The minimum level of uniform refinement
 unsigned Min_refinement_level=2;
 
 /// The additional levels of uniform refinement 
 unsigned Add_refinement_level=1;
 
 /// The number of adaptations allowed by the Newton solver
 unsigned N_adaptations=1;

 /// The choice of whether or not to use adaptation
 ///    0 = Uniform refinement
 ///    1 = Adaptive refinement
 unsigned Use_adaptation_flag=0;
 
 /// The choice of pre-smoother:
 ///    0 = Damped Jacobi
 ///    1 = Gauss-Seidel
 unsigned Pre_smoother_flag=0;
 
 /// The choice of post-smoother:
 ///    0 = Damped Jacobi
 ///    1 = Gauss-Seidel
 unsigned Post_smoother_flag=0;

 /// The choice of linear solver
 ///    0 = SuperLU
 ///    1 = Multigrid
 unsigned Linear_solver_flag=1;

 /// The MG solver allows for five different levels of output:
 ///    0 = Outputs everything
 ///    1 = Outputs everything and plots refinement and unrefinement patterns
 ///    2 = Outputs everything except the smoother timings 
 ///    3 = Outputs setup information but no V-cycle timings
 ///    4 = Suppresses all output
 /// Note: choosing '1' will also keep the coarser problems alive
 unsigned Output_management_flag=0;

 /// Variable used to decide whether or not convergence information
 /// is displayed:
 ///    0 = Don't display convergence information
 ///    1 = Display convergence information
 unsigned Doc_convergence_flag=0;
 
 /// DocInfo object used for documentation of the solution
 DocInfo Doc_info; 
} // End of namespace

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Returns a pointer to a smoother of the appropriate type
//========================================================================
namespace Smoother_Factory_Function_Helper
{
 /// Returns a pointer to a Smoother object which is to be used as
 /// the pre-smoother
 Smoother* set_pre_smoother()
 {
  // Set the pre-smoother
  if (Global_Parameters::Pre_smoother_flag==0)
  {
   // Create a new DampedJacobi object
   return new DampedJacobi<CRDoubleMatrix>;
  }
  else if (Global_Parameters::Pre_smoother_flag==1)
  {
   // Create a new GS object
   return new GS<CRDoubleMatrix>;
  }
  else
  {
   // The only choices for the parameter Pre_smoother_flag are 0 and 1 so
   // notify the user if anything else is being used
   std::string err_strng="Incorrect parameter value; presmoother can ";
   err_strng+="only be 0 or 1.";
   throw OomphLibError(err_strng,
		       OOMPH_CURRENT_FUNCTION,
		       OOMPH_EXCEPTION_LOCATION);
  }
 } // End of set_pre_smoother
 
 /// Returns a pointer to a Smoother object which is to be used as
 /// the post-smoother
 Smoother* set_post_smoother()
 {
  // Set the post-smoother
  if (Global_Parameters::Post_smoother_flag==0)
  {
   // Create a new DampedJacobi object
   return new DampedJacobi<CRDoubleMatrix>;
  }
  else if (Global_Parameters::Post_smoother_flag==1)
  {
   // Create a new GS object
   return new GS<CRDoubleMatrix>;
  }
  else
  {
   // The only choices for the parameter Post_smoother_flag are 0 and 1 so
   // notify the user if anything else is being used
   std::string err_strng="Incorrect parameter value; postsmoother can ";
   err_strng+="only be 0 or 1.";
   throw OomphLibError(err_strng,
		       OOMPH_CURRENT_FUNCTION,
		       OOMPH_EXCEPTION_LOCATION);
  }
 } // End of set_post_smoother
} // End of Smoother_Factory_Function_Helper

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//===== start_of_namespace================================================
/// Namespace for exact solution for Poisson equation
//========================================================================
namespace SolnForPoisson
{
 /// Parameter for steepness of "step"
 double Alpha=4.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }

 /// Flux required by the exact solution on a boundary on which x is fixed
 void prescribed_flux_on_fixed_x_boundary(const Vector<double>& x, 
                                          double& flux)
 {
  // The outer unit normal to the boundary is (1,0)
  double N[2] = {1.0,0.0};
  
  // The flux in terms of the normal is
  flux = -(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*TanPhi*N[0]+(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*N[1];
 } 
} // end of namespace

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//=================================start_of_problem_class=================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. Flux boundary conditions are applied
/// along boundary 1 (the boundary where x=L). The specific type of 
/// element and mesh is specified via the template parameter.
//========================================================================
template<class ELEMENT,class MESH> 
class FluxPoissonMGProblem : public MGProblem
{

public:

 /// Constructor: Pass pointer to source function
 FluxPoissonMGProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor
 ~FluxPoissonMGProblem();

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution();
  
 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();

 /// Pointer to the bulk mesh. Overloads the pure virtual function in
 /// the abstract base class MGProblem. Must be refineable to allow the
 /// use of refine_base_mesh_as_in_reference_mesh() in make_copy()
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
  {
   // Return the pointer to the bulk mesh
   return Bulk_mesh_pt;
  }

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
  {   
   // Document the solution
   doc_solution();
  }
  
 /// Build and set multgrid solver
 void set_multigrid_solver();

private:
 
 /// Access function for the mesh
 MESH* mesh_pt() 
  {
   return dynamic_cast<MESH*>(Problem::mesh_pt());
  }
  
 /// Return a pointer to a new instance of the same problem.
 MGProblem* make_new_problem()
  {
   // Make new problem of the FluxPoissonMGProblem class whose template
   // parameters are specified by the template parameters of the current
   // problem
   return new FluxPoissonMGProblem<ELEMENT,MESH>
    (&SolnForPoisson::source_function);
  }

 /// Create Poisson flux elements on the b-th boundary of the 
 /// problem's mesh
 void create_flux_elements(const unsigned &b,
			   MESH* const &bulk_mesh_pt,
			   Mesh* const &surface_mesh_pt);

 /// Delete Poisson flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);
  
 /// Pointer to the "bulk" mesh
 MESH* Bulk_mesh_pt;
  
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;
}; // End of FluxPoissonMGProblem class

//=======set_multigrid_solver=============================================
/// Build and set multigrid solver
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::set_multigrid_solver()
{
 // Make an object of the MGSolver class and get the pointer to it
 MGSolver<2>* mg_solver_pt=new MGSolver<2>(this);
 
 // Switch solver to MG
 linear_solver_pt()=mg_solver_pt;
  
 // Set the pre-smoother factory function
 mg_solver_pt->set_pre_smoother_factory_function
  (Smoother_Factory_Function_Helper::set_pre_smoother);

 // Set the post-smoother factory function
 mg_solver_pt->set_post_smoother_factory_function
  (Smoother_Factory_Function_Helper::set_post_smoother);
  
 // Suppress certain timings
 if (Global_Parameters::Output_management_flag==1)
 {
  mg_solver_pt->enable_doc_everything();
 }
 else if (Global_Parameters::Output_management_flag==2)
 {
  mg_solver_pt->disable_doc_time();
 }
 else if (Global_Parameters::Output_management_flag==3)
 {
  mg_solver_pt->disable_v_cycle_output();
 }
 else if (Global_Parameters::Output_management_flag==4)
 {
  mg_solver_pt->disable_output();
 }

 // If the user wishes to document the convergence information
 if (Global_Parameters::Doc_convergence_flag)
 {
  // Create a file to record the convergence history
  mg_solver_pt->open_convergence_history_file_stream("RESLT/conv.dat");
 }
} // End of set_multigrid_solver

//=======start_of_constructor=============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT,class MESH>
FluxPoissonMGProblem<ELEMENT,MESH>::
FluxPoissonMGProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{
 // If we're choosing to use multigrid as our linear solver
 if (Global_Parameters::Linear_solver_flag==1)
 {
  // Call the set_multigrid_solver function
  set_multigrid_solver();
 }

 // Build the mesh:
 //---------------- 
 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build "bulk" mesh
 Bulk_mesh_pt=new MESH(n_x,n_y,l_x,l_y);
 
 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Adjust error targets for adaptive refinement
 Bulk_mesh_pt->max_permitted_error()=0.001;
 Bulk_mesh_pt->min_permitted_error()=0.0001;
  
 // Create "surface mesh" that will contain only the prescribed-flux 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt=new Mesh;
  
 // Set boundary conditions:
 //-----------------------------------
 // Create prescribed-flux elements from all elements that are adjacent to
 // boundary 1, but add them to a separate mesh. Note that this is exactly
 // the same function as used in the single mesh version of the problem,
 // we merely pass different Mesh pointers
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);
 
 // Set the boundary conditions for this problem: All nodes are free by
 // default -- just pin the ones that have Dirichlet conditions here
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
 {
  // Leave nodes on boundary 1 free
  if(b!=1)
  {
   // Get the number of nodes on boundary b
   unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

   // Loop over the boundary nodes
   for (unsigned n=0;n<n_node;n++)
   {
    // Pin the boundary nodes
    Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
   }
  } // if(b!=1)
 } // for(unsigned b=0;b<n_bound;b++)

 // Build the global mesh from its submeshes:
 //------------------------------------------
 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
  
 // Rebuild the Problem's global mesh from its various sub-meshes
 build_global_mesh();
  
 // Complete the build of all elements so they are fully functional:
 //-----------------------------------------------------------------
 // Loop over the Poisson bulk elements to set up element-specific things
 // that cannot be handled by constructor: Pass pointer to source function
 for(unsigned e=0;e<Bulk_mesh_pt->nelement();e++)
 {
  // Upcast from GeneralisedElement to Poisson bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // Set the source function pointer
  el_pt->source_fct_pt()=Source_fct_pt;
 }
  
 // Setup equation numbering scheme
 assign_eqn_numbers(); 
} // end of constructor


//=======start_of_destructor==============================================
/// Destructor for Poisson problem.
//========================================================================
template<class ELEMENT,class MESH>
FluxPoissonMGProblem<ELEMENT,MESH>::~FluxPoissonMGProblem()
{
 // If we're using multigrid as the linear solver
 if (Global_Parameters::Linear_solver_flag==1)
 {
  // Delete the MG solver pointers
  delete linear_solver_pt();

  // Set the pointer to null
  linear_solver_pt()=0;
 }
   
 // Delete the error estimator (allocated as a 'new' object)
 delete Bulk_mesh_pt->spatial_error_estimator_pt();
 
 // Set the pointer to null
 Bulk_mesh_pt->spatial_error_estimator_pt()=0;

 // Delete the "bulk" mesh
 delete Bulk_mesh_pt;

 // Set the pointer to null
 Bulk_mesh_pt=0;
 
 // Delete the dynamically allocated surface 
 delete Surface_mesh_pt;
 
 // Set the pointer to null
 Surface_mesh_pt=0;  
} // end of destructor

//============start_of_create_flux_elements===============================
/// Create Poisson Flux Elements on the b-th boundary of the Mesh
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::
create_flux_elements(const unsigned& b, MESH* const& bulk_mesh_pt,
		     Mesh* const& surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_belement=bulk_mesh_pt->nboundary_element(b);
  
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_belement;e++)
 {
  // Get pointer to the bulk element that is adjacent to boundary b
  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
   (bulk_mesh_pt->boundary_element_pt(b,e));
    
  // Find the index of the face of the bulk element at the boundary
  int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
    
  // Build the corresponding prescribed-flux element
  PoissonFluxElement<ELEMENT>* flux_element_pt = new 
   PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);
    
  // Add the prescribed-flux element to the mesh
  surface_mesh_pt->add_element_pt(flux_element_pt);  
 } // end of loop over bulk elements adjacent to boundary b
 
 // Loop over the flux elements to pass pointer to prescribed flux function
 unsigned n_selement=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_selement;e++)
 {
  // Upcast from GeneralisedElement to Poisson flux element
  PoissonFluxElement<ELEMENT> *el_pt = 
   dynamic_cast< PoissonFluxElement<ELEMENT>*>(Surface_mesh_pt->element_pt(e));
   
  // Set the pointer to the prescribed flux function
  el_pt->flux_fct_pt()=&SolnForPoisson::prescribed_flux_on_fixed_x_boundary;
 } 
} // End of create_flux_elements


//============start_of_delete_flux_elements===============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
{
 // Find out how many surface elements there are in the surface mesh
 unsigned n_element=surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
 {
  // Kill surface element
  delete surface_mesh_pt->element_pt(e);
 }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();
} // End of delete_flux_elements

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_before_adapt()
{ 
 // Kill the flux elements and wipe surface mesh
 delete_flux_elements(Surface_mesh_pt);

 // Flush the submeshes
 flush_sub_meshes();

 // Add the bulk mesh back in
 add_sub_mesh(Bulk_mesh_pt);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
} // End of actions_before_adapt

//================================start_of_actions_after_adapt============
/// Actions after adapt: Rebuild the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1 and add them to surfac mesh
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the flux mesh back into the global mesh
 add_sub_mesh(Surface_mesh_pt);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();   
} // End of actions_after_adapt

//========================start_of_actions_before_newton_solve============
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 
 // Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
 {
  // Only update Dirichlet nodes
  if (i!=1)
  {
   // How many nodes are there on this boundary?
   unsigned n_node=Bulk_mesh_pt->nboundary_node(i);
     
   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
   {
    // Get pointer to node
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(i,n);
       
    // Extract nodal coordinates from node:
    Vector<double> x(2);
    x[0]=nod_pt->x(0);
    x[1]=nod_pt->x(1);
       
    // Compute the value of the exact solution at the nodal point
    Vector<double> u(1);
    SolnForPoisson::get_exact_u(x,u);
       
    // Assign the value to the one (and only) nodal value at this node
    nod_pt->set_value(0,u[0]);
   }
  } // if (i!=1)
 } // for(unsigned i=0;i<n_bound;i++)
} // End of actions_before_newton_solve

//==============================================start_of_doc==============
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT,class MESH>
void FluxPoissonMGProblem<ELEMENT,MESH>::doc_solution()
{
 // Output file stream
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",
	 Global_Parameters::Doc_info.directory().c_str(),
	 Global_Parameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",
	 Global_Parameters::Doc_info.directory().c_str(),
	 Global_Parameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,SolnForPoisson::get_exact_u); 
 some_file.close();
 
 // Doc error and return of the square of the L2 error. Note:
 //----------------------------------------------------------
 // this is only calculated over the bulk mesh since we cannot
 //-----------------------------------------------------------
 // properly deal with two different types of elements
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",
	 Global_Parameters::Doc_info.directory().c_str(),
	 Global_Parameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,SolnForPoisson::get_exact_u,error,norm); 
 some_file.close();

 // Document the L2 norm of the error and the L2 norm of the solution
 oomph_info << "Norm of error   : " << sqrt(error)
	    << "\nNorm of solution: " << sqrt(norm) << std::endl;
 
 // Increment the documentation counter
 Global_Parameters::Doc_info.number()++;
} // End of doc_solution

 

//===============================================start_of_main============
/// Demonstrate how to solve 2D Poisson problem with flux boundary 
/// conditions
//========================================================================
int main(int argc, char **argv)
{ 
 //------------------------
 // Command line arguments
 //------------------------
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Choose the number of nodes in one direction of an element;
 // Values of nnode_1d:
 //        2: Bilinear interpolation 
 //        3: Biquadratic interpolation 
 //        4: Bicubic interpolation 
 CommandLineArgs::specify_command_line_flag(
  "--nnode_1d",&Global_Parameters::Nnode_1d);

 // Choose the minimum level of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--min_ref",&Global_Parameters::Min_refinement_level);
 
 // Choose the additional levels of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--add_ref",&Global_Parameters::Add_refinement_level);
 
 // Choose the maximum number of adaptive refinements
 CommandLineArgs::specify_command_line_flag(
  "--n_adapt",&Global_Parameters::N_adaptations);
 
 // Choose how many additional levels of uniform refinement to use
 CommandLineArgs::specify_command_line_flag(
  "--use_adapt",&Global_Parameters::Use_adaptation_flag);
 
 // Choose the pre-smoother
 CommandLineArgs::specify_command_line_flag(
  "--presmoother",&Global_Parameters::Pre_smoother_flag);
 
 // Choose the post-smoother
 CommandLineArgs::specify_command_line_flag(
  "--postsmoother",&Global_Parameters::Post_smoother_flag);
 
 // Choose the linear solver
 CommandLineArgs::specify_command_line_flag(
  "--linear_solver",&Global_Parameters::Linear_solver_flag);
 
 // Decide whether or not to suppress all or some of the MG solver output
 CommandLineArgs::specify_command_line_flag(
  "--output_flag",&Global_Parameters::Output_management_flag);
     
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--conv_flag",&Global_Parameters::Doc_convergence_flag);
  
 // Parse command line
 CommandLineArgs::parse_and_assign();
  
 // Document what has been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 //--------------------------------
 // Set the documentation directory
 //--------------------------------
 // Set output directory
 Global_Parameters::Doc_info.set_directory("RESLT");
 
 //--------------------
 // Set up the problem
 //--------------------
 // Initialise a null pointer to the class Problem 
 Problem* problem_pt=0;
 
 // Set the problem pointer depending on the input (defaulted to nnode_1d=2)
 if (Global_Parameters::Nnode_1d==2)
 {
  // Using linear interpolation
  typedef RefineableQPoissonElement<2,2> ELEMENT;
  
  // Typedef the mesh and template it by the chosen element type
  typedef RefineableRectangularQuadMesh<ELEMENT> MESH;

  // Set the problem pointer
  problem_pt=new FluxPoissonMGProblem<ELEMENT,MESH>
   (&SolnForPoisson::source_function);
 }
 else if (Global_Parameters::Nnode_1d==3)
 {
  // Using quadratic interpolation
  typedef RefineableQPoissonElement<2,3> ELEMENT;
  
  // Typedef the mesh and template it by the chosen element type
  typedef RefineableRectangularQuadMesh<ELEMENT> MESH;

  // Set the problem pointer
  problem_pt=new FluxPoissonMGProblem<ELEMENT,MESH>
   (&SolnForPoisson::source_function);
 }  
 else if (Global_Parameters::Nnode_1d==4)
 {
  // Using cubic interpolation
  typedef RefineableQPoissonElement<2,4> ELEMENT;
  
  // Typedef the mesh and template it by the chosen element type
  typedef RefineableRectangularQuadMesh<ELEMENT> MESH;

  // Set the problem pointer
  problem_pt=new FluxPoissonMGProblem<ELEMENT,MESH>
   (&SolnForPoisson::source_function);
 }
 else
 {
  // Throw an error otherwise
  throw OomphLibError("nnode_1d can only be 2,3 or 4.",
		      OOMPH_CURRENT_FUNCTION,
		      OOMPH_EXCEPTION_LOCATION);
 }

 //-------------------
 // Solve the problem!
 //-------------------
 // If the user wishes to use adaptive refinement then we use the
 // Newton solver with a given argument to indicate how many
 // adaptations can be done
 if (Global_Parameters::Use_adaptation_flag)
 {
  // If the user did not choose to suppress everything
  if (Global_Parameters::Output_management_flag!=4)
  {
   oomph_info << "\n====================Initial Refinement===================="
	      << std::endl;
   
#ifndef PARANOID
   // If PARANOID isn't enabled we need an extra line
   oomph_info << std::endl;
#endif
  }

  // Refine once at least otherwise the V-cycle only uses SuperLU and converges
  // instantly
  problem_pt->refine_uniformly();
   
  // Solve the problem
  problem_pt->newton_solve();
   
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<Global_Parameters::N_adaptations;i++)
  { 
   // Adapt the problem
   problem_pt->adapt();
   
   // Solve the problem
   problem_pt->newton_solve();
  }
 }
 // If the user instead wishes to use uniform refinement
 else
 {
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<Global_Parameters::Min_refinement_level;i++)
  { 
   // If the user did not choose to suppress everything
   if (Global_Parameters::Output_management_flag!=4)
   {
    oomph_info << "\n====================Initial Refinement====================\n"
	       << std::endl;
   }

   // Add additional refinement
   problem_pt->refine_uniformly();
  }
 
  // Solve the problem
  problem_pt->newton_solve();
 
  // Refine and solve until the additional refinements have been completed
  for (unsigned i=0;i<Global_Parameters::Add_refinement_level;i++)
  {
   // If the user did not choose to suppress everything
   if (Global_Parameters::Output_management_flag!=4)
   {
    oomph_info << "===================Additional Refinement==================\n"
	       << std::endl;
   }
 
   // Add additional refinement
   problem_pt->refine_uniformly();
  
   // Solve the problem
   problem_pt->newton_solve();
  }
 } // if (Global_Parameters::Use_adaptation_flag)

 // Delete the problem pointer
 delete problem_pt;

 // Make it a null pointer
 problem_pt=0;
} // End of main

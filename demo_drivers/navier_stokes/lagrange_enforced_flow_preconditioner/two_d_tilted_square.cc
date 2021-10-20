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
//Driver for 2D tilted rectangle

#include <fenv.h>
#include <sstream>
#include <iomanip>
#include <ios>


// Generic includes
#include "generic.h"
#include "navier_stokes.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"


#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;

using namespace oomph;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{

 /// Problem dimension
 static const unsigned Dim = 2;

 /// Min and max x value respectively.
 static const double X_min = 0.0;
 static const double X_max = 1.0;

 /// Min and max y value respectively.
 static const double Y_min = 0.0;
 static const double Y_max = 1.0;

 /// The domain length in the x and y direction respectively.
 static const double Lx = X_max - X_min;
 static const double Ly = Y_max - Y_min;

 /// The viscous term
 unsigned Visc = 0;

 /// Reynolds number
 double Re = 100.0;

 /// Tilting angle of the domain with the x-axis
 double Ang_deg = 30.0;
 double Ang_rad = -1.0;

 /// Number of elements in 1D
 unsigned Noel = 4;

 /// Use trilinos?
 bool Use_trilinos = false;

 /// Use LSC preconditioner for the Navier-Stokes block?
 bool Use_lsc = false;

 /// Use Boomer AMG for the momentum block?
 bool Use_amg_for_f = false;

 /// Use Boomer AMG for the pressure block?
 bool Use_amg_for_p = false;

 /// Soln number (for doc_solution)
 unsigned Soln_num = 0;

 /// Convert degrees to radians
 inline double degtorad(const double& ang_deg)
 {
   return ang_deg * (MathematicalConstants::Pi / 180.0);
 }

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

} // namespace Global_Variables


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


namespace oomph
{
//========================================================================
///  A Sloping Mesh  class.
///
/// Derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi (in radians)
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly, const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
   {
    // Find out how many nodes there are
    const unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      const double x=nod_pt->x(0);
      const double y=nod_pt->x(1);

      // Set new nodal coordinates
      nod_pt->x(0)=x*cos(phi)-y*sin(phi);
      nod_pt->x(1)=x*sin(phi)+y*cos(phi);
     }
   }
 };
} // end of namespace oomph

//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

  ///  Constructor: Pass number of elements in x and y directions and
  /// lengths
  TiltedCavityProblem();

  ///  Update before Newton solve.
  void actions_before_newton_solve()
  {
    // Alias the namespace for convenience
    namespace GV = Global_Variables;
    GV::Iterations.clear();

    // Which is the inflow boundary?
    const unsigned in_bound = Left_bound;

    // The number of nodes on a boundary
    unsigned num_nod = mesh_pt()->nboundary_node(in_bound);
    for(unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=mesh_pt()->boundary_node_pt(in_bound,inod);

      // Pin both velocity components
      nod_pt->pin(0);
      nod_pt->pin(1);
    
      // Get the x and y cartesian coordinates
      double x0=nod_pt->x(0);
      double x1=nod_pt->x(1);

      // Tilt x1 by -Ang_rad, this will give us the original coordinate.
      double x1_old = x0*sin(-GV::Ang_rad) + x1*cos(-GV::Ang_rad);

      // Now calculate the parabolic inflow at this point.
      double u0_old = (x1_old - GV::Y_min)*(GV::Y_max - x1_old);
   
      // Now apply the rotation to u0_old, using rotation matrices.
      // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
      // velocity in the x direction only. There is no velocity
      // in the y direction.
      double u0=u0_old*cos(GV::Ang_rad);
      double u1=u0_old*sin(GV::Ang_rad);

      nod_pt->set_value(0,u0);
      nod_pt->set_value(1,u1);
    } // for  
  }

  ///  Update after Newton step - document the number of iterations 
  /// required for the iterative solver to converge.
  void actions_after_newton_step()
  {
    // Alias the namespace for convenience
    namespace GV = Global_Variables;

    unsigned iters = 0;
    // Get the iteration counts
#ifdef PARANOID
    IterativeLinearSolver* iterative_solver_pt
      = dynamic_cast<IterativeLinearSolver*>
      (this->linear_solver_pt());
    if(iterative_solver_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Cannot cast the solver pointer." << std::endl;

      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      iters = iterative_solver_pt->iterations();
      GV::Iterations.push_back(iters);
    }
#else
    iters = static_cast<IterativeLinearSolver*>
       (this->linear_solver_pt())->iterations();
    GV::Iterations.push_back(iters);
#endif
 }

 /// Doc the solution
 void doc_solution();

 ///  Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

private:


 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_P_pt;

 /// Preconditioner
 Preconditioner* Prec_pt;

 /// Preconditioner for the Navier-Stokes block
 Preconditioner* Navier_stokes_prec_pt;

 /// Preconditioner for the momentum block
 Preconditioner* F_preconditioner_pt;

 /// Preconditioner for the pressure block
 Preconditioner* P_preconditioner_pt;

 /// Iterative linear solver
 IterativeLinearSolver* Solver_pt;

 // Enumeration for the boundaries of the rectangular domain
 //         2
 //    -------------
 //    |           |
 //  3 |           | 1
 //    |           |
 //    |           |
 //    -------------
 //         0
 //
 
 /// Bottom boundary
 unsigned Bottom_bound; // 0

 /// Right boundary
 unsigned Right_bound;  // 1

 /// Top boundary
 unsigned Top_bound;    // 2

 /// Left boundary
 unsigned Left_bound;   // 3

};


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT>
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem()
{
  // Alias the namespace for convenience
  namespace GV = Global_Variables;

  // Assign the enumerations for the boundaries.
  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;
  
  // First we create the tilted cavity mesh.
  Bulk_mesh_pt = new SlopingQuadMesh<ELEMENT>(GV::Noel,GV::Noel,
                                              GV::Lx,GV::Ly,
                                              GV::Ang_rad); 

  // Set the boundary conditions, recall that the boundaries are:
  //
  //             2 non slip
  //         ----------
  //         |        |
  // 3 Inflow|        |1 P.O. (parallel outflow)
  //         |        |
  //         ----------
  //             0 non slip
  
  // Parallel outflow boundary is located at 1.
  const unsigned po_bound = Right_bound; 

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_P_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
  create_parall_outflow_lagrange_elements(po_bound,
                                          Bulk_mesh_pt,
                                          Surface_mesh_P_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_P_pt);

  // combine all sub-meshes into a single mesh.
  build_global_mesh();

  // All nodes are free by default
  // just pin the ones that have Dirichlet conditions
  const unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  { 
    if(ibound != po_bound)
    {
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
      {
        // Get node
        Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
 
        nod_pt->pin(0);
        nod_pt->pin(1);
        
        nod_pt->set_value(0,0);
        nod_pt->set_value(1,0);
      }
    }
  }
 
  //Complete the problem setup to make the elements fully functional

  //Loop over the elements
  unsigned n_el = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_el;e++)
  {
    //Cast to a fluid element
    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the Reynolds number, etc
    el_pt->re_pt() = &GV::Re;
  } // for(unsigned e=0;e<n_el;e++)

  //Assign equation numbers
  oomph_info << "\n equation numbers : "
             << assign_eqn_numbers() << std::endl;


  // We choose either trilinos or oomph-lib's GMRES as our linear solver
  if(GV::Use_trilinos)
  {
#ifdef OOMPH_HAS_TRILINOS
    // Create the trilinos solver.
    TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
    trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;

    // Store the solver pointer.
    Solver_pt = trilinos_solver_pt;
#endif
  }
  else
  {
    // Create oomph-lib iterative linear solver.
    IterativeLinearSolver* solver_pt = new GMRES<CRDoubleMatrix>;
    
    // We use RHS preconditioning. Note that by default,
    // left hand preconditioning is used.
    static_cast<GMRES<CRDoubleMatrix>*>(solver_pt)
      ->set_preconditioner_RHS();

    // Store the solver pointer.
    Solver_pt = solver_pt;
  }

  // Set up other solver parameters.
  Solver_pt->tolerance() = 1.0e-6;

  // Overview of solvers:
  //
  // The Jacobian takes the block form:
  // 
  // | F_ns | L^T |
  // |------------|
  // |   L  | 0   |
  //
  // where L correspond to the constrained block,
  // F_ns is the Navier-Stokes block with the following block structure
  //
  // |  F | B^T |
  // |----------|
  // |  B |  0  |
  //
  //
  // The lagrange enforced flow preconditioner takes the form:
  // F_ns + L^T inv(W) L |
  // --------------------------
  //                     | W
  //
  // where W=LL^T
  //
  // We use SuperLU to solve the W block (2,2)
  //
  // For the (1,1) block, we can use SuperLU or the LSC preconditioner.


  // Create the preconditioner
  LagrangeEnforcedFlowPreconditioner* lgr_prec_pt
    = new LagrangeEnforcedFlowPreconditioner;

  // Create the vector of mesh pointers!
  Vector<Mesh*> mesh_pt;
  mesh_pt.resize(2,0);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_P_pt;

  lgr_prec_pt->set_meshes(mesh_pt);
  
  NavierStokesSchurComplementPreconditioner* lsc_prec_pt = 0;
  if(GV::Use_lsc)
  {
    // Create the NS LSC preconditioner.
    lsc_prec_pt = new NavierStokesSchurComplementPreconditioner(this);
    lsc_prec_pt->set_navier_stokes_mesh(Bulk_mesh_pt);
    lsc_prec_pt->use_lsc();
    lgr_prec_pt->set_navier_stokes_preconditioner(lsc_prec_pt);

    if(GV::Use_amg_for_f)
    {
      if(GV::Visc == 0)
      {
        F_preconditioner_pt
          = Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper::
            boomer_amg_for_2D_momentum_simple_visc();
        lsc_prec_pt->set_f_preconditioner(F_preconditioner_pt);
      }
      else
      {
        F_preconditioner_pt
          = Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper::
            boomer_amg_for_2D_momentum_stressdiv_visc();
        lsc_prec_pt->set_f_preconditioner(F_preconditioner_pt);
      }
    }
 
    if(GV::Use_amg_for_p)
    {
      P_preconditioner_pt
        = Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper::
          boomer_amg_for_2D_poisson_problem();
      lsc_prec_pt->set_p_preconditioner(P_preconditioner_pt);
    }
  }
  else
  {
    lgr_prec_pt->set_superlu_for_navier_stokes_preconditioner();
  }

  // Store the preconditioner pointers.
  Navier_stokes_prec_pt = lsc_prec_pt;
  Prec_pt = lgr_prec_pt;

  // Pass the preconditioner to the solver.
  Solver_pt->preconditioner_pt() = lgr_prec_pt;

  // Pass the solver to the problem.
  this->linear_solver_pt() = Solver_pt;

  // Set the Newton solver tolerance.
  this->newton_solver_tolerance() = 1.0e-6;
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{    
  // Alias the namespace for convenience
  namespace GV = Global_Variables; 
 
  // Create the soln file name string.
  std::stringstream filename;
  filename <<"RESLT/soln"<<GV::Soln_num<<".dat";

  // Number of plot points
  unsigned npts=5;

  // Output solution
  std::ofstream some_file;
  some_file.open(filename.str().c_str());
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
}

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   // What is the index of the face of element e along boundary b?
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index,0);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//==start_of_main======================================================
/// Driver for Lagrange enforced flow preconditioner
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Alias the namespace for convenience.
  namespace GV = Global_Variables;

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // Store the viscous term, default is simple form.
  CommandLineArgs::specify_command_line_flag("--visc", &GV::Visc);

  // Default is 100.0
  CommandLineArgs::specify_command_line_flag("--re", &GV::Re);
 
  // Defaults to 30 degrees
  CommandLineArgs::specify_command_line_flag("--ang", &GV::Ang_deg);

  // Defaults to 4
  CommandLineArgs::specify_command_line_flag("--noel", &GV::Noel);

  // Use trilinos solver?
  CommandLineArgs::specify_command_line_flag("--use_trilinos");

  // Use the LSC preconditioner for the Navier-Stokes block?
  CommandLineArgs::specify_command_line_flag("--use_lsc");

  // Use the amg for the momentum block?
  CommandLineArgs::specify_command_line_flag("--use_amg_for_f");

  // Use the amg for the pressure block?
  CommandLineArgs::specify_command_line_flag("--use_amg_for_p");
  
  CommandLineArgs::specify_command_line_flag("--soln_num", &GV::Soln_num);

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  GV::Ang_rad = GV::degtorad(GV::Ang_deg);

  if(GV::Visc == 0)
  {
    oomph_info << "Doing simple form" << std::endl; 
    for (unsigned d = 0; d < GV::Dim; d++)
    {
      NavierStokesEquations<GV::Dim>::Gamma[d] = 0.0;
    }
  }
  else
  {
    oomph_info << "Doing stress divergence form" << std::endl; 
    for (unsigned d = 0; d < GV::Dim; d++)
    {
      NavierStokesEquations<GV::Dim>::Gamma[d] = 1.0;
    }
  }
 
  // Set the flag for trilinos solver
  if(CommandLineArgs::command_line_flag_has_been_set("--use_trilinos"))
  {
#ifdef OOMPH_HAS_TRILINOS
    GV::Use_trilinos = true;
#else
    std::ostringstream warning_stream;
    warning_stream << "WARNING: \n"
      << "No trilinos installed, using oomphlib's GMRES solver.\n";
    OomphLibWarning(warning_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);

    GV::Use_trilinos = false;
#endif
  }
  else
  {
    GV::Use_trilinos = false;
  }
 
  // Set the flag for lsc
  if(CommandLineArgs::command_line_flag_has_been_set("--use_lsc"))
  {
    GV::Use_lsc = true;
  }
  else
  {
    GV::Use_lsc = false;
  }

  // Set the flag for amg for the momentum block
  if(CommandLineArgs::command_line_flag_has_been_set("--use_amg_for_f"))
  {
    GV::Use_amg_for_f = true;
  }
  else
  {
    GV::Use_amg_for_f = false;
  }

  // Set the flag for amg for the pressure block
  if(CommandLineArgs::command_line_flag_has_been_set("--use_amg_for_p"))
  {
    GV::Use_amg_for_p = true;
  }
  else
  {
    GV::Use_amg_for_f = false;
  }

  TiltedCavityProblem< QTaylorHoodElement<GV::Dim> > problem;

  // Solve the problem
  problem.newton_solve();

  if(CommandLineArgs::command_line_flag_has_been_set("--soln_num"))
  {
    problem.doc_solution();

    // Create the soln file name string.
    std::stringstream filename;
    filename <<"RESLT/iter"<<GV::Soln_num<<".dat";
    std::ofstream some_file;
    some_file.open(filename.str().c_str());

    // Print out the iteration counts
    const unsigned num_newton_steps = GV::Iterations.size();
    some_file << num_newton_steps << "\n";
    for (unsigned stepi = 0; stepi < num_newton_steps; stepi++) 
    {
      some_file << GV::Iterations[stepi] << "\n";
    }
    some_file << "\n";
    some_file.close();

  }

  
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif  
 
} // end_of_main



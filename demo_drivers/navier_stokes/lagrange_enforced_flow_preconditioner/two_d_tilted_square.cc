//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
//LIC//=====================================================================
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

#ifdef OOMPH_HAS_HYPRE

//==========================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//==========================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}

#endif

 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
///
/// Namepace for all things related to the problem where the domain is a 
/// unit square.
///
////////////////////////////////////////////////////////////////////////////
namespace SquareLagrange
{
  // Prob id, set by main method
  const int* Prob_id_pt = 0;

  std::string Prob_str = "";
  std::string Ang_deg_str = "";
  std::string Noel_str = "";


  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  // Min and max x value respectively.
  static const double X_min = 0.0;
  static const double X_max = 1.0;

  // Min and max y value respectively.
  static const double Y_min = 0.0;
  static const double Y_max = 1.0;

  // The length in the x and y direction respectively.
  static const double Lx = X_max - X_min;
  static const double Ly = Y_max - Y_min;

////////////////////////////////////////////////////////////////////////////

  // CL - set directly from the command line.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.
  //
  // Problem parameter overview:
  //
  // // Solvers:
  //
  //
  // F_ns + L^T inv(W) L | B^T
  // --------------------------
  //                     | W
  //
  // W = 0 (SuperLU)
  // NS_solver = 0 (SuperLU) or 1 (LSC)
  // 
  // If NS_solver = 1, then we have:
  //
  // | F | B^T |   |
  // |----------   |
  // |   |-M_s |   |
  // |-------------|
  // |         | W |
  //
  // F_solver = 0 (SuperLU) or 1 (AMG)
  // P_solver = 0 (SuperLU) or 1 (AMG)
  // 

  // These are self explanatory:
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang = 0.0; //CL, Angle in degrees
  unsigned Noel = 4; //CL, Number of elements in 1D
  // the default is the norm of the momentum block.

  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

    CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

  inline void set_ang_str()
  {
    if(Prob_id_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Oh dear, Prob_id_pt is null. Please set this in main().\n"
        << "This should be stored in NSPP::Prob_id, and set by cmd via\n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // If this is the vanilla problem, we set the angle as -1 and set the
    // string as "A_". This would indicate that no angle is used.
    // Furthermore, we ensure that no --ang is set.
    if(Prob_str.compare("SqVa") == 0)
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "prob_id is 88, doing vanilla LSC with no tilt.\n"
          << "But you have set --ang, please do not set this."; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      Ang = -1.0;
      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A_";
      Ang_deg_str = strs.str();
    }
    else
    // This problem requires tilting, thus we set the Ang and Ang_deg_str.
    {
      // But first we check that --ang has been set.
      // Check that Ang has been set.
      if(!CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "Angle has not been set. Set (in degrees) with: \n"
          << "--ang \n"; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Now we need to convert Ang_deg into radians.
      Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A" << Ang_deg;
      Ang_deg_str = strs.str();
    }
  } // set_ang_str

  inline void set_noel_str()
  {
    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      Noel_str = strs.str();
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please supply the number of elements in 1D using --noel.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }

//  inline void generic_setup()
//  {
//    set_ang_str();
//    set_noel_str();
//  }

//  inline std::string prob_str()
//  {
//    set_prob_str();
//    return Prob_str;
//  }

//  inline std::string ang_deg_str()
//  {
//    set_ang_str();
//    return Ang_deg_str;
//  }

//  inline std::string noel_str()
//  {
//    set_noel_str();
//    return Noel_str;
//  }

//  inline std::string create_label()
//  {
//    std::string label = prob_str() + ang_deg_str() + noel_str();
//    return label; 
//  } // inlined function create_label

} // Namespace SquareLagrange


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{

 /// Enumeration for the problem ids
 enum {Driven_cavity, Through_flow};


 /// Reynolds number
 double Re=50.0;

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

 /// Storage for linear solver times during Newton steps 
 Vector<double> Linear_solver_time;

 /// Traction at the outflow boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double> &n,
                          Vector<double>& traction)
 {
  traction.resize(3);
  traction[0]=1.0;
  traction[1]=0.0;
  traction[2]=0.0;
 } 

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
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
    unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      double x=nod_pt->x(0);
      double y=nod_pt->x(1);

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

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 TiltedCavityProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
    // Initialise counters for each newton solve.
    Doc_linear_solver_info_pt->setup_new_time_step();
   }
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

 void actions_after_newton_step()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
//     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 void actions_before_distribute()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   namespace SL = SquareLagrange;

//   if(NSPP::Distribute_problem)
   {
//     if(NSPP::Prob_id == SL::PID_SQ_PO)
     {
//       GenericProblemSetup::delete_flux_elements(Surface_mesh_P_pt);

      rebuild_global_mesh();
     }
//     else
     {
   std::ostringstream err_msg;
   err_msg << "Please set up the distributed bit for problem id: "
//           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 }

 void actions_after_distribute()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   namespace SL = SquareLagrange;

//   if(NSPP::Distribute_problem)
   {
//   if(NSPP::Prob_id == SL::PID_SQ_PO)
   {
     create_parall_outflow_lagrange_elements(1,
         Bulk_mesh_pt,Surface_mesh_P_pt);
     rebuild_global_mesh();
   }
//   else
   {
   std::ostringstream err_msg;
   err_msg << "Please set up the distributed bit for problem id: "
//           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
   }
   }
 }

 /// Doc the solution
 void doc_solution();

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt);

 void set_inflow_BC(const unsigned &b,
                    Mesh* const &bulk_mesh_pt);
 void set_nonslip_BC(const unsigned &b,
                     Mesh* const &bulk_mesh_pt);

private:

 void set_mesh_bc_for_SqPo();
 void set_mesh_bc_for_SqTf();
 void set_mesh_bc_for_SqVa();

 /// Pointer to the "bulk" mesh
 //SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;
 Mesh* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 unsigned Right_bound;
 unsigned Left_bound;
 unsigned Top_bound;
 unsigned Bottom_bound;

};




//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================
class FpTestProblem : public Problem
{

public:


 /// Constructor
 FpTestProblem(const unsigned& n_element,
               const bool& use_tets,
               const bool& use_adaptivity, 
               const bool& use_lsc,
               const bool& use_hypre_for_pressure,
               const bool& use_block_diagonal_for_momentum,
               const bool& use_hypre_for_momentum_diagonals,
               const int& problem_id);

 
 /// Destructor: Cleanup
 ~FpTestProblem()
  {
   delete Solver_pt;
   delete Prec_pt;
   delete P_matrix_preconditioner_pt;
   delete F_matrix_preconditioner_pt;
  }

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 
 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<3>::
    unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());
   
   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 

 /// Actions after Newton step record number of iterations
 void actions_after_newton_step() 
  {                               
   Global_Variables::Iterations.push_back(
    dynamic_cast<IterativeLinearSolver*>
    (this->linear_solver_pt())->iterations());
   
   Global_Variables::Linear_solver_time.push_back(
    linear_solver_pt()->linear_solver_solution_time());
  }  

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {
  // Initialise counter for iterations
  Global_Variables::Iterations.clear();
  Global_Variables::Linear_solver_time.clear();

  // Driven cavity bcs
  if (Problem_id==Global_Variables::Driven_cavity)
   {
    // Setup tangential flow along driven boundary
    unsigned ibound=Driven_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      // Tangential flow
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,1.0);
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      
      // No penetration
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
     }
    
    // Overwrite with no flow along the other boundaries
    unsigned num_bound = Bulk_mesh_pt->nboundary();
    for(unsigned ibound=0;ibound<num_bound;ibound++)
     {
      if (ibound!=Driven_boundary)
       {
        unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
        for (unsigned inod=0;inod<num_nod;inod++)
         {
          for (unsigned i=0;i<3;i++)
           {
            Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(i,0.0);
           }
         }
       }
     }
   }
  // Inflow at left boundary
  else
   {
    // Inflow in upper half of inflow boundary
    unsigned ibound=Inflow_boundary; 
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      double y=nod_pt->x(1);
      double z=nod_pt->x(2);
      if ((y>0.5)&&(z>0.5))
       {
        double u=(y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
        nod_pt->set_value(0,u);
       }
      else
       {
        nod_pt->set_value(0,0.0);
       }
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
     }
   }
  
 } // end_of_actions_before_newton_solve

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Create traction elements on outflow boundary
 template<class ELEMENT>
 void create_traction_elements();

 /// Create refineable traction elements on outflow boundary
 template<class ELEMENT>
 void create_refineable_traction_elements();

 /// Validate fp
 template<class FP_ELEMENT>
 void validate_fp()
  {
   DocInfo my_doc_info;
   dynamic_cast<NavierStokesSchurComplementPreconditioner*>(
    Prec_pt)->template validate<FP_ELEMENT>(my_doc_info,this);
   pause("done validation");
  }

 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 /// ID of driven boundary
 unsigned Driven_boundary;

 /// ID of inflow boundary
 unsigned Inflow_boundary;

 /// ID of outflow boundary
 unsigned Outflow_boundary;

 /// Problem id
 unsigned Problem_id;

 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem()
{
  // Alias the namespace for convenience
//  namespace NSPP = NavierStokesProblemParameters;
//  namespace LPH = LagrangianPreconditionerHelpers;
//  namespace SL = SquareLagrange;

  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;
  
//  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

      /// Setup the mesh
    // # of elements in x-direction
    const unsigned nx=4;//SL::Noel;
    
    // # of elements in y-direction
    const unsigned ny=4;//SL::Noel;
    
    // Domain length in x-direction
    const double lx=1;//SL::Lx;

    // Domain length in y-direction
    const double ly=1;//SL::Ly;
  // First we set the mesh.
//  if((NSPP::Prob_id == SL::PID_SQ_TMP) ||
//     (NSPP::Prob_id == SL::PID_SQ_PO)  ||
//     (NSPP::Prob_id == SL::PID_SQ_TF)  ||
//     (NSPP::Prob_id == SL::PID_SQ_TFPO) )
  {
    // This is the tilted cavity mesh.
    Bulk_mesh_pt =
      new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,30); 
    // RAYRAY remember to convert
//      new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,SL::Ang);
  }
//  else if (NSPP::Prob_id == SL::PID_SQ_VA)
  {
    Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly);
  }
//  else
  {
//   std::ostringstream err_msg;
//   err_msg << "There is no mesh for the problem ID: "
//           << NSPP::Prob_id << ".\n"
//           << std::endl;
//
//   throw OomphLibError(err_msg.str(),
//       OOMPH_CURRENT_FUNCTION,
//       OOMPH_EXCEPTION_LOCATION);
  }

  // Set the boundary conditions
//  if(NSPP::Prob_id == SL::PID_SQ_PO)
  {
    set_mesh_bc_for_SqPo();
  }
//  else if(NSPP::Prob_id == SL::PID_SQ_TF)
  {
    set_mesh_bc_for_SqTf();
  }
//  else if(NSPP::Prob_id == SL::PID_SQ_VA)
  {
    set_mesh_bc_for_SqVa();
  }
//  else
  {
//    std::ostringstream err_msg;
//   err_msg << "There are no boundary conditions configured for prob_id: "
//           << NSPP::Prob_id << ".\n"
//           << std::endl;
//
//   throw OomphLibError(err_msg.str(),
//       OOMPH_CURRENT_FUNCTION,
//       OOMPH_EXCEPTION_LOCATION);
  }




// // Top boundary is slip.
// current_bound = 2;
// num_nod= mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   if(!nod_pt->is_on_boundary(3))
//   {
//     nod_pt->unpin(0);
//     nod_pt->pin(1);
//
//     nod_pt->set_value(1,0.0);
//   }
// }


 
 //set_nonslip_BC(0,Bulk_mesh_pt);
// set_nonslip_BC(2,Bulk_mesh_pt);

// set_inflow_BC(if_b,Bulk_mesh_pt);
 
 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
//   el_pt->re_pt() = &NSPP::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assign equation numbers
 oomph_info << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;
 

 // Only do this bit if we do NOT have a direct solver.
// if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
 {
   // Create the vector of mesh pointers!
   Vector<Mesh*> mesh_pt;
//   if(NSPP::Prob_id == SL::PID_SQ_PO)
   {
     mesh_pt.resize(2,0);
     mesh_pt[0] = Bulk_mesh_pt;
     mesh_pt[1] = Surface_mesh_P_pt;
   }
//   else if(NSPP::Prob_id == SL::PID_SQ_TF)
   {
//     mesh_pt.resize(2,0);
//     mesh_pt[0] = Bulk_mesh_pt;
//     mesh_pt[1] = Surface_mesh_T_pt;
   }
//   else if(NSPP::Prob_id == SL::PID_SQ_VA)
   {
//     mesh_pt.resize(1,0);
//     mesh_pt[0] = Bulk_mesh_pt;
   }

   // Quick check that the correct preconditioner is chosen.
//   if((NSPP::Prob_id == SL::PID_SQ_VA) && 
//       !CommandLineArgs::command_line_flag_has_been_set("--lsc_only"))
//   {
//     std::ostringstream err_msg;
//     err_msg << "You have requested Vanilla Navier-Stokes problem,\n"
//       << "NSPP::Prob_id is " << NSPP::Prob_id << "\n"
//       << "But you have not set the flag --lsc_only.\n" 
//       << "Please choose you preconditioner parameters again.\n"
//       << std::endl;
//
//     throw OomphLibError(err_msg.str(),
//         OOMPH_CURRENT_FUNCTION,
//         OOMPH_EXCEPTION_LOCATION);
//   }

//   LPH::Mesh_pt = mesh_pt;
//   LPH::Problem_pt = this;
//   Prec_pt = LPH::get_preconditioner();
 }
 const double solver_tol = 1.0e-6;
 const double newton_tol = 1.0e-6;
// GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
//                                   solver_tol,newton_tol,
//                                   NSPP::Solver_type,this,Prec_pt);
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{

//  namespace NSPP = NavierStokesProblemParameters;
  
//  std::ofstream some_file;
//  std::stringstream filename;
//  filename << NSPP::Soln_dir_str<<"/"<<NSPP::Label_str<<".dat";

  // Number of plot points
//  unsigned npts=5;

  // Output solution
//  some_file.open(filename.str().c_str());
//  Bulk_mesh_pt->output(some_file,npts);
//  some_file.close();
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqPo()
{
  // Alias the namespace for convenience
//  namespace SL = SquareLagrange;

  // Assign the boundaries:
  //             2 non slip
  //         ----------
  //         |        |
  // 3 Inflow|        |1 P.O.
  //         |        |
  //         ----------
  //             0 non slip
//  unsigned if_b = 3; // inflow
  const unsigned po_b = 1; // parallel outflow

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_P_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
  create_parall_outflow_lagrange_elements(po_b,
                                          Bulk_mesh_pt,Surface_mesh_P_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_P_pt);
  
  // combine all sub-meshes into a single mesh.
  build_global_mesh();
  const unsigned num_bound = mesh_pt()->nboundary();

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  { 
    if(ibound != po_b)
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

 // Which boundary are we dealing with?
 unsigned current_bound;
 
 // The number of nodes on a boundary.
 unsigned num_nod;

 // Inflow is at boundary 3
 current_bound = 3;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the x and y cartesian coordinates
   double x0=nod_pt->x(0);
   double x1=nod_pt->x(1);

   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
//   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);  RRR

   // Now calculate the parabolic inflow at this point.
//   double u0_old = (x1_old - SL::Y_min)*(SL::Y_max - x1_old); RRR
   
   // Now apply the rotation to u0_old, using rotation matrices.
   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
//   double u0=u0_old*cos(SL::Ang); RRR
//   double u1=u0_old*sin(SL::Ang); RRR

//   nod_pt->set_value(0,u0); RRR
//   nod_pt->set_value(1,u1); RRR
 }
} // set_mesh_bc_for_SqPo

////============RAYRAY===========
///// RAYRAY
////=======================================================================
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqTf()
//{
//  // Alias the namespace for convenience
//  namespace SL = SquareLagrange;
//
//  // Assign the boundaries:
//  //             2 slip bc
//  //         ----------
//  //         |        |
//  // 3 Inflow|        |1 Imposed outflow
//  //         |        |
//  //         ----------
//  //             0 non slip
////  const unsigned if_b = 3; // inflow
////  const unsigned po_b = 1; // parallel outflow
//  const unsigned slip_b = 2;
//
//  // Create a "surface mesh" that will contain only
//  // ImposeParallelOutflowElements in boundary 1
//  // The constructor just creates the mesh without
//  // giving it any elements, nodes, etc.
//  Surface_mesh_T_pt = new Mesh;
//
//  // Create ImposeParallelOutflowElement from all elements that are
//  // adjacent to the Neumann boundary.
//  create_impenetrable_lagrange_elements(slip_b,Bulk_mesh_pt,
//                                        Surface_mesh_T_pt);
////  create_parall_outflow_lagrange_elements(slip_b,
////                                          Bulk_mesh_pt,Surface_mesh_T_pt);
//
//  // Add the two meshes to the problem.
//  add_sub_mesh(Bulk_mesh_pt);
//  add_sub_mesh(Surface_mesh_T_pt);
//  
//  // combine all sub-meshes into a single mesh.
//  build_global_mesh();
////  const unsigned num_bound = mesh_pt()->nboundary();
//
//
// // Which boundary are we dealing with?
// unsigned current_bound;
// 
// // The number of nodes on a boundary.
// unsigned num_nod;
//
// // Inflow is at boundary 3
// current_bound = 3;
// num_nod= mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   // Pin both velocity components
//   nod_pt->pin(0);
//   nod_pt->pin(1);
//
//   // Get the x and y cartesian coordinates
//   double x0=nod_pt->x(0);
//   double x1=nod_pt->x(1);
//
//   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
//   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);
//
//   // Now calculate the parabolic inflow at this point.
//   double u0_old = (x1_old - SL::Y_min)*(2*SL::Y_max - x1_old);
//   
//   // Now apply the rotation to u0_old, using rotation matrices.
//   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
//   // velocity in the x direction only. There is no velocity
//   // in the y direction.
//   double u0=u0_old*cos(SL::Ang);
//   double u1=u0_old*sin(SL::Ang);
//
//   nod_pt->set_value(0,u0);
//   nod_pt->set_value(1,u1);
// }
//
// // Now do the outflow.
// // This is on boundary 1.
// current_bound = 1;
// num_nod= mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   // Pin both velocity components
//   nod_pt->pin(0);
//   nod_pt->pin(1);
//
//   // Get the x and y cartesian coordinates
//   double x0=nod_pt->x(0);
//   double x1=nod_pt->x(1);
//
//   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
//   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);
//
//   // Now calculate the parabolic inflow at this point.
//   double u0_old = (x1_old - SL::Y_min)*(2*SL::Y_max - x1_old);
//   
//   // Now apply the rotation to u0_old, using rotation matrices.
//   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
//   // velocity in the x direction only. There is no velocity
//   // in the y direction.
//   double u0=u0_old*cos(SL::Ang);
//   double u1=u0_old*sin(SL::Ang);
//
////   nod_pt->unpin(0);
//   nod_pt->set_value(0,u0);
//   nod_pt->set_value(1,u1);
// }
//
// // Now we do the bottom boundary, this is boundary 0
// current_bound = 0;
// num_nod= mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   // Pin both velocity components
//   nod_pt->pin(0);
//   nod_pt->pin(1);
//
//   nod_pt->set_value(0,0.0);
//   nod_pt->set_value(1,0.0);
// }
//
//} // set_mesh_bc_for_SqPo
//
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqVa()
//{
//  // Alias the namespace for convenience
//  namespace SL = SquareLagrange;
//
//  // Assign the boundaries:
//  //             2 non slip
//  //         ----------
//  //         |        |
//  // 3 Inflow|        |1 P.O.
//  //         |        |
//  //         ----------
//  //             0 non slip
////  unsigned if_b = 3; // inflow
//  unsigned po_b = 1; // parallel outflow
//
//  // Create a "surface mesh" that will contain only
//  // ImposeParallelOutflowElements in boundary 1
//  // The constructor just creates the mesh without
//  // giving it any elements, nodes, etc.
////  Surface_mesh_P_pt = new Mesh;
//
//  // Create ImposeParallelOutflowElement from all elements that are
//  // adjacent to the Neumann boundary.
////  create_parall_outflow_lagrange_elements(po_b,
////                                          Bulk_mesh_pt,Surface_mesh_P_pt);
//
//  // Add the two meshes to the problem.
//  add_sub_mesh(Bulk_mesh_pt);
////  add_sub_mesh(Surface_mesh_P_pt);
//  
//  // combine all sub-meshes into a single mesh.
//  build_global_mesh();
//  const unsigned num_bound = mesh_pt()->nboundary();
//
//  // Leave the x velocity on po_b to unpinned.
//  {
//    unsigned num_nod = mesh_pt()->nboundary_node(po_b);
//    for (unsigned inod = 0; inod < num_nod; inod++) 
//    {
//      Node* nod_pt = mesh_pt()->boundary_node_pt(po_b,inod);
//      // Unpin x
//      nod_pt->unpin(0);
//
//      // Pin y
//      nod_pt->pin(1);
//      // Set y value to zero.
//      nod_pt->set_value(1,0);
//    }
//  }
//
//
//  // Set the boundary conditions for this problem: All nodes are
//  // free by default -- just pin the ones that have Dirichlet conditions
//  // here.
//  for(unsigned ibound=0;ibound<num_bound;ibound++)
//  { 
//    if(ibound != po_b)
//    {
//      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
//      for (unsigned inod=0;inod<num_nod;inod++)
//      {
//        // Get node
//        Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
// 
//        nod_pt->pin(0);
//        nod_pt->pin(1);
//        
//        nod_pt->set_value(0,0);
//        nod_pt->set_value(1,0);
// 
//      }
//    }
//  }
//
// // Which boundary are we dealing with?
// unsigned current_bound;
// 
// // The number of nodes on a boundary.
// unsigned num_nod;
//
// // Inflow is at boundary 3
// current_bound = 3;
// num_nod=mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   // Pin both velocity components
//   nod_pt->pin(0);
//   nod_pt->pin(1);
//
//   // Get the y cartesian coordinates
//   double x1=nod_pt->x(1);
//
//   // Now calculate the parabolic inflow at this point.
//   double u0 = (x1 - SL::Y_min)*(SL::Y_max - x1);
//   
//   double u1=0.0;
//
//   nod_pt->set_value(0,u0);
//   nod_pt->set_value(1,u1);
// }
//} // set_mesh_bc_for_SqPo
//
////============RAYRAY===========
///// RAYRAY
////=======================================================================
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::
//set_nonslip_BC(const unsigned &b,
//               Mesh* const &bulk_mesh_pt)
//{
//  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
//  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
//   
//  for(unsigned inod=0;inod<num_nod;inod++)
//   {
//    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);
//    
//    // pin all velocity components and set the value to zero.
//    for (unsigned velo_i = 0; velo_i < dim; velo_i++) 
//    {
//      nod_pt->pin(velo_i);
//      nod_pt->set_value(velo_i,0);
//    }
//   }
//}
//
////============RAYRAY===========
///// RAYRAY
////=======================================================================
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::
//set_inflow_BC(const unsigned &b,
//              Mesh* const &bulk_mesh_pt)
//{
//
// // Alias the namespace for convenience
// namespace SL = SquareLagrange;
//
// // Check that the dimension is correct.
//#ifdef PARANOID
//  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
//  if(dim != 2)
//   {
//     std::ostringstream err_msg;
//     err_msg << "Inflow implemented for dim = 2 only." << std::endl;
//
//     throw OomphLibError(err_msg.str(),
//                         OOMPH_CURRENT_FUNCTION,
//                         OOMPH_EXCEPTION_LOCATION);
//   }
//#endif
//
//  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
//   
//  for(unsigned inod=0;inod<num_nod;inod++)
//   {
//    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);
//
//    // Pin both velocity components
//    nod_pt->pin(0);
//    nod_pt->pin(1);
//
//    // Get the x and y cartesian coordinates.
//    double x0 = nod_pt->x(0);
//    double x1 = nod_pt->x(1);
//
//    // Tilt x1 back the coordinate so we get the original coordinate.
//    double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);
//    
//    // Now calculate the parabolic inflow at this point
//    //double u0_old = (x1_old - SL::Y_min)*(SL::Y_max - x1_old);
//    double u0_old = (x1_old - SL::Y_min)*(2.0 - x1_old);
//
//    // Now apply the rotation to u0_old, using rotation matrices.
//    // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
//    // velocity in the x direction only. There is no velocity
//    // in the y direction.
//    double u0=u0_old*cos(SL::Ang);
//    double u1=u0_old*sin(SL::Ang);
//    
//    nod_pt->set_value(0,u0);
//    nod_pt->set_value(1,u1); 
//   }
//}


//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
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

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_impenetrable_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new
    ImposeImpenetrabilityElement<ELEMENT>(bulk_elem_pt,
                                          face_index);
//   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
//    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
//                                          face_index);

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 3 or 1?
     if ((nod_pt->is_on_boundary(3))||(nod_pt->is_on_boundary(1)))
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


//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
//========================================================================
FpTestProblem::FpTestProblem(const unsigned& n_el,
                             const bool& use_tets,
                             const bool& use_adaptivity, 
                             const bool& use_lsc,
                             const bool& use_hypre_for_pressure,
                             const bool& use_block_diagonal_for_momentum,
                             const bool& use_hypre_for_momentum_diagonals,
                             const int& problem_id)
{ 
 
 // Store flag
 Problem_id=problem_id;

 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x=n_el;
 
 // # of elements in y-direction
 unsigned n_y=n_el;

 // # of elements in z-direction
 unsigned n_z=n_el;
 
 // Domain length in x-direction
 double l_x=1.0;
 
 // Domain length in y-direction
 double l_y=1.0;
 
 // Domain length in y-direction
 double l_z=1.0;
 

 // Build and assign mesh
 if (use_tets)
  {
   Bulk_mesh_pt = new SimpleCubicTetMesh<TTaylorHoodElement<3> >
    (n_x,n_y,n_z,l_x,l_y,l_z);
     
   dynamic_cast<SimpleCubicTetMesh<TTaylorHoodElement<3> >*>(Bulk_mesh_pt)
    ->setup_boundary_element_info();

   Driven_boundary=4;
   Inflow_boundary=0;
   Outflow_boundary=1;
  }
 else
  {
   if (use_adaptivity)
    {
     Bulk_mesh_pt = 
      new RefineableSimpleCubicMesh<RefineableQTaylorHoodElement<3> >
      (n_x,n_y,n_z,l_x,l_y,l_z);
     
     Driven_boundary=0;
     Inflow_boundary=4;
     Outflow_boundary=2;     
    }
   else
    {
     Bulk_mesh_pt = 
      new SimpleCubicMesh<QTaylorHoodElement<3> >(n_x,n_y,n_z,l_x,l_y,l_z);
     
     Driven_boundary=0;
     Inflow_boundary=4;
     Outflow_boundary=2;
    }
  } 

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements.
 Surface_mesh_pt = new Mesh;

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // Build preconditoner
 NavierStokesSchurComplementPreconditioner* prec_pt = 
  new NavierStokesSchurComplementPreconditioner(this);
 Prec_pt=prec_pt;
   
 
 // By default, the Schur Complement Preconditioner uses SuperLU as
 // an exact preconditioner (i.e. a solver) for the
 // momentum and Schur complement blocks. 
 // Can overwrite this by passing pointers to 
 // other preconditioners that perform the (approximate)
 // solves of these blocks.


 // Create internal preconditioners used on Schur block
 P_matrix_preconditioner_pt=0;
 if (use_hypre_for_pressure)
  {
#ifdef OOMPH_HAS_HYPRE

   // Create preconditioner
   P_matrix_preconditioner_pt = new HyprePreconditioner;
   
   // Set parameters for use as preconditioner on Poisson-type problem
   Hypre_default_settings::set_defaults_for_2D_poisson_problem(
    static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
   
   // Use Hypre for the Schur complement block
   prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
   
   // Shut up!
   static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
    disable_doc_time();
   
#endif
  }
 
 // Create block-diagonal preconditioner used on momentum block
 F_matrix_preconditioner_pt=0;   
 if (use_block_diagonal_for_momentum)
  {
   
   F_matrix_preconditioner_pt = 
    new BlockDiagonalPreconditioner<CRDoubleMatrix>;
   
   // Use Hypre as block preconditioner
   if (use_hypre_for_pressure)
    {
#ifdef OOMPH_HAS_HYPRE
     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
      (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
#endif       
    }
   
   // Use Hypre for momentum block 
   prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
  }
 
 
   
   // Use LSC?
   if (use_lsc)
    {
     prec_pt->use_lsc();
    }
   else
    {
     prec_pt->use_fp();
    }
   
   // Set Navier Stokes mesh
   prec_pt->set_navier_stokes_mesh(Bulk_mesh_pt);    
   
   
   // Set the boundary conditions for this problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   unsigned num_bound = Bulk_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Loop over values (u, v and w velocities)
       for (unsigned i=0;i<3;i++)
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
        }
      }
    } // end loop over boundaries

 

 // In/outflow bcs
 if (Problem_id==Global_Variables::Through_flow)
  {
   unsigned ibound=Outflow_boundary;
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     // Only free if node is ONLY on a single boundary
     std::set<unsigned>* bnd_pt=0;
     nod_pt->get_boundaries_pt(bnd_pt);
     if (bnd_pt!=0)
      {
       if (bnd_pt->size()<2)
        {
         if (!(nod_pt->is_on_boundary(0)))
          {
           if ((nod_pt->x(1)<0.5)&&nod_pt->x(2)<0.5) nod_pt->unpin(0);
          }
        }
      }
    }
  }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   NavierStokesEquations<3>* el_pt = 
    dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Variables::Re;
  } // end loop over elements
 
 
 // Pin redundant pressure dofs
 if (use_adaptivity)
  {
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
  }

 // Now set the first pressure value in element 0 to 0.0
 if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

#ifdef OOMPH_HAS_TRILINOS

 // Build iterative linear solver
 oomph_info << "Using Trilinos GMRES\n"; 
 TrilinosAztecOOSolver* iterative_linear_solver_pt = 
  new TrilinosAztecOOSolver;

 Solver_pt=iterative_linear_solver_pt;

#else

 // Build solve and preconditioner
 Solver_pt = new GMRES<CRDoubleMatrix>;
 dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

#endif

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;
 
} // end_of_constructor


//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on outflow boundary
//=======================================================================
template<class ELEMENT>
void FpTestProblem::create_refineable_traction_elements()
{

 unsigned b=Outflow_boundary;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding prescribed-flux element
   RefineableNavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
    RefineableNavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   Surface_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set the pointer to the prescribed traction function
   flux_element_pt->traction_fct_pt() = &Global_Variables::prescribed_traction;
   
  } //end of loop over bulk elements adjacent to boundary b

 // Now rebuild the global mesh
 rebuild_global_mesh();

 // Reassign equation numbers
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of create_traction_elements






//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on outflow boundary
//=======================================================================
template<class ELEMENT>
void FpTestProblem::create_traction_elements()
{

 unsigned b=Outflow_boundary;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding prescribed-flux element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
      NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   Surface_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set the pointer to the prescribed traction function
   flux_element_pt->traction_fct_pt() = &Global_Variables::prescribed_traction;
   
  } //end of loop over bulk elements adjacent to boundary b

 // Now rebuild the global mesh
 rebuild_global_mesh();

 // Reassign equation numbers
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of create_traction_elements



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
void FpTestProblem::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_main======================================================
/// Driver for Fp preconditioner
//=====================================================================
int main(int argc, char **argv)
{


#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


 //Label for output
 DocInfo doc_info;
 
 //Set output directory
 doc_info.set_directory("RESLT");
 
 //Doc number of gmres iterations
 ofstream out_file;
 
 // Set flags
 bool use_hypre_for_pressure=true;
 bool use_block_diagonal_for_momentum=false;
 bool use_hypre_for_momentum_diagonals=false;
   
 //Loop over problems: Driven cavity and step
 for (unsigned problem_id=0;problem_id<2;problem_id++) 
  {
   
   if (problem_id==0)
    {
     out_file.open("three_d_iter_driven_cavity.dat");
    }
   else 
    {      
     out_file.open("three_d_iter_through_flow.dat");
    }
   
   out_file
    << "VARIABLES=\"nel_1d\","
    << "\"ndof\"," 
    << "\"Re\"," 
    << "\" Newton iteration\","
    << "\"GMRES iterations\","
    << "\"Linear solver time\","
    << "\"doc number\""
    << std::endl;
   
   std::string header1;
   std::string header2;
   
     
   // Loop over preconditioners: iprec=0: LSC
   //                            iprec=1: Fp
   bool use_lsc=true;
   //bool use_robin=true;
   for (unsigned iprec=0;iprec<2;iprec++)
    {
     
     // Loop over three cases (tets, non-refineable/refineable bricks)
     for (unsigned icase=0;icase<=2;icase++) 
      {
       bool use_tets=false;
       bool use_adaptivity=false;
       
       switch(icase)
        {
         
        case 0:
         
         header1=" Tets";
         use_tets=true;
         use_adaptivity=false;
         break;
         
        case 1:
         
         header1=" Bricks";
         use_tets=false;
         use_adaptivity=false;
         break;
         
        case 2:
         
         header1=" Refineable Bricks";
         use_tets=false;
         use_adaptivity=true;
         break;
         
        default:
         break;
        }
       
       oomph_info << "Doing it with " << header1 << std::endl;
       
       // Set preconditioner
       if (iprec==0)
        {
         use_lsc=true;
         header2=", LSC";
        }
       else if (iprec==1)
        {
         use_lsc=false;
         //use_robin=true;
         header2=", Fp";
        }   
       
        oomph_info << "Doing it with " << header2 << " preconditioner\n";

        // Write tecplot header
        string header="ZONE T=\""+header1+header2+"\"\n";
        out_file << header;
                 
        // Number of elements in x/y directions (reduced for validaton)
        unsigned max_nel_1d=4; 
        if (argc>1) max_nel_1d=2;
        for (unsigned nel_1d = 2; nel_1d <= max_nel_1d; nel_1d*=2) 
         {
           
          // Build the problem 
          FpTestProblem problem(
           nel_1d,use_tets,use_adaptivity,use_lsc,
           use_hypre_for_pressure,use_block_diagonal_for_momentum,
           use_hypre_for_momentum_diagonals,problem_id);
           
          
          
          // Refine a few times
          if (use_adaptivity)
           {             
            // Note: This manually chosen refinement pattern introduces
            // doubly hanging nodes on boundary and in the interior
            Vector<unsigned> elements_to_be_refined;
            unsigned nel=problem.bulk_mesh_pt()->nelement();
            unsigned e=0;
            while (e<nel)
             {
              elements_to_be_refined.push_back(e);
              e+=7;
             }
            problem.refine_selected_elements(0,elements_to_be_refined); 
            elements_to_be_refined.clear();
            elements_to_be_refined.push_back(1);
            problem.refine_selected_elements(0,elements_to_be_refined); 
            elements_to_be_refined.clear();
            elements_to_be_refined.push_back(4);
            problem.refine_selected_elements(0,elements_to_be_refined); 
           }
          // Attach traction elements now with the problem in its
          // most refined state
          if (use_tets)
           {
            problem.create_traction_elements<TTaylorHoodElement<3> >();
           }
          else 
           {
            if (use_adaptivity)
             {
              problem.create_refineable_traction_elements
               <RefineableQTaylorHoodElement<3> >();
             }
            else
             {
              problem.create_traction_elements<QTaylorHoodElement<3> >();
             }
           }
          
          
          
          // Loop over Reynolds numbers (limited during validation)
          double start_re = 0.0; 
          if (argc>1) start_re=50;
          double end_re = 50.0; 
          for (double re = start_re; re <= end_re; re+=50.0)
           {
            
            // Set Reynolds
            Global_Variables::Re=re;
            
            // Solve the problem 
            problem.newton_solve();
            
            // Doc solution
            problem.doc_solution(doc_info);
            doc_info.number()++;
            
            // Doc iteration counts for each Newton iteration
            unsigned ndof = problem.ndof();
            unsigned iter = Global_Variables::Iterations.size();

            // Doc for all Newton iterations or only the last one?
            unsigned j_lo=0;
            //j_lo=iter-1;
            for (unsigned j = j_lo; j < iter; j++)
             {
              out_file
               << nel_1d << " "
               << ndof << " "
               << re << " " 
               << j << " "
               << Global_Variables::Iterations[j] << " "
               << Global_Variables::Linear_solver_time[j] << " "
               << doc_info.number()-1 << " "
               << std::endl;
             }
            
           }
         }     
      }
     }
   
   out_file.close();
   
  }
 

#ifdef OOMPH_HAS_MPI
   MPI_Helpers::finalize();
#endif
   
   
} // end_of_main



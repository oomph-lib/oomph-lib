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
#include "navier_stokes.h"
#include "beam.h"
#include "multi_physics.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/collapsible_channel_mesh.h"

using namespace std;

using namespace oomph;

// Include the general-purpose fsi collapsible channel problem
#include "fsi_chan_problem.h"


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//====Namespace_for_flags================================
/// Extend namespace for flags
//======================================================
namespace Flags
{

 /// Solver flag [0: direct; 1: exact; 2: simple; 3: Schur/Superlu; 
 /// 4: Schur/Hypre]
 unsigned Solver_flag=0;

 /// Solver sub flag [0: diag; 1: retain fluid on solid; 
 /// 2: retain solid on fluid]
 unsigned Solver_sub_flag=0;

}



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
 



//====start_of_problem_class==========================================
/// Problem class
//====================================================================
template <class ELEMENT>
class PreconditionedFSICollapsibleChannelProblem : 
 public virtual FSICollapsibleChannelProblem<ELEMENT>
{

public :

/// Constructor: The arguments are the number of elements and
/// the lengths of the domain.
PreconditionedFSICollapsibleChannelProblem(const unsigned& nup, 
                                           const unsigned& ncollapsible,
                                           const unsigned& ndown,
                                           const unsigned& ny,
                                           const double& lup,
                                           const double& lcollapsible, 
                                           const double& ldown,
                                           const double& ly,
                                           const bool& displ_control,
                                           const bool& steady_flag,
                                           const unsigned& solver_flag,
                                           const unsigned& solver_sub_flag
                                           ) :
 FSICollapsibleChannelProblem<ELEMENT>(nup, 
                                       ncollapsible,
                                       ndown,
                                       ny,
                                       lup,
                                       lcollapsible, 
                                       ldown,
                                       ly,
                                       displ_control,
                                       steady_flag)
  {

   // Build iterative linear solver
   GMRES<CRDoubleMatrix>* iterative_linear_solver_pt = 
     new GMRES<CRDoubleMatrix>;
   
   // Set maximum number of iterations
   iterative_linear_solver_pt->max_iter() = 100;
   
   // Set tolerance
   iterative_linear_solver_pt->tolerance() = 1.0e-6;   

   // Choose solver: Direct vs. iterative
   if (solver_flag==0)
    {
     std::cout << "Using default direct solver." << std::endl;
    }
   else
    {
     std::cout << "Using iterative linear solver." << std::endl;
     this->linear_solver_pt()=iterative_linear_solver_pt;
    } //solver chosen
   
    
   // Choose solver/preconditioner
   switch (solver_flag)
    {

     // Direct solver -- no action required
     //====================================
    case 0:

     Flags::Run_identifier_string="Direct solver";

     break;


     // Exact preconditioner (re-arranged Jacobian)
     //============================================
    case 1:
     
     std::cout << "Using ExactFSIPreconditioner (re-arranged Jacobian)" 
               << std::endl;

     Flags::Run_identifier_string=
      "Exact preconditioner (re-arranged Jacobian)";
     
     {
      ExactBlockPreconditioner<CRDoubleMatrix>* prec_pt=
       new  ExactBlockPreconditioner<CRDoubleMatrix>;
     
      // With GeneralPurposeBlockPreconditioner, we set meshes with
      // add_mesh(...). The order of the meshes matter if we want to
      // use different subsidiary preconditioners for different blocks.

      // Set Navier Stokes mesh.
      prec_pt->add_mesh(this->bulk_mesh_pt()); 
      
      // Build  the mesh that contains all solid elements:
      
      // Create a vector of pointers to submeshes. Start with the solid
      // mesh itself.
      Vector<Mesh*> s_mesh_pt(1);
      s_mesh_pt[0]=this->wall_mesh_pt();
      
      // Add the displacement control mesh if required
      if (this->Displ_control) 
       s_mesh_pt.push_back(this->Displ_control_mesh_pt);
      
      // Build "combined" mesh from vector of solid submeshes
      Mesh* solid_mesh_pt = new Mesh(s_mesh_pt);

      // Set solid mesh with "true" to tolerate multiple element types
      // in the mesh.
      prec_pt->add_mesh(solid_mesh_pt,true);


      // Set preconditioner
      iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
     }
     
     break;

      // Simple preconditioner (block triangular Jacobian, solved exactly)
      //==================================================================
    case 2:

     std::cout << "Using SimpleFSIPreconditioner (block-triangular Jacobian)" 
               << std::endl;
     
     {
      SimpleFSIPreconditioner<CRDoubleMatrix>* prec_pt=
       new  SimpleFSIPreconditioner<CRDoubleMatrix>;
      
      // Set Navier Stokes mesh:
      prec_pt->set_navier_stokes_mesh(this->bulk_mesh_pt());
      
      // Build a compound mesh that contains all solid elements:
      
      // Create a vector of pointers to submeshes. Start with the solid
      // mesh itself.
      Vector<Mesh*> s_mesh_pt(1);
      s_mesh_pt[0]=this->wall_mesh_pt();
      
      // Add the displacement control mesh if required
      if (this->Displ_control) 
       s_mesh_pt.push_back(this->Displ_control_mesh_pt);
      
      // Build "combined" mesh from vector of solid submeshes
      Mesh* solid_mesh_pt = new Mesh(s_mesh_pt);
      
      // Set solid mesh with "true" to tolerate multiple element types in
      // the mesh.
      prec_pt->set_wall_mesh(solid_mesh_pt,true);
            
      
      switch (solver_sub_flag)
       {
        
        // Block diagonal
        //---------------
       case 0:
        
        std::cout << "(diagonal version)" << std::endl;
        Flags::Run_identifier_string=
         "Simple FSI preconditioner (diagonal)";
        prec_pt->use_block_diagonal_version();
        break;
        
        
        // Retain fluid onto solid terms
        //------------------------------
       case 1:
        
        std::cout << "(retaining fluid onto solid terms)" << std::endl;
        Flags::Run_identifier_string=
         "Simple FSI preconditioner (retaining fluid onto solid terms)";
        prec_pt->use_block_triangular_version_with_fluid_on_solid();
        break;
        
        // Retain solid onto fluid terms
        //------------------------------
       case 2:
        
        std::cout << "(retaining solid onto fluid terms)" << std::endl;

        Flags::Run_identifier_string=
         "Simple FSI preconditioner (retaining solid onto fluid terms)";
        prec_pt->use_block_triangular_version_with_solid_on_fluid();
        break;
        
        // Error
        //------
       default:
      
        std::ostringstream error_stream; 
        error_stream << "Error: Wrong solver_sub_flag: " 
                     << solver_sub_flag << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;
        
       }

      iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
     }
     break;
       

     // FSI preconditioner (block triangular Jacobian, solved with 
     //===========================================================
     // block decomposition, using least-squares-commutator (BFBt)
     //===========================================================
     // Navier-Stokes preconditioner for fluid block
     //=============================================
     
    case 3:

     std::cout << "Using FSIPreconditioner (block-triangular Jacobian + LSC)" 
               << std::endl;
     
     {

      // Create an instance of the FSI preconditioner -- pass the pointer
      // to the problem
      FSIPreconditioner* prec_pt=new FSIPreconditioner(this);

      // Set Navier Stokes mesh:
      prec_pt->set_navier_stokes_mesh(this->bulk_mesh_pt());
      
      // Build a compound mesh that contains all solid elements:
      
      // Create a vector of pointers to submeshes. Start with the solid
      // mesh itself.
      Vector<Mesh*> s_mesh_pt(1);
      s_mesh_pt[0]=this->wall_mesh_pt();
      
      // Add the displacement control mesh if required
      if (this->Displ_control) 
       {
        s_mesh_pt.push_back(this->Displ_control_mesh_pt);
       }

      // Build compound mesh from vector of solid submeshes
      Mesh* combined_solid_mesh_pt = new Mesh(s_mesh_pt);

      // Set solid mesh and tolerate multiple element types this is mesh.
      prec_pt->set_wall_mesh(combined_solid_mesh_pt,true);


      switch (solver_sub_flag)
       {
        
        // Block diagonal
        //---------------
       case 0:
        
        std::cout << "(diagonal version)" << std::endl;
        Flags::Run_identifier_string=
         "LSC FSI preconditioner (diagonal)";

        // Use block-diagonal preconditioner
        prec_pt->use_block_diagonal_version();
        break;
        
        
        // Retain fluid onto solid terms
        //------------------------------
       case 1:
        
        std::cout << "(retaining fluid onto solid terms)" << std::endl;
        Flags::Run_identifier_string=
         "LSC FSI preconditioner (retaining fluid onto solid terms)";

        // Choose preconditioner that retains fluid on solid terms
        prec_pt->use_block_triangular_version_with_fluid_on_solid();

        break;
        
        // Retain solid onto fluid terms
        //------------------------------
       case 2:
        
        std::cout << "(retaining solid onto fluid terms)" << std::endl;

        Flags::Run_identifier_string=
         "LSC FSI preconditioner (retaining solid onto fluid terms)";

        // Choose preconditioner that retains solid on fluid terms
        prec_pt->use_block_triangular_version_with_solid_on_fluid();
        break;
        
        // Error
        //------
       default:

        std::ostringstream error_stream; 
        error_stream << "Error: Wrong solver_sub_flag: " 
                     << solver_sub_flag << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;

       }

      // Pass preconditioner to iterative linear solver
      iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
     }

     break;
       
     // FSI preconditioner (block triangular Jacobian, solved with 
     //===========================================================
     // block decomposition, using least-squares-commutator (BFBt)
     //===========================================================
     // Navier-Stokes preconditioner for fluid block, Hypre AMG
     //========================================================
     // for P
     //======
     
    case 4:

     std::cout << "Using FSIPreconditioner (block-triangular Jacobian + LSC)" 
               << std::endl;
     
     {

      // Create an instance of the FSI preconditioner -- pass the pointer
      // to the problem
      FSIPreconditioner* prec_pt=new FSIPreconditioner(this);

      // Set Navier Stokes mesh:
      prec_pt->set_navier_stokes_mesh(this->bulk_mesh_pt());
      
      // Build a compound mesh that contains all solid elements:
      
      // Create a vector of pointers to submeshes. Start with the solid
      // mesh itself.
      Vector<Mesh*> s_mesh_pt(1);
      s_mesh_pt[0]=this->wall_mesh_pt();
      
      // Add the displacement control mesh if required
      if (this->Displ_control) 
       {
        s_mesh_pt.push_back(this->Displ_control_mesh_pt);
       }

      // Build compound mesh from vector of solid submeshes
      Mesh* combined_solid_mesh_pt = new Mesh(s_mesh_pt);

      // Set solid mesh with true to tolerate multiple element types in
      // the mesh.
      prec_pt->set_wall_mesh(combined_solid_mesh_pt,true);

#ifdef OOMPH_HAS_HYPRE
//If we are using MPI, then only use HYPRE if it has been initialised
#ifdef OOMPH_HAS_MPI
      if(MPI_Helpers::mpi_has_been_initialised())
#endif
       {
        
        // By default, the LSC Preconditioner uses SuperLU as
        // an exact preconditioner (i.e. a solver) for the
        // momentum and Schur complement blocks. 
        // Can overwrite this by passing pointers to 
        // other preconditioners that perform the (approximate)
        // solves of these blocks.
        
        // Create internal preconditioners used on Schur block
        HyprePreconditioner* P_matrix_preconditioner_pt = 
         new HyprePreconditioner;
        
        // Set defaults parameters for use as preconditioner on Poisson-type 
        // problem
        Hypre_default_settings::set_defaults_for_2D_poisson_problem(
         P_matrix_preconditioner_pt);
        
        // Use Hypre for the Schur complement block
        prec_pt->navier_stokes_preconditioner_pt()->
         set_p_preconditioner(P_matrix_preconditioner_pt);
        
        // Shut up
        P_matrix_preconditioner_pt->disable_doc_time();
       }
#endif // endif for we have hypre...

      switch (solver_sub_flag)
       {
        
        // Block diagonal
        //---------------
       case 0:
        
        std::cout << "(diagonal version)" << std::endl;
        Flags::Run_identifier_string=
         "LSC FSI preconditioner (diagonal)";

        // Use block-diagonal preconditioner
        prec_pt->use_block_diagonal_version();
        break;
        
        
        // Retain fluid onto solid terms
        //------------------------------
       case 1:
        
        std::cout << "(retaining fluid onto solid terms)" << std::endl;
        Flags::Run_identifier_string=
         "LSC FSI preconditioner (retaining fluid onto solid terms)";

        // Choose preconditioner that retains fluid on solid terms
        prec_pt->use_block_triangular_version_with_fluid_on_solid();

        break;
        
        // Retain solid onto fluid terms
        //------------------------------
       case 2:
        
        std::cout << "(retaining solid onto fluid terms)" << std::endl;

        Flags::Run_identifier_string=
         "LSC FSI preconditioner (retaining solid onto fluid terms)";

        // Choose preconditioner that retains solid on fluid terms
        prec_pt->use_block_triangular_version_with_solid_on_fluid();
        break;
        
        // Error
        //------
       default:

        std::ostringstream error_stream; 
        error_stream << "Error: Wrong solver_sub_flag: " 
                     << solver_sub_flag << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;

       }

      // Pass preconditioner to iterative linear solver
      iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
     }

     break;


     // Error
     //======
    default:
     
     std::ostringstream error_stream; 
     error_stream << "Error: Wrong solver_flag: " 
                  << solver_flag << std::endl;
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     break;
    }
   
  }

 
 /// Destructor (empty)
 ~PreconditionedFSICollapsibleChannelProblem(){}



 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape.
 void actions_before_newton_convergence_check()
  {
   // Update mesh
   this->Bulk_mesh_pt->node_update();

   // Try to cast to IterativeLinearSolver
   IterativeLinearSolver* it_lin_solver_pt=
    dynamic_cast<IterativeLinearSolver*>(this->linear_solver_pt());

   // Open convergence history file
   if (it_lin_solver_pt!=0)
    {
     // Close it first
     it_lin_solver_pt->close_convergence_history_file_stream();

     // File name
     std::ostringstream some_stream;
     some_stream << this->Doc_info.directory() << "/convergence_history" 
                 << this->Doc_info.number() << "-" << this->Newton_iter 
                 << ".dat";

     // Zone title
     std::ostringstream some_stream2;
     some_stream2 << "Step: " << this->Doc_info.number() 
                  << "; Newton iteration:" << this->Newton_iter 
                  << "; " << Flags::Run_identifier_string;

     it_lin_solver_pt->open_convergence_history_file_stream(
      some_stream.str(),some_stream2.str());
    }

   // Increment counter
   this->Newton_iter++;
  }


 /// Actions before solve. Reset counter for number of Newton iterations
 void actions_before_newton_solve()
  {
   this->Newton_iter=0;
  }


 /// Update the problem after solve: Close convergence file for 
 /// iterative linear solver
 void actions_after_newton_solve()
  {
   // Try to cast to IterativeLinearSolver
   IterativeLinearSolver* it_lin_solver_pt=
    dynamic_cast<IterativeLinearSolver*>(this->linear_solver_pt());

   // Close convergence history file
   if (it_lin_solver_pt!=0)
    {
     it_lin_solver_pt->close_convergence_history_file_stream();
    }
  }

};









/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//============start_of_main====================================================
/// Driver code for a collapsible channel problem with FSI.
/// Presence of command line arguments indicates validation run with 
/// coarse resolution and small number of steps.
//=============================================================================
int main(int argc, char* argv[])
{ 
//#ifdef OOMPH_HAS_MPI
//  MPI_Helpers::init(argc,argv);
//#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 if (CommandLineArgs::Argc==1)
  {
   std::cout << "Using default settings for flags" << std::endl;
  }
 else if ((CommandLineArgs::Argc==10)||(CommandLineArgs::Argc==12))
  {
   /// Resolution factor
   Flags::Resolution_factor=atoi(argv[1]);

   /// Use displacement control (1) or not (0)
   Flags::Use_displ_control=atoi(argv[2]);

   /// Min. y coordinate for parameter study with displacement control
   Global_Physical_Variables::Yprescr_min=double(atof(argv[3]));

   /// Steady (1) or unsteady (0) run
   Flags::Steady_flag=atoi(argv[4]);

   /// Number of steps
   Flags::Nsteps=atoi(argv[5]);

   // Solver flag [0: direct; 1: exact; 2: simple; 3: Schur/Superlu; 
   // 4: Schur/Hypre]
   Flags::Solver_flag=atoi(argv[6]);
   
   // Solver sub flag [0: diag; 1: retain fluid on solid; 
   // 2: retain solid on fluid]
   Flags::Solver_sub_flag=atoi(argv[7]);

   /// Reynolds number
   Global_Physical_Variables::Re=double(atof(argv[8]));

   /// Womersley number
   Global_Physical_Variables::ReSt=Global_Physical_Variables::Re;

   /// FSI parameter
   Global_Physical_Variables::Q=double(atof(argv[9]));

   // Restart?
   if (CommandLineArgs::Argc==12)
    {
     // Name of restart file
     Flags::Restart_file_name=argv[10];
     
     // Jump in pressure 
     Global_Physical_Variables::P_step=double(atof(argv[11]));
    }
  }
 else
  {
   std::cout 
    << "\n\n\n\n\n"
    << "Wrong number of command line args: Specify none, six or eight:"
    << std::endl
    << "- resolution factor [1,2,...]" << std::endl
    << "- use_displ_control [0/1]" << std::endl
    << "- min. y-coordinate of control point when using displ control" 
    <<    std::endl
    << "- steady_flag [0/1]" << std::endl
    << "- number of steps " << std::endl
    << "- solver flag [0: direct; 1: exact; 2: simple; 3: Schur/Superlu; 4: Schur/Hypre]" << std::endl
    << "- solver sub flag [0: diag; 1: retain fluid on solid; 2: retain solid on fluid]" << std::endl
    << "- Reynolds number" << std::endl
    << "- FSI parameter Q" << std::endl
    << "- restart file name [optional] " << std::endl
    << "- jump in pressure P_step [optional] " << std::endl
    << "You specified " << CommandLineArgs::Argc-1 << " command line arg[s]" 
    << "\n\n\n\n\n" << std::endl;
   abort();
  }
 Flags::doc_flags();


 // Number of elements in the domain
 unsigned nup=4*Flags::Resolution_factor;
 unsigned ncollapsible=20*Flags::Resolution_factor;
 unsigned ndown=40*Flags::Resolution_factor;
 unsigned ny=4*Flags::Resolution_factor;
  
 
 // Length of the domain
 double lup=1.0;
 double lcollapsible=5.0;
 double ldown=10.0;
 double ly=1.0;
 
 // Use displacement control?
 bool displ_control=false;
 if (Flags::Use_displ_control==1) displ_control=true;
 
 // Steady run?
 bool steady_flag=false;
 if (Flags::Steady_flag==1) steady_flag=true;

 // Build the problem with QTaylorHoodElements
 PreconditionedFSICollapsibleChannelProblem
  <AlgebraicElement<QTaylorHoodElement<2> > > 
 problem(nup, ncollapsible, ndown, ny, 
         lup, lcollapsible, ldown, ly, displ_control,
         steady_flag,Flags::Solver_flag,Flags::Solver_sub_flag);
  
 if (Flags::Steady_flag)
  {
   problem.steady_run();
  }
 else
  {
   problem.unsteady_run();
  }
 
//#ifdef OOMPH_HAS_MPI
// MPI_Helpers::finalize();
//#endif


}//end of main



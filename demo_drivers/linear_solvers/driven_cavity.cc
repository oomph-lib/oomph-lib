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
//Driver for 2D rectangular driven cavity

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100;

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}
#endif


//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class RectangularDrivenCavityProblem : public Problem
{

public:


 ///  Constructor: Specify multiplier for number of element 
 /// rows/columns and solver flag.
 RectangularDrivenCavityProblem(const unsigned& element_multiplier, 
                                const bool& use_iterative_solver,
                                const bool& use_hypre_for_pressure,
                                const bool& use_hypre_for_momentum,
                                const bool& use_block_diagonal_for_momentum);
 

 /// Destructor
 ~RectangularDrivenCavityProblem()
  {
   // Kill oomph-lib iterative linear solver
   if (Solver_pt!=0) delete Solver_pt;
   
   // Kill preconditioner
   if (Prec_pt!=0) delete Prec_pt;
   
   // Kill inexact solver for P block
   if  (P_matrix_preconditioner_pt!=0) delete P_matrix_preconditioner_pt;
   
   // Kill inexact solver for F block
   if  (F_matrix_preconditioner_pt!=0) delete F_matrix_preconditioner_pt;

   // Kill mesh
   delete mesh_pt();

  };

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}


 ///  Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
 {
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along the other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
 } // end_of_actions_before_newton_solve

 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }


 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


private:

 /// oomph-lib iterative linear solver
 IterativeLinearSolver* Solver_pt;
 
 /// Preconditioner
 NavierStokesSchurComplementPreconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RectangularDrivenCavity problem: Specify 
/// multiplier for number of element rows/columns and solver flag.
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::RectangularDrivenCavityProblem(
 const unsigned& element_multiplier, 
 const bool& use_iterative_solver,
 const bool& use_hypre_for_pressure,
 const bool& use_hypre_for_momentum,
 const bool& use_block_diagonal_for_momentum)
{ 

 // Initialise pointer to oomph-lib iterative linear solver
 Solver_pt=0;
 
 // Initialise pointer to Preconditioner
 Prec_pt=0;
 
 // Initialise pointer to inexact solver for P block
 P_matrix_preconditioner_pt=0;
 
 // Initialise pointer to inexact solver for F block
 F_matrix_preconditioner_pt=0;
   

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10*element_multiplier;

 // # of elements in y-direction
 unsigned n_y=10*element_multiplier;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);


 // Setup iterative linear solver if required
 if (use_iterative_solver)
  {
   // Create oomph-lib iterative linear solver
   Solver_pt=new GMRES<CRDoubleMatrix>;
   
   // Set linear solver
   linear_solver_pt() = Solver_pt;
   
   // Set preconditioner
   Prec_pt=new NavierStokesSchurComplementPreconditioner(this);
   Prec_pt->set_navier_stokes_mesh(this->mesh_pt());

   Solver_pt->preconditioner_pt()=Prec_pt;
   
   // By default, the LSC Preconditioner uses SuperLU as
   // an exact preconditioner (i.e. a solver) for the
   // momentum and Schur complement blocks. 
   // Can overwrite this by passing pointers to 
   // other preconditioners that perform the (approximate)
   // solves of these blocks.
   
   // Create internal preconditioners used on Schur block
   //-----------------------------------------------------
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
   if (use_hypre_for_pressure)
    {
     P_matrix_preconditioner_pt = new HyprePreconditioner;
     
     // Set parameters for use as preconditioner on Poisson-type problem
     Hypre_default_settings::set_defaults_for_2D_poisson_problem(
      static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
     
     // Use Hypre for the Schur complement block
     Prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
     
     // Shut up!
     static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
      disable_doc_time();
    }
#endif    
#endif

   // Create internal preconditioners used on momentum block
   //--------------------------------------------------------
   if (use_block_diagonal_for_momentum)
    {
     F_matrix_preconditioner_pt = 
      new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
     if (use_hypre_for_momentum)
      {
       dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
        (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
        (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
      }
#endif
#endif
       // Use Hypre for momentum block 
       Prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
    }
   else
    {
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
     if (use_hypre_for_momentum)
      {
       F_matrix_preconditioner_pt = new HyprePreconditioner;
       
       // Shut up!
       static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt)->
        disable_doc_time();
       
       // Set parameters for use as preconditioner in for momentum 
       // block in Navier-Stokes problem
       Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
        static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt));
       
       // Use Hypre for momentum block 
       Prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
      }
#endif
#endif
    }

  }

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element. Optional command line
/// args specify multiplier for number of element rows/columns
/// and flag to indicate if iterative solver is used.
/// Multiplier and flag both default to 1. 
//=====================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Default assignemts for flags
 unsigned element_multiplier=1;
 bool use_iterative_solver=true;
 bool use_hypre_for_pressure=true;
 bool use_block_diagonal_for_momentum=false;
 bool use_hypre_for_momentum=true;
 bool do_cr=true;
 bool do_th=true;
 oomph_info << CommandLineArgs::Argc << std::endl;

 if (CommandLineArgs::Argc==1)
  {
   oomph_info << "No command line arguments; running self-test" 
              << std::endl;
  }
 else if (CommandLineArgs::Argc==7)
  {
   element_multiplier=unsigned(atoi(CommandLineArgs::Argv[1]));
   if (atoi(CommandLineArgs::Argv[2])==0)
    {
     use_iterative_solver=false;
    }
   else if (atoi(CommandLineArgs::Argv[2])==1)
    {
     use_iterative_solver=true;
    }
   else
    {
     oomph_info << "Second command line arg must be 0 or 1 for\n"
                << "don't/do use iterative solver, but is" 
                << CommandLineArgs::Argv[2] << std::endl;
     abort();
    }
   if (atoi(CommandLineArgs::Argv[3])==0)
    {
     use_hypre_for_momentum=false;
    }
   else if (atoi(CommandLineArgs::Argv[3])==1)
    {
     use_hypre_for_momentum=true;
    }
   else
    {
     oomph_info << "Third command line arg must be 0 or 1 for\n"
                << "don't/do use hypre for momentum block solver, but is" 
                << CommandLineArgs::Argv[3] << std::endl;
     abort();
    }
   if (atoi(CommandLineArgs::Argv[4])==0)
    {
     use_hypre_for_pressure=false;
    }
   else if (atoi(CommandLineArgs::Argv[4])==1)
    {
     use_hypre_for_pressure=true;
    }
   else
    {
     oomph_info << "Fourth command line arg must be 0 or 1 for\n"
                << "don't/do use hypre for pressure block solver, but is" 
                << CommandLineArgs::Argv[4] << std::endl;
     abort();
    }
   if (atoi(CommandLineArgs::Argv[5])==0)
    {
     do_cr=true;
     do_th=true;
    }
   else if (atoi(CommandLineArgs::Argv[5])==1)
    {
     do_cr=true;
     do_th=false;
    }
   else if (atoi(CommandLineArgs::Argv[5])==2)
    {
     do_cr=false;
     do_th=true;
    }
   else
    {
     oomph_info << "Fifth command line arg must be 0, 1 or 2 for\n"
                << "run both (0) or only CR (1) or TH (2) elements, but is" 
                << CommandLineArgs::Argv[5] << std::endl;
     abort();
    }
   if (atoi(CommandLineArgs::Argv[6])==0)
    {
     use_block_diagonal_for_momentum=false;
    }
   else if (atoi(CommandLineArgs::Argv[6])==1)
    {
     use_block_diagonal_for_momentum=true;
    }
   else
    {
     oomph_info << "Sixth command line arg must be 0 or 1 for\n"
                << "don't/do use block diagonal for pressure solver, but is" 
                << CommandLineArgs::Argv[6] << std::endl;
     abort();
    }
  }
 else
  {
   oomph_info << "Wrong number of command line arguments" << std::endl;
   oomph_info << "Enter none (for default) or six" << std::endl;
   oomph_info << "- multiplier for number of element rows/columns" 
              << std::endl;
   oomph_info << "- flag (0/1) for (don't/do) use iterative solver" 
              << std::endl;
   oomph_info << "- flag (0/1) for (don't/do) use hypre for pressure block" 
              << std::endl;
   oomph_info << "- flag (0/1) for (don't/do) use hypre for momentum block" 
              << std::endl;
   oomph_info << "- flag (0/1/2) run both (0) or only CR (1) or TH (2) "
              << "elements" << std::endl;
   oomph_info << "- flag (0/1) for (don't/do) use block diagonal for "
              << "momentum block" << std::endl;
   abort();
  } 

 oomph_info << "Running with element multiplier: " 
            <<  element_multiplier << std::endl;

 if (use_iterative_solver)
  {
   oomph_info << "Using iterative solver" << std::endl;
  }
 else
  {
   oomph_info << "Using direct solver" << std::endl;
  }

 if (use_hypre_for_pressure)
  {
   oomph_info << "Using hypre for pressure block" << std::endl;
  }
 else
  {
   oomph_info << "Using SuperLU for pressure block" << std::endl;
  }

 if (use_hypre_for_momentum)
  {
   oomph_info << "Using hypre for momentum block" << std::endl;
  }
 else
  {
   oomph_info << "Using SuperLU for momentum block" << std::endl;
  }

 if (use_block_diagonal_for_momentum)
  {
   oomph_info << "Using block diagonal for momentum block" << std::endl;
  }
 
 // Set up doc info
 // ---------------

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;

 // ---------------
 // end of Set up doc info


 // Doing QCrouzeixRaviartElements
 if (do_cr)
 {
  // Build the problem with QCrouzeixRaviartElements
  RectangularDrivenCavityProblem<QCrouzeixRaviartElement<2> > 
   problem(element_multiplier,use_iterative_solver,
           use_hypre_for_pressure,use_hypre_for_momentum,
           use_block_diagonal_for_momentum);
  cout << "Doing QCrouzeixRaviartElement<2>" << std::endl;
  
  // Solve the problem
  problem.newton_solve();
  
  // Outpt the solution
  problem.doc_solution(doc_info);

  // Step number
  doc_info.number()++;

 } // end of QCrouzeixRaviartElements

 // Doing QTaylorHoodElements
 if (do_th)
  {
   
   // Build the problem with QTaylorHoodElements
   RectangularDrivenCavityProblem<QTaylorHoodElement<2> > 
    problem(element_multiplier,use_iterative_solver,
            use_hypre_for_pressure,use_hypre_for_momentum,
            use_block_diagonal_for_momentum);
   cout << "Doing QTaylorHoodElement<2>" << std::endl;
   
   // Solve the problem
   problem.newton_solve();
   
   // Outpt the solution
   problem.doc_solution(doc_info);
   
   // Step number
   doc_info.number()++;
   
  } // end of QTaylorHoodElements
 
} // end_of_main











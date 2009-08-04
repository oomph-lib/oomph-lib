//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#include "trilinos_solver.h"


namespace oomph
{


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// functions for TrilinosAztecOOSolver class
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=============================================================================
/// Function which uses problem_pt's get_jacobian(...) function to
/// generate a linear system which is then solved. This function deletes
/// any existing internal data and then generates a new AztecOO solver.
//=============================================================================
 void TrilinosAztecOOSolver::solve(Problem* const &problem_pt,
                                   DoubleVector &solution)
 { 
  // clean up from previous solve
  clean_up_memory();

  // if we were previously using a problem based solve then we should delete
  // the oomphlib matrix as it was created during the last solve
  if (Using_problem_based_solve)
   {
    delete Oomph_matrix_pt;
    Oomph_matrix_pt = NULL;
   }

  // this is a problem based solve
  Using_problem_based_solve = true;
  
  // store the problem_pt
  Problem_pt = problem_pt;
  
  //Get oomph-lib Jacobian matrix and residual vector
  
  // record the start time
  double start_t = TimingHelpers::timer();

  // create the residual
  DoubleVector residual;

  // create the jacobian
  CRDoubleMatrix* cr_matrix_pt = new CRDoubleMatrix;
  Oomph_matrix_pt = cr_matrix_pt;
  problem_pt->get_jacobian(residual,*cr_matrix_pt);
  Distribution_pt->rebuild(residual.distribution_pt());

  // record the end time and compute the matrix setup time
  double end_t = TimingHelpers::timer();
  Jacobian_setup_time = end_t-start_t;
  if (this->doc_time())
   {
    oomph_info << "Time to generate Jacobian [sec]    : "
               << Jacobian_setup_time << std::endl;
   }
  
  // store the distribution of the solution vector
  if (!solution.distribution_setup())
    {
      solution.build(Distribution_pt,0.0);
    }
  LinearAlgebraDistribution solution_dist(solution.distribution_pt());

  // redistribute the distribution
  solution.redistribute(Distribution_pt);

  // continue solving using matrix based solve function
  solve(Oomph_matrix_pt, residual, solution);

  // return to the original distribution
  solution.redistribute(&solution_dist);
}



//=============================================================================
/// Function to solve the linear system defined by matrix_pt and rhs.
/// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
/// DistributedCRDoubleMatrix.
/// \b NOTE 2. This function will delete any existing internal data and 
/// generate a new AztecOO solver.
//=============================================================================
void TrilinosAztecOOSolver::solve(DoubleMatrixBase* const& matrix_pt,
                                  const DoubleVector &rhs,
                                  DoubleVector &result)
{

 // start the timer
 double start_t = TimingHelpers::timer();

#ifdef PARANOID
 // check that the matrix is square
 if (matrix_pt->nrow() != matrix_pt->ncol())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The matrix at matrix_pt must be square.";
   throw OomphLibError(error_message_stream.str(),
                       "TrilinosAztecOOSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);    
  }
 // check that the matrix and the rhs vector have the same nrow()
 if (matrix_pt->nrow() != rhs.nrow())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The matrix and the rhs vector must have the same number of rows.";
   throw OomphLibError(error_message_stream.str(),
                       "TrilinosAztecOOSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // if the matrix is distributable then it too should have the same 
 // communicator as the rhs vector and should not be distributed
 CRDoubleMatrix* cr_matrix_pt =  dynamic_cast<CRDoubleMatrix*>(matrix_pt);
 if (cr_matrix_pt != 0)
  {
   OomphCommunicator temp_comm(*rhs.distribution_pt()->communicator_pt());
   if (!(temp_comm == *cr_matrix_pt->distribution_pt()->communicator_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The matrix matrix_pt must have the same communicator as the vectors"
      << " rhs and result must have the same communicator";
     throw OomphLibError(error_message_stream.str(),
                         "TrilinosAztecOOSolver::solve()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 else
  {
   throw OomphLibError("Matrix must be of type CRDoubleMatrix",
                       "TrilinosAztecOOSolver::factorise()",
                       OOMPH_EXCEPTION_LOCATION);   
  }

 // if the result vector is setup then check it is not distributed and has 
 // the same communicator as the rhs vector
 if (result.distribution_setup())
  {
   if (!(*result.distribution_pt() == *rhs.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The result vector distribution has been setup; it must have the "
      << "same distribution as the rhs vector.";
     throw OomphLibError(error_message_stream.str(),
                         "TrilinosAztecOOSolver::solve()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }  
#endif

 // setup the solver
 solver_setup(matrix_pt);

 // set up trilinos vectors
 Epetra_Vector* epetra_r_pt;
 Epetra_Vector* epetra_z_pt;

 // create Epetra version of r
 TrilinosHelpers::create_epetra_vector(rhs,Epetra_map_pt,epetra_r_pt);

 // create an empty Epetra vector for z
 TrilinosHelpers::create_epetra_vector(Epetra_map_pt,epetra_z_pt);

 // solve the system
 solve_using_AztecOO(epetra_r_pt,epetra_z_pt);            
 
 // Copy result to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,result);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // delete solver data if required
 if (!Enable_resolve)
 {
  clean_up_memory();
 }

 // stop timers and compute solve time
 double end_t = TimingHelpers::timer();
 Linear_solver_solution_time = end_t-start_t;

 // output timings and info
 if (this->doc_time())
 {
  oomph_info << "Time for solve                        : "
             << Linear_solver_solution_time
             << "s" << std::endl;
 } 
}


//=============================================================================
/// Helper function for setting up the solver. Converts the oomph-lib 
/// matrices to Epetra matrices, sets up the preconditioner, creates the 
/// Trilinos Aztec00 solver and passes in the matrix, preconditioner and 
/// parameters.
//=============================================================================
void TrilinosAztecOOSolver::solver_setup(DoubleMatrixBase* const& matrix_pt)
{

 // clean up the memory
 //  - delete all except Oomph_matrix_pt, which may have been set in the 
 //    problem based solve
 clean_up_memory();

 // cast to CRDoubleMatrix
 // note cast check performed in matrix based solve(...) method
 CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

 // store the distribution
 // distribution of preconditioner is same as matrix
 Distribution_pt->rebuild(cast_matrix_pt->distribution_pt());

 // create the new solver
 AztecOO_solver_pt = new AztecOO();

 // if the preconditioner is an oomph-lib preconditioner then we set it up
 TrilinosPreconditionerBase* trilinos_prec_pt = 
  dynamic_cast<TrilinosPreconditionerBase* >(Preconditioner_pt);
 if (trilinos_prec_pt == 0)
  {
   if (Setup_preconditioner_before_solve)
    {
     // setup the preconditioner
     // start of prec setup
     double prec_setup_start_t = TimingHelpers::timer();
     Preconditioner_pt->setup(Problem_pt,matrix_pt);
     // start of prec setup
     double prec_setup_finish_t = TimingHelpers::timer();
     if (Doc_time)
      {
       double t_prec_setup = prec_setup_finish_t - prec_setup_start_t;
       oomph_info << "Time for preconditioner setup [sec]: "
                  << t_prec_setup << std::endl;
      }
#ifdef PARANOID
     if (*Preconditioner_pt->distribution_pt() != *Distribution_pt)
      {
       std::ostringstream error_message;
       error_message << "The oomph-lib preconditioner and the solver must "
                     << "have the same distribution";
       throw OomphLibError(error_message.str(),
                           "TrilinosAztecOOSolver::solver_setup()",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

   // wrap the oomphlib preconditioner in the Epetra_Operator derived
   // OoomphLibPreconditionerEpetraOperator to allow it to be passed to the
   // trilinos preconditioner
   Epetra_preconditioner_pt = 
    new OomphLibPreconditionerEpetraOperator(Preconditioner_pt);
  }


 // assemble the map
#ifdef OOMPH_HAS_MPI
 // MPI version
 Epetra_comm_pt = 
  new Epetra_MpiComm(Distribution_pt->communicator_pt()->mpi_comm());
 TrilinosHelpers::create_epetra_map(Distribution_pt,
                                    Epetra_comm_pt,
                                    Epetra_map_pt,Epetra_global_rows);
#else
 // Serial version
 Epetra_comm_pt = new Epetra_SerialComm;
 TrilinosHelpers::create_epetra_map(cast_matrix_pt->distribution_pt(),
                                    Epetra_comm_pt,
                                    Epetra_map_pt);
#endif

 // create the matrices
 double start_t_matrix = TimingHelpers::timer();
 Epetra_col_map_pt = new Epetra_Map(matrix_pt->ncol(),matrix_pt->ncol(),
                                    0,*Epetra_comm_pt);
 TrilinosHelpers::create_epetra_matrix(matrix_pt,
                                       Epetra_map_pt,
                                       Epetra_col_map_pt,
                                       Epetra_matrix_pt);   

 // record the end time and compute the matrix setup time
 double end_t_matrix = TimingHelpers::timer();
 if (trilinos_prec_pt == 0)
  {
   if (Using_problem_based_solve)
    {
     dynamic_cast<CRDoubleMatrix*>(Oomph_matrix_pt)->clear();
     delete Oomph_matrix_pt;
     Oomph_matrix_pt = NULL;   
    }
   
   // delete Oomph-lib matrix if requested
   else if (Delete_matrix)
    {
     dynamic_cast<CRDoubleMatrix*>(matrix_pt)->clear();
    }
  }

 // output times
 if (Doc_time)
 {
  oomph_info << "Time to generate Trilinos matrix      : "
             << double(end_t_matrix-start_t_matrix)
             << "s" << std::endl;
 }

 // set the matrix
 AztecOO_solver_pt->SetUserMatrix(Epetra_matrix_pt);

 //set the preconditioner
 if (trilinos_prec_pt == 0)
  {
   AztecOO_solver_pt->
    SetPrecOperator(Epetra_preconditioner_pt);  
  }  

#ifdef PARANOID
 // paranoid check the preconditioner exists
 if (Preconditioner_pt == 0)
  {
     std::ostringstream error_message;
     error_message << "Preconditioner_pt == 0. (Remember default "
                   << "preconditioner is IdentityPreconditioner)";
     throw OomphLibError(error_message.str(),
                         "TrilinosAztecOOSolver::solver_setup()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

 // if the preconditioner is a trilinos preconditioner then setup the
 // preconditioner
 if (trilinos_prec_pt != 0)
  {
   // only setup the preconditioner if required
   if (Setup_preconditioner_before_solve)
    {
     // start of prec setup
     double prec_setup_start_t = TimingHelpers::timer();
     
     // setup the preconditioner
     trilinos_prec_pt->setup(Problem_pt,matrix_pt,Epetra_matrix_pt);
     
     // set the preconditioner
     AztecOO_solver_pt->
      SetPrecOperator(trilinos_prec_pt->epetra_operator_pt());
     
     // start of prec setup
     double prec_setup_finish_t = TimingHelpers::timer();
     if (Doc_time)
      {
       double t_prec_setup = prec_setup_finish_t - prec_setup_start_t;
       oomph_info << "Time for preconditioner setup [sec]: "
                  << t_prec_setup << std::endl;
      }
    }

   // delete the oomph-matrix if required
   if (Using_problem_based_solve)
    {
     dynamic_cast<CRDoubleMatrix*>(Oomph_matrix_pt)->clear();
     delete Oomph_matrix_pt;
     Oomph_matrix_pt = NULL;   
    }
   
   // delete Oomph-lib matrix if requested
   else if (Delete_matrix)
    {
     dynamic_cast<CRDoubleMatrix*>(matrix_pt)->clear();
    }
  }
 
 // set solver options
 AztecOO_solver_pt->SetAztecOption(AZ_output, AZ_warnings);
 AztecOO_solver_pt->SetAztecOption(AZ_kspace, Max_iter);
 
 // set solver type
 switch (Solver_type)
  {
  case CG:
   AztecOO_solver_pt->SetAztecOption(AZ_solver, AZ_cg);
   break;
   
  case GMRES:
   AztecOO_solver_pt->SetAztecOption(AZ_solver, AZ_gmres);
   break;
   
  case BiCGStab:
   AztecOO_solver_pt->SetAztecOption(AZ_solver, AZ_bicgstab);
   break;
   
  default:
   std::ostringstream error_message;
   error_message << "Solver_type set to incorrect value. " 
                 << "Acceptable values are "
                 << CG << ", " << GMRES << " and " << BiCGStab
                 << ". Current value is " << Solver_type << ".";
   throw OomphLibError(error_message.str(),
                       "TrilinosAztecOOSolver::solver_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//=============================================================================
/// \short Function to resolve a linear system using the existing solver
/// data, allowing a solve with a new right hand side vector. This
/// function must be used after a call to solve(...) with
/// enable_resolve set to true.
//=============================================================================
void TrilinosAztecOOSolver::resolve(const DoubleVector &rhs,
                                    DoubleVector &solution)
{
 // start the timer
 double start_t = TimingHelpers::timer();

#ifdef PARANOID
 if (Epetra_matrix_pt->NumGlobalRows() != static_cast<int>(rhs.nrow()))
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same number "
                 << "of rows.\n"
                 << "The rhs vector has " << rhs.nrow() << " rows.\n"
                 << "The matrix has " << Epetra_matrix_pt->NumGlobalRows() 
                 << " rows.\n";   
   throw OomphLibError(error_message.str(),
                       "TrilinosAztecOOSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // set up trilinos vectors
 Epetra_Vector* epetra_r_pt=0;
 Epetra_Vector* epetra_z_pt=0;

 // create Epetra version of r
 TrilinosHelpers::create_epetra_vector(rhs,Epetra_map_pt,epetra_r_pt);

 // create an empty Epetra vector for z
 TrilinosHelpers::create_epetra_vector(Epetra_map_pt,epetra_z_pt);

 // solve the system
 solve_using_AztecOO(epetra_r_pt,epetra_z_pt);            

 // Copy result to z
 if (!solution.distribution_setup())
  {
   solution.build(rhs.distribution_pt(),0.0);
  }
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 double end_t = TimingHelpers::timer();
 Linear_solver_solution_time = end_t-start_t;

 // output timings and info
 if (this->doc_time())
 {
  oomph_info << "Time for resolve                        : "
             << Linear_solver_solution_time
             << "s" << std::endl;
 } 
}


//=============================================================================
/// Helper function performs the actual solve once the AztecOO
/// solver is set up (i.e. solver_setup() is called)
//=============================================================================
void TrilinosAztecOOSolver::solve_using_AztecOO(Epetra_Vector* &rhs_pt, 
                                                Epetra_Vector* &soln_pt)
{
#ifdef PARANOID
 // check the matrix and rhs are of consistent sizes
 if (AztecOO_solver_pt == 0 )
  {
   std::ostringstream error_message;
   error_message << "Solver must be called with solve(...) "
                 << "before resolve(...) to set it up.\n";
    throw OomphLibError(error_message.str(),
                        "TrilinosAztecOOSolver::solve_using_AztecOO()",
                        OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // set the vectors
 AztecOO_solver_pt->SetLHS(soln_pt);
 AztecOO_solver_pt->SetRHS(rhs_pt);
 
 // perform solve
 AztecOO_solver_pt->Iterate(Max_iter, Tolerance);


 // output iterations and final norm
 Iterations = AztecOO_solver_pt->NumIters();
 double norm = AztecOO_solver_pt->TrueResidual();
 oomph_info << "Linear solver iterations    : "
            << Iterations << std::endl;
 oomph_info << "Final Relative Residual Norm: " << norm << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


}

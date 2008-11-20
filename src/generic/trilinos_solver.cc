//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
                                   Vector<double> &solution)
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
#ifdef OOMPH_HAS_MPI
  double start_t = MPI_Wtime();
#else
  clock_t start_t = clock();
#endif

  // get the jacobian
#ifdef OOMPH_HAS_MPI
  DistributedVector<double> residual_dist;
  Vector<double> residual;
  if (Assemble_serial_jacobian)
   {
    unsigned n_row = problem_pt->ndof();
    residual.resize(n_row);
    CRDoubleMatrix* matrix_pt = new CRDoubleMatrix;
    problem_pt->get_jacobian(residual,*matrix_pt);
    Oomph_matrix_pt = matrix_pt; 
   }
  else
   {
    DistributedCRDoubleMatrix* matrix_pt = new DistributedCRDoubleMatrix;
    problem_pt->get_jacobian(residual_dist,*matrix_pt);
    Oomph_matrix_pt = matrix_pt; 
   }
#else
  CRDoubleMatrix* matrix_pt = new CRDoubleMatrix;
  unsigned n_row = problem_pt->ndof();
  Vector<double> residual(n_row);
  problem_pt->get_jacobian(residual,*matrix_pt);
  Oomph_matrix_pt = matrix_pt;
#endif

  // record the end time and compute the matrix setup time
#ifdef OOMPH_HAS_MPI
  double end_t = MPI_Wtime();
  Jacobian_setup_time = end_t-start_t;
#else
  clock_t end_t = clock();
  Jacobian_setup_time = double(end_t-start_t)/CLOCKS_PER_SEC;
#endif
  
  if (this->doc_time())
   {
    oomph_info << "Time to generate Jacobian [sec]           : "
               << Jacobian_setup_time << std::endl;
   }
  
  // continue solving using matrix based solve function
#ifdef OOMPH_HAS_MPI
  if (Assemble_serial_jacobian)
   {
    solve(Oomph_matrix_pt, residual, solution);
   }
  else
   {
    solve(Oomph_matrix_pt, residual_dist, solution);
   }
#else
  solve(Oomph_matrix_pt, residual, solution);
#endif
}



//=============================================================================
/// Function to solve the linear system defined by matrix_pt and rhs.
/// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
/// DistributedCRDoubleMatrix.
/// \b NOTE 2. This function will delete any existing internal data and 
/// generate a new AztecOO solver.
//=============================================================================
void TrilinosAztecOOSolver::solve(DoubleMatrixBase* const& matrix_pt,
                                  const Vector<double> &rhs,
                                  Vector<double> &solution)
{

 // start the timer
#ifdef OOMPH_HAS_MPI
 double start_t = MPI_Wtime();
#else
 clock_t start_t = clock();
#endif

#ifdef PARANOID
 if (matrix_pt->nrow() != rhs.size())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same number "
                 << "of rows.\n"
                 << "The rhs vector has " << rhs.size() << " rows.\n"
                 << "The matrix has " << matrix_pt->nrow() << " rows.\n";   
   throw OomphLibError(error_message.str(),
                       "TrilinosAztecOOSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
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
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // delete solver data if required
 if (!Enable_resolve)
 {
  clean_up_memory();
 }

#ifdef OOMPH_HAS_MPI
  double end_t = MPI_Wtime();
  Linear_solver_solution_time = end_t-start_t;
#else
  clock_t end_t = clock();
  Linear_solver_solution_time = double(end_t-start_t)/CLOCKS_PER_SEC;
#endif

 // output timings and info
 if (this->doc_time())
 {
  oomph_info << "Time spent in solver                      : "
             << Linear_solver_solution_time
             << "s" << std::endl;
 } 
}

#ifdef OOMPH_HAS_MPI
//=============================================================================
/// Function to solve the linear system defined by matrix_pt and rhs\n.
/// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
/// DistributedCRDoubleMatrix.
/// \b NOTE 2. This function will delete any existing internal data and 
/// generate a new AztecOO solver.
/// \b NOTE 3. The vector rhs must have the same distribution as the matrix
//=============================================================================
void TrilinosAztecOOSolver::solve(DoubleMatrixBase* const& matrix_pt,
                                  const DistributedVector<double> &rhs,
                                  Vector<double> &solution)
{
 // start the timer
 double start_t = MPI_Wtime();

 // setup the solver
 solver_setup(matrix_pt);

#ifdef PARANOID
 if (Solver_distribution != rhs.distribution())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
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

 // copy the results to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // delete solver data if required
 if (!Enable_resolve)
 {
  clean_up_memory();
 }

 // end of solve time
 double end_t = MPI_Wtime();
 Linear_solver_solution_time = end_t-start_t;
 
 // output timings and info
 if (this->doc_time())
  {
   oomph_info << "Time for solve                            : "
              << Linear_solver_solution_time
              << "s" << std::endl;
  } 
}

//=============================================================================
/// Function to solve the linear system defined by matrix_pt and rhs\n.
/// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
/// DistributedCRDoubleMatrix.
/// \b NOTE 2. This function will delete any existing internal data and 
/// generate a new AztecOO solver. 
/// \b NOTE 3. The vectors rhs and solution must have the same distribution 
/// as the matrix
//=============================================================================
void TrilinosAztecOOSolver::solve(DoubleMatrixBase* const& matrix_pt,
                                  const DistributedVector<double> &rhs,
                                  DistributedVector<double> &solution)
{
 // start time
 double start_t = MPI_Wtime();

 // setup the solver
 solver_setup(matrix_pt);

#ifdef PARANOID
 if (Solver_distribution != rhs.distribution())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
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

 // copy the results to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,Solver_distribution,
                                          solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // delete solver data if required
 if (!Enable_resolve)
 {
  clean_up_memory();
 }

 // end of solve time
 double end_t = MPI_Wtime();
 Linear_solver_solution_time = end_t-start_t;
 
 // output timings and info
 if (this->doc_time())
  {
   oomph_info << "Time spent in linear solver               : "
              << Linear_solver_solution_time
              << "s" << std::endl;
  } 
}
#endif

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
 //=====================================================================
 clean_up_memory();

 // begin by setting up the Epetra Matrix
 //======================================

#ifdef OOMPH_HAS_MPI
 double start_t_setup = MPI_Wtime();
#else
 clock_t start_t_setup = clock();
#endif



 // if this is a Distributed build then attempt to build Trilinos Matrix from
 // a CRDoubleMatrix or a DistributedCRDoubleMatrix
#ifdef OOMPH_HAS_MPI
 bool cast_failed = true;

 // start by determining the distribution of this preconditioner
 // if DistributedCRDoubleMatrix then the distribution is the same as the 
 // preconditioner.
 // if CRDoubleMatrix then the distribution is the uniform distribution

 // first try DistributedCRDoubleMatrix
 DistributedCRDoubleMatrix* cast_dist_matrix_pt = 
  dynamic_cast<DistributedCRDoubleMatrix*>(matrix_pt);
 if (cast_dist_matrix_pt!=0)
  {
   #ifdef PARANOID
   // check the matrix is square
   if ( cast_dist_matrix_pt->nrow_global() != matrix_pt->ncol() )
    {
  std::ostringstream error_message;
  error_message << "TrilinosAztecOOSolver require a square matrix. "
                << "Matrix is " << matrix_pt->nrow()
                << " by " << matrix_pt->ncol() << std::endl;
  throw OomphLibError(error_message.str(),
                      "TrilinosAztecOOSolver::solver_setup()",
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif
   // the distribution of this preconditioner is the same as the matrix
   Solver_distribution = cast_dist_matrix_pt->distribution();
   cast_failed = false;
  }

 // if the cast failed try cast to CRDoubleMatrix
 if (cast_failed)
 {
  CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  if (cast_matrix_pt != 0)
  {
#ifdef PARANOID
   // check the matrix is square
   if ( matrix_pt->nrow() != matrix_pt->ncol() )
    {
     std::ostringstream error_message;
     error_message << "TrilinosAztecOOSolver require a square matrix. "
                   << "Matrix is " << matrix_pt->nrow()
                   << " by " << matrix_pt->ncol() << std::endl;
     throw OomphLibError(error_message.str(),
                         "TrilinosAztecOOSolver::solver_setup()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // use default distribution
   // NOTE: USE of MPI_COMM_WORLD
   Solver_distribution.distribute(MPI_COMM_WORLD,
                                  cast_matrix_pt->nrow());
   cast_failed = false;
  }
 }
 
#ifdef PARANOID
 // check the matrix is a 
 if (cast_failed)
  {
   std::ostringstream error_message;
   error_message << "TrilinosSolver only work with "
                 << "DistributedCRDoubleMatrix matrices" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosSolver::solver_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // assemble the map
 Epetra_comm_pt = new Epetra_MpiComm(MPI_COMM_WORLD);
 TrilinosHelpers::create_epetra_map(Solver_distribution,
                                    Epetra_comm_pt,
                                    Epetra_map_pt,Epetra_global_rows);

#else
 // Serial Version

#ifdef PARANOID
 // check the matrix is square
 if ( matrix_pt->nrow() != matrix_pt->ncol() )
 {
  std::ostringstream error_message;
  error_message << "TrilinosAztecOOSolver require a square matrix. "
                << "Matrix is " << matrix_pt->nrow()
                << " by " << matrix_pt->ncol() << std::endl;
  throw OomphLibError(error_message.str(),
                      "TrilinosAztecOOSolver::solver_setup()",
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // cast to CRDoubleMatrix*
 CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
#ifdef PARANOID
 // check the matrix is a CRDoubleMatrix
 if (cast_matrix_pt==0)
  {
   std::ostringstream error_message;
   error_message << "TrilinosSolver only work with "
                 << "CRDoubleMatrix matrices" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosSolver::solver_setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // setup the epetra matrix
 Epetra_comm_pt = new Epetra_SerialComm;
 TrilinosHelpers::create_epetra_map(cast_matrix_pt->nrow(),
                                    Epetra_comm_pt,
                                    Epetra_map_pt);
#endif
 
 // create the matrices
 TrilinosHelpers::create_epetra_matrix(matrix_pt,
                                       Epetra_map_pt,
                                       Epetra_matrix_pt);   

  // record the end time and compute the matrix setup time
#ifdef OOMPH_HAS_MPI
  double end_t_setup = MPI_Wtime();
#else
 clock_t end_t_setup = clock();
#endif

 // output times
 if (Doc_time)
 {
  oomph_info << "Time to generate Trilinos matrix          : "
#ifdef OOMPH_HAS_MPI
             << double(end_t_setup-start_t_setup)
#else
             << double(end_t_setup-start_t_setup)/CLOCKS_PER_SEC
#endif
             << "s" << std::endl;
 }

 // setup the Aztec00 solver and Preconditioner
 //============================================

 // create the new solver
 AztecOO_solver_pt = new AztecOO();

 // set the matrix
 AztecOO_solver_pt->SetUserMatrix(Epetra_matrix_pt);
 
 // if the preconditioner is a trilinos preconditioner then set the
 // Epetra_Operator
 TrilinosPreconditionerBase* trilinos_prec_pt = 
  dynamic_cast<TrilinosPreconditionerBase* >(Preconditioner_pt);
 if (trilinos_prec_pt != 0)
  {
   // setup the preconditioner
   trilinos_prec_pt->setup(Problem_pt,Oomph_matrix_pt,Epetra_matrix_pt);
   
   // set the preconditioner
   AztecOO_solver_pt->
    SetPrecOperator(trilinos_prec_pt->epetra_operator_pt());
  }
 
 // otherwise the preconditioner is an oomph-lib preconditioner and we wrap
 // it in the Epetra-Operator derived object 
 // OomphLibPreconditionerEpetraOperator
 else
  {
   // setup the preconditioner
   Preconditioner_pt->setup(Problem_pt,matrix_pt);
#ifdef OOMPH_HAS_MPI
   if (!Preconditioner_pt->distribution().setup())
    {
     Preconditioner_pt->distribution() = Solver_distribution;
    }
#ifdef PARANOID
   if (Preconditioner_pt->distribution() != Solver_distribution)
    {
     std::ostringstream error_message;
     error_message << "The oomph-lib preconditioner and the solver must "
                   << "have the same distribution";
     throw OomphLibError(error_message.str(),
                         "TrilinosAztecOOSolver::solver_setup()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
#endif
   
   // wrap the oomphlib preconditioner in the Epetra_Operator derived
   // OoomphLibPreconditionerEpetraOperator to allow it to be passed to the
   // trilinos preconditioner
#ifdef OOMPH_HAS_MPI
   Epetra_preconditioner_pt = 
    new OomphLibPreconditionerEpetraOperator(Preconditioner_pt);
#else
   Epetra_preconditioner_pt = 
    new OomphLibPreconditionerEpetraOperator(Preconditioner_pt,
                                             matrix_pt->nrow());
#endif
   AztecOO_solver_pt->
    SetPrecOperator(Epetra_preconditioner_pt);    
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
 
 
 // Delete the oomph-lib matrix if required
 //========================================

 // if the solve function was initially called via the Problem based interface
 // set flags so that the oomph-lib copy of the Jacobian matrix is deleted once
 // the Trilinos matrix is created

 if (Using_problem_based_solve)
  {
   dynamic_cast<CRDoubleMatrix*>(Oomph_matrix_pt)->clean_up_memory();
   delete Oomph_matrix_pt;
   Oomph_matrix_pt = NULL;   
  }
   
 // delete Oomph-lib matrix if requested
 else if (Delete_matrix)
  {
   dynamic_cast<CRDoubleMatrix*>(matrix_pt)->clean_up_memory();
  }
}

//=============================================================================
/// \short Function to resolve a linear system using the existing solver
/// data, allowing a solve with a new right hand side vector. This
/// function must be used after a call to solve(...) with
/// enable_resolve set to true.
//=============================================================================
void TrilinosAztecOOSolver::resolve(const Vector<double> &rhs,
                                    Vector<double> &solution)
{
 // start the timer
#ifdef OOMPH_HAS_MPI
 double start_t = MPI_Wtime();
#else
 clock_t start_t = clock();
#endif

#ifdef PARANOID
 if (Epetra_matrix_pt->NumGlobalRows() != static_cast<int>(rhs.size()))
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same number "
                 << "of rows.\n"
                 << "The rhs vector has " << rhs.size() << " rows.\n"
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
 
#ifdef OOMPH_HAS_MPI
 // Copy result to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);
#else
 // Copy result to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);
#endif

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

#ifdef OOMPH_HAS_MPI
  double end_t = MPI_Wtime();
  Linear_solver_solution_time = end_t-start_t;
#else
  clock_t end_t = clock();
  Linear_solver_solution_time = double(end_t-start_t)/CLOCKS_PER_SEC;
#endif

 // output timings and info
 if (this->doc_time())
 {
  oomph_info << "Time for resolve                            : "
             << Linear_solver_solution_time
             << "s" << std::endl;
 } 
}

#ifdef OOMPH_HAS_MPI
//=============================================================================
/// Function to resolve a linear system using the existing solver
/// data, allowing a solve with a new right hand side vector. This
/// function must be used after a call to solve(...) with
/// enable_resolve set to true.
//=============================================================================
void TrilinosAztecOOSolver::resolve(const DistributedVector<double> &rhs,
                                    Vector<double> &solution)
{
 // start the timer
 double start_t = MPI_Wtime();

#ifdef PARANOID
 if (Solver_distribution != rhs.distribution())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
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

 // copy the results to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // end of solve time
 double end_t = MPI_Wtime();
 Linear_solver_solution_time = end_t-start_t;
 
 // output timings and info
 if (this->doc_time())
  {
   oomph_info << "Time for resolve                            : "
              << Linear_solver_solution_time
              << "s" << std::endl;
  } 
}

//=============================================================================
/// \short Function to resolve a linear system using the existing solver
/// data, allowing a solve with a new right hand side vector. This
/// function must be used after a call to solve(...) with
/// enable_resolve set to true.
//=============================================================================
void TrilinosAztecOOSolver::resolve(const DistributedVector<double> &rhs,
                                    DistributedVector<double> &solution)
{
 // start time
 double start_t = MPI_Wtime();

#ifdef PARANOID
 if (Solver_distribution != rhs.distribution())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
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

 // copy the results to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,Solver_distribution,
                                          solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 // end of solve time
 double end_t = MPI_Wtime();
 Linear_solver_solution_time = end_t-start_t;
 
 // output timings and info
 if (this->doc_time())
  {
   oomph_info << "Time for resolve                            : "
              << Linear_solver_solution_time
              << "s" << std::endl;
  } 
}
#endif


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

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
   oomph_info << "RAYRAY: solve() begin" << std::endl; 
   


//   MemoryUsage::doc_memory_usage("start of TrilinosAztecOOSolver::solve");
//   MemoryUsage::insert_comment_to_continous_top(
//    "start of TrilinosAztecOOSolver::solve");

  // clean up from previous solve
  clean_up_memory();

  // if we were previously using a problem based solve then we should delete
  // the oomphlib matrix as it was created during the last solve
  oomph_info << "RAYRAY: Using_problem_based_solve: " << Using_problem_based_solve << std::endl; 
  
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

//   MemoryUsage::doc_memory_usage("start of get_jacobian()");
//   MemoryUsage::insert_comment_to_continous_top("start of get_jacobian()");

  // create the residual
  DoubleVector residual;

  // create the jacobian
  CRDoubleMatrix* cr_matrix_pt = new CRDoubleMatrix;
  Oomph_matrix_pt = cr_matrix_pt;
  oomph_info << "RAYRAY: About to get Jacobian" << std::endl;
  problem_pt->get_jacobian(residual,*cr_matrix_pt);
  oomph_info << "RAYRAY: Got Jacobian" << std::endl;

  unsigned rnrow = cr_matrix_pt->nrow();
  unsigned rncol = cr_matrix_pt->ncol();
  unsigned rnrow_local = cr_matrix_pt->nrow_local();
  unsigned rnnz = cr_matrix_pt->nnz();

  oomph_info << "rnrow = " << rnrow << std::endl; 
  oomph_info << "rncol = " << rncol << std::endl; 
  oomph_info << "rnrow_local = " << rnrow_local << std::endl; 
  oomph_info << "rnnz = " << rnnz << std::endl; 
  

  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
  // my rank and number of processors. 
  // This is used later for putting the data.
  const unsigned my_rank = comm_pt->my_rank();
  const unsigned nproc = comm_pt->nproc();

  if(Dump_matrices)
  {
    oomph_info << "About to output jac and res" << std::endl;

    std::ostringstream jac_ss;
    jac_ss << "raw_lin_system/jac" 
           << Tetgen_number << "_np" << nproc 
                            << "r" << my_rank;
    cr_matrix_pt->sparse_indexed_output(jac_ss.str(),15,true);

    std::ostringstream residual_ss;
    residual_ss << "raw_lin_system/residual" 
                << Tetgen_number << "_np" << nproc 
                                 << "r" << my_rank;
    residual.output_local_values_with_offset(residual_ss.str(),15);

    oomph_info << "Done output of jac and res" << std::endl;
  }
  
  this->build_distribution(residual.distribution_pt());

  // record the end time and compute the matrix setup time
  double end_t = TimingHelpers::timer();
  Jacobian_setup_time = end_t-start_t;
  if (this->Doc_time)
   {
    oomph_info << "Time to generate Jacobian [sec]    : "
               << Jacobian_setup_time << std::endl;
   }


  // RAYRAY new code, to read in the permuted matrix.
  // How to do this... the files are (on two cores)
  //
  // Jacobian:
  // (for rank 0)
  // jac13_np2r0_val
  // jac13_np2r0_row
  // jac13_np2r0_col
  //
  // (for rank 1)
  // jac13_np2r1_val
  // jac13_np2r1_row
  // jac13_np2r1_col
  //
  //
  // Residual:
  // res13_np2r0
  // res13_np2r0
  //
  // The jac13 is stored in Replaced_mat_str.
  // The res13 is stored in Replaced_res_str.
  // If Use_replaced_mat_res is true, then we assume that the files exist
  // and use them.
  if(Use_replacement_mat_res)
  {
    // Recall that we have my_rank and nproc
    std::ostringstream np_rank_oss;
    np_rank_oss << "np" << nproc << "r" << my_rank;

    ////////////////////////////////////////////////////////////////////////
    // create the val oss
    std::ostringstream val_fname_oss;
    val_fname_oss << "jac" << Tetgen_number << "_" 
                  << np_rank_oss.str() << "_val";

    // create the col oss
    std::ostringstream col_fname_oss;
    col_fname_oss << "jac" << Tetgen_number << "_"
                  << np_rank_oss.str() << "_col";

    // create the row oss
    std::ostringstream row_fname_oss;
    row_fname_oss << "jac" << Tetgen_number << "_" 
                  << np_rank_oss.str() << "_row";
    ////////////////////////////////////////////////////////////////////////
    std::string val_nlines_fname = val_fname_oss.str() + "_nlines";
    std::string col_nlines_fname = col_fname_oss.str() + "_nlines";
    std::string row_nlines_fname = row_fname_oss.str() + "_nlines";
    ////////////////////////////////////////////////////////////////////////

    // We now read in the nlines
    unsigned val_nlines = 0;
    unsigned col_nlines = 0;
    unsigned row_nlines = 0;

    // Get val_nlines //////////////////////////////////////////////////////
    std::string full_val_nlines_fname = "raw_lin_system/"+val_nlines_fname;
    std::ifstream val_nline_file(full_val_nlines_fname.c_str());
    if (!val_nline_file)
    {
      std::ostringstream error_message_stream;
      error_message_stream
       << "There was a problem opening file "
       << full_val_nlines_fname  << std::endl;
       throw OomphLibError(error_message_stream.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      val_nline_file >> val_nlines;
    }

    // Get col_nlines //////////////////////////////////////////////////////
    std::string full_col_nlines_fname = "raw_lin_system/"+col_nlines_fname;
    std::ifstream col_nline_file(full_col_nlines_fname.c_str());
    if (!col_nline_file)
    {
      std::ostringstream error_message_stream;
      error_message_stream
       << "There was a problem opening file "
       << full_col_nlines_fname  << std::endl;
       throw OomphLibError(error_message_stream.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      col_nline_file >> col_nlines;
    }

    // Get row_nlines //////////////////////////////////////////////////////
    std::string full_row_nlines_fname = "raw_lin_system/"+row_nlines_fname;
    std::ifstream row_nline_file(full_row_nlines_fname.c_str());
    if (!row_nline_file)
    {
      std::ostringstream error_message_stream;
      error_message_stream
       << "There was a problem opening file "
       << full_row_nlines_fname  << std::endl;
       throw OomphLibError(error_message_stream.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      row_nline_file >> row_nlines;
    }

    oomph_info << "val_nlines: " << val_nlines << std::endl; 
    oomph_info << "col_nlines: " << col_nlines << std::endl; 
    oomph_info << "row_nlines: " << row_nlines << std::endl; 
   
    // Now we load in the values in file jac1_np<NP>r<R>_val, this is in
    // e.g. jac1_np2r0_val. This is located in val_fname_oss 
    std::string full_val_fname = "raw_lin_system/"+val_fname_oss.str();
    std::ifstream val_file(full_val_fname.c_str());
    double* val_array = new double[val_nlines];
    for (unsigned val_i = 0; val_i < val_nlines; val_i++) 
    {
      val_file >> val_array[val_i];
    }

//    oomph_info << "read in values are:" << std::endl; 
//    
//    for (unsigned val_i = 0; val_i < val_nlines; val_i++) 
//    {
//      oomph_info << val_array[val_i] << std::endl; 
//    }
//    pause("done!"); 


    // Now we load in the values in file jac1_np<NP>r<R>_val, this is in
    // e.g. jac1_np2r0_val. This is located in col_fname_oss 
    std::string full_col_fname = "raw_lin_system/"+col_fname_oss.str();
    std::ifstream col_file(full_col_fname.c_str());
    int* col_array = new int[col_nlines];
    for (unsigned col_i = 0; col_i < col_nlines; col_i++) 
    {
      col_file >> col_array[col_i];
    }

//    oomph_info << "read in col are:" << std::endl; 
//    for (unsigned col_i = 0; col_i < col_nlines; col_i++) 
//    {
//      oomph_info << col_array[col_i] << std::endl; 
//    }
//    pause("done!"); 
    


    // Now we load in the values in file jac1_np<NP>r<R>_val, this is in
    // e.g. jac1_np2r0_val. This is located in row_fname_oss 
    std::string full_row_fname = "raw_lin_system/"+row_fname_oss.str();
    std::ifstream row_file(full_row_fname.c_str());
    int* row_array = new int[row_nlines];
    for (unsigned row_i = 0; row_i < row_nlines; row_i++) 
    {
      row_file >> row_array[row_i];
    }

//    oomph_info << "read in row are:" << std::endl; 
//    for (unsigned row_i = 0; row_i < row_nlines; row_i++) 
//    {
//      oomph_info << row_array[row_i] << std::endl; 
//    }
//    pause("done!"); 

    cr_matrix_pt->build_without_copy(rncol,val_nlines,
                                     val_array,col_array,row_array);

    



    // Now build the residual this is stored in
    // res<TETNUM>_np<NP>r<RANK>
    std::ostringstream res_fname_oss;
    res_fname_oss << "res" << Tetgen_number << "_"
                  << np_rank_oss.str();

    std::string full_res_fname = "raw_lin_system/"+res_fname_oss.str();
    std::ifstream res_file(full_res_fname.c_str());

    double* res_val_pt = residual.values_pt();
    for (unsigned res_i = 0; res_i < rnrow_local; res_i++)
    {
      res_file >> res_val_pt[res_i];
    }
    ///// re-output this matrix to make sure it's correct?
    if(Dump_replacement)
    {
    oomph_info << "About to output jac and res AGAIN" << std::endl;

    std::ostringstream jac_ss;
    jac_ss << "raw_lin_system/AAjac"
           << Tetgen_number << "_np" << nproc 
                            << "r" << my_rank;
    cr_matrix_pt->sparse_indexed_output(jac_ss.str(),15,true);

    std::ostringstream residual_ss;
    residual_ss << "raw_lin_system/AAresidual" 
                << Tetgen_number << "_np" << nproc 
                                 << "r" << my_rank;
    residual.output_local_values_with_offset(residual_ss.str(),15);

    oomph_info << "AADone output of jac and res" << std::endl;
    
    pause("Hello again!"); 
    }


  }


//   MemoryUsage::doc_memory_usage("after get_jacobian() in trilinos solver");
//   MemoryUsage::insert_comment_to_continous_top(
//    "after get_jacobian() in trilinos solver");

  // store the distribution of the solution vector
  if (!solution.built())
    {
     solution.build(this->distribution_pt(),0.0);
    }
  LinearAlgebraDistribution solution_dist(solution.distribution_pt());

  // redistribute the distribution
  solution.redistribute(this->distribution_pt());

//   MemoryUsage::doc_memory_usage("before trilinos solve");
//   MemoryUsage::insert_comment_to_continous_top("before trilinos solve ");

  oomph_info << "RAYRAY: About to solve" << std::endl; 
  
  // continue solving using matrix based solve function
  solve(Oomph_matrix_pt, residual, solution);
  oomph_info << "RAYRAY: End of solve" << std::endl; 
  


//   MemoryUsage::doc_memory_usage("after trilinos solve");
//   MemoryUsage::insert_comment_to_continous_top("after trilinos solve ");

  // return to the original distribution
  solution.redistribute(&solution_dist);

//   MemoryUsage::doc_memory_usage("end of TrilinosAztecOOSolver::solve");
//   MemoryUsage::insert_comment_to_continous_top(
//    "end of TrilinosAztecOOSolver::solve");


   oomph_info << "RAYRAY: solve() end" << std::endl; 
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
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that the matrix and the rhs vector have the same nrow()
 if (matrix_pt->nrow() != rhs.nrow())
  {
   std::ostringstream error_message_stream;
   error_message_stream
    << "The matrix and the rhs vector must have the same number of rows.";
   throw OomphLibError(error_message_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
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
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 else
  {
   throw OomphLibError("Matrix must be of type CRDoubleMatrix",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // if the result vector is setup then check it is not distributed and has
 // the same communicator as the rhs vector
 if (result.built())
  {
   if (!(*result.distribution_pt() == *rhs.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream
      << "The result vector distribution has been setup; it must have the "
      << "same distribution as the rhs vector.";
     throw OomphLibError(error_message_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // build the result if not built
 if (!result.built())
  {
   result.build(rhs.distribution_pt(),0.0);
  }

 // setup the solver
 solver_setup(matrix_pt);

 // create Epetra version of r
 Epetra_Vector* epetra_r_pt = TrilinosEpetraHelpers::
  create_distributed_epetra_vector(rhs);

 // create an empty Epetra vector for z
 Epetra_Vector* epetra_z_pt = TrilinosEpetraHelpers::
  create_distributed_epetra_vector(result);

 double start_t_trilinos = TimingHelpers::timer();

 // solve the system
 solve_using_AztecOO(epetra_r_pt,epetra_z_pt);

 double end_t_trilinos = TimingHelpers::timer();
 if (this->Doc_time)
 {
  oomph_info << "Time for trilinos solve itself                 : "
             << end_t_trilinos-start_t_trilinos
             << "s" << std::endl;
 }
 // Copy result to z
 TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_z_pt,result);

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
 if (this->Doc_time)
 {
  oomph_info << "Time for complete trilinos solve                  : "
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
 this->build_distribution(cast_matrix_pt->distribution_pt());

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
     Preconditioner_pt->setup(matrix_pt);
     // start of prec setup
     double prec_setup_finish_t = TimingHelpers::timer();
     if (Doc_time)
      {
       double t_prec_setup = prec_setup_finish_t - prec_setup_start_t;
       oomph_info << "Time for preconditioner setup [sec]: "
                  << t_prec_setup << std::endl;
      }
#ifdef PARANOID
     if (*Preconditioner_pt->distribution_pt() != *this->distribution_pt())
      {
       std::ostringstream error_message;
       error_message << "The oomph-lib preconditioner and the solver must "
                     << "have the same distribution";
       throw OomphLibError(error_message.str(),
                           OOMPH_CURRENT_FUNCTION,
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

 // create the matrix
 double start_t_matrix = TimingHelpers::timer();
 if (Use_aztecoo_workaround_for_epetra_matrix_setup)
  {
   Epetra_matrix_pt = TrilinosEpetraHelpers
    ::create_distributed_epetra_matrix_for_aztecoo(cast_matrix_pt);
  }
 else
  {
   Epetra_matrix_pt = TrilinosEpetraHelpers
    ::create_distributed_epetra_matrix(cast_matrix_pt,
                                       cast_matrix_pt->distribution_pt());
  }

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
                         OOMPH_CURRENT_FUNCTION,
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
     trilinos_prec_pt->set_matrix_pt(matrix_pt);
     trilinos_prec_pt->setup(Epetra_matrix_pt);


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
 if (Doc_time)
  {
   AztecOO_solver_pt->SetAztecOption(AZ_output, AZ_warnings);
  }
 else
  {
   AztecOO_solver_pt->SetAztecOption(AZ_output, AZ_none);
  }
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
                       OOMPH_CURRENT_FUNCTION,
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
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

// if the result vector is setup then check it is not distributed and has
 // the same communicator as the rhs vector
 if (solution.built())
  {
   if (!(*solution.distribution_pt() == *rhs.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream
      << "The result vector distribution has been setup; it must have the "
      << "same distribution as the rhs vector.";
     throw OomphLibError(error_message_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // build the result if not built
 if (!solution.built())
  {
   solution.build(rhs.distribution_pt(),0.0);
  }


 // create Epetra version of r
 Epetra_Vector* epetra_r_pt = TrilinosEpetraHelpers::
  create_distributed_epetra_vector(rhs);

 // create an empty Epetra vector for z
 Epetra_Vector* epetra_z_pt = TrilinosEpetraHelpers::
  create_distributed_epetra_vector(solution);

 // solve the system
 solve_using_AztecOO(epetra_r_pt,epetra_z_pt);

 // Copy result to z
 if (!solution.distribution_built())
  {
   solution.build(rhs.distribution_pt(),0.0);
  }
 TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_z_pt,solution);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;

 double end_t = TimingHelpers::timer();
 Linear_solver_solution_time = end_t-start_t;

 // output timings and info
 if (this->Doc_time)
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
                        OOMPH_CURRENT_FUNCTION,
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
  if (Doc_time)
   {
    double norm = AztecOO_solver_pt->TrueResidual();
    oomph_info << "Linear solver iterations    : "
               << Iterations << std::endl;
    oomph_info << "Final Relative Residual Norm: " << norm << std::endl;
   }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


}

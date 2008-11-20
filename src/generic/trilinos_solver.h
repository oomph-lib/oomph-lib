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
#ifndef OOMPH_TRILINOS_SOLVER_HEADER
#define OOMPH_TRILINOS_SOLVER_HEADER

#include "trilinos_preconditioners.h"
#include "iterative_linear_solver.h"

#include "AztecOO.h"

namespace oomph 
{

//=============================================================================
/// \short An Epetra_Operator class for oomph-lib preconditioners.
/// A helper class for TrilinosOomphLibPreconditioner to allow an oomph-lib
/// preconditioner (i.e. one derived from Preconditioner) to be used with 
/// a trilinos solver (TrilinosAztecOOSolver)
//=============================================================================
 class OomphLibPreconditionerEpetraOperator : public Epetra_Operator 
 {
  public:

#ifdef OOMPH_HAS_MPI
  /// \short Constructor - takes the pointer to the oomph-lib 
  /// preconditioner and the distribution of the preconditioner\n
  /// \b Note: the oomph-lib preconditioner must be setup
  OomphLibPreconditionerEpetraOperator(Preconditioner* preconditioner_pt)
   : Operator_comm(preconditioner_pt->distribution().communicator())
   {
    // set the ooomph-lib preconditioner
    Oomph_lib_preconditioner_pt = preconditioner_pt;
    
    // set the preconditioner label
    //Preconditioner_label = new char[26];
    Preconditioner_label = "oomph-lib Preconditioner";

    
    // setup the Epetra_map

    // number of local rows
    unsigned nrow_local = preconditioner_pt->distribution().nrow_local();

    // first row
    unsigned first_row = preconditioner_pt->distribution().first_row();
    
    // create the map
    My_global_rows = new int[nrow_local];
    for (unsigned i = 0; i < nrow_local; i++)
     {
      My_global_rows[i] = first_row + i;
     }
    
    // the number of global rows
    unsigned nrow_global = preconditioner_pt->distribution().nrow_global();
    
    // create the operator
    Operator_map_pt = new Epetra_Map(nrow_global,nrow_local,
                                     My_global_rows,0,Operator_comm);
   }
#else
  /// \short Constructor - takes the pointer to the oomph-lib 
  /// preconditioner and the number of rows\n
  /// \b Note: the oomph-lib preconditioner must be setup 
  OomphLibPreconditionerEpetraOperator(Preconditioner* preconditioner_pt,
                                       const unsigned n_rows)
   {
    // set the ooomph-lib preconditioner
    Oomph_lib_preconditioner_pt = preconditioner_pt;
    
    // set the preconditioner label
    Preconditioner_label = "oomph-lib Preconditioner";

    // setup the Epetra_map
    Operator_map_pt = new Epetra_Map(n_rows,0,Operator_comm);

    // store the number of rows
    Nrow = n_rows;
  }
#endif

  /// \short Destructor - deletes the Epetra_map and My_global_rows vector
  /// (if MPI)
  ~OomphLibPreconditionerEpetraOperator()
   {
#ifdef OOMPH_HAS_MPI
    delete[] My_global_rows;
    My_global_rows = 0;
#endif
    delete Operator_map_pt;
    Operator_map_pt = 0;
   }

  /// Broken copy constructor
  OomphLibPreconditionerEpetraOperator
   (const OomphLibPreconditionerEpetraOperator&)
#ifdef OOMPH_HAS_MPI
   // (NOTE: MPI_COMM_WORLD used just to construct Operator_comm)
   : Operator_comm(MPI_COMM_WORLD) 
#else
   : Operator_comm()
#endif
  {
   BrokenCopy::broken_copy("OomphLibPreconditionerEpetraOperator");
  }

  /// Broken assignment operator.
  void operator=(const OomphLibPreconditionerEpetraOperator&)
  {
   BrokenCopy::broken_assign("OomphLibPreconditionerEpetraOperator");
  }

  /// Broken Epetra_Operator member - SetUseTranspose
  int SetUseTranspose(bool UseTranspose)
   {
    std::ostringstream error_message;
    error_message << "SetUseTranspose() is a pure virtual Epetra_Operator "
                  << "member that is not required for a Preconditioner" 
                  << std::endl;
    throw OomphLibError(
     error_message.str(),
     "OomphLibPreconditionerEpetraOperator::SetUseTranspose()",
     OOMPH_EXCEPTION_LOCATION);
   }


  /// Broken Epetra_Operator member - Apply
  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
   {
    std::ostringstream error_message;
    error_message << "Apply() is a pure virtual Epetra_Operator member"
                  << "that is not required for a Preconditioner" << std::endl;
    throw OomphLibError(
     error_message.str(),
     "OomphLibPreconditionerEpetraOperator::Apply()",
     OOMPH_EXCEPTION_LOCATION);
   }


  /// \short applies the oomph-lib preconditioner. Converts the Epetra vector
  /// applys the preconditioner by calling the oomph-lib preconditioner's
  /// preconditioner_solve functionality. 
  /// NOTE : the oomph-lib preconditioner is setup prior to being passed to 
  /// this class
  int ApplyInverse(const Epetra_MultiVector &epetra_r, 
                   Epetra_MultiVector &epetra_z) const
   {
#ifdef OOMPH_HAS_MPI
    // oomph-lib vector for r
    DistributedVector<double> 
     oomph_r(Oomph_lib_preconditioner_pt->distribution());
    
    // copy the Epetra_MultiVector r into an oomph-lib vector
    double** r_pt;
    epetra_r.ExtractView(&r_pt);
    unsigned nrow_local =
     Oomph_lib_preconditioner_pt->distribution().nrow_local();
    for (unsigned i = 0; i < nrow_local; i++)
     {
      oomph_r[i] = r_pt[0][i];
     }
    
    // oomph-lib vector for Y
    DistributedVector<double> 
     oomph_z(Oomph_lib_preconditioner_pt->distribution());
      
    // apply the preconditioner
    Oomph_lib_preconditioner_pt->preconditioner_solve(oomph_r,oomph_z);

    // copy the oomph-lib vector oomph_Y back to Y
    for (unsigned i = 0; i < nrow_local; i++)
     {
      epetra_z.ReplaceMyValue(i,0,oomph_z[i]);
     }
  
    // return 0 to indicate success
    return 0;
#else 
    // oomph-lib vector for r
    Vector<double> oomph_r(Nrow);
    
    // copy the Epetra_MultiVector r into an oomph-lib vector
    double** r_pt;
    epetra_r.ExtractView(&r_pt);
    for (unsigned i = 0; i < Nrow; i++)
     {
      oomph_r[i] = r_pt[0][i];
     }
    
    // oomph-lib vector for Y
    Vector<double> oomph_z(Nrow);
      
    // apply the preconditioner
    Oomph_lib_preconditioner_pt->preconditioner_solve(oomph_r,oomph_z);

    // copy the oomph-lib vector oomph_Y back to Y
    for (unsigned i = 0; i < Nrow; i++)
     {
      epetra_z.ReplaceMyValue(i,0,oomph_z[i]);
     }
  
    // return 0 to indicate success
    return 0;
#endif
   }


  /// Broken Epetra_Operator member - NormInf
  double NormInf() const
   {
    std::ostringstream error_message;
    error_message << "NormInf() is a pure virtual Epetra_Operator member"
                  << "that is not required for a Preconditioner" << std::endl;
    throw OomphLibError(
     error_message.str(),
     "OomphLibPreconditionerEpetraOperator::NormInf()",
     OOMPH_EXCEPTION_LOCATION);
   }

  /// Epetra_Operator::Label - returns a string describing the operator
  const char* Label() const
   {
    return Preconditioner_label.c_str();
   }

  /// Broken Epetra_Operator member - UseTranspose
  bool UseTranspose() const
   {
    std::ostringstream error_message;
    error_message << "UseTranspose() is a pure virtual Epetra_Operator member "
                  << "that is not required for a Preconditioner" << std::endl;
    throw OomphLibError(
     error_message.str(),
     "OomphLibPreconditionerEpetraOperator::UseTranspose()",
     OOMPH_EXCEPTION_LOCATION);
   }

  /// Broken Epetra_Operator member - HasNormInf
  bool HasNormInf() const
   {
    std::ostringstream error_message;
    error_message << "HasNormInf() is a pure virtual Epetra_Operator member "
                  << "that is not required for a Preconditioner" << std::endl;
    throw OomphLibError(
     error_message.str(),
     "OomphLibPreconditionerEpetraOperator::HasNormInf()",
     OOMPH_EXCEPTION_LOCATION);
   }

  /// Returns the Epetra MPI_Comm object
  const Epetra_Comm& Comm() const
   {
    return Operator_comm;
   }

  /// Epetra_Operator member - OperatorDomainMap
  const Epetra_Map& OperatorDomainMap() const
   {
    return *Operator_map_pt;
   }

  /// Epetra_Operator member - OperatorRangeMap
  const Epetra_Map& OperatorRangeMap() const
   {
    return *Operator_map_pt;
   }

   private:

  /// A pointer to the oomph-lib preconditioner
  Preconditioner* Oomph_lib_preconditioner_pt;

#ifdef OOMPH_HAS_MPI
  /// An Epetra MPI Comm object 
  Epetra_MpiComm Operator_comm;
#else
  /// An Epetra Serial Comm object 
  Epetra_SerialComm Operator_comm;
#endif

  /// \short A pointer to an Epetra_Map object - describes distribution of the
  /// preconditioner, in this instance it is primarily used to prescribe the 
  /// distribution
  /// of the residual and solution vector
  Epetra_Map* Operator_map_pt;

#ifdef OOMPH_HAS_MPI
  /// \short the global of the preconditioner associated with the processor,
  /// for the trilinos Operator_map_pt
  int* My_global_rows;
#else
  /// number of row
  unsigned Nrow;
#endif

  /// a label for the preconditioner ( for Epetra_Operator::Label() )
  //char* Preconditioner_label;
  std::string Preconditioner_label;
  
 };


//=============================================================================
/// \short An interface to the Trilinos AztecOO classes allowing it to be used 
/// as an Oomph-lib LinearSolver. \n
/// The AztecOO solver is a Krylov Subspace solver; the solver type (either CG,
/// GMRES or BiCGStab) can be set using solver_type(). \n
/// This solver can be preconditioned with Trilinos Preconditioners (derived
/// from TrilinosPreconditionerBase) or Oomph-lib preconditioners (derived
/// from Preconditioner). Preconditioners are set using preconditioner_pt().
//=============================================================================
class TrilinosAztecOOSolver : public IterativeLinearSolver
{
 public:

 /// Constructor.
 TrilinosAztecOOSolver()
 {
  // set pointer to Null
  AztecOO_solver_pt = 0;

  // initially assume not problem based solve
  Using_problem_based_solve = false;

#ifdef OOMPH_HAS_MPI
  // if a problem based solve is performed then it should generate a 
  // serial matrix
  Assemble_serial_jacobian = false;
#endif

  // null the pts
  Problem_pt = 0;
  Epetra_matrix_pt = 0;
  Epetra_map_pt = 0;
  Epetra_comm_pt = 0;
  Oomph_matrix_pt = 0;
  Epetra_preconditioner_pt = 0;
#ifdef OOMPH_HAS_MPI
  Epetra_global_rows = 0;
#endif

  // set solver defaults
  Solver_type = GMRES;
  Tolerance = 1e-10;
  Max_iter = 1000;
 }

 /// Destructor - delete the solver and the matrices
 ~TrilinosAztecOOSolver()
 {
  // delete 
  clean_up_memory();

  // if Problem_based solve then the oomph matrix was generated by this class
  // and should be deleted
  if (Using_problem_based_solve)
   {
    delete Oomph_matrix_pt;
    Oomph_matrix_pt = 0; 
   }
 }

 /// Broken copy constructor.
 TrilinosAztecOOSolver(const TrilinosAztecOOSolver&)
 {
  BrokenCopy::broken_copy("TrilinosAztecOOSolver");
 }

 /// Broken assignment operator.
 void operator=(const TrilinosAztecOOSolver&)
 {
  BrokenCopy::broken_assign("TrilinosAztecOOSolver");
 }

 /// Clean up method - deletes the solver, the matrices and the preconditioner
 void clean_up_memory()
  {
   // delete the solver
   delete AztecOO_solver_pt;
   AztecOO_solver_pt = 0;

   // delete the matrices
   delete Epetra_matrix_pt;
   Epetra_matrix_pt = 0;

   // delete the Epetra_Map
   delete Epetra_map_pt;
   Epetra_map_pt = 0;

   // delete the Epetra_Map
   delete Epetra_comm_pt;
   Epetra_comm_pt = 0;

#ifdef OOMPH_HAS_MPI
   // delete my global rows
   delete[] Epetra_global_rows;
   Epetra_global_rows = 0;
#endif

   // delete the Epetra_Operator preconditioner (only if it is a wrapper to an 
   // oomphlib preconditioner in which case only the wrapper is deleted and 
   // not the actual preconditioner).
   // if the preconditioner is Trilinos preconditioner then the
   // Epetra_Operator is deleted when that preconditioner is deleted
   if (dynamic_cast<OomphLibPreconditionerEpetraOperator*>
       (Epetra_preconditioner_pt) != 0)
    {
     delete Epetra_preconditioner_pt;
     Epetra_preconditioner_pt = 0;
    }
  }

 /// \short Function which uses problem_pt's get_jacobian(...) function to
 /// generate a linear system which is then solved. This function deletes
 /// any existing internal data and then generates a new AztecOO solver.
 void solve(Problem* const &problem_pt,Vector<double> &solution);

 /// \short Function to solve the linear system defined by matrix_pt and rhs.
 /// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
 /// DistributedCRDoubleMatrix.
 /// \b NOTE 2. This function will delete any existing internal data and 
 /// generate a new AztecOO solver.
 void solve(DoubleMatrixBase* const &matrix_pt,
            const Vector<double> &rhs,
            Vector<double> &solution);

#ifdef OOMPH_HAS_MPI
 /// \short Function to solve the linear system defined by matrix_pt and rhs\n.
 /// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
 /// DistributedCRDoubleMatrix.
 /// \b NOTE 2. This function will delete any existing internal data and 
 /// generate a new AztecOO solver.
 /// \b NOTE 3. The vector rhs must have the same distribution as the matrix
 void solve(DoubleMatrixBase* const &matrix_pt,
            const DistributedVector<double> &rhs,
            Vector<double> &solution);

 /// \short Function to solve the linear system defined by matrix_pt and rhs\n.
 /// \b NOTE 1. The matrix has to be of type CRDoubleMatrix or 
 /// DistributedCRDoubleMatrix.
 /// \b NOTE 2. This function will delete any existing internal data and 
 /// generate a new AztecOO solver. 
 /// \b NOTE 3. The vectors rhs and solution must have the same distribution 
 /// as the matrix
 void solve(DoubleMatrixBase* const &matrix_pt,
            const DistributedVector<double> &rhs,
            DistributedVector<double> &solution);
#endif

 /// \short Function to resolve a linear system using the existing solver
 /// data, allowing a solve with a new right hand side vector. This
 /// function must be used after a call to solve(...) with
 /// enable_resolve set to true.
 void resolve(const Vector<double> &rhs,
              Vector<double> &solution);

#ifdef OOMPH_HAS_MPI
 /// \short Function to resolve a linear system using the existing solver
 /// data, allowing a solve with a new right hand side vector. This
 /// function must be used after a call to solve(...) with
 /// enable_resolve set to true.
 void resolve(const DistributedVector<double> &rhs,
              Vector<double> &solution);

 /// \short Function to resolve a linear system using the existing solver
 /// data, allowing a solve with a new right hand side vector. This
 /// function must be used after a call to solve(...) with
 /// enable_resolve set to true.
 void resolve(const DistributedVector<double> &rhs,
              DistributedVector<double> &solution);
#endif

 /// \short Disable resolve function (overloads the LinearSolver
 /// disable_resolve function).
 void disable_resolve()
 {
  Enable_resolve=false;
  clean_up_memory();
 }

 /// Access function to Delete_matrix
 bool& delete_matrix() {return Delete_matrix;}

 /// Access function to Max_iter
 unsigned& max_iter() {return Max_iter;}

 /// Acess function to Iterations
 unsigned iterations() {return Iterations;}

 /// Access function to Tolerance
 double& tolerance() {return Tolerance;}

 /// Access function to Solver_type
 unsigned& solver_type() {return Solver_type;}

 /// Function to return Jacobian_setup_time;
 double jacobian_setup_time() {return Jacobian_setup_time;}

 /// Function to return Linear_solver_solution_time
 double linear_solver_solution_time() {return Linear_solver_solution_time;}

#ifdef OOMPH_HAS_MPI
 /// Access function to Assemble serial jacobian
 bool& assemble_serial_jacobian() { return Assemble_serial_jacobian; }

 /// Access function to Assemble serial jacobian (const version)
 bool assemble_serial_jacobian() const { return Assemble_serial_jacobian; }
#endif

 /// \short Enumerated list to define which AztecOO solver is used
 enum AztecOO_solver_types{CG,
                           GMRES,
                           BiCGStab};

 protected:

 /// \short Helper function performs the actual solve once the AztecOO
 /// solver is set up
 void solve_using_AztecOO(Epetra_Vector* &rhs_pt, Epetra_Vector* &soln_pt);

 /// \short Helper function for setting up the solver. Converts the oomph-lib 
 /// matrices to Epetra matrices, sets up the preconditioner, creates the 
 /// Trilinos Aztec00 solver and passes in the matrix, preconditioner and 
 /// parameters.
 void solver_setup(DoubleMatrixBase* const& matrix_pt);

 /// Maximum number of iterations used in solver.
 unsigned Max_iter;

 /// Tolerance used to terminate solver.
 double Tolerance;

 /// Stores number of iterations used
 unsigned Iterations;

 /// Pointer to the AztecOO solver
 AztecOO* AztecOO_solver_pt;

 /// Stores set up time for Jacobian
 double Jacobian_setup_time;

 /// Stores time for the solution (excludes time to set up preconditioner)
 double Linear_solver_solution_time;

 /// \short Trilinos copies matrix data from oomph-lib's own CRDoubleMatrix
 /// or DistributedCRDoubleMatrix to Trilinos's Epetra format - the Trilinos
 /// solver no longer requires the oomph-lib matrices and therefore they could
 /// be deleted to save memory. This must be requested explicitly by setting
 /// this flag to TRUE. \n
 /// \b NOTE: The matrix is deleted after the preconditioner is setup.
 bool Delete_matrix;

#ifdef OOMPH_HAS_MPI
 /// \short If true, when performing a problem based solve a serial matrix
 /// will be requested from Problem::get_jacobian(...). Defaults to true
 bool Assemble_serial_jacobian;
#endif

 /// \short Defines which solver is set up - available types are
 /// defined in AztecOO_solver_types
 unsigned Solver_type;

 /// \short Helper flag keeping track of whether we called the 
 /// linear algebra or problem-based solve function.
 bool Using_problem_based_solve;
 
 /// \short A pointer for the linear system matrix in Epetra_CrsMatrix format 
 Epetra_CrsMatrix* Epetra_matrix_pt;

 /// \short A pointer to the Epetra_Operator for the preconditioner. This is 
 /// only used if the preconditioner NOT a Trilinos preconditioner.
 Epetra_Operator* Epetra_preconditioner_pt;
 
 /// \short A pointer to the Epetra_Map object for this solver, which is the
 /// same as the matrix
 Epetra_Map* Epetra_map_pt;

#ifdef OOMPH_HAS_MPI
   /// \short Global rows of the Epetra_matrix_pt - only used when this 
   /// preconditioner is setup using the oomph-lib interface (and MPI)
   int* Epetra_global_rows;
#endif

 /// Oomph lib matrix pointer
 DoubleMatrixBase* Oomph_matrix_pt;

 /// \short A pointer to the underlying problem (NULL if MATRIX based solve)\n
 /// The problem_pt is stored here in a problem based solve for the
 /// preconditioner
 Problem* Problem_pt;

#ifdef OOMPH_HAS_MPI
 /// \short the Distribution of the solver
 DistributionInfo Solver_distribution;

   /// \short Epetra communicator object (MPI version)
   Epetra_MpiComm* Epetra_comm_pt;
#else
   /// \short Epetra communicator object (serial version)
   Epetra_SerialComm* Epetra_comm_pt;
#endif
};

}
#endif

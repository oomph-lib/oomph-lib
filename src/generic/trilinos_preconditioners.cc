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
#include "trilinos_preconditioners.h"

namespace oomph
{


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// functions for the TrilinosPreconditionerBase class
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=============================================================================
/// \short Function to set up a preconditioner for the linear system 
/// defined by matrix_pt. This function must be called before using 
/// preconditioner_solve. \n
/// \b NOTE 1. matrix_pt must point to an object of class CRDoubleMatrix or 
/// DistributedCRDoubleMatrix\n
/// This method should be called by oomph-lib solvers and preconditioners
//=============================================================================
void TrilinosPreconditionerBase::setup(Problem* problem_pt, 
                                       DoubleMatrixBase* matrix_pt)
{
 //clean up the memory
 clean_up_memory();

#ifdef PARANOID
 // check the matrix is square
 if ( matrix_pt->nrow() != matrix_pt->ncol() )
 {
  std::ostringstream error_message;
  error_message << "Preconditioners require a square matrix. "
                << "Matrix is " << matrix_pt->nrow()
                << " by " << matrix_pt->ncol() << std::endl;
  throw OomphLibError(error_message.str(),
                      "TrilinosPreconditionerBase::setup()",
                      OOMPH_EXCEPTION_LOCATION);
 }
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
   // the distribution of this preconditioner is the same as the matrix
   Preconditioner_distribution = cast_dist_matrix_pt->distribution();
   cast_failed = false;
  }

 // if the cast failed try cast to CRDoubleMatrix
 if (cast_failed)
 {
  CRDoubleMatrix* cast_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  if (cast_matrix_pt != 0)
  {
   // use default distribution
   // NOTE: USE of MPI_COMM_WORLD
   Preconditioner_distribution.distribute(MPI_COMM_WORLD,
                                          cast_matrix_pt->nrow());
   cast_failed = false;
  }
 }
 
#ifdef PARANOID
 // check the matrix is a 
 if (cast_failed)
  {
   std::ostringstream error_message;
   error_message << "TrilinosPreconditionerBase only work with "
                 << "DistributedCRDoubleMatrix or CRDoubleMatrix matrices" 
                 << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosPreconditionerBase::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // assemble the map
 // NOTE use of MPI_COMM_WORLD
 Epetra_comm_pt = new Epetra_MpiComm(MPI_COMM_WORLD);
 TrilinosHelpers::create_epetra_map(Preconditioner_distribution,
                                    Epetra_comm_pt,
                                    Epetra_map_pt,Epetra_global_rows);

#else
 // Serial Version
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
                       "TrilinosSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // setup the epetra map
 Epetra_comm_pt = new Epetra_SerialComm;
 TrilinosHelpers::create_epetra_map(cast_matrix_pt->nrow(),
                                    Epetra_comm_pt,
                                    Epetra_map_pt);
#endif
 
 // create the matrices
 TrilinosHelpers::create_epetra_matrix(matrix_pt,
                                       Epetra_map_pt,
                                       Epetra_matrix_pt);   
 
 // set up preconditioner
 setup_trilinos_preconditioner(problem_pt,matrix_pt,Epetra_matrix_pt);
}
 
//===========================================================================
/// \short Function to setup a preconditioner for the linear system defined
/// by the oomph-lib oomph_matrix_pt and Epetra epetra_matrix_pt matrices.\n
/// This method is called by Trilinos solvers.
//===========================================================================
void TrilinosPreconditionerBase::setup(Problem* problem_pt, 
                                       DoubleMatrixBase* oomph_matrix_pt, 
                                       Epetra_CrsMatrix* epetra_matrix_pt)
{
// clean up old data
 clean_up_memory();

// if MPI then store the distribution of the matrix - this will be the 
// distribution of the preconditioner
#ifdef OOMPH_HAS_MPI
 bool cast_failed = true;

 // first try DistributedCRDoubleMatrix
 DistributedCRDoubleMatrix* cast_dist_matrix_pt = 
  dynamic_cast<DistributedCRDoubleMatrix*>(oomph_matrix_pt);
 if (cast_dist_matrix_pt!=0)
  {
   // the distribution of this preconditioner is the same as the matrix
   Preconditioner_distribution = cast_dist_matrix_pt->distribution();
   cast_failed = false;
  }

 // if the cast failed try cast to CRDoubleMatrix
 if (cast_failed)
 {
  CRDoubleMatrix* cast_matrix_pt = 
   dynamic_cast<CRDoubleMatrix*>(oomph_matrix_pt);
  if (cast_matrix_pt != 0)
  {
   // if a distribution for the preconditioner is provided, use this, otherwise
   // use default distribution
   // NOTE: USE of MPI_COMM_WORLD
   Preconditioner_distribution.distribute(MPI_COMM_WORLD,
                                          cast_matrix_pt->nrow());
   cast_failed = false;
  }
 }
 
#ifdef PARANOID
 // check the matrix is a DistributedCRDoubleMatrix
 if (cast_failed)
  {
   std::ostringstream error_message;
   error_message << "TrilinosSolver only work with "
                 << "DistributedCRDoubleMatrix matrices" << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosSolver::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
#endif 
 
// setup the specific preconditioner
 setup_trilinos_preconditioner(problem_pt,oomph_matrix_pt,epetra_matrix_pt);
}

#ifdef OOMPH_HAS_MPI
//=============================================================================
/// \short preconditioner_solve - applies the preconditioner to the vector r 
/// taking distributed oomph-lib vectors (DistributedVector<double>) as 
/// arguments. \n
//=============================================================================
void TrilinosPreconditionerBase::
preconditioner_solve(const DistributedVector<double> &r, 
                     DistributedVector<double> &z)
{
#ifdef PARANOID
 // check preconditioner data exists
 if (Epetra_preconditioner_pt == 0)
 {
  std::ostringstream error_message;
  error_message << "preconditioner_solve requires that solver data has "
                << "been set up" << std::endl;
  throw OomphLibError(error_message.str(),
                      "TrilinosPreconditionerBase::preconditioner_solve()",
                      OOMPH_EXCEPTION_LOCATION);
 }

 if (Preconditioner_distribution != r.distribution())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
   throw OomphLibError(error_message.str(),
                       "TrilinosPreconditionerBase::preconditioner_solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif


 // set up trilinos vectors
 Epetra_Vector* epetra_r_pt=0;
 Epetra_Vector* epetra_z_pt=0;

 // create Epetra version of r
 TrilinosHelpers::create_epetra_vector(r,Epetra_map_pt,epetra_r_pt);

 // create an empty Epetra vector for z
 TrilinosHelpers::create_epetra_vector(Epetra_map_pt,epetra_z_pt);
 
 // Apply preconditioner                                          
 Epetra_preconditioner_pt->ApplyInverse(*epetra_r_pt,*epetra_z_pt);

 // Copy result to z
 z.distribute(Preconditioner_distribution);
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,
                                          Preconditioner_distribution,z);

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;
}
#endif

//=============================================================================
/// preconditioner_solve preconditions vector r taking serial oomph-lib 
/// vectors (Vector<double>) as arguments. \n
/// \b if MPI then we convert a serial oomph-lib vector (Vector<double>) to a
/// distributed Epetra vector, apply the preconditioner (in parallel) before 
/// copying the results back to a serial oomph-lib vector
//=============================================================================
void TrilinosPreconditionerBase::preconditioner_solve(const Vector<double> &r, 
                                                      Vector<double> &z)
{

#ifdef PARANOID
 // check preconditioner data exists
 if (Epetra_preconditioner_pt == 0)
 {
  std::ostringstream error_message;
  error_message << "preconditioner_solve requires that solver data has "
                << "been set up" << std::endl;
  throw OomphLibError(error_message.str(),
                      "TrilinosPreconditionerBase::preconditioner_solve()",
                      OOMPH_EXCEPTION_LOCATION);
 }

 // check vector length is compatable with matrix
 if (Epetra_preconditioner_pt->OperatorDomainMap().NumGlobalElements() 
     != int(r.size()) )
  {
   std::ostringstream error_message;
   error_message 
    << "TrilinosPreconditioner require input vector length equals number of"
    << " rows/columns in the matrix. "
    << "Vector length is " << r.size()
    << " and number of rows/columns in matrix is "
    << Epetra_preconditioner_pt->OperatorDomainMap().NumGlobalElements()
    << std::endl;
   throw OomphLibError(error_message.str(),
                       "TrilinosPreconditioner::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // set up trilinos vectors
 Epetra_Vector* epetra_r_pt=0;
 Epetra_Vector* epetra_z_pt=0;

 // create Epetra version of r
 TrilinosHelpers::create_epetra_vector(r,Epetra_map_pt,epetra_r_pt);

 // create an empty Epetra vector for z
 TrilinosHelpers::create_epetra_vector(Epetra_map_pt,epetra_z_pt);

 // Apply preconditioner                                          
 Epetra_preconditioner_pt->ApplyInverse(*epetra_r_pt, *epetra_z_pt);

#ifdef OOMPH_HAS_MPI
 // Copy result to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,z);
#else
 // Copy result to z
 TrilinosHelpers::copy_to_oomphlib_vector(epetra_z_pt,z);
#endif

 // clean up memory
 delete epetra_r_pt;
 delete epetra_z_pt;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// function for the TrilinosML class
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=============================================================================
// Function to set up the ML preconditioner.
//=============================================================================
void TrilinosMLPreconditioner::
setup_trilinos_preconditioner(Problem* problem_pt, 
                              DoubleMatrixBase* oomph_matrix_pt, 
                              Epetra_CrsMatrix* epetra_matrix_pt)
{
 // create the preconditioner
 Epetra_preconditioner_pt =
  new ML_Epetra::MultiLevelPreconditioner(*epetra_matrix_pt,
                                          ML_parameters,
                                          true);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// function for the TrilinosIFPACK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=============================================================================
// Function to set up the IFPACK preconditioner.
//=============================================================================
void TrilinosIFPACKPreconditioner::
setup_trilinos_preconditioner(Problem* problem_pt, 
                              DoubleMatrixBase* oomph_matrix_pt, 
                              Epetra_CrsMatrix* epetra_matrix_pt)
{ 
 // set up a parameter list
 Teuchos::ParameterList ifpack_parameters;
 ifpack_parameters.set("fact: level-of-fill", ILU_fill_level);
 ifpack_parameters.set("fact: ilut level-of-fill", ILUT_fill_level);

 // create the preconditioner
 Ifpack ifpack_factory;
 Ifpack_Preconditioner* tmp_pt= ifpack_factory.Create(Preconditioner_type,
                                                      epetra_matrix_pt,
                                                      Overlap);
 tmp_pt->SetParameters(ifpack_parameters);
 tmp_pt->Initialize();
 tmp_pt->Compute();
 
 // Now store the pointer to the newly created Ifpack_Preconditioner
 Epetra_preconditioner_pt=tmp_pt;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


}

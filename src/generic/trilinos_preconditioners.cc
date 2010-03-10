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


 // get a pointer to the cr double matrix
 CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
 if (cr_matrix_pt == 0)
  {
   throw OomphLibError("Matrix must be of type CRDoubleMatrix",
                       "TrilinosAztecOOSolver::factorise()",
                       OOMPH_EXCEPTION_LOCATION);   
  }
#endif

 // build the distribution of this preconditioner
 this->build_distribution(cr_matrix_pt->distribution_pt());
 
 // create the matrices
 Epetra_matrix_pt = TrilinosEpetraHelpers::create_distributed_epetra_matrix
  (cr_matrix_pt,this->distribution_pt());

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

 // first try CRDoubleMatrix
 CRDoubleMatrix* cr_matrix_pt = 
  dynamic_cast<CRDoubleMatrix*>(oomph_matrix_pt);
 if (cr_matrix_pt==0)
  {
     std::ostringstream error_message;
     error_message << "TrilinosSolver only work with "
                   << "DistributedCRDoubleMatrix matrices" << std::endl;
     throw OomphLibError(error_message.str(),
                         "TrilinosSolver::solve()",
                         OOMPH_EXCEPTION_LOCATION);
  }
 
 // setup the specific preconditioner
 setup_trilinos_preconditioner(problem_pt,oomph_matrix_pt,epetra_matrix_pt);
}


//=============================================================================
/// \short preconditioner_solve - applies the preconditioner to the vector r 
/// taking distributed oomph-lib vectors (DistributedVector<double>) as 
/// arguments. \n
//=============================================================================
void TrilinosPreconditionerBase::
preconditioner_solve(const DoubleVector &r, DoubleVector &z)
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

 if (*this->distribution_pt() != r.distribution_pt())
  {
   std::ostringstream error_message;
   error_message << "The rhs vector and the matrix must have the same "
                 << "distribution.\n"; 
   throw OomphLibError(error_message.str(),
                       "TrilinosPreconditionerBase::preconditioner_solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // create Epetra version of r
 Epetra_Vector* epetra_r_pt = 
  TrilinosEpetraHelpers::create_distributed_epetra_vector(r);

 // create an empty Epetra vector for z
 Epetra_Vector* epetra_z_pt = 
  TrilinosEpetraHelpers::create_distributed_epetra_vector
  (this->distribution_pt());
 
 // Apply preconditioner                                          
 Epetra_preconditioner_pt->ApplyInverse(*epetra_r_pt,*epetra_z_pt);

 // Copy result to z
 z.build(this->distribution_pt(),0.0);
 TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_z_pt,z);

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

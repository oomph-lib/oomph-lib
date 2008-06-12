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
#include "navier_stokes_preconditioners.h"


namespace oomph
{


//===========================================================================
/// Identify the required blocks: Here we only need
/// the momentum, gradient and divergence blocks of the
/// 2x2 block-structured matrix -- this function can be overloaded in
/// derived preconditioners.
//===========================================================================
 void NavierStokesLSCPreconditioner::identify_required_blocks(
  DenseMatrix<bool>& required_blocks)
{
 // Only need momentum, gradient and divergence blocks
 required_blocks(0,0) = true;
 required_blocks(1,0) = true;
 required_blocks(0,1) = true;
 required_blocks(1,1) = false;
}




//===========================================================================
/// Setup the least-squares commutator Navier Stokes preconditioner. This
/// extracts blocks corresponding to the velocity and pressure unknowns,
/// creates the matrices actually needed in the application of the
/// preconditioner and deletes what can be deleted... Note that
/// this preconditioner needs a CRDoubleMatrix.
//============================================================================
 void NavierStokesLSCPreconditioner::
 setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // NOTE: In the interest of minimising memory usage, several containers
  //       are recycled, therefore their content/meaning changes
  //       throughout this function. The code is carefully annotated
  //       but you'll have to read it line by line!
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // make sure any old data is deleted
  clean_up_memory();

  // Get blocks
  // ----------

  // In comes the current Jacobian. Recast it to a CR double matrix;
  // shout if that can't be done.
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message << "NavierStokesLSCPreconditioner only works with "
                  << "CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                     	"NavierStokesLSCPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // Set up block look up schemes (done automatically in the
  // BlockPreconditioner base class, based on the information 
  // provided in the block-preconditionable elements in the problem)
  this->block_setup(problem_pt);

  // find number of block types
  unsigned n_block = this->nblock_types();

  // Matrix to indicate which blocks are required 
  DenseMatrix<bool> required_blocks(n_block,n_block);

  // Identify required blocks
  identify_required_blocks(required_blocks);

  // Get the blocks and store pointers to them in Block_matrix_pt
  this->get_blocks(cr_matrix_pt, required_blocks, Block_matrix_pt);
  

  // Generate pressure Poisson matrix
  // --------------------------------
  P_matrix_pt = new CRDoubleMatrix;

  // Vector to hold inverse of velocity mass matrix diagonal
  Vector<double> ivmm_diagonal;

  // If required,  multiply Block_matrix_pt(1,0) which currently
  // points to the divergence matrix, D, by the inverse diagonal mass matrix:
  if (P_matrix_using_scaling)
   {
    // Get the diagonal of the velocity mass matrix from helper
    // function that assembles it from the Navier-Stokes elements
    // in the mesh.
    assemble_velocity_mass_matrix_diagonal(ivmm_diagonal);

    // Invert diagonal mass matrix
    unsigned n_entry=ivmm_diagonal.size();
    for (unsigned i=0; i<n_entry; i++)
     {
      ivmm_diagonal[i] = 1.0/ivmm_diagonal[i];
     }

    // Now multiply the divergence matrix D (from the right) by the inverse
    // diagonal mass matrix: Block_matrix_pt(1,0) now contains
    // D Q^{-1}.
    mat_diag_multiply(ivmm_diagonal, *Block_matrix_pt(1,0));
   }

  // Set Block_matrix_pt(1,0) matrix-matrix multiplication method
  Block_matrix_pt(1,0)->matrix_matrix_multiply_method() = Mult_method;

  // Now multiply whatever's currently stored in  Block_matrix_pt(1,0)
  // (i.e. D or D Q^{-1}) from the right by G, which is stored
  // in Block_matrix_pt(0,1). Store the result in P_matrix_pt
  // which now contains DG or D Q^{-1} G as required.
  Block_matrix_pt(1,0)->multiply(*Block_matrix_pt(0,1), *P_matrix_pt);


  // Generate matrix E for multiplication in Schur complement approximation
  // ----------------------------------------------------------------------

  // Auxiliary matrix for intermediate results
  CRDoubleMatrix* aux_matrix_pt = new CRDoubleMatrix;

  // Multiply the momentum matrix, stored in Block_matrix_pt(0,0)
  // by the content of Block_matrix_pt(1,0) (either D or
  // D Q^{-1}) and store the result, namely B Q^{-1} F in aux_matrix_pt.
  Block_matrix_pt(1,0)->multiply(*Block_matrix_pt(0,0), *aux_matrix_pt);


  // We no longer need the divergence matrix -- kill it.
  delete Block_matrix_pt(1,0);
  Block_matrix_pt(1,0) = 0;


  // Multiply auxiliary matrix by diagonal of velocity mass matrix if required
  if (P_matrix_using_scaling)
   {
    // The auxiliary matrix now contains B Q^{-1} F Q^{-1}
    mat_diag_multiply(ivmm_diagonal, *aux_matrix_pt);
   }

  // Set aux_matrix_pt matrix-matrix multiplication method
  aux_matrix_pt->matrix_matrix_multiply_method() = Mult_method;

  // Now multiply the auxiliary matrix by the gradient matrix,
  // stored in Block_matrix_pt(0,1), to obtain either
  // E = B Q^{-1} F Q^{-1} G or  E = B  F  G, as required.
  aux_matrix_pt->multiply(*Block_matrix_pt(0,1),E_matrix);

  // Clean up memory
  delete aux_matrix_pt;
  ivmm_diagonal.clear();


  // Set up solver for the system involving the pressure Poisson matrix
  // ------------------------------------------------------------------
#ifdef PARANOID
  // check a solver has been set
  if (P_preconditioner_pt==0)
  {
   std::ostringstream error_message;
   error_message << "P_preconditioner_pt has not been set.";

   throw OomphLibError(error_message.str(),
                    	  "NavierStokesLSCPreconditioner::setup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Setup the preconditioner for the Pressure matrix
  P_preconditioner_pt->setup(problem_pt, P_matrix_pt);

  // Delete oomph-lib pressure Poisson matrix
  delete P_matrix_pt;
  P_matrix_pt=0;


  // Set up solver for solution of system with momentum matrix
  // ----------------------------------------------------------
#ifdef PARANOID
  // check a solver has been set
  if (F_preconditioner_pt==0)
  {
   std::ostringstream error_message;
   error_message << "F_preconditioner_pt has not been set.";

   throw OomphLibError(error_message.str(),
                   	  "NavierStokesLSCPreconditioner::setup()",
                      OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Setup the preconditioner the F matrix
  F_preconditioner_pt->setup(problem_pt, Block_matrix_pt(0,0));

  // Now we don't need the original momentum matrix any more. Kill it.
  delete Block_matrix_pt(0,0);
  Block_matrix_pt(0,0) = 0;

  // Remember that the preconditioner has been setup so
  // the stored information can be wiped when we
  // come here next...
  Preconditioner_has_been_setup = true;
 }



//=======================================================================
 /// Apply preconditioner to r.
//=======================================================================
 void NavierStokesLSCPreconditioner::
 preconditioner_solve(const Vector<double>&r, Vector<double> &z)
 {

#ifdef PARANOID
  if (Preconditioner_has_been_setup==false)
   {
    std::ostringstream error_message;
    error_message << "setup must be called before using preconditioner_solve";
    throw OomphLibError(
     error_message.str(),
     "NavierStokesLSCPreconditioner::preconditioner_solve()",
     OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Number of velocity and pressure degrees of freedom
  unsigned n_veloc_dof = this->block_dimension(0);
  unsigned n_press_dof = this->block_dimension(1);

  // Cache the number of rows into a local variable
  //const unsigned n_row = this->Nrow;


  // Step 1 - apply approximate Schur inverse to pressure unknowns (block 1)
  // -----------------------------------------------------------------------

  // Working vectors, big enough to hold all pressure dofs.
  Vector<double> temp_vec(n_press_dof);
  Vector<double> another_temp_vec(n_press_dof);

  // Copy pressure values from residual vector to temp_vec:
  // Loop over all entries in the global vector (this one
  // includes velocity and pressure dofs in some random fashion)
  this->get_block_vector(1,r,temp_vec);

  // NOTE: The vector temp_vec now contains the vector r_p.


  // Solve first pressure Poisson system
#ifdef PARANOID
  // check a solver has been set
  if (P_preconditioner_pt==0)
   {
    std::ostringstream error_message;
    error_message << "P_preconditioner_pt has not been set.";

    throw OomphLibError(error_message.str(),
                     	  "NavierStokesLSCPreconditioner::preconditioner_solve()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // use some Preconditioner's preconditioner_solve function
  P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);


  // NOTE: The vector another_temp_vec now contains the vector P^{-1} r_p

  // Multiply another_temp_vec by matrix E and stick the result into temp_vec
  temp_vec.clear();
  E_matrix.multiply(another_temp_vec, temp_vec);

  // NOTE: The vector temp_vec now contains E P^{-1} r_p

  // Solve second pressure Poisson system using preconditioner_solve
  P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);


  // NOTE: The vector another_temp_vec now contains  z_p = P^{-1} E P^{-1} r_p
  //       as required (apart from the sign which we'll fix in the
  //       next step.

  // Now copy another_temp_vec (i.e. z_p) back into the global vector z.
  // Loop over all entries in the global results vector z:
  temp_vec.resize(n_press_dof);
  for (unsigned i=0; i< n_press_dof; i++)
   {
    // copy the result across
    temp_vec[i] = -another_temp_vec[i];
   }
  return_block_vector(1,temp_vec,z);
  
  // Step 2 - apply preconditioner to velocity unknowns (block 0)
  // ------------------------------------------------------------
  
  // Recall that another_temp_vec (computed above) contains the
  // negative of the solution of the Schur complement systen, -z_p.
  // Multiply by G (stored in Block_matrix_pt(0,1) and store
  // result in temp_vec (vector resizes itself).
  temp_vec.clear();
  Block_matrix_pt(0,1)->multiply(another_temp_vec, temp_vec);

 // NOTE: temp_vec now contains -G z_p

 // The vector another_temp_vec is no longer needed -- re-use it to store
 // velocity quantities:
 another_temp_vec.resize(n_veloc_dof);

 // Loop over all enries in the global vector and find the
 // entries associated with the velocities:
 get_block_vector(0,r,another_temp_vec);
 for (unsigned i = 0; i < n_veloc_dof; i++)
  {
   another_temp_vec[i] += temp_vec[i];
  }
 
 // NOTE:  The vector another_temp_vec now contains r_u - G z_p

 // Solve momentum system
#ifdef PARANOID
 // check a solver has been set
 if (F_preconditioner_pt==0)
  {
   std::ostringstream error_message;
   error_message << "F_preconditioner_pt has not been set.";

   throw OomphLibError(error_message.str(),
                    	  "NavierStokesLSCPreconditioner::preconditioner_solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // use some Preconditioner's preconditioner solve
 F_preconditioner_pt->preconditioner_solve(another_temp_vec, temp_vec);


 // NOTE: The vector temp_vec now contains F^{-1}  ( r_u - G z_p ) = z_u,
 //       i.e. the required result.

 // Now copy back into the global vector
 return_block_vector(0,temp_vec,z);
}



//===================================================================
/// Helper function to multiply a CRDoubleMatrix by a
/// diagonal matrix held in diag_matrix. The input and output
/// matrices are held in matrix.
//===================================================================
 void NavierStokesLSCPreconditioner::
 mat_diag_multiply(const Vector<double>& diag_matrix,
                   CRDoubleMatrix& matrix)
 {

#ifdef PARANOID
  // check for compatible dimensions
  if (matrix.ncol() != diag_matrix.size() )
   {
    std::ostringstream error_message;
    error_message << "Incompatible matrix dimensions. "
                  << "CRDoubleMatrix has dimensions " << matrix.nrow()
                  << "x" << matrix.ncol()
                  << ". Diagonal matrix is " << diag_matrix.size()
                  << "x" << diag_matrix.size();

    throw OomphLibError(error_message.str(),
                     	"NavierStokesLSCPreconditioner::mat_diag_multiply()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif


  // get pointers to values etc in matrix
  const int* row_start = matrix.row_start();
  const int* column_index = matrix.column_index();
  double* value = matrix.value();

  // run through rows of matrix
  unsigned n_row=matrix.nrow();
  for (unsigned row=0; row<n_row; row++)
   {
    // run through entries in the row
    unsigned next_row_start=row_start[row+1];
    for (unsigned ptr=row_start[row];
    				 ptr<next_row_start;
         ptr++)
     {
      // find column index
      int col = column_index[ptr];

      // multiply value
      value[ptr] *= diag_matrix[col];
     }
   }
 }


//========================================================================
/// Helper function to assemble the diagonal of the velocity
/// mass matrix from the elemental contributions defined in
/// NavierStokesEquations<DIM>::get_velocity_mass_matrix_diagonal(...).
//========================================================================
 void NavierStokesLSCPreconditioner::
 assemble_velocity_mass_matrix_diagonal(Vector<double> &vmm_diagonal)
 {
  // find number of velocity degrees of freedom
  int n_dof = this->block_dimension(0);

  // find number of elements
  unsigned n_el = Mesh_pt[0]->nelement();

  // find the number of dimensions from the first element
  unsigned dim = dynamic_cast<FiniteElement*>(Mesh_pt[0]->element_pt(0) )->dim();

  // pointers to the only two possible element types
  // (required for accessing the function get_velocity_mass_matrix_diagonal)
  NavierStokesEquations<2>* el2d_pt=0;
  NavierStokesEquations<3>* el3d_pt=0;

  // Resize and initialise
  vmm_diagonal.assign(n_dof, 0.0);

  // loop over the elements
  for(unsigned e=0; e<n_el; e++)
   {
    // find number of degrees of freedom in the element
    // (this is slightly too big because it includes the
    // pressure dofs but this doesn't matter)
    unsigned el_dof = Mesh_pt[0]->element_pt(e)->ndof();

    // allocate local storage for the element's contribution to the
    // velocity mass matrix diagonal
    Vector<double> el_vmm_diagonal(el_dof);

    // get the element contribution
    if (dim==2)
    {
     el2d_pt = dynamic_cast< NavierStokesEquations<2>* >( Mesh_pt[0]->element_pt(e) );
     el2d_pt->get_velocity_mass_matrix_diagonal(el_vmm_diagonal);
    }
    else
    {
     el3d_pt = dynamic_cast< NavierStokesEquations<3>* >( Mesh_pt[0]->element_pt(e) );
     el3d_pt->get_velocity_mass_matrix_diagonal(el_vmm_diagonal);
    }

    // insert values into vmm_diagonal
    // - loop over local variables
    for(unsigned i=0; i<el_dof; i++)
     {
      //Get the equation number
      unsigned eqn_number = Mesh_pt[0]->element_pt(e)->eqn_number(i);

      // bypass pressure dof
      if ( this->block_number(eqn_number)==0 )
       {
        // add the contribution to the correct position based on block ordering
        vmm_diagonal[ this->index_in_block(eqn_number) ] += el_vmm_diagonal[i];
       }
     }
   }

 }




//=========================================================================
/// Helper function to delete preconditioner data.
//=========================================================================
 void NavierStokesLSCPreconditioner::clean_up_memory()
 {

  if (Preconditioner_has_been_setup)
   {
    // delete blocks
    unsigned nr = Block_matrix_pt.nrow();
    unsigned nc = Block_matrix_pt.ncol();
    for (unsigned i = 0 ; i < nr; i++)
     for (unsigned j = 0 ; j < nc; j++)
      {
       delete Block_matrix_pt(i,j);
       Block_matrix_pt(i,j) = 0;
      }

    // delete stuff from velocity solve
    if (F_preconditioner_pt!=0)
     {
      F_preconditioner_pt->clean_up_memory();
     }

    // delete stuff from Schur complement approx
    delete P_matrix_pt;
    P_matrix_pt = 0;
    if (P_preconditioner_pt!=0)
     {
      P_preconditioner_pt->clean_up_memory();
     }
   }
   
   // clean up memory in default preconditioners
   F_superlu_preconditioner.clean_up_memory();
   P_superlu_preconditioner.clean_up_memory();

 }
}

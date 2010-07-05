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

#include "pseudo_elastic_preconditioner.h"

namespace oomph
{

 //=============================================================================
 /// \short Functions to create instances of optimal subsidiary operators for
 /// the PseudoElasticPreconditioner
 //=============================================================================
 namespace Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper
 {
#ifdef HAVE_HYPRE
  /// \short AMG w/ GS smoothing for the augmented elastic subsidiary linear
  /// systems
  Preconditioner* get_elastic_preconditioner()
  {
   HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;
   hypre_preconditioner_pt->set_amg_iterations(2);
   hypre_preconditioner_pt->amg_using_simple_smoothing();
   hypre_preconditioner_pt->amg_simple_smoother() = 3;
   hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
   hypre_preconditioner_pt->amg_strength() = 0.25;
   hypre_preconditioner_pt->amg_damping() = 1.0;
   hypre_preconditioner_pt->amg_coarsening() = 6;
   return hypre_preconditioner_pt;
  }
#endif
  
#ifdef HAVE_TRILINOS
  /// \short CG with diagonal preconditioner for the lagrange multiplier
  /// subsidiary linear systems.
  Preconditioner* get_lagrange_multiplier_preconditioner()
  {
   InnerIterationPreconditioner
    <TrilinosAztecOOSolver,MatrixBasedDiagPreconditioner>* prec_pt = 
    new InnerIterationPreconditioner
    <TrilinosAztecOOSolver,MatrixBasedDiagPreconditioner>;
   // Note: This makes CG a proper "inner iteration" for
   // which GMRES (may) no longer converge. We should really
   // use FGMRES or GMRESR for this. However, here the solver
   // is so good that it'll converge very quickly anyway
   // so there isn't much to be gained by limiting the number
   // of iterations...
   // prec_pt->max_iter() = 4;
   prec_pt->solver_pt()->solver_type() = TrilinosAztecOOSolver::CG;
   prec_pt->solver_pt()->doc_time()=false;
   return prec_pt;
  }
#endif
 }


 
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////



 //=============================================================================
 // Setup method for the PseudoElasticPreconditioner.
 //=============================================================================
 void PseudoElasticPreconditioner::setup(Problem* problem_pt, 
                                         DoubleMatrixBase* matrix_pt)
 {
  // clean
  this->clean_up_memory();
  
#ifdef PARANOID
  // paranoid check that meshes have been set
  if (Elastic_mesh_pt==0)
   {
    std::ostringstream error_message;
    error_message << "The elastic mesh must be set.";
    throw OomphLibError(error_message.str(),
                        "PseudoElasticPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
    
   }
  if (Lagrange_multiplier_mesh_pt==0)
   {
    std::ostringstream error_message;
    error_message << "The Lagrange multiplier mesh must be set.";
    throw OomphLibError(error_message.str(),
                        "PseudoElasticPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // set the mesh
  unsigned n_solid_dof_types = 0;
  unsigned n_dof_types = 0;
  this->set_mesh(0,problem_pt,Elastic_mesh_pt);
  this->set_mesh(1,problem_pt,Lagrange_multiplier_mesh_pt);
  if (this->is_master_block_preconditioner())
   {
    
    // get the number of solid dof types from the first element
    n_solid_dof_types = this->ndof_types_in_mesh(0);
    
    // get the total number of dof types
    n_dof_types = n_solid_dof_types 
     + this->ndof_types_in_mesh(1);
   }
  else
   {
    n_dof_types = this->ndof_types();
    n_solid_dof_types = (int)(((double)2*n_dof_types)/3);
   }
#ifdef PARANOID
  if (n_dof_types%3 != 0)
   {
    std::ostringstream error_message;
    error_message << "This preconditioner requires DIM*3 types of DOF";
    throw OomphLibError(error_message.str(),
                        "PseudoElasticPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // determine the dimension
  Dim = n_dof_types/3;
  
  // Call block setup for this preconditioner
  this->block_setup(problem_pt,matrix_pt);

  // Recast Jacobian matrix to CRDoubleMatrix
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
  
#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message << "FSIPreconditioner only works with"
                  << " CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                        "PseudoElasticPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // compute the scaling (if required)
  if (Use_inf_norm_of_s_scaling)
   {
    Vector<unsigned> dof_list(n_solid_dof_types);
    for (unsigned i = 0; i < n_solid_dof_types; i++)
     {
      dof_list[i] = i;
     }
    PseudoElasticPreconditionerScalingHelper* helper_pt = 
     new PseudoElasticPreconditionerScalingHelper(problem_pt,
                                                  this,
                                                  cr_matrix_pt,
                                                  dof_list);
    Scaling = helper_pt->s_inf_norm();
    delete helper_pt;
   }
  else
   {
    Scaling = 1.0;
   }
  
  // setup the solid subsidiary preconditioner
  
  // this preconditioner uses the full S matrix
  if (E_preconditioner_type == Exact_block_preconditioner)
   {
    PseudoElasticPreconditionerSubsidiaryPreconditioner* s_prec_pt = 
     new PseudoElasticPreconditionerSubsidiaryPreconditioner;
    Vector<unsigned> dof_list(n_solid_dof_types);
    for (unsigned i = 0; i < n_solid_dof_types; i++)
     {
      dof_list[i] = i;
     }
    s_prec_pt->turn_into_subsidiary_block_preconditioner(this,dof_list);
    if (Elastic_subsidiary_preconditioner_function_pt != 0)
     {
      s_prec_pt->
       set_subsidiary_preconditioner_function
       (Elastic_subsidiary_preconditioner_function_pt);
     }
    s_prec_pt->scaling() = Scaling;
    s_prec_pt->setup(problem_pt,matrix_pt);
    Elastic_preconditioner_pt = s_prec_pt;
   }
  
  // otherwise it is a block based preconditioner
  else 
   {
    // create the preconditioner
    PseudoElasticPreconditionerSubsidiaryBlockPreconditioner* 
     s_prec_pt = 
     new PseudoElasticPreconditionerSubsidiaryBlockPreconditioner;
    Vector<unsigned> dof_list(n_solid_dof_types);
    for (unsigned i = 0; i < n_solid_dof_types; i++)
     {
      dof_list[i] = i;
     }
    s_prec_pt->turn_into_subsidiary_block_preconditioner(this,dof_list);
    
    // set the subsidiary solve method
    if (Elastic_subsidiary_preconditioner_function_pt != 0)
     {
      s_prec_pt->
       set_subsidiary_preconditioner_function
       (Elastic_subsidiary_preconditioner_function_pt);
     }
    
    // set the scaling
    s_prec_pt->scaling() = Scaling;
    
    // set the block preconditioning method
    switch (E_preconditioner_type)
     {
     case Block_diagonal_preconditioner:
      s_prec_pt->block_diagonal();
      break;
     case Block_upper_triangular_preconditioner:
      s_prec_pt->upper_triangular();
      break;
     case Block_lower_triangular_preconditioner:
      s_prec_pt->lower_triangular();
      break;
     default:
      break;
     }
    
    // setup
    s_prec_pt->setup(problem_pt,matrix_pt);
    Elastic_preconditioner_pt = s_prec_pt;
   }
  
  // next setup the lagrange multiplier preconditioners
  Lagrange_multiplier_preconditioner_pt.resize(Dim);
  for (unsigned d = 0; d < Dim; d++)
   {
    CRDoubleMatrix* b_pt = 0;
    this->get_block(2*Dim+d,Dim+d,cr_matrix_pt,b_pt);
    
    // if a non default preconditioner is specified create 
    // the preconditioners
    if (Lagrange_multiplier_subsidiary_preconditioner_function_pt != 0)
     {
      Lagrange_multiplier_preconditioner_pt[d] = 
       (*Lagrange_multiplier_subsidiary_preconditioner_function_pt)();
     }
    
    // else use default superlu preconditioner
    else
     {
      Lagrange_multiplier_preconditioner_pt[d] = new SuperLUPreconditioner;
     }
    
    // and setup
    Lagrange_multiplier_preconditioner_pt[d]->setup(problem_pt,b_pt);
    delete b_pt;
   }
 }
 
 //=============================================================================
 /// \short Apply the elastic subsidiary preconditioner.
 //=============================================================================
 void PseudoElasticPreconditioner::elastic_preconditioner_solve
 (const DoubleVector& r, DoubleVector& z)
 {
  // apply the solid preconditioner
  Elastic_preconditioner_pt->preconditioner_solve(r,z);
 }
 
 //=============================================================================
 /// \short Apply the lagrange multiplier subsidiary preconditioner.
 //=============================================================================
 void PseudoElasticPreconditioner::
 lagrange_multiplier_preconditioner_solve(const DoubleVector& r,
                                          DoubleVector& z)
 {
  // apply the lagrange multiplier preconditioner
  for (unsigned d = 0; d < Dim; d++)
   {
    DoubleVector x;
    this->get_block_vector(Dim*2+d,r,x);
    DoubleVector y;
    Lagrange_multiplier_preconditioner_pt[d]->preconditioner_solve(x,y);
    Lagrange_multiplier_preconditioner_pt[d]->preconditioner_solve(y,x);
    unsigned nrow_local = x.nrow_local();
    double* x_pt = x.values_pt();
    for (unsigned i = 0; i < nrow_local; i++)
     {
      x_pt[i] = x_pt[i] * Scaling;
     }     
    this->return_block_vector(Dim*2+d,x,z);
   }
 }

 //=============================================================================
 /// \short Clears the memory.
 //=============================================================================
 void PseudoElasticPreconditioner::clean_up_memory()
 {
  // clean the block preconditioner base class memory
  this->clear_block_preconditioner_base();
  
  // delete the solid preconditioner
  delete Elastic_preconditioner_pt;
  Elastic_preconditioner_pt = 0;
  
  // delete the lagrange multiplier preconditioner pt
  unsigned sz = Lagrange_multiplier_preconditioner_pt.size();
  for (unsigned i = 0; i < sz; i++)
   {
    delete Lagrange_multiplier_preconditioner_pt[i];
    Lagrange_multiplier_preconditioner_pt[i] = 0;
   }
 }



 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////



 //=============================================================================
 /// \short Setup the preconditioner
 //=============================================================================
 void PseudoElasticPreconditionerSubsidiaryPreconditioner::
 setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {
   // clean memory
   this->clean_up_memory();
     
#ifdef PARANOID
   // paranoid check that this preconditioner has an even number of DOF types
   if (this->ndof_types()%2 != 0)
    {
     std::ostringstream error_message;
     error_message
      << "This SUBSIDIARY preconditioner requires an even number of "
      << "types of DOF";
     throw OomphLibError(
      error_message.str(),
      "PrecribedBoundaryDisplacementSubisiaryPreconditioner::setup(...)",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
     
   // assemble dof_to_block_map
   unsigned ndof_types = this->ndof_types();
   Vector<unsigned> dof_to_block_map(ndof_types,0);
   for (unsigned i = ndof_types/2; i < ndof_types; i++)
    {
     dof_to_block_map[i] = 1;
    }
     
   // call block_setup(...)
   this->block_setup(problem_pt,matrix_pt,dof_to_block_map);

   // Recast Jacobian matrix to CRDoubleMatrix
   CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);
     
#ifdef PARANOID
   // paraoind check we could cast to CRDoubleMatrix
   if (cr_matrix_pt==0)
    {
     std::ostringstream error_message;
     error_message << "FSIPreconditioner only works with"
                   << " CRDoubleMatrix matrices" << std::endl;
     throw OomphLibError(error_message.str(),
                         "FSIPreconditioner::setup()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
     
   // get block 11
   CRDoubleMatrix* s11_pt = 0;
   this->get_block(1,1,cr_matrix_pt,s11_pt);

   // add the scaled identity matrix to block 11
   double* s11_values = s11_pt->value();
   int* s11_column_index = s11_pt->column_index();
   int* s11_row_start = s11_pt->row_start();
   int s11_nrow_local = s11_pt->nrow_local();
   int s11_first_row = s11_pt->first_row();
   for (int i = 0; i < s11_nrow_local; i++)
    {
     bool found = false;
     for (int j = s11_row_start[i]; 
          j < s11_row_start[i+1] && !found; j++)
      {
       if (s11_column_index[j] == i + s11_first_row)
        {
         s11_values[j] += Scaling;
         found = true;
        }
      }
    }

   // get the remaining block and build the preconditioner
   DenseMatrix<CRDoubleMatrix* > s_pt(2,2,0);
   this->get_block(0,0,cr_matrix_pt,s_pt(0,0));
   this->get_block(0,1,cr_matrix_pt,s_pt(0,1));
   this->get_block(1,0,cr_matrix_pt,s_pt(1,0));
   s_pt(1,1) = s11_pt;
   CRDoubleMatrix* s_prec_pt = 0;
   this->build_preconditioner_matrix(s_pt,s_prec_pt);
   delete s_pt(0,0);
   delete s_pt(0,1);
   delete s_pt(1,0);
   delete s_pt(1,1);
   
   // setup the preconditioner
   if (Subsidiary_preconditioner_function_pt != 0)
    {
     Preconditioner_pt = (*Subsidiary_preconditioner_function_pt)();
    }
   else
    {
     Preconditioner_pt = new SuperLUPreconditioner;
    }
   Preconditioner_pt->setup(problem_pt,s_prec_pt);
   delete s_prec_pt;
  }
   
 //=============================================================================
 /// \short Apply the preconditioner.
 //=============================================================================
 void PseudoElasticPreconditionerSubsidiaryPreconditioner::
 preconditioner_solve(const DoubleVector& r, DoubleVector& z)
 {
   DoubleVector x;
   this->get_block_ordered_preconditioner_vector(r,x);
   DoubleVector y;
   Preconditioner_pt->preconditioner_solve(x,y);
   this->return_block_ordered_preconditioner_vector(y,z);
  }



 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////


 
 //=============================================================================
 /// clean up the memory
 //=============================================================================
 void PseudoElasticPreconditionerSubsidiaryBlockPreconditioner::
 clean_up_memory()
 {
   //number of block types
  unsigned n_block = Diagonal_block_preconditioner_pt.size();
   
   //delete diagonal blocks
   for (unsigned i = 0 ; i < n_block; i++)
    {
     delete Diagonal_block_preconditioner_pt[i];
     Diagonal_block_preconditioner_pt[i] = 0;
     if (Method == 1)
      {
       for (unsigned j = i+1; j < n_block; j++)
        {
         delete Off_diagonal_matrix_vector_products(i,j);
         Off_diagonal_matrix_vector_products(i,j) = 0;
        }
      }
     else if (Method == 2)
      {
       for (unsigned j = 0; j < i; j++)
        {
         delete Off_diagonal_matrix_vector_products(i,j);
         Off_diagonal_matrix_vector_products(i,j) = 0;
        }
      }
    }
   
   // clean up the block preconditioner
   this->clear_block_preconditioner_base();
 }
 
 //=============================================================================
 /// \short Setup the preconditioner.
 //=============================================================================
 void PseudoElasticPreconditionerSubsidiaryBlockPreconditioner::
 setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {
  // clean the memory
  this->clean_up_memory();
  
  // determine the number of DOF types
  unsigned n_dof_types =  this->ndof_types();
  
#ifdef PARANOID
  // must be Dim*2 dof types
  if (n_dof_types%2 != 0)
   {
    std::ostringstream error_message;
    error_message << "This preconditioner requires DIM*3 types of DOF";
    throw OomphLibError(error_message.str(),
                        "PseudoElasticPreconditioner",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // store the dimension of the problem
  unsigned dim = n_dof_types/2;
  
  // assemble the dof to block lookup scheme
  Vector<unsigned> dof_to_block_map(n_dof_types,0);
  for (unsigned d = 0; d < dim; d++)
   {
    dof_to_block_map[d] = d;
    dof_to_block_map[d+dim] = d;
   }
  
  //setup the blocks look up schemes
  this->block_setup(problem_pt,matrix_pt,dof_to_block_map);

  // Need to recast here -- input type is determined by specs in
  // base class.
  CRDoubleMatrix* cast_matrix_pt=dynamic_cast<CRDoubleMatrix*>(matrix_pt);
#ifdef PARANOID
  if (cast_matrix_pt==0)
   {
    throw OomphLibError("Wasn't able to cast matrix to CRDoubleMatrix.",
                        "BlockTriangularPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Storage for the diagonal block preconditioners
  Diagonal_block_preconditioner_pt.resize(dim);
  
  // storage for the off diagonal matrix vector products
  Off_diagonal_matrix_vector_products.resize(dim,dim,0);
  
  // setup the subsidiary preconditioners
  for (unsigned d = 0; d < dim; d++)
   {
    Vector<unsigned> dof_list(2);
    dof_list[0]=d;
    dof_list[1]=d+dim;
    Diagonal_block_preconditioner_pt[d] = new 
     PseudoElasticPreconditionerSubsidiaryPreconditioner; 
    Diagonal_block_preconditioner_pt[d]->
     turn_into_subsidiary_block_preconditioner(this,dof_list);
    if (Subsidiary_preconditioner_function_pt != 0)
     {
      Diagonal_block_preconditioner_pt[d]->
       set_subsidiary_preconditioner_function
       (Subsidiary_preconditioner_function_pt);
     }
    Diagonal_block_preconditioner_pt[d]->scaling() = Scaling;
    Diagonal_block_preconditioner_pt[d]->setup(problem_pt,matrix_pt);
    
    // next setup the off diagonal mat vec operators if required
    if (Method == 1 || Method == 2)
     {
      unsigned l = d+1;
      unsigned u = dim;
      if (Method==2)
       {
        l = 0;
        u = d;
       }
      for (unsigned j = l; j < u; j++)
       {
        CRDoubleMatrix* block_matrix_pt = 0;
        this->get_block(d,j,cast_matrix_pt,block_matrix_pt);
        Off_diagonal_matrix_vector_products(d,j) 
         = new MatrixVectorProduct();
        Off_diagonal_matrix_vector_products(d,j)->setup(block_matrix_pt);
        delete block_matrix_pt;
       }
     }
   }
 }
 
 //=============================================================================
 /// Apply preconditioner to r
 //=============================================================================
 void PseudoElasticPreconditionerSubsidiaryBlockPreconditioner::
 preconditioner_solve(const DoubleVector &res, DoubleVector &z)
  {
   // copy r
   DoubleVector r(res);
   
   // Cache umber of block types (also the spatial DIM)
   const unsigned n_block = this->nblock_types();
     
   // loop parameters
   int start = n_block-1;
   int end = -1;
   int step = -1;
   if (Method != 1)
    {
     start = 0;
     end = n_block;
     step = 1;
    }

   // loop over the DIM
   for (int i = start; i != end; i+=step)
    {
       
     // solve
     Diagonal_block_preconditioner_pt[i]->preconditioner_solve(r,z);
       
     // if upper or lower triangular
     if (Method != 0)
      {

       // substitute
       for (int j = i + step; j !=end; j+=step)
        {
         DoubleVector x;
         this->get_block_vector(i,z,x);
         DoubleVector y;
         Off_diagonal_matrix_vector_products(j,i)->multiply(x,y);
         x.clear();
         this->get_block_vector(j,r,x);
         x -= y;
         this->return_block_vector(j,x,r);
        }
      }
    }
  }
}

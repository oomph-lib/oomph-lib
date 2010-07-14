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
#ifndef OOMPH_PSEUDO_ELASTIC_SUBSIDIARY_PRECONDITIONER
#define OOMPH_PSEUDO_ELASTIC_SUBSIDIARY_PRECONDITIONER

// includes
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"
#ifdef OOMPH_HAS_HYPRE
#include "../generic/hypre_solver.h"
#endif
#include "../generic/general_purpose_preconditioners.h"
#ifdef OOMPH_HAS_TRILINOS
#include "../generic/trilinos_solver.h"
#endif
namespace oomph
{

 //=============================================================================
 /// \short Functions to create instances of optimal subsidiary operators for
 /// the PseudoElasticPreconditioner
 //=============================================================================
 namespace Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper
 {
#ifdef OOMPH_HAS_HYPRE
  /// \short AMG w/ GS smoothing for the augmented elastic subsidiary linear
  /// systems
  Preconditioner* get_elastic_preconditioner();
#endif
  
#ifdef OOMPH_HAS_TRILINOS
  /// \short CG with diagonal preconditioner for the lagrange multiplier
  /// subsidiary linear systems.
  Preconditioner* get_lagrange_multiplier_preconditioner();
#endif
 }



 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////



//=============================================================================
/// \short A subsidiary preconditioner for the pseudo-elastic FSI
/// preconditioner. Also a stand-alone preconditioner for the problem of
/// non-linear elasticity subject to prescribed displacement by Lagrange
/// multiplier..\n\n
/// \b Enumeration of Elastic DOF types in the Pseudo-Elastic Elements \n
/// The method get_dof_types_for_unknowns() must be implemented such that
/// DOFs subject be Lagrange multiplier and DOFs NOT subject to Lagrange
/// multiplier have different labels. For example in a 3D problem there are
/// 6 DOF types and the following labelling must be implemented:
/// 0 - x displacement (without lagr mult traction)\n
/// 1 - y displacement (without lagr mult traction)\n
/// 2 - z displacement (without lagr mult traction)\n
/// 4 - x displacement (with lagr mult traction)\n
/// 5 - y displacement (with lagr mult traction)\n
/// 6 - z displacement (with lagr mult traction)\n
//=============================================================================
 class PseudoElasticPreconditioner 
  : public BlockPreconditioner<CRDoubleMatrix>
 {
  /// \short PseudoElasticFSIPreconditioner is a friend to access the private
  /// *_preconditioner_solve(...) method
  friend class PseudoElasticFSIPreconditioner;

   public:
  
  /// \short This preconditioner includes the option to use subsidiary 
  /// operators other than SuperLUPreconditioner for this problem. 
  /// This is the typedef of a function that should return an instance
  /// of a subsidiary preconditioning operator.  This preconditioner is 
  /// responsible for the destruction of the subsidiary preconditioners.
  typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();
  
  /// \short The augmented elasticity system can be preconditioned in one
  /// of four ways. \n
  /// 0 - Exact preconditioner \n
  /// 1 - Block diagonal preconditioning \n 
  /// 2 - Block upper triangular preconditioner \n
  /// 3 - Block lower triangular preconditioner \n
  /// We group together the different components of the displacement vector 
  /// field for the block decomposition.
  enum Elastic_preconditioner_type { Exact_block_preconditioner,
                                     Block_diagonal_preconditioner,
                                     Block_lower_triangular_preconditioner,
                                     Block_upper_triangular_preconditioner };
  
  /// \short Default (and only) constructor.
  PseudoElasticPreconditioner()
   {
    // null pointers
    Lagrange_multiplier_subsidiary_preconditioner_function_pt = 0;
    Elastic_subsidiary_preconditioner_function_pt = 0;
    Elastic_preconditioner_pt = 0;
    
    // set defaults
    Use_inf_norm_of_s_scaling = true;
    E_preconditioner_type = Exact_block_preconditioner;
    
    // resize the Mesh_pt
    this->set_nmesh(2);
    Lagrange_multiplier_mesh_pt = 0;
    Elastic_mesh_pt = 0;
   }
  
  /// destructor
  virtual ~PseudoElasticPreconditioner()
   {
    this->clean_up_memory();
   }
  
  /// Broken copy constructor
  PseudoElasticPreconditioner
   (const PseudoElasticPreconditioner&)
   { 
    BrokenCopy::broken_copy("PseudoElasticPreconditioner");
   } 
  
  /// Broken assignment operator
  void operator=
   (const PseudoElasticPreconditioner&) 
   {
    BrokenCopy::broken_assign(" PseudoElasticPreconditioner");
   }
  
  /// Setup method for the PseudoElasticPreconditioner.
  void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
  
  /// \short Apply the preconditioner. Method implemented in two
  /// other methods (elastic and lagrange multiplier subsidiary
  /// preocnditioner) for the PseudoElasticFSIPreconditioner
  void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
  {
   this->elastic_preconditioner_solve(r,z);
   this->lagrange_multiplier_preconditioner_solve(r,z);
  }
  
  /// \short Access function to mesh containing the block-preconditionable
  /// elastic elements
  void set_elastic_mesh(Mesh* mesh_pt) 
  {
   Elastic_mesh_pt = mesh_pt;
  }
  
  /// \short Access function to mesh containing the block-preconditionable
  /// lagrange multiplier elements 
  void set_lagrange_multiplier_mesh(Mesh* mesh_pt) 
  {
   Lagrange_multiplier_mesh_pt = mesh_pt;
  }
  
  /// \short access function to determine whether the inf norm of S should be
  /// used as scaling. Otherwise no scaling.
  bool& use_inf_norm_of_s_scaling() 
   {
    return Use_inf_norm_of_s_scaling;
   }
  
  /// \short By default the Lagrange multiplier subsidiary systems are 
  /// preconditioner with SuperLUPreconditioner. For a different 
  /// preconditioner, pass a function to this 
  /// method returning a different subsidiary operator.
  void set_lagrange_multiplier_subsidiary_preconditioner
   (SubsidiaryPreconditionerFctPt prec_fn)
  {
   Lagrange_multiplier_subsidiary_preconditioner_function_pt = prec_fn;
  }
  
  /// \short By default the elastic subsidiary systems are 
  /// preconditioner with SuperLUPreconditioner. For a different 
  /// preconditioner, pass a function to this 
  /// method returning a different subsidiary operator.
  void set_elastic_subsidiary_preconditioner
   (SubsidiaryPreconditionerFctPt prec_fn)
  {
   Elastic_subsidiary_preconditioner_function_pt = prec_fn;
  }
  
  /// \short Set the type of preconditioner applied to the elastic: \n
  /// 0 - Exact preconditioner \n
  /// 1 - Block diagonal preconditioning \n 
  /// 2 - Block upper triangular preconditioner \n
  /// 3 - Block lower triangular preconditioner \n
  /// We group together the different components of the displacement vector 
  /// field for the block decomposition.
  Elastic_preconditioner_type& elastic_preconditioner_type()
   {
    return E_preconditioner_type;
   }
  
  /// \short Clears the memory.
  void clean_up_memory();
  
   private:
  
  /// \short Apply the elastic subsidiary preconditioner.
  void elastic_preconditioner_solve(const DoubleVector& r, DoubleVector& z);
  
  /// \short  Apply the lagrange multiplier subsidiary preconditioner.
  void lagrange_multiplier_preconditioner_solve(const DoubleVector& r,
                                                DoubleVector& z);
  
  /// The scaling. Defaults to infinite norm of S.
  double Scaling;
  
  /// \short boolean indicating whether the inf-norm of S should be used as 
  /// scaling. Default = true;
  bool Use_inf_norm_of_s_scaling;
  
  /// \short an unsigned indicating which method should be used for 
  /// preconditioning the solid component.
  Elastic_preconditioner_type E_preconditioner_type;
  
  /// \short the dimension of the problem
  unsigned Dim;
  
  /// \short storage for the preconditioner for the solid system
  Preconditioner* Elastic_preconditioner_pt;
  
  /// \short lagrange multiplier preconditioner pt
  Vector<Preconditioner*> Lagrange_multiplier_preconditioner_pt;
  
  /// the lagrange multiplier subsidary preconditioner function pointer
  SubsidiaryPreconditionerFctPt 
   Lagrange_multiplier_subsidiary_preconditioner_function_pt;
  
  /// the solid subsidiary preconditioner function pointer
  SubsidiaryPreconditionerFctPt Elastic_subsidiary_preconditioner_function_pt;
  
  /// \short the Solid mesh pt
  Mesh* Elastic_mesh_pt;
  
  /// \short the Lagrange multiplier mesh pt
  Mesh* Lagrange_multiplier_mesh_pt;
 };
 


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
 


//=============================================================================
/// Subsidiary helper preconditioner for the PseudoElasticPreconditioner.
/// Required to construct the augmented elastic system prior to
/// preconditioning.\n
/// NOTE:\n
/// 1. This is only intended to be used as a subsidiary preconditioner within
/// the PseudoElasticPreconditioner.\n
/// 2. If this preconditioner has N DOF types then the first N/2 are assumed to
/// be ordinary solid DOF types, and the second N/2 are the solid DOF types 
/// with lagrange multiplier tractions applied.\n 
/// 3. By default this preconditioner uses a superlu preconditioner.
//=============================================================================
class PseudoElasticPreconditionerSubsidiaryPreconditioner 
: public BlockPreconditioner<CRDoubleMatrix>
{
   
  public:
   
 /// \short typedef for a function that allows other preconditioners to be
 /// emplyed to solve the subsidiary linear systems. \n
 /// The function should return a pointer to the requred subsidiary
 /// preconditioner generated using new. This preconditioner is responsible
 /// for the destruction of the subsidiary preconditioners.
 typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

 /// \short
 PseudoElasticPreconditionerSubsidiaryPreconditioner()
  {
   Scaling = 1.0;
   Preconditioner_pt = 0;
   Subsidiary_preconditioner_function_pt = 0;
  }
   
 /// destructor
 ~PseudoElasticPreconditionerSubsidiaryPreconditioner()
  {
   this->clean_up_memory();
  }
   
 /// Broken copy constructor
 PseudoElasticPreconditionerSubsidiaryPreconditioner 
  (const PseudoElasticPreconditionerSubsidiaryPreconditioner &)
  { 
   BrokenCopy::broken_copy
    ("PseudoElasticPreconditionerSubsidiaryPreconditioner ");
  } 
 
 /// Broken assignment operator
 void operator=(const PseudoElasticPreconditionerSubsidiaryPreconditioner&) 
  {
   BrokenCopy::broken_assign
   (" PseudoElasticPreconditionerSubsidiaryPreconditioner");
  }
 
 // setup the preconditioner
 void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
   
 // apply the preconditioner
 void preconditioner_solve(const DoubleVector& r, DoubleVector& z);
   
 /// \short Specify the scaling. Default is 1.0  Must be called before 
 /// setup(...).
 double& scaling()
  {
   return Scaling;
  }
   
 /// access function to set the subsidiary preconditioner function.
 void set_subsidiary_preconditioner_function
  (SubsidiaryPreconditionerFctPt sub_prec_fn)
  {
   Subsidiary_preconditioner_function_pt = sub_prec_fn;
  };

  private:
   
 /// clears the memory
 void clean_up_memory()
  {
   delete Preconditioner_pt;
   Preconditioner_pt = 0;
  }
   
 // the augmentation scaling
 double Scaling;
   
 /// the preconditioner pt 
 Preconditioner* Preconditioner_pt;

 /// the SubisidaryPreconditionerFctPt 
 SubsidiaryPreconditionerFctPt Subsidiary_preconditioner_function_pt;
}; // end of PseudoElasticPreconditionerSubsidiaryPreconditioner 



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



//=============================================================================
/// Subsidiary helper preconditioner for the PseudoElasticPreconditioner.
/// Required for block preconditioner of the augmented elastic subsidiary
/// problem.
/// NOTE:\n
/// 1. This is only intended to be used as a subsidiary preconditioner within
/// the PseudoElasticPreconditioner.\n
/// 2. If this preconditioner has N DOF types then the first N/2 are assumed to
/// be ordinary solid DOF types, and the second N/2 are the solid DOF types 
/// with lagrange multiplier tractions applied.\n 
/// 3. By default this preconditioner uses a superlu preconditioner.
//=============================================================================
 class PseudoElasticPreconditionerSubsidiaryBlockPreconditioner 
  : public BlockPreconditioner<CRDoubleMatrix>
 {
   public :
  
  /// \short This preconditioner includes the option to use subsidiary 
  /// operators other than SuperLUPreconditioner for this problem. 
  /// This is the typedef of a function that should return an instance
  /// of a subsidiary preconditioning operator.  This preconditioner is 
  /// responsible for the destruction of the subsidiary preconditioners.
  typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();
  
  /// Constructor. (By default this preconditioner is upper triangular).
  PseudoElasticPreconditionerSubsidiaryBlockPreconditioner ()
   : BlockPreconditioner<CRDoubleMatrix>()
   {
    // null the subsidiary preconditioner function pointer
    Subsidiary_preconditioner_function_pt = 0;
    
    // default to block diagonal
    Method = 0;
    
    // default scaling = 1.0
    Scaling = 1.0;
   };
  
  /// Destructor
  ~PseudoElasticPreconditionerSubsidiaryBlockPreconditioner ()
   {
    this->clean_up_memory();
   }
  
  /// Broken copy constructor
  PseudoElasticPreconditionerSubsidiaryBlockPreconditioner 
   (const PseudoElasticPreconditionerSubsidiaryBlockPreconditioner &)
   { 
    BrokenCopy::broken_copy
     ("PseudoElasticPreconditionerSubsidiaryBlockPreconditioner ");
   } 
  
  /// Broken assignment operator
  void operator=
   (const PseudoElasticPreconditionerSubsidiaryBlockPreconditioner&) 
   {
    BrokenCopy::broken_assign
    (" PseudoElasticPreconditionerSubsidiaryBlockPreconditioner");
   }
  
  /// clean up the memory
  void clean_up_memory();

  /// \short Setup the preconditioner 
  void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
  
  /// Apply preconditioner to r
  void preconditioner_solve(const DoubleVector &res, DoubleVector &z);

  /// access function to set the subsidiary preconditioner function.
  void set_subsidiary_preconditioner_function
   (SubsidiaryPreconditionerFctPt sub_prec_fn)
  {
   Subsidiary_preconditioner_function_pt = sub_prec_fn;
  };
  
  /// use as a block diagonal preconditioner
  void block_diagonal()
  {
   Method = 0;
  }
  
  /// Use as an upper triangular preconditioner
  void upper_triangular() 
  {
   Method = 1;
  }
  
  /// Use as a lower triangular preconditioner
  void lower_triangular() 
  {
   Method = 2;
  }
  
  /// \short specify the scaling. Default is 1.0  Must be called before 
  /// setup(...).
  double& scaling()
  {
   return Scaling;
  }
  
   private:
  
  /// \short Vector of SuperLU preconditioner pointers for storing the 
  /// preconditioners for each diagonal block
  Vector<PseudoElasticPreconditionerSubsidiaryPreconditioner*> 
   Diagonal_block_preconditioner_pt;   
  
  /// Matrix of matrix vector product operators for the off diagonals
  DenseMatrix<MatrixVectorProduct*> Off_diagonal_matrix_vector_products;
  
  /// the preconditioning method.\n
  /// 0 - block diagonal\n
  /// 1 - upper triangular\n
  /// 2 - lower triangular\n
  unsigned Method;
  
  /// the SubisidaryPreconditionerFctPt 
  SubsidiaryPreconditionerFctPt Subsidiary_preconditioner_function_pt;
  
  /// the scaling. default 1.0.
  double Scaling;
 };



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



//=============================================================================
/// \short A helper class for PseudoElasticPreconditioner.
/// \n NOT A PRECONDITIONER.
//=============================================================================
class PseudoElasticPreconditionerScalingHelper
: public BlockPreconditioner<CRDoubleMatrix>
{

  public:

 /// The constructor.\n
 /// NOTE: \n
 /// 1. master_prec_pt should point to the 
 /// PseudoElasticPreconditioner.
 /// 2. matrix_pt should point to the jacobian.
 /// 3. The vector dof_list should contain the full list of 
 /// DOFS associated with the solid subsidiary system.
 PseudoElasticPreconditionerScalingHelper
  (Problem* problem_pt, BlockPreconditioner<CRDoubleMatrix>* 
   master_prec_pt, CRDoubleMatrix* matrix_pt, Vector<unsigned>& dof_list)
  {
   // turn into a subisiary preconditioner
   this->turn_into_subsidiary_block_preconditioner(master_prec_pt,dof_list);
     
   // all dofs are of the same block type
   Vector<unsigned> dof_to_block_map(dof_list.size(),0);
     
   // call block_setup(...)
   this->block_setup(problem_pt,matrix_pt,dof_to_block_map);

   // store the matrix_pt
   Matrix_pt = matrix_pt;
  }

 /// Destructor. Does nothing.
 ~ PseudoElasticPreconditionerScalingHelper() 
  {
   this->clear_block_preconditioner_base();
  }

 /// Broken copy constructor
 PseudoElasticPreconditionerScalingHelper
  (const PseudoElasticPreconditionerScalingHelper&)
  {
   BrokenCopy::
    broken_copy("PseudoElasticPreconditionerScalingHelper");
  }

 /// Broken assignment operator
 void operator=(const PseudoElasticPreconditionerScalingHelper&)
  {
   BrokenCopy::
    broken_assign("PseudoElasticPreconditionerScalingHelper");
  }

 /// returns the infinite norm of S
 double s_inf_norm()
  {
   CRDoubleMatrix* m_pt = 0;
   this->get_block(0,0,Matrix_pt,m_pt);
   double s_inf_norm = m_pt->inf_norm();
   delete m_pt;
   return s_inf_norm;
  }
   
 // broken preconditioner setup
 void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
  {
   std::ostringstream error_message;
   error_message
    << "This method is intentionally broken. This is not a functioning "
    << "preconditioner.";
   throw OomphLibError(
    error_message.str(),
    "PrecribedBoundaryDisplacementSubidiaryPreconditioner::setup(...)",
    OOMPH_EXCEPTION_LOCATION);
  }
   
 // broken preconditioner solve
 void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
  {
   std::ostringstream error_message;
   error_message
    << "This method is intentionally broken. This is not a functioning "
    << "preconditioner.";
   throw OomphLibError(
    error_message.str(),
    "PseudoElasticPreconditionerScalingHelper::preconditioner_solve(...)",
    OOMPH_EXCEPTION_LOCATION);
  }
   
  private:
   
 /// pointer to the Jacobian
 CRDoubleMatrix* Matrix_pt;
}; // end of PseudoElasticPreconditionerScalingHelper
}
#endif

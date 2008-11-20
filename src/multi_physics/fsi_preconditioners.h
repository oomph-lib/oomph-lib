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
#ifndef OOMPH_FSI_PRECONDITIONERS_HEADER
#define OOMPH_FSI_PRECONDITIONERS_HEADER


namespace oomph
{


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



//============================================================================ 
/// \short FSI preconditioner. This extracts upper/lower triangular
/// blocks in the 3x3 overall block matrix structure arising from
/// the monolithic discretisation of FSI problems with algebraic
/// node updates. Dofs are decomposed into fluid velocity, pressure 
/// and solid unknowns. NavierStokesLSCPreconditioner is used
/// as the inexact solver for the fluid block; SuperLU (in 
/// its incarnation as an "exact" preconditioner) is used for
/// the solid block. By default we retain the fluid on solid off
/// diagonal blocks. 
//=============================================================================
class FSIPreconditioner : public virtual BlockPreconditioner<CRDoubleMatrix>
{
 
public :
 
 /// \short Constructor: By default use block triangular form with 
 /// retained fluid on solid terms.
 FSIPreconditioner()
  {
   // set the mesh pointers
   this->Mesh_pt.resize(2);
   this->Mesh_pt[0]=0;
   this->Mesh_pt[1]=0;

   // Default setting: Fluid onto solid as it this was shown to be
   // marginally faster than solid onto fluid; see Heil CMAME 193 (2004)
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=true;

   // Create the Navier Stokes LSC preconditioner
   Navier_stokes_preconditioner_pt = new NavierStokesLSCPreconditioner;

   // Create the Solid preconditioner
#ifdef OOMPH_HAS_MPI
   Solid_preconditioner_pt = new SuperLU_dist_Preconditioner;
#else
   Solid_preconditioner_pt = new SuperLU_Preconditioner;
#endif
   
   // Preconditioner hasn't been set up yet.
   Preconditioner_has_been_setup=false;

   // Null out pointers to additional block matrices not already
   // stored in Navier Stokes block preconditioner
   Block_matrix_0_2_pt=0;
   Block_matrix_1_2_pt=0;
   Block_matrix_2_2_pt=0;
   Block_matrix_2_0_pt=0;
   Block_matrix_2_1_pt=0;

   // set Doc_time to false
   Doc_time = false;
  }
 
 
 /// Destructor: Clean up.
 ~FSIPreconditioner()
  {
   // Do what it says
   clean_up_memory();

   //Delete the Navier-Stokes preconditioner (inexact solver)
   delete Navier_stokes_preconditioner_pt;

   //Delete the solid preconditioner (inexact solver)
   delete Solid_preconditioner_pt;
  }
 
 
 /// Broken copy constructor
 FSIPreconditioner(const FSIPreconditioner&)
  {
   BrokenCopy::broken_copy("FSIPreconditioner");
  }
 
 
 /// Broken assignment operator
 void operator=(const FSIPreconditioner&)
  {
   BrokenCopy::broken_assign("FSIPreconditioner");
  }


 /// Switch to block-diagonal preconditioner
 void use_block_diagonal_version()
  {
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=false;
  }

 /// \short Switch to block-triangular preconditioner in which
 /// action of fluid dofs onto solid equations is retained
 void use_block_triangular_version_with_fluid_on_solid()
  {
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=true;
  }
 
 /// \short Switch to block-triangular preconditioner in which
 /// action of solid dofs onto fluid equations is retained
 void use_block_triangular_version_with_solid_on_fluid()
  {
   Retain_solid_onto_fluid_terms=true;
   Retain_fluid_onto_solid_terms=false;
  } 
 
 /// \short Access function to mesh containing the block-preconditionable
 /// Navier-Stokes elements. 
 Mesh*& navier_stokes_mesh_pt() 
  {
   return Mesh_pt[0];
  }

 /// \short Access function to mesh containing the block-preconditionable
 /// FSI solid elements. 
 Mesh*& wall_mesh_pt() 
  {
   return Mesh_pt[1];
  }


 /// \short Setup the preconditioner
 void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
 
 
 /// \short Apply preconditioner to r
 void preconditioner_solve(const Vector<double>&r,
                           Vector<double> &z);

 /// Access function to the Navier Stokes preconditioner (inexact solver)
 NavierStokesLSCPreconditioner* navier_stokes_preconditioner_pt() const
  {
   return Navier_stokes_preconditioner_pt;
  }

 /// Access function for Doc_time
 bool& doc_time() {return Doc_time;}
 
private:

 /// Clean up memory, delete the blocks allocated in setup
 void clean_up_memory()
  {
   if (Preconditioner_has_been_setup)
    {
     if (Retain_solid_onto_fluid_terms)
      {
       delete Block_matrix_0_2_pt;
       Block_matrix_0_2_pt=0;
       delete Block_matrix_1_2_pt;
       Block_matrix_1_2_pt=0;
      }
     if (Retain_fluid_onto_solid_terms)
      {
       delete Block_matrix_2_0_pt;
       Block_matrix_2_0_pt=0;
       delete Block_matrix_2_1_pt;
       Block_matrix_2_1_pt=0;
      }
    }
  }

 /// Pointer the Navier Stokes preconditioner (inexact solver)
 NavierStokesLSCPreconditioner* Navier_stokes_preconditioner_pt;

 /// Pointer to the solid preconditioner  (inexact solver)
 Preconditioner* Solid_preconditioner_pt;

 /// Pointer to fluid velocity/solid interaction matrix
 CRDoubleMatrix* Block_matrix_0_2_pt;

 /// Pointer to fluid pressure/solid interaction matrix
 CRDoubleMatrix* Block_matrix_1_2_pt;

 /// Pointer to solid tangent stiffness matrix
 CRDoubleMatrix* Block_matrix_2_2_pt;

 /// Pointer to solid/fluid velocity solid interaction matrix
 CRDoubleMatrix* Block_matrix_2_0_pt;

 /// Pointer to solid/fluid pressure interaction matrix
 CRDoubleMatrix* Block_matrix_2_1_pt;

 /// Boolean indicating the preconditioner has been set up 
 bool Preconditioner_has_been_setup;

 /// \short Boolean flag used to indicate that the solid onto fluid
 /// interaction terms are to be retained
 bool Retain_solid_onto_fluid_terms;

 /// \short Boolean flag used to indicate that the fluid onto solid
 /// interaction terms are to be retained
 bool Retain_fluid_onto_solid_terms;

 /// Set Doc_time to true for outputting results of timings
 bool Doc_time;

 };


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// FSI preconditioner member functions
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




//=============================================================================
/// Setup the preconditioner. Note: Matrix must be a CRDoubleMatrix.
//=============================================================================
void FSIPreconditioner::setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {

  // Clean up any existing data
  this->clean_up_memory();

#ifdef PARANOID
  if (Mesh_pt.size()!=2)
   {
    std::ostringstream error_message;
    error_message << "FSIPreconditioner needs two meshes!\n"
                  << "We have: "  << Mesh_pt.size() << std::endl;
    throw OomphLibError(error_message.str(),
                     	"FSIPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  else
   {
    if (Mesh_pt[0]==0)
     {
      std::ostringstream error_message;
      error_message << "Pointer to fluid mesh hasn't been set!\n";
      throw OomphLibError(error_message.str(),
                          "FSIPreconditioner::setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    if (Mesh_pt[1]==0)
     {
      std::ostringstream error_message;
      error_message << "Pointer to solid mesh hasn't been set!\n";
      throw OomphLibError(error_message.str(),
                          "FSIPreconditioner::setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // Call block setup for this preconditioner
  this->block_setup(problem_pt);

  // Block mapping for the subsidiary Navier Stokes preconditioner:
  // blocks 0 and 1 in the FSI preconditioner are also blocks 0 and 1
  // in the subsidiary Navier Stokes one.
  Vector<unsigned> block_lookup(2);
  block_lookup[0] = 0;
  block_lookup[1] = 1;

  // Turn the Navier Stokes LSC preconditioner into a subsidiary preconditioner
  // of this preconditioner
  Navier_stokes_preconditioner_pt->
   turn_into_subsidiary_block_preconditioner(this,block_lookup);

  // Setup the navier stokes preconditioner: Tell it about the
  // Navier Stokes mesh and set it up.
  Navier_stokes_preconditioner_pt->navier_stokes_mesh_pt() = Mesh_pt[0];
  Navier_stokes_preconditioner_pt->setup(problem_pt,matrix_pt);


  // Recast Jacobian matrix to CRDoubleMatrix
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
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

  // Extract the additional blocks we need for FSI:
  
  // Solid tangent stiffness matrix
  this->get_block(2,2,cr_matrix_pt,Block_matrix_2_2_pt);

  // Solid on fluid terms (if needed)
  if (Retain_solid_onto_fluid_terms)
   {
    this->get_block(0,2,cr_matrix_pt,Block_matrix_0_2_pt);
    this->get_block(1,2,cr_matrix_pt,Block_matrix_1_2_pt);
   }
  
  // Fluid on solid terms (if needed)
  if (Retain_fluid_onto_solid_terms)
   {
    this->get_block(2,0,cr_matrix_pt,Block_matrix_2_0_pt);
    this->get_block(2,1,cr_matrix_pt,Block_matrix_2_1_pt);
   }
  

  // Setup the solid preconditioner (inexact solver)

#ifdef OOMPH_HAS_MPI   
  double t_start = MPI_Wtime();
#else
  clock_t t_start = clock();
#endif

  Solid_preconditioner_pt->setup(problem_pt,Block_matrix_2_2_pt);
 
#ifdef OOMPH_HAS_MPI   
  double t_end = MPI_Wtime();
  double setup_time= t_end-t_start;
#else
 clock_t t_end = clock();
 double setup_time=double(t_end-t_start)/CLOCKS_PER_SEC;
#endif

 // Delete the solid matrix
 delete Block_matrix_2_2_pt;
 Block_matrix_2_2_pt = 0;

 // Output times
 if(Doc_time)
  {
   oomph_info << "Solid sub-preconditioner setup time [sec]: "
              << setup_time << "\n";
   
   // Doc density of retained interaction block
   if (Retain_solid_onto_fluid_terms)
    {
     oomph_info << "Fill level of solid on fluid blocks (C_us and C_ps): " <<
      double(Block_matrix_0_2_pt->nnz())/
      double(Block_matrix_0_2_pt->nrow()*
             Block_matrix_0_2_pt->ncol())*100.0
                << "%  " <<
      double(Block_matrix_1_2_pt->nnz())/
      double(Block_matrix_1_2_pt->nrow()*
             Block_matrix_1_2_pt->ncol())*100.0
                << "%  " << std::endl;
    }
   
   // Doc density of retained interaction block
   if (Retain_fluid_onto_solid_terms)
    {
     oomph_info << "Fill level of fluid on solid blocks (C_su and C_sp): " <<
      double(Block_matrix_2_0_pt->nnz())/
      double(Block_matrix_2_0_pt->nrow()*
             Block_matrix_2_0_pt->ncol())*100.0
                << "%  " <<
      double(Block_matrix_2_1_pt->nnz())/
      double(Block_matrix_2_1_pt->nrow()*
             Block_matrix_2_1_pt->ncol())*100.0
                << "%  " << std::endl;
    }
  }
 
  // We're done (and we stored some data)
  Preconditioner_has_been_setup=true;
  
 }


//======================================================================
/// Apply preconditioner to Vector r
//======================================================================
void FSIPreconditioner::preconditioner_solve(const Vector<double>&r, 
                                             Vector<double> &z)
{
 // Block sizes
 unsigned n_veloc_dof= this->block_dimension(0);
 unsigned n_press_dof= this->block_dimension(1);
 unsigned n_solid_dof= this->block_dimension(2);
 

 // Make copy of residual vector (to overcome const-ness
 Vector<double> res(r);


 // Retain off-diagonals that represent effect of solid on fluid
 //-------------------------------------------------------------
 if (Retain_solid_onto_fluid_terms)
  {

   // Working vectors
   Vector<double> temp_solid_vec(n_solid_dof);
   Vector<double> temp_veloc_vec(n_veloc_dof);
   Vector<double> temp_press_vec(n_press_dof);

   // Copy solid values from residual to temp_vec:
   // Loop over all entries in the global vector (this one 
   // includes solid, velocity and pressure dofs in some random fashion)
   get_block_vector(2,res,temp_solid_vec);
    
   // Solve solid system by back-substitution
   // with LU-decomposed stiffness matrix   
   Solid_preconditioner_pt->preconditioner_solve(temp_solid_vec,
                                                 temp_solid_vec);
   
   // NOTE: temp_solid_vec now contains z_s = S^{-1} r_s

   // Multiply C_{us} by z_s
   Block_matrix_0_2_pt->multiply(temp_solid_vec,temp_veloc_vec);
      
   // Multiply C_{ps} by z_s
   Block_matrix_1_2_pt->multiply(temp_solid_vec,temp_press_vec);
   
   
   // Subtract from fluid residual vector for fluid solve
   // and copy solid result back into solution vector
   for (unsigned i=0; i<Nrow; i++)
    {
     unsigned i_block=this->block_number(i);
     
     switch (i_block)
      {
       // Subtract from velocity residual
      case 0:
       res[i]-=temp_veloc_vec[index_in_block(i)];
       break;
       
       // Subtract from pressure residual
      case 1:

       res[i]-=temp_press_vec[index_in_block(i)];
       break;

       // Copy into solid result:
      case 2:
       z[i] = temp_solid_vec[ this->index_in_block(i) ];
       break;

      default:
       std::ostringstream error_message;
       error_message << "Never get here: i_block= " << i_block << std::endl;
       throw OomphLibError(error_message.str(),
                           "FSIPreconditioner::preconditioner_solve()",
                           OOMPH_EXCEPTION_LOCATION);
       break;

      }
    }

   // now apply the navier stokes lsc preconditioner
   Navier_stokes_preconditioner_pt->preconditioner_solve(res,z);
  }
 

 // Retain off-diagonals that represent effect of fluid on solid
 //-------------------------------------------------------------
 // (or diagonal preconditioner)
 //-----------------------------
 else
  {

   // Call fluid preconditioner for fluid block
   Navier_stokes_preconditioner_pt->preconditioner_solve(res,z);

   // Working vectors
   Vector<double> temp_solid_vec(n_solid_dof);
   Vector<double> temp_veloc_vec;
   Vector<double> temp_press_vec;

   // Only allocate storage if needed
   if (Retain_fluid_onto_solid_terms)
    {
     temp_veloc_vec.resize(n_veloc_dof);
     temp_press_vec.resize(n_press_dof);
    }

   // Copy solid values from residual or solution vectors to temp_vecs:
   // Loop over all entries in the global vector (this one 
   // includes solid, velocity and pressure dofs in some random fashion)
   // Copy solid values from residual or solution vectors to temp_vecs:
   if (Retain_fluid_onto_solid_terms)
    {
     get_block_vector(0,z,temp_veloc_vec);
     get_block_vector(1,z,temp_press_vec);
    }
   get_block_vector(2,res,temp_solid_vec);

   // Do matrix vector products with fluid onto solid coupling matrices:
   if (Retain_fluid_onto_solid_terms)
    {
     // Auxiliary vector to hold the matrix vector product of the
     // fluid-onto-solid coupling matrices with the fluid solutions:
     Vector<double> aux_vec(n_solid_dof);
     
     // Multiply C_{su} by z_u
     Block_matrix_2_0_pt->multiply(temp_veloc_vec, aux_vec);
     
     // ...and subtract from r_s:
     for (unsigned i=0;i<n_solid_dof;i++)
      {
       temp_solid_vec[i]-=aux_vec[i];
      }
     
     
     // Multiply C_{sp} by z_p
     Block_matrix_2_1_pt->multiply(temp_press_vec, aux_vec);
     
     // ...and subtract from r_s:
     for (unsigned i=0;i<n_solid_dof;i++)
      {
       temp_solid_vec[i]-=aux_vec[i];
      }
    }
      
   // Solve solid system by back-substitution
   // with LU-decomposed stiffness matrix   
   Solid_preconditioner_pt->preconditioner_solve(temp_solid_vec,
                                                 temp_solid_vec);

   // Now copy result_vec (i.e. z_s) back into the global vector z.
   // Loop over all entries in the global results vector z:
   return_block_vector(2,temp_solid_vec,z);   
   for (unsigned i=0; i<Nrow; i++)
    {
     // If this entry is a solid dof
     if ( this->block_number(i) == 2 )
      {
       // copy the result across
       z[i] = temp_solid_vec[ this->index_in_block(i) ];
      }
    }   
  }
}



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//============================================================================ 
/// \short FSI preconditioner. This extracts upper/lower triangular
/// blocks in the 3x3 overall block matrix structure arising from
/// the monolithic discretisation of FSI problems with algebraic
/// node updates. Dofs are decomposed into fluid velocity, pressure 
/// and solid unknowns. Blocks are then re-assembled into one global
/// matrix and solved with a direct solver (SuperLU in its incarnation
/// as an exact preconditioner). By default we retain the fluid on solid off
/// diagonal blocks. 
//=============================================================================
template<typename MATRIX>
class SimpleFSIPreconditioner : public BlockPreconditioner<MATRIX>
{
 
public :
 
 /// Constructor. 
 SimpleFSIPreconditioner() : BlockPreconditioner<MATRIX>()
  {
   this->Mesh_pt.resize(2);
   this->Mesh_pt[0]=0;
   this->Mesh_pt[1]=0;

   // Default setting: Retain fluid on solid
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=true;

   // Null the preconditioner pointer (inexact solver)
   Preconditioner_pt = 0;
  }
 
 
 /// Destructor: Clean up
 ~SimpleFSIPreconditioner()
  {
   // Wiping preconditioner (inexact solver) wipes the stored 
   // LU decompositon
   if (Preconditioner_pt != 0)
    {
     delete Preconditioner_pt;
     Preconditioner_pt = 0;
    }
  }
 
 
 /// Broken copy constructor
 SimpleFSIPreconditioner(const SimpleFSIPreconditioner&)
  {
   BrokenCopy::broken_copy("SimpleFSIPreconditioner");
  }
 
 
 /// Broken assignment operator
 void operator=(const SimpleFSIPreconditioner&)
  {
   BrokenCopy::broken_assign("SimpleFSIPreconditioner");
  }
 
 
 /// \short Access function to mesh containing the block-preconditionable
 /// Navier-Stokes elements. 
 Mesh*& navier_stokes_mesh_pt() 
  {
   return this->Mesh_pt[0];
  }

 /// \short Access function to mesh containing the block-preconditionable
 /// FSI solid elements. 
 Mesh*& wall_mesh_pt() 
  {
   return this->Mesh_pt[1];
  }

 /// \short Setup the preconditioner
 void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt);
  
 /// \short Apply preconditioner to r
 void preconditioner_solve(const Vector<double>&r,
                           Vector<double> &z);

 /// Switch to block-diagonal preconditioner
 void use_block_diagonal_version()
  {
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=false;
  }

 /// \short Switch to block-triangular preconditioner in which
 /// action of fluid dofs onto solid equations is retained
 void use_block_triangular_version_with_fluid_on_solid()
  {
   Retain_solid_onto_fluid_terms=false;
   Retain_fluid_onto_solid_terms=true;
  }
 
 /// \short Switch to block-triangular preconditioner in which
 /// action of solid dofs onto fluid equations is retained
 void use_block_triangular_version_with_solid_on_fluid()
  {
   Retain_solid_onto_fluid_terms=true;
   Retain_fluid_onto_solid_terms=false;
  }
 
private:
 
 /// \short Preconditioner (inexact solver)
 Preconditioner* Preconditioner_pt;
 
 /// \short Boolean flag used to indicate that the solid onto fluid
 /// interaction terms are to be retained
 bool Retain_solid_onto_fluid_terms;

 /// \short Boolean flag used to indicate that the fluid onto solid
 /// interaction terms are to be retained
 bool Retain_fluid_onto_solid_terms;

 /// \short Identify the required blocks: Here we only need 
 /// the momentum, gradient and divergence blocks of the
 /// 2x2 block-structured fluid matrix, the 1x1 solid block
 /// and the selected FSI-off diagonals.
 virtual void identify_required_blocks(DenseMatrix<bool>& required_blocks);
 
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// FSI preconditioner member functions
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//===========================================================================
/// Identify the required blocks: Here we only need
/// the momentum, gradient and divergence blocks of the
/// 2x2 block-structured fluid matrix, the 1x1 solid block
/// and the selected FSI-off diagonals.
//===========================================================================
template<typename MATRIX>
void SimpleFSIPreconditioner<MATRIX>::identify_required_blocks(
  DenseMatrix<bool>& required_blocks)
{

 // find number of block types
 unsigned n_block = this->Nblock_types;
 
 // Initialise all blocks to false
 for (unsigned i=0;i<n_block;i++)
  {
   for (unsigned j=0;j<n_block;j++)
    {
     required_blocks(i,j)=false;
    }
  }
 
 // Fluid: Only need momentum, gradient and divergence blocks
 required_blocks(0,0) = true;
 required_blocks(1,0) = true;
 required_blocks(0,1) = true;
  
 // Always retain the solid block
 required_blocks(2,2) = true;

 // Switch on the required off-diagonals
 if (Retain_solid_onto_fluid_terms)
  {
   required_blocks(0,2)=true;
   required_blocks(1,2)=true;
  }
 if (Retain_fluid_onto_solid_terms)
  {
   required_blocks(2,0)=true;
   required_blocks(2,1)=true;
   if (Retain_solid_onto_fluid_terms)
    {
     std::ostringstream error_message;
     error_message << "Can't retain all off-diagonal blocks!\n";
     throw OomphLibError(error_message.str(),
                         "SimpleFSIPreconditioner::required_blocks()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 
}


//=============================================================================
/// Setup the preconditioner: Copy the upper/lower triangular 
/// block matrices back into a big matrix (with the entries 
/// re-ordered relative to the original Jacobian matrix). 
//=============================================================================
 template<typename MATRIX>
 void SimpleFSIPreconditioner<MATRIX>::
 setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
 {

  // Clean up memory
  if (Preconditioner_pt != 0)
   {
    delete Preconditioner_pt;
    Preconditioner_pt = 0;
   }
   
#ifdef PARANOID
  if (this->Mesh_pt.size()!=2)
   {
    std::ostringstream error_message;
    error_message << "SimpleFSIPreconditioner needs two meshes!\n"
                  << "We have: "  << this->Mesh_pt.size() << std::endl;
    throw OomphLibError(error_message.str(),
                     	"SimpleFSIPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  else
   {
    if (this->Mesh_pt[0]==0)
     {
      std::ostringstream error_message;
      error_message << "Pointer to fluid mesh hasn't been set!\n";
      throw OomphLibError(error_message.str(),
                          "SimpleFSIPreconditioner::setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    if (this->Mesh_pt[1]==0)
     {
      std::ostringstream error_message;
      error_message << "Pointer to solid mesh hasn't been set!\n";
      throw OomphLibError(error_message.str(),
                          "SimpleFSIPreconditioner::setup()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // Extract the upper (or lower) 3x3 blocks

  // Set up the blocks look up schemes
  this->block_setup(problem_pt);

  // find number of block types
  unsigned n_block = this->Nblock_types;

  // Create matrix that indicates which blocks are required
  DenseMatrix<bool> required_blocks(n_block,n_block);

  // Identify required blocks
  identify_required_blocks(required_blocks);

  // Get pointers to the blocks
  DenseMatrix<MATRIX*> block_matrix_pt;
  this->get_blocks(dynamic_cast<MATRIX*>(matrix_pt),
                   required_blocks,
                   block_matrix_pt);

  // Build a big matrix from the blocks
  MATRIX* P_matrix_pt=0;
  build_preconditioner_matrix(block_matrix_pt,P_matrix_pt);

  // Delete blocks -- they're no longer needed.
  for (unsigned i = 0 ; i < n_block; i++)
   {
    for (unsigned j = 0 ; j < n_block; j++)
     {
      delete block_matrix_pt(i,j);
      block_matrix_pt(i,j) = 0;
     }
   }

  // Setup preconditoner (i.e. inexact solver) -- does the LU decomposition
  Preconditioner_pt = new SuperLU_Preconditioner;
  Preconditioner_pt->setup(problem_pt,P_matrix_pt);

  // Now throw away the big matrix
  delete P_matrix_pt;
 
 }


//======================================================================
/// Apply preconditioner to Vector r
//======================================================================
template<typename MATRIX>
void SimpleFSIPreconditioner<MATRIX>::
preconditioner_solve(const Vector<double>&r, Vector<double> &z)
{ 
 // create a temporary vector to hold the result of preconditioning
 Vector<double> temp_vec(this->Nrow,0.0);
   
 // get the reordered vector
 this->get_block_ordered_preconditioner_vector(r,temp_vec);
 
 // apply preconditioner to z and store in r
 Preconditioner_pt->preconditioner_solve(temp_vec,temp_vec);
 
 // copy the solution back
 this->return_block_ordered_preconditioner_vector(temp_vec,z);
}



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


}

#endif

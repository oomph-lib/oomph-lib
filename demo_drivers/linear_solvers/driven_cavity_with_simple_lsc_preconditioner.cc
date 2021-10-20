//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
//Driver for 2D rectangular driven cavity

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=============================================================================
///  Simple implementation of LSC preconditioner.\n
/// This is demo code that shows how write preconditioners. There's
/// a more sophisticated version of this preconditioner in the 
/// library!
//=============================================================================
class SimpleLSCPreconditioner : public BlockPreconditioner<CRDoubleMatrix>
{
 
public :
 
 /// Constructor - sets defaults for control flags
 NavierStokesSchurComplementPreconditioner(Problem* problem_pt) :
  BlockPreconditioner<CRDoubleMatrix>(), Problem_pt(problem_pt)
  {

   // resize the mesh pt
   // note: meaningless if subsidiary preconditioner
   this->set_nmesh(1);
   Navier_stokes_mesh_pt = 0;
   
   // hierher inialise above
   P_preconditioner_pt = 0;
   F_preconditioner_pt = 0;
   Bt_mat_vec_pt = 0;
   F_mat_vec_pt = 0;
   
  }
 
 /// Destructor
 ~NavierStokesSchurComplementPreconditioner()
  {
   clean_up_memory();
  }
 
 /// Broken copy constructor
 NavierStokesSchurComplementPreconditioner(
  const NavierStokesSchurComplementPreconditioner&)
  {
   BrokenCopy::broken_copy("NavierStokesSchurComplementPreconditioner");
  }
 

 /// Broken assignment operator
 void operator=(const NavierStokesSchurComplementPreconditioner&)
  {
   BrokenCopy::broken_assign("NavierStokesSchurComplementPreconditioner");
  }
 
 /// Setup the preconditioner
 void setup();
 
 ///  for some reason we have to remind the compiler that there is a
 /// setup() function in Preconditioner base class.
 // hierher using Preconditioner::setup;
 
 /// Apply preconditioner to Vector r
 void preconditioner_solve(const DoubleVector&r, DoubleVector &z);
 
 /// Specify the mesh containing the block-preconditionable Navier-Stokes
 /// elements.
 void set_navier_stokes_mesh(Mesh* mesh_pt)
  { 
   Navier_stokes_mesh_pt = mesh_pt;
  }
 
 ///  Helper function to delete preconditioner data.
 void clean_up_memory();
 
 
private:
 
 // oomph-lib objects
 // -----------------
 
 // Pointers to preconditioner (=inexact solver) objects
 // -----------------------------------------------------
 /// Pointer to the 'preconditioner' for the pressure matrix
 Preconditioner* P_preconditioner_pt;
 
 /// Pointer to the 'preconditioner' for the F matrix
 Preconditioner* F_preconditioner_pt;
 
 
 // hierher
 // ///  Helper function to assemble the inverse diagonals of the pressure
 // /// and velocity mass matrices from the elemental contributions defined in
 // /// NavierStokesEquations<DIM>.
 // /// If do_both=true, both are computed, otherwise only the velocity
 // /// mass matrix (the LSC version of the preconditioner only needs
 // /// that one)
 // void assemble_inv_press_and_veloc_mass_matrix_diagonal(
 //  CRDoubleMatrix*& inv_p_mass_pt, 
 //  CRDoubleMatrix*& inv_v_mass_pt, 
 //  const bool& do_both);
   
 ///  Boolean indicating whether the momentum system preconditioner 
 /// is a block preconditioner
 bool F_preconditioner_is_block_preconditioner; // hierher
 
 /// MatrixVectorProduct operator for Qv^{-1} Bt 
 MatrixVectorProduct* QBt_mat_vec_pt;
 
 /// MatrixVectorProduct operator for Bt
 MatrixVectorProduct* Bt_mat_vec_pt;
 
 /// MatrixVectorProduct operator for F
 MatrixVectorProduct* F_mat_vec_pt;
 
 ///  the pointer to the mesh of block preconditionable Navier
 /// Stokes elements. 
 Mesh* Navier_stokes_mesh_pt;
 
};


//===========================================================================
/// Setup the least-squares commutator Navier Stokes preconditioner. This
/// extracts blocks corresponding to the velocity and pressure unknowns,
/// creates the matrices actually needed in the application of the
/// preconditioner and deletes what can be deleted... Note that
/// this preconditioner needs a CRDoubleMatrix.
//============================================================================
 void NavierStokesSchurComplementPreconditioner::setup()
 {

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // NOTE: In the interest of minimising memory usage, several containers
  //       are recycled, therefore their content/meaning changes
  //       throughout this function. The code is carefully annotated
  //       but you'll have to read it line by line!
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // make sure any old data is deleted
  clean_up_memory(); // hierher

#ifdef PARANOID
  // paranoid check that the navier stokes mesh pt has been set
  if (Navier_stokes_mesh_pt == 0)
   {
    std::ostringstream error_message;
    error_message << "The navier stokes elements mesh pointer must be set.\n"
                  << "Use method set_navier_stokes_mesh(...)";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // set the mesh
  this->set_nmesh(1);
  this->set_mesh(0,Navier_stokes_mesh_pt);
  
  // Get blocks
  // ----------

  // In comes the current Jacobian. Recast it to a CR double matrix;
  // shout if that can't be done.
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());


#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message 
     << "NavierStokesSchurComplementPreconditioner only works with "
     << "CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Set up block look up schemes 
  //-----------------------------

  // this preconditioner has two types of block:
  // type 0: velocity - corresponding to DOFs 0 to n-2
  // type 1: pressure - corresponding to DOF n-1
  unsigned ndof_types = 0;
  if (this->is_subsidiary_block_preconditioner())
   {
    ndof_types = this->ndof_types();
   }
  else
   {
    ndof_types = this->ndof_types_in_mesh(0);
   }

  // hierher
  Vector<unsigned> dof_to_block_map(ndof_types);
  dof_to_block_map[ndof_types-1]=1;
  this->block_setup(dof_to_block_map);
  
  // Get B (the divergence block)
  CRDoubleMatrix* b_pt = 0;
  this->get_block(1,0,b_pt);
  
  // // get the inverse velocity and pressure mass matrices
  // CRDoubleMatrix* inv_v_mass_pt = 0;
  // CRDoubleMatrix* inv_p_mass_pt = 0;
  // double ivmm_assembly_start_t = TimingHelpers::timer();
  // if (Use_LSC)
  //  {
  //   // We only need the velocity mass matrix
  //   bool do_both=false;
  //   assemble_inv_press_and_veloc_mass_matrix_diagonal(inv_p_mass_pt,
  //                                                     inv_v_mass_pt,
  //                                                     do_both);
  //  }
  // else
  //  {
  //   // We only need both mass matrices
  //   bool do_both=true;
  //   assemble_inv_press_and_veloc_mass_matrix_diagonal(inv_p_mass_pt,
  //                                                     inv_v_mass_pt,
  //                                                     do_both);
  //  }

  
  // Get gradient matrix Bt
  CRDoubleMatrix* bt_pt = 0;
  this->get_block(0,1,bt_pt);

  // Multiply inverse velocity mass matrix by gradient matrix B^T
  double t_QBt_matrix_start = TimingHelpers::timer();
  CRDoubleMatrix* qbt_pt = new CRDoubleMatrix;
  inv_v_mass_pt->multiply(*bt_pt, *qbt_pt);
  delete bt_pt;

  // Store product in bt_pt 
  bt_pt = qbt_pt;
  delete inv_v_mass_pt;
    
  // Multiply B from left by divergence matrix B and store result in 
  // pressure Poisson matrix.
  b_pt->multiply(*bt_pt, *p_matrix_pt);

  // Kill divergence matrix because we don't need it any more
  delete b_pt;
  
  // Build the matvec operator for QBt
  QBt_mat_vec_pt = new MatrixVectorProduct;
  //QBt_mat_vec_pt->setup(bt_pt);
  this->setup_matrix_vector_product(QBt_mat_vec_pt,
                                    bt_pt, 1);
  
  // Kill gradient matrix B^T (it's been overwritten anyway and
  // needs to be recomputed afresh below)
  delete bt_pt;
  
  // Get momentum block F
  CRDoubleMatrix* f_pt = 0;
  this->get_block(0,0,f_pt);
    
  // form the matrix vector product helper
  F_mat_vec_pt = new MatrixVectorProduct;
  //F_mat_vec_pt->setup(f_pt);
  this->setup_matrix_vector_product(F_mat_vec_pt,f_pt,0);
  
  // if F is a block preconditioner then we can delete the F matrix
  if (F_preconditioner_is_block_preconditioner)
   {
    delete f_pt;
   }
  
  // Rebuild Bt (remember that we temporarily overwrote
  // it by its product with the inverse velocity mass matrix)
  bt_pt = 0;
  this->get_block(0,1,bt_pt);
  
  // Form the matrix vector operator for Bt
  Bt_mat_vec_pt = new MatrixVectorProduct;
  //Bt_mat_vec_pt->setup(bt_pt);
  this->setup_matrix_vector_product(Bt_mat_vec_pt,bt_pt,1);
  delete bt_pt;
  
  // If the P preconditioner has not been setup
  P_preconditioner_pt = new SuperLUPreconditioner;
  P_preconditioner_pt->setup(p_matrix_pt,comm_pt());
  delete p_matrix_pt;
  p_matrix_pt=0;

  // Set up solver for solution of system with momentum matrix
  // ----------------------------------------------------------
  
  F_preconditioner_pt = new SuperLUPreconditioner;
  Using_default_f_preconditioner = true;

  
  // if F is a block preconditioner
  double t_f_prec_start = TimingHelpers::timer();
  if (F_preconditioner_is_block_preconditioner)
   {
    unsigned ndof_types = this->ndof_types();
    ndof_types--;
    Vector<unsigned> dof_map(ndof_types);
    for (unsigned i = 0; i < ndof_types; i++)
     {
      dof_map[i] = i;
     }
    F_block_preconditioner_pt->
     turn_into_subsidiary_block_preconditioner(this,dof_map);

    // Set the mesh in the subsidiary preconditioner. 
    // RAYRAY This is incorrect. We never set the mesh for a subsidiary 
    // block preconditioner. The self test still passes since this part of
    // the if-else condition is never executed.
    F_block_preconditioner_pt->set_nmesh(1);
    F_block_preconditioner_pt->set_mesh(0, Navier_stokes_mesh_pt);

    F_block_preconditioner_pt->setup(matrix_pt(),comm_pt());
   }
  // otherwise F is not a block preconditioner
  else
   {
    F_preconditioner_pt->setup(f_pt,comm_pt());
    delete f_pt;
   }
  double t_f_prec_finish = TimingHelpers::timer();
  if(Doc_time)
   {
    double t_f_prec_time = t_f_prec_finish - t_f_prec_start;
    oomph_info << "F sub-preconditioner setup time [sec]: "
               << t_f_prec_time << "\n";
   }
  
  // Remember that the preconditioner has been setup so
  // the stored information can be wiped when we
  // come here next...
  Preconditioner_has_been_setup = true;
 }



//=======================================================================
 /// Apply preconditioner to r.
//=======================================================================
 void NavierStokesSchurComplementPreconditioner:: 
 preconditioner_solve(const DoubleVector &r, DoubleVector &z)
 {

#ifdef PARANOID
  if (Preconditioner_has_been_setup==false)
   {
    std::ostringstream error_message;
    error_message << "setup must be called before using preconditioner_solve";
    throw OomphLibError(
     error_message.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
  if (z.built())
   {
    if (z.nrow() != r.nrow())
     {
      std::ostringstream error_message;
      error_message << "The vectors z and r must have the same number of "
                    << "of global rows";
      throw OomphLibError(
       error_message.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);      
     }
   }
#endif

  // if z is not setup then give it the same distribution
  if (!z.distribution_pt()->built())
   {
    z.build(r.distribution_pt(),0.0);
   }

  // Step 1 - apply approximate Schur inverse to pressure unknowns (block 1)
  // -----------------------------------------------------------------------

  // Working vectors
  DoubleVector temp_vec;
  DoubleVector another_temp_vec;
  DoubleVector yet_another_temp_vec;

  // Copy pressure values from residual vector to temp_vec:
  // Loop over all entries in the global vector (this one
  // includes velocity and pressure dofs in some random fashion)
  this->get_block_vector(1,r,temp_vec);

  // NOTE: The vector temp_vec now contains the vector r_p.


  // LSC version
  if (Use_LSC)
   {
    
    // Solve first pressure Poisson system
#ifdef PARANOID
    // check a solver has been set
    if (P_preconditioner_pt==0)
     {
      std::ostringstream error_message;
      error_message << "P_preconditioner_pt has not been set.";
      throw OomphLibError(
       error_message.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // use some Preconditioner's preconditioner_solve function
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains the vector P^{-1} r_p
    
    // Multiply another_temp_vec by matrix E and stick the result into temp_vec
    temp_vec.clear();  
    QBt_mat_vec_pt->multiply(another_temp_vec, temp_vec);
    another_temp_vec.clear();
    F_mat_vec_pt->multiply(temp_vec,another_temp_vec);
    temp_vec.clear();
    QBt_mat_vec_pt->multiply_transpose(another_temp_vec, temp_vec);
    
    
    // NOTE: The vector temp_vec now contains E P^{-1} r_p
    
    // Solve second pressure Poisson system using preconditioner_solve
    another_temp_vec.clear();
    P_preconditioner_pt->preconditioner_solve(temp_vec, another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains z_p = P^{-1} E P^{-1} r_p
    //       as required (apart from the sign which we'll fix in the
    //       next step.

   }
  // Fp version
  else
   {
      
    // Multiply temp_vec by matrix E and stick the result into 
    // yet_another_temp_vec
    E_mat_vec_pt->multiply(temp_vec,yet_another_temp_vec);
        
    // NOTE: The vector yet_another_temp_vec now contains Fp Qp^{-1} r_p

    // Solve pressure Poisson system
#ifdef PARANOID
    // check a solver has been set
    if (P_preconditioner_pt==0)
     {
      std::ostringstream error_message;
      error_message << "P_preconditioner_pt has not been set.";
      throw OomphLibError(
       error_message.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Solve second pressure Poisson system using preconditioner_solve
    another_temp_vec.clear();
    P_preconditioner_pt->preconditioner_solve(yet_another_temp_vec, 
                                              another_temp_vec);
    
    // NOTE: The vector another_temp_vec now contains 
    //       z_p = P^{-1} Fp Qp^{-1} r_p
    //       as required (apart from the sign which we'll fix in the
    //       next step.

   }

  // Now copy another_temp_vec (i.e. z_p) back into the global vector z.
  // Loop over all entries in the global results vector z:
  temp_vec.build(another_temp_vec.distribution_pt(),0.0);
  temp_vec -= another_temp_vec;
  return_block_vector(1,temp_vec,z);
  

    
  // Step 2 - apply preconditioner to velocity unknowns (block 0)
  // ------------------------------------------------------------
  
  // Recall that another_temp_vec (computed above) contains the
  // negative of the solution of the Schur complement systen, -z_p.
  // Multiply by G (stored in Block_matrix_pt(0,1) and store
  // result in temp_vec (vector resizes itself).
  temp_vec.clear();
  Bt_mat_vec_pt->multiply(another_temp_vec, temp_vec);

  // NOTE: temp_vec now contains -G z_p

  // The vector another_temp_vec is no longer needed -- re-use it to store
  // velocity quantities:
  another_temp_vec.clear();

  // Loop over all enries in the global vector and find the
  // entries associated with the velocities:
  get_block_vector(0,r,another_temp_vec);
  another_temp_vec += temp_vec;
 
  // NOTE:  The vector another_temp_vec now contains r_u - G z_p

  // Solve momentum system
#ifdef PARANOID
  // check a solver has been set
  if (F_preconditioner_pt==0)
   {
    std::ostringstream error_message;
    error_message << "F_preconditioner_pt has not been set."; 
    throw OomphLibError(
     error_message.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // use some Preconditioner's preconditioner solve
  // and return
  if (F_preconditioner_is_block_preconditioner)
   {
    return_block_vector(0,another_temp_vec,z);
    F_preconditioner_pt->preconditioner_solve(z,z);
   }
  else
   {
    F_preconditioner_pt->preconditioner_solve(another_temp_vec, temp_vec);
    return_block_vector(0,temp_vec,z);
   }
 }



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


// Namespace extension
namespace oomph
{

//===================================================================
/// Overload Navier-Stokes Taylor Hood element to make block 
/// preconditionable (again -- simply overwrites code that
/// already exists in library to make it more visible for 
/// documentation)
//===================================================================
template<unsigned DIM>
class MyTaylorHoodElement :  public virtual QTaylorHoodElement<DIM>
{
 
 ///  The number of "DOF types" that degrees of freedom in this element
 /// are sub-divided into: Velocity and pressure.
 unsigned ndof_types() const
  {
   return DIM+1;
  }
 
/// Create a list of pairs for all unknowns in this element,
/// so the first entry in each pair contains the global equation
/// number of the unknown, while the second one contains the number
/// of the "DOF type" that this unknown is associated with.
/// (Function can obviously only be called if the equation numbering
/// scheme has been set up.)
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // number of nodes
   unsigned n_node = this->nnode();
   
   // temporary pair (used to store dof lookup prior to being added to list)
   std::pair<unsigned,unsigned> dof_lookup;
   
   // loop over the nodes
   for (unsigned j=0;j<n_node;j++)
    {
     // Number of values
     unsigned nvalue=this->node_pt(j)->nvalue();

     //loop over these values
     for (unsigned i=0;i<nvalue;i++)
      {
       // determine local eqn number
       int local_eqn_number = this->nodal_local_eqn(j,i);
       
       // ignore pinned values - far away degrees of freedom resulting 
       // from hanging nodes can be ignored since these are be dealt
       // with by the element containing their master nodes
       if (local_eqn_number >= 0)
        {
         // store dof lookup in temporary pair: Global equation number
         // is the first entry in pair
         dof_lookup.first = this->eqn_number(local_eqn_number);
         
         // set dof numbers: Dof number is the second entry in pair
         dof_lookup.second = i;
         
         // add to list
         dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }

};

} // end namespace extension

/////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////// 



//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100;

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}
#endif


//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class RectangularDrivenCavityProblem : public Problem
{

public:


 ///  Constructor: Specify multiplier for number of element 
 /// rows/columns and solver flag.
 RectangularDrivenCavityProblem();
 
 /// Destructor
 ~RectangularDrivenCavityProblem()
  {
   // Kill oomph-lib iterative linear solver
   if (Solver_pt!=0) delete Solver_pt;
   
   // Kill preconditioner
   if (Prec_pt!=0) delete Prec_pt;
   
   // Kill inexact solver for P block
   if  (P_matrix_preconditioner_pt!=0) delete P_matrix_preconditioner_pt;
   
   // Kill inexact solver for F block
   if  (F_matrix_preconditioner_pt!=0) delete F_matrix_preconditioner_pt;

   // Kill mesh
   delete mesh_pt();

  };

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}


 ///  Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
 {
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along the other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
 } // end_of_actions_before_newton_solve

 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }


 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


private:

 /// oomph-lib iterative linear solver
 IterativeLinearSolver* Solver_pt;
 
 /// Preconditioner
 NavierStokesSchurComplementPreconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RectangularDrivenCavity problem: Specify 
/// multiplier for number of element rows/columns and solver flag.
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::RectangularDrivenCavityProblem()
{ 
 unsigned element_multiplier=1;
 bool use_iterative_solver=true;
 bool use_hypre_for_pressure=false;
 bool use_hypre_for_momentum=false;
 bool use_block_diagonal_for_momentum=false;

 // Initialise pointer to oomph-lib iterative linear solver
 Solver_pt=0;
 
 // Initialise pointer to Preconditioner
 Prec_pt=0;
 
 // Initialise pointer to inexact solver for P block
 P_matrix_preconditioner_pt=0;
 
 // Initialise pointer to inexact solver for F block
 F_matrix_preconditioner_pt=0;
   

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10*element_multiplier;

 // # of elements in y-direction
 unsigned n_y=10*element_multiplier;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);


 // Setup iterative linear solver if required
 if (use_iterative_solver)
  {
   // Create oomph-lib iterative linear solver
   Solver_pt=new GMRES<CRDoubleMatrix>;
   
   // Set linear solver
   linear_solver_pt() = Solver_pt;
   
   // Set preconditioner
   Prec_pt=new NavierStokesSchurComplementPreconditioner(this);
   Prec_pt->set_navier_stokes_mesh(this->mesh_pt());

//   Prec_pt->navier_stokes_mesh_pt() = this->mesh_pt();
   Solver_pt->preconditioner_pt()=Prec_pt;
   
   
   // By default, the LSC Preconditioner uses SuperLU as
   // an exact preconditioner (i.e. a solver) for the
   // momentum and Schur complement blocks. 
   // Can overwrite this by passing pointers to 
   // other preconditioners that perform the (approximate)
   // solves of these blocks.
   
   // Create internal preconditioners used on Schur block
   //-----------------------------------------------------
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
   if (use_hypre_for_pressure)
    {
     P_matrix_preconditioner_pt = new HyprePreconditioner;
     
     // Set parameters for use as preconditioner on Poisson-type problem
     Hypre_default_settings::set_defaults_for_2D_poisson_problem(
      static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
     
     // Use Hypre for the Schur complement block
     Prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
     
     // Shut up!
     static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
      disable_doc_time();
    }
#endif    
#endif

   // Create internal preconditioners used on momentum block
   //--------------------------------------------------------
   if (use_block_diagonal_for_momentum)
    {
     F_matrix_preconditioner_pt = 
      new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
     if (use_hypre_for_momentum)
      {
       dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
        (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
        (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
      }
#endif
#endif
       // Use Hypre for momentum block 
       Prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
    }
   else
    {
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
     if (use_hypre_for_momentum)
      {
       F_matrix_preconditioner_pt = new HyprePreconditioner;
       
       // Shut up!
       static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt)->
        disable_doc_time();
       
       // Set parameters for use as preconditioner in for momentum 
       // block in Navier-Stokes problem
       Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
        static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt));
       
       // Use Hypre for momentum block 
       Prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
      }
#endif
#endif
    }

  }

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element. Optional command line
/// args specify multiplier for number of element rows/columns
/// and flag to indicate if iterative solver is used.
/// Multiplier and flag both default to 1. 
//=====================================================================
int main(int argc, char **argv)
{

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
    
 // Build the problem with "my" TaylorHood elements
 RectangularDrivenCavityProblem<MyTaylorHoodElement<2> > problem;
 
 // Solve the problem
 problem.newton_solve();
 
 // Outpt the solution
 problem.doc_solution(doc_info);
 
} // end_of_main











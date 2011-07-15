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
#include "spherical_navier_stokes_preconditioners.h"


namespace oomph
{



//===========================================================================
/// Setup the least-squares commutator Navier Stokes preconditioner. This
/// extracts blocks corresponding to the velocity and pressure unknowns,
/// creates the matrices actually needed in the application of the
/// preconditioner and deletes what can be deleted... Note that
/// this preconditioner needs a CRDoubleMatrix.
//============================================================================
 void SphericalNavierStokesLSCPreconditioner::
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

#ifdef PARANOID
  // paranoid check that the navier stokes mesh pt has been set
  if (Navier_stokes_mesh_pt == 0)
   {
    std::ostringstream error_message;
    error_message << "The navier stokes elements mesh pointer must be set.\n"
                  << "Use method set_navier_stokes_mesh(...)";
    throw OomphLibError(error_message.str(),
                     	"NavierStokesLSCPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // set the mesh
  this->set_mesh(0,problem_pt,Navier_stokes_mesh_pt);

  // Get blocks
  // ----------

  // In comes the current Jacobian. Recast it to a CR double matrix;
  // shout if that can't be done.
  CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt);

#ifdef PARANOID
  if (cr_matrix_pt==0)
   {
    std::ostringstream error_message;
    error_message << "SphericalNavierStokesLSCPreconditioner only works with "
                  << "CRDoubleMatrix matrices" << std::endl;
    throw OomphLibError(error_message.str(),
                     	"SphericalNavierStokesLSCPreconditioner::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Set up block look up schemes (done automatically in the
  // BlockPreconditioner base class, based on the information 
  // provided in the block-preconditionable elements in the problem)

  // this preconditioner has two types of block:
  // type 0: velocity - corresponding to DOFs 0 to n-2
  // type 1: pressure - corresponding to DOF n-1
  unsigned ndof_types = this->ndof_types_in_mesh(0);
  Vector<unsigned> dof_to_block_map(ndof_types);
  dof_to_block_map[ndof_types-1]=1;
  this->block_setup(problem_pt,matrix_pt,dof_to_block_map);

  // determine whether the F preconditioner is a block preconditioner (and
  // therefore a subsidiary preconditioner)
  BlockPreconditioner<CRDoubleMatrix>* F_block_preconditioner_pt = 
   dynamic_cast<BlockPreconditioner<CRDoubleMatrix>* >(F_preconditioner_pt);
  F_preconditioner_is_block_preconditioner = true;
  if (F_block_preconditioner_pt == 0)
   {
    F_preconditioner_is_block_preconditioner = false;
   }
  
  // Get the blocks
  this->get_block(0,1,cr_matrix_pt,Block_matrix_0_1_pt);
  CRDoubleMatrix* block_matrix_1_0_pt = 0;
  this->get_block(1,0,cr_matrix_pt,block_matrix_1_0_pt);
  CRDoubleMatrix* block_matrix_0_0_pt = 0;
  this->get_block(0,0,cr_matrix_pt,block_matrix_0_0_pt);

  // Initialise timers 
  double t_start = TimingHelpers::timer();

  // Generate pressure Poisson matrix
  // --------------------------------
  P_matrix_pt = new CRDoubleMatrix;

  // Vector to hold inverse of velocity mass matrix diagonal
  Vector<double> ivmm_diagonal;

  // If required,  multiply Block_matrix_pt(1,0) which currently
  // points to the divergence matrix, D, by the inverse diagonal mass matrix:
  CRDoubleMatrix* ivmm_pt = 0;
  if (P_matrix_using_scaling)
   {
     ivmm_pt = assemble_velocity_mass_matrix_diagonal();
    CRDoubleMatrix* temp_matrix_pt = new CRDoubleMatrix;
    block_matrix_1_0_pt->multiply(*ivmm_pt,*temp_matrix_pt);
    delete block_matrix_1_0_pt;
    block_matrix_1_0_pt = temp_matrix_pt;
   }

  // Now multiply whatever's currently stored in  Block_matrix_pt(1,0)
  // (i.e. D or D Q^{-1}) from the right by G, which is stored
  // in Block_matrix_pt(0,1). Store the result in P_matrix_pt
  // which now contains DG or D Q^{-1} G as required.
  block_matrix_1_0_pt->multiply(*Block_matrix_0_1_pt, *P_matrix_pt);

  // Generate matrix E for multiplication in Schur complement approximation
  // ----------------------------------------------------------------------

  // Auxiliary matrix for intermediate results
  CRDoubleMatrix* aux_matrix_pt = new CRDoubleMatrix;

  // Multiply the momentum matrix, stored in Block_matrix_pt(0,0)
  // by the content of Block_matrix_pt(1,0) (either D or
  // D Q^{-1}) and store the result, namely B Q^{-1} F in aux_matrix_pt.
  block_matrix_1_0_pt->multiply(*block_matrix_0_0_pt, *aux_matrix_pt);

  // We no longer need the divergence matrix -- kill it.
  delete block_matrix_1_0_pt;

  // Multiply auxiliary matrix by diagonal of velocity mass matrix if required
  if (P_matrix_using_scaling)
   {
    CRDoubleMatrix* temp_matrix_pt = new CRDoubleMatrix;
    aux_matrix_pt->multiply(*ivmm_pt,*temp_matrix_pt);
    delete aux_matrix_pt;
    aux_matrix_pt = temp_matrix_pt;
   }

  // Now multiply the auxiliary matrix by the gradient matrix,
  // stored in Block_matrix_pt(0,1), to obtain either
  // E = B Q^{-1} F Q^{-1} G or  E = B  F  G, as required.
  E_matrix.clear();
  aux_matrix_pt->multiply(*Block_matrix_0_1_pt,E_matrix);

  // Clean up memory
  delete aux_matrix_pt;
  delete ivmm_pt;

  double t_end = TimingHelpers::timer();
  double setup_time= t_end-t_start;
  t_start = t_end;
 
 // Output times
 if(Doc_time)
  {
   oomph_info << "Time to generate additional LSC matrices [sec]: "
              << setup_time << "\n";
  }
 
  // if the P preconditioner has not been setup
  if (P_preconditioner_pt == 0)
   {
    P_preconditioner_pt = new SuperLUPreconditioner;
    Using_default_p_preconditioner = true;
   }

  // Setup the preconditioner for the Pressure matrix
  P_preconditioner_pt->setup(problem_pt, P_matrix_pt);

  // Delete oomph-lib pressure Poisson matrix
  delete P_matrix_pt;
  P_matrix_pt=0;
  t_end = TimingHelpers::timer();
  setup_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  t_start = t_end;
 
  // Output times
  if(Doc_time)
   {
    oomph_info << "P matrix sub-preconditioner setup time [sec]: "
               << setup_time << "\n";
   }
  
  // Set up solver for solution of system with momentum matrix
  // ----------------------------------------------------------

  // if the F preconditioner has not been setup
  if (F_preconditioner_pt == 0)
   {
    F_preconditioner_pt = new SuperLUPreconditioner;
    Using_default_f_preconditioner = true;
   }

  // if F is a block preconditioner
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
    F_block_preconditioner_pt->setup(problem_pt,matrix_pt);
   }
  // otherwise F is not a block preconditioner
  else
   {
    F_preconditioner_pt->setup(problem_pt,block_matrix_0_0_pt);
   }

  // can now delete block 0 0
  delete block_matrix_0_0_pt;

  // compute the setup time
  t_end = TimingHelpers::timer();
  setup_time=double(t_end-t_start)/CLOCKS_PER_SEC;
 
  // Output times
  if(Doc_time)
   {
    oomph_info << "F matrix sub-preconditioner setup time [sec]: "
               << setup_time << "\n";
   }
  
  // Remember that the preconditioner has been setup so
  // the stored information can be wiped when we
  // come here next...
  Preconditioner_has_been_setup = true;
 }



//=======================================================================
 /// Apply preconditioner to r.
//=======================================================================
 void SphericalNavierStokesLSCPreconditioner::
 preconditioner_solve(const DoubleVector &r, DoubleVector &z)
 {
#ifdef PARANOID
  if (Preconditioner_has_been_setup==false)
   {
    std::ostringstream error_message;
    error_message << "setup must be called before using preconditioner_solve";
    throw OomphLibError(
     error_message.str(),
     "SphericalNavierStokesLSCPreconditioner::preconditioner_solve()",
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
       "SphericalNavierStokesLSCPreconditioner::preconditioner_solve()",
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
    throw OomphLibError(
     error_message.str(),
     "SphericalNavierStokesLSCPreconditioner::preconditioner_solve()",
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
  Block_matrix_0_1_pt->multiply(another_temp_vec, temp_vec);

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
    "SphericalNavierStokesLSCPreconditioner::preconditioner_solve()",
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


//========================================================================
/// Helper function to assemble the diagonal of the velocity
/// mass matrix from the elemental contributions defined in
/// SphericalNavierStokesEquations<DIM>::get_velocity_mass_matrix_diagonal(...)
//========================================================================
 CRDoubleMatrix* 
 SphericalNavierStokesLSCPreconditioner
 ::assemble_velocity_mass_matrix_diagonal()
 {

  // determine the rows required by this processor
  unsigned first_row = this->block_distribution_pt(0)->first_row();
  unsigned nrow_local = this->block_distribution_pt(0)->nrow_local();
  unsigned nrow = this->block_distribution_pt(0)->nrow();
  
  // create storage for the diagonals
  double* m_values = new double[nrow_local];
  for (unsigned i = 0; i < nrow_local; i++)
   {
    m_values[i] = 0;
   }

#ifdef OOMPH_HAS_MPI
  // store the problem pt
  const Problem* problem_pt = this->problem_pt();
#endif  

  // if the problem is distributed
  bool distributed = false;
#ifdef OOMPH_HAS_MPI
  if (problem_pt->distributed() || 
      this->master_distribution_pt()->distributed())
   {
     distributed = true;
   }
#endif

  // next we get the diagonal velocity mass matrix data
  if (distributed)
   {
#ifdef OOMPH_HAS_MPI
    // the number of processors
    unsigned nproc = this->problem_pt()->communicator_pt()->nproc();

    // and my rank
    unsigned my_rank = this->problem_pt()->communicator_pt()->my_rank();

    // determine the rows for which we have lookup rows

    // if the problem is NOT distributed then we only classify global equation
    // on this processor to avoid duplication (as every processor holds 
    // elvery element)
    unsigned first_lookup_row = 0; 
    unsigned last_lookup_row = 0;
    if (!problem_pt->distributed())
     {
      first_lookup_row = this->master_distribution_pt()->first_row();
      last_lookup_row = first_lookup_row + 
       this->master_distribution_pt()->nrow_local() - 1;
     }
    // attempt to classify every degree of freedom
    else
     {
      first_lookup_row = this->first_row();
      last_lookup_row = first_lookup_row + this->nrow_local() - 1;
     }

    // find number of local elements
    unsigned n_el = Navier_stokes_mesh_pt->nelement();
    
    // the diagonal velocity mass matrix contributions that have been
    // classified and should be sent to another processor
    Vector<double>* classified_contributions_send 
     = new Vector<double>[nproc];

    // the corresponding block indices
    Vector<unsigned>* classified_indices_send
     = new Vector<unsigned>[nproc];

    // the maitrix contributions that cannot be classified by this processor
    // and therefore must be sent to another for classification
    Vector<double>* unclassified_contributions_send
     = new Vector<double>[nproc];

    // the corresponding global indices that require classification
    Vector<unsigned>* unclassified_indices_send
     = new Vector<unsigned>[nproc];

    // get the master distribution pt
    const LinearAlgebraDistribution* master_distribution_pt = 
     this->master_distribution_pt();

    // get the velocity distribution pt
    const LinearAlgebraDistribution* velocity_dist_pt 
     = this->block_distribution_pt(0);

    // get the contribution for each element
    for (unsigned e = 0; e < n_el; e++)
     {

      // check that the element is not halo d
      if (!Navier_stokes_mesh_pt->element_pt(e)->is_halo())
       {
        
        // find number of degrees of freedom in the element
        // (this is slightly too big because it includes the
        // pressure dofs but this doesn't matter)
        unsigned el_dof = Navier_stokes_mesh_pt->element_pt(e)->ndof();
        
        // allocate local storage for the element's contribution to the
        // velocity mass matrix diagonal
        Vector<double> el_vmm_diagonal(el_dof);
        dynamic_cast< SphericalNavierStokesEquations* >
         ( Navier_stokes_mesh_pt->element_pt(e) )
         ->get_velocity_mass_matrix_diagonal(el_vmm_diagonal);

        // get the contribution for each dof
        for (unsigned i = 0; i < el_dof; i++)
         {

          //Get the equation number
          unsigned eqn_number = Navier_stokes_mesh_pt
           ->element_pt(e)->eqn_number(i);

          // if I have lookup information on this processor
          if (eqn_number >= first_lookup_row && 
              eqn_number <= last_lookup_row)
           {
            
            // bypass non velocity DOFs
            if ( this->block_number(eqn_number)==0 )
             {
              
              // get the index in the block
              unsigned index = this->index_in_block(eqn_number);
             
              // determine which processor requires the block index
              for (unsigned p = 0; p < nproc; p++)
               {
                if (index >= velocity_dist_pt->first_row(p) &&
                    (index < 
                     (velocity_dist_pt->first_row(p)
                      +velocity_dist_pt->nrow_local(p))))
                 {
                  
                  // if it is required by this processor then add the 
                  // contribution
                  if (p == my_rank)
                   {
                    m_values[index-first_row] += el_vmm_diagonal[i];
                   }
                  // other wise store it for communication
                  else
                   {
                    classified_contributions_send[p]
                     .push_back(el_vmm_diagonal[i]);
                    classified_indices_send[p].push_back(index);
                   }
                 }
               }
             }
           }
        
        // if we do not have the lookup information on this processor
        // then we send the mass matrix contribution to a processor
        // which we know does have the lookup information
        // the assumption: the processor for which the master block
        // preconditioner distribution will definitely hold the lookup
        // data for eqn_number (although others may)
        else if (problem_pt->distributed())
         {
          
          // determine which processor requires the block index
          unsigned p = 0;
          while (!(eqn_number >= master_distribution_pt->first_row(p) &&
                   (eqn_number < (master_distribution_pt->first_row(p)
                                  +master_distribution_pt->nrow_local(p)))))
            {
             p++;
            }

            // store the data
            unclassified_contributions_send[p]
             .push_back(el_vmm_diagonal[i]);
            unclassified_indices_send[p].push_back(eqn_number);
           }
         }
       }
     }
   
    // next the unclassified contributions are communicated to 
    // processors that can classify them
    
    // first determine how many unclassified rows are to be sent to
    // each processor
    unsigned* n_unclassified_send = new unsigned[nproc];
    for (unsigned p = 0; p < nproc; p++)
     {
      if (p == my_rank)
       {
        n_unclassified_send[p] = 0;
       }
      else
       {
        n_unclassified_send[p] 
         = unclassified_contributions_send[p].size();
       }
     }
    
    // then all-to-all com number of unclassified to be sent / recv
    unsigned* n_unclassified_recv = new unsigned[nproc];
    MPI_Alltoall(n_unclassified_send,1,MPI_UNSIGNED,
                 n_unclassified_recv,1,MPI_UNSIGNED,
                 this->problem_pt()->communicator_pt()->mpi_comm());
    
    // the base displacement for the sends
    MPI_Aint base_displacement;
    MPI_Address(m_values,&base_displacement);
    
    // allocate storage for the data to be recieved
    // and post the sends and recvs
      Vector<double*> unclassified_contributions_recv(nproc);
      Vector<unsigned*> unclassified_indices_recv(nproc);
      Vector<MPI_Request> unclassified_recv_requests;
      Vector<MPI_Request> unclassified_send_requests;
      Vector<unsigned> unclassified_recv_proc;
      for (unsigned p = 0; p < nproc; p++)
       {
        if (p != my_rank)
         {
          // recv
          if (n_unclassified_recv[p] > 0)
           {
            unclassified_contributions_recv[p] 
             = new double[n_unclassified_recv[p]];
            unclassified_indices_recv[p] = new 
             unsigned[n_unclassified_recv[p]];
            
            // data for the struct data type
            MPI_Datatype recv_types[2];
            MPI_Aint recv_displacements[2];
            int recv_sz[2];
            
            // contributions
            MPI_Type_contiguous(n_unclassified_recv[p],MPI_DOUBLE,
                                &recv_types[0]);
            MPI_Type_commit(&recv_types[0]);
            MPI_Address(unclassified_contributions_recv[p],
                        &recv_displacements[0]);
            recv_displacements[0] -= base_displacement;
            recv_sz[0] = 1;
            
            // indices
            MPI_Type_contiguous(n_unclassified_recv[p],MPI_UNSIGNED,
                                &recv_types[1]);
            MPI_Type_commit(&recv_types[1]);
            MPI_Address(unclassified_indices_recv[p],
                        &recv_displacements[1]);
            recv_displacements[1] -= base_displacement;
            recv_sz[1] = 1;
            
            // build the final recv type
            MPI_Datatype final_recv_type;
            MPI_Type_struct(2,recv_sz,recv_displacements,recv_types,
                            &final_recv_type);
            MPI_Type_commit(&final_recv_type);
            
            // and recv
            MPI_Request req;
            MPI_Irecv(m_values,1,final_recv_type,p,0,
                      problem_pt->communicator_pt()->mpi_comm(),&req);
            unclassified_recv_requests.push_back(req);
            unclassified_recv_proc.push_back(p); 
            MPI_Type_free(&recv_types[0]);
            MPI_Type_free(&recv_types[1]);
            MPI_Type_free(&final_recv_type);
           }
          
          // send
          if (n_unclassified_send[p] > 0)
           {
            // data for the struct data type
            MPI_Datatype send_types[2];
            MPI_Aint send_displacements[2];
            int send_sz[2];
            
            // contributions
            MPI_Type_contiguous(n_unclassified_send[p],MPI_DOUBLE,
                                &send_types[0]);
            MPI_Type_commit(&send_types[0]);
            MPI_Address(&unclassified_contributions_send[p][0],
                        &send_displacements[0]);
            send_displacements[0] -= base_displacement;
            send_sz[0] = 1;
            
            // indices
            MPI_Type_contiguous(n_unclassified_send[p],MPI_UNSIGNED,
                                &send_types[1]);
            MPI_Type_commit(&send_types[1]);
            MPI_Address(&unclassified_indices_send[p][0],
                        &send_displacements[1]);
            send_displacements[1] -= base_displacement;
            send_sz[1] = 1;
            
            // build the final send type
            MPI_Datatype final_send_type;
            MPI_Type_struct(2,send_sz,send_displacements,send_types,
                            &final_send_type);
            MPI_Type_commit(&final_send_type);
            
            // and send
            MPI_Request req;
            MPI_Isend(m_values,1,final_send_type,p,0,
                      problem_pt->communicator_pt()->mpi_comm(),&req);
            unclassified_send_requests.push_back(req);
            MPI_Type_free(&send_types[0]);
            MPI_Type_free(&send_types[1]);
            MPI_Type_free(&final_send_type);
           }
         }
       }
   
      // next classify the data as it is received
      unsigned n_unclassified_recv_req = unclassified_recv_requests.size();
     while (n_unclassified_recv_req > 0)
      {
       // get the processor number and remove the completed request
       // for the vector of requests
       int req_num;
       MPI_Waitany(n_unclassified_recv_req,&unclassified_recv_requests[0],
                   &req_num,MPI_STATUS_IGNORE);
       unsigned p = unclassified_recv_proc[req_num];
       unclassified_recv_requests.erase(unclassified_recv_requests.begin()
                                        +req_num);    
       unclassified_recv_proc.erase(unclassified_recv_proc.begin()+req_num);
       n_unclassified_recv_req--;
       
       // next classify the dofs 
       // and store them for sending to other processors if required
       unsigned n_recv = n_unclassified_recv[p];
       for (unsigned i = 0; i < n_recv; i++)
        {
         unsigned eqn_number = unclassified_indices_recv[p][i];
         // bypass non velocity DOFs
         if ( this->block_number(eqn_number)==0 )
          {
           
           // get the index in the block
           unsigned index = this->index_in_block(eqn_number);
           
           // determine which processor requires the block index
           for (unsigned pp = 0; pp < nproc; pp++)
            {
             
             
             if (index >= velocity_dist_pt->first_row(pp) &&
                 (index < 
                  (velocity_dist_pt->first_row(pp)
                   +velocity_dist_pt->nrow_local(pp))))
              {
               
               // if it is required by this processor then add the 
               // contribution
               if (pp == my_rank)
                {
                 m_values[index-first_row]
                  += unclassified_contributions_recv[p][i];
                }
               // other wise store it for communication
               else
                {
                 double v = unclassified_contributions_recv[p][i];
                 classified_contributions_send[pp].push_back(v);
                 classified_indices_send[pp].push_back(index);
                }
              }
            }
          }
        }
       
       // clean up
       delete[] unclassified_contributions_recv[p];
       delete[] unclassified_indices_recv[p];
      }
     
     // now all indices have been classified
     
     // next the classified contributions are communicated to 
     // processors that require them
     
     // first determine how many classified rows are to be sent to
     // each processor
     unsigned* n_classified_send = new unsigned[nproc];
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p == my_rank)
        {
         n_classified_send[p] = 0;
        }
       else
        {
         n_classified_send[p] 
          = classified_contributions_send[p].size();
        }
      }
     
     // then all-to-all com number of classified to be sent / recv
     unsigned* n_classified_recv = new unsigned[nproc];
     MPI_Alltoall(n_classified_send,1,MPI_UNSIGNED,
                  n_classified_recv,1,MPI_UNSIGNED,
                  this->problem_pt()->communicator_pt()->mpi_comm());
     
     // allocate storage for the data to be recieved
     // and post the sends and recvs
     Vector<double*> classified_contributions_recv(nproc);
     Vector<unsigned*> classified_indices_recv(nproc);
     Vector<MPI_Request> classified_recv_requests;
     Vector<MPI_Request> classified_send_requests;
     Vector<unsigned> classified_recv_proc;
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p != my_rank)
        {
         // recv
         if (n_classified_recv[p] > 0)
          {
           classified_contributions_recv[p] 
            = new double[n_classified_recv[p]];
           classified_indices_recv[p] = new unsigned[n_classified_recv[p]];
           
           // data for the struct data type
           MPI_Datatype recv_types[2];
           MPI_Aint recv_displacements[2];
           int recv_sz[2];
           
           // contributions
           MPI_Type_contiguous(n_classified_recv[p],MPI_DOUBLE,
                               &recv_types[0]);
           MPI_Type_commit(&recv_types[0]);
           MPI_Address(classified_contributions_recv[p],
                       &recv_displacements[0]);
           recv_displacements[0] -= base_displacement;
           recv_sz[0] = 1;
           
           // indices
           MPI_Type_contiguous(n_classified_recv[p],MPI_UNSIGNED,
                               &recv_types[1]);
           MPI_Type_commit(&recv_types[1]);
           MPI_Address(classified_indices_recv[p],
                       &recv_displacements[1]);
           recv_displacements[1] -= base_displacement;
           recv_sz[1] = 1;
           
           // build the final recv type
           MPI_Datatype final_recv_type;
           MPI_Type_struct(2,recv_sz,recv_displacements,recv_types,
                           &final_recv_type);
           MPI_Type_commit(&final_recv_type);
           
           // and recv
           MPI_Request req;
           MPI_Irecv(m_values,1,final_recv_type,p,0,
                     problem_pt->communicator_pt()->mpi_comm(),&req);
           classified_recv_requests.push_back(req);
           classified_recv_proc.push_back(p);
           MPI_Type_free(&recv_types[0]);
           MPI_Type_free(&recv_types[1]);
           MPI_Type_free(&final_recv_type);
          }
         
         // send
         if (n_classified_send[p] > 0)
          {
           // data for the struct data type
           MPI_Datatype send_types[2];
           MPI_Aint send_displacements[2];
           int send_sz[2];
           
           // contributions
           MPI_Type_contiguous(n_classified_send[p],MPI_DOUBLE,
                               &send_types[0]);
           MPI_Type_commit(&send_types[0]);
           MPI_Address(&classified_contributions_send[p][0],
                       &send_displacements[0]);
           send_displacements[0] -= base_displacement;
           send_sz[0] = 1;
           
           // indices
           MPI_Type_contiguous(n_classified_send[p],MPI_UNSIGNED,
                               &send_types[1]);
           MPI_Type_commit(&send_types[1]);
           MPI_Address(&classified_indices_send[p][0],
                       &send_displacements[1]);
           send_displacements[1] -= base_displacement;
           send_sz[1] = 1;
           
           // build the final send type
           MPI_Datatype final_send_type;
           MPI_Type_struct(2,send_sz,send_displacements,send_types,
                           &final_send_type);
           MPI_Type_commit(&final_send_type);
           
           // and send
           MPI_Request req;
           MPI_Isend(m_values,1,final_send_type,p,0,
                     problem_pt->communicator_pt()->mpi_comm(),&req);
           classified_send_requests.push_back(req);
           MPI_Type_free(&send_types[0]);
           MPI_Type_free(&send_types[1]);
           MPI_Type_free(&final_send_type);
          }
        }
      }
     
     // next classify the data as it is received
     unsigned n_classified_recv_req = classified_recv_requests.size();
     while (n_classified_recv_req > 0)
      {
       // get the processor number and remove the completed request
       // for the vector of requests
       int req_num;
       MPI_Waitany(n_classified_recv_req,&classified_recv_requests[0],
                   &req_num,MPI_STATUS_IGNORE);
       unsigned p = classified_recv_proc[req_num];
       classified_recv_requests.erase(classified_recv_requests.begin()
                                      +req_num);    
       classified_recv_proc.erase(classified_recv_proc.begin()+req_num);
       n_classified_recv_req--;
       
       // next classify the dofs 
       // and store them for sending to other processors if required
       unsigned n_recv = n_classified_recv[p];
       for (unsigned i = 0; i < n_recv; i++)
        {
         m_values[classified_indices_recv[p][i]-first_row]
          += classified_contributions_recv[p][i];
        }
       
       // clean up
       delete[] classified_contributions_recv[p];
       delete[] classified_indices_recv[p];
      }
     
     // wait for the unclassified sends to complete
     unsigned n_unclassified_send_req = unclassified_send_requests.size();
     MPI_Waitall(n_unclassified_send_req,&classified_send_requests[0],
                 MPI_STATUS_IGNORE);
     delete[] unclassified_contributions_send;
     delete[] unclassified_indices_send;
     
     // wait for the classified sends to complete
     unsigned n_classified_send_req = classified_send_requests.size();
     MPI_Waitall(n_classified_send_req,&classified_send_requests[0],
                 MPI_STATUS_IGNORE);
     delete[] classified_indices_send;
     delete[] classified_contributions_send;
#endif
   }

  // or if the problem is not distributed
  else
   {

    // find number of elements
    unsigned n_el = Navier_stokes_mesh_pt->nelement();
    
    // get the contribution for each element
    for (unsigned e = 0; e < n_el; e++)
     {
      
      // find number of degrees of freedom in the element
      // (this is slightly too big because it includes the
      // pressure dofs but this doesn't matter)
      unsigned el_dof = Navier_stokes_mesh_pt->element_pt(e)->ndof();

      // allocate local storage for the element's contribution to the
      // velocity mass matrix diagonal
      Vector<double> el_vmm_diagonal(el_dof);
      dynamic_cast< SphericalNavierStokesEquations* >
       ( Navier_stokes_mesh_pt->element_pt(e) )
       ->get_velocity_mass_matrix_diagonal(el_vmm_diagonal);

      // get the contribution for each dof
      for (unsigned i = 0; i < el_dof; i++)
       {
        //Get the equation number
        unsigned eqn_number = Navier_stokes_mesh_pt
         ->element_pt(e)->eqn_number(i);
        
        // bypass non velocity DOFs
        if ( this->block_number(eqn_number)==0 )
         {

          // get the index in the block
          unsigned index = this->index_in_block(eqn_number);
          
          // if it is required on this processor
          if (index >= first_row &&
              index < first_row + nrow_local)
           {
            m_values[index-first_row] += el_vmm_diagonal[i];
           }
         }
       }
     } 
   }

  // create column index and row start
  int* m_column_index = new int[nrow_local];
  int* m_row_start = new int[nrow_local+1];
  for (unsigned i = 0; i < nrow_local; i++)
   {
    m_values[i] = 1 / m_values[i];
    m_column_index[i] = first_row + i;
    m_row_start[i] = i;
   }
  m_row_start[nrow_local] = nrow_local;
  
  // build the matrix
  CRDoubleMatrix* m_pt = new CRDoubleMatrix(this->block_distribution_pt(0));
  m_pt->build_without_copy(nrow,nrow_local,m_values,m_column_index,
                           m_row_start);
  
  // return the matrix;
  return m_pt;   
 }

//=========================================================================
/// Helper function to delete preconditioner data.
//=========================================================================
 void SphericalNavierStokesLSCPreconditioner::clean_up_memory()
 {
  if (Preconditioner_has_been_setup)
   {
    // delete blocks
    delete Block_matrix_0_1_pt;
    Block_matrix_0_1_pt = 0;  

    // delete stuff from velocity solve
    if (Using_default_f_preconditioner)
     {
      delete F_preconditioner_pt;
      F_preconditioner_pt = 0;
     }

    // delete stuff from Schur complement approx
    delete P_matrix_pt;
    P_matrix_pt = 0;
    if (Using_default_p_preconditioner)
     {
      delete P_preconditioner_pt;
      P_preconditioner_pt = 0;
     }
   }
 }
}

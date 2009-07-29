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
//Non-templated multi-domain functions which act on more than one mesh 
//and set up the storage and interaction between the two

//oomph-lib header
#include "multi_domain.h"
#include "multi_domain.template.cc"
#include "mesh.h"
#include "algebraic_elements.h"
#include "macro_element_node_update_element.h"
#include "Qelements.h"

namespace oomph
{

//======================================================================
// Namespace for "global" multi-domain functions
//======================================================================
namespace Multi_domain_functions
 {
  // Lookup scheme for whether an element's integration point
  // has had an external element assigned to it
  Vector<Vector<unsigned> > External_element_located;

  // List of all Vectors used by the external storage routines 
  
  /// \short Vector of (local) coordinates at integration points of
  /// elements on current processor
  Vector<double> Local_zetas;

  /// \short Vector of (local) coordinates at integration points of
  /// elements received from elsewhere
  Vector<double> Zetas;

  /// \short Vector of the dimension of the element on current processor
  Vector<unsigned> Local_zeta_dim;

  /// \short Vector of the dimension of the element received from elsewhere
  Vector<unsigned> Zeta_dim;

  /// \short Vector to indicate locally which processor a coordinate 
  /// has been located on
  Vector<int> Found_zeta;

  /// \short Vector of local coordinates within any elements found
  /// locally by current processor
  Vector<Vector<double> > Found_ss;

  /// \short Vector to indicate (on another process) whether a
  /// located element should be newly created (2), already exists (1), or
  /// is not on the current process at all (0)
  Vector<unsigned> Located_element;

  /// \short Vector of the local coordinates for each entry in Located_element
  Vector<double> Located_coord;

  /// \short Vector for current processor which indicates when an external
  /// halo element (and subsequent nodes) should be created
  Vector<unsigned> Located_zetas;

  /// \short Vector of doubles to be sent from another processor
  Vector<double> Double_values;

  /// \short Vector of unsigneds to be sent from another processor
  Vector<unsigned> Unsigned_values;

  // Counters for each of the above arrays

  /// \short Double_values
  unsigned Count_double_values;

  /// \short Unsigned_values
  unsigned Count_unsigned_values;

  /// \short Located_coord
  unsigned Count_located_coord;

  /// \short Local_zeta_dim
  unsigned Count_local_zeta_dim;

  /// \short Zeta_dim
  unsigned Count_zeta_dim;

  /// \short Local_zetas
  unsigned Count_local_zetas;
 
  /// \short Zetas
  unsigned Count_zetas;

  /// Default parameters for the binning method

  /// \short Bool to tell the MeshAsGeomObject whether to calculate
  /// the extreme coordinates of the bin structure
  bool Compute_extreme_bin_coordinates=true;

  /// \short Number of bins in the first dimension in binning method in
  /// setup_multi_domain_interaction(). Default value of 10.
  unsigned Nx_bin=10;

  /// \short Number of bins in the second dimension in binning method in
  /// setup_multi_domain_interaction(). Default value of 10.
  unsigned Ny_bin=10;

  /// \short Number of bins in the third dimension in binning method in
  /// setup_multi_domain_interaction(). Default value of 10.
  unsigned Nz_bin=10;

  /// \short (Measure of) the number of sampling points within the elements 
  /// when populating the bin
  unsigned Nsample_points=5;

  /// \short Minimum and maximum coordinates for
  /// each dimension of the external mesh used to "create" the bins in
  /// setup_multi_domain_interaction(). These can be set by user if they
  /// want to (otherwise the MeshAsGeomObject calculates these values based
  /// upon the mesh itself; see MeshAsGeomObject::get_max_and_min_coords(...))
  /// They default to "incorrect" values initially

  /// \short Minimum coordinate in first dimension
  double X_min=1.0;

  /// \short Maximum coordinate in first dimension
  double X_max=-1.0;

  /// \short Minimum coordinate in second dimension
  double Y_min=1.0;

  /// \short Maximum coordinate in second dimension
  double Y_max=-1.0;

  /// \short Minimum coordinate in third dimension
  double Z_min=1.0;

  /// \short Maximum coordinate in third dimension
  double Z_max=-1.0;

  /// \short Percentage offset to add to each extreme of the bin structure.
  /// Default value of 0.05.
  double Percentage_offset=0.05;

  /// \short Boolean to indicate when to use the bulk element as the
  /// external element.  Defaults to false, you must have set up FaceElements
  /// properly first in order for it to work
  bool Use_bulk_element_as_external=false;

  /// \short Boolean to indicate whether to doc timings or not.
  bool Doc_timings=false;

  /// \short Boolean to indicate whether to output basic info during
  ///        setup_multi_domain_interaction() routines
  bool Doc_stats=false;

  /// \short Boolean to indicate whether to output further info during
  ///        setup_multi_domain_interaction() routines
  bool Doc_full_stats=false;

#ifdef OOMPH_HAS_MPI

  /// \short Boolean to indicate when to check for duplicate data
  ///        between the external halo storage schemes
  bool Check_for_duplicates=true;

  // Functions for location method in multi-domain problems 

  //======================================================================
  /// Function which removes duplicate data that exist because
  /// they have been distinctly created by communications from different
  /// processors, whereas the data are in fact the same
  //======================================================================
  void remove_duplicate_data(Problem* problem_pt, Mesh* const &mesh_pt)
  {
   // Each individual container of external halo nodes has unique
   // nodes/equation numbers, but there may be some duplication between
   // two or more different containers; the following code checks for this
   // and removes the duplication by overwriting any data point with an already
   // existing eqn number with the original data point which had the eqn no.

   // Storage for existing global equation numbers for each node
   Vector<std::pair<Vector<int>,Node*> > existing_global_eqn_numbers;

   // Add all the global eqn numbers for external elements and nodes
   // that are stored locally first

   // Loop over external elements
   unsigned n_element=mesh_pt->nexternal_element();
   for (unsigned e_ext=0;e_ext<n_element;e_ext++)
    {
     FiniteElement* ext_el_pt=mesh_pt->external_element_pt(e_ext);
     // Loop over nodes
     unsigned n_node=ext_el_pt->nnode();
     for (unsigned j=0;j<n_node;j++)
      {
       Node* nod_pt=ext_el_pt->node_pt(j);
       unsigned n_val=nod_pt->nvalue();
       Vector<int> nodal_eqn_numbers(n_val);
       for (unsigned i_val=0;i_val<n_val;i_val++)
        {
         nodal_eqn_numbers[i_val]=nod_pt->eqn_number(i_val);
        }

       // All these nodes have unique equation numbers, so add all
       existing_global_eqn_numbers.push_back
        (make_pair(nodal_eqn_numbers,nod_pt));

       // If it's a SolidNode then add its extra equation numbers
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
       if (solid_nod_pt!=0)
        {
         unsigned n_val_solid=solid_nod_pt->variable_position_pt()->nvalue();
         Vector<int> solid_nodal_eqn_numbers(n_val_solid);
         for (unsigned i_val=0;i_val<n_val_solid;i_val++)
          {
           solid_nodal_eqn_numbers[i_val]=solid_nod_pt->
            variable_position_pt()->eqn_number(i_val);
          }
         // Add these equation numbers to the existing storage
         existing_global_eqn_numbers.push_back
          (make_pair(solid_nodal_eqn_numbers,solid_nod_pt));
        }

       // Take into account master nodes too
       if (dynamic_cast<RefineableElement*>(ext_el_pt)!=0)
        {
         int n_cont_int_values=dynamic_cast<RefineableElement*>
          (ext_el_pt)->ncont_interpolated_values();
         for (int i_cont=-1;i_cont<n_cont_int_values;i_cont++)
          {
           if (nod_pt->is_hanging(i_cont))
            {
             HangInfo* hang_pt=nod_pt->hanging_pt(i_cont);
             unsigned n_master=hang_pt->nmaster();
             for (unsigned m=0;m<n_master;m++)
              {
               Node* master_nod_pt=hang_pt->master_node_pt(m);
               unsigned n_val=master_nod_pt->nvalue();
               Vector<int> master_nodal_eqn_numbers(n_val);
               for (unsigned i_val=0;i_val<n_val;i_val++)
                {
                 master_nodal_eqn_numbers[i_val]=
                  master_nod_pt->eqn_number(i_val);
                }

               // Add these equation numbers to the existing storage
               existing_global_eqn_numbers.push_back
                (make_pair(master_nodal_eqn_numbers,master_nod_pt));

               // If this master is a SolidNode then add its extra eqn numbers
               SolidNode* master_solid_nod_pt=dynamic_cast<SolidNode*>
                (master_nod_pt);
               if (master_solid_nod_pt!=0)
                {
                 unsigned n_val_mst_solid=master_solid_nod_pt->
                  variable_position_pt()->nvalue();
                 Vector<int> master_solid_nodal_eqn_numbers(n_val_mst_solid);
                 for (unsigned i_val=0;i_val<n_val_mst_solid;i_val++)
                  {
                   master_solid_nodal_eqn_numbers[i_val]=master_solid_nod_pt->
                    variable_position_pt()->eqn_number(i_val);
                  }
                 // Add these equation numbers to the existing storage
                 existing_global_eqn_numbers.push_back
                  (make_pair(master_solid_nodal_eqn_numbers,
                             master_solid_nod_pt));
                }
              }
            }
          }
        }
      }
     // Internal data equation numbers do not need to be added since halo
     // elements cannot be external elements
    }

   // Now loop over the other processors from highest to lowest
   // (i.e. if there is a duplicate between these containers
   //  then this will use the node on the highest numbered processor)
   int n_proc=problem_pt->communicator_pt()->nproc();
   int my_rank=problem_pt->communicator_pt()->my_rank();
   for (int iproc=n_proc-1;iproc>=0;iproc--)
    {
     // Don't have external halo elements with yourself!
     if (iproc!=my_rank)
      {
       // Loop over external halo elements with iproc for internal data
       // to remove the duplicates in the external halo element storage
       unsigned n_element=mesh_pt->nexternal_halo_element(iproc);
       for (unsigned e_ext=0;e_ext<n_element;e_ext++)
        {
         FiniteElement* ext_el_pt=mesh_pt->
          external_halo_element_pt(iproc,e_ext);

         // Loop over nodes
         unsigned n_node=ext_el_pt->nnode();
         for (unsigned j=0;j<n_node;j++)
          {
           Node* nod_pt=ext_el_pt->node_pt(j);
           unsigned n_val=nod_pt->nvalue();
           Vector<int> nodal_eqn_numbers(n_val);
           for (unsigned i_val=0;i_val<n_val;i_val++)
            {
             nodal_eqn_numbers[i_val]=nod_pt->eqn_number(i_val);

            }
           // Check for duplicates with the existing set of global eqn numbers

           // Test to see if there is a duplicate with current
           unsigned n_exist_eqn=existing_global_eqn_numbers.size();
           bool is_a_duplicate=false;
           // Loop over the existing global equation numbers
           for (unsigned i=0;i<n_exist_eqn;i++)
            {
             // If the number of values is different from the size of the 
             // vector then they are already different, so only test if 
             // the sizes are the same
             if (n_val==(existing_global_eqn_numbers[i].first).size())
              {
               // Loop over values
               for (unsigned i_val=0;i_val<n_val;i_val++)
                {
                 // Make sure it isn't a pinned dof already
                 if (nodal_eqn_numbers[i_val]>=0)
                  {
                   // Test it against the equivalent entry in existing...
                   if (nodal_eqn_numbers[i_val]==
                       (existing_global_eqn_numbers[i].first)[i_val])
                    {
                     is_a_duplicate=true;
                     // It's a duplicate, so point the current node at the
                     // other position instead!
                     ext_el_pt->node_pt(j)=
                      existing_global_eqn_numbers[i].second;
                     break;
                    }
                  }
                }
              }
             // Break out of the loop if we have already found a duplicate
             if (is_a_duplicate)
              {
               break;
              }
            }

           // If it's not a duplicate, add it to the existing storage
           if (!is_a_duplicate)
            {
             existing_global_eqn_numbers.push_back
              (make_pair(nodal_eqn_numbers,nod_pt));
             // These external halo nodes need to be unclassified in order for
             // them to be bypassed in assign_(local)_eqn_numbers; they
             // receive the correct global equation numbers during the
             // synchronisation process in 
             // Problem::copy_external_haloed_eqn_numbers_helper(...)   
             for (unsigned i_val=0;i_val<n_val;i_val++)
              {
               nod_pt->eqn_number(i_val)=Data::Is_unclassified;
              }
            }

           // Do the same for any SolidNodes
           SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
           if (solid_nod_pt!=0)
            {
             unsigned n_val_solid=solid_nod_pt->
              variable_position_pt()->nvalue();
             Vector<int> solid_nodal_eqn_numbers(n_val_solid);
             for (unsigned i_val=0;i_val<n_val_solid;i_val++)
              {
               solid_nodal_eqn_numbers[i_val]=solid_nod_pt->
                variable_position_pt()->eqn_number(i_val);
              }

             // Check for duplicate equation numbers

             // Test to see if there is a duplicate with current
             unsigned n_exist_eqn=existing_global_eqn_numbers.size();
             bool is_a_duplicate=false;
             // Loop over the existing global equation numbers
             for (unsigned i=0;i<n_exist_eqn;i++)
              {
               // If the number of values is different from the size of the 
               // vector then they are already different, so only test if 
               // the sizes are the same
               if (n_val==(existing_global_eqn_numbers[i].first).size())
                {
                 // Loop over values
                 for (unsigned i_val=0;i_val<n_val;i_val++)
                  {
                   // Make sure it isn't a pinned dof already
                   if (solid_nodal_eqn_numbers[i_val]>=0)
                    {
                     // Test it against the equivalent entry in existing...
                     if (solid_nodal_eqn_numbers[i_val]==
                         (existing_global_eqn_numbers[i].first)[i_val])
                      {
                       is_a_duplicate=true;
                       // It's a duplicate, so point the current node at the
                       // other position instead!
                       solid_nod_pt=dynamic_cast<SolidNode*>
                        (existing_global_eqn_numbers[i].second);
                       break;
                      }
                    }
                  }
                }
               // Break out of the loop if we have already found a duplicate
               if (is_a_duplicate)
                {
                 break;
                }
              }

             // If it's not a duplicate, add it to the existing storage
             if (!is_a_duplicate)
              {
               existing_global_eqn_numbers.push_back
                (make_pair(solid_nodal_eqn_numbers,nod_pt));
               // These external halo nodes need to be unclassified in order
               // for them to be bypassed in assign_(local)_eqn_numbers; they
               // receive the correct global equation numbers during the
               // synchronisation process in 
               // Problem::copy_external_haloed_eqn_numbers_helper(...)   
               for (unsigned i_val=0;i_val<n_val;i_val++)
                {
                 solid_nod_pt->variable_position_pt()->
                  eqn_number(i_val)=Data::Is_unclassified;
                }
              }

            }

           // Do the same for any master nodes
           if (dynamic_cast<RefineableElement*>(ext_el_pt)!=0)
            {
             int n_cont_inter_values=dynamic_cast<RefineableElement*>
              (ext_el_pt)->ncont_interpolated_values();
             for (int i_cont=-1;i_cont<n_cont_inter_values;i_cont++)
              {
               if (nod_pt->is_hanging(i_cont))
                {
                 HangInfo* hang_pt=nod_pt->hanging_pt(i_cont);
                 unsigned n_master=hang_pt->nmaster();
                 for (unsigned m=0;m<n_master;m++)
                  {
                   Node* master_nod_pt=hang_pt->master_node_pt(m);
                   unsigned n_val=master_nod_pt->nvalue();
                   Vector<int> master_nodal_eqn_numbers(n_val);
                   for (unsigned i_val=0;i_val<n_val;i_val++)
                    {
                     master_nodal_eqn_numbers[i_val]=
                      master_nod_pt->eqn_number(i_val);
                    }

                   // Check for duplicate master eqn numbers

                   // Test to see if there is a duplicate with current
                   unsigned n_exist_eqn=existing_global_eqn_numbers.size();
                   bool is_a_duplicate=false;
                   // Loop over the existing global equation numbers
                   for (unsigned i=0;i<n_exist_eqn;i++)
                    {
                     // If the number of values is different from the size
                     // of the vector then they are already different,
                     // only test if the sizes are the same
                     if (n_val==
                         (existing_global_eqn_numbers[i].first).size())
                      {
                       // Loop over values
                       for (unsigned i_val=0;i_val<n_val;i_val++)
                        {
                         // Make sure it isn't a pinned dof already
                         if (master_nodal_eqn_numbers[i_val]>=0)
                          {
                           // Test it against the equivalent existing entry
                           if (master_nodal_eqn_numbers[i_val]==
                               (existing_global_eqn_numbers[i].first)
                               [i_val])
                            {
                             is_a_duplicate=true;
                             // It's a duplicate, so point this node's
                             // master at the original node instead!
                             // Need the weight of the original node
                             double m_weight=hang_pt->master_weight(m);
                             // Set master
                             ext_el_pt->node_pt(j)->hanging_pt(i_cont)
                              ->set_master_node_pt
                              (m,existing_global_eqn_numbers[i].second,
                               m_weight);
                             break;
                            }
                          }
                        }
                      }
                     // Break out of the loop if we have found a duplicate
                     if (is_a_duplicate)
                      {
                       break;
                      }
                    }

                   // If it's not a duplicate, add it to the existing storage
                   if (!is_a_duplicate)
                    {
                     existing_global_eqn_numbers.push_back
                      (make_pair(master_nodal_eqn_numbers,master_nod_pt));
                     // These external halo nodes need to be unclassified in
                     // order for them to be bypassed in assign_eqn_numbers;
                     // they receive the correct global equation numbers during
                     // the synchronisation process in
                     // Problem::copy_external_haloed_eqn_numbers_helper(...)
                     for (unsigned i_val=0;i_val<n_val;i_val++)
                      {
                       master_nod_pt->eqn_number(i_val)=Data::Is_unclassified;
                      }
                    }
           
                   // Do the same again for SolidNodes
                   SolidNode* solid_master_nod_pt=dynamic_cast<SolidNode*>
                    (master_nod_pt);
                   if (solid_master_nod_pt!=0)
                    {
                     unsigned n_val_solid=solid_master_nod_pt->
                      variable_position_pt()->nvalue();
                     Vector<int> sld_master_nodal_eqn_numbers(n_val_solid);
                     for (unsigned i_val=0;i_val<n_val_solid;i_val++)
                      {
                       sld_master_nodal_eqn_numbers[i_val]=solid_master_nod_pt
                        ->variable_position_pt()->eqn_number(i_val);
                      }

                     // Check for duplicate master eqn numbers

                     // Test to see if there is a duplicate with current
                     unsigned n_exist_eqn=existing_global_eqn_numbers.size();
                     bool is_a_duplicate=false;
                     // Loop over the existing global equation numbers
                     for (unsigned i=0;i<n_exist_eqn;i++)
                      {
                       // If the number of values is different from the size
                       // of the vector then they are already different,
                       // only test if the sizes are the same
                       if (n_val==
                           (existing_global_eqn_numbers[i].first).size())
                        {
                         // Loop over values
                         for (unsigned i_val=0;i_val<n_val;i_val++)
                          {
                           // Make sure it isn't a pinned dof already
                           if (sld_master_nodal_eqn_numbers[i_val]>=0)
                            {
                             // Test it against the equivalent existing entry
                             if (sld_master_nodal_eqn_numbers[i_val]==
                                 (existing_global_eqn_numbers[i].first)
                                 [i_val])
                              {
                               is_a_duplicate=true;
                               // It's a duplicate, so point the current master
                               // node at the original node instead!
                               // Need the weight of the original node
                               double m_weight=hang_pt->master_weight(m);
                               // Set master
                               ext_el_pt->node_pt(j)->hanging_pt(i_cont)
                                ->set_master_node_pt
                                (m,existing_global_eqn_numbers[i].second,
                                 m_weight);
                               break;
                              }
                            }
                          }
                        }
                       // Break out of the loop if we have found a duplicate
                       if (is_a_duplicate)
                        {
                         break;
                        }
                      }

                     // If it's not a duplicate then add it to existing storage
                     if (!is_a_duplicate)
                      {
                       existing_global_eqn_numbers.push_back
                        (make_pair(sld_master_nodal_eqn_numbers,
                                   master_nod_pt));
                       // These external halo nodes need to be unclassified in
                       // order for them to be bypassed in assign__eqn_numbers;
                       // they receive the correct global equation numbers
                       // during the synchronisation process in 
                       // Problem::copy_external_haloed_eqn_numbers_helper(...)
                       for (unsigned i_val=0;i_val<n_val;i_val++)
                        {
                         solid_master_nod_pt->variable_position_pt()->
                          eqn_number(i_val)=Data::Is_unclassified;
                        }
                      }


                    }

                  }
                }
              } // end hanging loop over continous interpolated variables
            }
          
          } // end loop over nodes on external halo elements
        
         // Reset any internal data to Is_unclassified as there
         // cannot be any duplication (a halo element cannot be an 
         // external halo element)
         unsigned n_int_data=ext_el_pt->ninternal_data();
         for (unsigned i=0;i<n_int_data;i++)
          {
           unsigned n_val=ext_el_pt->internal_data_pt(i)->nvalue();
           for (unsigned i_val=0;i_val<n_val;i_val++)
            {
             ext_el_pt->internal_data_pt(i)->eqn_number(i_val)=
              Data::Is_unclassified;
            }
          }
        } // end loop over external halo elements
      }
    } // end loop over processors

  }

//========================================================================
/// Send the zeta coordinates from the current process to 
/// the next process; receive from the previous process
//========================================================================
  void send_and_receive_missing_zetas(Problem* problem_pt)
  {
   // MPI info
   MPI_Status status;
   MPI_Request request;

   // Storage for number of processors, current process and communicator
   int n_proc=problem_pt->communicator_pt()->nproc();
   int my_rank=problem_pt->communicator_pt()->my_rank();
   OomphCommunicator* comm_pt=problem_pt->communicator_pt();

   // Work out processors to send and receive from
   int send_to_proc=my_rank+1;
   int recv_from_proc=my_rank-1;
   if (send_to_proc==n_proc) { send_to_proc=0; }
   if (recv_from_proc<0) { recv_from_proc=n_proc-1; }

   // Copy the "local" arrays in order to send to required process
   Count_zeta_dim=Count_local_zeta_dim;
   Zeta_dim.resize(Count_zeta_dim);
   for (unsigned i=0;i<Count_zeta_dim;i++)
    {
     Zeta_dim[i]=Local_zeta_dim[i];
    }
   Count_zetas=Count_local_zetas;
   Zetas.resize(Count_zetas);
   for (unsigned i=0;i<Count_zetas;i++)
    {
     Zetas[i]=Local_zetas[i];
    }

   // Send "local" values to next processor up
   MPI_Isend(&Count_local_zeta_dim,1,MPI_INT,
            send_to_proc,1,comm_pt->mpi_comm(),&request);
   // Receive from previous processor
   MPI_Recv(&Count_zeta_dim,1,MPI_INT,recv_from_proc,
            1,comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   // Now do the Zeta_dim array
   if (Count_local_zeta_dim!=0)
    {
     MPI_Isend(&Local_zeta_dim[0],Count_local_zeta_dim,MPI_INT,
              send_to_proc,3,comm_pt->mpi_comm(),&request);
    }
   if (Count_zeta_dim!=0)
    {
     Zeta_dim.resize(Count_zeta_dim);
     MPI_Recv(&Zeta_dim[0],Count_zeta_dim,MPI_INT,recv_from_proc,3,
              comm_pt->mpi_comm(),&status);
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Now do the Zetas array
   MPI_Isend(&Count_local_zetas,1,MPI_INT,send_to_proc,4,comm_pt->mpi_comm(),
             &request);
   MPI_Recv(&Count_zetas,1,MPI_INT,recv_from_proc,
            4,comm_pt->mpi_comm(),&status);

   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (Count_local_zetas!=0)
    {           
     MPI_Isend(&Local_zetas[0],Count_local_zetas,MPI_DOUBLE,
              send_to_proc,5,comm_pt->mpi_comm(),&request);
    }
   if (Count_zetas!=0)
    {       
     Zetas.resize(Count_zetas);
     MPI_Recv(&Zetas[0],Count_zetas,MPI_DOUBLE,recv_from_proc,5,
              comm_pt->mpi_comm(),&status);
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Now we should have the Zetas and Zeta_dim arrays set up correctly
   // for the next round of locations

  }

//========start of send_and_receive_located_info==========================
/// Send location information from current process; Received location
/// information from (current process + iproc) modulo (nproc)
//========================================================================
  void send_and_receive_located_info
  (int& iproc, Mesh* const &external_mesh_pt, Problem* problem_pt)
  {
   // Set MPI info
   MPI_Status status;
   MPI_Request request;

   // Storage for number of processors, current process and communicator
   OomphCommunicator* comm_pt=problem_pt->communicator_pt();
   int n_proc=comm_pt->nproc();
   int my_rank=comm_pt->my_rank();

   // Prepare vectors to receive information
   Vector<double> received_double_values;
   Vector<unsigned> received_unsigned_values;
   Vector<double> received_located_coord;
   Vector<int> received_found_zeta;

   // Communicate the located information back to the original process
   int orig_send_proc=my_rank-iproc;
   if (my_rank<iproc) { orig_send_proc=n_proc+orig_send_proc; }
   int orig_recv_proc=my_rank+iproc;
   if ((my_rank+iproc)>=n_proc) { orig_recv_proc=orig_recv_proc-n_proc; }

   // Store all this processor's Count_* variables to use in sends
   unsigned send_count_double_values=Count_double_values;
   unsigned send_count_zeta_dim=Count_zeta_dim;
   unsigned send_count_located_coord=Count_located_coord;
   unsigned send_count_unsigned_values=Count_unsigned_values;

   // Send the double values associated with external halos
   MPI_Isend(&send_count_double_values,1,MPI_INT,
            orig_send_proc,1,comm_pt->mpi_comm(),&request);
   MPI_Recv(&Count_double_values,1,MPI_INT,
            orig_recv_proc,1,comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_double_values!=0)
    {
     MPI_Isend(&Double_values[0],send_count_double_values,MPI_DOUBLE,
              orig_send_proc,2,comm_pt->mpi_comm(),&request);
    }

   // Receive the double values from the correct processor
   if (Count_double_values!=0)
    {
     received_double_values.resize(Count_double_values);
     MPI_Recv(&received_double_values[0],Count_double_values,
              MPI_DOUBLE,orig_recv_proc,2,comm_pt->mpi_comm(),&status);
    }

   if (send_count_double_values!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Now send unsigned values associated with external halos
   MPI_Isend(&send_count_unsigned_values,1,MPI_INT,
            orig_send_proc,14,comm_pt->mpi_comm(),&request);

   MPI_Recv(&Count_unsigned_values,1,MPI_INT,orig_recv_proc,14,
            comm_pt->mpi_comm(),&status);

   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_unsigned_values!=0)
    {     
     MPI_Isend(&Unsigned_values[0],send_count_unsigned_values,MPI_INT,
              orig_send_proc,15,comm_pt->mpi_comm(),&request);
    }

   // Receive the unsigned values from the correct processor
   if (Count_unsigned_values!=0)
    {
     received_unsigned_values.resize(Count_unsigned_values);
     MPI_Recv(&received_unsigned_values[0],Count_unsigned_values,MPI_INT,
              orig_recv_proc,15,comm_pt->mpi_comm(),&status);
    }

   if (send_count_unsigned_values!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Send and received the Located_element and Found_zeta arrays
   MPI_Isend(&send_count_zeta_dim,1,MPI_INT,orig_send_proc,
            20,comm_pt->mpi_comm(),&request);
   MPI_Recv(&Count_zeta_dim,1,MPI_INT,orig_recv_proc,
            20,comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_zeta_dim!=0)
    {
     MPI_Isend(&Located_element[0],send_count_zeta_dim,MPI_INT,
              orig_send_proc,3,comm_pt->mpi_comm(),&request);
    }
   if (Count_zeta_dim!=0)
    {
     Located_zetas.resize(Count_zeta_dim);
     MPI_Recv(&Located_zetas[0],Count_zeta_dim,MPI_INT,orig_recv_proc,3,
              comm_pt->mpi_comm(),&status);
    }
   if (send_count_zeta_dim!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }


   if (send_count_zeta_dim!=0)
    {
     MPI_Isend(&Found_zeta[0],send_count_zeta_dim,MPI_INT,
              orig_send_proc,13,comm_pt->mpi_comm(),&request);
    }
   if (Count_zeta_dim!=0)
    {
     received_found_zeta.resize(Count_zeta_dim);
     MPI_Recv(&received_found_zeta[0],Count_zeta_dim,MPI_INT,orig_recv_proc,13,
              comm_pt->mpi_comm(),&status);
    }
   if (send_count_zeta_dim!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }


   // And finally the Located_coord array
   MPI_Isend(&send_count_located_coord,1,MPI_INT,
            orig_send_proc,4,comm_pt->mpi_comm(),&request);
   MPI_Recv(&Count_located_coord,1,MPI_INT,orig_recv_proc,4,
            comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_located_coord!=0)
    {     
     MPI_Isend(&Located_coord[0],send_count_located_coord,MPI_DOUBLE,
              orig_send_proc,5,comm_pt->mpi_comm(),&request);
    }
   if (Count_located_coord!=0)
    {
     received_located_coord.resize(Count_located_coord);
     MPI_Recv(&received_located_coord[0],Count_located_coord,MPI_DOUBLE,
              orig_recv_proc,5,comm_pt->mpi_comm(),&status);
    }
   if (send_count_located_coord!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Fill in local arrays again: double and unsigned values
   Double_values.resize(Count_double_values);
   for (unsigned ii=0;ii<Count_double_values;ii++)
    {
     Double_values[ii]=received_double_values[ii];
    }
   Unsigned_values.resize(Count_unsigned_values);
   for (unsigned ii=0;ii<Count_unsigned_values;ii++)
    {
     Unsigned_values[ii]=received_unsigned_values[ii];
    }

   // Found_zeta and Located_coord
   Found_zeta.resize(Count_zeta_dim);
   for (unsigned ii=0;ii<Count_zeta_dim;ii++)
    {
     Found_zeta[ii]=received_found_zeta[ii];
    }
   Located_coord.resize(Count_located_coord);
   for (unsigned ii=0;ii<Count_located_coord;ii++)
    {
     Located_coord[ii]=received_located_coord[ii];
    }

  }


//========start of add_external_haloed_node_to_storage====================
/// Helper function to add external haloed nodes, including any masters
//========================================================================
 void add_external_haloed_node_to_storage(int& iproc, Node* nod_pt,
                                          Problem* problem_pt,
                                          Mesh* const &external_mesh_pt,
                                          int& n_cont_inter_values,
                                          FiniteElement* f_el_pt)
  {
   // Add the node if required
   add_external_haloed_node_helper(iproc,nod_pt,problem_pt,external_mesh_pt,
                                   n_cont_inter_values,f_el_pt);

   // Loop over continuously interpolated values and add masters
   for (int i_cont=-1;i_cont<n_cont_inter_values;i_cont++)
    {
     if (nod_pt->is_hanging(i_cont))
      {
       // Indicate that this node is a hanging node so the other
       // process knows to create HangInfo and masters, etc.
       Unsigned_values.push_back(1);
       Count_unsigned_values++;
       // If this is a hanging node then add all its masters as
       // external halo nodes if they have not yet been added
       HangInfo* hang_pt=nod_pt->hanging_pt(i_cont);
       // Loop over masters
       unsigned n_master=hang_pt->nmaster();
       // Indicate number of master nodes to add on other process
       Unsigned_values.push_back(n_master);
       Count_unsigned_values++;
       for (unsigned m=0;m<n_master;m++)
        {
         Node* master_nod_pt=hang_pt->master_node_pt(m);

         // Call the helper function for master nodes
         add_external_haloed_master_node_helper(iproc,master_nod_pt,
                                                problem_pt,
                                                external_mesh_pt,
                                                n_cont_inter_values);

         // Indicate the weight of this master
         Double_values.push_back(hang_pt->master_weight(m));
         Count_double_values++;
        }
      }
     else
      {
       // Indicate that it's not a hanging node in this variable
       Unsigned_values.push_back(0);
       Count_unsigned_values++;
      }
    } // end loop over continously interpolated values

  }

//==========start of add_external_haloed_node_helper======================
/// Helper to add external haloed node that is not a master
//========================================================================
 void add_external_haloed_node_helper(int& iproc, Node* nod_pt,
                                      Problem* problem_pt,
                                      Mesh* const &external_mesh_pt,
                                      int& n_cont_inter_values,
                                      FiniteElement* f_el_pt)
  {
   // Attempt to add this node as an external haloed node
   unsigned n_ext_haloed_nod=external_mesh_pt->nexternal_haloed_node(iproc);
   unsigned external_haloed_node_index;
   external_haloed_node_index=
    external_mesh_pt->add_external_haloed_node_pt(iproc,nod_pt);

   // If it was added then the new index should match the size of the storage
   if (external_haloed_node_index==n_ext_haloed_nod)
    {
     // Indicate that this node needs to be constructed on
     // the other process
     Unsigned_values.push_back(1);
     Count_unsigned_values++;

     // This helper function gets all the required information for the 
     // specified node and stores it into MPI-sendable information
     // so that a halo copy can be made on the receiving process
     get_required_nodal_information_helper(iproc,nod_pt,
                                           problem_pt,external_mesh_pt,
                                           n_cont_inter_values,f_el_pt);
    }
   else // It was already added
    {
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
     // This node is already an external haloed node, so tell
     // the other process its index in the equivalent external halo storage
     Unsigned_values.push_back(external_haloed_node_index);
     Count_unsigned_values++;
    }
  }


//==========start of add_external_haloed_master_node_helper===============
/// Helper function to add external haloed node that is a master
//========================================================================
 void add_external_haloed_master_node_helper(int& iproc, Node* master_nod_pt,
                                             Problem* problem_pt,
                                             Mesh* const &external_mesh_pt,
                                             int& n_cont_inter_values)
  {
   // Attempt to add node as an external haloed node
   unsigned n_ext_haloed_nod=external_mesh_pt->nexternal_haloed_node(iproc);
   unsigned external_haloed_node_index;
   external_haloed_node_index=
    external_mesh_pt->add_external_haloed_node_pt(iproc,master_nod_pt);

   // If it was added the returned index is the same as current storage size
   if (external_haloed_node_index==n_ext_haloed_nod)
    {
     // Indicate that this node needs to be constructed on
     // the other process
     Unsigned_values.push_back(1);
     Count_unsigned_values++;

     // This gets all the required information for the specified
     // master node and stores it into MPI-sendable information
     // so that a halo copy can be made on the receiving process
     get_required_master_nodal_information_helper(iproc,master_nod_pt,
                                                  problem_pt,
                                                  external_mesh_pt,
                                                  n_cont_inter_values);
    }
   else // It was already added
    {
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
     // This node is already an external haloed node, so tell
     // the other process its index in the equivalent external halo storage
     Unsigned_values.push_back(external_haloed_node_index);
     Count_unsigned_values++;
    }
  }

//========start of get_required_nodal_information_helper==================
/// Helper function to get the required nodal information from an
/// external haloed node so that a fully-functional external halo
/// node (and therefore element) can be created on the receiving process
//========================================================================
 void get_required_nodal_information_helper(int& iproc, Node* nod_pt,
                                            Problem* problem_pt,
                                            Mesh* const &external_mesh_pt,
                                            int& n_cont_inter_values,
                                            FiniteElement* f_el_pt)
  {
   unsigned n_dim=nod_pt->ndim();
   TimeStepper* time_stepper_pt=nod_pt->time_stepper_pt();

   // Find the timestepper in the list of problem timesteppers
   bool found_timestepper=false;
   unsigned time_stepper_index;
   unsigned n_time_steppers=problem_pt->ntime_stepper();
   for (unsigned i=0;i<n_time_steppers;i++)
    {
     if (time_stepper_pt==problem_pt->time_stepper_pt(i))
      {
       // Indicate the timestepper's index
       found_timestepper=true;
       time_stepper_index=i;
       break;
      }
    }

   if (found_timestepper)
    {
     Unsigned_values.push_back(1);
     Count_unsigned_values++;
     Unsigned_values.push_back(time_stepper_index);
     Count_unsigned_values++;
    }
   else
    {
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
    }

   // Default number of previous values to 1
   unsigned n_prev=1;
   if (time_stepper_pt!=0)
    {
     // Add number of history values to n_prev
     n_prev+=time_stepper_pt->nprev_values();
    }

   // Is the node on any boundaries?
   if (nod_pt->is_on_boundary())
    {
     Unsigned_values.push_back(1);
     Count_unsigned_values++;
     // Loop over the boundaries of the external mesh
     unsigned n_bnd=external_mesh_pt->nboundary();
     for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
      {
       // Which boundaries (could be more than one) is it on?
       if (nod_pt->is_on_boundary(i_bnd))
        {
         Unsigned_values.push_back(1);
         Count_unsigned_values++;
        }
       else
        {
         Unsigned_values.push_back(0);
         Count_unsigned_values++;
        }
      }
    }
   else
    {
     // Not on any boundary
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
    }

   // Is the Node algebraic?  If so, send its ref values and
   // an indication of its geometric objects if they are stored
   // in the algebraic mesh
   AlgebraicNode* alg_nod_pt=dynamic_cast<AlgebraicNode*>(nod_pt);
   if (alg_nod_pt!=0)
    {
     // The external mesh should be algebraic
     AlgebraicMesh* alg_mesh_pt=dynamic_cast<AlgebraicMesh*>
      (external_mesh_pt);

     // Get default node update function ID
     unsigned update_id=alg_nod_pt->node_update_fct_id();
     Unsigned_values.push_back(update_id);
     Count_unsigned_values++;

     // Get reference values at default...
     unsigned n_ref_val=alg_nod_pt->nref_value();
     Unsigned_values.push_back(n_ref_val);
     Count_unsigned_values++;
     for (unsigned i_ref_val=0;i_ref_val<n_ref_val;i_ref_val++)
      {
       Double_values.push_back(alg_nod_pt->ref_value(i_ref_val));
       Count_double_values++;
      }

     // Access geometric objects at default...
     unsigned n_geom_obj=alg_nod_pt->ngeom_object();
     Unsigned_values.push_back(n_geom_obj);
     Count_unsigned_values++;
     for (unsigned i_geom=0;i_geom<n_geom_obj;i_geom++)
      {
       GeomObject* geom_obj_pt=alg_nod_pt->geom_object_pt(i_geom);
       // Check this against the stored geometric objects in mesh
       unsigned n_geom_list=alg_mesh_pt->ngeom_object_list_pt();
       // Default found index to zero
       unsigned found_geom_object=0;
       for (unsigned i_list=0;i_list<n_geom_list;i_list++)
        {
         if (geom_obj_pt==alg_mesh_pt->geom_object_list_pt(i_list))
          {
           found_geom_object=i_list;
          }
        }
       Unsigned_values.push_back(found_geom_object);
       Count_unsigned_values++;
      }
    }

   // If it is a MacroElementNodeUpdateNode, everything has been
   // dealt with by the new element already

   // Is it a SolidNode?
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt!=0)
    {
     unsigned n_val=solid_nod_pt->variable_position_pt()->nvalue();
     for (unsigned i_val=0;i_val<n_val;i_val++)
      {
       for (unsigned t=0;t<n_prev;t++)
        {
         Double_values.push_back(solid_nod_pt->variable_position_pt()->
                          value(t,i_val));
         Count_double_values++;
        }
      }
    }

   // Finally copy info required for all node types

   // Halo copy needs to know all the history values
   unsigned n_val=nod_pt->nvalue();
   for (unsigned i_val=0;i_val<n_val;i_val++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Double_values.push_back(nod_pt->value(t,i_val));
       Count_double_values++;
      }
    }

   // Now do positions
   for (unsigned idim=0;idim<n_dim;idim++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Double_values.push_back(nod_pt->x(t,idim));
       Count_double_values++;
      }
    }


  }

//=========start of get_required_master_nodal_information_helper==========
/// Helper function to get the required master nodal information from an
/// external haloed master node so that a fully-functional external halo
/// master node (and possible element) can be created on the receiving process
//========================================================================
 void get_required_master_nodal_information_helper
 (int& iproc, Node* master_nod_pt, Problem* problem_pt,
  Mesh* const &external_mesh_pt, int& n_cont_inter_values)
  {
   // Need to send over dimension, position type and number of values
   Unsigned_values.push_back(master_nod_pt->ndim());
   Count_unsigned_values++;
   Unsigned_values.push_back(master_nod_pt->nposition_type());
   Count_unsigned_values++;
   Unsigned_values.push_back(master_nod_pt->nvalue());
   Count_unsigned_values++;
   
   // If it's a solid node, also need to send lagrangian dim and type
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(master_nod_pt);
   if (solid_nod_pt!=0)
    {
     Unsigned_values.push_back(solid_nod_pt->nlagrangian());
     Count_unsigned_values++;
     Unsigned_values.push_back(solid_nod_pt->nlagrangian_type());
     Count_unsigned_values++;
    }

   unsigned n_dim=master_nod_pt->ndim();
   TimeStepper* time_stepper_pt=master_nod_pt->time_stepper_pt();

   // Find the timestepper in the list of problem timesteppers
   bool found_timestepper=false;
   unsigned time_stepper_index;
   unsigned n_time_steppers=problem_pt->ntime_stepper();
   for (unsigned i=0;i<n_time_steppers;i++)
    {
     if (time_stepper_pt==problem_pt->time_stepper_pt(i))
      {
       // Indicate the timestepper's index
       // add 1 to the index so that 0 indicates no timestepper?
       found_timestepper=true;
       time_stepper_index=i;
       break;
      }
    }

   if (found_timestepper)
    {
     Unsigned_values.push_back(1);
     Count_unsigned_values++;
     Unsigned_values.push_back(time_stepper_index);
     Count_unsigned_values++;
    }
   else
    {
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
    }

   // Default number of previous values to 1
   unsigned n_prev=1;
   if (time_stepper_pt!=0)
    {
     // Add number of history values to n_prev
     n_prev+=time_stepper_pt->nprev_values();
    }

   // Is the node on any boundaries?
   if (master_nod_pt->is_on_boundary())
    {
     Unsigned_values.push_back(1);
     Count_unsigned_values++;
     // Loop over the boundaries of the external mesh
     unsigned n_bnd=external_mesh_pt->nboundary();
     for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
      {
       // Which boundaries (could be more than one) is it on?
       if (master_nod_pt->is_on_boundary(i_bnd))
        {
         Unsigned_values.push_back(1);
         Count_unsigned_values++;
        }
       else
        {
         Unsigned_values.push_back(0);
         Count_unsigned_values++;
        }
      }
    }
   else
    {
     // Not on any boundary
     Unsigned_values.push_back(0);
     Count_unsigned_values++;
    }

   // Is the Node algebraic?  If so, send its ref values and
   // an indication of its geometric objects if they are stored
   // in the algebraic mesh
   AlgebraicNode* alg_nod_pt=dynamic_cast<AlgebraicNode*>(master_nod_pt);
   if (alg_nod_pt!=0)
    {
     // The external mesh should be algebraic
     AlgebraicMesh* alg_mesh_pt=dynamic_cast<AlgebraicMesh*>
      (external_mesh_pt);

     // Get default node update function ID
     unsigned update_id=alg_nod_pt->node_update_fct_id();
     Unsigned_values.push_back(update_id);
     Count_unsigned_values++;

     // Get reference values at default...
     unsigned n_ref_val=alg_nod_pt->nref_value();
     Unsigned_values.push_back(n_ref_val);
     Count_unsigned_values++;
     for (unsigned i_ref_val=0;i_ref_val<n_ref_val;i_ref_val++)
      {
       Double_values.push_back(alg_nod_pt->ref_value(i_ref_val));
       Count_double_values++;
      }

     // Access geometric objects at default...
     unsigned n_geom_obj=alg_nod_pt->ngeom_object();
     Unsigned_values.push_back(n_geom_obj);
     Count_unsigned_values++;
     for (unsigned i_geom=0;i_geom<n_geom_obj;i_geom++)
      {
       GeomObject* geom_obj_pt=alg_nod_pt->geom_object_pt(i_geom);
       // Check this against the stored geometric objects in mesh
       unsigned n_geom_list=alg_mesh_pt->ngeom_object_list_pt();
       // Default found index to zero
       unsigned found_geom_object=0;
       for (unsigned i_list=0;i_list<n_geom_list;i_list++)
        {
         if (geom_obj_pt==alg_mesh_pt->geom_object_list_pt(i_list))
          {
           found_geom_object=i_list;
          }
        }
       Unsigned_values.push_back(found_geom_object);
       Count_unsigned_values++;
      }
    } // end AlgebraicNode check

   // Is it a MacroElementNodeUpdateNode?
   MacroElementNodeUpdateNode* macro_nod_pt=
    dynamic_cast<MacroElementNodeUpdateNode*>(master_nod_pt);
   if (macro_nod_pt!=0)
    {
     // Loop over current external haloed elements - has the element which
     // controls the node update for this node been added yet?
     FiniteElement* macro_node_update_el_pt=
      macro_nod_pt->node_update_element_pt();

     unsigned n_ext_haloed_el=external_mesh_pt->
      nexternal_haloed_element(iproc);
     unsigned external_haloed_el_index;
     external_haloed_el_index=external_mesh_pt->
      add_external_haloed_element_pt(iproc,macro_node_update_el_pt);

     // If it wasn't already added, we need to create a halo copy
     if (external_haloed_el_index==n_ext_haloed_el)
      {
       Unsigned_values.push_back(1);
       Count_unsigned_values++;

       // We're using macro elements to update...
       MacroElementNodeUpdateMesh* macro_mesh_pt=
        dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
       if (macro_mesh_pt!=0)
        {
         Unsigned_values.push_back(1);
         Count_unsigned_values++;

         // Need to send the macro element number in the mesh across
         MacroElement* macro_el_pt=macro_node_update_el_pt->macro_elem_pt();
         unsigned macro_el_num=macro_el_pt->macro_element_number();
         Unsigned_values.push_back(macro_el_num);
         Count_unsigned_values++;

         // Also need to send
         // the lower left and upper right coordinates of the macro element
         QElementBase* q_el_pt=dynamic_cast<QElementBase*>
          (macro_node_update_el_pt);
         if (q_el_pt!=0)
          {
           // The macro element needs to be set first before
           // its lower left and upper right coordinates can be accessed
           // Now send the lower left and upper right coordinates
           unsigned el_dim=q_el_pt->dim();
           for (unsigned i_dim=0;i_dim<el_dim;i_dim++)
            {
             Double_values.push_back(q_el_pt->s_macro_ll(i_dim));
             Count_double_values++;
             Double_values.push_back(q_el_pt->s_macro_ur(i_dim));
             Count_double_values++;
            }
          }
         else // Throw an error
          {
           std::ostringstream error_stream;
           error_stream << "You are using a MacroElement node update\n"
                        << "in a case with non-QElements. This has not\n"
                        << "yet been implemented.\n";
           throw OomphLibError
            (error_stream.str(),
             "Multi_domain_functions::get_required_master_nodal_...()",
             OOMPH_EXCEPTION_LOCATION);
          }

        }
       else // Not using macro elements for node update... umm, we're
            // already inside a loop over macro elements, so this
            // should never get here... an error should be thrown I suppose
        {
         Unsigned_values.push_back(0);
         Count_unsigned_values++;
        }

       // This element needs to be fully functioning on the other
       // process, so send all the information required to create it
       unsigned n_node=macro_node_update_el_pt->nnode();
       for (unsigned j=0;j<n_node;j++)
        {
         Node* new_nod_pt=macro_node_update_el_pt->node_pt(j);
         add_external_haloed_node_to_storage(iproc,new_nod_pt,
                                             problem_pt,
                                             external_mesh_pt,
                                             n_cont_inter_values,
                                             macro_node_update_el_pt);
        }
      }
     else // The external haloed element already exists
      {
       Unsigned_values.push_back(0);
       Count_unsigned_values++;

       Unsigned_values.push_back(external_haloed_el_index);
       Count_unsigned_values++;
      }

    } // end of MacroElementNodeUpdateNode check 

   // Is it a SolidNode?
   if (solid_nod_pt!=0)
    {
     unsigned n_val=solid_nod_pt->variable_position_pt()->nvalue();
     for (unsigned i_val=0;i_val<n_val;i_val++)
      {
       for (unsigned t=0;t<n_prev;t++)
        {
         Double_values.push_back(solid_nod_pt->variable_position_pt()->
                          value(t,i_val));
         Count_double_values++;
        }
      }
    }

   // Finally copy info required for all node types

   // Halo copy needs to know all the history values
   unsigned n_val=master_nod_pt->nvalue();
   for (unsigned i_val=0;i_val<n_val;i_val++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Double_values.push_back(master_nod_pt->value(t,i_val));
       Count_double_values++;
      }
    }

   // Now do positions
   for (unsigned idim=0;idim<n_dim;idim++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Double_values.push_back(master_nod_pt->x(t,idim));
       Count_double_values++;
      }
    }

  }

#endif

 /// Helper function that returns the dimension of the elements within
 /// each of the specified meshes (and checks they are the same)
 void get_dim_helper(Problem* problem_pt, Mesh* const &mesh_pt, 
                     Mesh* const &external_mesh_pt, unsigned& dim)
  {
#ifdef OOMPH_HAS_MPI
   // Storage for number of processors, current process and communicator
   OomphCommunicator* comm_pt=problem_pt->communicator_pt();
#endif

   // Extract the element dimensions from the first element of each mesh
   unsigned mesh_dim=0;
   if (mesh_pt->nelement() > 0)
    {
     mesh_dim=
      dynamic_cast<FiniteElement*>(mesh_pt->element_pt(0))->dim();
    }
   unsigned external_mesh_dim=0;
   if (external_mesh_pt->nelement() > 0)
    {
     external_mesh_dim=
      dynamic_cast<FiniteElement*>(external_mesh_pt->element_pt(0))->dim();
    }

   // Need to do an Allreduce
#ifdef OOMPH_HAS_MPI
   int n_proc=comm_pt->nproc();
   if (n_proc > 1)
    {
     unsigned mesh_dim_reduce;
     MPI_Allreduce(&mesh_dim,&mesh_dim_reduce,1,MPI_INT,
                   MPI_MAX,comm_pt->mpi_comm());
     mesh_dim=mesh_dim_reduce;

     unsigned external_mesh_dim_reduce;
     MPI_Allreduce(&external_mesh_dim,&external_mesh_dim_reduce,1,MPI_INT,
                   MPI_MAX,comm_pt->mpi_comm());
     external_mesh_dim=external_mesh_dim_reduce;
    }
#endif

   // Check the dimensions are the same!
   if (mesh_dim!=external_mesh_dim)
    {
     std::ostringstream error_stream;
     error_stream << "The elements within the two meshes do not\n"
                  << "have the same dimension, so the multi-domain\n"
                  << "method will not work.\n"
                  << "For the mesh, dim=" << mesh_dim 
                  << ", and the external mesh, dim=" << external_mesh_dim 
                  << "\n";
     throw OomphLibError
      (error_stream.str(),
       "Multi_domain_functions::get_dim_helper(...)",
       OOMPH_EXCEPTION_LOCATION);
    }

   // Set returned dimension
   dim=mesh_dim;
 }

 /// Helper function that clears all the information used
 /// during the external storage creation
 void clean_up()
  {
   // Clear every vector associated with the external storage creation
   Local_zetas.clear();
   Zetas.clear();
   Local_zeta_dim.clear();
   Zeta_dim.clear();

   Found_zeta.clear();
   Found_ss.clear();

   Located_element.clear();
   Located_zetas.clear();
   Located_coord.clear();

   Double_values.clear();
   Unsigned_values.clear();

   External_element_located.clear();
  }


}

}

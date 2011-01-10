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


  // Workspace for locate zeta methods
  //----------------------------------

  /// \short Dimension of zeta tuples (set by get_dim_helper) -- needed
  /// because we store the scalar coordinates in flat-packed form.
  unsigned Dim;

  /// \short Lookup scheme for whether a local element's integration point
  /// has had an external element assigned to it -- essentially boolean.
  /// External_element_located[e][ipt] = {0,1} if external element
  /// for ipt-th integration in local element e {has not, has} been found.
  /// Used locally to ensure that we're not searching for the same
  /// elements over and over again when we go around the spirals.
  Vector<Vector<unsigned> > External_element_located;
  
  /// \short Vector of flat-packed zeta coordinates for which the external
  /// element could not be found during current local search. These
  /// will be sent to the next processor in the ring-like parallel search.
  /// The zeta coordinates come in groups of Dim (scalar) coordinates.
  Vector<double> Flat_packed_zetas_not_found_locally;

  /// \short Vector of flat-packed zeta coordinates for which the external
  /// element could not be found on another processor and for which
  /// we're currently searching here. Whatever can't be found here,
  /// gets written into Flat_packed_zetas_not_found_locally and then
  /// passed on to the next processor during the ring-like parallel search.
  /// The zeta coordinates come in  groups of Dim (scalar) coordinates.
  Vector<double> Received_flat_packed_zetas_to_be_found;

  /// \short Proc_id_plus_one_of_external_element[i] contains the 
  /// processor id (plus one) of the processor
  /// on which the i-th zeta coordinate tuple received from elsewhere 
  /// (in the order in which these are stored in 
  /// Received_flat_packed_zetas_to_be_found) was located; it's zero if
  /// it wasn't found during the current stage of the ring-like parallel
  /// search.
  Vector<int> Proc_id_plus_one_of_external_element;

  /// \short Vector to indicate (to another processor) whether a
  /// located element (that will have to represented as an external
  /// halo element on that processor) should be newly created on that 
  /// processor (2), already exists on that processor (1), or
  /// is not on the current processor either (0).
  Vector<unsigned> Located_element_status;

  /// \short Vector of flat-packed local coordinates for zeta tuples
  /// that have been located
  Vector<double> Flat_packed_located_coordinates;

  /// \short Vector of flat-packed doubles to be communicated with
  /// other processors
  Vector<double> Flat_packed_doubles;

  /// \short Counter used when processing vector of flat-packed 
  /// doubles -- this is really "private" data, declared here
  /// to avoid having to pass it (and the associated array)
  /// between the various helper functions
  unsigned Counter_for_flat_packed_doubles;

  /// \short Vector of flat-packed unsigneds to be communicated with
  /// other processors -- this is really "private" data, declared here
  /// to avoid having to pass the array between the various helper 
  /// functions
  Vector<unsigned> Flat_packed_unsigneds;

  /// \short Counter used when processing vector of flat-packed 
  /// unsigneds -- this is really "private" data, declared here
  /// to avoid having to pass it (and the associated array)
  /// between the various helper functions
  unsigned Counter_for_flat_packed_unsigneds;


  // Default parameters for the binning method
  //------------------------------------------

  /// \short Bool to tell the MeshAsGeomObject whether to calculate
  /// the extreme coordinates of the bin structure
  bool Compute_extreme_bin_coordinates=true;

  /// \short Number of bins in the first dimension in binning method in
  /// setup_multi_domain_interaction(). 
  unsigned Nx_bin=1000;

  /// \short Number of bins in the second dimension in binning method in
  /// setup_multi_domain_interaction(). 
  unsigned Ny_bin=1000;

  /// \short Number of bins in the third dimension in binning method in
  /// setup_multi_domain_interaction(). 
  unsigned Nz_bin=1000;

  /// Number of spirals to be searched in one go
  unsigned N_spiral_chunk=1;

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


  // Other parameters
  //-----------------

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

   // Doc timings if required
   double t_start=0.0;
   if (Doc_timings)
    {
     t_start=TimingHelpers::timer();
    }

   // Loop over existing "normal" elements in mesh
   unsigned n_element=mesh_pt->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     FiniteElement* el_pt = mesh_pt->finite_element_pt(e);
     if (el_pt!=0)
      {
       // Loop over nodes
       unsigned n_node=el_pt->nnode();
       for (unsigned j=0;j<n_node;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);
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
         if (dynamic_cast<RefineableElement*>(el_pt)!=0)
          {
           int n_cont_int_values=dynamic_cast<RefineableElement*>
            (el_pt)->ncont_interpolated_values();
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
                 
                 // If this master is a SolidNode then add its extra 
                 // eqn numbers
                 SolidNode* master_solid_nod_pt=dynamic_cast<SolidNode*>
                  (master_nod_pt);
                 if (master_solid_nod_pt!=0)
                  {
                   unsigned n_val_mst_solid=master_solid_nod_pt->
                    variable_position_pt()->nvalue();
                   Vector<int> master_solid_nodal_eqn_numbers(n_val_mst_solid);
                   for (unsigned i_val=0;i_val<n_val_mst_solid;i_val++)
                    {
                     master_solid_nodal_eqn_numbers[i_val]=
                      master_solid_nod_pt->
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
      } //End of FiniteElement data

     // Internal data equation numbers do not need to be added since 
     // internal data cannot be shared between distinct elements, so 
     // internal data on locally-stored elements can never be halo.
    }
   

   // Doc timings if required
   if (Doc_timings)
    {
     double t_locate=TimingHelpers::timer();
     oomph_info 
      <<"CPU for assembly of existing global eqns in remove_duplicate_data(): "
      << t_locate-t_start << std::endl;
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
         GeneralisedElement* ext_el_pt=mesh_pt->
          external_halo_element_pt(iproc,e_ext);

         FiniteElement* finite_ext_el_pt = 
          dynamic_cast<FiniteElement*>(ext_el_pt);
         
         if(finite_ext_el_pt!=0)
          {
           // Loop over nodes
           unsigned n_node=finite_ext_el_pt->nnode();
           for (unsigned j=0;j<n_node;j++)
            {
             Node* nod_pt=finite_ext_el_pt->node_pt(j);
             unsigned n_val=nod_pt->nvalue();
             Vector<int> nodal_eqn_numbers(n_val);
             for (unsigned i_val=0;i_val<n_val;i_val++)
              {
               nodal_eqn_numbers[i_val]=nod_pt->eqn_number(i_val);
               
              }
             // Check for duplicates with the existing set of 
             //global eqn numbers
             
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
                       finite_ext_el_pt->node_pt(j)=
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
               // These external halo nodes need to be unclassified in 
               // order for
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
                 if (n_val_solid==
                     (existing_global_eqn_numbers[i].first).size())
                  {
                   // Loop over values
                   for (unsigned i_val=0;i_val<n_val_solid;i_val++)
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
                 for (unsigned i_val=0;i_val<n_val_solid;i_val++)
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
                               finite_ext_el_pt->node_pt(j)->hanging_pt(i_cont)
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
                       // they receive the correct global equation numbers 
                       // during
                       // the synchronisation process in
                       // Problem::copy_external_haloed_eqn_numbers_helper(...)
                       for (unsigned i_val=0;i_val<n_val;i_val++)
                        {
                         master_nod_pt->eqn_number(i_val)
                          =Data::Is_unclassified;
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
                         sld_master_nodal_eqn_numbers[i_val]=
                          solid_master_nod_pt
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
                         if (n_val_solid==
                             (existing_global_eqn_numbers[i].first).size())
                          {
                           // Loop over values
                           for (unsigned i_val=0;i_val<n_val_solid;i_val++)
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
                                 // It's a duplicate, 
                                 // so point the current master
                                 // node at the original node instead!
                                 // Need the weight of the original node
                                 double m_weight=hang_pt->master_weight(m);
                                 // Set master
                                 finite_ext_el_pt->node_pt(j)->
                                  hanging_pt(i_cont)
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
                       
                       // If it's not a duplicate then 
                       // add it to existing storage
                       if (!is_a_duplicate)
                        {
                         existing_global_eqn_numbers.push_back
                          (make_pair(sld_master_nodal_eqn_numbers,
                                     master_nod_pt));
                         // These external halo nodes need to be 
                         // unclassified in
                         // order for them to be bypassed in 
                         // assign__eqn_numbers;
                         // they receive the correct global equation numbers
                         // during the synchronisation process in 
                         // Problem::copy_external_haloed_eqn_numbers_helper(.)
                         for (unsigned i_val=0;i_val<n_val_solid;i_val++)
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
           
          } //End of check for finite element
 
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

   // Send the number  of flat-packed zetas that we couldn't find
   // locally to the next processor
   int n_missing_local_zetas=Flat_packed_zetas_not_found_locally.size();
   MPI_Isend(&n_missing_local_zetas,1,MPI_INT,send_to_proc,4,
             comm_pt->mpi_comm(),&request);

   // Receive the number of flat-packed zetas that couldn't be found
   // on the "previous" processor
   int count_zetas=0;
   MPI_Recv(&count_zetas,1,MPI_INT,recv_from_proc,
            4,comm_pt->mpi_comm(),&status);

   MPI_Wait(&request,MPI_STATUS_IGNORE);


   // Send the vector of flat-packed zetas that we couldn't find
   // locally to the next processor
   if (n_missing_local_zetas!=0)
    {           
     MPI_Isend(&Flat_packed_zetas_not_found_locally[0],
               n_missing_local_zetas,MPI_DOUBLE,
               send_to_proc,5,comm_pt->mpi_comm(),&request);
    }
   
   // Receive the vector of flat-packed zetas that couldn't be found
   // on the "previous" processor
   if (count_zetas!=0)
    {       
     Received_flat_packed_zetas_to_be_found.resize(count_zetas);
     MPI_Recv(&Received_flat_packed_zetas_to_be_found[0],
              count_zetas,MPI_DOUBLE,recv_from_proc,5,
              comm_pt->mpi_comm(),&status);
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
   
   // Now we should have the Zeta arrays set up correctly
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
   Vector<int> received_proc_id_plus_one_of_external_element;
   Vector<unsigned> received_located_element_status;

   // Communicate the located information back to the original process
   int orig_send_proc=my_rank-iproc;
   if (my_rank<iproc) { orig_send_proc=n_proc+orig_send_proc; }
   int orig_recv_proc=my_rank+iproc;
   if ((my_rank+iproc)>=n_proc) { orig_recv_proc=orig_recv_proc-n_proc; }



   // Send the double values associated with external halos
   //------------------------------------------------------
   unsigned send_count_double_values=Flat_packed_doubles.size(); 
   MPI_Isend(&send_count_double_values,1,MPI_INT,
            orig_send_proc,1,comm_pt->mpi_comm(),&request);
   int receive_count_double_values=0;
   MPI_Recv(&receive_count_double_values,1,MPI_INT,
            orig_recv_proc,1,comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_double_values!=0)
    {
     MPI_Isend(&Flat_packed_doubles[0],send_count_double_values,MPI_DOUBLE,
              orig_send_proc,2,comm_pt->mpi_comm(),&request);
    }
   if (receive_count_double_values!=0)
    {
     received_double_values.resize(receive_count_double_values);
     MPI_Recv(&received_double_values[0],receive_count_double_values,
              MPI_DOUBLE,orig_recv_proc,2,comm_pt->mpi_comm(),&status);
    }
   if (send_count_double_values!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Now send unsigned values associated with external halos
   //---------------------------------------------------------
   unsigned send_count_unsigned_values=Flat_packed_unsigneds.size(); 
   MPI_Isend(&send_count_unsigned_values,1,MPI_INT,
            orig_send_proc,14,comm_pt->mpi_comm(),&request);

   int receive_count_unsigned_values=0;
   MPI_Recv(&receive_count_unsigned_values,1,MPI_INT,orig_recv_proc,14,
            comm_pt->mpi_comm(),&status);

   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_unsigned_values!=0)
    {     
     MPI_Isend(&Flat_packed_unsigneds[0],send_count_unsigned_values,MPI_INT,
              orig_send_proc,15,comm_pt->mpi_comm(),&request);
    }
   if (receive_count_unsigned_values!=0)
    {
     received_unsigned_values.resize(receive_count_unsigned_values);
     MPI_Recv(&received_unsigned_values[0],receive_count_unsigned_values,
              MPI_INT,orig_recv_proc,15,comm_pt->mpi_comm(),&status);
    }

   if (send_count_unsigned_values!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   // Send and receive the Located_element_status
   //--------------------------------------------
   int send_count=Received_flat_packed_zetas_to_be_found.size()/Dim;
   MPI_Isend(&send_count,1,MPI_INT,orig_send_proc,
             20,comm_pt->mpi_comm(),&request);
   int receive_count=0;
   MPI_Recv(&receive_count,1,MPI_INT,orig_recv_proc,
            20,comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);
   
   if (send_count!=0)
    {
     MPI_Isend(&Located_element_status[0],send_count,MPI_INT,
               orig_send_proc,3,comm_pt->mpi_comm(),&request);
    }
   
   if (receive_count!=0)
    {
     received_located_element_status.resize(receive_count);
     MPI_Recv(&received_located_element_status[0],receive_count,
              MPI_INT,orig_recv_proc,3,
              comm_pt->mpi_comm(),&status);
    }
   if (send_count!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
   
   // Send and receive Proc_id_plus_one_of_external_element array
   //------------------------------------------------------------
   if (send_count!=0)
    {
     MPI_Isend(&Proc_id_plus_one_of_external_element[0],send_count,MPI_INT,
              orig_send_proc,13,comm_pt->mpi_comm(),&request);
    }
   if (receive_count!=0)
    {
     received_proc_id_plus_one_of_external_element.resize(receive_count);
     MPI_Recv(&received_proc_id_plus_one_of_external_element[0],
              receive_count,MPI_INT,orig_recv_proc,13,
              comm_pt->mpi_comm(),&status);
    }
   if (send_count!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

   
   // And finally the Flat_packed_located_coordinates array
   //------------------------------------------------------
   unsigned send_count_located_coord=Flat_packed_located_coordinates.size();
   MPI_Isend(&send_count_located_coord,1,MPI_INT,
            orig_send_proc,4,comm_pt->mpi_comm(),&request);
   int receive_count_located_coord=0;
   MPI_Recv(&receive_count_located_coord,1,MPI_INT,orig_recv_proc,4,
            comm_pt->mpi_comm(),&status);
   MPI_Wait(&request,MPI_STATUS_IGNORE);

   if (send_count_located_coord!=0)
    {     
     MPI_Isend(&Flat_packed_located_coordinates[0],
               send_count_located_coord,MPI_DOUBLE,
               orig_send_proc,5,comm_pt->mpi_comm(),&request);
    }
   if (receive_count_located_coord!=0)
    {
     received_located_coord.resize(receive_count_located_coord);
     MPI_Recv(&received_located_coord[0],receive_count_located_coord,MPI_DOUBLE,
              orig_recv_proc,5,comm_pt->mpi_comm(),&status);
    }
   if (send_count_located_coord!=0)
    {
     MPI_Wait(&request,MPI_STATUS_IGNORE);
    }



   // Copy across into original containers -- these can now
   //------------------------------------------------------
   // be processed by create_external_halo_elements() to generate
   //------------------------------------------------------------
   // external halo elements
   //------------------------
   Flat_packed_doubles.resize(receive_count_double_values);
   for (int ii=0;ii<receive_count_double_values;ii++)
    {
     Flat_packed_doubles[ii]=received_double_values[ii];
    }
   Flat_packed_unsigneds.resize(receive_count_unsigned_values);
   for (int ii=0;ii<receive_count_unsigned_values;ii++)
    {
     Flat_packed_unsigneds[ii]=received_unsigned_values[ii];
    }
   Proc_id_plus_one_of_external_element.resize(receive_count);
   Located_element_status.resize(receive_count);
   for (int ii=0;ii<receive_count;ii++)
    {
     Proc_id_plus_one_of_external_element[ii]=
      received_proc_id_plus_one_of_external_element[ii];
     Located_element_status[ii]=received_located_element_status[ii];
    }
   Flat_packed_located_coordinates.resize(receive_count_located_coord);
   for (int ii=0;ii<receive_count_located_coord;ii++)
    {
     Flat_packed_located_coordinates[ii]=received_located_coord[ii];
    }
   
  }


//========start of add_external_haloed_node_to_storage====================
/// Helper function to add external haloed nodes, including any masters
//========================================================================
 void add_external_haloed_node_to_storage(int& iproc, Node* nod_pt,
                                          Problem* problem_pt,
                                          Mesh* const &external_mesh_pt,
                                          int& n_cont_inter_values)
  {
   // Add the node if required
   add_external_haloed_node_helper(iproc,nod_pt,problem_pt,external_mesh_pt,
                                   n_cont_inter_values);

   // Loop over continuously interpolated values and add masters
   for (int i_cont=-1;i_cont<n_cont_inter_values;i_cont++)
    {
     if (nod_pt->is_hanging(i_cont))
      {
       // Indicate that this node is a hanging node so the other
       // process knows to create HangInfo and masters, etc.
       Flat_packed_unsigneds.push_back(1);

       // If this is a hanging node then add all its masters as
       // external halo nodes if they have not yet been added
       HangInfo* hang_pt=nod_pt->hanging_pt(i_cont);
       // Loop over masters
       unsigned n_master=hang_pt->nmaster();
       // Indicate number of master nodes to add on other process
       Flat_packed_unsigneds.push_back(n_master);
       for (unsigned m=0;m<n_master;m++)
        {
         Node* master_nod_pt=hang_pt->master_node_pt(m);

         // Call the helper function for master nodes
         add_external_haloed_master_node_helper(iproc,master_nod_pt,
                                                problem_pt,
                                                external_mesh_pt,
                                                n_cont_inter_values);

         // Indicate the weight of this master
         Flat_packed_doubles.push_back(hang_pt->master_weight(m));
        }
      }
     else
      {
       // Indicate that it's not a hanging node in this variable
       Flat_packed_unsigneds.push_back(0);
      }
    } // end loop over continously interpolated values

  }

//==========start of add_external_haloed_node_helper======================
/// Helper to add external haloed node that is not a master
//========================================================================
 void add_external_haloed_node_helper(int& iproc, Node* nod_pt,
                                      Problem* problem_pt,
                                      Mesh* const &external_mesh_pt,
                                      int& n_cont_inter_values)
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
     Flat_packed_unsigneds.push_back(1);

     // This helper function gets all the required information for the 
     // specified node and stores it into MPI-sendable information
     // so that a halo copy can be made on the receiving process
     get_required_nodal_information_helper(iproc,nod_pt,
                                           problem_pt,external_mesh_pt,
                                           n_cont_inter_values);
    }
   else // It was already added
    {
     Flat_packed_unsigneds.push_back(0);
     // This node is already an external haloed node, so tell
     // the other process its index in the equivalent external halo storage
     Flat_packed_unsigneds.push_back(external_haloed_node_index);
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
     Flat_packed_unsigneds.push_back(1);

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
     Flat_packed_unsigneds.push_back(0);

     // This node is already an external haloed node, so tell
     // the other process its index in the equivalent external halo storage
     Flat_packed_unsigneds.push_back(external_haloed_node_index);
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
                                            int& n_cont_inter_values)
  {
   // Tell the halo copy of this node how many values there are
   // [NB this may be different for nodes within the same element, e.g.
   //  when using Lagrange multipliers]
   unsigned n_val=nod_pt->nvalue();
   Flat_packed_unsigneds.push_back(n_val);

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
     Flat_packed_unsigneds.push_back(1);
     Flat_packed_unsigneds.push_back(time_stepper_index);
    }
   else
    {
     Flat_packed_unsigneds.push_back(0);
    }

   // Default number of previous values to 1
   unsigned n_prev=1;
   if (time_stepper_pt!=0)
    {
     // Add number of history values to n_prev
     n_prev=time_stepper_pt->ntstorage(); 
    }

   // Is the node on any boundaries?
   if (nod_pt->is_on_boundary())
    {
     Flat_packed_unsigneds.push_back(1);
     // Loop over the boundaries of the external mesh
     unsigned n_bnd=external_mesh_pt->nboundary();
     for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
      {
       // Which boundaries (could be more than one) is it on?
       if (nod_pt->is_on_boundary(i_bnd))
        {
         Flat_packed_unsigneds.push_back(1);
        }
       else
        {
         Flat_packed_unsigneds.push_back(0);
        }
      }
    }
   else
    {
     // Not on any boundary
     Flat_packed_unsigneds.push_back(0);
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
     Flat_packed_unsigneds.push_back(update_id);

     // Get reference values at default...
     unsigned n_ref_val=alg_nod_pt->nref_value();
     Flat_packed_unsigneds.push_back(n_ref_val);
     for (unsigned i_ref_val=0;i_ref_val<n_ref_val;i_ref_val++)
      {
       Flat_packed_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
      }

     // Access geometric objects at default...
     unsigned n_geom_obj=alg_nod_pt->ngeom_object();
     Flat_packed_unsigneds.push_back(n_geom_obj);
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
       Flat_packed_unsigneds.push_back(found_geom_object);
      }
    }

   // If it is a MacroElementNodeUpdateNode, everything has been
   // dealt with by the new element already

   // Is it a SolidNode?
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt!=0)
    {
     unsigned n_solid_val=solid_nod_pt->variable_position_pt()->nvalue();
     for (unsigned i_val=0;i_val<n_solid_val;i_val++)
      {
       for (unsigned t=0;t<n_prev;t++)
        {
         Flat_packed_doubles.push_back(solid_nod_pt->variable_position_pt()->
                          value(t,i_val));
        }
      }
    }

   // Finally copy info required for all node types
   for (unsigned i_val=0;i_val<n_val;i_val++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Flat_packed_doubles.push_back(nod_pt->value(t,i_val));
      }
    }

   // Now do positions
   for (unsigned idim=0;idim<n_dim;idim++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Flat_packed_doubles.push_back(nod_pt->x(t,idim));
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
   Flat_packed_unsigneds.push_back(master_nod_pt->ndim());
   Flat_packed_unsigneds.push_back(master_nod_pt->nposition_type());
   Flat_packed_unsigneds.push_back(master_nod_pt->nvalue());
   
   // If it's a solid node, also need to send lagrangian dim and type
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(master_nod_pt);
   if (solid_nod_pt!=0)
    {
     Flat_packed_unsigneds.push_back(solid_nod_pt->nlagrangian());
     Flat_packed_unsigneds.push_back(solid_nod_pt->nlagrangian_type());
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
     Flat_packed_unsigneds.push_back(1);
     Flat_packed_unsigneds.push_back(time_stepper_index);
    }
   else
    {
     Flat_packed_unsigneds.push_back(0);
    }

   // Default number of previous values to 1
   unsigned n_prev=1;
   if (time_stepper_pt!=0)
    {
     // Add number of history values to n_prev
     n_prev=time_stepper_pt->ntstorage();
    }

   // Is the node on any boundaries?
   if (master_nod_pt->is_on_boundary())
    {
     Flat_packed_unsigneds.push_back(1);
     // Loop over the boundaries of the external mesh
     unsigned n_bnd=external_mesh_pt->nboundary();
     for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
      {
       // Which boundaries (could be more than one) is it on?
       if (master_nod_pt->is_on_boundary(i_bnd))
        {
         Flat_packed_unsigneds.push_back(1);
        }
       else
        {
         Flat_packed_unsigneds.push_back(0);
        }
      }
    }
   else
    {
     // Not on any boundary
     Flat_packed_unsigneds.push_back(0);
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
     Flat_packed_unsigneds.push_back(update_id);

     // Get reference values at default...
     unsigned n_ref_val=alg_nod_pt->nref_value();
     Flat_packed_unsigneds.push_back(n_ref_val);
     for (unsigned i_ref_val=0;i_ref_val<n_ref_val;i_ref_val++)
      {
       Flat_packed_doubles.push_back(alg_nod_pt->ref_value(i_ref_val));
      }

     // Access geometric objects at default...
     unsigned n_geom_obj=alg_nod_pt->ngeom_object();
     Flat_packed_unsigneds.push_back(n_geom_obj);
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
       Flat_packed_unsigneds.push_back(found_geom_object);
      }
    } // end AlgebraicNode check

   // Is it a MacroElementNodeUpdateNode?
   MacroElementNodeUpdateNode* macro_nod_pt=
    dynamic_cast<MacroElementNodeUpdateNode*>(master_nod_pt);
   if (macro_nod_pt!=0)
    {
     // Loop over current external haloed elements - has the element which
     // controls the node update for this node been added yet?
     GeneralisedElement* macro_node_update_el_pt=
      macro_nod_pt->node_update_element_pt();

     unsigned n_ext_haloed_el=external_mesh_pt->
      nexternal_haloed_element(iproc);
     unsigned external_haloed_el_index;
     external_haloed_el_index=external_mesh_pt->
      add_external_haloed_element_pt(iproc,macro_node_update_el_pt);

     // If it wasn't already added, we need to create a halo copy
     if (external_haloed_el_index==n_ext_haloed_el)
      {
       Flat_packed_unsigneds.push_back(1);
       
       //Cast to a finite elemnet
       FiniteElement* macro_node_update_finite_el_pt 
        = dynamic_cast<FiniteElement*>(macro_node_update_el_pt);
       
       // We're using macro elements to update...
       MacroElementNodeUpdateMesh* macro_mesh_pt=
        dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
       if (macro_mesh_pt!=0)
        {
         Flat_packed_unsigneds.push_back(1);
         
         // Need to send the macro element number in the mesh across
         MacroElement* macro_el_pt= macro_node_update_finite_el_pt
          ->macro_elem_pt();
         unsigned macro_el_num=macro_el_pt->macro_element_number();
         Flat_packed_unsigneds.push_back(macro_el_num);

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
             Flat_packed_doubles.push_back(q_el_pt->s_macro_ll(i_dim));
             Flat_packed_doubles.push_back(q_el_pt->s_macro_ur(i_dim));
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
         Flat_packed_unsigneds.push_back(0);
        }

       // This element needs to be fully functioning on the other
       // process, so send all the information required to create it
       unsigned n_node=macro_node_update_finite_el_pt->nnode();
       for (unsigned j=0;j<n_node;j++)
        {
         Node* new_nod_pt=macro_node_update_finite_el_pt->node_pt(j);
         add_external_haloed_node_to_storage(iproc,new_nod_pt,
                                             problem_pt,
                                             external_mesh_pt,
                                             n_cont_inter_values);
        }
      }
     else // The external haloed element already exists
      {
       Flat_packed_unsigneds.push_back(0);
       Flat_packed_unsigneds.push_back(external_haloed_el_index);
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
         Flat_packed_doubles.push_back(solid_nod_pt->variable_position_pt()->
                          value(t,i_val));
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
       Flat_packed_doubles.push_back(master_nod_pt->value(t,i_val));
      }
    }

   // Now do positions
   for (unsigned idim=0;idim<n_dim;idim++)
    {
     for (unsigned t=0;t<n_prev;t++)
      {
       Flat_packed_doubles.push_back(master_nod_pt->x(t,idim));
      }
    }

  }







//=======start of add_external_halo_node_helper===========================
/// Helper functiono to add external halo node that is not a master
//========================================================================
  void add_external_halo_node_helper
  (Node* &new_nod_pt, Mesh* const &external_mesh_pt, unsigned& loc_p,
   unsigned& node_index, FiniteElement* const &new_el_pt, 
   int& n_cont_inter_values,
   Problem* problem_pt)
 {
  // Given the node and the external mesh, and received information
  // about them from process loc_p, construct them on the current process
  if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]==1)
   {
    // Increment counter
    Counter_for_flat_packed_unsigneds++;

    // Construct a new node based upon sent information
    construct_new_external_halo_node_helper(new_nod_pt,loc_p,
                                            node_index,new_el_pt,
                                            external_mesh_pt,problem_pt);
   }
  else
   {
    // Increment counter (node already exists)
    Counter_for_flat_packed_unsigneds++;

    // Copy node from received location
    new_nod_pt=external_mesh_pt->external_halo_node_pt
     (loc_p,Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]);
    new_el_pt->node_pt(node_index)=new_nod_pt;
    // Increment counter
    Counter_for_flat_packed_unsigneds++;
   }
 }




//========start of construct_new_external_halo_node_helper=================
/// Helper function which constructs a new external halo node (on new element)
/// with the required information sent from the haloed process
//========================================================================
 void construct_new_external_halo_node_helper
 (Node* &new_nod_pt, unsigned& loc_p, unsigned& node_index,
  FiniteElement* const &new_el_pt, 
  Mesh* const &external_mesh_pt, Problem* problem_pt)
 {
  // The first entry indicates the number of values at this new Node
  // (which may be different across the same element e.g. Lagrange multipliers)
  unsigned n_val=Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds];
  Counter_for_flat_packed_unsigneds++;

  // Null TimeStepper for now
  TimeStepper* time_stepper_pt=0;
  // Default number of previous values to 1
  unsigned n_prev=1;

  // The next entry in Flat_packed_unsigneds indicates
  // if a timestepper is required for this halo node
  if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]==1)
   {
    Counter_for_flat_packed_unsigneds++;
    // Index
    time_stepper_pt=problem_pt->time_stepper_pt
     (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]);
    Counter_for_flat_packed_unsigneds++;
    // Check whether number of prev values is "sent" across
    n_prev=time_stepper_pt->ntstorage(); 
   }
  else
   {
    // No timestepper, increment counter
    Counter_for_flat_packed_unsigneds++;
   }

  // If this node was on a boundary then it needs to
  // be on the same boundary here
  if (Flat_packed_unsigneds[Counter_for_flat_packed_unsigneds]==1)
   {
    Counter_for_flat_packed_unsigneds++;

    // Construct a new boundary node
    if (time_stepper_pt!=0)
     {
      new_nod_pt=new_el_pt->construct_boundary_node
       (node_index,time_stepper_pt);
     }
    else
     {
      new_nod_pt=new_el_pt->construct_boundary_node(node_index);
     }

    // How many boundaries on the external mesh?
    unsigned n_bnd=external_mesh_pt->nboundary();
    for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
     {
      if (Flat_packed_unsigneds
          [Counter_for_flat_packed_unsigneds]==1)
       {
        // Add to current boundary; increment counter
        external_mesh_pt->add_boundary_node(i_bnd,
                                            new_nod_pt);
        Counter_for_flat_packed_unsigneds++;
       }
      else
       {
        // Not on this boundary; increment counter
        Counter_for_flat_packed_unsigneds++;
       }
     }
   }
  else
   {
    // Not on boundary, increment counter
    Counter_for_flat_packed_unsigneds++;

    // Construct an ordinary (non-boundary) node
    if (time_stepper_pt!=0)
     {
      new_nod_pt=new_el_pt->construct_node
       (node_index,time_stepper_pt);
     }
    else
     {
      new_nod_pt=new_el_pt->construct_node(node_index);
     }
   }

  // Node constructed: add to external halo nodes
  external_mesh_pt->add_external_halo_node_pt(loc_p,new_nod_pt);

  // Is the new constructed node Algebraic?
  AlgebraicNode* new_alg_nod_pt=dynamic_cast<AlgebraicNode*>
   (new_nod_pt);

  // If it is algebraic, its node update functions will
  // not yet have been set up properly
  if (new_alg_nod_pt!=0)
   {
    // The AlgebraicMesh is the external mesh
    AlgebraicMesh* alg_mesh_pt=dynamic_cast<AlgebraicMesh*>
     (external_mesh_pt);

    /// The first entry of All_alg_nodal_info contains
    /// the default node update id
    /// e.g. for the quarter circle there are 
    /// "Upper_left_box", "Lower right box" etc...
    unsigned update_id=Flat_packed_unsigneds
     [Counter_for_flat_packed_unsigneds];
    Counter_for_flat_packed_unsigneds++;

    Vector<double> ref_value;

    // The size of this vector is in the next entry
    // of All_alg_nodal_info
    unsigned n_ref_val=Flat_packed_unsigneds
     [Counter_for_flat_packed_unsigneds];
    Counter_for_flat_packed_unsigneds++;

    // The reference values themselves are in
    // All_alg_ref_value
    ref_value.resize(n_ref_val);
    for (unsigned i_ref=0;i_ref<n_ref_val;i_ref++)
     {
      ref_value[i_ref]=Flat_packed_doubles
       [Counter_for_flat_packed_doubles];
      Counter_for_flat_packed_doubles++;
     }

    Vector<GeomObject*> geom_object_pt;
    /// again we need the size of this vector as it varies
    /// between meshes; we also need some indication
    /// as to which geometric object should be used...

    // The size of this vector is in the next entry
    // of All_alg_nodal_info
    unsigned n_geom_obj=Flat_packed_unsigneds
     [Counter_for_flat_packed_unsigneds];
    Counter_for_flat_packed_unsigneds++;

    // The remaining indices are in the rest of 
    // All_alg_nodal_info
    geom_object_pt.resize(n_geom_obj);
    for (unsigned i_geom=0;i_geom<n_geom_obj;i_geom++)
     {
      unsigned geom_index=Flat_packed_unsigneds
       [Counter_for_flat_packed_unsigneds];
      Counter_for_flat_packed_unsigneds++;
      // This index indicates which of the AlgebraicMesh's
      // stored geometric objects should be used
      // (0 is a null pointer; everything else should have
      //  been filled in by the specific Mesh).  If it
      // hasn't been filled in then the update_node_update
      // call should fix it
      geom_object_pt[i_geom]=alg_mesh_pt->
       geom_object_list_pt(geom_index);
     }

    /// For the received update_id, ref_value, geom_object
    /// call add_node_update_info
    new_alg_nod_pt->add_node_update_info
     (update_id,alg_mesh_pt,geom_object_pt,ref_value);

    /// Now call update_node_update
    alg_mesh_pt->update_node_update(new_alg_nod_pt);
   }

  // Is the node a MacroElementNodeUpdateNode?
  MacroElementNodeUpdateNode* macro_nod_pt=
   dynamic_cast<MacroElementNodeUpdateNode*>(new_nod_pt);

  if (macro_nod_pt!=0)
   {
    // Need to call set_node_update_info; this requires
    // a Vector<GeomObject*> (taken from the mesh)
    Vector<GeomObject*> geom_object_vector_pt;

    // Access the required geom objects from the
    // MacroElementNodeUpdateMesh
    MacroElementNodeUpdateMesh* macro_mesh_pt=
     dynamic_cast<MacroElementNodeUpdateMesh*>
     (external_mesh_pt);
    geom_object_vector_pt=
     macro_mesh_pt->geom_object_vector_pt();

    // Get local coordinate of node in new element
    Vector<double> s_in_macro_node_update_element;
    new_el_pt->local_coordinate_of_node
     (node_index,s_in_macro_node_update_element);

    // Set node update info for this node
    macro_nod_pt->set_node_update_info
     (new_el_pt,s_in_macro_node_update_element,
      geom_object_vector_pt);
   }

  // Is the new node a SolidNode? 
  SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(new_nod_pt);
  if (solid_nod_pt!=0)
   {
    unsigned n_solid_val=solid_nod_pt->variable_position_pt()->nvalue();
    for (unsigned i_val=0;i_val<n_solid_val;i_val++)
     {
      for (unsigned t=0;t<n_prev;t++)
       {
        solid_nod_pt->variable_position_pt()->
         set_value(t,i_val,
                   Flat_packed_doubles[Counter_for_flat_packed_doubles]);
        Counter_for_flat_packed_doubles++;
       }
     }
   }

  // If there are additional values, resize the node
  unsigned n_new_val=new_nod_pt->nvalue();
  if (n_val>n_new_val)
   {
    new_nod_pt->resize(n_val);
   }

  // Get copied history values
  //  unsigned n_val=new_nod_pt->nvalue();
  for (unsigned i_val=0;i_val<n_val;i_val++)
   {
    for (unsigned t=0;t<n_prev;t++)
     {
      new_nod_pt->set_value(t,i_val,Flat_packed_doubles
                            [Counter_for_flat_packed_doubles]);
      Counter_for_flat_packed_doubles++;
     }
   }

  // Get copied history values for positions
  unsigned n_dim=new_nod_pt->ndim();
  for (unsigned idim=0;idim<n_dim;idim++)
   {
    for (unsigned t=0;t<n_prev;t++)
     {
      // Copy to coordinate
      new_nod_pt->x(t,idim)=Flat_packed_doubles
       [Counter_for_flat_packed_doubles];
      Counter_for_flat_packed_doubles++;
     }
   }
 }

 
 //=====================================================================
 /// Locate zeta for current set of missing coordinates
 //=====================================================================
  void locate_zeta_for_missing_coordinates
  (int& iproc, Mesh* const &external_mesh_pt, Problem* problem_pt,
   MeshAsGeomObject* &mesh_geom_obj_pt)
  {
   // Storage for number of processors, current process and communicator
   OomphCommunicator* comm_pt=problem_pt->communicator_pt();
   int n_proc=comm_pt->nproc();
   int my_rank=comm_pt->my_rank();

   // Resize vectors containing data to be setn
   Flat_packed_doubles.resize(0);
   Flat_packed_unsigneds.resize(0);
   Flat_packed_located_coordinates.resize(0);

   // Flush storage for zetas not found locally (when
   // processing the zeta coordinates received from "previous"
   // processor)
   Flat_packed_zetas_not_found_locally.resize(0);
   
   // Number of zeta tuples to be dealt with
   unsigned n_zeta=Received_flat_packed_zetas_to_be_found.size()/Dim;
   
   // Create storage for the processor id (plus one) on which
   // the zetas stored in Flat_packed_zetas_not_found_locally[...]
   // were located 
   Proc_id_plus_one_of_external_element.resize(n_zeta);

   // Create storage for the status of the (external halo) element associated
   // the zetas stored in Flat_packed_zetas_not_found_locally[...].
   // It either hasn't been found, already exists on the processor
   // that needs it, or needs to be newly created.
   Located_element_status.resize(n_zeta);

   // Counter for flat-packed array of external zeta coordinates
   unsigned count=0;

   // Loop over the zeta tuples that we received from elsewhere and
   // are trying to find here
   for (unsigned i=0;i<n_zeta;i++)
    {
     // Storage for global coordinates to be located
     Vector<double> x_global(Dim);
     
     // Loop to fill in coordinates
     for (unsigned ii=0;ii<Dim;ii++)
      {
       x_global[ii]=Received_flat_packed_zetas_to_be_found[count];
       count++;
      }
     
     // Perform locate_zeta for these coordinates
     GeomObject *sub_geom_obj_pt;
     Vector<double> ss(Dim);
     bool called_within_spiral=true;
     mesh_geom_obj_pt->spiraling_locate_zeta(x_global,sub_geom_obj_pt,ss,
                                             called_within_spiral);
     
     // Did the locate method work?
     if (sub_geom_obj_pt!=0)
      {
       // Get the source element - bulk or not?
       GeneralisedElement *source_el_pt=0;
       if (!Use_bulk_element_as_external)
        {
         source_el_pt=dynamic_cast<FiniteElement*>(sub_geom_obj_pt);
        }
       else
        {
         FaceElement *face_el_pt=dynamic_cast<FaceElement*>(sub_geom_obj_pt);
         source_el_pt=dynamic_cast<FiniteElement*>(face_el_pt->
                                                   bulk_element_pt());
        }
        
       // Check if the returned element is halo
       if (!source_el_pt->is_halo())
        {
         // The correct non-halo element has been located; this will become
         // an external haloed element on the current process, and an
         // external halo copy needs to be created on the current process
         // minus wherever we are in the "ring-loop"
         int halo_copy_proc=my_rank-iproc;
         // If iproc is bigger than my_rank then we've "gone through" nproc-1
         if (my_rank<iproc) { halo_copy_proc=n_proc+halo_copy_proc; }
          
         // So, we found zeta on the current processor
         Proc_id_plus_one_of_external_element[i]=my_rank+1;
          
         // This source element is an external halo on process halo_copy_proc
         // but it should only be added to the storage if it hasn't
         // been added already, and this information also needs to be
         // communicated over to the other process
          
         unsigned n_extern_haloed=external_mesh_pt->
          nexternal_haloed_element(halo_copy_proc);
         unsigned external_haloed_el_index=
          external_mesh_pt->add_external_haloed_element_pt(halo_copy_proc,
                                                           source_el_pt);
          
         // If it was added to the storage then the returned index
         // will be the same as the (old) size of the storage
         if (external_haloed_el_index==n_extern_haloed)
          {
           // Set Located_element_status to say it 
           // should be newly created
           Located_element_status[i]=New;
            
           // How many continuously interpolated values are there?
           int n_cont_inter_values=-1;
           if (dynamic_cast<RefineableElement*>(source_el_pt)!=0)
            {
             n_cont_inter_values=dynamic_cast<RefineableElement*>
              (source_el_pt)->ncont_interpolated_values();
            }
            
           // Since it is (externally) haloed from the current process,
           // the info required to create a new element in the equivalent
           // external halo layer on process halo_copy_proc needs to be 
           // sent there
            
           // If we're using macro elements to update...
           MacroElementNodeUpdateMesh* macro_mesh_pt=
            dynamic_cast<MacroElementNodeUpdateMesh*>(external_mesh_pt);
           if (macro_mesh_pt!=0)
            {
             Flat_packed_unsigneds.push_back(1);
              
             //Cast to finite element... this must work because it's
             //a macroelement no update mesh
             FiniteElement* source_finite_el_pt 
              = dynamic_cast<FiniteElement*>(source_el_pt);
              
             MacroElement* macro_el_pt=source_finite_el_pt->macro_elem_pt();
             // Send the macro element number across
             unsigned macro_el_num=macro_el_pt->macro_element_number();
             Flat_packed_unsigneds.push_back(macro_el_num);
              
             // we need to send
             // the lower left and upper right coordinates of the macro
             QElementBase* q_el_pt=dynamic_cast<QElementBase*>(source_el_pt);
             if (q_el_pt!=0)
              {
               // The macro element needs to be set first before
               // its lower left and upper right coordinates can be accessed
               // Now send the lower left and upper right coordinates
               unsigned el_dim=q_el_pt->dim();
               for (unsigned i_dim=0;i_dim<el_dim;i_dim++)
                {
                 Flat_packed_doubles.push_back(q_el_pt->s_macro_ll(i_dim));
                 Flat_packed_doubles.push_back(q_el_pt->s_macro_ur(i_dim));
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
                 "Multi_domain_functions::locate_zeta_for_missing_coordinates()",
                 OOMPH_EXCEPTION_LOCATION);
              }
              
            }
           else // Not using macro elements to update
            {
             Flat_packed_unsigneds.push_back(0);
            }
            
            
           //Cast to finite element... this must work because it's
           //a macroelement no update mesh
           FiniteElement* source_finite_el_pt 
            = dynamic_cast<FiniteElement*>(source_el_pt);
#ifdef PARANOID
           if(source_finite_el_pt==0)
            {
             throw 
              OomphLibError(
               "Unable to cast source function to finite element\n",
               "Multi_domain_functions::locate_zeta_for_missing_coordinates()",
               OOMPH_EXCEPTION_LOCATION);
            }
#endif
            
            
           // Loop over the nodes of the new source element
           unsigned n_node=source_finite_el_pt->nnode();
           for (unsigned j=0;j<n_node;j++)
            {
             Node* nod_pt=source_finite_el_pt->node_pt(j);
              
             // Add the node to the storage; this routine
             // also takes care of any master nodes if the
             // node is hanging
             add_external_haloed_node_to_storage(halo_copy_proc,nod_pt,
                                                 problem_pt,
                                                 external_mesh_pt,
                                                 n_cont_inter_values);
            }            
          }
         else // it has already been added, so tell the other process
          {
           // Set Located_element_status to indicate an element has 
           // already been added
           Located_element_status[i]=Exists;
           Flat_packed_unsigneds.push_back(external_haloed_el_index);
          }
          
         // The coordinates returned by locate_zeta are also needed
         // in the setup of the source elements on the other process
         if (!Use_bulk_element_as_external)
          {
           for (unsigned ii=0;ii<Dim;ii++)
            {
             Flat_packed_located_coordinates.push_back(ss[ii]);
            }
          }
         else // translate the coordinates to the bulk element
          {
           // The translation is from Lagrangian to Eulerian
           FaceElement *face_el_pt=
            dynamic_cast<FaceElement*>(sub_geom_obj_pt);
           //Get the dimension of the BulkElement
           unsigned bulk_el_dim = 
            dynamic_cast<FiniteElement*>(source_el_pt)->dim();
           Vector<double> s_trans(bulk_el_dim);
           face_el_pt->get_local_coordinate_in_bulk(ss,s_trans);
           for (unsigned ii=0;ii<bulk_el_dim;ii++)
            {
             Flat_packed_located_coordinates.push_back(s_trans[ii]);
            }
          }
        }
       else // halo, so search again until non-halo equivalent is located
        {
         // Add required information to arrays (as below)
         for (unsigned ii=0;ii<Dim;ii++)
          {
           Flat_packed_zetas_not_found_locally.push_back(x_global[ii]);
          }
         // It wasn't found here
         Proc_id_plus_one_of_external_element[i]=0;

         // Set Located_element_status to indicate not found
         Located_element_status[i]=Not_found;
        }
      }
     else // not successful this time, so prepare for next process to try
      {
       // Add this global coordinate to the LOCAL zeta array
       for (unsigned ii=0;ii<Dim;ii++)
        {
         Flat_packed_zetas_not_found_locally.push_back(x_global[ii]);
        }
       // It wasn't found here
       Proc_id_plus_one_of_external_element[i]=0;

       // Set Located_element_status to indicate not found
       Located_element_status[i]=Not_found;
      }
    }
     
  }

#endif



 //=====================================================================
 /// locate zeta for current set of "local" coordinates
 //=====================================================================
  void locate_zeta_for_local_coordinates
  (Mesh* const &mesh_pt, Mesh* const &external_mesh_pt,
   MeshAsGeomObject* &mesh_geom_obj_pt,
   const unsigned& interaction_index)
  {
   // Number of local elements
   unsigned n_element=mesh_pt->nelement();

   // Flush storage for zetas not found locally 
   Flat_packed_zetas_not_found_locally.resize(0);

   // Loop over this processor's elements
   for (unsigned e=0;e<n_element;e++)
    {
     ElementWithExternalElement *el_pt=
      dynamic_cast<ElementWithExternalElement*>(mesh_pt->element_pt(e));
#ifdef OOMPH_HAS_MPI
     // Only visit non-halo elements
     if (!el_pt->is_halo())
#endif
      {
       // Find number of Gauss points and element dimension
       unsigned n_intpt=el_pt->integral_pt()->nweight();
       unsigned el_dim=el_pt->dim();
       // Set storage for local and global coordinates
       Vector<double> s_local(el_dim);
       Vector<double> x_global(el_dim);

       // Loop over integration points
       for (unsigned ipt=0;ipt<n_intpt;ipt++)
        {         
         // Has this integration point been done yet?
         if (External_element_located[e][ipt]==0)
          {
           // Get local coordinates
           for (unsigned i=0;i<el_dim;i++)
            {
             s_local[i]=el_pt->integral_pt()->knot(ipt,i);
            }
           // Interpolate to global coordinates
           el_pt->interpolated_zeta(s_local,x_global);

           // Storage for geometric object and its local coordinates
           GeomObject* sub_geom_obj_pt=0;
           Vector<double> s_ext(el_dim);

           // Perform locate_zeta locally for this coordinate
           bool called_within_spiral=true;
           mesh_geom_obj_pt->spiraling_locate_zeta(x_global,
                                                   sub_geom_obj_pt,s_ext,
                                                   called_within_spiral);

           // Has the required element been located?
           if (sub_geom_obj_pt!=0)
            {
             // The required element has been located
             // The located coordinates have the same dimension as the bulk
             GeneralisedElement* source_el_pt;
             Vector<double> s_source(el_dim);

             // Is the bulk element the actual external element?
             if (!Use_bulk_element_as_external)
              {
               // Use the object directly (it must be a finite element)
               source_el_pt=dynamic_cast<FiniteElement*>(sub_geom_obj_pt);
               s_source=s_ext;
              }
             else
              {
               // Cast to a FaceElement and use the bulk element
               FaceElement* face_el_pt=
                dynamic_cast<FaceElement*>(sub_geom_obj_pt);
               source_el_pt=face_el_pt->bulk_element_pt();

               //Need to resize the located coordinates to have the same
               //dimension as the bulk element
               s_source.resize(dynamic_cast<FiniteElement*>
                               (source_el_pt)->dim());

               // Translate the returned local coords into the bulk element
               face_el_pt->get_local_coordinate_in_bulk(s_ext,s_source);
              }

             // Check if it's a halo; if it is then the non-halo equivalent
             // needs to be located from another processor
#ifdef OOMPH_HAS_MPI
             if (!source_el_pt->is_halo())
#endif
              {
               //Need to cast to a FiniteElement
               FiniteElement* source_finite_el_pt = 
                dynamic_cast<FiniteElement*>(source_el_pt);

               // Set the external element pointer and local coordinates
               el_pt->external_element_pt(interaction_index,ipt)
                = source_finite_el_pt;
               el_pt->external_element_local_coord(interaction_index,ipt)
                =s_source;

               // Set the lookup array to 1/true 
               External_element_located[e][ipt]=1;
              }
#ifdef OOMPH_HAS_MPI
             else // elements can only be halo if MPI is turned on
              {
               // Add required information to arrays
               for (unsigned i=0;i<el_dim;i++)
                {
                 Flat_packed_zetas_not_found_locally.push_back(x_global[i]);
                }
              }
#endif
            }
           else
            {
             // If it has failed then add the required information to the
             // arrays which need to be sent to the other processors so that
             // they can perform the locate_zeta
             
             // Add this global coordinate to the LOCAL zeta array
             for (unsigned i=0;i<el_dim;i++)
              {
               Flat_packed_zetas_not_found_locally.push_back(x_global[i]);
              }
            }
          }
        } // end loop over integration points
      }
    } // end loop over local elements
   
  }




//=====================================================================
/// Helper function that computes the dimension of the elements within
/// each of the specified meshes (and checks they are the same).
/// Stores result in Dim.
//=====================================================================
 void get_dim_helper(Problem* problem_pt, Mesh* const &mesh_pt, 
                     Mesh* const &external_mesh_pt)
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

   // Set dimension
   Dim=mesh_dim;
 }


 //=================================================================
 /// Helper function that clears all the information used
 /// during the external storage creation
 //=================================================================
 void clean_up()
  {
   Flat_packed_zetas_not_found_locally.clear();
   Received_flat_packed_zetas_to_be_found.clear();
   Proc_id_plus_one_of_external_element.clear();
   Located_element_status.clear();
   Flat_packed_located_coordinates.clear();
   Flat_packed_doubles.clear();
   Flat_packed_unsigneds.clear();
   External_element_located.clear();
  }


}

}

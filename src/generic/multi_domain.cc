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
//Multi-domain functions which act on more than one mesh and set up the
//storage and interaction between the two

//oomph-lib header
#include "multi_domain.h"
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

  // List of all Vectors used by the external storage routines 
  
  /// \short Vector of (local) coordinates at integration points of
  /// elements on current processor
  Vector<double> Zetas;

  /// \short Vector of the dimension of the element on current processor
  Vector<unsigned> Zeta_dim;

  /// \short Vector to indicate which processor a coordinate is located on
  Vector<int> Found_zeta;

  /// \short Vector of any elements found locally by current processor
  Vector<FiniteElement*> Found_element;

  /// \short Vector of local coordinates within any elements found
  /// locally by current processor
  Vector<Vector<double> > Found_ss;

  /// \short Global vector which combines all Found_zeta vectors on
  /// every processor
  Vector<int> All_found_zeta;

  /// \short Vector to indicate (on another process) whether a
  /// located element should be newly created (2), already exists (1), or
  /// is not on the current process at all (0)
  Vector<unsigned> Located_element;

  /// \short Vector of the local coordinates for each entry in Located_element
  Vector<double> Located_coord;

  /// \short Vector of doubles to be sent from another processor
  Vector<double> Double_values;

  /// \short Vector of unsigneds to be sent from another processor
  Vector<unsigned> Unsigned_values;

  // Counters for arrays used in the external storage routines
  unsigned Count_double_values;
  unsigned Count_unsigned_values;
  unsigned Count_located_coord;
  unsigned Count_found_elements;
  unsigned Count_zeta_dim;
  unsigned Count_zetas;

  // The following Vector of Vectors (of size Nproc) store each required 
  // set of the above information from each process required to create 
  // external (halo) elements and nodes
  Vector<Vector<unsigned> > All_located_zetas;
  Vector<Vector<double> > All_located_coord;
  Vector<Vector<double> > All_double_values;
  Vector<Vector<unsigned> > All_unsigned_values;

  // These Vectors help to count through the above Vectors of Vectors
  Vector<unsigned> All_count_double_values;
  Vector<unsigned> All_count_located_coord;
  Vector<unsigned> All_count_unsigned_values;


  /// Default binning parameters

  /// \short Bool to tell the MeshAsGeomObject to setup bins again
  /// (should only be used in conjunction with a change from default N*_bin;
  ///  the minimum and maximum of the Mesh used must also be specified)
  bool Setup_bins_again=false;

  /// \short Number of bins along each dimension in binning method in
  /// set_external_storage(). Default value of 10.
  unsigned Nx_bin=10;
  unsigned Ny_bin=10;
  unsigned Nz_bin=10;

  /// \short Minimum and maximum coordinates for
  /// each dimension of the external mesh used to "create" the bins in
  /// set_external_storage(). No defaults; set by user if they want to
  /// (otherwise the MeshAsGeomObject calculates these values based
  ///  upon the mesh itself; see MeshAsGeomObject::get_max_and_min_coords(...))
  double X_min;
  double X_max;
  double Y_min;
  double Y_max;
  double Z_min;
  double Z_max;

  /// \short Percentage offset to add to each extreme of the bin structure.
  /// Default value of 0.05.
  double Percentage_offset=0.05;

  /// \short Boolean to indicate when to use the bulk element as the
  /// external element.  Defaults to false, you must have set up FaceElements
  /// properly first in order for it to work
  bool Use_bulk_element_as_external=false;

  /// \short Boolean to indicate whether to doc timings or not.
  bool Doc_timings=false;

  /// \short Boolean to indicate whether to output info during
  ///        set_external_storage routines
  bool Shut_up=true;

  // Functions for location method in multi-domain problems 

#ifdef OOMPH_HAS_MPI

  //======================================================================
  /// Function which removes duplicate data that exist because
  /// they have been distinctly created by communications from different
  /// processors, whereas the data are in fact the same
  //======================================================================
  void remove_duplicate_data(Mesh* const &mesh_pt)
  {
   // Each individual container of external halo elements has unique
   // nodes/equation numbers, but there may be some duplication between
   // two or more different containers; the following code checks for this
   // and removes the duplication by pointing the a data point with an already
   // existing eqn number to the original data point which had the eqn no.

   // Storage for existing global equation numbers
   Vector<std::pair<Vector<int>,Node*> > existing_global_eqn_numbers(0);

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
   for (int iproc=MPI_Helpers::Nproc-1;iproc>=0;iproc--)
    {
     // Don't have external halo elements with yourself!
     if (iproc!=MPI_Helpers::My_rank)
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

#endif


//======================================================================
/// Function which adds external data to current elements from the
/// external elements at each integration point of the current element
/// for the specified interaction index
//======================================================================
  void add_external_data_from_source_elements
  (Mesh* const &mesh_pt,const unsigned& interaction_index)
  {
   unsigned n_element=mesh_pt->nelement();
   // Loop over elements
   for (unsigned e=0;e<n_element;e++)
    {
     ElementWithExternalElement* el_pt=
      dynamic_cast<ElementWithExternalElement*>(mesh_pt->element_pt(e));
     // Only bother with non-halo elements
#ifdef OOMPH_HAS_MPI
     if (!el_pt->is_halo())
#endif
      {
       unsigned n_intpt=el_pt->integral_pt()->nweight();
       // Loop over integration points
       for (unsigned i=0;i<n_intpt;i++)
        {
         FiniteElement* source_el_pt=dynamic_cast<FiniteElement*>
          (el_pt->external_element_pt(interaction_index,i));
         if (source_el_pt==0)
          {
           std::ostringstream error_message;
           error_message << "Source element pointer not set for element " << e
                         << " at integration point " << i << ".\n";
           throw OomphLibError
            (error_message.str(),
             "Multi_domain_functions::add_external_data_from_source_elements",
             OOMPH_EXCEPTION_LOCATION);
          }
         unsigned n_node=source_el_pt->nnode();
         for (unsigned j=0;j<n_node;j++)
          {
           // Add the node as external data
           el_pt->add_external_data(source_el_pt->node_pt(j));
           // If this "source" node is hanging (in any variable) then its
           // master nodes also need to be added as external data
           if (dynamic_cast<RefineableElement*>(source_el_pt)!=0)
            {
             int n_cont_inter_values=dynamic_cast<RefineableElement*>
              (source_el_pt)->ncont_interpolated_values();
             for (int i_cont=-1;i_cont<n_cont_inter_values;i_cont++)
              {
               if (source_el_pt->node_pt(j)->is_hanging(i_cont))
                {
                 HangInfo* hang_pt=source_el_pt->node_pt(j)
                  ->hanging_pt(i_cont);
                 unsigned n_master=hang_pt->nmaster();
                 // Loop over master nodes
                 for (unsigned m=0; m<n_master; m++)
                  {
                   Node* master_nod_pt=hang_pt->master_node_pt(m);
                   // Add this node as external data
                   el_pt->add_external_data(master_nod_pt);
                  }
                }
              }
            }
          }
         // Add internal data points as external data
         unsigned n_int_data=source_el_pt->ninternal_data();
         for (unsigned i_int=0;i_int<n_int_data;i_int++)
          {
           el_pt->add_external_data(source_el_pt->internal_data_pt(i_int));
          }       

        } // end loop over integration points
      }
    }
  }



//========================================================================
/// Broadcast the zeta coordinates from the current process (iproc) to 
/// all other processes
//========================================================================
  void broadcast_local_zeta(int& iproc, Mesh* const &mesh_pt)
  {
#ifdef OOMPH_HAS_MPI
   MPI_Status status;
#endif
   // Resize arrays
   Zetas.resize(0);
   Zeta_dim.resize(0);
   // Initialise counters
   Count_zetas=0;
   Count_zeta_dim=0;

   // Number of elements
   unsigned n_element;
   // Current process?
   if (iproc==MPI_Helpers::My_rank)
    {
     // Find the "zetas" on this process and send them to all other processes
     n_element=mesh_pt->nelement();
     // If there are no elements then there's no need to do anything on any
     // process; "broadcast" the number of elements to every other process
     // "Broadcast" the number of elements to expect to every other process
#ifdef OOMPH_HAS_MPI
     for (int d=0;d<MPI_Helpers::Nproc;d++)
      {
       if (d!=iproc)
        {
         MPI_Send(&n_element,1,MPI_INT,d,0,MPI_COMM_WORLD);
        }
      }
#endif
     // Only need to do any work if this process has any elements
     if (n_element>0)
      {
       // Prepare vectors to send to other processes
       for (unsigned e=0;e<n_element;e++)
        {
         ElementWithExternalElement *el_pt=
          dynamic_cast<ElementWithExternalElement*>(mesh_pt->element_pt(e));
         // Only need to work on non-halo elements
#ifdef OOMPH_HAS_MPI
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
             // Get local coordinates
             for (unsigned i=0;i<el_dim;i++)
              {
               s_local[i]=el_pt->integral_pt()->knot(ipt,i);
              }
             // Interpolate to global coordinates
             el_pt->interpolated_zeta(s_local,x_global);
             // Add this global coordinate to the zeta array
             for (unsigned i=0;i<el_dim;i++)
              {
               Zetas.push_back(x_global[i]);
               Count_zetas++;
              }
             // Add the element dimension to the Zeta_dim array
             Zeta_dim.push_back(el_dim);
             Count_zeta_dim++;
            }
          }
        }

#ifdef OOMPH_HAS_MPI
       // Send the Zetas array to all other processes
       for (int d=0;d<MPI_Helpers::Nproc;d++)
        {
         // Don't send to yourself
         if (d!=iproc)
          {
           MPI_Send(&Count_zeta_dim,1,MPI_INT,d,1,MPI_COMM_WORLD);
           if (Count_zeta_dim!=0)
            {
             MPI_Send(&Zeta_dim[0],Count_zeta_dim,MPI_INT,
                      d,3,MPI_COMM_WORLD);
            }
           MPI_Send(&Count_zetas,1,MPI_INT,d,4,MPI_COMM_WORLD);
           if (Count_zetas!=0)
            {           
             MPI_Send(&Zetas[0],Count_zetas,MPI_DOUBLE,d,5,MPI_COMM_WORLD);
            }
          }
        }
#endif

      }  
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     MPI_Recv(&n_element,1,MPI_INT,iproc,0,MPI_COMM_WORLD,&status);

     if (n_element>0)
      {
       // Receive the zeta array from the loop process
       MPI_Recv(&Count_zeta_dim,1,MPI_INT,iproc,1,MPI_COMM_WORLD,&status);
       if (Count_zeta_dim!=0)
        {
         Zeta_dim.resize(Count_zeta_dim);
         MPI_Recv(&Zeta_dim[0],Count_zeta_dim,MPI_INT,iproc,3,
                  MPI_COMM_WORLD,&status);
        }

       MPI_Recv(&Count_zetas,1,MPI_INT,iproc,4,MPI_COMM_WORLD,&status);
       if (Count_zetas!=0)
        {       
         Zetas.resize(Count_zetas);
         MPI_Recv(&Zetas[0],Count_zetas,MPI_DOUBLE,iproc,5,
                  MPI_COMM_WORLD,&status);
        }
      }
#endif
    }
  }



#ifdef OOMPH_HAS_MPI
//========start of prepare_and_send_external_element_info=================
/// Prepares external element information from non-loop processor to
/// send to loop processor for the creation of external halo elements
//========================================================================
  void prepare_and_send_external_element_info
  (int& iproc, Mesh* const &external_mesh_pt, Problem* problem_pt)
  {
   MPI_Status status;

   // Initialise all the required arrays on the loop processor
   All_double_values.resize(MPI_Helpers::Nproc);
   All_located_coord.resize(MPI_Helpers::Nproc);
   All_located_zetas.resize(MPI_Helpers::Nproc);
   All_unsigned_values.resize(MPI_Helpers::Nproc);

   All_count_double_values.resize(MPI_Helpers::Nproc);
   All_count_located_coord.resize(MPI_Helpers::Nproc);
   All_count_unsigned_values.resize(MPI_Helpers::Nproc);
 
   // Are we on the "main" loop processor?
   if (iproc==MPI_Helpers::My_rank)
    {
     // Receive the sent information from the other processes
     // where locate_zeta succeeded and setup external halo structure
     
     // Loop over other processes
     for (int d=0;d<MPI_Helpers::Nproc;d++)
      {
       if (d!=iproc)
        {
         MPI_Recv(&Count_double_values,1,MPI_INT,d,1,MPI_COMM_WORLD,&status);
         All_count_double_values[d]=Count_double_values;

         if (Count_double_values!=0)
          {
           All_double_values[d].resize(Count_double_values);

           MPI_Recv(&All_double_values[d][0],Count_double_values,MPI_DOUBLE,d,
                    2,MPI_COMM_WORLD,&status);
          }

         if (Count_zeta_dim!=0)
          {
           All_located_zetas[d].resize(Count_zeta_dim);
           MPI_Recv(&All_located_zetas[d][0],Count_zeta_dim,MPI_INT,d,3,
                    MPI_COMM_WORLD,&status);
          }

         MPI_Recv(&Count_located_coord,1,MPI_INT,d,4,
                  MPI_COMM_WORLD,&status);

         All_count_located_coord[d]=Count_located_coord;
         if (Count_located_coord!=0)
          {
           All_located_coord[d].resize(Count_located_coord);
           MPI_Recv(&All_located_coord[d][0],Count_located_coord,MPI_DOUBLE,d,
                    5,MPI_COMM_WORLD,&status);
          }

         MPI_Recv(&Count_unsigned_values,1,MPI_INT,d,14,
                  MPI_COMM_WORLD,&status);

         All_count_unsigned_values[d]=Count_unsigned_values;
         if (Count_unsigned_values!=0)
          {
           All_unsigned_values[d].resize(Count_unsigned_values);
           MPI_Recv(&All_unsigned_values[d][0],Count_unsigned_values,MPI_INT,
                    d,15,MPI_COMM_WORLD,&status);
          }

        }
      }

    }
   else // iproc!=MPI_Helpers::My_rank i.e. not current process
    {
     // Now that locate_zeta ought to have worked for each integration
     // point on each (non-halo) element in iproc's mesh, communicate
     // the required element and node info back to process iproc
     // (Note that the order is preserved due to the found_zeta[..] array)

     MPI_Send(&Count_double_values,1,MPI_INT,iproc,1,MPI_COMM_WORLD);
     if (Count_double_values!=0)
      {
       MPI_Send(&Double_values[0],Count_double_values,MPI_DOUBLE,iproc,2,
                MPI_COMM_WORLD);
      }

     if (Count_zeta_dim!=0)
      {
       MPI_Send(&Located_element[0],Count_zeta_dim,MPI_INT,
                iproc,3,MPI_COMM_WORLD);
      }

     MPI_Send(&Count_located_coord,1,MPI_INT,iproc,4,MPI_COMM_WORLD);
     if (Count_located_coord!=0)
      {     
       MPI_Send(&Located_coord[0],Count_located_coord,MPI_DOUBLE,
                iproc,5,MPI_COMM_WORLD);
      }

     MPI_Send(&Count_unsigned_values,1,MPI_INT,iproc,14,MPI_COMM_WORLD);
     if (Count_unsigned_values!=0)
      {     
       MPI_Send(&Unsigned_values[0],Count_unsigned_values,MPI_INT,iproc,15,
                MPI_COMM_WORLD);
      }

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
                                                nod_pt,problem_pt,
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
   // Check whether node is currently an external haloed node
   bool is_an_external_haloed_node=false;
   unsigned external_haloed_node_index;
   unsigned n_ext_haloed_nod=external_mesh_pt->nexternal_haloed_node(iproc);
   for (unsigned k=0;k<n_ext_haloed_nod;k++)
    {
     if (nod_pt==
         external_mesh_pt->external_haloed_node_pt(iproc,k))
      {
       is_an_external_haloed_node=true;
       external_haloed_node_index=k;
       break;
      }
    }

   // If it's not then add it
   if (!is_an_external_haloed_node)
    {
     // Indicate that this node needs to be constructed on
     // the other process
     Unsigned_values.push_back(1);
     Count_unsigned_values++;
     // This gets all the required information for the specified
     // node and stores it into MPI-sendable information
     // so that a halo copy can be made on the receiving process
     get_required_nodal_information_helper(iproc,nod_pt,
                                           problem_pt,external_mesh_pt,
                                           n_cont_inter_values,f_el_pt);
     // Add it as an external haloed node
     external_mesh_pt->add_external_haloed_node_pt(iproc,nod_pt);
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
                                             Node* nod_pt, Problem* problem_pt,
                                             Mesh* const &external_mesh_pt,
                                             int& n_cont_inter_values)
  {
   // Check whether node is currently an external haloed node
   bool is_an_external_haloed_node=false;
   unsigned external_haloed_node_index;
   unsigned n_ext_haloed_nod=external_mesh_pt->nexternal_haloed_node(iproc);
   for (unsigned k=0;k<n_ext_haloed_nod;k++)
    {
     if (master_nod_pt==
         external_mesh_pt->external_haloed_node_pt(iproc,k))
      {
       is_an_external_haloed_node=true;
       external_haloed_node_index=k;
       break;
      }
    }

   // If it's not then add it
   if (!is_an_external_haloed_node)
    {
     // Add it as an external haloed node
     external_mesh_pt->add_external_haloed_node_pt(iproc,master_nod_pt);

     // Indicate that this node needs to be constructed on
     // the other process
     Unsigned_values.push_back(1);
     Count_unsigned_values++;

     // This gets all the required information for the specified
     // master node and stores it into MPI-sendable information
     // so that a halo copy can be made on the receiving process
     get_required_master_nodal_information_helper(iproc,master_nod_pt,
                                                  nod_pt,problem_pt,
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

   // Is it a MacroElementNodeUpdateNode?
   MacroElementNodeUpdateNode* macro_nod_pt=
    dynamic_cast<MacroElementNodeUpdateNode*>(nod_pt);
   if (macro_nod_pt!=0)
    {
     // It would appear as though everything required is taken care
     // of in the element anyway
    }

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
 (int& iproc, Node* master_nod_pt, Node* nod_pt, Problem* problem_pt,
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
     bool already_external_haloed_element=false;
     unsigned external_haloed_index;

     FiniteElement* macro_node_update_el_pt=
      macro_nod_pt->node_update_element_pt();
     unsigned n_ext_haloed_el=external_mesh_pt->
      nexternal_haloed_element(iproc);
     for (unsigned e=0;e<n_ext_haloed_el;e++)
      {
       if (macro_node_update_el_pt==
           (external_mesh_pt->external_haloed_element_pt(iproc,e)))
        {
         already_external_haloed_element=true;
         external_haloed_index=e;
         break;
        }
      }

     // If it's not already externally haloed then it needs to be
     if (!already_external_haloed_element)
      {
       Unsigned_values.push_back(1);
       Count_unsigned_values++;

       external_mesh_pt->add_external_haloed_element_pt
        (iproc,macro_node_update_el_pt);

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
     else // The external haloed node already exists
      {
       Unsigned_values.push_back(0);
       Count_unsigned_values++;

       Unsigned_values.push_back(external_haloed_index);
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

 /// Helper function that clears all the information used
 /// during the external storage creation
 void clean_up()
  {
   // Clear every vector associated with the external storage creation
   Zetas.clear();
   Zeta_dim.clear();
   Found_zeta.clear();
   Found_element.clear();
   Found_ss.clear();
   All_found_zeta.clear();

   Located_element.clear();
   Located_coord.clear();
   Double_values.clear();
   Unsigned_values.clear();

   // These Vector of Vectors (of size Nproc) stored each required set of 
   // information from each process to create external (halo) elements
   All_double_values.clear();
   All_located_coord.clear();
   All_located_zetas.clear();
   All_unsigned_values.clear();

   // Clear the corresponding counter arrays
   All_count_double_values.clear();
   All_count_located_coord.clear();
   All_count_unsigned_values.clear();
  }


}

}

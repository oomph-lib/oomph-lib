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
//Non-inline member functions for general mesh classes
#include<algorithm>

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

//oomph-lib headers
#include "oomph_utilities.h"
#include "mesh.h"
#include "problem.h"
#include "elastic_problems.h"
#include "refineable_mesh.h"

namespace oomph
{

//======================================================
/// The Steady Timestepper
//======================================================
Steady<0> Mesh::Default_TimeStepper;



//=======================================================================
/// Merge meshes.
/// Note: This simply merges the meshes' elements and nodes (ignoring
/// duplicates; no boundary information etc. is created). 
//=======================================================================
 void Mesh::merge_meshes(const Vector<Mesh*>& sub_mesh_pt)
 {
  // No boundary lookup scheme is set up for the combined mesh
   Lookup_for_elements_next_boundary_is_setup=false;

   //Number of submeshes
   unsigned nsub_mesh=sub_mesh_pt.size();
   
   // Initialise element, node and boundary counters for global mesh
   unsigned long n_element=0;
   unsigned long n_node=0;
   unsigned n_bound=0;
   
   // Loop over submeshes and get total number of elements, nodes and
   // boundaries
   for(unsigned imesh=0;imesh<nsub_mesh;imesh++)
    {
     n_element += sub_mesh_pt[imesh]->nelement();
     n_node += sub_mesh_pt[imesh]->nnode();
     n_bound += sub_mesh_pt[imesh]->nboundary();
    }  
   
   // Reserve storage for element and node pointers 
   Element_pt.clear();
   Element_pt.reserve(n_element);
   Node_pt.clear();
   Node_pt.reserve(n_node);

   //Resize vector of vectors of nodes
   Boundary_node_pt.clear();
   Boundary_node_pt.resize(n_bound);
   
   // Sets of pointers to elements and nodes (to exlude duplicates -- they
   // shouldn't occur anyway but if they do, they must only be added
   // once in the global mesh to avoid trouble in the timestepping)
   std::set<GeneralisedElement*> element_set_pt;
   std::set<Node*> node_set_pt;
   
   //Counter for total number of boundaries in all the submeshes
   unsigned ibound_global=0;   
   //Loop over the number of submeshes 
   for(unsigned imesh=0;imesh<nsub_mesh;imesh++)
    {
     //Loop over the elements of the submesh and add to vector
     //duplicates are ignored
     unsigned nel_before=0;
     unsigned long n_element=sub_mesh_pt[imesh]->nelement();
     for (unsigned long e=0;e<n_element;e++)
      {
       GeneralisedElement* el_pt=sub_mesh_pt[imesh]->element_pt(e);
       element_set_pt.insert(el_pt);
       // Was it a duplicate?
       unsigned nel_now=element_set_pt.size();
       if (nel_now==nel_before)
        {
         std::ostringstream warning_stream;
         warning_stream <<"WARNING: " << std::endl
                        <<"Element " << e << " in submesh " << imesh 
                        <<" is a duplicate \n and was ignored when assembling" 
                        <<" combined mesh." << std::endl;
         OomphLibWarning(warning_stream.str(),
                         "Mesh::Mesh(const Vector<Mesh*>&)",
                         OOMPH_EXCEPTION_LOCATION);
        }
       else
        {
         Element_pt.push_back(el_pt);
        }
       nel_before=nel_now;
      }
     
     //Loop over the nodes of the submesh and add to vector
     //duplicates are ignored
     unsigned nnod_before=0;
     unsigned long n_node=sub_mesh_pt[imesh]->nnode();
     for (unsigned long n=0;n<n_node;n++)
      {
       Node* nod_pt=sub_mesh_pt[imesh]->node_pt(n);
       node_set_pt.insert(nod_pt);
       // Was it a duplicate?
       unsigned nnod_now=node_set_pt.size();
       if (nnod_now==nnod_before)
        {
         std::ostringstream warning_stream;
         warning_stream<<"WARNING: " << std::endl
                       <<"Node " << n << " in submesh " << imesh 
                       <<" is a duplicate \n and was ignored when assembling " 
                       << "combined mesh." << std::endl;
         OomphLibWarning(warning_stream.str(),
                         "Mesh::Mesh(const Vector<Mesh*>&)",
                         OOMPH_EXCEPTION_LOCATION);
        }
       else
        {
         Node_pt.push_back(nod_pt);
        }
       nnod_before=nnod_now;
      }
     
     //Loop over the boundaries of the submesh
     unsigned n_bound=sub_mesh_pt[imesh]->nboundary();
     for (unsigned ibound=0;ibound<n_bound;ibound++)
      {
       //Loop over the number of nodes on the boundary and add to the 
       //global vector
       unsigned long n_bound_node=sub_mesh_pt[imesh]->nboundary_node(ibound);
       for (unsigned long n=0;n<n_bound_node;n++)
        {
         Boundary_node_pt[ibound_global].push_back(
          sub_mesh_pt[imesh]->boundary_node_pt(ibound,n));
        }
       //Increase the number of the global boundary counter
       ibound_global++;
      }
    } //End of loop over submeshes
   
  }
 

//========================================================
/// Remove the information about nodes stored on the 
/// b-th boundary of the mesh
//========================================================
void Mesh::remove_boundary_nodes(const unsigned &b)
{
 //Loop over all the nodes on the boundary and call 
 //their remove_from_boundary function
 unsigned n_boundary_node = Boundary_node_pt[b].size();
 for(unsigned n=0;n<n_boundary_node;n++)
  {
   boundary_node_pt(b,n)->remove_from_boundary(b);
  }
 //Clear the storage
 Boundary_node_pt[b].clear();
}

//=================================================================
/// Remove all information about mesh boundaries
//================================================================
void Mesh::remove_boundary_nodes()
{
 //Loop over each boundary call remove_boundary_nodes
 unsigned n_bound = Boundary_node_pt.size();
 for(unsigned b=0;b<n_bound;b++) {remove_boundary_nodes(b);}
 //Clear the storage
 Boundary_node_pt.clear();
}          

//============================================================
/// Remove the node node_pt from the b-th boundary of the mesh
/// This function also removes the information from the Node
/// itself
//===========================================================
void Mesh::remove_boundary_node(const unsigned &b, Node* const &node_pt)
{
 //Find the location of the node in the boundary
 Vector<Node*>::iterator it =
  std::find(Boundary_node_pt[b].begin(),
            Boundary_node_pt[b].end(),
            node_pt);
 //If the node is on this boundary
 if(it!=Boundary_node_pt[b].end())
  {
   //Remove the node from the mesh's list of boundary nodes
   Boundary_node_pt[b].erase(it);
  //Now remove the node's boundary information
   node_pt->remove_from_boundary(b);
  }
 //If not do nothing
}


//========================================================
/// Add the node node_pt to the b-th boundary of the mesh
/// This function also sets the boundary information in the
/// Node itself
//=========================================================
void Mesh::add_boundary_node(const unsigned &b, Node* const &node_pt)
{
 //Tell the node that it's on boundary b.
 //At this point, if the node is not a BoundaryNode, the function
 //should throw an exception.
 node_pt->add_to_boundary(b);
 
 //Get the size of the Boundary_node_pt vector
 unsigned nbound_node=Boundary_node_pt[b].size();
 bool node_already_on_this_boundary=false;
 //Loop over the vector
 for (unsigned n=0; n<nbound_node; n++)
  {
   // is the current node here already?
   if (node_pt==Boundary_node_pt[b][n])
    { 
     node_already_on_this_boundary=true;
    }
  }

 //Add the base node pointer to the vector if it's not there already
 if (!node_already_on_this_boundary)
  {
   Boundary_node_pt[b].push_back(node_pt); 
  }
 
}

//=======================================================
/// Update nodal positions in response to changes in the domain shape.
/// Uses the FiniteElement::get_x(...) function for FiniteElements
/// and doesn't do anything for other element types. \n\n 
/// If a MacroElement pointer has been set for a FiniteElement,
/// the MacroElement representation is used to update the
/// nodal positions; if not get_x(...) uses the FE interpolation
/// and thus leaves the nodal positions unchanged.
/// Virtual, so it can be overloaded by specific meshes,
/// such as AlgebraicMeshes or SpineMeshes. \n\n
/// Generally, this function updates the position of all nodes
/// in response to changes in the boundary position. For
/// SolidNodes it only applies the update to those SolidNodes
/// whose position is determined by the boundary position, unless
/// the bool flag is set to true.
//========================================================
void Mesh::node_update(const bool& update_all_solid_nodes)
{
 /// Local and global (Eulerian) coordinate
 Vector<double> s;
 Vector<double> r;

 // NB: This repeats nodes a lot - surely it would be 
 // quicker to modify it so that it only does each node once,
 // particularly in the update_all_solid_nodes=true case?

 // Loop over all elements 
 unsigned nel=nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt = dynamic_cast<FiniteElement*>(element_pt(e));
   
   // If it's a finite element we can proceed: FiniteElements have
   // nodes and a get_x() function
   if (el_pt!=0)
    {
     // Find out dimension of element = number of local coordinates
     unsigned ndim_el=el_pt->dim();
     s.resize(ndim_el);

     //Loop over nodal points
     unsigned n_node=el_pt->nnode();
     for (unsigned j=0;j<n_node;j++)
      {

       // Get pointer to node
       Node* nod_pt=el_pt->node_pt(j);

       // Get spatial dimension of node
       unsigned ndim_node=nod_pt->ndim();
       r.resize(ndim_node);

       // For non-hanging nodes
       if (!(nod_pt->is_hanging()))
        {
         //Get the position of the node
         el_pt->local_coordinate_of_node(j,s);

         // Get new position
         el_pt->get_x(s,r);

         // Try to cast to SolidNode
         SolidNode* solid_node_pt=dynamic_cast<SolidNode*>(nod_pt);

         // Loop over coordinate directions
         for (unsigned i=0;i<ndim_node;i++)
          {
           // It's a SolidNode:
           if (solid_node_pt!=0)
            {
             // Update nodal positon (only if it's pinned -- otherwise
             // the nodal position is determined by the equations of
             // solid mechanics)
             if (update_all_solid_nodes||
                 (solid_node_pt->position_is_pinned(i)&&
                  !solid_node_pt->is_hanging() ) )
              {
               solid_node_pt->x(i) = r[i];
              }
            }
           // Not a SolidNode: Update regardless
           else
            {
             nod_pt->x(i) = r[i];
            }
          }
        }
      }
    }
  }

 // Now loop over hanging nodes and adjust their position
 // in line with their hanging node constraints
 unsigned long n_node = nnode();
 for(unsigned long n=0;n<n_node;n++)
  {
   Node* nod_pt = Node_pt[n];
   if (nod_pt->is_hanging())
    {
     // Get spatial dimension of node
     unsigned ndim_node=nod_pt->ndim();

     // Initialise
     for (unsigned i=0;i<ndim_node;i++)
      {
       nod_pt->x(i)=0.0;
      }
     
     //Loop over master nodes
     unsigned nmaster=nod_pt->hanging_pt()->nmaster();
     for (unsigned imaster=0;imaster<nmaster;imaster++)
      {
       // Loop over directions
       for (unsigned i=0;i<ndim_node;i++)
        {
         nod_pt->x(i)+=nod_pt->hanging_pt()->
          master_node_pt(imaster)->x(i)*
          nod_pt->hanging_pt()->master_weight(imaster);
        }
      }
    }
  }
 
 // Loop over all nodes again and execute auxiliary node update
 // function
 for(unsigned long n=0;n<n_node;n++)
  {
   Node_pt[n]->perform_auxiliary_node_update_fct();
  }

}


//=======================================================
/// Reorder nodes in the order in which they are
/// encountered when stepping through the elements
//========================================================
void Mesh::reorder_nodes()
{

 // Setup map to check if nodes have been done yet
 std::map<Node*,bool> done;
 
 // Loop over all nodes
 unsigned nnod=nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   done[node_pt(j)]=false;
  }

 // Initialise counter for number of nodes
 unsigned long count=0;

 // Loop over all elements
 unsigned nel=nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Loop over nodes in element
   unsigned nnod=finite_element_pt(e)->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=finite_element_pt(e)->node_pt(j);
     // Has node been done yet?
     if (!done[nod_pt])
      {
       // Insert into node vector
       Node_pt[count]=nod_pt;
       done[nod_pt]=true;
       // Increase counter
       count++;
      }
    }
  }

 // Sanity check
 if (count!=nnod)
  {
   throw OomphLibError(
    "Trouble: Number of nodes hasn't stayed constant during reordering!\n",
    "Mesh::reorder_nodes()",
    OOMPH_EXCEPTION_LOCATION);
  }

}

//========================================================
/// Flush storage for elements and nodes by emptying the
/// vectors that store the pointers to them. This is
/// useful if a particular mesh is only built to generate
/// a small part of a bigger mesh. Once the elements and
/// nodes have been created, they are typically copied
/// into the new mesh and the auxiliary mesh can be
/// deleted. However, if we simply call the destructor
/// of the auxiliary mesh, it will also wipe out
/// the nodes and elements, because it still "thinks"
/// it's in charge of these...
//========================================================
void Mesh::flush_element_and_node_storage()
{
 //Clear vectors of pointers to the nodes and elements
 Node_pt.clear();
 Element_pt.clear();
}

//========================================================
/// Virtual Destructor to clean up all memory
//========================================================
Mesh::~Mesh()
{
 //Free the nodes first
 //Loop over the nodes in reverse order
 unsigned long Node_pt_range = Node_pt.size();
 for(unsigned long i=Node_pt_range;i>0;i--) 
  {
   delete Node_pt[i-1]; Node_pt[i-1] = 0;
  }
 //Free the elements
 //Loop over the elements in reverse order
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long i=Element_pt_range;i>0;i--) 
  {
   delete Element_pt[i-1]; Element_pt[i-1] = 0;
  }
}

//========================================================
/// Assign (global) equation numbers to the nodes
//========================================================
unsigned long Mesh::assign_global_eqn_numbers(Vector<double *> &Dof_pt)
{
 //Find out the current number of equations
 unsigned long equation_number=Dof_pt.size();

 //Loop over the nodes and call their assigment functions
 unsigned long Node_pt_range = Node_pt.size();
 for(unsigned long i=0;i<Node_pt_range;i++)
  {
   Node_pt[i]->assign_eqn_numbers(equation_number,Dof_pt);
  }

 //Loop over the elements and number their internals
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long i=0;i<Element_pt_range;i++)
  {  
   Element_pt[i]->assign_internal_eqn_numbers(equation_number,Dof_pt);
  }

 //Return the total number of equations
 return(equation_number);
}

//========================================================
/// Assign local equation numbers in all elements
//========================================================
void Mesh::assign_local_eqn_numbers()
{ 
 //Now loop over the elements and assign local equation numbers
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long i=0;i<Element_pt_range;i++)
  {
   Element_pt[i]->assign_local_eqn_numbers();
  }
}

//========================================================
/// Self-test: Check elements and nodes. Return 0 for OK
//========================================================
unsigned Mesh::self_test()
{ 

 // Initialise
 bool passed=true;

 // Check the mesh for repeated nodes (issues its own error message)
 if (0!=check_for_repeated_nodes()) passed=false;

 //Loop over the elements, check for duplicates and do self test
 std::set<GeneralisedElement*> element_set_pt;
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long i=0;i<Element_pt_range;i++)
  {  
   if (Element_pt[i]->self_test()!=0)
    {
     passed=false;
     oomph_info << "\n ERROR: Failed Element::self_test() for element i=" 
               << i << std::endl;
    }
   // Add to set (which ignores duplicates):
   element_set_pt.insert(Element_pt[i]);
  }

 // Check for duplicates:
 if (element_set_pt.size()!=Element_pt_range)
  {
   oomph_info << "ERROR:  " << Element_pt_range-element_set_pt.size() 
             << " duplicate elements were encountered in mesh!" << std::endl;
   passed=false;
  }
 
 
 //Loop over the nodes, check for duplicates and do self test
 std::set<Node*> node_set_pt;
 unsigned long Node_pt_range = Node_pt.size();
 for(unsigned long i=0;i<Node_pt_range;i++)
  {
   if (Node_pt[i]->self_test()!=0)
    {
     passed=false;
     oomph_info << "\n ERROR: Failed Node::self_test() for node i=" 
               << i << std::endl;
    }
   // Add to set (which ignores duplicates):
   node_set_pt.insert(Node_pt[i]);
  }

 // Check for duplicates:
 if (node_set_pt.size()!=Node_pt_range)
  {
   oomph_info << "ERROR:  " << Node_pt_range-node_set_pt.size() 
             << " duplicate nodes were encountered in mesh!" << std::endl;
   passed=false;
  }

 // Return verdict
 if (passed) {return 0;}
 else {return 1;}

}


//========================================================
/// Nodes that have been marked as obsolete are removed
/// from the mesh and the its boundaries
//========================================================
void Mesh::prune_dead_nodes()
{
 // Only copy the 'live' nodes across to new mesh
 //----------------------------------------------
 // New Vector of pointers to nodes
 Vector<Node *> new_node_pt;

 // Loop over all nodes in mesh 
 unsigned long n_node = nnode();
 for(unsigned long n=0;n<n_node;n++)
  {
   // If the node still exists: Copy across
   if(!(Node_pt[n]->is_obsolete())) {new_node_pt.push_back(Node_pt[n]);}
   // Otherwise the Node is gone: 
   // Delete it for good if it does not lie on a boundary
   else {if(Node_pt[n]->is_on_boundary()==false) {delete Node_pt[n];}}
  }
 
 // Now update old vector by setting it equal to the new vector
 Node_pt = new_node_pt;
 
 //---BOUNDARIES
 Vector<double> bnd_coords;
 // Only copy the 'live' nodes into new boundary node arrays
 //---------------------------------------------------------
 //Loop over the boundaries
 unsigned num_bound = nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // New Vector of pointers to existent boundary nodes 
   Vector<Node *> new_boundary_node_pt;

   //Loop over the boundary nodes
   unsigned long Nboundary_node = Boundary_node_pt[ibound].size();

   // Reserve contiguous memory for new vector of pointers
   // Must be equal in size to the number of nodes or less
   new_boundary_node_pt.reserve(Nboundary_node);

   for(unsigned long n=0;n<Nboundary_node;n++)
    {
     // If node still exists: Copy across
     if (!(Boundary_node_pt[ibound][n]->is_obsolete()))
      {new_boundary_node_pt.push_back(Boundary_node_pt[ibound][n]);}
       // Otherwise Node is gone: Delete it for good
     else
      {
       //The node may lie on multiple boundaries, so remove the node 
       //from the current boundary
       //Remove the node from the bounadry
       Boundary_node_pt[ibound][n]->remove_from_boundary(ibound);
       //Now if the node is no longer on any boundaries, delete it
       if(!Boundary_node_pt[ibound][n]->is_on_boundary())
        {
         delete Boundary_node_pt[ibound][n];
        }
      }
    }

   //Update the old vector by setting it equal to the new vector
   Boundary_node_pt[ibound] = new_boundary_node_pt;
  
  } //End of loop over boundaries
}




//========================================================
/// Output function for the mesh boundaries
///
/// Loop over all boundaries and dump out the coordinates
/// of the points on the boundary (in individual tecplot
/// zones) 
//========================================================
void Mesh::output_boundaries(std::ostream &outfile)
{
 //Loop over the boundaries
 unsigned num_bound = nboundary();
 for(unsigned long ibound=0;ibound<num_bound;ibound++)
  {
   unsigned nnod=Boundary_node_pt[ibound].size();
   if (nnod>0)
    {
     outfile << "ZONE T=\"boundary" << ibound << "\"\n";
     
     for (unsigned inod=0;inod<nnod;inod++)
      {
       Boundary_node_pt[ibound][inod]->output(outfile);
      }
    }
  }
}



//===================================================================
/// Dump function for the mesh class.
/// Loop over all nodes and elements and dump them
//===================================================================
void Mesh::dump(std::ofstream &dump_file)
{
 // Find number of nodes
 unsigned long Node_pt_range = this->nnode(); 
   
 // Doc # of nodes
 dump_file << Node_pt_range << " # number of nodes " << std::endl;
   
 //Loop over all the nodes and dump their data
 for(unsigned long n=0;n<Node_pt_range;n++)
  {
   this->node_pt(n)->dump(dump_file);
  }
 
 // Loop over elements and deal with internal data
 unsigned n_element = this->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   GeneralisedElement* el_pt = this->element_pt(e);
   unsigned n_internal = el_pt->ninternal_data();
   if(n_internal > 0)
    {
     dump_file << n_internal 
               << " # number of internal Data items in element " 
               << e << std::endl;
     for(unsigned i=0;i<n_internal;i++)
      {
       el_pt->internal_data_pt(i)->dump(dump_file);
      }
    }
  }
}


//=======================================================
/// Read solution from restart file
//=======================================================
void Mesh::read(std::ifstream &restart_file)
{
 std::string input_string;
 
 //Read nodes
 
 // Find number of nodes
 unsigned long n_node = this->nnode(); 
 
 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 
 // Check # of nodes:
 unsigned long check_n_node=atoi(input_string.c_str());
 if (check_n_node!=n_node)
  {
   std::ostringstream error_stream;
   error_stream << "The number of nodes allocated " << n_node 
                << " is not the same as specified in the restart file "
                << check_n_node << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "Mesh::read()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 //Loop over the nodes
 for(unsigned long n=0;n<n_node;n++)
  {     
   /// Try to cast to elastic node 
   SolidNode* el_node_pt=dynamic_cast<SolidNode*>(
    this->node_pt(n));
   if (el_node_pt!=0)
    {
     el_node_pt->read(restart_file);
    }
   else
    {
     this->node_pt(n)->read(restart_file);
    }
  }
 
 // Read internal data of elements:
 //--------------------------------
 // Loop over elements and deal with internal data
 unsigned n_element = this->nelement();
 for (unsigned e=0;e<n_element;e++)
  {
   GeneralisedElement* el_pt = this->element_pt(e);
   unsigned n_internal=el_pt->ninternal_data();
   if (n_internal>0)
    {
     // Read line up to termination sign
     getline(restart_file,input_string,'#');
     
     // Ignore rest of line
     restart_file.ignore(80,'\n');
     
     // Check # of internals :
     unsigned long check_n_internal=atoi(input_string.c_str());
     if (check_n_internal!=n_internal)
      {
       std::ostringstream error_stream;
       error_stream << "The number of internal data  " << n_internal 
                    << " is not the same as specified in the restart file "
                    << check_n_internal << std::endl;
   
       throw OomphLibError(error_stream.str(),
                           "Mesh::read()",
                           OOMPH_EXCEPTION_LOCATION);
      }
     
     for (unsigned i=0;i<n_internal;i++)
      {
       el_pt->internal_data_pt(i)->read(restart_file);
      }
    }
  }
}



//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function)
//========================================================
void Mesh::output(std::ostream &outfile)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output(outfile);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output(outfile);
        }
      }
#endif
    }
  }
}

//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function). Use
/// n_plot plot points in each coordinate direction.
//========================================================
void Mesh::output(std::ostream &outfile, const unsigned &n_plot)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();

 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output(outfile,n_plot);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output(outfile,n_plot);
        }
      }
#endif
    }
  }
}





//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function)
/// (C style output)
//========================================================
void Mesh::output(FILE* file_pt)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output(file_pt);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output(file_pt);
        }
      }
#endif
    }
  }
}

//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function). Use
/// n_plot plot points in each coordinate direction.
/// (C style output)
//========================================================
void Mesh::output(FILE* file_pt, const unsigned &n_plot)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output(file_pt,n_plot);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output(file_pt,n_plot);
        }
      }
#endif
    }
  }
}


//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function). Use
/// n_plot plot points in each coordinate direction.
//========================================================
void Mesh::output_fct(std::ostream &outfile, const unsigned &n_plot, 
                      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output_fct(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output_fct(outfile,n_plot,exact_soln_pt);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output_fct(outfile,n_plot,exact_soln_pt);
        }
      }
#endif
    }
  }
}

//========================================================
/// Output function for the mesh class
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function) at time t. Use
/// n_plot plot points in each coordinate direction.
//========================================================
void Mesh::output_fct(std::ostream &outfile, const unsigned &n_plot, 
                      const double& time, 
                      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
{
 //Loop over the elements and call their output functions
 //Assign Element_pt_range
 unsigned long Element_pt_range = Element_pt.size();
 for(unsigned long e=0;e<Element_pt_range;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     oomph_info << "Can't execute output_fct(...) for non FiniteElements" 
               << std::endl;
    }
   else
    {
#ifdef OOMPH_HAS_MPI
     if (Output_halo_elements)
#endif
      {
       el_pt->output_fct(outfile,n_plot,time,exact_soln_pt);
      }
#ifdef OOMPH_HAS_MPI
     else
      {
       if (!el_pt->is_halo())
        {
         el_pt->output_fct(outfile,n_plot,time,exact_soln_pt);
        }
      }
#endif
    }
  }
}

//==================================================================
/// Assign the initial values for an impulsive start, which is 
/// acheived by looping over all data in the mesh (internal element
/// data and data stored at nodes) and setting the calling the 
/// assign_initial_values_impulsive() function for each data's 
/// timestepper
//=================================================================
void Mesh::assign_initial_values_impulsive()
{
 //Loop over the elements
 unsigned long Nelement=nelement();
 for(unsigned long e=0;e<Nelement;e++)
  {
   //Find the number of internal dofs
   unsigned Ninternal = element_pt(e)->ninternal_data();
   //Loop over internal dofs and shift the time values
   //using the internals data's timestepper
   for(unsigned j=0;j<Ninternal;j++)
    {
     element_pt(e)->internal_data_pt(j)->time_stepper_pt()->
      assign_initial_values_impulsive(element_pt(e)->internal_data_pt(j));
    }
  }
 
 //Loop over the nodes
 unsigned long n_node=nnode();
 for (unsigned long n=0;n<n_node;n++)
  {
   // Assign initial values using the Node's timestepper
   Node_pt[n]->time_stepper_pt()->
    assign_initial_values_impulsive(Node_pt[n]);
   // Assign initial positions using the Node's timestepper
   Node_pt[n]->position_time_stepper_pt()->
    assign_initial_positions_impulsive(Node_pt[n]);
  }

}

//===============================================================
/// Shift time-dependent data along for next timestep:
/// Again this is achieved by looping over all data and calling
/// the functions defined in each data object's timestepper.
//==============================================================
void Mesh::shift_time_values()
{
 // Loop over the elements which shift their internal data
 // via their own timesteppers
 unsigned long Nelement=nelement();
 for (unsigned long e=0;e<Nelement;e++)
  {
   //Find the number of internal dofs
   unsigned Ninternal = element_pt(e)->ninternal_data();
   //Loop over internal dofs and shift the time values
   //using the internals data's timestepper
   for(unsigned j=0;j<Ninternal;j++)
    {
     element_pt(e)->internal_data_pt(j)->time_stepper_pt()->
      shift_time_values(element_pt(e)->internal_data_pt(j));
    }
  }
 
 //Loop over the nodes
 unsigned long n_node=nnode();
 for (unsigned long n=0;n<n_node;n++)
  {
   // Shift the Data associated with the nodes with the Node's own 
   // timestepper
   Node_pt[n]->time_stepper_pt()->shift_time_values(Node_pt[n]);
   // Push history of nodal positions back
   Node_pt[n]->position_time_stepper_pt()->shift_time_positions(Node_pt[n]);
  }

}

//=========================================================================
/// Calculate predictions for all Data and positions associated 
/// with the mesh. This is usually only used for adaptive time-stepping
/// when the comparison between a predicted value and the actual value
/// is usually used to determine the change in step size. Again the
/// loop is over all data in the mesh and individual timestepper functions
/// for each data value are called.
//=========================================================================
void Mesh::calculate_predictions()
{
 // Loop over the elements which shift their internal data
 // via their own timesteppers
 unsigned long Nelement=nelement();
 for (unsigned long e=0;e<Nelement;e++)
  {
   //Find the number of internal dofs
   unsigned Ninternal = element_pt(e)->ninternal_data();
   //Loop over internal dofs and calculate predicted positions
   //using the internals data's timestepper
   for(unsigned j=0;j<Ninternal;j++)
    {
     element_pt(e)->internal_data_pt(j)->time_stepper_pt()->
      calculate_predicted_values(element_pt(e)->internal_data_pt(j));
    }
  }
 
 //Loop over the nodes
 unsigned long n_node=nnode();
 for (unsigned long n=0;n<n_node;n++)
  {
   // Calculate the predicted values at the nodes
   Node_pt[n]->time_stepper_pt()->calculate_predicted_values(Node_pt[n]);
   //Calculate the predicted positions
   Node_pt[n]->position_time_stepper_pt()->
    calculate_predicted_positions(Node_pt[n]);
  }
}

//========================================================================
/// A function that upgrades an ordinary node to a boundary node.
/// All pointers to the node from the mesh's elements are found.
/// and replaced by pointers to the new boundary node. If the node
/// is present in the mesh's list of nodes, that pointer is also
/// replaced. Finally, the pointer argument node_pt addresses the new
/// node on return from the function.
/// We shouldn't ever really use this, but it does make life that
/// bit easier for the lazy mesh writer.
//=======================================================================
void Mesh::convert_to_boundary_node(Node* &node_pt)
{

 //If the node is already a boundary node, then return straight away,
 //we don't need to do anything
 if (dynamic_cast<BoundaryNodeBase*>(node_pt)!=0)
  {
   return;
  }

 //Loop over all the elements in the mesh and find all those in which
 //the present node is referenced and the corresponding local node number
 //in those elements.
 
 //Storage for elements and local node number
 std::list<std::pair<unsigned long, int> > 
  list_of_elements_and_local_node_numbers;
 
 //Loop over all elements
 unsigned long n_element=this->nelement();
 for(unsigned long e=0;e<n_element;e++)
  {
   //Buffer the case when we have not yet filled up the element array
   //Unfortunately, we should not assume that the array has been filled
   //in a linear order, so we can't break out early.
   if(Element_pt[e]!=0)
    {
     //Find the local node number of the passed node
     int node_number = finite_element_pt(e)->get_node_number(node_pt);
     //If the node is present in the element, add it to our list and 
     //NULL out the local element entries
     if(node_number!=-1)
      {
       list_of_elements_and_local_node_numbers.insert(
        list_of_elements_and_local_node_numbers.end(),
        std::make_pair(e,node_number));
       //Null it out
       finite_element_pt(e)->node_pt(node_number)=0;
      }
    }
  } //End of loop over elements

 //If there are no entries in the list we are in real trouble
 if(list_of_elements_and_local_node_numbers.empty())
  {
   std::ostringstream error_stream;
   error_stream << "Node " << node_pt 
                << " is not contained in any elements in the Mesh."
                << std::endl
                << "How was it created then?" << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "Mesh::convert_to_boundary_node()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 
 //Create temporary storage for a pointer to the old node.
 //This is required because if we have passed a reference to the
 //first element that we find, constructing the new node
 //will over-write our pointer and we'll get segmentation faults.
 Node* old_node_pt = node_pt;
 
 //We now create the new node by using the first element in the list
 std::list<std::pair<unsigned long,int> >::iterator list_it
  = list_of_elements_and_local_node_numbers.begin();
 
 //Create a new boundary node, using the timestepper from the
 //original node
 Node* new_node_pt = finite_element_pt(list_it->first)
  ->construct_boundary_node(list_it->second,node_pt->time_stepper_pt());

  //Now copy all the information accross from the old node
 
 //Can we cast the node to a solid node
 SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(new_node_pt);
 //If it's a solid node, do the casting
 if(solid_node_pt!=0)
  {
   solid_node_pt->copy(dynamic_cast<SolidNode*>(old_node_pt));
  }
 else
  {
   new_node_pt->copy(old_node_pt);
  }

 //Loop over all other elements in the list and set their pointers
 //to the new node
 for(++list_it; //Increment the iterator
     list_it!=list_of_elements_and_local_node_numbers.end();
     ++list_it)
  {
   finite_element_pt(list_it->first)->node_pt(list_it->second)
    = new_node_pt;
  }
 
 //Finally, find the position of the node in the global mesh
 Vector<Node*>::iterator it=
  std::find(Node_pt.begin(),Node_pt.end(),old_node_pt);
 
 //If it is in the mesh, update the pointer
 if(it!=Node_pt.end()) {*it = new_node_pt;}
   
 //Can now delete the old node
 delete old_node_pt;
  
 //Replace the passed pointer by a pointer to the new node
 //Note that in most cases, this will be wasted work because the node
 //pointer will either the pointer in the mesh's or an element's 
 //node_pt vector. Still assignment is quicker than an if to check this.
 node_pt = new_node_pt;

}



#ifdef OOMPH_HAS_MPI


//========================================================================
/// Classify all halo and haloed information in the mesh
//========================================================================
void Mesh::classify_halo_and_haloed_nodes(OomphCommunicator* comm_pt,
                                          DocInfo& doc_info,
                                          const bool& report_stats)
{
 //Wipe existing storage schemes for halo(ed) nodes
 Halo_node_pt.clear();
 Haloed_node_pt.clear();

 // Storage for number of processors and current processor
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();
 MPI_Status status;
 
 // Determine which processors the nodes are associated with
 // and hence who's in charge
 std::map<Data*,std::set<unsigned> > processors_associated_with_data;
 std::map<Data*,unsigned> processor_in_charge;
 
 // Loop over all processors and associate any nodes on the halo
 // elements involved with that processor
 for (int domain=0;domain<n_proc;domain++) 
  {   
   // Get vector of halo elements by copy operation
   Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(domain));
   
   // Loop over halo elements associated with this adjacent domain
   unsigned nelem=halo_elem_pt.size();
   for (unsigned e=0;e<nelem;e++)
    {
     // Get element
     FiniteElement* el_pt=halo_elem_pt[e];
     
     //Loop over nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=el_pt->node_pt(j);
       // Associate node with this domain
       processors_associated_with_data[nod_pt].insert(domain);

       // Do the same if the node is solid
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
       if (solid_nod_pt!=0)
        {
         processors_associated_with_data[solid_nod_pt->variable_position_pt()].
          insert(domain);
        }
      }
    }
  }
 
 
 // Loop over all [non-halo] elements and associate their nodes
 // with current procesor
 unsigned nelem=this->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=this->finite_element_pt(e);
   
   // Only visit non-halos
   if (!el_pt->is_halo())
    {
     // Loop over nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=el_pt->node_pt(j);

       // Associate this node with current processor
       processors_associated_with_data[nod_pt].insert(my_rank);

       // do the same if we have a SolidNode
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
       if (solid_nod_pt!=0)
        {
         processors_associated_with_data
          [solid_nod_pt->variable_position_pt()].insert(my_rank);
        }
      }
    }
  }

 // At this point we need to "syncrhonise" the nodes on halo(ed) elements
 // so that the processors_associated_with_data agrees for the same node
 // on all processors [this is only a problem in 3D, where hanging nodes
 // on an edge which is at a junction between at least 3 processors do not
 // necessarily know that they should be associated with all of the processors
 // at the junction, on all of the processors]
 
 // Loop over all domains
 for (int d=0;d<n_proc;d++)
  {
   // Prepare vector to send/receive
   Vector<unsigned> processors_associated_with_data_on_other_proc;

   if (d!=my_rank)
    {
     // Communicate the processors associated with data on haloed elements

     // Get haloed elements
     Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(d));

     // Initialise counter for this haloed layer
     unsigned count_data=0;

     // Loop over haloed elements
     unsigned n_haloed_elem=haloed_elem_pt.size();
     for (unsigned e=0;e<n_haloed_elem;e++)
      {
       FiniteElement* haloed_el_pt=haloed_elem_pt[e];
       // Loop over nodes
       unsigned n_node=haloed_el_pt->nnode();
       for (unsigned j=0;j<n_node;j++)
        {
         Node* nod_pt=haloed_el_pt->node_pt(j);

         // Number of processors associated with this node
         unsigned n_assoc=processors_associated_with_data[nod_pt].size();
       
         // This number needs to be sent
         processors_associated_with_data_on_other_proc.push_back(n_assoc);
         count_data++;

         // Now add the process IDs associated to the vector to be sent
         std::set<unsigned> procs_set=processors_associated_with_data[nod_pt];
         for (std::set<unsigned>::iterator it=procs_set.begin();
              it!=procs_set.end();it++)
          {
           processors_associated_with_data_on_other_proc.push_back(*it);
           count_data++;
          }
        }
      }

     // Send the information
     MPI_Send(&count_data,1,MPI_INT,d,0,comm_pt->mpi_comm());
     if (count_data!=0)
      {
       MPI_Send(&processors_associated_with_data_on_other_proc[0],count_data,
                MPI_INT,d,1,comm_pt->mpi_comm());
      }
    }
   else
    {
     // Receive the processors associated with data onto halo elements
     for (int dd=0;dd<n_proc;dd++)
      {
       if (dd!=my_rank) // (my_rank=d)
        {
         // We will be looping over the halo elements with process dd
         Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(dd));
         unsigned n_halo_elem=halo_elem_pt.size();

         unsigned count_data=0;
         MPI_Recv(&count_data,1,MPI_INT,dd,0,comm_pt->mpi_comm(),&status);

         if (count_data!=0)
          {
           processors_associated_with_data_on_other_proc.resize(count_data);

           MPI_Recv(&processors_associated_with_data_on_other_proc[0],
                    count_data,MPI_INT,dd,1,comm_pt->mpi_comm(),&status);

           // Reset counter and loop through nodes on halo elements
           count_data=0;
           for (unsigned e=0;e<n_halo_elem;e++)
            {
             FiniteElement* halo_el_pt=halo_elem_pt[e];
             unsigned n_node=halo_el_pt->nnode();
             for (unsigned j=0;j<n_node;j++)
              {
               Node* nod_pt=halo_el_pt->node_pt(j);

               // Get number of processors associated with data that was sent
               unsigned n_assoc=
                processors_associated_with_data_on_other_proc[count_data];
               count_data++;

               for (unsigned i_assoc=0;i_assoc<n_assoc;i_assoc++)
                {
                 // Get the process ID
                 unsigned sent_domain=
                  processors_associated_with_data_on_other_proc[count_data];
                 count_data++;

                 // Add it to this processor's list of IDs
                 processors_associated_with_data[nod_pt].insert(sent_domain);

                 // If the node is solid then add the ID to the solid data
                 SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
                 if (solid_nod_pt!=0)
                  {
                   processors_associated_with_data
                    [solid_nod_pt->variable_position_pt()].insert(sent_domain);
                  }
                }
              }
            }
          }
        }
      }
    }
  }


 // Loop over all nodes on the present processor and put the highest-numbered
 // processor associated with each node "in charge" of the node
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=this->node_pt(j);

   // Reset halo status of node to false
   nod_pt->is_halo()=false;

   // If it's a SolidNode then the halo status of the data 
   // associated with that must also be reset to false
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt!=0)
    {
     solid_nod_pt->variable_position_pt()->is_halo()=false;
    }

   // Now put the highest-numbered one in charge
   unsigned proc_max=0;
   std::set<unsigned> procs_set=processors_associated_with_data[nod_pt];
   for (std::set<unsigned>::iterator it=procs_set.begin();
        it!=procs_set.end();it++)
    {
     if (*it>proc_max) proc_max=*it;
    }
   processor_in_charge[nod_pt]=proc_max;

   // Do the same if we have a SolidNode
   if (solid_nod_pt!=0)
    {
     // Now put the highest-numbered one in charge
     unsigned proc_max_solid=0;
     std::set<unsigned> procs_set_solid=processors_associated_with_data
      [solid_nod_pt->variable_position_pt()];
     for (std::set<unsigned>::iterator it=procs_set_solid.begin();
          it!=procs_set_solid.end();it++)
      {
       if (*it>proc_max_solid) proc_max_solid=*it;
      }
     processor_in_charge[solid_nod_pt->variable_position_pt()]=proc_max_solid;
    }
  }

 // Determine halo nodes. They are located on the halo
 // elements and the processor in charge differs from the
 // current processor

 // Only count nodes once (map is initialised to 0 = false)
 std::map<Node*,bool> done;

 // Loop over all processors
 for (int domain=0;domain<n_proc;domain++) 
  {
   // Get vector of halo elements by copy operation
   Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(domain));

   // Loop over halo elements associated with this adjacent domain
   unsigned nelem=halo_elem_pt.size();

   for (unsigned e=0;e<nelem;e++)
    {
     // Get element
     FiniteElement* el_pt=halo_elem_pt[e];

     //Loop over nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=el_pt->node_pt(j);
       
       // Have we done this node already?
       if (!done[nod_pt])
        {
         // Is the other processor/domain in charge of this node?
         int proc_in_charge=processor_in_charge[nod_pt];

         // To keep the order of the nodes consistent with that
         // in the haloed node lookup scheme, only 
         // allow it to be added when the current domain is in charge
         if (proc_in_charge==int(domain))
          {
           if (proc_in_charge!=my_rank)
            {
             // Add it as being halo node whose non-halo counterpart
             // is located on processor domain
             this->add_halo_node_pt(proc_in_charge,nod_pt);

             // The node itself needs to know it is a halo
             nod_pt->is_halo()=true;

             // If it's a SolidNode then the data associated with that
             // must also be halo
             SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
             if (solid_nod_pt!=0)
              {
               solid_nod_pt->variable_position_pt()->is_halo()=true;
              }

             // We're done with this node
             done[nod_pt]=true;
            }
          }

        }

      }
     // Now make sure internal data on halo elements is also halo
     unsigned nintern_data = el_pt->ninternal_data();
     for (unsigned iintern=0;iintern<nintern_data;iintern++)
      {
       el_pt->internal_data_pt(iintern)->is_halo()=true;
      }
    }
  }

 // Determine haloed nodes. They are located on the haloed
 // elements and the processor in charge is the current processor
 
 // Loop over processors that share haloes with the current one
 for (int domain=0;domain<n_proc;domain++) 
  {
   // Only count nodes once (map is initialised to 0 = false)
   std::map<Node*,bool> node_done;
  
   // Get vector of haloed elements by copy operation
   Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(domain));

   // Loop over haloed elements associated with this adjacent domain
   unsigned nelem=haloed_elem_pt.size();

   for (unsigned e=0;e<nelem;e++)
    {
     // Get element
     FiniteElement* el_pt=haloed_elem_pt[e];

     //Loop over nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=el_pt->node_pt(j);
       
       // Have we done this node already?
       if (!node_done[nod_pt])
        {
         // Is the current processor/domain in charge of this node?
         int proc_in_charge=processor_in_charge[nod_pt];

         if (proc_in_charge==my_rank)
          {
           // Add it as being haloed from specified domain
           this->add_haloed_node_pt(domain,nod_pt);
           // We're done with this node
           node_done[nod_pt]=true;
          }
        }

      }
    }
  }

 // Doc stats
 if (report_stats)
  {
   // Report total number of halo(ed) and shared nodes for this process
   oomph_info << "Processor " << my_rank 
              << " holds " << this->nnode() 
              << " nodes of which " << this->nhalo_node()
              << " are halo nodes \n while " << this->nhaloed_node()
              << " are haloed nodes, and " << this->nshared_node()
              << " are shared nodes." << std::endl;   

   // Report number of halo(ed) and shared nodes with each domain 
   // from the current process
   for (int iproc=0;iproc<n_proc;iproc++)
    {
     oomph_info << "With process " << iproc << ", there are " 
                << this->nhalo_node(iproc) << " halo nodes, and " << std::endl
                << this->nhaloed_node(iproc) << " haloed nodes, and " 
                << this->nshared_node(iproc) << " shared nodes" << std::endl; 
    }
  }

}





//========================================================================
/// Get halo node stats for this distributed mesh:
/// Average/max/min number of halo nodes over all processors.
/// \b Careful: Involves MPI Broadcasts and must therefore
/// be called on all processors!
//========================================================================
void Mesh::get_halo_node_stats(OomphCommunicator* comm_pt,
                               double& av_number,
                               unsigned& max_number,
                               unsigned& min_number)
{
 // Storage for number of processors and current processor
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();

 // Create vector to hold number of halo nodes
 Vector<int> nhalo_nodes(n_proc);
 
 // Stick own number of halo nodes into appropriate entry 
 nhalo_nodes[my_rank]=nhalo_node();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of dofs on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhalo_nodes[my_rank],1,MPI_INT,
            &nhalo_nodes[0],1, MPI_INT,
            0,comm_pt->mpi_comm());

 // Initialise stats
 av_number=0.0;
 int max=-1;
 int min=1000000000;

 if (my_rank==0)
  {
   for (int i=0;i<n_proc;i++)
    {
     av_number+=double(nhalo_nodes[i]);
     if (int(nhalo_nodes[i])>max) max=nhalo_nodes[i];
     if (int(nhalo_nodes[i])<min) min=nhalo_nodes[i];

    }  
   av_number/=double(n_proc); 
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_INT,0,comm_pt->mpi_comm());
 MPI_Bcast(&min,1,MPI_INT,0,comm_pt->mpi_comm());
 MPI_Bcast(&av_number,1,MPI_DOUBLE,0,comm_pt->mpi_comm());
 
 max_number=max;
 min_number=min;  
}


//========================================================================
/// Get haloed node stats for this distributed mesh:
/// Average/max/min number of haloed nodes over all processors.
/// \b Careful: Involves MPI Broadcasts and must therefore
/// be called on all processors!
//========================================================================
 void  Mesh::get_haloed_node_stats(OomphCommunicator* comm_pt,
                                   double& av_number,
                                   unsigned& max_number,
                                   unsigned& min_number)
{
 // Storage for number of processors and current processor
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();

 // Create vector to hold number of haloed nodes
 Vector<int> nhaloed_nodes(n_proc);
 
 // Stick own number of haloed nodes into appropriate entry
 nhaloed_nodes[my_rank]=nhaloed_node();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of dofs on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhaloed_nodes[my_rank],1,MPI_INT,
            &nhaloed_nodes[0],1, MPI_INT,
            0,comm_pt->mpi_comm());

 // Initialise stats
 av_number=0.0;
 int max=-1;
 int min=1000000000;

 if (my_rank==0)
  {
   for (int i=0;i<n_proc;i++)
    {
     av_number+=double(nhaloed_nodes[i]);
     if (int(nhaloed_nodes[i])>max) max=nhaloed_nodes[i];
     if (int(nhaloed_nodes[i])<min) min=nhaloed_nodes[i];

    }  
   av_number/=double(n_proc); 
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_INT,0,comm_pt->mpi_comm());
 MPI_Bcast(&min,1,MPI_INT,0,comm_pt->mpi_comm());
 MPI_Bcast(&av_number,1,MPI_DOUBLE,0,comm_pt->mpi_comm());
 
 max_number=max;
 min_number=min;  
}

//========================================================================
/// Distribute the mesh
//========================================================================
void Mesh::distribute(OomphCommunicator* comm_pt,
                      const Vector<unsigned>& element_domain,
                      DocInfo& doc_info,
                      const bool& report_stats)
{ 
 // Storage for number of processors and current processor
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();

 // Storage for number of elements and number of nodes on this mesh
 unsigned nelem=this->nelement();
 unsigned nnod=this->nnode();

 char filename[100];

 // Doc the partitioning (only on processor 0) 
 //-------------------------------------------
 if (doc_info.doc_flag())
  {
   if (my_rank==0)
    {
     // Open files for doc of element partitioning
     Vector<std::ofstream*> domain_file(n_proc);
     for (int d=0;d<n_proc;d++)
      {
       // Note: doc_info.number() was set in Problem::distribute(...) to
       // reflect the submesh number
       sprintf(filename,"%s/domain%i-%i.dat",doc_info.directory().c_str(),
               d,doc_info.number());
       domain_file[d]=new std::ofstream(filename);
      }
     
     // Doc
     for (unsigned e=0;e<nelem;e++)
      {
       this->finite_element_pt(e)->
        output(*domain_file[element_domain[e]],5);
      }
     
     for (int d=0;d<n_proc;d++)
      {
       domain_file[d]->close();
      }
    }
  }

 // Loop over all elements, associate all 
 //--------------------------------------
 // nodes with the highest-numbered processor and record all
 //---------------------------------------------------------
 // processors the node is associated with
 //---------------------------------------

 // Storage for processors in charge and processors associated with data
 std::map<Data*,std::set<unsigned> > processors_associated_with_data;
 std::map<Data*,unsigned> processor_in_charge;

 // For all nodes set the processor in charge to zero
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=this->node_pt(j);
   processor_in_charge[nod_pt]=0;
  }

 // Loop over elements
 for (unsigned e=0;e<nelem;e++)
  {
   // Get element and its domain
   FiniteElement* el_pt=this->finite_element_pt(e);
   unsigned el_domain=element_domain[e];

   // Associate nodes with highest numbered processor
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=el_pt->node_pt(j); 

     // processor in charge was initialised to 0 above
     if (el_domain>processor_in_charge[nod_pt])
      {
       processor_in_charge[nod_pt]=el_domain;
      }
     processors_associated_with_data[nod_pt].insert(el_domain);
    }
  }

 // Doc the partitioning (only on processor 0) 
 //-------------------------------------------
 if (doc_info.doc_flag())
  {
   if (my_rank==0)
    {
     // Open files for doc of node partitioning
     Vector<std::ofstream*> node_file(n_proc);
     for (int d=0;d<n_proc;d++)
      {
       // Note: doc_info.number() was set in Problem::distribute(...) to
       // reflect the submesh number...
       sprintf(filename,"%s/node%i-%i.dat",doc_info.directory().c_str(),
               d,doc_info.number());
       node_file[d]=new std::ofstream(filename);
      }
     
     // Doc
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=this->node_pt(j);
       *node_file[processor_in_charge[nod_pt]]
        << nod_pt->x(0) << " " 
        << nod_pt->x(1) << std::endl;
      }
     for (int d=0;d<n_proc;d++)
      {
       node_file[d]->close();
      }
    }
  }

 // Declare all nodes as obsolete. We'll
 // change this setting for all nodes that must be retained
 // further down
 for (unsigned j=0;j<nnod;j++)
  {
   this->node_pt(j)->set_obsolete();
  }


 // Backup old mesh data and flush mesh
 //-------------------------------------

 // Backup pointers to elements in this mesh
 Vector<FiniteElement*> backed_up_el_pt(nelem);
 for (unsigned e=0;e<nelem;e++)
  {
   backed_up_el_pt[e]=this->finite_element_pt(e);
  }

 // Backup pointers to nodes in this mesh
 Vector<Node*> backed_up_nod_pt(nnod);
 for (unsigned j=0;j<nnod;j++)
  {
   backed_up_nod_pt[j]=this->node_pt(j);
  }

 // Flush the mesh storage
 this->flush_element_and_node_storage();

 // Flush any storage of external elements and nodes
 this->flush_all_external_storage();

 // Now loop over the (backed up) elements and identify the ones
 //--------------------------------------------------------------
 // that must be retained. 
 //-----------------------

 // Boolean to indicate which element is to be retained on 
 // which processor. This is needed because elements have
 // to be added into the various halo/haloed lookup schemes in the
 // same order and we base this on the order in which
 // they were in the original mesh.
 Vector<std::vector<bool> > tmp_element_retained;
 tmp_element_retained.resize(n_proc);
 nelem=backed_up_el_pt.size();
 for (int i=0;i<n_proc;i++)
  {
   tmp_element_retained[i].resize(nelem,false);
  }
 
 // Temporary storage for root halo elements on the various
 // processors. Needed to figure out haloed lookup schemes:
 // When setting these up on any given processor we have to know
 // which elements will (have) become halo elements on other processors.
 Vector<Vector<Vector<FiniteElement*> > > tmp_root_halo_element_pt;
 tmp_root_halo_element_pt.resize(n_proc);
 for (int i=0;i<n_proc;i++)
  {
   tmp_root_halo_element_pt[i].resize(n_proc);
  }

 // Determine which elements are going to end up on which processor
 //----------------------------------------------------------------

 // This procedure needs to be repeated to catch elements which may
 // be missed the first time round but which contain nodes from this process

 unsigned elements_retained=true;
 int myi=1;
 while (elements_retained) 
  {
   Vector<unsigned> number_of_retained_elements(n_proc,0);
   int number_of_retained_halo_elements=0; // not dependent on dummy_my_rank

   for (int dummy_my_rank=0;dummy_my_rank<n_proc;dummy_my_rank++)
    {
     // Loop over all backed up elements
     nelem=backed_up_el_pt.size();
     for (unsigned e=0;e<nelem;e++)
      {
       // Get element and its domain
       FiniteElement* el_pt=backed_up_el_pt[e];
       unsigned el_domain=element_domain[e];
     
       // If element is located on current processor add it back to the mesh
       if (el_domain==unsigned(dummy_my_rank))
        {
         // Add element to current processor 
         tmp_element_retained[dummy_my_rank][e]=true;
         number_of_retained_elements[dummy_my_rank]++;
        }
       // Otherwise we may still need it if it's a halo element:
       else
        {       
         // Is one of the nodes associated with the current processor?
         unsigned nnod=el_pt->nnode();
         for (unsigned j=0;j<nnod;j++)
          {
           Node* nod_pt=el_pt->node_pt(j); 

           // Keep element?
           unsigned keep_it=false;
           for (std::set<unsigned>::iterator 
                 it=processors_associated_with_data[nod_pt].begin();
                it!=processors_associated_with_data[nod_pt].end();
                it++)
            {
             if (*it==unsigned(dummy_my_rank))
              {
               keep_it=true;
               break;
              }
            }
         
           // Add a root halo element either if keep_it=true OR this 
           // current mesh has been told to keep all elements as halos,
           // OR the element itself knows that it must be kept
           if ((keep_it) || (keep_all_elements_as_halos())
               || (el_pt->must_be_kept_as_halo()))
            {
             // Add as root halo element whose non-halo counterpart is
             // located on processor el_domain
             tmp_root_halo_element_pt[dummy_my_rank][el_domain].
              push_back(el_pt);
             tmp_element_retained[dummy_my_rank][e]=true;
             number_of_retained_elements[dummy_my_rank]++;
             break;
            }
          }
        }
      }


     // Now loop over all halo elements to check if any of their
     // nodes are located on a higher-numbered processor.
     // The elements associated with these must be added as halo elements, too
     std::map<Node*,bool> higher_numbered_node;
         
     // Loop over all domains
     for (int d=0;d<n_proc;d++)
      {  
       // Loop over root halo elements 
       unsigned nelem=tmp_root_halo_element_pt[dummy_my_rank][d].size();     
       for (unsigned e=0;e<nelem;e++)
        {
         FiniteElement* el_pt=tmp_root_halo_element_pt[dummy_my_rank][d][e];
       
         // Loop over its nodes
         unsigned nnod=el_pt->nnode();
         for (unsigned j=0;j<nnod;j++)
          {
           Node* nod_pt=el_pt->node_pt(j);
           int proc_in_charge=processor_in_charge[nod_pt];
           if (proc_in_charge>d) // too many halos if > changed to >=
            {
             higher_numbered_node[nod_pt]=true;
            }
           else
            {
             higher_numbered_node[nod_pt]=false;
            }
          }
        }
      }

     // Now loop over all the original elements again
     nelem=backed_up_el_pt.size();
     for (unsigned e=0;e<nelem;e++)
      {
       // Get element and its domain
       FiniteElement* el_pt=backed_up_el_pt[e];
       unsigned el_domain=element_domain[e];
     
       // By default, don't keep it
       bool keep_it=false;
     
       // Check if it's already retained
       if (!tmp_element_retained[dummy_my_rank][e])
        {
         // Loop over its nodes
         unsigned nnod=el_pt->nnode();
         for (unsigned j=0;j<nnod;j++)
          {
           Node* nod_pt=el_pt->node_pt(j);
           if (higher_numbered_node[nod_pt]&&
               (element_domain[e]==processor_in_charge[nod_pt]))
            {
             keep_it=true; 
             break; // doesn't help if the break is removed
            }
          }
         if (keep_it)
          {
           tmp_root_halo_element_pt[dummy_my_rank][el_domain].push_back(el_pt);
           tmp_element_retained[dummy_my_rank][e]=true;
           number_of_retained_elements[dummy_my_rank]++;
           number_of_retained_halo_elements++; // need to count these
          }
        }

      }

     if (report_stats)
      {
       // Check number of retained halo elements on this process
       if (number_of_retained_elements[dummy_my_rank]!=0)
        {
         oomph_info << "Percentage of extra halo elements retained: "
                    << 100.0*double(number_of_retained_halo_elements)/
          double(number_of_retained_elements[dummy_my_rank])
                    << " on process " << dummy_my_rank 
                    << " in loop number " << myi << std::endl;
        }
       else // Dummy output in case a process has no retained elements
            // (relevant in some multi-mesh problems)
        {
         oomph_info << "Percentage of extra halo elements retained: "
                    << 0.0 << " on process " << dummy_my_rank
                    << " in loop number " << myi << std::endl;
        }
      }

    } // end of loop over all "processors"; we've now established the
   // elements and the root halo elements for all processors

   int total_number_of_retained_halo_elements=0;
   
   // Sum values over all processes 
   // - must be zero retained in order to continue
   MPI_Allreduce(&number_of_retained_halo_elements, 
                 &total_number_of_retained_halo_elements,1,MPI_INT,
                 MPI_SUM,comm_pt->mpi_comm());

   if (report_stats)
    {
     oomph_info << "Total number of extra halo elements retained: " 
                << total_number_of_retained_halo_elements
                << " in loop: " << myi << std::endl;
    }

   if (total_number_of_retained_halo_elements==0) {elements_retained=false;} 
   
   myi++;

  } // end of while(elements_retained)

 // Copy the elements associated with the actual
 // current processor into its own permanent storage.
 // Do it in the order in which the elements appeared originally
 nelem=backed_up_el_pt.size();
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=backed_up_el_pt[e];
   if (tmp_element_retained[my_rank][e])
    {
     this->add_element_pt(el_pt);
    }
   else
    {
     // Flush the object attached to the tree for this element?
     RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(el_pt);
     if (ref_el_pt!=0)
      {
       ref_el_pt->tree_pt()->flush_object();
      }
     // Delete the element
     delete el_pt;
    }
  }
 
 // Copy the root halo elements associated with the actual
 // current processor into its own permanent storage; the order
 // here is somewhat random but we compensate for that by
 // ensuring that the corresponding haloed elements are 
 // added in the same order below
 for (int d=0;d<n_proc;d++)
  {  
   nelem=tmp_root_halo_element_pt[my_rank][d].size();
   for (unsigned e=0;e<nelem;e++)
    {
     this->add_root_halo_element_pt(d,
      tmp_root_halo_element_pt[my_rank][d][e]);
    }
  }
  
//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Determine root haloed elements....";

 // Determine root haloed elements
 //-------------------------------

 // Loop over all other processors
 for (int d=0;d<n_proc;d++)
  {
   if (d!=my_rank)
    {
     // Loop over root halo elements that are held on that processor
     unsigned nelem_other=tmp_root_halo_element_pt[d][my_rank].size();
     for (unsigned e2=0;e2<nelem_other;e2++)
      {
       // Loop over all elements on current processor.
       // We loop over these in the inner loop) to ensure that they are 
       // added to the haloed lookup scheme in the same
       // order in which elements were added to the
       // corresponding halo lookup scheme.
       nelem=this->nelement();
       for (unsigned e=0;e<nelem;e++)
        {
         // Get pointer to element
         FiniteElement* el_pt=this->finite_element_pt(e);
         
         // Halo elements can't be haloed themselves
         if (!el_pt->is_halo())
          {
           if (el_pt==tmp_root_halo_element_pt[d][my_rank][e2])
            {
             // Current element is haloed by other processor
             this->add_root_haloed_element_pt(d,el_pt);
             break;
            }
          }  
        }
      }
    }
  }

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;


 // Doc stats
 if (report_stats)
  {
   oomph_info << "Processor " << my_rank 
              << " holds " << this->nelement() 
              << " elements of which " << this->nroot_halo_element()
              << " are root halo elements \n while " 
              << this->nroot_haloed_element()
              << " are root haloed elements" << std::endl;
  }

//  oomph_info << "Retain nodes....";
 
 // Loop over all retained elements and mark their nodes
 //-----------------------------------------------------
 // as to be retained too (some double counting going on here)
 //-----------------------------------------------------------
 nelem=this->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=this->finite_element_pt(e);

   // Loop over nodes
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=el_pt->node_pt(j);
     nod_pt->set_non_obsolete();
    }
  }


 // Complete rebuild of mesh by adding retained nodes
 // Note that they are added in the order in which they 
 // occured in the original mesh as this guarantees the
 // synchronisity between the serialised access to halo
 // and haloed nodes from different processors.
 nnod=backed_up_nod_pt.size();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=backed_up_nod_pt[j];
   if(!nod_pt->is_obsolete())
    {
     // Not obsolete so add it back to the mesh
     this->add_node_pt(nod_pt);
    }
  }

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;


 // Prune and rebuild mesh
 //-----------------------

//  oomph_info << "Pruning dead nodes...";

 // Now remove the pruned nodes from the boundary lookup scheme
 this->prune_dead_nodes();

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Setup boundary info....";

 // And finally re-setup the boundary lookup scheme for elements
 this->setup_boundary_element_info();

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Recreate forest....";

 // Re-setup tree forest if needed
 if (this->nelement()>0)
  {
   RefineableMeshBase* ref_mesh_pt=dynamic_cast<RefineableMeshBase*>(this);
   if (ref_mesh_pt!=0)
    {
     ref_mesh_pt->setup_tree_forest();
    }
  }

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Classify nodes....";

 // Classify nodes 
 classify_halo_and_haloed_nodes(comm_pt,doc_info,report_stats);

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

 // Mesh has now been distributed 
 // (required for Z2ErrorEstimator::get_element_errors)
 Mesh_has_been_distributed=true;

 // Doc?
 //-----
 if (doc_info.doc_flag())
  {
   doc_mesh_distribution(comm_pt,doc_info);
  }


}


//========================================================================
/// (Irreversibly) prune halo(ed) elements and nodes, usually
/// after another round of refinement, to get rid of
/// excessively wide halo layers. Note that the current
/// mesh will be now regarded as the base mesh and no unrefinement
/// relative to it will be possible once this function 
/// has been called.
//========================================================================
void Mesh::prune_halo_elements_and_nodes(OomphCommunicator* comm_pt,
                                         DocInfo& doc_info,
                                         const bool& report_stats)
{

 RefineableMeshBase* ref_mesh_pt=dynamic_cast<RefineableMeshBase*>(this);
 if (ref_mesh_pt!=0)
  {

#ifdef OOMPH_HAS_MPI
   // Flush any external element storage before performing the redistribution
   // (in particular, external halo nodes that are on mesh boundaries)
   this->flush_all_external_storage();
#endif
   
   // Storage for number of processors and current processor
   int n_proc=comm_pt->nproc();
   int my_rank=comm_pt->my_rank();
   
   // Doc stats
   if (report_stats)
    {
     oomph_info << "Before pruning: Processor " << my_rank 
                << " holds " << this->nelement() 
                << " elements of which " << this->nroot_halo_element()
                << " are root halo elements \n while " 
                << this->nroot_haloed_element()
                << " are root haloed elements" << std::endl;
    }
   
   // Declare all nodes as obsolete. We'll
   // change this setting for all nodes that must be retained
   // further down
   unsigned nnod=this->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     this->node_pt(j)->set_obsolete();
    }
   
   // Backup pointers to elements in this mesh
   unsigned nelem=this->nelement();
   Vector<FiniteElement*> backed_up_el_pt(nelem);
   std::map<FiniteElement*,bool> keep_element;
   for (unsigned e=0;e<nelem;e++)
    {
     FiniteElement* el_pt=this->finite_element_pt(e);
     backed_up_el_pt[e]=el_pt;
    }
   
   // Get the min and max refinement level, and current refinement pattern
   unsigned min_ref=0;
   unsigned max_ref=0;
   
   // Skip this first bit if you have no elements 
   if (nelem>0)
    {
     // Get min and max refinement level
     ref_mesh_pt->get_refinement_levels(min_ref,max_ref);

     // Get refinement pattern
     Vector<Vector<unsigned> > current_refined;
     ref_mesh_pt->get_refinement_pattern(current_refined);

     // get_refinement_pattern refers to the elements at each level
     // that were refined when proceeding to the next level
     unsigned n_ref=current_refined.size();

     // Loop over all elements; keep those on the min refinement level
     // Need to go back to the level indicated by min_ref
     unsigned base_level=n_ref-(max_ref-min_ref);

#ifdef PARANOID
     if (base_level<0)
      {
       std::ostringstream error_stream;
       error_stream  << "Error: the base level of refinement " 
                     << "is negative; this cannot be the case" << std::endl;
       throw OomphLibError(error_stream.str(),
                           "Mesh::prune_halo_elements_and_nodes(...)",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // Get the elements at the specified "base" refinement level
     Vector<RefineableElement*> base_level_elements_pt;
     ref_mesh_pt->get_elements_at_refinement_level(base_level,
                                                   base_level_elements_pt);
     unsigned n_base_el=base_level_elements_pt.size();

     // Loop over the elements at this level
     for (unsigned e=0;e<n_base_el;e++)
      {
       // Extract correct element...
       RefineableElement* ref_el_pt=base_level_elements_pt[e];

       // Check it exists
       if (ref_el_pt!=0)
        {
         // Keep all non-halo elements, remove excess halos
         if (!ref_el_pt->is_halo())
          {
           keep_element[ref_el_pt]=true;

           // Loop over this non-halo element's nodes and retain them too
           unsigned nnod=ref_el_pt->nnode();
           for (unsigned j=0;j<nnod;j++)
            {
             ref_el_pt->node_pt(j)->set_non_obsolete();
            }
          }
        }
      } // end loop over base level elements
    }

   // Now work on which "root" halo elements to keep at this level
   // Can't use the current set directly; however, 
   // we know the refinement level of the current halo, so
   // it is possible to go from that backwards to find the "father
   // halo element" necessary to complete this step

   // Temp map of vectors holding the pointers to the root halo elements
   std::map<unsigned, Vector<FiniteElement*> > tmp_root_halo_element_pt;

   // Temp map of vectors holding the pointers to the root haloed elements
   std::map<unsigned, Vector<FiniteElement*> > tmp_root_haloed_element_pt;
 
   // Map to store if a halo element survives
   std::map<FiniteElement*,bool> halo_element_is_retained;

   for (int domain=0;domain<n_proc;domain++)
    {
     // Get vector of halo elements with processor domain by copy operation
     Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(domain));
     
     // Loop over halo elements associated with this adjacent domain
     unsigned nelem=halo_elem_pt.size();
     for (unsigned e=0;e<nelem;e++)
      {
       // Get element
       RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>
        (halo_elem_pt[e]);
       
       // An element should only be kept if its refinement
       // level is the same as the minimum refinement level
       unsigned halo_el_level=ref_el_pt->refinement_level();

       RefineableElement* el_pt;
       if (halo_el_level==min_ref)
        {
         // Already at the correct level
         el_pt=ref_el_pt;
        }
       else
        {
         // Need to go up the tree to the father element at min_ref
         RefineableElement* father_el_pt;
         ref_el_pt->get_father_at_refinement_level(min_ref,father_el_pt);
         el_pt=father_el_pt;
        }

       //Loop over nodes
       unsigned nnod=el_pt->nnode();
       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);
         if (!nod_pt->is_obsolete())
          {
           // Keep element and add it to preliminary storage for
           // halo elements associated with current neighbouring domain
           if (!halo_element_is_retained[el_pt])
            {
             keep_element[el_pt]=true;
             tmp_root_halo_element_pt[domain].push_back(el_pt);
             halo_element_is_retained[el_pt]=true;
             break;
            }
          }
        }       
      }

    }

   // Make sure everybody finishes this part
   MPI_Barrier(comm_pt->mpi_comm());

   // Now all processors have decided (independently) which of their
   // (to-be root) halo elements they wish to retain. Now we need to figure out
   // which of their elements are haloed and add them in the appropriate
   // order into the haloed element scheme. For this we exploit that
   // the halo and haloed elements are accessed in the same order on
   // all processors!
   
   // Identify haloed elements on domain d
   for (int d=0;d<n_proc;d++)
    {
     // Loop over domains that halo this domain
     for (int dd=0;dd<n_proc;dd++)
      {       
       // Dont't talk to yourself
       if (d!=dd)
        {
         // If we're identifying my haloed elements:
         if (d==my_rank)
          {
           // Get vector all elements that are currently haloed by domain dd
           Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(dd));
           // Create a vector of ints to indicate if the halo element
           // on processor dd processor was kept
           unsigned nelem=haloed_elem_pt.size();
           Vector<int> halo_kept(nelem);
           
           // Receive this vector from processor dd 
           if (nelem!=0)
            {
             MPI_Status status;
             MPI_Recv(&halo_kept[0],nelem,MPI_INT,dd,0,comm_pt->mpi_comm(),
                      &status);
           
             // Classify haloed element accordingly
             for (unsigned e=0;e<nelem;e++)
              {
               RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>
                (haloed_elem_pt[e]);

               // An element should only be kept if its refinement
               // level is the same as the minimum refinement level
               unsigned haloed_el_level=ref_el_pt->refinement_level();

               // Go up the tree to the correct level
               RefineableElement* el_pt;

               if (haloed_el_level==min_ref)
                {
                 // Already at the correct level
                 el_pt=ref_el_pt;
                }
               else
                {
                 // Need to go up the tree to the father element at min_ref
                 RefineableElement* father_el_pt;
                 ref_el_pt->get_father_at_refinement_level
                  (min_ref,father_el_pt);
                 el_pt=father_el_pt;
                }

               if (halo_kept[e]==1)
                {
                 // I am being haloed by processor dd
                 // Only keep it if it's not already in the storage
                 unsigned n_root_haloed=tmp_root_haloed_element_pt[dd].size();
                 bool already_root_haloed=false;
                 for (unsigned e_root=0;e_root<n_root_haloed;e_root++)
                  {
                   if (el_pt==tmp_root_haloed_element_pt[dd][e_root])
                    {
                     already_root_haloed=true;
                     break;
                    }
                  }
                 if (!already_root_haloed)
                  {
                   tmp_root_haloed_element_pt[dd].push_back(el_pt);
                  }
                }
              }
            }
          }
         else
          {
           // If we're dealing with my halo elements:
           if (dd==my_rank)
            {
             // Find (current) halo elements on processor dd whose non-halo is 
             // on processor d
             Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(d));
             
             // Create a vector of ints to indicate if the halo 
             // element was kept
             unsigned nelem=halo_elem_pt.size();
             Vector<int> halo_kept(nelem,0);
             for (unsigned e=0;e<nelem;e++)
              {
               RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>
                (halo_elem_pt[e]);

               // An element should only be kept if its refinement
               // level is the same as the minimum refinement level
               unsigned halo_el_level=ref_el_pt->refinement_level();

               // Go up the tree to the correct level
               RefineableElement* el_pt;
               if (halo_el_level==min_ref)
                {
                 // Already at the correct level
                 el_pt=ref_el_pt;
                }
               else
                {
                 // Need to go up the tree to the father element at min_ref
                 RefineableElement* father_el_pt;
                 ref_el_pt->get_father_at_refinement_level
                  (min_ref,father_el_pt);
                 el_pt=father_el_pt;
                }

               if (halo_element_is_retained[el_pt])
                {
                 halo_kept[e]=1;
                }
              }
             
             // Now send this vector to processor d to tell it which of
             // the haloed elements (which are listed in the same order)
             // are to be retained as haloed elements.
             if (nelem!=0)
              {
               MPI_Send(&halo_kept[0],nelem,MPI_INT,d,0,comm_pt->mpi_comm());
              }
            }
          }
        }
      }
    }
 
   // Backup pointers to nodes in this mesh
   nnod=this->nnode();
   Vector<Node*> backed_up_nod_pt(nnod);
   for (unsigned j=0;j<nnod;j++)
    {
     backed_up_nod_pt[j]=this->node_pt(j);
    }
   
   // Flush the mesh storage
   this->flush_element_and_node_storage();

   // Loop over all backed-up elements
   nelem=backed_up_el_pt.size();
   for (unsigned e=0;e<nelem;e++)
    {
     RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>
      (backed_up_el_pt[e]);

     // Get refinement level
     unsigned level=ref_el_pt->refinement_level();

     // Go up the tree to the correct level
     RefineableElement* el_pt;

     if (level==min_ref)
      {
       // Already at the correct level
       el_pt=ref_el_pt;
      }
     else
      {
       // Need to go up the tree to the father element at min_ref
       RefineableElement* father_el_pt;
       ref_el_pt->get_father_at_refinement_level
        (min_ref,father_el_pt);
       el_pt=father_el_pt;
      }     

     // If the base element is going to be kept, then add the current element
     // to the "new" mesh
     if (keep_element[el_pt])
      {
       this->add_element_pt(ref_el_pt);
      }
     else
      {
       // Flush the object attached to the tree for this element?
       RefineableElement* my_el_pt=dynamic_cast<RefineableElement*>(ref_el_pt);
       if (my_el_pt!=0)
        {
         my_el_pt->tree_pt()->flush_object();
        }

       // Delete the element
       delete ref_el_pt;
      }
    }

   // Wipe the storage scheme for (root) halo(ed) elements and then re-assign
   Root_haloed_element_pt.clear();
   Root_halo_element_pt.clear();     
   for (int domain=0;domain<n_proc;domain++)
    {
     
     unsigned nelem=tmp_root_halo_element_pt[domain].size();
     for (unsigned e=0;e<nelem;e++)
      {
       Root_halo_element_pt[domain].push_back(
        tmp_root_halo_element_pt[domain][e]);
      }
    
     nelem=tmp_root_haloed_element_pt[domain].size();
     for (unsigned e=0;e<nelem;e++)
      {
       Root_haloed_element_pt[domain].push_back(
        tmp_root_haloed_element_pt[domain][e]);
      }
    }
   
   // Doc stats
   if (report_stats)
    {
     oomph_info << "AFTER pruning:  Processor " << my_rank 
                << " holds " << this->nelement() 
                << " elements of which " << this->nroot_halo_element()
                << " are root halo elements \n while " 
                << this->nroot_haloed_element()
                << " are root haloed elements" << std::endl;
    }

   // Loop over all retained elements at this level and mark their nodes
   //-------------------------------------------------------------------
   // as to be retained too (some double counting going on here)
   //-----------------------------------------------------------
   nelem=this->nelement();
   for (unsigned e=0;e<nelem;e++)
    {
     FiniteElement* el_pt=this->finite_element_pt(e);
     
     // Loop over nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=el_pt->node_pt(j);
       nod_pt->set_non_obsolete();
      }
    }
   
   // Complete rebuild of mesh by adding retained nodes
   // Note that they are added in the order in which they 
   // occured in the original mesh as this guarantees the
   // synchronisity between the serialised access to halo
   // and haloed nodes from different processors.
   nnod=backed_up_nod_pt.size();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=backed_up_nod_pt[j];
     if(!nod_pt->is_obsolete())
      {
       // Not obsolete so add it back to the mesh
       this->add_node_pt(nod_pt);
      }
    }
   
   // Prune and rebuild mesh
   //-----------------------
   
   // Now remove the pruned nodes from the boundary lookup scheme
   this->prune_dead_nodes();
    
   // And finally re-setup the boundary lookup scheme for elements
   this->setup_boundary_element_info();

   // Re-setup tree forest if needed
   if (this->nelement()>0)
    {
     RefineableMeshBase* ref_mesh_pt=dynamic_cast<RefineableMeshBase*>(this);
     if (ref_mesh_pt!=0)
      {
       ref_mesh_pt->setup_tree_forest();
      }
    }

   // Classify nodes 
   classify_halo_and_haloed_nodes(comm_pt,doc_info,report_stats);
   
   // Doc?
   //-----
   if (doc_info.doc_flag())
    {
     doc_mesh_distribution(comm_pt,doc_info);
    }

  }

}






//========================================================================
///  Get efficiency of mesh distribution: In an ideal distribution
/// without halo overhead, each processor would only hold its own
/// elements. Efficieny per processor =  (number of non-halo elements)/
/// (total number of elements). 
//========================================================================
void Mesh::get_efficiency_of_mesh_distribution(OomphCommunicator* comm_pt,
                                               double& av_efficiency,
                                               double& max_efficiency,
                                               double& min_efficiency)
{
 // Storage for number of processors and current processor
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();

 // Create vector to hold number of elements and halo elements
 Vector<int> nhalo_elements(n_proc);
 Vector<int> n_elements(n_proc);

 // Count total number of halo elements
 unsigned count=0;
 for (int d=0;d<n_proc;d++)
  {
   Vector<FiniteElement*> halo_elem_pt(halo_element_pt(d));
   count+=halo_elem_pt.size();
  }

 // Stick own number into appropriate entry
 nhalo_elements[my_rank]=count;
 n_elements[my_rank]=nelement();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of elements on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhalo_elements[my_rank],1,MPI_INT,
            &nhalo_elements[0],1, MPI_INT,
            0,comm_pt->mpi_comm());
 MPI_Gather(&n_elements[my_rank],1,MPI_INT,
            &n_elements[0],1, MPI_INT,
            0,comm_pt->mpi_comm());

 // Initialise stats
 av_efficiency=0.0;
 double max=-1.0;
 double min=1000000000.0;

 if (my_rank==0)
  {
   for (int i=0;i<n_proc;i++)
    {
     double eff=double(n_elements[i]-nhalo_elements[i])/double(n_elements[i]);
     av_efficiency+=eff;
     if (eff>max) max=eff;
     if (eff<min) min=eff;
       
    }
   av_efficiency/=double(n_proc);
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_DOUBLE,0,comm_pt->mpi_comm());
 MPI_Bcast(&min,1,MPI_DOUBLE,0,comm_pt->mpi_comm());
 MPI_Bcast(&av_efficiency,1,MPI_DOUBLE,0,comm_pt->mpi_comm());

 max_efficiency=max;
 min_efficiency=min;

}



//========================================================================
/// Doc the mesh distribution
//========================================================================
void Mesh::doc_mesh_distribution(OomphCommunicator* comm_pt,DocInfo& doc_info)
{ 
 // Storage for current processor and number of processors
 int my_rank=comm_pt->my_rank();
 int n_proc=comm_pt->nproc();

 char filename[100];
 std::ofstream some_file;

 // Doc elements on this processor
 sprintf(filename,"%s/elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 this->output(some_file,5);
 some_file.close();

 // Doc non-halo elements on this processor
 sprintf(filename,"%s/non_halo_elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);

 // Get to elements on processor
 unsigned nelem=this->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=this->finite_element_pt(e);

   if (!el_pt->is_halo()) // output if non-halo
    {
     el_pt->output(some_file,5);
    }
  }

 some_file.close();
 
 // Doc halo elements on this processor
 sprintf(filename,"%s/halo_elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 for (int domain=0; domain<n_proc; domain++)
  {
   // Get vector of halo elements by copy operation
   Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(domain));
   unsigned nelem=halo_elem_pt.size();
//   oomph_info
//    << "Processor " << my_rank << " holds " << nelem 
//    << " halo elem whose non-halo counterparts are located on domain "
//    << domain << std::endl;
   for (unsigned e=0;e<nelem;e++)
    {
     halo_elem_pt[e]->output(some_file,5);
    }
  }
 some_file.close();
 
 
 // Doc haloed elements on this processor
 sprintf(filename,"%s/haloed_elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 for (int domain=0; domain<n_proc; domain++)
  {
   // Get vector of haloed elements by copy operation
   Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(domain));
   unsigned nelem=haloed_elem_pt.size();
//   oomph_info
//    << "Processor " << my_rank << " holds " << nelem
//    << " haloed elem whose halo counterparts are located on domain "
//    << domain << std::endl;
   for (unsigned e=0;e<nelem;e++)
    {
     haloed_elem_pt[e]->output(some_file,5);
    }
  }
 some_file.close();
 
 
 // Doc nodes on this processor
 sprintf(filename,"%s/nodes_on_proc%i_%i.dat",doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=this->node_pt(j);
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt==0) // not a SolidNode (see comment below)
    {
     unsigned ndim=nod_pt->ndim();
     for (unsigned i=0;i<ndim;i++)
      {
       some_file << nod_pt->x(i) << " " ;
      }
     some_file
//      << nod_pt->processor_in_charge() << " " 
      << nod_pt->is_halo() << " " 
      << nod_pt->eqn_number(0) << " "   // these two won't work for SolidNodes
      << nod_pt->is_pinned(0) << " "    // with eqn numbers for position only
      << std::endl;
    }
  }
 some_file.close();
 
 // Doc solid nodes on this processor
 sprintf(filename,"%s/solid_nodes_on_proc%i_%i.dat",
         doc_info.directory().c_str(),my_rank,doc_info.number());
 some_file.open(filename);
 unsigned nsnod=this->nnode();
 for (unsigned j=0;j<nsnod;j++)
  {
   Node* nod_pt=this->node_pt(j);
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt!=0)
    {
     unsigned ndim=solid_nod_pt->ndim();
     for (unsigned i=0;i<ndim;i++)
      {
       some_file << nod_pt->x(i) << " " ;
      }
     some_file
//      << solid_nod_pt->processor_in_charge() << " " 
      << solid_nod_pt->is_halo() << " ";
     unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
     for (unsigned ival=0;ival<nval;ival++)
      {
       some_file
        << solid_nod_pt->variable_position_pt()->eqn_number(ival) << " "
        << solid_nod_pt->variable_position_pt()->is_pinned(ival) << " ";
      }
     some_file << std::endl;
    }
  }
 some_file.close();
 
 
 // Doc halo nodes on this processor
 sprintf(filename,"%s/halo_nodes_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 for (int domain=0; domain<n_proc; domain++)
  {
   unsigned nnod=this->nhalo_node(domain);
//   oomph_info << "Processor " << my_rank 
//             << " holds " << nnod << " halo nodes assoc with domain "
//              << domain << std::endl;
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=this->halo_node_pt(domain,j);
     unsigned ndim=nod_pt->ndim();
     for (unsigned i=0;i<ndim;i++)
      {
       some_file << nod_pt->x(i) << " " ;
      }
     some_file << domain << std::endl;
    }
  }
 some_file.close();
 
 
 
 // Doc haloed nodes on this processor
 sprintf(filename,"%s/haloed_nodes_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 for (int domain=0; domain<n_proc; domain++)
  {
   unsigned nnod=this->nhaloed_node(domain);
//   oomph_info << "Processor " << my_rank 
//              << " holds " << nnod << " haloed nodes assoc with domain "
//              << domain << std::endl;
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=this->haloed_node_pt(domain,j);
     unsigned ndim=nod_pt->ndim();
     for (unsigned i=0;i<ndim;i++)
      {
       some_file << nod_pt->x(i) << " " ;
      }
     some_file << domain << std::endl;
    }
  }
 some_file.close();
 
 
 // Doc mesh
 sprintf(filename,"%s/mesh%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 this->output(some_file,5);
 some_file.close();
 
 
 // Doc boundary scheme
 sprintf(filename,"%s/boundaries%i_%i.dat",
         doc_info.directory().c_str(),
         my_rank,doc_info.number());
 some_file.open(filename);
 this->output_boundaries(some_file);
 some_file.close();
 
 
 // Doc elements next to boundaries scheme
 
 // How many finite elements are adjacent to boundary b?
 unsigned nbound=this->nboundary();
 for (unsigned b=0;b<nbound;b++)
  {
   sprintf(filename,"%s/boundary_elements%i_%i_%i.dat",
           doc_info.directory().c_str(),
           my_rank,b,doc_info.number());
   some_file.open(filename);
   unsigned nelem=this->nboundary_element(b);
   for (unsigned e=0;e<nelem;e++)
    {
     this->boundary_element_pt(b,e)->output(some_file,5);
    }
   some_file.close();
  }
 
}


//========================================================================
/// Check the halo/haloed/shared node/element schemes on the Mesh
//========================================================================
void Mesh::check_halo_schemes(OomphCommunicator* comm_pt, DocInfo& doc_info,
                              double& max_permitted_error_for_halo_check)
{
 // Moved this from the Problem class so that it would work better
 // in multiple mesh problems; there remains a simple "wrapper"
 // function in the Problem class that calls this for each (sub)mesh.

 MPI_Status status;
 char filename[100];
 std::ofstream shared_file;
 std::ofstream halo_file;
 std::ofstream haloed_file;

 // Storage for current processor and number of processors
 int n_proc=comm_pt->nproc();
 int my_rank=comm_pt->my_rank();

 // Check the shared node scheme first: if this is incorrect then
 // the halo(ed) node scheme is likely to be wrong too

 // Doc shared nodes lookup schemes
 //-------------------------------------
 if (doc_info.doc_flag())
  {
   // Loop over domains for shared nodes
   for (int dd=0;dd<n_proc;dd++)
    {   
     sprintf(filename,"%s/shared_node_check%i_%i.dat",
             doc_info.directory().c_str(),my_rank,dd);
     shared_file.open(filename);
     shared_file << "ZONE " << std::endl;
     
     unsigned nnod=nshared_node(dd);
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=shared_node_pt(dd,j);
       unsigned ndim=nod_pt->ndim();
       for (unsigned i=0;i<ndim;i++)
        {
         shared_file << nod_pt->position(i) << " ";
        }
       shared_file << std::endl;
      }
     // Dummy output for processor that doesn't share nodes
     // (needed for tecplot)
     if ((nnod==0) && (nelement()!=0))
      {
       unsigned ndim=finite_element_pt(0)->node_pt(0)->ndim();
       if (ndim==2)
        {
         shared_file   << " 1.0 1.1 " << std::endl;
        }
       else
        {
         shared_file   << " 1.0 1.1 1.1" << std::endl;
        }
      }
     shared_file.close(); 
    }

  }

 // Check shared nodes lookup schemes
 //---------------------------------------
 double max_error=0.0;

 // Loop over domains for shared nodes
 for (int d=0;d<n_proc;d++)
  {
   // Are my shared nodes being checked?
   if (d==my_rank)
    {
     // Loop over domains for shared nodes
     for (int dd=0;dd<n_proc;dd++)
      {
       // Don't talk to yourself
       if (dd!=d)
        {
         // How many of my nodes are shared nodes with processor dd?
         int nnod_shared=nshared_node(dd);

         if (nnod_shared!=0)
          {         
           // Receive from processor dd how many of his nodes are shared
           // with this processor
           int nnod_share=0;
           MPI_Recv(&nnod_share,1, MPI_INT,dd,0,comm_pt->mpi_comm(),&status);
         
           if (nnod_shared!=nnod_share)
            {
             std::ostringstream error_message;
           
             error_message
              << "Clash in numbers of shared nodes! " 
              << std::endl;
             error_message 
              << "# of shared nodes on proc "
              << dd << ": " << nnod_shared << std::endl;
             error_message
              << "# of shared nodes on proc "
              << d << ": " << nnod_share << std::endl;
             error_message 
              << "(Re-)run Problem::check_halo_schemes() with DocInfo object"
              << std::endl;
             error_message 
              << "to identify the problem" << std::endl;
             throw OomphLibError(error_message.str(),
                                 "Mesh::check_halo_schemes()",
                                 OOMPH_EXCEPTION_LOCATION);
            }


           unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
         
           // Get strung-together nodal positions from other processor
           Vector<double> other_nodal_positions(nod_dim*nnod_share);
           MPI_Recv(&other_nodal_positions[0],nod_dim*nnod_share,MPI_DOUBLE,dd,
                    0,comm_pt->mpi_comm(),&status);

           // Check
           unsigned count=0;
           for (int j=0;j<nnod_share;j++)
            {
             double x_shared=shared_node_pt(dd,j)->position(0);
             double y_shared=shared_node_pt(dd,j)->position(1);
             double z_shared=0.0;
             if (nod_dim==3)
              {
               z_shared=shared_node_pt(dd,j)->position(2);
              }
             double x_share=other_nodal_positions[count];
             count++;
             double y_share=other_nodal_positions[count];
             count++;
             double z_share=0.0;
             if (nod_dim==3)
              {
               z_share=other_nodal_positions[count];
               count++;
              }
             double error=sqrt( pow(x_shared-x_share,2)+
                                pow(y_shared-y_share,2)+
                                pow(z_shared-z_share,2));
             if (fabs(error)>max_error)
              {
//              std::cout << "ZONE" << std::endl;
//              std::cout << x_halo << " " 
//                        << y_halo << " " 
//                        << y_halo << " " 
//                        << d << " " << dd 
//                        << std::endl;
//              std::cout << x_haloed << " " 
//                        << y_haloed << " " 
//                        << y_haloed << " "
//                        << d << " " << dd  
//                        << std::endl;
//              std::cout << std::endl;
               max_error=fabs(error);       
              }
            }
          }
        }
      }
    }
   // My shared nodes are not being checked: Send my shared nodes
   // to the other processor
   else
    {
     int nnod_share=nshared_node(d);
     
     if (nnod_share!=0)
      {
       // Send it across to the processor whose shared nodes are being checked
       MPI_Send(&nnod_share,1,MPI_INT,d,0,comm_pt->mpi_comm());

       unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
         
       // Now string together the nodal positions of all shared nodes
       Vector<double> nodal_positions(nod_dim*nnod_share);
       unsigned count=0;
       for (int j=0;j<nnod_share;j++)
        {
         nodal_positions[count]=shared_node_pt(d,j)->position(0);
         count++;
         nodal_positions[count]=shared_node_pt(d,j)->position(1);
         count++;
         if (nod_dim==3)
          {
           nodal_positions[count]=shared_node_pt(d,j)->position(2);
           count++;
          }
        }
       // Send it across to the processor whose shared nodes are being checked
       MPI_Send(&nodal_positions[0],nod_dim*nnod_share,MPI_DOUBLE,d,0,
                comm_pt->mpi_comm());
      }
    }
  }

 oomph_info << "Max. error for shared nodes " << max_error
            << std::endl;

 if (max_error>max_permitted_error_for_halo_check)
  {         
   std::ostringstream error_message;
   error_message
    << "This is bigger than the permitted threshold "
    << max_permitted_error_for_halo_check << std::endl;
   error_message
    << "If you believe this to be acceptable for your problem\n"
    << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
   throw OomphLibError(error_message.str(),
                       "Mesh::check_halo_schemes()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Now check the halo/haloed element lookup scheme

 // Doc halo/haoloed element lookup schemes
 //-----------------------------------------
 if (doc_info.doc_flag())
  {
   // Loop over domains for halo elements
   for (int dd=0;dd<n_proc;dd++)
    {
     sprintf(filename,"%s/halo_element_check%i_%i_mesh_%i.dat",
             doc_info.directory().c_str(),my_rank,dd,doc_info.number());
     halo_file.open(filename);
     
     // Get vectors of halo/haloed elements by copy operation
     Vector<FiniteElement*> 
      halo_elem_pt(halo_element_pt(dd));
     
     unsigned nelem=halo_elem_pt.size();

     for (unsigned e=0;e<nelem;e++)
      {
       halo_file << "ZONE " << std::endl;
       unsigned nnod=halo_elem_pt[e]->nnode();
       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=halo_elem_pt[e]->node_pt(j);
         unsigned ndim=nod_pt->ndim();
         for (unsigned i=0;i<ndim;i++)
          {
           halo_file << nod_pt->position(i) << " ";
          }
         halo_file << std::endl;
        }
      }
     halo_file.close(); 
    }
      
   // Loop over domains for halo elements
   for (int d=0;d<n_proc;d++)
    {
     sprintf(filename,"%s/haloed_element_check%i_%i_mesh_%i.dat",
             doc_info.directory().c_str(),d,my_rank,doc_info.number());
     haloed_file.open(filename);
     
     // Get vectors of halo/haloed elements by copy operation
     Vector<FiniteElement*> 
      haloed_elem_pt(haloed_element_pt(d));
     
     unsigned nelem2=haloed_elem_pt.size(); 
     for (unsigned e=0;e<nelem2;e++)
      {
       haloed_file << "ZONE " << std::endl;
       unsigned nnod2=haloed_elem_pt[e]->nnode();
       for (unsigned j=0;j<nnod2;j++)
        {
         Node* nod_pt=haloed_elem_pt[e]->node_pt(j);
         unsigned ndim=nod_pt->ndim();
         for (unsigned i=0;i<ndim;i++)
          {
           haloed_file << nod_pt->position(i) << " ";
          }
         haloed_file << std::endl;
        }
      }
     haloed_file.close(); 
    }
  }

 // Check halo/haloed element lookup schemes
 //-----------------------------------------
 max_error=0.0;

 // Loop over domains for haloed elements
 for (int d=0;d<n_proc;d++)
  {
   // Are my haloed elements being checked?
   if (d==my_rank)
    {
     // Loop over domains for halo elements
     for (int dd=0;dd<n_proc;dd++)
      {
       // Don't talk to yourself
       if (dd!=d)
        {
         // Get vectors of haloed elements by copy operation
         Vector<FiniteElement*> 
          haloed_elem_pt(haloed_element_pt(dd));
         
         // How many of my elements are haloed elements whose halo
         // counterpart is located on processor dd?
         int nelem_haloed=haloed_elem_pt.size();

         if (nelem_haloed!=0)
          {
           // Receive from processor dd how many of his elements are halo
           // nodes whose non-halo counterparts are located here
           int nelem_halo=0;
           MPI_Recv(&nelem_halo,1,MPI_INT,dd,0,comm_pt->mpi_comm(),&status);
           if (nelem_halo!=nelem_haloed)
            {
             std::ostringstream error_message;
             error_message 
              << "Clash in numbers of halo and haloed elements! " 
              << std::endl;           
             error_message 
              << "# of haloed elements whose halo counterpart lives on proc "
              << dd << ": " << nelem_haloed << std::endl;
             error_message
              << "# of halo elements whose non-halo counterpart lives on proc "
              << d << ": " << nelem_halo << std::endl;
             error_message 
              << "(Re-)run Problem::check_halo_schemes() with DocInfo object"
              << std::endl;
             error_message 
              << "to identify the problem" << std::endl;
             throw OomphLibError(error_message.str(),
                                 "Mesh::check_halo_schemes()",
                                 OOMPH_EXCEPTION_LOCATION);
            }


           // Get strung-together elemental nodal positions 
           // from other processor
           unsigned nnod_per_el=finite_element_pt(0)->nnode();
           unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
           Vector<double> other_nodal_positions
            (nod_dim*nnod_per_el*nelem_halo);
           Vector<int> other_nodal_hangings(nnod_per_el*nelem_halo);
           MPI_Recv(&other_nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
                    MPI_DOUBLE,dd,0,comm_pt->mpi_comm(),&status);
           MPI_Recv(&other_nodal_hangings[0],nnod_per_el*nelem_halo,MPI_INT,dd,
                    1,comm_pt->mpi_comm(),&status);

//         oomph_info << "Received from process " << dd 
//            << ", with size=" << nod_dim*nnod_per_el*nelem_halo << std::endl;
         
           sprintf(filename,"%s/error_haloed_check%i_%i.dat",
                   doc_info.directory().c_str(),dd,my_rank);
           haloed_file.open(filename);
           sprintf(filename,"%s/error_halo_check%i_%i.dat",
                   doc_info.directory().c_str(),dd,my_rank);
           halo_file.open(filename);
     
           unsigned count=0;         
           unsigned count_hanging=0;
           for (int e=0;e<nelem_haloed;e++)
            {   
             for (unsigned j=0;j<nnod_per_el;j++)
              {
               // Testing POSITIONS, not x location 
               // (cf hanging nodes, nodes.h)
               double x_haloed=haloed_elem_pt[e]->node_pt(j)->position(0);
               double y_haloed=haloed_elem_pt[e]->node_pt(j)->position(1);
               double z_haloed=0.0;
               if (nod_dim==3)
                {
                 z_haloed=haloed_elem_pt[e]->node_pt(j)->position(2);
                }
               double x_halo=other_nodal_positions[count];
               count++;
               double y_halo=other_nodal_positions[count];
               count++;
               int other_hanging=other_nodal_hangings[count_hanging];
               count_hanging++;
               double z_halo=0.0;
               if (nod_dim==3)
                {
                 z_halo=other_nodal_positions[count];
                 count++;
                }
               double error=sqrt( pow(x_haloed-x_halo,2)+ 
                                  pow(y_haloed-y_halo,2)+
                                  pow(z_haloed-z_halo,2));
               if (fabs(error)>max_error) max_error=fabs(error);
               if (fabs(error)>0.0)
                {
                 // Report error. NOTE: ERROR IS THROWN BELOW ONCE 
                 // ALL THIS HAS BEEN PROCESSED.
                 oomph_info
                  << "Discrepancy between nodal coordinates of halo(ed)"
                  << "element.  Error: " << error << std::endl;
                 oomph_info
                  << "Domain with non-halo (i.e. haloed) elem: " 
                  << dd << std::endl;
                 oomph_info
                  << "Domain with    halo                elem: " << d
                  << std::endl;
                 oomph_info
                  << "Current processor is " << my_rank 
                  << std::endl
                  << "Nodal positions: " << x_halo << " " << y_halo 
                  << std::endl
                  << "and haloed: " << x_haloed << " " << y_haloed << std::endl
                  << "Node pointer: " << haloed_elem_pt[e]->node_pt(j)
                  << std::endl;
//               oomph_info << "Haloed: " << x_haloed << " " << y_haloed << " "
//                          << error << " " << my_rank << " "
//                          << dd << std::endl;
//               oomph_info << "Halo: " << x_halo << " " << y_halo << " "
//                          << error << " " << my_rank << " "
//                          << dd << std::endl;
                 haloed_file << x_haloed << " " << y_haloed << " "
                             << error << " " << my_rank << " " 
                             << dd << " "
                             << haloed_elem_pt[e]->node_pt(j)->is_hanging() 
                             << std::endl;
                 halo_file << x_halo << " " << y_halo << " "
                           << error << " " << my_rank << " " 
                           << dd << " " 
                           << other_hanging << std::endl; 
                 // (communicated is_hanging value)
                }
              } // j<nnod_per_el
            } // e<nelem_haloed
//         oomph_info << "Check count (receive)... " << count << std::endl;
           haloed_file.close();
           halo_file.close();  
          }
        }
      }
    }
   // My haloed elements are not being checked: Send my halo elements
   // whose non-halo counterparts are located on processor d
   else
    {

     // Get vectors of halo elements by copy operation
     Vector<FiniteElement*> 
      halo_elem_pt(halo_element_pt(d));
     
     // How many of my elements are halo elements whose non-halo
     // counterpart is located on processor d?
     unsigned nelem_halo=halo_elem_pt.size();

     if (nelem_halo!=0)
      {          
       // Send it across to the processor whose haloed nodes are being checked
       MPI_Send(&nelem_halo,1,MPI_INT,d,0,comm_pt->mpi_comm());


       // Now string together the nodal positions of all halo nodes
       unsigned nnod_per_el=finite_element_pt(0)->nnode();
       unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
       Vector<double> nodal_positions(nod_dim*nnod_per_el*nelem_halo); 
       Vector<int> nodal_hangings(nnod_per_el*nelem_halo);
       unsigned count=0;
       unsigned count_hanging=0;
       for (unsigned e=0;e<nelem_halo;e++)
        {
         FiniteElement* el_pt= halo_elem_pt[e];
         for (unsigned j=0;j<nnod_per_el;j++)
          {
           // Testing POSITIONS, not x location (cf hanging nodes, nodes.h)
           nodal_positions[count]=el_pt->node_pt(j)->position(0);
           count++;
           nodal_positions[count]=el_pt->node_pt(j)->position(1);
           count++;
           if (el_pt->node_pt(j)->is_hanging())
            {
             nodal_hangings[count_hanging]=1;
            }
           else
            {
             nodal_hangings[count_hanging]=0;
            }
           count_hanging++;
           if (nod_dim==3)
            {
             nodal_positions[count]=el_pt->node_pt(j)->position(2);
             count++;
            }
          }
        }
       // Send it across to the processor whose haloed elements are being 
       // checked

       MPI_Send(&nodal_positions[0],nod_dim*nnod_per_el*nelem_halo,
                MPI_DOUBLE,d,0,comm_pt->mpi_comm());
       MPI_Send(&nodal_hangings[0],nnod_per_el*nelem_halo,
                MPI_INT,d,1,comm_pt->mpi_comm());
      }
    }
  }

 oomph_info << "Max. error for halo/haloed elements " << max_error
            << std::endl;

 if (max_error>max_permitted_error_for_halo_check)
  {         
   std::ostringstream error_message;
   error_message
    << "This is bigger than the permitted threshold "
    << max_permitted_error_for_halo_check << std::endl;
   error_message
    << "If you believe this to be acceptable for your problem\n"
    << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
   throw OomphLibError(error_message.str(),
                       "Mesh::check_halo_schemes()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Now check the halo/haloed nodes lookup schemes
 
 // Doc halo/haloed nodes lookup schemes
 //-------------------------------------
 if (doc_info.doc_flag())
  {
   // Loop over domains for halo nodes
   for (int dd=0;dd<n_proc;dd++)
    {   
     sprintf(filename,"%s/halo_node_check%i_%i_mesh_%i.dat",
             doc_info.directory().c_str(),my_rank,dd,doc_info.number());
     halo_file.open(filename);
     halo_file << "ZONE " << std::endl;
     
     unsigned nnod=nhalo_node(dd);
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=halo_node_pt(dd,j);
       unsigned ndim=nod_pt->ndim();
       for (unsigned i=0;i<ndim;i++)
        {
         halo_file << nod_pt->position(i) << " ";
        }
       halo_file << nod_pt->is_hanging() << std::endl;
      }
     // Dummy output for processor that doesn't share halo nodes
     // (needed for tecplot)
     // (This will only work if there are elements on this processor...)
     if ((nnod==0) && (nelement()!=0))
      {
       unsigned ndim=finite_element_pt(0)->node_pt(0)->ndim();
       if (ndim==2)
        {
         halo_file   << " 1.0 1.1 " << std::endl;
        }
       else
        {
         halo_file   << " 1.0 1.1 1.1" << std::endl;
        }
      }
     halo_file.close(); 
    }
   
   
   // Loop over domains for haloed nodes
   for (int d=0;d<n_proc;d++)
    {
     sprintf(filename,"%s/haloed_node_check%i_%i_mesh_%i.dat",
             doc_info.directory().c_str(),d,my_rank,doc_info.number());
     haloed_file.open(filename);
     haloed_file << "ZONE " << std::endl;
     
     unsigned nnod=nhaloed_node(d);
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=haloed_node_pt(d,j);
       unsigned ndim=nod_pt->ndim();
       for (unsigned i=0;i<ndim;i++)
        {
         haloed_file << nod_pt->position(i) << " ";
        }
       haloed_file << nod_pt->is_hanging() << std::endl;
      }
     // Dummy output for processor that doesn't share halo nodes
     // (needed for tecplot)
     if ((nnod==0) && (nelement()!=0))
      {
       unsigned ndim=finite_element_pt(0)->node_pt(0)->ndim();
       if (ndim==2)
        {
         halo_file   << " 1.0 1.1 " << std::endl;
        }
       else
        {
         halo_file   << " 1.0 1.1 1.1" << std::endl;
        }
      }
     haloed_file.close(); 
    }
  }

 // Check halo/haloed nodes lookup schemes
 //---------------------------------------
 max_error=0.0;

 // Loop over domains for haloed nodes
 for (int d=0;d<n_proc;d++)
  {
   // Are my haloed nodes being checked?
   if (d==my_rank)
    {
     // Loop over domains for halo nodes
     for (int dd=0;dd<n_proc;dd++)
      {
       // Don't talk to yourself
       if (dd!=d)
        {
         // How many of my nodes are haloed nodes whose halo
         // counterpart is located on processor dd?
         int nnod_haloed=nhaloed_node(dd);

         if (nnod_haloed!=0)
          {         
           // Receive from processor dd how many of his nodes are halo
           // nodes whose non-halo counterparts are located here
           int nnod_halo=0;
           MPI_Recv(&nnod_halo,1,MPI_INT,dd,0,comm_pt->mpi_comm(),&status);
         
           if (nnod_haloed!=nnod_halo)
            {
             std::ostringstream error_message;
           
             error_message
              << "Clash in numbers of halo and haloed nodes! " 
              << std::endl;
             error_message 
              << "# of haloed nodes whose halo counterpart lives on proc "
              << dd << ": " << nnod_haloed << std::endl;
             error_message
              << "# of halo nodes whose non-halo counterpart lives on proc "
              << d << ": " << nnod_halo << std::endl;
             error_message 
              << "(Re-)run Mesh::check_halo_schemes() with DocInfo object"
              << std::endl;
             error_message 
              << "to identify the problem" << std::endl;
             throw OomphLibError(error_message.str(),
                                 "Mesh::check_halo_schemes()",
                                 OOMPH_EXCEPTION_LOCATION);
            }


           unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
         
           // Get strung-together nodal positions from other processor
           Vector<double> other_nodal_positions(nod_dim*nnod_halo);
           MPI_Recv(&other_nodal_positions[0],nod_dim*nnod_halo,MPI_DOUBLE,dd,
                    0,comm_pt->mpi_comm(),&status);

           // Check
           unsigned count=0;
           for (int j=0;j<nnod_halo;j++)
            {
             double x_haloed=haloed_node_pt(dd,j)->position(0);
             double y_haloed=haloed_node_pt(dd,j)->position(1);
             double z_haloed=0.0;
             if (nod_dim==3)
              {
               z_haloed=haloed_node_pt(dd,j)->position(2);
              }
             double x_halo=other_nodal_positions[count];
             count++;
             double y_halo=other_nodal_positions[count];
             count++;
             double z_halo=0.0;
             if (nod_dim==3)
              {
               z_halo=other_nodal_positions[count];
               count++;
              }
             double error=sqrt( pow(x_haloed-x_halo,2)+
                                pow(y_haloed-y_halo,2)+
                                pow(z_haloed-z_halo,2));
             if (fabs(error)>max_error)
              {
//              std::cout << "ZONE" << std::endl;
//              std::cout << x_halo << " " 
//                        << y_halo << " " 
//                        << y_halo << " " 
//                        << d << " " << dd 
//                        << std::endl;
//              std::cout << x_haloed << " " 
//                        << y_haloed << " " 
//                        << y_haloed << " "
//                        << d << " " << dd  
//                        << std::endl;
//              std::cout << std::endl;
               max_error=fabs(error);       
              }
            }
          }
        }
      }
    }
   // My haloed nodes are not being checked: Send my halo nodes
   // whose non-halo counterparts are located on processor d
   else
    {
     int nnod_halo=nhalo_node(d);

     if (nnod_halo!=0)
      {     
       // Send it across to the processor whose haloed nodes are being checked
       MPI_Send(&nnod_halo,1,MPI_INT,d,0,comm_pt->mpi_comm());

       unsigned nod_dim=finite_element_pt(0)->node_pt(0)->ndim();
         
       // Now string together the nodal positions of all halo nodes
       Vector<double> nodal_positions(nod_dim*nnod_halo);
       unsigned count=0;
       for (int j=0;j<nnod_halo;j++)
        {
         nodal_positions[count]=halo_node_pt(d,j)->position(0);
         count++;
         nodal_positions[count]=halo_node_pt(d,j)->position(1);
         count++;
         if (nod_dim==3)
          {
           nodal_positions[count]=halo_node_pt(d,j)->position(2);
           count++;
          }
        }
       // Send it across to the processor whose haloed nodes are being checked
       MPI_Send(&nodal_positions[0],nod_dim*nnod_halo,MPI_DOUBLE,d,0,
                comm_pt->mpi_comm());
      }
    }
  }

 oomph_info << "Max. error for halo/haloed nodes " << max_error
            << std::endl;

 if (max_error>max_permitted_error_for_halo_check)
  {         
   std::ostringstream error_message;
   error_message
    << "This is bigger than the permitted threshold "
    << max_permitted_error_for_halo_check << std::endl;
   error_message
    << "If you believe this to be acceptable for your problem\n"
    << "increase Problem::Max_permitted_error_for_halo_check and re-run \n";
   throw OomphLibError(error_message.str(),
                       "Mesh::check_halo_schemes()",
                       OOMPH_EXCEPTION_LOCATION);
  }

}


#endif

//========================================================================
/// Wipe the storage for all externally-based elements and delete halos
//========================================================================
void Mesh::flush_all_external_storage()
{
 // Clear the local storage for elements and nodes
 External_element_pt.clear();
 External_node_pt.clear();

#ifdef OOMPH_HAS_MPI
 // Storage for number of processors - use size of external halo node array
 int n_proc=External_halo_node_pt.size();

 // Careful: some of the external halo nodes are also in boundary
 //          node storage and should be removed from this first
 for (int d=0;d<n_proc;d++)
  {
   // How many external haloes with this process?
   unsigned n_ext_halo_nod=nexternal_halo_node(d);
   for (unsigned j=0;j<n_ext_halo_nod;j++)
    {
     Node* ext_halo_nod_pt=external_halo_node_pt(d,j);
     unsigned n_bnd=nboundary();
     for (unsigned i_bnd=0;i_bnd<n_bnd;i_bnd++)
      {
       // Call this for all boundaries; it will do nothing
       // if the node is not on the current boundary
       remove_boundary_node(i_bnd,ext_halo_nod_pt);
      }
    }
  }

 // A loop to delete external halo nodes
 for (int d=0;d<n_proc;d++)
  {
   unsigned n_ext_halo_nod=nexternal_halo_node(d);
   for (unsigned j=0;j<n_ext_halo_nod;j++)
    {
     // Only delete if it's not a node stored in the current mesh
     bool is_a_mesh_node=false;
     unsigned n_node=nnode();
     for (unsigned jj=0;jj<n_node;jj++)
      {
       if (Node_pt[jj]==External_halo_node_pt[d][j])
        {
         is_a_mesh_node=true;
        }
      }

     // There will also be duplications between multiple processors,
     // so make sure that we don't try to delete these twice
     if (!is_a_mesh_node)
      {
       // Loop over all other higher-numbered processors and check
       // for duplicated external halo nodes
       // (The highest numbered processor should delete all its ext halos)
       for (int dd=d+1;dd<n_proc;dd++)
        {
         unsigned n_ext_halo=nexternal_halo_node(dd);
         for (unsigned jjj=0;jjj<n_ext_halo;jjj++)
          {
           if (External_halo_node_pt[dd][jjj]==External_halo_node_pt[d][j])
            {
             is_a_mesh_node=true;
            }
          }
        }
      }

     // Only now if no duplicates exist can the node be safely deleted
     if (!is_a_mesh_node)
      {
       delete External_halo_node_pt[d][j];
//       External_halo_node_pt[d][j]=0;
      }
    }
  }

 // Another loop to delete external halo elements (which are distinct)
 for (int d=0;d<n_proc;d++)
  {
   unsigned n_ext_halo_el=nexternal_halo_element(d);
   for (unsigned e=0;e<n_ext_halo_el;e++)
    {
     delete External_halo_element_pt[d][e];
//     External_halo_element_pt[d][e]=0;
    }
  }

 // Now we are okay to clear the external halo node storage
 External_halo_node_pt.clear();
 External_halo_element_pt.clear();

 // External haloed nodes and elements are actual members
 // of the external mesh and should not be deleted
 External_haloed_node_pt.clear();
 External_haloed_element_pt.clear();
#endif
}


//========================================================================
/// \short Add external element to this Mesh.
//========================================================================
bool Mesh::add_external_element_pt(FiniteElement*& el_pt)
{
 // Only add the element to the external element storage if
 // it's not already there
 bool already_external_element=false;
 unsigned n_ext_el=nexternal_element();
 for (unsigned e_ext=0;e_ext<n_ext_el;e_ext++)
  {
   if (el_pt==External_element_pt[e_ext])
    {
     already_external_element=true;
     break;
    }
  }
 if (!already_external_element)
  {
   External_element_pt.push_back(el_pt);
   return true;
  }
 else
  {
   return false;
  }
}

//========================================================================
/// \short Add external node to this Mesh.
//========================================================================
bool Mesh::add_external_node_pt(Node*& nod_pt)
{
 // Only add the node if it's not already there
 bool node_is_external=false;
 unsigned n_ext_nod=nexternal_node();
 for (unsigned j_ext=0;j_ext<n_ext_nod;j_ext++)
  {
   if (nod_pt==External_node_pt[j_ext])
    {
     node_is_external=true;
     break;
    }
  }
 if (!node_is_external)
  {
   External_node_pt.push_back(nod_pt);
   return true;
  }
 else
  {
   return false;
  }
}

#ifdef OOMPH_HAS_MPI

// NOTE: the add_external_haloed_node_pt and add_external_haloed_element_pt
//       functions need to check whether the Node/FiniteElement argument
//       has been added to the storage already; this is not the case
//       for the add_external_halo_node_pt and add_external_halo_element_pt
//       functions as these are newly-created elements that are created and
//       added to the storage based on the knowledge of when their haloed 
//       counterparts were created and whether they were newly added

//========================================================================
/// \short Add external haloed element whose non-halo counterpart is held 
/// on processor p to the storage scheme for external haloed elements.
/// If the element is already in the storage scheme then return its index
//========================================================================
unsigned Mesh::add_external_haloed_element_pt(const unsigned& p, 
                                              FiniteElement*& el_pt)
{
 // Loop over current storage
 unsigned n_extern_haloed=nexternal_haloed_element(p);

 // Is this already an external haloed element?
 bool already_external_haloed_element=false;
 unsigned external_haloed_el_index=0;
 for (unsigned eh=0;eh<n_extern_haloed;eh++)
  {
   if (el_pt==External_haloed_element_pt[p][eh])
    {
     // It's already there, so...
     already_external_haloed_element=true;
     // ...set the index of this element
     external_haloed_el_index=eh;
     break;
    }
  }

 // Has it been found?
 if (!already_external_haloed_element)
  {
   // Not found, so add it:
   External_haloed_element_pt[p].push_back(el_pt);
   // Return the index where it's just been added
   return n_extern_haloed;
  }
 else
  {
   // Return the index where it was found
   return external_haloed_el_index;
  }
}

//========================================================================
/// \short Add external haloed node whose halo (external) counterpart
/// is held on processor p to the storage scheme for external haloed nodes.
/// If the node is already in the storage scheme then return its index
//========================================================================
unsigned Mesh::add_external_haloed_node_pt(const unsigned& p, Node*& nod_pt)
{
 // Loop over current storage
 unsigned n_ext_haloed_nod=nexternal_haloed_node(p);

 // Is this already an external haloed node?
 bool is_an_external_haloed_node=false;
 unsigned external_haloed_node_index=0;
 for (unsigned k=0;k<n_ext_haloed_nod;k++)
  {
   if (nod_pt==External_haloed_node_pt[p][k])
    {
     is_an_external_haloed_node=true;
     external_haloed_node_index=k;
     break;
    }
  }

 // Has it been found?
 if (!is_an_external_haloed_node)
  {
   // Not found, so add it
   External_haloed_node_pt[p].push_back(nod_pt);
   // Return the index where it's just been added
   return n_ext_haloed_nod;
  }
 else
  {
   // Return the index where it was found
   return external_haloed_node_index;
  }
}

#endif


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Functions for solid meshes
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//========================================================================
/// Make the current configuration the undeformed one by
/// setting the nodal Lagrangian coordinates to their current
/// Eulerian ones
//========================================================================
void SolidMesh::set_lagrangian_nodal_coordinates()
{ 

 //Find out how many nodes there are
 unsigned long n_node = nnode();
 
 //Loop over all the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Cast node to solid node (can safely be done because
   // SolidMeshes consist of SolidNodes
   SolidNode* node_pt = static_cast<SolidNode*>(Node_pt[n]);
   
   // Number of Lagrangian coordinates
   unsigned n_lagrangian = node_pt->nlagrangian();

   // Number of generalised Lagrangian coordinates
   unsigned n_lagrangian_type = node_pt->nlagrangian_type();

   //The assumption here is that there must be fewer lagrangian coordinates
   //than eulerian (which must be true?)

   // Set (generalised) Lagrangian coords = (generalised) Eulerian coords
   for(unsigned k=0;k<n_lagrangian_type;k++)
    {
     // Loop over lagrangian coordinates and set their values
     for(unsigned j=0;j<n_lagrangian;j++)
      {
       node_pt->xi_gen(k,j)=node_pt->x_gen(k,j);
      }
    }
  }
}


//=======================================================================
/// Static problem that can be used to assign initial conditions
/// on a given mesh.
//=======================================================================
SolidICProblem SolidMesh::Solid_IC_problem;


}

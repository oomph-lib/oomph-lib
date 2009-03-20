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
#include "partitioning.h"
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
 
 //The above call is fine for all nodes, since the boundary gets added
 //to Boundaries_pt, which is a set and therefore has no duplication.
 //However, this is not true for the Boundary_node_pt vector, and this
 //is causing a fault in the parallel version of the static fish problem,
 //where Mesh::prune_dead_nodes attempts to remove the same point twice,
 //but obviously can't the second time and therefore causes a seg fault.
 //The new code checks whether the node has already been added to the
 //Boundary_node_pt vector for this particular boundary.

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
   /// Try to cast to elastic node \todo there's got to be a better way
   /// but making Problem::mesh_pt() virtual doesn't do the right thing...
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
     el_pt->output(outfile);
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
     el_pt->output(outfile,n_plot);
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
     el_pt->output(file_pt);
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
     el_pt->output(file_pt,n_plot);
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
     el_pt->output_fct(outfile,n_plot,exact_soln_pt);
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
     el_pt->output_fct(outfile,n_plot,time,exact_soln_pt);
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
/// Classify the halo and haloed nodes in the mesh
/// TO DO add internal and external data and then possibly rename
//========================================================================
void Mesh::classify_halo_and_haloed_nodes(DocInfo& doc_info,
                                          const bool& report_stats)
{
 //Wipe existing storage schemes for halo(ed) and shared nodes
 Halo_node_pt.clear();
 Haloed_node_pt.clear();
 Shared_node_pt.clear();
 
 // Determine which processors the nodes are associated with
 // and hence who's in charge; only needs to be done for
 // nodes that aren't classified yet as their affiliation never
 // changes
 
 // Loop over processors that share haloes with the current one
 // and associate nnodes with that processor
// for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
//      it!=Processors_that_share_halos.end();it++)
 for (int domain=0;domain<MPI_Helpers::Nproc;domain++) // try looping over every processor
  {   
   // Which adjacent domain/processor are we dealing with?
//   unsigned domain=*it;

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
       if (nod_pt->processor_in_charge()==
           Data::Not_associated_with_any_processor)
        {
         nod_pt->processors_associated_with_data().insert(domain);
        }
       // do the same if the node is solid
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
       if (solid_nod_pt!=0)
        {
         if (solid_nod_pt->variable_position_pt()->processor_in_charge()==
           Data::Not_associated_with_any_processor)
          {
           solid_nod_pt->variable_position_pt()->
            processors_associated_with_data().insert(domain);
          }
        }
      }
    }
  }
 
 
 // Loop over all non-halo elements and associate their nodes
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
       if (nod_pt->processor_in_charge()==
           Data::Not_associated_with_any_processor)
        {
         // Associate node with current processor
         nod_pt->processors_associated_with_data().
          insert(MPI_Helpers::My_rank);
        }
       // do the same if we have a SolidNode
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
       if (solid_nod_pt!=0)
        {
         if (solid_nod_pt->variable_position_pt()->processor_in_charge()==
           Data::Not_associated_with_any_processor)
          {
           // Associate node with current processor
           solid_nod_pt->variable_position_pt()->
            processors_associated_with_data().insert(MPI_Helpers::My_rank);
          }
        }
      }
    }
  }



 // Loop over all nodes on the present processor and put the highest-numbered
 // one in charge of the node
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=this->node_pt(j);
   if (nod_pt->processor_in_charge()==
       Data::Not_associated_with_any_processor)
    {
     // Now put the highest-numbered one in charge
     unsigned proc_max=0;
     std::set<unsigned> procs_set=nod_pt->processors_associated_with_data();
     for (std::set<unsigned>::iterator it=procs_set.begin();
          it!=procs_set.end();it++)
      {
       if (*it>proc_max) proc_max=*it;
      }
     nod_pt->processor_in_charge()=proc_max;
     
    }
   // Do the same if we have a SolidNode
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);
   if (solid_nod_pt!=0)
    {
     if (solid_nod_pt->variable_position_pt()->processor_in_charge()==
         Data::Not_associated_with_any_processor)
      {
       // Now put the highest-numbered one in charge
       unsigned proc_max_solid=0;
     std::set<unsigned> procs_set_solid=solid_nod_pt->
           variable_position_pt()->processors_associated_with_data();
     for (std::set<unsigned>::iterator it=procs_set_solid.begin();
          it!=procs_set_solid.end();it++)
      {
       if (*it>proc_max_solid) proc_max_solid=*it;
      }
     solid_nod_pt->variable_position_pt()->
          processor_in_charge()=proc_max_solid;
      }
    }
  }

 // Determine halo nodes. They are located on the halo
 // elements and the processor in charge differs from the
 // current processor

 // Only count nodes once (map is initialised to 0 = false)
 std::map<Node*,bool> done;

 // Loop over processors that share haloes with the current one
// for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
//      it!=Processors_that_share_halos.end();it++)
 for (int domain=0;domain<MPI_Helpers::Nproc;domain++) 
// try looping over every processor?
  {
   // Which adjacent domain/processor are we dealing with?
//   unsigned domain=*it;

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
         int proc_in_charge=nod_pt->processor_in_charge();

         // To keep the order of the nodes consistent with that
         // in the haloed node lookup scheme, only 
         // allow it to be added when the current domain is in charge
         if (proc_in_charge==int(domain))
          {
           if (proc_in_charge!=MPI_Helpers::My_rank)
            {
             // Add it as being halo node whose non-halo counterpart
             // is located on processor domain
             this->add_halo_node_pt(proc_in_charge,nod_pt);

             // We're done with this node
             done[nod_pt]=true;
            }
          }

        }
       // It is not necessary to do this for SolidNodes
       // as the halo/haloed structure is taken care of by using
       // dynamic casting in Problem::synchronise_*

      }
     // Now make sure internal data on halo elements is also halo
     unsigned nintern_data = el_pt->ninternal_data();
     for (unsigned iintern=0;iintern<nintern_data;iintern++)
      {
       // no need to go in and pin eqn number (done in 
       // Data::assign_eqn_numbers); just need to set as halo
       // by telling it that the "other processor" is in charge
       // (see Data::is_halo() ...)
       //
       // domain != MPI_Helpers::My_rank (i.e. the current process)
       // since a process cannot share halo elements with itself
       //
       el_pt->internal_data_pt(iintern)->processor_in_charge()=domain;

      }
    }

  }



 // Determine haloed nodes. They are located on the haloed
 // elements and the processor in charge is the current processor
 
 // Loop over processors that share haloes with the current one
// for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
//      it!=Processors_that_share_halos.end();it++)
 for (int domain=0;domain<MPI_Helpers::Nproc;domain++) // try looping over every processor
  {
   // Which adjacent domain/processor are we dealing with?
//   unsigned domain=*it;

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
         int proc_in_charge=nod_pt->processor_in_charge();

         if (proc_in_charge==MPI_Helpers::My_rank)
          {
           // Add it as being haloed from specified domain
           this->add_haloed_node_pt(domain,nod_pt);
           
          }

         // We're done with this node
         node_done[nod_pt]=true;
          
        }

      }
    }

  }

 // Determine shared nodes - located on the halo(ed) elements
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {
   // map of bools for whether the node has been shared,
   // initialised to 0 (false) for each d
   std::map<Node*,bool> node_shared;

   if (d<MPI_Helpers::My_rank)
    {
     // Get the nodes from the halo elements first
     Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(d));
     unsigned nhalo_elem=halo_elem_pt.size();

     for (unsigned e=0;e<nhalo_elem;e++)
      {
       // Get element
       FiniteElement* el_pt=halo_elem_pt[e];
       unsigned nnod=el_pt->nnode();

       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);

         // Add it as a shared node from current domain (?)
         if (!node_shared[nod_pt])
          {
           this->add_shared_node_pt(d,nod_pt);
           node_shared[nod_pt]=true;
          }

        } // loop over nodes

      } // loop over elements

     // Now get the nodes from any haloed elements left over
     Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(d));
     unsigned nhaloed_elem=haloed_elem_pt.size();

     for (unsigned e=0;e<nhaloed_elem;e++)
      {
       // Get element
       FiniteElement* el_pt=haloed_elem_pt[e];
       unsigned nnod=el_pt->nnode();

       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);

         // Add it as a shared node from current domain (?)
         if (!node_shared[nod_pt])
          {
           this->add_shared_node_pt(d,nod_pt);
           node_shared[nod_pt]=true;
          }

        } // loop over nodes

      } // loop over elements

    }
   if (d>MPI_Helpers::My_rank)
    {
     // Get the nodes from the haloed elements first
     Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(d));
     unsigned nhaloed_elem=haloed_elem_pt.size();

     for (unsigned e=0;e<nhaloed_elem;e++)
      {
       // Get element
       FiniteElement* el_pt=haloed_elem_pt[e];
       unsigned nnod=el_pt->nnode();

       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);

         // Add it as a shared node from current domain (?)
         if (!node_shared[nod_pt])
          {
           this->add_shared_node_pt(d,nod_pt);
           node_shared[nod_pt]=true;
          }

        } // loop over nodes

      } // loop over elements

     // Now get the nodes from any halo elements left over
     Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(d));
     unsigned nhalo_elem=halo_elem_pt.size();

     for (unsigned e=0;e<nhalo_elem;e++)
      {
       // Get element
       FiniteElement* el_pt=halo_elem_pt[e];
       unsigned nnod=el_pt->nnode();

       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt=el_pt->node_pt(j);

         // Add it as a shared node from current domain (?)
         if (!node_shared[nod_pt])
          {
           this->add_shared_node_pt(d,nod_pt);
           node_shared[nod_pt]=true;
          }

        } // loop over nodes

      } // loop over elements

    } // end if (d ...)
    
  } // loop over processes

 // Doc stats
// if (report_stats)
//  {
   oomph_info << "Processor " << MPI_Helpers::My_rank 
              << " holds " << this->nnode() 
              << " nodes of which " << this->nhalo_node()
              << " are halo nodes \n while " << this->nhaloed_node()
              << " are haloed nodes, and " << this->nshared_node()
              << " are shared nodes." << std::endl;
   
//  }

}





//========================================================================
/// Get halo node stats for this distributed mesh:
/// Average/max/min number of halo nodes over all processors.
/// \b Careful: Involves MPI Broadcasts and must therefore
/// be called on all processors!
//========================================================================
void Mesh::get_halo_node_stats(double& av_number, 
                               unsigned& max_number,
                               unsigned& min_number)
{
 // Create vector to hold number of halo nodes
 Vector<int> nhalo_nodes(MPI_Helpers::Nproc);
 
 // Stick own number of halo nodes into appropriate entry 
 nhalo_nodes[MPI_Helpers::My_rank]=nhalo_node();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of dofs on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhalo_nodes[MPI_Helpers::My_rank],1,MPI_INT,
            &nhalo_nodes[0],1, MPI_INT,
            0,MPI_COMM_WORLD);

 // Initialise stats
 av_number=0.0;
 int max=-1;
 int min=1000000000;

 if (MPI_Helpers::My_rank==0)
  {
   for (int i=0;i<MPI_Helpers::Nproc;i++)
    {
     av_number+=double(nhalo_nodes[i]);
     if (int(nhalo_nodes[i])>max) max=nhalo_nodes[i];
     if (int(nhalo_nodes[i])<min) min=nhalo_nodes[i];

    }  
   av_number/=double(MPI_Helpers::Nproc); 
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);
 MPI_Bcast(&min,1,MPI_INT,0,MPI_COMM_WORLD);
 MPI_Bcast(&av_number,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 
 max_number=max;
 min_number=min;  
}


//========================================================================
/// Get haloed node stats for this distributed mesh:
/// Average/max/min number of haloed nodes over all processors.
/// \b Careful: Involves MPI Broadcasts and must therefore
/// be called on all processors!
//========================================================================
 void  Mesh::get_haloed_node_stats(double& av_number, 
                                   unsigned& max_number,
                                   unsigned& min_number)
{
 // Create vector to hold number of haloed nodes
 Vector<int> nhaloed_nodes(MPI_Helpers::Nproc);
 
 // Stick own number of haloed nodes into appropriate entry
 nhaloed_nodes[MPI_Helpers::My_rank]=nhaloed_node();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of dofs on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhaloed_nodes[MPI_Helpers::My_rank],1,MPI_INT,
            &nhaloed_nodes[0],1, MPI_INT,
            0,MPI_COMM_WORLD);

 // Initialise stats
 av_number=0.0;
 int max=-1;
 int min=1000000000;

 if (MPI_Helpers::My_rank==0)
  {
   for (int i=0;i<MPI_Helpers::Nproc;i++)
    {
     av_number+=double(nhaloed_nodes[i]);
     if (int(nhaloed_nodes[i])>max) max=nhaloed_nodes[i];
     if (int(nhaloed_nodes[i])<min) min=nhaloed_nodes[i];

    }  
   av_number/=double(MPI_Helpers::Nproc); 
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_INT,0,MPI_COMM_WORLD);
 MPI_Bcast(&min,1,MPI_INT,0,MPI_COMM_WORLD);
 MPI_Bcast(&av_number,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 
 max_number=max;
 min_number=min;  
}

//========================================================================
/// Distribute the mesh
//========================================================================
void Mesh::distribute(DocInfo& doc_info, 
                      const bool& report_stats)
{ 

 char filename[100];
 std::ofstream some_file;


 // Doc the original mesh on proc 0
 //--------------------------------
 if (doc_info.doc_flag())
  {
   if (MPI_Helpers::My_rank==0)
    {
     sprintf(filename,"%s/complete_mesh%i.dat",doc_info.directory().c_str(),
             doc_info.number());
     this->output(filename,5);
    }
  }

//  oomph_info << "Partioning mesh....";
//  clock_t t_start=clock();

 // Partition the mesh
 //-------------------
 // Objective (0: minimise edge cut; 1: minimise total comm volume)
 unsigned objective=0;

 // Record the original total number of elements in the mesh
 // to be able to assess the efficiency of the distribution
 unsigned orig_nelem=this->nelement();

 // Vector listing the affiliation of each element
 unsigned nelem=orig_nelem;
 Vector<unsigned> element_domain(nelem);
  
 // Do the partitioning
 METIS::partition_mesh(this,MPI_Helpers::Nproc,objective,element_domain);

 // On very coarse meshes with larger numbers of processors, METIS occasionally
 // returns an element_domain Vector for which a particular processor has no
 // elements affiliated to it; the following takes care of this...

 // Convert element_domain to integer storage
 Vector<int> int_element_domain(nelem);
 for (unsigned e=0;e<nelem;e++)
  {
   int_element_domain[e]=element_domain[e];
  }

 // Global storage for number of elements on each process
 int my_number_of_elements=0;
 Vector<int> number_of_elements(MPI_Helpers::Nproc,0);

 for (unsigned e=0;e<nelem;e++)
  {
   if (int_element_domain[e]==MPI_Helpers::My_rank)
    {
     my_number_of_elements++;
    }
  }

 // Communicate the correct value for each single process into
 // the global storage vector
 MPI_Allgather(&my_number_of_elements,1,MPI_INT,
               &number_of_elements[0],1,MPI_INT,MPI_COMM_WORLD);

 // If a process has no elements then switch an element with the
 // process with the largest number of elements, assuming
 // that it still has enough elements left to share
 int max_number_of_elements=0;
 int process_with_max_elements=0;
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {
   if (number_of_elements[d]==0)
    {
     // Find the process with maximum number of elements, if necessary
     if (max_number_of_elements<=1)
      {
       for (int dd=0;dd<MPI_Helpers::Nproc;dd++)
        {
         if (number_of_elements[dd]>max_number_of_elements)
          {
           max_number_of_elements=number_of_elements[dd];
           process_with_max_elements=dd;
          }
        }
      }

     // Check that this number of elements is okay for sharing...
     if (max_number_of_elements<=1)
      {
       // Throw error; we shouldn't arrive here, but...
       std::ostringstream error_stream;
       error_stream << "No process has more than 1 element, and\n"
                    << "at least one process has no elements!\n"
                    << "Suggest rerunning with more refinement.\n"
                    << std::endl;
       throw OomphLibError(error_stream.str(),
                           "Mesh::distribute()",
                           OOMPH_EXCEPTION_LOCATION);

      }

     // Loop over the element domain vector and switch
     // one value for process "process_with_max_elements" with d
     for (unsigned e=0;e<nelem;e++)
      {
       if (int_element_domain[e]==process_with_max_elements)
        {
         int_element_domain[e]=d;
         // Change the numbers associated with these processes
         number_of_elements[d]++;
         number_of_elements[process_with_max_elements]--;
         // Reduce the number of elements available on "max" process
         max_number_of_elements--;
         // Inform the user that a switch has taken place
         oomph_info << "INFO: Switched element domain at position " << e 
                    << std::endl 
                    << "from process " << process_with_max_elements 
                    << " to process " << d 
                    << std::endl
                    << "which was given no elements by METIS partition" 
                    << std::endl;           
         // Only need to do this once for this element loop, otherwise
         // this will take all the elements from "max" process and put them
         // in process d, thus leaving essentially the same bug!
         break;
        }
      }
    }

  }

 // Reassign new values to the element_domain vector
 for (unsigned e=0;e<nelem;e++)
  {
   element_domain[e]=int_element_domain[e];
  }

 // Doc the partitioning (only on processor 0) 
 //-------------------------------------------
 if (doc_info.doc_flag())
  {
   if (MPI_Helpers::My_rank==0)
    {
     // Open files for doc of element partitioning
     Vector<std::ofstream*> domain_file(MPI_Helpers::Nproc);
     for (int d=0;d<MPI_Helpers::Nproc;d++)
      {
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
     
     for (int d=0;d<MPI_Helpers::Nproc;d++)
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

     // Recall that processor in charge is initialised to -1
     if (int(el_domain)>nod_pt->processor_in_charge())
      {
       nod_pt->processor_in_charge()=el_domain;
      }
     nod_pt->processors_associated_with_data().insert(el_domain);
    }
  }

 // Doc the partitioning (only on processor 0) 
 //-------------------------------------------
 if (doc_info.doc_flag())
  {
   if (MPI_Helpers::My_rank==0)
    {
     // Open files for doc of node partitioning
     Vector<std::ofstream*> node_file(MPI_Helpers::Nproc);
     for (int d=0;d<MPI_Helpers::Nproc;d++)
      {
       sprintf(filename,"%s/node%i-%i.dat",doc_info.directory().c_str(),
               d,doc_info.number());
       node_file[d]=new std::ofstream(filename);
      }
     
     // Doc
     unsigned nnod=this->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=this->node_pt(j);
       *node_file[nod_pt->processor_in_charge()]
        << nod_pt->x(0) << " " 
        << nod_pt->x(1) << std::endl;
      }
     for (int d=0;d<MPI_Helpers::Nproc;d++)
      {
       node_file[d]->close();
      }
    }
  }

 // Declare all nodes as obsolete. We'll
 // change this setting for all nodes that must be retained
 // further down
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   this->node_pt(j)->set_obsolete();
  }

 // Backup old mesh data and flush mesh
 //-------------------------------------

 // Backup pointers to elements in this mesh
 nelem=this->nelement();
 Vector<FiniteElement*> backed_up_el_pt(nelem);
 for (unsigned e=0;e<nelem;e++)
  {
   backed_up_el_pt[e]=this->finite_element_pt(e);
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
 tmp_element_retained.resize(MPI_Helpers::Nproc);
 nelem=backed_up_el_pt.size();
 for (int i=0;i<MPI_Helpers::Nproc;i++)
  {
   tmp_element_retained[i].resize(nelem,false);
  }
 
 // Temporary storage for root halo elements on the various
 // processors. Needed to figure out haloed lookup schemes:
 // When setting these up on any given processor we have to know
 // which elements will (have) become halo elements on other processors.
 Vector<Vector<Vector<FiniteElement*> > > tmp_root_halo_element_pt;
 tmp_root_halo_element_pt.resize(MPI_Helpers::Nproc);
 for (int i=0;i<MPI_Helpers::Nproc;i++)
  {
   tmp_root_halo_element_pt[i].resize(MPI_Helpers::Nproc);
  }

 // Determine which elements are going to end up on which processor
 //----------------------------------------------------------------

 // This procedure needs to be repeated to catch elements which may
 // be missed the first time round but which contain nodes from this process

 unsigned elements_retained=true;
 int myi=1;
 while (elements_retained) 
  {
 Vector<unsigned> number_of_retained_elements(MPI_Helpers::Nproc,0);
 int number_of_retained_halo_elements=0; // not dependent on dummy_my_rank

 for (int dummy_my_rank=0;dummy_my_rank<MPI_Helpers::Nproc;
      dummy_my_rank++)
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
               it=nod_pt->processors_associated_with_data().begin();
              it!=nod_pt->processors_associated_with_data().end();
              it++)
          {
           if (*it==unsigned(dummy_my_rank))
            {
             keep_it=true;
             break;
            }
          }
         
         if (keep_it)
          {
           // Add as root halo element whose non-halo counterpart is
           // located on processor el_domain
           tmp_root_halo_element_pt[dummy_my_rank][el_domain].push_back(el_pt);
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
   for (int d=0;d<MPI_Helpers::Nproc;d++)
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
         int proc_in_charge=nod_pt->processor_in_charge();
         if (proc_in_charge>d)
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
             (int(element_domain[e])==nod_pt->processor_in_charge()))
          {
           keep_it=true; 
           break;
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

   // Check number of retained halo elements on this process
   oomph_info << "percentage of extra halo elements retained: " 
              << 100.0*double(number_of_retained_halo_elements)/
                  double(number_of_retained_elements[dummy_my_rank])
              << " on process " << dummy_my_rank 
              << " in loop number " << myi << std::endl;

  } // end of loop over all "processors"; we've now established the
    // elements and the root halo elements for all processors

   int total_number_of_retained_halo_elements;
   
   // Sum values over all processes 
   // - must be zero retained in order to continue
   MPI_Allreduce(&number_of_retained_halo_elements, 
                &total_number_of_retained_halo_elements, 1, MPI_INT, 
                 MPI_SUM, MPI_COMM_WORLD);

   oomph_info << "total number of extra halo elements retained: " 
              << total_number_of_retained_halo_elements
              << " in loop: " << myi << std::endl;

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
   if (tmp_element_retained[MPI_Helpers::My_rank][e])
    {
     this->add_element_pt(el_pt);
    }
   else
    {
     delete el_pt;
    }
  }
 
 // Copy the root halo elements associated with the actual
 // current processor into its own permanent storage; the order
 // here is somewhat random but we compensate for that by
 // ensuring that the corresponding haloed elements are 
 // added in the same order below
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {  
   nelem=tmp_root_halo_element_pt[MPI_Helpers::My_rank][d].size();
   for (unsigned e=0;e<nelem;e++)
    {
     this->add_root_halo_element_pt(d,
      tmp_root_halo_element_pt[MPI_Helpers::My_rank][d][e]);
    }
  }
  
//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Determine root haloed elements....";

 // Determine root haloed elements
 //-------------------------------

 // Loop over all other processors
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {
   if (d!=MPI_Helpers::My_rank)
    {
     // Loop over root halo elements that are held on that processor
     unsigned nelem_other=
      tmp_root_halo_element_pt[d][MPI_Helpers::My_rank].size();
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
           if (el_pt==
               tmp_root_halo_element_pt[d][MPI_Helpers::My_rank][e2])
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
   oomph_info << "Processor " << MPI_Helpers::My_rank 
              << " holds " << this->nelement() 
              << " elements of which " << this->nroot_halo_element()
              << " are root halo elements \n while " 
              << this->nroot_haloed_element()
              << " are root haloed elements" << std::endl;
  }

// Return an error if any process has no elements!
 if (this->nelement()==0)
  {
   throw OomphLibError(
    "Process has no elements; re-run with more uniform refinement!\n",
    "Mesh::distribute()",
    OOMPH_EXCEPTION_LOCATION);
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
 RefineableMeshBase* ref_mesh_pt=dynamic_cast<RefineableMeshBase*>(this);
 if (ref_mesh_pt!=0)
  {
   ref_mesh_pt->setup_tree_forest();
  }

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

//  oomph_info << "Classify nodes....";

 // Classify nodes 
 classify_halo_and_haloed_nodes(doc_info,report_stats);

//  elapsed_time=double(clock()-t_start)/CLOCKS_PER_SEC; t_start=clock();
//  oomph_info << "....done " << elapsed_time << std::endl;

 // Mesh has now been distributed 
 // (required for Z2ErrorEstimator::get_element_errors)
 Mesh_has_been_distributed=true;

 // Doc?
 //-----
 if (doc_info.doc_flag())
  {
   doc_mesh_distribution(doc_info);
  }


}


//========================================================================
/// (Irreversibly) redistribute elements and nodes, usually
/// after another round of refinement, to get rid of
/// excessively wide halo layers. Note that the current
/// mesh will be now regarded as the base mesh and no unrefinement
/// relative to it will be possible once this function 
/// has been called.
//========================================================================
void Mesh::redistribute(DocInfo& doc_info, 
                        const bool& report_stats)
{
 // Doc stats
 if (report_stats)
  {
   oomph_info << "Before redistribution: Processor " << MPI_Helpers::My_rank 
              << " holds " << this->nelement() 
              << " elements of which " << this->nroot_halo_element()
              << " are root halo elements \n while " 
              << this->nroot_haloed_element()
              << " are root haloed elements" << std::endl;
  }
   

 // Declare all nodes as obsolete. We'll
 // change this setting for all nodes that must to be retained
 // further down
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   this->node_pt(j)->set_obsolete();
  }
 

 // Backup old mesh data
 //---------------------
   
 // Backup pointers to elements in this mesh
 unsigned nelem=this->nelement();
 Vector<FiniteElement*> backed_up_el_pt(nelem);
 std::map<FiniteElement*,bool> keep_element;
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=this->finite_element_pt(e);
   backed_up_el_pt[e]=el_pt;

   // Definitely retain the element if it's not a halo
   if (!el_pt->is_halo())
    {
     keep_element[el_pt]=true;

     // Loop over the element's nodes and retain them too
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       el_pt->node_pt(j)->set_non_obsolete();
      }
    }
  }
   

 // Now loop over all halo elements and check if they
 // have at least one node that's to be retained, i.e.
 // check if the halo element is directly connected
 // with the bulk of the mesh containing the non-halo elements.
         
 // Temp map of vectors holding the pointers to the root halo elements
 std::map<unsigned, Vector<FiniteElement*> > tmp_root_halo_element_pt;

 // Temp map of vectors holding the pointers to the root haloed elements
 std::map<unsigned, Vector<FiniteElement*> > tmp_root_haloed_element_pt;
 
 // Map to store if a halo element survives
 std::map<FiniteElement*,bool> halo_element_is_retained;

//    sprintf(filename,"%s/retained_halo_elements%i_on_proc%i.dat",
//            doc_info.directory().c_str(),
//            doc_info.number(),MPI_Helpers::My_rank);
//    some_file.open(filename);

  
 for (int domain=0;domain<MPI_Helpers::Nproc;domain++)
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
       if (!nod_pt->is_obsolete())
        {
         // Keep element and add it to preliminary storage for
         // halo elements associated with current neighbouring domain
         keep_element[el_pt]=true;
//             el_pt->output(some_file,5);
         tmp_root_halo_element_pt[domain].push_back(el_pt);
         halo_element_is_retained[el_pt]=true;
         break;
        }
      }       
    }
  }
  

//   some_file.close();


 // Make sure everybody finishes this part
 MPI_Barrier(MPI_COMM_WORLD);


//    sprintf(filename,"%s/retained_haloed_elements%i_on_proc%i.dat",
//            doc_info.directory().c_str(),
//            doc_info.number(),MPI_Helpers::My_rank);
//    some_file.open(filename);


 // Now all processors have decided (independently) which of their
 // (to-be root) halo elements they wish to retain. Now we need to figure out
 // which of their elements are haloed and add them in the appropriate
 // order into the haloed element scheme. For this we exploit that
 // the halo and haloed elements are accessed in the same order on
 // all processors!
   
 // Identify haloed elements on domain d
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {
   // Loop over domains that halo this domain
   for (int dd=0;dd<MPI_Helpers::Nproc;dd++)
    {       
     // Dont't talk to yourself
     if (d!=dd)
      {

       // If we're identifying my haloed elements:
       if (d==MPI_Helpers::My_rank)
        {
         // Get vector all elements that are currently haloed by domain dd
         Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(dd));         
         // Create a vector of ints to indicate if the halo element
         // on processor dd processor was kept
         unsigned nelem=haloed_elem_pt.size();
         Vector<int> halo_kept(nelem);
           
         // Receive this vector from processor dd 
         MPI_Status status;
         MPI_Recv(&halo_kept[0],nelem,MPI_INT,dd,0,MPI_COMM_WORLD,&status);
           
         // Classify haloed element accordingly
         for (unsigned e=0;e<nelem;e++)
          {
           FiniteElement* el_pt=haloed_elem_pt[e];
           if (halo_kept[e]==1)
            {
             // I am being haloed by processor dd
             tmp_root_haloed_element_pt[dd].push_back(el_pt);
//               el_pt->output(some_file,5);
            }
          }
        }
       else
        {
         // If we're dealing with my halo elements:
         if (dd==MPI_Helpers::My_rank)
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
             FiniteElement* el_pt=halo_elem_pt[e];
             if (halo_element_is_retained[el_pt])
              {
               halo_kept[e]=1;
              }
            }
             
           // Now send this vector to processor d to tell it which of
           // the haloed elements (which are listed in the same order)
           // are to be retained as haloed elements.
           MPI_Send(&halo_kept[0],nelem,MPI_INT,d,0,MPI_COMM_WORLD);
          }
        }
      }
    }
  }
 
//   some_file.close();
  
 // Backup pointers to nodes in this mesh
 nnod=this->nnode();
 Vector<Node*> backed_up_nod_pt(nnod);
 for (unsigned j=0;j<nnod;j++)
  {
   backed_up_nod_pt[j]=this->node_pt(j);
  }
   
 // Flush the mesh storage
 this->flush_element_and_node_storage();
   
 // Loop over all backed up elements
 nelem=backed_up_el_pt.size();
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=backed_up_el_pt[e];
   if (keep_element[el_pt])
    {
     this->add_element_pt(el_pt);
    }
   else
    {
     delete el_pt;
    }
  }

 // Wipe the storage scheme for halo(ed) elements and then re-assign
 Root_haloed_element_pt.clear();
 Root_halo_element_pt.clear();     
 for (int domain=0;domain<MPI_Helpers::Nproc;domain++)
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
   oomph_info << "AFTER:  Processor " << MPI_Helpers::My_rank 
              << " holds " << this->nelement() 
              << " elements of which " << this->nroot_halo_element()
              << " are root halo elements \n while " 
              << this->nroot_haloed_element()
              << " are root haloed elements" << std::endl;
  }
   

      
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
   
 // Prune and rebuild mesh
 //-----------------------
   
 // Now remove the pruned nodes from the boundary lookup scheme
 this->prune_dead_nodes();
    
 // And finally re-setup the boundary lookup scheme for elements
 this->setup_boundary_element_info();
      
 // Re-setup tree forest if needed
 RefineableMeshBase* ref_mesh_pt=dynamic_cast<RefineableMeshBase*>(this);
 if (ref_mesh_pt!=0)
  {
   ref_mesh_pt->setup_tree_forest();
  }
   
 // Classify nodes 
 classify_halo_and_haloed_nodes(doc_info,report_stats);
   
 // Doc?
 //-----
 if (doc_info.doc_flag())
  {
   doc_mesh_distribution(doc_info);
  }
   
   
}






//========================================================================
///  Get efficiency of mesh distribution: In an ideal distribution
/// without halo overhead, each processor would only hold its own
/// elements. Efficieny per processor =  (number of non-halo elements)/
/// (total number of elements). 
//========================================================================
void Mesh:: get_efficiency_of_mesh_distribution(double& av_efficiency, 
                                                double& max_efficiency,
                                                double& min_efficiency)
{
 // Create vector to hold number of elements and halo elements
 Vector<int> nhalo_elements(MPI_Helpers::Nproc);
 Vector<int> n_elements(MPI_Helpers::Nproc);

 // Count total number of halo elements
 unsigned count=0;
 for (int d=0;d<MPI_Helpers::Nproc;d++)
  {
   Vector<FiniteElement*> halo_elem_pt(halo_element_pt(d));
   count+=halo_elem_pt.size();
  }

 // Stick own number into appropriate entry
 nhalo_elements[MPI_Helpers::My_rank]=count;
 n_elements[MPI_Helpers::My_rank]=nelement();

 // Gather information on root processor: First argument group
 // specifies what is to be sent (one int from each procssor, indicating
 // the number of elements on it), the second group indicates where
 // the results are to be gathered (in rank order) on root processor.
 MPI_Gather(&nhalo_elements[MPI_Helpers::My_rank],1,MPI_INT,
            &nhalo_elements[0],1, MPI_INT,
            0,MPI_COMM_WORLD);
 MPI_Gather(&n_elements[MPI_Helpers::My_rank],1,MPI_INT,
            &n_elements[0],1, MPI_INT,
            0,MPI_COMM_WORLD);

 // Initialise stats
 av_efficiency=0.0;
 double max=-1.0;
 double min=1000000000.0;

 if (MPI_Helpers::My_rank==0)
  {
   for (int i=0;i<MPI_Helpers::Nproc;i++)
    {
     double eff=double(n_elements[i]-nhalo_elements[i])/double(n_elements[i]);
     av_efficiency+=eff;
     if (eff>max) max=eff;
     if (eff<min) min=eff;
       
    }
   av_efficiency/=double(MPI_Helpers::Nproc);
  }

 // Now broadcast the result back out
 MPI_Bcast(&max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 MPI_Bcast(&min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 MPI_Bcast(&av_efficiency,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

 max_efficiency=max;
 min_efficiency=min;

}



//========================================================================
/// Doc the mesh distribution -- hierher extend to multiple values
/// per Data object.
//========================================================================
void Mesh::doc_mesh_distribution(DocInfo& doc_info)
{ 

 char filename[100];
 std::ofstream some_file;

 // Doc elements on this processor
 sprintf(filename,"%s/elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 this->output(some_file,5);
 some_file.close();

 // Doc non-halo elements on this processor
 sprintf(filename,"%s/non_halo_elements_on_proc%i_%i.dat",
         doc_info.directory().c_str(),
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);

   // Get to elements on processor
   unsigned nelem=this->nelement();
//   oomph_info
//    << "Processor " << MPI_Helpers::My_rank << " holds " << nelem 
//    << "non-halo elements " << std::endl;
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
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
      it!=Processors_that_share_halos.end();it++)
  {
   unsigned domain=*it;

   // Get vector of halo elements by copy operation
   Vector<FiniteElement*> halo_elem_pt(this->halo_element_pt(domain));
   unsigned nelem=halo_elem_pt.size();
//   oomph_info
//    << "Processor " << MPI_Helpers::My_rank << " holds " << nelem 
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
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
      it!=Processors_that_share_halos.end();it++)
  {
   unsigned domain=*it;
   // Get vector of haloed elements by copy operation
   Vector<FiniteElement*> haloed_elem_pt(this->haloed_element_pt(domain));
   unsigned nelem=haloed_elem_pt.size();
//   oomph_info
//    << "Processor " << MPI_Helpers::My_rank << " holds " << nelem
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
         MPI_Helpers::My_rank,doc_info.number());
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
      << nod_pt->processor_in_charge() << " " 
      << nod_pt->is_halo() << " " 
      << nod_pt->eqn_number(0) << " "   // these two won't work for SolidNodes
      << nod_pt->is_pinned(0) << " "    // with eqn numbers for position only
      << std::endl;
    }
  }
 some_file.close();
 
 // Doc solid nodes on this processor
 sprintf(filename,"%s/solid_nodes_on_proc%i_%i.dat",
    doc_info.directory().c_str(),MPI_Helpers::My_rank,doc_info.number());
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
      << solid_nod_pt->processor_in_charge() << " " 
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
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
      it!=Processors_that_share_halos.end();it++)
  {
   unsigned domain=*it;
   unsigned nnod=this->nhalo_node(domain);
//   oomph_info << "Processor " << MPI_Helpers::My_rank 
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
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 for (std::set<unsigned>::iterator it=Processors_that_share_halos.begin();
      it!=Processors_that_share_halos.end();it++)
  {
   unsigned domain=*it;
   unsigned nnod=this->nhaloed_node(domain);
//   oomph_info << "Processor " << MPI_Helpers::My_rank 
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
         MPI_Helpers::My_rank,doc_info.number());
 some_file.open(filename);
 this->output(some_file,5);
 some_file.close();
 
 
 // Doc boundary scheme
 sprintf(filename,"%s/boundaries%i_%i.dat",
         doc_info.directory().c_str(),
         MPI_Helpers::My_rank,doc_info.number());
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
           MPI_Helpers::My_rank,b,doc_info.number());
   some_file.open(filename);
   unsigned nelem=this->nboundary_element(b);
   for (unsigned e=0;e<nelem;e++)
    {
     this->boundary_element_pt(b,e)->output(some_file,5);
    }
   some_file.close();
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

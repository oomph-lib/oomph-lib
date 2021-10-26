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
// Non-inline functions for the PerturbedSpineNode, PerturbedSpineElement
// and PerturbedSpineMesh classes

// oomph-lib headers
#include "generic.h"
#include <cstdlib>

namespace oomph
{



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
// Functions for the PerturbedSpineNode class
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=======================================================================
/// Update function: Call the update function in the Node's
/// PerturbedSpineMesh
//=======================================================================
void PerturbedSpineNode::node_update(
 const bool& update_all_time_levels_for_new_node)
{
 PerturbedSpine_mesh_pt->perturbed_spine_node_update(this);

 // Perform any auxiliary updates (i.e. reseting boundary conditions)
 if(Aux_node_update_fct_pt != 0)
  {
   Aux_node_update_fct_pt(this);
  }
}



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
// Functions for the PerturbedSpineElement class
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=======================================================================
/// Construct and fill the node_update_data vector
//=======================================================================
template<class ELEMENT>
void PerturbedSpineElement<ELEMENT>::complete_setup_of_dependencies()
{
 // Call function of underlying element
 ElementWithSpecificMovingNodes<ELEMENT,PerturbedSpineNode>::
  complete_setup_of_dependencies();

 // Sort out the perturbed spine index stuff
 // Find the data that correspond to spine heights
 {
  const unsigned n_node = this->nnode();

  // Allocate memory
  if(PerturbedSpine_geometric_index)
   {
    delete[] PerturbedSpine_geometric_index;
   }
  PerturbedSpine_geometric_index = new unsigned[n_node];
  
  // Now loop over data and find out where it fits
  for(unsigned n=0;n<n_node;n++)
   {
    // Find pointer to the perturbed spine
    PerturbedSpine* const perturbedspine_pt =
     static_cast<PerturbedSpineNode*>(this->node_pt(n))->perturbed_spine_pt();

    // If there is a spine then find the pointer to the data
    if(perturbedspine_pt)
     {
      // Find the pointer to the data
      Data* perturbed_spine_height_data_pt =
       perturbedspine_pt->height_pt();
      
      // Now find the index of the corresponding perturbed spine
      const unsigned n_node_update_data = this->ngeom_data();
      for(unsigned i=0;i<n_node_update_data;i++)
       {
        if(this->Geom_data_pt[i]==perturbed_spine_height_data_pt)
         {
          PerturbedSpine_geometric_index[n] = i;
          break;
         }
       }
     }
    // Otherwise issue a warning
    else
     {
      // Set the spine_geometric_index out of range,
      // which will cause the spine_local_eqn to return a pinned value
      PerturbedSpine_geometric_index[n] = this->ngeom_data();
     }
   }
 }
 
} // End of complete_setup_of_dependencies



/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
// Functions for the PerturbedSpineMesh class
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//=======================================================================
/// Update function to update all nodes of mesh.
/// [Doesn't make sense to use this mesh with SolidElements anyway,
/// so we buffer the case if update_all_solid_nodes (which defaults
/// to false) is set to true.]
//=======================================================================
void PerturbedSpineMesh::node_update(const bool& update_all_solid_nodes)
{
#ifdef PARANOID
 if (update_all_solid_nodes)
  {
   std::string error_message =
    "Doesn't make sense to use an PerturbedSpineMesh with\n";
   error_message +=
    "SolidElements so specifying update_all_solid_nodes=true\n";
   error_message += "doesn't make sense either\n";
   
   throw OomphLibError(error_message,
                       "PerturbedSpineMesh::node_update()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Loop over all the nodes
 unsigned long Node_pt_range = Node_pt.size();
 for(unsigned long l=0;l<Node_pt_range;l++)
  {
#ifdef PARANOID
   if(!dynamic_cast<PerturbedSpineNode*>(Node_pt[l]))
    {
     std::ostringstream error_stream;
     error_stream  
      << "Error: Node " << l << "is a " << typeid(Node_pt[l]).name() 
      << ", not a PerturbedSpineNode" << std::endl;
     throw OomphLibError(error_stream.str(),
                         "PerturbedSpineMesh::node_update()",
                         OOMPH_EXCEPTION_LOCATION);
    } 
#endif
   
   // Need to cast to spine node to get to update function
   dynamic_cast<PerturbedSpineNode*>(Node_pt[l])->node_update();
  }
}



//=======================================================================
/// Assign (global) equation numbers to spines, nodes and elements
//=======================================================================
unsigned long PerturbedSpineMesh::assign_global_eqn_numbers
(Vector<double *> &Dof_pt)
{
 // First call the general mesh assign_eqn_numbers
 unsigned long equation_number = Mesh::assign_global_eqn_numbers(Dof_pt);

 // Loop over spines and set global equation numbers for the spine heights
 // (they are the only Data items whose global eqn numbers are assigned
 // here)
 unsigned long PerturbedSpine_pt_range = PerturbedSpine_pt.size();
 for(unsigned long i=0;i<PerturbedSpine_pt_range;i++)
  {
   PerturbedSpine_pt[i]->height_pt()->assign_eqn_numbers(equation_number,
                                                         Dof_pt);
  }

 // Return total number of equations
 return(equation_number);
}



//=======================================================================
/// Overload the dump function so that the spine data is also dumped
//=======================================================================
  void PerturbedSpineMesh::dump(std::ofstream &dump_file, 
				const bool &use_old_ordering) const
{
 // Call the standard mesh dump function
  Mesh::dump(dump_file,use_old_ordering);

 // Now loop over the spine data and dump the spine height data
 // The ASSUMPTION is that the geometric data is stored elsewhere and will
 // be dumped elsewhere

 // Find the number of perturbed spines
 const unsigned long n_spine = nspine();

 // Doc number of spines
 dump_file << n_spine << " # number of spines " << std::endl;

 // Loop over the perturbed spines
 for(unsigned long s=0;s<n_spine;s++)
  {
   perturbed_spine_pt(s)->height_pt()->dump(dump_file);
  }
}



//========================================================================
/// Overload the read function so that the spine data is also read
//========================================================================
void PerturbedSpineMesh::read(std::ifstream &restart_file)
{
 // Call the standard mesh read function
 Mesh::read(restart_file);

 // Now loop over the spine data and dump the spine height data
 // The ASSUMPTION is that the geometric data is stored elsewhere and will
 // be dumped elsewhere

 // Get the number of perturbed spines
 const unsigned long n_spine = nspine();

 std::string input_string;
 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 // Ignore the restr of the line
 restart_file.ignore(80,'\n');
 
 // Check the number of spines
 const unsigned long check_n_spine = atoi(input_string.c_str());

 if(check_n_spine != n_spine)
  {
   std::ostringstream error_stream;
    error_stream 
     << "Number of spines in the restart file, " << check_n_spine 
     << std::endl << "does not equal the number of spines in the mesh " 
     << n_spine << std::endl;

    throw OomphLibError(error_stream.str(),
                        "PerturbedSpineMesh::read()",
                        OOMPH_EXCEPTION_LOCATION);
  }

 // Loop over the spines and read the data
 for(unsigned long s=0;s<n_spine;s++)
  {
   perturbed_spine_pt(s)->height_pt()->read(restart_file);
  }
}


} // End of oomph namespace

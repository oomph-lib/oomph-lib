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
//Demo driver for adaptive bifurcation detection
//based on flow of a fluid in a porous channel.
//There is a pitchfork bifurcation at Re 9.3
//Care must be taken to ensure that the adapted eigenfunction has
//the correct symmetries otherwise the problem converges to the wrong
//solution (or fails to converge entirely).

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/rectangular_quadmesh.h"

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
 // Reynolds number
 double Re=9.3;
 
 // Length of the channel
 double L = 8.0;
} // end of GPV namespace


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//============================================================================
/// A class to do comparison of the nodes based on
/// lexicographical ordering in the Cartesian coordinates.
/// This is required to symmetrise eigenfunction (which is a bit of a hack)
//============================================================================
 class CompareNodes
 {
 public:
  
  /// Define the comparison operator using the round brackets
  int operator() (Node* const &node1_pt,
                  Node* const &node2_pt)
   {
    //Find the dimension of the node
     unsigned n_dim = node1_pt->ndim();
     //Check that the dimensions of both nodes are the same
     {
      unsigned n_dim2 = node2_pt->ndim();
      //If not issue a warning
      if(n_dim != n_dim2)
       {
        std::ostringstream warn_message;
        warn_message
         << "Warning: Two nodess do not have the same dimension"
         << n_dim << " and " << n_dim2 << ". They cannot be compared\n";
        OomphLibWarning(warn_message.str(),
                        "CompareNodes::()",
                        OOMPH_EXCEPTION_LOCATION);
       }
     }
     
     //Get the global coordinates of the nodes
     Vector<double> x1(n_dim), x2(n_dim);
     for(unsigned i=0;i<n_dim;i++)
      {
       x1[i] = node1_pt->x(i);
       x2[i] = node2_pt->x(i);
      }

     //This is the Stroustrup-approved way to do lexicographical ordering
     //Loop over the components until they are not equal
     //to within a given tolerance
     {
      unsigned i=0; double tol = 1.0e-14;
      while(i!=n_dim && (std::abs(x1[i]-x2[i]) < tol)){ ++i;}
      //If we've reached the end, the coordinates are equal, return false
      if(i==n_dim) {return 0;}
      //Otherwise, return the ordering on the final component
      return x1[i] < x2[i];
     }
    }
  };
 

//==start_of_problem_class============================================
/// Porous channel flow on a refineable mesh
//====================================================================
template<class ELEMENT>
class RefineablePorousChannelProblem : public Problem
{

public:

 /// Make a copy of the problem for using in adaptive bifurcation tracking
 Problem* make_copy()
  {
   //Make a copy based on the current parameters
   return(new RefineablePorousChannelProblem<ELEMENT>());
  }

 /// Constructor
 RefineablePorousChannelProblem();

 /// Destructor to clean up memory
 ~RefineablePorousChannelProblem();

 /// Set the boundary conditions
 void set_boundary_conditions();

 /// No actions are required after the change in bifurcation parameter.
 /// This overloads the default which calls actions before and after 
 /// newton solve.
 void actions_after_change_in_bifurcation_parameter() {}
 
 /// Update the boundary conditions before next timestep: 
 void actions_before_implicit_timestep() {set_boundary_conditions();}

 /// Hacky function to symmetrise the problem
 /// The idea is to sort all nodes in the mesh lexicographically
 /// Then we can loop over all nodes with the same x value and
 /// apply the appropriate symmetry conditions.
 void symmetrise_eigenfunction_for_adaptive_pitchfork_tracking()
  {
   //Find the number of nodes
   const unsigned n_node = this->mesh_pt()->nnode();
   //Allocate storage for the nodes
   Vector<Node*> local_nodes_pt(n_node);
   //Load the nodes into the local copy
   for(unsigned n=0;n<n_node;n++)
    {
     local_nodes_pt[n] = this->mesh_pt()->node_pt(n);
    }

   //Now let's sort the nodes lexicographically
   std::sort(local_nodes_pt.begin(),local_nodes_pt.end(),
             CompareNodes());

   //The nodes are now sorted, we proceed to find all nodes 
   //in a column (same x location)
   //Storage for start and end of the column
   unsigned column_start=0, column_end = 0;
   //Find the x-location of the first node
   double x = local_nodes_pt[0]->x(0);
   //Find the number of values stored at the nodes 
   //(ASSUMED TO BE THE SAME FOR ALL NODES)
   unsigned n_value = local_nodes_pt[0]->nvalue();
   //Specify the symmetries of the u,v and p components
   int symm[3] = {-1,1,-1};

   //Loop over all other nodes
   for(unsigned n=1;n<n_node;n++)
    {
     //If the x-location has changed, then we
     //know all nodes in the column
     if(std::abs(local_nodes_pt[n]->x(0)-x) > 1.0e-14)
      {
       //Set the end of the previous column
       column_end = n-1;
       //Find the number of entries in the column
       unsigned n_entries = (column_end - column_start + 1);
       
       //Storage for "half" the number of nodes
       unsigned n_half=0;
       //Boolean to indicate whether there is an odd (or even) number of nodes
       bool odd = false;

       //Is the number of entries odd  (remainder mod 2 is non-zero)
       if(n_entries%2)
        {
         //Find the number of nodes in the half-domain
         n_half = (n_entries - 1)/2;
         odd = true;
        }
       //Otherwise the number is even
       else
        {
         n_half = n_entries/2;
        }

       //Loop over the variables stored at every node
       for(unsigned i=0;i<n_value;i++)
        {
         //If the symmetry is odd
         //The idea is to take the average of the sum of the absolute values 
         //located symmetrically about the centreline and set the new values
         //to be +average and -average, preserving the signs.
         if(symm[i]==-1)
          {
           //Set the middle to zero if we have an odd number of nodes
           if(odd) {local_nodes_pt[column_start+n_half]->set_value(i,0.0);}
           //Loop over the half-width
           for(unsigned n=0;n<n_half;n++)
            {
             //Get the values of the node
             Node* node1_pt = local_nodes_pt[column_start+n];
             Node* node2_pt = local_nodes_pt[column_end-n];
             double val1 = node1_pt->value(i);
             double val2 = node2_pt->value(i);

             //If the "lower value" is negative,
             //the "upper value" should be positive
             if(val1 < 0)
              {
               //Find the average of the absolute values
               //and set new values preserving the sign
               double average = 0.5*(val2-val1);
               val1 = -average; val2 = average;
              }
             //Otherwise the "lower value" is positive
             //the "upper value" is negative
             else if(val1 > 0)
              {
               //Find the average of the absolute values
               //and set new values preserving the sign
               double average = 0.5*(val1 - val2);
               val1 = average; val2 = -average;
              }
             //Otherwise both values should be zero
             else
              {
               val1 = val2 = 0.0;
              }

             //Set the values stored at the nodes
             node1_pt->set_value(i,val1);
             node2_pt->set_value(i,val2);
            }
          }
         //Otherwise the required symmetry is even
         //The idea is to take the simple average of the two values
         //located symmetrically about the midplane
         else
          {
           //Loop over the half-width
           for(unsigned n=0;n<n_half;n++)
            {
             //Get the values stored at the node
             Node* node1_pt = local_nodes_pt[column_start+n];
             Node* node2_pt = local_nodes_pt[column_end-n];
             double val1 = node1_pt->value(i);
             double val2 = node2_pt->value(i);
             //Find the average of the two values
             double average = 0.5*(val1+val2);
             //Set the new values to the average
             val1 = val2 = average;
             //Now actually set the values stored at the node
             node1_pt->set_value(i,val1);
             node2_pt->set_value(i,val2);
            }
          }
        }
       
       //Finally, set the new column start
       column_start = n; column_end=n;
       //Set the new x-location of the column
       x = local_nodes_pt[column_start]->x(0);
      }
    }
  }


 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   //Reset the boundary conditions
   set_boundary_conditions();
  }
 
 
 // Access function for the specific mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RefineablePorousChannel problem 
//========================================================================
template<class ELEMENT>
RefineablePorousChannelProblem<ELEMENT>::RefineablePorousChannelProblem()
{
 //Set the maximum residuals to be large so that the initial 
 //(coarse) problem converges
 Max_residuals = 1000.0;
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timesteps. 
 add_time_stepper_pt(new BDF<2>);
 
 //Number of elements in one unit
 unsigned N=2;
 
 // Domain length in r-direction
 double l_x=Global_Physical_Variables::L;
 
 //Number of elements across
 unsigned n_y = 2*N;

 // # of elements in x-direction
 unsigned n_x= static_cast<unsigned>(l_x)*N;
 
 // Build and assign rectangular mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,0.0,l_x,-1.0,1.0, 
                                             time_stepper_pt());
 

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 //Let this problem be conventional form by setting gamma to zero
 ELEMENT::Gamma[0] = 0.0; //x-momentum
 ELEMENT::Gamma[1] = 0.0; //y-momentum
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 
 unsigned num_bound = mesh_pt()->nboundary();
 
 // Pin all two velocities on boundaries 0 and 2
 for(unsigned ibound=0;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u/v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries 0 and 2
  
 // Now pin the x velocity on boundaries 3
 {
  unsigned ibound = 3;
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Loop over the theta- and phi-velocities
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
   }
 } 
 // end of set boundary conditions

 // Complete the build of all elements so they are fully functional
 //================================================================

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
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;
   
   //The Mesh does not move
   el_pt->disable_ALE();
  } // end loop over elements

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=========start of actions_before_implicit_timestep======================
/// Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void RefineablePorousChannelProblem<ELEMENT>::set_boundary_conditions()
{
 // Get current time
 const double time=time_pt()->time();


 //Find number of nodes on top row
 unsigned n_node = mesh_pt()->nboundary_node(2);

 //Set the value of a perturbation
 double epsilon = 0.1;

 //Loop over the nodes on the top row and set their u and v velocities
 //We are going to add a small time-dependent perturbation
 //to the velocity and also allow the posibility for an initial 
 //symmetry-breaking perturbation
 for(unsigned i=0;i<n_node;i++)
  {
   //Top row
   double u = 1.0 + epsilon*time*exp(-time);
   
   mesh_pt()->boundary_node_pt(2,i)->set_value(0,0.0);
   mesh_pt()->boundary_node_pt(2,i)->set_value(1,u);
  }

 //Find number of nodes on bottom row
 n_node = mesh_pt()->nboundary_node(0);
 
 for(unsigned i=0;i<n_node;i++)
  {
   //Bottom row
   double u =  -1.0 + epsilon*time*exp(-time);

   mesh_pt()->boundary_node_pt(0,i)->set_value(0,0.0);
   mesh_pt()->boundary_node_pt(0,i)->set_value(1,u);
  }

} // end of actions_before_implicit_timestep

//==start_of_destructor===================================================
/// Destructor for RefineablePorousChannel problem 
//========================================================================
template<class ELEMENT>
RefineablePorousChannelProblem<ELEMENT>::~RefineablePorousChannelProblem()
{ 
 //Kill the error estimator
 delete mesh_pt()->spatial_error_estimator_pt();
 //Kill the mesh
 delete mesh_pt();
 //Kill the timestepper
 delete time_stepper_pt();
} // end_of_destructor


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for RefineablePorousChannel test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{
 // Build the problem with QTaylorHoodElements
 RefineablePorousChannelProblem<RefineableQTaylorHoodElement<2> > 
  problem;

 problem.set_boundary_conditions();
 problem.steady_newton_solve();

 
 //Assign memory for the eigenvalues and eigenvectors
 Vector<std::complex<double> > eigenvalues;
 Vector<DoubleVector> eigenvectors;
 //Set the eigen solver to the LAPACK version which is
 //in every distributioin of oomph-lib
 problem.eigen_solver_pt() = new LAPACK_QZ;
 //Solve the eigenproblem
 problem.solve_eigenproblem(4,eigenvalues,eigenvectors);  

 //Find the eigenvalue with greatest real part
 unsigned ev_crit_index=0;  double ev_crit_value=0.0;
 
 //Loop over the values until we get a finite one
 for(unsigned ev=0;ev<eigenvalues.size();ev++)
  {
   double ev_crit_value = eigenvalues[ev].real();
   if(std::isfinite(ev_crit_value)) {ev_crit_index = ev; break;}
  }
 
 //Now loop over the other values and find the eigenvalue with 
 //greatest real part
 for(unsigned ev=ev_crit_index+1;ev<eigenvalues.size();ev++)
  {
   //Get the current value
   double real_part = eigenvalues[ev].real();
   //If it's finite do the comparison
   if(std::isfinite(real_part))
    {
     if(real_part > ev_crit_value)
      {
       ev_crit_value = real_part;
       ev_crit_index = ev;
      }
    }
  }

 problem.activate_pitchfork_tracking(&Global_Physical_Variables::Re,
                                     eigenvectors[ev_crit_index]);

 //Solve with two rounds of adaptivity
 problem.steady_newton_solve(2);

 //Output the solution at the bifurcation point
 problem.mesh_pt()->output("bif_soln.dat",5);

 //Report the final answer
 std::cout << "Pitchfork found at " << Global_Physical_Variables::Re  
           << "\n";
 std::cout << "The slack parameter is " << problem.dof(problem.ndof()-1) 
           << "\n";

 //Write the final answer into a file
 std::ofstream trace("trace.dat");
 trace << Global_Physical_Variables::Re << " "
       << problem.dof(problem.ndof()-1) << "\n";
 trace.close();

 //Deactivate the bifurcation tracking (for neatness)
 problem.deactivate_bifurcation_tracking();

 //Delete the eigensolver that we allocated
 delete problem.eigen_solver_pt();
 
} // end_of_main

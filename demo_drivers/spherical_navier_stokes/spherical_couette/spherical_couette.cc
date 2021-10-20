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
//Driver for flow between two concentric spheres, with one of the 
//spheres rotating. The Reynolds number is based on the radius
//of the outer sphere, which is therefore fixed at a radius of one.
//The torque on the inner sphere is calculated and is in good agreement
//with published data of Ni and Nigro (1994) among others.


//Generic includes
#include "generic.h"
#include "spherical_navier_stokes.h"

//Standard rectangular mesh
#include "meshes/rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 // Reynolds number (default is zero --- Stokes flow)
 double Re=0;
} // end of GPV namespace


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=====================================================================
/// Class of face elements whose sole raison d'etre is to calculate the
/// torque on the boundary of the selected sphere
//=====================================================================
template <class ELEMENT>
class TorqueCalculationElement : public virtual FaceGeometry<ELEMENT>, 
 public virtual FaceElement 
{

public:

 ///Constructor, which takes a "bulk" element and the value of the index
 ///and its limit
 TorqueCalculationElement(FiniteElement* const &element_pt, 
                          const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
  }

 /// Destructor should not delete anything
 ~TorqueCalculationElement() { }

 /// Return the torque calculation
 double calculate_torque()
  {
   //Set the value of n_intpt
   const unsigned n_intpt = integral_pt()->nweight();

   //Get the dimension of the element
   const unsigned dim = this->dim();

   //Get the storage for the local coordinate
   Vector<double> s(dim);

   //Get the position in the bulk
   Vector<double> s_bulk(dim+1);

   //Get the outer unit normal
   Vector<double> N(2,0.0);

   //Need to find position 
   Vector<double> x(dim+1);

   //Need the traction as well
   Vector<double> traction(3);

   double sum = 0.0;

   ELEMENT* const bulk_elem_pt = 
    dynamic_cast<ELEMENT*>(this->bulk_element_pt());

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the local coordinates
     for(unsigned i=0;i<dim;i++) {s[i] = integral_pt()->knot(ipt,i);}

     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Find the shape and test functions and return the Jacobian
     //of the mapping
     double J = J_eulerian(s);

     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Now get the position in the bulk
     this->get_local_coordinate_in_bulk(s,s_bulk);

     //Get the position
     this->interpolated_x(s,x);

     //Get the normal
     this->outer_unit_normal(s,N);

     //Now get the traction from the bulk element
     bulk_elem_pt->get_traction(s_bulk,N,traction);

     //Now sort out the magnitude of the torque the phi component
     //of the stress about the distance from the axis of rotation,
     // r*sin(theta)
     sum += x[0]*sin(x[1])*traction[2]*W*x[0]*x[0]*sin(x[1]);
    }

   return sum;
  }

}; 




//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain - time dependent version
//====================================================================
template<class ELEMENT>
class RefineableSphericalCouetteProblem : public Problem
{

public:


 /// Constructor
 RefineableSphericalCouetteProblem();

 /// Destructor to clean up memory
 ~RefineableSphericalCouetteProblem();


 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure



 ///  Set the boundary conditions
 void set_boundary_conditions();

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableSphericalNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableSphericalNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   //Reset the boundary conditions
   set_boundary_conditions();
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);
  }
 
 
 // Access function for the specific mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Document the solution
 void doc_solution(DocInfo& doc_info, std::ofstream&);

 /// Compute the torque on the sphere
 double compute_torque();

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RefineableSphericalCouette problem 
//========================================================================
template<class ELEMENT>
RefineableSphericalCouetteProblem<ELEMENT>::RefineableSphericalCouetteProblem()
{ 
 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------

 // pi definition
 const double pi = MathematicalConstants::Pi;
     
 // # of elements in r-direction
 unsigned n_r=4;

 // # of elements in theta-direction
 unsigned n_theta=4;

 // Radius of inner sphere
 double R_inner = 0.5;

 // Radius of outer sphere
 double R_outer=1.0;
 
  // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_r,n_theta,R_inner,R_outer,
                                             0.0,pi);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 
 // Pin all three velocities on boundaries 1 and 3
 for(unsigned ibound=1;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u/v/w velocities)
     for (unsigned i=0;i<3;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries 1 and 3
  
  // Now pin the theta and phi velocities on boundaries 0 and 2
   for(unsigned ibound=0;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over the theta- and phi-velocities
     for (unsigned i=1; i<3; i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
      }
    }
  } // end loop over boundaries 0 and 2
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

   //Disable ALE
   el_pt->disable_ALE();
   
  } // end loop over elements

 // Pin redudant pressure dofs
 RefineableSphericalNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
} // end_of_constructor


//=========start of set_boundary_conditions===============================
///  Set the boundary conditions so that the inner sphere has
/// a constant angular rotation of angular velocity one.
//========================================================================
template<class ELEMENT>
void RefineableSphericalCouetteProblem<ELEMENT>::set_boundary_conditions()
{
 //Setting for boundary 0 - zero theta and phi velocities
 unsigned ibound=0;
  
  // Loop over the nodes on boundary
  unsigned num_nod=mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     nod_pt->set_value(1,0.0);
     nod_pt->set_value(2,0.0);
    }
  
  
  //Set velocity for boundary 1 - outer wall (zero)
  ibound=1;
 
   // Loop over the nodes on boundary
   num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     // Get current values of the boundary conditions from the
     // exact solution
     nod_pt->set_value(0,0.0);
     nod_pt->set_value(1,0.0);
     nod_pt->set_value(2,0.0);
    }
    
    
  //Setting for boundary 2 - zero theta and phi velocities
  ibound=2;
  
   // Loop over the nodes on boundary
   num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     nod_pt->set_value(1,0.0);
     nod_pt->set_value(2,0.0);
    }
    
  
 //Setting for boundary 3 (inner sphere that is driven)
 ibound=3;
  
   // Loop over the nodes on boundary
   num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     nod_pt->set_value(0,0.0);
     nod_pt->set_value(1,0.0);
     //The velocity is sin(theta)
     nod_pt->set_value(2,x[0]*sin(x[1]));
    }
  
  
} // end of actions_before_implicit_timestep


//==start_of_destructor===================================================
/// Destructor for RefineableSphericalCouette problem 
//========================================================================
template<class ELEMENT>
RefineableSphericalCouetteProblem<ELEMENT>::~RefineableSphericalCouetteProblem()
{ 
 //Delete the error estimator
 delete mesh_pt()->spatial_error_estimator_pt();

 //Clean up the memory allocated for the mesh
 delete Problem::mesh_pt();

} // end_of_destructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableSphericalCouetteProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream&)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 
 
 
 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for Re =" << Global_Physical_Variables::Re 
      << std::endl;
 cout << "=================================================" << std::endl;


 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%g.dat",doc_info.directory().c_str(),
         Global_Physical_Variables::Re);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 
 some_file.close();
} // end_of_doc_solution



 /// Compute the torque on the sphere
template<class ELEMENT>
double RefineableSphericalCouetteProblem<ELEMENT>::compute_torque()
{
 unsigned bound = 3;
 //Loop over the elements adjacent to the boundary 3 and make face elements
 unsigned n_bound_element = this->mesh_pt()->nboundary_element(bound);

 double torque = 0.0;

 for(unsigned e=0;e<n_bound_element;e++)
  {
   //Get pointer to the bulk element
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    this->mesh_pt()->boundary_element_pt(bound,e));

   //FInd the face index
   int face_index = this->mesh_pt()->face_index_at_boundary(bound,e);
   
   //Build the flux element
   TorqueCalculationElement<ELEMENT>* torque_element_pt = new
    TorqueCalculationElement<ELEMENT>(bulk_elem_pt,face_index);

   //Now calculate the torque
   torque += torque_element_pt->calculate_torque();

   //Delete our element (it's work is done)
   delete torque_element_pt;
  }

 return torque;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for RefineableSphericalCouette test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");

 // ---------------
 // end of Set up doc info
 

 // Doing QCrouzeixRaviartElements
 {
  // Build the problem with QCrouzeixRaviartElements
  RefineableSphericalCouetteProblem<
   RefineableQSphericalCrouzeixRaviartElement > 
   problem;
  cout << "Doing QCrouzeixRaviartElement" << std::endl;
  
  // Open a trace file
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);
  
  // Over-ride the maximum and minimum permitted errors
  problem.mesh_pt()->max_permitted_error() = 1.0e-3; //Default = 1.0e-3
  problem.mesh_pt()->min_permitted_error() = 1.0e-5; //Default = 1.0e-5
  
  // Over-ride the maximum and minimum permitted refinement levels
  problem.mesh_pt()->max_refinement_level() = 4;//maximum_ref_level;
  problem.mesh_pt()->min_refinement_level() = 1;//minimum_ref_level;
  
  //Set the boundary conditions
  problem.set_boundary_conditions();
  
  //Set the maximum adaptation
  unsigned max_adapt = 10;

  for(unsigned i=0;i<5;i++)
   {
    // Solve the problem
    problem.steady_newton_solve(max_adapt);

    //Scale the torque to be consistent with the non-dimensionalisation
    //in Dennis, Singh & Ingham (1980) and many others
    double torque = 3.0*problem.compute_torque()/4.0;
    std::cout << "Torque is " 
              << torque
              << "\n";
    
    problem.doc_solution(doc_info,trace_file);
    
    trace_file << Global_Physical_Variables::Re << " " 
               << " " << torque
               << std::endl;

      Global_Physical_Variables::Re += 25.0;
   }
  
  // Close trace file
  trace_file.close();
  
 } // end of QCrouzeixRaviartElements


} // end_of_main

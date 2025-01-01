//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Driver for flow past odd-shaped obstacle -- domain meshed with triangle.
// This is a warm-up problem for an unstructured fsi problem.

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/triangle_mesh.h"


using namespace std;

using namespace oomph; 


//===========start_mesh====================================================
/// Triangle-based mesh upgraded to become a (pseudo-) solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticTriangleMesh : public virtual TriangleMesh<ELEMENT>,
                            public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 ElasticTriangleMesh(const std::string& node_file_name,
                     const std::string& element_file_name,
                     const std::string& poly_file_name,
                     TimeStepper* time_stepper_pt=
                     &Mesh::Default_TimeStepper) : 
 TriangleMesh<ELEMENT>(node_file_name,element_file_name,
                       poly_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
   
   //Identify the domain boundaries
   this->identify_boundaries();
  }


 /// Function used to identify the domain boundaries
 void identify_boundaries()
  {
   // We will have three boundaries
   this->set_nboundary(3);

   unsigned n_node=this->nnode();
   for (unsigned j=0;j<n_node;j++)
    {
     Node* nod_pt=this->node_pt(j);

     // Boundary 1 is left (inflow) boundary
     if (nod_pt->x(0)<0.226)
      {
       // If it's not on the upper or lower channel walls remove it
       // from boundary 0
       if ((nod_pt->x(1)<4.08)&&(nod_pt->x(1)>0.113))
        {
         this->remove_boundary_node(0,nod_pt);
        }
       this->add_boundary_node(1,nod_pt);
      }

     // Boundary 2 is right (outflow) boundary
     if (nod_pt->x(0)>8.28)
      {
       // If it's not on the upper or lower channel walls remove it
       // from boundary 0
       if ((nod_pt->x(1)<4.08)&&(nod_pt->x(1)>0.113))
        {
         this->remove_boundary_node(0,nod_pt);
        }
       this->add_boundary_node(2,nod_pt);
      }
    }
   this->setup_boundary_element_info();
  }

 /// Empty Destructor
 virtual ~ElasticTriangleMesh() { }


};



//==start_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=0.0; 

 /// Pseudo-solid Poisson ratio
 double Nu=0.3;

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt=0;

} // end_of_namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Unstructured fluid problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem : public Problem
{

public:

 /// Constructor 
 UnstructuredFluidProblem();
 
 /// Destructor (empty)
 ~UnstructuredFluidProblem(){}

 /// Access function for the fluid mesh
 ElasticTriangleMesh<ELEMENT>*& fluid_mesh_pt() 
  {
   return Fluid_mesh_pt;
  }

 /// Set the boundary conditions
 void set_boundary_conditions();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Fluid mesh
 ElasticTriangleMesh<ELEMENT>* Fluid_mesh_pt;


}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor 
//========================================================================
template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
 
 //Create fluid mesh
 string fluid_node_file_name="fluid.fig.1.node";
 string fluid_element_file_name="fluid.fig.1.ele";
 string fluid_poly_file_name="fluid.fig.1.poly"; 
 Fluid_mesh_pt = new ElasticTriangleMesh<ELEMENT>(fluid_node_file_name,
                                                  fluid_element_file_name,
                                                  fluid_poly_file_name);

 //Set the boundary conditions
 this->set_boundary_conditions();
 
 // Add fluid mesh
 add_sub_mesh(fluid_mesh_pt());
 
 // Build global mesh
 build_global_mesh();
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor

 

//============start_of_set_boundary_conditions=============================
/// Set the boundary conditions for the problem and document the boundary
/// nodes
//==========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::set_boundary_conditions()
{ 
 // Doc pinned nodes
 std::ofstream solid_bc_file("pinned_solid_nodes.dat");
 std::ofstream u_bc_file("pinned_u_nodes.dat");
 std::ofstream v_bc_file("pinned_v_nodes.dat");

 // Set the boundary conditions for fluid problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);

     // Pin velocity everywhere apart from outlet where we
     // have parallel outflow
     if (ibound!=2)
      {
       // Pin it
       nod_pt->pin(0);
       // Doc it as pinned
       u_bc_file << nod_pt->x(0) << " " << nod_pt->x(1) << std::endl;
      }
     // Pin it
     nod_pt->pin(1);
     // Doc it as pinned
     v_bc_file << nod_pt->x(0) << " " << nod_pt->x(1) << std::endl;
     
     // Pin pseudo-solid positions everywhere
     for(unsigned i=0;i<2;i++)
      {
       // Pin the node
       SolidNode* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       nod_pt->pin_position(i);
       
       // Doc it as pinned
       solid_bc_file << nod_pt->x(i) << " ";
      }
     solid_bc_file << std::endl;    
    }
  } // end loop over boundaries

 // Close
 solid_bc_file.close();
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = fluid_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(fluid_mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
  } 


 // Apply fluid boundary conditions: Poiseuille at inflow
 //------------------------------------------------------

 // Find max. and min y-coordinate at inflow
 unsigned ibound=1;
 //Initialise y_min and y_max to y-coordinate of first node on boundary
 double y_min=fluid_mesh_pt()->boundary_node_pt(ibound,0)->x(1);
 double y_max=y_min;
 //Loop over the rest of the nodes
 unsigned num_nod= fluid_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=1;inod<num_nod;inod++)
  {
   double y=fluid_mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
   if (y>y_max)
    {
     y_max=y;
    }
   if (y<y_min)
    {
     y_min=y;
    }
  }
 double y_mid=0.5*(y_min+y_max);
 
 // Loop over all boundaries
 for (unsigned ibound=0;ibound<fluid_mesh_pt()->nboundary();ibound++)
  {
   unsigned num_nod= fluid_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parabolic inflow at the left boundary (boundary 1)
     if (ibound==1)
      {
       double y=fluid_mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
       double veloc=1.5/(y_max-y_min)*
        (y-y_min)*(y_max-y)/((y_mid-y_min)*(y_max-y_mid));
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,veloc);
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
     // Zero flow elsewhere apart from x-velocity on outflow boundary
     else
      {
       if(ibound!=2)
        {
         fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
    }
  }
} // end_of_set_boundary_conditions




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 fluid_mesh_pt()->output(some_file,npts);
 some_file.close();
} //end_of_doc_solution


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for unstructured fluid test problem
//=====================================================================
int main()
{
 // Label for output
 DocInfo doc_info;
 
 //Set the constitutive law for the pseudo-elasticity
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);

 //Taylor--Hood formulation
 {
  // Set output directory
  doc_info.set_directory("RESLT_TH");
  
  // Step number
  doc_info.number()=0;
  
  // Build the problem with TTaylorHoodElements
  UnstructuredFluidProblem<PseudoSolidNodeUpdateElement<
   TTaylorHoodElement<2>, 
   TPVDElement<2,3> > > problem;
  
  // Output boundaries 
  problem.fluid_mesh_pt()->output_boundaries("RESLT_TH/boundaries.dat");
  
  // Output the initial guess for the solution
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Parameter study
  double re_increment=5.0;
  unsigned nstep=2; // 10;
  for (unsigned i=0;i<nstep;i++)
   {
    // Solve the problem
    problem.newton_solve();
    
    // Output the solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    
    // Bump up Re
    Global_Physical_Variables::Re+=re_increment;
   }  //end_of_parameter_study
 }



 //Crouzeix--Raviart formulation
 {
 // Set output directory
 doc_info.set_directory("RESLT_CR");
 
 // Step number
 doc_info.number()=0;

 // Reset Reynolds number
 Global_Physical_Variables::Re = 0.0;
 
 // Build the problem with TCrouzeixRaviartElements
 UnstructuredFluidProblem<PseudoSolidNodeUpdateElement<
  TCrouzeixRaviartElement<2>, 
  TPVDBubbleEnrichedElement<2,3> > > problem;
 
 // Output boundaries 
  problem.fluid_mesh_pt()->output_boundaries("RESLT_CR/boundaries.dat");
  
  // Output the initial guess for the solution
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Parameter study
  double re_increment=5.0;
  unsigned nstep=2; // 10;
  for (unsigned i=0;i<nstep;i++)
   {
    // Solve the problem
    problem.newton_solve();
    
    // Output the solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    
    // Bump up Re
    Global_Physical_Variables::Re+=re_increment;
   }
 }


} // end_of_main











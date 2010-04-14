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
// Driver code for a simple unstructured fsi problem using meshes
// generated with VMTK.


//Generic libraries
#include "generic.h"
#include "solid.h"
#include "constitutive.h"
#include "navier_stokes.h"
#include "poisson.h"
#include "linear_elasticity.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;
using namespace oomph;

//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 
 double Lambda = 0.7;

 double Mu = 1.0;

 /// The elasticity tensor
 IsotropicElasticityTensor E(Lambda,Mu);
}

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// hierher
//================================================================
namespace hierher_namespace
{

 /// Poisson's ratio
 double Nu=0.0;


 /// Young's modulus
 double E=0.00001;

 /// Generalised Mooney Rivling coefficient
 double C1=0.01;



 StrainEnergyFunction* strain_energy_pt=new GeneralisedMooneyRivlin(&Nu,&C1,&E);


 /// Create constitutive law
 //ConstitutiveLaw* Constitutive_law_pt=new GeneralisedHookean(&Nu,&E);
 ConstitutiveLaw* Constitutive_law_pt=
  new IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_pt);

}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//====================================================================
/// Auxiliary Problem to smooth a SolidMesh by adjusting the internal
/// nodal positions by solving a LINEAR solid mechanics problem for the
/// nodal displacements between the specified displacements of certain
/// pinned nodes (usually located on boundaries). 
/// hierher turn into functor?
/// hierher specify faster Poisson solver?
//====================================================================
class LinearSmoothMeshProblem : public Problem
{

public:

 /// \short Constructor: Specify SolidMesh whose nodal positions are to 
 /// be adjusted, and set of nodes in that mesh whose position
 /// are to remain fixed.
 LinearSmoothMeshProblem(SolidMesh* orig_mesh_pt,
                         std::set<Node*> pinned_nodes)
  {
   // Create new mesh and read out node/element numbers from old one
   mesh_pt()=new Mesh;
   unsigned nelem=orig_mesh_pt->nelement();
   unsigned nnode=orig_mesh_pt->nnode();
   
   // Have we already created that node?
   std::map<Node*,Node*> new_node;
   
   // Create new elements
   for (unsigned e=0;e<nelem;e++)
    {
     
     // Make/add new element
     TLinearElasticityElement<3,3>* el_pt=new TLinearElasticityElement<3,3>;
     mesh_pt()->add_element_pt(el_pt);
     
     //Set the Reynolds number, etc
     el_pt->elasticity_tensor_pt() = &Global_Parameters::E;
     
     // Find corresponding original element
     SolidFiniteElement* orig_elem_pt=
      dynamic_cast<SolidFiniteElement*>(orig_mesh_pt->finite_element_pt(e));
     unsigned nnod=orig_elem_pt->nnode();
     
     // Create nodes
     for (unsigned j=0;j<nnod;j++)
      {
       // Does it not exist yet?
       if (new_node[orig_elem_pt->node_pt(j)]==0)
        {
         Node* new_nod_pt=mesh_pt()->finite_element_pt(e)->construct_node(j);
         new_node[orig_elem_pt->node_pt(j)]=new_nod_pt;
         mesh_pt()->add_node_pt(new_nod_pt);
         for (unsigned i=0;i<3;i++)
          {
           // Set new nodal position to be the old one in the
           // SolidMesh (assumed to contain no inverted elements)
           new_nod_pt->x(i)=
            dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))->xi(i);
          }         
        }
       // It already exists -- copy across
       else
        {
         mesh_pt()->finite_element_pt(e)->node_pt(j)=
          new_node[orig_elem_pt->node_pt(j)];
        }
      }
    }
   
   
   mesh_pt()->output("old_mesh.dat");
      
   // Loop over pinned nodes -- pin their positions and assign updated nodal 
   // positions
   double scale=1.0;
   for (std::set<Node*>::iterator it=pinned_nodes.begin();
        it!=pinned_nodes.end();it++)
    {
     for (unsigned i=0;i<3;i++)
      {
       new_node[*it]->pin(i);
       new_node[*it]->set_value(i,scale*
                                (dynamic_cast<SolidNode*>(*it)->x(i)-
                                 dynamic_cast<SolidNode*>(*it)->xi(i)));
      }
    }


   mesh_pt()->output("new_mesh.dat");

   oomph_info << "Number of equations for smoothing problem: " 
              << assign_eqn_numbers() << std::endl;
   
   
   // Solve
   double backup=Problem::Newton_solver_tolerance;
   //Newton_solver_tolerance=1.0e29;
   newton_solve();
   Problem::Newton_solver_tolerance=backup;
 
   mesh_pt()->output("linear_soln.dat");

   // Loop over nodes and assign displacement difference
   for (unsigned j=0;j<nnode;j++)
    {
     // Get nodes
     SolidNode* orig_node_pt=orig_mesh_pt->node_pt(j);
     Node* new_node_pt=new_node[orig_node_pt];
     
     // Assign displacement difference
     for (unsigned i=0;i<3;i++)
      {
       orig_node_pt->x(i)=orig_node_pt->xi(i)+new_node_pt->value(i);
      }
    }     
     
   // Now re-assign undeformed position
   orig_mesh_pt->set_lagrangian_nodal_coordinates();
   
   // Clean up -- mesh deletes nodes and elements
   delete mesh_pt();
  }
 
 /// Destructor (empty)
 ~LinearSmoothMeshProblem(){}

};


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//====================================================================
/// Auxiliary Problem to smooth a SolidMesh by adjusting the internal
/// nodal positions by solving a solid mechanics problem for the
/// nodal displacements between the specified displacements of certain
/// pinned nodes (usually located on boundaries). Assumptions is
/// that the Lagrangian coordinates in the SolidMesh still reflect
/// the original nodal positions before the boundary nodes were
/// moved. 
/// hierher turn into functor?
/// hierher specify faster Poisson solver?
//====================================================================
class SmoothMeshProblem : public Problem
{

public:

 /// \short Constructor: Specify SolidMesh whose nodal positions are to 
 /// be adjusted, and set of nodes in that mesh whose position
 /// are to remain fixed.
 SmoothMeshProblem(SolidMesh* orig_mesh_pt,
                   std::set<Node*> pinned_nodes)
  {
   // Create new mesh and read out node/element numbers from old one
   mesh_pt()=new SolidMesh;
   unsigned nelem=orig_mesh_pt->nelement();
   unsigned nnode=orig_mesh_pt->nnode();
   
   // Have we already created that node?
   std::map<SolidNode*,SolidNode*> new_node;
   
   // Create new elements
   for (unsigned e=0;e<nelem;e++)
    {
     
     // Make/add new element
     TPVDElement<3,3>* el_pt=new TPVDElement<3,3>;
     mesh_pt()->add_element_pt(el_pt);
     
     // Set the constitutive law   
     el_pt->constitutive_law_pt() =
      hierher_namespace::Constitutive_law_pt;
     
     // Find corresponding original element
     SolidFiniteElement* orig_elem_pt=
      dynamic_cast<SolidFiniteElement*>(orig_mesh_pt->finite_element_pt(e));
     unsigned nnod=orig_elem_pt->nnode();
     
     // Create nodes
     for (unsigned j=0;j<nnod;j++)
      {
       // Does it not exist yet?
       if (new_node[dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))]==0)
        {
         SolidNode* new_nod_pt=
          dynamic_cast<SolidNode*>(mesh_pt()->finite_element_pt(e)->
                                   construct_node(j));
         new_node[dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))]=
          new_nod_pt;
         mesh_pt()->add_node_pt(new_nod_pt);
         for (unsigned i=0;i<3;i++)
          {
           // Set new nodal position to be the old one in the
           // SolidMesh (assumed to contain no inverted elements)
           new_nod_pt->x(i)=
            dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))->xi(i);
          }         
        }
       // It already exists -- copy across
       else
        {
         mesh_pt()->finite_element_pt(e)->node_pt(j)=
          new_node[dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))];
        }
      }
    }
   

   mesh_pt()->output("old_mesh.dat");

   // Now assign undeformed position
   dynamic_cast<SolidMesh*>(mesh_pt())->set_lagrangian_nodal_coordinates();
      
   // Loop over pinned nodes -- pin their positions and assign updated nodal 
   // positions
   double scale=0.001;
   for (std::set<Node*>::iterator it=pinned_nodes.begin();
        it!=pinned_nodes.end();it++)
    {
     new_node[dynamic_cast<SolidNode*>(*it)]->pin_position(0);
     new_node[dynamic_cast<SolidNode*>(*it)]->pin_position(1);
     new_node[dynamic_cast<SolidNode*>(*it)]->pin_position(2);

     new_node[dynamic_cast<SolidNode*>(*it)]->x(0)=
      dynamic_cast<SolidNode*>(*it)->xi(0)+scale*
      (dynamic_cast<SolidNode*>(*it)->x(0)-
       dynamic_cast<SolidNode*>(*it)->xi(0));

     new_node[dynamic_cast<SolidNode*>(*it)]->x(1)=
      dynamic_cast<SolidNode*>(*it)->xi(1)+scale*
      (dynamic_cast<SolidNode*>(*it)->x(1)-
       dynamic_cast<SolidNode*>(*it)->xi(1));

     new_node[dynamic_cast<SolidNode*>(*it)]->x(2)=
      dynamic_cast<SolidNode*>(*it)->xi(2)+scale*
      (dynamic_cast<SolidNode*>(*it)->x(2)-
       dynamic_cast<SolidNode*>(*it)->xi(2));
    }
   

   mesh_pt()->output("new_mesh.dat");


   oomph_info << "Number of equations for smoothing problem: " 
              << assign_eqn_numbers() << std::endl;
   
   
   // Solve
   double backup=Problem::Newton_solver_tolerance;
   //Newton_solver_tolerance=1.0e29;
   newton_solve();
   Problem::Newton_solver_tolerance=backup;
 

   // Loop over nodes and assign displacement difference
   for (unsigned j=0;j<nnode;j++)
    {
     // Get nodes
     SolidNode* orig_node_pt=orig_mesh_pt->node_pt(j);
     SolidNode* new_node_pt=new_node[orig_node_pt];
     
     // Assign displacement difference
     for (unsigned i=0;i<3;i++)
      {
       orig_node_pt->x(i)=new_node_pt->x(i);
      }
    }     
     
   // Now re-assign undeformed position
   orig_mesh_pt->set_lagrangian_nodal_coordinates();
   
   // Clean up -- mesh deletes nodes and elements
   delete mesh_pt();
  }
 
 /// Destructor (empty)
 ~SmoothMeshProblem(){}

};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//====================================================================
/// Functor to smooth a SolidMesh by adjusting the internal
/// nodal positions by solving a Poisson problem for the
/// nodal displacements in the interior. The displacements of the specified
/// pinned nodes (usually located on boundaries) remain fixed (their
/// displacements are computed from the difference between their
/// Lagrangian and Eulerian coordinates. The assumptions is
/// that the Lagrangian coordinates in the SolidMesh still reflect
/// the original nodal positions before the boundary nodes were
/// moved. 
/// hierher specify faster Poisson solver?
//====================================================================
template<class POISSON_ELEMENT>
class PoissonSmoothMesh : public Problem
{

public:

 /// \short Functor: Specify SolidMesh whose nodal positions are to 
 /// be adjusted, and set of nodes in that mesh whose position
 /// are to remain fixed.
 void operator()(SolidMesh* orig_mesh_pt,
                 std::set<Node*> pinned_nodes)  
  {
   // Create new mesh and read out node/element numbers from old one
   mesh_pt()=new Mesh;
   unsigned nelem=orig_mesh_pt->nelement();
   unsigned nnode=orig_mesh_pt->nnode();

   // Have we already created that node?
   std::map<Node*,Node*> new_node;

   // Create new elements
   for (unsigned e=0;e<nelem;e++)
    {
     mesh_pt()->add_element_pt(new POISSON_ELEMENT);

     // Find corresponding original element
     SolidFiniteElement* orig_elem_pt=
      dynamic_cast<SolidFiniteElement*>(orig_mesh_pt->finite_element_pt(e));
     unsigned nnod=orig_elem_pt->nnode();

     // Create nodes
     for (unsigned j=0;j<nnod;j++)
      {
       // Does it not exist yet?
       if (new_node[orig_elem_pt->node_pt(j)]==0)
        {
         Node* new_nod_pt=mesh_pt()->finite_element_pt(e)->construct_node(j);
         new_node[orig_elem_pt->node_pt(j)]=new_nod_pt;
         mesh_pt()->add_node_pt(new_nod_pt);
         for (unsigned i=0;i<3;i++)
          {
           // Set new nodal position to be the old one in the
           // SolidMesh (assumed to contain no inverted elements)
           new_nod_pt->x(i)=
            dynamic_cast<SolidNode*>(orig_elem_pt->node_pt(j))->xi(i);
          }         
        }
       // It already exists -- copy across
       else
        {
         mesh_pt()->finite_element_pt(e)->node_pt(j)=
          new_node[orig_elem_pt->node_pt(j)];
        }
      }
    }
   

   // Loop over pinned nodes
   for (std::set<Node*>::iterator it=pinned_nodes.begin();
        it!=pinned_nodes.end();it++)
    {
     new_node[*it]->pin(0);
    }

   oomph_info << "Number of equations for Poisson displacement smoothing: " 
              << assign_eqn_numbers() << std::endl;

   // Solve separate displacement problems
   for (unsigned i=0;i<3;i++)
    {
     // Loop over nodes and assign displacement difference
     for (unsigned j=0;j<nnode;j++)
      {
       // Get nodes
       SolidNode* orig_node_pt=orig_mesh_pt->node_pt(j);
       Node* new_node_pt=new_node[orig_node_pt];
       
       // Assign displacement difference
       new_node_pt->set_value(0,orig_node_pt->x(i)-orig_node_pt->xi(i));
      }

     // Solve
     newton_solve();
     
     // Loop over nodes and assign displacement difference
     for (unsigned j=0;j<nnode;j++)
      {
       // Get nodes
       SolidNode* orig_node_pt=orig_mesh_pt->node_pt(j);
       Node* new_node_pt=new_node[orig_node_pt];
       
       // Assign displacement difference 
       orig_node_pt->x(i)=orig_node_pt->xi(i)+new_node_pt->value(0);
      }     
    }

   // Now re-assign undeformed position
   orig_mesh_pt->set_lagrangian_nodal_coordinates();
    
   // Clean up -- mesh deletes nodes and elements
   delete mesh_pt();
   
  }
 
};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//==========start_solid_mesh===============================================
/// Tetgen-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class SolidTetMesh : public virtual TetgenMesh<ELEMENT>, 
                     public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 SolidTetMesh(const std::string& node_file_name,
                const std::string& element_file_name,
                const std::string& face_file_name,
                TimeStepper* time_stepper_pt=
                &Mesh::Default_TimeStepper) : 
  TetgenMesh<ELEMENT>(node_file_name, element_file_name,
                      face_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();

   // Find elements next to boundaries
   setup_boundary_element_info();

   // Setup boundary coordinates for all boundaries
   char filename[100];
   ofstream some_file;
   unsigned nb=this->nboundary();
   for (unsigned b=0;b<nb;b++)
    {
     sprintf(filename,"RESLT/solid_boundary_test%i.dat",b);
     some_file.open(filename);
     this->setup_boundary_coordinates(b,some_file);
     some_file.close();
    }

  }

 /// Scale all nodal coordinates by given factor
 void scale_mesh(const double& factor)
  {
   unsigned nnod=this->nnode();
   unsigned dim=this->node_pt(0)->ndim();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=this->node_pt(j);
     for (unsigned i=0;i<dim;i++)
      {
       nod_pt->x(i)*=factor;
      }
    }

   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~SolidTetMesh() { }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//==============start_fluid_mesh===========================================
/// Tetgen-based mesh upgraded to become a (pseudo-) solid mesh
//=========================================================================
template<class ELEMENT>
class FluidTetMesh : public virtual TetgenMesh<ELEMENT>,
                     public virtual SolidMesh 
{
 
public:
 
 /// \short Constructor: 
 FluidTetMesh(const std::string& node_file_name,
              const std::string& element_file_name,
              const std::string& face_file_name,
              const bool& split_corner_elements,
              TimeStepper* time_stepper_pt=
              &Mesh::Default_TimeStepper) : 
  TetgenMesh<ELEMENT>(node_file_name, element_file_name,
                      face_file_name, split_corner_elements, 
                      time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
   
   // Find out elements next to boundary
   setup_boundary_element_info();

   // Setup boundary coordinates for boundary.
   // To be consistent with the boundary coordinates generated
   // in the solid, we switch the direction of the normal.
   // (Both meshes are generated from the same polygonal facets
   // at the FSI interface).
   bool switch_normal=true;

   // Setup boundary coordinates for all boundaries
   char filename[100];
   ofstream some_file;
   unsigned nb=this->nboundary();
   for (unsigned b=0;b<nb;b++) 
    {
     sprintf(filename,"RESLT/fluid_boundary_test%i.dat",b);
     some_file.open(filename);
     this->setup_boundary_coordinates(b,switch_normal,some_file);
     some_file.close();
    }
  }


 /// Scale all nodal coordinates by given factor
 void scale_mesh(const double& factor)
  {
   unsigned nnod=this->nnode();
   unsigned dim=this->node_pt(0)->ndim();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=this->node_pt(j);
     for (unsigned i=0;i<dim;i++)
      {
       nod_pt->x(i)*=factor;
      }
    }

   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();

  }

 /// Empty Destructor
 virtual ~FluidTetMesh() { }

};
 

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=======start_of_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Default Reynolds number
 double Re=10.0; 

 /// Default FSI parameter
 double Q=0.0;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Poisson's ratio for generalised Hookean constitutive equation
 double Nu=0.3;
 
 /// Timescale ratio (non-dim density) for solid
 double Lambda_sq=0.0;

 /// Boolean to choose if inflow pressure (true) or volume flux are prescribed
 bool Prescribe_pressure=false;

 /// Max. inflow pressure (during periodic variation)
 double P_in_max=0.0025;

 /// Period of periodic variation in inflow pressure
 double Period=0.1;

 /// Fluid pressure on inflow boundary
 double P_in=0.0; //0.25;

 /// Fluid pressure on outflow boundary
 double P_out=-0.25; 

 /// Peak prescribed flow rate
 double Peak_prescribed_flow_rate=-1.0; 

 /// Prescribed flow rate
 double Prescribed_flow_rate=Peak_prescribed_flow_rate;

 /// \short Pointer to Data object whose one and only value stores the
 /// upstream pressure used to enforce flow rate
 Data* P_in_data_pt=0; 

 /// \short IDs for the two types of Lagrange multipliers used
 /// in this problem
 enum{Parallel_flow_lagrange_multiplier_id, 
      FSI_interface_displacement_lagrange_multiplier_id};
 
} //end_of_namespace






//===============start_of_problem_class===============================
/// Unstructured 3D FSI problem
//====================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
class UnstructuredFSIProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredFSIProblem();

 /// Destructor (empty)
 ~UnstructuredFSIProblem(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Create Lagrange multiplier elements that impose parallel flow
 void create_parallel_flow_lagrange_elements();

 /// Create FSI traction elements
 void create_fsi_traction_elements();

 /// \short Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 /// Actions before timestep 
 void actions_before_implicit_timestep()
  {
   // Get the current time
   double time=time_pt()->time();

   // Change inflow conditions
   if (Global_Parameters::Prescribe_pressure)
    {
     Global_Parameters::P_in=Global_Parameters::P_in_max*
      cos(2.0*MathematicalConstants::Pi/Global_Parameters::Period*time);
    }
   else
    {
     // hierher

     Global_Parameters::Prescribed_flow_rate=
      Global_Parameters::Peak_prescribed_flow_rate*
      cos(2.0*MathematicalConstants::Pi/Global_Parameters::Period*time);  
     cout << "New imposed prescribed flow rate: " 
          << Global_Parameters::Prescribed_flow_rate << std::endl;
    }
  }


 /// Update no slip bc after update of unknowns
 void actions_before_newton_convergence_check()
  {
   // hierher Fluid_mesh_pt->node_update();


   // The velocity of the fluid nodes on the fsi boundary
   // is set by the wall motion -- hence the no-slip condition must be
   // re-applied whenever a node update is performed for these nodes. 
   // Such tasks may be performed automatically by the auxiliary node update 
   // function specified by a function pointer:
   unsigned n=nfluid_fsi_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=Fluid_fsi_boundary_id[i];
     
     // How many nodes on boundary b?
     unsigned num_nod =  Fluid_mesh_pt->nboundary_node(b);     
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       FSI_functions::apply_no_slip_on_moving_wall(
        Fluid_mesh_pt->boundary_node_pt(b,inod));
      }
    }
   

  }

 /// Doc parameters
 void doc_parameters()
  {
   std::cout << "\n\n=================================================\n";
   std::cout << "Average wall thickness: " << Wall_thickness 
             << std::endl << std::endl;
   double speed_index=
    sqrt(Global_Parameters::Q/(0.5*Wall_thickness*Global_Parameters::Re));
   std::cout << "Q                     : " << Global_Parameters::Q << std::endl;
   std::cout << "Re                    : " << Global_Parameters::Re 
             << std::endl;
   std::cout << "Period                : " << Global_Parameters::Period 
             <<std::endl;
   std::cout << "Speed index           : " << speed_index  << std::endl;
   std::cout << "Wavelength            : " 
             << Global_Parameters::Period/speed_index << std::endl;
   double alpha_sq=Global_Parameters::Re/Global_Parameters::Period;
   std::cout << "Square of Womersley   : " 
             << alpha_sq << std::endl;
   std::cout << "=================================================\n\n";
  }
 
 
private:
 
 /// Sanity check: Doc boundary coordinates on i-th solid FSI interface
 void doc_solid_boundary_coordinates(const unsigned& i);
 
 /// \short Return total number of mesh boundaries that make up the inflow 
 /// boundary
 unsigned nfluid_inflow_traction_boundary()
  {return Inflow_boundary_id.size();}

 ///  \short Return total number of mesh boundaries that make up the outflow 
 /// boundary
 unsigned nfluid_outflow_traction_boundary()
  {return Outflow_boundary_id.size();}

 /// \short Return total number of mesh boundaries that make up the 
 /// in- and outflow boundaries where a traction has to be applied
 unsigned nfluid_traction_boundary()
  {return Inflow_boundary_id.size()+Outflow_boundary_id.size();}

 /// \short Return total number of mesh boundaries in the solid mesh that
 /// make up the FSI interface
 unsigned nsolid_fsi_boundary()
  {return Solid_fsi_boundary_id.size();}

 /// \short Return total number of mesh boundaries in the fluid mesh that
 /// make up the FSI interface
 unsigned nfluid_fsi_boundary()
  {return Fluid_fsi_boundary_id.size();}

 /// \short Return total number of mesh boundaries in the solid mesh 
 /// where the position is pinned.
 unsigned npinned_solid_boundary()
  {return Pinned_solid_boundary_id.size();} 
  //end npinned_solid_boundary

 /// Bulk solid mesh
 SolidTetMesh<SOLID_ELEMENT>* Solid_mesh_pt;

 /// Meshes of FSI traction elements
 Vector<SolidMesh*> Solid_fsi_traction_mesh_pt;

 /// Bulk fluid mesh
 FluidTetMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

 /// Mesh containing the FaceElements that monitor the inflow
 Mesh* Inflow_flux_monitor_mesh_pt;

 /// Mesh of FaceElements that control the inflow
 Mesh* Inflow_flux_control_mesh_pt;

 /// Mesh containing the master flux control element
 Mesh* Inflow_flux_control_master_mesh_pt;

 /// Meshes of Lagrange multiplier elements that impose parallel flow
 Vector<Mesh*> Parallel_flow_lagrange_multiplier_mesh_pt;

 /// Meshes of Lagrange multiplier elements
 Vector<SolidMesh*> Lagrange_multiplier_mesh_pt;

 /// GeomObject incarnations of the FSI boundary in the solid mesh
 Vector<MeshAsGeomObject<2,3,FSISolidTractionElement<SOLID_ELEMENT,3> >*>
 Solid_fsi_boundary_pt;

 /// IDs of solid mesh boundaries where displacements are pinned
 Vector<unsigned> Pinned_solid_boundary_id;
  
 /// \short IDs of solid mesh boundaries which make up the FSI interface
 Vector<unsigned> Solid_fsi_boundary_id;

 /// \short IDs of solid mesh boundaries which make up the outer surface
 Vector<unsigned> Solid_outer_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Outflow_boundary_id;

 /// \short IDs of fluid mesh boundaries which make up the FSI interface
 Vector<unsigned> Fluid_fsi_boundary_id;
 
 /// \short Average wall thickness (computed from area of fsi interface and
 /// volume of solid mesh.
 double Wall_thickness;
 
};



//==========start_of_constructor==========================================
/// Constructor for unstructured 3D FSI problem
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::UnstructuredFSIProblem()
{ 

 // Allocate the timestepper for the Navier-Stokes equations
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;

 // Add the fluid timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(fluid_time_stepper_pt);
 
 // Create a Steady timestepper that stores two history values
 Steady<2>* wall_time_stepper_pt = new Steady<2>; 

 // Add the wall timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(wall_time_stepper_pt);

 
 // Read boundary enumeration
 string input_string;
 
 // Open input file
 ifstream* input_file_pt=new ifstream("boundary_enumeration.dat"); 
 
 // Check if it's been opened succesfully
 if (input_file_pt==0)
  {
   throw OomphLibError(
   	"ERROR while trying to open boundary enumeration file ",
        "UnstructuredFSIProblem::UnstructuredFSIProblem()",
	 OOMPH_EXCEPTION_LOCATION);
  }
   
 std::ifstream boundary_enumeration("boundary_enumeration.dat");

 getline(*input_file_pt,input_string,' ');
 unsigned fluid_fsi_lo=atoi(input_string.c_str());
 getline(*input_file_pt,input_string,'#');
 unsigned fluid_fsi_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');

 getline(*input_file_pt,input_string,' ');
 unsigned fluid_in_lo=atoi(input_string.c_str());
 getline(*input_file_pt,input_string,'#');
 unsigned fluid_in_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');

 getline(*input_file_pt,input_string,' ');
 unsigned fluid_out_lo=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');
 getline(*input_file_pt,input_string,' ');
 getline(*input_file_pt,input_string,'#');
 unsigned fluid_out_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');
 
 getline(*input_file_pt,input_string,' ');
 unsigned solid_fsi_lo=atoi(input_string.c_str());
 getline(*input_file_pt,input_string,'#');
 unsigned solid_fsi_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');

 getline(*input_file_pt,input_string,' ');
 unsigned solid_pin_lo=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');
 input_file_pt->ignore(80,'\n');
 getline(*input_file_pt,input_string,' ');
 getline(*input_file_pt,input_string,'#');
 unsigned solid_pin_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');

 getline(*input_file_pt,input_string,' ');
 unsigned solid_outer_lo=atoi(input_string.c_str());
 getline(*input_file_pt,input_string,'#');
 unsigned solid_outer_hi=atoi(input_string.c_str());
 input_file_pt->ignore(80,'\n');


 boundary_enumeration.close();

 // Define fluid mesh and its distinguished boundaries
 //---------------------------------------------------
 
 //Create fluid bulk mesh, sub-dividing "corner" elements
 string node_file_name="fluid.1.node";
 string element_file_name="fluid.1.ele";
 string face_file_name="fluid.1.face";
 bool split_corner_elements=true;
 Fluid_mesh_pt =  new FluidTetMesh<FLUID_ELEMENT>(node_file_name,
                                                  element_file_name,
                                                  face_file_name,
                                                  split_corner_elements,
                                                  fluid_time_stepper_pt);
 

 // The following corresponds to the boundaries as specified by
 // facets in the xda input:

 // Fluid mesh inflow boundaries
 unsigned nfluid_in=fluid_in_hi-fluid_in_lo+1;
 Inflow_boundary_id.resize(nfluid_in);
 for(unsigned i=0; i<nfluid_in; i++)
  {
   Inflow_boundary_id[i]=fluid_fsi_hi+1+i;
  }
 
 // Fluid mesh outflow boundaries
 unsigned nfluid_out=fluid_out_hi-fluid_out_lo+1;
 Outflow_boundary_id.resize(nfluid_out);
 for(unsigned i=0; i<nfluid_out; i++)
  {
   Outflow_boundary_id[i]=fluid_out_lo+i;
  }

 // The FSI boundaries :
 unsigned nfluid_fsi=fluid_fsi_hi-fluid_fsi_lo+1;
 Fluid_fsi_boundary_id.resize(nfluid_fsi);
 for(unsigned i=0; i<nfluid_fsi; i++)
  {
   Fluid_fsi_boundary_id[i]=fluid_fsi_lo+i;
  }


 // Define solid mesh and its distinguished boundaries
 //---------------------------------------------------
   //Create solid bulk mesh
 node_file_name="solid.1.node";
 element_file_name="solid.1.ele";
 face_file_name="solid.1.face";
 Solid_mesh_pt =  new SolidTetMesh<SOLID_ELEMENT>(node_file_name,
                                                  element_file_name,
                                                  face_file_name,
                                                  wall_time_stepper_pt);
 
 // The following corresponds to the boundaries as specified by
 // facets in the Tetgen input:
 
 /// IDs of solid mesh boundaries where displacements are pinned
 unsigned nsolid_pin=solid_pin_hi-solid_pin_lo+1;
 Pinned_solid_boundary_id.resize(nsolid_pin);
 for(unsigned i=0; i<nsolid_pin; i++)
  {
   Pinned_solid_boundary_id[i]=solid_pin_lo+i;
  }

  // The solid and fluid fsi boundaries are numbered in the same way.
 unsigned nsolid_fsi=solid_fsi_hi-solid_fsi_lo+1;
 Solid_fsi_boundary_id.resize(nsolid_fsi);
 for(unsigned i=0; i<nsolid_fsi; i++)
  {
   Solid_fsi_boundary_id[i]=solid_fsi_lo+i;
  }

  // Outer solid boundary
 unsigned nsolid_outer=solid_outer_hi-solid_outer_lo+1;
 Solid_outer_boundary_id.resize(nsolid_outer);
 for(unsigned i=0; i<nsolid_outer; i++)
  {
   Solid_outer_boundary_id[i]=solid_outer_lo+i;
  }


 // Map to quadratic?
 bool fluid_quadr=true;
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   if (arg=="coarse_test")
    {
     oomph_info << "Not doing adjustment to quadratic geometry\n";
     fluid_quadr=false;
    }
  }

 if (fluid_quadr)
  {

   // Snap curved fluid FSI boundaries onto quadratic surface
   //--------------------------------------------------------

   
   // Check for inverted elements before quadratic trafo
   bool mesh_has_inverted_elements;
   std::ofstream inverted_fluid_elements;
   Fluid_mesh_pt->output("fluid_mesh_before_snap.dat");
   inverted_fluid_elements.open(
    "inverted_fluid_elements_before_snap.dat");
   Fluid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_fluid_elements);  
   inverted_fluid_elements.close(); 
   cout << "Before quadratic snapping fluid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";
   

   // Snap fluid FSI interface to quadratic surface
   DocInfo doc_info;
   doc_info.set_directory("RESLT");
   doc_info.label()="fluid_fsi";
   bool switch_normal=false;
   Fluid_mesh_pt->snap_to_quadratic_surface(Fluid_fsi_boundary_id,
                                            "quadratic_fsi_boundary.dat",
                                            switch_normal,
                                            doc_info);

   // Doc
   Fluid_mesh_pt->output("fluid_mesh_after_snap.dat");
   inverted_fluid_elements.open(
    "inverted_fluid_elements_after_snap.dat");
   Fluid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_fluid_elements);  
   inverted_fluid_elements.close(); 
   cout << "After quadratic snapping fluid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";
   
   
   // Now smooth nodal displacements by Poisson solves
   
   // Create set containing the nodes whose displacements are contrained
   // (nodes on FSI boundary)
   std::set<Node*> pinned_nodes;
   unsigned nbound_pinned=nfluid_fsi_boundary();
   for(unsigned i=0;i<nbound_pinned;i++)
    {
     //Get the mesh boundary
     unsigned b = Fluid_fsi_boundary_id[i];
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       pinned_nodes.insert(Fluid_mesh_pt->boundary_node_pt(b,inod));
      }
    }
      
   // Smooth
   bool use_poisson=false;
   if (use_poisson)
    {
     PoissonSmoothMesh<TPoissonElement<3,3> >()(Fluid_mesh_pt,pinned_nodes);
    }
   else
    {
     LinearSmoothMeshProblem(Fluid_mesh_pt,pinned_nodes);
    }

   // Doc
   Fluid_mesh_pt->output("fluid_mesh_after_poisson_smooth.dat");
   inverted_fluid_elements.open(
    "inverted_fluid_elements_after_poisson_smooth.dat");
   Fluid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_fluid_elements);   
   inverted_fluid_elements.close();
   cout << "After Poisson smoothing fluid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";
   
  }




 // Map to quadratic?
 bool solid_quadr=true;
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   if (arg=="coarse_test")
    {
     oomph_info << "Not doing adjustment to quadratic geometry\n";
     solid_quadr=false;
    }
  }
 if (solid_quadr)
  {
   
   
   // Snap curved solid boundaries (fsi & outer boundary) onto quadratic surface
   //---------------------------------------------------------------------------
   
   // Check for inverted elements before quadratic trafo
   bool mesh_has_inverted_elements;
   std::ofstream inverted_solid_elements;
   Solid_mesh_pt->output("solid_mesh_before_snap.dat");
   inverted_solid_elements.open(
    "inverted_solid_elements_before_snap.dat");
   Solid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_solid_elements);  
   inverted_solid_elements.close(); 
   cout << "Before quadratic snapping solid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";
   
   
   // Snap solid FSI interface to quadratic surface
   DocInfo doc_info;
   doc_info.set_directory("RESLT");
   doc_info.label()="solid_fsi";
   bool switch_normal=false;
   Solid_mesh_pt->snap_to_quadratic_surface(Solid_fsi_boundary_id,
                                            "quadratic_fsi_boundary.dat",
                                            switch_normal,
                                            doc_info);
   
   // Snap outer solid interface to quadratic surface
   doc_info.set_directory("RESLT");
   doc_info.label()="solid_outer_boundary";
   switch_normal=true;
   Solid_mesh_pt->snap_to_quadratic_surface(
    Solid_outer_boundary_id,
    "quadratic_outer_solid_boundary.dat",
    switch_normal,
    doc_info);
   

   // Doc outcome of snapping
   Solid_mesh_pt->output("solid_mesh_after_snap.dat");
   inverted_solid_elements.open(
    "inverted_solid_elements_after_snap.dat");
   Solid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_solid_elements);  
   inverted_solid_elements.close(); 
   cout << "After quadratic snapping solid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";
   
   // Loop over nodes on the FSI boundary in the solid mesh
   std::set<Node*> pinned_nodes;
   unsigned nbound_pinned=nsolid_fsi_boundary();
   for(unsigned i=0;i<nbound_pinned;i++)
    {
     //Get the mesh boundary
     unsigned b = Solid_fsi_boundary_id[i];
     unsigned num_nod=Solid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       pinned_nodes.insert(Solid_mesh_pt->boundary_node_pt(b,inod));
      }
    }
   
   // Loop over nodes on outer boundary of solid mesh
   nbound_pinned=Solid_outer_boundary_id.size();
   for(unsigned i=0;i<nbound_pinned;i++)
    {
     //Get the mesh boundary
     unsigned b = Solid_outer_boundary_id[i];
     unsigned num_nod=Solid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node     
       pinned_nodes.insert(Solid_mesh_pt->boundary_node_pt(b,inod));
      }
    }
   

   bool use_poisson=false;
   if (use_poisson)
    {
     PoissonSmoothMesh<TPoissonElement<3,3> >()(Solid_mesh_pt,pinned_nodes);
    }
   else
    {
     LinearSmoothMeshProblem(Solid_mesh_pt,pinned_nodes);
    }

   
   // Doc
   Solid_mesh_pt->output("solid_mesh_after_poisson_smooth.dat");
   inverted_solid_elements.open(
    "inverted_solid_elements_after_poisson_smooth.dat");
   Solid_mesh_pt->check_inverted_elements(mesh_has_inverted_elements,
                                          inverted_solid_elements);   
   inverted_solid_elements.close();
   cout << "After Poisson smoothing solid mesh does ";
   if (!mesh_has_inverted_elements) cout << "not ";
   cout << "have inverted elements. \n";

   
  }

 // Scale 
 double factor=0.0666;
 Fluid_mesh_pt->scale_mesh(factor);
 Solid_mesh_pt->scale_mesh(factor);


 // Create (empty) meshes of lagrange multiplier elements at inflow/outflow
 //------------------------------------------------------------------------
 
 // Create the meshes
 unsigned n=nfluid_traction_boundary();
 Parallel_flow_lagrange_multiplier_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Parallel_flow_lagrange_multiplier_mesh_pt[i]=new Mesh;
  } 
 

 // Meshes of inflow elements
 Inflow_flux_monitor_mesh_pt=new Mesh;
 Inflow_flux_control_mesh_pt=new Mesh;
 Inflow_flux_control_master_mesh_pt=new Mesh;

 // Populate them with elements
 create_parallel_flow_lagrange_elements();

 cout << "Number of inflow flux monitor elements: " 
      << Inflow_flux_monitor_mesh_pt->nelement() << std::endl;

// Create FSI Traction elements
//-----------------------------
 
// Create (empty) meshes of FSI traction elements
 n=nsolid_fsi_boundary();
 Solid_fsi_traction_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Solid_fsi_traction_mesh_pt[i]=new SolidMesh;
  }
 
 // Build the FSI traction elements
 create_fsi_traction_elements();


 // Compute average wall thickness
 //-------------------------------
 
 // Surface area of fsi interface
 double fsi_area=0.0;
 n=nsolid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   unsigned nel=Solid_fsi_traction_mesh_pt[i]->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     fsi_area+=Solid_fsi_traction_mesh_pt[i]->finite_element_pt(e)->size();
    }
  }
 
 // Volume of solid mesh
 double solid_vol=0.0;
 unsigned nel=Solid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   solid_vol+=Solid_mesh_pt->finite_element_pt(e)->size();
  }
 Wall_thickness=solid_vol/fsi_area;  
 std::cout << "Surface area          : " << fsi_area << std::endl;
 std::cout << "Solid volume          : " << solid_vol << std::endl;
 
 
 // Create Lagrange multiplier mesh for boundary motion of fluid mesh
 //------------------------------------------------------------------
 
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 n=nfluid_fsi_boundary();
 Lagrange_multiplier_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Lagrange_multiplier_mesh_pt[i]=new SolidMesh;
  }
 
 // Create elements
 create_lagrange_multiplier_elements();
 
 

 // Combine the lot
 //----------------
 
 // Add sub meshes:

 // The solid bulk mesh
 add_sub_mesh(Solid_mesh_pt);

 // Fluid bulk mesh
 add_sub_mesh(Fluid_mesh_pt);
 
 // The fluid traction meshes
 n=nfluid_traction_boundary();
 for (unsigned i=0;i<n;i++)
  { 
   add_sub_mesh(Parallel_flow_lagrange_multiplier_mesh_pt[i]);
  }
 
 // The solid fsi traction meshes
 n=nsolid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   add_sub_mesh(Solid_fsi_traction_mesh_pt[i]);
  }
 
 // The Lagrange multiplier meshes for the fluid
 n=nfluid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {   
   add_sub_mesh(Lagrange_multiplier_mesh_pt[i]);
  }

 // Add flux control meshes if required
 if (!Global_Parameters::Prescribe_pressure)
  {
   add_sub_mesh(Inflow_flux_control_mesh_pt);
   add_sub_mesh(Inflow_flux_control_master_mesh_pt);
  }

 // Build global mesh
 build_global_mesh();



 
 // Apply BCs for fluid and Lagrange multiplier elements
 //-----------------------------------------------------
  
 // Doc position of pinned pseudo solid nodes
 std::ofstream pseudo_solid_bc_file("RESLT/pinned_pseudo_solid_nodes.dat");
 
 // Loop over inflow/outflow boundaries to pin pseudo-solid displacements
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over in/outflow boundaries
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {

     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }

     // Number of nodes on that boundary
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get the node
       SolidNode* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,inod);
       
       // Pin the nodal (pseudo-solid) displacements
       for(unsigned i=0;i<3;i++)
        {         
         nod_pt->pin_position(i);         
         
         // Doc it as pinned
         pseudo_solid_bc_file << nod_pt->x(i) << " ";
        }
      }
    }
  }
 
 // Close
 pseudo_solid_bc_file.close();
 
// Doc bcs for Lagrange multipliers
 ofstream pinned_file("RESLT/pinned_lagrange_multiplier_nodes.dat");
 
 // Loop over nodes on the FSI boundary in the fluid mesh
 unsigned nbound=nfluid_fsi_boundary();
 for(unsigned i=0;i<nbound;i++)
  {
   //Get the mesh boundary
   unsigned b = Fluid_fsi_boundary_id[i];
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt= Fluid_mesh_pt->boundary_node_pt(b,inod);
     
     // Pin all velocities
     nod_pt->pin(0); 
     nod_pt->pin(1); 
     nod_pt->pin(2); 
     
     // Find out whether node is also on in/outflow
     bool is_in_or_outflow_node=false;
     unsigned n=nfluid_inflow_traction_boundary();
     for (unsigned k=0;k<n;k++)
      {
       if (nod_pt->is_on_boundary(Inflow_boundary_id[k]))
        {
         is_in_or_outflow_node=true;
         break;
        }
      }
     if (!is_in_or_outflow_node)
      {
       unsigned n=nfluid_outflow_traction_boundary();
       for (unsigned k=0;k<n;k++)
        {
         if (nod_pt->is_on_boundary(Outflow_boundary_id[k]))
          {
           is_in_or_outflow_node=true;
           break;
          }
        }
      } // ...now we know if the node is on an in- or outflow boundary
     
     // Pin the Lagrange multipliers for the imposition of
     // parallel flow if the nodes is also on the in/outflow boundaries
     if(is_in_or_outflow_node)
      {
       //Cast to a boundary node
       BoundaryNode<SolidNode> *bnod_pt = 
        dynamic_cast<BoundaryNode<SolidNode>*>
        ( Fluid_mesh_pt->boundary_node_pt(b,inod) );
       
       // Get the index of the first Lagrange multiplier
       unsigned first_index=bnod_pt->
        index_of_first_value_assigned_by_face_element(
         Global_Parameters::Parallel_flow_lagrange_multiplier_id);

       //Pin the Lagrange multipliers (as the velocity is already
       //determined via the no slip condition on the fsi boundary
       for (unsigned l=0;l<2;l++)
        {
         nod_pt->pin(first_index+l);
        }

       
       // Get the first index of the second Lagrange multiplier 
       first_index=bnod_pt->index_of_first_value_assigned_by_face_element(
        Global_Parameters::FSI_interface_displacement_lagrange_multiplier_id);

       // Loop over the Lagrange multipliers that deform the FSI boundary
       // of the pseudo-solid fluid mesh.
       for (unsigned l=0;l<3;l++)
        {
         // Pin the Lagrange multipliers that impose the displacement
         // because the positon of the fluid nodes at the in/outflow
         // is already determined. 
         nod_pt->pin(first_index+l);
        }

       // Doc that we've pinned the Lagrange multipliers at this node
       pinned_file << nod_pt->x(0) << " "
                   << nod_pt->x(1) << " "
                   << nod_pt->x(2) << endl;
      }
    }
  } // end of BC for fluid mesh

 // Done pinning Lagrange nultipliers
 pinned_file.close();
  
 
 // Complete the build of the fluid elements so they are fully functional
 //----------------------------------------------------------------------
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   FLUID_ELEMENT* el_pt = 
    dynamic_cast<FLUID_ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;
   el_pt->re_st_pt() = &Global_Parameters::Re; 
   
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pt;
   
  } // end loop over elements



 // Apply BCs for solid
 //--------------------
 
 // Doc pinned solid nodes
 std::ofstream bc_file("RESLT/pinned_solid_nodes.dat");
 
 // Pin positions at inflow boundary (boundaries 0 and 1)
 n=npinned_solid_boundary();
 for (unsigned i=0;i<n;i++)
  {
   // Get boundary ID
   unsigned b=Pinned_solid_boundary_id[i];
   unsigned num_nod= Solid_mesh_pt->nboundary_node(b);  
   for (unsigned inod=0;inod<num_nod;inod++)
    {    
     // Get node
     SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(b,inod);
     
     // Pin all directions
     for (unsigned i=0;i<3;i++)
      {
       nod_pt->pin_position(i);
       
       // ...and doc it as pinned
       bc_file << nod_pt->x(i) << " ";
      }
     
     bc_file << std::endl;
    }
  }
 bc_file.close();
 
 
 
 // Complete the build of Solid elements so they are fully functional
 //----------------------------------------------------------------
 n_element = Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   SOLID_ELEMENT *el_pt = dynamic_cast<SOLID_ELEMENT*>(
    Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law   
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pt;

   // Set the timescale ratio 
   el_pt->lambda_sq_pt() =
    &Global_Parameters::Lambda_sq;
  }


 // Setup FSI
 //----------
     
 // The velocity of the fluid nodes on the fsi boundary
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 n=nfluid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   // Get boundary ID
   unsigned b=Fluid_fsi_boundary_id[i];
   
   // How many nodes on boundary b?
   unsigned num_nod =  Fluid_mesh_pt->nboundary_node(b);
   
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Fluid_mesh_pt->boundary_node_pt(b,inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  }
 

 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. 
 n=nsolid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   // Sanity check: Doc boundary coordinates from solid side
   doc_solid_boundary_coordinates(i);
   
   //Doc boundary coordinates in fluid
   char filename[100];
   sprintf(filename,"RESLT/fluid_boundary_coordinates%i.dat",i);
   FSI_functions::Doc_boundary_coordinate_file.open(filename);
   
   // Setup FSI: Pass ID of fluid FSI boundary and associated
   // mesh of solid fsi traction elements.
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,3>
    (this,Fluid_fsi_boundary_id[i],Fluid_mesh_pt,Solid_fsi_traction_mesh_pt[i]);
   
   // Close the doc file
   FSI_functions::Doc_boundary_coordinate_file.close();
  } 
 
 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
}


//============start_of_create_fsi_traction_elements======================
/// Create FSI traction elements 
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_fsi_traction_elements()
{

 // Loop over FSI boundaries in solid
 unsigned n=nsolid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   // Get boundary ID
   unsigned b=Solid_fsi_boundary_id[i];
   
   // How many bulk elements are adjacent to boundary b?
   unsigned n_element = Solid_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
      Solid_mesh_pt->boundary_element_pt(b,e));
     
     //What is the index of the face of the element e along boundary b
     int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element 
     FSISolidTractionElement<SOLID_ELEMENT,3>* el_pt=
      new FSISolidTractionElement<SOLID_ELEMENT,3>(bulk_elem_pt,face_index);
     
     // Add it to the mesh
     Solid_fsi_traction_mesh_pt[i]->add_element_pt(el_pt);
     
     // Specify boundary number
     el_pt->set_boundary_number_in_bulk_mesh(b);
     
     // Function that specifies the load ratios
     el_pt->q_pt() = &Global_Parameters::Q; 
    }
  }
 
} // end of create_fsi_traction_elements


//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Make space
 unsigned n=nfluid_fsi_boundary();
 Solid_fsi_boundary_pt.resize(n);
    
 // Loop over FSI interfaces in fluid
 for (unsigned i=0;i<n;i++)
  {   
   // Get boundary ID
   unsigned b=Fluid_fsi_boundary_id[i];
   
   // Create  GeomObject incarnation of fsi boundary in solid mesh
   Solid_fsi_boundary_pt[i]=
    new MeshAsGeomObject<2,3,FSISolidTractionElement<SOLID_ELEMENT,3> >
    (Solid_fsi_traction_mesh_pt[i]);
   
   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is adjacent to boundary b
     FLUID_ELEMENT* bulk_elem_pt = dynamic_cast<FLUID_ELEMENT*>(
      Fluid_mesh_pt->boundary_element_pt(b,e));

     //Find the index of the face of element e along boundary b
     int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     ImposeDisplacementByLagrangeMultiplierElement<FLUID_ELEMENT>* el_pt =
      new  ImposeDisplacementByLagrangeMultiplierElement<FLUID_ELEMENT>(
       bulk_elem_pt,face_index, 
       Global_Parameters::FSI_interface_displacement_lagrange_multiplier_id);   
     
     // Add it to the mesh
     Lagrange_multiplier_mesh_pt[i]->add_element_pt(el_pt);
     
     // Set the GeomObject that defines the boundary shape and set
     // which bulk boundary we are attached to (needed to extract
     // the boundary coordinate from the bulk nodes)
     el_pt->set_boundary_shape_geom_object_pt(Solid_fsi_boundary_pt[i],b);
    }
  }

} // end of create_lagrange_multiplier_elements



//============start_of_create_parallel_flow_lagrange_elements============
/// Create Lagrange multiplier elements  that impose parallel flow
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_parallel_flow_lagrange_elements()
{
 // Counter for number of Lagrange multiplier elements
 unsigned count=0;

 // Loop over inflow/outflow boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over boundaries with fluid traction elements
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }
     
     // How many bulk elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       FLUID_ELEMENT* bulk_elem_pt = dynamic_cast<FLUID_ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //What is the index of the face of the element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Build the corresponding Lagrange multiplier element
       ImposeParallelOutflowElement<FLUID_ELEMENT>* el_pt 
        = new ImposeParallelOutflowElement<FLUID_ELEMENT>
        (bulk_elem_pt, face_index, 
         Global_Parameters::Parallel_flow_lagrange_multiplier_id);

       // Add it to the mesh
       Parallel_flow_lagrange_multiplier_mesh_pt[count]->add_element_pt(el_pt);
       
       // Set the pointer to the prescribed pressure
       if (in_out==0)
        {

         // Attach elements that monitor volume flux (etc).
         NavierStokesSurfacePowerElement<FLUID_ELEMENT>* power_element_pt = new 
          NavierStokesSurfacePowerElement<FLUID_ELEMENT>(
           bulk_elem_pt,face_index);
         
         //Add the elements to the mesh
         Inflow_flux_monitor_mesh_pt->add_element_pt(power_element_pt);     
         
         // Pressure or flux control
         if (Global_Parameters::Prescribe_pressure)
          {
           // Set inflow pressure
           el_pt->pressure_pt()= &Global_Parameters::P_in;
          }
         else
          {
           // Build the flux control element
           NavierStokesFluxControlElement<FLUID_ELEMENT>* flux_element_pt = new 
            NavierStokesFluxControlElement<FLUID_ELEMENT>(bulk_elem_pt, 
                                                          face_index);
           
           //Add the new element to its mesh
           Inflow_flux_control_mesh_pt->add_element_pt(flux_element_pt);  
          }
        }
       else
        {
         el_pt->pressure_pt()= &Global_Parameters::P_out;
        }
       //end of element setup

       
      } 
     // Bump up counter
     count++;
    }
  }

 // Build master flux control element
 if (!Global_Parameters::Prescribe_pressure)
  {
   // Build master element
   NetFluxControlElement* flux_control_el_pt = 
    new NetFluxControlElement(Inflow_flux_control_mesh_pt,
                              &Global_Parameters::Prescribed_flow_rate);   

   // Add NetFluxControlElement to its own mesh
   Inflow_flux_control_master_mesh_pt->add_element_pt(flux_control_el_pt);

   // Get pointer to the inflow pressure data
   Global_Parameters::P_in_data_pt = flux_control_el_pt->pressure_data_pt();
   
   // Loop over the elements in the sub mesh and add P_in_data_pt
   // to external data
   unsigned n_el = Inflow_flux_control_mesh_pt->nelement();
   for (unsigned e=0; e< n_el; e++)
    {
     // Get pointer to the element
     GeneralisedElement* el_pt = Inflow_flux_control_mesh_pt->element_pt(e);

     // Dynamic cast
     dynamic_cast<NavierStokesFluxControlElement<FLUID_ELEMENT>* >(el_pt)->   
      add_pressure_data(Global_Parameters::P_in_data_pt);
    }
   
  }
  
}  // end of create_parallel_flow_lagrange_elements



//============start_doc_solid_zeta=======================================
/// Doc boundary coordinates of i-th solid FSI boundary.
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
doc_solid_boundary_coordinates(const unsigned& i)
{
 
 //Doc boundary coordinates in fluid
 char filename[100];
 sprintf(filename,"RESLT/solid_boundary_coordinates%i.dat",i);
 std::ofstream the_file(filename);
 
 // Loop over traction elements
 unsigned n_face_element = Solid_fsi_traction_mesh_pt[i]->nelement();
 for(unsigned e=0;e<n_face_element;e++)
  {
   //Cast the element pointer
   FSISolidTractionElement<SOLID_ELEMENT,3>* el_pt=
    dynamic_cast<FSISolidTractionElement<SOLID_ELEMENT,3>*>
    (Solid_fsi_traction_mesh_pt[i]->element_pt(e));
   
   // Doc boundary coordinate
   Vector<double> s(2);
   Vector<double> zeta(2);
   Vector<double> x(3);
   unsigned n_plot=3;
   the_file << el_pt->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {         
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,n_plot,s);         
     el_pt->interpolated_zeta(s,zeta);
     el_pt->interpolated_x(s,x);
     for (unsigned i=0;i<3;i++)
      {
       the_file << x[i] << " ";
      }
     for (unsigned i=0;i<2;i++)
      {
       the_file << zeta[i] << " ";
      }

     the_file << std::endl;
    }
   el_pt->write_tecplot_zone_footer(the_file,n_plot);
  } 

 // Close doc file
 the_file.close();

} // end doc_solid_zeta


//========start_of_doc_solution===========================================
/// Doc the solution
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 


 // hierher Elements can still be inverted near boundary
 // even if they're non-inverted at the 3D Gauss points
 // thus code can die here when evaluating the fluid traction
 // Deal with this later...
 FiniteElement::Accept_negative_jacobian=true; 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
 
 // Output solid boundaries
 //------------------------
 sprintf(filename,"%s/solid_boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output_boundaries(some_file);
 some_file.close();
 
 
 // Output solid solution
 //-----------------------
 sprintf(filename,"%s/solid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 
 // Output fluid boundaries
 //------------------------
 sprintf(filename,"%s/fluid_boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output_boundaries(some_file);
 some_file.close();
 
 
 // Output fluid solution
 //-----------------------
 sprintf(filename,"%s/fluid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();
  
   
 // Output fsi traction
 //--------------------
 sprintf(filename,"%s/fsi_traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned n=nsolid_fsi_boundary();
 for (unsigned i=0;i<n;i++)
  {
   Solid_fsi_traction_mesh_pt[i]->output(some_file,npts);
  }
 some_file.close();



 // Document the volume flux at inflow
 double inflow_volume_flux = 0.0;
 unsigned nel = Inflow_flux_monitor_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   NavierStokesSurfacePowerElement<FLUID_ELEMENT>* el_pt=dynamic_cast<
    NavierStokesSurfacePowerElement<FLUID_ELEMENT>*>(
     Inflow_flux_monitor_mesh_pt->element_pt(e));
   inflow_volume_flux += el_pt->get_volume_flux();
  }
 cout << "Volume flux at inflow: " << inflow_volume_flux << std::endl;


 // Output nodal positions (only) if test of quadratic snapping
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   if (arg=="quadratic_test")
    {
     
     sprintf(filename,"%s/solid_nodes%i.dat",doc_info.directory().c_str(),
             doc_info.number());
     some_file.open(filename);
     unsigned nnod=Solid_mesh_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       some_file << Solid_mesh_pt->node_pt(j)->x(0) << " "
                 << Solid_mesh_pt->node_pt(j)->x(1) << " "
                 << Solid_mesh_pt->node_pt(j)->x(2) << std::endl;
      }
     some_file.close();

     sprintf(filename,"%s/fluid_nodes%i.dat",doc_info.directory().c_str(),
             doc_info.number());
     some_file.open(filename);
     nnod=Fluid_mesh_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       some_file << Fluid_mesh_pt->node_pt(j)->x(0) << " "
                 << Fluid_mesh_pt->node_pt(j)->x(1) << " "
                 << Fluid_mesh_pt->node_pt(j)->x(2) << std::endl;
      }
     some_file.close();
    }
  }



 // hierher Elements can still be inverted near boundary
 // even if they're non-inverted at the 3D Gauss points
 // thus code can die here when evaluating the fluid traction
 // Deal with this later...
 FiniteElement::Accept_negative_jacobian=false; 

} // end_of_doc





//========================= start_of_main=================================
/// Demonstrate how to solve an unstructured 3D FSI problem
//========================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // hierher
 TetMeshBase::Tolerance_for_boundary_finding = 1e-10;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 // Create generalised Hookean constitutive equations
 Global_Parameters::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Parameters::Nu);
 
 //Set up the problem
 UnstructuredFSIProblem<
 PseudoSolidNodeUpdateElement<TTaylorHoodElement<3>, TPVDElement<3,3> >,
  TPVDElement<3,3> > problem;

 //Output initial configuration
 problem.doc_solution(doc_info);
 doc_info.number()++;   

 // Stop if quadratic test
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   cout << "Arg: " << arg << std::endl;
   if (arg=="quadratic_test")
    {
     oomph_info << "Only testing the quadratic mapping; stopping now.\n";
     exit(0);
    }
  }

 //Set value of dt
 double nsteps_per_period=40;
 double dt = Global_Parameters::Period/double(nsteps_per_period);
 problem.initialise_dt(dt);

 // Increment FSI parameter during steady runs
 unsigned nq=4;
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   if (arg=="coarse_test")
    {
     oomph_info << "Smaller number of steps for test\n";
     nq=1;
    }
  }
 double q_increment=2.0e-6;   
 for (unsigned iq=0;iq<nq;iq++)
  {
   // Bump up fsi parameter
   Global_Parameters::Q+=q_increment;
   
   // Doc parameters
   problem.doc_parameters();
   
   cout << "Steady solve...\n";

   // Solve the steady problem first
   problem.steady_newton_solve();
   
   //Output solution
   problem.doc_solution(doc_info);
   doc_info.number()++;
  }

 // Do timestepping
 unsigned nperiod=20;
 unsigned nstep=nperiod*nsteps_per_period;
 if (CommandLineArgs::Argc==2)
  {
   string arg=CommandLineArgs::Argv[1];
   if (arg=="coarse_test")
    {
     oomph_info << "Smaller number of steps for test\n";
     nstep=2;
    }
  }
 for (unsigned istep=0;istep<nstep;istep++)
  {

   cout << "Unsteady solve...\n";

   // Solve the problem
   problem.unsteady_newton_solve(dt);
   
   //Output solution
   problem.doc_solution(doc_info);
   doc_info.number()++;
  }
 
} // end_of_main

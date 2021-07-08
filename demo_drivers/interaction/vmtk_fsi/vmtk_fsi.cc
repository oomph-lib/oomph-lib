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
// Driver code for a simple unstructured fsi problem using meshes
// generated with VMTK.


//Generic libraries
#include "generic.h"
#include "solid.h"
#include "constitutive.h"
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;
using namespace oomph;



//==========start_solid_mesh===============================================
/// Tetgen-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class MySolidTetgenMesh : public virtual TetgenMesh<ELEMENT>, 
                       public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 MySolidTetgenMesh(const std::string& node_file_name,
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
     this->template setup_boundary_coordinates<ELEMENT>(b,some_file);
     some_file.close();
    }

  }

 /// Empty Destructor
 virtual ~MySolidTetgenMesh() { }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//==============start_fluid_mesh===========================================
/// Tetgen-based mesh upgraded to become a (pseudo-) solid mesh
//=========================================================================
template<class ELEMENT>
class MyFluidTetMesh : public virtual TetgenMesh<ELEMENT>,
                       public virtual SolidMesh 
{
 
public:
 
 /// \short Constructor: 
 MyFluidTetMesh(const std::string& node_file_name,
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
     this->template setup_boundary_coordinates<ELEMENT>
      (b,switch_normal,some_file);
     some_file.close();
    }
  }

 /// Empty Destructor
 virtual ~MyFluidTetMesh() { }

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
 double Re=50.0; 

 /// Default FSI parameter
 double Q=0.0;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Poisson's ratio for generalised Hookean constitutive equation
 double Nu=0.3;
 
 /// Fluid pressure on inflow boundary
 double P_in=0.25;

 /// Fluid pressure on outflow boundary
 double P_out=-0.25; 


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
 MySolidTetgenMesh<SOLID_ELEMENT>* Solid_mesh_pt;

 /// Meshes of FSI traction elements
 Vector<SolidMesh*> Solid_fsi_traction_mesh_pt;

 /// Bulk fluid mesh
 MyFluidTetMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

 /// Meshes of Lagrange multiplier elements that impose parallel flow
 Vector<Mesh*> Parallel_flow_lagrange_multiplier_mesh_pt;

 /// Meshes of Lagrange multiplier elements
 Vector<SolidMesh*> Lagrange_multiplier_mesh_pt;

 /// GeomObject incarnations of the FSI boundary in the solid mesh
 Vector<MeshAsGeomObject*>
 Solid_fsi_boundary_pt;

 /// IDs of solid mesh boundaries where displacements are pinned
 Vector<unsigned> Pinned_solid_boundary_id;
  
 /// \short IDs of solid mesh boundaries which make up the FSI interface
 Vector<unsigned> Solid_fsi_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Outflow_boundary_id;

 /// \short IDs of fluid mesh boundaries which make up the FSI interface
 Vector<unsigned> Fluid_fsi_boundary_id;

};



//==========start_of_constructor==========================================
/// Constructor for unstructured 3D FSI problem
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::UnstructuredFSIProblem()
{ 

 // We have a large number of sub-meshes with very few elements
 // in this problem: Reduce the number of bins in the MeshAsGeomObject
 // representations of these meshes to avoid memory problems
 NonRefineableBinArray::Default_n_bin_1d=1;


 // Define fluid mesh and its distinguished boundaries
 //---------------------------------------------------
 
   //Create fluid bulk mesh, sub-dividing "corner" elements
 string node_file_name="fluid_iliac.1.node";
 string element_file_name="fluid_iliac.1.ele";
 string face_file_name="fluid_iliac.1.face";
 bool split_corner_elements=true;
 Fluid_mesh_pt =  new MyFluidTetMesh<FLUID_ELEMENT>(node_file_name,
                                                    element_file_name,
                                                    face_file_name,
                                                    split_corner_elements);
 

 // The following corresponds to the boundaries as specified by
 // facets in the xda input:

 // Fluid mesh inflow boundaries
 Inflow_boundary_id.resize(22);
 for(unsigned i=0; i<22; i++)
  {
   Inflow_boundary_id[i]=215+i;
  }
 
 // Fluid mesh outflow boundaries
 Outflow_boundary_id.resize(11);
 for(unsigned i=0; i<11; i++)
  {
   Outflow_boundary_id[i]=237+i;
  }

 // The FSI boundaries :
 Fluid_fsi_boundary_id.resize(215);
 for(unsigned i=0; i<215; i++)
  {
   Fluid_fsi_boundary_id[i]=i;
  }

 
 // Define solid mesh and its distinguished boundaries
 //---------------------------------------------------
   //Create solid bulk mesh
 node_file_name="solid_iliac.1.node";
 element_file_name="solid_iliac.1.ele";
 face_file_name="solid_iliac.1.face";
 Solid_mesh_pt =  new MySolidTetgenMesh<SOLID_ELEMENT>(node_file_name,
                                                    element_file_name,
                                                    face_file_name);
 
 // The following corresponds to the boundaries as specified by
 // facets in the Tetgen input:
 
 /// IDs of solid mesh boundaries where displacements are pinned
 Pinned_solid_boundary_id.resize(42);
 for(unsigned i=0; i<42; i++)
  {
   Pinned_solid_boundary_id[i]=215+i;
  }

  // The solid and fluid fsi boundaries are numbered int he same way.

 Solid_fsi_boundary_id.resize(215);
 for(unsigned i=0; i<215; i++)
  {
   Solid_fsi_boundary_id[i]=i;
  }



 // Create (empty) meshes of lagrange multiplier elements at inflow/outflow
 //------------------------------------------------------------------------
 
 // Create the meshes
 unsigned n=nfluid_traction_boundary();
 Parallel_flow_lagrange_multiplier_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Parallel_flow_lagrange_multiplier_mesh_pt[i]=new Mesh;
  } 
 
 // Populate them with elements
 create_parallel_flow_lagrange_elements();


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
  }

 // Setup FSI
 //----------
     
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
   Multi_domain_functions::Doc_boundary_coordinate_file.open(filename);
   
   // Setup FSI: Pass ID of fluid FSI boundary and associated
   // mesh of solid fsi traction elements.
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,3>
    (this,Fluid_fsi_boundary_id[i],Fluid_mesh_pt,Solid_fsi_traction_mesh_pt[i]);
   
   // Close the doc file
   Multi_domain_functions::Doc_boundary_coordinate_file.close();
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
    new MeshAsGeomObject
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
         el_pt->pressure_pt()= &Global_Parameters::P_in;
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

} // end_of_doc





//========================= start_of_main=================================
/// Demonstrate how to solve an unstructured 3D FSI problem
//========================================================================
int main(int argc, char **argv)
{

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

 // Parameter study
 unsigned nstep=2;

 // Increment in FSI parameter
 double q_increment=5.0e-2;

 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution(doc_info);
   doc_info.number()++;

   // Bump up FSI parameter
   Global_Parameters::Q+=q_increment;   
  }

} // end_of_main

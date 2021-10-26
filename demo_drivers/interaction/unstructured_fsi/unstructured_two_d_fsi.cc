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
//Driver for fsi on unstructured mesh

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"

// The meshes
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph; 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re=0.0; 

 /// FSI parameter
 double Q=0.0;

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }

 /// Pseudo-solid Poisson ratio
 double Nu=0.3;

 /// Constitutive law for the solid (and pseudo-solid) mechanics
 ConstitutiveLaw *Constitutive_law_pt=0;

 /// Boolean to identify if node is on fsi boundary
 bool is_on_fsi_boundary(Node* nod_pt)
 {
  if (
   (
    // Is it a boundary node?
    dynamic_cast<BoundaryNodeBase*>(nod_pt)!=0)&&
   (
    // Horizontal extent of main immersed obstacle
    (  (nod_pt->x(0)>1.6)&&(nod_pt->x(0)<4.75)&&
       // Vertical extent of main immersed obstacle
       (nod_pt->x(1)>0.1125)&&(nod_pt->x(1)<2.8) ) ||
    // Two nodes on the bottom wall are below y=0.3
    (  (nod_pt->x(1)<0.3)&&
       // ...and bracketed in these two x-ranges
       (  ( (nod_pt->x(0)>3.0)&&(nod_pt->x(0)<3.1) ) ||
          ( (nod_pt->x(0)<4.6)&&(nod_pt->x(0)>4.5) )   ) 
     )
    )
   )
   {
    return true;
   }
  else
   {
    return false;
   }
 }


} // end_of_namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//========start_solid_mesh=================================================
/// Triangle-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class MySolidTriangleMesh : public virtual TriangleMesh<ELEMENT>, 
                          public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 MySolidTriangleMesh(const std::string& node_file_name,
                     const std::string& element_file_name,
                     const std::string& poly_file_name,
                     TimeStepper* time_stepper_pt=
                     &Mesh::Default_TimeStepper) :
  TriangleMesh<ELEMENT>(node_file_name,element_file_name,
                        poly_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();

   // Identify special boundaries
   set_nboundary(3);

   unsigned n_node=this->nnode();
   for (unsigned j=0;j<n_node;j++)
    {
     Node* nod_pt=this->node_pt(j);

     // Boundary 1 is lower boundary
     if (nod_pt->x(1)<0.15)
      {
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(1,nod_pt);
      }

     // Boundary 2 is FSI interface
     if (Global_Parameters::is_on_fsi_boundary(nod_pt))
      {
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(2,nod_pt);
      }
    }// done boundary assignment

   // Identify the elements next to the newly created boundaries
   TriangleMesh<ELEMENT>::setup_boundary_element_info();
   
   // Setup boundary coordinates for boundaries 1 and 2
   this->template setup_boundary_coordinates<ELEMENT>(1);
   this->template setup_boundary_coordinates<ELEMENT>(2);
  }
 
 /// Empty Destructor
 virtual ~MySolidTriangleMesh() { }

};



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//======================start_fluid_mesh===================================
/// Triangle-based mesh upgraded to become a pseudo-solid mesh
//=========================================================================
template<class ELEMENT>
class FluidTriangleMesh : public virtual TriangleMesh<ELEMENT>,
                          public virtual SolidMesh 
{
 
public:
 
 /// Constructor
 FluidTriangleMesh(const std::string& node_file_name,
                   const std::string& element_file_name,
                   const std::string& poly_file_name,
                   TimeStepper* time_stepper_pt=
                   &Mesh::Default_TimeStepper) :
  TriangleMesh<ELEMENT>(node_file_name,element_file_name,
                        poly_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();

   // Identify special boundaries
   this->set_nboundary(4);

   unsigned n_node=this->nnode();
   for (unsigned j=0;j<n_node;j++)
    {
     Node* nod_pt=this->node_pt(j);

     // Boundary 1 is left (inflow) boundary
     if (nod_pt->x(0)<0.226)
      {
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(1,nod_pt);

       // Add overlapping nodes back to boundary 0
       if (nod_pt->x(1)<0.2) this->add_boundary_node(0,nod_pt);
       if (nod_pt->x(1)>4.06) this->add_boundary_node(0,nod_pt);
      }

     // Boundary 2 is right (outflow) boundary
     if (nod_pt->x(0)>8.28)
      {
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(2,nod_pt);

       // Add overlapping nodes back to boundary 0
       if (nod_pt->x(1)<0.2) this->add_boundary_node(0,nod_pt);
       if (nod_pt->x(1)>4.06) this->add_boundary_node(0,nod_pt);

      }

     // Boundary 3 is FSI boundary
     if (Global_Parameters::is_on_fsi_boundary(nod_pt))
      {
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(3,nod_pt);
       
       //If it's below y=0.2 it's also on boundary 0 so stick it back on
       if (nod_pt->x(1)<0.2) this->add_boundary_node(0,nod_pt);
      }
    }
   TriangleMesh<ELEMENT>::setup_boundary_element_info();
   
   // Open a file to doc the FaceElements that are used to
   // create the boundary coordinates. The elements must
   // form a simply-connected line. This may not work
   // if the mesh is too coarse so that, e.g. an element
   // that goes through the interior has endpoints that are
   // both located on the same boundary. Outputting the
   // FaceElements can help identify such cases.
   ofstream some_file("RESLT/boundary_generation_test.dat");
   
   // Setup boundary coordinates for boundaries 1, 2 and 3
   this->template setup_boundary_coordinates<ELEMENT>(1);
   this->template setup_boundary_coordinates<ELEMENT>(2);
   this->template setup_boundary_coordinates<ELEMENT>(3,some_file);

   // Close it again
   some_file.close();
  }
 
 /// Empty Destructor
 virtual ~FluidTriangleMesh() { }

};



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Unstructured FSI Problem
//====================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
class UnstructuredFSIProblem : public Problem
{

public:

 /// Constructor
 UnstructuredFSIProblem();

 /// Destructor (empty)
 ~UnstructuredFSIProblem(){}

 /// Access function for the fluid mesh
 FluidTriangleMesh<FLUID_ELEMENT>*& fluid_mesh_pt() 
  {
   return Fluid_mesh_pt;
  }

 /// Access function for the solid mesh
 MySolidTriangleMesh<SOLID_ELEMENT>*& solid_mesh_pt() 
  {
   return Solid_mesh_pt;
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Create FSI traction elements
 void create_fsi_traction_elements();

 /// Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 /// Sanity check: Doc boundary coordinates from solid side
 void doc_solid_boundary_coordinates();
 
 /// Fluid mesh
 FluidTriangleMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

 /// Solid mesh
 MySolidTriangleMesh<SOLID_ELEMENT>* Solid_mesh_pt;

 /// Pointers to mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// Vector of pointers to mesh of FSI traction elements
 SolidMesh* Traction_mesh_pt;

 /// GeomObject incarnation of fsi boundary in solid mesh
 MeshAsGeomObject* 
 Solid_fsi_boundary_pt;

}; 



//==start_of_constructor==================================================
/// Constructor for unstructured FSI problem.
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::UnstructuredFSIProblem()
{ 

 // Fluid mesh
 //-----------

 //Create fluid mesh
 string fluid_node_file_name="fluid.fig.1.node";
 string fluid_element_file_name="fluid.fig.1.ele";
 string fluid_poly_file_name="fluid.fig.1.poly"; 
 Fluid_mesh_pt = new FluidTriangleMesh<FLUID_ELEMENT>(fluid_node_file_name,
                                                      fluid_element_file_name,
                                                      fluid_poly_file_name);
 
 // Doc pinned solid nodes
 std::ofstream pseudo_solid_bc_file("pinned_pseudo_solid_nodes.dat");

 // Set the boundary conditions for fluid problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin velocity everywhere apart from outlet where we
     // have parallel outflow
     if (ibound!=2)
      {
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }
     Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(1); 

     // Pin pseudo-solid positions everywhere apart from boundary 3, 
     // the fsi boundary 
     if ((ibound==0)||(ibound==1)||(ibound==2))
      {
       for(unsigned i=0;i<2;i++)
        {
         // Pin the node
         SolidNode* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
         nod_pt->pin_position(i);

         // Doc it as pinned
         pseudo_solid_bc_file << nod_pt->x(i) << " ";
        }
       pseudo_solid_bc_file << std::endl;
      }
    }
  } // end loop over boundaries

 // Close
 pseudo_solid_bc_file.close();

 
 // Complete the build of the fluid elements so they are fully functional
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


 // Apply fluid boundary conditions: Poiseuille at inflow

 // Find max. and min y-coordinate at inflow
 unsigned ibound=1;
 //Initialise both to the y-coordinate of the first boundary node
 double y_min=fluid_mesh_pt()->boundary_node_pt(ibound,0)->x(1);;
 double y_max=y_min;

 //Loop over the rest of the boundary nodes
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
 const unsigned n_boundary = fluid_mesh_pt()->nboundary();
 for (unsigned ibound=0;ibound<n_boundary;ibound++)
  {
   const unsigned num_nod= fluid_mesh_pt()->nboundary_node(ibound);
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
     // Zero flow elsewhere 
     else
      {
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       fluid_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
    }
  } // end Poiseuille
 
 
 // Solid mesh
 //-----------
 
 //Create solid mesh
 string solid_node_file_name="solid.fig.1.node";
 string solid_element_file_name="solid.fig.1.ele";
 string solid_poly_file_name="solid.fig.1.poly"; 
 Solid_mesh_pt = new MySolidTriangleMesh<SOLID_ELEMENT>(solid_node_file_name,
                                                     solid_element_file_name,
                                                     solid_poly_file_name);

 
 // Complete the build of all solid elements so they are fully functional
 n_element = Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   SOLID_ELEMENT *el_pt = 
    dynamic_cast<SOLID_ELEMENT*>(Solid_mesh_pt->element_pt(i));

   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Parameters::Constitutive_law_pt;
   
   //Set the body force
   el_pt->body_force_fct_pt() = Global_Parameters::gravity;
  }
 
 // Pin both positions at lower boundary of solid mesh (boundary 1)
 ibound=1;
 num_nod=Solid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {    
   Solid_mesh_pt->boundary_node_pt(ibound,inod)->pin_position(0);
   Solid_mesh_pt->boundary_node_pt(ibound,inod)->pin_position(1);
  }


 // Create FSI Traction elements
 //-----------------------------

 // Now construct the (empty) traction element mesh
 Traction_mesh_pt=new SolidMesh;

  // Build the FSI traction elements and add them to the traction mesh
 create_fsi_traction_elements();

 // Create Lagrange multiplier mesh for boundary motion
 //----------------------------------------------------
 
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();


 // Combine meshes
 //---------------

 // Add sub meshes
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Solid_mesh_pt);
 add_sub_mesh(Traction_mesh_pt);
 add_sub_mesh(Lagrange_multiplier_mesh_pt);
 
 // Build global mesh
 build_global_mesh();


 // Setup FSI
 //----------

 // Document the boundary coordinate  along the FSI interface 
 // of the fluid mesh during call to
 // setup_fluid_load_info_for_solid_elements()
 Multi_domain_functions::Doc_boundary_coordinate_file.open(
  "fluid_boundary_test.dat");

 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 3
 // of the 2D fluid mesh.
 FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
  (this,3,Fluid_mesh_pt,Traction_mesh_pt);
 
 // Close the doc file
 Multi_domain_functions::Doc_boundary_coordinate_file.close();

 // Sanity check: Doc boundary coordinates from solid side
 doc_solid_boundary_coordinates();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor



//============start_doc_solid_zeta=======================================
/// Doc boundary coordinates in solid and plot GeomObject representation
/// of FSI boundary.
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
doc_solid_boundary_coordinates()
{

 // Doc boundary coordinates for fsi boundary in solid mesh
 std::ofstream the_file("solid_boundary_test.dat");

 // Initialise max/min boundary coordinate
 double zeta_min= 1.0e40;
 double zeta_max=-1.0e40;

 // Loop over FSI traction elements
 unsigned n_face_element = Traction_mesh_pt->nelement();
 for(unsigned e=0;e<n_face_element;e++)
  {
   //Cast the element pointer
   FSISolidTractionElement<SOLID_ELEMENT,2>* el_pt=
    dynamic_cast<FSISolidTractionElement<SOLID_ELEMENT,2>*>
    (Traction_mesh_pt->element_pt(e));

   // Doc boundary coordinate
   Vector<double> s(1);
   Vector<double> zeta(1);
   Vector<double> x(2);
   unsigned n_plot=5;
   the_file << el_pt->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {         
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,n_plot,s);         
     el_pt->interpolated_zeta(s,zeta);
     el_pt->interpolated_x(s,x);
     for (unsigned i=0;i<2;i++)
      {
       the_file << x[i] << " ";
      }
     the_file << zeta[0] << " ";

     // Update max/min boundary coordinate
     if (zeta[0]<zeta_min) zeta_min=zeta[0];
     if (zeta[0]>zeta_max) zeta_max=zeta[0];

     the_file << std::endl;
    }
  } 

 // Close doc file
 the_file.close();


 // Doc compound GeomObject
 the_file.open("fsi_geom_object.dat");
 unsigned nplot=10000;
 Vector<double> zeta(1);
 Vector<double> r(2);
 for (unsigned i=0;i<nplot;i++)
  {
   zeta[0]=zeta_min+(zeta_max-zeta_min)*double(i)/double(nplot-1);
   Solid_fsi_boundary_pt->position(zeta,r);
   the_file << r[0] << " " << r[1] << std::endl;
  }
 the_file.close();
 
} //end doc_solid_zeta




//============start_of_create_traction_elements==========================
/// Create FSI traction elements 
//=======================================================================
template<class FLUID_ELEMENT,class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_fsi_traction_elements()
{
 // Traction elements are located on boundary 2 of solid bulk mesh
 unsigned b=2;
 
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //What is the index of the face of the element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
  
   // Create new element 
   FSISolidTractionElement<SOLID_ELEMENT,2>* el_pt=
    new FSISolidTractionElement<SOLID_ELEMENT,2>(bulk_elem_pt,face_index);
   
   // Add it to the mesh
   Traction_mesh_pt->add_element_pt(el_pt);
   
   // Specify boundary number
   el_pt->set_boundary_number_in_bulk_mesh(b);
   
   // Function that specifies the load ratios
   el_pt->q_pt() = &Global_Parameters::Q; 
  } 

 } // end of create_traction_elements



//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
void UnstructuredFSIProblem<FLUID_ELEMENT,SOLID_ELEMENT>::
create_lagrange_multiplier_elements()
{

 // Create  GeomObject incarnation of fsi boundary in solid mesh
 Solid_fsi_boundary_pt=
  new MeshAsGeomObject
  (Traction_mesh_pt);
 
 // Lagrange multiplier elements are located on boundary 3 of the fluid mesh
 unsigned b=3;

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
    new ImposeDisplacementByLagrangeMultiplierElement<FLUID_ELEMENT>(
     bulk_elem_pt,face_index);   
   
   // Add it to the mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);

   // Set the GeomObject that defines the boundary shape and set
   // which bulk boundary we are attached to (needed to extract
   // the boundary coordinate from the bulk nodes)
   el_pt->set_boundary_shape_geom_object_pt(Solid_fsi_boundary_pt,b);
   
   // Loop over the nodes to apply boundary conditions
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = el_pt->node_pt(j);
     
     // Is the node also on boundary 0? 
     if (nod_pt->is_on_boundary(0))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=el_pt->nbulk_value(j);
       
        // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
 
} // end of create_lagrange_multiplier_elements


//==start_of_doc_solution=================================================
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

 // Output fluid solution 
 sprintf(filename,"%s/fluid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output solid solution 
 sprintf(filename,"%s/solid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();


} // end_of_doc_solution




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for unstructured fsi problem
//=====================================================================
int main()
{
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 //Create the constitutive law
 Global_Parameters::Constitutive_law_pt = new GeneralisedHookean(
  &Global_Parameters::Nu);
 
 // Build the problem with triangular Taylor Hood for fluid and solid
 UnstructuredFSIProblem<
 PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, TPVDElement<2,3> >, 
  TPVDElement<2,3> > problem;

 // Output boundaries 
 problem.fluid_mesh_pt()->output_boundaries("RESLT/fluid_boundaries.dat");
 problem.solid_mesh_pt()->output_boundaries("RESLT/solid_boundaries.dat");
 
 // Output the initial guess for the solution
 problem.doc_solution(doc_info);
 doc_info.number()++;
 
 // Parameter study
 Global_Parameters::Gravity=2.0e-4;
 double q_increment=1.0e-6;

 // Solve the problem at zero Re and Q
 problem.newton_solve();
 
 // Output the solution
 problem.doc_solution(doc_info);
 doc_info.number()++;
 
 // Bump up Re
 Global_Parameters::Re=10.0;

 // Now do proper parameter study with increase in Q
 unsigned nstep=2; // 10;
 for (unsigned i=0;i<nstep;i++)
  {
   // Solve the problem
   problem.newton_solve();
   
   // Output the solution
   problem.doc_solution(doc_info);
   doc_info.number()++;

   // Bump up Q
   Global_Parameters::Q+=q_increment;
  }
 

} // end_of_main











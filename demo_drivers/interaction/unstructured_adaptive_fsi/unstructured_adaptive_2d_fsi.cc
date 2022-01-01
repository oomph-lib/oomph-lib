//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Driver code for a simple unstructured FSI problem using a mesh
// generated from within the driver.
// The problem is the flow past a (thin) solid leaflet.
// The fluid mesh is of the form:
//
//                                 5
//    (x_i,H) -----------------------------------------------------(x_i+L,H)
//            |                   1                               |
//            |       (x_l-w/2,h)----(x_l+w/2,h)                  |
//          6 |                  |  |                             |4
//            |                0 |  | 2                           |
//            |                  |  |                             |
//     (x_i,0)-------------------|  |------------------------------(x_i+L,0)
//                   7 (x_l-w/2,0)  (x_l+w/2,0)       3
//
// where the key variables are x_i = x_inlet, H = channel_height, 
// L = channel_length, x_l = x_leaflet, w = leaflet_width, h = leaflet_height
// and the boundaries are labelled as shown.
//
// The solid domain is (initially) a simple rectangle with the following
// boundary labels and positions:
//                        1
//           (x_l-w/2,h)------(x_l+w/2,h)
//                      |    |
//                    0 |    |2
//                      |    |
//                      |    |
//           (x_l-w/2,0)------(x_l+w/2,0)
//                        3
// For convenience, the coincident solid and fluid boundaries have the same
// boundary labels in each mesh, but this is not required.

//Generic routines
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"


// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re = 0.0;

 /// FSI parameter
 double Q = 0.0;

 /// Poisson's ratio
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;
 
 /// Mesh poisson ratio
 double Mesh_Nu = 0.1;

 /// Pointer to constitutive law for the mesh
 ConstitutiveLaw* Mesh_constitutive_law_pt=0;

} //end namespace



//==============start_problem=========================================
/// Unstructured FSI problem
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

 /// Actions before adapt
 void actions_before_adapt()
  {
   //Delete the boundary meshes
   this->delete_lagrange_multiplier_elements();
   this->delete_fsi_traction_elements();

   //Rebuild the global mesh
   this->rebuild_global_mesh();
  }


 /// Actions after adapt
 void actions_after_adapt()
  {
   //Ensure that the lagrangian coordinates of the mesh are set to be
   //the same as the eulerian
   Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

   //Apply boundary conditions again
   
   // Pin both positions at lower boundary (boundary 3)
   unsigned ibound=3;
   unsigned num_nod= Solid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {  
     
     // Get node
     SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin both directions
     for (unsigned i=0;i<2;i++)
      {
       nod_pt->pin_position(i);
      }
    }
   
   // Complete the build of all elements so they are fully functional
   unsigned n_element = Solid_mesh_pt->nelement();
   for(unsigned i=0;i<n_element;i++)
    {
     //Cast to a solid element
     SOLID_ELEMENT *el_pt = 
      dynamic_cast<SOLID_ELEMENT*>(Solid_mesh_pt->element_pt(i));
     
     // Set the constitutive law
     el_pt->constitutive_law_pt() =
      Global_Physical_Variables::Constitutive_law_pt;
    } // end complete solid build

   
   // Set the boundary conditions for fluid problem: All nodes are
   // free by default 
   // --- just pin the ones that have Dirichlet conditions here. 
   
   //Pin velocity everywhere apart from parallel outflow (boundary 4)
   unsigned nbound=Fluid_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<nbound;ibound++)
    {
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Pin velocity everywhere apart from outlet where we
       // have parallel outflow
       if (ibound!=4)
        {
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
        }
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(1); 
       
       // Pin pseudo-solid positions everywhere apart from boundaries 0, 1, 2 
       // the fsi boundaries
       if(ibound > 2)
        {
         for(unsigned i=0;i<2;i++)
          {
           // Pin the node
           SolidNode* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
           nod_pt->pin_position(i);
          }
        }
      }
    } // end loop over boundaries
   
   
   // Complete the build of the fluid elements so they are fully functional
   n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     FLUID_ELEMENT* el_pt = 
      dynamic_cast<FLUID_ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     
     //Set the Reynolds number
     el_pt->re_pt() = &Global_Physical_Variables::Re;
     
     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt() =
      Global_Physical_Variables::Mesh_constitutive_law_pt;
     
    } // end loop over elements
   
   
   // Apply fluid boundary conditions: Poiseuille at inflow
   const unsigned n_boundary = Fluid_mesh_pt->nboundary();
   for (unsigned ibound=0;ibound<n_boundary;ibound++)
    {
     const unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Parabolic inflow at the left boundary (boundary 6)
       if(ibound==6)
        {
         double y=Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
         double veloc = y*(1.0-y);
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,veloc);
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }
       // Zero flow elsewhere
       else 
        {
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }
      }
    } // end Poiseuille
   
   //Recreate the boundary elements
   this->create_fsi_traction_elements();
   this->create_lagrange_multiplier_elements();
   
   //Rebuild the global mesh
   this->rebuild_global_mesh();

 // Setup FSI (again)
 //------------------
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary are boundaries 0, 1 and 2
 // of the 2D fluid mesh.
 for(unsigned b=0;b<3;b++)
  {
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,b,Fluid_mesh_pt,Traction_mesh_pt[b]);
  }
 
  }

 /// Output function to compute the strain energy in the solid and the 
 /// dissipation in the fluid and write to the output stream trace
 void output_strain_and_dissipation(std::ostream &trace)
  {
   double strain_energy = this->get_solid_strain_energy();
   double dissipation = this->get_fluid_dissipation();

   trace << Global_Physical_Variables::Q <<
    " " << dissipation << " " << strain_energy << std::endl;
  }


private:

 /// Create the traction element
 void create_fsi_traction_elements()
  {
   // Traction elements are located on boundaries 0 1 and 2 of solid bulk mesh
   for(unsigned b=0;b<3;++b)
    {
     // How many bulk elements are adjacent to boundary b?
     const unsigned n_element = Solid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       SOLID_ELEMENT* bulk_elem_pt = dynamic_cast<SOLID_ELEMENT*>(
        Solid_mesh_pt->boundary_element_pt(b,e));
       
       //What is the index of the face of the element e along boundary b
       int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element 
       FSISolidTractionElement<SOLID_ELEMENT,2>* el_pt=
        new FSISolidTractionElement<SOLID_ELEMENT,2>(bulk_elem_pt,face_index);
   
       // Add it to the mesh
       Traction_mesh_pt[b]->add_element_pt(el_pt);
   
       // Specify boundary number
       el_pt->set_boundary_number_in_bulk_mesh(b);
   
       // Function that specifies the load ratios
       el_pt->q_pt() = &Global_Physical_Variables::Q; 
      } 
    }
  }

 /// Delete the traction elements
 void delete_fsi_traction_elements()
  {
   //There are three traction element boundaries
   for(unsigned b=0;b<3;++b)
    {
     const unsigned n_element = Traction_mesh_pt[b]->nelement();
     //Delete the elements
     for(unsigned e=0;e<n_element;e++)
      {
       delete Traction_mesh_pt[b]->element_pt(e);
      } 
     //Wipe the mesh
     Traction_mesh_pt[b]->flush_element_and_node_storage();
     //Also wipe out the Mesh as Geometric objects
     delete Solid_fsi_boundary_pt[b]; Solid_fsi_boundary_pt[b] = 0;
    }
  }

/// Create the multipliers that add lagrange multipliers to the fluid
/// elements that apply the solid displacement conditions
 void create_lagrange_multiplier_elements() 
  {
   //Loop over the FSI boundaries
   for(unsigned b=0;b<3;b++)
    {
     // Create GeomObject incarnation of fsi boundary in solid mesh
     Solid_fsi_boundary_pt[b] = new MeshAsGeomObject(Traction_mesh_pt[b]);
     
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
       Lagrange_multiplier_mesh_pt[b]->add_element_pt(el_pt);
       
       // Set the GeomObject that defines the boundary shape and set
       // which bulk boundary we are attached to (needed to extract
       // the boundary coordinate from the bulk nodes)
       el_pt->set_boundary_shape_geom_object_pt(Solid_fsi_boundary_pt[b],b);
       
       // Loop over the nodes to apply boundary conditions
       unsigned nnod=el_pt->nnode();
       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt = el_pt->node_pt(j);
         
         // How many nodal values were used by the "bulk" element
         // that originally created this node?
         unsigned n_bulk_value=el_pt->nbulk_value(j);
         
         // The remaining ones are Lagrange multipliers
         unsigned nval=nod_pt->nvalue();
         //If we have more than two Lagrange multipliers, pin the rest
         if(nval > n_bulk_value + 2)
          {
           for (unsigned j=n_bulk_value+2;j<nval;j++)
            {
             nod_pt->pin(j);
            }
          }

         //If I'm also on one of the base boundaries, pin the Lagrange
         //Multipliers
         if((nod_pt->is_on_boundary(7)) || (nod_pt->is_on_boundary(3)))
          {
           for(unsigned j=n_bulk_value;j<nval;j++)
            {
             nod_pt->pin(j);
            }
          }

        }
      }
    }
  } // end of create_lagrange_multiplier_elements


 /// Delete the traction elements
 void delete_lagrange_multiplier_elements()
  {
   //There are three Lagrange-multiplier element boundaries
   for(unsigned b=0;b<3;++b)
    {
     const unsigned n_element = Lagrange_multiplier_mesh_pt[b]->nelement();
     //Delete the elements
     for(unsigned e=0;e<n_element;e++)
      {
       delete Lagrange_multiplier_mesh_pt[b]->element_pt(e);
      } 
     //Wipe the mesh
     Lagrange_multiplier_mesh_pt[b]->flush_element_and_node_storage();
    }
  }



 /// Calculate the strain energy of the solid
 double get_solid_strain_energy()
  {
   double strain_energy=0.0;
   const unsigned n_element = Solid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a solid element
     SOLID_ELEMENT *el_pt = 
      dynamic_cast<SOLID_ELEMENT*>(Solid_mesh_pt->element_pt(e));
     
     double pot_en, kin_en;
     el_pt->get_energy(pot_en,kin_en);
     strain_energy += pot_en;
    }
   return strain_energy;
  }

 /// Calculate the fluid dissipation
 double get_fluid_dissipation()
  {
   double dissipation=0.0;
   const unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     FLUID_ELEMENT *el_pt = 
      dynamic_cast<FLUID_ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     //Add to the dissipation
     dissipation += el_pt->dissipation(); 
    }
   return dissipation;
  }

 /// Bulk solid mesh
 RefineableSolidTriangleMesh<SOLID_ELEMENT>* Solid_mesh_pt;

public: 
 /// Bulk fluid mesh
 RefineableSolidTriangleMesh<FLUID_ELEMENT>* Fluid_mesh_pt;

private:

 /// Vector of pointers to mesh of Lagrange multiplier elements
 Vector<SolidMesh*> Lagrange_multiplier_mesh_pt;

 /// Vectors of pointers to mesh of traction elements
 Vector<SolidMesh*> Traction_mesh_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Solid_outer_boundary_polyline_pt; 

  /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Fluid_outer_boundary_polyline_pt; 

 // GeomObject incarnation of fsi boundaries in solid mesh
 Vector<MeshAsGeomObject*> Solid_fsi_boundary_pt;

};



//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class FLUID_ELEMENT, class SOLID_ELEMENT>
UnstructuredFSIProblem<FLUID_ELEMENT, SOLID_ELEMENT>::UnstructuredFSIProblem() 
{  

 //Some geometric parameters
 double x_inlet = 0.0;
 double channel_height = 1.0;
 double channel_length = 4.0;
 double x_leaflet = 1.0;
 double leaflet_width = 0.2;
 double leaflet_height = 0.5;

 // Solid Mesh
 //---------------

 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshCurveSection*> solid_boundary_segment_pt(4);
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(2);
 for(unsigned i=0;i<2;i++) {bound_seg[i].resize(2);}
 
 // First boundary segment
 bound_seg[0][0]=x_leaflet - 0.5*leaflet_width;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=x_leaflet - 0.5*leaflet_width;
 bound_seg[1][1]=leaflet_height;
 
 // Specify 1st boundary id
 unsigned bound_id = 0;

 // Build the 1st boundary segment
 solid_boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Second boundary segment
 bound_seg[0][0]=x_leaflet - 0.5*leaflet_width;
 bound_seg[0][1]=leaflet_height;
 bound_seg[1][0]=x_leaflet + 0.5*leaflet_width;
 bound_seg[1][1]=leaflet_height;

 // Specify 2nd boundary id
 bound_id = 1;

 // Build the 2nd boundary segment
 solid_boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Third boundary segment
 bound_seg[0][0]=x_leaflet + 0.5*leaflet_width;
 bound_seg[0][1]=leaflet_height;
 bound_seg[1][0]=x_leaflet + 0.5*leaflet_width;
 bound_seg[1][1]=0.0;

 // Specify 3rd boundary id
 bound_id = 2;

 // Build the 3rd boundary segment
 solid_boundary_segment_pt[2] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fourth boundary segment
 bound_seg[0][0]=x_leaflet + 0.5*leaflet_width;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=x_leaflet - 0.5*leaflet_width;
 bound_seg[1][1]=0.0;

 // Specify 4th boundary id
 bound_id = 3;

 // Build the 4th boundary segment
 solid_boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Solid_outer_boundary_polyline_pt = 
  new TriangleMeshPolygon(solid_boundary_segment_pt);

 // There are no holes
 //-------------------------------
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------
 double uniform_element_area= leaflet_width*leaflet_height/20.0;

 TriangleMeshClosedCurve* solid_closed_curve_pt=
  Solid_outer_boundary_polyline_pt;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters_solid(
   solid_closed_curve_pt);

 // Define the maximum element area
 triangle_mesh_parameters_solid.element_area() =
   uniform_element_area;

 // Create the mesh
 Solid_mesh_pt =
   new RefineableSolidTriangleMesh<SOLID_ELEMENT>(
     triangle_mesh_parameters_solid);

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Solid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Solid_mesh_pt->max_permitted_error()=0.0001;
 Solid_mesh_pt->min_permitted_error()=0.001; 
 Solid_mesh_pt->max_element_size()=0.2;
 Solid_mesh_pt->min_element_size()=0.001; 
   
 // Output boundary and mesh
 this->Solid_mesh_pt->output_boundaries("solid_boundaries.dat");
 this->Solid_mesh_pt->output("solid_mesh.dat");
 
 // Pin both positions at lower boundary (boundary 3)
 unsigned ibound=3;
 unsigned num_nod= Solid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {  

   // Get node
   SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
   
   // Pin both directions
   for (unsigned i=0;i<2;i++)
    {
     nod_pt->pin_position(i);
    }
  } // end_solid_boundary_conditions

 // Complete the build of all elements so they are fully functional
 unsigned n_element = Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   SOLID_ELEMENT *el_pt = 
    dynamic_cast<SOLID_ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
  }


 // Fluid Mesh
 //--------------

 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshCurveSection*> fluid_boundary_segment_pt(8);
 
 //The first three boundaries should be in common with the solid
 for(unsigned b=0;b<3;b++)
  {
   fluid_boundary_segment_pt[b] = solid_boundary_segment_pt[b];
  }

 //Now fill in the rest 
 // Fourth boundary segment
 bound_seg[0][0]=x_leaflet + 0.5*leaflet_width;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=x_inlet + channel_length;
 bound_seg[1][1]=0.0;
 
 // Specify 4th boundary id
 bound_id = 3;

 // Build the 4th boundary segment
 fluid_boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fifth boundary segment
 bound_seg[0][0]=x_inlet + channel_length;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=x_inlet + channel_length;
 bound_seg[1][1]=channel_height;
 
 // Specify 5th boundary id
 bound_id = 4;

 // Build the 4th boundary segment
 fluid_boundary_segment_pt[4] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Sixth boundary segment
 bound_seg[0][0]=x_inlet + channel_length;
 bound_seg[0][1]=channel_height;
 bound_seg[1][0]=x_inlet;
 bound_seg[1][1]=channel_height;
 
 // Specify 6th boundary id
 bound_id = 5;
 
 // Build the 6th boundary segment
 fluid_boundary_segment_pt[5] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Seventh boundary segment
 bound_seg[0][0]=x_inlet;
 bound_seg[0][1]=channel_height;
 bound_seg[1][0]=x_inlet;
 bound_seg[1][1]=0.0;
 
 // Specify 7th boundary id
 bound_id = 6;

 // Build the 7th boundary segment
 fluid_boundary_segment_pt[6] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Eighth boundary segment
 bound_seg[0][0]=x_inlet;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=x_leaflet - 0.5*leaflet_width;
 bound_seg[1][1]=0.0;
 
 // Specify 8th boundary id
 bound_id = 7;

 // Build the 8th boundary segment
 fluid_boundary_segment_pt[7] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Fluid_outer_boundary_polyline_pt = 
  new TriangleMeshPolygon(fluid_boundary_segment_pt);

 // There are no holes
 //-------------------------------
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------
 uniform_element_area= channel_length*channel_height/40.0;;

 TriangleMeshClosedCurve* fluid_closed_curve_pt=
  Fluid_outer_boundary_polyline_pt;
 
 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters_fluid(
   fluid_closed_curve_pt);

 // Define the maximum element area
 triangle_mesh_parameters_fluid.element_area() =
   uniform_element_area;

 // Create the mesh
 Fluid_mesh_pt =
   new RefineableSolidTriangleMesh<FLUID_ELEMENT>(
     triangle_mesh_parameters_fluid);

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* fluid_error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=fluid_error_estimator_pt;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=0.0001;
 Fluid_mesh_pt->min_permitted_error()=0.001; 
 Fluid_mesh_pt->max_element_size()=0.2;
 Fluid_mesh_pt->min_element_size()=0.001; 
   
 // Output boundary and mesh
 this->Fluid_mesh_pt->output_boundaries("fluid_boundaries.dat");
 this->Fluid_mesh_pt->output("fluid_mesh.dat");

 // Set the boundary conditions for fluid problem: All nodes are
 // free by default 
 // --- just pin the ones that have Dirichlet conditions here. 

 //Pin velocity everywhere apart from parallel outflow (boundary 4)
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin velocity everywhere apart from outlet where we
     // have parallel outflow
     if (ibound!=4)
      {
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }
     Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(1); 

     // Pin pseudo-solid positions everywhere apart from boundaries 0, 1, 2 
     // the fsi boundaries
     if(ibound > 2)
      {
       for(unsigned i=0;i<2;i++)
        {
         // Pin the node
         SolidNode* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
         nod_pt->pin_position(i);
        }
      }
    }
  } // end loop over boundaries

 
 // Complete the build of the fluid elements so they are fully functional
 n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   FLUID_ELEMENT* el_pt = 
    dynamic_cast<FLUID_ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Mesh_constitutive_law_pt;
   
  } // end loop over elements


 // Apply fluid boundary conditions: Poiseuille at inflow
 const unsigned n_boundary = Fluid_mesh_pt->nboundary();
 for (unsigned ibound=0;ibound<n_boundary;ibound++)
  {
   const unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parabolic inflow at the left boundary (boundary 6)
     if(ibound==6)
      {
       double y=Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
       double veloc = y*(1.0-y);
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,veloc);
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
     // Zero flow elsewhere
     else 
      {
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
    }
  } // end Poiseuille
 

 // Make traction mesh 
 //(This must be done first because the resulting meshes are used
 // as the geometric objects that set the boundary locations of the fluid
 // mesh, as enforced by the Lagrange multipliers)
 Traction_mesh_pt.resize(3);
 for(unsigned m=0;m<3;m++) {Traction_mesh_pt[m] = new SolidMesh;}
 this->create_fsi_traction_elements();
 
 //Make the Lagrange multiplier mesh
 Lagrange_multiplier_mesh_pt.resize(3);
 Solid_fsi_boundary_pt.resize(3);
 for(unsigned m=0;m<3;m++) {Lagrange_multiplier_mesh_pt[m] = new SolidMesh;}
 this->create_lagrange_multiplier_elements();

 // Add sub meshes
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Solid_mesh_pt);
 for(unsigned m=0;m<3;m++)
  {
   add_sub_mesh(Traction_mesh_pt[m]);
   add_sub_mesh(Lagrange_multiplier_mesh_pt[m]);
  }
 
 // Build global mesh
 build_global_mesh();

 // Setup FSI
 //----------
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary are boundaries 0, 1 and 2
 // of the 2D fluid mesh.
 for(unsigned b=0;b<3;b++)
  {
   FSI_functions::setup_fluid_load_info_for_solid_elements<FLUID_ELEMENT,2>
    (this,b,Fluid_mesh_pt,Traction_mesh_pt[b]);
  }
   
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} //end constructor


//========================================================================
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

 // Output solution
 //----------------
 sprintf(filename,"%s/solid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 //----------------
 sprintf(filename,"%s/fluid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();

}

 



//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");

 //Create a trace file
 std::ofstream trace("RESLT/trace.dat");
 
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);

 // Create generalised Hookean constitutive equations for the mesh as well
 Global_Physical_Variables::Mesh_constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Mesh_Nu);
 
 //Set up the problem
 UnstructuredFSIProblem<
 ProjectableTaylorHoodElement<
 PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, TPVDElement<2,3> > >, 
  ProjectablePVDElement<TPVDElement<2,3> > > problem;

//Output initial configuration
problem.doc_solution(doc_info);
doc_info.number()++;

// Solve the problem
problem.newton_solve();

//Output solution
problem.doc_solution(doc_info);
doc_info.number()++;

//Calculate the strain energy of the solid and dissipation in the
//fluid as global output measures of the solution for validation purposes
problem.output_strain_and_dissipation(trace);

//Number of steps to be taken
unsigned n_step = 10;
//Reduce the number of steps if validating
if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
 {
  n_step=3;
 }

//Now Crank up interaction
for(unsigned i=0;i<n_step;i++)
 {
  Global_Physical_Variables::Q += 1.0e-4;
  problem.newton_solve(1);

  //Reset the lagrangian nodal coordinates in the fluid mesh
  //(Obviously we shouldn't do this in the solid mesh)
  problem.Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  //Output solution
  problem.doc_solution(doc_info);
  doc_info.number()++;

  //Calculate the strain energy of the solid and dissipation in the
  //fluid as global output measures of the solution for validation purposes
  problem.output_strain_and_dissipation(trace);
 }

trace.close();
} // end main




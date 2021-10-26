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

//Generic routines
#include "generic.h"


// The equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
 

namespace oomph
{


//==============================================================
/// Overload CrouzeixRaviart element to modify output
//==============================================================
 class MyCrouzeixRaviartElement : 
  public virtual PseudoSolidNodeUpdateElement<TCrouzeixRaviartElement<2>, 
  TPVDBubbleEnrichedElement<2,3> >
 {
  
 private:
  
  /// Storage for elemental error estimate -- used for post-processing
  double Error;

 public:

  /// Constructor initialise error
  MyCrouzeixRaviartElement() : PseudoSolidNodeUpdateElement<TCrouzeixRaviartElement<2>, 
                                                            TPVDBubbleEnrichedElement<2,3> >()
   {
    Error=0.0;
   }

  /// Set error value for post-processing
  void set_error(const double& error){Error=error;}
  
  /// Return variable identifier
  std::string variable_identifier()
   {
    std::string txt="VARIABLES=";
    txt+="\"x\",";
    txt+="\"y\",";
    txt+="\"u\",";
    txt+="\"v\",";
    txt+="\"p\",";   
    txt+="\"du/dt\",";
    txt+="\"dv/dt\",";
    txt+="\"u_m\",";   
    txt+="\"v_m\",";
    txt+="\"x_h1\",";
    txt+="\"y_h1\",";   
    txt+="\"x_h2\",";
    txt+="\"y_h2\",";   
    txt+="\"u_h1\",";
    txt+="\"v_h1\",";   
    txt+="\"u_h2\",";
    txt+="\"v_h2\",";   
    txt+="\"error\",";   
    txt+="\"size\",";   
    txt+="\n";
    return txt;
   }

  
  /// Overload output function
  void output(std::ostream &outfile, 
              const unsigned &nplot)
   {
    
    // Assign dimension 
    unsigned el_dim=2;
    
    // Vector of local coordinates
    Vector<double> s(el_dim);
    
    // Acceleration
    Vector<double> dudt(el_dim);
    
    // Mesh elocity
    Vector<double> mesh_veloc(el_dim,0.0);
   
    // Tecplot header info
    outfile << tecplot_zone_string(nplot);
   
    // Find out how many nodes there are
    unsigned n_node = nnode();
   
    //Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node,el_dim);
   
    // Loop over plot points
    unsigned num_plot_points=nplot_points(nplot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
     {
     
      // Get local coordinates of plot point
      get_s_plot(iplot,nplot,s);
     
      //Call the derivatives of the shape and test functions
      dshape_eulerian(s,psif,dpsifdx);
     
      //Allocate storage
      Vector<double> mesh_veloc(el_dim);
      Vector<double> dudt(el_dim);
      Vector<double> dudt_ALE(el_dim);
      DenseMatrix<double> interpolated_dudx(el_dim,el_dim);
     
      //Initialise everything to zero
      for(unsigned i=0;i<el_dim;i++)
       {
        mesh_veloc[i]=0.0;
        dudt[i]=0.0;
        dudt_ALE[i]=0.0;
        for(unsigned j=0;j<el_dim;j++)
         {
          interpolated_dudx(i,j) = 0.0;
         }
       }
     
      //Calculate velocities and derivatives

      //Loop over directions
      for(unsigned i=0;i<el_dim;i++)
       {
        //Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over nodes
        for(unsigned l=0;l<n_node;l++) 
         {
          dudt[i]+=du_dt_nst(l,u_nodal_index)*psif[l];
          mesh_veloc[i]+=dnodal_position_dt(l,i)*psif[l];
          
          //Loop over derivative directions for velocity gradients
          for(unsigned j=0;j<el_dim;j++)
           {                               
            interpolated_dudx(i,j) += nodal_value(l,u_nodal_index)*
             dpsifdx(l,j);
           }
         }
       }
     
     
      // Get dudt in ALE form (incl mesh veloc)
      for(unsigned i=0;i<el_dim;i++)
       {
        dudt_ALE[i]=dudt[i];
        for (unsigned k=0;k<el_dim;k++)
         {
          dudt_ALE[i]-=mesh_veloc[k]*interpolated_dudx(i,k);
         }
       }
     
     
      // Coordinates
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_x(s,i) << " ";
       }
     
      // Velocities
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_u_nst(s,i) << " ";
       }
     
      // Pressure
      outfile << interpolated_p_nst(s)  << " ";
     
      // Accelerations
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << dudt_ALE[i] << " ";
       }
     
      // Mesh velocity
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << mesh_veloc[i] << " ";
       }
     
      // History values of coordinates
      unsigned n_prev=node_pt(0)->position_time_stepper_pt()->ntstorage();
      for (unsigned t=1;t<n_prev;t++)
       {
        for(unsigned i=0;i<el_dim;i++) 
         {
          outfile << interpolated_x(t,s,i) << " ";
         }
       }
     
      // History values of velocities
      n_prev=node_pt(0)->time_stepper_pt()->ntstorage();
      for (unsigned t=1;t<n_prev;t++)
       {
        for(unsigned i=0;i<el_dim;i++) 
         {
          outfile << interpolated_u_nst(t,s,i) << " ";
         }
       }

      outfile << Error << " " 
              << size() << std::endl;        
     }
    
    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile,nplot); 
    }

 };


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<MyCrouzeixRaviartElement>
  : public virtual SolidTElement<1,3> 
 {
 public:
  FaceGeometry() : SolidTElement<1,3>() {}
 };


//=======================================================================
/// Face geometry of Face geometry for element is the same 
/// as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<FaceGeometry<MyCrouzeixRaviartElement> >
  : public virtual SolidPointElement 
 {
 public:
  FaceGeometry() : SolidPointElement() {}
 };


} //End of namespace extension



/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////


//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    
  /// Doc info
  DocInfo Doc_info;
  
  /// Reynolds number
  double Re=5.0;

  /// Strouhal number
  double St = 1.0;

  /// Womersley number (Reynolds x Strouhal)
  double ReSt = 5.0;
  
  /// Product of Reynolds number and inverse of Froude number
  double ReInvFr = 5.0;
  
  /// Ratio of viscosity in upper fluid to viscosity in lower
  /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
  double Viscosity_Ratio = 0.1;
  
  /// Ratio of density in upper fluid to density in lower
  /// fluid. Reynolds number etc. is based on density in lower fluid.
  double Density_Ratio = 0.5;
  
  /// Capillary number
  double Ca = 0.01; 

  /// Direction of gravity
  Vector<double> G(2);
  
  /// Pseudo-solid Poisson ratio
  double Nu=0.1;

  /// Length of the channel
  double Length = 1.0;

  // Height of the channel
  double Height = 2.0;
  
  // Relative Height of the interface
  double Interface_Height = 0.5;
  
  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=0;

  /// Trace file
  ofstream Trace_file;

  /// File to document the norm of the solution (for validation 
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

 } // end_of_namespace



//==start_of_problem_class============================================
/// Problem class to simulate viscous drop propagating along 2D channel
//====================================================================
template<class ELEMENT>
class TwoLayerInterfaceProblem : public Problem
{

public:

 /// Constructor
 TwoLayerInterfaceProblem();
 
 /// Destructor
 ~TwoLayerInterfaceProblem()
  {
   // Fluid timestepper
   delete this->time_stepper_pt(0);

   // Kill data associated with outer boundary
   unsigned n=Outer_boundary_polyline_pt->npolyline();
   for (unsigned j=0;j<n;j++)
    {
     delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
   delete Outer_boundary_polyline_pt;
   
   // Flush element of free surface elements
   delete_free_surface_elements();
   delete Free_surface_mesh_pt;

   // Delete error estimator
   delete Fluid_mesh_pt->spatial_error_estimator_pt();

   // Delete fluid mesh
   delete Fluid_mesh_pt;

   // Delete the global pressure drop data
   delete Drop_pressure_data_pt;

   // Kill const eqn
   delete Problem_Parameter::Constitutive_law_pt;

  }
 

 /// Actions before adapt: Wipe the mesh of free surface elements
 void actions_before_adapt()
  {
   // Kill the  elements and wipe surface mesh
   delete_free_surface_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  
  }// end of actions_before_adapt

 
 /// Actions after adapt: Rebuild the mesh of free surface elements
 void actions_after_adapt()
  {
   // Create the elements that impose the displacement constraint 
   create_free_surface_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
   
   // Setup the problem again -- remember that fluid mesh has been
   // completely rebuilt and its element's don't have any
   // pointers to Re etc. yet
   complete_problem_setup();
   
  }// end of actions_after_adapt

 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve 
 void actions_before_newton_solve()
  {
   //Reset the Lagrangian coordinates of the nodes to be the current
   //Eulerian coordinates -- this makes the current configuration
   //stress free
   Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  }
 
 /// Set boundary conditions and complete the build of all elements
 void complete_problem_setup();

 /// Doc the solution
 void doc_solution(const std::string& comment="");
 
 /// Compute the error estimates and assign to elements for plotting
 void compute_error_estimate(double& max_err,
                             double& min_err);

/// Set the initial conditions
void set_initial_condition()
{
 // Determine number of nodes in mesh
 const unsigned n_node = Fluid_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Loop over the two velocity components
   for(unsigned i=0;i<2;i++)
    {
     // Set velocity component i of node n to zero
     Fluid_mesh_pt->node_pt(n)->set_value(i,0.0);
    }
  }
 
 // Initialise the previous velocity values for timestepping
 // corresponding to an impulsive start
 assign_initial_values_impulsive();
 
} // End of set_initial_condition



 
 /// Function to deform the interface
void deform_interface(const double &epsilon,const unsigned &n_periods)
  {
   // Determine number of nodes in the "bulk" mesh
   const unsigned n_node = Fluid_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Determine eulerian position of node
   const double current_x_pos = Fluid_mesh_pt->node_pt(n)->x(0);
   const double current_y_pos = Fluid_mesh_pt->node_pt(n)->x(1);
   
   // Determine new vertical position of node
   const double new_y_pos = current_y_pos
    + (1.0-fabs(1.0-current_y_pos))*epsilon
    *(cos(2.0*n_periods*MathematicalConstants::Pi*current_x_pos/Problem_Parameter::Length));
   
   // Set new position
   Fluid_mesh_pt->node_pt(n)->x(1) = new_y_pos;
  }
} // End of deform_interface


 
 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e,
                   const unsigned &pdof, 
                   const double &pvalue)
  {
   // Fix the pressure at that element
   dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }


private:
 

 /// Create free surface elements
 void create_free_surface_elements();

 /// Delete free surface elements 
 void delete_free_surface_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Free_surface_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Free_surface_mesh_pt->element_pt(e);
    }
  

   // Wipe the mesh
   Free_surface_mesh_pt->flush_element_and_node_storage();
   

  } // end of delete_free_surface_elements
 
 
 /// Pointers to mesh of free surface elements
 Mesh* Free_surface_mesh_pt;
 

 /// Pointer to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;

 /// Vector storing pointer to the drop polygons
 Vector<TriangleMeshPolygon*> Drop_polygon_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 

 /// Pointer to a global drop pressure datum
 Data* Drop_pressure_data_pt;

 /// Backed up drop pressure between adaptations
 double Initial_value_for_drop_pressure;

 /// Pointer to hijacked element
 ELEMENT* Hijacked_element_pt;
 
 /// Bool to indicate if volume constraint is applied (only for steady run)
 bool Use_volume_constraint;

 /// Enumeration of channel boundaries
 enum 
 {
  Inflow_boundary_id=0,
  Upper_wall_boundary_id=1,
  Outflow_boundary_id=2,
  Bottom_wall_boundary_id=3,
  Interface_boundary_id=4,
 };
 

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
TwoLayerInterfaceProblem<ELEMENT>::TwoLayerInterfaceProblem()
{ 
 // Output directory
 Problem_Parameter::Doc_info.set_directory("RESLT");
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 this->add_time_stepper_pt(new BDF<2>);
 
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polylines
 //------------------------
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 
 // Each polyline only has three vertices -- provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord(3);
 for(unsigned i=0;i<3;i++)
  {
   vertex_coord[i].resize(2);
  }

 const double H = Problem_Parameter::Height;
 const double h = Problem_Parameter::Interface_Height*H;
 
 // First polyline: Inflow
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]= h;
 vertex_coord[2][0]=0.0;
 vertex_coord[2][1]= H;
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,
                                                   Inflow_boundary_id);
 
 // Second boundary polyline: Upper wall
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=H;
 vertex_coord[1][0]=0.5*Problem_Parameter::Length;
 vertex_coord[1][1]=H;
 vertex_coord[2][0]=Problem_Parameter::Length;
 vertex_coord[2][1]=H;

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,
                                                   Upper_wall_boundary_id);

 // Third boundary polyline: Outflow
 vertex_coord[0][0]=Problem_Parameter::Length;
 vertex_coord[0][1]=H;
 vertex_coord[1][0]=Problem_Parameter::Length;
 vertex_coord[1][1]=h;
 vertex_coord[2][0]=Problem_Parameter::Length;
 vertex_coord[2][1]=0.0;

 // Build the 3rd boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,
                                                   Outflow_boundary_id);

 // Fourth boundary polyline: Bottom wall
 vertex_coord[0][0]=Problem_Parameter::Length;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=0.5*Problem_Parameter::Length;
 vertex_coord[1][1]=0.0;
 vertex_coord[2][0]=0.0;
 vertex_coord[2][1]=0.0;

 // Build the 4th boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,
                                                    Bottom_wall_boundary_id);
 
 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);
 
 //Here we need to put the dividing internal line in
 Vector<TriangleMeshOpenCurve *> interface_pt(1);
 //Set the vertex coordinates
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=h;
 vertex_coord[1][0]=0.5*Problem_Parameter::Length;
 vertex_coord[1][1]=h;
 vertex_coord[2][0]=Problem_Parameter::Length;
 vertex_coord[2][1]=h;

//Create the internal line
  TriangleMeshPolyLine* interface_polyline_pt =
   new TriangleMeshPolyLine(vertex_coord,
                            Interface_boundary_id);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the inflow boundary
  interface_polyline_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[0]),1);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the outflow boundary
  interface_polyline_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[2]),1);

  Vector<TriangleMeshCurveSection*> interface_curve_pt(1);
  interface_curve_pt[0] = interface_polyline_pt;
  
  interface_pt[0] = new TriangleMeshOpenCurve(interface_curve_pt);
  
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------

 // Convert to "closed curve" objects
 TriangleMeshClosedCurve* outer_closed_curve_pt=Outer_boundary_polyline_pt;

 unsigned n_internal_closed_boundaries = 0;
 Vector<TriangleMeshClosedCurve *>
  inner_boundaries_pt(n_internal_closed_boundaries);
 
 // Target area for initial mesh
 double uniform_element_area=0.01;

 // Use the TriangleMeshParameter object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   outer_closed_curve_pt);

 //Define the inner boundaries
 triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundaries_pt;
 
 // Define the holes on the boundary
 triangle_mesh_parameters.internal_open_curves_pt() = interface_pt;

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 Vector<double> lower_region(2);
 lower_region[0] = 0.5*Problem_Parameter::Length;
 lower_region[1] = 0.5*Problem_Parameter::Interface_Height;
 
 // Define the region
 triangle_mesh_parameters.add_region_coordinates(1, lower_region);
 
 // Create the mesh
 Fluid_mesh_pt =
   new RefineableSolidTriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=0.005;
 Fluid_mesh_pt->min_permitted_error()=0.001;  
 Fluid_mesh_pt->max_element_size()=0.2;
 Fluid_mesh_pt->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   Fluid_mesh_pt->min_element_size()=0.01; 
  }

 // Output boundary and mesh initial mesh for information
 this->Fluid_mesh_pt->output_boundaries("boundaries.dat");
 this->Fluid_mesh_pt->output("mesh.dat");

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();
 
 // Construct the mesh of free surface elements
 Free_surface_mesh_pt=new Mesh;
 create_free_surface_elements();

 // Combine meshes
 //---------------

 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Free_surface sub meshes
 this->add_sub_mesh(this->Free_surface_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor


//============start_of_create_free_surface_elements===============
/// Create elements that impose the kinematic and dynamic bcs
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void TwoLayerInterfaceProblem<ELEMENT>::create_free_surface_elements()
{ 

 //Loop over the free surface boundaries
 unsigned nb=Fluid_mesh_pt->nboundary();
 for(unsigned b=4;b<nb;b++)
  {
   // Note: region is important
   // How many bulk fluid elements are adjacent to boundary b in region 0?
   unsigned n_element = Fluid_mesh_pt->nboundary_element_in_region(b,0);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b in region 0
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Fluid_mesh_pt->boundary_element_in_region_pt(b,0,e));
     
     //Find the index of the face of element e along boundary b in region 0
     int face_index = Fluid_mesh_pt->face_index_at_boundary_in_region(b,0,e);
     
     // Create new element
     ElasticLineFluidInterfaceElement<ELEMENT>* el_pt =
      new ElasticLineFluidInterfaceElement<ELEMENT>(
       bulk_elem_pt,face_index);   
     
     // Add it to the mesh
     Free_surface_mesh_pt->add_element_pt(el_pt);
     
     //Add the appropriate boundary number
     el_pt->set_boundary_number_in_bulk_mesh(b);

     //Specify the Strouhal number
     el_pt->st_pt() = &Problem_Parameter::St;
     
     //Specify the capillary number
     el_pt->ca_pt() = &Problem_Parameter::Ca;
    }
  }
}
// end of create_free_surface_elements




//==start_of_complete_problem_setup=======================================
/// Set boundary conditions and complete the build of all elements
//========================================================================
template<class ELEMENT>
void TwoLayerInterfaceProblem<ELEMENT>::complete_problem_setup()
  {      
   // Re-set the boundary conditions for fluid problem: All nodes are
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
      
       //Pin both velocities on side boundaries (1 and 3)
       if((ibound==1) || (ibound==3))
        {
         nod_pt->pin(0);
         nod_pt->pin(1);
        }
      
       //If it's the inflow or outflow pin only the horizontal velocity
       if((ibound==0) || (ibound==2)) {nod_pt->pin(0);}
       
       // Pin horizontal pseudo-solid positions
       SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
       solid_node_pt->pin_position(0);

       //Pin the vertical position on the upper and lower walls
       if((ibound==1) || (ibound==3))
        {
         solid_node_pt->pin_position(1);
        }
       else
        {
         solid_node_pt->unpin_position(1);
        }
      } //End of loop over nodes
    } // end loop over boundaries
   
   // Complete the build of all elements so they are fully functional
   // Remember that adaptation for triangle meshes involves a complete
   // regneration of the mesh (rather than splitting as in tree-based
   // meshes where such parameters can be passed down from the father
   // element!)
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     
     // Set the Reynolds number
     el_pt->re_pt() = &Problem_Parameter::Re;

     // Set the Womersley number (same as Re since St=1)
     el_pt->re_st_pt() = &Problem_Parameter::ReSt;
     
     // Set the product of the Reynolds number and the inverse of the
     // Froude number
     el_pt->re_invfr_pt() = &Problem_Parameter::ReInvFr;
     
     // Set the direction of gravity
     el_pt->g_pt() = &Problem_Parameter::G;
     
     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;
    }

   //For the elements in the upper region (region 0),
   //set the viscosity ratio
   n_element = Fluid_mesh_pt->nregion_element(1);
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = 
      dynamic_cast<ELEMENT*>(Fluid_mesh_pt->region_element_pt(1,e));
     
     el_pt->viscosity_ratio_pt() = &Problem_Parameter::Viscosity_Ratio;

     el_pt->density_ratio_pt() = &Problem_Parameter::Density_Ratio;
    }

   
   // Re-apply boundary values on Dirichlet boundary conditions 
   // (Boundary conditions are ignored when the solution is transferred
   // from the old to the new mesh by projection; this leads to a slight
   // change in the boundary values (which are, of course, never changed,
   // unlike the actual unknowns for which the projected values only
   // serve as an initial guess)

   // Set velocity and history values of velocity on walls
   nbound=this->Fluid_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<nbound;++ibound)
    {
     // Loop over nodes on this boundary
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();

       //If we are on the upper or lower walls
       // Velocity is and was zero at all previous times
       if((ibound==1) || (ibound==3))
        {
         //Loop over time history values
         for (unsigned t=0;t<=n_prev;t++)
          {
           nod_pt->set_value(t,0,0.0); 
           nod_pt->set_value(t,1,0.0);
           
           // Nodes have always been there...
           nod_pt->x(t,0)=nod_pt->x(0,0);
           nod_pt->x(t,1)=nod_pt->x(0,1);
          }
        }
       //If we are on the side wals there is no horizontal velocity or
       //change in horizontal position
       if((ibound==0) || (ibound==2))
        {
         for (unsigned t=0;t<=n_prev;t++)
          {
           nod_pt->set_value(t,0,0.0); 
           // Nodes have always been there...
           nod_pt->x(t,0)=nod_pt->x(0,0);
          }
        }
       //But there is a change in vertical position!
      }
    } //End of loop over boundaries
   
   fix_pressure(0,0,0.0);

  }



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TwoLayerInterfaceProblem<ELEMENT>::doc_solution(const std::string& comment)
{

 // Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 double min_x_coordinate = 2.0;
 unsigned min_boundary_node = 0;
 //Find the left-most node on the boundary
 unsigned n_bound = Fluid_mesh_pt->nboundary_node(4);
 for(unsigned n=0;n<n_bound;++n)
  {
   Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(4,n);
   if(nod_pt->x(0) < min_x_coordinate)
    {
     min_x_coordinate = nod_pt->x(0);
     min_boundary_node = n;
    }
  }
 
 // Document time and vertical position of left hand side of interface
 // in trace file
 Problem_Parameter::Trace_file
  << time_pt()->time() << " "
  << Fluid_mesh_pt->boundary_node_pt(4,min_boundary_node)->x(1) << std::endl;
 
 ofstream some_file;
 char filename[100];
 
 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;
 
 // Open solution output file
 sprintf(filename,"%s/soln%i.dat",
         Problem_Parameter::Doc_info.directory().c_str(),
         Problem_Parameter::Doc_info.number());
 some_file.open(filename);

 // Output solution to file
 Fluid_mesh_pt->output(some_file,npts);

 // Close solution output file
 some_file.close();

 // Open interface solution output file
 sprintf(filename,"%s/interface_soln%i.dat",
         Problem_Parameter::Doc_info.directory().c_str(),
         Problem_Parameter::Doc_info.number());
 some_file.open(filename);
 
 // Output solution to file
 Free_surface_mesh_pt->output(some_file,npts);
 
 // Close solution output file
 some_file.close();

 // Increment the doc_info number
 Problem_Parameter::Doc_info.number()++;
}

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void TwoLayerInterfaceProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                             double& min_err)
{ 
 // Get error estimator
 ErrorEstimator* err_est_pt=Fluid_mesh_pt->spatial_error_estimator_pt();
 
 // Get/output error estimates
 unsigned nel=Fluid_mesh_pt->nelement();
 Vector<double> elemental_error(nel);
 
 // We need a dynamic cast, get_element_errors from the Fluid_mesh_pt
 // Dynamic cast is used because get_element_errors require a Mesh* ans
 // not a SolidMesh*
 Mesh* fluid_mesh_pt=dynamic_cast<Mesh*>(Fluid_mesh_pt);
 err_est_pt->get_element_errors(fluid_mesh_pt,
                                elemental_error);

 // Set errors for post-processing and find extrema
 max_err=0.0;
 min_err=DBL_MAX;
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<MyCrouzeixRaviartElement*>(Fluid_mesh_pt->element_pt(e))->
    set_error(elemental_error[e]);

   max_err=std::max(max_err,elemental_error[e]);
   min_err=std::min(min_err,elemental_error[e]);
  }
  
}


//============================================================
/// Driver code for moving block problem
//============================================================
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

 // Create generalised Hookean constitutive equations
 Problem_Parameter::Constitutive_law_pt = 
  new GeneralisedHookean(&Problem_Parameter::Nu);

 Problem_Parameter::G[0] = 0.0;
 Problem_Parameter::G[1] = -1.0;
 
 // Create problem in initial configuration
 TwoLayerInterfaceProblem<Hijacked<ProjectableCrouzeixRaviartElement<
 MyCrouzeixRaviartElement> > > problem;  

 // Set value of epsilon
 const double epsilon = 0.1;

 // Set number of periods for cosine term
 const unsigned n_periods = 1;

 const double dt = 0.0025;

 double t_max = 0.6;
 
 // Deform the mesh/free surface
 problem.deform_interface(epsilon,n_periods);
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Problem_Parameter::Doc_info.directory().c_str());
 Problem_Parameter::Trace_file.open(filename);

 // Initialise trace file
 Problem_Parameter::Trace_file << "time, interface height" << std::endl;

 // Initialise timestep
 problem.initialise_dt(dt);

 // Set initial condition
 problem.set_initial_condition();

 // Maximum number of spatial adaptations per timestep
 unsigned max_adapt = 2;

 // Doc initial solution
 problem.doc_solution();

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);

 // Are we on the first timestep? At this point, yes!
 bool first_timestep = true;

 // Timestepping loop
 for(unsigned t=1;t<=n_timestep;t++)
  {
   // Output current timestep to screen
   cout << "\nTimestep " << t << " of " << n_timestep << std::endl;
   
   // Take one fixed timestep with spatial adaptivity
   problem.unsteady_newton_solve(dt,max_adapt,first_timestep);

   // No longer on first timestep, so set first_timestep flag to false
   first_timestep = false; 

   // Reset maximum number of adaptations for all future timesteps
   max_adapt = 1;

   // Doc solution
   problem.doc_solution();

  } // End of timestepping loop


} //End of main

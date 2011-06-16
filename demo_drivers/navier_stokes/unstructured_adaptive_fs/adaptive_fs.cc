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
#include <fenv.h> 

//Generic routines
#include "generic.h"


// The equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Elements that impose the volume constraint for discretisations in which
// the fluid mesh is updated by pseudo-elasticity 
// hierher shouldn't they be moved elsewhere?
#include "fix_vol_int_elastic_elements.h"

using namespace std;
using namespace oomph;


namespace oomph
{
//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    
  /// Doc info
  DocInfo Doc_info;
  
  /// Reynolds number
  double Re=0.0;

  /// Capillary number
  double Ca = 10.0;

  /// Pseudo-solid Poisson ratio
  double Nu=0.3;

  /// Initial radius of bubble
  double Radius = 0.25;

  /// Volume of the interface
  double Volume = MathematicalConstants::Pi*Radius*Radius;

  /// \short Scaling factor for inflow velocity (allows it to be switched off
  /// to do hydrostatics)
  double Inflow_veloc_magnitude = 0.0;

  /// \short Length of the channel
  double Length = 3.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=0;

  /// Trace file
  ofstream Trace_file;

  /// \short File to document the norm of the solution (for validation 
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

 } // end_of_namespace
 

//==============================================================
/// Overload TaylorHood element to modify output
//==============================================================
 class MyTaylorHoodElement : 
  public virtual PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, 
  TPVDElement<2,3> >
 {
  
 private:
  
  /// Storage for elemental error estimate -- used for post-processing
  double Error;

 public:

  /// Constructor initialise error
  MyTaylorHoodElement()
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


  /// Get square of L2 norm of velocity 
  double square_of_l2_norm()
   {

    // Assign dimension 
    unsigned el_dim=2; 
    // Initalise
    double sum=0.0;
    
    //Find out how many nodes there are
    unsigned n_node = nnode();
    
    //Find the indices at which the local velocities are stored
    unsigned u_nodal_index[el_dim];
    for(unsigned i=0;i<el_dim;i++) {u_nodal_index[i] = u_index_nst(i);}
    
    //Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node,el_dim);
    
    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(el_dim);
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt,psif,dpsidx);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      //Calculate velocities 
      Vector<double> interpolated_u(el_dim,0.0);      
      
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
       {
        //Loop over directions
        for(unsigned i=0;i<el_dim;i++)
         {
          //Get the nodal value
          double u_value = raw_nodal_value(l,u_nodal_index[i]);
          interpolated_u[i] += u_value*psif[l];
         }
       }

      //Assemble square of L2 norm
      for(unsigned i=0;i<el_dim;i++)
       {
        sum+=interpolated_u[i]*interpolated_u[i]*W;
       }           
     }

    return sum;

   }

 };


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<MyTaylorHoodElement>
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
 class FaceGeometry<FaceGeometry<MyTaylorHoodElement> >
  : public virtual SolidPointElement 
 {
 public:
  FaceGeometry() : SolidPointElement() {}
 };


} //End of namespace extension



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Problem class to simulate inviscid bubble propagating along 2D channel
//====================================================================
template<class ELEMENT>
class BubbleInChannelProblem : public Problem
{

public:

 /// Constructor
 BubbleInChannelProblem();
 
 /// Destructor
 ~BubbleInChannelProblem()
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

   //Kill data associated with bubbles
   unsigned n_bubble = Bubble_polygon_pt.size();
   for(unsigned ibubble=0;ibubble<n_bubble;ibubble++)
    {
     unsigned n=Bubble_polygon_pt[ibubble]->npolyline();
     for (unsigned j=0;j<n;j++)
      {
       delete Bubble_polygon_pt[ibubble]->polyline_pt(j);
      }
     delete Bubble_polygon_pt[ibubble];
    }
   
   // Flush element of free surface elements
   delete_free_surface_elements();
   delete Free_surface_mesh_pt;

   // Delete error estimator
   delete Fluid_mesh_pt->spatial_error_estimator_pt();

   // Delete fluid mesh
   delete Fluid_mesh_pt;

   // Delete the volume constraint mesh
   delete Volume_constraint_mesh_pt;

   // Delete the global pressure bubble data
   delete Bubble_pressure_data_pt;

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

   // Output solution after adaptation/projection
   //doc_solution("new mesh with projected solution");
   
  }// end of actions_after_adapt

 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
 
 
 /// \short Set boundary conditions and complete the build of all elements
 void complete_problem_setup()
  {      
   // Map to record if a given boundary is on a bubble or not
   map<unsigned,bool> is_on_bubble_bound;
   
   // Loop over the bubbles 
   unsigned nbubble=Bubble_polygon_pt.size();
   for(unsigned ibubble=0;ibubble<nbubble;ibubble++)
    {
     // Get the vector all boundary IDs associated with the polylines that
     // make up the closed polygon
     Vector<unsigned> bubble_bound_id=this->Bubble_polygon_pt[ibubble]->
      polygon_boundary_id();
     
     // Get the number of boundary
     unsigned nbound=bubble_bound_id.size();
     
     // Fill in the map
     for(unsigned ibound=0;ibound<nbound;ibound++)
      {
       // This boundary...
       unsigned bound_id=bubble_bound_id[ibound];
       
       // ...is on the bubble
       is_on_bubble_bound[bound_id]=true;
      }
    }
   
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
       
       //Pin both velocities on inflow (0) and side boundaries (1 and 3)
       if((ibound==0) || (ibound==1) || (ibound==3))
        {
         nod_pt->pin(0);
         nod_pt->pin(1);
        }
       
       //If it's the outflow pin only the vertical velocity
       if(ibound==2) {nod_pt->pin(1);}
       
       // Pin pseudo-solid positions apart from bubble boundary which
       // we allow to move
       SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
       if(is_on_bubble_bound[ibound])
        {
         solid_node_pt->unpin_position(0);
         solid_node_pt->unpin_position(1);
        }
       else
        {
         solid_node_pt->pin_position(0);
         solid_node_pt->pin_position(1);
        }
      }
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
     
     // Set pointer to continous time
     el_pt->time_pt()=this->time_pt();
     
     // Set the Reynolds number
     el_pt->re_pt() = &Problem_Parameter::Re;

     // Set the Womersley number (same as Re since St=1)
     el_pt->re_st_pt() = &Problem_Parameter::Re;
     
     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;
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
     if ((ibound==Upper_wall_boundary_id)||
         (ibound==Bottom_wall_boundary_id)||
         (ibound==Outflow_boundary_id))
      {
       // Loop over nodes on this boundary
       unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
       for (unsigned inod=0;inod<num_nod;inod++)
        {
         // Get node
         Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
         
         // Get number of previous (history) values
         unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
         
         // Velocity is and was zero at all previous times
         for (unsigned t=0;t<=n_prev;t++)
          {
           // Parallel outflow
           if (ibound!=Outflow_boundary_id)
            {
             nod_pt->set_value(t,0,0.0); 
            }
           nod_pt->set_value(t,1,0.0);
          }
        }
      }
    }

   // Re-assign prescribed inflow velocity at inlet
   unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(Inflow_boundary_id);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(Inflow_boundary_id,
                                                        inod);
     
     //Now set the boundary velocity
     double y = nod_pt->x(1);
     nod_pt->set_value(0,Problem_Parameter::Inflow_veloc_magnitude*y*(1-y));
    }
  }
     
 /// Doc the solution
 void doc_solution(const std::string& comment="");
 
 /// Compute the error estimates and assign to elements for plotting
 void compute_error_estimate(double& max_err,
                             double& min_err);
 
private:
 

 /// \short Create free surface elements
 void create_free_surface_elements();

 /// \short Delete free surface elements 
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
 
 /// Pointer to mesh containing single element that imposes volume constraint
 Mesh* Volume_constraint_mesh_pt;

 /// Pointer to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
 
 /// Vector storing pointer to the bubble polygons
 Vector<TriangleMeshInternalPolygon*> Bubble_polygon_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 
 
 /// Pointer to a global bubble pressure datum
 Data* Bubble_pressure_data_pt;

 /// Enumeration of channel boundaries
 enum 
 {
  Inflow_boundary_id=0,
  Upper_wall_boundary_id=1,
  Outflow_boundary_id=2,
  Bottom_wall_boundary_id=3
 };
 

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
BubbleInChannelProblem<ELEMENT>::BubbleInChannelProblem()
{ 
 // Output directory
 Problem_Parameter::Doc_info.set_directory("RESLT");
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 this->add_time_stepper_pt(new BDF<2>);

  // Create bubble pressure as global Data
 Bubble_pressure_data_pt =  new Data(1);
 this->add_global_data(Bubble_pressure_data_pt);
 
 //Provide a reasonable initial guess (hydrostatics)
 Bubble_pressure_data_pt->set_value(0,
                                    Problem_Parameter::Ca/
                                    Problem_Parameter::Radius);
 

 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polylines
 //------------------------
 Vector<TriangleMeshPolyLine*> boundary_polyline_pt(4);
 
 // Each polyline only has two vertices -- provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord(2);
 for(unsigned i=0;i<2;i++)
  {
   vertex_coord[i].resize(2);
  }
 
 // First polyline: Inflow
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=1.0;
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,
                                                   Inflow_boundary_id);
 
 // Second boundary polyline: Upper wall
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=1.0;
 vertex_coord[1][0]=Problem_Parameter::Length;
 vertex_coord[1][1]=1.0;

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,
                                                   Upper_wall_boundary_id);

 // Third boundary polyline: Outflow
 vertex_coord[0][0]=Problem_Parameter::Length;
 vertex_coord[0][1]=1.0;
 vertex_coord[1][0]=Problem_Parameter::Length;
 vertex_coord[1][1]=0.0;

 // Build the 3rd boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,
                                                   Outflow_boundary_id);

 // Fourth boundary polyline: Bottom wall
 vertex_coord[0][0]=Problem_Parameter::Length;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=0.0;

 // Build the 4th boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,
                                                    Bottom_wall_boundary_id);
 
 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);
 

 // Now define initial shape of bubble(s) with polygon
 //---------------------------------------------------

 // Counter for number of boundaries (bit over the top but future
 // proof for multiple bubbles...)
 unsigned boundary_id=Bottom_wall_boundary_id+1;

 // We have one bubble
 Bubble_polygon_pt.resize(1);

 // Place it smack in the middle of the channel
 double x_center = 0.5*Problem_Parameter::Length;
 double y_center = 0.5;
 Ellipse * bubble_pt = new Ellipse(Problem_Parameter::Radius,
                                       Problem_Parameter::Radius);
 
 // Intrinsic coordinate along GeomObject defining the bubble
 Vector<double> zeta(1);
 
 // Position vector to GeomObject defining the bubble
 Vector<double> coord(2);
 
 // Number of points defining bubble
 unsigned npoints = 16; 
 double unit_zeta = MathematicalConstants::Pi/double(npoints-1);
 
 // This bubble is bounded by two distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshPolyLine*> bubble_polyline_pt(2);
 
 // Vertex coordinates
 Vector<Vector<double> > bubble_vertex(npoints);
 
 // Create points on boundary
 for(unsigned ipoint=0; ipoint<npoints;ipoint++)
  {
   // Resize the vector 
   bubble_vertex[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=unit_zeta*double(ipoint);
   bubble_pt->position(zeta,coord);

   // Shift
   bubble_vertex[ipoint][0]=coord[0]+x_center;
   bubble_vertex[ipoint][1]=coord[1]+y_center;
  }
 
 // Build the 1st bubble polyline
 bubble_polyline_pt[0] = new TriangleMeshPolyLine(bubble_vertex,boundary_id++);

 // Second boundary of bubble
 for(unsigned ipoint=0; ipoint<npoints;ipoint++)
  {
   // Resize the vector 
   bubble_vertex[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=(unit_zeta*double(ipoint))+MathematicalConstants::Pi;
   bubble_pt->position(zeta,coord);

   // Shift
   bubble_vertex[ipoint][0]=coord[0]+x_center;
   bubble_vertex[ipoint][1]=coord[1]+y_center;
  }

 // Build the 2nd bubble polyline
 bubble_polyline_pt[1] = new TriangleMeshPolyLine(bubble_vertex,boundary_id++);


 // Define coordinates of a point inside the bubble
 Vector<double> bubble_center(2);
 bubble_center[0]=x_center;
 bubble_center[1]=y_center;
 
 
 // Create closed polygon from two polylines
 Bubble_polygon_pt[0] = new TriangleMeshInternalPolygon(
  bubble_center,bubble_polyline_pt);
 

 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------

 // Convert to "closed curve" objects
 TriangleMeshClosedCurve* outer_closed_curve_pt=Outer_boundary_polyline_pt;
 unsigned nb=Bubble_polygon_pt.size();
 Vector<TriangleMeshInternalClosedCurve*> bubble_closed_curve_pt(nb);
 for (unsigned i=0;i<nb;i++)
  {
   bubble_closed_curve_pt[i]=Bubble_polygon_pt[i];
  }

 // Target area for initial mesh
 double uniform_element_area=0.2;
 Fluid_mesh_pt = 
  new RefineableSolidTriangleMesh<ELEMENT>(outer_closed_curve_pt, 
                                           bubble_closed_curve_pt,
                                           uniform_element_area,
                                           this->time_stepper_pt());
 
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

 // Set the problem pointer (required for mesh adaptation/projection)
 Fluid_mesh_pt->problem_pt()=this;
   
 // Output boundary and mesh initial mesh for information
 this->Fluid_mesh_pt->output_boundaries("boundaries.dat");
 this->Fluid_mesh_pt->output("mesh.dat");
 
 // Set boundary condition and complete the build of all elements
 complete_problem_setup();
 
 // Construct the mesh of free surface elements
 Free_surface_mesh_pt=new Mesh;
 create_free_surface_elements();

 //Create volume constraint element and store in its own mesh
 Volume_constraint_mesh_pt = new Mesh;
 Volume_constraint_mesh_pt->add_element_pt(new ElasticVolumeConstraintElement);

 // Specify constrained pressure hierher
 dynamic_cast<ElasticVolumeConstraintElement*>(
  Volume_constraint_mesh_pt->element_pt(0))
  ->set_traded_pressure_data(this->global_data_pt(0));

 // Specify target volume for bubble
 dynamic_cast<ElasticVolumeConstraintElement*>(
  Volume_constraint_mesh_pt->element_pt(0))
  ->volume_pt() = &Problem_Parameter::Volume;


 // Combine meshes
 //---------------
 
 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Free_surface sub meshes
 this->add_sub_mesh(this->Free_surface_mesh_pt);
 
 // Add volume constraint sub mesh
 this->add_sub_mesh(this->Volume_constraint_mesh_pt);

 // Build global mesh
 this->build_global_mesh();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor


//============start_of_create_free_surface_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_free_surface_elements()
{ 
 //Loop over the free surface boundaries
 unsigned nb=Fluid_mesh_pt->nboundary();
 for(unsigned b=Bottom_wall_boundary_id+1;b<nb;b++)
  {
   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Fluid_mesh_pt->boundary_element_pt(b,e));
     
     //Find the index of the face of element e along boundary b
     int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     FixedVolumeElasticLineFluidInterfaceElement<ELEMENT>* el_pt =
      new FixedVolumeElasticLineFluidInterfaceElement<ELEMENT>(
       bulk_elem_pt,face_index);   
     
     // Add it to the mesh
     Free_surface_mesh_pt->add_element_pt(el_pt);
     
     //Add the appropriate boundary number
     el_pt->set_boundary_number_in_bulk_mesh(b);
     
     //Specify the capillary number
     el_pt->ca_pt() = &Problem_Parameter::Ca;

     //Specify the external pressure hierher ???
     el_pt->set_external_pressure_data(this->global_data_pt(0));

     //Set the traded pressure
     el_pt->set_traded_pressure_data(this->global_data_pt(0));

    } 
  }
}
// end of create_free_surface_elements


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::doc_solution(const std::string& comment)
{ 

 oomph_info << "Docing step: " << Problem_Parameter::Doc_info.number()
            << std::endl;

 ofstream some_file;
 char filename[100];
 sprintf(filename,"%s/soln%i.dat",
         Problem_Parameter::Doc_info.directory().c_str(),
         Problem_Parameter::Doc_info.number());

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Compute errors and assign to each element for plotting
 double max_err;
 double min_err;
 compute_error_estimate(max_err,min_err);
 
 // Assemble square of L2 norm 
 double square_of_l2_norm=0.0;
 unsigned nel=Fluid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   square_of_l2_norm+=
    dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))->
    square_of_l2_norm();
  }
 Problem_Parameter::Norm_file << sqrt(square_of_l2_norm) << std::endl;
 

 some_file.open(filename);
 some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
  ->variable_identifier();
 this->Fluid_mesh_pt->output(some_file,npts);   
 some_file << "TEXT X = 25, Y = 78, CS=FRAME T = \"Global Step " 
           << Problem_Parameter::Doc_info.number() << "  " 
           << comment << "\"\n";
 some_file.close();

 // Get max/min area
 double max_area;
 double min_area;
 Fluid_mesh_pt->max_and_min_area(max_area, min_area);

 // Write trace file
 Problem_Parameter::Trace_file 
  << this->time_pt()->time() << " " 
  << Fluid_mesh_pt->nelement() << " "
  << max_area << " "
  << min_area << " "
  << max_err << " "
  << min_err << " "
  << sqrt(square_of_l2_norm) << " "
  << std::endl;

 // Increment the doc_info number
 Problem_Parameter::Doc_info.number()++;


}

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::compute_error_estimate(double& max_err,
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
 err_est_pt->get_element_errors(this->communicator_pt(),
                                fluid_mesh_pt,
                                elemental_error);

 // Set errors for post-processing and find extrema
 max_err=0.0;
 min_err=DBL_MAX;
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<MyTaylorHoodElement*>(Fluid_mesh_pt->element_pt(e))->
    set_error(elemental_error[e]);

   max_err=std::max(max_err,elemental_error[e]);
   min_err=std::min(min_err,elemental_error[e]);
  }
  
}


//============================================================
///Driver code for moving block problem
//============================================================
int main(int argc, char **argv)
{
 
 // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
 
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
 
 // Open trace file
 Problem_Parameter::Trace_file.open("RESLT/trace.dat");
 
 // Open norm file
 Problem_Parameter::Norm_file.open("RESLT/norm.dat");
 

 // Create problem in initial configuration
 BubbleInChannelProblem<ProjectableTaylorHoodElement<MyTaylorHoodElement> > 
  problem;  


 // Output the problem's state with the bubble in its
 // initial polygonal representation
 problem.doc_solution();
 
 // Before starting the time-integration we want to "inflate" it to form 
 // a proper circular bubble. We do this by setting the inflow to zero
 // and doing a steady solve (with one adaptation)
 Problem_Parameter::Inflow_veloc_magnitude=0.0;
 problem.steady_newton_solve(1);

 // If all went well, this should show us a nice circular bubble
 // in a stationary fluid
 problem.doc_solution();

  // Initialise timestepper
 double dt=0.025;
 problem.initialise_dt(dt);
 
 // Perform impulsive start from current state
 problem.assign_initial_values_impulsive();

 // Now switch on the inflow and re-assign the boundary conditions
 // (Call to complete_problem_setup() is a bit expensive given that we
 // we only want to set the inflow velocity but who cares -- it's just
 // a one off.
 Problem_Parameter::Inflow_veloc_magnitude=1.0;
 problem.complete_problem_setup();
 
 // Solve problem on fixed mesh
 unsigned nstep=2;
 for (unsigned i=0;i<nstep;i++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt);    
   problem.doc_solution();
  }

 // Now do a proper loop, doing nstep timesteps before adapting/remeshing
 // and repeating the lot ncycle times
 unsigned ncycle=1000;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   ncycle=1;
   oomph_info << "Only doing one cycle during validation\n";
  }

 // Do the cycles
 for(unsigned j=0;j<ncycle;j++)
  {       
   // Adapt
   problem.adapt();
   
   //Solve problem a few times
   for (unsigned i=0;i<nstep;i++)
    {     
     // Solve the problem
     problem.unsteady_newton_solve(dt); 
     
     // Build the label for doc
     std::stringstream label;
     label << "Adaptation " <<j << " Step "<< i;
     problem.doc_solution(label.str());
    }
  }


} //End of main

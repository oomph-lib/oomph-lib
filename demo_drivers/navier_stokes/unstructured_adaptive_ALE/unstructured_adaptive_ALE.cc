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
#include <fenv.h> 

//Generic routines
#include "generic.h"


// The equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "rigid_body.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;


namespace oomph
{
//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    
  /// Block velocity
  double Block_x_velocity=0.5;
  
  /// Block velocity
  double Block_y_velocity=0.3;

  /// Block velocity
  double Block_rotation_velocity=1.0;
 
  /// Vector of Data object containing the three values (x-displacement,
  /// y-displacement and rotation angle) for the polygons' centreline
  /// motion.
  Vector<Data*> Centre_displacement_data_pt;
  
  /// Doc info
  DocInfo Doc_info;
  
  /// Reynolds number
  double Re=50.0;

  /// Pseudo-solid Poisson ratio
  double Nu=0.3;

  ///  Pseudo solid "density" -- set to zero because we don't want
  /// inertia in the node update!
  double Lambda_sq=0.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=0;

  /// Trace file
  ofstream Trace_file;

  /// File to document the norm of the solution (for validation purposes)
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

} //End of namespace extension



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Unstructured Navier-Stokes ALE Problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem : public Problem
{

public:

 /// Constructor
 UnstructuredFluidProblem();
 
 /// Destructor
 ~UnstructuredFluidProblem()
  {
   // Fluid timestepper
   delete this->time_stepper_pt(0);
   
   // Kill hole displacement data
   for(unsigned ihole=0;ihole<2;ihole++)
    {
     delete Problem_Parameter::Centre_displacement_data_pt[ihole];
    }

   // Kill data associated with outer boundary
   unsigned n=Outer_boundary_polyline_pt->npolyline();
   for (unsigned j=0;j<n;j++)
    {
     delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
   delete Outer_boundary_polyline_pt;


   //Kill data associated with inner holes
   for(unsigned ihole=0;ihole<2;ihole++)
    {
     unsigned n=Inner_hole_pt[ihole]->npolyline();
     for (unsigned j=0;j<n;j++)
      {
       delete Inner_hole_pt[ihole]->polyline_pt(j);
      }
     delete Inner_hole_pt[ihole];
    }
   
   // Flush Lagrange multiplier mesh
   delete_lagrange_multiplier_elements();
   delete Lagrange_multiplier_mesh_pt;

   // Delete error estimator
   delete Fluid_mesh_pt->spatial_error_estimator_pt();

   // Delete fluid mesh
   delete Fluid_mesh_pt;

   // Kill const eqn
   delete Problem_Parameter::Constitutive_law_pt;

  }
 
 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt()
  {
   // Kill the  elements and wipe surface mesh
   delete_lagrange_multiplier_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  
  }// end of actions_before_adapt

 
 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt()
  {
   //Reassign the lagrangian coordinates (an updated lagrangian approach)
   Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

   // Create the elements that impose the displacement constraint 
   create_lagrange_multiplier_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
   
   // Setup the problem again -- remember that fluid mesh has been
   // completely rebuilt and its element's don't have any
   // pointers to Re etc. yet
   complete_problem_setup();

   // Output solution after adaptation/projection
   bool doc_projection=true;
   doc_solution("new mesh with projected solution",doc_projection);
   
  }// end of actions_after_adapt


 ///  Update before implicit timestep: Move hole boundary
 void actions_before_implicit_timestep()
  {
   double time=this->time_pt()->time();
   double ampl=0.5*(exp(-10.0*time*time)-
                    cos(2.0*MathematicalConstants::Pi*time));

   // Assign the x increment of the holes using velocity and timestep 
   double dx=Problem_Parameter::
    Block_x_velocity*ampl*this->time_pt()->dt(); 
   
   // Assign the y increment of the hole using velocity and timestep 
   double dy=Problem_Parameter::
    Block_y_velocity*ampl*this->time_pt()->dt(); 

   // Assign the rotation of the hole using velocity and timestep 
   double ampl2=1.0-exp(-10.0*time*time);
   double drotation=Problem_Parameter::
    Block_rotation_velocity*ampl2*this->time_pt()->dt(); 
      
   // Update the geom object position   
   unsigned nhole=this->inner_hole_pt().size();
   for(unsigned ihole=0;ihole<nhole;ihole++)
    {
     //Shift the values for the centre displacement
     Problem_Parameter::Centre_displacement_data_pt[ihole]->time_stepper_pt()
      ->shift_time_values(
       Problem_Parameter::Centre_displacement_data_pt[ihole]);

     //Now backup
     double old_x=
      Problem_Parameter::Centre_displacement_data_pt[ihole]->value(0);     
     double old_y=
      Problem_Parameter::Centre_displacement_data_pt[ihole]->value(1);
     double old_angle=
      Problem_Parameter::Centre_displacement_data_pt[ihole]->value(2);

     
     // Update values 
     Problem_Parameter::Centre_displacement_data_pt[ihole]->
      set_value(1,old_y+dy);
     
     // Rotate one hole clockwise the other anti-clockwise
     if(ihole==0)
      {
       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(0,old_x-dx);

       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(1,old_y+dy);
       
       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(2,old_angle+drotation);
      }
     else
      {
       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(0,old_x+dx);

       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(1,old_y-dy);
       
       Problem_Parameter::Centre_displacement_data_pt[ihole]->
        set_value(2,old_angle-drotation);
      }
    }
  }
 
 ///  Re-apply the no slip condition (imposed indirectly via dependent
 /// velocities)
 void actions_before_newton_convergence_check()
  {
   // Update mesh -- this applies the auxiliary node update function
   Fluid_mesh_pt->node_update();
  }
 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
 
 
 ///  Set boundary condition, assign auxiliary node update fct.
 /// Complete the build of all elements, attach power elements that allow
 /// computation of drag vector
 void complete_problem_setup()
  {   
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
       
       // Pin everywhere apart from in/outflow (boundaries 0 and 2)
       // where we only impose parallel flow)
       if ((ibound!=0)&&(ibound!=2))
        {
         nod_pt->pin(0);
        }
       nod_pt->pin(1);
       
       // Pin pseudo-solid positions apart from hole boundary we want to move
       SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
       
       // Unpin the position of all the nodes on hole boundaries
       // since they will be moved using Lagrange Multiplier
       if(ibound > 3)
        {
         solid_node_pt->unpin_position(0);
         solid_node_pt->unpin_position(1);
         
         // Assign auxiliary node update fct if we're dealing with a 
         // hole boundary
         nod_pt->set_auxiliary_node_update_fct_pt(
          FSI_functions::apply_no_slip_on_moving_wall); 
        }
       else
        {
         solid_node_pt->pin_position(0);
         solid_node_pt->pin_position(1);
        }
      }

    } // end loop over boundaries
   
   // Complete the build of all elements so they are fully functional
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
     // Set the Reynolds number
     el_pt->re_pt() = &Problem_Parameter::Re;

     // Set the Womersley number (same as Re since St=1)
     el_pt->re_st_pt() = &Problem_Parameter::Re;

     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;

     // Set the "density" for pseudo-elastic mesh deformation
     el_pt->lambda_sq_pt()=&Problem_Parameter::Lambda_sq;
    }

   // Re-apply Dirichlet boundary conditions (projection ignores
   // boundary conditions!)
   
   // Zero velocity and history values of velocity on walls 
   // (boundaries 1 and 3)
   nbound=this->Fluid_mesh_pt->nboundary();
   for(unsigned ibound=1;ibound<4;ibound=ibound+2)
    {
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       
       // Zero all current and previous veloc values
       for (unsigned t=0;t<=n_prev;t++)
        {
         nod_pt->set_value(t,0,0.0);
         nod_pt->set_value(t,1,0.0);
        }

      }
    } 
  }
 
 /// Doc the solution
 void doc_solution(const std::string& comment="", const bool& project=false);
 
 /// Compute the error estimates and assign to elements for plotting
 void compute_error_estimate(double& max_err,
                             double& min_err);
  
 /// Sanity check: Doc boundary coordinates from mesh and GeomObject
 void doc_boundary_coordinates();
 
 /// Get the TriangleMeshPolygon objects
 Vector<TriangleMeshPolygon*>& inner_hole_pt()
  {return Inner_hole_pt;}
 
private:
 

 ///  Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 ///  Delete elements that impose the prescribed boundary displacement
 /// and wipe the associated mesh
 void delete_lagrange_multiplier_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Lagrange_multiplier_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();
   
  } // end of delete_lagrange_multiplier_elements
 
 /// Pointers to mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;
 
 /// Pointers to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
 
 /// Vector storing pointer to the hole polygon
 Vector<TriangleMeshPolygon*> Inner_hole_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 


}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor: build the first mesh with TriangleMeshPolygon and
///              TriangleMeshPolygon object
//========================================================================
template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
 // Allow for rough startup
 this->Problem::Max_residuals=1000.0;

 // Output directory
 Problem_Parameter::Doc_info.set_directory("RESLT");

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 this->add_time_stepper_pt(new BDF<2>);

 // Define the boundaries: Polyline with 4 different
 // boundaries for the outer boundary and 2 internal holes, 
 // egg shaped, with 2 boundaries each
 
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshCurveSection*> boundary_segment_pt(4);
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(2);
 for(unsigned i=0;i<2;i++)
  {
   bound_seg[i].resize(2);
  }
 
 // First boundary segment
 bound_seg[0][0]=0.0;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=0.0;
 bound_seg[1][1]=1.0;
 
 // Specify 1st boundary id
 unsigned bound_id = 0;

 // Build the 1st boundary segment
 boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Second boundary segment
 bound_seg[0][0]=0.0;
 bound_seg[0][1]=1.0;
 bound_seg[1][0]=2.0;
 bound_seg[1][1]=1.0;

 // Specify 2nd boundary id
 bound_id = 1;

 // Build the 2nd boundary segment
 boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Third boundary segment
 bound_seg[0][0]=2.0;
 bound_seg[0][1]=1.0;
 bound_seg[1][0]=2.0;
 bound_seg[1][1]=0.0;

 // Specify 3rd boundary id
 bound_id = 2;

 // Build the 3rd boundary segment
 boundary_segment_pt[2] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fourth boundary segment
 bound_seg[0][0]=2.0;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=0.0;
 bound_seg[1][1]=0.0;

 // Specify 4th boundary id
 bound_id = 3;

 // Build the 4th boundary segment
 boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_segment_pt);




 // Now deal with the moving holes
 //-------------------------------

 // We have two holes
 Inner_hole_pt.resize(2);
 
 // Create Data objects that contains the three values (x-displacement,
 // y-displacement and rotation angle) for the two polygons' centreline
 // motion.
 Problem_Parameter::Centre_displacement_data_pt.resize(2);
 for(unsigned ihole=0;ihole<2;ihole++)
  {
   // Create time-dependent Data with three values
   Problem_Parameter::Centre_displacement_data_pt[ihole]=
    new Data(this->time_stepper_pt(),3);
   
   // For now pin all three values
   Problem_Parameter::Centre_displacement_data_pt[ihole]->pin(0);
   Problem_Parameter::Centre_displacement_data_pt[ihole]->pin(1);
   Problem_Parameter::Centre_displacement_data_pt[ihole]->pin(2);
  }
 
 
 // Build first hole
 //-----------------
 double x_center = 1.3;
 double y_center = 0.5;
 double A = 0.1;
 double B = 0.25;
 Ellipse * egg_hole_pt = new Ellipse(A,B);
 
 // Define the vector of angle value to build the hole
 Vector<double> zeta(1);
 
 // Initialize the vector of coordinates
 Vector<double> coord(2);
 
 // Number of points defining hole
 unsigned ppoints = 8; 
 double unit_zeta = MathematicalConstants::Pi/double(ppoints-1);
 
 // This hole is bounded by two distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> hole_segment_pt(2);
 
 // Vertex coordinates
 Vector<Vector<double> > bound_hole(ppoints);
 
 // Create points on boundary
 for(unsigned ipoint=0; ipoint<ppoints;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=unit_zeta*double(ipoint);
   egg_hole_pt->position(zeta,coord);
   bound_hole[ipoint][0]=coord[0]+x_center;
   bound_hole[ipoint][1]=coord[1]+y_center;
  }
 
 // Inner hole center coordinates
 Vector<double> hole_center(2);
 hole_center[0]=x_center;
 hole_center[1]=y_center;
 
 // Specify the hole boundary id
 unsigned hole_id = 4;

 // Build the 1st hole polyline
 hole_segment_pt[0] = new TriangleMeshPolyLine(bound_hole,hole_id);

 // Second boundary of hole
 for(unsigned ipoint=0; ipoint<ppoints;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=(unit_zeta*double(ipoint))+MathematicalConstants::Pi;
   egg_hole_pt->position(zeta,coord);
   bound_hole[ipoint][0]=coord[0]+x_center;
   bound_hole[ipoint][1]=coord[1]+y_center;
  }

 // Specify the hole boundary id
 hole_id=5;

 // Build the 2nd hole polyline
 hole_segment_pt[1] = new TriangleMeshPolyLine(bound_hole,hole_id);

 
 // Fill in the vector of holes. Specify data that define centre's
 // displacement
 Inner_hole_pt[0] = new ImmersedRigidBodyTriangleMeshPolygon(
  hole_center,hole_segment_pt,this->time_stepper_pt(),
  Problem_Parameter::Centre_displacement_data_pt[0]);
 
 // Build the second hole
 //----------------------
 x_center = 0.5;
 y_center = 0.6;
 A = 0.1;
 B = 0.2;
 delete egg_hole_pt;
 egg_hole_pt = new Ellipse(A,B);
 
 // Number of points defining the hole
 ppoints = 8;
 unit_zeta = MathematicalConstants::Pi/double(ppoints-1);

 // Create points on the boundary
 for(unsigned ipoint=0; ipoint<ppoints;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=unit_zeta*double(ipoint);
   egg_hole_pt->position(zeta,coord);
   bound_hole[ipoint][0]=coord[0]+x_center;
   bound_hole[ipoint][1]=coord[1]+y_center;
  }
 
 // Specify the hole boundary id
 hole_id=6;
 
 // Build the 1st hole polyline
 hole_segment_pt[0] = new TriangleMeshPolyLine(bound_hole,hole_id);
 
 // Create points on second boundary
 for(unsigned ipoint=0; ipoint<ppoints;ipoint++)
  {
   // Resize the vector 
   bound_hole[ipoint].resize(2);
   
   // Get the coordinates
   zeta[0]=(unit_zeta*double(ipoint))+MathematicalConstants::Pi;
   egg_hole_pt->position(zeta,coord);
   bound_hole[ipoint][0]=coord[0]+x_center;
   bound_hole[ipoint][1]=coord[1]+y_center;
  }
 
 // Inner hole center coordinates
 hole_center[0]=x_center;
 hole_center[1]=y_center;
 
 // Clean up
 delete egg_hole_pt;

 // Specify the hole boundary id
 hole_id=7;

 // Build the 2nd hole polyline
 hole_segment_pt[1] = new TriangleMeshPolyLine(bound_hole,hole_id);
 
 // Fill in the second hole. Specify data that define centre's
 // displacement
 Inner_hole_pt[1] = new ImmersedRigidBodyTriangleMeshPolygon(
  hole_center,hole_segment_pt,this->time_stepper_pt(),
  Problem_Parameter::Centre_displacement_data_pt[1]);
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------
 double uniform_element_area=0.2;

 TriangleMeshClosedCurve* closed_curve_pt=Outer_boundary_polyline_pt;
 unsigned nh=Inner_hole_pt.size();
 Vector<TriangleMeshClosedCurve*> hole_pt(nh);
 for (unsigned i=0;i<nh;i++)
  {
   hole_pt[i]=Inner_hole_pt[i];
  }
 
 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   closed_curve_pt);

 // Define the holes on the domain
 triangle_mesh_parameters.internal_closed_curve_pt() =
   hole_pt;

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

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
   
 // Output boundary and mesh
 this->Fluid_mesh_pt->output_boundaries("boundaries.dat");
 this->Fluid_mesh_pt->output("mesh.dat");
   
 // Set boundary condition, assign auxiliary node update fct,
 // complete the build of all elements, attach power elements that allow
 // computation of drag vector
 complete_problem_setup();
 
 // Create Lagrange multiplier mesh for boundary motion
 //----------------------------------------------------
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();
 
 // Combine meshes
 //---------------
 
 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Lagrange_multiplier sub meshes
 this->add_sub_mesh(this->Lagrange_multiplier_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();
  
 // Sanity check: Doc boundary coordinates from mesh and GeomObject
 doc_boundary_coordinates();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor




//============start_doc_solid_zeta=======================================
/// Doc boundary coordinates in mesh and plot GeomObject representation
/// of inner boundary.
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_boundary_coordinates()
{

 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npoints = 5;
 
 // Output solution and projection files
 sprintf(filename,"RESLT/inner_hole_boundary_from_geom_obj.dat");
 some_file.open(filename);

 //Initialize zeta and r
 Vector<double> zeta(1);
 zeta[0]=0;
 
 Vector<double> r(2);
 r[0]=0;
 r[1]=0;
   
 // Loop over the boundaries of the fluid mesh and print out geometric
 // objects associated with each
 unsigned n_boundary = Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary;b++)
  {
   //Get geometric object for b-th bounday
   GeomObject* boundary_geom_obj_pt = 
    Fluid_mesh_pt->boundary_geom_object_pt(b);
   
   //Only bother to do anything if there is a geometric object
   if(boundary_geom_obj_pt!=0)
    {
     // Zone label 
     some_file << "ZONE T=boundary" << b << std::endl;
     
     //Get the coordinate limits
     Vector<double> zeta_limits = 
      Fluid_mesh_pt->boundary_coordinate_limits(b);
     
     //Set the increment
     double zeta_inc = 
      (zeta_limits[1] - zeta_limits[0])/(double)(npoints-1);
     
     for(unsigned i=0;i<npoints;i++)
      {
       // Get coordinate
       zeta[0] = zeta_limits[0] + i*zeta_inc;
       boundary_geom_obj_pt->position(zeta,r);
       
       // Print it
       some_file <<r[0]<<" "<<r[1]<<" "<<zeta[0]<<std::endl;  
      }
    }
  }
 some_file.close();
 
 // Doc boundary coordinates using Lagrange_multiplier_mesh_pt
 std::ofstream the_file("RESLT/inner_hole_boundary_from_mesh.dat");
 
 // Initialise max/min boundary coordinate
 double zeta_min= DBL_MAX;
 double zeta_max=-DBL_MAX;

 // Loop over Lagrange_multiplier elements
 unsigned n_face_element = this->Lagrange_multiplier_mesh_pt->nelement();
 
 for(unsigned e=0;e<n_face_element;e++)
  {
   
   //Cast the element pointer
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt=
    dynamic_cast< ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(e));
   
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
 
 
} //end doc_solid_zeta

//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_lagrange_multiplier_elements()
{
 // The idea is to apply Lagrange multipliers to the boundaries in the mesh
 // that have associated geometric objects
 
 //Find the number of boundaries
 unsigned n_boundary = Fluid_mesh_pt->nboundary();
 
 // Loop over the boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   //Get the geometric object associated with the boundary
   GeomObject* boundary_geom_obj_pt = 
    Fluid_mesh_pt->boundary_geom_object_pt(b);

   //Only bother to do anything if there is a geometric object
   if(boundary_geom_obj_pt!=0)
    {
     // How many bulk fluid elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk fluid elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk fluid element that is 
       // adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element. Note that we use different Lagrange
       // multiplier fields for each distinct boundary (here indicated
       // by b.
       ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt =
        new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
         bulk_elem_pt,face_index,b);   
       
       // Add it to the mesh
       Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);
       
       // Set the GeomObject that defines the boundary shape and set
       // which bulk boundary we are attached to (needed to extract
       // the boundary coordinate from the bulk nodes)
       el_pt->set_boundary_shape_geom_object_pt(
        boundary_geom_obj_pt,b);
       
       // Loop over the nodes to pin Lagrange multiplier
       unsigned nnod=el_pt->nnode();
       for(unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt = el_pt->node_pt(j);

         // How many nodal values were used by the "bulk" element
         // that originally created this node?
         unsigned n_bulk_value=el_pt->nbulk_value(j);
         
         // Pin two of the four Lagrange multipliers at vertices
         // This is not totally robust, but will work in this application
         unsigned nval=nod_pt->nvalue();
         if (nval==7)
          {
           for (unsigned i=0;i<2;i++) 
            { 
             // Pin lagrangian multipliers
             nod_pt->pin(n_bulk_value+2+i);
            }
          }
        }
      } // end loop over the element
    }  //End of case if there is a geometric object
  } //End of loop over boundaries
 
}
// end of create_lagrange_multiplier_elements


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_solution(
 const std::string& comment,
 const bool& project)
{ 

 oomph_info << "Docing step: " << Problem_Parameter::Doc_info.number()
            << std::endl;


 oomph_info 
  << "Docing step: " << Problem_Parameter::Doc_info.number() 
  << "\n cenrelinedispl = "
  <<Problem_Parameter::Centre_displacement_data_pt[0]->value(0)<< "\n"
  << "centre of mass: " 
  << dynamic_cast<ImmersedRigidBodyTriangleMeshPolygon*>(
   Inner_hole_pt[0])->initial_centre_of_mass(0) 
  << std::endl;


 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 


 // Compute errors and assign to each element for plotting
 double max_err;
 double min_err;
 compute_error_estimate(max_err,min_err);
 
 // Output solution and projection files
 if(!project)
  {
   sprintf(filename,"%s/soln%i.dat",
           Problem_Parameter::Doc_info.directory().c_str(),
           Problem_Parameter::Doc_info.number());
  }
 else
  {
   sprintf(filename,"%s/proj%i.dat",
           Problem_Parameter::Doc_info.directory().c_str(),
           Problem_Parameter::Doc_info.number()-1);
  }


 // Assemble square of L2 norm 
 double square_of_l2_norm=0.0;
 unsigned nel=Fluid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   square_of_l2_norm+=
    dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))->
    square_of_l2_norm();
  }
 Problem_Parameter::Norm_file << sqrt(square_of_l2_norm) << "\n";
 

 some_file.open(filename);
 some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
  ->variable_identifier();
 this->Fluid_mesh_pt->output(some_file,npts);   
 some_file << "TEXT X = 25, Y = 78, CS=FRAME T = \"Global Step " 
           << Problem_Parameter::Doc_info.number() << "  " 
           << comment << "\"\n";
 some_file.close();

 // No trace file writing after projection
 if(project) return;

 // Get max/min area
 double max_area;
 double min_area;
 Fluid_mesh_pt->max_and_min_element_size(max_area, min_area);

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
void UnstructuredFluidProblem<ELEMENT>::compute_error_estimate(double& max_err,
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
 UnstructuredFluidProblem<ProjectableTaylorHoodElement<MyTaylorHoodElement> > 
  problem;  
 
 // Initialise timestepper
 double dt=0.025;
 problem.initialise_dt(dt);
 
 // Perform impulsive start
 problem.assign_initial_values_impulsive();

 // Output initial conditions
 problem.doc_solution();

 // Solve problem a few times on given mesh
 unsigned nstep=3;
 for (unsigned i=0;i<nstep;i++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt);    
   problem.doc_solution();
  }

 // Now do a couple of adaptations
 unsigned ncycle=100;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   ncycle=1;
   oomph_info << "Only doing one cycle during validation\n";
  }
 for (unsigned j=0;j<ncycle;j++)
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
     label << "Cycle " <<j << " Step "<< i;
     problem.doc_solution(label.str());
    }
  }


} //End of main

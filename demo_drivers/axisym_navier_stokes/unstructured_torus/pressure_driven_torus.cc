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
//
//Flow in a torus computed by using a cylindrical polar coordinate system
//and assuming axisymmetry.

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

//=================================================================
//A namespace for the physical parameters in the problem
//=================================================================
namespace Global_Physical_Variables
{
 /// The Reynolds number
 double Re = 1.0;

 /// The radius of the torus
 double Radius = 1.0;

 /// The curvature of the torus
 double Delta = 0.1;

 /// The Dean number
 double Dean = 100.0;

 /// Pseudo-solid (mesh) Poisson ratio
 double Nu=0.3;
 
 /// Pseudo-solid (mesh) "density" 
 /// Set to zero because we don't want inertia in the node update!
 double Lambda_sq=0.0;
 
 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt=   
  new GeneralisedHookean(&Nu);

 /// Parameter that switches between circle (1.0) and square (0.0) 
 /// cross section
 double Mu = 0.0;


/// A function to specify a constant axial body force
void axial_pressure_gradient(const double &time,
                             const Vector<double> &x,
                             Vector<double> &result)
{
 using namespace Global_Physical_Variables;

 result[0] = 0.0;
 result[1] = 0.0;
 //Axial pressure gradient in siggers waters scaling is our original
 //Reynolds number
 result[2] = Dean/(x[0]*Delta*sqrt(2.0*Delta)); //1.0;
}

};

namespace oomph
{

//==============================================================
/// Overload TaylorHood element to modify output
//==============================================================
 class MyTaylorHoodElement : 
  public virtual PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement, 
  TPVDElement<2,3> >
 {
 public:

  /// Constructor initialise error
  MyTaylorHoodElement() {}
  
  /// Overload output function
  void output(std::ostream &outfile, 
              const unsigned &nplot)
   {
    
    // Assign dimension 
    unsigned el_dim=2;
    
    // Vector of local coordinates
    Vector<double> s(el_dim);
   
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
     
      // Coordinates
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_x(s,i) << " ";
       }

      for(unsigned i=0;i<el_dim;i++)
       {
        outfile << interpolated_xi(s,i) << " ";
       }
     
      // Velocities
      for(unsigned i=0;i<3;i++) 
       {
        outfile << interpolated_u_axi_nst(s,i) << " ";
       }
     
      // Pressure
      outfile << interpolated_p_axi_nst(s)  << " ";
     
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
        for(unsigned i=0;i<3;i++) 
         {
          outfile << interpolated_u_axi_nst(t,s,i) << " ";
         }
       }

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
 class FaceGeometry<MyTaylorHoodElement>
  : public virtual SolidTElement<1,3> 
 {
 public:
  FaceGeometry() : SolidTElement<1,3>() {}
 };

//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<FaceGeometry<MyTaylorHoodElement> >
  : public virtual PointElement 
 {
 public:
  FaceGeometry() : PointElement() {}
 };


}



//=start_of_FillCircle=============================================
//A geometric object that represents the geometry of the domain
//a circle of given centre and radius. One could use a non-linear
//stretch in r (xi[0]) to shift the elements towards the edge
//(boundary layer).
//=================================================================
class GeneralCircle : public GeomObject
{
public:

 /// Constructor that takes the centre position and radius of the circle
 /// as its arguments
 GeneralCircle(const double &centre_y,
               const double &radius) :
  GeomObject(1,2), Centre_y(centre_y), Radius(radius) { }
 
/// Destructor
virtual ~GeneralCircle(){}

/// Lagrangian coordinate xi
void position (const Vector<double>& xi, Vector<double>& r) const
{
 Vector<double> r_circle(2);

 r_circle[0] = 1.0/Global_Physical_Variables::Delta + Radius*cos(xi[0]);
 r_circle[1] = Centre_y + Radius*sin(xi[0]);


 Vector<double> r_square(2);
 const double pi = MathematicalConstants::Pi;
 if(std::fabs(xi[0]) < 0.25*pi)
  {
   r_square[0] = 1.0/Global_Physical_Variables::Delta + Radius;
   r_square[1] = Centre_y + Radius*tan(xi[0]);
  }
 else if((xi[0] >= 0.25*pi) && (xi[0] < 0.75*pi))
  {
   r_square[0] = 1.0/Global_Physical_Variables::Delta 
    + Radius*tan(0.5*pi - xi[0]);
   r_square[1] = Centre_y + Radius;
  }
 else if(std::fabs(xi[0]) >= 0.75*pi)
  {
   r_square[0] = 1.0/Global_Physical_Variables::Delta - Radius;
   r_square[1] = Centre_y + Radius*tan(pi - xi[0]);
  }
 else
  {
   r_square[0] = 1.0/Global_Physical_Variables::Delta 
    + Radius*tan(0.5*pi + xi[0]);
   r_square[1] = Centre_y - Radius;
  }


 //Now smoothly match between the two
 for(unsigned i=0;i<2;i++)
  {
   r[i] = Global_Physical_Variables::Mu*r_circle[i] + 
    (1.0 - Global_Physical_Variables::Mu)*r_square[i];
  }


}


 //Do not interpolated history values (Really this is only true when continuing)
 //This is going to be a total pain in general....
 bool interpolated_history(const unsigned &t)
  {
   if(t==0) {return true;}
   else {return false;}
  }

/// Return the position of the circle as a function of time 
/// (doesn't move as a function of time)
void position(const unsigned& t, 
              const Vector<double>& xi, Vector<double>& r) const
  {
   position(xi,r);
  }

 unsigned ngeom_data() const {return 0;}

 /// Return pointer to the j-th (only) Data item that the object's 
 /// shape depends on.
 Data* geom_data_pt(const unsigned& j) 
  {return 0;}
 


private:

 /// Storage for the y-coordinate of the centre
 double Centre_y;

 /// Storage for the radius of the circle
 double Radius;

};

//==========================================================================
/// Solve the Axisymmetric Navier Stokes equations in a torus
//==========================================================================
template<class ELEMENT>
class UnstructuredTorusProblem : public Problem
{
public:
 /// Constructor taking the maximum refinement level and
 /// the minimum and maximum error targets.
 UnstructuredTorusProblem(
              const double &min_error_target, 
              const double &max_error_target);

 /// Calculate the square of the l2 norm
 double calculate_square_of_l2_norm()
  {
   //Initialise
   double sum = 0.0;

   //Now loop over all elements and add the contributions to the 
   //components of the norm
   const unsigned n_element = this->Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     sum += dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))
      ->square_of_l2_norm();
    }
   return sum;
  }


 /// Calculate the cross-sectional area of the domain
 double calculate_area()
  {
   //Initialise
   double sum = 0.0;

   //Now loop over all elements and add the contributions to the 
   //components of the norm
   const unsigned n_element = this->Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     sum += this->Fluid_mesh_pt->finite_element_pt(e)->size();
    }
   return sum;
  }


/// Set the initial conditions: all nodes have zero velocity
void set_initial_condition() 
  {
   const unsigned n_node = Fluid_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     for(unsigned i=0;i<3;i++)
      {
       Fluid_mesh_pt->node_pt(n)->set_value(i,0.0);
      }
    }
  }



 /// Set boundary conditions on the walls
 void set_boundary_conditions(const double &time);

 /// Function that is used to run the parameter study
 void solve_system(const double &dt, const unsigned &nstep,
                   const std::string &directory);

 /// Update the problem specs before next timestep: 
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time());}

 void actions_before_adapt()
  {
   using namespace Global_Physical_Variables;
   /*char filename[100];
   sprintf(filename,"./pre_adapt_%g_%g_%g.dat", Re, Delta,Mu);
   std::ofstream out_file(filename);
   Fluid_mesh_pt->output(out_file,5);
   out_file.close();*/


   const unsigned n_boundary = this->Fluid_mesh_pt->nboundary();
   //Back up the surface meshes
   for(unsigned b=0;b<n_boundary;++b)
    {
     Backed_up_surface_mesh_pt[b] = 
      new BackupMeshForProjection<TElement<1,3> >(
       Lagrange_multiplier_mesh_pt[b],b,b);
    }
   
   // Kill the  elements and wipe surface mesh
   delete_lagrange_multiplier_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  }

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   using namespace Global_Physical_Variables;
   //Reset the lagrangian coordinates for the solid mechanics
   //an updated lagrangian approach
   //Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
   
   // Create the elements that impose the displacement constraint 
   create_lagrange_multiplier_elements();
   
   //Now to the projection
   const unsigned n_boundary = this->Fluid_mesh_pt->nboundary();
   for(unsigned b=0;b<n_boundary;++b)
    {
     Backed_up_surface_mesh_pt[b]->
      project_onto_new_mesh(this->Lagrange_multiplier_mesh_pt[b]);
    }

   this->rebuild_global_mesh();

   //Loop over all the (fluid) elements 
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to the particular element type, this is necessary because
     //the base elements don't have the member functions that we're about
     //to call!
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     
     //There is no need for ALE
     el_pt->disable_ALE();
     
     //Set the Reynolds number for each element 
     //(yes we could have different Reynolds number in each element!!)
     el_pt->re_pt() = &Global_Physical_Variables::Re;
     //Set the product of Reynolds and Strouhal numbers
     el_pt->re_st_pt() = &Global_Physical_Variables::Re;
     //Set the body force
     el_pt->axi_nst_body_force_fct_pt() = &Global_Physical_Variables::
      axial_pressure_gradient;
     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt()=
      Global_Physical_Variables::Constitutive_law_pt;
     // Set the "density" for pseudo-elastic mesh deformation
     el_pt->lambda_sq_pt()=&Global_Physical_Variables::Lambda_sq;
    }
   
   //Repin the boundary nodes
   for(unsigned b=0;b<n_boundary;b++)
    {
     unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;++n)
      {
       Node* nod_pt= Fluid_mesh_pt->boundary_node_pt(b,n);
       //Repin all the nodes
       for(unsigned i=0;i<3;i++) 
        {nod_pt->pin(i);}
      }
    }

   //Set the boundary conditions
   this->set_boundary_conditions(this->time());

   //Pin a single pressure value 
   dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(0))->fix_pressure(0,0.0);

   //Kill the backed up mesh
   for(unsigned b=0;b<n_boundary;++b)
    {
     delete Backed_up_surface_mesh_pt[b];
     Backed_up_surface_mesh_pt[b] = 0;
    }

   //Dump the output
   /*char filename[100];
   sprintf(filename,"./post_adapt_%g_%g_%g.dat", Re, Delta,Mu);
   std::ofstream out_file(filename);
   Fluid_mesh_pt->output(out_file,5);
   out_file.close();*/
  }

 
 /// Pointer to the Backedup Surface mesh
 Vector<BackupMeshForProjection<TElement<1,3> >*> Backed_up_surface_mesh_pt;

 /// Pointers to mesh of Lagrange multiplier elements
 Vector<SolidMesh*> Lagrange_multiplier_mesh_pt;

 /// Pointer to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;


//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
void create_lagrange_multiplier_elements()
{ 
 // The idea is to apply Lagrange multipliers to the boundaries in 
 // the mesh that have associated geometric objects

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
     
     // Loop over the bulk fluid elements adjacent to boundary b?
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
       Lagrange_multiplier_mesh_pt[b]->add_element_pt(el_pt);
       
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
         if (nval==8)
          {
           for (unsigned i=0;i<2;i++) 
            { 
             // Pin lagrangian multipliers
             nod_pt->pin(n_bulk_value+2+i);
            }
          }
        }
      } // end loop over the element
    } //End of  case if there is a geometric object
  } //End of loop over boundaries
}
// end of create_lagrange_multiplier_elements


//===============start_delete_lagrange_multiplier_elements==================
/// Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//==========================================================================
void delete_lagrange_multiplier_elements()
{
 unsigned n_bound = this->Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Lagrange_multiplier_mesh_pt[b]->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Lagrange_multiplier_mesh_pt[b]->element_pt(e);
    }
   
   // Wipe the mesh
   Lagrange_multiplier_mesh_pt[b]->flush_element_and_node_storage();
   
  } // end of delete_lagrange_multiplier_elements
} 

};

//============================================================================
/// Constructor: specify the maximum refinement level, the minimum and
/// maximum error targets
//============================================================================
template<class ELEMENT>
UnstructuredTorusProblem<ELEMENT>::UnstructuredTorusProblem(
 const double &min_error_target,
 const double &max_error_target)
{
 Use_continuation_timestepper=true;

 Max_residuals = 1e10;

 using namespace Global_Physical_Variables;

 //Create a timestepper
 add_time_stepper_pt(new BDF<2>);

 //Create the domain for the mesh, which consists of a circle of
 //radius Radius and centred at y=0 
 GeomObject* area_pt = new GeneralCircle(0.0,Radius);

 // No holes
 Vector<TriangleMeshClosedCurve*> Inner_hole_pt;
 
 //Now create the mesh
 double uniform_element_area = 0.2;

 // Build the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);
 
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned nsegment=8; 
 unsigned boundary_id=0; 
 curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  area_pt,zeta_start,zeta_end, 
  nsegment,boundary_id);
 
 zeta_start=-MathematicalConstants::Pi;
 zeta_end=0.0; //2.0*MathematicalConstants::Pi;
 nsegment=8; 
 boundary_id=1; 
 curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
  area_pt,zeta_start,zeta_end, 
  nsegment,boundary_id);
 
 // Combine to hole
 Vector<double> hole_coords(2);
 hole_coords[0]=0.0;
 hole_coords[1]=0.0;
 TriangleMeshClosedCurve* curvilinear_outer_boundary_pt=
  new TriangleMeshClosedCurve(curvilinear_boundary_pt);

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   curvilinear_outer_boundary_pt);

 // Take the holes into the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() =
   Inner_hole_pt;

 // Take the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 // Create the mesh
 this->Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
  triangle_mesh_parameters, this->time_stepper_pt());
 
 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Error targets for adaptive refinement
 Fluid_mesh_pt->max_permitted_error() = max_error_target; 
 Fluid_mesh_pt->min_permitted_error() = min_error_target; 
 
 //Loop over all the (fluid) elements 
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to the particular element type, this is necessary because
   //the base elements don't have the member functions that we're about
   //to call!
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

   //There is no need for ALE
   el_pt->disable_ALE();

   //Set the Reynolds number for each element 
   //(yes we could have different Reynolds number in each element!!)
   el_pt->re_pt() = &Re;
   //Set the product of Reynolds and Strouhal numbers
   el_pt->re_st_pt() = &Re;
   //Set the body force
   el_pt->axi_nst_body_force_fct_pt() = &Global_Physical_Variables::
    axial_pressure_gradient;
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt()=
    Global_Physical_Variables::Constitutive_law_pt;
   // Set the "density" for pseudo-elastic mesh deformation
   el_pt->lambda_sq_pt()=&Global_Physical_Variables::Lambda_sq;
  }

 //Let this problem be conventional form by setting gamma to zero
 ELEMENT::Gamma[0] = 0.0; //r-momentum
 ELEMENT::Gamma[1] = 0.0; //z-momentum

 
 
 //Set the boundary conditions (no slip on the torus walls)
 //Loop over the nodes on the (only) mesh boundary
 const unsigned n_boundary = Fluid_mesh_pt->nboundary();
 //Pin the boundary nodes
 for(unsigned b=0;b<n_boundary;b++)
  {
   unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
   for(unsigned n=0;n<n_boundary_node;++n)
    {
     //Repin all the nodes
     for(unsigned i=0;i<3;i++) 
      {Fluid_mesh_pt->boundary_node_pt(b,n)->pin(i);}
    }
  }
 
 //Pin a single pressure value 
 dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(0))->fix_pressure(0,0.0);

 // Create Lagrange multiplier mesh for boundary motion
 //----------------------------------------------------
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 Backed_up_surface_mesh_pt.resize(n_boundary);
 Lagrange_multiplier_mesh_pt.resize(n_boundary);
 for(unsigned b=0;b<n_boundary;++b)
  {
   Lagrange_multiplier_mesh_pt[b]=new SolidMesh;
  }
 create_lagrange_multiplier_elements();

 // Combine meshes
 //---------------
 
 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Lagrange_multiplier sub meshes
 for(unsigned b=0;b<n_boundary;++b)
  {
   this->add_sub_mesh(this->Lagrange_multiplier_mesh_pt[b]);
  }

 // Build global mesh
 this->build_global_mesh();
    
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 }

//========================================================================
/// Set the boundary conditions as a function of time we are going
/// to spin up the torus
//========================================================================
template<class ELEMENT>
void UnstructuredTorusProblem<ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //NOTE: The default value of all parameters is zero, so we need only 
 //set the values that are non-zero on the boundary, i.e. the swirl

 const unsigned n_boundary = Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary;++b)
  {
   //Loop over the nodes on the boundary
   unsigned n_boundary_node = Fluid_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     //Get the radial values
     double r = Fluid_mesh_pt->boundary_node_pt(b,n)->x(0);
     //Set the value of the u-, v- and w-velocity
     //Fast(ish) spin-up
     Fluid_mesh_pt->boundary_node_pt(b,n)->set_value(0,0.0);
     Fluid_mesh_pt->boundary_node_pt(b,n)->set_value(1,0.0);
     Fluid_mesh_pt->boundary_node_pt(b,n)->set_value(2,r*(1.0 - exp(-20.0*time)));
    }
  }
}


//==========================================================================
/// Solve the system for a number of different values of the Reynolds number
//==========================================================================
template<class ELEMENT>
void UnstructuredTorusProblem<ELEMENT>::solve_system(const double &dt, 
                                                     const unsigned &nstep,
                                                     const std::string 
                                                     &directory)
{
 using namespace Global_Physical_Variables;

 //Standard newton solve
 this->steady_newton_solve(2);
 
 char filename[100];
 std::ofstream file;
 //Output data after the first timestep
 //Create the filename, including the array index
 sprintf(filename,"%s/soln_De%g_Delta%g_Mu%g.dat",directory.c_str(),Dean,
         Delta,Mu);
 //Actually, write the data
 file.open(filename);
 Fluid_mesh_pt->output(file,5);
 file.close();

 Bifurcation_detection = true;

 //Crank up Dean number
 double ds = 50.0;
 for(unsigned i=0;i<50;i++)
  {
   //Delta += 0.1;
   ds = this->arc_length_step_solve(&Dean,ds,1);
   //this->steady_newton_solve(1);
   this->Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

   sprintf(filename,"%s/soln_De%g_Delta%g_Mu%g.dat",directory.c_str(),Dean,
           Delta,Mu);
   //Actually, write the data
   file.open(filename);
   Fluid_mesh_pt->output(file,5);
   file.close();
  }

 //Now move shapes from circle to square
 this->reset_arc_length_parameters();
 ds = +0.1;
 bool exit_flag = false;
 for(unsigned i=0;i<10;i++)
  {
   if(Mu > 1.0) 
    {
     Mu = 1.0;
     this->steady_newton_solve(1);
     exit_flag = true;
    }
   else
    {
     //Mu -= 0.1;
     //this->steady_newton_solve(1);
     ds = this->arc_length_step_solve(&Mu,ds,1);
    }

   this->Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
   
   sprintf(filename,"%s/soln_De%g_Delta%g_Mu%g.dat",directory.c_str(),Dean,
           Delta,Mu);
   //Actually, write the data
   file.open(filename);
   Fluid_mesh_pt->output(file,5);
   file.close();
   
   if(exit_flag) {break;}
  }

}

//Main driver loop
int main()
{

 double max_error = 1.0e-3;
 double min_error = 1.0e-5;

 {
  UnstructuredTorusProblem<
   ProjectableAxisymmetricTaylorHoodElement<MyTaylorHoodElement> >
    /* PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement, 
       TPVDElement<2,3> > > >   */
   problem(min_error,max_error);

  //Now timestep
  problem.solve_system(0.01,2,"RESLT_TH");
 }


 }









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
//
//Flow in a torus computed by using a cylindrical polar coordinate system
//and assuming axisymmetry.

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

//=================================================================
//A namespace for the physical parameters in the problem
//=================================================================
namespace Global_Physical_Variables
{
 /// The Reynolds number
 double Re = 0.0;

 /// The radius of the torus
 double Radius = 1.0;

 /// The curvature of the torus
 double Delta = 0.128;
};



//=start_of_FillCircle=============================================
//A geometric object that represents the geometry of the domain
//a circle of given centre and radius. One could use a non-linear
//stretch in r (xi[0]) to shift the elements towards the edge
//(boundary layer).
//=================================================================
class GeneralCircle : public GeomObject
{
public:

 /// Constructor that takes the centre position and raidus of the circle
 /// as its arguments
 GeneralCircle(const double &centre_x, const double &centre_y,
                  const double &radius) :
  GeomObject(1,2), Centre_x(centre_x), Centre_y(centre_y), Radius(radius) { }
 
/// Destructor
virtual ~GeneralCircle(){}

/// Lagrangian coordinate xi
void position (const Vector<double>& xi, Vector<double>& r) const
{
 r[0] = Centre_x + Radius*cos(xi[0]);
 r[1] = Centre_y + Radius*sin(xi[0]);
}


/// Return the position of the circle as a function of time 
/// (doesn't move as a function of time)
void position(const unsigned& t, 
              const Vector<double>& xi, Vector<double>& r) const
  {
   position(xi,r);
  }

private:

 /// Storage for the x-coordinate of the centre
 double Centre_x;
 
 /// Storage for the y-coordinate of the centre
 double Centre_y;

 /// Storage for the radius of the circle
 double Radius;

};


namespace oomph
{

//==============================================================
/// Overload Element to allow calculation of the flux
//==============================================================
template<class ELEMENT>
class MyAxisymmetricFluidElement : public ELEMENT
 {
 public:

  /// Empty constructor
  MyAxisymmetricFluidElement(){ }

  /// Get square of L2 norm of velocity components
  /// and add to the entries in the vector norm
  double square_of_l2_norm()
   {
    //Initialise the sum to zero
    double sum = 0.0;
    
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();
    
    //Find the indices at which the local velocities are stored
    unsigned u_nodal_index[3];
    for(unsigned i=0;i<3;i++) {u_nodal_index[i] = this->u_index_axi_nst(i);}
    
    //Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node,2);
    
    //Number of integration points
    const unsigned n_intpt = this->integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(2);
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<2;i++) s[i] = this->integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = this->integral_pt()->weight(ipt);
      
      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt,psif,dpsidx);


      // Calculate position
      Vector<double> interpolated_x(2,0.0);
      //Calculate velocities 
      Vector<double> interpolated_u(3,0.0);      
      
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
       {
        const double psif_ = psif[l];

        //Loop over physical coordinates
        for(unsigned i=0;i<2;i++)
         {
          interpolated_x[i] += this->raw_nodal_position(l,i)*psif_;
         }

        //Loop over velocity components
        for(unsigned i=0;i<3;i++)
         {
          interpolated_u[i] += this->raw_nodal_value(l,u_nodal_index[i])*psif_;
         }
       }

      
      //Premultiply the weights and the Jacobian
      double W = interpolated_x[0]*w*J;

      //Assemble square of L2 norm
      for(unsigned i=0;i<3;i++)
       {
        sum +=interpolated_u[i]*interpolated_u[i]*W;
       }           
     }

    return sum;
   }

 };


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<class ELEMENT>
 class FaceGeometry<MyAxisymmetricFluidElement<ELEMENT> >
  : public virtual FaceGeometry<ELEMENT> 
 {
 public:
  FaceGeometry() : FaceGeometry<ELEMENT>() {}
 };


} //End of namespace extension



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
   const unsigned n_element = this->mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     sum += dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(e))
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
   const unsigned n_element = this->mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     sum += this->mesh_pt()->finite_element_pt(e)->size();
    }
   return sum;
  }


/// Set the initial conditions: all nodes have zero velocity
void set_initial_condition() 
  {
   const unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     for(unsigned i=0;i<3;i++)
      {
       mesh_pt()->node_pt(n)->set_value(i,0.0);
      }
    }
  }

 /// Set boundary conditions on the walls
 void set_boundary_conditions(const double &time);

 /// Function that is used to run the parameter study
 void solve_system(const double &dt, const unsigned &nstep,
                   const std::string &directory);
 
 /// Return a pointer to the specific mesh used
 RefineableTriangleMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());}
 /// Update the problem specs before next timestep: 
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time());}

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   //Loop over all the (fluid) elements 
   unsigned n_element = mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to the particular element type, this is necessary because
     //the base elements don't have the member functions that we're about
     //to call!
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
     
     //There is no need for ALE
     el_pt->disable_ALE();
     
     //Set the Reynolds number for each element 
     //(yes we could have different Reynolds number in each element!!)
     el_pt->re_pt() = &Global_Physical_Variables::Re;
     //Set the product of Reynolds and Strouhal numbers
     el_pt->re_st_pt() = &Global_Physical_Variables::Re;
    }
   
   const unsigned n_boundary = mesh_pt()->nboundary();
   //Repin the boundary nodes
   for(unsigned b=0;b<n_boundary;b++)
    {
     unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
     Vector<double> new_x(2);
     Vector<double> b_coord(1);
     bool geom_object=false;
     if(this->mesh_pt()->boundary_geom_object_pt(b)!=0) 
      {geom_object=true;}
     for(unsigned n=0;n<n_boundary_node;++n)
      {
       //Now move the boundary nodes exactly onto the geometric object
       Node* const nod_pt = this->mesh_pt()->boundary_node_pt(b,n);
       if(geom_object)
        {
         nod_pt->get_coordinates_on_boundary(b,b_coord);
         //Let's find boundary coordinates of the new node
         this->mesh_pt()->boundary_geom_object_pt(b)->position(b_coord,new_x);
         //Snap to the boundary
         for(unsigned i=0;i<2;i++)
          {
           nod_pt->x(i) = new_x[i];
          }
        }
       //Repin all the nodes
       for(unsigned i=0;i<3;i++) 
        {nod_pt->pin(i);}
      }
    }

   //Set the boundary conditions
   this->set_boundary_conditions(this->time());

   //Pin a single pressure value 
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
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
 
 using namespace Global_Physical_Variables;

 //Create a timestepper
 add_time_stepper_pt(new BDF<2>);

 //Create the domain for the mesh, which consists of a circle of
 //radius Radius and centred at (1/Delta, 0) 
 GeomObject* area_pt = new GeneralCircle(1.0/Delta,0.0,Radius);

 //Specify the number of boundary points
 unsigned n_points = 11;
 
 //Set the range of the coordinates
 double xi_frac = 4.0*atan(1.0)/((double)(n_points-1));
 
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshPolyLine*> boundary_segment_pt(2);
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(n_points);
 for(unsigned i=0;i<n_points;i++)
  {bound_seg[i].resize(2);}
 
 //Position vector required to build the hole
 Vector<double> zeta(1);
 // Initialize the vector of coordinates
 Vector<double> coord(2); 
 
//Create the boundary points
 for(unsigned p=0;p<n_points;p++)
  {
   double zeta_spec = 0.0 + p*xi_frac;
   
   zeta[0] = zeta_spec;
   area_pt->position(zeta,coord);
   // First boundary segment
   bound_seg[p][0]=coord[0];
   bound_seg[p][1]=coord[1];
  }
 
 // Specify 1st boundary id
 unsigned bound_id = 0;
 
 // Build the 1st boundary segment
 boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 

//Create the next set of boundary points
 for(unsigned p=0;p<n_points;p++)
  {
   double zeta_spec = 4.0*atan(1.0) + p*xi_frac;
   
   zeta[0] = zeta_spec;
   area_pt->position(zeta,coord);
   // First boundary segment
   bound_seg[p][0]=coord[0];
   bound_seg[p][1]=coord[1];
  }
 
 // Specify 2nd boundary id
 bound_id = 1;
 
 // Build the 2nd boundary segment
 boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 
//  // Create the triangle mesh polygon for outer boundary using boundary segment
//  TriangleMeshPolygon* 
//   Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_segment_pt);

 // No holes
 Vector<TriangleMeshClosedCurve*> Inner_hole_pt;
 
 //Now create the mesh
 double uniform_element_area = 0.2;
 //Problem::mesh_pt() = new RefineableTriangleMesh<ELEMENT>(
 // Outer_boundary_polyline_pt, 
 // Inner_hole_pt,
 // uniform_element_area,
 // this->time_stepper_pt());
 




 // Build the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);
 
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned nsegment=8; 
 unsigned boundary_id=0; 
 curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  area_pt,zeta_start,zeta_end, 
  nsegment,boundary_id);
 
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
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
 

//  //Set the boundaries
//  Vector<Vector<double> > outer_split_coord(2);
//  outer_split_coord[0].resize(2);
//  outer_split_coord[0][0] = 0.0; outer_split_coord[0][1] = 4.0*atan(1.0);
//  outer_split_coord[1].resize(2);
//  outer_split_coord[1][0] = outer_split_coord[0][1]; 
//  outer_split_coord[1][1] = 8.0*atan(1.0);

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
 Problem::mesh_pt() = new RefineableTriangleMesh<ELEMENT>(
   triangle_mesh_parameters, this->time_stepper_pt());

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error() = max_error_target; 
 mesh_pt()->min_permitted_error() = min_error_target; 
 
 //Loop over all the (fluid) elements 
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to the particular element type, this is necessary because
   //the base elements don't have the member functions that we're about
   //to call!
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //There is no need for ALE
   el_pt->disable_ALE();

   //Set the Reynolds number for each element 
   //(yes we could have different Reynolds number in each element!!)
   el_pt->re_pt() = &Re;
   //Set the product of Reynolds and Strouhal numbers
   el_pt->re_st_pt() = &Re;
  }

 //Let this problem be conventional form by setting gamma to zero
 ELEMENT::Gamma[0] = 0.0; //r-momentum
 ELEMENT::Gamma[1] = 0.0; //z-momentum
 
 //Set the boundary conditions (no slip on the torus walls)
 //Loop over the nodes on the (only) mesh boundary
 const unsigned n_boundary = mesh_pt()->nboundary();
 //Pin the boundary nodes
 for(unsigned b=0;b<n_boundary;b++)
  {
   unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   Vector<double> new_x(2);
   Vector<double> b_coord(1);
   for(unsigned n=0;n<n_boundary_node;++n)
    {
     //Now move the boundary nodes exactly onto the geometric object
     Node* const nod_pt = this->mesh_pt()->boundary_node_pt(b,n);
     nod_pt->get_coordinates_on_boundary(b,b_coord);
     //Let's find boundary coordinates of the new node
     this->mesh_pt()->boundary_geom_object_pt(b)->position(b_coord,new_x);
     //Snap to the boundary
     for(unsigned i=0;i<2;i++)
      {
       nod_pt->x(i) = new_x[i];
      }
     
     //Repin all the nodes
     for(unsigned i=0;i<3;i++) 
      {mesh_pt()->boundary_node_pt(b,n)->pin(i);}
    }
  }
 
 //Pin a single pressure value 
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
 
 //Setup all the equation numbering and look-up schemes 
 std::cout << assign_eqn_numbers() << std::endl; 
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

 const unsigned n_boundary = mesh_pt()->nboundary();
 for(unsigned b=0;b<n_boundary;++b)
  {
   //Loop over the nodes on the boundary
   unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     //Get the radial values
     double r = mesh_pt()->boundary_node_pt(b,n)->x(0);
     //Set the value of the u-, v- and w-velocity
     //Fast(ish) spin-up
     mesh_pt()->boundary_node_pt(b,n)->set_value(0,0.0);
     mesh_pt()->boundary_node_pt(b,n)->set_value(1,0.0);
     mesh_pt()->boundary_node_pt(b,n)->set_value(2,r*(1.0 - exp(-20.0*time)));
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

 //Open a trace file
 std::stringstream trace_name;
 trace_name << directory << "/time_trace.dat";
 ofstream trace(trace_name.str().c_str());

 //Define a string that we can set to be the name of the output file
 char filename[100];
 //Define an output filestream
 ofstream file;

 //Set the Reynolds number
 Re = 10.0;//1000.0;

 //Set an impulsive start from rest
 assign_initial_values_impulsive(dt);

 //Calculate the l2 norm
 double l2_norm = calculate_square_of_l2_norm();
 double area = calculate_area();

 //Output intital data
 trace << time() << " " << area << " " << l2_norm << std::endl;

 //Increase the maximum value of the residuals to get
 //past the first few steps
 Max_residuals = 50.0;

 //Now perform the first timestep with 2 steps of refinement
 unsteady_newton_solve(dt,2,true);

 //Calculate the l2 norm and area
 l2_norm = calculate_square_of_l2_norm();
 area = calculate_area();

 //Output intital data
 trace << time() << " " << area << " " << l2_norm << std::endl;
 
 //Output data after the first timestep
 //Create the filename, including the array index
 sprintf(filename,"%s/soln_Re%g_t%g.dat",directory.c_str(),Re,time());
 //Actually, write the data
 file.open(filename);
 mesh_pt()->output(file,5);
 file.close();
 
 //Now do the other steps with only one adaptation per step
 for(unsigned n=1;n<nstep;n++)
  {
   //Solve the problem
   unsteady_newton_solve(dt,1,false);
   
   //Calculate the l2 norm and area
   l2_norm = calculate_square_of_l2_norm();
   area = calculate_area();
   
   //Output intital data
   trace << time() << " " << area << " " << l2_norm << std::endl;
   
   //Output data at each step
   //Create the filename, including the array index
   sprintf(filename,"%s/soln_Re%g_t%g.dat",directory.c_str(),Re,time());
   //Actually, write the data
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();
  }

 //Close the trace file
 trace.close();
}

//Main driver loop
int main()
{

 //Construct and solve the problem
 //This maximum refinement level means that we fit (easily) into a 
 //1G machine, but it could probably go higher if you start to 
 //see refinement being overruled.
 //If you have this too high you get ridiculous refinement at early 
 //times that isn't really necessary.
 double max_error = 1.0e-3;
 double min_error = 1.0e-5;


{
  UnstructuredTorusProblem<MyAxisymmetricFluidElement<
   ProjectableAxisymmetricCrouzeixRaviartElement<
   AxisymmetricTCrouzeixRaviartElement> > >  problem(
    min_error,max_error);
  
  //Now timestep
  problem.solve_system(0.01,2,"RESLT_CR");
 }

 {
  UnstructuredTorusProblem<MyAxisymmetricFluidElement<
   ProjectableAxisymmetricTaylorHoodElement<
   AxisymmetricTTaylorHoodElement> > >  problem(
    min_error,max_error);
  
  //Now timestep
  problem.solve_system(0.01,2,"RESLT_TH");
 }


 }









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

//Convecting-heated sphere

// Generic oomph-lib stuff
#include "generic.h"

//Specific mesh
#include "half_rectangle_with_hole_mesh.h"

// Navier Stokes
#include "axisym_buoyant_navier_stokes.h"

using namespace std;

using namespace oomph;



//===============================================
/// Global parameters
//===============================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re = 1.0;

 /// Prandtl number
 double Pr = 0.73;

 /// Peclet number
 double Pe = 1.0;

 /// Rayleigh number
 double Ra = 10.0;
 
 /// Gravity
 Vector<double> G;

 /// Location of the centre of the sphere on the axis
 double Sphere_centre_z= 0.0;

}


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

namespace StokesFlowExactWind
{

 /// Wind that represents a constantly translating sphere
 void get_wind(const Vector<double>& x, Vector<double>& wind)
 {
  double x_sc = x[0] - 0.0;
  double z_sc = x[1] - Global_Parameters::Sphere_centre_z;
  double r = sqrt(x_sc*x_sc + z_sc*z_sc);
  double theta = atan2(x_sc,z_sc);
  double sp_r=0.5;
  double r_sc = r/sp_r;

  double u_r_sph = -(-1.0  + 1.5/(r_sc) - 0.5/(r_sc*r_sc*r_sc))*cos(theta);
  double u_theta_sph = -(1.0 - 0.75/(r_sc) - 0.25/(r_sc*r_sc*r_sc))*sin(theta);

  wind[0]= u_r_sph*sin(theta) + u_theta_sph*cos(theta);
  wind[1]= u_r_sph*cos(theta) - u_theta_sph*sin(theta);
 }

}



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//=====================================================================
/// Class of face elements whose sole raison d'etre is to calculate the
/// drag and mass transfter on the boundary of the selected sphere
//=====================================================================
template <class ELEMENT>
class DragNusseltCalculationElement : 
 public virtual FaceGeometry<ELEMENT>, 
 public virtual FaceElement 
{
 
public:

 /// Constructor, which takes a "bulk" element and the value of the index
 /// and its limit
 DragNusseltCalculationElement(FiniteElement* const &element_pt, 
                                      const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
  }

 /// Destructor should not delete anything
 ~DragNusseltCalculationElement() { }

 /// Return the torque calculation
 void calculate_drag_nusselt(Vector<double> &drag, 
                             double &nusselt, double &area)
  {
   //Set the value of n_intpt
   const unsigned n_intpt = integral_pt()->nweight();

   //Get the dimension of the element
   const unsigned dim = this->dim();

   //Get the storage for the local coordinate
   Vector<double> s(dim);

   //Get the position in the bulk
   Vector<double> s_bulk(dim+1);

   //Get the outer unit normal
   Vector<double> N(2,0.0);

   //Need to find position 
   Vector<double> x(dim+1);

   //Need the traction as well
   Vector<double> traction(3,0.0);

   //Need the flux
   Vector<double> flux(2,0.0);

   double sum[4] = {0.0,0.0,0.0,0.0};

   ELEMENT* const bulk_elem_pt = 
    dynamic_cast<ELEMENT*>(this->bulk_element_pt());

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the local coordinates
     for(unsigned i=0;i<dim;i++) {s[i] = integral_pt()->knot(ipt,i);}

     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Find the shape and test functions and return the Jacobian
     //of the mapping
     double J = J_eulerian(s);

     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Now get the position in the bulk
     this->get_local_coordinate_in_bulk(s,s_bulk);

     //Get the position
     this->interpolated_x(s,x);

     //Get the normal
     //Resize to two to stop PARANOID compaints
     N.resize(2);
     this->outer_unit_normal(s,N);

     //Now get the traction from the bulk element
     //Resize normal to three to stop PARANOID complaints
     N.resize(3,0.0);
     bulk_elem_pt->get_traction(s_bulk,N,traction);

     //Get the concentration flux from the bulk element
     bulk_elem_pt->get_flux(s_bulk,flux);

     //Now work out the radial drag, the z-component of the traction 
     sum[0] += traction[0]*x[0]*W;
     //Axial drag
     sum[1] += traction[1]*x[0]*W;

     //Also sort out the normal mass transfer which is 
     sum[2] += (N[0]*flux[0] + N[1]*flux[1])*x[0]*W;

     //Area
     sum[3] += x[0]*W;
    }

   drag[0] = sum[0]; drag[1] = sum[1]; nusselt = sum[2]; area = sum[3];
  }

}; 






//===================================================================
/// Flow around a cylinder in rectangular domain
//===================================================================
template<class ELEMENT>
class FlowAroundHalfCylinderProblem : public Problem
{

public:

 /// Constructor: Pass geometric object that represents
 /// central cylinder, and length and height of domain.
 FlowAroundHalfCylinderProblem(GeomObject* cylinder_pt, 
                               const double &radius, 
                               const double &length);
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty; all prescribed
 /// velocities are constant along their respective boundares, therefore
 /// their FE interpolation onto the newly created nodes is sufficiently
 /// accurate)
 void actions_before_newton_solve() {}

 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
   // Pin redundant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());

   //Pin all swirl velocities to zero
   unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     mesh_pt()->node_pt(n)->pin(2);
     mesh_pt()->node_pt(n)->set_value(2,0.0);
    }


   //Plug flow everywhere
   for(unsigned ibound=0;ibound<5;++ibound)
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //No radial flow on symmetry boundary
       if(ibound==3) 
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }

       //No radial flow or concentration on outlet
       if(ibound==2) 
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }

       //No flow and concentration on side boundary
       if(ibound==1)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,0.0);
        }

       //Plug flow and no concentration on inlet only
       if(ibound==0) 
        {
         //Specify the exact flow at inlet
         /*Vector<double> x(2);
         x[0] = mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
         x[1] = mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
         Vector<double> wind(2);
         StokesFlowExactWind::get_wind(x,wind);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,wind[0]);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,wind[1]);*/
        

         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,0.0);
        }

       //No slip on sphere and fixed concentration on sphere
       if(ibound==4)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,1.0);
       }
      }
    }

  } // end_of_actions_after_adapt


 /// Access function for the specific mesh
 RefineableHalfRectangleWithHoleMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableHalfRectangleWithHoleMesh<ELEMENT>*>
    (Problem::mesh_pt());
  }


 void compute_drag_nusselt(Vector<double> &drag, double &nusselt, double &area)
  {
   unsigned bound = 4;
   //Loop over the elements adjacent to the boundary 4 and make face elements
   unsigned n_bound_element = this->mesh_pt()->nboundary_element(bound);

   drag[0] = 0.0; drag[1] = 0.0; nusselt = 0.0; area = 0.0;
   
   for(unsigned e=0;e<n_bound_element;e++)
    {
     //Get pointer to the bulk element
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      this->mesh_pt()->boundary_element_pt(bound,e));
     
     //FInd the face index
     int face_index = this->mesh_pt()->face_index_at_boundary(bound,e);
     
     //Build the flux element
     DragNusseltCalculationElement<ELEMENT>* drag_element_pt = new
      DragNusseltCalculationElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Now calculate the torque
     Vector<double> el_drag(2);
     double  el_nusselt, el_area;
     drag_element_pt->calculate_drag_nusselt(el_drag,el_nusselt,el_area);
     //Delete our element (it's work is done)
     delete drag_element_pt;
     
     //Add elemental contribution to total
     drag[0] += el_drag[0]; drag[1] += el_drag[1]; 
     nusselt += el_nusselt; area += el_area;
    }
   
   //Need to multiply the drag and area by the aziumthal component
   double two_pi = 2.0*MathematicalConstants::Pi;
   
   //Multiply the drag by the surface area of the sphere
   drag[0] *= two_pi; drag[1] *= two_pi; area *= two_pi;
  }


 private:
 
  /// Height of the domain
  double Domain_radius;

 /// Length of the domain
 double Domain_length;
 
};




//========================================================================
/// Constructor 
//========================================================================
template<class ELEMENT>
FlowAroundHalfCylinderProblem<ELEMENT>::FlowAroundHalfCylinderProblem(
 GeomObject* cylinder_pt, const double &radius, const double &length) 
{ 
 Domain_radius = radius;
 Domain_length = length;
  
 // Build mesh
 Problem::mesh_pt()=
  new RefineableHalfRectangleWithHoleMesh<ELEMENT>(cylinder_pt,radius,length,
                                                 4.0,2,4.0,2,1.0,2);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 //Allow a fair amount of refinement
 //mesh_pt()->max_refinement_level() = 10;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 

 //Pin all swirl velocities
 unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   mesh_pt()->node_pt(n)->pin(2);
  }
 
 //Set boundary conditions
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Pin boundary conditions on cylinder, inlet and side
     if((ibound==0) || (ibound==1) || (ibound==4))
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(3);
      }
     //Otherwise pin radial flow
     else 
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }

     //Pin outlet concentration
     /*if(ibound==2)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(3);
       }*/
    }
  }
    
  // Pin redudant pressure dofs
  RefineableAxisymmetricNavierStokesEquations::
   pin_redundant_nodal_pressures(mesh_pt()->element_pt());
  
  // Pass pointer to Reynolds number to elements
  unsigned nelem=mesh_pt()->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

    // Set the Peclet number
    el_pt->pe_pt() = &Global_Parameters::Pe;
    
    // Set the Peclet number multiplied by the Strouhal number
    el_pt->pe_st_pt() =&Global_Parameters::Pe;
    
    // Set the Reynolds number (1/Pr in our non-dimensionalisation)
    el_pt->re_pt() = &Global_Parameters::Re;
    
   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Parameters::Re;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Parameters::Ra;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Parameters::G;
    
   }
  
  //Attach the boundary conditions to the mesh
  cout <<"Number of equations: " << assign_eqn_numbers() << endl; 


  //Set swirl velocity to zero
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->node_pt(n)->set_value(2,0.0);
   }

  /*ofstream boundary[5];
  boundary[0].open("bound0.dat");
  boundary[1].open("bound1.dat");
  boundary[2].open("bound2.dat");
  boundary[3].open("bound3.dat");
  boundary[4].open("bound4.dat");

  ofstream output("initial.dat");
  mesh_pt()->output(output,3);
  output.close();*/
 
  //Plug flow everywhere
  for(unsigned ibound=0;ibound<5;++ibound)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      /*boundary[ibound] << mesh_pt()->boundary_node_pt(ibound,inod)->x(0) 
                       << " " 
                       << mesh_pt()->boundary_node_pt(ibound,inod)->x(1) 
                       << "\n";*/

      //No radial flow on symmetry boundary
      if(ibound==3) 
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       }
      
      //No radial flow or concentration on outlet
      if(ibound==2) 
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        // mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,0.0);
       }
      
      //No radial flow on side boundary
      if(ibound==1)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,0.0);
       }
      

      //Zero flow and zero concentration on inlet 
      if(ibound==0)
       {
        //Specify the exact flow at inlet
        /*Vector<double> x(2);
        x[0] = mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
        x[1] = mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
        Vector<double> wind(2);
        StokesFlowExactWind::get_wind(x,wind);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,wind[0]);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,wind[1]);*/

        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,0.0);
       }

      //No slip on sphere and fixed concentration
      if(ibound==4)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(3,1.0);
       }
     }
   }
}

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=====================================================================
/// Driver
//=====================================================================
int main()
{
 Global_Parameters::G.resize(3);
 Global_Parameters::G[0] = 0.0;
 Global_Parameters::G[1] = -1.0;
 Global_Parameters::G[2] = 0.0;
 
 // radius and Length of domain
 double radius=5.0;
 double length=10.0;
 Global_Parameters::Sphere_centre_z = 5.0;

 //Create a new ellipse object as the central cylinder
 HalfEllipse* cylinder_pt = 
  new HalfEllipse(Global_Parameters::Sphere_centre_z,0.5,0.5);

 // Create Problem
 FlowAroundHalfCylinderProblem 
 <RefineableBuoyantQAxisymCrouzeixRaviartElement> 
  problem(cylinder_pt,radius,length);

 // Solve adaptively with up to max_adapt rounds of refinement
 unsigned max_adapt=1;

 // Output filename
 char filename[100];
 // Trace file
 std::ofstream trace("trace.dat");

 //Step up in the Rayleigh number
 for(unsigned i=0;i<2;i++)
  {
   problem.newton_solve(max_adapt);

   //Open an output file
   sprintf(filename,"soln_Re%g_Ra%g.dat", Global_Parameters::Re,
           Global_Parameters::Ra);
   //Doc result
   ofstream outfile(filename);
   problem.mesh_pt()->output(outfile,5);
   outfile.close();

   //Compute drag and mass transfer
   Vector<double> drag(2); double nusselt = 0.0; double area = 0.0;
   problem.compute_drag_nusselt(drag,nusselt, area);
   
   trace << Global_Parameters::Re << " " << Global_Parameters::Pr
         << " " << Global_Parameters::Pe << " " 
         << drag[0] << " " << drag[1] << " " << 2.0*nusselt << std::endl;

   //Global_Parameters::Re += 10.0;
   //Global_Parameters::Pe = Global_Parameters::Pr*Global_Parameters::Re;
   Global_Parameters::Ra += 40.0;
  }

 trace.close();
}

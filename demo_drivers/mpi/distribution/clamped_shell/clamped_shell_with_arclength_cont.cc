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
//Driver function for a simple test shell problem:
//Calculate the deformation of an elastic tube approximated
//using Kirchoff--Love shell theory

//Standard system includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio>

//Include files from the finite-element library
#include "generic.h"
#include "shell.h"
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//======================================================================== 
/// Global variables that represent physical properties
//======================================================================== 
namespace Global_Physical_Variables
{

 /// Prescribed position of control point
 double Prescribed_y = 1.0;

 ///  Pointer to pressure load (stored in Data so it can 
 /// become an unknown in the problem when displacement control is used
 Data* Pext_data_pt;

 ///  Return a reference to the external pressure 
 /// load on the elastic tube.
 double external_pressure() 
  {return (*Pext_data_pt->value_pt(0))*pow(0.05,3)/12.0;}


 /// Load function, normal pressure loading
 void press_load(const Vector<double> &xi,
                 const Vector<double> &x,
                 const Vector<double> &N,
                 Vector<double>& load)
 {
  for(unsigned i=0;i<3;i++) 
   {
    load[i] = external_pressure()*N[i];
   }
 }

}

//======================================================================== 
/// A 2D Mesh class. The tube wall is represented by two Lagrangian 
/// coordinates that correspond to z and theta in cylindrical polars. 
/// The required mesh is therefore a  2D mesh and is therefore inherited 
/// from the generic RectangularQuadMesh 
//=======================================================================
template <class ELEMENT>
class ShellMesh : public virtual RectangularQuadMesh<ELEMENT>,
                  public virtual SolidMesh
{
public:
 
 ///Constructor for the mesh
 ShellMesh(const unsigned &nx, const unsigned &ny, 
           const double &lx, const double &ly);
 
 ///  In all elastic problems, the nodes must be assigned an undeformed,
 /// or reference, position, corresponding to the stress-free state
 /// of the elastic body. This function assigns the undeformed position 
 /// for the nodes on the elastic tube
 void assign_undeformed_positions(GeomObject* const &undeformed_midplane_pt);

};





//=======================================================================
/// Mesh constructor
/// Argument list:
/// nx  : number of elements in the axial direction
/// ny : number of elements in the azimuthal direction
/// lx  : length in the axial direction
/// ly  : length in theta direction
//=======================================================================
template <class ELEMENT>
ShellMesh<ELEMENT>::ShellMesh(const unsigned &nx, 
                              const unsigned &ny,
                              const double &lx, 
                              const double &ly) :
 RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();
 
 //Now in this case it is the Lagrangian coordinates that we want to set,
 //so we have to loop over all nodes and set them to the Eulerian
 //coordinates that are set by the generic mesh generator
 for(unsigned i=0;i<n_node;i++)
  {
   node_pt(i)->xi(0) = node_pt(i)->x(0);
   node_pt(i)->xi(1) = node_pt(i)->x(1);
  }
 

 //Assign gradients, etc for the Lagrangian coordinates of 
 //hermite-type elements
 
 //Read out number of position dofs
 unsigned n_position_type = finite_element_pt(0)->nnodal_position_type();

 //If this is greater than 1 set the slopes, which are the distances between 
 //nodes. If the spacing were non-uniform, this part would be more difficult
 if(n_position_type > 1)
  {
   double xstep = (this->Xmax - this->Xmin)/((this->Np-1)*this->Nx);
   double ystep = (this->Ymax - this->Ymin)/((this->Np-1)*this->Ny);
   for(unsigned n=0;n<n_node;n++)
    {
     //The factor 0.5 is because our reference element has length 2.0
     node_pt(n)->xi_gen(1,0) = 0.5*xstep;
     node_pt(n)->xi_gen(2,1) = 0.5*ystep;
    }
  }
}


//=======================================================================
/// Set the undeformed coordinates of the nodes
//=======================================================================
template <class ELEMENT>
void ShellMesh<ELEMENT>::assign_undeformed_positions(
 GeomObject* const &undeformed_midplane_pt)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Loop over all the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Get the Lagrangian coordinates
   Vector<double> xi(2);
   xi[0] = node_pt(n)->xi(0); 
   xi[1] = node_pt(n)->xi(1);
   
   //Assign memory for values of derivatives, etc
   Vector<double> R(3);
   DenseMatrix<double> a(2,3); 
   RankThreeTensor<double>  dadxi(2,2,3);
   
   //Get the geometrical information from the geometric object
   undeformed_midplane_pt->d2position(xi,R,a,dadxi);
   
   //Loop over coordinate directions
   for(unsigned i=0;i<3;i++)
    {
     //Set the position
     node_pt(n)->x_gen(0,i) = R[i];

     //Set the derivative wrt Lagrangian coordinates
     //Note that we need to scale by the length of each element here!!
     node_pt(n)->x_gen(1,i) = 0.5*a(0,i)*((this->Xmax - this->Xmin)/this->Nx);
     node_pt(n)->x_gen(2,i) = 0.5*a(1,i)*((this->Ymax - this->Ymin)/this->Ny);

     //Set the mixed derivative 
     //(symmetric so doesn't matter which one we use)
     node_pt(n)->x_gen(3,i) = 0.25*dadxi(0,1,i);
    }
  }
}



// //======================================================================== 
// /// A 2D mesh class for the problem.
// /// The tube is represented by two Lagrangian coordinates that correspond
// /// to z and theta in cylindrical polars. The required mesh is therefore a 
// /// 2D mesh and is therefore inherited from the generic RectangularQuadMesh 
// //=======================================================================
// template <class ELEMENT>
// class ShellMesh : public virtual RectangularQuadMesh<ELEMENT>,
//                   public virtual SolidMesh
// {
// public:
// //Constructor for the mesh
//  ShellMesh(const unsigned &nx, const unsigned &ny, 
//            const double &lx, const double &ly);
 
//  //In all elastic problems, the nodes must be assigned an undeformed,
//  //or reference, position, corresponding to the stress-free state
//  //of the elastic body. This function assigns the undeformed position 
//  //for the nodes on the elastic tube
//  void assign_undeformed_positions(GeomObject* const &undeformed_midplane_pt);
// };

// //Define the mesh constructor
// //Argument list:
// // nx  : number of elements in the axial direction
// // ny : number of elements in the azimuthal direction
// // lx  : length in the axial direction
// // ly  : length in theta direction
// template <class ELEMENT>
// ShellMesh<ELEMENT>::ShellMesh(const unsigned &nx, 
//                               const unsigned &ny,
//                               const double &lx, const double &ly) :
// RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
// {
//  //Find out how many nodes there are
//  unsigned Nnode = nnode();

//  //Now in this case it is the Lagrangian coordinates that we want to set,
//  //so we have to loop over all nodes and set them to the Eulerian
//  //coordinates that are set by the generic mesh generator
//  for(unsigned i=0;i<Nnode;i++)
//   {
//    node_pt(i)->xi(0) = node_pt(i)->x(0);
//    node_pt(i)->xi(1) = node_pt(i)->x(1);
//   }

//  //Assign gradients, etc for the Lagrangian coordinates of 
//  //hermite-type elements
 
//  //Read out number of position dofs
//  unsigned Nposition_type = finite_element_pt(0)->nnodal_position_type();
//  std::cout << Nposition_type << std::endl;
//  //If this is greater than 1 set the slopes, which are the distances between 
//  //nodes. If the spacing were non-uniform, this part would be more difficult
//  if(Nposition_type > 1)
//   {
//    double xstep = (this->Xmax - this->Xmin)/((this->Np-1)*this->Nx);
//    double ystep = (this->Ymax - this->Ymin)/((this->Np-1)*this->Ny);
//    for(unsigned n=0;n<Nnode;n++)
//     {
//      //The factor 0.5 is because our reference element has length 2.0
//      node_pt(n)->xi_gen(1,0) = 0.5*xstep;
//      node_pt(n)->xi_gen(2,1) = 0.5*ystep;
//     }
//   }
// }

// //Set the undeformed values of the nodes
// template <class ELEMENT>
// void ShellMesh<ELEMENT>::
// assign_undeformed_positions(GeomObject* const &undeformed_midplane_pt)
// {
//  //Find out how many nodes there are
//  unsigned n_node = nnode();
//  //Loop over all the nodes
//  for(unsigned n=0;n<n_node;n++)
//   {
//    //Get the lagrangian coordinates
//    Vector<double> xi(2);
//    xi[0] = node_pt(n)->xi(0); xi[1] = node_pt(n)->xi(1);

//    //Assign memory for values of derivatives, etc
//    Vector<double> R(3);
//    DenseMatrix<double> a(2,3); 
//    RankThreeTensor<double>  dadxi(2,2,3);

//    //Get the geometrical information from the geometric object
//    undeformed_midplane_pt->d2position(xi,R,a,dadxi);

//    //Loop over coordinate directions
//    for(unsigned i=0;i<3;i++)
//     {
//      //Set the position
//      node_pt(n)->x_gen(0,i) = R[i];
//      //Set the derivative wrt Lagrangian coordinates
//      //Note that we need to scale by the length of each element here!!
//      node_pt(n)->x_gen(1,i) = 0.5*a(0,i)*((this->Xmax - this->Xmin)/this->Nx);
//      node_pt(n)->x_gen(2,i) = 0.5*a(1,i)*((this->Ymax - this->Ymin)/this->Ny);
//      //Set the mixed derivative 
//      //(symmetric so doesn't matter which one we use)
//      node_pt(n)->x_gen(3,i) = 0.25*dadxi(0,1,i);
//     }
//   }
// }



//======================================================================
//Problem class to solve the deformation of an elastic tube
//=====================================================================
template<class ELEMENT>
class ShellProblem : public Problem
{

public:

 /// Constructor
 ShellProblem(const unsigned &nx, const unsigned &ny, 
              const double &lx, const double &ly);

 /// Destructor: delete mesh, geometric object
 ~ShellProblem()
  {
   delete Displacement_control_mesh_pt;
   delete Solid_mesh_pt;
   delete Undeformed_midplane_pt;
  }


 /// Overload Access function for the mesh
 ShellMesh<ELEMENT>* solid_mesh_pt() 
  {return dynamic_cast<ShellMesh<ELEMENT>*>(this->Solid_mesh_pt);}

 /// Actions after solve empty
 void actions_after_newton_solve() {}

 /// Actions before solve empty
 void actions_before_newton_solve() {}
 
 ///  Actions after problem distribution.
 /// Need to reset the pointer to stored shape functions on all processors
 void actions_after_distribute()
  {
   //If there are any solid elements calculate and store the 
   //shape functions in the first element
   unsigned n_solid_element = this->solid_mesh_pt()->nelement();

   //Explicit pointer to first element in the mesh if there is one
   ELEMENT* first_el_pt =0;
   if(n_solid_element > 0) 
    {
     first_el_pt = 
      dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(0));
    }
   
   //Loop over all elements
   for(unsigned e=0;e<n_solid_element;e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(e));
     
     //Pre-compute the second derivatives wrt Lagrangian coordinates
     //for the first element only
     if(e==0)
      {
       el_pt->pre_compute_d2shape_lagrangian_at_knots();
      }
     //Otherwise set the values to be the same as those in the first element
     //this is OK because the Lagrangian mesh is uniform.
       //
       else
        {
         el_pt->set_dshape_lagrangian_stored_from_element(first_el_pt);
        }
    }
  }

 //A self_test function
 void solve();

private:

 /// Pointer to GeomObject that specifies the undeformed midplane
 GeomObject* Undeformed_midplane_pt;

 /// First trace node
 Node* Trace_node_pt;

 /// Second trace node
 Node* Trace_node2_pt;

 /// Pointer to the solid mesh
 Mesh* Solid_mesh_pt;

 /// Pointer to the mesh that contains the displacement control element
 Mesh* Displacement_control_mesh_pt;

};



//======================================================================
/// Constructor
//======================================================================
template<class ELEMENT>
ShellProblem<ELEMENT>::ShellProblem(const unsigned &nx, const unsigned &ny, 
                                    const double &lx, const double &ly)
{
 Mesh::Suppress_warning_about_empty_mesh_level_time_stepper_function=true;

 Use_continuation_timestepper=true;
 //Create the undeformed midplane object
 Undeformed_midplane_pt = new EllipticalTube(1.0,1.0);

 //Now create the mesh
 Solid_mesh_pt = new ShellMesh<ELEMENT>(nx,ny,lx,ly); 

 //Set the undeformed positions in the mesh
 solid_mesh_pt()->assign_undeformed_positions(Undeformed_midplane_pt);

 //Reorder the elements, since I know what's best for them....
 solid_mesh_pt()->element_reorder();

 //Apply boundary conditions to the ends of the tube
 unsigned n_ends = solid_mesh_pt()->nboundary_node(1);
 //Loop over the node
 for(unsigned i=0;i<n_ends;i++)
  {
   //Pin in the axial direction (prevents rigid body motions)
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(2);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(2);
   //Derived conditions
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(2,2);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(2,2);

   //------------------CLAMPING CONDITIONS----------------------
   //------Pin positions in the transverse directions-----------
   // Comment these out to get the ring case
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(0);
   //Derived conditions
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(2,0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(2,0);

   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(1);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1);
   //Derived conditions
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(2,1);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(2,1);
   //----------------------------------------------------------

   // Set the axial gradients of the transverse coordinates to be
   // zero --- need to be enforced for ring or tube buckling
   //Pin dx/dz and dy/dz
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(1,0);
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(1,1);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1,0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1,1);
   //Derived conditions
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(3,0);
   solid_mesh_pt()->boundary_node_pt(1,i)->pin_position(3,1);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(3,0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(3,1);
 }

 //Now loop over the sides and apply symmetry conditions
 unsigned n_side = solid_mesh_pt()->nboundary_node(0);
 for(unsigned i=0;i<n_side;i++)
  {
   //At the side where theta is 0, pin in the y direction
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(1);
   //Derived condition
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(1,1);
   //Pin dx/dtheta and dz/dtheta
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(2,0);
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(2,2);
   //Pin the mixed derivative
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(3,0);
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(3,2);

   //At the side when theta is 0.5pi  pin in the x direction
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(0);
   //Derived condition
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(1,0);
   //Pin dy/dtheta and dz/dtheta
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(2,1);
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(2,2);
   //Pin the mixed derivative
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(3,1);
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(3,2);

   //Set an initial kick to make sure that we hop onto the
   //non-axisymmetric branch
   if((i>1) && (i<n_side-1))
    {
     solid_mesh_pt()->boundary_node_pt(0,i)->x(0) += 0.05;
     solid_mesh_pt()->boundary_node_pt(2,i)->x(1) -= 0.1;
    }
  }


 // Setup displacement control
 //---------------------------



//  //Setup displacement control
//  //Fix the displacement at the mid-point of the tube in the "vertical"
//  //(y) direction.
//  //Set the displacement control element (located halfway along the tube)
// Disp_ctl_element_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(3*Ny-1));
//  //The midpoint of the tube is located exactly half-way along the element
//  Vector<double> s(2);  s[0] = 1.0; s[1] = 0.0; //s[1] = 0.5
//  //Fix the displacement at this point in the y (1) direction
//  Disp_ctl_element_pt->fix_displacement_for_displacement_control(s,1);
//  //Set the pointer to the prescribed position
//  Disp_ctl_element_pt->prescribed_position_pt() = &Prescribed_y; 
 


 // Choose element in which displacement control is applied: This
 // one is located halfway along the tube)
 SolidFiniteElement* controlled_element_pt=
  dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(3*ny-1));
 
 // Fix the displacement in the y (1) direction...
 unsigned controlled_direction=1;

 // Local coordinate of the control point within the controlled element
 Vector<double> s_displ_control(2);
 s_displ_control[0]=1.0;
 s_displ_control[1]=0.0;
 
 // Pointer to displacement control element
 DisplacementControlElement* displ_control_el_pt;
 
 // Build displacement control element
 displ_control_el_pt=
  new DisplacementControlElement(controlled_element_pt,
                                 s_displ_control,
                                 controlled_direction,
                                 &Global_Physical_Variables::Prescribed_y);
 
 // The constructor of the  DisplacementControlElement has created
 // a new Data object whose one-and-only value contains the
 // adjustable load: Use this Data object in the load function:
 Global_Physical_Variables::Pext_data_pt=displ_control_el_pt->
  displacement_control_load_pt();
 
 // Add the displacement-control element to the mesh
 Displacement_control_mesh_pt = new Mesh;
 Displacement_control_mesh_pt->add_element_pt(displ_control_el_pt); 

 //This mesh must be retained as halo on all processors
#ifdef OOMPH_HAS_MPI
 Displacement_control_mesh_pt->set_keep_all_elements_as_halos();
#endif


 // Complete build of shell elements
 //---------------------------------

 //Find number of shell elements in mesh
 unsigned n_element = nx*ny;

 //Explicit pointer to first element in the mesh
 ELEMENT* first_el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(0));
 
 //Loop over the elements 
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to a shell element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(e));

   //Set the load function
   el_pt->load_vector_fct_pt() = & Global_Physical_Variables::press_load;

   //Set the undeformed surface
   el_pt->undeformed_midplane_pt() = Undeformed_midplane_pt;

   //The external pressure is external data for all elements
   el_pt->add_external_data(Global_Physical_Variables::Pext_data_pt);
   
   //Pre-compute the second derivatives wrt Lagrangian coordinates
   //for the first element only
   if(e==0)
    {
     el_pt->pre_compute_d2shape_lagrangian_at_knots();
    }
   //Otherwise set the values to be the same as those in the first element
   //this is OK because the Lagrangian mesh is uniform.
   //
   else
    {
     el_pt->set_dshape_lagrangian_stored_from_element(first_el_pt);
    }
  }

 //Set pointers to two trace nodes, used for output
 //in the middle of the shell so available on both processors!
 Trace_node_pt = solid_mesh_pt()->finite_element_pt(2*ny-1)->node_pt(3);
 Trace_node2_pt = solid_mesh_pt()->finite_element_pt(ny)->node_pt(1);

 //Now add the two submeshes
 this->add_sub_mesh(Solid_mesh_pt);
 this->add_sub_mesh(Displacement_control_mesh_pt);
 this->build_global_mesh();

 // Do equation numbering
 cout << std::endl;
 cout << "------------------DISPLACEMENT CONTROL--------------------" 
      << std::endl;
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 cout << "----------------------------------------------------------" 
      << std::endl;
 cout << std::endl;

}


//================================================================
// /Define the solve function, disp ctl and then continuation
//================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::solve()
{

 //Increase the maximum number of Newton iterations.
 //Finding the first buckled solution requires a large(ish) number
 //of Newton steps -- shells are just a bit twitchy
 Max_newton_iterations = 40;
 
 //Open an output trace file
 char filename[100];
 sprintf(filename,"trace_on_proc%i.dat",this->communicator_pt()->my_rank());
 ofstream trace(filename);

 //Change in displacemenet
 double disp_incr = 0.05;

 //Gradually compress the tube by decreasing the value of the prescribed
 //position from 1 to zero in steps of 0.05 initially and then 0.1
 for(unsigned i=1;i<13;i++)
  {
   //By the time we reach the second time round increase the incremenet
   if(i==3) {disp_incr = 0.1;}
   //Reduce prescribed y by our chosen increment
   Global_Physical_Variables::Prescribed_y -= disp_incr;

   cout << std::endl << "Increasing displacement: Prescribed_y is " 
        << Global_Physical_Variables::Prescribed_y << std::endl;

   // Solve
   newton_solve();
   
   //Output the pressure
   trace << Global_Physical_Variables::external_pressure()/(pow(0.05,3)/12.0)
         << " "
         //Position of first trace node
         << Trace_node_pt->x(0) << " " << Trace_node_pt->x(1) << " " 
          //Position of second trace node
         << Trace_node2_pt->x(0) << " " << Trace_node2_pt->x(1) << std::endl;
  }

 //Close the trace file
 trace.close();
 
 //Output the tube shape in the most strongly collapsed configuration
 sprintf(filename,"final_shape_on_proc%i.dat",
         this->communicator_pt()->my_rank());
 ofstream file(filename);
 solid_mesh_pt()->output(file,5);
 file.close();


 //Switch from displacement control to arc-length continuation and
 //trace back up the solution branch

 //Now pin the external pressure
 Global_Physical_Variables::Pext_data_pt->pin(0);

 //Re-assign the equation numbers
 cout << std::endl;
 cout << "-----------------ARC-LENGTH CONTINUATION --------------" 
      << std::endl;
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 cout << "-------------------------------------------------------" 
      << std::endl;
 cout << std::endl;

 //Set the maximum number of Newton iterations to something more reasonable
 Max_newton_iterations=6;

 //Set the desired number of Newton iterations per arc-length step
 Desired_newton_iterations_ds=3;

 //Set the proportion of the arc length
 Desired_proportion_of_arc_length = 0.2;

 //Set an initial value for the step size
 double ds = -0.5;

 //Open a different trace file
 sprintf(filename,"trace_disp_on_proc%i.dat",
         this->communicator_pt()->my_rank());
 trace.open(filename);
 //Take fifteen continuation steps
 for(unsigned i=0;i<15;i++)
   {
    ds = arc_length_step_solve(
     Global_Physical_Variables::Pext_data_pt,0,ds);
    
    //Output the pressure
    trace << Global_Physical_Variables::external_pressure()/(pow(0.05,3)/12.0)
          << " "
          //Position of first trace node
          << Trace_node_pt->x(0) << " " << Trace_node_pt->x(1) << " " 
          //Position of second trace node
          << Trace_node2_pt->x(0) << " " << Trace_node2_pt->x(1) << std::endl;
   }

 //Close the trace file
 trace.close();

}


//====================================================================
/// Driver
//====================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif
 { 
  //SuperLUSolver::Suppress_incorrect_rhs_distribution_warning_in_resolve=true;;

 //Length of domain
 double L = 10.0;
 double L_phi=0.5*MathematicalConstants::Pi;

 //Set up the problem
 ShellProblem<StorableShapeSolidElement<DiagHermiteShellElement> > 
  problem(5,5,L,L_phi);

 //Let's just be crazy and distribut it
#ifdef OOMPH_HAS_MPI
  //Set up a dummy partition
  unsigned n_element = problem.mesh_pt()->nelement();
  Vector<unsigned> element_partition(n_element);
  for(unsigned e=0;e<n_element/2;e++) {element_partition[e]=0;}
  for(unsigned e=n_element/2;e<n_element;e++) {element_partition[e]=1;}
  
  DocInfo mesh_doc_info;
  bool report_stats=true;
  mesh_doc_info.set_directory("RESLT_MESH");
  problem.distribute(element_partition,mesh_doc_info,report_stats);
  problem.check_halo_schemes(mesh_doc_info);
#endif

 //Solve the problem
 problem.solve();
 }
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
}







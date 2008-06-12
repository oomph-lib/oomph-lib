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
//Driver function for a simple test shell problem:
//Calculate the deformation of an elastic tube approximated
//using Kirchoff--Love shell theory using various boundary conditions


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

 /// \short Pointer to pressure load (stored in Data so it can 
 /// become an unknown in the problem when displacement control is used
 Data* Pext_data_pt;

 /// Perturbation pressure
 double Pcos=1.0;


 /// \short Return a reference to the external pressure 
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
    load[i] = (external_pressure() - 
               Pcos*pow(0.05,3)/12.0*cos(2.0*xi[1]))*N[i];
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
 
 /// \short In all elastic problems, the nodes must be assigned an undeformed,
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


//======================================================================
//Problem class to solve the deformation of an elastic tube
//=====================================================================
template<class ELEMENT>
class ShellProblem : public Problem
{

public:

 /// Constructor
 ShellProblem(const unsigned &nx, const unsigned &ny, 
              const double &lx, const double &ly,
              const unsigned& problem_id);

 /// Overload Access function for the mesh
 ShellMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<ShellMesh<ELEMENT>*>(Problem::mesh_pt());}

 /// Actions after solve empty
 void actions_after_newton_solve() {}

 /// Actions before solve empty
 void actions_before_newton_solve() {}
 
 /// Create clamping bc face elements on the b-th boundary of the Mesh
 void create_bc_elements(const unsigned& b);

 /// Do parameter study
 void solve(string& dir_name);

private:

 /// Pointer to GeomObject that specifies the undeformed midplane
 GeomObject* Undeformed_midplane_pt;

 /// First trace node
 Node* Trace_node_pt;

 /// Second trace node
 Node* Trace_node2_pt;

 /// Number of shell elements
 unsigned Nshell;

};



//======================================================================
/// Constructor
//======================================================================
template<class ELEMENT>
ShellProblem<ELEMENT>::ShellProblem(const unsigned &nx, const unsigned &ny, 
                                    const double &lx, const double &ly,
                                    const unsigned& problem_id)
{
 //Create the undeformed midplane object
 Undeformed_midplane_pt = new EllipticalTube(1.0,1.0);

 //Now create the mesh
 Problem::mesh_pt() = new ShellMesh<ELEMENT>(nx,ny,lx,ly); 

 // Store number of genuine shell elements
 Nshell=mesh_pt()->nelement();

 //Set the undeformed positions in the mesh
 mesh_pt()->assign_undeformed_positions(Undeformed_midplane_pt);

 //Reorder the elements, since I know what's best for them....
 mesh_pt()->element_reorder();


 
 //Loop over the nodes at the ends of the tube and pin the positions
 //-----------------------------------------------------------------
 // (All constraints apply at all nodes on boundaries 1 and 3)
 //-----------------------------------------------------------
 unsigned n_node_at_ends = mesh_pt()->nboundary_node(1);
 for(unsigned j=0;j<n_node_at_ends;j++)
  {

   // Loop over all three displacement components
   for (unsigned i=0;i<3;i++)
    {
     // Pin positions for all s_1 (along tube circumference)
     mesh_pt()->boundary_node_pt(1,j)->pin_position(i);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(i);

     // ...implying that dr/ds_1=0, too:  [pin_position() takes type (2) and 
     // direction (i)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(2,i);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(2,i);
    }
  
   
   // Use symmetry conditions at ends?
   //--------------------------------
   if (problem_id==2)
    {
     // Pin the derivative (into the shell) of the x-position  for all 
     // s_1 (along tube circumference)[pin_position() takes type (1) 
     // and direction(0)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(1,0);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(1,0);
     
     // ...implying that d^2x/ds_0 ds_1=0, too:  [pin_position() takes 
     // type (3) and direction (0)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(3,0);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(3,0);
     
     // Pin the derivative (into the shell) of the y-position  for all 
     // s_1 (along tube circumference)[pin_position() takes type (1) 
     // and direction(1)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(1,1);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(1,1);
     
     // ...implying that d^2y/ds_0 ds_1=0, too:  [pin_position() takes 
     // type (3) and direction (1)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(3,1);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(3,1);
    }
 }


 //Now loop over the sides and apply symmetry conditions
 //-----------------------------------------------------
 unsigned n_node_at_side = mesh_pt()->nboundary_node(0);
 for(unsigned j=0;j<n_node_at_side;j++)
  {

   //At the side where theta is 0 (boundary 0), pin in the y displacement
   //--------------------------------------------------------------------

   // y=0 for all s_0 (along the tube)...
   mesh_pt()->boundary_node_pt(0,j)->pin_position(1);

   // ...implying that dy/ds_0=0 too: [pin_position() takes type (1) and 
   // direction(1)]
   mesh_pt()->boundary_node_pt(0,j)->pin_position(1,1);


   //At the side where theta is 0 (boundary 0), apply symm cond for x and z
   //----------------------------------------------------------------------

   // d{x,z}/ds_1=0 for all s_0 (along the tube) because of symmetry...
   // [pin_position() takes type (2) and direction {0,2}]
   mesh_pt()->boundary_node_pt(0,j)->pin_position(2,0);
   mesh_pt()->boundary_node_pt(0,j)->pin_position(2,2);

   // ...implying that d^2{x,z}/ds_0 ds_1=0 too: [pin_position() takes 
   // type (3) and direction {0,2}]
   mesh_pt()->boundary_node_pt(0,j)->pin_position(3,0);
   mesh_pt()->boundary_node_pt(0,j)->pin_position(3,2);


   //At the side where theta is pi/2 (boundary 2), pin in the x displacement
   //--------------------------------------------------------------------

   // x=0 for all s_0 (along the tube)...
   mesh_pt()->boundary_node_pt(2,j)->pin_position(0);

   // ...implying that dx/ds_0=0 too: [pin_position() takes type (1) and 
   // direction(0)]
   mesh_pt()->boundary_node_pt(2,j)->pin_position(1,0);


   //At the side where theta is pi/2 (boundary 2), apply symm cond for y and z
   //----------------------------------------------------------------------

   // d{y,z}/ds_1=0 for all s_0 (along the tube) because of symmetry...
   // [pin_position() takes type (2) and direction {1,2}]
   mesh_pt()->boundary_node_pt(2,j)->pin_position(2,1);
   mesh_pt()->boundary_node_pt(2,j)->pin_position(2,2);

   // ...implying that d^2{y,z}/ds_0 ds_1=0 too: [pin_position() takes 
   // type (3) and direction {1,2}]
   mesh_pt()->boundary_node_pt(2,j)->pin_position(3,1);
   mesh_pt()->boundary_node_pt(2,j)->pin_position(3,2);

  }


 // Setup displacement control
 //---------------------------

 // Choose element in which displacement control is applied: This
 // one is located about halfway along the tube -- remember that
 // we've renumbered the elements!
 unsigned nel_ctrl=0;
 Vector<double> s_displ_control(2);

 // Even/odd number of elements in axial direction
 if (nx%2==1)
  {
   nel_ctrl=unsigned(floor(0.5*double(nx))+1.0)*ny-1;
   s_displ_control[0]=0.0;
   s_displ_control[1]=1.0;
  }
 else
  {
   nel_ctrl=unsigned(floor(0.5*double(nx))+1.0)*ny-1;
   s_displ_control[0]=-1.0;
   s_displ_control[1]=1.0;
  }

 // Controlled element
 SolidFiniteElement* controlled_element_pt=
  dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(nel_ctrl));
 

 // Fix the displacement in the y (1) direction...
 unsigned controlled_direction=1;

 // Pointer to displacement control element
 DisplacementControlElement* displ_control_el_pt;
 
 // Build displacement control element
 displ_control_el_pt=
  new DisplacementControlElement(controlled_element_pt,
                                 s_displ_control,
                                 controlled_direction,
                                 &Global_Physical_Variables::Prescribed_y);
 

 // Doc control point
 Vector<double> xi(2);
 Vector<double> x(3);
 controlled_element_pt->interpolated_xi(s_displ_control,xi);
 controlled_element_pt->interpolated_x(s_displ_control,x);
 std::cout << std::endl;
 std::cout << "Controlled element: " << nel_ctrl << std::endl;
 std::cout << "Displacement control applied at xi = (" 
           << xi[0] << ", " << xi[1] << ")" << std::endl;
 std::cout << "Corresponding to                x  = (" 
           << x[0] << ", " << x[1] << ", " << x[2] << ")" << std::endl;


 // The constructor of the  DisplacementControlElement has created
 // a new Data object whose one-and-only value contains the
 // adjustable load: Use this Data object in the load function:
 Global_Physical_Variables::Pext_data_pt=displ_control_el_pt->
  displacement_control_load_pt();
 
 // Add the displacement-control element to the mesh
 mesh_pt()->add_element_pt(displ_control_el_pt); 

 

 // Complete build of shell elements
 //---------------------------------

 //Find number of shell elements in mesh
 unsigned n_element = nx*ny;

 //Explicit pointer to first element in the mesh
 ELEMENT* first_el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));
 
 //Loop over the elements 
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to a shell element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

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
   else
    {
     el_pt->set_dshape_lagrangian_stored_from_element(first_el_pt);
    }
  }

 //Set pointers to two trace nodes, used for output
 Trace_node_pt = mesh_pt()->finite_element_pt(2*ny-1)->node_pt(3);
 Trace_node2_pt = mesh_pt()->finite_element_pt(ny)->node_pt(1);

 // Attach boundary condition elements that enforce clamping condition
 if (problem_id==0)
  {
   create_bc_elements(1);
   create_bc_elements(3);
  }


 // Do equation numbering
 cout << std::endl;
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 cout << std::endl;

}


//=======================================================================
/// Create clamping bc face elements on the b-th boundary of the Mesh
//=======================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::create_bc_elements(const unsigned& b)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = mesh_pt()->nboundary_element(b);

 std::cout << "Creating " << n_element 
           << " boundary condition elements"  << std::endl;


 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    mesh_pt()->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = mesh_pt()->face_index_at_boundary(b,e);

   // Build the corresponding bc element
   ClampedHermiteShellBoundaryConditionElement* flux_element_pt = new 
   ClampedHermiteShellBoundaryConditionElement(bulk_elem_pt,face_index);

   //Add the bc element to the mesh
   mesh_pt()->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

}


//================================================================
// /Define the solve function, disp ctl and then continuation
//================================================================
template<class ELEMENT>
void ShellProblem<ELEMENT>::solve(string& dir_name)
{
 ofstream some_file;
 char filename[100];
 
 //Increase the maximum number of Newton iterations.
 //Finding the first buckled solution requires a large(ish) number
 //of Newton steps -- shells are just a bit twitchy
 Max_newton_iterations = 40;
 Max_residuals=1.0e6;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory(dir_name);

  //Open an output trace file
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 ofstream trace_file(filename);
 trace_file << "VARIABLES=\"p_e_x_t\",\"R_1\",\"R_2\"" << std::endl;
 trace_file << "ZONE" << std::endl;

 //Gradually compress the tube by decreasing the value of the prescribed
 //position 
 for(unsigned i=1;i<11;i++)
  {

   Global_Physical_Variables::Prescribed_y -= 0.1;

   // Solve
   newton_solve();   

   //Output the pressure (on the bending scale)
   trace_file 
    << Global_Physical_Variables::external_pressure()/(pow(0.05,3)/12.0) << " "
    << Trace_node_pt->x(1) << " " 
    << Trace_node2_pt->x(0) << std::endl;
   
   
   //Output the tube shape 
   sprintf(filename,"%s/shell%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   for (unsigned e=0;e<Nshell;e++)
    {
     mesh_pt()->finite_element_pt(e)->output(some_file,15);
    }
   some_file.close();

   //Output the Lagrange multipliers
   sprintf(filename,"%s/lagrange_multiplier%i.dat",
           doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   unsigned n_el=mesh_pt()->nelement();
   for (unsigned e=Nshell+1;e<n_el;e++)
    {
     mesh_pt()->finite_element_pt(e)->output(some_file,15);
    }
   some_file.close();

   // Increment counter for output files
   doc_info.number()++;

   // Reset perturbation
   Global_Physical_Variables::Pcos=0.0;
  }

 //Close the trace file
 trace_file.close();
 
}


//====================================================================
/// Driver
//====================================================================
int main(int argc, char* argv[])
{

 //Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // hierher
 //FiniteElement::Accept_negative_jacobian=true;

 //Length of domain
 double L = 10.0;
 double L_phi=0.5*MathematicalConstants::Pi;

 // Choose problem type: 0: Clamped; 1: Pinned; 2: Pinned & symmetry
 unsigned problem_id=0;
 string dir_name;



 if (CommandLineArgs::Argc==2)
  {
   if (atoi(CommandLineArgs::Argv[1])==0)
    {
     std::cout << "Doing clamped case " << std::endl;
     problem_id=0;
     dir_name="RESLT_clamped";
    }
   else if (atoi(CommandLineArgs::Argv[1])==1)
    {
     std::cout << "Doing pinned case " << std::endl;
     problem_id=1;
     dir_name="RESLT_pinned";
    }
   else if (atoi(CommandLineArgs::Argv[1])==2)
    {
     std::cout << "Doing pinned periodic case " << std::endl;
     problem_id=2;
     dir_name="RESLT_pinned_periodic";
    }
   else
    {
     std::cout << "Wrong command line argument -- using default (0) instead" 
               << std::endl;
     problem_id=0;
     dir_name="RESLT_clamped";
    }
  }
 else
  {
   std::cout 
    << "No (or more than one) input argument(s) -- using default (0) instead" 
    << std::endl;  
     problem_id=0;
     dir_name="RESLT_clamped";  
  }


 // Run with pinned boundary conditions
 //------------------------------------
 if (problem_id==1)
 {
  //Set up the problem
  ShellProblem<StorableShapeSolidElement<DiagHermiteShellElement> >
   problem(5,3,L,L_phi,problem_id);
  
  //Solve the problem
  problem.solve(dir_name);
 }
 // Run with proper clamped boundary conditions
 //--------------------------------------------
 else if (problem_id==0)
  {
   {
    //Set up the problem
    ShellProblem<StorableShapeSolidElement<DiagHermiteShellElement> >
     problem(5,3,L,L_phi,problem_id);
   
    //Solve the problem
    problem.solve(dir_name);   
   }
//    {    
//     //Set up the problem
//     ShellProblem<StorableShapeSolidElement<DiagHermiteShellElement> >
//      problem(5,3,L,L_phi,problem_id);
    
//     //Solve the problem
//     problem.solve(dir_name);   
//    }
  }
 // Run with pinned periodic boundary conditions
 //---------------------------------------------
 else if (problem_id==2)
  {
   //Set up the problem
   ShellProblem<StorableShapeSolidElement<DiagHermiteShellElement> >
    problem(5,3,L,L_phi,problem_id);
   //Solve the problem
   problem.solve(dir_name);
   
  }
 else
  {
   std::cout  << "Wrong problem id " << problem_id << std::endl;
  }



}







//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
//Driver function for a simple test plate problem:
//Calculate the deformation of a flat pre-stressed approximated
//using Kirchoff--Love shell theory


//Include files from the finite-element library
#include "generic.h"
#include "shell.h"
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//===========================================================
/// Flat plate in x-y plane
//===========================================================
class FlatPlate : public GeomObject
{
public:

 /// Constructor: 
 FlatPlate() :  GeomObject(2,3) {}
 
 /// Broken copy constructor
 FlatPlate(const FlatPlate& node) 
  { 
   BrokenCopy::broken_copy("FlatPlate");
  } 
 
 /// Broken assignment operator
 void operator=(const FlatPlate&) 
  {
   BrokenCopy::broken_assign("FlatPlate");
  }

 /// Position vector
 void position(const Vector<double>& zeta, Vector<double>& r)const
  {
   r[0]=zeta[0];
   r[1]=zeta[1];
   r[2]=0.0;
  }


 /// Position vector (dummy unsteady version returns steady version)
 void position(const unsigned& t, 
               const Vector<double>& zeta,
               Vector<double>& r)const
  {
   position(zeta,r);
  }

 /// How many items of Data does the shape of the object depend on?
 virtual unsigned ngeom_data() const
  {
   return 0;
  }

 /// Position Vector and 1st and 2nd derivs w.r.t. zeta.
 void d2position(const Vector<double>& zeta, Vector<double>& r,
                 DenseMatrix<double> &drdzeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   //Flat plate
   r[0] = zeta[0];
   r[1] = zeta[1];
   r[2] = 0.0;

   //Do the zeta0 derivatives
   drdzeta(0,0) = 1.0;
   drdzeta(0,1) = 0.0;
   drdzeta(0,2) = 0.0;
   
   //Do the zeta 1 derivatives
   drdzeta(1,0) = 0.0;
   drdzeta(1,1) = 1.0;
   drdzeta(1,2) = 0.0;
   
   //Now let's do the second derivatives
   ddrdzeta(0,0,0) = 0.0;
   ddrdzeta(0,0,1) = 0.0;
   ddrdzeta(0,0,2) = 0.0;
   
   ddrdzeta(1,1,0) = 0.0;
   ddrdzeta(1,1,1) = 0.0;
   ddrdzeta(1,1,2) = 0.0;
   
   //Mixed derivatives
   ddrdzeta(0,1,0) = 0.0;
   ddrdzeta(0,1,1) = 0.0;
   ddrdzeta(0,1,2) = 0.0;
  }

};


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//========================================================================
/// A 2D solid mesh for (topologically) flat 2D plates.
/// The plate is represented by two Lagrangian coordinates that correspond
/// to x and y in the undeformed position. The required mesh is therefore a
/// 2D mesh and is therefore inherited from the generic RectangularQuadMesh
//=======================================================================
template <class ELEMENT>
class FlatPlateMesh : 
 public virtual RectangularQuadMesh<ELEMENT>,
 public virtual SolidMesh
{
public:

 /// Constructor for the mesh
 FlatPlateMesh(const unsigned &nx,
               const unsigned &ny,
               const double &lx,
               const double &ly,
               TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper);
 
 /// In all elastic problems, the nodes must be assigned an undeformed,
 /// or reference, position, corresponding to the stress-free state
 /// of the elastic body. This function assigns the undeformed position
 /// for the nodes on the flat plate
 void assign_undeformed_positions(GeomObject* const &undeformed_midplane_pt);

};


//=======================================================================
/// Mesh constructor
/// Argument list:
/// nx  : number of elements in the axial direction
/// ny  : number of elements in the azimuthal direction
/// lx  : length in the x direction
/// ly  : length in the y direction
//=======================================================================
template <class ELEMENT>
FlatPlateMesh<ELEMENT>::FlatPlateMesh(
 const unsigned &nx,
 const unsigned &ny,
 const double &lx,
 const double &ly,
 TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt)
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
 //Hermite-type elements

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
void FlatPlateMesh<ELEMENT>::assign_undeformed_positions(
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


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////



//======================================================================== 
/// Global variables that represent physical properties
//======================================================================== 
namespace Global_Physical_Variables
{

 /// Prescribed position of control point
 double Prescribed_z = 0.0;

 /// Pointer to pressure load (stored in Data so it can 
 /// become an unknown in the problem when displacement control is used
 Data* Pext_data_pt;

 /// 2nd Piola Kirchhoff pre-stress
 double Sigma0=0.1;

 /// Poisson's ratio
 double Nu=0.3;

 /// Wall thickness
 double H=0.01;

 /// Return a reference to the external pressure 
 /// load on the elastic tube.
 double external_pressure() 
 {return (*Pext_data_pt->value_pt(0)); } // *pow(0.05,3)/12.0;}


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


//======================================================================
//Problem class to solve the deformation of an elastic plate
//=====================================================================
template<class ELEMENT>
class PlateProblem : public Problem
{

public:

 /// Constructor
 PlateProblem(const unsigned &nx, const unsigned &ny, 
              const double &lx, const double &ly);

 /// Overload Access function for the mesh
 FlatPlateMesh<ELEMENT>* mesh_pt() 
  {return 
    dynamic_cast<FlatPlateMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Actions after solve empty
 void actions_after_newton_solve() {}

 /// Actions before solve empty
 void actions_before_newton_solve() {}
 
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
PlateProblem<ELEMENT>::PlateProblem(const unsigned &nx, const unsigned &ny, 
                                    const double &lx, const double &ly)
{
 //Create the undeformed midplane object
 Undeformed_midplane_pt = new FlatPlate;

 //Now create the mesh
 Problem::mesh_pt() = new FlatPlateMesh<ELEMENT>(nx,ny,lx,ly); 

 // Store number of genuine shell elements
 Nshell=mesh_pt()->nelement();

 //Set the undeformed positions in the mesh
 mesh_pt()->assign_undeformed_positions(Undeformed_midplane_pt);

 //Reorder the elements, since I know what's best for them....
 mesh_pt()->element_reorder();

 // Loop over all nodes and suppress displacements in y (=1) direction
 //-------------------------------------------------------------------
 unsigned n_node = mesh_pt()->nnode();
 for (unsigned j=0;j<n_node;j++)
  {
   SolidNode* nod_pt=mesh_pt()->node_pt(j);

   //pin y-displacement
   nod_pt->pin_position(1);

   // ...for all s_0:
   nod_pt->pin_position(1,1);

   // ...for all s_1:
   nod_pt->pin_position(2,1);

   // ...for all s_0, s_1:
   nod_pt->pin_position(3,1);
  }


 //Loop over the nodes at the ends of the plate and pin the positions
 //-----------------------------------------------------------------
 // (All constraints apply at all nodes on boundaries 1 and 3)
 //-----------------------------------------------------------
 unsigned n_node_at_ends = mesh_pt()->nboundary_node(1);
 for(unsigned j=0;j<n_node_at_ends;j++)
  {

   // Loop over all three displacement components
   for (unsigned i=0;i<3;i++)
    {
     // Pin positions for all s_1 
     mesh_pt()->boundary_node_pt(1,j)->pin_position(i);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(i);

     // ...implying that dr/ds_1=0, too:  [pin_position() takes type (2) and 
     // direction (i)]
     mesh_pt()->boundary_node_pt(1,j)->pin_position(2,i);
     mesh_pt()->boundary_node_pt(3,j)->pin_position(2,i);
    }
   
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
 

 // Fix the displacement in the z (2) direction...
 unsigned controlled_direction=2;

 // Pointer to displacement control element
 DisplacementControlElement* displ_control_el_pt;
 
 // Build displacement control element
 displ_control_el_pt=
  new DisplacementControlElement(controlled_element_pt,
                                 s_displ_control,
                                 controlled_direction,
                                 &Global_Physical_Variables::Prescribed_z);
 

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

   // Set Nu
   el_pt->nu_pt() = &Global_Physical_Variables::Nu;

   // Set H
   el_pt->h_pt() = &Global_Physical_Variables::H;

   //The external pressure is external data for all elements
   el_pt->add_external_data(Global_Physical_Variables::Pext_data_pt);
   
   // set pre-stress
   el_pt->set_prestress_pt(0,0,&Global_Physical_Variables::Sigma0); 

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

 // Do equation numbering
 cout << std::endl;
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 cout << std::endl;

}




//====================================================================
/// Driver
//====================================================================
int main(int argc, char* argv[])
{
 
 //Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Length of domain
 double L_x=10.0;
 double L_y=1.0;

 //Set up the problem
 PlateProblem<StorableShapeSolidElement<DiagHermiteShellElement> >
  problem(10,3,L_x,L_y);

 ofstream some_file;
 char filename[100];
 
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 //Open an output trace file
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 ofstream trace_file(filename);
 trace_file << "VARIABLES=\"p_e_x_t\",\"z_ctrl\"" << std::endl;
 trace_file << "ZONE" << std::endl;
 
 //Gradually deflect the plate by decreasing the value of the prescribed
 //position 
 unsigned nstep=2;
 for(unsigned i=0;i<nstep;i++)
  {
   
   // Increase displacement
   Global_Physical_Variables::Prescribed_z -= 0.1;
   
   // Solve
   problem.newton_solve();   
   
   
   // Calculate exact solution for `string under tension' (applicable for
   // small wall thickness and pinned ends)
   
   // The tangent of the angle beta
   double tanbeta =-2.0*Global_Physical_Variables::Prescribed_z/L_x;
   
   double exact_pressure = 0.0;
   //If the plate has deformed, calculate the pressure required
   if(tanbeta!=0)
    {
     
     //Calculate the opening angle alpha
     double alpha = 2.0*atan(2.0*tanbeta/(1.0-tanbeta*tanbeta));

      // Jump back onto the main branch if alpha>180 degrees
     if (alpha<0) alpha+=2.0*MathematicalConstants::Pi;
     
     // Green strain:
     double gamma=0.5*(0.25*alpha*alpha/(sin(0.5*alpha)*sin(0.5*alpha))-1.0);
     
     //Calculate the exact pressure
     exact_pressure=Global_Physical_Variables::H*
      (Global_Physical_Variables::Sigma0+gamma)*alpha/L_x;
    } 
   

   //Output the pressure 
   trace_file 
    << Global_Physical_Variables::external_pressure() << " " ///(pow(0.05,3)/12.0) << " "
    << Global_Physical_Variables::Prescribed_z  << " " 
    << exact_pressure << " "
    << std::endl;
   
   
   //Output the shape 
   sprintf(filename,"%s/plate%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   problem.mesh_pt()->output(some_file,15);
//    for (unsigned e=0;e<nelem;e++)
//     {
//      mesh_pt()->finite_element_pt(e)->output(some_file,15);
//     }
   some_file.close();
   
   // Increment counter for output files
   doc_info.number()++;
   
  }

 //Close the trace file
 trace_file.close();
 
}







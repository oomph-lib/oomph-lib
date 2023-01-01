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
// Driver for 2-D Channel with changing width, using Taylor Hood 
// and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=========================================================================
/// Geometric object representing a wavy wall, parametrised by
///  \f[ x = \zeta \f]
///  \f[ y = A (1 - \cos(2\pi\zeta/L) )\f]
//=========================================================================
class WavyWall : public GeomObject
{

public:

 /// Constructor:  Pass wavelength and amplitude
 WavyWall(const double& l, const double& amplitude)
  : GeomObject(1,2), L(l), A(amplitude) {}

 /// Position vector to wavy wall.
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   r[0]=zeta[0];
   r[1]=A*(1.0-cos(2.0*MathematicalConstants::Pi*zeta[0]/L));
  }

protected:

 /// Wavelength
 double L;

 /// Amplitude 
 double A;


};


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//======================================================================
/// Simple spine mesh class derived from standard 2D mesh.
/// Vertical lines of nodes are located on spines.
//======================================================================
template <class ELEMENT>
class SimpleSpineMesh : public SimpleRectangularQuadMesh<ELEMENT >, 
                        public SpineMesh
{

public:

 /// Constructor: Pass number of elements in x-direction,
 /// number of elements in y-direction, length in x direction, initial 
 /// height of mesh, and pointer to timestepper (defaults to 
 /// Steady timestepper)
 SimpleSpineMesh(const unsigned &nx, 
                 const unsigned &ny,
                 const double &lx, 
                 const double &h, 
                 GeomObject* substrate_pt,
                 TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper);


 /// Function pointer to  function, h(x), that may be used 
 /// to specify the "height" of the domain, by assigning the function
 /// values to the spine heights.
 typedef double (*HeightFctPt)(const double& x);

 /// Access function: Pointer to height function
 HeightFctPt& height_fct_pt() 
  {
   return Height_fct_pt;
  }


 /// Height function -- this is called by update_spine_heights()
 /// when spine heights are assigned.
 double height_fct(const double& x) 
  {
   // Resolve function pointer if non-NULL
   if (Height_fct_pt==0)
    {
     return Default_height;
    }
   else
   {
    return Height_fct_pt(x);
   }
  }


 /// Update the spine heights according to the function specified
 /// by height_fct_pt().
 void update_spine_heights()
  {
   unsigned n_spine=nspine();
   for (unsigned i=0;i<n_spine;i++)
    {
     // Zero-th geometric parameter of the spine is the x-coordinate
     // of its basis
     double x=spine_pt(i)->geom_parameter(0);

     // Set spine height
     spine_pt(i)->height()=height_fct(x);
    }
  }


 /// General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific node update
 /// actions: Nodes are located along vertical "spines" that emanate
 /// from the "substrate" (the lower wall) specified by a
 /// GeomObject.
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {

   //Get fraction along the spine
   double w = spine_node_pt->fraction();
   
   // Get height of spine
   double h = spine_node_pt->spine_pt()->height();

   // Get x-position of spine origin (as vector)
   Vector<double> x_spine_origin(1);
   x_spine_origin[0]=spine_node_pt->spine_pt()->geom_parameter(0);
   
   // Get position vector to substrate (lower wall) for this spine
   Vector<double> spine_base(2);
   spine_node_pt->spine_pt()->geom_object_pt(0)->
    position(x_spine_origin,spine_base);
     
   //Update the nodal position
   spine_node_pt->x(0) = spine_base[0];
   spine_node_pt->x(1) = spine_base[1]+w*h;
  }


private:

 /// Default height
 double Default_height;

 /// Pointer to height function
 HeightFctPt Height_fct_pt;

 /// Pointer to GeomObject that specifies the "substrate" (the lower wall)
 GeomObject* Substrate_pt;

};

 

//===========================================================================
/// Constructor for spine 2D mesh: Pass number of elements in x-direction, 
/// number of elements in y-direction, axial length and height of layer, 
/// pointer to geometric object that specifies the substrate (the lower wall)
/// and pointer to timestepper (defaults to Static timestepper).
///
/// The mesh contains a layer of spine-ified fluid elements (of type ELEMENT;
/// e.g  SpineElement<QCrouzeixRaviartElement<2>)
//===========================================================================
template<class ELEMENT>
SimpleSpineMesh<ELEMENT>::SimpleSpineMesh(const unsigned &nx, 
                                          const unsigned &ny,
                                          const double &lx, 
                                          const double &h, 
                                          GeomObject* substrate_pt,
                                          TimeStepper* time_stepper_pt) :
 SimpleRectangularQuadMesh<ELEMENT>(nx,ny,lx,h,time_stepper_pt),
 Default_height(h), Height_fct_pt(0), Substrate_pt(substrate_pt)

{

 // We've already called the constructor for the SimpleRectangularQuadMesh
 // which sets up nodal positons/topology etc. Now attach the spine
 // information 


 // Allocate memory for the spines and fractions along spines

 // Read out number of linear points in the element
 unsigned np = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();
 
 // Number of spines
 unsigned nspine = (np-1)*nx+1;
 Spine_pt.reserve(nspine);

 // Set up x-increments between spine origins
 double x_lo=0.0;
 double dx=lx/double(nx);

 //FIRST SPINE
 //-----------

 //Element 0
 //Node 0
 //Assign the new spine with specified height
 Spine* new_spine_pt=new Spine(h);
 Spine_pt.push_back(new_spine_pt);

 // Get pointer to node
 SpineNode* nod_pt=element_node_pt(0,0);

 //Pass the pointer to the spine to the node
 nod_pt->spine_pt() = new_spine_pt;

 //Set the node's fraction along the spine
 nod_pt->fraction() = 0.0;

 // Set the pointer to the mesh that implements the update fct for this node
 nod_pt->spine_mesh_pt() = this; 

 // Set update fct id (not really needed here...)
 nod_pt->node_update_fct_id()=0;


 // Create vector that stores additional geometric parameters that
 // are required during the execution of the remesh function:
 // For our remesh function, we need the x-coordinate of the 
 // spine origin. 

 // Vector with enough storage for one geometic parameter
  Vector<double> parameters(1);

  // Store the x-coordinate in it
  parameters[0]=x_lo;

  // Pass the parameter to the spine
  new_spine_pt->set_geom_parameter(parameters);
     

 
  // The remesh function involving this spine only involves
  // a single geometric object: The substrate
  Vector<GeomObject*> geom_object_pt(1);
  geom_object_pt[0] = Substrate_pt;
  
  // Pass geom object(s) to spine
  new_spine_pt->set_geom_object_pt(geom_object_pt);
  

 //Loop vertically along the first spine

 //Loop over the elements 
 for(unsigned i=0;i<ny;i++)
  {
   //Loop over the vertical nodes, apart from the first
   for(unsigned l1=1;l1<np;l1++)
    {
     // Get pointer to node
     SpineNode* nod_pt=element_node_pt(i*nx,l1*np);

     //Pass the pointer to the spine to the node
     nod_pt->spine_pt() = new_spine_pt;

     //Set the fraction
     nod_pt->fraction()=(double(i)+double(l1)/double(np-1))/double(ny);

     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 

     // Set update fct id (not really needed here)
     nod_pt->node_update_fct_id()=0;

    }
  } // end loop over elements


 //LOOP OVER OTHER SPINES
 //----------------------

 //Now loop over the elements horizontally
 for(unsigned j=0;j<nx;j++)
  {
   //Loop over the nodes in the elements horizontally, ignoring 
   //the first column
   unsigned npmax=np;
   for(unsigned l2=1;l2<npmax;l2++)
    {
     //Assign the new spine with specified height
     new_spine_pt=new Spine(h);

     // Add to collection of spines
     Spine_pt.push_back(new_spine_pt);
     
     // Get the node
     SpineNode* nod_pt=element_node_pt(j,l2);

     //Pass the pointer to the spine to the node
     nod_pt->spine_pt() = new_spine_pt;

     //Set the fraction
     nod_pt->fraction() = 0.0;

     // Pointer to the mesh that implements the update fct
     nod_pt->spine_mesh_pt() = this; 

     // Set update fct id (not really needed here)
     nod_pt->node_update_fct_id()=0;

     // Create vector that stores additional geometric parameters that
     // are required during the execution of the remesh function:
     // For our remesh function, we need the x-coordinate of the 
     // spine origin. 
     
     // Vector with enough storage for one geometic parameter
     Vector<double> parameters(1);
     
     // Store the x-coordinate in it
     parameters[0]=x_lo+double(j)*dx + double(l2)*dx/double(np-1);
     
     // Pass the parameter to the spine
     new_spine_pt->set_geom_parameter(parameters);
     

     // The remesh function involving this spine only involves
     // a single geometric object: The substrate
     Vector<GeomObject*> geom_object_pt(1);
     geom_object_pt[0] = Substrate_pt;
     
     // Pass geom object(s) to spine
     new_spine_pt->set_geom_object_pt(geom_object_pt);
     
     //Loop vertically along the spine
     //Loop over the elements 
     for(unsigned i=0;i<ny;i++)
      {
       //Loop over the vertical nodes, apart from the first
       for(unsigned l1=1;l1<np;l1++)
        {
         // Get the node
         SpineNode* nod_pt=element_node_pt(i*nx+j,l1*np+l2);

         //Set the pointer to the spine
         nod_pt->spine_pt() = new_spine_pt;

         //Set the fraction
         nod_pt->fraction()=(double(i)+double(l1)/double(np-1))/double(ny);

         // Pointer to the mesh that implements the update fct
         nod_pt->spine_mesh_pt() = this; 

         // Set update fct id (not really needed here)
         nod_pt->node_update_fct_id()=0;
        }  
      }
    }
  }

} // end of constructor





/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100;

 /// Start of indented region
 double X_indent_start=0.5;

 /// Length of indented region
 double L=1.0;

 /// Total length of domain
 double L_total=4.0;

 /// Undeformed height of domain
 double H=1.0;

 /// Amplitude of indentation
 double A=-0.6;

 /// Height of domain
 double height(const double& x)
 {
  if ((x>X_indent_start)&&(x<(X_indent_start+L)))
   {
    return H + A*sin(MathematicalConstants::Pi*(x-X_indent_start)/L);
   }
  else
   {
    return H;
   }
 }

} // end_of_namespace





/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Channel flow through a non-uniform channel whose geometry is defined 
/// by a spine mesh.
//====================================================================
template<class ELEMENT>
class ChannelSpineFlowProblem : public Problem
{

public:

 /// Constructor
 ChannelSpineFlowProblem();
 
 /// Destructor: (empty)
 ~ChannelSpineFlowProblem(){}

 /// Update the problem specs before solve. Update the nodal
 /// positions
 void actions_before_newton_solve()
  { 
   // Update the mesh
   mesh_pt()->node_update();

  } // end_of_actions_before_newton_solve

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Overload access to mesh
 SimpleSpineMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<SimpleSpineMesh<ELEMENT>*>(Problem::mesh_pt());
  }


private:
 
 /// Width of channel
 double Ly;

 
}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for ChannelSpineFlow problem 
//========================================================================
template<class ELEMENT>
ChannelSpineFlowProblem<ELEMENT>::ChannelSpineFlowProblem()
{ 

 // Setup mesh

 // Domain length in x-direction
 double Lx=Global_Physical_Variables::L_total;
 
 // Domain length in y-direction
 Ly=Global_Physical_Variables::H;

 // # of elements in x-direction
 unsigned Nx=40;

 // # of elements in y-direction
 unsigned Ny=10;
 
 // Substrate (lower wall): A wavy wall
 GeomObject* substrate_pt=
  new WavyWall(Global_Physical_Variables::L_total,-0.2);

 // Build and assign mesh -- pass pointer to geometric object
 // that represents the sinusoidal bump on the upper wall
 Problem::mesh_pt() = new SimpleSpineMesh<ELEMENT>(Nx,Ny,Lx,Ly,substrate_pt);

 // Set function pointer for height function
 mesh_pt()->height_fct_pt()=&Global_Physical_Variables::height;

 // Update spine heights according to the function just specified
 mesh_pt()->update_spine_heights();

 // Update nodal positions
 mesh_pt()->node_update();

 // Pin all spine heights
 unsigned nspine=mesh_pt()->nspine();
 for (unsigned i=0;i<nspine;i++)
  {
   mesh_pt()->spine_pt(i)->spine_height_pt()->pin(0);
  }


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries, except on boundary 1
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     if (ibound!=1)
      {
       // Loop over values (u and v velocities)
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
        }
      }
     else
      {
       // Parallel outflow ==> no-slip
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
      }
    }
  } // end loop over boundaries

 
 // No slip on stationary upper and lower walls (boundaries 0 and 2) 
 // and parallel outflow (boundary 1) 
 for (unsigned ibound=0;ibound<num_bound-1;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     if (ibound!=1)
      {
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
        }
      }
     else
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
    }
  }


 // Setup parabolic inflow along boundary 3:
 unsigned ibound=3; 
 unsigned num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
   // Parallel, parabolic inflow
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,y*(Ly-y));
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
  }


 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass 
 // pointer to Reynolds number
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } 

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ChannelSpineFlowProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution



//==start_of_main======================================================
/// Driver for channel flow problem with spine mesh.
//=====================================================================
int main()
{

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 doc_info.number()=0;
 
 // Solve problem with Taylor Hood elements
 //---------------------------------------
 {
  //Build problem
  ChannelSpineFlowProblem<SpineElement<QTaylorHoodElement<2> > > problem;
  
  // Solve the problem
  problem.newton_solve();
  
  //Output solution
  problem.doc_solution(doc_info);

  // Step number
  doc_info.number()++;

 } // end of Taylor Hood elements
 
 
 // Solve problem with Crouzeix Raviart elements
 //--------------------------------------------
 {
  // Build problem
  ChannelSpineFlowProblem<SpineElement<QCrouzeixRaviartElement<2> > >
   problem;
  
  // Solve the problem with automatic adaptation
  problem.newton_solve();
  
  //Output solution
  problem.doc_solution(doc_info);
  // Step number
  doc_info.number()++;
  
 } // end of Crouzeix Raviart elements
      

     
} // end_of_main




// //====================================================================
// /// Create the files to illustrate the sparse node update operation
// //====================================================================
// void doc_sparse_node_update()
// {

//  // Set output directory
//  DocInfo doc_info;
//  doc_info.set_directory("RESLT");

//  // Setup mesh
 
//  // # of elements in x-direction
//  unsigned Nx=5;
 
//  // # of elements in y-direction
//  unsigned Ny=5;
 
//  // Domain length in x-direction
//  double Lx=5.0;
 
//  // Domain length in y-direction
//  double Ly=1.0;
 
//  // Build and assign mesh
//  SimpleSpineMesh<SpineElement<QTaylorHoodElement<2> > >* mesh_pt =
//   new SimpleSpineMesh<SpineElement<QTaylorHoodElement<2> > >(Nx,Ny,Lx,Ly);
 
//  // Update *all* nodal positions
//  mesh_pt->node_update();
 
//  unsigned count=0;
//  ofstream some_file;
//  char filename[100];
 
//  // Number of plot points
//  unsigned npts=5; 
 
//  // Output solution 
//  sprintf(filename,"%s/mesh_update%i.dat",doc_info.directory().c_str(),
//          count);
//  count++;
 
//  // Loop over spines
//  unsigned n_node = mesh_pt->nnode();
//  for (unsigned inode=0;inode<n_node;inode++)
//   {
//    SpineNode* node_pt=dynamic_cast<SpineNode*>(mesh_pt->node_pt(inode));
//    node_pt->node_update();
//    // Output solution 
//    some_file.open(filename);
//    sprintf(filename,"%s/mesh_update%i.dat",doc_info.directory().c_str(),
//            count);
//    count++;
//    mesh_pt->output(some_file,npts);
//    some_file.close();
//   }
// }

     













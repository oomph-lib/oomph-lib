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
// Driver code for a simple unstructured fluid problem using a mesh
// generated from an input file generated by the 3d mesh generator
// tetgen


//Generic routines
#include "generic.h"
#include "constitutive.h"
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;
using namespace oomph;


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Default Reynolds number
 double Re=100.0;

 /// Fluid pressure on inflow boundary
 double P_in=0.5;

 /// Applied traction on fluid at the inflow boundary
 void prescribed_inflow_traction(const double& t,
                                 const Vector<double>& x,
                                 const Vector<double>& n,
                                 Vector<double>& traction)
 {
  traction[0]=0.0;
  traction[1]=0.0;
  traction[2]=P_in;
 } 


 /// Fluid pressure on outflow boundary
 double P_out=-0.5; 

 /// Applied traction on fluid at the inflow boundary
 void prescribed_outflow_traction(const double& t,
                                  const Vector<double>& x,
                                  const Vector<double>& n,
                                  Vector<double>& traction)
 {
  traction[0]=0.0;
  traction[1]=0.0;
  traction[2]=-P_out;
 } 
 
} //end namespace






//======start_problem_class===========================================
/// Unstructured fluid problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredFluidProblem();

 /// Destructor (empty)
 ~UnstructuredFluidProblem(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Return total number of fluid inflow traction boundaries
 unsigned nfluid_inflow_traction_boundary()
  {
   return Inflow_boundary_id.size();
  }

 /// Return total number of fluid outflow traction boundaries
 unsigned nfluid_outflow_traction_boundary()
  {
   return Outflow_boundary_id.size();
  }

 /// Return total number of fluid outflow traction boundaries
 unsigned nfluid_traction_boundary()
  {
   return Inflow_boundary_id.size()+Outflow_boundary_id.size();
  }

 //private:

 /// Create fluid traction elements at inflow
 void create_fluid_traction_elements();

 /// Bulk fluid mesh
 TetgenMesh<ELEMENT>* Fluid_mesh_pt;

 /// Meshes of fluid traction elements that apply pressure at in/outflow
 Vector<Mesh*> Fluid_traction_mesh_pt;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Outflow_boundary_id;

};



//==========start_constructor=============================================
/// Constructor for unstructured 3D fluid problem
//========================================================================
template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
 
 //Create fluid bulk mesh, sub-dividing "corner" elements
 string node_file_name="fsi_bifurcation_fluid.1.node";
 string element_file_name="fsi_bifurcation_fluid.1.ele";
 string face_file_name="fsi_bifurcation_fluid.1.face";
 bool split_corner_elements=true;
 Fluid_mesh_pt =  new TetgenMesh<ELEMENT>(node_file_name,
                                          element_file_name,
                                          face_file_name,
                                          split_corner_elements);
 
 // Find elements next to boundaries
 //Fluid_mesh_pt->setup_boundary_element_info();

 // The following corresponds to the boundaries as specified by
 // facets in the tetgen input:

 // Fluid mesh has one inflow boundary: Boundary 0
 Inflow_boundary_id.resize(1);
 Inflow_boundary_id[0]=0;
 
 // Fluid mesh has two outflow boundaries: Boundaries 1 and 2
 Outflow_boundary_id.resize(2);
 Outflow_boundary_id[0]=1;
 Outflow_boundary_id[1]=2;
 
 // Apply BCs
 //----------
 
 // Map to indicate which boundary has been done
 std::map<unsigned,bool> done; 
  
 // Loop over inflow/outflow boundaries to impose parallel flow
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over in/outflow boundaries
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }

     // Number of nodes on that boundary
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get the node
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,inod);
       
       // Pin transverse velocities
       nod_pt->pin(0);
       nod_pt->pin(1);
      }
     
     // Done!
     done[b]=true;
    }

  } // done in and outflow
 
 
 
 // Loop over all fluid mesh boundaries and pin velocities
 // of nodes that haven't been dealt with yet
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<nbound;b++)
  {

   // Has the boundary been done yet?
   if (!done[b])
    {
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt= Fluid_mesh_pt->boundary_node_pt(b,inod);
       
       // Pin all velocities
       nod_pt->pin(0); 
       nod_pt->pin(1); 
       nod_pt->pin(2); 
      }
    }

  } // done no slip elsewhere 
 
 
 // Complete the build of the fluid elements so they are fully functional
 //----------------------------------------------------------------------
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {

   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;   

  } 
 
 
 // Create meshes of fluid traction elements at inflow/outflow
 //-----------------------------------------------------------
 
 // Create the meshes
 unsigned n=nfluid_traction_boundary();
 Fluid_traction_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Fluid_traction_mesh_pt[i]=new Mesh;
  } 
 
 // Populate them with elements
 create_fluid_traction_elements();
 
 
 // Combine the lot
 //----------------
 
 // Add sub meshes:

 // Fluid bulk mesh
 add_sub_mesh(Fluid_mesh_pt);
 
 // The fluid traction meshes
 n=nfluid_traction_boundary();
 for (unsigned i=0;i<n;i++)
  { 
   add_sub_mesh(Fluid_traction_mesh_pt[i]);
  }
 
 // Build global mesh
 build_global_mesh();

 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end constructor



//============start_of_fluid_traction_elements==============================
/// Create fluid traction elements 
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_fluid_traction_elements()
{

 // Counter for number of fluid traction meshes
 unsigned count=0;

 // Loop over inflow/outflow boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over boundaries with fluid traction elements
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }
     
     // How many bulk elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //What is the index of the face of the element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element 
       NavierStokesTractionElement<ELEMENT>* el_pt=
        new NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,
                                                       face_index);
       
       // Add it to the mesh
       Fluid_traction_mesh_pt[count]->add_element_pt(el_pt);
       
       // Set the pointer to the prescribed traction function
       if (in_out==0)
        {
         el_pt->traction_fct_pt() = 
          &Global_Parameters::prescribed_inflow_traction;
        }
       else
        {
         el_pt->traction_fct_pt() = 
          &Global_Parameters::prescribed_outflow_traction;
        }
      }
     // Bump up counter
     count++;
    }
  }
 
 } // end of create_traction_elements



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
  
 
 // Output fluid solution
 sprintf(filename,"%s/fluid_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();
  
}





//=============start_main=================================================
/// Demonstrate how to solve an unstructured 3D fluids problem
//========================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Label for output
 DocInfo doc_info;
 
 // Parameter study
 double Re_increment=100.0;
 unsigned nstep=4;
 if (CommandLineArgs::Argc==2)
  {
   std::cout << "Validation -- only doing two steps" << std::endl;
   nstep=2;
  }
 
 
 //Taylor--Hood
 {
  // Output directory
  doc_info.set_directory("RESLT_TH");
  
  //Set up the problem
  UnstructuredFluidProblem<TTaylorHoodElement<3> > problem;
  
  //Output initial guess
  problem.doc_solution(doc_info);
  doc_info.number()++;   
  
  // Parameter study: Crank up the pressure drop along the vessel
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Solve the problem
    problem.newton_solve();
    
    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    
    // Bump up Reynolds number (equivalent to increasing the imposed pressure
    // drop)
    Global_Parameters::Re+=Re_increment;   
   }
 }

 //Crouzeix Raviart
 {
  //Reset to default Reynolds number
  Global_Parameters::Re = 100.0;

  //Reset doc info number
  doc_info.number()=0;   

  // Output directory
  doc_info.set_directory("RESLT_CR");
  
  //Set up the problem
  UnstructuredFluidProblem<TCrouzeixRaviartElement<3> > problem;

  //Output initial guess
  problem.doc_solution(doc_info);
  doc_info.number()++;   
  
  // Parameter study: Crank up the pressure drop along the vessel
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Solve the problem
    problem.newton_solve();
    
    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    
    // Bump up Reynolds number (equivalent to increasing the imposed pressure
    // drop)
    Global_Parameters::Re+=Re_increment;   
   }
 }
 
} // end_of_main





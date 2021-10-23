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
//Generic routines
#include "generic.h"
#include "constitutive.h"
#include "navier_stokes.h"
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
 double Re=10.0;

 /// Fluid pressure on inflow boundary
 double P_in=0.5;


 /// Fluid pressure on outflow boundary
 double P_out=-0.5; 
 
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

 /// Update the problem specs before solve: empty
 void actions_before_newton_solve() {}

 /// Update the problem specs before solve: empty
 void actions_after_newton_solve() {}

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

private:
 
 /// Create Lagrange multiplier elements that enforce parallel outflow
 void create_parallel_outflow_lagrange_elements();

 /// Bulk fluid mesh
 TetgenMesh<ELEMENT>* Fluid_mesh_pt;

 /// Meshes of FaceElements imposing parallel outflow 
 /// and a pressure at the in/outflow
 Vector<Mesh*> Parallel_outflow_lagrange_multiplier_mesh_pt;

 /// IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

 /// IDs of fluid mesh boundaries along which inflow boundary conditions
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
 string node_file_name="fluid_iliac.1.node";
 string element_file_name="fluid_iliac.1.ele";
 string face_file_name="fluid_iliac.1.face";
 bool split_corner_elements=true;

 Fluid_mesh_pt =  new TetgenMesh<ELEMENT>(node_file_name,
                                            element_file_name,
                                            face_file_name,
                                          split_corner_elements);
 
 // Find elements next to boundaries
 Fluid_mesh_pt->setup_boundary_element_info();

 // The following corresponds to the boundaries as specified by
 // facets in the '.poly' input file:

 // Fluid mesh inflow boundaries
 Inflow_boundary_id.resize(22);
 for(unsigned i=0; i<22; i++)
  {
   Inflow_boundary_id[i]=215+i;
  }
 
 // Fluid mesh outflow boundaries
 Outflow_boundary_id.resize(11);
 for(unsigned i=0; i<11; i++)
  {

   Outflow_boundary_id[i]=237+i;

  } // done outflow boundaries


 // Create meshes of Lagrange multiplier elements at inflow/outflow
 //----------------------------------------------------------------

 // Create the meshes
 unsigned n=nfluid_traction_boundary();
 Parallel_outflow_lagrange_multiplier_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Parallel_outflow_lagrange_multiplier_mesh_pt[i]=new Mesh;
  } 
 
 // Populate them with elements
 create_parallel_outflow_lagrange_elements();

 
 // Combine the lot
 //----------------
 
 // Add sub meshes:

 // Fluid bulk mesh
 add_sub_mesh(Fluid_mesh_pt);
 
 // The fluid traction meshes
 n=nfluid_traction_boundary();
 for (unsigned i=0;i<n;i++)
  { 
   add_sub_mesh(Parallel_outflow_lagrange_multiplier_mesh_pt[i]);
  }
 
 // Build global mesh
 build_global_mesh();


 // Apply BCs
 //----------
 unsigned nbound=Fluid_mesh_pt->nboundary();

 // Vector indicating the boundaries where we have no slip
 std::vector<bool> pin_velocity(nbound, true);
 
 // Loop over inflow/outflow boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over in/outflow boundaries
   n=nfluid_inflow_traction_boundary();
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
     
     pin_velocity[b]=false;
    }

  } // done identification of boundaries where velocities are pinned


 // Loop over all boundaries to apply no slip where required
 for(unsigned b=0;b<nbound;b++)
  {
   if(pin_velocity[b])
    {
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,inod);

       // Pin all velocities
       nod_pt->pin(0); 
       nod_pt->pin(1); 
       nod_pt->pin(2); 
       
       // Find out if the node is also located on an in- or outflow
       // boundary
       bool is_in_or_outflow_node=false;
       for (unsigned in_out=0;in_out<2;in_out++)
        {
         // Loop over boundaries with Lagrange multiplier elements
         n=nfluid_inflow_traction_boundary();
         if (in_out==1) n=nfluid_outflow_traction_boundary();
         for (unsigned i=0;i<n;i++)
          {
           // Get boundary ID
           unsigned bb=0;
           if (in_out==0)
            {
             bb=Inflow_boundary_id[i];
            }
           else
            {
             bb=Outflow_boundary_id[i];
            }
           
           if(nod_pt->is_on_boundary(bb))
            {
             is_in_or_outflow_node=true;
            }
          }
        } // now we know if it's on the an in- or outflow boundary...


       // If its on an in- or outflow boundary pin the Lagrange multipliers
       if(is_in_or_outflow_node)
        {
         //Cast to a boundary node
         BoundaryNodeBase *bnod_pt = 
          dynamic_cast<BoundaryNodeBase*>
          (Fluid_mesh_pt->boundary_node_pt(b,inod) );
         
         // What's the index of the first Lagrange multiplier
         // in the node's values? 
         unsigned first_index=bnod_pt->index_of_first_value_assigned_by_face_element();
         
         // Pin the lagrange multiplier components 
         // in the out/in_flow boundaries
         for (unsigned l=0;l<2;l++)
          {
           nod_pt->pin(first_index+l);
          }
        }
      }
    }
  } // end of BC 

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
 
 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end constructor


//============start_of_lagrange_multiplier_elements======================
/// Create Lagrange multiplier elements that impose parallel outflow
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::
create_parallel_outflow_lagrange_elements()
{
 // Counter for number of Lagrange multiplier meshes
 unsigned count=0;

 // Loop over inflow/outflow boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over boundaries with Lagrange multiplier elements
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
       
       // Build the corresponding lagrange element
       ImposeParallelOutflowElement<ELEMENT>* el_pt = new 
        ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,face_index);
       
       // Add it to the mesh
       Parallel_outflow_lagrange_multiplier_mesh_pt[count]->
        add_element_pt(el_pt);
       
       // Set the pointer to the prescribed pressure
       if (in_out==0)
        {
         el_pt->pressure_pt()= &Global_Parameters::P_in;
        }
       else
        {
         el_pt->pressure_pt()= &Global_Parameters::P_out;
        }
      }
     // Bump up counter
     count++;
    }
  }

}  // done


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
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 //Set up the problem
 UnstructuredFluidProblem<TTaylorHoodElement<3> > problem;
 
 //Output initial guess
 problem.doc_solution(doc_info);
 doc_info.number()++;   

 // Parameter study
 double Re_increment=100.0;
 unsigned nstep=2;
 
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
 
} // end_of_main





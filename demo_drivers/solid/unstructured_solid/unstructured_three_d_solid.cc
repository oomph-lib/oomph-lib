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
// Driver code for a simple unstructured solid problem using a mesh
// generated from an input file generated by the 3d mesh generator
// tetgen

//Generic routines
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;

using namespace oomph;



//=========================================================================
/// Triangle-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticTetMesh : public virtual TetgenMesh<ELEMENT>, 
                       public virtual SolidMesh 
{
 
public:
 
 /// \short Constructor: 
 ElasticTetMesh(const std::string& node_file_name,
                const std::string& element_file_name,
                const std::string& poly_file_name,
                TimeStepper* time_stepper_pt=
                &Mesh::Default_TimeStepper) : 
  TetgenMesh<ELEMENT>(node_file_name, element_file_name,
                      poly_file_name, time_stepper_pt)
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
      
   // Identify special boundaries but remember that
   // outer boundary is already set to boundary 0; inner
   // (hole) boundary is already boundary 1.
   set_nboundary(4);

   unsigned n_node=this->nnode();
   for (unsigned j=0;j<n_node;j++)
    {
     Node* nod_pt=this->node_pt(j);

     // Boundary 2 is right boundary
     if (nod_pt->x(1)>2.99)
      {
       this->convert_to_boundary_node(nod_pt);
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(2,nod_pt);
      }

     // Boundary 3 is upper boundary
     if (nod_pt->x(2)>2.99)
      {
       this->convert_to_boundary_node(nod_pt);
       this->remove_boundary_node(0,nod_pt);
       this->add_boundary_node(3,nod_pt);
      }

    }
   TetgenMesh<ELEMENT>::setup_boundary_element_info();

  }

 /// Empty Destructor
 virtual ~ElasticTetMesh() { }


};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Poisson's ratio
 double Nu=0.3;

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=0.0;
  b[2]=-Gravity;
 }
 
 /// Uniform pressure
 double P = 0.0;

 /// \short Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
 void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 } // end traction


} //end namespace






//====================================================================
/// Unstructured solid problem
//====================================================================
template<class ELEMENT> 
class UnstructuredSolidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredSolidProblem();

 /// Destructor (empty)
 ~UnstructuredSolidProblem(){}

 /// Update the problem specs before solve: empty
 void actions_before_newton_solve() {}

 /// Update the problem specs before solve: empty
 void actions_after_newton_solve() {}

  /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Bulk mesh
 ElasticTetMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

};



//========================================================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT>
UnstructuredSolidProblem<ELEMENT>::
UnstructuredSolidProblem()
{ 

 //Create bulk mesh
 string node_file_name="cube_hole.1.node";
 string element_file_name="cube_hole.1.ele";
 string face_file_name="cube_hole.1.face";
 Solid_mesh_pt =  new ElasticTetMesh<ELEMENT>(node_file_name,
                                              element_file_name,
                                              face_file_name);
 
  // Traction elements are located on boundary 3:
  unsigned b=3;
  
  // Make traction mesh
  Traction_mesh_pt=new SolidMesh;
  
  // How many bulk elements are adjacent to boundary b?
  unsigned n_element = Solid_mesh_pt->nboundary_element(b);
  
  // Loop over the bulk elements adjacent to boundary b
  for(unsigned e=0;e<n_element;e++)
   {
    // Get pointer to the bulk element that is adjacent to boundary b
    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
     Solid_mesh_pt->boundary_element_pt(b,e));
    
    //Find the index of the face of element e along boundary b
    int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
    
    //Create solid traction element
    SolidTractionElement<ELEMENT> *el_pt = 
     new SolidTractionElement<ELEMENT>(bulk_elem_pt,face_index);   
    
    // Add to mesh
    Traction_mesh_pt->add_element_pt(el_pt);
    
    //Set the traction function
    el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
    
   }  
  
  
  // Add sub meshes
  add_sub_mesh(Solid_mesh_pt);
  add_sub_mesh(Traction_mesh_pt);
  
  // Build global mesh
  build_global_mesh();
  

  // Doc pinned solid nodes
  std::ofstream bc_file("pinned_nodes.dat");

  // Pin positions at right boundary (boundary 2)
  unsigned ibound=2;
  unsigned num_nod= Solid_mesh_pt->nboundary_node(ibound);  
  for (unsigned inod=0;inod<num_nod;inod++)
   {    
    // Get node
    SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
    
    // Pin all directions
    for (unsigned i=0;i<3;i++)
     {
      nod_pt->pin_position(i);
      
      // ...and doc it as pinned
      bc_file << nod_pt->x(i) << " ";
     }
    
    bc_file << std::endl;
   }
  
  // Complete the build of all elements so they are fully functional
  n_element = Solid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
   {
    //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   //Set the body force
   el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
   }
  
  
  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  
}



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output traction
 //----------------
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Traction_mesh_pt->output(some_file,npts);
 some_file.close();

}

 



//========================================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main()
{

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 //Set up the problem
 UnstructuredSolidProblem<TPVDElement<3,3> > problem;
 
 //Output initial configuration
 problem.doc_solution(doc_info);
 doc_info.number()++;

 // Parameter study
 Global_Physical_Variables::Gravity=12.0e-3;
 Global_Physical_Variables::P=0.0;
 double pressure_increment=-8.0e-3;
 
 unsigned nstep=2; // 10;
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution(doc_info);
   doc_info.number()++;
   
   // Bump up suction
   Global_Physical_Variables::P+=pressure_increment;
  }

}




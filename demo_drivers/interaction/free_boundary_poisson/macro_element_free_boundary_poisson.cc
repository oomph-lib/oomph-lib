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
// Driver for solution of "free boundary" 2D Poisson equation in 
// fish-shaped domain with adaptivity

 
// Generic oomph-lib headers
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The fish mesh 
#include "meshes/fish_mesh.h"

// Circle as generalised element:
#include "circle_as_generalised_element.h"

using namespace std;

using namespace oomph;


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=============start_of_namespace=====================================
/// Namespace for const source term in Poisson equation
//====================================================================
namespace ConstSourceForPoisson
{

/// Const source function
 void get_source(const Vector<double>& x, double& source)
 {
  source = -1.0;
 }
 
} // end of namespace







////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// MacroElementNodeUpdate-version of RefineableFishMesh
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//==========start_of_mesh=================================================
/// Refineable, fish-shaped mesh with MacroElement-based node update.
//========================================================================
template<class ELEMENT>
class MyMacroElementNodeUpdateRefineableFishMesh : 
 public virtual RefineableFishMesh<ELEMENT>,
 public virtual MacroElementNodeUpdateMesh
{

public: 

 ///  Constructor: Pass pointer to GeomObject that defines
 /// the fish's back and pointer to timestepper
 /// (defaults to (Steady) default timestepper defined in the Mesh
 /// base class).
 MyMacroElementNodeUpdateRefineableFishMesh(GeomObject* back_pt, 
                   TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper) :  
  FishMesh<ELEMENT>(back_pt,time_stepper_pt),
  RefineableFishMesh<ELEMENT>(time_stepper_pt)
  {
   // Set up all the information that's required for MacroElement-based
   // node update: Tell the elements that their geometry depends on the
   // fishback geometric object. 
   unsigned n_element = this->nelement();
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from FiniteElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(this->element_pt(i));

     // There's just one GeomObject
     Vector<GeomObject*> geom_object_pt(1);
     geom_object_pt[0] = back_pt;
     
     // Tell the element which geom objects its macro-element-based
     // node update depends on     
     el_pt->set_node_update_info(geom_object_pt);
    }

  } //end of constructor

 ///  Destructor: empty
 virtual ~MyMacroElementNodeUpdateRefineableFishMesh(){}

 ///  Resolve mesh update: Node update current nodal
 /// positions via sparse MacroElement-based update.
 //void node_update()
 // {
 //  MacroElementNodeUpdateMesh::node_update();
 // }

}; // end of mesh class



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==========start_of_problem_class====================================
/// Refineable "free-boundary" Poisson problem in deformable 
/// fish-shaped domain. Template parameter identifies the element.
//====================================================================
template<class ELEMENT>
class FreeBoundaryPoissonProblem : public Problem
{

public:

 ///   Constructor
 FreeBoundaryPoissonProblem();

 /// Destructor (empty)
 virtual ~FreeBoundaryPoissonProblem(){};

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Access function for the fish mesh
 MyMacroElementNodeUpdateRefineableFishMesh<ELEMENT>* fish_mesh_pt() 
  {
   return Fish_mesh_pt;
  }

 /// Doc the solution
 void doc_solution();

 ///  Before checking the new residuals in Newton's method
 /// we have to update nodal positions in response to possible 
 /// changes in the position of the domain boundary
 void actions_before_newton_convergence_check()
  {
   fish_mesh_pt()->node_update();
  }

private:

 /// Pointer to fish mesh
 MyMacroElementNodeUpdateRefineableFishMesh<ELEMENT>* Fish_mesh_pt;

 /// Pointer to single-element mesh that stores the GeneralisedElement
 /// that represents the fish's back
 Mesh* Fish_back_mesh_pt;

}; // end of problem class




//=========start_of_constructor===========================================
/// Constructor for adaptive free-boundary Poisson problem in 
/// deformable fish-shaped domain. 
//========================================================================
template<class ELEMENT>
FreeBoundaryPoissonProblem<ELEMENT>::FreeBoundaryPoissonProblem()
{ 

 // Set coordinates and radius for the circle that will become the fish back
 double x_c=0.5;
 double y_c=0.0;
 double r_back=1.0;

 // Build geometric object that will become the fish back
 ElasticallySupportedRingElement* fish_back_pt=
  new ElasticallySupportedRingElement(x_c,y_c,r_back);

 // Build fish mesh with geometric object that specifies the fish back 
 Fish_mesh_pt=new 
  MyMacroElementNodeUpdateRefineableFishMesh<ELEMENT>(fish_back_pt);

 // Add the fish mesh to the problem's collection of submeshes:
 add_sub_mesh(Fish_mesh_pt);

 // Create/set error estimator for the fish mesh
 fish_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Build mesh that will store only the geometric wall element
 Fish_back_mesh_pt=new Mesh;

 // So far, the mesh is completely empty. Let's add the 
 // GeneralisedElement that represents the shape
 // of the fish's back to it:
 Fish_back_mesh_pt->add_element_pt(fish_back_pt);

 // Add the fish back mesh to the problem's collection of submeshes:
 add_sub_mesh(Fish_back_mesh_pt);

 // Now build global mesh from the submeshes
 build_global_mesh();
 
 // Choose a control node: We'll use the
 // central node that is shared by all four elements in
 // the base mesh because it exists at all refinement levels.
 
 // How many nodes does element 0 have?
 unsigned nnod=fish_mesh_pt()->finite_element_pt(0)->nnode();

 // The central node is the last node in element 0:
 Node* control_node_pt=fish_mesh_pt()->finite_element_pt(0)->node_pt(nnod-1);
 
 // Use the solution (value 0) at the control node as the load
 // that acts on the ring. [Note: Node == Data by inheritance]
 dynamic_cast<ElasticallySupportedRingElement*>(Fish_mesh_pt->fish_back_pt())->
  set_load_pt(control_node_pt);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.  Set homogeneous boundary conditions everywhere
 unsigned num_bound = fish_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= fish_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     fish_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
     fish_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);

    }
  }
 
 /// Loop over elements and set pointers to source function
 unsigned n_element = fish_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(fish_mesh_pt()->element_pt(i));
   
   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::get_source;
  }

 // Do equation numbering
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//============start_of_doc================================================
/// Doc the solution in tecplot format.
//========================================================================
template<class ELEMENT>
void FreeBoundaryPoissonProblem<ELEMENT>::doc_solution()
{ 

 // Number of plot points in each coordinate direction.
 unsigned npts=5;

 // Output solution 
 ofstream some_file("RESLT/soln0.dat");
 fish_mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc

 



//==================start_of_main=========================================
/// Driver for "free-boundary" fish poisson solver with adaptation.
//========================================================================
int main()
{

 // Shorthand for element type
 typedef MacroElementNodeUpdateElement<RefineableQPoissonElement<2,3> > 
  ELEMENT;

 // Build problem
 FreeBoundaryPoissonProblem<ELEMENT> problem;

 // Do some uniform mesh refinement first
 problem.refine_uniformly();
 problem.refine_uniformly();
 
 // Solve/doc fully coupled problem, allowing for up to two spatial
 // adaptations. 
 unsigned max_solve=2; 
 problem.newton_solve(max_solve);
 problem.doc_solution();

} // end of main



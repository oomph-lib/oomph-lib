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
// (non-refineable) fish-shaped domain. This driver tests the non-refineable
// version of PoissonEquations<DIM>::get_dresidual_dnodal_coordinates(...)
// and so it evaluates the derivatives of the residual equation w.r.t.
// the nodal coordinates analytically.

 
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
// MacroElementNodeUpdate-version of FishMesh
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//==========start_of_mesh=================================================
/// Fish-shaped mesh with MacroElement-based node update.
//========================================================================
template<class ELEMENT>
class MyMacroElementNodeUpdateFishMesh : 
 public virtual FishMesh<ELEMENT>,
 public virtual MacroElementNodeUpdateMesh
{

public: 

 /// \short Constructor: Pass pointer to GeomObject that defines
 /// the fish's back and pointer to timestepper
 /// (defaults to (Steady) default timestepper defined in the Mesh
 /// base class).
 MyMacroElementNodeUpdateFishMesh(
  GeomObject* back_pt, 
  TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper) :  
  FishMesh<ELEMENT>(back_pt,time_stepper_pt)
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

 /// \short Destructor: empty
 virtual ~MyMacroElementNodeUpdateFishMesh(){}

}; // end of mesh class



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==========start_of_problem_class====================================
/// Non-refineable "free-boundary" Poisson problem in deformable 
/// fish-shaped domain. Template parameter identifies the element.
//====================================================================
template<class ELEMENT>
class FreeBoundaryPoissonProblem : public Problem
{

public:

 /// \short  Constructor
 FreeBoundaryPoissonProblem();

 /// Destructor (empty)
 virtual ~FreeBoundaryPoissonProblem(){};

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Access function for the fish mesh
 MyMacroElementNodeUpdateFishMesh<ELEMENT>* fish_mesh_pt() 
  {
   return Fish_mesh_pt;
  }

 /// Doc the solution
 void doc_solution();

 /// \short Before checking the new residuals in Newton's method
 /// we have to update nodal positions in response to possible 
 /// changes in the position of the domain boundary
 void actions_before_newton_convergence_check()
  {
   fish_mesh_pt()->node_update();
  }

private:

 /// Pointer to fish mesh
 MyMacroElementNodeUpdateFishMesh<ELEMENT>* Fish_mesh_pt;

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
  MyMacroElementNodeUpdateFishMesh<ELEMENT>(fish_back_pt);

 // Add the fish mesh to the problem's collection of submeshes:
 add_sub_mesh(Fish_mesh_pt);

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
 
 // Choose a control node:
 
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

   // Evaluate derivatives of residual equation w.r.t. nodal coordinates
   // analytically

   // It's broken but let's call it anyway to keep self-test alive
   bool i_know_what_i_am_doing=true;
   el_pt->evaluate_shape_derivs_by_chain_rule(i_know_what_i_am_doing);
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
 typedef MacroElementNodeUpdateElement<QPoissonElement<2,2> > 
  ELEMENT;

 // Build problem
 FreeBoundaryPoissonProblem<ELEMENT> problem;

 // Solve/doc fully coupled problem
 problem.newton_solve();
 problem.doc_solution();

} // end of main



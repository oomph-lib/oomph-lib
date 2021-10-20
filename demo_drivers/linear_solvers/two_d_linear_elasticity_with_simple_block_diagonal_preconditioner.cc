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
// Driver for a periodically loaded elastic body

// The oomphlib headers
#include "generic.h"
#include "linear_elasticity.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

// The simple block preconditioners
#include "simple_block_preconditioners.h"

using namespace std;

using namespace oomph;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


// Namespace extension
namespace oomph
{

//==start_of_mylinearelasticityelement===============================
/// Wrapper to make quadratic linear elasticity element block
/// preconditionable 
//===================================================================
template<unsigned DIM>
class MyLinearElasticityElement : public virtual QLinearElasticityElement<DIM,3>
{
 
public: 

 ///  The number of "DOF types" that degrees of freedom in this element
 /// are sub-divided into: The displacement components
 unsigned ndof_types() const
  {
   return DIM;
  }
 
/// Create a list of pairs for all unknowns in this element,
/// so the first entry in each pair contains the global equation
/// number of the unknown, while the second one contains the number
/// of the "DOF type" that this unknown is associated with.
/// (Function can obviously only be called if the equation numbering
/// scheme has been set up.)
/// 
/// The dof type enumeration (in 3D) is as follows:
/// S_x = 0
/// S_y = 1
/// S_z = 2
/// 
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // number of nodes
   unsigned n_node = this->nnode();
   
   // temporary pair (used to store dof lookup prior to being added to list)
   std::pair<unsigned,unsigned> dof_lookup;
   
   // loop over the nodes
   for (unsigned j=0;j<n_node;j++)
    {
     //loop over displacement components
     for (unsigned i=0;i<DIM;i++)
      {
       // determine local eqn number
       int local_eqn_number = this->nodal_local_eqn(j,i);
       
       // ignore pinned values - far away degrees of freedom resulting 
       // from hanging nodes can be ignored since these are be dealt
       // with by the element containing their master nodes
       if (local_eqn_number >= 0)
        {
         // store dof lookup in temporary pair: Global equation number
         // is the first entry in pair
         dof_lookup.first = this->eqn_number(local_eqn_number);
         
         // set dof numbers: Dof number is the second entry in pair
         dof_lookup.second = i;
         
         // add to list
         dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }

};


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
template<unsigned DIM>
class FaceGeometry<MyLinearElasticityElement<DIM> >
 : public virtual QElement<DIM-1,3> 
 {
 public:
  FaceGeometry() : QElement<DIM-1,3>() {}
 };

} // end namespace extension


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Amplitude of traction applied
 double Amplitude = 1.0;

 ///  Specify problem to be solved (boundary conditons for finite or
 /// infinite domain).
 bool Finite=false;

 /// Define Poisson coefficient Nu
 double Nu = 0.3;

 /// Length of domain in x direction
 double Lx = 1.0;

 /// Length of domain in y direction
 double Ly = 2.0;

 /// The elasticity tensor
 IsotropicElasticityTensor E(Nu);

 /// The exact solution for infinite depth case
 void exact_solution(const Vector<double> &x,
                     Vector<double> &u)
 {
  u[0] = -Amplitude*cos(2.0*MathematicalConstants::Pi*x[0]/Lx)*
	exp(2.0*MathematicalConstants::Pi*(x[1]-Ly))/
        (2.0/(1.0+Nu)*MathematicalConstants::Pi);
  u[1] = -Amplitude*sin(2.0*MathematicalConstants::Pi*x[0]/Lx)*
	exp(2.0*MathematicalConstants::Pi*(x[1]-Ly))/
        (2.0/(1.0+Nu)*MathematicalConstants::Pi);
 }

 /// The traction function
void periodic_traction(const double &time,
                       const Vector<double> &x,
                       const Vector<double> &n,
                       Vector<double> &result)
 {
  result[0] = -Amplitude*cos(2.0*MathematicalConstants::Pi*x[0]/Lx);
  result[1] = -Amplitude*sin(2.0*MathematicalConstants::Pi*x[0]/Lx);
 }
} // end_of_namespace


//===start_of_problem_class=============================================
/// Periodic loading problem
//======================================================================
template<class ELEMENT>
class PeriodicLoadProblem : public Problem
{
public:

 ///  Constructor: Pass number of elements in x and y directions 
 /// and lengths
 PeriodicLoadProblem(const unsigned &nx, const unsigned &ny, 
                     const double &lx, const double &ly);

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Allocate traction elements on the top surface
 void assign_traction_elements();
 
 /// Pointer to the bulk mesh
 Mesh* Bulk_mesh_pt;

 /// Pointer to the mesh of traction elements
 Mesh* Surface_mesh_pt;

 /// Solver
 IterativeLinearSolver* Solver_pt; 

 /// Preconditioner
 SimpleBlockDiagonalPreconditioner<CRDoubleMatrix>* Prec_pt;

}; // end_of_problem_class


//===start_of_constructor=============================================
/// Problem constructor: Pass number of elements in coordinate
/// directions and size of domain.
//====================================================================
template<class ELEMENT>
PeriodicLoadProblem<ELEMENT>::PeriodicLoadProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& ly)
{
 //Now create the mesh with periodic boundary conditions in x direction
 bool periodic_in_x=true;
 Bulk_mesh_pt = 
  new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,periodic_in_x);

 //Create the surface mesh of traction elements
 Surface_mesh_pt=new Mesh;
 assign_traction_elements();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin & set the ones that have Dirichlet 
 // conditions here
 unsigned ibound=0;
 unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get pointer to node
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

   // Pinned in x & y at the bottom and set value
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Check which boundary conditions to set and set them
   if (Global_Parameters::Finite)
     {
      // Set the displacements to zero
      nod_pt->set_value(0,0);
      nod_pt->set_value(1,0);
     }
   else
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);

      // Compute the value of the exact solution at the nodal point
      Vector<double> u(2);
      Global_Parameters::exact_solution(x,u);

      // Assign these values to the nodal values at this node
      nod_pt->set_value(0,u[0]);
      nod_pt->set_value(1,u[1]);
     }
  } // end_loop_over_boundary_nodes

 // Complete the problem setup to make the elements fully functional

 // Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Cast to a bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the elasticity tensor
   el_pt->elasticity_tensor_pt() = &Global_Parameters::E;
  }// end loop over elements

 // Loop over the traction elements
 unsigned n_traction =  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_traction;e++)
  {
   // Cast to a surface element
   LinearElasticityTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<LinearElasticityTractionElement<ELEMENT>* >
    (Surface_mesh_pt->element_pt(e));
   
   // Set the applied traction
   el_pt->traction_fct_pt() = &Global_Parameters::periodic_traction;
  }// end loop over traction elements

 // Add the submeshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Now build the global mesh
 build_global_mesh();

 // Assign equation numbers
 cout << assign_eqn_numbers() << " equations assigned" << std::endl; 

 // Create the solver.
 Solver_pt = new GMRES<CRDoubleMatrix>;

 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

 // Set linear solver
 linear_solver_pt() = Solver_pt;

 // Doc convergence
 Solver_pt->open_convergence_history_file_stream
  ("RESLT/iterative_solver_convergence.dat");

 // Create the preconditioner
 Prec_pt=new SimpleBlockDiagonalPreconditioner<CRDoubleMatrix>;

 // Block preconditioner can work with just the bulk mesh
 // since its elements contain all the degrees of freedom that
 // need to be classified. 
 Prec_pt->add_mesh(Bulk_mesh_pt);

 // Set the preconditioner
 Solver_pt->preconditioner_pt()=Prec_pt;


} // end of constructor


//===start_of_traction===============================================
/// Make traction elements along the top boundary of the bulk mesh
//===================================================================
template<class ELEMENT>
void PeriodicLoadProblem<ELEMENT>::assign_traction_elements()
{
// return;

 // How many bulk elements are next to boundary 2 (the top boundary)?
 unsigned bound=2;
 unsigned n_neigh = Bulk_mesh_pt->nboundary_element(bound); 
 
 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   FiniteElement *traction_element_pt 
    = new LinearElasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
 
   // Add to mesh
   Surface_mesh_pt->add_element_pt(traction_element_pt);
  }
 
} // end of assign_traction_elements

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PeriodicLoadProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 sprintf(filename,"%s/exact_soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,
                          Global_Parameters::exact_solution); 
 some_file.close();

 // Doc error
 double error=0.0;
 double norm=0.0;
 sprintf(filename,"%s/error.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             Global_Parameters::exact_solution, 
                             error,norm);
 some_file.close();

// Doc error norm:
 cout << "\nNorm of error    " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;


} // end_of_doc_solution   


//===start_of_main======================================================
/// Driver code for PeriodicLoad linearly elastic problem
//======================================================================
int main(int argc, char* argv[]) 
{
#ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Number of elements in x-direction
 unsigned nx=5;
 
 // Number of elements in y-direction (for (approximately) square elements)
 unsigned ny=unsigned(double(nx)*Global_Parameters::Ly/Global_Parameters::Lx);
 
 // Set up doc info
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 //Build the problem
 PeriodicLoadProblem<MyLinearElasticityElement<2> > 
  problem(nx,ny,Global_Parameters::Lx, Global_Parameters::Ly);
 
 // Solve
 problem.newton_solve();
 
 // Output the solution
 problem.doc_solution(doc_info);

#ifdef OOMPH_HAS_MPI
 // finalize MPI
 MPI_Helpers::finalize();
#endif

 return(EXIT_SUCCESS);  
} // end_of_main

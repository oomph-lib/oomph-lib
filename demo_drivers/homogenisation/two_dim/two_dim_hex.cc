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
// Driver code for a unit (hexagonal) cell homogenisation problem, solved in
// Willoughby et al (2012) IJSS, 49 pp 1421-1422

//Generic routines
#include "generic.h"


// The equations
#include "./Thomo_lin_elasticity_elements.h"
//#include "../../src/linear_elasticity/elasticity_tensor.cc"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
 

/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////


//A Comparison operator for the boundary nodes
class CompareNodeCoordinatesX
{
public:
/// The actual comparison operator
 int operator() (Node* const &node1_pt,
                 Node* const &node2_pt)
  {
   unsigned n_dim = node1_pt->ndim();
   if(n_dim != node2_pt->ndim())
    {
     throw OomphLibError("Can't compare two nodes of different dimension",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   //Make sure to handle the finite precision problems
   return node1_pt->x(0) < node2_pt->x(0);
  }
};

//A Comparison operator for the boundary nodes
class CompareNodeCoordinatesY
{
public:
/// The actual comparison operator
 int operator() (Node* const &node1_pt,
                 Node* const &node2_pt)
  {
   unsigned n_dim = node1_pt->ndim();
   if(n_dim != node2_pt->ndim())
    {
     throw OomphLibError("Can't compare two nodes of different dimension",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   //Make sure to handle the finite precision problems
   return node1_pt->x(1) < node2_pt->x(1);
  }
};




//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    
  //The lame coefficients for the bulk
  double P0 = 3.0;
  double M0 = 1.3;
  
  //The lame coefficients for the fibre
  double Pr = 23.5;
  double Mr = 30.0;

  //Pointer to the tensor for the bulk
  ElasticityTensor *E_bulk_pt;

  //Pointer to the tensor for the fibre
  ElasticityTensor *E_fibre_pt;

  /// Function that returns a pointer to the elasticity tensor
  /// associated with the bulk at position x
  void bulk_elasticity_tensor_pt(const Vector<double> &x, 
                                 ElasticityTensor* &E_pt)
  {
   E_pt = E_bulk_pt;
  }

  /// /Function that returns a pointer to the elasticity tensor
  /// associated with the fibre at position x
  void fibre_elasticity_tensor_pt(const Vector<double> &x, 
                                  ElasticityTensor* &E_pt)
  {
   E_pt = E_fibre_pt;
  }

  
 } // end_of_namespace


//==start_of_problem_class============================================
/// Problem class to simulate viscous inclusion propagating along 2D channel
//====================================================================
template<class ELEMENT>
class HomogenisationProblem : public Problem
{
 //Integers that are used to specify which sub-problem is being solved
 unsigned M, P;

 //Storage for the coefficients of the effective modulus
 Vector<Vector<DenseMatrix<double> > > C_eff;

 //Storage for the Internal Circular boundaries
 Vector<GeomObject*> Internal_circle_pt;

public:

 /// Constructor
 HomogenisationProblem();
 
 /// Destructor
 ~HomogenisationProblem()
  {
   // Kill data associated with outer boundary
   unsigned n= this->Outer_boundary_polyline_pt->npolyline();
   for (unsigned j=0;j<n;j++)
    {
     delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
   delete Outer_boundary_polyline_pt;

   //Kill data associated with inclusions
   unsigned n_inclusion = Inclusion_polygon_pt.size();
   for(unsigned iinclusion=0;iinclusion<n_inclusion;iinclusion++)
    {
     unsigned n=Inclusion_polygon_pt[iinclusion]->npolyline();
     for (unsigned j=0;j<n;j++)
      {
       delete Inclusion_polygon_pt[iinclusion]->polyline_pt(j);
      }
     delete Inclusion_polygon_pt[iinclusion];
    }
   // Delete solid mesh
   delete Bulk_mesh_pt;
  }

 /// Finish the steup of the problem
 void complete_problem_setup();
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve: 
 void actions_before_newton_solve() {}

 /// Calculate the values of the effective modulus by
 /// integrating over each element
 void calculate_coefficients()
  {
   //Loop over all elements
   unsigned n_element = this->Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     ELEMENT* el_pt = 
      dynamic_cast<ELEMENT*>(this->Bulk_mesh_pt->element_pt(e));
     
     //Add the contribution to the effective modulus
     el_pt->calculate_effective_modulus(C_eff[P][M]);
    }
  }
 
 /// Solve the sub-problems
 /// This will only solve the problem if the first_solve flag is true
 void sub_solve(const unsigned &n_dofs, DoubleVector &dx, DoubleVector &res,
                const bool &first_solve)
  {
   //Update anything that needs updating
   actions_before_newton_solve();
   //Do any updates that are required 
   actions_before_newton_step();
   
   //Now do the linear solve
   if(first_solve)
    {
     linear_solver_pt()->solve(this,dx);
    }
   else
    {
     //Get the new residuals
     get_residuals(res);
     //Now do the linear solve
     linear_solver_pt()->resolve(res,dx);
    }
   
   //Subtract the new values from the true dofs
   for(unsigned l=0;l<n_dofs;l++)
    { 
     // This is needed during parallel runs when dofs that are not
     // held on the current processor are nulled out. Can change
     // this once/if the Dof_pt vector is distributed too. 
     if (Dof_pt[l]!=0) *Dof_pt[l] -= dx[l];
    }

        
#ifdef OOMPH_HAS_MPI
     // Synchronise the solution on different processors
   this->synchronise_all_dofs();
#endif

     // Do any updates that are required 
     actions_after_newton_step();
     actions_before_newton_convergence_check();
     
     // Maximum residuals
     double maxres=0.0;
     //Calculate the new residuals
     get_residuals(dx);
     //Get the maximum residuals
     maxres = dx.max();

     oomph_info << "Final  residuals " << maxres << std::endl;
     
     //Now update anything that needs updating
     actions_after_newton_solve();
     
     //Calculate the current contribution to the coefficients
     calculate_coefficients();
  }


 /// Make our own solve function
 void solve() 
  {
   //Enable the resolve
   linear_solver_pt()->enable_resolve();

   //Find total number of dofs
   unsigned long n_dofs = ndof();
   
   //Set up the Vector to hold the solution
   //Vector<double> dx(n_dofs,0.0);
   DoubleVector dx;
   //Set up a vector to hold the residuals
   //Vector<double> res(n_dofs);
   DoubleVector res;

   //Only need to loop over the upper "half" of the matrix
   for(P=0;P<3;P++)
    {
     for(M=P;M<3;M++)
      {
       if(M==0)
        {
         //Do the solve once only
         sub_solve(n_dofs,dx,res,true);
        }
       //Otherwise it's resolves
       else
        {
         sub_solve(n_dofs,dx,res,false);
        }
      }
    }

   //Disable the resolve
   linear_solver_pt()->disable_resolve();
  }  


 //Document the current solution
 void doc_solution();
 
private:
 
 /// Pointer to Bulk_mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Vector storing pointer to the inclusion polygons
 Vector<TriangleMeshPolygon*> Inclusion_polygon_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
HomogenisationProblem<ELEMENT>::HomogenisationProblem() :  M(0), P(0)
{ 
 // Allocate the timestepper(a Steady default)
 this->add_time_stepper_pt(new Steady<0>);

 //Set the number of elements on each side
 unsigned n_element_on_side = 10;

 //The length of the side of a hexagon with unit area is
 const double l = sqrt(2.0)/(pow(3.0,0.75));

 //The height of the equilateral triangle adjacent to each side is the
 const double l_h = 1.0/(pow(12.0,0.25));

 //Set the volume fraction of the (single) fibre
 double volume_fraction = 0.15;


 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // six separate polylines
 //------------------------
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(6);
 
 //Each polyline has n_element_on_side + 1 vertices
 unsigned n_vertex= n_element_on_side + 1;
 Vector<Vector<double> > vertex_coord(n_vertex);
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i].resize(2);
  }

 //Set the start and end positions
 Vector<double> start(2), end(2);
 start[0] = -0.5*l;
 start[1] = -l_h;
 end[0] = -l;
 end[1] = 0.0;

 // First polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,5);

 //Set new start and end positions
 start = end;
 end[0] = -0.5*l;
 end[1] = l_h;

 // Second polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,4);

 
 //Set new start and end positions
 start = end;
 end[0] = 0.5*l;
 end[1] = l_h;

 // Second polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,3);
 

 //Set new start and end positions
 start = end;
 end[0] = l;
 end[1] = 0.0;

 // Second polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,2);


 //Set new start and end positions
 start = end;
 end[0] = 0.5*l;
 end[1] = -l_h;

 // Second polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[4] = new TriangleMeshPolyLine(vertex_coord,1);

 //Set new start and end positions
 start = end;
 end[0] = -0.5*l;
 end[1] = -l_h;

 // Second polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   for(unsigned j=0;j<2;j++)
    {
     vertex_coord[i][j] = 
      start[j] + i*(end[j] - start[j])/(double)(n_vertex-1);
    }
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[5] = new TriangleMeshPolyLine(vertex_coord,0);


 
 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);

 //Storage for the outer boundary
 TriangleMeshClosedCurve* outer_boundary_closed_curve_pt = 
  Outer_boundary_polyline_pt;

 // Now define initial shape of inclusion with a geometric object
 //---------------------------------------------------------------
 // Build an inclusion 
 Internal_circle_pt.resize(2);
 double rad = std::sqrt(volume_fraction/(4.0*atan(1.0)));
 double Radius[2] = {0.5*rad,rad};

 Vector<TriangleMeshClosedCurve*> curvilinear_inclusion_pt(2);

 for(unsigned h=0;h<2;h++)
  {
   Internal_circle_pt[h] = 
    new Circle(0.0,0.0,Radius[h],this->time_stepper_pt());
   
   // Build the two parts of the curvilinear boundary
   //Note that there could well be a memory leak here owing to some stupid 
   //choices
   Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);
   
   // First part of curvilinear boundary
   //-----------------------------------
   double zeta_start=0.0;
   double zeta_end=MathematicalConstants::Pi;
   unsigned nsegment=10;
   unsigned boundary_id=6+2*h;
   curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
    Internal_circle_pt[h],zeta_start,zeta_end, 
    nsegment,boundary_id);
   
   // Second part of curvilinear boundary
   //-------------------------------------
   zeta_start=MathematicalConstants::Pi;
   zeta_end=2.0*MathematicalConstants::Pi;
   nsegment=10;
   boundary_id=7 + 2*h;
   curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
    Internal_circle_pt[h],zeta_start,zeta_end, 
    nsegment,boundary_id);

   curvilinear_inclusion_pt[h]=
    new TriangleMeshClosedCurve(curvilinear_boundary_pt);

   //Delete the curvilines not possible because there are some memory
   //issues to consider
   //delete curvilinear_boundary_pt[1];
   //delete curvilinear_boundary_pt[0];
  }

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   outer_boundary_closed_curve_pt);

 // Define the holes on the boundary
 triangle_mesh_parameters.internal_closed_curve_pt() =
   curvilinear_inclusion_pt;


 double uniform_element_area = 0.001;

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 //Coordinates 
 Vector<double> region1(2,0.0);

 // Define the region (the internal fibre)
 triangle_mesh_parameters.add_region_coordinates(1, region1);
 
  Vector<double> region2(2);
  region2[0] = 0.0;
  region2[1] = 0.5*(Radius[1] + Radius[0]);

 //Set the second region (the annulus)
 triangle_mesh_parameters.add_region_coordinates(2, region2);

 //Prevent refinement on the boundaries so that we can apply periodic 
 //conditions
 triangle_mesh_parameters.disable_boundary_refinement();
 //Disable internal boundary refinement as well for precise control
 //triangle_mesh_parameters.disable_internal_boundary_refinement();

 // Create the mesh
 Bulk_mesh_pt =
   new TriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());

 //Output the regions
 unsigned n_region = this->Bulk_mesh_pt->nregion();
 for(unsigned i=0;i<n_region;i++)
  {
   unsigned n_element = Bulk_mesh_pt->nregion_element(i);

   std::stringstream output_string;
   output_string <<  "Region" << i << ".dat";
   std::ofstream output_file(output_string.str().c_str());
   
   for(unsigned e=0;e<n_element;e++)
    {
     Bulk_mesh_pt->region_element_pt(i,e)->output(output_file,5);
    }
   output_file.close();
  }

 // Output boundary and mesh initial mesh for information
 this->Bulk_mesh_pt->output_boundaries("boundaries.dat");
 this->Bulk_mesh_pt->output("mesh.dat");
 
 // Set boundary condition and complete the build of all elements
 this->complete_problem_setup();

 // Add meshes to the problem
 //----------------------------
 
 // Add Bulk_mesh_pt sub meshes
 this->add_sub_mesh(Bulk_mesh_pt);

 // Build global mesh
 this->build_global_mesh();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
 //Build storage for the effective modulus
 C_eff.resize(3);
 for(unsigned i=0;i<3;i++)
  {
   C_eff[i].resize(3);
   for(unsigned j=0;j<3;j++) 
    {
     C_eff[i][j].resize(3,3,0.0);
    }
  }

 //Solve the problem
 solve();

 //The final effective tensor is actually C_pmij because of the way I have
 //set up the storage. It shouldn't matter because of the symmetries of the
 //stress tensor, however.

 std::cout << "Output effective coefficients\n";

 std::cout << "C11 " << C_eff[0][0](0,0) << "\n";
 std::cout << "C12 " << C_eff[0][0](1,1) << "\n";
 std::cout << "C13 " << C_eff[0][0](2,2) << "\n";
 std::cout << "C16 " << C_eff[0][0](0,1) << "\n";
 std::cout << "C22 " << C_eff[1][1](1,1) << "\n";
 std::cout << "C23 " << C_eff[1][1](2,2) << "\n";
 std::cout << "C33 " << C_eff[2][2](2,2) << "\n";
 std::cout << "C36 " << C_eff[2][2](0,1) << "\n";
 std::cout << "C44 " << C_eff[1][2](1,2) << "\n";
 std::cout << "C45 " << C_eff[1][2](0,2) << "\n";
 std::cout << "C55 " << C_eff[0][2](0,2) << "\n";
 std::cout << "C66 " << C_eff[0][1](0,1) << "\n";

 //Output to a file
 std::ofstream output("C_eff.dat");

 output << "C11 " << C_eff[0][0](0,0) << "\n";
 output << "C12 " << C_eff[0][0](1,1) << "\n";
 output << "C13 " << C_eff[0][0](2,2) << "\n";
 output << "C16 " << C_eff[0][0](0,1) << "\n";
 output << "C22 " << C_eff[1][1](1,1) << "\n";
 output << "C23 " << C_eff[1][1](2,2) << "\n";
 output << "C33 " << C_eff[2][2](2,2) << "\n";
 output << "C36 " << C_eff[2][2](0,1) << "\n";
 output << "C44 " << C_eff[1][2](1,2) << "\n";
 output << "C45 " << C_eff[1][2](0,2) << "\n";
 output << "C55 " << C_eff[0][2](0,2) << "\n";
 output << "C66 " << C_eff[0][1](0,1) << std::endl;

 //Close the output file
 output.close();

} // end_of_constructor


//==start_of_complete_problem_setup=======================================
/// Set boundary conditions and complete the build of all elements
//========================================================================
template<class ELEMENT>
void HomogenisationProblem<ELEMENT>::complete_problem_setup()
  {      
   //Collect vectors of the nodes on the mesh boundaries
   Vector<Vector<Node*> > boundary_nodes_pt(6);
   
   for(unsigned b=0;b<6;b++)
    {
     unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
     boundary_nodes_pt[b].resize(n_node);
     for(unsigned n=0;n<n_node;n++)
      {
       boundary_nodes_pt[b][n] = this->Bulk_mesh_pt->boundary_node_pt(b,n);
      }
    }
   
   //Now let's sort each one
   std::sort(boundary_nodes_pt[0].begin(), boundary_nodes_pt[0].end(),
             CompareNodeCoordinatesX());
   std::sort(boundary_nodes_pt[1].begin(), boundary_nodes_pt[1].end(),
             CompareNodeCoordinatesY());
   std::sort(boundary_nodes_pt[2].begin(), boundary_nodes_pt[2].end(),
             CompareNodeCoordinatesY());
   std::sort(boundary_nodes_pt[3].begin(), boundary_nodes_pt[3].end(),
             CompareNodeCoordinatesX());
   std::sort(boundary_nodes_pt[4].begin(), boundary_nodes_pt[4].end(),
             CompareNodeCoordinatesY());
   std::sort(boundary_nodes_pt[5].begin(), boundary_nodes_pt[5].end(),
             CompareNodeCoordinatesY());

   //Now we can make it periodic
   double tol = 1.0e-7;
   unsigned n_node = boundary_nodes_pt[0].size();
   for(unsigned n=1;n<n_node-1;n++)
    {
     //If we have different x coordinates complain
     if(std::abs(boundary_nodes_pt[0][n]->x(0) 
                 - boundary_nodes_pt[3][n]->x(0)) > tol)
       {
        std::ostringstream error_stream;
        error_stream << 
         "Trying to make periodic nodes across the top/bottom boundary,\n"
                     << "but the nodes have x-coordinates "
                     << boundary_nodes_pt[0][n]->x(0) << " and "
                     << boundary_nodes_pt[3][n]->x(0);
                
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }

     //Make it periodic
     boundary_nodes_pt[0][n]->make_periodic(boundary_nodes_pt[3][n]);
    }
   
   n_node = boundary_nodes_pt[1].size();
   for(unsigned n=1;n<n_node-1;n++)
    {
     //Test?

     //Make it periodic
     boundary_nodes_pt[1][n]->make_periodic(boundary_nodes_pt[4][n]);
    }

   n_node = boundary_nodes_pt[2].size();
   for(unsigned n=1;n<n_node-1;n++)
    {
     //Test?

     //Make it periodic
     boundary_nodes_pt[2][n]->make_periodic(boundary_nodes_pt[5][n]);
    }


    //Vector of corner nodes
   Vector<Node*> corner_nodes_pt(3);
   corner_nodes_pt[0] = boundary_nodes_pt[0][0];
   corner_nodes_pt[1] = boundary_nodes_pt[1][n_node-1];
   corner_nodes_pt[2] = boundary_nodes_pt[3][0];
   
   //PUT IN A TEST TO CHECK THAT THE CORNERS ARE REALLY THE CORNERS!
   
   //Make these nodes periodic from the bottom node
   corner_nodes_pt[0]->make_periodic_nodes(corner_nodes_pt);
   
   //Vector of other corner nodes
   Vector<Node*> corner_nodes1_pt(3);
   corner_nodes1_pt[0] = boundary_nodes_pt[0][n_node-1];
   corner_nodes1_pt[1] = boundary_nodes_pt[4][0];
   corner_nodes1_pt[2] = boundary_nodes_pt[3][n_node-1];
   
   //PUT IN A TEST TO CHECK THAT THE CORNERS ARE REALLY THE CORNERS!
   
   //Make these nodes periodic from the bottom node
   corner_nodes1_pt[0]->make_periodic_nodes(corner_nodes1_pt);

   //Pin the corner to suppress rigid body motions
   corner_nodes_pt[0]->pin(0);
   corner_nodes_pt[0]->pin(1);
   corner_nodes_pt[0]->pin(2);
   
   // Complete the build of all elements so they are fully functional
   // Remember that adaptation for triangle meshes involves a complete
   // regneration of the mesh (rather than splitting as in tree-based
   // meshes where such parameters can be passed down from the father
   // element!)
   unsigned n_element = this->Bulk_mesh_pt->nregion_element(0);
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = 
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->region_element_pt(0,e));
     
     // Set pointer to the elasticity tensor function
     el_pt->elasticity_tensor_fct_pt() = 
      &Problem_Parameter::bulk_elasticity_tensor_pt;
     
     // Set the values of m and p
     el_pt->m_pt() = &(this->M);
     el_pt->p_pt() = &(this->P);
    }

   //For the elements within the inclusion (regions 1 and 2
   //set the viscosity ratio
   for(unsigned r=1;r<3;r++)
    {
     n_element = Bulk_mesh_pt->nregion_element(r);
     for(unsigned e=0;e<n_element;e++)
      {
       // Upcast from GeneralisedElement to the present element
       ELEMENT* el_pt = 
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->region_element_pt(r,e));

       el_pt->elasticity_tensor_fct_pt() = 
      &Problem_Parameter::fibre_elasticity_tensor_pt;
     
     // Set the values of m and p
     el_pt->m_pt() = &(this->M);
     el_pt->p_pt() = &(this->P);
      }
     
    }

  }



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void HomogenisationProblem<ELEMENT>::doc_solution()
{ 
 this->Bulk_mesh_pt->output("soln.dat",5);
}

//============================================================
/// Driver code for moving block problem
//============================================================
int main(int argc, char **argv)
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 {
  using namespace Problem_Parameter;
  
  Vector<double> lame_bulk(2);
  lame_bulk[0] = P0/M0; lame_bulk[1] = 1.0;
  E_bulk_pt = new IsotropicElasticityTensor(lame_bulk);
  Vector<double> lame_fibre(2);
  lame_fibre[0] = Pr/M0; lame_fibre[1] = Mr/M0;
  E_fibre_pt = new IsotropicElasticityTensor(lame_fibre);
  }


 
 // Create problem in initial configuration
 // which solves the required problem
 HomogenisationProblem<THomogenisedLinearElasticityElement<2,3> > problem;  

 problem.doc_solution();

} //End of main

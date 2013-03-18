//LIC// ====================================================================Ela
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


//Generic routines
#include "generic.h"


// The equations
#include "./homo_lin_elasticity_elements.h"
#include "./Thomo_lin_elasticity_elements.h"
#include "./homo_lin_elasticity_traction_elements.h"
#include "../../src/linear_elasticity/elasticity_tensor.cc"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
 

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//A Comparison operator for the boundary nodes
class CompareNodeCoordinatesX
{
public:
///The actual comparison operator
 int operator() (Node* const &node1_pt,
                 Node* const &node2_pt)
  {
   unsigned n_dim = node1_pt->ndim();
   if(n_dim != node2_pt->ndim())
    {
     throw OomphLibError("Can't compare two nodes of different dimension",
                         "CompareNodeCoordinates::operator()",
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
///The actual comparison operator
 int operator() (Node* const &node1_pt,
                 Node* const &node2_pt)
  {
   unsigned n_dim = node1_pt->ndim();
   if(n_dim != node2_pt->ndim())
    {
     throw OomphLibError("Can't compare two nodes of different dimension",
                         "CompareNodeCoordinates::operator()",
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
  double P0 = 3.0;
  double M0 = 1.3;
  
  double Pr = 23.5;
  double Mr = 30.0;

  ElasticityTensor *E_bulk_pt;

  ElasticityTensor *E_fibre_pt;

  void bulk_elasticity_tensor_pt(const Vector<double> &x, 
                                 ElasticityTensor* &E_pt)
  {
   E_pt = E_bulk_pt;
  }

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

 //Storage for the coefficients H
 Vector<Vector<Vector<DenseMatrix<double> > > > H;

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
   
   // Flush element of free surface elements
   delete_traction_elements();
   unsigned n_traction_mesh = Traction_mesh_pt.size();
   for(unsigned i=0;i<n_traction_mesh;i++) {delete Traction_mesh_pt[i];}

   // Delete error estimator
   //delete Bulk_mesh_pt->spatial_error_estimator_pt();

   // Delete fluid mesh
   delete Bulk_mesh_pt;
  }

 void complete_problem_setup();
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve: 
 void actions_before_newton_solve() {}

 /// \short Calculate the values of H by 
 /// performing the line integration around each hole
 void calculate_coefficients()
  {
   //Loop over each boundary
   unsigned nfaces = Traction_mesh_pt.size();
   
   //Loop over the faces
   for(unsigned b=0;b<nfaces;b++)
    {
     //Zero the matrix
     for(unsigned i=0;i<3;i++)
      {
       for(unsigned j=0;j<2;j++)
        {
         H[b][P][M](i,j) = 0.0;
        }
      }
     //Calculate the components
     unsigned nel = Traction_mesh_pt[b]->nelement();
     for(unsigned e=0;e<nel;e++)
      {
       HomogenisedLinearElasticityTractionElement<ELEMENT>* traction_el_pt = 
       dynamic_cast<HomogenisedLinearElasticityTractionElement<ELEMENT>*> 
        (Traction_mesh_pt[b]->element_pt(e));
       //Call the integration function
       traction_el_pt->calculate_H(H[b][P][M]);
      }
    }
  }

 ///Solve the sub-problems
 ///This will only solve the problem if the first_solve flag is true
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

     //calculate_coefficients();
  
     //Print the output
     /*ofstream mesh;
       char filename[100];
       sprintf(filename,"output%i%i.dat",P,M);
       mesh.open(filename);
       //Only output the bulk elements
       unsigned n_element = mesh_pt()->nelement();
       for(unsigned e=0;e<n_element;e++)
       {
       FiniteElement* elem_pt = mesh_pt()->finite_element_pt(e);
       //Only output if it's 2D (not a face element)
       if(elem_pt->dim() == 2) {elem_pt->output(mesh,5);}
       }
       //mesh_pt()->output(mesh,5);
       mesh.close();*/
  }


 /// \short Make our own solve function
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

 /// Actions before adapt: Wipe the mesh of free surface elements
 void actions_before_adapt()
  {
   delete_traction_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  
  }// end of actions_before_adapt

 
 /// Actions after adapt: Rebuild the mesh of free surface elements
 void actions_after_adapt()
  {
   // Create the elements that impose the displacement constraint 
   //create_traction_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
   
   // Setup the problem again -- remember that fluid mesh has been
   // completely rebuilt and its element's don't have any
   // pointers to Re etc. yet
   this->complete_problem_setup();
   
  }// end of actions_after_adapt

 void doc_solution();

private:
 

 /// \short Create free surface elements
 void create_traction_elements();

 /// \short Delete free surface elements 
 void delete_traction_elements()
  {
   //How many traction meshes
   unsigned n_traction_mesh = Traction_mesh_pt.size();
   for(unsigned m=0;m<n_traction_mesh;m++)
    {
     // How many surface elements are in the surface mesh
     unsigned n_element = Traction_mesh_pt[m]->nelement();
     
     // Loop over the surface elements
     for(unsigned e=0;e<n_element;e++)
      {
       // Kill surface element
       delete Traction_mesh_pt[m]->element_pt(e);
      }
     
     // Wipe the mesh
     Traction_mesh_pt[m]->flush_element_and_node_storage();
    }
  } // end of delete_traction_elements
 
 /// Pointers to mesh of traction elements
 Vector<Mesh*> Traction_mesh_pt;

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
 unsigned n_element_on_side = 40;
 //Each side has uniform length so calculate the spacing between vertices
 double vertex_spacing = 1.0 / (double)(n_element_on_side); 

 double volume_fraction = 0.05;

 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polylines
 //------------------------
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 
 //Each polyline has n_element_on_side + 1 vertices
 unsigned n_vertex= n_element_on_side + 1;
 Vector<Vector<double> > vertex_coord(n_vertex);
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i].resize(2);
  }

 // First polyline: Left-hand edge, set the nodal positions
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i][0]=0.0;
   vertex_coord[i][1] = i*vertex_spacing;
  }
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,3);
 
 // Second boundary polyline: Top edge
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i][0]= i*vertex_spacing;
   vertex_coord[i][1]=1.0;
  }

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord,2);

 // Third boundary polyline: Right-hand edge
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i][0]=1.0;
   vertex_coord[i][1]=1.0 - i*vertex_spacing;
  }

 // Build the 3rd boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,1);

 // Fourth boundary polyline: Bottom edge
 for(unsigned i=0;i<n_vertex;i++)
  {
   vertex_coord[i][0] = 1.0 - i*vertex_spacing;
   vertex_coord[i][1] = 0.0;
  }

 // Build the 4th boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord,0);
 
 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);

 //Storage for the outer boundary
 TriangleMeshClosedCurve* outer_boundary_closed_curve_pt = 
  Outer_boundary_polyline_pt;

 // Now define initial shape of inclusion with a geometric object
 //---------------------------------------------------------------
 // Build an inclusion 
 Internal_circle_pt.resize(1);
 double rad = std::sqrt(volume_fraction/(4.0*atan(1.0)));
 double Radius[2] = {rad,rad}; //{0.5*rad,rad};

 Vector<TriangleMeshClosedCurve*> curvilinear_inclusion_pt(1);

 for(unsigned h=0;h<1;h++)
  {
   Internal_circle_pt[h] = 
    new Circle(0.5,0.5,Radius[h],this->time_stepper_pt());
   
   // Build the two parts of the curvilinear boundary
   //Note that there could well be a memory leak here owing to some stupid 
   //choices
   Vector<TriangleMeshCurveSection*> curvilinear_boundary_pt(2);
   
   // First part of curvilinear boundary
   //-----------------------------------
   double zeta_start=0.0;
   double zeta_end=MathematicalConstants::Pi;
   unsigned nsegment=20;
   unsigned boundary_id=4+2*h;
   curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
    Internal_circle_pt[h],zeta_start,zeta_end, 
    nsegment,boundary_id);
   
   // Second part of curvilinear boundary
   //-------------------------------------
   zeta_start=MathematicalConstants::Pi;
   zeta_end=2.0*MathematicalConstants::Pi;
   nsegment=20;
   boundary_id=5 + 2*h;
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
 Vector<double> region1(2,0.5);

 // Define the region (the internal fibre)
 triangle_mesh_parameters.add_region_coordinates(1, region1);
 
 //Vector<double> region2(2);
 //region2[0] = 0.5;
 //region2[1] = 0.5 + 0.5*(Radius[1] + Radius[0]);

 //Set the second region (the annulus)
 //triangle_mesh_parameters.add_region_coordinates(2, region2);

 // Establish the use of regions when setting use attributes = true
 triangle_mesh_parameters.enable_use_attributes();

 //Prevent refinement on the boundaries
 triangle_mesh_parameters.disable_boundary_refinement();

 // Create the mesh
 Bulk_mesh_pt =
   new TriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());

 //Output the inner region
 for(unsigned i=0;i<3;i++)
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
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 //Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 //Bulk_mesh_pt->max_permitted_error()=0.005;
 //Bulk_mesh_pt->min_permitted_error()=0.001; 
 //Bulk_mesh_pt->max_element_size()=0.2;
 //Bulk_mesh_pt->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   //Bulk_mesh_pt->min_element_size()=0.01; 
  }

 // Output boundary and mesh initial mesh for information
 this->Bulk_mesh_pt->output_boundaries("boundaries.dat");
 this->Bulk_mesh_pt->output("mesh.dat");
 
 // Set boundary condition and complete the build of all elements
 this->complete_problem_setup();
 
 // Construct the mesh of free surface elements
 Traction_mesh_pt.resize(1);
 Traction_mesh_pt[0]=new Mesh;
 this->create_traction_elements();

 // Combine meshes
 //---------------
 
 // Add Bulk_mesh_pt sub meshes
 this->add_sub_mesh(Bulk_mesh_pt);

 // Add Traction sub meshes
 //this->add_sub_mesh(this->Traction_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 



 //Build storage for the answer
 H.resize(1);
 for(unsigned i=0;i<1;i++)
  {
   H[i].resize(3);
   for(unsigned p=0;p<3;p++) 
    {
     H[i][p].resize(3);
     for(unsigned m=0;m<3;m++) 
      {
       H[i][p][m].resize(3,2,0.0);
      }
    }
  }

 solve();



} // end_of_constructor


//============start_of_create_traction_elements===============
/// Create elements that impose the kinematic and dynamic bcs
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void HomogenisationProblem<ELEMENT>::create_traction_elements()
{ 
 //Loop over and assign from the fibre side

 for(unsigned b=4;b<6;b++)
  {
   // Note: region is important
   // How many bulk fluid elements are adjacent to boundary b in region 0?
   unsigned n_element = Bulk_mesh_pt->nboundary_element_in_region(b,2);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b in region 2
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_in_region_pt(b,2,e));
     
     //Find the index of the face of element e along boundary b in region 2
     int face_index = Bulk_mesh_pt->face_index_at_boundary_in_region(b,2,e);
     
     // Create new element
     HomogenisedLinearElasticityTractionElement<ELEMENT>* el_pt =
      new HomogenisedLinearElasticityTractionElement<ELEMENT>(
       bulk_elem_pt,face_index);   
     
     // Add it to the mesh
     Traction_mesh_pt[0]->add_element_pt(el_pt);
     
     //Add the appropriate boundary number
     el_pt->set_boundary_number_in_bulk_mesh(b);

    }
  }
}
// end of create_traction_elements




//==start_of_complete_problem_setup=======================================
/// Set boundary conditions and complete the build of all elements
//========================================================================
template<class ELEMENT>
void HomogenisationProblem<ELEMENT>::complete_problem_setup()
  {      
   //Collect vectors of the nodes on the mesh boundaries
   Vector<Vector<Node*> > boundary_nodes_pt(4);
   
   for(unsigned b=0;b<4;b++)
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
             CompareNodeCoordinatesX());
   std::sort(boundary_nodes_pt[3].begin(), boundary_nodes_pt[3].end(),
             CompareNodeCoordinatesY());


   //Now we can make it periodic
   double tol = 1.0e-7;
   unsigned n_node = boundary_nodes_pt[0].size();
   for(unsigned n=1;n<n_node-1;n++)
    {
     //If we have different x coordinates complain
     if(std::abs(boundary_nodes_pt[0][n]->x(0) 
                 - boundary_nodes_pt[2][n]->x(0)) > tol)
       {
        std::ostringstream error_stream;
        error_stream << 
         "Trying to make periodic nodes across the top/bottom boundary,\n"
                     << "but the nodes have x-coordinates "
                     << boundary_nodes_pt[0][n]->x(0) << " and "
                     << boundary_nodes_pt[2][n]->x(0);
                
        throw OomphLibError(error_stream.str(),
                            "HomogenisationProblem",
                            OOMPH_EXCEPTION_LOCATION);
       }

     //Make it periodic
     boundary_nodes_pt[0][n]->make_periodic(boundary_nodes_pt[2][n]);
    }
   
   n_node = boundary_nodes_pt[1].size();
   for(unsigned n=1;n<n_node-1;n++)
    {
     //If we have different y coordinates complain
     if(std::abs(boundary_nodes_pt[1][n]->x(1) - boundary_nodes_pt[3][n]->x(1)) > tol)
       {
        std::ostringstream error_stream;
        error_stream << 
         "Trying to make periodic nodes across the left/right boundary,\n"
                     << "but the nodes have y-coordinates "
                     << boundary_nodes_pt[1][n]->x(1) << " and "
                     << boundary_nodes_pt[3][n]->x(1);
                
        throw OomphLibError(error_stream.str(),
                            "HomogenisationProblem",
                            OOMPH_EXCEPTION_LOCATION);
       }

     //Make it periodic
     boundary_nodes_pt[1][n]->make_periodic(boundary_nodes_pt[3][n]);
    }


    //Vector of corner nodes
   Vector<Node*> corner_nodes_pt(4);
   corner_nodes_pt[0] = boundary_nodes_pt[1][0];
   corner_nodes_pt[1] = boundary_nodes_pt[1][n_node-1];
   corner_nodes_pt[2] = boundary_nodes_pt[3][0];
   corner_nodes_pt[3] = boundary_nodes_pt[3][n_node-1];
   
    //PUT IN A TEST TO CHECK THAT THE CORNERS ARE REALLY THE CORNERS!

    //Make these nodes periodic from the bottom node
    corner_nodes_pt[0]->make_periodic_nodes(corner_nodes_pt);

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
///Driver code for moving block problem
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
 HomogenisationProblem<THomogenisedLinearElasticityElement<2,3> > problem;  

 problem.doc_solution();
 

} //End of main

//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Demonstration code to solve a one-dimensional Poisson equation
//using the OOMPH-LIB base classes only

//Include std C++ IO library
#include <iostream>

//Include the OOMPH-LIB base classes
#include "generic.h"

using namespace std;

using namespace oomph;

//------------------------GEOMETRIC ELEMENT----------------------------------

//===========================================================================
/// A "Geometric" Element consists of Nodes only and contains all the 
/// information necessary to find the value of one, or more, global 
/// coordinates at the nodes and at any point within the element. It should 
/// also specify a spatial integration scheme that can be used to 
/// calculate the size of the element. 
///
/// A linear, one dimensional geometric element consists of two Nodes and 
/// uses linear interpolation between the Nodes. The integration scheme is 
/// a two-point Gaussian quadrature.
//===========================================================================
class TwoNodeGeometricElement : public FiniteElement
{
private:
 /// Integration scheme that will be used to integrate over the element.
 /// Simple Gaussian quadrature in one dimension, with two Gauss points.
 static Gauss<1,2> Default_spatial_integration_scheme;

public:

 /// Element constructor; resize Node_pt to the number of nodes (two)
 /// and set the integration scheme. 
 TwoNodeGeometricElement()
  {
   //Linear interpolation requires two Nodes per element.
   //In fact, calling this function merely provides storage for 
   //the pointers to the Nodes and initialises the pointers to NULL. 
   //The Nodes themselves are created during the mesh generation 
   //process by the functions FiniteElement::construct_node(...) 
   //which stores the pointers to the newly created Nodes 
   //in the element's own internal storage.
   this->set_n_node(2);

   //The element is one-dimensional 
   this->set_dimension(1);

   //Set the pointer to the spatial integration scheme
   set_integration_scheme(&Default_spatial_integration_scheme);
  }

 /// Define the shape functions for interpolation across the element.
 void shape(const Vector<double> &s, Shape &psi) const
  {
   //There are two shape functions (one per node)
   //In terms of the local coordinate, s[0]:
   //Node 0 is at s[0] = -1.0
   //Node 1 is at s[0] =  1.0
   
   //psi[0] takes the value one at node 0 and zero at node 1
   psi[0] = 0.5*(1.0 - s[0]);

   //psi[1] takes the value one at node 1 and zero at node 0
   psi[1] = 0.5*(1.0 + s[0]);
  }

 /// Define the derivatives of the shape functions w.r.t. the local coordinate
 void dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsids) const
  {
   //Call the shape functions
   shape(s,psi);

   //The derivative of psi[0] wrt s[0] is -0.5
   dpsids(0,0) = -0.5;

   //The derivative of psi[0] wrt s[0] is 0.5
   dpsids(1,0) = 0.5;
  }
 
 /// The element consists of two nodes along its one (and only) edge
 unsigned nnode_1d() const {return 2;}

};

//-----------------------SPECIFIC FINITE ELEMENT---------------------------

//=========================================================================
/// A specific Element needs to implement the functions that compute 
/// the elemental contributions to the global residual vector and the 
/// global Jacobian matrix formed when using Newton's method. 
/// By using calls to the interpolation scheme defined in an
/// an underlying "Geometric" Element, the equations can be written 
/// generically and Elements that solve the same equations, 
/// but using different orders of interpolation may be easily 
/// constructed, using a different "geometric" base class.
///
/// A linear, Poisson element solves a one-dimensional Poisson equation
/// using the integration and interpolation scheme defined in the
/// above "Geometric" element.
//=========================================================================
class TwoNodePoissonElement : public TwoNodeGeometricElement
 {

 public:

  /// Constructor: Initialise the sign of the source function to +1
  TwoNodePoissonElement()
   {
    // Initialise sign of source function
    Sign=1;
   }

  /// Define a specific source function. For greater generality, and
  /// inclusion in a library, this could be defined as a source function 
  /// pointer, that would then be set externally.
  double f(const double &x) const 
   {
    return double(Sign)*30.0*sin(sqrt(30.0)*x);
   }

  /// Access function to the sign in the source function
  int& sign() {return Sign;}

  /// Define an access function to the first data value stored 
  /// at each node. In a more general "Equation" element, such abstraction
  /// is essential, because different Elements will store the same variables
  /// in different locations.
  double u(const unsigned &n) {return node_pt(n)->value(0);}

  /// For the Poisson equation, only one value is stored at each node
  unsigned required_nvalue(const unsigned &n) const {return 1;}

  /// Calculate the elemental contributions to the global 
  /// residual vector for the weak form of the Poisson equation
  void get_residuals(Vector<double> &residuals)
   {
    //Find the number of degrees of freedom (unpinned values) in the element
    unsigned n_dof = ndof();

    //Initialise all the residuals to zero
    for(unsigned i=0;i<n_dof;i++) {residuals[i] = 0.0;}
    
    //Find the number of nodes in the element
    unsigned n_node = nnode();

    //Allocate memory for shape functions and their derivatives:
    // There's one shape function for each node:
    Shape psi(n_node);

    // Each of the n_node shape functions has one derivative with 
    // respect to the single local coordinate:
    DShape dpsidx(n_node,1);

    //Storage for the single local coordinate
    Vector<double> s(1);
    
    //Find the number of integration points in the underlying 
    //geometric element's integration scheme 
    unsigned n_intpt = integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Set the value of the local coordinate to be the integration 
      //scheme's knot point
      s[0] = integral_pt()->knot(ipt,0);

      //Find the weight of the integration scheme at this knot point
      double w = integral_pt()->weight(ipt);

      //Find the shape functions and their derivatives at the knot point. 
      //This function is implemented in FiniteElement.
      //It also returns the Jacobian of the mapping from local to 
      //global coordinates.
      double J = dshape_eulerian(s,psi,dpsidx);

      //Premultiply the weight and the Jacobian
      double W = w*J;
      
      //Allocate storage for the value of the field variable u,
      //its derivative and the global position at the knot point.
      //Initialise them all to zero.
      double interpolated_x=0.0, interpolated_u=0.0, interpolated_dudx=0.0;

      //Calculate the interpolated values by  looping over the shape 
      //functions and summing the appropriate contributions
      for(unsigned n=0;n<n_node;n++) 
       {
        interpolated_x += nodal_position(n,0)*psi[n];
        interpolated_u += u(n)*psi[n];
        interpolated_dudx += u(n)*dpsidx(n,0);
       }
      
      // Evaluate the source function
      double source=f(interpolated_x);

      //ASSEMBLE THE RESIDUALS
      
      //Loop over the test functions (same as the shape functions
      //since we're implementing an isoparametric element)
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the local equation number
        //The variable is the first (only) value stored at the nodes
        int local_eqn_number = nodal_local_eqn(l,0);

        //If the equation is not a Dirichlet boundary condition
        if(local_eqn_number >= 0)
         {
          //Add body force/source term here 
          residuals[local_eqn_number] += source*psi[l]*W;

          //Add the Poisson bit itself
          residuals[local_eqn_number] += interpolated_dudx*dpsidx(l,0)*W;
         }
       }

     } //End of loop over the integration points

   } //End of function


  /// Calculate the elemental contribution to the global residual
  /// vector and to the Jacobian matrix dR_{i}/du_{j} used in the Newton method
  void get_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian)
   {
    //First, calculate the residuals
    get_residuals(residuals);

    //Find the number of degrees of freedom (unpinned values) in the element
    unsigned n_dof = ndof();

    //Initialise all entries of the Jacobian matrix to zero
    for(unsigned i=0;i<n_dof;i++) 
     {
      for(unsigned j=0;j<n_dof;j++) {jacobian(i,j) = 0.0;}
     }
    
    //Find the number of nodes in the element
    unsigned n_node = nnode();
    //Allocate memory for shape functions and their derivatives
    Shape psi(n_node);
    DShape dpsidx(n_node,1);

    //Storage for the local coordinate
    Vector<double> s(1);
    
    //Find the number of integration points in the underlying
    //geometric element's integration scheme 
    unsigned n_intpt = integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Set the value of the local coordinate to be the integration 
      //scheme's knot point
      s[0] = integral_pt()->knot(ipt,0);

      //Find the weight of the integration scheme at this knot point
      double w = integral_pt()->weight(ipt);

      //Find the shape functions and their derivatives at the knot point. 
      //This function is implemented in FiniteElement.
      //It also returns the Jacobian of the mapping from local to 
      //global coordinates.
      double J = dshape_eulerian(s,psi,dpsidx);

      //Premultiply the weight and the Jacobian
      double W = w*J;
            
      //ASSEMBLE THE JACOBIAN TERMS
      
      //Loop over the test (shape) functions
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the local equation number
        //The variable is the first (only) value stored at the nodes
        int local_eqn_number = nodal_local_eqn(l,0);

        //If the equation is not a Dirichlet boundary condition
        if(local_eqn_number >= 0)
         {
          //Loop over the degrees of freedom
          for(unsigned l2=0;l2<n_node;l2++)
           {
            //Get the local degree of freedom number
            //The variable is the first (only) value stored at the nodes
            int local_dof_number = nodal_local_eqn(l2,0);

            //If the degree of freedom is not pinned
            if(local_dof_number >= 0)
             {
              //Add the contribution to the Jacobian
              jacobian(local_eqn_number,local_dof_number) += 
               dpsidx(l2,0)*dpsidx(l,0)*W;
             }
           }
         }
       }
     } //End of loop over the integration points

   } //End of function
  
  //Define an output function for the element 
  void output(ostream &output) 
   {
    //Read out the number of nodes in the element   
    unsigned n_node = nnode();

    //Loop over the nodes and print out the global coordinate 
    //and value of the field variable, u, at each node
    for(unsigned n=0;n<n_node;n++)
     {
      output << nodal_position(n,0) << " " << u(n) << std::endl;
     }
   } //End of function



  /// Self test function: The sign in the source function should
  /// only have the values +/- 1. Following the general oomph-lib convention,
  /// the self_test() returns 0 for success, and 1 for failure:
  unsigned self_test()
   {
    // Initialise success flag
    unsigned success=0;

    // Run the generic FiniteElement self test
    success=FiniteElement::self_test();

    // Do additional test for this function
    if ((Sign!=1)&&(Sign!=-1))
     {
      cout << "Sign of source function should be +/- 1," << std::endl;
      cout << "but it is:" << Sign << std::endl;
      success=1;
     }

    // Return success flag
    return success;

   } // End of self test


 private:

  /// The sign of the source function
  int Sign;

}; //End of the class

/// Define the static spatial integration scheme
Gauss<1,2> TwoNodeGeometricElement::Default_spatial_integration_scheme;

//----------------------ONE DIMENSIONAL (LINE) MESH---------------------------

//============================================================================
/// A simple one dimensional mesh: uniformly spaced nodes in the domain x=0,1
//============================================================================
template<class ELEMENT>
class OneDimMesh : public Mesh
{

public:

 /// Mesh Constructor. The argument is the desired number of elements
 OneDimMesh(const unsigned &n_element)
 {
  //Resize the vector of pointers to elements: there are n_element elements
  Element_pt.resize(n_element); 

  //Construct the first element (Note the use of the template parameter)
  Element_pt[0] = new ELEMENT;

  //Construct the first node and add it to the Mesh::Node_pt vector
  //Note: The FiniteElement::construct_boundary_node(j) function
  //builds the element's j-th local node, and provides the functionality
  //that allows it to be located on a Mesh boundary -- essentially this
  //involves allocating additional storage to the Node.
  //The function obtains the Node's
  //characteristics (e.g. its spatial dimension, the number of
  //values to be stored, etc) from various virtual FiniteElement
  //member functions, such as FiniteElement::required_nvalue(). 
  //FiniteElement::construct_boundary_node(...) also
  //stores a pointer to the newly created Node in the element's own
  //Node_pt vector.
  //Finally, the function returns the pointer to the
  //newly created Node, so that it can be stored in the Mesh's Node_pt
  //vector, as done here:
  Node_pt.push_back(finite_element_pt(0)->construct_boundary_node(0));

  //Find the number of nodes per element (N.B. all elements are identical
  //so we can determine this value once and for all). 
  unsigned n_node = finite_element_pt(0)->nnode();

  //Loop over the remaning nodes of the first element
  for(unsigned n=1;n<n_node;n++)
   {
    //Construct the next node and add it to the Mesh::Node_pt vector
    //Note that these interior nodes need not (and should not)
    //be boundary nodes, so they are created using the construct_node
    //function, which has the same interface as
    //construct_boundary_node()
    Node_pt.push_back(finite_element_pt(0)->construct_node(n));
   }

  //Loop over the remaining elements apart from the last
  for(unsigned e=1;e<(n_element-1);e++)
   {
    //Construct the e-th element
    Element_pt[e] = new ELEMENT;

    //The first local node of the e-th element is the last local node
    //of the (e-1)-th element. We MUST NOT construct the node twice.
    //Instead, we set the pointer in the e-th element to point to the
    //previously created node in the (e-1)-th element.
    finite_element_pt(e)->node_pt(0) = 
      finite_element_pt(e-1)->node_pt(n_node-1);

    //Loop over the remaining nodes of the e-th element
    for(unsigned n=1;n<n_node;n++)
     {
      //Construct the next node and add it to the Mesh::Node_pt vector
      //Note that these interior nodes need not (and should not)
      //be boundary nodes, so they are created using the construct_node
      //function, which has the same interface as
      //construct_boundary_node()
      Node_pt.push_back(finite_element_pt(e)->construct_node(n));
     }
   } //End of loop over elements
  
  
  //Construct the final element
  Element_pt[n_element-1] = new ELEMENT;
  
  //The first local node of the final element is the last local node
  //of the penultimate element. We MUST NOT construct the node twice.
  //Instead, we set the pointer in the final element to point to the
  //previously created node in the penultimate element.
  finite_element_pt(n_element-1)->node_pt(0) = 
   finite_element_pt(n_element-2)->node_pt(n_node-1);

  //Loop over the remaining central nodes of the final element
  for(unsigned n=1;n<(n_node-1);n++)
   {
    //Construct the next node and add it to the Mesh::Node_pt vector
    //Note that these interior nodes need not (and should not)
    //be boundary nodes, so they are created using the construct_node
    //function()
    Node_pt.push_back(finite_element_pt(n_element-1)->construct_node(n));
   }

  //Construct the final node and add it to the Mesh::Node_pt vector.
  //This node will be located on a boundary, and hence we use
  //the construct_boundary_node function.
  Node_pt.push_back(finite_element_pt(n_element-1)
                    ->construct_boundary_node(n_node-1));

  //We've now created all the nodes -- let's set their positions:

  //Find the total number of nodes
  unsigned n_global_node = nnode();

  //Loop over all nodes
  for(unsigned n=0;n<n_global_node;n++)
   {
    //Set the position of the node (equally spaced through the unit interval)
    Node_pt[n]->x(0) = double(n)/double(n_global_node-1);
   }
 
  //Set the boundary data:

  //There are two boundaries in this mesh
  set_nboundary(2);

  //Boundary 0 contains the first node in the mesh:
  add_boundary_node(0,Node_pt[0]);

  //Boundary 1 contains the final node in the mesh:
  add_boundary_node(1,Node_pt[n_global_node-1]); 
  
 } // End of constructor

}; // End of OneDimMesh class.

//-----------------------PROBLEM CLASS-------------------------------------

//==========================================================================
/// Define the Problem which creates the mesh, applies the
/// boundary conditions, and assigns equation numbers.
//==========================================================================
class DemoPoissonProblem : public Problem 
 {
  public:

  /// Problem constructor: Pass the sign of the source function (default
  /// is +1)
  DemoPoissonProblem(const int& sign=1) : Sign(sign)
   {
    //Create a OneDimMesh Mesh object and set it to be the problem's mesh.
    //The element type, TwoNodePoissonElement, is passed  as a template 
    //parameter to the mesh. The argument to the constructor indicates
    //the number of elements in the mesh.
    Problem::mesh_pt() = new OneDimMesh<TwoNodePoissonElement>(10);

    //Pin the unknowns at the ends of the 1D domain:
    //The 1D mesh has 2 boundaries, each of which contains a single node;
    //the nodes on the boundary are available from Mesh::boundary_node_pt(...)

    //Pin the single nodal value at the single node on mesh 
    //boundary 0 (= the left domain boundary at x=0)
    mesh_pt()->boundary_node_pt(0,0)->pin(0);

    //Pin the single nodal value at the single node on mesh 
    //boundary 1 (= the right domain boundary at x=1)
    mesh_pt()->boundary_node_pt(1,0)->pin(0);

    // All values are initialised to zero. This is consistent with the
    // boundary condition at x=0 and no further action is required
    // at that node.

    // Apply the boundary condition at x=1: u(x=1)=-/+1
    mesh_pt()->boundary_node_pt(1,0)->set_value(0,-double(Sign));

    
    // Finish problem setup: Set the sign for the source function
    // in all elements
    
    //Find number of elements in mesh
    unsigned n_element = mesh_pt()->nelement();
    
    // Loop over the elements 
    for(unsigned i=0;i<n_element;i++)
     {
      // The sign() member function is defined in the 
      // TwoNodePoissonElement class not the base GeneralisedElement class.
      // In order to use it, we must
      // upcast from GeneralisedElement to the specific element type,
      // which is achieved by a C++ dynamic_cast.
      TwoNodePoissonElement *specific_element_pt 
       = dynamic_cast<TwoNodePoissonElement*>(mesh_pt()->element_pt(i));
      
      // Set the sign of the source function
      specific_element_pt->sign() = Sign;
     }

    //Assign the global and local equations numbers for the problem
    cout << "Number of equations is " << assign_eqn_numbers() << std::endl;
   }

  
  /// Check that everything has been set up properly
  void actions_before_newton_solve()
   {
    if (0==self_test())
     {
      cout << "Problem has been set up correctly and can be solved." 
           << std::endl;
     }
    else
     {
      throw 
       OomphLibError("Trouble! Check error messages and fix the problems.\n",
                     "DemoPoissonProblem::actions_before_newton_solve()",
                     OOMPH_EXCEPTION_LOCATION);
     }
   }

  
  /// Print out the result after the solve
  void actions_after_newton_solve() 
   {
    ofstream file("result.dat");
    mesh_pt()->output(file);
   }

 private:

  /// The sign of the source function
  int Sign;


}; //End of problem definition

//----------------------------MAIN FUNCTION-------------------------------

int main()
 {
  //Build the problem 
  DemoPoissonProblem problem;
  
  //Solve the problem, using Newton's method
  problem.newton_solve();

 }

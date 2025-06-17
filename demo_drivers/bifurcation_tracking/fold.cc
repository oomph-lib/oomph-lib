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
//Driver code for a simple fold bifurcation detection and tracking problem

//Standard system includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio>

//Include files from the finite-element library
#include "generic.h"
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

//======================================================================== 
/// Global variables that represent physical properties
//======================================================================== 
namespace Global_Physical_Variables
{
 double *Lambda_pt;
 double *Mu_pt;
}

template<unsigned NNODE_1D>
class GelfandBratuElement : public QElement<1,NNODE_1D> 
{
 double *Lambda_pt;

 double *Mu_pt;
 
 public:

 GelfandBratuElement() {}

 //Interface to the parameter
 const double &lambda() const {return *Lambda_pt;}

 const double &mu() const {return *Mu_pt;}

 //Set the pointer
 double* &lambda_pt() {return Lambda_pt;}

 double* &mu_pt() {return Mu_pt;}

 /// For the Equation, only one value is stored at each node
 unsigned required_nvalue(const unsigned &n) const {return 1;}

  /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }
 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,jacobian,1);
  }
 

 
 /// Calculate the elemental contributions to the global 
 /// residual vector for the weak form of the Gelfand-Bratu equation
 void fill_in_generic_residual_contribution(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian,
                                            unsigned flag)
  {
    //Find the number of nodes in the element
    unsigned n_node = this->nnode();

    //Allocate memory for shape functions and their derivatives:
    // There's one shape function for each node:
    Shape psi(n_node);
    // Each of the n_node shape functions has one derivative with 
    // respect to the single local coordinate:
    DShape dpsidx(n_node,1);

    //Get the value of the parameter
    double lam = lambda();
    double m = mu();

    //Find the number of integration points in the underlying 
    //geometric element's integration scheme 
    unsigned n_intpt = this->integral_pt()->nweight();
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Set the value of the local coordinate to be the integration 
      //scheme's knot point
      //Find the weight of the integration scheme at this knot point
      double w = this->integral_pt()->weight(ipt);
      //Find the shape functions and their derivatives at the knot point. 
      //This function is implemented in FiniteElement.
      //It also returns the Jacobian of the mapping from local to 
      //global coordinates.
      double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
      //Premultiply the weight and the Jacobian
      double W = w*J;
      
      //Allocate storage for the value of the field variable u,
      //its derivative and the global position at the knot point.
      //Initialise them all to zero.
      double interpolated_u=0.0, interpolated_dudx=0.0;
      //Calculate the interpolated values by  looping over the shape 
      //functions and summing the appropriate contributions
      for(unsigned n=0;n<n_node;n++) 
       {
        double value = this->nodal_value(n,0);
        interpolated_u += value*psi[n];
        interpolated_dudx += value*dpsidx(n,0);
       }
   
      //Get the source function
      double source = (lam*lam + m*m)*exp(interpolated_u);
   
      //ASSEMBLE THE RESIDUALS
      
      //Loop over the test functions (same as the shape functions
      //since we're implementing an isoparametric element)
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the local equation number
        //The variable is the first (only) value stored at the nodes
        int local_eqn_number = this->nodal_local_eqn(l,0);
        //If the equation is not a Dirichlet boundary condition
        if(local_eqn_number >= 0)
         {
          //Add body force/source term here 
          residuals[local_eqn_number] -= source*psi[l]*W;
          //Add the Poisson bit itself
          residuals[local_eqn_number] += interpolated_dudx*dpsidx(l,0)*W;

          //If we are doing the jacobian terms
          if(flag)
           {
            for(unsigned l2=0;l2<n_node;l2++)
             {
              int local_unknown = this->nodal_local_eqn(l2,0);
              if(local_unknown >= 0)
               {
                jacobian(local_eqn_number,local_unknown) -=
                 source*psi[l2]*psi[l]*W;
                jacobian(local_eqn_number,local_unknown) +=
                 dpsidx(l2,0)*dpsidx(l,0)*W;
               }
             }
           }
         }
       }
     } //End of loop over the integration points
   } //End of function

 /// Add the element's contribution to the derivatives of 
 /// its residual vector with respect to a parameter (wrapper)
 void fill_in_contribution_to_dresiduals_dparameter(
  double* const &parameter_pt, Vector<double> &dres_dparam)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_dresidual_contribution(
    parameter_pt,dres_dparam,GeneralisedElement::Dummy_matrix,0);
  }
 
 /// Add the element's contribution to the derivaives of its 
 /// residual vector and element Jacobian matrix with respect to a parameter
 /// (wrapper)
 void fill_in_contribution_to_djacobian_dparameter(
  double* const &parameter_pt, Vector<double> &dres_dparam,
  DenseMatrix<double> &djac_dparam)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_dresidual_contribution(parameter_pt,
                                          dres_dparam,djac_dparam,1);
  }

 
 /// Calculate the elemental contributions to the derivatives of
 /// the global residual vector and jacobian for the weak form of the 
 /// Gelfand-Bratu equation with respect to the passed parameter.
 void fill_in_generic_dresidual_contribution(
  double* const &parameter_pt,
  Vector<double> &dres_dparam,
  DenseMatrix<double> &djac_dparam,
  unsigned flag)
  {
   //There are only two parameters, if it's not either of them, then
   //die
   if((parameter_pt!=Lambda_pt) && (parameter_pt!=Mu_pt))
    {
     std::ostringstream error_stream;
     error_stream << 
      "Cannot compute analytic jacobian for parameter addressed by " 
                  << parameter_pt << "\n";
     error_stream << "Can only compute derivatives wrt lambda ("
                  << Lambda_pt << ") or mu (" << Mu_pt << ")\n";
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }

   //Find the number of nodes in the element
   unsigned n_node = this->nnode();
   
   //Allocate memory for shape functions and their derivatives:
   // There's one shape function for each node:
   Shape psi(n_node);
   // Each of the n_node shape functions has one derivative with 
   // respect to the single local coordinate:
   DShape dpsidx(n_node,1);
   
   //Get the value of the parameter
   double lam = lambda();
   double m = mu();

    //Find the number of integration points in the underlying 
    //geometric element's integration scheme 
    unsigned n_intpt = this->integral_pt()->nweight();
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Set the value of the local coordinate to be the integration 
      //scheme's knot point
      //Find the weight of the integration scheme at this knot point
      double w = this->integral_pt()->weight(ipt);
      //Find the shape functions and their derivatives at the knot point. 
      //This function is implemented in FiniteElement.
      //It also returns the Jacobian of the mapping from local to 
      //global coordinates.
      double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
      //Premultiply the weight and the Jacobian
      double W = w*J;
      
      //Allocate storage for the value of the field variable u,
      //its derivative and the global position at the knot point.
      //Initialise them all to zero.
      double interpolated_u=0.0, interpolated_dudx=0.0;
      //Calculate the interpolated values by  looping over the shape 
      //functions and summing the appropriate contributions
      for(unsigned n=0;n<n_node;n++) 
       {
        double value = this->nodal_value(n,0);
        interpolated_u += value*psi[n];
        interpolated_dudx += value*dpsidx(n,0);
       }
   
      //Get the derivative of the source function with respect to the 
      //parameter
      double source = exp(interpolated_u);

      //If the parameter is lambda
      if(parameter_pt==Lambda_pt) {source *= 2.0*lam;} 
      //Only other possibility is mu 
      else {source *= 2.0*m;}
   
      //ASSEMBLE THE DERIVATIVE OF THE RESIDUALS
      
      //Loop over the test functions (same as the shape functions
      //since we're implementing an isoparametric element)
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the local equation number
        //The variable is the first (only) value stored at the nodes
        int local_eqn_number = this->nodal_local_eqn(l,0);
        //If the equation is not a Dirichlet boundary condition
        if(local_eqn_number >= 0)
         {
          //Add body force/source term here 
          dres_dparam[local_eqn_number] -= source*psi[l]*W;

          //If we are doing the jacobian terms
          if(flag)
           {
            for(unsigned l2=0;l2<n_node;l2++)
             {
              int local_unknown = this->nodal_local_eqn(l2,0);
              if(local_unknown >= 0)
               {
                djac_dparam(local_eqn_number,local_unknown) -=
                 source*psi[l2]*psi[l]*W;
               }
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
    unsigned n_node = this->nnode();
    //Loop over the nodes and print out the global coordinate 
    //and value of the field variable, u, at each node
    for(unsigned n=0;n<n_node;n++)
     {
      output << this->nodal_position(n,0) << " " 
             << this->nodal_value(n,0) << std::endl;
     }
   } //End of function

}; //End of the class


//======================================================================
//Problem class to solve the deformation of an elastic tube
//=====================================================================
template<class ELEMENT>
class BratuProblem : public Problem
{
private:

 Node* Trace_node_pt;

public:

 //Constructor
 BratuProblem(const unsigned &Nx);

 //Overload Access function for the mesh
 OneDMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<OneDMesh<ELEMENT>*>(Problem::mesh_pt());}

 //Update functions are both empty
 void actions_after_newton_solve() {}
 void actions_before_newton_solve() {}

 void solve();
 
};

//Constructor
template<class ELEMENT>
BratuProblem<ELEMENT>::BratuProblem
(const unsigned &Nx)
{
 //Now create the mesh
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(Nx,1.0); 

 //Pin the ends
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 mesh_pt()->boundary_node_pt(1,0)->pin(0);

 //Assign the physical parameters
 using namespace Global_Physical_Variables;
 //Set the pointer to the external pressure 
 Lambda_pt = new double(0.0);
 Mu_pt = new double(0.1);
 
 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();
 //Loop over the elements 
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to a shell element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the load function
   el_pt->lambda_pt() = Lambda_pt;
   el_pt->mu_pt() = Mu_pt;
  }

 Trace_node_pt = mesh_pt()->finite_element_pt(n_element/2)->node_pt(0);

 cout << assign_eqn_numbers() << " Equation numbers assigned" << std::endl;
}

//Define the solve function, disp ctl and then continuation
template<class ELEMENT>
void BratuProblem<ELEMENT>::solve()
{
 //Set the desired number of Newton iterations
 Desired_newton_iterations_ds = 2;
 Desired_proportion_of_arc_length = 0.5;
 //Open an output trace file
 ofstream trace("trace.dat");
 double ds = 0.05;
 //Let's have a look at the solution
 for(unsigned i=0;i<5;i++)
  {

   ds = arc_length_step_solve(Global_Physical_Variables::Lambda_pt,ds);
   
   //Output the pressure
   trace << *Global_Physical_Variables::Lambda_pt  << " "
    //Position of first trace node
         << Trace_node_pt->x(0) << " " << Trace_node_pt->value(0) << " " 
         << std::endl;
  }
 trace.close();

 trace.open("trace_mu.dat");

 using namespace Global_Physical_Variables;
 *Mu_pt = 2.0*atan(1.0);





 activate_fold_tracking(Global_Physical_Variables::Lambda_pt);





 newton_solve();




 std::cout << "Fold at " << *Global_Physical_Variables::Lambda_pt << std::endl;
 trace << *Mu_pt << " " << *Lambda_pt << " "
       << Trace_node_pt->x(0) << " " << Trace_node_pt->value(0) << std::endl;
 



 reset_arc_length_parameters();
  
 Desired_proportion_of_arc_length = 0.9;
 Desired_newton_iterations_ds  = 2;

 ds = 0.01;
 for(unsigned i=0;i<15;i++)
  {

   ds = arc_length_step_solve(Mu_pt,ds);
   trace << *Mu_pt << " " << *Lambda_pt << " "
         << Trace_node_pt->x(0) << " " << Trace_node_pt->value(0) << std::endl;
  }

 trace.close();

 deactivate_bifurcation_tracking();
}

//Set up and solve the problem
int main()
{
 //Length of domain
 //Set up the problem
 BratuProblem<GelfandBratuElement<3> > 
  problem(100);

 //Compute derivatives with respect to parameters mu and lambda
 //analytically
 problem.set_analytic_dparameter(Global_Physical_Variables::Mu_pt);
 problem.set_analytic_dparameter(Global_Physical_Variables::Lambda_pt);
 //Solve the problem
 problem.solve();
}







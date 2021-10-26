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
//Driver code to detect and track a Hopf bifurcation problem

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
#include <cmath>

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

class PredatorPreyElement : public GeneralisedElement 
{
 double *Lambda_pt;

 double *Mu_pt;
 
 unsigned Internal_index;

 public:

 PredatorPreyElement() 
  {
   Internal_index = this->add_internal_data(new Data(3));
  }

 //Interface to the parameter
 const double &lambda() const {return *Lambda_pt;}

 const double &mu() const {return *Mu_pt;}

 //Set the pointer
 double* &lambda_pt() {return Lambda_pt;}

 double* &mu_pt() {return Mu_pt;}

  /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix arguments
   fill_in_generic_residual_contribution(
    residuals,GeneralisedElement::Dummy_matrix,
    GeneralisedElement::Dummy_matrix,0);
  }
 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,
                                         jacobian,
                                         GeneralisedElement::Dummy_matrix,1);
  }
 

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the generic routine with the flag set to 2
   fill_in_generic_residual_contribution(residuals,jacobian,mass_matrix,2);
  }
 
 /// Calculate the elemental contributions to the global 
 /// residual vector for the weak form of the Poisson equation
 void fill_in_generic_residual_contribution(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian,
                                            DenseMatrix<double> &mass_matrix,
                                            unsigned flag)
  {
   //Set the mass matrix
   if(flag==2)
    {
     for(unsigned i=0;i<3;i++) {mass_matrix(i,i) = 1.0;}
    }

    //Get the value of the parameter
    double lam = lambda();
    double m = mu();

    //Get the values of the unknowns
    double u[3];
    for(unsigned i=0;i<3;i++) 
     {u[i] = internal_data_pt(Internal_index)->value(i);}

    unsigned local_eqn = internal_local_eqn(Internal_index,0);
    residuals[local_eqn] = u[0]*(1.0-u[0]) - m*u[0]*u[1];
    
    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,0);
      jacobian(local_eqn,local_unknown) = 1.0 - 2.0*u[0] - m*u[1];
      
      local_unknown = internal_local_eqn(Internal_index,1);
      jacobian(local_eqn,local_unknown) = -m*u[0];
     }


    local_eqn = internal_local_eqn(Internal_index,1);
    residuals[local_eqn] = 
     -0.25*u[1] + m*u[0]*u[1] - 3.0*u[1]*u[2] - lam*(1.0 - exp(-5.0*u[1]));

    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,0);
      jacobian(local_eqn,local_unknown) =
       m*u[1];

      local_unknown = internal_local_eqn(Internal_index,1);
      jacobian(local_eqn,local_unknown) =
       -0.25 + m*u[0] - 3.0*u[2] - lam*5.0*exp(-5.0*u[1]);

      local_unknown = internal_local_eqn(Internal_index,2);
      jacobian(local_eqn,local_unknown) = -3.0*u[1];
     }

    local_eqn = internal_local_eqn(Internal_index,2);
    residuals[local_eqn]  =
     -0.5*u[2] + 3.0*u[1]*u[2];

    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,1);
      jacobian(local_eqn,local_unknown) = 3.0*u[2];

      local_unknown = internal_local_eqn(Internal_index,2);
      jacobian(local_eqn,local_unknown) = -0.5 + 3.0*u[1];
     }

   } //End of function

  //Define an output function for the element 
  void output(ostream &output) 
   {
    for(unsigned n=0;n<3;n++)
     {
      output << this->internal_data_pt(Internal_index)->value(n) << " ";
     }
    output << std::endl;
   } //End of function

}; //End of the class


//======================================================================
//Problem class to solve the deformation of an elastic tube
//=====================================================================
template<class ELEMENT>
class PredPreyProblem : public Problem
{

public:
 //Constructor
 PredPreyProblem();

 //Update functions are both empty
 void actions_after_newton_solve() {}
 void actions_before_newton_solve() {}
 void actions_before_newton_convergence_check() {}

 void solve();
 
};

//Constructor
template<class ELEMENT>
PredPreyProblem<ELEMENT>::PredPreyProblem()
{
 eigen_solver_pt() = new LAPACK_QZ();

 //Now create the mesh
 Problem::mesh_pt() = new Mesh;

 //Single element
 ELEMENT* elem_pt = new ELEMENT;
 mesh_pt()->add_element_pt(elem_pt);

 //Assign the physical parameters
 using namespace Global_Physical_Variables;
 //Set the pointer to the external pressure 
 Lambda_pt = new double(0.5); //0.7
 Mu_pt = new double(4.0);//(5.0);

 //Set the load function
 elem_pt->lambda_pt() = Lambda_pt;
 elem_pt->mu_pt() = Mu_pt;

 //elem_pt->internal_data_pt(0)->set_value(0,0.5);
 //elem_pt->internal_data_pt(0)->set_value(1,0.125);
 //elem_pt->internal_data_pt(0)->set_value(2,0.5);

 elem_pt->internal_data_pt(0)->set_value(0,0.25);
 elem_pt->internal_data_pt(0)->set_value(1,0.125);
 elem_pt->internal_data_pt(0)->set_value(2,-0.1);


 cout << assign_eqn_numbers() << " Equation numbers assigned" << std::endl;
}

//Define the solve function, disp ctl and then continuation
template<class ELEMENT>
void PredPreyProblem<ELEMENT>::solve()
{
 //Assign memory for the eigenvalues and eigenvectors
 Vector<complex<double> > eigenvalues;
 Vector< DoubleVector > eigenvectors;

 Desired_newton_iterations_ds = 2;
 Desired_proportion_of_arc_length = 0.5;
 //Open an output trace file
 ofstream trace("trace.dat");
 double ds = -0.05;
 //Let's have a look at the solution
 for(unsigned i=0;i<4/*15*/;i++)
  {
   //*Global_Physical_Variables::Lambda_pt += 0.1;
   //newton_solve();
   /*ds = */arc_length_step_solve(Global_Physical_Variables::Lambda_pt,ds);

   //Output the pressure
   trace << *Global_Physical_Variables::Lambda_pt  << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
         << " " 
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
         << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
         << std::endl;
   
   solve_eigenproblem(3,eigenvalues,eigenvectors);
   
   for(unsigned e=0;e<eigenvalues.size();e++)
    {
     std::cout << eigenvalues[e] << std::endl;
    }
  }
 
 trace.close();

 Max_newton_iterations = 10;
 for(unsigned i=0;i<3;i++) {std::cout << eigenvectors[0][i] << " "	
        << eigenvectors[1][i] << std::endl;}
 //activate_hopf_tracking(Global_Physical_Variables::Lambda_pt);
 activate_hopf_tracking(Global_Physical_Variables::Lambda_pt,
                        eigenvalues[0].imag(),
                        eigenvectors[0],eigenvectors[1]);



 newton_solve();

 std::cout << "Hopf bifurcation found at " 
           << *Global_Physical_Variables::Lambda_pt  << " "
           << dof(ndof()-1) << " "
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
           << " " 
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
           << " "
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
           << std::endl;

 trace.open("trace_hopf.dat");
 //Output the details
 trace << *Global_Physical_Variables::Lambda_pt  << " "
       << *Global_Physical_Variables::Mu_pt << " "
       << dof(ndof() -1) << " "
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
       << " " 
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
       << " "
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
       << std::endl;

 //Let's track the old hopfer
 ds = -0.01;
 //Let's have a look at the solution
 for(unsigned i=0;i<90;i++)
  {
   /*ds = */arc_length_step_solve(Global_Physical_Variables::Mu_pt,ds);

   //Output the pressure
   trace << *Global_Physical_Variables::Lambda_pt  << " "
         << *Global_Physical_Variables::Mu_pt << " "
         << dof(ndof() -1) << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
         << " " 
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
         << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
         << std::endl;
  }
 
 trace.close();
 
}

//Set up and solve the problem
int main()
{
 //Length of domain
 //Set up the problem
 PredPreyProblem<PredatorPreyElement> problem;
 //Solve the problem
 problem.solve();
}







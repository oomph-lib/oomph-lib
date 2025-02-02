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
//Driver function for simple periodic orbit testing
//based on the A-> B -> C chemical reactions in a tank

//Standard system includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio>

//Include files from the finite-element library
#include "generic.h"
#include <cmath>

using namespace std;

using namespace oomph;

//======================================================================== 
/// Global variables that represent physical properties
//======================================================================== 
namespace Global_Physical_Variables
{
 //Parameters shamelessly stolen from Goevearts et al
 Vector<double> P(5);
}


class ABCElement : public virtual GeneralisedElement, 
public virtual PeriodicOrbitBaseElement 
{
 //Pointer to vector of physical variables
 Vector<double> *P_pt;
 
 unsigned Internal_index;

 public:

 ABCElement() 
  {
   //Internal_index = this->add_internal_data(new Data(3));
  }

 //Switch for previous data


 //Contruct the internal data
 void construct_internal_data(TimeStepper* const &time_stepper_pt)
  {
   Internal_index = this->add_internal_data(new Data(time_stepper_pt,3));
  }

 //Interface to the parameter
 inline const double &p(const unsigned &i) const {return (*P_pt)[i];}

 //Set the pointer
 Vector<double>* &p_pt() {return P_pt;}

 /// Switch


 /// Interface to get the current value of all (internal and shared) unknowns
 void get_non_external_dofs(Vector<double> &u)
  {
   Vector<double> val(3);
   internal_data_pt(Internal_index)->time_stepper_pt()->
    time_derivative(0,internal_data_pt(Internal_index),val);

   for(unsigned i=0;i<3;i++)
    {
     unsigned local_eqn = internal_local_eqn(Internal_index,i);
     u[local_eqn] = val[i];
    }
  }


 /// Interface to get the current value of the time derivative of
 /// all (internal and shared) unknowns
 void get_non_external_ddofs_dt(Vector<double> &du_dt)
  {
   Vector<double> dval_dt(3);
   internal_data_pt(Internal_index)->time_stepper_pt()->
    time_derivative(1,internal_data_pt(Internal_index),dval_dt);
   
   for(unsigned i=0;i<3;i++)
    {
     unsigned local_eqn = internal_local_eqn(Internal_index,i);
     du_dt[local_eqn] = dval_dt[i];
    }
  }

 /// Get the inner product matrix
 void get_inner_product_matrix(DenseMatrix<double> &inner_product)
  {
   inner_product.initialise(0.0);
   for(unsigned i=0;i<3;i++) {inner_product(i,i) = 1.0;}
  }


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

    //Get the value of the parameters
   double p0 = p(0);
   double p1 = p(1);
   double p2 = p(2);
   double p3 = p(3);
   double p4 = p(4);

    //Get the values of the unknowns at the current time
    Vector<double> u(3);
    internal_data_pt(Internal_index)
     ->time_stepper_pt()->time_derivative(0,
                                          internal_data_pt(Internal_index),u);
 
    //Get the values of the derivatives of the unknowns at the current time
    Vector<double> du_dt(3);
    internal_data_pt(Internal_index)
     ->time_stepper_pt()->time_derivative(1,
                                          internal_data_pt(Internal_index),
                                          du_dt);
    

    //  for(unsigned i=0;i<3;i++) 
    // {u[i] = internal_data_pt(Internal_index)->value(i);}
    double exp_u2 = exp(u[2]);

    unsigned local_eqn = internal_local_eqn(Internal_index,0);
    residuals[local_eqn] = -du_dt[0] -u[0] + p0*(1.0 - u[0])*exp_u2;
    
    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,0);

      //Don't return time derivative terms if also asked for the mass matrix
      if(flag==2)
       {
        jacobian(local_eqn,local_unknown) =
         -1.0 - p0*exp_u2;
       }
      else
       {
        jacobian(local_eqn,local_unknown) =
         - internal_data_pt(Internal_index)->time_stepper_pt()->weight(1,0)
         -1.0 - p0*exp_u2;
       }
      
      local_unknown = internal_local_eqn(Internal_index,2);
      jacobian(local_eqn,local_unknown) = p0*(1.0 - u[0])*exp_u2;
     }


    local_eqn = internal_local_eqn(Internal_index,1);
    residuals[local_eqn] = -du_dt[1] - u[1] + 
     p0*(1.0 - u[0] - p4*u[1])*exp_u2;

    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,0);
      jacobian(local_eqn,local_unknown) = - p0*exp_u2;

      local_unknown = internal_local_eqn(Internal_index,1);
      //Don't return time derivative terms if also asked for the mass matrix
      if(flag==2)
       {
        jacobian(local_eqn,local_unknown) =
         -1.0 - p0*p4*exp_u2;
       }
      else
       {
        jacobian(local_eqn,local_unknown) =
       - internal_data_pt(Internal_index)->time_stepper_pt()->weight(1,0)
         -1.0 - p0*p4*exp_u2;
       }

      local_unknown = internal_local_eqn(Internal_index,2);
      jacobian(local_eqn,local_unknown) = 
       p0*(1.0 - u[0] - p4*u[1])*exp_u2;
     }

    local_eqn = internal_local_eqn(Internal_index,2);
    residuals[local_eqn]  = -du_dt[2]
     -u[2] - p2*u[2] + p0*p3*(1.0 - u[0] + p1*p4*u[1])*exp_u2;

    if(flag)
     {
      unsigned local_unknown = internal_local_eqn(Internal_index,0);
      jacobian(local_eqn,local_unknown) = -p0*p3*exp_u2;

      local_unknown = internal_local_eqn(Internal_index,1);
      jacobian(local_eqn,local_unknown) = p0*p1*p3*p4*exp_u2;

      local_unknown = internal_local_eqn(Internal_index,2);
      //Don't add the time derivative terms if also asked for the mass
      //matrix
      if(flag==2)
       {
        jacobian(local_eqn,local_unknown) = 
         -1.0 - p2 
         + p0*p3*(1.0 - u[0] + p1*p4*u[1])*exp_u2;
       }
      else
       {
        jacobian(local_eqn,local_unknown) = 
         - internal_data_pt(Internal_index)->time_stepper_pt()->weight(1,0)
         -1.0 - p2 
         + p0*p3*(1.0 - u[0] + p1*p4*u[1])*exp_u2;
       }
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

 //Ouput for spacetime problem
 void spacetime_output(std::ostream &outfile, const unsigned &Nplot,
                       const double &time=0.0)
  {
   outfile << time << " ";

   Vector<double> u(3), du_dt(3);
   this->internal_data_pt(Internal_index)->time_stepper_pt()
    ->time_derivative(0,internal_data_pt(Internal_index),u);
   this->internal_data_pt(Internal_index)->time_stepper_pt()
    ->time_derivative(1,internal_data_pt(Internal_index),du_dt);

   for(unsigned i=0;i<3;i++) {outfile << u[i] << " ";}
   for(unsigned i=0;i<3;i++) {outfile << du_dt[i] << " ";}
   
   
   outfile << std::endl;
  }

}; //End of the class


//======================================================================
//Problem class to solve the A->B->C reaction equations
//=====================================================================
template<class ELEMENT,class TIMESTEPPERT>
class ABCProblem : public Problem
{

 unsigned count;
public:
 //Constructor
 ABCProblem();

 /// Make a copy for using in bifurcation tracking
 Problem* make_copy()
  {
   //Make a copy based on the current parameters
   return(new ABCProblem());
  }

 //Update functions are both empty
 void actions_after_newton_solve() {}
 void actions_before_newton_solve() {}
 void actions_before_newton_convergence_check() {}

 void solve();
 
};

//Constructor
template<class ELEMENT, class TIMESTEPPER>
ABCProblem<ELEMENT,TIMESTEPPER>::ABCProblem() 
{
 count=0;
 //Asssign the timestepper
 add_time_stepper_pt(new TIMESTEPPER);

 eigen_solver_pt() = new LAPACK_QZ();

 //Now create the mesh
 Problem::mesh_pt() = new Mesh;

 //Single element
 ELEMENT* elem_pt = new ELEMENT;
 //Create the internal data
 elem_pt->construct_internal_data(this->time_stepper_pt());
 //Add the element to the mesh
 mesh_pt()->add_element_pt(elem_pt);

 //Assign the physical parameters
 using namespace Global_Physical_Variables;

 //Set the load function
 elem_pt->p_pt() = &P;

 elem_pt->internal_data_pt(0)->set_value(0,0.5);
 elem_pt->internal_data_pt(0)->set_value(1,0.5);
 elem_pt->internal_data_pt(0)->set_value(2,2.0);


 cout << assign_eqn_numbers() << " Equation numbers assigned" << std::endl;
}

//Define the solve function, disp ctl and then continuation
template<class ELEMENT, class TIMESTEPPER>
void ABCProblem<ELEMENT,TIMESTEPPER>::solve()
{
 //Assign memory for the eigenvalues and eigenvectors
 Vector<std::complex<double> > eigenvalues;
 Vector<DoubleVector> eigenvector_real;
 Vector<DoubleVector> eigenvector_imag;

 Desired_newton_iterations_ds = 2;
 Desired_proportion_of_arc_length = 0.5;
 //Open an output trace file
 ofstream trace("trace.dat");
 double ds = 0.01;
 //Let's have a look at the solution
 for(unsigned i=0;i<13;i++)
  {
   //Take some arc-length steps, but ignore the return value
   /*ds = */arc_length_step_solve(&Global_Physical_Variables::P[0],ds);

   //Output the pressure
   trace << Global_Physical_Variables::P[0]  << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
         << " " 
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
         << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
         << std::endl;
   
   this->solve_eigenproblem(3,eigenvalues,eigenvector_real,eigenvector_imag);
   
   for(unsigned e=0;e<eigenvalues.size();e++)
    {
     std::cout << eigenvalues[e] << std::endl;
    }
  }
 
 trace.close();

 //Now let's track the Hopf 
 activate_hopf_tracking(&Global_Physical_Variables::P[0],
                        eigenvalues[0].imag(),
                        eigenvector_real[0],eigenvector_imag[0]);

 this->steady_newton_solve();

 std::cout << "Hopf bifurcation found at " 
           << Global_Physical_Variables::P[0]  << " "
           << dof(ndof()-1) << " "
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
           << " " 
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
           << " "
           << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
           << std::endl;

 //Store the frequency
 double omega = dof(ndof()-1);

 //Go back to normal problem
 this->deactivate_bifurcation_tracking();

 //Increment the parameter a little bit
 Global_Physical_Variables::P[0] += 0.001;

 //Solve for the steady state
 this->steady_newton_solve();
 
 unsigned n_time_element=50;//200;

 //Find out how many time points we'd want
 Vector<double> time_point;
 {
  Mesh* temp_mesh_pt = 
   new OneDMesh<QSpectralElement<1,7> >(n_time_element,1.0);
  const unsigned n_node = temp_mesh_pt->nnode();
  time_point.resize(n_node);
  for(unsigned n=0;n<n_node;n++)
   {
    time_point[n] = temp_mesh_pt->node_pt(n)->x(0);
   }
 }


 //Now assume we've been where we are for all time
 double period = (8.0*atan(1.0))/omega;
 double dt = 0.05*period;
 assign_initial_values_impulsive(dt);

 //Add a tiny kick
 double u = mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0);
 mesh_pt()->element_pt(0)->internal_data_pt(0)->set_value(0,u*1.01);
 

 trace.open("time_trace.dat");
 //Output the details
 trace << time() << " "
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
       << " " 
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
       << " "
       << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
       << std::endl;


 //Integrate over 10 putative periods and then start storing stuff
 
 //Let's have a look at the solution
 for(unsigned i=0;i<200;i++)
  {
   unsteady_newton_solve(dt);

   //Output the pressure
   trace << time() << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0) 
         << " " 
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(1) 
         << " "
         << mesh_pt()->element_pt(0)->internal_data_pt(0)->value(2) 
         << std::endl;
  }
 
 trace.close();
 
 //OK let's do this
 const unsigned n_time = time_point.size();
 const unsigned n_dof = this->ndof();
 DenseMatrix<double> initial_guess(n_time,n_dof);
 for(unsigned n=0;n<n_dof;n++)
  {
   initial_guess(0,n) = this->dof(n);
  }
 

 //Take a timestep
 //period = 8.4;
 for(unsigned t=1;t<n_time;t++)
  {
   dt = period*(time_point[t] - time_point[t-1]);
   unsteady_newton_solve(dt);
   for(unsigned n=0;n<n_dof;n++)
    {
     initial_guess(t,n) = this->dof(n);
    }
  }
 
 //Have now setup the initial guess time to try the timestepping again
 
 //Pepare a periodic orbit handler
 assembly_handler_pt() = new PeriodicOrbitAssemblyHandler<7>(this,
                                                             n_time_element,
                                                             initial_guess,
                                                             8.0*atan(1.0)/
                                                             period);

 //Now let's setup the initial values
 this->newton_solve();

 Max_newton_iterations = 100;

 std::cout << "Orbit found with period " << 1.0/
  this->dof(this->ndof()-1) << "\n";

 //Output
 std::ofstream junk("first_orbit.dat");
 dynamic_cast<PeriodicOrbitAssemblyHandler<7>*>(assembly_handler_pt())
  ->orbit_output(junk,5);
 junk.close();

 dynamic_cast<PeriodicOrbitAssemblyHandler<7>*>(
  assembly_handler_pt())->set_previous_dofs_to_current_dofs();

 std::ofstream orbit_trace("orbit_trace.dat");

 Max_newton_iterations = 10;  

 Desired_proportion_of_arc_length = 0.01;
 //Now we can continue the orbit!
 char filename[100];
 ds = 0.0005;
 unsigned count=0;
 for(unsigned i=0;i<2;i++)
  {
   ++count;
   std::cout << "Taking ds " << ds << "\n";
   ds = this->arc_length_step_solve(&Global_Physical_Variables::P[0],ds);

   dynamic_cast<PeriodicOrbitAssemblyHandler<7>*>(
    assembly_handler_pt())->set_previous_dofs_to_current_dofs();
   sprintf(filename,"orbit%g_%g.dat",Global_Physical_Variables::P[0],
           1.0/this->dof(this->ndof()-1));
   std::ofstream crap(filename);
   dynamic_cast<PeriodicOrbitAssemblyHandler<7>*>(assembly_handler_pt())
    ->orbit_output(crap,5);
   crap.close();
   orbit_trace << Global_Physical_Variables::P[0] << " "
               << this->dof(this->ndof()-1) << " " 
               << this->mesh_pt()->element_pt(0)->internal_data_pt(0)->value(10,0) << " "
               << this->mesh_pt()->element_pt(0)->internal_data_pt(0)->value(10,1) << " "
               << this->mesh_pt()->element_pt(0)->internal_data_pt(0)->value(10,2) << std::endl;
  }

 orbit_trace.close();
}



//Set up and solve the problem
int main()
{
 //Set the physical parameters
 using namespace Global_Physical_Variables;
 P[0] = 0.1;
 P[1] = 1.0;
 P[2] = 1.5;
 P[3] = 8.0;
 P[4] = 0.04;

 //Length of domain
 //Set up the problem
 ABCProblem<ABCElement,BDF<2> > problem;
 //Solve the problem
 problem.steady_newton_solve();
 problem.solve();
}







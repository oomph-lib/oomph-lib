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
// Driver 
#include <fenv.h> 

// The oomphlib headers
#include "generic.h"
#include "axisym_linear_elasticity.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===start_of_Global_Parameters_namespace===============================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Define Poisson's ratio Nu
 double Nu = 0.3;

 /// Define the non-dimensional Young's modulus
 double E = 1.0;

 /// Lame parameters
 double Lambda = E*Nu/(1.0+Nu)/(1.0-2.0*Nu);
 double Mu = E/2.0/(1.0+Nu);

 /// Square of the frequency of the time dependence
 double Omega_sq = 0.5;

 /// Number of elements in r-direction
 unsigned Nr = 5;

 /// Number of elements in z-direction
 unsigned Nz = 10;

 /// Length of domain in r direction
 double Lr = 1.0;

 /// Length of domain in z-direction
 double Lz = 2.0;

 /// Set up min r coordinate
 double Rmin = 0.1;

 /// Set up min z coordinate
 double Zmin = 0.3;

 /// Set up max r coordinate
 double Rmax = Rmin+Lr;

 /// Set up max z coordinate
 double Zmax = Zmin+Lz;
 
 /// The traction function at r=Rmin: (t_r, t_z, t_theta)
 void boundary_traction(const double &time,
                        const Vector<double> &x,
                        const Vector<double> &n,
                        Vector<double> &result)
 {
  result[0] = cos(time)*(-6.0*pow(x[0],2)*Mu*cos(x[1])-
   Lambda*(4.0*pow(x[0],2)+pow(x[0],3))*cos(x[1]));
  result[1] = cos(time)*(-Mu*(3.0*pow(x[0],2)-pow(x[0],3))*sin(x[1]));
  result[2] = cos(time)*(-Mu*pow(x[0],2)*(2*pow(x[1],3)));
 }
 
 ///  The body force function; returns vector of doubles
 /// in the order (b_r, b_z, b_theta)
 void body_force(const double &time,
                 const Vector<double> &x,
                 Vector<double> &result)
 {
  result[0] = cos(time)*(
   x[0]*(-cos(x[1])*
         (Lambda*(8.0+3.0*x[0])-
          Mu*(-16.0+x[0]*(x[0]-3.0))+pow(x[0],2)*Omega_sq)));
  result[1] = cos(time)*(
   x[0]*sin(x[1])*(Mu*(-9.0)+
                   4.0*x[0]*(Lambda+Mu)+pow(x[0],2)*
                   (Lambda+2.0*Mu-Omega_sq)));
  result[2] = cos(time)*(
   -x[0]*(8.0*Mu*pow(x[1],3)+pow(x[0],2)*(pow(x[1],3)*Omega_sq+6.0*Mu*x[1])));
 } // end of body force
 
 ///  Helper function - spatial components of the exact solution in a
 /// vector. This is necessary because we need to multiply this by different
 /// things to obtain the velocity and acceleration
 /// 0: u_r, 1: u_z, 2: u_theta
 void exact_solution_th(const Vector<double> &x,
                        Vector<double> &u)
 {
  u[0] = pow(x[0],3)*cos(x[1]);
  u[1] = pow(x[0],3)*sin(x[1]);
  u[2] = pow(x[0],3)*pow(x[1],3);
 }


 //Displacement, velocity and acceleration functions which are used to provide
 //initial values to the timestepper. We must provide a separate function for
 //each component in each case.

 /// Calculate the time dependent form of the r-component of displacement
 double u_r(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return cos(time)*displ[0];
  } // end_of_u_r

 /// Calculate the time dependent form of the z-component of displacement
 double u_z(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return cos(time)*displ[1];
  }

 /// Calculate the time dependent form of the theta-component of displacement
 double u_theta(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return cos(time)*displ[2];
  }

 /// Calculate the time dependent form of the r-component of velocity
 double d_u_r_dt(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -sin(time)*displ[0];
  }

 /// Calculate the time dependent form of the z-component of velocity
 double d_u_z_dt(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -sin(time)*displ[1];
  }

 /// Calculate the time dependent form of the theta-component of velocity
 double d_u_theta_dt(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -sin(time)*displ[2];
  }

 /// Calculate the time dependent form of the r-component of acceleration
 double d2_u_r_dt2(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -cos(time)*displ[0];
  }

 /// Calculate the time dependent form of the z-component of acceleration
 double d2_u_z_dt2(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -cos(time)*displ[1];
  }
 
 /// Calculate the time dependent form of the theta-component of acceleration
 double d2_u_theta_dt2(const double &time, const Vector<double> &x)
  {
   Vector<double> displ(3);
   exact_solution_th(x,displ);
   return -cos(time)*displ[2];
  }

 /// The exact solution in a vector:
 /// 0: u_r, 1: u_z, 2: u_theta and their 1st and 2nd derivs
 void exact_solution(const double &time,
                     const Vector<double> &x,
                     Vector<double> &u)
 {
  u[0]=u_r(time,x);
  u[1]=u_z(time,x);
  u[2]=u_theta(time,x);

  u[3]=d_u_r_dt(time,x);
  u[4]=d_u_z_dt(time,x);
  u[5]=d_u_theta_dt(time,x);

  u[6]=d2_u_r_dt2(time,x);
  u[7]=d2_u_z_dt2(time,x);
  u[8]=d2_u_theta_dt2(time,x);
 } // end_of_exact_solution

} // end_of_Global_Parameters_namespace

//===start_of_problem_class=============================================
/// Class to validate time harmonic linear elasticity (Fourier 
/// decomposed)
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class AxisymmetricLinearElasticityProblem : public Problem
{
public:

 ///  Constructor: Pass number of elements in r and z directions,
 /// boundary locations and whether we are doing an impulsive start or not
 AxisymmetricLinearElasticityProblem();

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 /// Actions before implicit timestep
 void actions_before_implicit_timestep()
  {
   // Just need to update the boundary conditions
   set_boundary_conditions();
  }

 ///  Set the initial conditions, either for an impulsive start or
 /// with history values for the time stepper
 void set_initial_conditions();

 /// Set the boundary conditions
 void set_boundary_conditions();
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:
 
 /// Allocate traction elements on the bottom surface
 void assign_traction_elements();
 
 /// Pointer to the bulk mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the mesh of traction elements
 Mesh* Surface_mesh_pt;
}; // end_of_problem_class


//===start_of_constructor=============================================
/// Problem constructor: Pass number of elements in coordinate
/// directions and size of domain.
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
AxisymmetricLinearElasticityProblem<ELEMENT, TIMESTEPPER>::
AxisymmetricLinearElasticityProblem()
{
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER());

 //Now create the mesh
 Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(
   Global_Parameters::Nr,
   Global_Parameters::Nz,
   Global_Parameters::Rmin,
   Global_Parameters::Rmax,
   Global_Parameters::Zmin,
   Global_Parameters::Zmax,
   time_stepper_pt());

 //Create the surface mesh of traction elements
 Surface_mesh_pt=new Mesh;
 assign_traction_elements();
 
 //Set the boundary conditions
 set_boundary_conditions();

 // Complete the problem setup to make the elements fully functional

 // Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Cast to a bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the body force
   el_pt->body_force_fct_pt() = &Global_Parameters::body_force;

   // Set the pointer to Poisson's ratio
   el_pt->nu_pt() = &Global_Parameters::Nu;

   // Set the pointer to non-dim Young's modulus
   el_pt->youngs_modulus_pt() = &Global_Parameters::E;

   // Set the pointer to the Lambda parameter
   el_pt->lambda_sq_pt() = &Global_Parameters::Omega_sq;

  }// end_loop_over_elements

 // Loop over the traction elements
 unsigned n_traction =  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_traction;e++)
  {
   // Cast to a surface element
   AxisymmetricLinearElasticityTractionElement<ELEMENT>*
    el_pt = 
    dynamic_cast<AxisymmetricLinearElasticityTractionElement
    <ELEMENT>* >(Surface_mesh_pt->element_pt(e));
   
   // Set the applied traction
   el_pt->traction_fct_pt() = &Global_Parameters::boundary_traction;
   
  }// end_loop_over_traction_elements
 
 // Add the submeshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Now build the global mesh
 build_global_mesh();

 // Assign equation numbers
 cout << assign_eqn_numbers() << " equations assigned" << std::endl; 


 //{ 
 // // Create animation
 // unsigned ntime=20;
 // double displ_r=0;
 // char filename[30];
 // unsigned n_node=Bulk_mesh_pt->nnode();
 // Node* nod_pt=0;

 // for (unsigned it=0;it<ntime;it++)
 //  {
 //   sprintf(filename,"animation%i.dat",it);
 //   Global_Parameters::Output_stream.open(filename);
 //   double t=(2*MathematicalConstants::Pi)*(double(it)/(ntime-1));
 //   std::cout << t << std::endl;
 //   
 //   // Loop over nodes of the mesh
 //   for(unsigned j=0;j<n_node;j++)
 //    {
 //     nod_pt=Bulk_mesh_pt->node_pt(j);
 //     Vector<double> x(2);
 //     x[0]=nod_pt->x(0);
 //     x[1]=nod_pt->x(1);

 //     displ_r = Global_Parameters::u_r(t,x);
 //     Global_Parameters::Output_stream << x[0] << ' ';
 //     Global_Parameters::Output_stream << x[1] << ' ';
 //     Global_Parameters::Output_stream << displ_r << std::endl;
 //    }
 //   Global_Parameters::Output_stream.close();
 //  }
 //} //end_of_animation

} // end_of_constructor


//===start_of_traction===============================================
/// Make traction elements along the boundary r=Rmin
//===================================================================
template<class ELEMENT, class TIMESTEPPER>
void AxisymmetricLinearElasticityProblem<ELEMENT, TIMESTEPPER>::
assign_traction_elements()
{
 unsigned bound, n_neigh;

 // How many bulk elements are next to boundary 3
 bound=3;
 n_neigh = Bulk_mesh_pt->nboundary_element(bound); 

 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   FiniteElement *traction_element_pt 
    = new AxisymmetricLinearElasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
 
   // Add to mesh
   Surface_mesh_pt->add_element_pt(traction_element_pt);
  }

} // end of assign_traction_elements

//===start_of_set_initial_conditions=================================
/// Set the initial conditions (history values)
//===================================================================
template<class ELEMENT, class TIMESTEPPER>
void AxisymmetricLinearElasticityProblem<ELEMENT, TIMESTEPPER>::
set_initial_conditions()
{
 // Upcast the timestepper to the specific type we have
 TIMESTEPPER* timestepper_pt =
  dynamic_cast<TIMESTEPPER*>(time_stepper_pt());

 // By default do a non-impulsive start and provide initial conditions
 bool impulsive_start=false;

 if(impulsive_start)
  {
   // Number of nodes in the bulk mesh
   unsigned n_node = Bulk_mesh_pt->nnode();

   // Loop over all nodes in the bulk mesh
   for(unsigned inod=0;inod<n_node;inod++)
    {
     // Pointer to node
     Node* nod_pt = Bulk_mesh_pt->node_pt(inod);

     // Get nodal coordinates
     Vector<double> x(2);
     x[0] = nod_pt->x(0);
     x[1] = nod_pt->x(1);

     // Assign zero solution at t=0
     nod_pt->set_value(0,0);
     nod_pt->set_value(1,0);
     nod_pt->set_value(2,0);

     // Set the impulsive initial values in the timestepper
     timestepper_pt->assign_initial_values_impulsive(nod_pt);
    }
  } // end_of_impulsive_start
 else // Smooth start
  {
   // Storage for pointers to the functions defining the displacement,
   // velocity and acceleration components
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
    initial_value_fct(3);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
    initial_veloc_fct(3);
   Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
    initial_accel_fct(3);

   // Set the displacement function pointers
   initial_value_fct[0]=&Global_Parameters::u_r;
   initial_value_fct[1]=&Global_Parameters::u_z;
   initial_value_fct[2]=&Global_Parameters::u_theta;

   // Set the velocity function pointers
   initial_veloc_fct[0]=&Global_Parameters::d_u_r_dt;
   initial_veloc_fct[1]=&Global_Parameters::d_u_z_dt;
   initial_veloc_fct[2]=&Global_Parameters::d_u_theta_dt;

   // Set the acceleration function pointers
   initial_accel_fct[0]=&Global_Parameters::d2_u_r_dt2;
   initial_accel_fct[1]=&Global_Parameters::d2_u_z_dt2;
   initial_accel_fct[2]=&Global_Parameters::d2_u_theta_dt2;

   // Number of nodes in the bulk mesh
   unsigned n_node = Bulk_mesh_pt->nnode();

   // Loop over all nodes in bulk mesh
   for(unsigned inod=0;inod<n_node;inod++)
    {
     // Pointer to node
     Node* nod_pt = Bulk_mesh_pt->node_pt(inod);

     // Assign the history values
     timestepper_pt->assign_initial_data_values(nod_pt,
                                                initial_value_fct,
                                                initial_veloc_fct,
                                                initial_accel_fct);
    } // end_of_loop_over_nodes

   // Paranoid checking of history values
   double err_max=0.0;

   // Loop over all nodes in bulk mesh
   for(unsigned jnod=0;jnod<n_node;jnod++)
    {
     // Pointer to node
     Node* nod_pt=Bulk_mesh_pt->node_pt(jnod);

     // Get nodal coordinates
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Get exact displacements
     double u_r_exact=
      Global_Parameters::u_r(time_pt()->time(),x);
     double u_z_exact=
      Global_Parameters::u_z(time_pt()->time(),x);
     double u_theta_exact=
      Global_Parameters::u_theta(time_pt()->time(),x);

     // Get exact velocities
     double d_u_r_dt_exact=
      Global_Parameters::d_u_r_dt(time_pt()->time(),x);
     double d_u_z_dt_exact=
      Global_Parameters::d_u_z_dt(time_pt()->time(),x);
     double d_u_theta_dt_exact=
      Global_Parameters::d_u_theta_dt(time_pt()->time(),x);

     // Get exact accelerations
     double d2_u_r_dt2_exact=
      Global_Parameters::d2_u_r_dt2(time_pt()->time(),x);
     double d2_u_z_dt2_exact=
      Global_Parameters::d2_u_z_dt2(time_pt()->time(),x);
     double d2_u_theta_dt2_exact=
      Global_Parameters::d2_u_theta_dt2(time_pt()->time(),x);

     // Get Newmark approximations for:
     // zero-th time derivatives
     double u_r_fe=timestepper_pt->time_derivative(0,nod_pt,0);
     double u_z_fe=timestepper_pt->time_derivative(0,nod_pt,1);
     double u_theta_fe=timestepper_pt->time_derivative(0,nod_pt,2);

     // first time derivatives
     double d_u_r_dt_fe=timestepper_pt->time_derivative(1,nod_pt,0);
     double d_u_z_dt_fe=timestepper_pt->time_derivative(1,nod_pt,1);
     double d_u_theta_dt_fe=timestepper_pt->time_derivative(1,nod_pt,2);

     // second time derivatives
     double d2_u_r_dt2_fe=timestepper_pt->time_derivative(2,nod_pt,0);
     double d2_u_z_dt2_fe=timestepper_pt->time_derivative(2,nod_pt,1);
     double d2_u_theta_dt2_fe=timestepper_pt->time_derivative(2,nod_pt,2);

     // Calculate the error as the norm of all the differences between the
     // Newmark approximations and the 'exact' (numerical) expressions
     double error=sqrt(pow(u_r_exact-u_r_fe,2)+
                       pow(u_z_exact-u_z_fe,2)+
                       pow(u_theta_exact-u_theta_fe,2)+
                       pow(d_u_r_dt_exact-d_u_r_dt_fe,2)+
                       pow(d_u_z_dt_exact-d_u_z_dt_fe,2)+
                       pow(d_u_theta_dt_exact-d_u_theta_dt_fe,2)+
                       pow(d2_u_r_dt2_exact-d2_u_r_dt2_fe,2)+
                       pow(d2_u_z_dt2_exact-d2_u_z_dt2_fe,2)+
                       pow(d2_u_theta_dt2_exact-d2_u_theta_dt2_fe,2));

     // If there is an error greater than one previously seen, keep hold of it
     if(error>err_max)
      {
       err_max=error;
      }
    } // end of loop over nodes
   std::cout << "Max error in assignment of initial conditions "
             << err_max << std::endl;
  }
} // end_of_set_initial_conditions

//==start_of_set_boundary_conditions======================================
/// Set the boundary conditions
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void AxisymmetricLinearElasticityProblem<ELEMENT, TIMESTEPPER>::
set_boundary_conditions()
{
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin & set the ones that have Dirichlet 
 // conditions here
 
 // storage for nodal position
 Vector<double> x(2);

 // Storage for prescribed displacements
 Vector<double> u(9);

 // Storage for pointers to the functions defining the displacement,
 // velocity and acceleration components
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  initial_value_fct(3);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  initial_veloc_fct(3);
 Vector<typename TIMESTEPPER::NodeInitialConditionFctPt>
  initial_accel_fct(3);
 
 // Set the displacement function pointers
 initial_value_fct[0]=&Global_Parameters::u_r;
 initial_value_fct[1]=&Global_Parameters::u_z;
 initial_value_fct[2]=&Global_Parameters::u_theta;
 
 // Set the velocity function pointers
 initial_veloc_fct[0]=&Global_Parameters::d_u_r_dt;
 initial_veloc_fct[1]=&Global_Parameters::d_u_z_dt;
 initial_veloc_fct[2]=&Global_Parameters::d_u_theta_dt;
 
 // Set the acceleration function pointers
 initial_accel_fct[0]=&Global_Parameters::d2_u_r_dt2;
 initial_accel_fct[1]=&Global_Parameters::d2_u_z_dt2;
 initial_accel_fct[2]=&Global_Parameters::d2_u_theta_dt2;
 

 // Now set displacements on boundaries 0 (z=Zmin),
 //------------------------------------------------
 // 1 (r=Rmax) and 2 (z=Zmax)
 //--------------------------
 for (unsigned ibound=0;ibound<=2;ibound++)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++) 
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pinned in r, z and theta
     nod_pt->pin(0);nod_pt->pin(1);nod_pt->pin(2);
     
     // Direct assigment of just the current values...
     bool use_direct_assigment=true;
     if (use_direct_assigment)
      {
       // get r and z coordinates
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       
       // Compute the value of the exact solution at the nodal point
       Global_Parameters::exact_solution(time_pt()->time(),x,u);
       
       // Set the displacements
       nod_pt->set_value(0,u[0]);
       nod_pt->set_value(1,u[1]);
       nod_pt->set_value(2,u[2]);
      }
     // ...or the history values too:
     else
      {
       // Upcast the timestepper to the specific type we have
       TIMESTEPPER* timestepper_pt =
        dynamic_cast<TIMESTEPPER*>(time_stepper_pt());

       // Assign the history values
       timestepper_pt->assign_initial_data_values(nod_pt,
                                                  initial_value_fct,
                                                  initial_veloc_fct,
                                                  initial_accel_fct);
      }
    }
  } // end_of_loop_over_boundary_nodes
} // end_of_set_boundary_conditions

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void AxisymmetricLinearElasticityProblem<ELEMENT, TIMESTEPPER>::
doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts=10; 
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
   doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
   doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                          Global_Parameters::exact_solution);
 some_file.close();

 // Doc error
 double error=0.0;
 double norm=0.0;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
   doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             Global_Parameters::exact_solution, 
                             time_pt()->time(),
                             error,norm);
 some_file.close();

 // Doc error norm:
 cout << "\nNorm of error:    " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;
} // end_of_doc_solution   

//======================================================================
// End of Axisymmetric problem class
//======================================================================



//===start_of_main======================================================
/// Driver code 
//======================================================================
int main(int argc, char* argv[])
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

 // Set up doc info
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Time dependent problem instance
 AxisymmetricLinearElasticityProblem
  <QAxisymmetricLinearElasticityElement<3>, Newmark<1> > problem;

 // Set the initial time to t=0
 problem.time_pt()->time()=0.0;

 // Set and initialise timestep
 double dt=0;

 // If we're validating, use a larger timestep so that we can do fewer steps
 // before reaching interesting behaviour
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   dt=0.1;
  }
 else // Otherwise use a small timestep
  {
   dt=0.01;
  }
 problem.time_pt()->initialise_dt(dt);

 // Set the initial conditions
 problem.set_initial_conditions();

 // Doc the initial conditions and increment the doc_info number
 problem.doc_solution(doc_info);
 doc_info.number()++;

 // Find the number of timesteps to perform
 unsigned nstep=0;

 // If we're validating, only do a few timesteps; otherwise do a whole period
 if(CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   nstep=5;
  }
 else // Otherwise calculate based on timestep
  {
   // Solve for one full period
   double t_max=2*MathematicalConstants::Pi;

   nstep=unsigned(t_max/dt);
  } //end_of_calculate_number_of_timesteps

 // Do the timestepping
 for(unsigned istep=0;istep<nstep;istep++)
  {
   // Solve for this timestep
   problem.unsteady_newton_solve(dt);

   // Doc the solution and increment doc_info number
   problem.doc_solution(doc_info);
   doc_info.number()++;
  }
} // end_of_main


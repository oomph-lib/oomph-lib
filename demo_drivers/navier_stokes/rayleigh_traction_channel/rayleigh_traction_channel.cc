//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver for "Rayleigh channel" driven by an applied traction

// The oomphlib headers
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===start_of_namespace=================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re;

 /// Womersley = Reynolds times Strouhal
 double ReSt;

 /// Flag for long/short run: Default =  perform long run
 unsigned Long_run_flag=1;

 /// Flag for impulsive start: Default = start from exact
 /// time-periodic solution. 
 unsigned Impulsive_start_flag=0;

} // end of namespace


//==start_of_exact_solution=============================================
/// Namespace for exact solution
//======================================================================
namespace ExactSoln
{
 
 /// Exact solution of the problem as a vector
 void get_exact_u(const double& t, const Vector<double>& x, Vector<double>& u)
 {
  double y=x[1];
  // I=sqrt(-1)
  complex<double> I(0.0,1.0);
  // omega
  double omega=2.0*MathematicalConstants::Pi;
  // lambda
  complex<double> lambda(0.0,omega*Global_Parameters::ReSt);
  lambda = I*sqrt(lambda);

  // Solution
  complex<double> sol(
   exp(complex<double>(0.0,omega*t)) * 
   (exp(lambda*complex<double>(0.0,y))-exp(lambda*complex<double>(0.0,-y)))
   /(exp(I*lambda)-exp(-I*lambda)) );
  
  // Extract real part and stick into vector
  u.resize(2);
  u[0]=real(sol);
  u[1]=0.0;

 }

 /// Traction required at the top boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double> &n,
                          Vector<double>& traction)
 {
  double y=x[1];
  // I=sqrt(-1)
  complex<double> I(0.0,1.0);
  // omega
  double omega=2.0*MathematicalConstants::Pi;
  // lambda
  complex<double> lambda(0.0,omega*Global_Parameters::ReSt);
  lambda = I*sqrt(lambda);

  // Solution
  complex<double> sol(
   exp(complex<double>(0.0,omega*t)) * 
   (exp(lambda*complex<double>(0.0,y))+exp(lambda*complex<double>(0.0,-y)))
   *I*lambda/(exp(I*lambda)-exp(-I*lambda)) );
   
  //Extract real part and stick into vector
  traction.resize(2);
  traction[0]=real(sol);
  traction[1]=0.0;
 } 

}  // end of exact_solution


//===start_of_problem_class=============================================
/// Rayleigh-type problem: 2D channel flow driven by traction bc
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class RayleighTractionProblem : public Problem
{
public:

 /// Constructor: Pass number of elements in x and y directions and 
 /// lengths
 RayleighTractionProblem(const unsigned &nx, const unsigned &ny, 
                 const double &lx, const double &ly);

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 /// Actions before timestep (empty)
 void actions_before_implicit_timestep() {}
   
 /// Run an unsteady simulation
 void unsteady_run(DocInfo& doc_info); 
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Set initial condition (incl previous timesteps) according
 /// to specified function. 
 void set_initial_condition();

 /// Create traction elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_traction_elements(const unsigned &b, 
                               Mesh* const &bulk_mesh_pt,
                               Mesh* const &surface_mesh_pt);

private:

 /// Pointer to the "bulk" mesh
 RectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem_class



//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT,class TIMESTEPPER>
RayleighTractionProblem<ELEMENT,TIMESTEPPER>::RayleighTractionProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& ly)
{
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the mesh with periodic boundary conditions in x direction
 bool periodic_in_x=true;
 Bulk_mesh_pt = 
  new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,periodic_in_x,
                                   time_stepper_pt());

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Create prescribed-traction elements from all elements that are 
 // adjacent to boundary 2, but add them to a separate mesh.
 create_traction_elements(2,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 unsigned num_bound=Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // No slip on bottom
     // (DO NOT PIN TOP BOUNDARY, SINCE TRACTION DEFINES ITS VELOCITY!)
     if (ibound==0)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
      }
     // Horizontal outflow on the left (and right -- right bc not
     // strictly necessary because of symmetry)
     else if ((ibound==1) || (ibound==3))
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  } // end loop over boundaries

 //Complete the problem setup to make the elements fully functional
 
 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Parameters::Re;
   el_pt->re_st_pt() = &Global_Parameters::ReSt;
  }

 // Loop over the flux elements to pass pointer to prescribed traction function
 // and pointer to global time object
 n_el=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Upcast from GeneralisedElement to traction element
   NavierStokesTractionElement<ELEMENT> *el_pt = 
    dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed traction function
   el_pt->traction_fct_pt() = 
    &ExactSoln::prescribed_traction;
  }

 //Assgn equation numbers
 cout << assign_eqn_numbers() << std::endl; 

} // end of constructor

//======================start_of_set_initial_condition====================
/// Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighTractionProblem<ELEMENT,TIMESTEPPER>::set_initial_condition()
{ 
 // Backup time in global Time object
 double backed_up_time=time_pt()->time();
         
 // Past history needs to be established for t=time0-deltat, ...
 // Then provide current values (at t=time0) which will also form
 // the initial guess for the first solve at t=time0+deltat
 
 // Vector of exact solution value
 Vector<double> soln(2);
 Vector<double> x(2);

 //Find number of nodes in mesh
 unsigned num_nod = mesh_pt()->nnode();

 // Set continuous times at previous timesteps:
 // How many previous timesteps does the timestepper use?
 int nprev_steps=time_stepper_pt()->nprev_values();
 Vector<double> prev_time(nprev_steps+1);
 for (int t=nprev_steps;t>=0;t--)
  {
   prev_time[t]=time_pt()->time(unsigned(t));
  } 

 // Loop over current & previous timesteps
 for (int t=nprev_steps;t>=0;t--)
  {
   // Continuous time
   double time=prev_time[t];
   cout << "setting IC at time =" << time << std::endl;
   
   // Loop over the nodes to set initial guess everywhere
   for (unsigned n=0;n<num_nod;n++)
    {
     // Get nodal coordinates
     x[0]=mesh_pt()->node_pt(n)->x(0);
     x[1]=mesh_pt()->node_pt(n)->x(1);

     // Get exact solution at previous time
     ExactSoln::get_exact_u(time,x,soln);
     
     // Assign solution
     mesh_pt()->node_pt(n)->set_value(t,0,soln[0]);
     mesh_pt()->node_pt(n)->set_value(t,1,soln[1]);
     
     // Loop over coordinate directions: Mesh doesn't move, so 
     // previous position = present position
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->node_pt(n)->x(t,i)=x[i];
      }
    } 
  }

 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;

} // end of set_initial_condition




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighTractionProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);

 // Write file as a tecplot text object
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 // ...and draw a horizontal line whose length is proportional
 // to the elapsed time
 some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 some_file << "1" << std::endl;
 some_file << "2" << std::endl;
 some_file << " 0 0" << std::endl;
 some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 some_file.close();
 
 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                          ExactSoln::get_exact_u); 
 some_file.close();

 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             ExactSoln::get_exact_u,
                             time_pt()->time(),
                             error,norm); 
 some_file.close();

 // Doc solution and error
 //-----------------------
 cout << "error: " << error << std::endl; 
 cout << "norm : " << norm << std::endl << std::endl;

 // Get time, position and exact soln
 Vector<double> x(2), u(2);
 double time=time_pt()->time();
 Node* node_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(37))->node_pt(1);
 x[0] = node_pt->x(0);
 x[1] = node_pt->x(1);
 ExactSoln::get_exact_u(time,x,u);

 // doc
 Trace_file << time << " " 
            << x[0] << " "
            << x[1] << " "
            << node_pt->value(0) << " "
            << node_pt->value(1) << " "
            << u[0] << " "
            << u[1] << " " 
            << abs(u[0]-node_pt->value(0)) << " "
            << abs(u[1]-node_pt->value(1)) << " " 
            << error << " "
            << norm << " "
            << std::endl;

} // end_of_doc_solution   

//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the 
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighTractionProblem<ELEMENT,TIMESTEPPER>::
create_traction_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                         Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
    NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_traction_elements



//===start_of_unsteady_run=====================================================
/// Unsteady run...
//=============================================================================
template<class ELEMENT,class TIMESTEPPER>
void RayleighTractionProblem<ELEMENT,TIMESTEPPER>::unsteady_run(DocInfo& doc_info)
{

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // initialise Trace file
 Trace_file << "time" << ",      " 
            << "x" << ",         "
            << "y" << ",         "
            << "u_1" << ",       "
            << "u_2" << ",       "
            << "u_exact_1" << ", "
            << "u_exact_2" << ", "
            << "error_1" << ",   "
            << "error_2" << ",   " 
            << "L2 error" << ",  "
            << "L2 norm" << ",   " << std::endl;

 //Set value of dt
 double dt = 0.025;

 if (Global_Parameters::Impulsive_start_flag==1)
  {
   // Initialise all history values for an impulsive start
   assign_initial_values_impulsive(dt);
   cout << "IC = impulsive start" << std::endl;
  }
 else
  {
   // Initialise timestep
   initialise_dt(dt);
   // Set initial conditions.
   set_initial_condition();
   cout << "IC = exact solution" << std::endl;
  } 

 //Now do many timesteps
 unsigned ntsteps=500;

 // If validation run only do 5 timesteps
 if (Global_Parameters::Long_run_flag==0)
  { 
   ntsteps=5; 
   cout << "validation run" << std::endl;
  }

 // Doc initial condition
 doc_solution(doc_info);
 
 // increment counter
 doc_info.number()++;

 //Loop over the timesteps
 for(unsigned t=1;t<=ntsteps;t++)
  {
   cout << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   unsteady_newton_solve(dt);

   //Output the time
   cout << "Time is now " << time_pt()->time() << std::endl;

   // Doc solution
   doc_solution(doc_info);

   // increment counter
   doc_info.number()++;
  }

} // end of unsteady run


//===start_of_main======================================================
/// Driver code for Rayleigh-type problem
//======================================================================
int main(int argc, char* argv[]) 
{


 /// Convert command line arguments (if any) into flags:
 if (argc==1)
  {
   cout << "No command line arguments specified -- using defaults." << std::endl;
  }
 else if (argc==3)
  {
   cout << "Two command line arguments specified:" << std::endl;
   // Flag for long run
   Global_Parameters::Long_run_flag=atoi(argv[1]);
   /// Flag for impulsive start
   Global_Parameters::Impulsive_start_flag=atoi(argv[2]);
  }
 else
  {
   std::string error_message = 
    "Wrong number of command line arguments. Specify none or two.\n";
   error_message +=
    "Arg1: Long_run_flag [0/1]\n";
   error_message +=
    "Arg2: Impulsive_start_flag [0/1]\n";

   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 cout << "Long run flag: " 
      <<  Global_Parameters::Long_run_flag << std::endl;
 cout << "Impulsive start flag: " 
      <<  Global_Parameters::Impulsive_start_flag << std::endl;


 // Set physical parameters:

 // Womersly number = Reynolds number (St = 1)
 Global_Parameters::ReSt = 10.0;
 Global_Parameters::Re = Global_Parameters::ReSt;

 //Horizontal length of domain
 double lx = 1.0;

 //Vertical length of domain
 double ly = 1.0;

 // Number of elements in x-direction
 unsigned nx=5;

 // Number of elements in y-direction
 unsigned ny=10;

 // Solve with Crouzeix-Raviart elements
 {
  // Set up doc info
  DocInfo doc_info;
  doc_info.number()=0;
  doc_info.set_directory("RESLT_CR");
  
  //Set up problem
  RayleighTractionProblem<QCrouzeixRaviartElement<2>,BDF<2> > 
   problem(nx,ny,lx,ly);
  
  // Run the unsteady simulation
  problem.unsteady_run(doc_info);
 }


 // Solve with Taylor-Hood elements
 {
  // Set up doc info
  DocInfo doc_info;
  doc_info.number()=0;
  doc_info.set_directory("RESLT_TH");

  //Set up problem
  RayleighTractionProblem<QTaylorHoodElement<2>,BDF<2> > 
   problem(nx,ny,lx,ly);
  
  // Run the unsteady simulation
  problem.unsteady_run(doc_info);
 }

} // end of main

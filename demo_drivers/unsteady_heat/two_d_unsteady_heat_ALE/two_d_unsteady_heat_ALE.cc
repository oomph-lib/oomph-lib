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
//Driver for adaptive 2D unsteady heat problem in moving domain

//Generic routines
#include "generic.h"

// The unsteady heat equations
#include "unsteady_heat.h"

// Mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;


//============start_of_MyEllipse===========================================
/// Oscillating ellipse
/// \f[ x = (a + \widehat{a} \sin(2\Pi t/T)) \cos(\xi)  \f]
/// \f[ y = (b + \widehat{b} \sin(2\Pi t/T)) \sin(\xi)  \f]
//=========================================================================
class MyEllipse : public GeomObject
{

public:

 /// Constructor:  Pass half axes, amplitudes of their variation, period
 /// of oscillation and pointer to time object.
 MyEllipse(const double& a, const double& b, 
           const double& a_hat, const double& b_hat, 
           const double& period, Time* time_pt) : 
  GeomObject(1,2), A(a), B(b), A_hat(a_hat), B_hat(b_hat), 
  T(period), Time_pt(time_pt) {}

 /// Destructor: Empty
 virtual ~MyEllipse() {}

 /// Current position vector to material point at 
 /// Lagrangian coordinate xi 
 void position(const Vector<double>& xi, Vector<double>& r) const
  {
   // Get current time:
   double time=Time_pt->time();

   // Position vector
   r[0] = (A+A_hat*sin(2.0*MathematicalConstants::Pi*time/T))*cos(xi[0]);
   r[1] = (B+B_hat*sin(2.0*MathematicalConstants::Pi*time/T))*sin(xi[0]);

  } // end of position(...)



 /// Parametrised position on object: r(xi). Evaluated at
 /// previous time level. t=0: current time; t>0: previous
 /// time level.
 void position(const unsigned& t, const Vector<double>& xi,
               Vector<double>& r) const
  {
   // Get current time:
   double time=Time_pt->time(t);
   
   // Position vector
   r[0] = (A+A_hat*sin(2.0*MathematicalConstants::Pi*time/T))*cos(xi[0]);
   r[1] = (B+B_hat*sin(2.0*MathematicalConstants::Pi*time/T))*sin(xi[0]);

  } // end of position(...)


protected:

 /// x-half axis
 double A;

 /// y-half axis
 double B;

 /// Amplitude of variation in x-half axis
 double A_hat;

 /// Amplitude of variation in y-half axis
 double B_hat;

 /// Period of oscillation
 double T;

 /// Pointer to time object
 Time* Time_pt;

}; // end of MyEllipse

/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_TanhSolnForUnsteadyHeat==============================
/// Namespace for exact solution of unsteady heat equation 
/// with sharp step 
//====================================================================
namespace TanhSolnForUnsteadyHeat
{

 /// Parameter for steepness of step
 double Alpha;

 /// Parameter for amplitude of step translation
 double Beta;

 /// Parameter for timescale of step translation
 double Gamma;

 /// Parameter for angle of step
 double TanPhi;

 /// Position of step (x-axis intercept)
 double step_position(const double& time)
 {
  return Beta*tanh(Gamma*cos(0.2E1*MathematicalConstants::Pi*time));
 }

 /// Exact solution as a Vector
 void get_exact_u(const double& time, const Vector<double>& x, 
                  Vector<double>& u)
 {
  double X=x[0];
  double Y=x[1];
  u[0]=tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
       MathematicalConstants::Pi*time)))-Y));
 }

 /// Exact solution as a scalar
 void get_exact_u(const double& time, const Vector<double>& x, double& u)
 {
  double X=x[0];
  double Y=x[1];
  u=tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
         MathematicalConstants::Pi*time)))-Y));
 }

 /// Source function to make it an exact solution 
 void get_source(const double& time, const Vector<double>& x, double& source)
 {
  double X=x[0];
  double Y=x[1];
  source =  -2.0*tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
MathematicalConstants::Pi*time)))-Y))*(1.0-pow(tanh(0.1E1-Alpha*(TanPhi*(X-
Beta*tanh(Gamma*cos(0.2E1*MathematicalConstants::Pi*time)))-Y)),2.0))*Alpha*
Alpha*TanPhi*TanPhi-2.0*tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
MathematicalConstants::Pi*time)))-Y))*(1.0-pow(tanh(0.1E1-Alpha*(TanPhi*(X-
Beta*tanh(Gamma*cos(0.2E1*MathematicalConstants::Pi*time)))-Y)),2.0))*Alpha*
Alpha+0.2E1*(1.0-pow(tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
MathematicalConstants::Pi*time)))-Y)),2.0))*Alpha*TanPhi*Beta*(1.0-pow(tanh(
Gamma*cos(0.2E1*MathematicalConstants::Pi*time)),2.0))*Gamma*sin(0.2E1*
MathematicalConstants::Pi*time)*MathematicalConstants::Pi;
 }

 /// Flux required by the exact solution on a boundary on which y is fixed
 void prescribed_flux_on_fixed_y_boundary(const double& time,
                                          const Vector<double>& x, 
                                          double& flux)
 {
  double X=x[0];
  double Y=x[1];

  //The outer unit normal to the boundary is (0,-1)
  double NX =  0.0;
  double NY = -1.0;

  //The flux in terms of the normal is
  flux =  -(1.0-pow(tanh(0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*
MathematicalConstants::Pi*time)))-Y)),2.0))*Alpha*TanPhi*NX+(1.0-pow(tanh(
0.1E1-Alpha*(TanPhi*(X-Beta*tanh(Gamma*cos(0.2E1*MathematicalConstants::Pi*
time)))-Y)),2.0))*Alpha*NY;
 }

} // end of TanhSolnForUnsteadyHeat


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// Unsteady heat problem in deformable ellipse domain.
//====================================================================
template<class ELEMENT>
class RefineableUnsteadyHeatProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 RefineableUnsteadyHeatProblem(UnsteadyHeatEquations<2>::
                             UnsteadyHeatSourceFctPt source_fct_pt);

 /// Destructor: Close trace file
 ~RefineableUnsteadyHeatProblem();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after timestep (empty)
 void actions_after_implicit_timestep(){}

 /// Update the problem specs before next timestep: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_implicit_timestep();
 
 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();

 /// Set initial condition (incl previous timesteps) according
 /// to specified function.  Note that his overloads the virtual
 /// function in the Problem base class and is therefore executed 
 /// automatically to re-assign the initial conditions during the 
 /// spatially adaptive solution at the first timestep.
 void set_initial_condition();

 /// Create UnsteadyHeat flux elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const &surface_mesh_pt);

 /// Delete UnsteadyHeat flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 /// Doc the solution
 void doc_solution();

 /// Dump problem data to allow for later restart
 void dump_it(ofstream& dump_file);

 /// Read problem data for restart
 void restart(ifstream& restart_file);

 /// Pointer to bulk mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* bulk_mesh_pt()
  {
   return Bulk_mesh_pt;
  }


private:

 /// Pointer to GeomObject that specifies the domain bondary
 GeomObject* Boundary_pt;

 /// Pointer to source function
 UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt Source_fct_pt;

 /// Pointer to the "bulk" mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to central node (exists at all refinement levels) for doc
 Node* Doc_node_pt;

 /// Doc info object
 DocInfo Doc_info;

 /// Trace file
 ofstream Trace_file;

}; // end of problem_class

//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem on deformable ellipse domain.
/// Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineableUnsteadyHeatProblem<ELEMENT>::RefineableUnsteadyHeatProblem(
   UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt) : 
         Source_fct_pt(source_fct_pt)
{ 


 // Setup labels for output
 //------------------------

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0; 
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 
 Trace_file << "VARIABLES=\"time t\",\"u<SUB>FE</SUB>\",\"u<SUB>exact</SUB>\","
            << "\"A\","
            << "\"X<SUB>step</SUB>\","
            << "\"N<SUB>element</SUB>\","
            << "\"N<SUB>refined</SUB>\","
            << "\"N<SUB>unrefined</SUB>\","
            << "\"norm of error\","
            << "\"norm of solution\""
            << std::endl;


 // Setup parameters for tanh solution
 // ----------------------------------

 // Steepness of step
 TanhSolnForUnsteadyHeat::Alpha=10.0;

 // Orientation of step
 TanhSolnForUnsteadyHeat::TanPhi=1.0;

 // Amplitude for movement of step
 TanhSolnForUnsteadyHeat::Beta=0.3; 

 // Parameter for time-dependence of step movement
 TanhSolnForUnsteadyHeat::Gamma=5.0;



 //Allocate the timestepper -- This constructs the time object as well
 add_time_stepper_pt(new BDF<2>());


 // Setup mesh
 //-----------

 // Build geometric object that forms the curvilinear domain boundary:
 // an oscillating ellipse

 // Half axes
 double a=1.0;
 double b=1.0;

 // Variations of half axes
 double a_hat= 0.1;
 double b_hat=-0.1;

 // Period of the oscillation
 double period=1.0;

 // Create GeomObject
 Boundary_pt=new MyEllipse(a,b,a_hat,b_hat,period,Problem::time_pt()); 

 // Start and end coordinates of curvilinear domain boundary on ellipse
 double xi_lo=0.0;
 double xi_hi=MathematicalConstants::Pi/2.0;

 // Now create the bulk mesh. Separating line between the two 
 // elements next to the curvilinear boundary is located half-way
 // along the boundary.
 double fract_mid=0.5;
 Bulk_mesh_pt = new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Boundary_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());

 // Create the surface mesh as an empty mesh
 Surface_mesh_pt=new Mesh;

 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 0 (the horizontal lower boundary), and add them 
 // to the (so far empty) surface mesh.
 create_flux_elements(0,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   // Leave nodes on boundary 0 free -- this is where we apply the flux
   // boundary condition
   if (b!=0)
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    }
  }

 // Extract pointer to the central node (this exists at all refinement levels)
 // for doc of solution
 FiniteElement* el0_pt=Bulk_mesh_pt->finite_element_pt(0);
 unsigned nnod=el0_pt->nnode();
 Doc_node_pt=el0_pt->node_pt(nnod-1);


 // Complete the build of all elements so they are fully functional
 //----------------------------------------------------------------

 // Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Loop over the flux elements to pass pointer to prescribed flux function
 n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to UnsteadyHeat flux element
   UnsteadyHeatFluxElement<ELEMENT> *el_pt = 
    dynamic_cast<UnsteadyHeatFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForUnsteadyHeat::prescribed_flux_on_fixed_y_boundary;
  }

 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//======start_of_destructor===============================================
/// Destructor: Close trace file
//========================================================================
template<class ELEMENT>
RefineableUnsteadyHeatProblem<ELEMENT>::~RefineableUnsteadyHeatProblem()
{ 
 // Close trace file
 Trace_file.close();

} // end of destructor


//=========start of actions_before_implicit_timestep===============================
/// Actions before timestep: Update the domain shape, then set the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Update the domain shape
 Bulk_mesh_pt->node_update();
 
 // Get current time
 double time=time_pt()->time();
 
 //Loop over all boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<num_bound;b++)
  {
   // Exclude flux boundary
   if (b!=0)
    {
     // Loop over the nodes on boundary
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(b);
     for (unsigned j=0;j<num_nod;j++)
      {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
       double u;
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       TanhSolnForUnsteadyHeat::get_exact_u(time,x,u);
       nod_pt->set_value(0,u);
      }
    }
  }
} // end of actions_before_implicit_timestep


//=========start_of_actions_before_adapt==================================
/// Actions before adapt: Wipe the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::actions_before_adapt()
{

 // Kill the flux elements and wipe surface mesh
 delete_flux_elements(Surface_mesh_pt);
 
 // Rebuild the global mesh. 
 rebuild_global_mesh();

} // end of actions_before_adapt


//==========start_of_actions_after_adapt==================================
/// Actions after adapt: Rebuild the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 0 and add them to surface mesh
 create_flux_elements(0,Bulk_mesh_pt,Surface_mesh_pt);

 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Loop over the flux elements to pass pointer to prescribed flux function
 unsigned n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to UnsteadyHeat flux element
   UnsteadyHeatFluxElement<ELEMENT> *el_pt = 
    dynamic_cast<UnsteadyHeatFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForUnsteadyHeat::prescribed_flux_on_fixed_y_boundary;
  }
} // end of actions_after_adapt


//======================start_of_set_initial_condition====================
/// Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::set_initial_condition()
{ 

 // Pointer to restart file
 ifstream* restart_file_pt=0;

 // Restart?
 //---------
 // Restart file specified via command line [all programs have at least
 // a single command line argument: their name. Ignore this here.]
 if (CommandLineArgs::Argc==1)
  {
   cout << "No restart -- setting IC from exact solution" << std::endl;
  }
 else if (CommandLineArgs::Argc==2)
  {
   // Open restart file
   restart_file_pt= new ifstream(CommandLineArgs::Argv[1],ios_base::in);
   if (restart_file_pt!=0)
    {
     cout << "Have opened " << CommandLineArgs::Argv[1] << 
      " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream 
      << "ERROR while trying to open " << CommandLineArgs::Argv[1] << 
      " for restart." << std::endl;
     
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
  }
 // More than one command line argument?
 else 
  {
   std::ostringstream error_stream;
   error_stream << "Can only specify one input file\n" 
                << "You specified the following command line arguments:\n";
   //Fix this
   CommandLineArgs::output();
   
   throw OomphLibError( 
    error_stream.str(),
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }
 
 
 // Read restart data:
 //-------------------
 if (restart_file_pt!=0)
  {
   // Read the data from restart file and find out if the restart file
   // was from an unsteady run
   restart(*restart_file_pt);
  }
 // Assign initial condition from exact solution
 //---------------------------------------------
 else
  {
   // Choose initial timestep
   double dt0=0.01;

   // Initialise timestep -- also sets the weights for all timesteppers
   // in the problem.
   initialise_dt(dt0);

   // Backup time in global timestepper
   double backed_up_time=time_pt()->time();
   
   // Past history for velocities must be established for t=time0-deltat, ...
   // Then provide current values (at t=time0) which will also form
   // the initial guess for first solve at t=time0+deltat
   
   // Vector of exact solution value
   Vector<double> soln(1);
   Vector<double> x(2);
   
   //Find number of nodes in mesh
   unsigned num_nod = Bulk_mesh_pt->nnode();
   
   // Get continuous times at previous timesteps
   int nprev_steps=time_stepper_pt()->nprev_values();
   Vector<double> prev_time(nprev_steps+1);
   for (int itime=nprev_steps;itime>=0;itime--)
    {
     prev_time[itime]=time_pt()->time(unsigned(itime));
    }
   
   // Loop over current & previous timesteps (in outer loop because
   // the mesh also moves!)
   for (int itime=nprev_steps;itime>=0;itime--)
    {
     double time=prev_time[itime];
     
     // Set global time (because this is how the geometric object refers 
     // to continous time 
     time_pt()->time()=time;
     
     cout << "setting IC at time =" << time << std::endl;
     
     // Update the mesh for this value of the continuous time
     // (The wall object reads the continous time from the same
     // global time object)
     Bulk_mesh_pt->node_update(); 
     
     // Loop over the nodes to set initial guess everywhere
     for (unsigned jnod=0;jnod<num_nod;jnod++)
      {
       // Get nodal coordinates
       x[0]=Bulk_mesh_pt->node_pt(jnod)->x(0);
       x[1]=Bulk_mesh_pt->node_pt(jnod)->x(1);
       
       // Get intial solution
       TanhSolnForUnsteadyHeat::get_exact_u(time,x,soln);
       
       // Assign solution
       Bulk_mesh_pt->node_pt(jnod)->set_value(itime,0,soln[0]);
       
       // Loop over coordinate directions
       for (unsigned i=0;i<2;i++)
        {
         Bulk_mesh_pt->node_pt(jnod)->x(itime,i)=x[i];
        }
      } 
    } // end of loop over previous timesteps
   
   // Reset backed up time for global timestepper
   time_pt()->time()=backed_up_time;
  }

} // end of set_initial_condition


//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::doc_solution()
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for t=" << time_pt()->time() << std::endl;
 cout << "=================================================" << std::endl;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 some_file << "1" << std::endl;
 some_file << "2" << std::endl;
 some_file << " 0 0" << std::endl;
 some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 -1.2" << std::endl;
 some_file << "1.3 0.0 -1.2" << std::endl;
 some_file << "0.0 1.3 -1.2" << std::endl;
 some_file << "1.3 1.3 -1.2" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 1.2" << std::endl;
 some_file << "1.3 0.0 1.2" << std::endl;
 some_file << "0.0 1.3 1.2" << std::endl;
 some_file << "1.3 1.3 1.2" << std::endl;

 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                       TanhSolnForUnsteadyHeat::get_exact_u); 

 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 -1.2" << std::endl;
 some_file << "1.3 0.0 -1.2" << std::endl;
 some_file << "0.0 1.3 -1.2" << std::endl;
 some_file << "1.3 1.3 -1.2" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "0.0 0.0 1.2" << std::endl;
 some_file << "1.3 0.0 1.2" << std::endl;
 some_file << "0.0 1.3 1.2" << std::endl;
 some_file << "1.3 1.3 1.2" << std::endl;

 some_file.close();


 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                          TanhSolnForUnsteadyHeat::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();



 // Doc error and write trace file
 //--------------------------------
 cout << "error: " << error << std::endl; 
 cout << "norm : " << norm << std::endl << std::endl;

 Vector<double> x(2);
 x[0]=Doc_node_pt->x(0);
 x[1]=Doc_node_pt->x(1);
 double u_exact;
 TanhSolnForUnsteadyHeat::get_exact_u(time_pt()->time(),x,u_exact);
 Vector<double > xi_wall(1);
 Vector<double > r_wall(2);
 xi_wall[0]=0.0;
 Boundary_pt->position(xi_wall,r_wall);
 Trace_file << time_pt()->time() 
            << " " << Doc_node_pt->value(0)
            << " " << u_exact
            << " " << r_wall[0]
            << " " << TanhSolnForUnsteadyHeat::step_position(time_pt()->time())
            << " " << Bulk_mesh_pt->nelement() 
            << " " << Bulk_mesh_pt->nrefined()
            << " " << Bulk_mesh_pt->nunrefined() 
            << " " << error 
            << " " << norm << std::endl;

 // Plot wall posn
 //---------------
 sprintf(filename,"%s/Wall%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 
 unsigned nplot=100;
 for (unsigned iplot=0;iplot<nplot;iplot++)
  {
   xi_wall[0]=0.5*Pi*double(iplot)/double(nplot-1);
   Boundary_pt->position(xi_wall,r_wall);
   some_file << r_wall[0] << " " << r_wall[1] << std::endl;
  }
 some_file.close();
 
 // Write restart file
 sprintf(filename,"%s/restart%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 dump_it(some_file);
 some_file.close();

 // Increment number of doc
 Doc_info.number()++;


} // end of doc_solution


//============start_of_create_flux_elements==============================
/// Create UnsteadyHeat Flux Elements on the b-th boundary of the Mesh object
/// pointed to by bulk_mesh_pt and add the elements to the Mesh object
/// pointed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::
create_flux_elements(const unsigned& b, Mesh* const &bulk_mesh_pt,
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
   
   //What is the face index of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   UnsteadyHeatFluxElement<ELEMENT>* flux_element_pt = new 
   UnsteadyHeatFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements


//============start_of_delete_flux_elements==============================
/// Delete UnsteadyHeat Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements


//=======start_of_dump_it=================================================
/// Dump the solution to disk
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::dump_it(ofstream& dump_file)
{

 // Write step number
 dump_file << Doc_info.number() << " # step number" << std::endl;

 // Dump the refinement pattern and the generic problem data
 Problem::dump(dump_file);
  
} // end of dump_it

//=========start_of_restart===============================================
/// Read solution from disk
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::restart(ifstream& restart_file)
{
 // Read line up to termination sign
 string input_string;
 getline(restart_file,input_string,'#');

 // Ignore rest of line
 restart_file.ignore(80,'\n');

 // Read in step number
 Doc_info.number()=unsigned(atof(input_string.c_str()));

 // Refine the mesh and read in the generic problem data
 Problem::read(restart_file);

} // end of restart



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======start_of_main=====================================================
/// Demonstrate how to solve an unsteady heat problem in deformable domain
/// with mesh adaptation. Command line arguments specify 
/// the name of the restart file. 
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Build problem
 RefineableUnsteadyHeatProblem<RefineableQUnsteadyHeatElement<2,3> >
  problem(&TanhSolnForUnsteadyHeat::get_source);
   
 // Specify duration of the simulation
// double t_max=3.0;

 // Set targets for spatial adaptivity
 problem.bulk_mesh_pt()->max_permitted_error()=0.001;
 problem.bulk_mesh_pt()->min_permitted_error()=0.0001;

 // First timestep?
 bool first=true;
 
 // Max. number of spatial adaptations per timestep. Allow plenty
 // of adaptations at first timestep as the initial conditions
 // can be reset "exactly" from without any interpolation error.
 unsigned max_adapt=10; 
 
 // Set IC
 problem.set_initial_condition();

 // Initial timestep: Use the one used when setting up the initial
 // condition
 double dt=problem.time_pt()->dt();

 // If restart: The first step isn't really the first step,
 // i.e. initial condition should not be re-set when 
 // adaptive refinement has been performed. Also, limit
 // the max. number of refinements per timestep to the
 // normal value straightaway.
 if (CommandLineArgs::Argc==2)
  {
   first=false;
   max_adapt=1;
  }
 // If no restart, refine mesh uniformly before we get started
 else
  {
   problem.refine_uniformly();
   problem.refine_uniformly();
  }

 //Output solution
 problem.doc_solution();
 
 // find number of steps
 unsigned nstep = 6; //unsigned(t_max/dt);

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Take timestep 
   problem.unsteady_newton_solve(dt,max_adapt,first);
   
   // Now we've done the first timestep -- don't re-set the IC
   // in subsequent steps
   first=false;
   
   // Reduce the number of spatial adaptations to one per 
   // timestep
   max_adapt=1;
   
   //Output solution
   problem.doc_solution();

  }
 

}; // end of main

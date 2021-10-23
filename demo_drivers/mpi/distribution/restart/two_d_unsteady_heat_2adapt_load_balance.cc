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
// Driver for doubly adaptive 2D unsteady heat problem in moving domain
// with restart and load balancing and pruning all in one code

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


//======start_of_GlobalParameters==============================
/// Namespace for exact solution of unsteady heat equation 
/// with sharp step 
//====================================================================
namespace GlobalParameters
{

 /// Name of restart file
 std::string Restart_file="";

 /// Name of file specifying the partitioning of the problem
 std::string Partitioning_file="";

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

} // end of GlobalParameters


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


 /// Constructor: Pass pointer to source function and timestep
 RefineableUnsteadyHeatProblem(
  UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt);

 /// Destructor: Close trace file
 ~RefineableUnsteadyHeatProblem();

 /// Build meshes (helper fct accessed from constructor and
 /// from the load balancing routines)
 void build_mesh();

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
 void actions_before_adapt()
  {
   generic_actions_before();
  }
 
 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt()
  {
   generic_actions_after();
  }

 /// Actions before distribute: Wipe the mesh of prescribed flux elements
 void actions_before_distribute()
  {
   generic_actions_before();
  }

 /// Actions after distribute: Rebuild the mesh of prescribed flux elements
 void actions_after_distribute()
  {
   generic_actions_after();
  }

 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();

 /// Set initial condition (incl previous timesteps) according
 /// to specified function. 
 void set_initial_condition();

 /// Restart
 void restart();

 /// Create UnsteadyHeat flux elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const &surface_mesh_pt);

 /// Delete UnsteadyHeat flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 /// Doc the solution
 void doc_solution(const std::string& comment="");

 /// Return DocInfo object
 DocInfo& doc_info(){return Doc_info;}

 /// Dump problem data to allow for later restart
 void dump_it(ofstream& dump_file);

 /// Read problem data for restart
 void restart(ifstream& restart_file);

 /// Pointer to bulk mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* bulk_mesh_pt()
  {
   return Bulk_mesh_pt;
  }

 /// Target error for adaptive timestepping
 double& epsilon_t() {return Epsilon_t;}
 
 /// Write header for trace file
 void write_trace_file_header();

 /// Suggestion for next timestep (stored to allow it to be written
 /// to (or read from) restart file)
 double& next_dt(){return Next_dt;}

private:

 /// Actions before adapt/distribute: Wipe the mesh of prescribed flux elements
 void generic_actions_before();
 
 /// Actions after adapt/distribute: Rebuild the mesh of prescribed 
 /// flux elements
 void generic_actions_after();

 /// Doc info object
 DocInfo Doc_info;

 /// Suggestion for next timestep (stored to allow it to be written
 /// to (or read from) restart file)
 double Next_dt;
 
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

 /// Trace file
 ofstream Trace_file;

 /// Target error for adaptive timestepping
 double Epsilon_t;
 
}; // end of problem_class

//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem on deformable ellipse domain.
/// Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineableUnsteadyHeatProblem<ELEMENT>::RefineableUnsteadyHeatProblem(
 UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt)
 : Source_fct_pt(source_fct_pt)
{ 

 // Setup labels for output
 //------------------------

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0; 


 // Setup timestepping
 //-------------------

 // Allocate the timestepper -- this constructs the time object as well.
 // The boolean flag to the constructor enables adaptivity.
 add_time_stepper_pt(new BDF<2>(true));

 // Set target error for adaptive timestepping
 Epsilon_t=1.0e-2;

 // Setup parameters for tanh solution
 // ----------------------------------

 // Steepness of step
 GlobalParameters::Alpha=10.0;

 // Orientation of step
 GlobalParameters::TanPhi=1.0;

 // Amplitude for movement of step
 GlobalParameters::Beta=0.3; 

 // Parameter for time-dependence of step movement
 GlobalParameters::Gamma=5.0;

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

 // Now build the mesh
 build_mesh();

 // Do equation numbering
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

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


//========start_of_buil_mesh==============================================
/// Build mesh function
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::build_mesh()
{ 

 // Start and end coordinates of curvilinear domain boundary on ellipse
 double xi_lo=0.0;
 double xi_hi=MathematicalConstants::Pi/2.0;

 // Now create the bulk mesh. Separating line between the two 
 // elements next to the curvilinear boundary is located half-way
 // along the boundary.
 double fract_mid=0.5;
 Bulk_mesh_pt = new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Boundary_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());

 // Note: Mesh is completely rebuilt in here
 // so target errors etc need to be reassigned

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=0.001;
 Bulk_mesh_pt->min_permitted_error()=0.0001;

 // Need to do this in build_mesh() because it's supposed to
 // get the problem into the state it was when it was distributed.
 Bulk_mesh_pt->refine_uniformly();
 Bulk_mesh_pt->refine_uniformly();

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

 // In parallel make sure that we retain this element on all processors
 // so the node remains accessible for post-processing
 el0_pt->must_be_kept_as_halo();

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
    &GlobalParameters::prescribed_flux_on_fixed_y_boundary;
  }

} // end of build_mesh


//=========start of actions_before_implicit_timestep======================
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
       GlobalParameters::get_exact_u(time,x,u);
       nod_pt->set_value(0,u);
      }
    }
  }
} // end of actions_before_implicit_timestep


//=======start_of_generic_actions_before==================================
/// Actions before adapt/distribute: Wipe the mesh of prescribed flux 
/// elements
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::generic_actions_before()
{

 // Kill the flux elements and wipe surface mesh
 delete_flux_elements(Surface_mesh_pt);
 
 // Rebuild the global mesh. 
 rebuild_global_mesh();

} // end of generic_actions_before


//==========start_of_generic_actions_after==================================
/// Actions after adapt/distribute: Rebuild the mesh of prescribed flux
/// elements
//==========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::generic_actions_after()
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
    &GlobalParameters::prescribed_flux_on_fixed_y_boundary;
  }

} // end of generic_actions_after


//======================start_of_set_initial_condition====================
/// Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::set_initial_condition()
{ 

 // Choose initial timestep
 double dt0=0.005;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 initialise_dt(dt0);
 
 // Use this as the first "suggested" timestep
 Next_dt=dt0;
 
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
   
   oomph_info << "setting IC at time =" << time << std::endl;
   
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
     GlobalParameters::get_exact_u(time,x,soln);
     
     // Assign solution
     Bulk_mesh_pt->node_pt(jnod)->set_value(itime,0,soln[0]);
     
     // Loop over coordinate directions
     for (unsigned i=0;i<2;i++)
      {
       Bulk_mesh_pt->node_pt(jnod)->x(itime,i)=x[i];
      }
    } 
  }
 
 // Reset backed up time for global timestepper
 time_pt()->time()=backed_up_time;

} // end of set_initial_condition


//======================start_of_restart==================================
/// Restart
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::restart()
{ 

 // Pointer to restart file
 ifstream* restart_file_pt=0;

 // Open restart file from stem
 char filename[100];
 sprintf(filename,"%s_on_proc%i.dat",GlobalParameters::Restart_file.c_str(),
         this->communicator_pt()->my_rank());
 restart_file_pt= new ifstream(filename,ios_base::in);
 if (restart_file_pt!=0)
  {
   oomph_info << "Have opened " << filename << 
    " for restart. " << std::endl;
  }
 else
  {
   std::ostringstream error_stream;
   error_stream
    << "ERROR while trying to open " << filename
    << " for restart." << std::endl;
   
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

} // end of restart




//=======start_of_write_trace_file_header=================================
/// Write the tecplot header for the trace file
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::write_trace_file_header() 
{
 // Open trace file on root only
 if (communicator_pt()->my_rank()==0)
  {
   char filename[100];   
   sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
   Trace_file.open(filename);
   
   Trace_file << "VARIABLES=\"time t\","
              << "\"u<SUB>FE</SUB>\","
              << "\"u<SUB>exact</SUB>\","
              << "\"A\","
              << "\"dt\","
              << "\"global temporal error norm\","
              << "\"X<SUB>step</SUB>\","
              << "\"N<SUB>element</SUB>\","
              << "\"N<SUB>refined</SUB>\","
              << "\"N<SUB>unrefined</SUB>\","
              << "\"norm of error\","
              << "\"norm of solution\""
              << std::endl;
   Trace_file << "ZONE T=\"" << time_stepper_pt()->type() 
              << time_stepper_pt()->nprev_values()
              << "; initial dt="
              <<  time_pt()->dt() << "; ";
   if (time_stepper_pt()->adaptive_flag())
    {
     Trace_file << "adaptive; target error " << Epsilon_t 
                << " \"" << std::endl;
    }
   else
    {
     Trace_file << "(fixed timestep)\"" << std::endl;
    }
  }

} // end of write_trace_file_header



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::doc_solution(const 
                                                          std::string& comment)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 oomph_info << std::endl;
 oomph_info << "=================================================" 
            << std::endl;
 oomph_info << "Docing solution " << Doc_info.number() << " for t=" 
      << time_pt()->time() << " [ndof=" << ndof() << "]\n";
 oomph_info << "=================================================" 
            << std::endl;

 // Doc halo schemes 
 //-----------------
 Bulk_mesh_pt->doc_mesh_distribution(Doc_info);

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "; "<< comment << "\"";
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
 sprintf(filename,"%s/exact_soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,time_pt()->time(),
                       GlobalParameters::get_exact_u); 

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

 sprintf(filename,"%s/error%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                          GlobalParameters::get_exact_u,
                          time_pt()->time(),
                          error,norm); 
 some_file.close();


 // Collect error and norm from all processors
 Vector<double> local_stuff(2);
 local_stuff[0]=error;
 local_stuff[1]=norm;
 Vector<double> global_stuff(2);
#ifdef OOMPH_HAS_MPI
 MPI_Allreduce(&local_stuff[0],&global_stuff[0],2,MPI_DOUBLE,MPI_SUM,
               communicator_pt()->mpi_comm());
#else
 global_stuff[0]=error;
 global_stuff[1]=norm;
#endif

 // Doc solution and error
 //-----------------------
 oomph_info << "error: " << global_stuff[0] << std::endl; 
 oomph_info << "norm : " << global_stuff[1] << std::endl << std::endl;
 
 // Get global error norm (on all processors as it involves comms)
 double temp_error_norm=global_temporal_error_norm();
 
 // Gather refinement stats
 Vector<int> ref_stats(3);
 ref_stats[0]=Bulk_mesh_pt->nnon_halo_element();
 ref_stats[1]=Bulk_mesh_pt->nrefined();
 ref_stats[2]=Bulk_mesh_pt->nunrefined();
 Vector<int> global_ref_stats(3);
 MPI_Allreduce(&ref_stats[0],&global_ref_stats[0],3,MPI_INT,MPI_SUM,
               communicator_pt()->mpi_comm());
 
 // Write trace file on root only
 if (Communicator_pt->my_rank()==0)
  {
   Vector<double> x(2,0.0);
   double u_exact;
   GlobalParameters::get_exact_u(time_pt()->time(),x,u_exact);
   Vector<double > xi_wall(1);
   Vector<double > r_wall(2);
   xi_wall[0]=0.0;
   Boundary_pt->position(xi_wall,r_wall);
   Trace_file << time_pt()->time() 
              << " " << u_exact
              << " " << r_wall[0]
              << " " << time_pt()->dt()
              << " " << temp_error_norm
              << " " << GlobalParameters::step_position(time_pt()->
                                                               time())
              << " " << global_ref_stats[0] 
              << " " << global_ref_stats[1]
              << " " << global_ref_stats[2] 
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
  }


 // Write restart file
 sprintf(filename,"%s/restart%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 some_file.precision(20); 
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
   
   //What is the face index of element e on boundary b
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


//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template<class ELEMENT>
double RefineableUnsteadyHeatProblem<ELEMENT>::global_temporal_error_norm()
{

#ifdef OOMPH_HAS_MPI

 double global_error = 0.0;
   
 //Find out how many nodes there are in the problem
 unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Loop over the nodes and calculate the estimated error in the values
 // for non-haloes
 int count=0;
 for(unsigned i=0;i<n_node;i++)
  {
   Node* nod_pt=Bulk_mesh_pt->node_pt(i);
   if (!(nod_pt->is_halo()))
    {
     // Get error in solution: Difference between predicted and actual
     // value for nodal value 0
     double error = nod_pt->time_stepper_pt()->
      temporal_error_in_value(nod_pt,0);
     
     //Add the square of the individual error to the global error
     global_error += error*error;
     count++;
    }
  }
 
 // Accumulate
 int n_node_local=count;
 int n_node_total=0;
 MPI_Allreduce(&n_node_local,&n_node_total,1,MPI_INT,MPI_SUM,
               communicator_pt()->mpi_comm());
 
 double global_error_total=0.0;
 MPI_Allreduce(&global_error,&global_error_total,1,MPI_DOUBLE,MPI_SUM,
               communicator_pt()->mpi_comm());
 
 // Divide by the number of nodes
 global_error_total /= double(n_node_total);
 
 // Return square root...
 oomph_info << "Total error " << n_node_total << " " 
            <<  sqrt(global_error_total) << std::endl;
 return sqrt(global_error_total);
 
#else
 
 double global_error = 0.0;
 
 //Find out how many nodes there are in the problem
 unsigned n_node = Bulk_mesh_pt->nnode();

 //Loop over the nodes and calculate the errors in the values
 for(unsigned i=0;i<n_node;i++)
  {
   // Get error in solution: Difference between predicted and actual
   // value for nodal value 0
   double error = 
    Bulk_mesh_pt->node_pt(i)->time_stepper_pt()->
    temporal_error_in_value(Bulk_mesh_pt->node_pt(i),0);
   
   //Add the square of the individual error to the global error
   global_error += error*error;
  }
    
   //Now the global error must be divided by the number of nodes
 global_error /= double(n_node);

 //Return the square root of the error
 return sqrt(global_error);

#endif


} // end of global_temporal_error_norm


//=======start_of_dump_it=================================================
/// Dump the solution to disk
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::dump_it(ofstream& dump_file)
{

 // Write step number
 dump_file << Doc_info.number() << " # step number" << std::endl;

 // Write suggested next timestep
 dump_file << Next_dt << " # suggested next timestep" << std::endl;

 // Dump the refinement pattern and the generic problem data
 Problem::dump(dump_file);
  
} // end of dump_it

//=========start_of_restart===============================================
/// Read solution from disk
//========================================================================
template<class ELEMENT>
void RefineableUnsteadyHeatProblem<ELEMENT>::restart(ifstream& restart_file)
{

 double local_next_dt=0.0;
 unsigned local_doc_info_number=0;
 
 if (restart_file.is_open())
  {
   oomph_info <<"restart file exists"<<endl;
     
   // Read line up to termination sign
   string input_string;
   getline(restart_file,input_string,'#');
   
   // Ignore rest of line
   restart_file.ignore(80,'\n');
   
   // Read in step number
   local_doc_info_number=unsigned(atof(input_string.c_str()));
   
   // Read line up to termination sign
   getline(restart_file,input_string,'#');
   
   // Ignore rest of line
   restart_file.ignore(80,'\n');
   
   // Read suggested next timestep
   local_next_dt=double(atof(input_string.c_str()));
   
  }
 else
  {
   oomph_info <<"restart file does not exist"<<endl;
  }

 // Get next dt
 double next_dt=0.0;
 MPI_Allreduce(&local_next_dt,&next_dt,1,MPI_DOUBLE,MPI_MAX,
               communicator_pt()->mpi_comm());
 Next_dt=next_dt;
 oomph_info << "Next dt: " << Next_dt << std::endl;

 // Get doc info
 unsigned doc_info_number=0; 
 MPI_Allreduce(&local_doc_info_number,&doc_info_number,1,MPI_UNSIGNED,MPI_MAX,
               communicator_pt()->mpi_comm());
 Doc_info.number()=doc_info_number;
 oomph_info << "Restarted doc_info.number(): " << doc_info_number << std::endl;
 

 // Refine the mesh and read in the generic problem data
 Problem::read(restart_file);
 
 // Output restarted solution
 std::stringstream comment_stream;
 comment_stream << "Restarted; Step " 
                << doc_info().number();
 doc_solution(comment_stream.str());


} // end of restart


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main=====================================================
/// Demonstrate how to solve an unsteady heat problem in deformable domain
/// with mesh adaptation.  
//========================================================================
int main(int argc, char* argv[])
{

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Name of restart file
 CommandLineArgs::specify_command_line_flag("--restart_file",
                                            &GlobalParameters::Restart_file);
 
 // Name of file that specifies the problem partitioning
 CommandLineArgs::specify_command_line_flag(
  "--partitioning_file",&GlobalParameters::Partitioning_file); 

 // Flag to indicate that we're doing a validation run
 CommandLineArgs::specify_command_line_flag("--validation_run");

 // Flag to indicate that we do load_balance first (before prune)
 CommandLineArgs::specify_command_line_flag("--load_balance_first");
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Build problem: Pass pointer to source function and initial timestep
 RefineableUnsteadyHeatProblem<RefineableQUnsteadyHeatElement<2,3> >
  problem(&GlobalParameters::get_source);

 // Work out doc_info directory
 char doc_info_directory[100];
 if(CommandLineArgs::command_line_flag_has_been_set("--validation_run"))
  {
   //Validation run
   if(CommandLineArgs::command_line_flag_has_been_set("--load_balance_first"))
    {
     //Load balancing
     if(!CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
      {
       //First run
       sprintf(doc_info_directory,"RESLT_load_balance_first_for_restart");
      }
     else
      {
       //Restarting
       //cout << GlobalParameters::Restart_file.c_str() << endl;
       char filename[100];
       sprintf(filename,"%s",GlobalParameters::Restart_file.c_str());
       char* step;
       step = strtok(filename,"/");
       step = strtok(NULL,"/");
       sprintf(doc_info_directory,"RESLT_load_balance_first_restarted_from_step_%s",step);
       //cout << doc_info_directory << endl;
       //exit(1);
      }
    }
   else
    {
     //Pruning
     if(!CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
      {
       //First run
       sprintf(doc_info_directory,"RESLT_prune_first_for_restart");
      }
     else
      {
       //Restarting
       char filename[100];
       sprintf(filename,"%s",GlobalParameters::Restart_file.c_str());
       char* step;
       step = strtok(filename,"/");
       step = strtok(NULL,"/");
       sprintf(doc_info_directory,"RESLT_prune_first_restarted_from_step_%s",step);
      }
    }
  }
 else
  {
   sprintf(doc_info_directory,"RESLT");
  }
 
 // Set doc_info directory
 problem.doc_info().set_directory(doc_info_directory);

  // Switch off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;

 // Define processor-labeled output file for all on-screen stuff
 std::ofstream output_stream;
 char filename[100];
 sprintf(filename,"%s/OUTPUT.%i",
         problem.doc_info().directory().c_str(),
         MPI_Helpers::communicator_pt()->my_rank());
 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);  

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
   
 // First timestep?
 bool first=true;
 
 // Max. number of spatial adaptations per timestep. Allow plenty
 // of adaptations at first timestep as the initial conditions
 // can be reset "exactly" from without any interpolation error.
 unsigned max_adapt=10; 

 // Set IC from analytical soln
 problem.set_initial_condition();
 
 // If restart: The first step isn't really the first step,
 // i.e. initial condition should not be re-set when 
 // adaptive refinement has been performed. Also, limit
 // the max. number of refinements per timestep to the
 // normal value straightaway.
 if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
   first=false;
   max_adapt=1;
  }
 

 // Prune first or load balance first?
 unsigned first_load_balance=9;
 unsigned second_load_balance=18;
 unsigned first_prune=4;
 unsigned second_prune=13;

 // Create storage for partitioning
 Vector<unsigned> element_partition;
 
 // No partitioning specified -- simply distribute problem and
 // record distribution in file to allow restart
 if (!CommandLineArgs::command_line_flag_has_been_set("--partitioning_file"))
  {
   oomph_info << "Distributing problem with METIS \n"
              << "Documenting partitioning in RESLT/partitioning.dat.\n";
   
   // Distribute
   element_partition=problem.distribute();
   
   // Write partition to disk
   std::ofstream output_file;
   char filename[100];
   sprintf(filename,"%s/partitioning.dat",problem.doc_info().directory().c_str());
   output_file.open(filename);
   unsigned n=element_partition.size();
   output_file << element_partition.size() << std::endl;
   for (unsigned e=0;e<n;e++)
    {
     output_file << element_partition[e] << std::endl;
    }
   output_file.close();
  }
 // Read in partitioning from disk using the specified file
 else 
  {
   
   oomph_info << "Distributing problem with partitioning \n"
              << "specified in: " << GlobalParameters::Partitioning_file 
              << std::endl;
   
   // Read in partitioning from disk: Name of partitioning file specified
   // as command line argument
   std::ifstream input_file;
   input_file.open(GlobalParameters::Partitioning_file.c_str());
   if (!input_file.is_open())
    {
     oomph_info << "Error opening input file\n";
     assert(false);
    }
   std::string input_string;
   unsigned n=0;
   if (getline(input_file,input_string,'\n'))
    {
     n=atoi(input_string.c_str());
    }
   else
    {
     oomph_info << "Reached end of file when reading partitioning file\n";
     assert(false);
    }
   element_partition.resize(n);
   for (unsigned e=0;e<n;e++)
    {
     if (getline(input_file,input_string,'\n'))
      {
       element_partition[e]=atoi(input_string.c_str());
      }
     else
      {
       oomph_info << "Reached end of file when reading partitioning file\n";
       assert(false);
      }
    }
   
   // Create storage for actually used partitioning
   Vector<unsigned> used_element_partition;
 
   // Now perform the distribution 
   used_element_partition=problem.distribute(element_partition);


   // Write used partition to disk
   std::ofstream output_file;
   char filename[100];
   sprintf(filename,"%s/partitioning.dat",problem.doc_info().directory().c_str());
   output_file.open(filename);
   unsigned n_used=used_element_partition.size();
   output_file << n_used << std::endl;
   for (unsigned e=0;e<n_used;e++)
    {
     output_file << used_element_partition[e] << std::endl;
    }
   output_file.close();

  }


 // Doc
 std::stringstream comment_stream0;
 comment_stream0 << "Distributed; Step " 
                 << problem.doc_info().number();
 problem.doc_solution(comment_stream0.str());
 oomph_info << "Finished distribution\n";
 
// Restart file specified via command line 
 if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
   problem.restart();
   oomph_info << "Finished restart\n";
  }

 // Initial timestep: Use the one used when setting up the initial
 // condition or the "suggested next dt" from the restarted run
 double dt=problem.next_dt();
 
 // Write header for trace file
 problem.write_trace_file_header();
 
 // Timestepping loop
 if (CommandLineArgs::command_line_flag_has_been_set("--validation_run"))
  {
   // Fake but repeatable load balancing for self-test
   problem.set_default_partition_in_load_balance();
  }


 // Keep solving...
 while (problem.doc_info().number()<22)  
  {
   
   oomph_info << "Doing solve\n";
   
   // Take timestep with temporal and spatial adaptivity
   double dt_new=
    problem.doubly_adaptive_unsteady_newton_solve(dt,problem.epsilon_t(),
                                                  max_adapt,first);
   oomph_info << "Suggested new dt: " << dt_new << std::endl;
   dt=dt_new;
   
   // Store for restart
   problem.next_dt()=dt;
   
   // Now we've done the first timestep -- don't re-set the IC
   // in subsequent steps
   first=false;
   
   // Reduce the number of spatial adaptations to one per 
   // timestep -- too scared that the interpolation error will 
   // wipe out any further gains...
   max_adapt=1;
   
   //Output solution
   std::stringstream comment_stream;
   comment_stream << "Step " << problem.doc_info().number();
   problem.doc_solution(comment_stream.str());
   
   // Prune
   if ((problem.doc_info().number()==first_prune)||
       (problem.doc_info().number()==second_prune))
    {
     std::stringstream comment_stream1;
     if ((problem.doc_info().number()==second_prune)||
         ((problem.doc_info().number()==first_prune)&&
          (!CommandLineArgs::command_line_flag_has_been_set(
            "--load_balance_first"))))
      {
       oomph_info << "Refining uniformly before pruning\n";
       problem.refine_uniformly();       
       comment_stream1 << "Uniform refinement; Step " 
                       << problem.doc_info().number();
      }
     else
      {
       comment_stream1 << "Skipped uniform refinement; Step " 
                       << problem.doc_info().number();
      }
     problem.doc_solution(comment_stream1.str());

     std::stringstream comment_stream2;
     if ((problem.doc_info().number()==second_prune+1)||
         ((problem.doc_info().number()==first_prune+1)&&
          (!CommandLineArgs::command_line_flag_has_been_set(
            "--load_balance_first"))))
      {
       oomph_info << "Pruning\n";
       problem.prune_halo_elements_and_nodes();
       comment_stream2 << "Pruned; Step " 
                       << problem.doc_info().number();
      }
     else
      {
       comment_stream2 << "Skipped prune; Step " 
                       << problem.doc_info().number();
      }
     problem.doc_solution(comment_stream2.str());
    }
   
   
   
   // Load balance?
   if ((problem.doc_info().number()==first_load_balance)||
       (problem.doc_info().number()==second_load_balance))
    {
     oomph_info << "\n\n\n LOAD BALANCING \n\n\n";
     
     // Load balance
     problem.load_balance(); 
     
     std::stringstream comment_stream3;
     comment_stream3 << "Load balanced; Step " 
                     << problem.doc_info().number();
     problem.doc_solution(comment_stream3.str());
     
     
     oomph_info << "\n\n\n DONE LOAD BALANCING \n\n\n";
    }
  }

 oomph_info << "Done\n";
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

}; // end of main

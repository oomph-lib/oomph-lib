//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Driver for Axisymmetic spin up in a cylinder using black
// box adaptation, using Axisymmetric Taylor Hood and Crouzeix 
// Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "navier_stokes.h"

// Axisym Navier Stokes headers
#include "axisym_navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re = 5.0;

 /// Womersley number
 double ReSt = 5.0;

} // End of namespace



//==start_of_problem_class===============================================
/// \short Refineable rotating cylinder problem in rectangular 
/// axisymmetric domain.
//=======================================================================
template<class ELEMENT, class TIMESTEPPER>
class RefineableRotatingCylinderProblem : public Problem
{

public:

 /// Constructor
 RefineableRotatingCylinderProblem(const unsigned& n_r,
                                   const unsigned& n_z,
                                   const double& l_r,
                                   const double& l_z);

 /// Destructor (empty)
 ~RefineableRotatingCylinderProblem() {}

 /// Update before solve (empty)
 void actions_before_newton_solve() {}
  
 /// Update after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Actions before timestep: (Re-)set velocity boundary conditions
 void actions_before_implicit_timestep()
  {
   // Determine number of mesh boundaries
   const unsigned n_boundary = mesh_pt()->nboundary();
   
   // Loop over mesh boundaries
   for(unsigned b=0;b<n_boundary;b++)
    {
     // Determine number of nodes on boundary b
     const unsigned n_node = mesh_pt()->nboundary_node(b);
     
     // Loop over nodes on boundary b
     for(unsigned n=0;n<n_node;n++)
      {
       // Loop over the three velocity components
       for(unsigned i=0;i<3;i++)
	{
         // For the solid boundaries (boundaries 0,1,2)
         if(b<3)
          {
           // Get the radial value
           double r = mesh_pt()->boundary_node_pt(b,n)->x(0);
           
           // Set all velocity components to no flow along boundary
           switch(i)
            {
            case 2: // Azimuthal velocity
             mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,r);
             break;
            case 1: // Axial velocity
             mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
             break;
            case 0: // Radial velocity
             mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
             break;
            }
          }
         
         // For the symmetry boundary (boundary 3)
         if(b==3)
          {
           // Set only the radial (i=0) and azimuthal (i=2) velocity
           // components to no flow along boundary
           // (axial component is unconstrained)
           if(i!=1)
            {
             mesh_pt()->boundary_node_pt(b,n)->set_value(0,i,0.0);
            }
          }
        } // End of loop over velocity components
      } // End of loop over nodes on boundary b
    } // End of loop over mesh boundaries
  } // End of actions_before_implicit_timestep
 
 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);

  } // End of actions_after_adapt

 /// Set initial conditions: Set all nodal velocities to zero and
 /// initialise the previous velocities to correspond to an impulsive
 /// start
 void set_initial_condition()
  {
   // Determine number of nodes in mesh
   const unsigned n_node = mesh_pt()->nnode();

   // Loop over all nodes in mesh
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       // Set velocity component i of node n to zero
       mesh_pt()->node_pt(n)->set_value(i,0.0);
      }
    }

   // Initialise the previous velocity values for timestepping
   // corresponding to an impulsive start
   assign_initial_values_impulsive();

  } // End of set_initial_condition
 
 /// \short Access function for the specific mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double& t_max, const double& dt, string dir_name, 
                   const unsigned& maximum_ref_level,
                   const unsigned& minimum_ref_level);

private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned& e, const unsigned& pdof, 
                   const double& pvalue)
  {
   // Fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }

 /// Trace file
 ofstream Trace_file;

}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for refineable rotating cylinder problem
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::
RefineableRotatingCylinderProblem(const unsigned& n_r,
                                  const unsigned& n_z,
                                  const double& l_r,
                                  const double& l_z)
{
 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER);

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_r,n_z,l_r,l_z,
                                             time_stepper_pt());

 // Set error estimator for spatial adaptivity
 Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;

 // Pass error estimator to the mesh
 mesh_pt()->spatial_error_estimator_pt() = error_estimator_pt;
 
 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here

 // Determine number of mesh boundaries
 const unsigned n_boundary = mesh_pt()->nboundary();

 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = mesh_pt()->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Pin values for radial velocity on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(0);

     // Pin values for axial velocity on all SOLID boundaries
     // (i.e. b = 0,1,2)
     if(b!=3) { mesh_pt()->boundary_node_pt(b,n)->pin(1); }

     // Pin values for azimuthal velocity on all boundaries
     mesh_pt()->boundary_node_pt(b,n)->pin(2);

    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------

 // Determine number of elements in mesh
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Assign the time pointer
   el_pt->time_pt() = time_pt();

   // The mesh remains fixed
   el_pt->disable_ALE();

  } // End of loop over elements
 
 // Pin redundant pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the pressure in first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // End of constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo& doc_info)
{ 
 // Document in trace file
 Trace_file << time_pt()->time() << " " 
            << mesh_pt()->max_permitted_error() << " "
            << mesh_pt()->min_permitted_error() << " "
            << mesh_pt()->max_error() << " "
            << mesh_pt()->min_error() << " " << std::endl;

 ofstream some_file;
 char filename[100];

 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;

 // Open solution output file
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Output solution to file
 mesh_pt()->output(some_file,npts);

 // Write file as a tecplot text object...
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 // ...and draw a horizontal line whose length is proportional
 // to the elapsed time
 some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 some_file << "1" << std::endl;
 some_file << "2" << std::endl;
 some_file << " 0 0" << std::endl;
 some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "-0.05 -0.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "1.05 -0.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "-0.05 1.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "1.05 1.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "-0.05 -0.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "1.05 -0.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "-0.05 1.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "1.05 1.05 0.007 0.007 1.05 2.55" << std::endl;

 // Close solution output file
 some_file.close();
 
} // End of doc_solution



//==start_of_unsteady_run=================================================
/// Perform run up to specified time t_max with given timestep dt
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::unsteady_run(
 const double& t_max, const double& dt, string dir_name, 
 const unsigned& maximum_ref_level, const unsigned& minimum_ref_level)
{
 // Create DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory(dir_name); 

 // Set step number to zero
 doc_info.number()=0;

 // Open trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise trace file
 Trace_file << "time" << ", " 
            << "max permitted error" << ", "
            << "min permitted error" << ", "
            << "max error" << ", "
            << "min error" << ", " << std::endl;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial conditions
 set_initial_condition();

 // Override the maximum and minimum permitted errors
 mesh_pt()->max_permitted_error() = 1.0e-2; // Default = 1.0e-3
 mesh_pt()->min_permitted_error() = 1.0e-3; // Default = 1.0e-5
 
 // Override the maximum and minimum permitted refinement levels
 mesh_pt()->max_refinement_level() = maximum_ref_level;
 mesh_pt()->min_refinement_level() = minimum_ref_level;

 // Maximum number of spatial adaptations per timestep
 unsigned max_adapt = maximum_ref_level - minimum_ref_level;

 // Initial refinement level
 for(unsigned i=0;i<minimum_ref_level;i++) { refine_uniformly(); }
 
 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);
 
 // Doc initial solution
 doc_solution(doc_info);

 // Increment counter for solutions 
 doc_info.number()++;

 // Are we on the first timestep? At this point, yes!
 bool first_timestep = true;
 
 // Specify normalising factor explicitly
 Z2ErrorEstimator* error_pt=dynamic_cast<Z2ErrorEstimator*>(
  mesh_pt()->spatial_error_estimator_pt());
 error_pt->reference_flux_norm() = 0.01;
 
 // Timestepping loop
 for(unsigned t=1;t<=n_timestep;t++)
  {
   // Output current timestep to screen
   cout << "\nTimestep " << t << " of " << n_timestep << std::endl;

   // Take fixed timestep with spatial adaptivity
   unsteady_newton_solve(dt,max_adapt,first_timestep);
   
   // No longer on first timestep, so set first_timestep flag to false
   first_timestep = false; 

   // Reset maximum number of adaptations for all future timesteps
   max_adapt = 1;
   
   // Doc solution
   doc_solution(doc_info);
   
   // Increment counter for solutions 
   doc_info.number()++;
   
  } // End of timestepping loop
 
} // End of unsteady_run


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for RotatingCylinderProblem
//=====================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 /// Maximum time
 double t_max = 1.0;

 /// Duration of timestep
 const double dt = 0.01;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::Argc>1) { t_max = 0.02; }
 
 // Number of elements in radial (r) direction
 const unsigned n_r = 2;

 // Number of elements in axial (z) direction
 const unsigned n_z = 2;

 // Length in radial (r) direction
 const double l_r = 1.0;

 // Length in axial (z) direction
 const double l_z = 1.3;

 /// Maximum refinement level
 const unsigned maximum_ref_level = 4;

 /// Minimum refinement level
 const unsigned minimum_ref_level = 1;

 // Doing RefineableAxisymmetricQTaylorHoodElements
 {
  // Build the problem with RefineableAxisymmetricQTaylorHoodElements
  RefineableRotatingCylinderProblem<
   RefineableAxisymmetricQTaylorHoodElement, BDF<2> > 
   problem(n_r,n_z,l_r,l_z);

  cout << "Doing RefineableAxisymmetricQTaylorHoodElement" << std::endl;

  // Solve the problem and output the solution
  problem.unsteady_run(t_max,dt,"RESLT_TH",
                       maximum_ref_level,minimum_ref_level);

 } // End of RefineableAxisymmetricQTaylorHoodElements
  
 // Doing RefineableAxisymmetricQCrouzeixRaviartElements
 {
  // Build the problem with RefineableAxisymmetricQCrouzeixRaviartElements
  RefineableRotatingCylinderProblem<
   RefineableAxisymmetricQCrouzeixRaviartElement, BDF<2> > 
   problem(n_r,n_z,l_r,l_z);

  cout << "Doing RefineableAxisymmetricQCrouzeixRaviartElement" << std::endl;
  
  // Solve the problem and output the solution
  problem.unsteady_run(t_max,dt,"RESLT_CR",
                       maximum_ref_level,minimum_ref_level);

 } // End of RefineableAxisymmetricQCrouzeixRaviartElements
 
} // End of main







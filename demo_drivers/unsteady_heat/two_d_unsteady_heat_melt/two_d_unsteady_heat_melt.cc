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
//Driver for 2D unsteady heat problem with "melting" constraint
#include <fenv.h> 

//Generic routines
#include "generic.h"

// The unsteady heat equations
#include "unsteady_heat.h"

// Mesh
#include "meshes/rectangular_quadmesh.h"

// The unsteady heat equations
#include "unsteady_heat_flux_melt_elements.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;




//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{

 /// Get flux applied along boundary x=0.
 void flux(const double& time, const Vector<double>& x, double& flux)
 {
  flux = 20.0*sin(2.0*4.0*MathematicalConstants::Pi*time)*
   x[1]*(1.0-x[1]);
 }

 /// Melt-temperature
 double Melt_temperature=0.0;


} // end of ProblemParameters

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// UnsteadyHeat problem 
//====================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public Problem
{

public:

  /// Constructor
 UnsteadyHeatProblem();

 /// Destructor (empty)
 ~UnsteadyHeatProblem(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Bulk mesh
 RectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem()
{ 
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Setup mesh
 //-----------

 // Number of elements in x and y directions
 unsigned nx=5;
 unsigned ny=5;
 
 // Lengths in x and y directions
 double lx=1.0;
 double ly=1.0;

 // Build mesh
 Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt());

 // Flux boundary
 unsigned b_flux=3;

 // Create the surface mesh as an empty mesh
 Surface_mesh_pt=new Mesh;

 // How many bulk elements are adjacent to boundary b?
 unsigned b=b_flux; 
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the face index of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   UnsteadyHeatFluxMeltElement<ELEMENT>* flux_element_pt = new 
    UnsteadyHeatFluxMeltElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   Surface_mesh_pt->add_element_pt(flux_element_pt);
   
  } //end of loop over bulk elements adjacent to boundary b
 

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 
 // Combine all submeshes into a single global Mesh
 build_global_mesh();
 
 
// hierher
 // Set the initial conditions
 // --------------------------
 unsigned nnod = Bulk_mesh_pt->nnode();
 for(unsigned j=0;j<nnod;j++)
  {
   Bulk_mesh_pt->node_pt(j)->set_value(0,-0.1); 
  } 
 

 // Set the boundary conditions for this problem: 
 // ---------------------------------------------

 // All nodes are free by default -- just pin the ones that have 
 // Dirichlet conditions here. 
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   if (b!=b_flux) 
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    }
  } // end of set boundary conditions



 // Complete the build of all elements so they are fully functional
 //----------------------------------------------------------------

 // Find number of elements in mesh
 n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   // Set pointer to continous time
   el_pt->time_pt()=time_pt();
  }

 // Loop over the flux elements to pass pointer to prescribed flux function
 n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to UnsteadyHeat flux element
   UnsteadyHeatFluxMeltElement<ELEMENT> *el_pt = 
    dynamic_cast<UnsteadyHeatFluxMeltElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = &ProblemParameters::flux;

   // Set pointer to continous time
   el_pt->time_pt()=time_pt();

   // Set melt temperature
   el_pt->melt_temperature_pt()=&ProblemParameters::Melt_temperature;

   // Suppress melting?
   if (CommandLineArgs::command_line_flag_has_been_set("--disable_melting"))
    {
     el_pt->disable_melting();
    }
   else
    {
     // Disable suppression of refreezing?
     if (CommandLineArgs::command_line_flag_has_been_set(
          "--disable_suppression_of_refreezing"))
      {
       el_pt->disable_suppression_of_refreezing();
      }
    }
  }


 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output flux
 sprintf(filename,"%s/flux%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

} // end of doc_solution



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// \short Driver code for unsteady heat equation
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress melting?
 CommandLineArgs::specify_command_line_flag("--disable_melting");

 // Disable suppression of refreezing?
 CommandLineArgs::specify_command_line_flag(
  "--disable_suppression_of_refreezing");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Build problem
 UnsteadyHeatProblem<QUnsteadyHeatElement<2,3> > problem;

 // Setup labels for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");
 
 // Output number
 doc_info.number()=0;
 
 // Choose simulation interval and timestep
 double t_max=0.5;
 double dt=0.01;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);
 
 // Set impulsive IC
 problem.assign_initial_values_impulsive(dt);
 
 //Output initial condition
 problem.doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Find number of steps
 unsigned nstep = unsigned(t_max/dt);

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   cout << " Timestep " << istep << std::endl;
   
   // Take timestep
   problem.unsteady_newton_solve(dt);
   
   //Output solution
   problem.doc_solution(doc_info);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }
 

}; // end of main

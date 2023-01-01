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
//Driver for axisymmetric FvK equations
#include <fenv.h> 

// Generic oomph-lib routines
#include "generic.h"

// Include axisymmetric fvk elements/equations
#include "axisym_foeppl_von_karman.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;


//==================================================================
/// Namespace for pressure solution
//==================================================================
namespace AxisymFvKParameters
{
 /// Constant pressure
 double Pressure=0.0;
 
 /// Function that computes the pressure at radius r
 void pressure_function(const double& r,double& pressure ) 
 {
  pressure = Pressure;
 }
 
 /// FvK parameter
 double Eta = 2.39e6;
 
 /// Function to get the exact solution for the pure bending model
 void get_exact_u(const Vector<double>& r, Vector<double>& u )
 {
  u[0]=AxisymFvKParameters::Pressure*(r[0]*r[0]-1.0)*(r[0]*r[0]-1.0)/64.0;
 }
 
 /// Directory
 string Directory = "RESLT";
 }


//==start_of_problem_class============================================
/// 
//====================================================================
template<class ELEMENT> 
class AxisymFvKProblem : public Problem
{
public:
 
 /// Constructor: Pass number of elements
 AxisymFvKProblem(const unsigned& n_element);
 
 /// Destructor
 ~AxisymFvKProblem()
  {
   delete mesh_pt();
  }
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution.
 void doc_solution();
 
private:
 
 /// Doc info object for labeling output
 DocInfo Doc_info;

}; // end of problem class



//=====start_of_constructor===============================================
/// Constructor
//========================================================================
template<class ELEMENT>
AxisymFvKProblem<ELEMENT>::
AxisymFvKProblem(const unsigned& n_element) 
 
{ 
 // Set domain length 
 double L=1.0;
 
 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);
 
 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 
 
 // Do-nothing on laplacian_w in order to pin dwdr at r=0
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at r=1)
 mesh_pt()->boundary_node_pt(1,0)->pin(0);
 mesh_pt()->boundary_node_pt(1,0)->pin(2);
 
 // Complete the setup
 // Loop over elements and set pointers to pressure function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the pressure function pointer
   elem_pt->pressure_fct_pt() =  &AxisymFvKParameters::pressure_function;
   
   //Set the pointer to eta
   elem_pt->eta_pt() = &AxisymFvKParameters::Eta;
   
   // Choose between pure bending or fvk model
   if (CommandLineArgs::command_line_flag_has_been_set("--use_linear_eqns"))
    {
     elem_pt->use_linear_bending_model();
    }
  }
 
 // Setup equation numbering scheme
 assign_eqn_numbers();
 
} // end of constructor



//===start_of_doc=========================================================
/// Doc the solution in tecplot format.
//========================================================================
template<class ELEMENT>
void AxisymFvKProblem<ELEMENT>::doc_solution()
{ 
 //  Number of plot points
 unsigned npts;
 npts=5; 
 
 // Output solution with specified number of plot points per element
 char filename[100];
 sprintf(filename, "%s/sol_%i.dat",
         AxisymFvKParameters::Directory.c_str(),Doc_info.number());
 ofstream solution_file(filename,ios::app);
 mesh_pt()->output(solution_file,npts);
 solution_file.close();
 
 
 //  Output exact solution
 sprintf(filename, "%s/exact_sol_%i.dat",
         AxisymFvKParameters::Directory.c_str(),Doc_info.number());
 ofstream exact_file(filename,ios::app); 
 mesh_pt()->output_fct(exact_file,npts,AxisymFvKParameters::get_exact_u); 
 exact_file.close();
 
 
 // Output solution at the centre (r=0) as function of the pressure
 sprintf(filename, "%s/w_centre.dat",AxisymFvKParameters::Directory.c_str());
 ofstream w_centre_file(filename,ios::app); 
 w_centre_file << AxisymFvKParameters::Pressure << " " 
               << mesh_pt()->node_pt(0)->value(0) << std::endl;
 w_centre_file.close();


 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver 
//=====================================================================
int main(int argc, char **argv)
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Value of eta
 CommandLineArgs::specify_command_line_flag("--eta",
                                            &AxisymFvKParameters::Eta);
 
 // Choose between pure bending or fvk model
 CommandLineArgs::specify_command_line_flag("--use_linear_eqns");
 
 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &AxisymFvKParameters::Directory);
 
 // Number of elements in the mesh
 unsigned n_element=100;
 CommandLineArgs::specify_command_line_flag("--n_element", &n_element);
 
 // Increment displacements
 unsigned n_step=10; // hierher 10;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);
 double dp=10;
 CommandLineArgs::specify_command_line_flag("--dp", &dp);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Set up the problem: 
 AxisymFvKProblem<AxisymFoepplvonKarmanElement<3> >  
  problem(n_element);
 
 // Set initial value for pressure 
 AxisymFvKParameters::Pressure=0.0;
 
 for (unsigned i=0;i<n_step;i++)
  {
   // Solve the problem
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution();
   
   // Increment pressure
   AxisymFvKParameters::Pressure+=dp; 
  }

 cout << "\n Done \n" << std::endl;
 
} // end of main

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
// Time-harmonic deformation of an elastic annulus subject to 
// a nonuniform pressure load.

#include<fenv.h>

#include "math.h"
#include <complex>

//Oomph-lib includes
#include "generic.h"
#include "pml_time_harmonic_linear_elasticity.h"
#include "oomph_crbond_bessel.h"

//The meshes
#include "meshes/triangle_mesh.h"

// The meshes needed in the PML constructions
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 ///  Body force
 void body_force(const Vector<double>& x,
                 Vector<std::complex<double> >& b)
 {
  /// Define the magnitude of the forcing
  double magnitude=1.0;

  /// Define the decay rate
  double alpha=1.0;

  /// The forcing used is of damped exponential type
  /// to generate axisymmetric results after post-processing
  double source_fct=magnitude*exp(-alpha*sqrt(x[0]*x[0]+x[1]*x[1]));
 
  /// Assign the source in each dimension
  b[0]=source_fct*x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
  b[1]=source_fct*x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);
 }

 /// helpers to time the code
 double T_start = 0.0;
 double T_end = 0.0;

 /// PML width in elements for the right layer
 unsigned N_pml_multiplier = 1;
 double L_pml_multiplier = 1.0;

 /// PML width in elements for the right layer
 unsigned N_x_right_pml = 8;

 /// PML width in elements for the top layer
 unsigned N_y_top_pml = 8;

 /// PML width in elements for the left layer
 unsigned N_x_left_pml = 8; 

 /// PML width in elements for the left layer
 unsigned N_y_bottom_pml = 8; 

 // Outer physical length of the PML layers
 // defaults to 0.2, so 10% of the size of the 
 // physical domain
 double Width_x_right_pml  = 2.0;
 double Width_y_top_pml    = 2.0;
 double Width_x_left_pml   = 2.0;
 double Width_y_bottom_pml = 2.0;

 /// Function to compute dependent parameters
 void compute_dependent_parameters()
 {
  /// Adjust number of PML elements, set to be equal for all layers
  N_x_right_pml  = N_x_right_pml  * N_pml_multiplier;
  N_x_left_pml   = N_x_left_pml   * N_pml_multiplier;
  N_y_top_pml    = N_y_top_pml    * N_pml_multiplier;
  N_y_bottom_pml = N_y_bottom_pml * N_pml_multiplier;

  ///Adjust physical size of PML layers, set to be equal for all layers
  Width_x_right_pml  = Width_x_right_pml  * L_pml_multiplier;
  Width_x_left_pml   = Width_x_left_pml   * L_pml_multiplier;
  Width_y_top_pml    = Width_y_top_pml    * L_pml_multiplier;
  Width_y_bottom_pml = Width_y_bottom_pml * L_pml_multiplier;
 }

 /// Poisson's ratio
 double Nu=0.3;

 /// Square of non-dim frequency 
 double Omega_sq=10.0; 
  
 /// The elasticity tensor
 PMLTimeHarmonicIsotropicElasticityTensor* E_pt;

 /// Output directory
 string Directory="RESLT";
} //end_namespace



//=============begin_problem============================================ 
/// Annular disk
//====================================================================== 
template<class ELASTICITY_ELEMENT>
class ElasticAnnulusProblem : public Problem
{

public:

 /// Constructor:
 ElasticAnnulusProblem();

 /// Destructor (empty)
 ~ElasticAnnulusProblem(){}
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Create PML meshes
 void create_pml_meshes();

 /// Actions before adapt: Wipe the mesh of traction elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of traction elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Helper function to complete problem setup
 void complete_problem_setup();

#ifdef ADAPTIVE

 /// Pointer to refineable solid mesh
 /// Adaptivity is not verified in this particular case
 /// The computation will run as a non-adaptive solution
 RectangularQuadMesh<ELASTICITY_ELEMENT>* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 RectangularQuadMesh<ELASTICITY_ELEMENT>* Solid_mesh_pt;

#endif

 /// Pointer to the right PML mesh
 Mesh* PML_right_mesh_pt;

 /// Pointer to the top PML mesh
 Mesh* PML_top_mesh_pt;

 /// Pointer to the left PML mesh
 Mesh* PML_left_mesh_pt;

 /// Pointer to the bottom PML mesh
 Mesh* PML_bottom_mesh_pt;

 /// Pointer to the top right corner PML mesh
 Mesh* PML_top_right_mesh_pt;

 /// Pointer to the top left corner PML mesh
 Mesh* PML_top_left_mesh_pt;

 /// Pointer to the bottom right corner PML mesh
 Mesh* PML_bottom_right_mesh_pt;

 /// Pointer to the bottom left corner PML mesh
 Mesh* PML_bottom_left_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;

 /// Boundary ID of upper inner boundary
 unsigned Upper_inner_boundary_id;
 
 /// Boundary ID of upper outer boundary
 unsigned Upper_outer_boundary_id;

 /// Boundary ID of lower inner boundary
 unsigned Lower_inner_boundary_id;
 
 /// Boundary ID of lower outer boundary
 unsigned Lower_outer_boundary_id;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELASTICITY_ELEMENT>
ElasticAnnulusProblem<ELASTICITY_ELEMENT>::ElasticAnnulusProblem() 
{
 
 // Solid mesh
 //-----------

 /// Number of elements in x direction in the bulk mesh
 unsigned n_bulk_x = 20;

 /// Number of elements in y direction in the bulk mesh
 unsigned n_bulk_y = 20;
 
 /// Start and end spatial coordinates of the bulk mesh
 /// in x direction
 double x_bulk_start = -2.0;
 double x_bulk_end   =  2.0;

 /// Start and end spatial coordinates of the bulk mesh
 /// in y direction
 double y_bulk_start = -2.0;
 double y_bulk_end   =  2.0;

#ifdef ADAPTIVE

 // Build the mesh
 Solid_mesh_pt = new RectangularQuadMesh<ELASTICITY_ELEMENT>(n_bulk_x,n_bulk_y,
                                 x_bulk_start,x_bulk_end,
                                 y_bulk_start,y_bulk_end);

#else

 // Build the mesh
 Solid_mesh_pt = new RectangularQuadMesh<ELASTICITY_ELEMENT>(n_bulk_x,n_bulk_y,
                                 x_bulk_start,x_bulk_end,
                                 y_bulk_start,y_bulk_end);

#endif

 Solid_mesh_pt->output("bulk_mesh.dat");
 Solid_mesh_pt->output_boundaries("bulk_mesh_boundary.dat");

 // Create the main triangular mesh
 add_sub_mesh(Solid_mesh_pt);

 // Create PML meshes and add them to the global mesh
 create_pml_meshes();

 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Let's have a look where the boundaries are
 this->mesh_pt()->output("global_mesh.dat");
 this->mesh_pt()->output_boundaries("global_mesh_boundary.dat");
 
 // Complete problem setup
 complete_problem_setup();

 //Assign equation numbers
 cout << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory(Global_Parameters::Directory);

} //end_of_constructor




//=====================start_of_complete_problem_setup====================
/// Complete problem setup
//========================================================================
template<class ELASTICITY_ELEMENT>
void ElasticAnnulusProblem<ELASTICITY_ELEMENT>::complete_problem_setup()
{

#ifdef ADAPTIVE

 // Min element size allowed during adaptation
 if (!CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {   
   Solid_mesh_pt->min_element_size()=1.0e-5; 
  }

#endif

 //Assign the physical properties to the elements
 //----------------------------------------------
 unsigned nel=this->mesh_pt()->nelement();
 for (unsigned e=0;e<nel;e++)
  {     
   ///  Upcast from GeneralisedElement to time harmonic 
   /// linear elasticity bulk element
   PMLTimeHarmonicLinearElasticityEquations<2> *el_pt = 
    dynamic_cast<PMLTimeHarmonicLinearElasticityEquations<2>*>
    (mesh_pt()->element_pt(e));
   
   // Set the constitutive law
   el_pt->elasticity_tensor_pt() = Global_Parameters::E_pt;
   
   // Square of non-dim frequency
   el_pt->omega_sq_pt() = &Global_Parameters::Omega_sq;

   // Set the body force
   el_pt->body_force_fct_pt() = &Global_Parameters::body_force;
  }                        

}



//=====================start_of_actions_before_adapt======================
/// Actions before adapt: empty
//========================================================================
template<class ELASTICITY_ELEMENT>
void ElasticAnnulusProblem<ELASTICITY_ELEMENT>::actions_before_adapt()
{
 // Before adapting the added PML meshes must be removed
 // as they are not refineable and are to be rebuilt from the
 // newly refined triangular mesh
 delete PML_right_mesh_pt;
 PML_right_mesh_pt=0;
 delete PML_top_mesh_pt;
 PML_top_mesh_pt=0;
 delete PML_left_mesh_pt;
 PML_left_mesh_pt=0;
 delete PML_bottom_mesh_pt;
 PML_bottom_mesh_pt=0;
 delete PML_top_right_mesh_pt;
 PML_top_right_mesh_pt=0;
 delete PML_top_left_mesh_pt;
 PML_top_left_mesh_pt=0;
 delete PML_bottom_right_mesh_pt;
 PML_bottom_right_mesh_pt=0;
 delete PML_bottom_left_mesh_pt;
 PML_bottom_left_mesh_pt=0;

 // Rebuild the Problem's global mesh from its various sub-meshes
 // but first flush all its submeshes
 flush_sub_meshes();
 
 // Then add the triangular mesh back
 add_sub_mesh(Solid_mesh_pt);

 // Rebuild the global mesh
 rebuild_global_mesh();
}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: complete problem setup
//========================================================================
template<class ELASTICITY_ELEMENT>
void ElasticAnnulusProblem<ELASTICITY_ELEMENT>::actions_after_adapt()
{
 
 // Build PML meshes  and add them to the global mesh
 create_pml_meshes();

 // Build the entire mesh from its submeshes
 rebuild_global_mesh();

 // Complete problem setup
 complete_problem_setup();   
 
}// end of actions_after_adapt

//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELASTICITY_ELEMENT>
void ElasticAnnulusProblem<ELASTICITY_ELEMENT>::doc_solution(DocInfo& doc_info)
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot=5; 

 // Output displacement field
 //--------------------------
 sprintf(filename,"%s/soln_bulk%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_right%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_right_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_top%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_top_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_left%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_left_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_bottom%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_bottom_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_top_right%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_top_right_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_top_left%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_top_left_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_bottom_right%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_bottom_right_mesh_pt->output(some_file,n_plot);
 some_file.close();

 sprintf(filename,"%s/soln_pml_bottom_left%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 PML_bottom_left_mesh_pt->output(some_file,n_plot);
 some_file.close();

 // Output runtime (wall clock time) in s in a file
 sprintf(filename,"%s/wall_clock_time%i.dat",Doc_info.directory().c_str(),
         Doc_info.number()); 
 some_file.open(filename);
 some_file << Global_Parameters::T_end-Global_Parameters::T_start << std::endl;
 some_file.close();

 // Output number of degrees of freedom in a file
 sprintf(filename,"%s/ndof%i.dat",Doc_info.directory().c_str(),
         Doc_info.number()); 
 some_file.open(filename);
 some_file << this->ndof() << std::endl;
 some_file.close();

 // Output norm of solution 
 sprintf(filename,"%s/elast_soln_norm%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());   
 some_file.open(filename);   
 double norm_soln=0.0;
 this->mesh_pt()->compute_norm(norm_soln);  
 some_file << sqrt(norm_soln) << std::endl;
 cout << "Norm of computed solution: "   << sqrt(norm_soln)  << endl;

 // Increment label for output files
 Doc_info.number()++;

} //end doc

//============start_of_create_pml_meshes======================================
/// Create PML meshes and add them to the problem's sub-meshes
//============================================================================
template<class ELASTICITY_ELEMENT>
void ElasticAnnulusProblem<ELASTICITY_ELEMENT>::create_pml_meshes()
{

 // Bulk mesh left boundary id
 unsigned int left_boundary_id = 3;

 // Bulk mesh top boundary id
 unsigned int top_boundary_id = 2;

 // Bulk mesh right boundary id
 unsigned int right_boundary_id = 1;
 
 // Bulk mesh bottom boundary id
 unsigned int bottom_boundary_id = 0;
 
 // Build the PML meshes based on the new adapted triangular mesh
 PML_right_mesh_pt = 
  TwoDimensionalPMLHelper::create_right_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (Solid_mesh_pt, right_boundary_id, 
  Global_Parameters::N_x_right_pml, 
  Global_Parameters::Width_x_right_pml);

 PML_top_mesh_pt   = 
  TwoDimensionalPMLHelper::create_top_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (Solid_mesh_pt, top_boundary_id, 
   Global_Parameters::N_y_top_pml, 
   Global_Parameters::Width_y_top_pml);

 PML_left_mesh_pt  = 
  TwoDimensionalPMLHelper::create_left_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (Solid_mesh_pt, left_boundary_id, 
   Global_Parameters::N_x_left_pml, 
   Global_Parameters::Width_x_left_pml);

 PML_bottom_mesh_pt= 
  TwoDimensionalPMLHelper::create_bottom_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (Solid_mesh_pt, bottom_boundary_id, 
   Global_Parameters::N_y_bottom_pml, 
   Global_Parameters::Width_y_bottom_pml);
 
 // Add submeshes to the global mesh
 add_sub_mesh(PML_right_mesh_pt);
 add_sub_mesh(PML_top_mesh_pt);
 add_sub_mesh(PML_left_mesh_pt);
 add_sub_mesh(PML_bottom_mesh_pt);
 
// Rebuild corner PML meshes
 PML_top_right_mesh_pt    = 
  TwoDimensionalPMLHelper::create_top_right_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (PML_right_mesh_pt, PML_top_mesh_pt, 
   Solid_mesh_pt, right_boundary_id);

 PML_bottom_right_mesh_pt = 
  TwoDimensionalPMLHelper::create_bottom_right_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (PML_right_mesh_pt,PML_bottom_mesh_pt, 
   Solid_mesh_pt, right_boundary_id);
  
 PML_top_left_mesh_pt     = 
  TwoDimensionalPMLHelper::create_top_left_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (PML_left_mesh_pt, PML_top_mesh_pt, 
   Solid_mesh_pt, left_boundary_id);
 
 PML_bottom_left_mesh_pt  = 
  TwoDimensionalPMLHelper::create_bottom_left_pml_mesh
  <PMLLayerElement<ELASTICITY_ELEMENT> >
  (PML_left_mesh_pt, PML_bottom_mesh_pt, 
   Solid_mesh_pt, left_boundary_id);
 
 // Add submeshes to the global mesh
 add_sub_mesh(PML_top_right_mesh_pt);
 add_sub_mesh(PML_bottom_right_mesh_pt);
 add_sub_mesh(PML_top_left_mesh_pt);
 add_sub_mesh(PML_bottom_left_mesh_pt);
 
} // end of create_pml_meshes



//=======start_of_main==================================================
/// Driver for annular disk loaded by pressure
//======================================================================
int main(int argc, char **argv)
{

 // Start timing of the code
 Global_Parameters::T_start=TimingHelpers::timer();

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified 

 // Over-write PML layers element number in each dimension
 CommandLineArgs::specify_command_line_flag("--n_pml",
                 &Global_Parameters::N_pml_multiplier);

 // Over-write PML layers physical length in each dimension
 CommandLineArgs::specify_command_line_flag("--l_pml",
                 &Global_Parameters::L_pml_multiplier);

 // Output directory
 CommandLineArgs::specify_command_line_flag(
  "--dir", &Global_Parameters::Directory);

  
#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=0;
 CommandLineArgs::specify_command_line_flag("--max_adapt",&max_adapt);

#endif

 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Validation run?
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   oomph_info << "Using coarser resolution for self-test\n";
   
   /// Number of elements for each layer
   Global_Parameters::N_x_right_pml  = 2;
   Global_Parameters::N_y_top_pml    = 2;
   Global_Parameters::N_x_left_pml   = 2; 
   Global_Parameters::N_y_bottom_pml = 2; 

   /// Thickness of each layer
   Global_Parameters::Width_x_right_pml  = 1.0;
   Global_Parameters::Width_y_top_pml    = 1.0;
   Global_Parameters::Width_x_left_pml   = 1.0;
   Global_Parameters::Width_y_bottom_pml = 1.0;
  }

 /// Update dependent parameters
 Global_Parameters::compute_dependent_parameters();

 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory(Global_Parameters::Directory);

 // Build elasticity tensor
 Global_Parameters::E_pt=new PMLTimeHarmonicIsotropicElasticityTensor(
    Global_Parameters::Nu);


#ifdef ADAPTIVE

 //Set up the problem
 ElasticAnnulusProblem<
 ProjectablePMLTimeHarmonicLinearElasticityElement
  <QPMLTimeHarmonicLinearElasticityElement<2,3> > 
 > problem;

#else

 //Set up the problem
 ElasticAnnulusProblem<QPMLTimeHarmonicLinearElasticityElement<2,3> > 
  problem;

#endif

 
#ifdef ADAPTIVE
 
 // Solve the problem with Newton's method, allowing
 // up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);
#else
 
 // Solve the problem using Newton's method
 problem.newton_solve();
#endif

 // End timing of the code
 Global_Parameters::T_end=TimingHelpers::timer();

 // Doc solution
 problem.doc_solution(doc_info);
 
} //end of main









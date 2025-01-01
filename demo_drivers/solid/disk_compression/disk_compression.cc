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
//Driver function for a simple test elasticity problem: the
//compression of an axisymmetric disk. We also demonstrate how
//to incorporate isotropic growth into the model and how to 
//switch between different constitutive equations.
#include <iostream>
#include <fstream>
#include <cmath>

//My own includes
#include "generic.h"
#include "solid.h"

//Need to instantiate templated mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;


//============namespace_for_problem_parameters=====================
/// Global variables
//=================================================================
namespace Global_Physical_Variables
{
 /// Pointer to strain energy function
 StrainEnergyFunction* Strain_energy_function_pt;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0; 

 /// Poisson's ratio
 double Nu=0.3;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C1=1.3;

 /// Uniform pressure
 double P = 0.0;

 /// Constant pressure load
 void constant_pressure(const Vector<double> &xi,const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 } // end of pressure load


 /// Uniform volumetric expansion
 double Uniform_gamma=1.1;

 /// Growth function
 void growth_function(const Vector<double>& xi, double& gamma)
 {
  gamma = Uniform_gamma;
 }
 
} // end namespace


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//==========start_mesh=================================================
/// Elastic quarter circle sector mesh with functionality to
/// attach traction elements to the curved surface. We "upgrade"
/// the RefineableQuarterCircleSectorMesh to become an
/// SolidMesh and equate the Eulerian and Lagrangian coordinates,
/// thus making the domain represented by the mesh the stress-free 
/// configuration. 
/// \n\n
/// The member function \c make_traction_element_mesh() creates
/// a separate mesh of SolidTractionElements that are attached to the
/// mesh's curved boundary (boundary 1). 
//=====================================================================
template <class ELEMENT>
class ElasticRefineableQuarterCircleSectorMesh :
 public virtual RefineableQuarterCircleSectorMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
 /// ones so that the initial configuration is the stress-free one.
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(GeomObject* wall_pt,
                                         const double& xi_lo,
                                         const double& fract_mid,
                                         const double& xi_hi,
                                         TimeStepper* time_stepper_pt=
                                         &Mesh::Default_TimeStepper) :
  RefineableQuarterCircleSectorMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                             time_stepper_pt)
  {
   /// Make the current configuration the undeformed one by
   /// setting the nodal Lagrangian coordinates to their current
   /// Eulerian ones
   set_lagrangian_nodal_coordinates();
  }


 /// Function to create mesh made of traction elements
 void make_traction_element_mesh(SolidMesh*& traction_mesh_pt)
  {

   // Make new mesh
   traction_mesh_pt = new SolidMesh;

   // Loop over all elements on boundary 1:
   unsigned b=1;
   unsigned n_element = this->nboundary_element(b);
   for (unsigned e=0;e<n_element;e++)
    {
     // The element itself:
     FiniteElement* fe_pt = this->boundary_element_pt(b,e);
     
     // Find the index of the face of element e along boundary b
     int face_index = this->face_index_at_boundary(b,e);
     
     // Create new element
     traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
                                      (fe_pt,face_index));
    }
  }

};


//====================================================================== 
/// Uniform compression of a circular disk in a state of plane strain,
/// subject to uniform growth. 
//====================================================================== 
template<class ELEMENT>
class StaticDiskCompressionProblem : public Problem
{

public:

 /// Constructor:
 StaticDiskCompressionProblem();

 /// Run simulation: Pass case number to label output files
 void parameter_study(const unsigned& case_number);
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

private:

 /// Trace file
 ofstream Trace_file;
 
 /// Vector of pointers to nodes whose position we're tracing
 Vector<Node*> Trace_node_pt;

 /// Pointer to solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

};

//====================================================================== 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
StaticDiskCompressionProblem<ELEMENT>::StaticDiskCompressionProblem() 
{
 // Build the geometric object that describes the curvilinear
 // boundary of the quarter circle domain
 Ellipse* curved_boundary_pt = new Ellipse(1.0,1.0);

 // The curved boundary of the mesh is defined by the geometric object
 // What follows are the start and end coordinates on the geometric object:
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 // Fraction along geometric object at which the radial dividing line
 // is placed
 double fract_mid=0.5;

 //Now create the mesh using the geometric object
 Solid_mesh_pt = new ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(
  curved_boundary_pt,xi_lo,fract_mid,xi_hi);

 // Setup trace nodes as the nodes on boundary  1 (=curved boundary) 
 // in the original mesh.
 unsigned n_boundary_node = Solid_mesh_pt->nboundary_node(1);
 Trace_node_pt.resize(n_boundary_node);
 for(unsigned j=0;j<n_boundary_node;j++)
  {Trace_node_pt[j]=Solid_mesh_pt->boundary_node_pt(1,j);}

 // Refine the mesh uniformly
 Solid_mesh_pt->refine_uniformly();

 // Now construct the traction element mesh
 Solid_mesh_pt->make_traction_element_mesh(Traction_mesh_pt);
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(Solid_mesh_pt);

 // Traction mesh is second sub-mesh
 add_sub_mesh(Traction_mesh_pt);

 // Build combined "global" mesh
 build_global_mesh();


 // Pin the left edge in the horizontal direction
 unsigned n_side = mesh_pt()->nboundary_node(2);
 for(unsigned i=0;i<n_side;i++)
  {Solid_mesh_pt->boundary_node_pt(2,i)->pin_position(0);}

 // Pin the bottom in the vertical direction
 unsigned n_bottom = mesh_pt()->nboundary_node(0);
 for(unsigned i=0;i<n_bottom;i++)
  {Solid_mesh_pt->boundary_node_pt(0,i)->pin_position(1);}

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  Solid_mesh_pt->element_pt());

 //Complete the build process for elements in "bulk" solid mesh
 unsigned n_element =Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   // Set the isotropic growth function pointer
   el_pt->isotropic_growth_fct_pt()=Global_Physical_Variables::growth_function;
  }

 // Complete build process for SolidTractionElements
 n_element=Traction_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<SolidTractionElement<ELEMENT>*>
    (Traction_mesh_pt->element_pt(i));

   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }

 //Set up equation numbering scheme
 cout << "Number of equations: " <<  assign_eqn_numbers() << std::endl; 
}


//==================================================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void StaticDiskCompressionProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts = 5; 

 // Output shape of deformed body
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 //Find number of solid elements
 unsigned nelement = Solid_mesh_pt->nelement();

 // Work out volume
 double volume = 0.0;
 for(unsigned e=0;e<nelement;e++) 
  {volume+= Solid_mesh_pt->finite_element_pt(e)->size();}
 
 // Exact outer radius for linear elasticity
 double nu=Global_Physical_Variables::Nu;
 double exact_r=sqrt(Global_Physical_Variables::Uniform_gamma)*
  (1.0-Global_Physical_Variables::P/Global_Physical_Variables::E
   *((1.0+nu)*(1.0-2.0*nu)));


 // Write trace file: Problem parameters
 Trace_file << Global_Physical_Variables::P  << " " 
            << Global_Physical_Variables::Uniform_gamma << " " 
            << volume << " " 
            << exact_r << " ";
   
 // Write radii of trace nodes
 unsigned ntrace_node=Trace_node_pt.size();
 for (unsigned j=0;j<ntrace_node;j++)
  {
   Trace_file << sqrt(pow(Trace_node_pt[j]->x(0),2)+
                      pow(Trace_node_pt[j]->x(1),2)) << " ";
  }
 Trace_file << std::endl;

} // end of doc_solution


//==================================================================
/// Run the paramter study
//==================================================================
template<class ELEMENT>
void StaticDiskCompressionProblem<ELEMENT>::parameter_study(
 const unsigned& case_number)
{
 // Output
 DocInfo doc_info;

 char dirname[100];
 sprintf(dirname,"RESLT%i",case_number);

 // Set output directory
 doc_info.set_directory(dirname);

 // Step number
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 
 //Parameter incrementation
 double delta_p=0.0125;
 unsigned nstep=21; 

 // Perform fewer steps if run as self-test (indicated by nonzero number
 // of command line arguments)
 if (CommandLineArgs::Argc!=1)
  {
   nstep=3;
  }

 // Offset external pressure so that the computation sweeps 
 // over a range of positive and negative pressures
 Global_Physical_Variables::P =-delta_p*double(nstep-1)*0.5; 

 // Do the parameter study
 for(unsigned i=0;i<nstep;i++)
  {
   //Solve the problem for current load
   newton_solve();

   // Doc solution
   doc_solution(doc_info);
   doc_info.number()++;

   // Increment pressure load
   Global_Physical_Variables::P += delta_p;   
  }

} // end of parameter study


//=====start_of_main====================================================
/// Driver code for disk-compression
//======================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Define a strain energy function: Generalised Mooney Rivlin
 Global_Physical_Variables::Strain_energy_function_pt = 
  new GeneralisedMooneyRivlin(&Global_Physical_Variables::Nu,
                              &Global_Physical_Variables::C1,
                              &Global_Physical_Variables::E);
 
 // Define a constitutive law (based on strain energy function)
 Global_Physical_Variables::Constitutive_law_pt = 
  new IsotropicStrainEnergyFunctionConstitutiveLaw(
   Global_Physical_Variables::Strain_energy_function_pt);
 
 // Case 0: No pressure, generalised Mooney Rivlin
 //-----------------------------------------------
 {
    //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElement<2,3> > problem;
  
  cout << "gen. M.R.: RefineableQPVDElement<2,3>" << std::endl;
  
  //Run the simulation
  problem.parameter_study(0);
  
 } // done case 0
 
 
 // Case 1: Continuous pressure formulation with generalised Mooney Rivlin
 //------------------------------------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<
   RefineableQPVDElementWithContinuousPressure<2> > problem;
  
  cout << "gen. M.R.: RefineableQPVDElementWithContinuousPressure<2> " 
       << std::endl;
  
  //Run the simulation
  problem.parameter_study(1);

 }  // done case 1
 
 
 
 // Case 2: Discontinuous pressure formulation with generalised Mooney Rivlin
 //--------------------------------------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElementWithPressure<2> > 
   problem;
  
  cout << "gen. M.R.: RefineableQPVDElementWithPressure<2>" << std::endl;
  
  //Run the simulation
  problem.parameter_study(2);
  
 }  // done case 2
 
 
 // Change the consitutive law: Delete the old one
 delete Global_Physical_Variables::Constitutive_law_pt;
 
 // Create oomph-lib's generalised Hooke's law constitutive equation
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu,
                         &Global_Physical_Variables::E);
 
 // Case 3: No pressure, generalised Hooke's law
 //----------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElement<2,3> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElement<2,3> " << std::endl;
  
  //Run the simulation
  problem.parameter_study(3);  

 }  // done case 3
 
 // Case 4: Continuous pressure formulation with generalised Hooke's law
 //---------------------------------------------------------------------
 {
  
  //Set up the problem
  StaticDiskCompressionProblem<
   RefineableQPVDElementWithContinuousPressure<2> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElementWithContinuousPressure<2> " 
       << std::endl;
  
  //Run the simulation
  problem.parameter_study(4);
  
 }  // done case 4
 
 
 // Case 5:  Discontinous pressure formulation with generalised Hooke's law
 //------------------------------------------------------------------------
 {
  
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElementWithPressure<2> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElementWithPressure<2> " << std::endl;
  
  //Run the simulation
  problem.parameter_study(5);
    
 }  // done case 5
 
 // Clean up 
 delete Global_Physical_Variables::Constitutive_law_pt;
 Global_Physical_Variables::Constitutive_law_pt=0;
 
} // end of main









//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//Driver function for a simple test elasticity problem
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

// hierher: Should really use an example with non-uniform pressure
// add body force for self-gravity

//================================================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{
 /// Pointer to strain energy function
 StrainEnergyFunction*Strain_energy_function_pt;

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

 /// Uniform volumetric expansion
 double Uniform_gamma=1.1;

 /// Constant pressure load
 void constant_pressure(const Vector<double> &xi,const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++) {traction[i] = -P*n[i];}
 }


 /// Growth function: gamma=sqrt(g)=sqrt(g_ij)=const
 void growth_function(const Vector<double>& xi, double& gamma)
 {
  gamma = Uniform_gamma;
 }
 
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//================================================================
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
//================================================================
template <class ELEMENT>
class ElasticRefineableQuarterCircleSectorMesh :
 public virtual RefineableQuarterCircleSectorMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// \short Constructor: Build mesh and copy Eulerian coords to Lagrangian
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
 void run(const unsigned& case_number);
 
 /// Access function for the solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Access function for the mesh of surface traction elements
 SolidMesh*& traction_mesh_pt() {return Traction_mesh_pt;} 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

private:

 /// Geometric object that defines the boundary of the undeformed disk
 Ellipse* Curved_boundary_pt;

 /// Trace file
 ofstream Trace_file;
 
 /// Vector of pointers to nodes whose position we're tracing
 Vector<Node*> Trace_node_pt;

 /// Pointer to solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of traction elements
 SolidMesh*  Traction_mesh_pt;

};

//====================================================================== 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
StaticDiskCompressionProblem<ELEMENT>::StaticDiskCompressionProblem() 
{
 // Build the geometric object that describes the outer wall
 Curved_boundary_pt = new Ellipse(1.0,1.0);

 // The curved boundary of the mesh is defined by the geometric object
 // What follows are the start and end coordinates on the geometric object:
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 // Fraction along geometric object at which the radial dividing line
 // is placed
 double fract_mid=0.5;

 //Now create the mesh using the geometric object
 solid_mesh_pt() = new ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(
  Curved_boundary_pt,xi_lo,fract_mid,xi_hi);

 // Setup trace nodes as the nodes on boundary  1 (=curved boundary) 
 // in the original mesh (they exist under any refinement!) 
 unsigned n_boundary_node = solid_mesh_pt()->nboundary_node(1);
 Trace_node_pt.resize(n_boundary_node);
 for(unsigned n=0;n<n_boundary_node;n++)
  {Trace_node_pt[n]=solid_mesh_pt()->boundary_node_pt(1,n);}

 // Refine the mesh uniformly
 solid_mesh_pt()->refine_uniformly();

 // Mesh has been adapted: Need to setup boundary info again
 solid_mesh_pt()->setup_boundary_element_info();

 // Now construct the traction element mesh
 solid_mesh_pt()->make_traction_element_mesh(traction_mesh_pt());
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Traction mesh is first sub-mesh
 add_sub_mesh(traction_mesh_pt());

 // Build combined "global" mesh
 build_global_mesh();

 // Pin the bottom in the vertical direction
 unsigned n_bottom = mesh_pt()->nboundary_node(0);
 //Loop over the nodes
 for(unsigned i=0;i<n_bottom;i++)
  {solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(1);}

 // Pin the left edge in the horizontal direction
 unsigned n_side = mesh_pt()->nboundary_node(2);
 //Loop over the node
 for(unsigned i=0;i<n_side;i++)
  {solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(0);}


 //Find number of elements in solid mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 //Loop over the elements in the main mesh
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   // Set the isotropic growth function pointer
   el_pt->isotropic_growth_fct_pt()=Global_Physical_Variables::growth_function;
  }

 // Pin the redundant solid pressures
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 //Find number of elements in traction mesh
 n_element=traction_mesh_pt()->nelement();
 //Loop over the elements in the traction element mesh
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<SolidTractionElement<ELEMENT>*>
    (traction_mesh_pt()->element_pt(i));

   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 
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
 solid_mesh_pt()->output(some_file,npts);
 some_file.close();

 //Find number of solid elements
 unsigned nelement = solid_mesh_pt()->nelement();

 // Work out volume
 double volume = 0.0;
 for(unsigned e=0;e<nelement;e++) 
  {volume+= solid_mesh_pt()->finite_element_pt(e)->size();}
 
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


 cout << "Doced solution for step " 
      << doc_info.number() 
      << std::endl << std::endl << std::endl;
}


//==================================================================
/// Run the problem
//==================================================================
template<class ELEMENT>
void StaticDiskCompressionProblem<ELEMENT>::run(const unsigned& case_number)
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
 
 // Initial parameter values
 Global_Physical_Variables::Uniform_gamma=1.1; 
 
 //Parameter incrementation
 double delta_p=0.025;
 unsigned nstep=21; 
 if (CommandLineArgs::Argc!=1)
  {
   nstep=3;
  }
 Global_Physical_Variables::P =-delta_p*double(nstep-1)*0.5; 
 for(unsigned i=0;i<nstep;i++)
  {
   //Solve the problem with Newton's method
   newton_solve();

   // Doc solution
   doc_solution(doc_info);
   doc_info.number()++;

   // Increment pressure load
   Global_Physical_Variables::P += delta_p;   
  }

}

//======================================================================
/// Driver for simple elastic problem
//======================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Initialise physical parameters
 Global_Physical_Variables::E = 2.1;  // ADJUST
 Global_Physical_Variables::Nu = 0.3; // ADJUST
 Global_Physical_Variables::C1 = 1.3; // ADJUST
 
 
 
 // Define a strain energy function: Generalised Mooney Rivlin
 Global_Physical_Variables::Strain_energy_function_pt = 
  new GeneralisedMooneyRivlin(Global_Physical_Variables::Nu,
                              Global_Physical_Variables::C1,
                              Global_Physical_Variables::E);
 
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
  problem.run(0);
  
 }
 
 
 // Case 1: Continuous pressure formulation with generalised Mooney Rivlin
 //------------------------------------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<
   RefineableQPVDElementWithContinuousPressure<2> > problem;
  
  cout << "gen. M.R.: RefineableQPVDElementWithContinuousPressure<2> " 
       << std::endl;
  
  //Run the simulation
  problem.run(1);
 }
 
 
 
 // Case 2: Discontinuous pressure formulation with generalised Mooney Rivlin
 //--------------------------------------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElementWithPressure<2> > 
   problem;
  
  cout << "gen. M.R.: RefineableQPVDElementWithPressure<2>" << std::endl;
  
  //Run the simulation
  problem.run(2);
  
 }
 
 
 // Change the consitutive law
 delete Global_Physical_Variables::Constitutive_law_pt;
 
 // "Big G" Linear constitutive equations:
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(Global_Physical_Variables::Nu,
                         Global_Physical_Variables::E);
 
 // Case 3: No pressure, generalised Hooke's law
 //----------------------------------------------
 {
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElement<2,3> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElement<2,3> " << std::endl;
  
  //Run the simulation
  problem.run(3);  
 }
 
 // Case 4: Continuous pressure formulation with generalised Hooke's law
 //---------------------------------------------------------------------
 {
  
  //Set up the problem
  StaticDiskCompressionProblem<
   RefineableQPVDElementWithContinuousPressure<2> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElementWithContinuousPressure<2> " 
       << std::endl;
  
  //Run the simulation
  problem.run(4);
  
 }
 
 
 // Case 5:  Discontinous pressure formulation with generalised Hooke's law
 //------------------------------------------------------------------------
 {
  
  //Set up the problem
  StaticDiskCompressionProblem<RefineableQPVDElementWithPressure<2> > problem;
  
  cout << "gen. Hooke: RefineableQPVDElementWithPressure<2> " << std::endl;
  
  //Run the simulation
  problem.run(5);
  
  
 }
 
 // Clean up 
 delete Global_Physical_Variables::Constitutive_law_pt;
 Global_Physical_Variables::Constitutive_law_pt=0;
 
}









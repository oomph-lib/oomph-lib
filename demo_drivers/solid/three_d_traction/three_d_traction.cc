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
// Driver for elastic deformation of a cuboidal domain
 
 
// Generic oomph-lib headers
#include "generic.h"

// Solid mechanics
#include "solid.h"

// The mesh 
#include "meshes/simple_cubic_mesh.h"

using namespace std;

using namespace oomph;


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


#ifdef NOREFINE

//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticCubicMesh : public virtual SimpleCubicMesh<ELEMENT>, 
                         public virtual SolidMesh
{

public:

 ///  Constructor: 
 ElasticCubicMesh(const unsigned &nx, const unsigned &ny, 
                  const unsigned &nz,
                  const double &a, const double &b, 
                  const double &c,
                  TimeStepper* time_stepper_pt = 
                  &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,-a,a,-b,b,-c,c,time_stepper_pt),
  SolidMesh()
  {
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~ElasticCubicMesh() { }

};


#else

//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class RefineableElasticCubicMesh : public virtual SimpleCubicMesh<ELEMENT>, 
                                  public virtual RefineableBrickMesh<ELEMENT>,
                                  public virtual SolidMesh
{

public:

 ///  Constructor: 
 RefineableElasticCubicMesh(const unsigned &nx, const unsigned &ny, 
                            const unsigned &nz,
                            const double &a, const double &b, 
                            const double &c,
                            TimeStepper* time_stepper_pt = 
                            &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,-a,a,-b,b,-c,c,time_stepper_pt),
  RefineableBrickMesh<ELEMENT>(), SolidMesh()
  {
   
   this->setup_octree_forest();

   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~RefineableElasticCubicMesh() { }

};

#endif


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




//================================================================
/// Global variables
//================================================================
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

}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//====================================================================== 
/// Deformation of elastic block
//====================================================================== 
template<class ELEMENT>
class BlockCompressionProblem : public Problem
{

public:

 /// Constructor:
 BlockCompressionProblem();

 /// Run simulation.
 void run(const std::string &dirname);
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Update function (empty)
 void actions_after_newton_solve() {}

#ifndef NOREFINE

 /// Actions before adapt
 void actions_before_adapt() 
  {
   // Kill the traction elements and wipe surface mesh
   delete_traction_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }


 /// Actions after adapt
 void actions_after_adapt() 
  {
   // Create traction elements
   create_traction_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
   
   // Pin the redundant solid pressures
   PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
    Solid_mesh_pt->element_pt());
  }

#endif

 /// Delete traction elements and wipe the  traction meshes
 void delete_traction_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Traction_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Traction_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Traction_mesh_pt->flush_element_and_node_storage();   
  }
   



 /// Create traction elements
 void create_traction_elements()
  {

   unsigned b=5;
   unsigned n_element = Solid_mesh_pt->nboundary_element(b);
   for (unsigned e=0;e<n_element;e++)
    {
     // The element itself:
     FiniteElement* fe_pt = Solid_mesh_pt->boundary_element_pt(b,e);
     
     // Find the index of the face of element e along boundary b
     int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
     
#ifdef NOREFINE

     // Create new element
     Traction_mesh_pt->add_element_pt(
      new SolidTractionElement<ELEMENT>
      (fe_pt,face_index));

#else

     // Create new element
     Traction_mesh_pt->add_element_pt(
      new RefineableSolidTractionElement<ELEMENT>
      (fe_pt,face_index));

#endif

    }

   // Complete build process for SolidTractionElements
   n_element=Traction_mesh_pt->nelement();
   for(unsigned i=0;i<n_element;i++)
    {

#ifdef NOREFINE

     //Cast to a solid traction element
     SolidTractionElement<ELEMENT> *el_pt = 
      dynamic_cast<SolidTractionElement<ELEMENT>*>
      (Traction_mesh_pt->element_pt(i));

#else

     //Cast to a solid traction element
     RefineableSolidTractionElement<ELEMENT> *el_pt = 
      dynamic_cast<RefineableSolidTractionElement<ELEMENT>*>
      (Traction_mesh_pt->element_pt(i));

#endif

     
     //Set the traction function
     el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
    }
   
  }

 ///  Update before solve: Empty
 void actions_before_newton_solve(){}

private:

#ifdef NOREFINE

 /// Pointer to solid mesh
 ElasticCubicMesh<ELEMENT>* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 RefineableElasticCubicMesh<ELEMENT>* Solid_mesh_pt;

#endif

 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

};

//====================================================================== 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
BlockCompressionProblem<ELEMENT>::BlockCompressionProblem() 
{
 double a = 1.0, b = 1.0, c = 1.0;
 unsigned n_x = 2, n_y = 2, n_z = 2;


#ifdef NOREFINE

 //Now create the mesh 
 Solid_mesh_pt = new ElasticCubicMesh<ELEMENT>(n_x,n_y,n_z,a,b,c);

#else

 //Now create the mesh 
 Solid_mesh_pt = new RefineableElasticCubicMesh<ELEMENT>(n_x,n_y,n_z,a,b,c);

 // Set error estimator
 Solid_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Solid_mesh_pt->max_permitted_error()=0.02;
 Solid_mesh_pt->min_permitted_error()=0.001;

#endif


 // Make new mesh
 Traction_mesh_pt = new SolidMesh;
 create_traction_elements();

 // Solid mesh is first sub-mesh
 add_sub_mesh(Solid_mesh_pt);
 
 // Traction mesh is second sub-mesh
 add_sub_mesh(Traction_mesh_pt);

 // Build combined "global" mesh
 build_global_mesh();

 //Loop over the elements in the mesh to set parameters/function pointers
 unsigned  n_element=Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
  }
 


 // Setup boundary conditions
 
 //Loop over nodes in the bottom boundary
 {
  unsigned b=0;
  unsigned n_node = Solid_mesh_pt->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    //Pin all nodes 
    for(unsigned i=0;i<3;i++)
     {
      Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(i); 
     }
   }
 }
 
 // Pin the redundant solid pressures
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  Solid_mesh_pt->element_pt());
 
 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 

} 


//==================================================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void BlockCompressionProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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


    
 // Output traction
 //----------------
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Traction_mesh_pt->output(some_file,npts);
 some_file.close();

}
 

//==================================================================
/// Run the problem
//==================================================================
template<class ELEMENT>
void BlockCompressionProblem<ELEMENT>::run(const std::string &dirname)
{

 // Output
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory(dirname);

 // Step number
 doc_info.number()=0;
  
 // Doc initial configuration
 doc_solution(doc_info);
 doc_info.number()++;

 // Initial parameter values

 // Gravity and pressure
 Global_Physical_Variables::P=0.05; 


#ifndef NOREFINE

 //Refine according to a pattern
 Vector<unsigned> refine_pattern(2);
 refine_pattern[0] = 0; refine_pattern[1] = 7;
 refine_selected_elements(0,refine_pattern);
 
 #endif

 //Parameter incrementation
 unsigned nstep=10; 
 if (CommandLineArgs::Argc!=1)
  {
   std::cout << "Validation -- only doing one step" << std::endl;
   nstep=1;
  }

 for(unsigned i=0;i<nstep;i++)
  {


#ifdef NOREFINE

   // Solve the problem with Newton's method
   newton_solve();

#else

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   unsigned max_adapt=1;
   newton_solve(max_adapt);

#endif

   // Doc solution
   doc_solution(doc_info);
   doc_info.number()++;

   //Increase the body force and pressure
   Global_Physical_Variables::P += 0.05;
  }

}


//======================================================================
/// Driver for simple elastic problem
//======================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 //Initialise physical parameters
 Global_Physical_Variables::E  = 2.1; 
 Global_Physical_Variables::Nu = 0.4; 
 Global_Physical_Variables::C1 = 1.3; 
 
 
 // Define a strain energy function: Generalised Mooney Rivlin
 Global_Physical_Variables::Strain_energy_function_pt = 
  new GeneralisedMooneyRivlin(&Global_Physical_Variables::Nu,
                              &Global_Physical_Variables::C1,
                              &Global_Physical_Variables::E);
 
 // Define a constitutive law (based on strain energy function)
 Global_Physical_Variables::Constitutive_law_pt = 
  new IsotropicStrainEnergyFunctionConstitutiveLaw(
   Global_Physical_Variables::Strain_energy_function_pt);
 
#ifdef NOREFINE

 //Set up the problem with pure displacement formulation
 BlockCompressionProblem<QPVDElement<3,3> > problem;
 problem.run("RESLT");
 
#else

 //Set up the problem with pure displacement formulation
 BlockCompressionProblem<RefineableQPVDElement<3,3> > problem;
 problem.run("RESLT");

#endif

 
}






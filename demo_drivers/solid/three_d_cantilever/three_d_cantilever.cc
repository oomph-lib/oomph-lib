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
//Driver for 3D Airy cantilever beam problem

//#include <fenv.h>

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// The mesh 
#include "meshes/simple_cubic_mesh.template.h"

// The mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;
using namespace oomph;


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//======================start_mesh=========================================
/// Simple quarter tube mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class RefineableElasticQuarterTubeMesh : 
 public virtual RefineableQuarterTubeMesh<ELEMENT>, 
 public virtual SolidMesh 
{
 
public:
 
 ///  Constructor: 
 RefineableElasticQuarterTubeMesh(GeomObject* wall_pt,
                                  const Vector<double>& xi_lo,
                                  const double& fract_mid,
                                  const Vector<double>& xi_hi,
                                  const unsigned& nlayer,
                                  TimeStepper* time_stepper_pt=
                                  &Mesh::Default_TimeStepper) :
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt),
  RefineableQuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                     nlayer,time_stepper_pt) 
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }
 
 /// Empty Destructor
 virtual ~RefineableElasticQuarterTubeMesh() { }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=========================================================================
/// Simple quarter tube mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticQuarterTubeMesh : public virtual QuarterTubeMesh<ELEMENT>, 
                               public virtual SolidMesh 
{
 
public:
 
 ///  Constructor: 
 ElasticQuarterTubeMesh(GeomObject* wall_pt,
                 const Vector<double>& xi_lo,
                 const double& fract_mid,
                 const Vector<double>& xi_hi,
                 const unsigned& nlayer,
                 TimeStepper* time_stepper_pt=
                 &Mesh::Default_TimeStepper) :
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt) 
  {
   //Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~ElasticQuarterTubeMesh() { }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Length of beam
 double L=10.0;

 /// Pointer to strain energy function
 StrainEnergyFunction* Strain_energy_function_pt=0;

 /// First "Mooney Rivlin" coefficient 
 double C1=1.3;

 /// Second "Mooney Rivlin" coefficient 
 double C2=1.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
  b[2]=0.0;
 }
 
} //end namespace





//================================================================
/// Extension of global variables for self tests
//================================================================
namespace Global_Physical_Variables
{

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;

} 



//=============begin_problem============================================ 
/// Problem class for the 3D cantilever "beam" structure.
//====================================================================== 
template<class ELEMENT>
class CantileverProblem : public Problem
{

public:

 /// Constructor:
 CantileverProblem(); 
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Actions before adapt. Empty
 void actions_before_adapt(){}

 /// Actions after adapt
 void actions_after_adapt()
  {
   // Pin the redundant solid pressures (if any)
   PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
    mesh_pt()->element_pt());
  }

 /// Doc the solution
 void doc_solution();

#ifdef REFINE

 /// Access function for the mesh
 RefineableElasticQuarterTubeMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

#else

 /// Access function for the mesh
 ElasticQuarterTubeMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<ElasticQuarterTubeMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }
 
#endif

 /// Run extended tests -- doc in RESLTi_case
 void run_tests(const unsigned& i_case,
                const bool& incompress, 
                const bool& use_fd);

private:

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
CantileverProblem<ELEMENT>::CantileverProblem() 
{

 // Create geometric object that defines curvilinear boundary of
 // beam: Elliptical tube with half axes = radius = 1.0
 double radius=1.0;
 GeomObject* wall_pt=new EllipticalTube(radius,radius);
 
 // Bounding Lagrangian coordinates
 Vector<double> xi_lo(2);
 xi_lo[0]=0.0;
 xi_lo[1]=0.0;

 Vector<double> xi_hi(2);
 xi_hi[0]= Global_Physical_Variables::L;
 xi_hi[1]=2.0*atan(1.0);


#ifdef REFINE

 // # of layers
 unsigned nlayer=6;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;

 //Now create the mesh 
 Problem::mesh_pt() = new RefineableElasticQuarterTubeMesh<ELEMENT>
  (wall_pt,xi_lo,frac_mid,xi_hi,nlayer);
 
 // Set error estimator
 dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>* >(
  mesh_pt())->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error()=0.05;
 mesh_pt()->min_permitted_error()=0.005;

#else

 // # of layers
 unsigned nlayer=6;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;

 //Now create the mesh 
 Problem::mesh_pt() = new ElasticQuarterTubeMesh<ELEMENT>
  (wall_pt,xi_lo,frac_mid,xi_hi,nlayer);

#endif
 

 // Complete build of elements
 unsigned n_element=mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;

   // Set the body force
   el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;

   // Material is incompressble: Use incompressible displacement/pressure
   // formulation (if the element is pressure based, that is!)
   PVDEquationsWithPressure<3>* cast_el_pt = 
    dynamic_cast<PVDEquationsWithPressure<3>*>(mesh_pt()->element_pt(i));
   if (cast_el_pt!=0)
    {
     cast_el_pt->set_incompressible();
    }

  } // done build of elements

 
 // Pin the left boundary (boundary 0) in all directions
 unsigned b=0; 
 unsigned n_side = mesh_pt()->nboundary_node(b);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
   mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
   mesh_pt()->boundary_node_pt(b,i)->pin_position(2);
  }

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  mesh_pt()->element_pt());

 //Assign equation numbers
 assign_eqn_numbers();

 // Prepare output directory
 Doc_info.set_directory("RESLT");
  
} //end of constructor



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of deformed body
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Increment label for output files
 Doc_info.number()++;

} //end doc




//==============start_run_tests========================================
/// Run tests
//==================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::run_tests(const unsigned& i_case,
                                           const bool& incompress,
                                           const bool& use_fd)
{

 // Set output directory
 char dirname[100];   

#ifdef REFINE
 sprintf(dirname,"RESLT_refine%i",i_case);
#else
 sprintf(dirname,"RESLT_norefine%i",i_case);
#endif

 // Prepare output
 Doc_info.set_directory(dirname);

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element=mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Get Jacobian by FD?
   if(use_fd)
    {
     el_pt->enable_evaluate_jacobian_by_fd();
    }
   else
    {
     el_pt->disable_evaluate_jacobian_by_fd();
    }
   
   // Is the material actually not incompressible?
   if (!incompress)
    {  
     PVDEquationsWithPressure<3>* cast_el_pt = 
      dynamic_cast<PVDEquationsWithPressure<3>*>(
       mesh_pt()->element_pt(i));
     if (cast_el_pt!=0)
      {
       cast_el_pt->set_compressible();
      }
    }
  }


 // Doc solution
 doc_solution();

 // Initial values for parameter values
 Global_Physical_Variables::Gravity=0.0;
 
 //Parameter incrementation
 unsigned nstep=1; 

 double g_increment=1.0e-5;   
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment load
   Global_Physical_Variables::Gravity+=g_increment;

#ifdef REFINE

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   unsigned max_adapt=1;   
   newton_solve(max_adapt);

#else

   // Solve it
   newton_solve();

#endif

   // Doc solution
   doc_solution();
    
  }

}


//=======start_of_main==================================================
/// Driver for 3D cantilever beam loaded by gravity
//======================================================================
int main(int argc, char* argv[])
{

 // Run main demo code if no command line arguments are specified
 if (argc==1)
  {
   
   // Create incompressible Mooney Rivlin strain energy function
   Global_Physical_Variables::Strain_energy_function_pt = 
    new MooneyRivlin(&Global_Physical_Variables::C1,
                     &Global_Physical_Variables::C2);
   
   // Define a constitutive law (based on strain energy function)
   Global_Physical_Variables::Constitutive_law_pt = 
    new IsotropicStrainEnergyFunctionConstitutiveLaw(
     Global_Physical_Variables::Strain_energy_function_pt);
   
   //Set up the problem with continous pressure/displacement
   CantileverProblem<RefineableQPVDElementWithContinuousPressure<3> > problem;
   
   // Doc solution
   problem.doc_solution();
   
   // Initial values for parameter values
   Global_Physical_Variables::Gravity=0.0;
   
   //Parameter incrementation
   unsigned nstep=10; 
   
   double g_increment=5.0e-4;   
   for(unsigned i=0;i<nstep;i++)
    {
     // Increment load
     Global_Physical_Variables::Gravity+=g_increment;
     
     // Solve the problem with Newton's method, allowing
     // up to max_adapt mesh adaptations after every solve.
     unsigned max_adapt=1;   
     problem.newton_solve(max_adapt);
     
     // Doc solution
     problem.doc_solution();    
    }

  } // end main demo code

 // Run self-tests
 else
  {
   
   // Additional test cases combining adaptive/non-adaptive
   // elements with displacement/displacement-pressure
   // discretisation and various constitutive equations
   
   
   // Initialise flag for FD evaluation
   bool use_fd=false; 
   
   // Number of cases per implementation
   unsigned ncase=5;
   
   // Is the material incomressible?
   bool incompress=true;
   
   // Loop over fd and analytical implementation
   for (unsigned i=0;i<2;i++) 
    {
     
     // Generalised Hookean constitutive equations
     //-------------------------------------------
     {
      Global_Physical_Variables::Constitutive_law_pt = 
       new GeneralisedHookean(&Global_Physical_Variables::Nu,
                              &Global_Physical_Variables::E);
      
      incompress=false;
      
#ifdef REFINE
      {
       //Set up the problem with pure displacement based elements
       CantileverProblem<RefineableQPVDElement<3,3> > problem; 
       problem.run_tests(0+i*ncase,incompress,use_fd);
      }
#else   
      {
       //Set up the problem with pure displacement based elements
       CantileverProblem<QPVDElement<3,3> > problem;
       problem.run_tests(0+i*ncase,incompress,use_fd);
      }
#endif
      
      
#ifdef REFINE
      {
       //Set up the problem with continous pressure/displacement
       CantileverProblem<RefineableQPVDElementWithContinuousPressure<3> > problem;
       problem.run_tests(1+i*ncase,incompress,use_fd);
      }
#else
      {
       //Set up the problem with continous pressure/displacement
       CantileverProblem<QPVDElementWithContinuousPressure<3> > problem;
       problem.run_tests(1+i*ncase,incompress,use_fd);
      }
#endif
      
      
#ifdef REFINE    
      {
       //Set up the problem with discontinous pressure/displacement
       CantileverProblem<RefineableQPVDElementWithPressure<3> > problem;
       problem.run_tests(2+i*ncase,incompress,use_fd);
     }
#else
      {
       //Set up the problem with discontinous pressure/displacement
       CantileverProblem<QPVDElementWithPressure<3> > problem;
       problem.run_tests(2+i*ncase,incompress,use_fd);
      }
#endif
      
      delete Global_Physical_Variables::Constitutive_law_pt;
      Global_Physical_Variables::Constitutive_law_pt=0;
     }
     
     
     
     // Incompressible Mooney Rivlin
     //-----------------------------
     {
      Global_Physical_Variables::Strain_energy_function_pt = 
       new MooneyRivlin(&Global_Physical_Variables::C1,
                        &Global_Physical_Variables::C2);
    
      // Define a constitutive law (based on strain energy function)
      Global_Physical_Variables::Constitutive_law_pt = 
       new IsotropicStrainEnergyFunctionConstitutiveLaw(
        Global_Physical_Variables::Strain_energy_function_pt);
    
      incompress=true;


#ifdef REFINE
      {
       //Set up the problem with continous pressure/displacement
       CantileverProblem<RefineableQPVDElementWithContinuousPressure<3> > problem;
       problem.run_tests(3+i*ncase,incompress,use_fd);
      }
#else
      {
       //Set up the problem with continous pressure/displacement
       CantileverProblem<QPVDElementWithContinuousPressure<3> > problem;     
       problem.run_tests(3+i*ncase,incompress,use_fd);
      }
#endif


    
#ifdef REFINE
      {
       //Set up the problem with discontinous pressure/displacement
       CantileverProblem<RefineableQPVDElementWithPressure<3> > problem;
       problem.run_tests(4+i*ncase,incompress,use_fd);
      }
#else
      {
       //Set up the problem with discontinous pressure/displacement
       CantileverProblem<QPVDElementWithPressure<3> > problem;
       problem.run_tests(4+i*ncase,incompress,use_fd);
      }
#endif
    
      delete  Global_Physical_Variables::Strain_energy_function_pt;
      Global_Physical_Variables::Strain_energy_function_pt=0;
    
      delete Global_Physical_Variables::Constitutive_law_pt;
      Global_Physical_Variables::Constitutive_law_pt=0;
     }

 
     use_fd=true;
     std::cout << "\n\n\n CR Total fill_in... : bla \n\n\n " << std::endl;
   
    }
  }


} //end of main







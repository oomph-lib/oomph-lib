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
//Driver for compressed square

//Oomph-lib includes
#include "generic.h"
#include "solid.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//#define REFINE

namespace oomph
{

//=================start_wrapper==================================
/// Wrapper class for solid elements to modify their output 
/// functions.
//================================================================
template <class ELEMENT>
class MySolidElement : public virtual ELEMENT
{

public:

 /// Constructor: Call constructor of underlying element
 MySolidElement() : ELEMENT() {};

 /// Overload output function:
 void output(std::ostream &outfile, const unsigned &n_plot)
  {

   // Element dimension
   unsigned el_dim = this->dim();

   Vector<double> s(el_dim);
   Vector<double> x(el_dim);
   Vector<double> xi(el_dim);
   DenseMatrix<double> sigma(el_dim,el_dim);
   
   switch(el_dim)
    {
     
    case 2:

     //Tecplot header info 
     outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
     
     //Loop over element nodes
     for(unsigned l2=0;l2<n_plot;l2++)
      {
       s[1] = -1.0 + l2*2.0/(n_plot-1);
       for(unsigned l1=0;l1<n_plot;l1++)
        {
         s[0] = -1.0 + l1*2.0/(n_plot-1);
         
         // Get Eulerian coordinates and stress
         this->interpolated_x(s,x);
         this->interpolated_xi(s,xi);

         this->get_stress(s,sigma);

         //Output the x,y,..
         for(unsigned i=0;i<el_dim;i++) 
          {outfile << x[i] << " ";}

         for(unsigned i=0;i<el_dim;i++) 
          {outfile << x[i]-xi[i] << " ";}


         // hierher change to physical stress:
         // sigma_phys_up(i,j)=sigma_up(i,j) sqrt(g_down(ii))  sqrt(g_down(jj))
 
         // Output stress
         outfile << sigma(0,0) << " "
                 << sigma(1,0) << " "
                 << sigma(1,1) << " "
                 << std::endl;
        }
      }

     break;
     
    default:

     std::ostringstream error_message;
     error_message << "Output for dim !=2 not implemented" << std::endl;
     throw OomphLibError(error_message.str(),"MySolidElement::output()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  
  }

};



//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<MySolidElement<ELEMENT> > :
 public virtual FaceGeometry<ELEMENT>
{
};


}






///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Pointer to strain energy function
 StrainEnergyFunction*Strain_energy_function_pt;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C1=1.3;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C2=1.3;

 /// Height of domain
 double H=1.0;

 /// Length of domain
 double L=1.0;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.4999;

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi,
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }
 
} //end namespace



//=============begin_problem============================================ 
/// Problem class for the cantilever "beam" structure.
//====================================================================== 
template<class ELEMENT>
class CompressedSquareProblem : public Problem
{

public:

 /// Constructor:
 CompressedSquareProblem(const bool& incompress);
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}


 /// Actions before adapt
 void actions_before_adapt(){}

 /// Actions after adapt
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

 /// Run the job -- doc in RESLTi_case
 void run_it(const unsigned& i_case);

private:

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
CompressedSquareProblem<ELEMENT>::CompressedSquareProblem(
 const bool& incompress)
{

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=5;

 // # of elements in y-direction
 unsigned n_y=5;

 // Domain length in x-direction
 double l_x=Global_Physical_Variables::L;

 // Domain length in y-direction
 double l_y=Global_Physical_Variables::H;

 // Shift mesh downwards so that centreline is at y=0:
 Vector<double> origin(2);
 origin[0]=0.0;
 origin[1]=-0.5*l_y;

#ifdef REFINE

 //Now create the mesh 
 mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y,origin);
 
 // Set error estimator
  dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>*>(
   mesh_pt())->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
#else
 
 //Now create the mesh 
 mesh_pt() = new ElasticRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y,origin);

#endif
 

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element=mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;

   //Set the body force
   el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
  
   // Is it incompressible -- hierher switch to const eqn
   if (incompress)
    {     
     PVDEquationsWithPressure<2>* test_pt = 
      dynamic_cast<PVDEquationsWithPressure<2>*>(mesh_pt()->element_pt(i));
     if (test_pt!=0)
      {
       test_pt->incompressible()=true;
      }
    }
  }


 // Choose a control node: The last node in the solid mesh
 unsigned nnod=mesh_pt()->nnode();
 Trace_node_pt=mesh_pt()->node_pt(nnod-1);

#ifdef REFINE
 
 // Do a nonuniform refinement
Vector<unsigned> elements_to_be_refined(3);
elements_to_be_refined[0]=1;
elements_to_be_refined[1]=3;
elements_to_be_refined[2]=5;
 dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>*>(
  mesh_pt())->refine_selected_elements(elements_to_be_refined);

#endif
 
 
// Pin the left and right boundaries (1 and 2) in the horizontal directions
for (unsigned b=1;b<4;b+=2)
{
 unsigned n_side = mesh_pt()->nboundary_node(b);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   dynamic_cast<SolidNode*>(mesh_pt()->boundary_node_pt(b,i))->pin_position(0);
  }
}


// Pin the bottom boundary (0) in the vertical direction
unsigned b=0;
{
 unsigned n_side = mesh_pt()->nboundary_node(b);
 
//Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   dynamic_cast<SolidNode*>(mesh_pt()->boundary_node_pt(b,i))->pin_position(1);
  }
}

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  mesh_pt()->element_pt());

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 

 
} //end of constructor




//=====================start_of_actions_after_adapt=======================
///  Actions after adapt
//========================================================================
template<class ELEMENT>
void CompressedSquareProblem<ELEMENT>::actions_after_adapt()
{
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
mesh_pt()->element_pt());
 
}// end of actions_after_adapt



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void CompressedSquareProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of and stress in deformed body
 //--------------------------------------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Write trace file: Load/displacement characteristics
 Trace_file << Global_Physical_Variables::Gravity  << " " 
            << Trace_node_pt->x(0) << " " 
            << Trace_node_pt->x(1) << " " 
            << std::endl;

 // Increment label for output files
 Doc_info.number()++;

} //end doc





//==============start_run_it========================================
/// Run it
//==================================================================
template<class ELEMENT>
void CompressedSquareProblem<ELEMENT>::run_it(const unsigned& i_case)
{

#ifdef TIME_SOLID_JAC
  PVDEquationsBase<2>::Solid_timer.reset();
#endif

 // Set output directory
 char dirname[100];   

#ifdef REFINE
 sprintf(dirname,"RESLT_refine%i",i_case);
#else
 sprintf(dirname,"RESLT_norefine%i",i_case);
#endif

 Doc_info.set_directory(dirname);

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);


 // Doc solution
 doc_solution();

 // Initial values for parameter values
 Global_Physical_Variables::Gravity=0.0;
 
 //Parameter incrementation
 unsigned nstep=1; 
 double g_increment=1.0e2;   
 for(unsigned i=0;i<nstep;i++)
  {
   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.

   // Initial values for parameter values
   Global_Physical_Variables::Gravity+=g_increment;
 
#ifdef REFINE

   // Max. number of adaptations per solve
   unsigned max_adapt=1;
   
   newton_solve(max_adapt);

#else
 
   newton_solve();

#endif

   // Doc solution
   doc_solution();
   
  }

}


//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main()
{
 
 // Is the material incomressible hierher switch to const eqn
 bool incompress=false;

 // Nu values
 Vector<double> nu_value(3);
 nu_value[0]=0.45;
 nu_value[1]=0.499999;
 nu_value[2]=0.5;
 
 unsigned case_number=0;
 
 // Generalised Hookean constitutive equations
 //-------------------------------------------
 {
  Global_Physical_Variables::Constitutive_law_pt = 
   new GeneralisedHookean(Global_Physical_Variables::Nu,
                          Global_Physical_Variables::E);
  
  
  // Loop over nu values
  for (unsigned i=0;i<3;i++)
   {
    // Set Poisson ratio
    Global_Physical_Variables::Nu=nu_value[i];
    
    std::cout << "===========================================" << std::endl;
    std::cout << "Doing Nu=" <<  Global_Physical_Variables::Nu << std::endl;
    std::cout << "===========================================" << std::endl;
    
    
    incompress=false;
    
    // nly do displacement-based elements if nu is not 1/2
    if (i!=2) 
     {
      
      std::cout 
       << "Doing Generalised Hookean with displacement formulation: Case "
       << case_number << std::endl;
      
#ifdef REFINE
      {
       //Set up the problem with pure displacement based elements
       CompressedSquareProblem<MySolidElement<RefineableQPVDElement<2,3> > > 
        problem(incompress);
       
       problem.run_it(case_number);
      }
#else   
      {
       //Set up the problem with pure displacement based elements
       CompressedSquareProblem<MySolidElement<QPVDElement<2,3> > > 
        problem(incompress);
       
       problem.run_it(case_number);
      }
#endif
     }
    
    
    case_number++;
    std::cout 
     << "Doing Generalised Hookean with displacement/cont pressure formulation: Case "
     << case_number << std::endl;

#ifdef REFINE
    {
     //Set up the problem with continous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      RefineableQPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     problem.run_it(case_number);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     problem.run_it(case_number);
    }
#endif
    
   
    case_number++;
    std::cout 
     << "Doing Generalised Hookean with displacement/discont pressure formulation: Case "
     << case_number << std::endl; 

#ifdef REFINE    
    {
     //Set up the problem with discontinous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      RefineableQPVDElementWithPressure<2> > > problem(incompress); 
     
     problem.run_it(case_number);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithPressure<2> > > problem(incompress); 
     
     problem.run_it(case_number);
    }
#endif


    incompress=true;

    case_number++;
    std::cout 
     << "Doing Generalised Hookean with displacement/cont pressure, incompressible formulation: Case "
     << case_number << std::endl;

#ifdef REFINE
    {
     //Set up the problem with continous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      RefineableQPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     problem.run_it(case_number);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     problem.run_it(case_number);
    }
#endif
   
    case_number++;
    std::cout 
     << "Doing Generalised Hookean with displacement/discont pressure, incompressible formulation: Case "
     << case_number << std::endl; 

#ifdef REFINE    
    {
     //Set up the problem with discontinous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      RefineableQPVDElementWithPressure<2> > > problem(incompress); 
     
     problem.run_it(case_number);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithPressure<2> > > problem(incompress); 
     
     problem.run_it(case_number);
    }
#endif
    


   } // end of loop over Poisson ratios
  
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;
 }
   
 // Incompressible Mooney Rivlin
 //-----------------------------
 {
  Global_Physical_Variables::Strain_energy_function_pt = 
   new MooneyRivlin(Global_Physical_Variables::C1,
                    Global_Physical_Variables::C2);
  
  // Define a constitutive law (based on strain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
   new IsotropicStrainEnergyFunctionConstitutiveLaw(
    Global_Physical_Variables::Strain_energy_function_pt);
    
  incompress=true;
    
  case_number++;
  std::cout 
   << "Doing Mooney Rivlin with cont pressure formulation: Case "
   << case_number << std::endl; 

    
#ifdef REFINE
  {
   //Set up the problem with continous pressure/displacement
   CompressedSquareProblem<MySolidElement<
    RefineableQPVDElementWithContinuousPressure<2> > > 
    problem(incompress); 
     
   problem.run_it(case_number);
  }
#else
  {
   //Set up the problem with continous pressure/displacement
   CompressedSquareProblem<MySolidElement<
    QPVDElementWithContinuousPressure<2> > > 
    problem(incompress); 
     
   problem.run_it(case_number);
  }
#endif
    
    
  case_number++;
  std::cout 
   << "Doing Mooney Rivlin with discont pressure formulation: Case "
   << case_number << std::endl; 

    
#ifdef REFINE
  {
   //Set up the problem with discontinous pressure/displacement
   CompressedSquareProblem<MySolidElement<
    RefineableQPVDElementWithPressure<2> > > 
    problem(incompress); 
     
   problem.run_it(case_number);
  }
#else
  {
   //Set up the problem with discontinous pressure/displacement
   CompressedSquareProblem<MySolidElement<
    QPVDElementWithPressure<2> > > 
    problem(incompress); 
     
   problem.run_it(case_number);
  }
#endif
    
  delete  Global_Physical_Variables::Strain_energy_function_pt;
  Global_Physical_Variables::Strain_energy_function_pt=0;
    
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;
 }
   
  

} //end of main









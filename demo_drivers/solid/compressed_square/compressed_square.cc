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
//Driver for compressed square

//Oomph-lib includes
#include "generic.h"
#include "solid.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

namespace oomph
{

//=================start_wrapper==================================
/// Wrapper class for solid element to modify the output 
//================================================================
template <class ELEMENT>
class MySolidElement : public virtual ELEMENT
{

public:

 /// Constructor: Call constructor of underlying element
 MySolidElement() : ELEMENT() {}

 /// Overload output function
 void output(std::ostream &outfile, const unsigned &n_plot)
  {
   Vector<double> s(2);
   Vector<double> x(2);
   Vector<double> xi(2);
   DenseMatrix<double> sigma(2,2);
   
   //Tecplot header info 
   outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
   
   //Loop over plot points
   for(unsigned l2=0;l2<n_plot;l2++)
    {
     s[1] = -1.0 + l2*2.0/(n_plot-1);
     for(unsigned l1=0;l1<n_plot;l1++)
      {
       s[0] = -1.0 + l1*2.0/(n_plot-1);
       
       // Get Eulerian and Lagrangian coordinates and the stress
       this->interpolated_x(s,x);
       this->interpolated_xi(s,xi);
       this->get_stress(s,sigma);
       
       //Output the x,y coordinates
       for(unsigned i=0;i<2;i++) 
        {outfile << x[i] << " ";}
       
       // Output displacements, the difference between Eulerian and Lagrangian
       // coordinates
       for(unsigned i=0;i<2;i++) 
        {outfile << x[i]-xi[i] << " ";}
       
       //Output stress
       outfile << sigma(0,0) << " "
               << sigma(1,0) << " "
               << sigma(1,1) << " "
               << std::endl;
      }
    }  
  }
};

} //end namespace extension



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Poisson's ratio for Hooke's law
 double Nu=0.45;

 /// Pointer to strain energy function 
 StrainEnergyFunction* Strain_energy_function_pt=0;

 /// First "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C1=1.3;

 /// Second "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C2=1.3;

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
/// Problem class
//====================================================================== 
template<class ELEMENT>
class CompressedSquareProblem : public Problem
{

public:

 /// Constructor: Pass flag that determines if we want to use
 /// a true incompressible formulation
 CompressedSquareProblem(const bool& incompress);
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Doc the solution & exact (linear) solution for compressible 
 /// or incompressible materials
 void doc_solution(const bool& incompress);

 /// Run the job -- doc in RESLTi_case
 void run_it(const int& i_case,const bool& incompress);

private:

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: Pass flag that determines if we want to enforce
/// incompressibility
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
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;
 
 //Now create the mesh 
 mesh_pt() = new ElasticRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 //Assign the physical properties to the elements
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


   // Is the element based on the pressure/displacement formulation?
   PVDEquationsWithPressure<2>* test_pt = 
    dynamic_cast<PVDEquationsWithPressure<2>*>(mesh_pt()->element_pt(i));
   if (test_pt!=0)
    {
     // Do we want true incompressibility (formulation III in the 
     // associated tutorial) or not (formulation II)
     if (incompress)
      {     
       test_pt->set_incompressible();
      }
     else
      {
       // Note that this assignment isn't strictly necessary as it's the
       // default setting, but it doesn't do any harm to be explicit.
       test_pt->set_compressible();
      }
    }
  } // end compressibility
 
 
 // Choose a control node: The last node in the solid mesh
 unsigned nnod=mesh_pt()->nnode();
 Trace_node_pt=mesh_pt()->node_pt(nnod-1);
 
// Pin the left and right boundaries (1 and 2) in the horizontal directions
 for (unsigned b=1;b<4;b+=2)
  {
   unsigned nnod = mesh_pt()->nboundary_node(b);
   for(unsigned i=0;i<nnod;i++)
    {
     dynamic_cast<SolidNode*>(
      mesh_pt()->boundary_node_pt(b,i))->pin_position(0);
    }
  }
 
// Pin the bottom boundary (0) in the vertical direction
 unsigned b=0;
 {
  unsigned nnod= mesh_pt()->nboundary_node(b);
  for(unsigned i=0;i<nnod;i++)
   {
    dynamic_cast<SolidNode*>(
     mesh_pt()->boundary_node_pt(b,i))->pin_position(1);
   }
 }
 
 //Assign equation numbers
 assign_eqn_numbers();

} //end of constructor



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void CompressedSquareProblem<ELEMENT>::doc_solution(const bool& incompress)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of and stress in deformed body
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


 // Output exact solution for linear elasticity
 // -------------------------------------------
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nelem=mesh_pt()->nelement(); 
 Vector<double> s(2);
 Vector<double> x(2);
 DenseMatrix<double> sigma(2,2);
 
 // Poisson's ratio
 double nu=Global_Physical_Variables::Nu;
 if (incompress) nu=0.5;

 // Loop over all elements
 for (unsigned e=0;e<nelem;e++)
  {  
   //Cast to a solid element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
    
   //Tecplot header info 
   some_file << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
   
   //Loop over plot points
   for(unsigned l2=0;l2<n_plot;l2++)
    {
     s[1] = -1.0 + l2*2.0/(n_plot-1);
     for(unsigned l1=0;l1<n_plot;l1++)
      {
       s[0] = -1.0 + l1*2.0/(n_plot-1);
       
       // Get Lagrangian coordinates
         el_pt->interpolated_x(s,x);

         // Output the x,y,..
         for(unsigned i=0;i<2;i++) 
          {some_file << x[i] << " ";}

         // Exact vertical displacement
         double v_exact=Global_Physical_Variables::Gravity*
          (1.0+nu)*(1.0-2*nu)/(1.0-nu)*(0.5*x[1]*x[1]-x[1]);
         
         // x and y displacement
         some_file << "0.0 " <<   v_exact << " ";

         // Stresses
         sigma(0,0)=nu/(1.0-nu)*Global_Physical_Variables::Gravity*(x[1]-1.0);
         sigma(1,0)=0.0;
         sigma(1,1)=Global_Physical_Variables::Gravity*(x[1]-1.0);
          
         // Output linear stress tensor
         some_file << sigma(0,0) << " "
                   << sigma(1,0) << " "
                   << sigma(1,1) << " "
                   << std::endl;
      }
    }
  }
 some_file.close();

 // Increment label for output files
 Doc_info.number()++;

} //end doc





//==============start_run_it========================================
/// Run it
//==================================================================
template<class ELEMENT>
void CompressedSquareProblem<ELEMENT>::run_it(const int& i_case,
                                              const bool& incompress)
{

 // Set output directory
 char dirname[100];   
 sprintf(dirname,"RESLT%i",i_case);
 Doc_info.set_directory(dirname);

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);


 // Doc initial configuration
 doc_solution(incompress);

 // Initial values for parameter values
 Global_Physical_Variables::Gravity=0.0;
 
 //Parameter incrementation
 unsigned nstep=1; 
 double g_increment=1.0e-2;   
 for(unsigned i=0;i<nstep;i++)
  {
   // Bump up gravity
   Global_Physical_Variables::Gravity+=g_increment;
   
   // Solve
   newton_solve();

   // Doc solution
   doc_solution(incompress);
  }

}


//=======start_of_main==================================================
/// Driver for compressed square
//======================================================================
int main()
{
 
 //Flag to indicate if we want the material to be incompressible 
 bool incompress=false;

 // Label for different cases
 int case_number=-1;
 
 // Generalised Hookean constitutive equations
 //===========================================
 {
  // Nu values
  Vector<double> nu_value(3);
  nu_value[0]=0.45;
  nu_value[1]=0.499999;
  nu_value[2]=0.5;
  
  // Loop over nu values
  for (unsigned i=0;i<3;i++)
   {
    // Set Poisson ratio
    Global_Physical_Variables::Nu=nu_value[i];
  
    std::cout << "===========================================\n";
    std::cout << "Doing Nu=" <<  Global_Physical_Variables::Nu 
              << std::endl;
    std::cout << "===========================================\n";

    // Create constitutive equation
    Global_Physical_Variables::Constitutive_law_pt = 
     new GeneralisedHookean(&Global_Physical_Variables::Nu);
    
    // Displacement-based formulation
    //--------------------------------
    // (not for nu=1/2, obviously)
    if (i!=2) 
     {
      case_number++;
      std::cout 
       << "Doing Generalised Hookean with displacement formulation: Case "
       << case_number << std::endl;
      
      //Set up the problem with pure displacement based elements
      incompress=false;
      CompressedSquareProblem<MySolidElement<QPVDElement<2,3> > > 
       problem(incompress);
      
      // Run it
      problem.run_it(case_number,incompress);
     }
    

    // Generalised Hookean with continuous pressure, compressible
    //-----------------------------------------------------------
    {
     case_number++;
     std::cout 
      << "Doing Generalised Hookean with displacement/cont pressure "
      << "formulation: Case " << case_number << std::endl;     
     
     //Set up the problem with continous pressure/displacement
     incompress=false;
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     // Run it
     problem.run_it(case_number,incompress);
    }
     
    
    // Generalised Hookean with discontinuous pressure, compressible
    //--------------------------------------------------------------
    {
     case_number++;
     std::cout 
      << "Doing Generalised Hookean with displacement/discont pressure "
      << "formulation: Case " << case_number << std::endl; 

     //Set up the problem with discontinous pressure/displacement
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithPressure<2> > > problem(incompress); 
     
     // Run it
     problem.run_it(case_number,incompress);
    }



    // Generalised Hookean with continous pressure/displacement; 
    //----------------------------------------------------------
    // incompressible
    //---------------
    {
     case_number++;
     std::cout 
      << "Doing Generalised Hookean with displacement/cont pressure, "
      << "incompressible formulation: Case " << case_number << std::endl;
     
     //Set up the problem with continous pressure/displacement
     incompress=true;
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress); 
     
     // Run it
     problem.run_it(case_number,incompress);
    }


    // Generalised Hookean with discontinous pressure/displacement; 
    //-------------------------------------------------------------
    // incompressible
    //---------------
    {
     case_number++;
     std::cout 
      << "Doing Generalised Hookean with displacement/discont pressure, "
      << "incompressible formulation: Case " << case_number << std::endl; 
     
     //Set up the problem with discontinous pressure/displacement
     incompress=true;
     CompressedSquareProblem<MySolidElement<
      QPVDElementWithPressure<2> > > problem(incompress); 
     
     // Run it
     problem.run_it(case_number,incompress);
    }

    // Clean up
    delete Global_Physical_Variables::Constitutive_law_pt;
    Global_Physical_Variables::Constitutive_law_pt=0;
    
   } 
 
 } // end generalised Hookean

 
 // Incompressible Mooney Rivlin
 //=============================
 {

  // Create strain energy function
  Global_Physical_Variables::Strain_energy_function_pt = 
   new MooneyRivlin(&Global_Physical_Variables::C1,
                    &Global_Physical_Variables::C2);
  
  // Define a constitutive law (based on strain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
   new IsotropicStrainEnergyFunctionConstitutiveLaw(
    Global_Physical_Variables::Strain_energy_function_pt);
    

  // Mooney-Rivlin with continous pressure/displacement; 
  //----------------------------------------------------
  // incompressible
  //---------------
  {
   case_number++;
   std::cout 
    << "Doing Mooney Rivlin with cont pressure formulation: Case "
    << case_number << std::endl; 
   
   //Set up the problem with continous pressure/displacement
   incompress=true;
   CompressedSquareProblem<MySolidElement<
    QPVDElementWithContinuousPressure<2> > > 
    problem(incompress); 
     
   // Run it
   problem.run_it(case_number,incompress);
  }
    

  // Mooney-Rivlin with discontinous pressure/displacement; 
  //-------------------------------------------------------
  // incompressible
  //---------------
  {
   case_number++;
   std::cout 
    << "Doing Mooney Rivlin with discont pressure formulation: Case "
    << case_number << std::endl; 
   
   //Set up the problem with discontinous pressure/displacement
   incompress=true;
   CompressedSquareProblem<MySolidElement<
    QPVDElementWithPressure<2> > > 
    problem(incompress); 
     
   // Run it
   problem.run_it(case_number,incompress);
  }


  // Clean up
  delete  Global_Physical_Variables::Strain_energy_function_pt;
  Global_Physical_Variables::Strain_energy_function_pt=0;
  
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;
  
 }
   
} //end of main









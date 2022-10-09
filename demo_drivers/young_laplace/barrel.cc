//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Driver code for a 2D YoungLaplace problem

// Generic routines
#include "generic.h"

// The YoungLaplace equations
#include "young_laplace.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;
using namespace oomph;

//===== start_of_namespace========================================
/// Namespace for "global" problem parameters
//================================================================
namespace GlobalParameters
{

 // Displacement control:
 //----------------------

 /// Height control value
 double Controlled_height = 0.0;

 /// Exact kappa 
 double get_exact_kappa()
 { 

  // Mean curvature of barrel-shaped meniscus
  return 2.0*Controlled_height/
   (Controlled_height*Controlled_height+1.0/4.0);

 } //end exact kappa


 // Spine basis
 //------------

 /// Spine basis: The position vector to the basis of the spine
 /// as a function of the two coordinates x_1 and x_2, and its
 /// derivatives w.r.t. to these coordinates. 
 /// dspine_B[i][j] = d spine_B[j] / dx_i
 /// Spines start in the (x_1,x_2) plane at (x_1,x_2).
 void spine_base_function(const Vector<double>& x, 
                          Vector<double>& spine_B, 
                          Vector< Vector<double> >& dspine_B)
 {
  
  // Bspines and derivatives 
  spine_B[0]     = x[0];
  spine_B[1]     = x[1];
  spine_B[2]     = 0.0 ;
  dspine_B[0][0] = 1.0 ;
  dspine_B[1][0] = 0.0 ;
  dspine_B[0][1] = 0.0 ; 
  dspine_B[1][1] = 1.0 ;
  dspine_B[0][2] = 0.0 ;
  dspine_B[1][2] = 0.0 ;
  
 } // End of bspine functions
 
 

 // Spines rotate in the y-direction
 //---------------------------------

 /// Min. spine angle against horizontal plane
 double Alpha_min = MathematicalConstants::Pi/2.0*1.5;

 /// Max. spine angle against horizontal plane
 double Alpha_max = MathematicalConstants::Pi/2.0*0.5;

 /// Spine: The spine vector field as a function of the two 
 /// coordinates x_1 and x_2, and its derivatives w.r.t. to these coordinates:
 /// dspine[i][j] = d spine[j] / dx_i
 void spine_function(const Vector<double>& x, 
                     Vector<double>& spine, 
                     Vector< Vector<double> >& dspine)
 {
  
  /// Spines (and derivatives)  are independent of x[0] and rotate 
  /// in the x[1]-direction
  spine[0]=0.0;
  dspine[0][0]=0.0; 
  dspine[1][0]=0.0; 
  
  spine[1]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]); 
  dspine[0][1]=0.0;                                   
  dspine[1][1]=-sin(Alpha_min+(Alpha_max-Alpha_min)*x[1])
   *(Alpha_max-Alpha_min);            
  
  spine[2]=sin(Alpha_min+(Alpha_max-Alpha_min)*x[1]);
  dspine[0][2]=0.0;                                  
  dspine[1][2]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]) 
   *(Alpha_max-Alpha_min);            

 } // End spine function


} // end of namespace




//====== start_of_problem_class=======================================
/// 2D YoungLaplace problem on rectangular domain, discretised with
/// 2D QYoungLaplace elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class YoungLaplaceProblem : public Problem
{

public:

 /// Constructor: 
 YoungLaplaceProblem();

 /// Destructor (empty)
 ~YoungLaplaceProblem(){}

 /// Update the problem before solve 
 void actions_before_newton_solve()
  {
   // This only has an effect if displacement control is disabled
   double new_kappa=Kappa_pt->value(0)-0.5;
   Kappa_pt->set_value(0,new_kappa);
  }

 /// Update the problem after solve: Empty
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to and the trace file
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

private:

 /// Node at which the height (displacement along spine) is controlled/doced
 Node* Control_node_pt;

 /// Pointer to Data object that stores the prescribed curvature
 Data* Kappa_pt;

}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for YoungLaplace problem
//========================================================================
template<class ELEMENT>
YoungLaplaceProblem<ELEMENT>::YoungLaplaceProblem()
{ 

 // Setup mesh
 //-----------

 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;
 
 // Build and assign mesh
 Problem::mesh_pt()=new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 // Check that we've got an even number of elements otherwise
 // out counting doesn't work...
 if ((n_x%2!=0)||(n_y%2!=0))
  {
   cout << "n_x n_y should be even" << endl;
   abort();
  }
  
 //  This is the element that contains the central node:
 ELEMENT* prescribed_height_element_pt= dynamic_cast<ELEMENT*>(
  mesh_pt()->element_pt(n_y*n_x/2+n_x/2));
 
 // The central node is node 0 in that element
 Control_node_pt= static_cast<Node*>(prescribed_height_element_pt->node_pt(0));

 std::cout << "Controlling height at (x,y) : (" << Control_node_pt->x(0) 
           << "," << Control_node_pt->x(1)  << ")" << "\n" << endl;


 // Create a height control element
 HeightControlElement* height_control_element_pt=new HeightControlElement(
  Control_node_pt,&GlobalParameters::Controlled_height);
 
 // Store pointer to kappa data
 Kappa_pt=height_control_element_pt->kappa_pt();


 // Comment out the previous two commands and uncomment the following two
 // to prescribe the pressure drop (the curvature) directly
 //Kappa_pt=new Data(1);
 //Kappa_pt->pin(0);


 // Boundary conditions
 //--------------------

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = mesh_pt()->nboundary(); 
 for(unsigned b=0;b<n_bound;b++)
  {

   // Pin meniscus displacement at all nodes boundaries  0 and 2
   if ((b==0)||(b==2))
    {
     unsigned n_node = mesh_pt()->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); 
      }
    }

  } // end bc
 
 // Complete build of elements
 //---------------------------

 // Complete the build of all elements so they are fully functional 
 unsigned nelement = mesh_pt()->nelement();
 for(unsigned i=0;i<nelement;i++)
  {
   // Upcast from GeneralsedElement to YoungLaplace element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the spine function pointers
   el_pt->spine_base_fct_pt() = GlobalParameters::spine_base_function;
   el_pt->spine_fct_pt() =  GlobalParameters::spine_function;
  
   // Set the curvature data for the element
   el_pt->set_kappa(Kappa_pt); 
  }

 // Add the height control element to mesh (comment this out
 // if you're not using displacement control)
 mesh_pt()->add_element_pt(height_control_element_pt);
 
 // Setup equation numbering scheme
 cout <<"\nNumber of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor




//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void YoungLaplaceProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                                ofstream& trace_file)
{ 

 // Output kappa vs height of the apex 
 //------------------------------------
 trace_file << -1.0*Kappa_pt->value(0) << " ";
 trace_file << GlobalParameters::get_exact_kappa()  << " ";
 trace_file << Control_node_pt->value(0) ;
 trace_file << endl; 
 
 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output full solution 
 //---------------------
 ofstream some_file;
 char filename[100];
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc


//===================start_of_main========================================
/// Driver code
//========================================================================
int main()
{
 
 // Create label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 //Open a trace file
 ofstream trace_file;
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 // Write kappa, exact kappa and height values
 trace_file
  << "VARIABLES=\"<GREEK>k</GREEK>\",\"<GREEK>k</GREEK>_{ex}\",\"h\"" 
  << std::endl;
 trace_file << "ZONE" << std::endl;
  
 // Create the problem
 //-------------------
 
 // Create the problem with 2D nine-node elements from the
 // QYoungLaplaceElement family. 
 YoungLaplaceProblem<QYoungLaplaceElement<3> > problem;

 //Output the solution
 problem.doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
 doc_info.number()++;

 
 // Parameter incrementation
 //------------------------- 
 double increment=0.1;

 // Loop over steps
 unsigned nstep=2; // 10;
 for (unsigned istep=0;istep<nstep;istep++)
  {
   
   // Increment prescribed height value
   GlobalParameters::Controlled_height+=increment;
   
   // Solve the problem
   problem.newton_solve();
   
   //Output the solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
   
  }
 
 // Close output file
 trace_file.close();

} // end of main







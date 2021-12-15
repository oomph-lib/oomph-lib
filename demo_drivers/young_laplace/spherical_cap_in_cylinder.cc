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
//Driver for refineable 2D Young Laplace problem on a circle sector

//Generic routines
#include "generic.h"

// The YoungLaplace equations
#include "young_laplace.h"


// The mesh
#include "meshes/quarter_circle_sector_mesh.h"


using namespace std;
using namespace oomph;


// Namespace (shared with other driver codes)
#include "common_young_laplace_stuff.h"


//====== start_of_problem_class=======================================
/// 2D RefineableYoungLaplace problem on a circle sector, discretised with
/// 2D QRefineableYoungLaplace elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableYoungLaplaceProblem : public Problem
{

public:

 /// Constructor: 
 RefineableYoungLaplaceProblem();

 /// Destructor (empty)
 ~RefineableYoungLaplaceProblem(){};

 /// Update the problem specs before solve: Empty
 void actions_before_solve(){};

 /// Update the problem after solve: Empty
 void actions_after_solve(){};

 /// Increment problem parameters
 void increment_parameters();

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to and the trace file
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

private:

 /// Pointer to GeomObject that specifies the domain bondary
 GeomObject* Boundary_pt;

 /// Pointer to the "bulk" mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to mesh containing the height control element
 Mesh* Height_control_mesh_pt;

 /// Pointer to height control element
 HeightControlElement* Height_control_element_pt;

 /// Node at which the height (displacement along spine) is controlled/doced
 Node* Control_node_pt;

}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for RefineableYoungLaplace problem
//========================================================================
template<class ELEMENT>
RefineableYoungLaplaceProblem<ELEMENT>::RefineableYoungLaplaceProblem()
{ 

 // Setup dependent parameters in namespace
 GlobalParameters::setup_dependent_parameters_and_sanity_check();

 // Setup bulk mesh
 //----------------

 // Build geometric object that forms the curvilinear domain boundary:
 // a unit circle 

 // Create GeomObject
 Boundary_pt=new Circle(0.0,0.0,1.0);

 // Start and end coordinates of curvilinear domain boundary on circle  
 double xi_lo=0.0;
 double xi_hi=MathematicalConstants::Pi/2.0;
 
 // Now create the bulk mesh. Separating line between the two 
 // elements next to the curvilinear boundary is located half-way
 // along the boundary.
 double fract_mid=0.5;
 Bulk_mesh_pt = new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Boundary_pt,xi_lo,fract_mid,xi_hi);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=1.0e-4;
 Bulk_mesh_pt->min_permitted_error()=1.0e-6;

 // Add bulk mesh to the global mesh
 add_sub_mesh(Bulk_mesh_pt);


 // Prescribed height?
 //-------------------


 // Which element are we using for displacement control?
 GlobalParameters::Control_element=0;

 //  Choose the prescribed height element
 ELEMENT* prescribed_height_element_pt= dynamic_cast<ELEMENT*>(
  Bulk_mesh_pt->element_pt(GlobalParameters::Control_element));
 
 // ...and the associated control node (node 0 in that element)
 // (we're storing this node even if there's no height-control, for
 // output purposes...)
 Control_node_pt= static_cast<Node*>(
  prescribed_height_element_pt->node_pt(0));

 cout << "Controlling height at (x,y) : (" << Control_node_pt->x(0) 
      << "," << Control_node_pt->x(1)  << ")" << endl;

 // If needed, create a height control element and store the
 // pointer to the Kappa Data created by this object
 Height_control_element_pt=0;
 Height_control_mesh_pt=0;
 if (GlobalParameters::Use_height_control)
  {
   Height_control_element_pt=new HeightControlElement(
    Control_node_pt,&GlobalParameters::Controlled_height);

   GlobalParameters::Kappa_pt=Height_control_element_pt->kappa_pt();

   // Add to mesh
   Height_control_mesh_pt = new Mesh;
   Height_control_mesh_pt->add_element_pt(Height_control_element_pt);

   // Add height control mesh to the global mesh
   add_sub_mesh(Height_control_mesh_pt);

  }
 //...otherwise create a kappa data item from scratch and pin its 
 // single unknown value
 else
  {
   GlobalParameters::Kappa_pt=new Data(1);
   GlobalParameters::Kappa_pt->set_value(0,GlobalParameters::Kappa_initial);
   GlobalParameters::Kappa_pt->pin(0);
  }
  
 // Build global mesh
 //------------------
 build_global_mesh();


 // Boundary conditions
 //--------------------

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_node = Bulk_mesh_pt->nboundary_node(1);
 for (unsigned n=0;n<n_node;n++)
  {
   Bulk_mesh_pt->boundary_node_pt(1,n)->pin(0); 
  }
 
 
 // Complete build of elements
 //---------------------------
 
 // Complete the build of all elements so they are fully functional 
 unsigned n_bulk=Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_bulk;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
   
   if ( GlobalParameters::Use_spines )
    {
     //Set the spine function pointers
     el_pt->spine_base_fct_pt() = GlobalParameters::spine_base_function;
     el_pt->spine_fct_pt() =  GlobalParameters::spine_function;
    }
   
   // Set the curvature data for the element
   el_pt->set_kappa(GlobalParameters::Kappa_pt);  
   
  }
 
 // Setup equation numbering scheme
 cout <<"\nNumber of equations: " << assign_eqn_numbers() << endl; 
 cout << "\n********************************************\n" <<  endl;
 
} // end of constructor



//===============start_of_update_parameters==============================
/// Update (increase/decrease) parameters
//=======================================================================
template<class ELEMENT>
void RefineableYoungLaplaceProblem<ELEMENT>::increment_parameters()
{ 
 
 GlobalParameters::Controlled_height+=
  GlobalParameters::Controlled_height_increment;
 
 cout << "Solving for Prescribed Height Value = " ;
 cout << GlobalParameters::Controlled_height << "\n" << endl;
 
}



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void RefineableYoungLaplaceProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                              ofstream& trace_file)
{ 

 // Output kappa vs height
 //-----------------------
 trace_file << -1.0*GlobalParameters::Kappa_pt->value(0) << " "; 
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
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

}



//===== start_of_main=====================================================
/// Driver code for 2D RefineableYoungLaplace problem. Input arguments: none
/// (for validation) or number of steps.
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // No command line args: Running with limited number of steps
 if (CommandLineArgs::Argc==1)
  {
   std::cout 
    << "Running with limited number of steps for validation" 
    << std::endl;

   // Number of steps
   GlobalParameters::Nsteps=2;
  } 
 else 
  {
   // Number of steps
   GlobalParameters::Nsteps=atoi(argv[1]);
  }
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set outputs
 //------------ 
 
 // Trace file
 ofstream trace_file;

// Set output directory
 doc_info.set_directory("RESLT_adapt_pinned_spherical_cap_in_cylinder");
 
//Open a trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 trace_file 
  << "VARIABLES=\"<GREEK>k</GREEK>\",\"<GREEK>k</GREEK>_{exact}\",\"h\""
  << std::endl;
 trace_file << "ZONE" << std::endl;
 

 // Set case
 GlobalParameters::Case=GlobalParameters::Spherical_cap_in_cylinder_pinned;
 
 // Run with spines
 GlobalParameters::Use_spines=true; 


 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // RefineableQYoungLaplaceElement family. 
 RefineableYoungLaplaceProblem<RefineableQYoungLaplaceElement<3> > problem;
 
 //Output the solution
 problem.doc_solution(doc_info,trace_file);
  
 //Increment counter for solutions 
 doc_info.number()++;
 
 // Parameter incrementation
 //------------------------- 
 
 // Loop over steps
 for (unsigned istep=0;istep<GlobalParameters::Nsteps;istep++)
  {
   // Solve the problem
   unsigned max_adapt=1;
   problem.newton_solve(max_adapt);
   
   //Output the solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
   
   // Increase the parameters
   problem.increment_parameters();
  }
 
 // Close output file
 trace_file.close();
 






} //end of main



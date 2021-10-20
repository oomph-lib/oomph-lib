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
// Driver code for a 2D YoungLaplace problem

// Generic routines
#include "generic.h"

// The YoungLaplace equations
#include "young_laplace.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;
using namespace oomph;


// Namespace (shared with refineable version)
#include "common_young_laplace_stuff.h"



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
 ~YoungLaplaceProblem(){};

 ///  Update the problem specs before solve
 void actions_before_newton_solve();

 /// Update the problem after solve: Empty
 void actions_after_newton_solve(){};

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to and the trace file
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

private:

 ///  Create YoungLaplace contact angle elements on the 
 /// b-th boundary of the problem's mesh and add them to mesh
 void create_contact_angle_elements(const unsigned& b);
 
 ///  Number of YoungLaplace "bulk" elements (We're attaching the  
 /// contact angle elements to the bulk mesh --> only the first 
 /// N_bulk_elements elements in the mesh are bulk elements!)
 unsigned N_bulk_elements;

 /// Number of last FaceElement on boundary 1
  unsigned Last_element_on_boundary1; 
 
 /// Number of last FaceElement on boundary 3
  unsigned Last_element_on_boundary3; 

 /// Node at which the height (displacement along spine) is controlled/doced
 Node* Control_node_pt;

}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for YoungLaplace problem
//========================================================================
template<class ELEMENT>
YoungLaplaceProblem<ELEMENT>::YoungLaplaceProblem()
{ 

 // Setup dependent parameters in namespace
 GlobalParameters::setup_dependent_parameters_and_sanity_check();
 
 // Setup mesh
 //-----------

 // # of elements in x-direction
 unsigned n_x=GlobalParameters::N_x;

 // # of elements in y-direction
 unsigned n_y=GlobalParameters::N_y;

 // Domain length in x-direction
 double l_x=GlobalParameters::L_x;

 // Domain length in y-direction
 double l_y=GlobalParameters::L_y;
 
 // Print Size of the mesh
 cout <<  "Lx = " << l_x << " and Ly = " << l_y << endl;

 // Build and assign mesh
 Problem::mesh_pt()=new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 
 //  Choose the prescribed height element
 ELEMENT* prescribed_height_element_pt= dynamic_cast<ELEMENT*>(
  mesh_pt()->element_pt(GlobalParameters::Control_element));
 
 // ...and the associated control node (node 0 in that element)
 // (we're storing this node even if there's no height-control, for
 // output purposes...)
 Control_node_pt= static_cast<Node*>(
  prescribed_height_element_pt->node_pt(0));

 // If needed, create a height control element and store the
 // pointer to the Kappa Data created by this object
 HeightControlElement* height_control_element_pt=0;
 if (GlobalParameters::Use_height_control)
  {
   cout << "Controlling height at (x,y) : (" << Control_node_pt->x(0) 
	<< "," << Control_node_pt->x(1)  << ")" << "\n" << endl;

   height_control_element_pt=new HeightControlElement(
	 Control_node_pt,&GlobalParameters::Controlled_height);

   GlobalParameters::Kappa_pt=height_control_element_pt->kappa_pt();

   // Assign initial value for kappa value
   height_control_element_pt->kappa_pt()
     ->set_value(0,GlobalParameters::Kappa_initial);
  }
 //...otherwise create a kappa data item from scratch and pin its 
 // single unknown value
 else
  {
   GlobalParameters::Kappa_pt=new Data(1);
   GlobalParameters::Kappa_pt->set_value(0,GlobalParameters::Kappa_initial);
   GlobalParameters::Kappa_pt->pin(0);
  }
  
 // Number of elements in the bulk mesh
 N_bulk_elements = mesh_pt()->nelement();
 
 // Create contact angle elements from all elements that are 
 // adjacent to boundary 1 and 3 and add them to the (single) global mesh
 if (GlobalParameters::Case==
     GlobalParameters::T_junction_with_nonzero_contact_angle)
   {
     create_contact_angle_elements(1);
     Last_element_on_boundary1=mesh_pt()->nelement();
     
     create_contact_angle_elements(3);
     Last_element_on_boundary3=mesh_pt()->nelement();
   }
 else
   {
     Last_element_on_boundary1=0;
     Last_element_on_boundary3=0;
   }

 // Boundary conditions
 //--------------------

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = mesh_pt()->nboundary(); 
 for(unsigned b=0;b<n_bound;b++)
  {
   // Pin all boundaries for two cases and only boundaries
   // 0 and 2 in all others:
   if ((GlobalParameters::Case==GlobalParameters::All_pinned)||
       (b==0)||
       (b==2))
    {
     unsigned n_node = mesh_pt()->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); 
      }
    }
  }

 // Complete build of elements
 //---------------------------

 // Complete the build of all elements so they are fully functional 
 for(unsigned i=0;i<N_bulk_elements;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   if ( GlobalParameters::Use_spines )
     {
       //Set the spine function pointers
       el_pt->spine_base_fct_pt() = GlobalParameters::spine_base_function;
       el_pt->spine_fct_pt() =  GlobalParameters::spine_function;
     }

   // Set the curvature data for the element
   el_pt->set_kappa(GlobalParameters::Kappa_pt);

  }

 // Set function pointers for contact-angle elements
 if (GlobalParameters::Case==
     GlobalParameters::T_junction_with_nonzero_contact_angle)
  {
   // Loop over the contact-angle elements (located at the "end" of the
   // mesh) to pass function pointer to contact angle function.
   for (unsigned e=N_bulk_elements;e<Last_element_on_boundary1;e++)
    {
     // Upcast from GeneralisedElement to YoungLaplace contact angle element
     YoungLaplaceContactAngleElement<ELEMENT>* el_pt = 
      dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
	  mesh_pt()->element_pt(e));
     
     // Set the pointer to the prescribed contact angle
     el_pt->prescribed_cos_gamma_pt() = &GlobalParameters::Cos_gamma;
    }
   for (unsigned e=Last_element_on_boundary1;e<Last_element_on_boundary3;e++)
    {
     // Upcast from GeneralisedElement to YoungLaplace contact angle element
     YoungLaplaceContactAngleElement<ELEMENT> *el_pt = 
      dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
       mesh_pt()->element_pt(e));
     
     // Set the pointer to the prescribed contact angle
     el_pt->prescribed_cos_gamma_pt() = &GlobalParameters::Cos_gamma;
    }
  }

 // Height Element setup
 //---------------------

 /// Add height control element to mesh at the very end
 if (GlobalParameters::Use_height_control)
  {
   mesh_pt()->add_element_pt(height_control_element_pt);
  }


 // Setup equation numbering scheme
 cout <<"\nNumber of equations: " << assign_eqn_numbers() << endl; 
 cout << "\n********************************************\n" <<  endl;

} // end of constructor



//============start_of_create_contact_angle_elements=====================
/// Create YoungLaplace contact angle elements on the b-th boundary of the Mesh
//=======================================================================
template<class ELEMENT>
void YoungLaplaceProblem<ELEMENT>::create_contact_angle_elements(
 const unsigned &b)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = mesh_pt()->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    mesh_pt()->boundary_element_pt(b,e));
   
   // What is the index of the face of the bulk element at the boundary
   int face_index = mesh_pt()->face_index_at_boundary(b,e);

   // Build the corresponding prescribed contact angle element
   YoungLaplaceContactAngleElement<ELEMENT>* contact_angle_element_pt = new 
   YoungLaplaceContactAngleElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed contact angle element to the mesh
   mesh_pt()->add_element_pt(contact_angle_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_contact_angle_elements




//========================================start_of_actions_before_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void YoungLaplaceProblem<ELEMENT>::actions_before_newton_solve()
{
  
 // Increment kappa or height value
 if (!GlobalParameters::Use_height_control)
  {
   double kappa=GlobalParameters::Kappa_pt->value(0);
   kappa+=GlobalParameters::Kappa_increment;
   GlobalParameters::Kappa_pt->set_value(0,kappa);
   
   cout << "Solving for Prescribed KAPPA Value = " ;
   cout << GlobalParameters::Kappa_pt->value(0) << "\n" << endl;
  }
 else
  {
   GlobalParameters::Controlled_height+=
    GlobalParameters::Controlled_height_increment;
   
   cout << "Solving for Prescribed HEIGHT Value = " ;
   cout << GlobalParameters::Controlled_height << "\n" << endl;
  }
 
}  // end of actions before solve




//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void YoungLaplaceProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                                ofstream& trace_file)
{ 

 // Output kappa vs height of the apex 
 //------------------------------------
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
 for(unsigned i=0;i<N_bulk_elements;i++)
  {
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   el_pt->output(some_file,npts);
  }
 some_file.close();

 // Output contact angle 
 //---------------------
 
 //Doc contact angle stuff
 if (GlobalParameters::Case==
     GlobalParameters::T_junction_with_nonzero_contact_angle)
  {
   
   ofstream tangent_file;
   sprintf(filename,"%s/tangent_to_contact_line%i.dat",
           doc_info.directory().c_str(),
           doc_info.number());
   tangent_file.open(filename);
   
   ofstream normal_file;
   sprintf(filename,"%s/normal_to_contact_line%i.dat",
           doc_info.directory().c_str(),
           doc_info.number());
   normal_file.open(filename);
   
   
   ofstream contact_angle_file;
   sprintf(filename,"%s/contact_angle%i.dat",
           doc_info.directory().c_str(),
           doc_info.number());
   contact_angle_file.open(filename);
   
   // Tangent and normal vectors to contact line
   Vector<double> tangent(3);
   Vector<double> normal(3);
   Vector<double> r_contact(3);
   
   
   // Loop over both boundaries
   for (unsigned bnd=0;bnd<2;bnd++)
    {
     unsigned el_lo=N_bulk_elements;
     unsigned el_hi=Last_element_on_boundary1;
     if (1==bnd)
      {
       el_lo=Last_element_on_boundary1;
       el_hi=Last_element_on_boundary3;
      }
     
     // Loop over the contact-angle elements (located at the "end" of the
     // mesh)
     for (unsigned e=el_lo;e<el_hi;e++)
      {
       
       tangent_file << "ZONE" << std::endl;
       normal_file << "ZONE" << std::endl;
       contact_angle_file << "ZONE" << std::endl;
       
       // Upcast from GeneralisedElement to YoungLaplace contact angle element
       YoungLaplaceContactAngleElement<ELEMENT>* el_pt = 
        dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
         mesh_pt()->element_pt(e));
       
       // Loop over a few points in the contact angle element
       Vector<double> s(1);
       for (unsigned i=0;i<npts;i++)
        {
         s[0]=-1.0+2.0*double(i)/double(npts-1);
         
         dynamic_cast<ELEMENT*>(el_pt->bulk_element_pt())->
          position(el_pt->local_coordinate_in_bulk(s),r_contact);
         
         el_pt->contact_line_vectors(s,tangent,normal);
         tangent_file << r_contact[0] << " " 
                      << r_contact[1] << " " 
                      << r_contact[2] << " " 
                      << tangent[0] << " " 
                      << tangent[1] << " " 
                      << tangent[2] << " "  << std::endl;
         
         normal_file << r_contact[0] << " " 
                     << r_contact[1] << " " 
                     << r_contact[2] << " " 
                     << normal[0] << " " 
                     << normal[1] << " " 
                     << normal[2] << " "  << std::endl;
        

         contact_angle_file << r_contact[1] << " " 
                            << el_pt->actual_cos_contact_angle(s)
                            << std::endl;
        }
      }
     
    } // end of loop over both boundaries
   
   tangent_file.close();
   normal_file.close();
   contact_angle_file.close();
   
  }

cout << "\n********************************************" << endl <<  endl;

} // end of doc

//========================================================================
/// Run code for current setting of parameter values -- specify name
/// of output directory
//========================================================================
void run_it(const string& output_directory)
{
 
 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set outputs
 //------------ 

 // Trace file
 ofstream trace_file;
 
 // Set output directory
 doc_info.set_directory(output_directory);
 
 //Open a trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 // Write kappa, exact kappa if exists and height values
 trace_file
  << "VARIABLES=\"<GREEK>k</GREEK>\",\"<GREEK>k</GREEK>_{ex}\",\"h\"" 
  << std::endl;
 trace_file << "ZONE" << std::endl;
  
 //Set up the problem
 //------------------
 
 // Create the problem with 2D nine-node elements from the
 // QYoungLaplaceElement family. 
 YoungLaplaceProblem<QYoungLaplaceElement<3> > problem;

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
   problem.newton_solve();
   
   //Output the solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
   
  }
 
 // Close output file
 trace_file.close();

} //end of run_it()



//===== start_of_main=====================================================
/// Driver code for 2D YoungLaplace problem. Input arguments: none
/// (for validation) or case (0,1,2 for all pinned, barrel and T junction)
/// and number of steps.
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Cases to run (By default (validation) run all
 unsigned case_lo=0;
 unsigned case_hi=2;

 // No command line args: Running every case with 
 // limited number of steps
 if (CommandLineArgs::Argc==1)
  {
   std::cout 
    << "Running every case with limited number of steps for validation" 
    << std::endl;
   
   // Number of steps
   GlobalParameters::Nsteps=2;
  }
 else 
  {
   // Which case to run?
   case_lo=atoi(argv[1]);
   case_hi=atoi(argv[1]);

   // Number of steps
   GlobalParameters::Nsteps=atoi(argv[2]);
  }


 // Loop over chosen case(s)
 //-------------------------
 for (unsigned my_case=case_lo;my_case<=case_hi;my_case++)
  {
   
   // Choose
   switch (my_case)
    {
     
    case 0:

     cout << endl << endl 
          << "//////////////////////////////////////////////////////////\n"
          << "All pinned solution \n" 
          << "//////////////////////////////////////////////////////////\n\n";
     
     GlobalParameters::Case=GlobalParameters::All_pinned;
      
     // Run with spines
     GlobalParameters::Use_spines=true; 
     run_it("RESLT_all_pinned");
       
     break;

    case 1:
          
     cout << endl << endl 
          << "//////////////////////////////////////////////////////////\n"
          << "Barrel-shaped solution \n" 
          << "/////////////////////////////////////////////////////////\n\n";
     
     GlobalParameters::Case= GlobalParameters::Barrel_shape;
      
     // Run with spines
     GlobalParameters::Use_spines=true; 
     run_it("RESLT_barrel_shape");
          
     break;
     
    case 2:
     
     cout << endl << endl 
          << "//////////////////////////////////////////////////////////\n"
          << "T-junction solution \n" 
          << "//////////////////////////////////////////////////////////\n\n";
     
     GlobalParameters::Case= 
      GlobalParameters::T_junction_with_nonzero_contact_angle;

     // Adjust spine orientation
     GlobalParameters::Controlled_height_increment=0.05;
     GlobalParameters::Gamma=MathematicalConstants::Pi/6.0;
 
     // Run with spines
     GlobalParameters::Use_spines=true;
     run_it("RESLT_T_junction");
         
     break;

    default:
     
     std::cout << "Wrong case! Options are:\n"
               << "0: All pinned\n" 
               << "1: Barrel \n"
               << "2: T_junction\n"
	       << std::endl;
     assert(false);
     
    }
   
  }
 
 
} //end of main







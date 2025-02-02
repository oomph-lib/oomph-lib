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
//Driver for refineable Young Laplace problem

//Generic routines
#include "generic.h"

// The YoungLaplace equations
#include "young_laplace.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

// Namespaces
using namespace std;
using namespace oomph;



//======start_of_namespace========================================
/// Namespace for "global" problem parameters
//================================================================
namespace GlobalParameters
{

 /// Cos of contact angle 
 double Cos_gamma=cos(MathematicalConstants::Pi/6.0);

 /// Height control value for displacement control
 double Controlled_height = 0.0;

 /// Length of domain
 double L_x = 1.0; 

 /// Width of domain
 double L_y = 5.0; 

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
  
  spine[1]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]/L_y); 
  dspine[0][1]=0.0;                                   
  dspine[1][1]=-sin(Alpha_min+(Alpha_max-Alpha_min)*x[1]/L_y)
   *(Alpha_max-Alpha_min)/L_y;            
  
  spine[2]=sin(Alpha_min+(Alpha_max-Alpha_min)*x[1]/L_y);
  dspine[0][2]=0.0;                                  
  dspine[1][2]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]/L_y) 
   *(Alpha_max-Alpha_min)/L_y;            

 } // End spine function


} // end of namespace






//====== start_of_problem_class=======================================
/// 2D RefineableYoungLaplace problem on rectangular domain, discretised with
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
 void actions_before_newton_solve(){};

 /// Update the problem after solve: Empty
 void actions_after_newton_solve(){};

 /// Actions before adapt: Wipe the mesh of contact angle elements
 void actions_before_adapt()
  {
   // Kill the contact angle elements and wipe contact angle mesh
   if (Contact_angle_mesh_pt!=0) delete_contact_angle_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }


 ///  Actions after adapt: Rebuild the mesh of contact angle elements
 void actions_after_adapt()
  {
   create_contact_angle_elements(1);
   create_contact_angle_elements(3);
   
   // Set function pointers for contact-angle elements
   unsigned nel=Contact_angle_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     // Upcast from GeneralisedElement to YoungLaplace contact angle
     // element
     YoungLaplaceContactAngleElement<ELEMENT> *el_pt = 
      dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
       Contact_angle_mesh_pt->element_pt(e));
     
     // Set the pointer to the prescribed contact angle
     el_pt->prescribed_cos_gamma_pt() = &GlobalParameters::Cos_gamma;
    }
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
   
  }

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to and the trace file
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

private:

 /// Create YoungLaplace contact angle elements on the 
 /// b-th boundary of the bulk mesh and add them to contact angle mesh
 void create_contact_angle_elements(const unsigned& b);
 
 /// Delete contact angle elements 
 void delete_contact_angle_elements();

 /// Pointer to the "bulk" mesh
 RefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the contact angle mesh
 Mesh* Contact_angle_mesh_pt;

 /// Pointer to mesh containing the height control element
 Mesh* Height_control_mesh_pt;

 /// Node at which the height (displacement along spine) is controlled/doced
 Node* Control_node_pt;

 /// Pointer to Data object that stores the prescribed curvature
 Data* Kappa_pt;

}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for RefineableYoungLaplace problem
//========================================================================
template<class ELEMENT>
RefineableYoungLaplaceProblem<ELEMENT>::RefineableYoungLaplaceProblem()
{ 

 // Setup bulk mesh
 //----------------

 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=GlobalParameters::L_x;

 // Domain length in y-direction
 double l_y=GlobalParameters::L_y;
 
 // Build and assign mesh
 Bulk_mesh_pt=new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=1.0e-4;
 Bulk_mesh_pt->min_permitted_error()=1.0e-6;

 // Check that we've got an even number of elements otherwise
 // out counting doesn't work...
 if ((n_x%2!=0)||(n_y%2!=0))
  {
   cout << "n_x n_y should be even" << endl;
   abort();
  }
  
 //  This is the element that contains the central node:
 ELEMENT* prescribed_height_element_pt= dynamic_cast<ELEMENT*>(
  Bulk_mesh_pt->element_pt(n_y*n_x/2+n_x/2));
 
 // The central node is node 0 in that element
 Control_node_pt= static_cast<Node*>(prescribed_height_element_pt->node_pt(0));

 std::cout << "Controlling height at (x,y) : (" << Control_node_pt->x(0) 
           << "," << Control_node_pt->x(1)  << ")" << "\n" << endl;

 // Create a height control element and store the
 // pointer to the Kappa Data created by this object
 HeightControlElement* height_control_element_pt=new HeightControlElement(
  Control_node_pt,&GlobalParameters::Controlled_height);
 
 // Add to mesh
 Height_control_mesh_pt = new Mesh;
 Height_control_mesh_pt->add_element_pt(height_control_element_pt);

 // Store curvature data
 Kappa_pt=height_control_element_pt->kappa_pt();


 // Contact angle elements
 //-----------------------

 // Create prescribed-contact-angle elements from all elements that are 
 // adjacent to boundary 1 and 3 and add them to their own mesh

 // set up new mesh
 Contact_angle_mesh_pt=new Mesh;
 
 // creation of contact angle elements
 create_contact_angle_elements(1);
 create_contact_angle_elements(3);


 // Add various meshes and build the global mesh
 //----------------------------------------------
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Height_control_mesh_pt);
 add_sub_mesh(Contact_angle_mesh_pt);
 build_global_mesh();


 // Boundary conditions
 //--------------------

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = Bulk_mesh_pt->nboundary(); 
 for(unsigned b=0;b<n_bound;b++)
  {
   // Pin all boundaries for three cases and only boundaries
   // 0 and 2 in all others:
   if ((b==0)||(b==2))
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    }
  } // end bcs
 
 // Complete build of elements
 //---------------------------

 // Complete the build of all elements so they are fully functional 
 unsigned n_bulk=Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_bulk;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   //Set the spine function pointers
   el_pt->spine_base_fct_pt() = GlobalParameters::spine_base_function;
   el_pt->spine_fct_pt() =  GlobalParameters::spine_function;
  
   // Set the curvature data for the element
   el_pt->set_kappa(Kappa_pt); 
  }
 
 // Set function pointers for contact-angle elements
 unsigned nel=Contact_angle_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   // Upcast from GeneralisedElement to YoungLaplace contact angle 
   // element
   YoungLaplaceContactAngleElement<ELEMENT> *el_pt = 
    dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
     Contact_angle_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed contact angle
   el_pt->prescribed_cos_gamma_pt() = &GlobalParameters::Cos_gamma;
  }


 // Setup equation numbering scheme
 cout <<"\nNumber of equations: " << assign_eqn_numbers() << endl; 
 cout << "\n********************************************\n" <<  endl;

} // end of constructor


//============start_of_create_contact_angle_elements=====================
/// Create YoungLaplace contact angle elements on the b-th boundary of the 
/// bulk mesh and add them to the contact angle mesh
//=======================================================================
template<class ELEMENT>
void RefineableYoungLaplaceProblem<ELEMENT>::create_contact_angle_elements(
 const unsigned &b)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));

   // What is the index of the face of the bulk element at the boundary
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding contact angle element
   YoungLaplaceContactAngleElement<ELEMENT>* contact_angle_element_pt = new 
   YoungLaplaceContactAngleElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the contact angle element to the contact angle mesh
   Contact_angle_mesh_pt->add_element_pt(contact_angle_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_contact_angle_elements


//============start_of_delete_contact_angle_elements=====================
/// Delete YoungLaplace contact angle elements
//=======================================================================
template<class ELEMENT>
void RefineableYoungLaplaceProblem<ELEMENT>::delete_contact_angle_elements()
{

 // How many contact angle elements are there?
 unsigned n_element = Contact_angle_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Contact_angle_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 Contact_angle_mesh_pt->flush_element_and_node_storage();


} // end of delete_contact_angle_elements



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void RefineableYoungLaplaceProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                              ofstream& trace_file)
{ 

 // Output kappa vs height
 //-----------------------
 trace_file << -1.0*Kappa_pt->value(0) << " ";
 trace_file << Control_node_pt->value(0) ;
 trace_file << endl;
  
 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output full solution 
 //---------------------
 ofstream some_file;
 char filename[100];
 //YoungLaplaceEquations::Output_meniscus_and_spines=false;
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output contact angle 
 //---------------------

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
 
 // How many contact angle elements are there?
 unsigned n_element = Contact_angle_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   
   tangent_file << "ZONE" << std::endl;
   normal_file << "ZONE" << std::endl;
   contact_angle_file << "ZONE" << std::endl;
   
   // Upcast from GeneralisedElement to YoungLaplace contact angle element
   YoungLaplaceContactAngleElement<ELEMENT>* el_pt = 
    dynamic_cast<YoungLaplaceContactAngleElement<ELEMENT>*>(
     Contact_angle_mesh_pt->element_pt(e));
   
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
   
   
  } // end of loop over both boundaries
 
 tangent_file.close();
 normal_file.close();
 contact_angle_file.close();
 

cout << "\n********************************************" << endl <<  endl;

} // end of doc



//===============start_of_main============================================
/// Drive code
//========================================================================
int main()
{

 // Create label for output
 DocInfo doc_info;

 // Trace file
 ofstream trace_file;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Open a trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);

 // Tecplot header for trace file: kappa and height value
 trace_file << "VARIABLES=\"<GREEK>k</GREEK>\",\"h\"" << std::endl;
 trace_file << "ZONE" << std::endl;


 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // RefineableQYoungLaplaceElement family. 
 RefineableYoungLaplaceProblem<RefineableQYoungLaplaceElement<3> > problem;

 // Perform one uniform refinement
 problem.refine_uniformly();

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
   GlobalParameters::Controlled_height+=increment;

   // Solve the problem  
   unsigned max_adapt=1;
   problem.newton_solve(max_adapt);
 
   //Output the solution
   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }
 
 // Close output file
 trace_file.close();
 
} //end of main



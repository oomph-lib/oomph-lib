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
//Driver function for a simple test ring problem with 
//non-axis-aligned symmetry boundary conditions

//Generic oomph-lib sources
#include "generic.h"

// Beam sources
#include "beam.h"

// The mesh
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;

using namespace oomph;


/// ////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////



//=======start_of_namespace=========================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Nondim thickness 
 double H=0.05;
 
 /// Prescribed position (only used for displacement control)
 double Xprescr = 1.0;

 /// Perturbation pressure
 double Pcos=0.0;

 /// Pointer to pressure load
 Data* Pext_data_pt;

 /// Buckling wavenumber
 unsigned Nbuckl=3;

 /// Load function: Constant external pressure with cos variation to
 /// induce buckling in desired mode
 void press_load(const Vector<double>& xi,
                 const Vector<double> &x,
                 const Vector<double>& N,
                 Vector<double>& load)
 {
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = (Pext_data_pt->value(0)+Pcos*cos(double(Nbuckl)*xi[0]))*N[i];
   }
 }
 
 /// Return a reference to the external pressure 
 /// load on the elastic ring.
 /// A reference is obtained by de-referencing the pointer to the
 /// data value that contains the external load
 double &external_pressure() 
  {return *Pext_data_pt->value_pt(0);}


} // end of namespace



/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//========start_of_problem_class========================================
/// Ring problem
//======================================================================
template<class ELEMENT>
class ElasticRingProblem : public Problem
{

public:

 /// Constructor: Pass umber of elements
 ElasticRingProblem(const unsigned &n_element);

 /// Access function for the specific mesh
 OneDLagrangianMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<OneDLagrangianMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Update function is empty 
 void actions_after_newton_solve() {}

 /// Update function is empty 
 void actions_before_newton_solve() {}

 /// Doc solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);

 /// Perform the parameter study
 void parameter_study(DocInfo& doc_info);

private:

 /// Pointer to geometric object that represents the undeformed shape
 GeomObject* Undef_geom_pt;

 /// Number of elements in the beam mesh
 unsigned Nbeam_element;

}; // end of problem class







//======start_of_constructor============================================
/// Constructor for elastic ring problem
//======================================================================
template<class ELEMENT>
ElasticRingProblem<ELEMENT>::ElasticRingProblem(const unsigned& n_element)
{

 // Undeformed beam is an elliptical ring 
 Undef_geom_pt=new Ellipse(1.0,1.0); 

 // Angle
 double phi=MathematicalConstants::Pi/
  double(Global_Physical_Variables::Nbuckl);

 // Length of the doamin (in terms of the Lagrangian coordinates)
 double length=phi;

 // Actual number of elements
 Nbeam_element=n_element;

 //Now create the (Lagrangian!) mesh
 Problem::mesh_pt() = 
  new OneDLagrangianMesh<ELEMENT>(Nbeam_element,length,Undef_geom_pt); 

 // "Bulk" beam element that the boundary condition element is attached to
 int face_index=1;
 ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
  mesh_pt()->element_pt(n_element-1));

 // Create boundary condition element
 ClampedSlidingHermiteBeamBoundaryConditionElement*
  bc_el_pt = new
  ClampedSlidingHermiteBeamBoundaryConditionElement(
   bulk_el_pt,face_index);
 
 // Set vectors that define the symmetry line
 Vector<double> vector_to_symmetry_line(2);
 vector_to_symmetry_line[0]=0.0;
 vector_to_symmetry_line[1]=0.0;
 
 Vector<double> normal_to_symmetry_line(2);
 normal_to_symmetry_line[0]=-sin(phi);
 normal_to_symmetry_line[1]= cos(phi);
 
 bc_el_pt->set_symmetry_line(vector_to_symmetry_line,
                             normal_to_symmetry_line);
 
 // Add boundary condition element to mesh
 mesh_pt()->add_element_pt(bc_el_pt);
 

 // Boundary condition: Only at bottom
 //-----------------------------------

 // Bottom: 
 unsigned ibound=0;

 // No vertical displacement
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(1); 

 // Infinite slope: Pin type 1 (slope) dof for displacement direction 0 
 mesh_pt()->boundary_node_pt(ibound,0)->pin_position(1,0);


 // Displacement control
 //---------------------

 // Choose element in which displacement control is applied: the first one
 SolidFiniteElement* controlled_element_pt=
  dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));
 
 // Fix the displacement in the horizontal direction...
 unsigned controlled_direction=0;
 
 //... at left end of the control element
 Vector<double> s_displ_control(1);
 s_displ_control[0]=0.0;
 
 // Pointer to displacement control element
 DisplacementControlElement* displ_control_el_pt;
 
 
 // Build displacement control element
 displ_control_el_pt=
  new DisplacementControlElement(controlled_element_pt,
                                 s_displ_control,
                                 controlled_direction,
                                 &Global_Physical_Variables::Xprescr);
 
 // The constructor of the  DisplacementControlElement has created
 // a new Data object whose one-and-only value contains the
 // adjustable load: Use this Data object in the load function:
 Global_Physical_Variables::Pext_data_pt=displ_control_el_pt->
  displacement_control_load_pt();


 // Add the displacement-control element to the mesh
 mesh_pt()->add_element_pt(displ_control_el_pt); 

 //Loop over the elements to set physical parameters etc.
 for(unsigned i=0;i<Nbeam_element;i++)
  {
   //Cast to proper element type
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   // Set wall thickness
   elem_pt->h_pt() = &Global_Physical_Variables::H;
   
   // Function that specifies load Vector
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::press_load;

   //Assign the undeformed beam shape
   elem_pt->undeformed_beam_pt() = Undef_geom_pt;

   // Displacement control: The load on *all* elements
   // is affected by an unknown -- the external pressure, stored
   // as the one-and-only value in a Data object: Add it to the
   // elements' external Data.
   elem_pt->add_external_data(Global_Physical_Variables::Pext_data_pt);
  } 
 
 
 
 // Do equation numbering
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor


//=======start_of_doc=====================================================
/// Document solution
//========================================================================
template<class ELEMENT>
void ElasticRingProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                               ofstream& trace_file)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts=5;
 
 // Output solution 
 sprintf(filename,"%s/ring%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for(unsigned i=0;i<Nbeam_element;i++)
  {
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   el_pt->output(some_file,npts);
  }
 some_file.close();

 // Local coordinates of plot points at left and right end of domain
 Vector<double> s_left(1);
 s_left[0]=-1.0;
 Vector<double> s_right(1);
 s_right[0]=1.0;

 // Write trace file: Pressure, two radii
 trace_file 
  << Global_Physical_Variables::Pext_data_pt->value(0)/
  (pow(Global_Physical_Variables::H,3)/12.0)  
  << " " 
  << sqrt(
   pow(dynamic_cast<ELEMENT*>(mesh_pt()->
                              element_pt(0))->interpolated_x(s_left,0),2)+
   pow(dynamic_cast<ELEMENT*>(mesh_pt()->
                              element_pt(0))->interpolated_x(s_left,1),2))
  << " " 
  << sqrt(
   pow(dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(Nbeam_element-1))->
       interpolated_x(s_right,0),2)+
   pow(dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(Nbeam_element-1))->
       interpolated_x(s_right,1),2))
  << std::endl;

  
} // end of doc



//========start_of_run=====================================================
/// Solver loop to perform parameter study
//=========================================================================
template<class ELEMENT>
void ElasticRingProblem<ELEMENT>::parameter_study(DocInfo& doc_info)
{

 //Open a trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 ofstream trace_file(filename);
 trace_file << "VARIABLES=\"p_e_x_t\",\"R_1\",\"R_2\"" << std::endl;
 trace_file << "ZONE" << std::endl;

 //Output initial data
 doc_solution(doc_info,trace_file);

 // Number of steps
 unsigned nstep= 11; 
 
 // Increment and initial value (which will immediately be decremented)
 double displ_increment=1.0/double(nstep-1);
 Global_Physical_Variables::Xprescr=1.0+displ_increment;


 // Downward loop over parameter incrementation with pcos>0 
 //--------------------------------------------------------
 
 /// Perturbation pressure
 Global_Physical_Variables::Pcos=1.0e-4; 


 // Downward loop over parameter incrementation
 //---------------------------------------------
 for(unsigned i=1;i<=nstep;i++)
  {
   // Increment control displacement
   Global_Physical_Variables::Xprescr-=displ_increment;

   // Solve
   newton_solve();
   
   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info,trace_file);
   
  } // end of downward loop
 
 // Reset perturbation pressure
 //----------------------------
 Global_Physical_Variables::Pcos=0.0;
 
 // Set initial values for parameters that are to be incremented
 Global_Physical_Variables::Xprescr-=displ_increment;

 // Start new zone for tecplot
 trace_file << "ZONE" << std::endl;

 // Upward loop over parameter incrementation
 //------------------------------------------
 for(unsigned i=nstep;i<2*nstep;i++)
  {
   // Increment control displacement
   Global_Physical_Variables::Xprescr+=displ_increment;

   // Solve
   newton_solve();
   
   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info,trace_file);
   
  } 
 
} // end of run



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=======start_of_main=================================================
/// Driver for ring test problem 
//=====================================================================
int main()
{
 // Number of elements
 unsigned n_element = 13;
 
 // Label for output
 DocInfo doc_info;
 
 //Set up the problem
 ElasticRingProblem<HermiteBeamElement> problem(n_element);
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 // Do static run
 problem.parameter_study(doc_info);


} // end of main







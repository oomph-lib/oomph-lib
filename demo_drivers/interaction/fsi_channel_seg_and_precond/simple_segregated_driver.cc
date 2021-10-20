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
// Generic oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "beam.h"
#include "multi_physics.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/collapsible_channel_mesh.h"


using namespace std;

using namespace oomph;


// Include the general-purpose fsi collapsible channel problem
#include "fsi_chan_problem.h"

//====Namespace_for_flags================================
/// Extend namespace for control flags
//======================================================
namespace Flags
{

 /// Use Newton solver (0) or segregated solver (1)?
 unsigned Use_segregated_solver=1;

 /// Use pointwise Aitken extrapolation (1) or not (0)
 unsigned Use_pointwise_aitken=0;

 /// Under-relaxation parameter (1.0: no under-relaxation; 0.0: freeze) 
 double Omega_under_relax=1.0;

 /// Use Irons and Tuck extrapolation (1) or not (0)
 unsigned Use_irons_and_tuck_extrapolation=0;

 /// Convergence criterion: 0: global resmax; 1: abs. change; 2: rel. change
 unsigned Convergence_criterion=0;

 /// Convergence tolerance
 double Convergence_tolerance=1.0e-8;

}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//====start_of_problem_class==========================================
/// Problem class -- add segregated solver capability to an existing
/// problem.
//====================================================================
template <class ELEMENT>
class SegregatedFSICollapsibleChannelProblem : 
 public virtual FSICollapsibleChannelProblem<ELEMENT>,
 public virtual SegregatableFSIProblem
{

public :
 
 ///  Constructor: The arguments are the same as the original
 /// (non-segregated) problem, namely, numbers of elements and lengths
 /// of different sections of the domain.
 SegregatedFSICollapsibleChannelProblem(const unsigned& nup, 
                                        const unsigned& ncollapsible,
                                        const unsigned& ndown,
                                        const unsigned& ny,
                                        const double& lup,
                                        const double& lcollapsible, 
                                        const double& ldown,
                                        const double& ly,
                                        const bool& displ_control,
                                        const bool& steady_flag);
 
 /// Empty Destructor 
 ~SegregatedFSICollapsibleChannelProblem(){}


 ///  Identify the fluid and solid Data and meshes that
 /// contain only elements involved in the respective sub-problems. 
 /// This is a specific implementation of a pure virtual function in the 
 /// SegregatableFSIProblem base class.
 void identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                                    Vector<Data*>& solid_data_pt,
                                    Mesh*& fluid_mesh_pt,
                                    Mesh*& solid_mesh_pt);

 // start_of_convergence_checks

 ///  Update nodal positions in the fluid mesh in
 /// response to changes in the wall displacement field after every
 /// Newton step in a monolithic or segregated solid solve. Note
 /// the use of the (protected) flag Solve_type, which can take the
 /// values Full_solve, Fluid_solve or Solid_solve. This flag is used
 /// to allow specification of different actions depending on the
 /// precise solve taking place.
 void actions_before_newton_convergence_check()
  {
   //For a "true" segregated solver, we would not do this in fluid or solid
   //solves, but adding the bulk node update to the solid solve phase aids
   //convergence and makes it possible for larger values of Q. Of course,
   //there is a small cost associated with doing this.
   if(Solve_type!=Fluid_solve) {this->Bulk_mesh_pt->node_update();}
  }


 /// Update nodal positions in the fluid mesh
 /// in response to any changes in the wall displacement field after every 
 /// segregated solve. This is not strictly necessary because we
 /// do the solid solve last, which performs its own node update before the 
 /// convergence check of the sub problem. It remains here because if we
 /// were solving in a completely segregated fashion a node update would be 
 /// required for the fluid mesh in the final converged solution to be
 /// consistent with the solid positions.
 void actions_before_segregated_convergence_check()
  { 
   this->Bulk_mesh_pt->node_update();
  } 

 // end_of_convergence_checks

 /// Document the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Perform a steady run
 void steady_run();

};


//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
SegregatedFSICollapsibleChannelProblem< ELEMENT>::
SegregatedFSICollapsibleChannelProblem(const unsigned& nup, 
                                       const unsigned& ncollapsible,
                                       const unsigned& ndown,
                                       const unsigned& ny,
                                       const double& lup,
                                       const double& lcollapsible, 
                                       const double& ldown,
                                       const double& ly,
                                       const bool& displ_control,
                                       const bool& steady_flag) :
 FSICollapsibleChannelProblem<ELEMENT>(nup, 
                                       ncollapsible,
                                       ndown,
                                       ny,
                                       lup,
                                       lcollapsible, 
                                       ldown,
                                       ly,
                                       displ_control,
                                       steady_flag) 
{
 // Choose convergence criterion based on Flag::Convergence criterion
 // with tolerance given by Flag::Convergence_tolerance
 if (Flags::Convergence_criterion==0)
  {
   assess_convergence_based_on_max_global_residual(
    Flags::Convergence_tolerance);
  }
 else if (Flags::Convergence_criterion==1)
  {
   assess_convergence_based_on_absolute_solid_change(
    Flags::Convergence_tolerance);
  }
 else if (Flags::Convergence_criterion==2)
  {
   assess_convergence_based_on_relative_solid_change(
    Flags::Convergence_tolerance);
  }
 
 //Select a convergence-acceleration technique based on control flags

 // Pointwise Aitken extrapolation
 if(Flags::Use_pointwise_aitken)
  {
   this->enable_pointwise_aitken();
  }
 else
  {
   this->disable_pointwise_aitken();
  }

 // Under-relaxation
 this->enable_under_relaxation(Flags::Omega_under_relax);

 // Irons and Tuck's extrapolation
 if(Flags::Use_irons_and_tuck_extrapolation)
  {
   this->enable_irons_and_tuck_extrapolation();
  }
 else
  {
   this->disable_irons_and_tuck_extrapolation();
  }

} //end_of_constructor



//=====start_of_identify_fluid_and_solid======================================
/// Identify the fluid and solid Data and the meshes that
/// contain only elements that are involved in the respective sub-problems. 
/// This implements a pure virtual function in the 
/// SegregatableFSIProblem base class.
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::
identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                              Vector<Data*>& solid_data_pt,
                              Mesh*& fluid_mesh_pt,
                              Mesh*& solid_mesh_pt)
{

 //FLUID DATA: 
 //All fluid elements are stored in the Mesh addressed by bulk_mesh_pt() 

 //Reset the storage
 fluid_data_pt.clear();

 //Find number of fluid elements
 unsigned n_fluid_elem=this->bulk_mesh_pt()->nelement();
 //Loop over fluid elements and add internal data to fluid_data_ptt
 for(unsigned e=0;e<n_fluid_elem;e++)
  {
   GeneralisedElement* el_pt=this->bulk_mesh_pt()->element_pt(e);
   unsigned n_internal=el_pt->ninternal_data();
   for(unsigned i=0;i<n_internal;i++)
    {
     fluid_data_pt.push_back(el_pt->internal_data_pt(i));
    }
  }
 
 //Find number of nodes in fluid mesh
 unsigned n_fluid_node=this->bulk_mesh_pt()->nnode();
 //Loop over nodes and add the nodal data to fluid_data_pt
 for (unsigned n=0;n<n_fluid_node;n++)
  {
   fluid_data_pt.push_back(this->bulk_mesh_pt()->node_pt(n));
  }
  
 // The bulk_mesh_pt() is a mesh that contains only fluid elements
 fluid_mesh_pt = this->bulk_mesh_pt(); 
 

 //SOLID DATA
 //All solid elements are stored in the Mesh addressed by wall_mesh_pt()

 //Reset the storage
 solid_data_pt.clear();

 //Find number of nodes in the solid mesh
 unsigned n_solid_node=this->wall_mesh_pt()->nnode();
 //Loop over nodes and add nodal position data to solid_data_pt
 for(unsigned n=0;n<n_solid_node;n++)
  {
   solid_data_pt.push_back(
    this->wall_mesh_pt()->node_pt(n)->variable_position_pt());
  }
   
 //If we are using displacement control then the displacement control element
 //and external pressure degree of freedom should be treated as part
 //of the solid problem

 //We will assemble a single solid mesh from a vector of pointers to meshes
 Vector<Mesh*> s_mesh_pt(1);
 //The wall_mesh_pt() contains all solid elements and is the first
 //entry in our vector
 s_mesh_pt[0]=this->wall_mesh_pt();
  
 //If we are using displacement control
 if (this->Displ_control)
  {
   //Add the external pressure data to solid_data_pt
   solid_data_pt.push_back(Global_Physical_Variables::P_ext_data_pt);
   //Add a pointer to a Mesh containing the displacement control element
   //to the vector of pointers to meshes
   s_mesh_pt.push_back(this->Displ_control_mesh_pt);
  } 

 // Build "combined" mesh from our vector of solid meshes
 solid_mesh_pt = new Mesh(s_mesh_pt);

} //end_of_identify_fluid_and_solid


//====start_of_doc_solution============================================
/// Document the solution
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>:: doc_solution(
 DocInfo& doc_info)
{ 
 //Output stream filenames
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5;

 // Output fluid solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 this->bulk_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Document the wall shape
 sprintf(filename,"%s/beam%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 this->wall_mesh_pt()->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution



//====steady_run==============================================================
/// Perform a steady run in which the external pressure 
/// (or presribed displacement) is varied causing the channel to collapse.
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::steady_run()
{ 

 // Set initial value for external pressure (on the wall stiffness scale). 
 // This can be overwritten in set_initial_condition.
 Global_Physical_Variables::P_ext_data_pt->
  set_value(0,Global_Physical_Variables::Pmin);
 
 // Apply initial condition
 set_initial_condition();
 
 //Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Output the initial solution
 doc_solution(doc_info);
 
 // Increment step number
 doc_info.number()++;
  
 // Increment for external pressure (on the wall stiffness scale)
 double delta_p=(Global_Physical_Variables::Pmax-
                 Global_Physical_Variables::Pmin)/double(Flags::Nsteps-1);
 
 // Initial and final values for control position
 Global_Physical_Variables::Yprescr=1.0;
 
 // Final value of prescribed displacement 
 double y_min=0.65;
 //Change in y-coordinate to attain the minimum position in a given
 //number of steps
 double delta_y=(y_min-Global_Physical_Variables::Yprescr)/
  double(Flags::Nsteps-1);

 // Parameter study (loop over the number of steps)
 for (unsigned istep=0;istep<Flags::Nsteps;istep++)
  {
   // Setup segregated solver 
   //(Default behaviour will identify the fluid and solid dofs and
   // allocate memory, etc every time. This is a bit inefficient in 
   // this case, but it is safe and will always work)
   setup_segregated_solver();

   // SEGREGATED SOLVER
   if(Flags::Use_segregated_solver)
    {
     //Set the maximum number of Picard steps
     Max_picard =50;
     
     // Solve ignoring return type (convergence data)
     (void)steady_segregated_solve();
    }
   // NEWTON SOLVER
   else
    {
     //Explit call to the steady Newton solve.
     steady_newton_solve();
    }
   
   // Output the solution
   doc_solution(doc_info);
   
   //Increase the Step number
   doc_info.number()++;
   
   // Adjust control parameters
   //If displacment control increment position
   if (this->Displ_control)
    {
     Global_Physical_Variables::Yprescr+=delta_y;
    }
   //Otherwise increment external pressure
   else
    {
     double old_p=Global_Physical_Variables::P_ext_data_pt->value(0);
     Global_Physical_Variables::P_ext_data_pt->set_value(0,old_p+delta_p);
    }

  } // End of parameter study

}

//============start_of_main====================================================
/// Driver code for a segregated collapsible channel problem with FSI.
//=============================================================================
int main()
{
 // Number of elements in the domain
 unsigned nup=4*Flags::Resolution_factor;
 unsigned ncollapsible=20*Flags::Resolution_factor;
 unsigned ndown=40*Flags::Resolution_factor;
 unsigned ny=4*Flags::Resolution_factor;
  
 
 // Geometry of the domain
 double lup=1.0;
 double lcollapsible=5.0;
 double ldown=10.0;
 double ly=1.0;
 
 // Steady run by default
 bool steady_flag=true;
 // with displacement control
 bool displ_control=true;

 // Build the problem with QTaylorHoodElements
 SegregatedFSICollapsibleChannelProblem
  <AlgebraicElement<QTaylorHoodElement<2> > > 
  problem(nup, ncollapsible, ndown, ny, 
          lup, lcollapsible, ldown, ly, displ_control,
          steady_flag);
 
 //Perform a steady run
 problem.steady_run();
 
}//end of main

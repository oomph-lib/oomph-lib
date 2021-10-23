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
// Include the general-purpose fsi driven cavity problem
#include "fsi_driven_cavity_problem.h"


#include "multi_physics.h"



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//====start_of_problem_class==========================================
/// Problem class -- add segregated solver capability to existing
/// problem.
//====================================================================
template <class ELEMENT>
class SegregatedFSIDrivenCavityProblem : 
 public virtual FSIDrivenCavityProblem<ELEMENT>,
 public virtual SegregatableFSIProblem
{

public :



 /// Constructor: The arguments are the number of elements,
 /// the lengths of the domain, the fractional height of the gap
 /// next to the moving lid and the period of the lid's oscillation 
 SegregatedFSIDrivenCavityProblem(const unsigned& nx, 
                                  const unsigned& ny,
                                  const double& lx,
                                  const double& ly,
                                  const double& gap_fraction,
                                  const double& period);
 
 
 /// Destructor (empty)
 ~SegregatedFSIDrivenCavityProblem(){}

 /// Identify the fluid and solid Data and the meshes that
 /// contain only elements that are involved in the respective sub-problems. 
 /// This implements a pure pure virtual function in the 
 /// SegregatableFSIProblem base class.
 void identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                                    Vector<Data*>& solid_data_pt,
                                    Mesh*& fluid_mesh_pt,
                                    Mesh*& solid_mesh_pt);
 
 /// Overload empty virtual function that is called before header for
 /// convergence history is written. This can be overloaded to insert
 /// zone information that can be used to identify 
 /// the problem parameters for this solve.
 void write_zone_info_for_convergence_history(ofstream& convergence_file);


 /// empty final overload
 void actions_after_newton_solve()  {}


 /// Initialise timer and reset counter for Newton iterations
 /// if monolithic solver is used.
 void actions_before_newton_solve()
  {
//    if (!Flags::Use_segregated_solver)
//     {
//      // Initialise counter for Newton iteration
//      Newton_iter=0;
     
//      // Initialise timer that allows doc of iteration/cpu time
//      reset_timer();
//     }
  }

 /// Overload actions before Newton step: Update nodal 
 /// positions in the fluid mesh in response to any changes in 
 /// the wall displacement field if segregated solver is used.
 void actions_before_newton_step()
  {
   this->Bulk_mesh_pt->node_update();
  }

 /// Overload actions after Newton step: Update nodal 
 /// positions in the fluid mesh
 /// in response to any changes in the wall displacement field. If 
 /// monolithic Newton solver is used, doc progress of Newton iteration, 
 /// using the same output as during Picard iteration.
 void actions_before_newton_convergence_check();
 
private:

 /// Output stream to document the convergence history
 ofstream Convergence_file;

};








//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
SegregatedFSIDrivenCavityProblem< ELEMENT>::
 SegregatedFSIDrivenCavityProblem(const unsigned& nx, 
                                  const unsigned& ny,
                                  const double& lx,
                                  const double& ly,
                                  const double& gap_fraction,
                                  const double& period) :
 FSIDrivenCavityProblem<ELEMENT>(nx, ny, lx, ly, gap_fraction, period)
{
 Max_picard = 50;
//  // Choose convergence criterion
//  if (Flags::Convergence_criterion==0)
//   {
//    assess_convergence_based_on_max_global_residual();
//   }
//  else if (Flags::Convergence_criterion==1)
//   {
//    assess_convergence_based_on_absolute_solid_change();
//   }
//  else if (Flags::Convergence_criterion==2)
//   {
//    assess_convergence_based_on_relative_solid_change();
//   }
 
 
 // Convergence criterion
 assess_convergence_based_on_absolute_solid_change(1.0e-7);


 // Pointwise Aitken extrapolation?
 this->disable_pointwise_aitken(); // Flags::Use_pointwise_aitken;

 // Use under-relaxation?
 //enable_under_relaxation(Flags::Omega_under_relax);
 this->enable_under_relaxation(1.0e-2);

 // Use Irons and Tuck's extrapolation
 //this->enable_irons_and_tuck_extrapolation();

 // Doc max. global residual during Picard iteration
 Doc_max_global_residual=true;
 
 // Number of wall control points
 unsigned n_control=10;
 Vector<std::pair<GeomObject* ,Vector<double> > >
  control_point_pairs(n_control);

 // Get wall control points
 for (unsigned i=0;i<n_control;i++)
  {
   // Get pointer to/local coordinate in wall element that contains 
   // control node 
   Vector<double> zeta_ctrl(1);
   zeta_ctrl[0]=double(i+1)/double(n_control+1);
   control_point_pairs[i].second.resize(1);
   this->Wall_geom_object_pt->locate_zeta(zeta_ctrl,control_point_pairs[i].first,
                                    control_point_pairs[i].second);
 
   Vector<double> r_ctrl(2);
   control_point_pairs[i].first->position(control_point_pairs[i].second,
                                          r_ctrl);
  }

 // Identify control points for doc of picard convergence
 //identify_solid_control_points(control_point_pairs);
 
}


//============================================================================
/// Actions after Newton step: Update nodal positions in the fluid mesh
/// in response to any changes in the wall displacement field. If 
/// monolithic Newton solver is used, doc progress of Newton iteration, 
/// using the same output as during Picard iteration.
//============================================================================
template <class ELEMENT>
void SegregatedFSIDrivenCavityProblem<ELEMENT>::
actions_before_newton_convergence_check()
{ 

 this->Bulk_mesh_pt->node_update();

// // if (!Flags::Use_segregated_solver)
//   {
//    // Halt timer
//    halt_timer();
    
//    double rms_change;
//    double max_change;
//    double rms_norm;
//    double max_res=0.0;
//    get_solid_change(rms_change,max_change,rms_norm);

//    Vector<double> residual(this->ndof());
//    this->get_residuals(residual);
//    max_res=std::abs(*std::max_element(residual.begin(),
//                                       residual.end(),          
//                                       AbsCmp<double>()));    
   
//    std::cout << "==================================================\n";
//    std::cout <<   "Iteration             : " 
//              << Newton_iter << std::endl;
//    std::cout <<   "RMS  change           :       "
//              << rms_change << std::endl;
//    std::cout <<   "Max. change           :       "
//              << max_change << std::endl;
//    std::cout <<   "RMS norm              :       "
//              << rms_norm   << std::endl;
//    std::cout << "Max. global residual  :       "
//              << max_res   << std::endl;
//    std::cout << "==================================================\n\n";

//    bool get_max_global_residual=true;
//    write_convergence_history(Newton_iter,
//                              rms_change,
//                              max_change,
//                              rms_norm,
//                              max_res,
//                              Convergence_file, 
//                              get_max_global_residual);
  
//    // Store the current values of the solid dofs as reference values
//    store_solid_dofs();
   
//    // Increment counter
//    Newton_iter++;

//    // Restart timer
//    restart_timer();
//   }
}


//============================================================================
/// Overload empty virtual function that is called before header for
/// convergence history is written. Overloaded to insert
/// zone information that can be used to identify 
/// the problem parameters for this solve.
//============================================================================
template<class ELEMENT>
void SegregatedFSIDrivenCavityProblem<ELEMENT>::
write_zone_info_for_convergence_history(ofstream& convergence_file)
{

 convergence_file << "ZONE T=\"";
 convergence_file << "Re=" << Global_Physical_Variables::Re << ", ";
 convergence_file << "Q=" << Global_Physical_Variables::Q << " \"" << std::endl;
// convergence_file << "resolution factor: " << Flags::Resolution_factor << " ";
//  if (Flags::Use_segregated_solver)
//   {
//    convergence_file << "Picard";
//   }
//  else
//   {
//    convergence_file << "Newton";
//   }
//  convergence_file << "\""<< std::endl;

}


//=====start_of_identify_fluid_and_solid======================================
/// Identify the fluid and solid Data and the meshes that
/// contain only elements that are involved in the respective sub-problems. 
/// This implements a pure pure virtual function in the 
/// SegregatableFSIProblem base class.
//============================================================================
template <class ELEMENT>
void SegregatedFSIDrivenCavityProblem<ELEMENT>::
identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                              Vector<Data*>& solid_data_pt,
                              Mesh*& fluid_mesh_pt,
                              Mesh*& solid_mesh_pt)
{

 // Fluid
 //------

 // Reset
 fluid_data_pt.clear();
 
 // Get internal data of fluid elements in bulk
 unsigned n_fluid_elem=this->bulk_mesh_pt()->nelement();
 for (unsigned e=0;e<n_fluid_elem;e++)
  {
   GeneralisedElement* el_pt=this->bulk_mesh_pt()->element_pt(e);
   unsigned n_internal=el_pt->ninternal_data();
   for (unsigned j=0;j<n_internal;j++)
    {
     fluid_data_pt.push_back(el_pt->internal_data_pt(j));
    }
  }
 
 // Loop over nodes in fluid bulk mesh
 unsigned n_fluid_node=this->bulk_mesh_pt()->nnode();
 for (unsigned j=0;j<n_fluid_node;j++)
  {
   fluid_data_pt.push_back(this->bulk_mesh_pt()->node_pt(j));
  }
  
 // Here's the mesh that contains only fluid elements:
 fluid_mesh_pt = this->bulk_mesh_pt(); 
  

 // Solid
 //------

 // Reset
 solid_data_pt.clear();

 // Loop over nodes in solid mesh 
 unsigned n_solid_node=this->wall_mesh_pt()->nnode();
 for (unsigned j=0;j<n_solid_node;j++)
  {
   solid_data_pt.push_back(
    this->wall_mesh_pt()->node_pt(j)->variable_position_pt());
  }
   

 // Build  the mesh that contains only solid elements:
 Vector<Mesh*> s_mesh_pt(1);
 s_mesh_pt[0]=this->wall_mesh_pt();
  
 // Build "combined" mesh
 solid_mesh_pt = new Mesh(s_mesh_pt);

}









///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////





//============start_of_main====================================================
/// Driver code for a driven cavity problem with FSI.
/// Presence of command line arguments indicates validation run with 
/// coarse resolution and small number of steps.
//=============================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Number of elements in the domain (W. Wall standard resolution)
 unsigned nx=16; //32;
 unsigned ny=16; //32;

 // Length of the domain
 double lx=1.0;
 double ly=1.0;

 // Fractional height of the gap next to the moving lid
 double gap_fraction=1.0/16.0;
   
 // Period of lid oscillation
 double period=5.0;

 // Build the problem
 SegregatedFSIDrivenCavityProblem<AlgebraicElement<QTaylorHoodElement<2> > > 
  problem(nx, ny, lx, ly, gap_fraction, period);

 // Timestep. 
 double dt=period/40.0; 
 
 // Maximum time for simulation
 double t_max=50.0; 

 // Initialise timestep 
 problem.initialise_dt(dt);

 // Apply initial condition
 problem.set_initial_condition();

 //Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 
 // Output the initial solution
 problem.doc_solution(doc_info, trace_file);

 // Increment step number
 doc_info.number()++;

 // Find number of timesteps (reduced for validation)
 unsigned nstep = unsigned(t_max/dt);
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }
 
 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   
   // Object that stores the convergence data for Picard iteration
   PicardConvergenceData conv_data;
   
   // Number iterations taken
   //unsigned niter=0; 
   
   // Setup segregated solver in either case
   problem.setup_segregated_solver();
   
   // Solve the problem
   //------------------

   conv_data=problem.unsteady_segregated_solve(dt);

   //niter=conv_data.niter();
   if (conv_data.cpu_total()!=0.0)
    {
     std::cout 
      << std::endl << std::endl << std::endl 
      << "PERCENTAGE OF CPU TIME SPENT IN COMPUTATION OF GLOBAL RESIDUAL="
      << double(conv_data.cpu_for_global_residual()/
                conv_data.cpu_total())*100.0
      << std::endl << std::endl << std::endl <<std::endl;
     
     std::cout 
      <<"PERCENTAGE OF CPU TIME SPENT IN COMPUTATION OF NON-ESSENTIAL BITS="
      << double((conv_data.cpu_total()-conv_data.essential_cpu_total())/
                conv_data.cpu_total())*100.0 << " " 
      << conv_data.cpu_total() << " " 
      << conv_data.essential_cpu_total() << " "
      << std::endl << std::endl << std::endl <<std::endl;
    }
  
   

   // Solve the problem
   // problem.unsteady_newton_solve(dt);
   
   // Outpt the solution
   problem.doc_solution(doc_info, trace_file);
   
   // Step number
   doc_info.number()++;
  }


 // Close trace file.
 trace_file.close();



}//end of main


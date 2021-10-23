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


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//====Namespace_for_flags================================
/// Extend namespace for flags
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
 double Convergence_tolerance=1.0e-9;

 /// Timestep
 double Dt=0.1;

 /// Doc flags for extended namespace
 void segregated_doc_flags()
 {
  std::cout << "\nFlags:\n"
            <<   "======\n"; 

  std::cout << "-- Resolution factor: " << Resolution_factor << std::endl;


  if (Use_segregated_solver)
   {
    std::cout << "-- Using segregated solver " << std::endl;
   }
  else
   {
    std::cout << "-- Using monolithic Newton solver " << std::endl;
   }


  if (Use_pointwise_aitken)
   {
    std::cout << "-- Using pointwise Aitken extrapolation " << std::endl;
   }
  else
   {
    std::cout << "-- Not using pointwise Aitken extrapolation " << std::endl;
   }



  if (Use_irons_and_tuck_extrapolation)
   {
    std::cout << "-- Using Irons and Tuck's extrapolation " << std::endl;
   }
  else
   {
    std::cout << "-- Not using Irons and Tuck's extrapolation " << std::endl;
   }

  
  if (Convergence_criterion==0)
   {
    std::cout << "-- Convergence based on max. global residual " << std::endl;
   }
  else if (Convergence_criterion==1)
   {
    std::cout << "-- Convergence based on max. abs. change " << std::endl;
   }
  else if (Convergence_criterion==2)
   {
    std::cout << "-- Convergence based on max. rel. change " << std::endl;
   }
  else
   {
    std::cout << "-- Wrong convergence criterion [0,1,2]: " 
              << Convergence_criterion << std::endl;
    exit(1);
   }


  std::cout << "-- Convergence tolerance: " 
            << Convergence_tolerance << std::endl;
   

  std::cout << "-- Under-relaxation parameter: " 
            << Omega_under_relax << std::endl;
   
  if (Steady_flag)
   { 
    std::cout << "-- Steady run " << std::endl;
    if (Use_displ_control)
     {
      std::cout << "-- Using displacement control " << std::endl;
     }
    else
     {
      std::cout << "-- Not using displacement control " << std::endl;
     }
   }
  else
   {
    std::cout << "-- Unsteady run " << std::endl;
    if (Use_displ_control)
     {
      std::cout << "-- Not using displacement control (command line flag\n"
                << "   overwritten because it's an unsteady run!) " 
                << std::endl;
     }
   }


  std::cout << "-- Reynolds number: " 
            << Global_Physical_Variables::Re  << std::endl;


  std::cout << "-- FSI parameter Q: " 
            << Global_Physical_Variables::Q  << std::endl;


  std::cout << "-- Timestep dt    : "  << Flags::Dt << std::endl;


  if (Restart_file_name!="")
   {
    std::cout << "-- Performing restart from: " << Restart_file_name 
              << std::endl;
    std::cout << "-- Jump in pressure: " << Global_Physical_Variables::P_step
              << std::endl;

   }
  else
   {
    std::cout << "-- No restart " << std::endl;
   }
  std::cout << std::endl;

 }

}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//====start_of_problem_class==========================================
/// Problem class -- add segregated solver capability to existing
/// problem.
//====================================================================
template <class ELEMENT>
class SegregatedFSICollapsibleChannelProblem : 
 public virtual FSICollapsibleChannelProblem<ELEMENT>,
 public virtual SegregatableFSIProblem
{
 //Count the number of Picard iterations
 int Picard_iter;
 
 /// Vector of pairs of pointers to wall elements 
 /// (in their incarnation as 
 /// GeomObjects) that contains the control points and the local coordinate
 /// in those objects
 Vector<std::pair<GeomObject*, Vector<double> > > Control_point_pair;
 

public :

 /// Constructor: The arguments are the same as the original
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
 
 /// Empty destructor
 ~SegregatedFSICollapsibleChannelProblem(){}


 /// Identify the fluid and solid Data and meshes that
 /// contain only elements involved in the respective sub-problems. 
 /// This is a specific implementation of a pure virtual function in the 
 /// SegregatableFSIProblem base class.
 void identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                                    Vector<Data*>& solid_data_pt,
                                    Mesh*& fluid_mesh_pt,
                                    Mesh*& solid_mesh_pt);
 
//============================================================================
/// Write header for convergence history file
//============================================================================
 void write_convergence_history_header(std::ofstream& convergence_file, 
                                       const bool& doc_max_global_residual)
  {
   convergence_file 
    << "VARIABLES=\"iteration\",\"CPU time (excluding doc etc.)\",";
   
   unsigned n_control=Control_point_pair.size();
   for (unsigned i=0;i<n_control;i++)
    {
     convergence_file << "\"control point" << i << "\",";
    }
   convergence_file << "\"rms change\","
                    << "\"max change\",\"rms norm\"";
   if (doc_max_global_residual)
    {
     convergence_file << ",\"max global residual\"";
    }
   convergence_file << std::endl;
 
   // Fill in zone header as defined in the derived problem
   // (or call empty version defined in here....)
   write_zone_info_for_convergence_history(convergence_file);

  }


//============================================================================
/// Write to convergence history file
//============================================================================
 void write_convergence_history(
  const unsigned& iter,
  const double& rms_change,
  const double& max_change,
  const double& rms_norm,
  const double& max_res,
  std::ofstream& convergence_file, 
  const bool& doc_max_global_residual)
  {
 
   // Doc iteration number
   convergence_file << iter << " ";
 
   // Doc elapsed cpu time
   convergence_file << t_spent_on_actual_solve()
                    << " ";
 
   // Get wall control points
   unsigned n_control=Control_point_pair.size();
   for (unsigned i=0;i<n_control;i++)
    {
     // Get position
     Vector<double> r_ctrl(2);
     Control_point_pair[i].first->position(Control_point_pair[i].second,
                                           r_ctrl);
     convergence_file << r_ctrl[1] << " ";
    }
 
   // Max. change:
   convergence_file << rms_change << " "
                    << max_change << " "
                    << rms_norm   << " " ;
   if (doc_max_global_residual) convergence_file << max_res << " ";
   convergence_file << std::endl;

  }




 /// Overload empty virtual function that is called before header for
 /// convergence history is written. This can be overloaded to insert
 /// zone information that can be used to identify 
 /// the problem parameters for this solve.
 void write_zone_info_for_convergence_history(ofstream& convergence_file);

 /// Initialise timer and reset counter for Newton iterations
 /// if monolithic solver is used.
 void actions_before_newton_solve()
  {
   if (!Flags::Use_segregated_solver)
    {
     // Initialise counter for Newton iteration
     this->Newton_iter=0;
     
     // Initialise timer that allows doc of iteration/cpu time
     reset_timer();
    }
  }

 void actions_before_segregated_solve()
  {
   //Initialise the counter for the Picard iteration
   Picard_iter = -1;
  }


 /// Overload actions before Newton step: Update nodal 
 /// positions in the fluid mesh in response to any changes in 
 /// the wall displacement field if segregated solver is used.
 void actions_before_newton_step() {}

 /// Overload actions after Newton step: Update nodal 
 /// positions in the fluid mesh
 /// in response to any changes in the wall displacement field. If 
 /// monolithic Newton solver is used, doc progress of Newton iteration, 
 /// using the same output as during Picard iteration.
 void actions_before_newton_convergence_check();

/// Overload actions before the segregated_convergence_check
 /// for documentation purposes
 void actions_before_segregated_convergence_check();

 /// Perform output of wall shape  during after the iter-th 
 /// Picard iteration. This re-implements an empty virtual fct
 /// in the SegregatableFSIproblem base class.
 void solid_output_during_picard(const unsigned& iter,
                                 DocInfo& doc_info);
 
 /// Perform output of fluid flow during after the iter-th 
 /// Picard iteration. This re-implements an empty virtual fct
 /// in the SegregatableFSIproblem base class.
 void fluid_output_during_picard(const unsigned& iter,
                                 DocInfo& doc_info);

 /// Doc the steady solution
 void doc_solution_steady(DocInfo& doc_info,ofstream& trace_file, 
                          const double& cpu, const unsigned& niter);

 /// Doc the unsteady solution
 void doc_solution_unsteady(DocInfo& doc_info,ofstream& trace_file, 
                            const double& cpu, const unsigned& niter);
 
 /// Steady run
 void steady_run();

 /// Unsteady run
  void unsteady_run(const double &dummy_dt=0.1);

private:

 /// Output stream to document the convergence history
 ofstream Convergence_file;

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
 
 //Selet a convergence-acceleration technique based on control flags

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

 //Initialise the number of Picard iterations
 Picard_iter = -1;

 // Doc max. global residual during Picard iteration
 Doc_max_global_residual = true;

 // Number of wall control points
 unsigned n_control=10;
 Control_point_pair.resize(n_control);

 // Get wall control points
 for (unsigned i=0;i<n_control;i++)
  {
   // Get pointer to/local coordinate in wall element that contains 
   // control node 
   Vector<double> zeta_ctrl(1);
   zeta_ctrl[0]=this->Lcollapsible*double(i+1)/double(n_control+1);
   Control_point_pair[i].second.resize(1);
   this->Wall_geom_object_pt->locate_zeta(zeta_ctrl,
                                          Control_point_pair[i].first,
                                          Control_point_pair[i].second);
 
   Vector<double> r_ctrl(2);
   Control_point_pair[i].first->position(Control_point_pair[i].second,
                                         r_ctrl);
  }

} //end_of_constructor





//============================================================================
/// Actions after Newton step: Update nodal positions in the fluid mesh
/// in response to any changes in the wall displacement field. If 
/// monolithic Newton solver is used, doc progress of Newton iteration, 
/// using the same output as during Picard iteration.
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::
actions_before_newton_convergence_check()
{
 //For a "true" segregated solver, we would not do this in fluid or solid
 //solves, but adding the bulk node update to the solid solve phase aids
 //convergence and makes it possible for larger values of Q. Of course,
 //there is a small cost associated with doing this.
 if(Solve_type!=Fluid_solve) {this->Bulk_mesh_pt->node_update();}

 //If we are solving the full problem using a Newton method
 //The default Solve_type is Full_solve, so we will have updated the nodes
 //in the above call
 if (!Flags::Use_segregated_solver)
  {
   //Update the positions
   //this->Bulk_mesh_pt->node_update();

   // Halt timer
   halt_timer();
    
   double rms_change;
   double max_change;
   double rms_norm;
   double max_res=0.0;
   get_solid_change(rms_change,max_change,rms_norm);

   DoubleVector residual;
   get_residuals(residual);
   max_res=residual.max();
   
   std::cout << "==================================================\n";
   std::cout <<   "Iteration             : " 
             << this->Newton_iter << std::endl;
   std::cout <<   "RMS  change           :       "
             << rms_change << std::endl;
   std::cout <<   "Max. change           :       "
             << max_change << std::endl;
   std::cout <<   "RMS norm              :       "
             << rms_norm   << std::endl;
   std::cout << "Max. global residual  :       "
             << max_res   << std::endl;
   std::cout << "==================================================\n\n";

   bool get_max_global_residual=true;
   write_convergence_history(this->Newton_iter,
                             rms_change,
                             max_change,
                             rms_norm,
                             max_res,
                             Convergence_file, 
                             get_max_global_residual);
  
   // Store the current values of the solid dofs as reference values
   store_solid_dofs();
   
   // Increment counter
   this->Newton_iter++;

   // Restart timer
   restart_timer();
  }
}


//============================================================================
/// Actions after Segregated step: Update nodal positions in the fluid mesh
/// in response to any changes in the wall displacement field. If 
/// monolithic Newton solver is used, doc progress of Newton iteration, 
/// using the same output as during Picard iteration.
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::
actions_before_segregated_convergence_check()
{ 
 //Update the nodal positions (not necesary because we do the solid solve
 //last, but doesn't really cost us much.
 this->Bulk_mesh_pt->node_update();
 
 //If we are NOT going to do the Aitken acceleration, then increase
 //the Picard iteration
 if(!Recheck_convergence_after_pointwise_aitken) {++Picard_iter;}

 //Stop timing for the documentation
 halt_timer();
 double rms_change;
 double max_change;
 double rms_norm;
 double max_res=0.0;
 get_solid_change(rms_change,max_change,rms_norm);
 
 bool get_max_global_residual = Doc_max_global_residual;

 if (get_max_global_residual) 
  {
   restore_solid_dofs();
   restore_fluid_dofs(); 
   rebuild_monolithic_mesh();
   assign_eqn_numbers();
   DoubleVector residual;
   get_residuals(residual);
   
   //Get maximum residuals, using our own abscmp function
   max_res =  residual.max();
  }
 
 // Write
 write_convergence_history(Picard_iter,
                           rms_change,
                           max_change,
                           rms_norm,
                           max_res,
                           Convergence_file, 
                           get_max_global_residual);
 
 
 // Restart timer
 restart_timer();
}


//============================================================================
/// Overload empty virtual function that is called before header for
/// convergence history is written. Overloaded to insert
/// zone information that can be used to identify 
/// the problem parameters for this solve.
//============================================================================
template<class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::
write_zone_info_for_convergence_history(ofstream& convergence_file)
{

 convergence_file << "ZONE T=\"";
 convergence_file << "Re=" << Global_Physical_Variables::Re << ", ";
 convergence_file << "Q=" << Global_Physical_Variables::Q << ", ";
 convergence_file << "resolution factor: " << Flags::Resolution_factor << ", ";
 if (Flags::Use_segregated_solver)
  {
   convergence_file << "Picard";
   if (Flags::Use_irons_and_tuck_extrapolation)
    {
     convergence_file << " with Irons and Tuck";
     convergence_file << " with initial under-relaxation factor of "
                      << Flags::Omega_under_relax;
    }
   else
    {
     convergence_file << " with fixed under-relaxation factor of "
                      << Flags::Omega_under_relax;
    }
   convergence_file << "; Tol = "
                    << Flags::Convergence_tolerance;
  }
 else
  {
   convergence_file << "Newton";
  }
 convergence_file << "\""<< std::endl;

}

//=====start_of_identify_fluid_and_solid======================================
/// Identify the fluid and solid Data and the meshes that
/// contain only elements that are involved in the respective sub-problems. 
/// This implements a pure pure virtual function in the 
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


//====start_of_doc_solution_steady============================================
/// Doc the solution for a steady run
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>:: doc_solution_steady(
 DocInfo& doc_info, 
 ofstream& trace_file,
 const double& cpu,
 const unsigned& niter)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

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
 
// Write restart file
 sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 some_file.precision(16);                          
 this->dump_it(some_file);
 some_file.close();

 // Write trace file 
 trace_file << Global_Physical_Variables::P_ext_data_pt->value(0)  << " ";
 trace_file << Global_Physical_Variables::Yprescr  << " ";

 // Get vertical position of control points
 unsigned n_control= Control_point_pair.size();

 // Get wall control points
 for (unsigned i=0;i<n_control;i++)
  {
   std::pair<GeomObject*,Vector<double> > point=Control_point_pair[i];
    
   // Get position
   Vector<double> r_ctrl(2);
   point.first->position(point.second,r_ctrl);
   trace_file << r_ctrl[1] << " ";
  }

 // Write trace file 
 trace_file << this->Left_node_pt->value(0) << " "
            << this->Right_node_pt->value(0) << " "
            << Flags::Use_segregated_solver << " " 
            << Flags::Use_pointwise_aitken << " " 
            << Flags::Omega_under_relax << " " 
            << Flags::Use_irons_and_tuck_extrapolation << " " 
            << cpu << " "
            << niter << " " 
            << std::endl; 

 
} // end_of_doc_solution_steady








//====start_of_doc_solution_unsteady==========================================
/// Doc the solution for an unstady run
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::doc_solution_unsteady(
 DocInfo& doc_info, 
 ofstream& trace_file,
 const double& cpu,
 const unsigned& niter)
{ 

 std::cout << "Doc-ing " << doc_info.number() << std::endl;

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

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
 
// Write restart file
 sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 this->dump_it(some_file);
 some_file.close();

 // Write trace file 
 trace_file << time_pt()->time() << " ";

 // Get/doc y-coordinate of control point
 Vector<double> r(2);
 this->Ctrl_geom_obj_pt->position(this->S_displ_ctrl,r);
 trace_file << r[1]  << " ";

 // Get vertical position of control points
 unsigned n_control= Control_point_pair.size();

 // Get wall control points
 for (unsigned i=0;i<n_control;i++)
  {
   std::pair<GeomObject*,Vector<double> > point=Control_point_pair[i];
    
   // Get position
   Vector<double> r_ctrl(2);
   point.first->position(point.second,r_ctrl);
   trace_file << r_ctrl[1] << " ";
  }

 // Write trace file 
 trace_file << this->Left_node_pt->value(0) << " "
            << this->Right_node_pt->value(0) << " "
            << Flags::Use_segregated_solver << " " 
            << Flags::Use_pointwise_aitken << " " 
            << Flags::Omega_under_relax << " " 
            << Flags::Use_irons_and_tuck_extrapolation << " " 
            << cpu << " "
            << niter << " " 
            << std::endl; 

 
} // end_of_doc_solution_steady






//====steady_run==============================================================
/// Steady run
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
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);



 // Write trace file header
 trace_file << "VARIABLES=\"p<sub>ext</sub>\","
            << "\"y<sub>ctrl</sub>\",";
 unsigned n_control= Control_point_pair.size();
 for (unsigned i=0;i<n_control;i++)
  {
   trace_file << "\"y_" << i << "\",";
  }
 trace_file << "\"u_1\","
            << "\"u_2\","
            << "\"segregated solver flag\","
            << "\"Aitken flag\","
            << "\"Under-relaxation factor\","
            << "\"Irons and Tuck flag\","
            << "\"CPU time for solve\","
            << "\"Number of iterations\"" 
            << std::endl;

 trace_file << "ZONE T=\"";
 trace_file << "Re=" << Global_Physical_Variables::Re << ", ";
 trace_file << "Q=" << Global_Physical_Variables::Q << ", ";
 trace_file << "resolution factor: " << Flags::Resolution_factor << ", ";
 if (Flags::Use_segregated_solver)
  {
   trace_file << "Picard";
   if (Flags::Use_irons_and_tuck_extrapolation)
    {
     trace_file << " with Irons and Tuck";
     trace_file << " with initial under-relaxation factor of "
                << Flags::Omega_under_relax;
    }
   else
    {
     trace_file << " with fixed under-relaxation factor of "
                << Flags::Omega_under_relax;
    }
   trace_file << "; Tol = "
              << Flags::Convergence_tolerance;
  }
 else
  {
   trace_file << "Newton";
  }
 trace_file << "\"" << std::endl;

 // Output the initial solution (dummies for CPU time and # of picard iters
 doc_solution_steady(doc_info, trace_file,0.0,0);
 
 // Increment step number
 doc_info.number()++;
 
 
 // Increment for external pressure (on the wall stiffness scale)
 double delta_p=(Global_Physical_Variables::Pmax-
                 Global_Physical_Variables::Pmin)/double(Flags::Nsteps-1);
 
 // Initial and final values for control position
 Global_Physical_Variables::Yprescr=1.0;
 
 // Steady run: Go down to 0.6. Preparation for unsteady: Go down to 0.65
 double y_min=0.65; // 0.65; // 0.6;
 double delta_y=(y_min-Global_Physical_Variables::Yprescr)/
  double(Flags::Nsteps-1);
 
 // Boolean flag used to specify whether a full setup of solid and fluid dofs
 // is required
 bool full_setup = true;
 
 // Parameter study
 //----------------
 for (unsigned istep=0;istep<Flags::Nsteps;istep++)
  {
   
   // Displacement control?
   if (this->Displ_control)
    {
     std::cout << "Solving for control displ = " 
               << Global_Physical_Variables::Yprescr 
               << std::endl;
    }
   else
    {
     std::cout << "Solving for p_ext = " 
               << Global_Physical_Variables::P_ext_data_pt->value(0) 
               << std::endl;
    }
   
   // Object that stores the convergence data for Picard iteration
   PicardConvergenceData conv_data;

   // Number of iterations taken
   unsigned niter=0; 
   
   // Setup segregated solver in either case
   setup_segregated_solver(full_setup);

   // Solve the problem
   //------------------
   clock_t t_start = clock();


   // SEGREGATED SOLVER
   //==================
   if (Flags::Use_segregated_solver)
    {
     Max_picard =50;
     // Doc convergence history?
     if (doc_info.is_doc_enabled())
      { 
       // Initialise timer that allows doc of iteration/cpu time
       reset_timer();

       // Halt timer before output of reference data
       halt_timer();
   
       char filename[100];
       sprintf(filename,"%s/picard_convergence%i.dat",
               doc_info.directory().c_str(),
               doc_info.number());
       Convergence_file.open(filename);

       bool get_max_global_residual= Doc_max_global_residual;
       write_convergence_history_header(Convergence_file, 
                                        get_max_global_residual);
       // Restart timer after doc
       restart_timer();
      }
     
     // Solve and return converence data
     conv_data=steady_segregated_solve();
       
     niter=conv_data.niter();
     
     Convergence_file.close();

     if (conv_data.cpu_total()!=0.0)
      {
       std::cout 
        << std::endl << std::endl << std::endl 
        << "PERCENTAGE OF CPU TIME SPENT IN COMPUTATION OF GLOBAL RESIDUAL="
        << double(conv_data.cpu_for_global_residual()/
                  conv_data.cpu_total())*100.0
        << std::endl << std::endl << std::endl <<std::endl;
  
       std::cout 
        << "PERCENTAGE OF CPU TIME SPENT IN COMPUTATION OF NON-ESSENTIAL BITS="
        << double((conv_data.cpu_total()-conv_data.essential_cpu_total())/
                  conv_data.cpu_total())*100.0 << " " 
        << conv_data.cpu_total() << " " 
        << conv_data.essential_cpu_total() << " "
        << std::endl << std::endl << std::endl <<std::endl;
      }
       
     
    }

   // NEWTON SOLVER
   //==============
   else
    {
     // Doc convergence history?
     if (doc_info.is_doc_enabled())
      { 
       // Initialise timer that allows doc of iteration/cpu time
       reset_timer();

       // Halt timer before output of reference data
       halt_timer();

       // Store the current values of the solid dofs as reference values
       store_solid_dofs();
   
       char filename[100];
       sprintf(filename,"%s/newton_convergence%i.dat",
               doc_info.directory().c_str(),
               doc_info.number());
       Convergence_file.open(filename);

       bool get_max_global_residual=true;
       write_convergence_history_header(Convergence_file, 
                                        get_max_global_residual);

       // Restart timer after doc
       restart_timer();
      }
     
     // Explit call to the steady Newton solve.
     steady_newton_solve();

     // Number of iterations
     niter=this->Newton_iter;

     Convergence_file.close();
    }
   clock_t t_end= clock();
   double cpu=double(t_end-t_start)/CLOCKS_PER_SEC;
   
   
   
   // Outpt the solution
   //-------------------
   doc_solution_steady(doc_info,trace_file,cpu,niter);
   
   // Step number
   doc_info.number()++;
   
   // Adjust control parameter
   if (this->Displ_control)
    {
     // Increment control position
     Global_Physical_Variables::Yprescr+=delta_y;
    }
   else
    {
     // Increment external pressure 
     double old_p=Global_Physical_Variables::P_ext_data_pt->value(0);
     Global_Physical_Variables::P_ext_data_pt->set_value(0,old_p+delta_p);
    }

   //We no longer need a full setup of the dofs
   full_setup = false;
   
  }
 
 // Close trace file.
 trace_file.close();

}






//====unsteady_run============================================================
/// Unsteady run, dummy argument to be consistent with 
/// definition in FSICollapsibleChannelProblem
//============================================================================
template <class ELEMENT>
void SegregatedFSICollapsibleChannelProblem<ELEMENT>::unsteady_run(
								   const double &dummy_dt)
{ 
  
  // Set initial value for external pressure (on the wall stiffness scale). 
  // Will be overwritten by restart data if a restart file (and pressure
  // jump) are specified
  Global_Physical_Variables::P_ext_data_pt->
   set_value(0,Global_Physical_Variables::Pmax);

  // Timestep
  double dt=Flags::Dt;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  initialise_dt(dt);

  std::cout << "Pressure before set initial: " 
            << Global_Physical_Variables::P_ext_data_pt->value(0)
            << std::endl;

  // Apply initial condition
  set_initial_condition();

  std::cout << "Pressure after set initial: " 
            << Global_Physical_Variables::P_ext_data_pt->value(0)
            << std::endl;

  //Set output directory
  DocInfo doc_info;
  doc_info.set_directory("RESLT");

  // Open a trace file 
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);


 // Write trace file header
 trace_file << "VARIABLES=\"time\","
            << "\"y<sub>ctrl</sub>\",";
 unsigned n_control=Control_point_pair.size();
 for (unsigned i=0;i<n_control;i++)
  {
   trace_file << "\"y_" << i << "\",";
  }
 trace_file << "\"u_1\","
            << "\"u_2\","
            << "\"segregated solver flag\","
            << "\"Aitken flag\","
            << "\"Under-relaxation factor\","
            << "\"Irons and Tuck flag\","
            << "\"CPU time for solve\","
            << "\"Number of iterations\"" 
            << std::endl;

  trace_file << "ZONE T=\"";
  trace_file << "Re=" << Global_Physical_Variables::Re << ", ";
  trace_file << "Q=" << Global_Physical_Variables::Q << ", ";
  trace_file << "resolution factor: " << Flags::Resolution_factor << ", ";
  if (Flags::Use_segregated_solver)
   {
    trace_file << "Picard";
    if (Flags::Use_irons_and_tuck_extrapolation)
     {
      trace_file << " with Irons and Tuck";
      trace_file << " with initial under-relaxation factor of "
                 << Flags::Omega_under_relax;
     }
    else
     {
      trace_file << " with fixed under-relaxation factor of "
                 << Flags::Omega_under_relax;
     }
    trace_file << "; Tol = "
               << Flags::Convergence_tolerance;
   }
  else
   {
    trace_file << "Newton";
   }
  trace_file << "\"" << std::endl;
  
  
  // Output the initial solution (dummies for CPU time and # of picard iters
  doc_solution_unsteady(doc_info, trace_file,0.0,0);

  // Increment step number
  doc_info.number()++;

  // Boolean flag used to specify whether a full setup of solid and fluid dofs
  // is required
  bool full_setup = true;

  // Timestepping loop
  //------------------
  for (unsigned istep=0;istep<Flags::Nsteps;istep++)
   {
 
    // Object that stores the convergence data for Picard iteration
    PicardConvergenceData conv_data;

    // Number iterations taken
    unsigned niter=0; 

   // Setup segregated solver in either case
   setup_segregated_solver(full_setup);
 
    // Solve the problem
    //------------------
    clock_t t_start = clock();
    if (Flags::Use_segregated_solver)
     {
      Max_picard = 50;
      // Doc convergence history?
      if (doc_info.is_doc_enabled())
       { 
        // Initialise timer that allows doc of iteration/cpu time
        reset_timer();
        
        // Halt timer before setup of reference data/doc
        halt_timer();
        
        // Store the current values of the solid dofs as reference values
        //store_solid_dofs();
        
        char filename[100];
        sprintf(filename,"%s/picard_convergence%i.dat",
                doc_info.directory().c_str(),
                doc_info.number());
        Convergence_file.open(filename);
        
        bool get_max_global_residual=Doc_max_global_residual;
        write_convergence_history_header(Convergence_file, 
                                         get_max_global_residual);
        
        // Restart timer after setup/doc of reference data
        restart_timer();
        }
      
      // Solve and return convergence data
      conv_data=unsteady_segregated_solve(dt);

      niter=conv_data.niter();
      
      Convergence_file.close();
      
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
     }
    else
     {
      
      // Doc convergence history?
      if (doc_info.is_doc_enabled())
       { 
        // Initialise timer that allows doc of iteration/cpu time
        reset_timer();
        
        // Halt timer before setup of reference data/doc
        halt_timer();
        
        // Store the current values of the solid dofs as reference values
        store_solid_dofs();
        
        char filename[100];
        sprintf(filename,"%s/newton_convergence%i.dat",
                doc_info.directory().c_str(),
                doc_info.number());
        Convergence_file.open(filename);
        
        bool get_max_global_residual=true;
        write_convergence_history_header(Convergence_file, 
                                         get_max_global_residual);
        
        // Restart timer after setup/doc of reference data
        restart_timer();
       }
      
      // Explit call to the unsteady Newton solve.
      unsteady_newton_solve(dt);
      
      // Number of Newton iterations
      niter=this->Newton_iter;
      
      Convergence_file.close();
      
     }
    
    clock_t t_end= clock();
    double cpu=double(t_end-t_start)/CLOCKS_PER_SEC;
    
 
    // Output the solution
    //--------------------
    doc_solution_unsteady(doc_info,trace_file,cpu,niter);
 
    // Step number
    doc_info.number()++;
 
    // Only keep the last ten solutions unless we're using the
    // Newton solver
    if (Flags::Use_segregated_solver==1)
     {
      if (doc_info.number()==10) doc_info.number()-=10;
      std::cout << "Resetting doc file numbers" << std::endl;
     }
    
    
    //We no longer need a full setup of the dofs
    full_setup = false;
   
   }

  // Close trace file.
  trace_file.close();

}



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////





//============start_of_main====================================================
/// Driver code for a collapsible channel problem with FSI.
/// Presence of command line arguments indicates validation run with 
/// coarse resolution and small number of steps.
//=============================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 if (CommandLineArgs::Argc==1)
  {
   std::cout << "Using default settings for flags" << std::endl;
  }
 else if ((CommandLineArgs::Argc==15)||(CommandLineArgs::Argc==17))
  {
   /// Resolution factor
   Flags::Resolution_factor=atoi(argv[1]);

   /// Use Newton solver (0) or segregated solver (1)?
   Flags::Use_segregated_solver=atoi(argv[2]);

   /// Use pointwise Aitken extrapolation (1) or not (0)
   Flags::Use_pointwise_aitken=atoi(argv[3]);

   /// Under-relaxation parameter
   Flags::Omega_under_relax=atof(argv[4]);

   /// Use Irons and Tuck extrapolation
   Flags::Use_irons_and_tuck_extrapolation=atoi(argv[5]);

   /// Use displacement control (1) or not (0)
   Flags::Use_displ_control=atoi(argv[6]);

   /// Min. y coordinate for parameter study with displacement control
   Global_Physical_Variables::Yprescr_min=double(atof(argv[7]));

   /// Steady (1) or unsteady (0) run
   Flags::Steady_flag=atoi(argv[8]);

   /// Convergence criterion: 0: global resmax; 1: abs. change; 2: rel. change
   Flags::Convergence_criterion=atoi(argv[9]);

   /// Convergence tolerance
   Flags::Convergence_tolerance=atof(argv[10]);

   /// Number of steps
   Flags::Nsteps=atoi(argv[11]);

   /// Reynolds number
   Global_Physical_Variables::Re=double(atof(argv[12]));

   /// Womersley number
   Global_Physical_Variables::ReSt=Global_Physical_Variables::Re;

   /// FSI parameter
   Global_Physical_Variables::Q=double(atof(argv[13]));

   /// Timestep
   Flags::Dt=double(atof(argv[14]));

   // Restart?
   if (CommandLineArgs::Argc==17)
    {
     // Name of restart file
     Flags::Restart_file_name=argv[15];
     
     // Jump in pressure 
     Global_Physical_Variables::P_step=double(atof(argv[16]));
    }
  }
 else
  {
   std::cout 
    << "\n\n\n\n\n"
    << "Wrong number of command line args: Specify none, fourteen or sixteen:"
    << std::endl
    << "- resolution factor [1,2,...]" << std::endl
    << "- use_segregated_solver [0/1]"<< std::endl
    << "- use_pointwise_aitken [0/1]" << std::endl
    << "- under-relaxation parameter" << std::endl
    << "- use_irons_and_tuck_extrapolation [0/1]" << std::endl
    << "- use_displ_control [0/1]" << std::endl
    << "- min. y-coordinate of control point when using displ control" 
    <<    std::endl
    << "- steady_flag [0/1]" << std::endl
    << "- convergence criterion [0/1/2]" << std::endl
    << "- convergence tolerance " << std::endl
    << "- number of steps " << std::endl
    << "- Reynolds number" << std::endl
    << "- FSI parameter Q" << std::endl
    << "- Timestep dt" << std::endl
    << "- restart file name [optional] " << std::endl
    << "- jump in pressure P_step [optional] " << std::endl
    << "You specified " << CommandLineArgs::Argc-1 << " command line arg[s]" 
    << "\n\n\n\n\n" << std::endl;
   abort();
  }
 Flags::doc_flags();
 Flags::segregated_doc_flags();

 
 // Number of elements in the domain
 unsigned nup=4*Flags::Resolution_factor;
 unsigned ncollapsible=20*Flags::Resolution_factor;
 unsigned ndown=40*Flags::Resolution_factor;
 unsigned ny=4*Flags::Resolution_factor;
  
 
 // Length of the domain
 double lup=1.0;
 double lcollapsible=5.0;
 double ldown=10.0;
 double ly=1.0;
 

 // Use displacement control?
 bool displ_control=false;
 if (Flags::Use_displ_control==1) displ_control=true;
 
 // Steady run?
 bool steady_flag=false;
 if (Flags::Steady_flag==1) steady_flag=true;

 // Build the problem with QTaylorHoodElements
 SegregatedFSICollapsibleChannelProblem
  <AlgebraicElement<QTaylorHoodElement<2> > > 
 problem(nup, ncollapsible, ndown, ny, 
         lup, lcollapsible, ldown, ly, displ_control,
         steady_flag);
  
 if (Flags::Steady_flag)
  {
   problem.steady_run();
  }
 else
  {
   problem.unsteady_run();
  }
 
}//end of main

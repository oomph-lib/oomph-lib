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
#include <fenv.h> 

//Generic routines
#include "generic.h"

// The equations
#include "generalised_newtonian_axisym_navier_stokes.h" 
//#include "num_rec.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;


namespace oomph
{
//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    

  /// Which bdf timestepper do we use?
  unsigned BDF_type=2;

  /// Number of previous timesteps to be used for
  /// extrapolation of strain rate (cannot be bigger than  BDF_type)
  unsigned Nprev_for_extrapolation_of_strain_rate=2;

  /// Uniform initial element area
  double Uniform_element_area=0.004;

  /// Output directory
  std::string Directory = "RESLT";
  
  /// Name of restart file
  std::string Restart_file="";

  /// Doc info trace file
  DocInfo Doc_info_trace;

  /// Doc info solutions
  DocInfo Doc_info_soln;

  /// Set external pressure
  double External_pressure=0.0;

  /// Reynolds number
  double Re=1.967e-1;

  // Reynolds inverse Froude number
  double ReInvFr=1.0;

  /// Capillary number
  double Ca = 42.19;

  /// Strouhal number
  double St = 1.440e0;

  /// Reynolds multiplied Strouhal
  double Re_St = Re*St;

  /// Maximum force of the applied oscillation (non-dimensionalised)
  double F_max = 40.24;

  /// Body force multiplier
  double Oscillating_body_force_multiplier=1.0;

  /// Time after which oscillation starts
  double Initial_settling_time = 0.25;

  /// Minimum allowable timestep
  double Dt_min=1.0e-6;

  /// Function describing the oscillation
  void oscillation(const double& time, const Vector<double>& x, 
                   Vector<double>& force)
  {
   force[0]=0.0;
   
   // Start the oscillation after initial settling time
   if(time >= Initial_settling_time)
    {
     force[1]=-1.0 +  Oscillating_body_force_multiplier*
      F_max*sin(2.0*MathematicalConstants::Pi*(time-Initial_settling_time));
    }
   else
    {
     force[1]=-1.0;
    }

   // Horizontal body force is always zero
   force[2]=0.0;

  }

  // Free surface profile
  double free_surface_profile(const double r)
  {
   return -0.05*(cos(MathematicalConstants::Pi*r)+1.0);
  }

  /// Pseudo-solid Poisson ratio
  double Nu=0.3;

  /// Pseudo-solid Mooney-Rivlin parameter
  double C1=1.0;

  /// Pseudo-solid Young's modulus
  double E=2.2;

  /// Radius of the shell
  double Shell_radius=1.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=0;

  /// have a generic yield stress, required for the output function
  double Yield_stress = 2.861e-1;

  // Casson model
  //-------------------------------------------------------------

  /// Yield stress
  double C_yield_stress=2.526e-1;

  /// Reynolds number
  double C_Re=7.827e-1;

  /// Strouhal number
  double C_St = 2.966e-1;

  /// Reynolds multiplied Strouhal
  double C_Re_St = C_Re*C_St;
 
  //-------------------------------------------------------------

  // Herschel-Bulkley model
  //-------------------------------------------------------------

  /// Yield stress
  double HB_yield_stress=2.861e-1;
 
  /// Power law exponent
  double HB_flow_index = 0.828;

  /// Reynolds number
  double HB_Re=2.430e-1;

  /// Strouhal number
  double HB_St =1.295e0;

  /// Reynolds multiplied Strouhal
  double HB_Re_St = HB_Re*HB_St;
 
  //-------------------------------------------------------------

  // Sisko model
  //-------------------------------------------------------------

  /// Pre-factor in Sisko model K / ( mu_inf T^(n-1) )
  double S_alpha=3.181e-1;

  /// Flow index in Sisko model
  double S_flow_index=4.89e-6;

  /// Reynolds number
  double S_Re=1.967e-1;

  /// Strouhal number
  double S_St = 1.440e0;

  /// Reynolds multiplied Strouhal
  double S_Re_St = S_Re*S_St;

  // hierher rename to cutoff for second invariant below which the
  // viscosity remains constant. 

  /// Critical strain rate
  double Critical_strain_rate=1.4e-14;

  /// Fluid constitutive equation
  GeneralisedNewtonianConstitutiveEquation<3>* Const_eqn_pt=0;

  /// Trace file
  ofstream Trace_file;

  /// File to document the norm of the solution (for validation 
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

 } // end_of_namespace


//==============================================================
/// Overload TaylorHood element to modify output
//==============================================================
 class MyTaylorHoodElement : 
  public virtual PseudoSolidNodeUpdateElement
  <GeneralisedNewtonianAxisymmetricTTaylorHoodElement, TPVDElement<2,3> >
 {

#include "overloaded_element_body.h"
  
 };



//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<MyTaylorHoodElement>
  : public virtual SolidTElement<1,3> 
 {
 public:
  FaceGeometry() : SolidTElement<1,3>() {}


 };


//=======================================================================
/// Face geometry of Face geometry for element is the same 
/// as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<FaceGeometry<MyTaylorHoodElement> >
  : public virtual SolidPointElement 
 {
 public:
  FaceGeometry() : SolidPointElement() {}
 };

} //End of namespace extension


/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Problem class to simulate the settling of a viscous drop
//====================================================================
template<class ELEMENT>
class AxisymmetricVibratingShellProblem : public Problem
{

public:

 /// Constructor
 AxisymmetricVibratingShellProblem();
 
 /// Destructor
 ~AxisymmetricVibratingShellProblem();
 
 /// Actions before adapt: Wipe the mesh of free surface elements
 void actions_before_adapt()
  {
   // Kill the  elements and wipe surface mesh
   delete_free_surface_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  
  }// end of actions_before_adapt

 
 /// Actions after adapt: Rebuild the mesh of free surface elements
 void actions_after_adapt()
  {

   // Create the elements that impose the displacement constraint 
   create_free_surface_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();

   // Setup the problem again -- remember that fluid mesh has been
   // completely rebuilt and its element's don't have any
   // pointers to Re etc. yet
   complete_problem_setup();

  }// end of actions_after_adapt

 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve
 void actions_before_newton_solve()
  {
   //Reset the Lagrangian coordinates of the nodes to be the current
   //Eulerian coordinates -- this makes the current configuration
   //stress free
   Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  }


 /// Allow inverted elements during Newton iteration
 void actions_before_newton_convergence_check()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

     // accept negative jacobian in solution process
     el_pt->Accept_negative_jacobian=true;
    }

  }

 /// Don't allow inverted  elements in general
 void actions_after_newton_convergence_check()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

     el_pt->Accept_negative_jacobian=false;
    }

  }

 /// Set boundary conditions and complete the build of all elements
 void complete_problem_setup();

 /// Doc the solution
 void doc_solution(const std::string& comment="");

 /// Compute the error estimates and assign to elements for plotting
 void compute_error_estimate(double& max_err,double& min_err);

 /// Error norm to determine the next time step
 double global_temporal_error_norm();
 
 /// Access the next suggested timestep
 double& next_dt() {return Next_dt;}
 
 /// Get the height at the centre node
 double height_central_node()
  {
   return Central_node_on_free_surface_pt->x(1);
  }

 /// Update latest guess for strain rate
 void update_latest_fixed_point_iteration_guess_for_strain_rate_for_all_elements()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     el_pt->update_latest_fixed_point_iteration_guess_for_strain_rate();
    }
  }

 /// Enable use of fixed point iteration
 /// for all elements
  void enable_fixed_point_iteration_for_strain_rate_for_all_elements()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     el_pt->enable_fixed_point_iteration_for_strain_rate();
   }
  }

  /// Disable use of fixed point iteration
  void disable_fixed_point_iteration_for_strain_rate_for_all_elements()
   {
    unsigned n_element = Fluid_mesh_pt->nelement();
    for(unsigned e=0;e<n_element;e++)
     {
      // Upcast from GeneralisedElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
      el_pt->disable_fixed_point_iteration_for_strain_rate();
     }
   }

 /// Enable use of Aitken extrapolation for all elements
 void enable_aitken_extrapolation_for_all_elements()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     el_pt->enable_aitken_extrapolation();
   }
  }

 /// Disable use of Aitken extrapolation
 void disable_aitken_extrapolation_for_all_elements()
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     el_pt->disable_aitken_extrapolation();
    }
  }

 // Set up extrapolation stuff
 void set_nprev_for_extrapolation_of_strain_rate_for_all_elements(
  const unsigned& nprev)
  {
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
     
     // Use extrapolated strain rate when determining viscosity?  
     if (!CommandLineArgs::command_line_flag_has_been_set
         ("--use_current_strainrate_for_viscosity"))
      {
       el_pt->use_extrapolated_strainrate_to_compute_second_invariant();
       
       // Set extrapolation order
       el_pt->nprev_for_extrapolation_of_strain_rate()=nprev;
      }
     else
      {
       el_pt->use_current_strainrate_to_compute_second_invariant();
      }
    }
  }

 /// get the error of the fixed point iteration
 double calculate_error_of_fixed_point_iteration()
  {
   // Get elemental max/min/norms
   double norm_squared=0.0;
   double latest_guess_norm_squared=0.0;
   double error_norm_squared=0.0;

   unsigned nel=Fluid_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     // Get element
     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
      Fluid_mesh_pt->element_pt(e));

     // Get norms of invariant
     double el_norm_squared=0.0;
     double el_latest_guess_norm_squared=0.0;
     double el_error_norm_squared=0.0;
     el_pt->square_of_norm_of_fixed_point(
      el_norm_squared,
      el_latest_guess_norm_squared,
      el_error_norm_squared);

     // Add it...
     norm_squared+=el_norm_squared;
     latest_guess_norm_squared+=el_latest_guess_norm_squared;
     error_norm_squared+=el_error_norm_squared;
    }
   
   oomph_info << "Norm of current strain rate invariant: "
              << sqrt(norm_squared) << std::endl;
   oomph_info << "Norm of latest fixed point iteration guess for "
              << "strain rate invariant: "
              << sqrt(latest_guess_norm_squared) << std::endl;
   oomph_info << "Norm of error in fixed point iteration "
              << "strain rate invariant: "
              << sqrt(error_norm_squared) << " equivalent to " ;
   if (sqrt(norm_squared)!=0.0)
    {
     oomph_info << sqrt(error_norm_squared)/sqrt(norm_squared)*100.0 << " %";
    }
   oomph_info << std::endl;

   if (sqrt(norm_squared)!=0.0)
    {
     return sqrt(error_norm_squared)/sqrt(norm_squared)*100.0;
    }
   else
    {
     return 0.0;
    }

  }

 // error in extrapolation
 double strainrate_norm()
  {
   // Compute current volume and error in extrapolation from fluid elements
   double norm_squared=0.0;
   double extrapolated_norm_squared=0.0;
   double error_norm_squared=0.0;
   unsigned nel=Fluid_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     // Get element
     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
      Fluid_mesh_pt->element_pt(e));

     // Get norms of invariant
     double el_norm_squared=0.0;
     double el_extrapolated_norm_squared=0.0;
     double el_error_norm_squared=0.0;
     double test_size=0.0;
     test_size=el_pt->square_of_norm_of_strain_invariant(
      el_norm_squared,
      el_extrapolated_norm_squared,
      el_error_norm_squared);
     
     // Add it...
     norm_squared+=el_norm_squared;
     extrapolated_norm_squared+=el_extrapolated_norm_squared;
     error_norm_squared+=el_error_norm_squared;
    }

   return sqrt(norm_squared);
  }

 /// Dump problem data to allow for later restart
 void dump_it(ofstream& dump_file)
  {

   // Write doc numbers
   dump_file << Problem_Parameter::Doc_info_trace.number() 
             << " # current doc number for trace" << std::endl;
   dump_file << Problem_Parameter::Doc_info_soln.number()  
             << " # current doc number for soln" << std::endl;

   // Next timestep (required for restart from temporally adaptive
   // run
   dump_file << Next_dt << " # next timestep " << std::endl;

   // Dump the refinement pattern and the generic problem data
   Problem::dump(dump_file);
  }


 /// Restart
 void restart()
  {
   // Pointer to restart file
   ifstream* restart_file_pt=0;
   
   // Open restart file from stem
   restart_file_pt=new ifstream(Problem_Parameter::Restart_file.c_str(),
                                ios_base::in);
   if (restart_file_pt!=0)
    {
     oomph_info << "Have opened "
                << Problem_Parameter::Restart_file.c_str() 
                << " for restart. " << std::endl;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream
      << "ERROR while trying to open " 
      << Problem_Parameter::Restart_file.c_str()
      << " for restart." << std::endl;
     
     throw OomphLibError(
      error_stream.str(),
      "restart()",
      OOMPH_EXCEPTION_LOCATION);
    }
   
   
   // Read restart data:
   //-------------------
   if (restart_file_pt!=0)
    { 

     // Doc number for trace
     //---------------------
     // Read line up to termination sign
     string input_string;
     getline(*restart_file_pt,input_string,'#');
     
     // Ignore rest of line
     restart_file_pt->ignore(80,'\n');
     
     // Doc number
     Problem_Parameter::Doc_info_trace.number()=unsigned(atoi(input_string.c_str()));


     // Doc number for solution
     //---------------------
     // Read line up to termination sign
     getline(*restart_file_pt,input_string,'#');
     
     // Ignore rest of line
     restart_file_pt->ignore(80,'\n');
     
     // Doc number
     Problem_Parameter::Doc_info_soln.number()=unsigned(atoi(input_string.c_str()));


     // Next timestep (required for restart from temporally adaptive run
     //-----------------------------------------------------------------

     // Read line up to termination sign
     getline(*restart_file_pt,input_string,'#');
     
     // Ignore rest of line
     restart_file_pt->ignore(80,'\n');
     
     // Next timestep
     Next_dt=double(atof(input_string.c_str()));
     

     // Refine the mesh and read in the generic problem data
     Problem::read(*restart_file_pt);
    }
  }


private:

 /// Pointer to free surface node on symmetry line
 Node* Central_node_on_free_surface_pt;

 /// Suggestion for the next timestep (allows it to be written to
 /// or read from a restart file)
 double Next_dt;

 /// Create free surface elements
 void create_free_surface_elements();

 /// Delete free surface elements 
 void delete_free_surface_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Free_surface_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Free_surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Free_surface_mesh_pt->flush_element_and_node_storage();
   
  } // end of delete_free_surface_elements

 /// Pointers to mesh of free surface elements
 Mesh* Free_surface_mesh_pt;

 /// Pointer to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;

 /// Pointer to a global external pressure datum
 Data* External_pressure_data_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshClosedCurve* Outer_boundary_polyline_pt; 

 Vector<TriangleMeshOpenCurve* > Internal_open_boundary_pt;

 /// Enumeration of boundaries
 enum 
 {
  Internal_boundary1_id=0,
  Free_surface_boundary_id=1,
  Symmetry_line_id=2,
  Shell_wall_boundary_id=3
 };
 
 /// External pressure
 double Pext;

}; // end_of_problem_class

//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
AxisymmetricVibratingShellProblem<ELEMENT>::AxisymmetricVibratingShellProblem()
{ 

 Always_take_one_newton_step = true;
 Minimum_dt_but_still_proceed=Problem_Parameter::Dt_min;
 Max_residuals=10.0;
 Max_newton_iterations=25;

 //enable_globally_convergent_newton_method();

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps.
 switch (Problem_Parameter::BDF_type)
  {
  case 1:
   this->add_time_stepper_pt(new BDF<1>(true));
   oomph_info << "Using BDF1\n";
   break;

  case 2: 
   this->add_time_stepper_pt(new BDF<2>(true));
   oomph_info << "Using BDF2\n";
   break;

  case 4:
   this->add_time_stepper_pt(new BDF<4>(true));
   oomph_info << "Using BDF4\n";
   break;

  default:
   oomph_info << "Wrong BDF type: " << Problem_Parameter::BDF_type 
              << std::endl;
   break;
  }

  // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // 2 separate polylines
 //------------------------

 Vector<TriangleMeshCurveSection*> boundary_curve_section_pt(3);
 Vector<TriangleMeshPolyLine*> internal_polyline_pt(2);
 
 // Each polyline only has two vertices -- provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord;

 Ellipse* Shell_wall_pt = new Ellipse(Problem_Parameter::Shell_radius,
                                      Problem_Parameter::Shell_radius);

 
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi/2.0;

 // number of points along the shell wall
 unsigned npoints_wall=80;

 // get the increment in boundary coordinate
 double unit_zeta = (zeta_end-zeta_start)/double(npoints_wall-1);

 // resize the vector storing the vertices
 vertex_coord.resize(npoints_wall);

 // boundary coordinate
 Vector<double> zeta(1,0.0);

 // coordinates of point on the boundary
 Vector<double> coord(2,0.0);

 // Create points on boundary
 for(unsigned ipoint=0; ipoint<npoints_wall;ipoint++)
  {
   // Get the coordinates
   zeta[0]=zeta_start+unit_zeta*double(ipoint);
   Shell_wall_pt->position(zeta,coord);

   // resize vertex coordinate
   vertex_coord[ipoint].resize(2);

   // Shift
   vertex_coord[ipoint][0]=coord[0];
   //std::cout<<vertex_coord[ipoint][0]<< " ";
   vertex_coord[ipoint][1]=coord[1];
   //std::cout<<vertex_coord[ipoint][1]<< std::endl;
  }

 boundary_curve_section_pt[0]= 
  new TriangleMeshPolyLine(vertex_coord, Shell_wall_boundary_id);

 ofstream file;
 file.open((Problem_Parameter::Directory+"/boundary_section0.dat").c_str());
 boundary_curve_section_pt[0]->output(file);
 file.close();

 // Resize the vertex coordinates vector
 vertex_coord.resize(2);

 for(unsigned i=0;i<2;i++)
  {
   vertex_coord[i].resize(2);
  }

 // Build the symmetry polyline
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=Problem_Parameter::Shell_radius;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=Problem_Parameter::free_surface_profile(0.0); // obacht
 boundary_curve_section_pt[1]= new TriangleMeshPolyLine(vertex_coord,
                                                        Symmetry_line_id);

 file.open((Problem_Parameter::Directory+"/boundary_section1.dat").c_str());
 boundary_curve_section_pt[1]->output(file);
 file.close();

 // number of points along the shell wall
 unsigned npoints_fs=80;

 // get the increment in boundary coordinate
 double increment = 1.0/double(npoints_fs-1);

 // resize the vector storing the vertices
 vertex_coord.resize(npoints_fs);

 // Create points on boundary
 for(unsigned ipoint=0; ipoint<npoints_fs;ipoint++)
  {
   // resize vertex coordinate
   vertex_coord[ipoint].resize(2);

   if(ipoint==0)
    {
     vertex_coord[ipoint][0]=0.0;
     //std::cout<<vertex_coord[ipoint][0]<< " ";
     vertex_coord[ipoint][1]=
      Problem_Parameter::free_surface_profile(0.0); // obacht
     //std::cout<<vertex_coord[ipoint][1]<< std::endl;
     continue;
    }
   else if(ipoint==npoints_fs-1)
    {
     vertex_coord[ipoint][0]=Problem_Parameter::Shell_radius;
     //std::cout<<vertex_coord[ipoint][0]<< " ";
     vertex_coord[ipoint][1]=0.0;
     //std::cout<<vertex_coord[ipoint][1]<< std::endl;
     continue;
    }

   // Get the coordinates
   coord[0]=double(ipoint)*increment;
   coord[1]=Problem_Parameter::free_surface_profile(coord[0]);
   //-1.0*(-0.05*tanh(10.0*coord[0]-7.5)+0.05); // obacht

   // Shift
   vertex_coord[ipoint][0]=coord[0];
   //std::cout<<vertex_coord[ipoint][0]<< " ";
   vertex_coord[ipoint][1]=coord[1];
   //std::cout<<vertex_coord[ipoint][1]<< std::endl;
  }

 // Build the boundary polyline
 boundary_curve_section_pt[2]= 
  new TriangleMeshPolyLine(vertex_coord, Free_surface_boundary_id);

 file.open((Problem_Parameter::Directory+"/boundary_section2.dat").c_str());
 boundary_curve_section_pt[2]->output(file);
 file.close();

 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = 
  new TriangleMeshClosedCurve(boundary_curve_section_pt);

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   Outer_boundary_polyline_pt);

 // Define the maximum element areas
 triangle_mesh_parameters.element_area() =
  Problem_Parameter::Uniform_element_area;

 // enable the boundary refinement
 triangle_mesh_parameters.enable_boundary_refinement();

 // Create the mesh
 Fluid_mesh_pt =
   new RefineableSolidTriangleMesh<ELEMENT>(
    triangle_mesh_parameters,this->time_stepper_pt());


 // Increase number of bins to reduce bias/drift in area targets
 Fluid_mesh_pt->nbin_x_for_area_transfer()=500;
 Fluid_mesh_pt->nbin_y_for_area_transfer()=500;

// Fluid_mesh_pt->output((Problem_Parameter::Directory+"/Mesh_before_snapping.dat").c_str());

 //----------------------------------------------------------------
 // Snap nodes manually onto the curved boundary
 //----------------------------------------------------------------

 // loop over the nodes on the wall boundary
 unsigned n_boundary_node = Fluid_mesh_pt->
  nboundary_node(Shell_wall_boundary_id);

 for(unsigned n=0;n<n_boundary_node;n++)
  {
   /// Get pointer to the free surface node
   Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(Shell_wall_boundary_id,n);

   /// Get the node's boundary coordinate
   Vector<double> zeta(1);
   nod_pt->get_coordinates_on_boundary(Shell_wall_boundary_id, zeta);

   /// Get the node's coordinate in the spline representation
   Vector<double> new_x(2,0.0);
   // updating zeta
   zeta[0]=zeta_start+zeta[0]*(zeta_end-zeta_start);
   Shell_wall_pt->position(zeta,new_x);
   nod_pt->x(0) = new_x[0];
   nod_pt->x(1) = new_x[1];

  }

 n_boundary_node=Fluid_mesh_pt->nboundary_node(Free_surface_boundary_id);

 for(unsigned n=0;n<n_boundary_node;n++)
  {
   /// Get pointer to the free surface node
   Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id,n);

   if(nod_pt->is_on_boundary(Shell_wall_boundary_id) ||
      nod_pt->is_on_boundary(Symmetry_line_id))
    {
     continue;
    }

   /// Get the node's r-coordinate
   double r=nod_pt->x(0);
   
   // calculate node's z-coordinate
   double new_z=Problem_Parameter::free_surface_profile(r); // obacht

   // move node
   nod_pt->x(1) = new_z;

  }

 //-----------------------------------------------------------------
 //   End of snapping
 //-----------------------------------------------------------------

 //Fluid_mesh_pt->output((Problem_Parameter::Directory+"/Mesh_after_snapping.dat").c_str());
 Fluid_mesh_pt->output_boundaries((Problem_Parameter::Directory+"/boundaries.dat").c_str());

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=5.0e-3; // was 0.005; 0.0005
 Fluid_mesh_pt->min_permitted_error()=1.0e-3; // was 0.001; 0.0001
 Fluid_mesh_pt->max_element_size()=Problem_Parameter::Uniform_element_area;
 Fluid_mesh_pt->min_element_size()=1.0e-5; // was 0.001 or 0.00001


 /// Set external pressure to zero
 Pext=0.0;

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Pext);

 // The external pressure is pinned -- the external pressure
 // sets the pressure throughout the domain -- we do not have
 // the liberty to fix another pressure value!
 External_pressure_data_pt->pin(0);

 // Construct the mesh of free surface elements
 Free_surface_mesh_pt=new Mesh;
 create_free_surface_elements();

 // Set boundary condition and complete the build of all bulk elements
 complete_problem_setup();

 // Combine meshes
 //---------------

 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Free_surface sub meshes
 this->add_sub_mesh(this->Free_surface_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();

 // Set lagrangian coordinates for pseudo-solid
 Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

 // Use mumps
 //linear_solver_pt()=new MumpsSolver;

 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
     
 
} // end_of_constructor


//==start_destructor======================================================
/// Destructor
//========================================================================
template<class ELEMENT>
AxisymmetricVibratingShellProblem<ELEMENT>::~AxisymmetricVibratingShellProblem()
{
  // Fluid timestepper
  delete this->time_stepper_pt(0);
  
  // Kill data associated with outer boundary
  unsigned n=Outer_boundary_polyline_pt->ncurve_section();
  for (unsigned j=0;j<n;j++)
   {
    delete Outer_boundary_polyline_pt->curve_section_pt(j);
   }
  delete Outer_boundary_polyline_pt;
  
  

  // Flush element of free surface elements
  delete_free_surface_elements();
  delete Free_surface_mesh_pt;
  
  // Delete error estimator
  delete Fluid_mesh_pt->spatial_error_estimator_pt();
  
  // Delete fluid mesh
  delete Fluid_mesh_pt;

  // Delete pressure data
  delete External_pressure_data_pt;
  
  // Kill const eqn
  delete Problem_Parameter::Constitutive_law_pt;

  delete Problem_Parameter::Const_eqn_pt;
  
} // end_of_destructor


//=============start_of_complete_problem_setup===========================
/// Set up the problem: apply BC and make bulk elements fully functional
//=======================================================================
template<class ELEMENT>
void AxisymmetricVibratingShellProblem<ELEMENT>::complete_problem_setup()
{ 
 // Re-set the boundary conditions for fluid problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {  
     // Get node
     Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
     
     //Pin both velocities at the shell wall boundary
     if(ibound==Shell_wall_boundary_id)
      {
       nod_pt->pin(0);
       nod_pt->pin(1);
      }
     
     // pin horizontal velocity at symmetry boundary
     if(ibound==Symmetry_line_id)
      {
       nod_pt->pin(0);
      }
     
     // pin Lagrange multiplier at the intersection of the shell wall
     // boundary and the free surface
     if( (nod_pt->is_on_boundary(Shell_wall_boundary_id)) &&
         (nod_pt->is_on_boundary(Free_surface_boundary_id)) )
      {
       // Get the number of values at this node
       unsigned n_value=nod_pt->nvalue();

       // check that it is the corner node with 5 values
       // (u,v,w,p,delta)
       if(n_value != 5) 
        {
         oomph_info <<" Here!\n";
         abort();
        }

       nod_pt->pin(n_value-1);
      }
     
     // Pin pseudo-solid positions apart from free surface boundary which
     // we allow to move
     SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
     if(ibound==Shell_wall_boundary_id)
      {
       solid_node_pt->pin_position(0);
       solid_node_pt->pin_position(1);
      }
     else if(ibound==Symmetry_line_id)
      {
       solid_node_pt->pin_position(0);
      }

    }
   
  } // end loop over boundaries
 
 // Complete the build of all elements so they are fully functional
 // Remember that adaptation for triangle meshes involves a complete
 // regneration of the mesh (rather than splitting as in tree-based
 // meshes where such parameters can be passed down from the father
 // element!)
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

   // Modify the tolerance for when the mapping is considered singular
   el_pt->Tolerance_for_singular_jacobian=1.0e-26;

   // Set the Reynolds number
   el_pt->re_pt() = &Problem_Parameter::Re;
   
   // Set the Womersley number
   el_pt->re_st_pt() = &Problem_Parameter::Re_St;
   
   // Set the body force
   el_pt->axi_nst_body_force_fct_pt() = &Problem_Parameter::oscillation;

   // Set the product of Reynolds number and inverse Froude number
   el_pt->re_invfr_pt() = &Problem_Parameter::ReInvFr;
   
   // Set the constitutive law for pseudo-elastic mesh deformation
   el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;

   // Assign Constitutive equation
   el_pt->constitutive_eqn_pt() = Problem_Parameter::Const_eqn_pt;

   // Use extrapolated strain rate when determining viscosity
   el_pt->use_extrapolated_strainrate_to_compute_second_invariant();

   // for all nodes, pin azimuthal velocity (fewer dofs)
   unsigned n_node=el_pt->nnode();
   for(unsigned inod=0;inod<n_node;inod++)
    {
     el_pt->node_pt(inod)->pin(2);
     el_pt->node_pt(inod)->set_value(2,0.0);
    }

  } // end of functional elements
 

 // Setup extrapolation
 set_nprev_for_extrapolation_of_strain_rate_for_all_elements(
  Problem_Parameter::Nprev_for_extrapolation_of_strain_rate);
      
 
 // Re-apply boundary values on Dirichlet boundary conditions 
 // (Boundary conditions are ignored when the solution is transferred
 // from the old to the new mesh by projection; this leads to a slight
 // change in the boundary values (which are, of course, never changed,
 // unlike the actual unknowns for which the projected values only
 // serve as an initial guess)
 
 // Set velocity and history values of velocity on walls
 nbound=this->Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;++ibound)
  {
   if( (ibound==Shell_wall_boundary_id) ||
       (ibound==Symmetry_line_id) )
    {
     // Loop over nodes on this boundary
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       
       // Velocity is and was zero at all previous times
       for (unsigned t=0;t<=n_prev;t++)
        {
         if(ibound==Shell_wall_boundary_id)
          {
           nod_pt->set_value(t,1,0.0);
           nod_pt->set_value(t,0,0.0);
          }
         else
          {
           nod_pt->set_value(t,0,0.0);
          }

         // Move nodes on symmetry line exactly onto r=0 (for all times)
         if(ibound==Symmetry_line_id)
          {
           nod_pt->x(t,0)=0.0;
          }
        }
      }
    }
  }
 

 // Update pointer to central node
 oomph_info << "Updating central node" << std::endl;
 Central_node_on_free_surface_pt=0;
 const unsigned n_boundary_node = Fluid_mesh_pt->
  nboundary_node(Free_surface_boundary_id);
 for(unsigned n=0;n<n_boundary_node;n++)
  {
   /// Get pointer to the free surface node
   Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id,n);
   
   // Is this the one?
   if (nod_pt->is_on_boundary(Symmetry_line_id))
    {
     if (Central_node_on_free_surface_pt!=0)
      {
       oomph_info << "Odd -- more than one free surface node on sym line?\n";
       abort();
      }
     Central_node_on_free_surface_pt=nod_pt;
    }
  }
 
 if (Central_node_on_free_surface_pt==0)
  {
   oomph_info << "Odd -- not found the free surface node on sym line...\n";
   abort();
  }
 else
  {
   oomph_info << "Updated central node" << std::endl;
   oomph_info << "Central node now at: " 
              << Central_node_on_free_surface_pt->x(0) << " " 
              << Central_node_on_free_surface_pt->x(1) << " " 
              << std::endl;
  }

} // end of complete_problem_setup


//============start_of_create_free_surface_elements===============
/// Create elements that impose the kinematic and dynamic bcs
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void AxisymmetricVibratingShellProblem<ELEMENT>::create_free_surface_elements()
{ 
 // How many bulk fluid elements are adjacent to boundary b?
 unsigned n_element = 
  Fluid_mesh_pt->nboundary_element(Free_surface_boundary_id);
 
 // Loop over the bulk fluid elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk fluid element that is 
   // adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Fluid_mesh_pt->boundary_element_pt(Free_surface_boundary_id,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = 
    Fluid_mesh_pt->face_index_at_boundary(Free_surface_boundary_id,e);
   
   // Create new element
   ElasticAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt =
    new ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
     bulk_elem_pt,face_index); 
   
   // Add it to the mesh
   Free_surface_mesh_pt->add_element_pt(el_pt);
   
   //Add the appropriate boundary number
   el_pt->set_boundary_number_in_bulk_mesh(Free_surface_boundary_id);
   
   //Specify the capillary number
   el_pt->ca_pt() = &Problem_Parameter::Ca;

   //Specify the Strouhal number
   el_pt->st_pt() = &Problem_Parameter::St;
   
   // Specify the bubble pressure (pointer to Data object and 
   // index of value within that Data object that corresponds
   // to the traded pressure
   el_pt->set_external_pressure_data(External_pressure_data_pt); 
   
  } 
 
} // end of create_free_surface_elements


//==============start of global_temporal_error_norm=======================
/// Calculate the global temporal error norm
//========================================================================
template<class ELEMENT>
double AxisymmetricVibratingShellProblem<ELEMENT>::global_temporal_error_norm()
{

 //oomph_info << "hierher bypassed global temporal error norm\n";
 //return 0.0;

 // Initialise
 double global_error = 0.0;
 
 //Find out how many nodes there are in the problem
 const unsigned n_node = Fluid_mesh_pt->nnode();
 

 oomph_info << "in here with node = " << n_node << std::endl;

 //Loop over the nodes and calculate the errors in the positions
 for(unsigned n=0;n<n_node;n++)
  {
   //Find number of dimensions of the node
   const unsigned n_dim = Fluid_mesh_pt->node_pt(n)->ndim();
   //Set the position error to zero
   double node_position_error = 0.0;
   //Loop over the dimensions
   for(unsigned i=0;i<n_dim;i++)
    {
     //Get position error
     double error = 
      Fluid_mesh_pt->node_pt(n)->position_time_stepper_pt()->
      temporal_error_in_position(Fluid_mesh_pt->node_pt(n),i);
     
     //Add the square of the individual error to the position error
     node_position_error += error*error;
    }
   
   //Divide the position error by the number of dimensions
   node_position_error /= n_dim;
 
  //Now add to the global error
   global_error += node_position_error;
  }
 
 //Now the global error must be divided by the number of nodes
 global_error /= n_node;

 oomph_info << "done global error = " << global_error << std::endl;

 //Return the square root of the errr
 return sqrt(global_error);
 

}


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void AxisymmetricVibratingShellProblem<ELEMENT>::
doc_solution(const std::string& comment)
{ 

 ofstream some_file;
 char filename[1000];
 
 oomph_info << "Docing trace step: " 
            << Problem_Parameter::Doc_info_trace.number() << std::endl;

 // Compute errors and assign to each element for plotting
 double max_err=0.0;
 double min_err=0.0;
 compute_error_estimate(max_err,min_err);
 
 // Write restart file
 if ( !CommandLineArgs::command_line_flag_has_been_set
     ("--suppress_restart_files") &&
     !CommandLineArgs::command_line_flag_has_been_set
     ("--validation") )
  {
   sprintf(filename,"%s/restart%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   ofstream dump_file;
   dump_file.open(filename);
   dump_file.precision(20); 
   dump_it(dump_file);
   dump_file.close();
  }
 
 // Get body force
 Vector<double> body_force(3,0.0);
 Vector<double> x(3,0.0);
 Problem_Parameter::oscillation(this->time_pt()->time(),x,body_force);

 // only output the actual solution when it's the right time
 // hierher if(this->time_pt()->time() >= double(Count_doc)*Dt_doc)
  {
   oomph_info << "Docing soln step: " 
            << Problem_Parameter::Doc_info_soln.number() << std::endl;

   // Number of plot points
   unsigned npts=5;

   
   // Actual solution
   sprintf(filename,"%s/soln%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   some_file.open(filename);
   some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
    ->variable_identifier();
   this->Fluid_mesh_pt->output(some_file,npts);
   some_file.close();
   

   // Actual solution on the free surface
   sprintf(filename,"%s/free_surface_soln%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   some_file.open(filename);
   unsigned nel=Free_surface_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<ElasticAxisymmetricFluidInterfaceElement<ELEMENT>*>(
      Free_surface_mesh_pt->element_pt(e))->output(some_file,npts);
    }
   some_file.close();
   

   // Output free surface
   sprintf(filename,"%s/free_surface%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   some_file.open(filename);
   unsigned num_nod=Fluid_mesh_pt->nboundary_node(Free_surface_boundary_id);
   if (num_nod>0)
    {
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Fluid_mesh_pt->boundary_node_pt(Free_surface_boundary_id, inod)
        ->output(some_file);
      }
    }
   some_file.close();


   // Coarse solution (mesh)
   unsigned npts_coarse=2;
   sprintf(filename,"%s/coarse_soln%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   some_file.open(filename);
   some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
    ->variable_identifier();
   this->Fluid_mesh_pt->output(some_file,npts_coarse);
   some_file.close();

   // Doc body force
   sprintf(filename,"%s/body_force%i.dat",
           Problem_Parameter::Doc_info_soln.directory().c_str(),
           Problem_Parameter::Doc_info_soln.number());
   some_file.open(filename);
   some_file << "-0.5 0.0 0.0 "
             << body_force[0] << " " 
             << body_force[1] << " " 
             << body_force[2] << " "
             << std::endl;
   some_file.close();
   

   // Increment the doc_info number
   Problem_Parameter::Doc_info_soln.number()++;

  }

 // Assemble square of L2 norm 
 double square_of_l2_norm=0.0;
 unsigned nel=Fluid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   square_of_l2_norm+=
    dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))->
    square_of_l2_norm();
  }
 Problem_Parameter::Norm_file << sqrt(square_of_l2_norm) << std::endl;


 // Output boundaries
 sprintf(filename,"%s/boundaries%i.dat",
         Problem_Parameter::Doc_info_soln.directory().c_str(),
         Problem_Parameter::Doc_info_soln.number());
 some_file.open(filename);
 this->Fluid_mesh_pt->output_boundaries(some_file);
 some_file.close();
 
 // Get max/min area
 double max_area=0.0;
 double min_area=0.0;
 Fluid_mesh_pt->max_and_min_element_size(max_area, min_area);


 // Compute current volume and error in extrapolation from fluid elements
 double current_vol=0.0;
 double norm_squared=0.0;
 double visc_norm_squared=0.0;
 double extrapolated_norm_squared=0.0;
 double error_norm_squared=0.0;
 double min_invariant=DBL_MAX;
 double max_invariant=-DBL_MAX;
 double min_viscosity=DBL_MAX;
 double max_viscosity=-DBL_MAX;
 for (unsigned e=0;e<nel;e++)
  {
   // Get element
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
    Fluid_mesh_pt->element_pt(e));

   // Add to physical size (actual volume)
   current_vol+=el_pt->compute_physical_size();

   // Check size (in computational coordinates
   double size=el_pt->size();

   // Get norms of invariant
   double el_norm_squared=0.0;
   double el_extrapolated_norm_squared=0.0;
   double el_error_norm_squared=0.0;
   double test_size=0.0;
   test_size=el_pt->square_of_norm_of_strain_invariant(
    el_norm_squared,
    el_extrapolated_norm_squared,
    el_error_norm_squared);

   if (std::fabs(test_size-size)>1.0e-10)
    {
     oomph_info << "Trouble: " 
                << test_size << " " 
                << size << std::endl;
    }
   
   // Get norms of viscosity
   double el_visc_norm_squared=0.0;
   double el_visc_extrapolated_norm_squared=0.0;
   double el_visc_error_norm_squared=0.0;
   el_pt->square_of_norm_of_viscosity(
    el_visc_norm_squared,
    el_visc_extrapolated_norm_squared,
    el_visc_error_norm_squared);

   // Add it...
   norm_squared+=el_norm_squared;
   visc_norm_squared+=el_visc_norm_squared;
   extrapolated_norm_squared+=el_extrapolated_norm_squared;
   error_norm_squared+=el_error_norm_squared;

   // Get viscosity extrema
   double el_min_invariant=0.0;
   double el_max_invariant=0.0;
   double el_min_viscosity=0.0;
   double el_max_viscosity=0.0;
   el_pt->max_and_min_invariant_and_viscosity(el_min_invariant,
                                              el_max_invariant,
                                              el_min_viscosity,
                                              el_max_viscosity);

   // Update overall extrema
   min_invariant=std::min(min_invariant,el_min_invariant);
   max_invariant=std::max(max_invariant,el_max_invariant);
   min_viscosity=std::min(min_viscosity,el_min_viscosity);
   max_viscosity=std::max(max_viscosity,el_max_viscosity);
  }
   
 oomph_info << "Norm of strain rate invariant: "
            << sqrt(norm_squared) << std::endl;
 oomph_info << "Norm of extrapolated strain rate invariant: "
            << sqrt(extrapolated_norm_squared) << std::endl;
 oomph_info << "Norm of error in extrapolated strain rate invariant: "
            << sqrt(error_norm_squared) << " equivalent to " ;
 if (sqrt(norm_squared)!=0.0)
  {
   oomph_info << sqrt(error_norm_squared)/sqrt(norm_squared)*100.0 << " %";
  }
 oomph_info << std::endl;


 oomph_info << "min_invariant = " << min_invariant << "\n"
            << "max_invariant = " << max_invariant << "\n"
            << "min_viscosity = " << min_viscosity << "\n"
            << "max_viscosity = " << max_viscosity 
            << std::endl;

 // Write trace file
 Problem_Parameter::Trace_file 
  << this->time_pt()->time() << " "  // 1
  << Central_node_on_free_surface_pt->x(1) << " " // 2
  << body_force[1] << " " // 3
  << current_vol << " " // 4
  << Fluid_mesh_pt->nelement() << " " // 5 
  << max_area << " " // 6
  << min_area << " " // 7
  << max_err << " " // spatial error 8
  << min_err << " " // spatial error 9
  << sqrt(norm_squared) << " "  // strain invariant 10
  << sqrt(extrapolated_norm_squared) << " "  // strain invariant 11
  << sqrt(error_norm_squared) << " " // strain invariant  12
  << min_invariant << " " // 13
  << max_invariant << " " // 14
  << min_viscosity << " " // 15
  << max_viscosity << " " // 16
  << max_viscosity << " " // 17
  << this->time_pt()->dt() << " " // 18
  << global_temporal_error_norm() << " " // temporal error measure 19
  << Problem_Parameter::Doc_info_trace.number() << " "  // 20
  << std::endl; 

 // Increment the doc_info number
 Problem_Parameter::Doc_info_trace.number()++;

} // end_of_doc_full_solution


//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void AxisymmetricVibratingShellProblem<ELEMENT>::
compute_error_estimate(double& max_err,double& min_err)
{ 
 // Get error estimator
 ErrorEstimator* err_est_pt=Fluid_mesh_pt->spatial_error_estimator_pt();
 
 // Get/output error estimates
 unsigned nel=Fluid_mesh_pt->nelement();
 Vector<double> elemental_error(nel);
 
 // We need a dynamic cast, get_element_errors from the Fluid_mesh_pt
 // Dynamic cast is used because get_element_errors require a Mesh* ans
 // not a SolidMesh*
 Mesh* fluid_mesh_pt=dynamic_cast<Mesh*>(Fluid_mesh_pt);
 err_est_pt->get_element_errors(fluid_mesh_pt,
                                elemental_error);

 // Set errors for post-processing and find extrema
 max_err=0.0;
 min_err=DBL_MAX;
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<MyTaylorHoodElement*>(Fluid_mesh_pt->element_pt(e))->
    set_error(elemental_error[e]);

   max_err=std::max(max_err,elemental_error[e]);
   min_err=std::min(min_err,elemental_error[e]);
  }
  
}


//====start_of_main===========================================
/// Driver code for vibrating drop problem
//============================================================
int main(int argc, char **argv)
{

 ToleranceForVertexMismatchInPolygons::Tolerable_error=1e-3;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress writing of restart files (they're big!)
 CommandLineArgs::specify_command_line_flag("--suppress_restart_files");
 

 // Use extrapolated strain rate when determining viscosity?
 CommandLineArgs::specify_command_line_flag
  ("--use_current_strainrate_for_viscosity");
 
 // Name of restart file
 CommandLineArgs::specify_command_line_flag
  ("--restart_file",
   &Problem_Parameter::Restart_file);
  
 // Output directory
 CommandLineArgs::specify_command_line_flag(
   "--dir",
   &Problem_Parameter::Directory);

 // Multiplier for oscillatory body force (hierher MH hack)
 CommandLineArgs::specify_command_line_flag(
  "--osc_body_force_multiplier",
  &Problem_Parameter::Oscillating_body_force_multiplier);

 // BDF TYPE (1,2 or 4)
 CommandLineArgs::specify_command_line_flag(
   "--bdf_type",
   &Problem_Parameter::BDF_type);

 // Number of previous values to be used for extrapolation
 // of strain rate (<= Problem_Parameter::BDF_type)
 CommandLineArgs::specify_command_line_flag(
   "--nprev_for_extrapolation_of_strain_rate",
   &Problem_Parameter::Nprev_for_extrapolation_of_strain_rate);

 // Min number of steps period (as timestep limitation)
 unsigned min_ntstep_per_period=40;
 CommandLineArgs::specify_command_line_flag(
   "--min_ntstep_per_period",&min_ntstep_per_period);

 // Initial number of steps period 
 unsigned initial_ntstep_per_period=40;
 CommandLineArgs::specify_command_line_flag(
   "--initial_ntstep_per_period",&initial_ntstep_per_period);

 // Manual adaptation?
 unsigned nadapt_interval=1; 
 CommandLineArgs::specify_command_line_flag(
  "--manual_adapt",&nadapt_interval);

 // Uniform element area
 CommandLineArgs::specify_command_line_flag(
  "--uniform_element_area",
  &Problem_Parameter::Uniform_element_area);

 // Use Newtonian constitutive equation (and specify its 
 // viscosity ratio)
 double newtonian_viscosity_ratio=1.0;
 CommandLineArgs::specify_command_line_flag(
  "--use_newtonian_constitutive_equation",
  &newtonian_viscosity_ratio);

 // Use Casson const eqn
 CommandLineArgs::specify_command_line_flag(
  "--use_casson_constitutive_equation",
  &Problem_Parameter::Critical_strain_rate);

 // Use Herschel Bulkley const eqn
 CommandLineArgs::specify_command_line_flag(
  "--use_herschel_bulkley_constitutive_equation",
  &Problem_Parameter::Critical_strain_rate);

 // Use Sisko const eqn
 CommandLineArgs::specify_command_line_flag(
  "--use_sisko_constitutive_equation",
  &Problem_Parameter::Critical_strain_rate);

 // Temporal accuray
 double epsilon_t=2.0e-7;
 CommandLineArgs::specify_command_line_flag(
  "--epsilon_t",&epsilon_t);

 // Flag indicating whether we use fixed point iteration or not 
 // Number of fixed point iterations before we use Aitken extrapolation
 unsigned fixed_point_steps_before_aitken_extrapolation=0;
 CommandLineArgs::specify_command_line_flag(
  "--use_fixed_point_iteration",
  &fixed_point_steps_before_aitken_extrapolation);

 // Validation command line flag
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Choose constitutive equation
 if (CommandLineArgs::command_line_flag_has_been_set("--use_newtonian_constitutive_equation"))
  {
   oomph_info << "Using Newtonian constitutive equation with viscosity ratio: "
              << newtonian_viscosity_ratio << std::endl;
   Problem_Parameter::Const_eqn_pt = 
    new NewtonianConstitutiveEquation<3>(newtonian_viscosity_ratio);
  }
 else if (CommandLineArgs::command_line_flag_has_been_set
          ("--use_casson_constitutive_equation"))
  {
   oomph_info << "Using Casson constitutive equation"
              << " with cutoff strain rate: "
              << Problem_Parameter::Critical_strain_rate 
              << std::endl;
   Problem_Parameter::Const_eqn_pt = 
    new CassonTanMilRegWithBlendingConstitutiveEquation<3>
    (&Problem_Parameter::C_yield_stress,
     &Problem_Parameter::Critical_strain_rate);

   Problem_Parameter::Re=Problem_Parameter::C_Re;
   Problem_Parameter::St=Problem_Parameter::C_St;
   Problem_Parameter::Re_St=Problem_Parameter::C_Re_St;
   Problem_Parameter::Yield_stress=Problem_Parameter::C_yield_stress;
  }
 else if (CommandLineArgs::command_line_flag_has_been_set
          ("--use_herschel_bulkley_constitutive_equation"))
  {
   oomph_info << "Using Herschel Bulkley constitutive equation"
              << " with cutoff strain rate: "
              << Problem_Parameter::Critical_strain_rate 
              << std::endl;
   Problem_Parameter::Const_eqn_pt = 
    new HerschelBulkleyTanMilRegWithBlendingConstitutiveEquation<3>
    (&Problem_Parameter::HB_yield_stress,
     &Problem_Parameter::HB_flow_index,
     &Problem_Parameter::Critical_strain_rate);

   Problem_Parameter::Re=Problem_Parameter::HB_Re;
   Problem_Parameter::St=Problem_Parameter::HB_St;
   Problem_Parameter::Re_St=Problem_Parameter::HB_Re_St;
   Problem_Parameter::Yield_stress=Problem_Parameter::HB_yield_stress;
  }
 else if( (CommandLineArgs::command_line_flag_has_been_set
          ("--use_sisko_constitutive_equation")) ||
          (CommandLineArgs::command_line_flag_has_been_set
           ("--validation")) )
  {
   oomph_info << "Using Sisko constitutive equation"
              << " with cutoff strain rate: "
              << Problem_Parameter::Critical_strain_rate 
              << std::endl;
   Problem_Parameter::Const_eqn_pt = 
    new SiskoTanMilRegWithBlendingConstitutiveEquation<3>
    (&Problem_Parameter::S_alpha,
     &Problem_Parameter::S_flow_index,
     &Problem_Parameter::Critical_strain_rate);

   Problem_Parameter::Re=Problem_Parameter::S_Re;
   Problem_Parameter::St=Problem_Parameter::S_St;
   Problem_Parameter::Re_St=Problem_Parameter::S_Re_St;
   Problem_Parameter::Yield_stress=0.0;
  }
 else
  {
   oomph_info << "Please select constitutive equation!\n";
   abort();
  }


 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   Problem_Parameter::Initial_settling_time=0.1;
   nadapt_interval=5;
   fixed_point_steps_before_aitken_extrapolation=500;
   Problem_Parameter::Directory="RESLT/";
   Problem_Parameter::Critical_strain_rate=1.0e-12;
  }


 // Create generalised Hookean constitutive equations
 Problem_Parameter::Constitutive_law_pt = 
  new GeneralisedHookean(&Problem_Parameter::Nu);
 
 // Open trace file
 Problem_Parameter::Trace_file.open((Problem_Parameter::Directory+"/trace.dat").c_str());
 
 // Write trace file
 Problem_Parameter::Trace_file  
  << "VARIABLES="
  << "\"time\", "  // 1
  << "\"drop height\"," // 2 
  << "\"body force\","  // 3
  << "\"volume\"," // 4
  << "\"number of elements\","// 5 
  << "\"max element area \"," // 6 
  << "\"min element area \"," // 7 
  << "\"max spatial error\"," // 8
  << "\"min spatial error\"," // 9
  << "\"norm of strain invariant\","  // 10
  << "\"norm of extrapolated strain invariant\"," // 11
  << "\"norm of error in extrapolated strain invariant\"," // 12
  << "\"min_invariant \"," // 13
  << "\"max_invariant \"," // 14
  << "\"min_viscosity \"," // 15
  << "\"max_viscosity \"," // 16
  << "\"norm of viscosity \"," // 17
  << "\"dt\"," // 18
  << "\"global temporal error norm\"," // 19
  << "\"doc number\"" // 20
  << std::endl; 


 // Open norm file
 Problem_Parameter::Norm_file.open((Problem_Parameter::Directory+"/norm.dat").c_str()); // hierher needed? Maybe keep 
 // alive for distant future when this becomes a demo driver code...

 // Output directory
 Problem_Parameter::Doc_info_trace.set_directory(Problem_Parameter::Directory);
 Problem_Parameter::Doc_info_soln.set_directory(Problem_Parameter::Directory);


 // Doc constitutive equation
 ofstream constitutive_file;
 constitutive_file.open((Problem_Parameter::Directory+"/constitutive.dat").c_str());
 double invariant_min=1.0e-15;
 double invariant_max=1.0e2;
 unsigned nstep_per_decade=10;
 double ratio=pow(0.1,double(1.0/nstep_per_decade));
 double invariant=invariant_max;
 while (invariant>invariant_min)
  {
   constitutive_file << invariant << " " 
                     << Problem_Parameter::Const_eqn_pt->viscosity(invariant) 
                     << std::endl;
   invariant*=ratio;
  }
 constitutive_file.close();

 // Create problem in initial configuration
 AxisymmetricVibratingShellProblem<
 GeneralisedNewtonianProjectableAxisymmetricTaylorHoodElement
  <MyTaylorHoodElement> > problem;

 // Timestepping parameters:
 //-------------------------

 oomph_info << "Tolerance for adaptive timestepping: epsilon_t = "
            << epsilon_t << std::endl;

 // Initial tiemstep
 double dt_new=1.0/double(initial_ntstep_per_period); 

 // max. permitted timestep
 double dt_max=1.0/double(min_ntstep_per_period);

 // Do timestepping until tmax
 double t_max=500.25; //2.0

 // Tolerance for fixed point iteration in %
 double fixed_point_tol=1.0; //1.0e-6;

 // Restart?
 if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
   problem.restart();
   oomph_info << "Done restart\n";
   
   // Have read in next timestep
   dt_new=problem.next_dt();
   
   // Doc the read-in initial conditions 
   problem.doc_solution();
  }
 else
  {
   oomph_info << "Not doing restart\n";
   
   // Initialise timestepper
   problem.next_dt()=dt_new;
   problem.initialise_dt(problem.next_dt());
   
   // Perform impulsive start from current state
   problem.assign_initial_values_impulsive();

   oomph_info << "zeroth order extrapol (one prev value)" << std::endl;
   problem.set_nprev_for_extrapolation_of_strain_rate_for_all_elements(1);

   // Take one timestep and doc 
   problem.unsteady_newton_solve(dt_new);
   problem.doc_solution();
   
   // Take a few timesteps with fixed timestep and without adaptation
   // to get some history values under our belt
   unsigned nstep=unsigned(t_max/problem.next_dt()); // hierher //5; 
   oomph_info << "Doing " << nstep 
              << " steps with constant timestep and mesh " << std::endl;
   //for (unsigned i=0;i<nstep;i++)
   unsigned count = 0;
   while (problem.time_pt()->time()<t_max)
    {
     // Adapt?
     if( (CommandLineArgs::command_line_flag_has_been_set("--manual_adapt") ||
          CommandLineArgs::command_line_flag_has_been_set("--validation") ) &&
         ((Problem_Parameter::Doc_info_trace.number())%nadapt_interval==0))
      {
       oomph_info << "hierher: time for another adaptation at doc number ="
                  << Problem_Parameter::Doc_info_trace.number() << std::endl;
       problem.adapt();
      }

     if (count==3)
      {
       oomph_info << "next order extrapol (two prev values)" << std::endl;
       problem.set_nprev_for_extrapolation_of_strain_rate_for_all_elements(2);
      }

     oomph_info << "HIERHER: SKIPPING FOURTH ORDER STUFF!\n";
     // if (i==7)
     //  {
     //   oomph_info << "next order extrapol (four prev values)" << std::endl;
     //   problem.set_nprev_for_extrapolation_of_strain_rate_for_all_elements(4);
     //  }


     // problem.unsteady_newton_solve(dt_new);     

     // Adaptive timestep with specified temporal tolerance epsilon_t
     dt_new=problem.adaptive_unsteady_newton_solve(problem.next_dt(),epsilon_t);

     // Next timestep 
     cout << "Suggested new dt: " << dt_new << std::endl;
     if(dt_new > dt_max)
      {
       dt_new=dt_max;
       cout<<"dt limited (from above) to: "<<dt_new<<endl;
      }
     if(dt_new < Problem_Parameter::Dt_min)
      {
       dt_new=Problem_Parameter::Dt_min;
       cout<<"dt limited (from below) to: "<<Problem_Parameter::Dt_min<<endl;
      }   
     problem.next_dt()=dt_new;

     //problem.doc_solution();

     //if(i<3) continue;

     // Fixed point iteration?
     if ( CommandLineArgs::command_line_flag_has_been_set("--use_fixed_point_iteration")||
          CommandLineArgs::command_line_flag_has_been_set("--validation") )
      {
       oomph_info << "DOING FIXED POINT ITERATION\n";
       
       problem.enable_fixed_point_iteration_for_strain_rate_for_all_elements();

       unsigned i_fp=0;
       double error_fixed_point_iteration = DBL_MAX;

       //for (unsigned i_fp=0;i_fp<nfixed_point;i_fp++)
       while(error_fixed_point_iteration > fixed_point_tol)
        {
         if(i_fp == fixed_point_steps_before_aitken_extrapolation)
          {
           oomph_info << "DOING AITKEN EXTRAPOLATION\n";
           problem.enable_aitken_extrapolation_for_all_elements();
          }
         oomph_info << "FIXED POINT RE-SOLVE " << i_fp << " FOR t = " 
                    << problem.time_pt()->time() << "\n";
         problem.newton_solve();
         //problem.doc_solution();
         error_fixed_point_iteration=
          problem.calculate_error_of_fixed_point_iteration();
         
         // if( error_fixed_point_iteration < fixed_point_tol)
         //  {   
         //   oomph_info << "Error in fixed point iteration is less than the "
         //              << "specified tolerance of "<<fixed_point_tol<<"%\n";
           
         //   oomph_info << "\nSTOPPING FIXED POINT ITERATION\n";
         //   problem.doc_solution();
         //   break;
         //  }
         
         // if(i_fp==nfixed_point-1)
         //  {
         //   oomph_info << "Maximum number of fixed point iterations ("
         //              <<nfixed_point<<") reached.\n";

         //   oomph_info << "\nSTOPPING FIXED POINT ITERATION\n";
         //   problem.doc_solution();
         //   break;
         //  }

         // Update latest guess
         problem.update_latest_fixed_point_iteration_guess_for_strain_rate_for_all_elements();

         i_fp++;
        }

       problem.doc_solution();

       // Switch off
       problem.disable_fixed_point_iteration_for_strain_rate_for_all_elements();
      }
     else
      {
       problem.doc_solution();
      }
     
     count++;     
     
     if( (count==20) &&
         CommandLineArgs::command_line_flag_has_been_set("--validation") )
      {
       break;
      }
     

    }
   oomph_info << "hierher bailing out after const timesteps\n";
   exit(0);
  }
 
 
 // Now do temporally adaptive solves, explicitly adaptiving every specified number
 //--------------------------------------------------------------------------------
 // of steps
 //---------
 if (CommandLineArgs::command_line_flag_has_been_set("--manual_adapt"))
  {
 while (problem.time_pt()->time()<t_max)
  {
     // Adaptive timestep with specified temporal tolerance epsilon_t
     dt_new=problem.adaptive_unsteady_newton_solve(problem.next_dt(),epsilon_t);

   // Next timestep 
   cout << "Suggested new dt: " << dt_new << std::endl;
   if(dt_new > dt_max)
    {
     dt_new=dt_max;
     cout<<"dt limited (from above) to: "<<dt_new<<endl;
    }
   if(dt_new < Problem_Parameter::Dt_min)
    {
     dt_new=Problem_Parameter::Dt_min;
     cout<<"dt limited from below to: "<<Problem_Parameter::Dt_min<<endl;
    }   
   problem.next_dt()=dt_new;
   
   // Doc it...
   problem.doc_solution();

   // Adapt?
   if ((Problem_Parameter::Doc_info_trace.number())%nadapt_interval==0)
    {
     oomph_info << "hierher: time for another adaptation at doc number ="
                << Problem_Parameter::Doc_info_trace.number() << std::endl;
     problem.adapt();
    }
    }
  }
 // Doubly adaptive run
 //--------------------
 else
  {
   unsigned max_adapt=1;
   bool first=false;
   dt_new=problem.doubly_adaptive_unsteady_newton_solve(problem.next_dt(),epsilon_t,
                                                        max_adapt,first);

   // Next timestep 
   cout << "Suggested new dt: " << dt_new << std::endl;
   if(dt_new > dt_max)
    {
     dt_new=dt_max;
     cout<<"dt limited (from above) to: "<<dt_new<<endl;
    }
   if(dt_new < Problem_Parameter::Dt_min)
    {
     dt_new=Problem_Parameter::Dt_min;
     cout<<"dt limited from below to: "<<Problem_Parameter::Dt_min<<endl;
    }   
   problem.next_dt()=dt_new;
   
   // Doc it...
   problem.doc_solution();
  }








 exit(0);

} //End of main

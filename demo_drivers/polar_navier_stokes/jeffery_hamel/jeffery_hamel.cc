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
/// ////////////////////////////////////////////////////////////////////////////
// What it does:                                                             //
/// //////////////////////////////////////////////////////////////////////////// 

// Created 15/07/09
// (by copying jh_streamfunction/jh_bifurcated_solutions.cc)

// Steps up along bifurcated branch monitoring Jacobian sign for 
// subsequent bifurcations


/// ////////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////////

// Standard includes
#include <algorithm>

//Generic includes
#include "generic.h"
#include "polar_navier_stokes.h"


using namespace std;
using namespace oomph;
 
//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
  // Wedge dimensions
  double R_l=0.01;
  double R_r=1.; 

  // Initial number of elements
 int xmesh=30;//60;
 int ymesh=15;//30;
  
  // Initial value for the Reynolds number
  double Re=34.;
  // Initial value for alpha
  double Alpha = 0.1;

  // Stepping parameters for Re
  double Rstep_prestart=30.0;
  double Rmax_prestart=94.;
  double Rstep=0.1;
  double Rmax=100.;

  // Step beyond pitchfork to find bifurcated state
  double epsilon=0.1;

  // Inlet traction
  bool inlet_traction=false;
  double eta_inlet=1.0;
  // Outlet traction
  bool outlet_traction=true; 
  double eta_outlet=0.0;

  // Pin v at inlet switch
  bool pinv=true;

  // Use stokes flow as Dirichlet?
  bool stokes=false;

  //Log mesh (with or without outlet)
  bool log_mesh=true;
  bool new_outlet_region=true;

  /// Unused (but assigned) function to specify tractions
  void traction_function(const double &time, 
                         const Vector<double> &x, 
                         Vector<double> &traction)
  {
   //Impose traction
   traction[0] = 0.0;
   traction[1] = 0.0;
  }

  //Storage for number of adapts so far
  int uniform = 0; int adaptive = 0;

} // end_of_namespace


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
///  Further includes:
#include "meshes/rectangular_quadmesh.h"
#include "./refineable_r_mesh.h"
#include "./streamfunction_include.h"


// Also include the FluxConstraint and jh_mesh classes to save space
#include "./jh_includes.h"

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class PolarNSProblem : public Problem
{
private:

 /// Data object whose single value stores the external pressure
 Data* External_pressure_data_pt;

public:

 /// Constructor
 PolarNSProblem();

 /// Destructor to clean up memory
 ~PolarNSProblem();

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned int &e, const unsigned int &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_solve(){}

 /// Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_solve()
 {
  unsigned long start=1,end=4;
  if(Global_Physical_Variables::inlet_traction) end=2;
  if(Global_Physical_Variables::outlet_traction) start=3;

  double A = Global_Physical_Variables::Alpha;

  for(unsigned long ibound=start;ibound<end;ibound+=2)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      //Get r
      double r = mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
      double phi = mesh_pt()->boundary_node_pt(ibound,inod)->x(1);

      //Use exact stokes solution if necessary
      double exact_u;
      if(Global_Physical_Variables::stokes)
	{
	 double m=2.*A;
	 double mu=(1./(sin(m)-m*cos(m)));
	 exact_u=(mu/r)*(cos(m*phi)-cos(m));
	}
      else{exact_u=0.75*(1.-pow(phi,2.))/(A*r);}
      
      mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,exact_u);
      mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
     }
   }//End of dirichlet inlet/outlet in u

  if(Global_Physical_Variables::pinv)
   {
    //If we haven't just pinned v at the inlet (end=2)
    //then pin it here
    if(end==2)
     {
      //Pin v at inlet 
      unsigned long ibound=3;
      unsigned num_nod= mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
       {     
	// No flow in azimuthal direction:
	mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
       }
     }
   }//End of pinv

  // Loop over remaining boundaries (side walls)
  unsigned num_bound = mesh_pt()->nboundary();

  // No penetration and no slip condition on remaining boundaries
  for(unsigned long ibound=0;ibound<num_bound;ibound=ibound+2)
    {
      unsigned num_nod= mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
	{
	  mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.);
	  //No flow in theta direction
	  mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.);
        }
    }//End of setting side walls to zero

 } // end_of_actions_before_solve

 /// Before adaptation: 
 void actions_before_adapt()
  {
   //Strip off all surface elements
   mesh_pt()->remove_traction_elements();
   mesh_pt()->remove_shear_elements();

   //Reset the external pressure data pointer if there is one
   //using namespace Global_Physical_Variables;
   //if(inlet_traction && outlet_traction) 
   //{delete External_pressure_data_pt;External_pressure_data_pt=0;}

   //Remove global data from problem
   flush_global_data();

  }

 /// After adaptation: Unpin all pressures and then pin redudant pressure dofs.
 void actions_after_adapt()
  {
   //Reassign fluid elements to storage
   mesh_pt()->assign_fluid_element_vector();

   //Reattach surface elements
   if(Global_Physical_Variables::inlet_traction) mesh_pt()->make_traction_elements(false);
   if(Global_Physical_Variables::outlet_traction) mesh_pt()->make_traction_elements(true);
   if(Global_Physical_Variables::inlet_traction && Global_Physical_Variables::outlet_traction) mesh_pt()->make_flux_element();

   //Reattach shear elements
   mesh_pt()->make_shear_elements();

   // Unpin all pressure dofs
   RefineablePolarNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->fluid_elt_vector());
  
    // Pin redundant pressure dofs
   RefineablePolarNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->fluid_elt_vector());

   //Determine whether we have traction elements at both ends
   if(Global_Physical_Variables::inlet_traction && Global_Physical_Variables::outlet_traction) {assign_external_pressure();}

   // Now set the first pressure value in element 0 to 0.0
   // But only if all B.C.s are Dirichlet
   if(!Global_Physical_Variables::inlet_traction && !Global_Physical_Variables::outlet_traction) fix_pressure(0,0,0.0);

  } // end_of_actions_after_adapt

 /// Pin boundaries
 void pin_boundaries()
  {
   // Inlet/Outlet perturbation  
   //Which ends do we need to pin both u and v at?
   unsigned long start=1,end=4;
   if(Global_Physical_Variables::inlet_traction) end=2;
   if(Global_Physical_Variables::outlet_traction) start=3;

   for(unsigned long ibound=start;ibound<end;ibound+=2)
    {
     cout << "Pinning u on boundary: " << ibound << endl;
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //pin u velocity
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       //pin v velocity
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
      }
    } // End pinning of inlet and outlet

   if(Global_Physical_Variables::pinv)
    {
     //If we haven't just pinned v at the inlet (end=2)
     //then pin it here
     if(end==2)
      {
       unsigned long ibound=3;
       cout << "Pinning v at inlet. " << endl;
       unsigned num_nod= mesh_pt()->nboundary_node(ibound);
       for (unsigned inod=0;inod<num_nod;inod++)
	{
	 //pin u velocity
	 mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
	}
      }
    } // end pin v at inlet

   //Pin velocities on side walls to zero
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound+=2)
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Loop over values (u and v velocities)
       for (unsigned i=0;i<2;i++)
	{
	 mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
	}
      }
    } // end loop over side walls

  } // end_of_pin_boundaries

 /// Set up external pressure and pass to traction elements
 void setup_external_pressure()
  {
    //Create a Data object whose single value stores the
   //external pressure
   External_pressure_data_pt = new Data(1);

   // Set external pressure
   External_pressure_data_pt->set_value(0,0.0);

   //Now assign it to the problem
   assign_external_pressure();

  } // end of setup_external_pressure

 /// Set up external pressure and pass to traction elements
 void assign_external_pressure()
  {
   // Regard the external pressure as an unknown and add
   // it to the problem's global data.
   add_global_data(External_pressure_data_pt);

   //Can pin our external pressure if needed
   //External_pressure_data_pt->pin(0);

   //pass flux element the external pressure data
   FluxConstraint* flux_el_pt = this->mesh_pt()->flux_constraint_elt_pt();

   flux_el_pt->set_pressure_data(External_pressure_data_pt);
 
   //Need to pass pointer to external pressure to inlet traction elements
   unsigned n_inlet = mesh_pt()->inlet_traction_elt_length();
   for(unsigned e=0;e<n_inlet;e++)
    {
     this->mesh_pt()->inlet_traction_elt_pt(e)->set_external_pressure_data(External_pressure_data_pt);
    }

  } // end of assign_external_pressure

 // Access function for the specific mesh
 jh_mesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<jh_mesh<ELEMENT>*>(Problem::mesh_pt());
  }

 // Experiment! - Still necessary as Alpha included in boundary conditions (if not mesh itself anymore)
 void actions_before_newton_convergence_check()
 {
  actions_before_solve();
 }// End of actions_before_newton_convergence_check

 //Output the solution at nplot points in each direction across an individual element
 void element_output(std::ostream& outfile,unsigned e,const unsigned& nplot);

 // Return the value of the external pressure
 double get_pext();

 // Set the value of the external pressure
 void set_pext(double p);

 // A method for computing the shear stress along both walls
 void get_shear_stress(Vector<double> &shear_stress);

 // Pin nodes on boundary to zero for outputting eigenfunction
 void pin_boundaries_to_zero();

 // Create and solve a streamfunction problem
 void output_streamfunction(DocInfo &doc_info,bool eigen);

 // Function to return a vector specifying the symmetry being broken
 void get_symmetry(DoubleVector &symmetry,double pressure_norm);

 // Function to return sign of the Jacobian
 int get_Jacobian_sign();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 // Header for doc_solution
 void header( ofstream &some_file );

}; // end_of_problem_class

//========================================================================
/// Output an individual element 
//========================================================================
template<class ELEMENT>
void PolarNSProblem<ELEMENT>::element_output(std::ostream& outfile,unsigned e,const unsigned& nplot)
{ 
  // Upcast from GeneralisedElement to the present element
  ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

  el_pt->output(outfile,nplot);

} // end_of_element_output

//========================================================================
// Function to return the value of pext
//========================================================================

template<class ELEMENT> 
double PolarNSProblem<ELEMENT>::get_pext()
 // Return the value of the external pressure
{
 return External_pressure_data_pt->value(0);
}

//========================================================================
// Function to set the value of pext
//========================================================================

template<class ELEMENT> 
void PolarNSProblem<ELEMENT>::set_pext(double p)
 // Return the value of the external pressure
{
 External_pressure_data_pt->set_value(0,p);
}

//========================================================================
// My function for calculating shear stress along walls
//========================================================================
template<class ELEMENT> 
void PolarNSProblem<ELEMENT>::get_shear_stress(Vector<double> &shear_stress)
{
 //Find the number of shear elements
 unsigned Nshear_lower = mesh_pt()->lower_stress_integral_elt_length();
 unsigned Nshear_upper = mesh_pt()->upper_stress_integral_elt_length();

 //Storage for element pointer
 PolarStressIntegralElement<ELEMENT>* el_pt;
 
 //Storage for answer
 double lower=0.0,upper=0.0;

 for(unsigned e=0;e<Nshear_lower;e++)
  {
   el_pt = dynamic_cast<PolarStressIntegralElement<ELEMENT>*>(this->mesh_pt()->lower_stress_integral_elt_pt(e));
   lower+= el_pt->get_shear_stress();
  }
 for(unsigned e=0;e<Nshear_upper;e++)
  {
   el_pt = dynamic_cast<PolarStressIntegralElement<ELEMENT>*>(this->mesh_pt()->upper_stress_integral_elt_pt(e));
   upper+= el_pt->get_shear_stress();
  }

 shear_stress[0]=lower;
 shear_stress[1]=upper;

} //End of get shear_stress

//========================================================================
// My function for calculating shear stress along walls
//========================================================================
template<class ELEMENT> 
void PolarNSProblem<ELEMENT>::pin_boundaries_to_zero()
{
  //Initially set to zero both inlet and outlet
  unsigned long start=1,end=4;
  //But which ends have we pinned u at?
  if(Global_Physical_Variables::inlet_traction) end=2;
  if(Global_Physical_Variables::outlet_traction) start=3;
  for(unsigned long ibound=start;ibound<end;ibound+=2)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      //Pin to zero
      mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
     }
   }//End of dirichlet inlet/outlet in u

} //End of pin_boundaries_to_zero

//========================================================================
// My function for outputting the streamfunction
//========================================================================
template<class ELEMENT> 
void PolarNSProblem<ELEMENT>::output_streamfunction(DocInfo &doc_info,bool eigen)
{
 //Create a streamfunction problem 
 // This needs to know if its solving for the flow or the eigenfunction
 StreamfunctionProblem stream_problem(eigen); 

 // Document the refines to date
 Vector<Vector<unsigned> > refinement_pattern;
 cout << "Obtaining refinement pattern of mesh" << endl;
 this->mesh_pt()->get_refinement_pattern(refinement_pattern);
 cout << "Refining streamfunction mesh to same level" << endl;
 stream_problem.actions_before_adapt();
 stream_problem.mesh_pt()->refine_base_mesh(refinement_pattern);
 stream_problem.actions_after_adapt();
 cout << "Streamfunction mesh refined to same level" << endl;
 cout << "Number of (streamfunction) equations:" << stream_problem.assign_eqn_numbers() << endl;
 cout << "Assigning velocities to streamfunction problem" << endl;
 clock_t start=clock();

 unsigned num_nod1 = this->mesh_pt()->nnode();
 unsigned num_nod2 = stream_problem.mesh_pt()->nnode();
 cout << "num_nod1: " << num_nod1 << endl;
 cout << "num_nod2: " << num_nod2 << endl;
 for(unsigned i=0;i<num_nod2;i++)
   {
    double r2 = stream_problem.mesh_pt()->node_pt(i)->x(0);
    double phi2 = stream_problem.mesh_pt()->node_pt(i)->x(1);

    int found=-1;int j=0;
    double u=0.;double v=0.;
    while(found<0)
      {
	double r = this->mesh_pt()->node_pt(j)->x(0);
	double phi = this->mesh_pt()->node_pt(j)->x(1);
	if((abs(r-r2)<1.e-8) && (abs(phi-phi2)<1.e-8))
	{
	 u = this->mesh_pt()->node_pt(j)->value(0);
	 v = this->mesh_pt()->node_pt(j)->value(1);
	 found=1;
	}
	j+=1;
      }

    // Transfer values across from old node to new node
    stream_problem.mesh_pt()->node_pt(i)->set_value(1,u);
    stream_problem.mesh_pt()->node_pt(i)->set_value(2,v);
   }

 clock_t end=clock();
 cout << "This took  "<<(double(end-start))/CLOCKS_PER_SEC<<" seconds."<<endl; 

 stream_problem.newton_solve();
 stream_problem.doc_solution(doc_info);
 //Select output file 
 char file_name[100];
 sprintf(file_name,"my_streamfunction%i.dat",doc_info.number());
 stream_problem.my_output(201,81,false,file_name);

} //End of output_streamfunction


//========================================================================
// Function to output vector describing the symmetry being broken
//========================================================================
template<class ELEMENT> 
void PolarNSProblem<ELEMENT>::get_symmetry(DoubleVector &symmetry,
                                           double pressure_norm)
{
 // First store the current dof values
 this->store_current_dof_values();

 //Loop over nodes and set nodal values according to symmetry
 unsigned num_nod= mesh_pt()->nnode();
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get coordinates of node
   double r = mesh_pt()->node_pt(inod)->x(0);
   double phi = mesh_pt()->node_pt(inod)->x(1);
   
   // Anti-symmetric in u, symmetric in v
   mesh_pt()->node_pt(inod)->set_value(0,r*phi*(1.-phi*phi));
   mesh_pt()->node_pt(inod)->set_value(1,(1.-phi*phi));
  }

 //Loop over fluid elements to set pressure dofs
 unsigned n_fluid = mesh_pt()->fluid_elt_length();
 for(unsigned e=0;e<n_fluid;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   Vector<double> s(2,0.0);
   double r = el_pt->interpolated_x(s,0);
   double phi = el_pt->interpolated_x(s,1);
   s[0]=-1.;double r_min = el_pt->interpolated_x(s,0);
   s[0]=1.;double r_max = el_pt->interpolated_x(s,0);
   s[1]=-1.;double phi_min = el_pt->interpolated_x(s,1);
   s[1]=1.;double phi_max = el_pt->interpolated_x(s,1);

   double a = phi*(1.-phi*phi)*r*r;
   double b = phi*(1.-phi*phi)*2.*r*(r_max-r_min)/2.;
   double c = (1.-3.*phi*phi)*r*r*(phi_max-phi_min)/2.;

   //Set the pressure dofs to give an antisymmetric pressure
   *el_pt->internal_data_pt(0)->value_pt(0) = a*pressure_norm;
   *el_pt->internal_data_pt(0)->value_pt(1) = b*pressure_norm;
   *el_pt->internal_data_pt(0)->value_pt(2) = c*pressure_norm;
   
  } // end loop over fluid elements

 //Finally set the value of pext (if appropriate)
 if(Global_Physical_Variables::inlet_traction 
    && Global_Physical_Variables::outlet_traction) this->set_pext(0.0);

 //Get the dofs
 this->get_dofs(symmetry);

 // Finally restore the dof values
 this->restore_dof_values(); 

} //End of get_symmetry

//========================================================================
// My new function to return the sign of the Jacobian
//========================================================================

template<class ELEMENT> 
int PolarNSProblem<ELEMENT>::get_Jacobian_sign()
{
  int sign=Sign_of_jacobian;
  return sign;
} //End of get_Jacobian_sign

//==start_of_constructor==================================================
/// Constructor for PolarNS problem 
//========================================================================
template<class ELEMENT>
PolarNSProblem<ELEMENT>::PolarNSProblem()
{ 
 // Increace max residuals slightly
 Max_residuals=20000;
 Max_newton_iterations=20;

 // Setup mesh

 // # of elements in x-direction
 unsigned Nx=Global_Physical_Variables::xmesh;

 // # of elements in y-direction
 unsigned Ny=Global_Physical_Variables::ymesh;

 // Build and assign mesh
 Problem::mesh_pt() = new jh_mesh<ELEMENT>(Nx,Ny);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 //dynamic_cast<Refineable_r_mesh<ELEMENT>*>(mesh_pt())->
 //spatial_error_estimator_pt()=error_estimator_pt;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 //pin the necessary boundaries
 pin_boundaries();

 // Complete the build of all elements so they are fully functional

 //Find number of fluid elements in mesh
 unsigned n_fluid = mesh_pt()->fluid_elt_length();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_fluid;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   //Set alpha
   el_pt->alpha_pt() = &Global_Physical_Variables::Alpha;
   //Set the re_st pointer
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;
   
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineablePolarNavierStokesEquations::
   pin_redundant_nodal_pressures(mesh_pt()->fluid_elt_vector());

 //Determine whether we have traction elements at both ends
 if(Global_Physical_Variables::inlet_traction && Global_Physical_Variables::outlet_traction) {setup_external_pressure();}
 
 // Now set the first pressure value in element 0 to 0.0
 // But only if all B.C.s are Dirichlet
 if(!Global_Physical_Variables::inlet_traction && !Global_Physical_Variables::outlet_traction) fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << endl; 

} // end_of_constructor


//==start_of_destructor===================================================
/// Destructor for PolarNS problem 
//========================================================================
template<class ELEMENT>
PolarNSProblem<ELEMENT>::~PolarNSProblem()
{ 

 // Mesh gets killed in general problem destructor

} // end_of_destructor

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PolarNSProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points 
 unsigned npts;
 npts=3; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 header( some_file );
 unsigned n_fluid = mesh_pt()->fluid_elt_length();
 for(unsigned e=0;e<n_fluid;e++)
  {
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   el_pt->output(some_file,npts);
  }
 some_file.close();

} // end_of_doc_solution

template<class ELEMENT>
void PolarNSProblem<ELEMENT>::header( ofstream &some_file )
 {
  using namespace Global_Physical_Variables;
  some_file << "# Refineable mesh" << "\n";
  some_file << "# Re = " << Re << " Alpha = " << Alpha << " R_l = " << R_l << "\n";
  some_file << "# Initial xmesh = " << xmesh << " Initial ymesh = " << ymesh << "\n";
  some_file << "# Uniform refinements: " << uniform << " Adaptive refinements: " << adaptive << "\n";
  some_file << "# inlet_traction = " << inlet_traction << " eta_inlet = " << eta_inlet << "\n";
  some_file << "# outlet_traction = " << outlet_traction << " eta_outlet = " << eta_outlet << "\n";
  some_file << "# log_mesh = " << log_mesh << " new_outlet_region = " << new_outlet_region << "\n";
  some_file << "# pinv = " << pinv << "\n";
  some_file << "# Total equations = " << this->ndof() << "\n";
  some_file << "\n";
 }

//==start_of_main======================================================
/// Driver for PolarNS test problem -- test drive
/// with two different types of element.
//=====================================================================
int main()
{ 
 // Set up doc info
 // ---------------
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 doc_info.number()=0; 
 // And again for eigenmode
 DocInfo eigenmode;
 eigenmode.set_directory("RESLT");
 eigenmode.number() = 100;
 // ---------------

 using namespace Global_Physical_Variables;
 Vector<double> shear_stress(2,0.0);

 // Build the problem
 PolarNSProblem<RefineablePolarCrouzeixRaviartElement > problem;
 //PolarNSProblem<RefineablePolarTaylorHoodElement > problem;
  
 cout << endl << "# xmesh: " << xmesh << " ymesh: " << ymesh << endl;

 //Increace the Reynolds number in large increments 
 for(;Re<Rmax_prestart;Re+=Rstep_prestart)
   {
     cout << endl << "Prestart solving for Re: " << Re << endl << endl; 

     problem.newton_solve(); 
     problem.doc_solution(doc_info);
     doc_info.number()+=1;
     //problem.output_streamfunction(doc_info,false);
   }

 exit(1);

 //Always start from same Reynolds number
 Re=Rmax_prestart;
 
 problem.newton_solve();
 problem.doc_solution(doc_info);
 doc_info.number()+=1;

 bool sym=true;

 //Storage for eigenvalues and eigenvectors 
 Vector<complex<double> > eigenvalues;
 Vector<DoubleVector> eigenvectors;
 //Desired number eigenvalues
 unsigned n_eval=6;
 //Set shift
 static_cast<ARPACK*>( problem.eigen_solver_pt() ) -> set_shift( 0.4 );

 //Solve the eigenproblem
 cout << "facking" << endl;
 problem.solve_eigenproblem(n_eval,eigenvalues,eigenvectors);
 for(unsigned k=0;k<n_eval;k++)  cout << "eigenvalue: " << eigenvalues[k] << endl;
 cout << "arpack" << endl;
 //Activate pitchfork tracking
 DoubleVector eigenvector(eigenvectors[0]);

 //Compute symmetry vector
 DoubleVector symmetry;
 double norm = 1.e-2;
 problem.get_symmetry(symmetry,norm);

 ofstream out("scrap_up.dat");
 out << "# Output { Re, Alpha, Int du_dphi(-1), Int du_dphi(1), Delta_Shear_Stress, Jacobian sign }" << endl;
 out << "# Radius ratio: " << R_l << " xmesh: " << xmesh << " ymesh " << ymesh << endl;
 out << "# Alpha = " << Alpha << endl; 
 out << "# New outlet region: " << new_outlet_region << endl; 
 out << "# Inlet traction: " << inlet_traction << " Outlet traction: " << outlet_traction << endl;
 out << "# eta_inlet (Homotopy parameter): " << eta_inlet << " eta_outlet: " << eta_outlet << endl;
 out << "# Pin v at inlet: " << pinv << endl;
 out << "# Using symmetry for pitchfork solve: " << sym << endl << endl;

 //Activate pitchfork tracking
 if(sym)
  {
   cout << endl << "Using symmetry vector from get_symmetry routine.  Re = " << Re << endl << endl;
   problem.activate_pitchfork_tracking(&Global_Physical_Variables::Re,symmetry);
  }
 else
  {
   cout << endl << "Using eigenvector for pitchfork detection.  Re = " << Re << endl << endl;
   problem.activate_pitchfork_tracking(&Global_Physical_Variables::Re,eigenvector);
  }

 //Find an intitial pitchfork
 problem.newton_solve();

 //Check on slack parameter
 unsigned long n_dof = problem.ndof();cout << "ndof: " << n_dof << endl;
 DoubleVector soln;
 problem.get_dofs(soln);
 unsigned long half=(n_dof/2)-1;
 cout << "According to dofs Re is: " << problem.dof(half) 
      << " and slack parameter is: " << problem.dof(n_dof-1) << endl;

 double slack=problem.dof(n_dof-1);

 problem.deactivate_bifurcation_tracking();

 doc_info.number()=10;
 problem.doc_solution(doc_info);
 problem.get_shear_stress(shear_stress);
 cout << "Shear stress is: Lower = " << shear_stress[0] << " Upper = " << shear_stress[1] 
      << " Delta = " << (shear_stress[0]+shear_stress[1]) << endl;
 
 cout << endl << "Pitchfork detected at Re = " << Re << " and Alpha = " << Alpha << endl;
 cout << "========================================================" << endl;

 out << Re << " " << Alpha << " " << shear_stress[0] << " " << shear_stress[1] << " " 
     << (shear_stress[0]+shear_stress[1]) << " " << 0 << " " << slack << endl;


 cout << "========================================================" << endl;
 cout << endl << "Computing eigenvalue at pitchfork.  Re = " << Re << endl << endl;
 //Solve the eigenproblem
 problem.solve_eigenproblem(n_eval,eigenvalues,eigenvectors);
 for(unsigned k=0;k<n_eval;k++)  cout << "eigenvalue: " << eigenvalues[k] << endl;

 // OUTPUT THE EIGENMODE
 eigenmode.number()=100;
 problem.store_current_dof_values();
 problem.pin_boundaries_to_zero();
 problem.assign_eigenvector_to_dofs( eigenvectors[0] );
 problem.doc_solution( eigenmode ); 
 problem.restore_dof_values(); 

 //Store solution at bifurcation
 n_dof = problem.ndof(); cout << "ndof: " << n_dof << endl;
 DoubleVector soln_at_bifurcation;
 problem.get_dofs(soln_at_bifurcation);
 DoubleVector critical_eigenvector(eigenvectors[0]);

 // part one ---------------------------------------------------------------------------------------

 double multiplier=-5.0;
 int sign=0;

 cout << endl << "Adding eigenvector[0] with multiplier " << multiplier << " to degrees of freedom" << endl; 
 for(unsigned long d=0;d<n_dof;d++)
  {
   problem.dof(d)=soln_at_bifurcation[d]+multiplier*critical_eigenvector[d];
  }

 //Now we increace/decreace the Reynolds number slightly and compute a solution
 Re+=epsilon;
 cout << "Adjusting Reynolds number by " << epsilon << " to Re: " << Re
      << " and newton solving" << endl << endl;
 problem.newton_solve();
 problem.get_shear_stress(shear_stress);
 sign=problem.get_Jacobian_sign();
 cout << "Shear stress is: Lower = " << shear_stress[0] << " Upper = " << shear_stress[1] 
      << " Delta = " << (shear_stress[0]+shear_stress[1]) << endl;

 out << Re << " " << Alpha << " " << shear_stress[0] << " " << shear_stress[1] << " " << (shear_stress[0]+shear_stress[1]) << " " << sign << endl;

 //Now we step in Re and try to continue the bifurcated state
 double ds=Rstep;
 Re-=ds*0.9;

 ds=problem.arc_length_step_solve( &Re, ds );
 cout << endl << "Re: " << Re << endl;
 problem.get_shear_stress(shear_stress);
 sign=problem.get_Jacobian_sign();
 cout << "Shear stress is: Lower = " << shear_stress[0] << " Upper = " << shear_stress[1] 
      << " Delta = " << (shear_stress[0]+shear_stress[1]) << endl;
 if(inlet_traction && outlet_traction) cout << " pext: " << problem.get_pext() << endl;
 out << Re << " " << Alpha << " " << shear_stress[0] << " " << shear_stress[1] << " " << (shear_stress[0]+shear_stress[1]) << " " << sign << endl;

 for(;Re<Rmax;)
  {
   problem.arc_length_step_solve( &Re, ds );
   cout << endl << "Re: " << Re << endl;
   problem.get_shear_stress(shear_stress);
 sign=problem.get_Jacobian_sign();
   cout << "Shear stress is: Lower = " << shear_stress[0] << " Upper = " << shear_stress[1] 
	<< " Delta = " << (shear_stress[0]+shear_stress[1]) << endl;
   if(inlet_traction && outlet_traction) cout << " pext: " << problem.get_pext() << endl;
   out << Re << " " << Alpha << " " << shear_stress[0] << " " << shear_stress[1] << " " << (shear_stress[0]+shear_stress[1]) << " " << sign << endl;
  }

 //Compute the picture at Re=Rmax/Re=Rmin
 Re=Rmax;
 cout << endl << "Computing picture at Re: " << Re << endl << endl;
 problem.newton_solve(2);
 doc_info.number()=20;
 problem.doc_solution(doc_info);
 problem.output_streamfunction(doc_info,false);
 doc_info.number()+=1;
 problem.get_shear_stress(shear_stress);
 sign=problem.get_Jacobian_sign();
 out << Re << " " << Alpha << " " << shear_stress[0] << " " << shear_stress[1] << " " << (shear_stress[0]+shear_stress[1]) << " " << sign << endl;

 out.close();

} // end_of_main


/*
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver for an adaptive 2D advection diffusion problem with flux boundary 
// conditions using two separate meshes for the bulk and surface meshes*/

//Generic routines
#include "generic.h"

// The Advection Diffusion equations
//#include "advection_diffusion.h"
#include "advection_diffusion/refineable_gen_advection_diffusion_elements_time_wind.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"



using namespace std;

using namespace oomph;
/* template<unsigned DIM>
class GeneralisedAdvectionDiffusionEquations_time_wind : public GeneralisedAdvectionDiffusionEquations<DIM>
{
	public:
	
	 GeneralisedAdvectionDiffusionEquations_time_wind()
      : Source_fct_pt(0),
        Wind_fct_pt(0),
        Conserved_wind_fct_pt(0),
        Diff_fct_pt(0),
        ALE_is_disabled(false)
    {
      // Set Peclet number to default
      Pe_pt = &Default_peclet_number;
      // Set Peclet Strouhal number to default
      PeSt_pt = &Default_peclet_number;
    }
	
	
	typedef void (*GeneralisedAdvectionDiffusionWindFctPt)(
      const Vector<double>& x, double& tval, Vector<double>& wind);
	
	GeneralisedAdvectionDiffusionWindFctPt& wind_fct_pt()
    {
      return Wind_fct_pt;
    }

	 void get_wind_cons_adv_diff(const unsigned& ipt,
                                               const Vector<double>& s,
                                               const Vector<double>& x,
                                               Vector<double>& wind) const;
											   
   protected:
	GeneralisedAdvectionDiffusionWindFctPt Wind_fct_pt;
};


template<unsigned DIM>
void GeneralisedAdvectionDiffusionEquations_time_wind<DIM>::get_wind_cons_adv_diff(const unsigned& ipt,
                                               const Vector<double>& s,
                                               const Vector<double>& x,
                                               Vector<double>& wind) const
											   {
												   double tval=this->time_pt()->time();
												   (*Wind_fct_pt)(x,tval,wind);
											   }
 */


//======start_of_namespace============================================
/// Namespace for AdvectionDiffusion equation 
//====================================================================
namespace ParameterSpace
{
 /// Peclet number
 double Peclet=50.0;
 
 /// Peclet-Strouhal
 double Peclet_St = 1.0;
 
  ///'swimming Peclet' number
 double Beta=10.0;
 
 /// Channel length
 double Length=1200.0;
 
 /// Flow direction sign, changes both the flow and shear profiles
 double flow_sign=1.0;
 
 /// Wind function, the fluid flow
/*  void wind_function(const Vector<double>& x, Vector<double>& wind)
 {
  wind[0]=0.0;
  wind[1]=flow_sign*2.0*(1-(x[0]-1)*(x[0]-1));
 } */
 
 
void wind_function(const Vector<double>& x, double& tval, Vector<double>& wind)
 {
  wind[0]=1/tval;
  wind[1]=flow_sign*2.0*(1-(x[0]-1)*(x[0]-1));
 }
 
 /// Conservative wind function, i.e. the swimming function 
 /// Taken from Bearon et al.
void swimming(const Vector<double> &x, double &tval, Vector<double> &swim)
{
  //local shear	
   double sigma=-flow_sign*(1.0-x[0])*2.0*Peclet/(Beta*Beta);	

  //Radial direction
  double a0=2.05e-001;
  double a2= 1.86e-002;
  double b2= 1.74e-001;
  double b4= 1.27e-002;
 	
  swim[0] = -Beta*sigma*(a0+a2*sigma*sigma)/(1+b2*sigma*sigma+b4*sigma*sigma*sigma*sigma);
	
  //Vertical direction
  a0=5.70e-001;
  a2=3.66e-002;
  b2=1.75e-001;
  b4=1.25e-002;

  swim[1] = -Beta*(a0+a2*sigma*sigma)/(1+b2*sigma*sigma+b4*sigma*sigma*sigma*sigma);
 
  //swim[0]=0.0;
  //swim[1]=0.1;
 
  //Azimuthal direction (ignored)
  //swim[2] = 0.0;
 }
 
 /// Source function: No effect on the system, need to keep it because
 /// I can't skip inputs in the constructor function.
  void source_function(const Vector<double>& x_vect, double& source)
 {
  source=0.0;
 }
 
/// Diffusion function. Taken from Bearon et al.
 void diff_function(const Vector<double> &x, DenseMatrix<double> &D)
{

 //local shear	  
  double sigma=-flow_sign*(1.0-x[0])*2.0*Peclet/(Beta*Beta);	

  //Radial component
  double a0= 9.30e-002;
  double a2= 1.11e-004;
  double b2= 1.19e-001;
  double b4= 1.63e-004; 

  D(0,0) = (a0+a2*sigma*sigma)/(1+b2*sigma*sigma+b4*sigma*sigma*sigma*sigma);
 
  //Radial-Axial component
  a0= 9.17e-002;
  a2= 1.56e-004;
  b2= 2.81e-001;
  b4= 2.62e-002; 

  D(0,1) =  -sigma*(a0+a2*sigma*sigma)/(1+b2*sigma*sigma+b4*sigma*sigma*sigma*sigma);

  //Axial-Radial component
  D(1,0) =  D(0,1);

  //Axial component
  a0=5.00e-002 ;
  a2= 1.11e-001;
  double a4= 3.71e-005; 
  b2= 1.01e-001;
  b4= 1.86e-002 ; 

  D(1,1) =  (a0+a2*sigma*sigma+a4*sigma*sigma*sigma*sigma)/(1+b2*sigma*sigma+b4*sigma*sigma*sigma*sigma);
 
  //D(0,0) = 1.0;
  //D(0,1) = 1.0;
  //D(1,0) = 1.0;
  //D(1,1) = 1.0;
 
 
 /* 2nd row components (commented out)
  //These will not affect the equation
  //D(0,2) = 0.0;
  //D(1,2) = 0.0;
  //D(2,0) = 0.0;
  //D(2,1) = 0.0;
  //D(2,2) = 0.0;*/
 } 

} //end of namespace

//========= start_of_problem_class=====================================
/// 2D AdvectionDiffusion problem on rectangular domain, discretised with
/// 2D QAdvectionDiffusion elements. The specific type of 
/// element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class BiasedActiveMatterDispersionProblem : public Problem
{

public:

 /// Constructor: Pass pointer to wind, conservative wind and diffusion
 BiasedActiveMatterDispersionProblem(RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionSourceFctPt source_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt conserved_wind_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt);

 /// Destructor (empty)
 ~BiasedActiveMatterDispersionProblem(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
 void set_initial_condition();
 
  double global_temporal_error_norm();

private:

 /// No action before solving (empty)
 void actions_before_newton_solve();

 /// Update the problem specs after solve
 void actions_after_newton_solve(){};

 /// Actions before adapt (empty)
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh
 void actions_after_adapt();
 
 /// Actions before implicit timestep (redefining the wind function)
 void actions_before_implicit_timestep(){};
 
  /// Pointer to source function (set to zero, no effect)
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionSourceFctPt Source_fct_pt;

 /// Pointer to the "bulk" mesh
 RefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to wind function
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt Wind_fct_pt;
 
 /// Pointer to swimming function
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt Conserved_wind_fct_pt;
 
 /// Diffusivity
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionDiffFctPt Diff_fct_pt;
}; // end of problem class

//=======start_of_constructor=============================================
/// Constructor for AdvectionDiffusion problem: Pass pointer to wind and
/// diffusion functions.
//========================================================================
template<class ELEMENT>
BiasedActiveMatterDispersionProblem<ELEMENT>::
BiasedActiveMatterDispersionProblem(
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionSourceFctPt source_fct_pt,
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt, 
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionWindFctPt conserved_wind_fct_pt,
 RefineableGeneralisedAdvectionDiffusionEquations_time_wind<2>::GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt)
 : Source_fct_pt(source_fct_pt), Wind_fct_pt(wind_fct_pt), Conserved_wind_fct_pt(conserved_wind_fct_pt), Diff_fct_pt(diff_fct_pt)
{ 

// Put this above as the final input to include diffusion
// GeneralisedAdvectionDiffusionEquations<2>::GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt)

 // Add the time stepper
 this->add_time_stepper_pt(new BDF<2>(true));
 
 // Setup "bulk" mesh

 // # of elements in x-direction
 unsigned n_x=15;

 // # of elements in y-direction
 unsigned n_y=ParameterSpace::Length/2;
 
 // Domain length in x-direction
 // Setting this l_x=2.0 so that x can vary between -1 and 1 when
 // I subtract 1 in vorticity calculations
 double l_x=2.0;

 // Domain length in y-direction
 double l_y=ParameterSpace::Length;

 // Build "bulk" mesh
 Bulk_mesh_pt=new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,this->time_stepper_pt());

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
    Bulk_mesh_pt->max_permitted_error() = 1.0e-3;
  Bulk_mesh_pt->min_permitted_error() = 1.0e-5;
  Bulk_mesh_pt->max_refinement_level() = 10;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 
 // Hakan - Periodic boundary addition
   unsigned nx=n_x;
 unsigned ny=n_y;
 
 unsigned n_node = Bulk_mesh_pt->nboundary_node(0);
 for(unsigned n=0;n<n_node;n++)
  {
   Bulk_mesh_pt->boundary_node_pt(0,n)
    ->make_periodic(Bulk_mesh_pt->boundary_node_pt(2,n));
  }
  
   // Get pointers to tree roots associated with elements on the 
 // left and right boundaries
  Vector<TreeRoot*> left_root_pt(ny);
  Vector<TreeRoot*> right_root_pt(ny);
  for(unsigned i=0;i<ny;i++) 
   {
    left_root_pt[i] = 
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i*nx))->
     tree_pt()->root_pt();
    right_root_pt[i] = 
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(nx-1+i*nx))->
     tree_pt()->root_pt();
   }
   
     // Switch on QuadTreeNames for enumeration of directions
   using namespace QuadTreeNames;
 
  //Set the neighbour and periodicity
  for(unsigned i=0;i<ny;i++) 
   {
    // The western neighbours of the elements on the left
    // boundary are those on the right
    left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
    left_root_pt[i]->set_neighbour_periodic(W); 
    
    // The eastern neighbours of the elements on the right
    // boundary are those on the left
    right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
    right_root_pt[i]->set_neighbour_periodic(E);     
   } // done
     
 // Hakan - Periodic boundary addition ends
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);

 // Combine all submeshes into a single Mesh
  build_global_mesh();

 // Loop over the AdvectionDiffusion bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wind and diffusion function
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to AdvectionDiffusion bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;

   //Set the wind function pointer
   el_pt->wind_fct_pt() = &ParameterSpace::wind_function;
   
   // Swimming component
   el_pt->conserved_wind_fct_pt() =  &ParameterSpace::swimming;
   
   // Diffusivity
   el_pt->diff_fct_pt() = &ParameterSpace::diff_function;

   // Set the Peclet number
   el_pt->pe_pt() = &ParameterSpace::Peclet;
   
  //Set the Peclet Strouhal number
   el_pt->pe_st_pt() = &ParameterSpace::Peclet_St;

   // Set the Peclet number
   el_pt->pe_pt() = &ParameterSpace::Peclet;
   
   std::cout << el_pt->time_fct_pt();
   
   // Set the time pointer
   //el_pt->time_fct_pt() = &this->time_pt()->time();
  }
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/* template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_before_implicit_timestep()
{
	 //Find out how many nodes there are in the problem
 unsigned n_node = Bulk_mesh_pt->nnode();
 
 //Loop over the nodes and calculate the estimated error in the values
 for(unsigned i=0;i<n_node;i++)
  {
	  Node* current_Node=Bulk_mesh_pt->node_pt(i);
	  double u = current_Node->value(0);
	  current_Node->set_value(0,std::max(0.0,u));
  }
 	/*unsigned n_element = Bulk_mesh_pt->nelement();
	
	
 for(unsigned e=0;e<n_element;e++)
  {
	  ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
	  inline Vector<double> windy(int t) {return time_pt()->time()+2.0;}
	  el_pt->wind_fct_pt() = &windy(0);
  }
//}


/* template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_after_newton_solve()
{
 //Find out how many nodes there are in the problem
 unsigned n_node = Bulk_mesh_pt->nnode();
 
 //Loop over the nodes and calculate the estimated error in the values
 for(unsigned i=0;i<n_node;i++)
  {
	  Node* current_Node=Bulk_mesh_pt->node_pt(i);
	  double u = current_Node->value(0);
	  current_Node->set_value(0,std::max(0.0,u));
  }
}
 */

//==========================Temporal error norm function===========================

template<class ELEMENT>
double BiasedActiveMatterDispersionProblem<ELEMENT>::global_temporal_error_norm()
{
 double global_error = 0.0;
   
 //Find out how many nodes there are in the problem
 unsigned n_node = Bulk_mesh_pt->nnode();
 
 //Loop over the nodes and calculate the estimated error in the values
 for(unsigned i=0;i<n_node;i++)
  {
   // Get error in solution: Difference between predicted and actual
   // value for nodal value 0
   double error = Bulk_mesh_pt->node_pt(i)->time_stepper_pt()->
    temporal_error_in_value(Bulk_mesh_pt->node_pt(i),0);
   
   //Add the square of the individual error to the global error
   global_error += error*error;
  }
    
 // Divide by the number of nodes
 global_error /= double(n_node);
 
 // Return square root...
 return sqrt(global_error);
 
} // end of global_temporal_error_norm




//====================start_of_actions_before_newton_solve=======================
/// Update the problem specs before solve (empty)
//========================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_before_newton_solve()
{
} // end of actions before solve

//========================set_initial_condition===========================
// New addition: Initial concentration distribution - 24/06/24
// From Bearon et al. (2012)
//========================================================================
 template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::
set_initial_condition()
{

 unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
  
   Node* nod_pt = mesh_pt()->node_pt(n);
   //Lets have a little blob of swimmers along the axis
   //Exponentially decaying in z direction; uniform in r
   double x = nod_pt->x(0);
   double y = nod_pt->x(1);
 
  
    /* double u = exp(-4*r*r)*
     exp( -(z-0.1*ParameterSpace::Length)*
                  (z-0.1*ParameterSpace::Length)/(0.01*0.01*ParameterSpace::Length*ParameterSpace::Length));
     */

	//double u = 1*exp(-0.1*y);
	
  double u=1.0; // Uniform initial distribution
 
   nod_pt->set_value(0,u);
  }
} 


//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 

 // Doc refinement levels in bulk mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 Bulk_mesh_pt->get_refinement_levels(min_refinement_level,
                                     max_refinement_level); 
 cout << "Ultimate min/max. refinement levels in bulk mesh : " 
      << min_refinement_level << " " 
      << max_refinement_level << std::endl;

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 /* // Output exact solution 
 //----------------------
  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,ParameterSpace::get_exact_u); 
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,ParameterSpace::get_exact_u,
                               error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
  cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;
 */

} // end of doc
 

//============start_of_actions_before_adapt==============================
// Remove the flux elements from the mesh.
//=======================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_before_adapt()
{
 // Rebuild the global mesh
 // 01/07/24: Probably not needed (didn't change the solution, will remove it in the future)
 //rebuild_global_mesh();
} // end of actions_before_adapt

//============start_of_actions_after_adapt==============================
// Attach flux elements to the mesh.
//=======================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_after_adapt()
{
	
 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Doc refinement levels in bulk mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 Bulk_mesh_pt->get_refinement_levels(min_refinement_level,
                                     max_refinement_level); 
 cout << "Min/max. refinement levels in bulk mesh: " 
      << min_refinement_level << " " 
      << max_refinement_level << std::endl;

} // end of actions_after_adapt

//==========start_of_main=================================================
/// Demonstrate how to solve 2D AdvectionDiffusion problem with flux boundary 
/// conditions, using two meshes.
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 //Set up the problem with 2D nine-node elements from the
 //RefineableQuadAdvectionDiffusionElement family. 
 BiasedActiveMatterDispersionProblem<RefineableQGeneralisedAdvectionDiffusionElement_time_wind<2,
  3> > problem(&ParameterSpace::source_function,
			   &ParameterSpace::wind_function,
			   &ParameterSpace::swimming,
			   &ParameterSpace::diff_function);
 
 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 //Output solution
 problem.doc_solution(doc_info);

 double dt = 0.01;
 double t_max=20;
  problem.initialise_dt(dt);
 // Set initial condition: From Bearon - 24/06/24
 problem.set_initial_condition();
 problem.doc_solution(doc_info);
 doc_info.number()++;
 //problem.assign_initial_values_impulsive(dt);

 bool first = true;
 unsigned max_adapt=4;
  double epsilon_t=1.0e-4;
 
  unsigned nstep = unsigned(t_max/dt);
unsigned istep=0;
 while (problem.time_pt()->time()<t_max)
  {
	// Solve the problem, allowing for up to four adaptations
	istep=istep+1;
	double dt_next=problem.adaptive_unsteady_newton_solve(dt,epsilon_t);
	dt=dt_next;
	if(problem.time_pt()->time()==0.0) 
	{
		first=false;
		//Output solution
		problem.doc_solution(doc_info);
 
		//Increment counter for solutions 
		doc_info.number()++; 
	}
	
	if (istep>0)
	{
		if (istep % 10 == 0)
		{
			//Output solution
			problem.doc_solution(doc_info);
		 
			//Increment counter for solutions 
			doc_info.number()++;
		}
	}
  } 
} //end of main

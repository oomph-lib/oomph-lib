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
//
//Flow in a torus computed by using a cylindrical polar coordinate system
//and assuming axisymmetry.

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "navier_stokes.h" // for preconditioner
#include "meshes/full_circle_mesh.h"

using namespace std;
using namespace oomph;

//=================================================================
//A namespace for the physical parameters in the problem
//=================================================================
namespace Global_Physical_Variables
{
 /// The Reynolds number
 double Re = 0.0;

 /// The radius of the torus
 double Radius = 1.0;

 /// The curvature of the torus
 double Delta = 0.128;
};



//=start_of_FillCircle=============================================
//A geometric object that represents the geometry of the domain
//a circle of given centre and radius. One could use a non-linear
//stretch in r (xi[0]) to shift the elements towards the edge
//(boundary layer).
//=================================================================
class FilledCircle : public GeomObject
{
public:

 /// Constructor that takes the centre position and raidus of the circle
 /// as its arguments
 FilledCircle(const double &centre_x, const double &centre_y,
              const double &radius) :
  GeomObject(2,2), Centre_x(centre_x), Centre_y(centre_y), Radius(radius) { }
 
/// Destructor
virtual ~FilledCircle(){}

///Lagrangian coordinate xi
void position (const Vector<double>& xi, Vector<double>& r) const
{
 r[0] = Centre_x + Radius*xi[0]*cos(xi[1]);
 r[1] = Centre_y + Radius*xi[0]*sin(xi[1]);
}


/// Return the position of the circle as a function of time 
/// (doesn't move as a function of time)
void position(const unsigned& t, 
              const Vector<double>& xi, Vector<double>& r) const
  {
   position(xi,r);
  }

private:

 ///Storage for the x-coordinate of the centre
 double Centre_x;
 
 ///Storage for the y-coordinate of the centre
 double Centre_y;

 ///Storage for the radius of the circle
 double Radius;

};


#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}
#endif


//==========================================================================
/// Solve the Axisymmetric Navier Stokes equations in a torus
//==========================================================================
template<class ELEMENT>
class TorusProblem : public Problem
{
public:
 /// Constructor taking the maximum refinement level and
 /// the minimum and maximum error targets.
 TorusProblem(const unsigned &max_refinement_level,
              const double &min_error_target, 
              const double &max_error_target);

/// Set the initial conditions: all nodes have zero velocity
void set_initial_condition() 
  {
   const unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     for(unsigned i=0;i<3;i++)
      {
       mesh_pt()->node_pt(n)->set_value(i,0.0);
      }
    }
  }

 /// Set boundary conditions on the walls
 void set_boundary_conditions(const double &time);

 /// Function that is used to run the parameter study
 void solve_system(const double &dt, const unsigned &nstep);
 
 /// Return a pointer to the specific mesh used
 RefineableFullCircleMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<RefineableFullCircleMesh<ELEMENT>*>(Problem::mesh_pt());}
 /// Update the problem specs before next timestep: 
 void actions_before_implicit_timestep() 
  {set_boundary_conditions(time());}

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   //Pin a single pressure value 
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
  }

};

//============================================================================
/// Constructor: specify the maximum refinement level, the minimum and
/// maximum error targets
//============================================================================
template<class ELEMENT>
TorusProblem<ELEMENT>::TorusProblem(const unsigned &max_refinement_level,
                                    const double &min_error_target,
                                    const double &max_error_target)
{
 
 using namespace Global_Physical_Variables;

 //Create a timestepper
 add_time_stepper_pt(new BDF<2>);

 //Create the domain for the mesh, which consists of a circle of
 //radius Radius and centred at (1/Delta, 0) 
 GeomObject* area_pt = new FilledCircle(1.0/Delta,0.0,Radius);

 //Define pi
 const double pi = MathematicalConstants::Pi;
 
 //Set the positions of the angles that divide the outer ring
 //These must be in the range -pi,pi, ordered from smallest to
 //largest
 Vector<double> theta_positions(4);
 theta_positions[0] = -0.75*pi;
 theta_positions[1] = -0.25*pi;
 theta_positions[2] = 0.25*pi;
 theta_positions[3] = 0.75*pi;

 //Define the radial fraction of the central box (always halfway
 //along the radius)
 Vector<double> radial_frac(4,0.5);

 //Now create the mesh
 Problem::mesh_pt() = new RefineableFullCircleMesh<ELEMENT>(
  area_pt,theta_positions,radial_frac,time_stepper_pt());

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Maximum number of refinements (increase this if you want a finer mesh)
 mesh_pt()->max_refinement_level() = max_refinement_level;

 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error() = max_error_target; 
 mesh_pt()->min_permitted_error() = min_error_target; 

 //Build iterative linear solver
 GMRES<CRDoubleMatrix>* iterative_linear_solver_pt = new
  GMRES<CRDoubleMatrix>;
 
 // Set maximum number of iterations
 iterative_linear_solver_pt->max_iter() = 100;
 
 // Set tolerance
 iterative_linear_solver_pt->tolerance() = 1.0e-8;   
 
 NavierStokesSchurComplementPreconditioner* prec_pt = 
  new NavierStokesSchurComplementPreconditioner(this);
 //Set the mesh
 prec_pt->set_navier_stokes_mesh(this->mesh_pt());

 //Set the preconditioner
 iterative_linear_solver_pt->preconditioner_pt() = prec_pt;
 
 //Set the linear solver
 this->linear_solver_pt() = iterative_linear_solver_pt;

 //Precondition the pressure block with AMG if we have Hypre
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when OOMPH_HAS_MPI, but we
//run in serial
#ifndef OOMPH_HAS_MPI

 //Set up the internal preconditioners
 Preconditioner* P_matrix_preconditioner_pt = new HyprePreconditioner;
 
 // Set parameters for use as preconditioner on Poisson-type problem
 Hypre_default_settings::set_defaults_for_2D_poisson_problem(
  static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
 
 // Use Hypre for the Schur complement block
 prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
 
 // Shut up!
 static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
  disable_doc_time();
#endif
#endif
 
//Set Block diagonal preconditioner for the momentum block
BlockDiagonalPreconditioner<CRDoubleMatrix>* F_matrix_preconditioner_pt = 
  new BlockDiagonalPreconditioner<CRDoubleMatrix>;

// Set mesh
F_matrix_preconditioner_pt->add_mesh(this->mesh_pt());

//Use AMG on the momentum blocks if we have Hypre 
#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
 dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
  (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
  (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
#endif
#endif
  
 // Use Hypre for momentum block 
 prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);

 //Loop over all the (fluid) elements 
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to the particular element type, this is necessary because
   //the base elements don't have the member functions that we're about
   //to call!
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //There is no need for ALE
   el_pt->disable_ALE();

   //Set the Reynolds number for each element 
   //(yes we could have different Reynolds number in each element!!)
   el_pt->re_pt() = &Re;
   //Set the product of Reynolds and Strouhal numbers
   el_pt->re_st_pt() = &Re;
  }

 //Let this problem be conventional form by setting gamma to zero
 ELEMENT::Gamma[0] = 0.0; //r-momentum
 ELEMENT::Gamma[1] = 0.0; //z-momentum
 
 //Set the boundary conditions (no slip on the torus walls)
 //Loop over the nodes on the (only) mesh boundary
 unsigned n_boundary_node = mesh_pt()->nboundary_node(0);
 for(unsigned n=0;n<n_boundary_node;n++)
  {
   //Pin all three velocity components on the wall
   for(unsigned i=0;i<3;i++)
    {
     mesh_pt()->boundary_node_pt(0,n)->pin(i);
    }
  }
 
 // Pin redudant pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 //Pin a single pressure value 
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
 
 //Setup all the equation numbering and look-up schemes 
 std::cout << assign_eqn_numbers() << std::endl; 
}

//========================================================================
/// Set the boundary conditions as a function of time we are going
/// to spin up the torus
//========================================================================
template<class ELEMENT>
void TorusProblem<ELEMENT>::set_boundary_conditions(const double &time)
{
 //NOTE: The default value of all parameters is zero, so we need only 
 //set the values that are non-zero on the boundary, i.e. the swirl

 //Loop over the nodes on the boundary
 unsigned n_boundary_node = mesh_pt()->nboundary_node(0);
 //Loop over the nodes on the boundary
 for(unsigned n=0;n<n_boundary_node;n++)
  {
   //Get the radial values
   double r = mesh_pt()->boundary_node_pt(0,n)->x(0);
   //Set the value of the w-velocity
   //Fast(ish) spin-up
   mesh_pt()->boundary_node_pt(0,n)->set_value(2,r*(1.0 - exp(-20.0*time)));
  }
}


//==========================================================================
///Solve the system for a number of different values of the Reynolds number
//==========================================================================
template<class ELEMENT>
void TorusProblem<ELEMENT>::solve_system(const double &dt, 
                                         const unsigned &nstep)
{
 using namespace Global_Physical_Variables;

 //Open a trace file
 ofstream trace("time_trace.dat");

 //Define a string that we can set to be the name of the output file
 char filename[100];
 //Define an output filestream
 ofstream file;

 //Set the Reynolds number
 Re = 1000.0;

 //Set an impulsive start from rest
 assign_initial_values_impulsive(dt);

 //Output intital data
 trace << time() << " " << mesh_pt()->boundary_node_pt(0,0)->value(2) 
       << " " << mesh_pt()->node_pt(0)->value(0) 
       << " " << mesh_pt()->node_pt(0)->value(1) 
       << " " << mesh_pt()->node_pt(0)->value(2) << std::endl;

 //Increase the maximum value of the residuals to get
 //past the first few steps
 Max_residuals = 50.0;

 //Now perform the first timestep with 2 steps of refinement
 unsteady_newton_solve(dt,2,true);

 trace << time() << " " << mesh_pt()->boundary_node_pt(0,0)->value(2) 
       << " " << mesh_pt()->node_pt(0)->value(0) 
       << " " << mesh_pt()->node_pt(0)->value(1) 
       << " " << mesh_pt()->node_pt(0)->value(2) << std::endl;
 
 //Output data after the first timestep
 //Create the filename, including the array index
 sprintf(filename,"soln_Re%g_t%g.dat",Re,time());
 //Actually, write the data
 file.open(filename);
 mesh_pt()->output(file,5);
 file.close();
 
 //Now do the other steps with only one adaptation per step
 for(unsigned n=1;n<nstep;n++)
  {
   //Solve the problem
   unsteady_newton_solve(dt,1,false);
 
   trace << time() << " " << mesh_pt()->boundary_node_pt(0,0)->value(2) 
         << " " << mesh_pt()->node_pt(0)->value(0) 
         << " " << mesh_pt()->node_pt(0)->value(1) 
         << " " << mesh_pt()->node_pt(0)->value(2) << std::endl;
  
   //Output data at each step
   //Create the filename, including the array index
   sprintf(filename,"soln_Re%g_t%g.dat",Re,time());
   //Actually, write the data
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();
  }

 //Close the trace file
 trace.close();
}

//Main driver loop
int main()
{
 //Construct and solve the problem
 //This maximum refinement level means that we fit (easily) into a 
 //1G machine, but it could probably go higher if you start to 
 //see refinement being overruled.
 //If you have this too high you get ridiculous refinement at early 
 //times that isn't really necessary.
 unsigned max_refinement_level = 6;
 double max_error = 1.0e-3;
 double min_error = 1.0e-5;
 TorusProblem<RefineableAxisymmetricQTaylorHoodElement> problem(
  max_refinement_level,min_error,max_error);

 //Refine uniformly once
 problem.refine_uniformly();

 //Now timestep
 problem.solve_system(0.01,2);
}









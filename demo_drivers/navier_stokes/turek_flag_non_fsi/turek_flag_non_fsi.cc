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
//------------------------------------------------------------------------
// Flow past cylinder with flag -- first draft developed by 
// Stefan Kollmannsberger and students Iason Papaioannou and 
// Orkun Oezkan Doenmez.
//
// http://www.inf.bauwesen.tu-muenchen.de/~kollmannsberger/
//
// Completed by Floraine Cordier.
//------------------------------------------------------------------------

// Generic oomph-lib stuff
#include "generic.h"

// Navier Stokes
#include "navier_stokes.h"

// The mesh
#include "meshes/cylinder_with_flag_mesh.h"

using namespace std;

using namespace oomph;

//#define USE_MACRO_ELEMENTS


//====start_of_global_parameters==========================================
/// Global parameters
//========================================================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re=100.0;
}



//====start_of_flag_definition===========================================
/// Namespace for definition of flag boundaries
//=======================================================================
namespace Flag_definition
{

 /// Period of prescribed flag oscillation
 double Period=10.0;

 /// Height of flag
 double H=0.2;

 /// Length of flag
 double L=3.5;
   
 /// x position of centre of cylinder
 double Centre_x=2.0;

 /// y position of centre of cylinder
 double Centre_y=2.0;

 /// Radius of cylinder
 double Radius=0.5;

 /// Amplitude of tip deflection
 double Amplitude=0.33;

 ///Pointer to the global time object
 Time* Time_pt=0;
 


 /// Time-dependent vector to upper tip of the "flag"
 Vector<double> upper_tip(const double& t)
  {
   double tmp_ampl=Amplitude;
   if (t<=0.0) tmp_ampl=0.0;
   Vector<double> uppertip(2);
   uppertip[0]= Centre_x+Radius*sqrt(1.0-H*H/(4.0*Radius*Radius))+L;
   uppertip[1]= Centre_y+0.5*H-
    tmp_ampl*sin(2.0*MathematicalConstants::Pi*t/Period);
   return uppertip;
  }  
 
 /// Time-dependent vector to bottom tip of the "flag"
 Vector<double> lower_tip(const double& t)
  {
   double tmp_ampl=Amplitude;
   if (t<=0.0) tmp_ampl=0.0;
   Vector<double> lowertip(2);
   lowertip[0]= Centre_x+Radius*sqrt(1.0-H*H/(4.0*Radius*Radius))+L;
   lowertip[1]= Centre_y-0.5*H-
    tmp_ampl*sin(2.0*MathematicalConstants::Pi*t/Period);
   return lowertip;
  } // end of bottom tip





//-----start_of_top_of_flag--------------------------------------
/// GeomObject that defines the upper boundary of the flag
//---------------------------------------------------------------
 class TopOfFlag : public GeomObject
 {
  
 public:
  
  ///Constructor: It's a 1D object in 2D space 
  TopOfFlag() : GeomObject(1,2) {}

  ///Destructor (empty)
  ~TopOfFlag(){}
  
  /// Return the position along the top of the flag (xi[0] varies 
  /// between 0 and Lx)
  void position(const unsigned& t,const Vector<double> &xi, Vector<double> &r)
   const
   {
    // Compute position of fixed point on the cylinder
    Vector<double> point_on_circle(2);     
    point_on_circle[0]=Centre_x+Radius*sqrt(1.0-H*H/(4.0*Radius*Radius));
    point_on_circle[1]=Centre_y+H/2.0;

    r[0] = point_on_circle[0]+xi[0]/L*(upper_tip(Time_pt->time(t))[0]-
                                        point_on_circle[0]);
    
    r[1] = point_on_circle[1]+xi[0]/L*(upper_tip(Time_pt->time(t))[1]-
                                        point_on_circle[1])+
     1.0/3.0*sin((r[0]-point_on_circle[0])/
                 (upper_tip(Time_pt->time(t))[0]-
                  point_on_circle[0])*MathematicalConstants::Pi)
     *sin(2.0* MathematicalConstants::Pi*Time_pt->time(t)/Period);
    
   }
 
  
  /// Current position
  void position(const Vector<double> &xi, Vector<double> &r) const
   {
    return position(0,xi,r);
   }

  /// Number of geometric Data in GeomObject: None.
  unsigned ngeom_data() const {return 0;} 
  
};



//-----start_of_bottom_of_flag-----------------------------------
/// GeomObject that defines the lower boundary of the flag
//---------------------------------------------------------------
 class BottomOfFlag : public GeomObject
 {
  
 public:
  
  ///Constructor: 
  BottomOfFlag() : GeomObject(1,2) {}
  
  ///Destructor (empty)
  ~BottomOfFlag(){}
  

  /// Return the position along the bottom of the flag (xi[0] varies 
  /// between 0 and Lx)
  void position(const unsigned& t,const Vector<double> &xi, Vector<double> &r)
   const
   {
    // Compute position of fixed point on the cylinder
    Vector<double> point_on_circle(2);     
    point_on_circle[0]=Centre_x+Radius*sqrt(1.0-H*H/(4.0*Radius*Radius));
    point_on_circle[1]=Centre_y-H/2.0;

    r[0] = point_on_circle[0]+ xi[0]/L*(lower_tip(Time_pt->time(t))[0]-
                                        point_on_circle[0]);
    
    r[1] = point_on_circle[1]+ xi[0]/L*(lower_tip(Time_pt->time(t))[1]-
                                        point_on_circle[1])+
     1.0/3.0*sin((r[0]-point_on_circle[0])/
                 (lower_tip(Time_pt->time(t))[0]-
                  point_on_circle[0])*MathematicalConstants::Pi)
     *sin(2.0*MathematicalConstants::Pi*Time_pt->time(t)/Period);
   }
 
  
  /// Current position
  void position(const Vector<double> &xi, Vector<double> &r) const
   {
    return position(0,xi,r);
   }

  /// Number of geometric Data in GeomObject: None.
  unsigned ngeom_data() const {return 0;} 
    
};




 //-----start_of_tip_of_flag--------------------------------------
 /// GeomObject that defines the tip of the flag
 //---------------------------------------------------------------
 class TipOfFlag : public GeomObject
  { 
   
  public:
   
   ///Constructor
   TipOfFlag() : GeomObject(1,2) {}
   
   ///Destructor (empty)
   ~TipOfFlag(){}
   
   /// Return the position
   /// This object links the tips of the top and bottom by a straight line
   /// whilst xi[0] goes from -H/2 to H/2.
   void position(const unsigned& t,const Vector<double> &xi, Vector<double> &r)
    const
    {
     Vector<double> point_top(upper_tip(Time_pt->time(t)));
     Vector<double> point_bottom(lower_tip(Time_pt->time(t)));
     
     r[1]= point_bottom[1]+(xi[0]+H/2.0)/H*(point_top[1]-point_bottom[1]);
     r[0]= point_bottom[0]+(xi[0]+H/2.0)/H*(point_top[0]-point_bottom[0]);
    }
   
   /// Current position.
   void position(const Vector<double> &xi, Vector<double> &r) const
    {
     return position(0,xi,r);
    }
   
   /// Number of geometric Data in GeomObject: None.
   unsigned ngeom_data() const {return 0;}
   
 };


 /// Pointer to GeomObject that bounds the upper edge of the flag
 TopOfFlag* Top_flag_pt=0;

 /// Pointer to GeomObject that bounds the bottom edge of the flag
 BottomOfFlag* Bottom_flag_pt=0;

 /// Pointer to GeomObject that bounds the tip edge of the flag
 TipOfFlag* Tip_flag_pt=0;

 /// Pointer to GeomObject of type Circle that defines the
 /// central cylinder.
 Circle* Cylinder_pt=0;

 /// Create all GeomObjects needed to define the cylinder and the flag
 void setup(Time* time_pt)
 {
  // Assign pointer to global time object
  Time_pt=time_pt;
  
  // Create GeomObject of type Circle that defines the
  // central cylinder.
  Cylinder_pt=new Circle(Centre_x,Centre_y,Radius);
  
  /// Create GeomObject that bounds the upper edge of the flag
  Top_flag_pt=new TopOfFlag;
  
  /// Create GeomObject that bounds the bottom edge of the flag
  Bottom_flag_pt=new BottomOfFlag;
  
  /// Create GeomObject that bounds the tip edge of the flag
  Tip_flag_pt=new TipOfFlag;
 
 }

} // end of namespace that specifies the flag



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//======= start_of_problem_class=====================================
/// Flow around a cylinder with flag
//===================================================================
template<class ELEMENT>
class TurekNonFSIProblem : public Problem
{

public:

 /// Constructor: Pass length and height of domain.
 TurekNonFSIProblem(const double &length, 
                         const double &height);
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// After adaptation: Unpin pressures and pin redudant pressure dofs.
 void actions_after_adapt();

 /// Update the velocity boundary condition on the flag
 void actions_before_implicit_timestep();


#ifdef USE_MACRO_ELEMENTS

 /// Access function for the specific mesh
 RefineableCylinderWithFlagMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableCylinderWithFlagMesh<ELEMENT>*>
    (Problem::mesh_pt());
  }

#else

 /// Access function for the specific mesh
 RefineableAlgebraicCylinderWithFlagMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableAlgebraicCylinderWithFlagMesh<ELEMENT>*>
    (Problem::mesh_pt());
  }
 
#endif

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


private:

 /// Height of channel
 double Height;

};//end_of_problem_class




//=====start_of_constructor===============================================
/// Constructor: Pass length and height of domain.
//========================================================================
template<class ELEMENT>
TurekNonFSIProblem<ELEMENT>::TurekNonFSIProblem(
 const double &length, const double &height) : Height(height)
{ 

 // Bump up Newton solver parameters to allow for crappy initial guesses
 Max_residuals=100.0;;
 Max_newton_iterations=50;

 // Allocate the timestepper
 add_time_stepper_pt(new BDF<2>);

 // Setup flag/cylinder geometry
 Flag_definition::setup(time_pt());

#ifdef USE_MACRO_ELEMENTS



 // Build mesh (with Domain/MacroElement-based node update)
 Problem::mesh_pt()=new RefineableCylinderWithFlagMesh<ELEMENT>
  (Flag_definition::Cylinder_pt,
   Flag_definition::Top_flag_pt, 
   Flag_definition::Bottom_flag_pt, 
   Flag_definition::Tip_flag_pt,
   length, height,
   Flag_definition::L,
   Flag_definition::H, 
   Flag_definition::Centre_x,
   Flag_definition::Centre_y,
   Flag_definition::Radius,
   time_stepper_pt());

 std::cout << "Using Domain/MacroElement-based node update" << std::endl;

#else

 // Build mesh (with AlgebraicMesh-based node update)
 Problem::mesh_pt()=new RefineableAlgebraicCylinderWithFlagMesh<ELEMENT>
  (Flag_definition::Cylinder_pt,
   Flag_definition::Top_flag_pt, 
   Flag_definition::Bottom_flag_pt, 
   Flag_definition::Tip_flag_pt,
   length, height, 
   Flag_definition::L,
   Flag_definition::H, 
   Flag_definition::Centre_x,
   Flag_definition::Centre_y,
   Flag_definition::Radius,
   time_stepper_pt());

 std::cout << "Using AlgebraicMesh-based node update" << std::endl;

#endif

 // Refine uniformly
 mesh_pt()->refine_uniformly();
 mesh_pt()->refine_uniformly();

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 
 //Pin both velocities at all boundaries apart from outflow
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parallel, axially traction free outflow at downstream end
     if (ibound != 1)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     else
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  } // done bc
  
  
  // Pass pointer to Reynolds number to elements
  unsigned nelem=mesh_pt()->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
    
    //Set the Reynolds number
    el_pt->re_pt() = &Global_Parameters::Re;   
    
    //Set the Womersley number (assuming St=1)
    el_pt->re_st_pt() = &Global_Parameters::Re;
   }
  

  // Setup Poiseuille flow along boundary 3
  unsigned ibound=3;
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    double ycoord = mesh_pt()->boundary_node_pt(ibound,inod)->x(1);

    // Set Poiseuille velocity
    double uy = 6.0*ycoord/Height*(1.0-ycoord/Height);
    
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,uy);
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
   }   

  
  // Pin redudant pressure dofs
  RefineableNavierStokesEquations<2>::
   pin_redundant_nodal_pressures(mesh_pt()->element_pt());
  
  // Assign equations numbers
  cout <<"Number of equations: " << assign_eqn_numbers() << endl;  
 
} // end of constructor



//=====actions_after_adapt================================================
/// Actions after adapt
//========================================================================
template <class ELEMENT> 
void TurekNonFSIProblem<ELEMENT>::actions_after_adapt()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(mesh_pt()->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
} // end_of_actions_after_adapt


//==== start_of_actions_before_implicit_timestep==========================
/// Actions before implicit timestep: Update velocity boundary conditions
//========================================================================
template <class ELEMENT> 
void TurekNonFSIProblem<ELEMENT>::actions_before_implicit_timestep()
{

 // Update the domain shape
 mesh_pt()->node_update();

 // Moving leaflet: No slip; this implies that the velocity needs
 // to be updated in response to flag motion
 for( unsigned ibound=5;ibound<8;ibound++)
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Which node are we dealing with?
     Node* node_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     
     // Apply no slip
     FSI_functions::apply_no_slip_on_moving_wall(node_pt);
    }
  }
} //end_of_actions_before_implicit_timestep


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TurekNonFSIProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver code -- pass a command line argument if you want to run
/// the code in validation mode where it only performs a few steps
//=====================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Set up doc info
 DocInfo doc_info; 
 doc_info.set_directory("RESLT");
 doc_info.number()=0;
 
 // Length and height of domain
 double length=25.0;
 double height=4.1;

#ifdef USE_MACRO_ELEMENTS

 // Create Problem with Domain/MacroElement-based node update
 TurekNonFSIProblem<RefineableQCrouzeixRaviartElement<2> >  
  problem(length,height);

#else

 // Create Problem with AlgebraicMesh-based node update
 TurekNonFSIProblem
  <AlgebraicElement<RefineableQCrouzeixRaviartElement<2> > > 
  problem(length, height);

#endif

 
// Number of timesteps per period
 unsigned nsteps_per_period=40;

 // Number of periods
 unsigned nperiod=3; 

 // Number of timesteps  (reduced for validation)
 unsigned nstep=nsteps_per_period*nperiod;
 if (CommandLineArgs::Argc>1)
  {
   nstep=2;
   // Also reduce the Reynolds number to reduce the mesh refinement
   Global_Parameters::Re=10.0;
  }

 //Timestep: 
 double dt=Flag_definition::Period/double(nsteps_per_period);
 
 /// Initialise timestep 
 problem.initialise_dt(dt);

 // Solve adaptively with up to max_adapt rounds of refinement (fewer if 
 // run during self-test)
 unsigned max_adapt=3;
 if (CommandLineArgs::Argc>1)
  {
   max_adapt=1;
  } 

 /// Output intial guess for steady Newton solve
 problem.doc_solution(doc_info);
 doc_info.number()++;

 // Do steady solve first -- this also sets the history values
 // to those corresponding to an impulsive start from the
 // steady solution
 problem.steady_newton_solve(max_adapt);
 
 /// Output steady solution = initial condition for subsequent unsteady solve
 problem.doc_solution(doc_info);
 doc_info.number()++;

 /// Reduce the max number of adaptations for time-dependent simulation
 max_adapt=1;

 // We don't want to re-assign the initial condition after the mesh
 // adaptation
 bool first=false;
 
// Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  { 
   // Solve the problem
   problem.unsteady_newton_solve(dt,max_adapt,first);
   
   // Output the solution
   problem.doc_solution(doc_info);
   
   // Step number
   doc_info.number()++;
  }
 
}//end of main




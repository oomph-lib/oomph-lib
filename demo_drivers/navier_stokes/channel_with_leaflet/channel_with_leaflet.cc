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

// Generic oomph-lib includes
#include "navier_stokes.h"

//Include the mesh
#include "meshes/channel_with_leaflet_mesh.h"

using namespace std;
using namespace oomph;


//===start_of_leaflet_class===========================================
/// GeomObject representing a vertical leaflet that performs
/// bending and stretching oscillations. 
//====================================================================
class Leaflet : public GeomObject
{

public:

 /// Constructor: Pass length (in Lagrangian coordinates),
 /// the amplitude of the horizontal and vertical deflection of the tip,
 /// the x-coordinate of the origin and the period of the oscillation.
 /// Passes the number of Lagrangian and Eulerian coordinates to the
 /// constructor of the GeomObject base class.
 Leaflet(const double& length, const double& d_x,const double& d_y,
         const double& x_0, const double& period, Time* time_pt)
  : GeomObject(1,2), Length(length), D_x(d_x), D_y(d_y), X_0(x_0),
   T(period),Time_pt(time_pt) {}


 /// Destructor -- emtpy
 virtual ~Leaflet(){}

 /// Position vector, r, to the point identified by  
 /// its 1D Lagrangian coordinate, xi (passed as a 1D Vector) at discrete time
 /// level t (t=0: present; t>0: previous).
 void position(const unsigned& t, const Vector<double>& xi, 
               Vector<double>& r) const
  {
    using namespace MathematicalConstants;

   //Position
    r[0] =  X_0 + D_x*xi[0]*xi[0]/Length/Length*sin(2.0*Pi*Time_pt->time(t)/T);
    r[1] = xi[0]*(1.0+D_y/Length*0.5*(1.0-cos(4.0*Pi*Time_pt->time(t)/T)));
  }
 
 /// Steady version: Get current shape
 void position(const Vector<double>& xi, Vector<double>& r) const
  {
   position(0,xi,r); 
  }

 /// Number of geometric Data in GeomObject: None.
 unsigned ngeom_data() const {return 0;}  

 /// Length of the leaflet
 double length() { return Length; }

 /// Amplitude of horizontal tip displacement
 double& d_x() {return D_x;}

 /// Amplitude of vertical tip displacement
 double d_y() {return D_y;}

 /// x-coordinate of leaflet origin
 double x_0() {return X_0;}
 
private :
 
 /// Length in terms of Lagrangian coordinates
 double Length;
 
 ///Horizontal displacement of tip
 double D_x;
 
 ///Vertical displacement of tip
 double D_y;
 
///Origin
 double X_0;
 
 ///Period of the oscillations
 double T;

 ///Pointer to the global time object
 Time* Time_pt;

}; //end_of_the_GeomObject



///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//==start_of_global_parameters=======================================
/// Global parameters
//===================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=20.0;

} // end_of_namespace



///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//==start_of_problem_class===========================================
/// Problem class
//===================================================================
template<class ELEMENT>
class ChannelWithLeafletProblem : public Problem
{

public:

 /// Constructor: Pass the length of the domain at the left
 /// of the leaflet lleft,the length of the domain at the right of the
 /// leaflet lright,the height of the leaflet hleaflet, the total height
 /// of the domain htot, the number of macro-elements at the left of the
 /// leaflet nleft, the number of macro-elements at the right of the
 /// leaflet nright, the number of macro-elements under hleaflet ny1,
 /// the number of macro-elements above hleaflet ny2,the x-displacement
 /// of the leaflet d_x,the y-displacement of the leaflet d_y,the abscissa 
 /// of the origin of the leaflet x_0, the period of the moving leaflet.
 ChannelWithLeafletProblem(const double& l_left,
                           const double& l_right, const double& h_leaflet,
                           const double& h_tot,
                           const unsigned& nleft, const unsigned& nright,
                           const unsigned& ny1, const unsigned&  ny2,
                           const double& d_x,const double& d_y,
                           const double& x_0, const double& period);  

 /// Destructor (empty)
 ~ChannelWithLeafletProblem(){}

 /// Overloaded access function to specific mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve(){}

 /// Update before solve (empty)
 void actions_before_newton_solve(){}

 /// Actions after adaptation: Pin redundant pressure dofs
 void actions_after_adapt();

 /// Update the velocity boundary condition on the moving leaflet
 void actions_before_implicit_timestep();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to the GeomObject
 GeomObject* Leaflet_pt;
 
};




//==start_of_constructor=================================================
/// Constructor
//=======================================================================
template <class ELEMENT>
ChannelWithLeafletProblem<ELEMENT>::ChannelWithLeafletProblem(
 const double& l_left,
 const double& l_right, const double& h_leaflet,
 const double& h_tot,
 const unsigned& nleft, const unsigned& nright,
 const unsigned& ny1, const unsigned&  ny2,
 const double& d_x,const double& d_y,
 const double& x_0, const double& period)
{
 // Allocate the timestepper
 add_time_stepper_pt(new BDF<2>);
 
 //Create the geometric object that represents the leaflet
 Leaflet_pt = new Leaflet(h_leaflet, d_x, d_y, x_0, period, time_pt()) ;

//Build the mesh
Problem::mesh_pt()=new RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>(
 Leaflet_pt,
 l_left, l_right,
 h_leaflet,
 h_tot,nleft,
 nright,ny1,ny2,
 time_stepper_pt());     


 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;

 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  
   // Set the Womersley number (product of Reynolds and Strouhal).
   // We're assuming a Strouhal number of one, corresponding to
   // a non-dimensionalisation of time on the flow's natural timescale.
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;
  
  } // end loop over elements


 //Pin the boundary nodes
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      //do not pin the x velocity of the outflow
      if( ibound != 1)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
      }     
    }
  }
 
 // Setup parabolic flow along the inflow boundary 3
 unsigned ibound=3; 
 unsigned num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = mesh_pt()->boundary_node_pt(ibound,inod)->x(1); 
   double uy = 6.0*(ycoord/h_tot)*(1.0-(ycoord/h_tot));

   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,uy);
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }// end of setup boundary condition
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Problem::mesh_pt()->element_pt());
 
 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 
 
}//end of constructor




//=====start_of_actions_before_implicit_timestep=========================
/// Actions before implicit timestep: Update domain shape and
/// the velocity boundary conditions
//=======================================================================
template <class ELEMENT> 
void ChannelWithLeafletProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Update the domain shape
 mesh_pt()->node_update();

 // Moving leaflet: No slip; this implies that the velocity needs
 // to be updated in response to leaflet motion
 for( unsigned ibound=4;ibound<6;ibound++)
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



//==========start_of_actions_after_adaptation============================
// Actions after adaptation: Pin redundant pressure dofs
//=======================================================================
template<class ELEMENT>
void ChannelWithLeafletProblem<ELEMENT>::actions_after_adapt()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(mesh_pt()->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
} // end_of_actions_after_adapt



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ChannelWithLeafletProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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

 // Output boundaries
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();

} // end_of_doc_solution


///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
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
 
 // Parameters for the leaflet
 //----------------------------
 
 // Height
 double h_leaflet = 0.5;


 // Tip deflection
 double d_x = 0.25;
 double d_y = -0.05;

 // x-positon of root
 double x_0 = 3.0; 

 // Period of the oscillation on the natural timescale of the flow
 double period = 20.0;


 //Parameters for the domain
 //-------------------------

 // Length of the mesh to right and left of the leaflet
 double l_left =2.0;
 double l_right= 3.0;

 // Total height of domain (unity because lengths have been scaled on it)
 double h_tot=1.0;

 // Initial number of element rows/columns in various mesh regions
 unsigned nleft=8;
 unsigned nright=12;
 unsigned ny1=2;
 unsigned ny2=2; 
 
 //Build the problem
 ChannelWithLeafletProblem<AlgebraicElement<RefineableQTaylorHoodElement<2> > >
  problem(l_left, l_right, h_leaflet,
          h_tot,nleft, nright,ny1,ny2,
          d_x, d_y, x_0,
          period);
 
 
 // Number of timesteps per period
 unsigned nsteps_per_period=40;

 // Number of periods
 unsigned nperiod=3; 

 // Number of timesteps (reduced for validation)
 unsigned nstep=nsteps_per_period*nperiod;
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }

 //Timestep: 
 double dt=period/double(nsteps_per_period);
 
 /// Initialise timestep 
 problem.initialise_dt(dt);


 /// Set max. number of adaptations (reduced for validation)
 unsigned max_adapt=5;
 if (CommandLineArgs::Argc>1)
  {
   max_adapt=2;
  }

 // Do steady solve first -- this also sets the history values
 // to those corresponding to an impulsive start from the
 // steady solution
 problem.steady_newton_solve(max_adapt);
 
 /// Output steady solution
 problem.doc_solution(doc_info);
 doc_info.number()++;


 /// Reduce the max number of adaptations for time-dependent simulation
 max_adapt=1;

 // We don't want to re-assign the initial condition 
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



//A demo driver that solves the classic fluid flow problem of flow
//of a fluid film along an inclined plane. Stability analysis performed by 
//Yih (1963), Benjamin (1957), etc.

//This is an example of the subtleties involved in even a seemingly simple
//free surface problem.

//Standard C++ library includes
#include <iostream>
#include <fstream>
#include <cmath>

//Finite-Element library routines
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "fluid_interface.h"
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//The global physical variables
namespace Global_Physical_Variables
{
 ///Reynolds number, based on the velocity of the free-surface of a vertical 
 ///film
 double Re=0.0;

 ///The product of Reynolds number and inverse Froude number 
 ///is set to two in this problem, a natural velocity 
 ///scale.
 double ReInvFr=2.0;

 ///Angle of incline of the slope (45 degrees)
 double Alpha = 1.0*atan(1.0);

 ///The Vector direction of gravity, set to zero by default
 Vector<double> G(2,0.0);

 ///The Capillary number
 double Ca= 2.0;

 ///Direction of the wall normal vector
 Vector<double> Wall_normal;

 /// \short Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }

 ///Set the wavenumber
 double K = 0.5;

 ///Set the number of waves desired in the domain
 double N_wave = 2;

 ///The length of the domain to fit the desired number of waves
 double Length = 2*N_wave*4.0*atan(1.0)/K;

 /// Function that prescribes the hydrostatic pressure field at the outlet
 void hydrostatic_pressure_outlet(const double& time, const Vector<double> &x, 
                           Vector<double> &traction)
 {
  traction[0] = ReInvFr*G[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet(const double& time, const Vector<double> &x, 
                                 Vector<double> &traction)
 {
  traction[0] = -ReInvFr*G[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }
 //end of traction functions

 ///Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 // MH: No longer needed
 /// Timescale ratio (density) for (pseudo-)solid
 //double Lambda_sq=0.0;

}

//=====================================================================
/// A simple class that adds additional functionality to the basic
/// rectangular mesh class to create the appropriate FaceElements
/// for this problem.
///====================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class InclinedPlaneMesh : 
 public SimpleRectangularQuadMesh<ELEMENT>
{
 public:
 ///Constructor, do nothing
 InclinedPlaneMesh(const unsigned &nx, const unsigned &ny,
                   const double &lx, const double &ly,
                   TimeStepper* time_stepper_pt) :
  SimpleRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt)
  { } //end of constructor

 ///\short Function to add the traction boundary elements to boundaries
 /// 3(inlet) and 1(outlet) of the mesh
 void make_traction_elements()
  {
   //Inlet boundary conditions (boundary 3)
   {
    unsigned b = 3;
    //Find the number of elements adjacent to mesh boundary
    unsigned n_boundary_element = this->nboundary_element(b);
    //Loop over these elements and create the traction elements
    for(unsigned e=0;e<n_boundary_element;e++)
     {
      NavierStokesTractionElement<ELEMENT> *surface_element_pt =
       new NavierStokesTractionElement<ELEMENT>
       (this->boundary_element_pt(b,e),this->face_index_at_boundary(b,e));
      //Add to mesh's elements
      this->Element_pt.push_back(surface_element_pt);
      //Set the traction function
      surface_element_pt->traction_fct_pt() = 
      &Global_Physical_Variables::hydrostatic_pressure_inlet;
     }
   }

   //Outlet boundary conditions (boundary 1)
   {
    unsigned b=1;
    //Find the number of elements adjacent to mesh boundary
    unsigned n_boundary_element = this->nboundary_element(b);
    //Loop over these elements and create the traction elements
    for(unsigned e=0;e<n_boundary_element;e++)
     {
      NavierStokesTractionElement<ELEMENT> *surface_element_pt =
       new NavierStokesTractionElement<ELEMENT>
       (this->boundary_element_pt(b,e),this->face_index_at_boundary(b,e));
      //Add to mesh's elements
      this->Element_pt.push_back(surface_element_pt);
      //Set the traction function
      surface_element_pt->traction_fct_pt() = 
       &Global_Physical_Variables::hydrostatic_pressure_outlet;
     }
   }
  } //end of make_traction_elements
  
 //Make the free surface elements on the top surface
 void make_free_surface_elements()
  {
   //The free surface is on the boundary 2
   unsigned b = 2;
   unsigned n_boundary_element = this->nboundary_element(b);
   for(unsigned e=0;e<n_boundary_element;e++)
    {
     INTERFACE_ELEMENT *surface_element_pt =
      new INTERFACE_ELEMENT
      (this->boundary_element_pt(b,e),this->face_index_at_boundary(b,e));
     //Add to mesh's elements
     this->Element_pt.push_back(surface_element_pt);
     //Assign the capillary number to the free surface
     surface_element_pt->ca_pt() = 
      &Global_Physical_Variables::Ca;
     
     //Make a point element from left-hand side of the 
     //first surface element
     if(e==0)
      {
       FluidInterfaceBoundingElement* point_element_pt =
        surface_element_pt->make_bounding_element(-1);
       //Add ot mesh's elements
       this->Element_pt.push_back(point_element_pt);
       //Set the capillary number
       point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;
       point_element_pt->wall_unit_normal_fct_pt() = 
        &Global_Physical_Variables::wall_unit_normal_fct;
      }
      
     //Make another point element from the right-hand side of the 
     //last surface element
     if(e==n_boundary_element-1)
      {
       FluidInterfaceBoundingElement* point_element_pt =
        surface_element_pt->make_bounding_element(1);

       //Add to the mesh's elements
       this->Element_pt.push_back(point_element_pt);

       //Set the capillary number
       point_element_pt->ca_pt() = &Global_Physical_Variables::Ca;

       // Set the function that specifies the wall normal
       point_element_pt->wall_unit_normal_fct_pt() = 
        &Global_Physical_Variables::wall_unit_normal_fct;
      }
    }
  } //end of make_free_surface_elements
};


//=====================================================================
///Generic problem class
//====================================================================
template<class ELEMENT>
class InclinedPlaneProblem : public Problem
{
protected:

 ///Prefix for output files
 std::string Output_prefix;

public:

 ///Constructor
 InclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length) :
  Output_prefix("Unset") { }
 
 ///Solve the steady problem
 void solve_steady();
 
 ///Take n_tsteps timesteps of size dt
 void timestep(const double &dt, const unsigned &n_tsteps);

 ///\short Complete the build of the problem setting all standard
 ///parameters and boundary conditions
 void complete_build(const unsigned &nx, const unsigned &ny)
  {
   using namespace Global_Physical_Variables;
   
   //Complete the build of the fluid elements by passing physical parameters
   //Find the number of elements
   unsigned long n_element = nx*ny;
   //Loop over all the fluid elements 
   for(unsigned long e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
     
     //Set the Reynolds number, etc
     temp_pt->re_pt() = &Re;
     temp_pt->re_st_pt() = &Re;
     temp_pt->re_invfr_pt() = &ReInvFr;
     temp_pt->g_pt() = &G;

     //Assign the global time
     temp_pt->time_pt() = time_pt();
    }
   
   //------------Set the boundary conditions for this problem----------

   {
    //Determine whether we are solving an elastic problem
    bool elastic = false;
    if(dynamic_cast<SolidNode*>(mesh_pt()->node_pt(0))) {elastic=true;}

    //Loop over the bottom of the mesh (the wall of the channel)
    unsigned n_node = mesh_pt()->nboundary_node(0);
    for(unsigned j=0;j<n_node;j++)
     {
      //Pin the u- and v- velocities
      mesh_pt()->boundary_node_pt(0,j)->pin(0);
      mesh_pt()->boundary_node_pt(0,j)->pin(1);

      //If we are formulating the elastic problem pin both positions
      //of nodes
      if(elastic)
       {
        static_cast<SolidNode*>(mesh_pt()->boundary_node_pt(0,j))
         ->pin_position(0);
        static_cast<SolidNode*>(mesh_pt()->boundary_node_pt(0,j))
         ->pin_position(1);
       }
     }
    
    //Loop over the inlet and set the Dirichlet condition
    //of no vertical velocity
    n_node = mesh_pt()->nboundary_node(3);
    for(unsigned j=0;j<n_node;j++)
     {
      mesh_pt()->boundary_node_pt(3,j)->pin(1);

      //If elastic pin horizontal position
      if(elastic)
       { 
        static_cast<SolidNode*>(mesh_pt()->boundary_node_pt(3,j))
         ->pin_position(0);
       }
     }
    
    
    //Loop over the outlet and set the Dirichlet condition
    //of no vertical velocity
    n_node = mesh_pt()->nboundary_node(1);
    for(unsigned j=0;j<n_node;j++)
     {
      mesh_pt()->boundary_node_pt(1,j)->pin(1);

      //If elastic pin horizontal position
      if(elastic)
       { 
        static_cast<SolidNode*>(mesh_pt()->boundary_node_pt(1,j))
         ->pin_position(0);
       }
     }
   }
 
   //Attach the boundary conditions to the mesh
   cout << assign_eqn_numbers() << " in the main problem" << std::endl; 
  } //end of complete_build
 
};


//-------------------------------------------------------------------------
///Solve the steady problem
//-------------------------------------------------------------------------
template<class ELEMENT>
void InclinedPlaneProblem<ELEMENT>::solve_steady()
{
 //Load the namespace
 using namespace Global_Physical_Variables;
 
 //Set all nodes to the Nusselt flat-film solution
 {
  unsigned n_node = mesh_pt()->nnode();
  for(unsigned long i=0;i<n_node;i++)
   {
    double y = mesh_pt()->node_pt(i)->x(1);
    //Top row
    mesh_pt()->node_pt(i)->set_value(0,0.5*ReInvFr*sin(Alpha)*(2.0*y - y*y));
   }
 }
 
 //Do one steady solve
 steady_newton_solve();
 
 //Output the full flow field
 char filename[100];
 sprintf(filename,"%s_output.dat",Output_prefix.c_str());
 ofstream file(filename);
 mesh_pt()->output(file,5);
 file.close();
} //end of solve_steady


//----------------------------------------------------------------------
/// Perform n_tsteps timesteps of size dt
//----------------------------------------------------------------------
template<class ELEMENT>
void InclinedPlaneProblem<ELEMENT>::
timestep(const double &dt, const unsigned &n_tsteps)
{
 //Need to use the Global variables here
 using namespace Global_Physical_Variables;
 
 //Open an output file
 char filename[100];
 sprintf(filename,"%s_time_trace.dat",Output_prefix.c_str());
 ofstream trace(filename); 
 //Counter that will be used to output the full flowfield
 //at certain timesteps
 int counter=0; 
 
 //Initial output of the time and the value of the vertical position at the
 //left and right-hand end of the free surface
 trace << time_pt()->time() << " " 
       << mesh_pt()->boundary_node_pt(2,0)->value(1) 
       << " "
       <<  mesh_pt()->boundary_node_pt(2,mesh_pt()->nboundary_node(2)-1)->x(1) 
       << " "
       << std::endl;
 
 //Loop over the desired number of timesteps
 for(unsigned t=1;t<=n_tsteps;t++)
  {
   //Increase the counter
   counter++;
   cout << std::endl;
   cout << "--------------TIMESTEP " << t<< " ------------------" << std::endl;
   
   //Take a timestep of size dt
   unsteady_newton_solve(dt);
   
   //Uncomment to get full solution output
   if(counter==2) //Change this number to get output every n steps
    {
     ofstream file;
     char filename[100];
     sprintf(filename, "%s_step%g_%i.dat",Output_prefix.c_str(),Re,t);
     
     file.open(filename);
     mesh_pt()->output(file,5);
     file.close();
     
     counter=0;
    }
   
   //Output the time and value of the vertical position of the free surface
   //at the left- and right-hand ends
   trace << time_pt()->time() << " "
         << mesh_pt()->boundary_node_pt(2,0)->value(1) << " "
         << 
    mesh_pt()->boundary_node_pt(2,mesh_pt()->nboundary_node(2)-1)->x(1) << " "
         << std::endl;

  }
} //end of timestep






//======================================================================
/// Convert the basic inclined plane mesh into an spine mesh
//======================================================================
template <class ELEMENT>
class SpineInclinedPlaneMesh : 
 public InclinedPlaneMesh<ELEMENT,SpineLineFluidInterfaceElement<ELEMENT> >, 
 public SpineMesh
{
public:
 SpineInclinedPlaneMesh(const unsigned &nx, const unsigned &ny,
                        const double &lx, const double &ly,
                        TimeStepper* time_stepper_pt) :
  InclinedPlaneMesh<ELEMENT,SpineLineFluidInterfaceElement<ELEMENT> >
 (nx,ny,lx,ly,time_stepper_pt), SpineMesh()
  {
   //Find the number of linear points in the element
   unsigned n_p =  dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();
   //Reserve storage for the number of spines
   Spine_pt.reserve((n_p-1)*nx + 1);
   
   //Create single pointer to a spine
   Spine* new_spine_pt=0;

   //Now loop over the elements horizontally
   for(unsigned long j=0;j<nx;j++)
    {
     //In most elements, we don't assign a spine to the last column,
     //beacuse that will be done by the next element
     unsigned n_pmax = n_p-1;
     //In the last element, however, we must assign the final spine
     if(j==nx-1) {n_pmax = n_p;}

     //Loop over all nodes horizontally
     for(unsigned l2=0;l2<n_pmax;l2++)
      {
       //Create a new spine with unit height and add to the mesh
       new_spine_pt=new Spine(1.0);
       Spine_pt.push_back(new_spine_pt);
       
       // Get the node
       SpineNode* nod_pt=element_node_pt(j,l2);
       //Set the pointer to spine
       nod_pt->spine_pt() = new_spine_pt;
       //Set the fraction
       nod_pt->fraction() = 0.0;
       // Pointer to the mesh that implements the update fct
       nod_pt->spine_mesh_pt() = this; 
       
       //Loop vertically along the spine
       //Loop over the elements 
       for(unsigned long i=0;i<ny;i++)
        {
         //Loop over the vertical nodes, apart from the first
         for(unsigned l1=1;l1<n_p;l1++)
          {
           // Get the node
           SpineNode* nod_pt=element_node_pt(i*nx+j,l1*n_p+l2);
           //Set the pointer to the spine
           nod_pt->spine_pt() = new_spine_pt;
           //Set the fraction
           nod_pt->fraction()=(double(i)+double(l1)/double(n_p-1))/double(ny);
           // Pointer to the mesh that implements the update fct
           nod_pt->spine_mesh_pt() = this; 
          }  
        }
      }
    } //End of horizontal loop over elements  
   
   //Now we can make the traction elements
   this->make_traction_elements();

   //Now we can make the free surface elements
   this->make_free_surface_elements();
  } //end of constructor

 /// \short General node update function implements pure virtual function 
 /// defined in SpineMesh base class and performs specific node update
 /// actions:  along vertical spines
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of y
   spine_node_pt->x(1) = W*H;
  }
};



//============================================================================
//Specific class for inclined plane problem using spines
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class SpineInclinedPlaneProblem : 
 public InclinedPlaneProblem<ELEMENT>
{
public:

 //Constructor
 SpineInclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length): 
  InclinedPlaneProblem<ELEMENT>(nx,ny,length) 
  {
   //Set the name
   this->Output_prefix = "spine";

   //Create our one and only timestepper, with adaptive timestepping
   add_time_stepper_pt(new TIMESTEPPER);

   Problem::mesh_pt() = new  SpineInclinedPlaneMesh<ELEMENT>(
    nx,ny,length,1.0,this->time_stepper_pt());

   this->complete_build(nx,ny);
  }
 
 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check()
  {this->mesh_pt()->node_update();}

 //Destructor, destory mesh and timestepper
 ~SpineInclinedPlaneProblem()
  {
   delete this->mesh_pt();
   delete this->time_stepper_pt();
  }
};


//====================================================================
// Elastic formulation of the problem
//===================================================================


//======================================================================
/// Convert the basic inclined plane mesh into an Elastic mesh
//======================================================================
template <class ELEMENT>
class ElasticInclinedPlaneMesh : 
 public InclinedPlaneMesh<ELEMENT,ElasticLineFluidInterfaceElement<ELEMENT> >, 
 public SolidMesh
{
 //Public functions
 public:
 ElasticInclinedPlaneMesh(const unsigned &nx, const unsigned &ny,
                          const double &lx, const double &ly,
                          TimeStepper* time_stepper_pt) :
  InclinedPlaneMesh<ELEMENT,ElasticLineFluidInterfaceElement<ELEMENT> >
 (nx,ny,lx,ly,time_stepper_pt), SolidMesh()
  {
   //Now we can make the traction elements
   this->make_traction_elements();

   //Now we can make the free surface elements
   this->make_free_surface_elements();

   //Make the current configuration the undeformed one
   set_lagrangian_nodal_coordinates();
  }
};




//============================================================================
//Specific class for inclined plane problem using pseudo-elastic formulation
//============================================================================
template<class ELEMENT, class TIMESTEPPER>
class ElasticInclinedPlaneProblem : 
 public InclinedPlaneProblem<ELEMENT>
{
public:
 //Constructor
 ElasticInclinedPlaneProblem(const unsigned &nx, const unsigned &ny,
                      const double &length) :
  InclinedPlaneProblem<ELEMENT>(nx,ny,length) 
  {
   //Set the name
   this->Output_prefix = "elastic";
   
   //Create our one and only timestepper, with adaptive timestepping
   add_time_stepper_pt(new TIMESTEPPER);

   Problem::mesh_pt() = new  ElasticInclinedPlaneMesh<ELEMENT>(
    nx,ny,length,1.0,this->time_stepper_pt());

   //Set the consititutive law for the elements
   unsigned long n_element = nx*ny;
   //Loop over all the fluid elements 
   for(unsigned long e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(e));
     temp_pt->constitutive_law_pt() = 
      Global_Physical_Variables::Constitutive_law_pt;

     // MH: No longer needed
     //temp_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;
    }

   //Complete the rest of the build
   this->complete_build(nx,ny);
  } //end of constructor

 ///\short Update Lagrangian positions after each timestep 
 ///(updated-lagrangian approach)
 void actions_after_implicit_timestep()
  {
   //Now loop over all the nodes and reset their Lagrangian coordinates
   unsigned n_node = this->mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     //Cast node to an elastic node
     SolidNode* temp_pt = static_cast<SolidNode*>(this->mesh_pt()->node_pt(n));
     for(unsigned j=0;j<2;j++) {temp_pt->xi(j) = temp_pt->x(j);}
    }
  } //end of actions_after_implicit_timestep

 //Destructor, destory mesh and timestepper
 ~ElasticInclinedPlaneProblem()
  {
   delete this->mesh_pt();
   delete this->time_stepper_pt();
  }
};


//start of main
int main(int argc, char **argv)
{

 using namespace Global_Physical_Variables;
 
 //Set the constitutive law for the mesh deformation
 Constitutive_law_pt  = new GeneralisedHookean(&Global_Physical_Variables::Nu);


#ifdef CR_ELEMENT
#define FLUID_ELEMENT QCrouzeixRaviartElement<2>
#else
#define FLUID_ELEMENT QTaylorHoodElement<2>
#endif 

 //Initialise physical parameters
 Re = 4.0/sin(Alpha);
 G[0] = sin(Alpha);
 G[1] = -cos(Alpha);
 
 Wall_normal.resize(2);
 Wall_normal[0] = 1.0;
 Wall_normal[1] = 0.0;

 //Spine problem
 {
  //Create the problem
  SpineInclinedPlaneProblem<SpineElement<FLUID_ELEMENT >, 
   BDF<2> > 
   problem(10,4,Global_Physical_Variables::Length);
  
  //Solve the steady problem
  problem.solve_steady();
  
  //Prepare the problem for timestepping 
  //(assume that it's been at the flat-film solution for all previous time)
  double dt = 0.1;
  problem.assign_initial_values_impulsive(dt);
  
  //Apply a small sinusoidal perturbation to the flat film
 {
  const unsigned n_node = problem.mesh_pt()->nboundary_node(2);
  for(unsigned n=0;n<n_node;n++)
   {
    Node* nod_pt = problem.mesh_pt()->boundary_node_pt(2,n);
    double x = nod_pt->x(0);
    double arg = Global_Physical_Variables::K*x; 
    static_cast<SpineNode*>(problem.mesh_pt()->boundary_node_pt(2,n))->
     spine_pt()->height() 
     += 0.01*sin(arg);
   }

  //Remember to update the nodal positions
  problem.mesh_pt()->node_update();
 } 
 
 //Timestep it
 problem.timestep(dt,2);
 } //End of spine problem

 //Elastic problem
 {
  //Create the problem
  ElasticInclinedPlaneProblem<
   PseudoSolidNodeUpdateElement<FLUID_ELEMENT,
   QPVDElement<2,3> >, BDF<2> > 
     problem(10,4,Global_Physical_Variables::Length);

  //Solve the steady problem
  problem.solve_steady();
  
  //Prepare the problem for timestepping 
  //(assume that it's been at the flat-film solution for all previous time)
  double dt = 0.1;
  problem.assign_initial_values_impulsive(dt);
  
  //Apply a small sinusoidal perturbation to the flat film
  {
   const unsigned n_node = problem.mesh_pt()->nboundary_node(2);
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = problem.mesh_pt()->boundary_node_pt(2,n);
     double x = nod_pt->x(0);
     double arg = Global_Physical_Variables::K*x; 
     problem.mesh_pt()->boundary_node_pt(2,n)->x(1) += 0.01*sin(arg);
     }
  }
  
  //Timestep it
  problem.timestep(dt,2);
 } //End of elastic problem
}








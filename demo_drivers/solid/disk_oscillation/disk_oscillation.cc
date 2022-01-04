//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
//Driver function for a simple test elasticity problem

// oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "oomph_crbond_bessel.h"

//Need to instantiate templated mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{
 /// Poisson's ratio
 double Nu=0.3;

 /// Timescale ratio 
 double Lambda_sq=(1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));
 
 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Multiplier for inertia terms (needed for consistent assignment 
 /// of initial conditions in Newmark scheme)
 double multiplier(const Vector<double>& xi)
 {
  return Global_Physical_Variables::Lambda_sq;
 }


} // end namespace


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
// Axisymmetrically oscillating disk
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//================disk_as_geom_object======================================
/// Axisymmetrially oscillating disk with displacement
/// field according to linear elasticity.
//=========================================================================
class AxisymOscillatingDisk : public GeomObject
{

public:

 /// Constructor: 2 Lagrangian coordinate, 2 Eulerian coords. Pass 
 /// amplitude of oscillation, Poisson ratio nu, and pointer to 
 /// global timestepper.
 AxisymOscillatingDisk(const double& ampl, const double& nu, 
                       TimeStepper* time_stepper_pt);

 /// Destructor (empty)
 ~AxisymOscillatingDisk(){}

 /// Position vector at Lagrangian coordinate xi at present
 /// time
 void position(const Vector<double>& xi, Vector<double>& r) const;

 /// Parametrised velocity on object at current time: veloc = d r(xi)/dt.
 void veloc(const Vector<double>& xi, Vector<double>& veloc);

 /// Parametrised acceleration on object at current time: 
 /// accel = d^2 r(xi)/dt^2.
 void accel(const Vector<double>& xi, Vector<double>& accel);
 
 /// Parametrised j-th time-derivative on object at current time: 
 /// \f$ \frac{d^{j} r(\zeta)}{dt^j} \f$.
 void dposition_dt(const Vector<double>& xi, const unsigned& j, 
                  Vector<double>& drdt)
  {
   switch (j)
    {
     // Current position
    case 0:
     position(xi,drdt);
     break;
     
     // Velocity:
    case 1:
     veloc(xi,drdt);
     break;

     // Acceleration:
    case 2:
     accel(xi,drdt);
     break;
     
    default:
     std::ostringstream error_message;
     error_message << j << "th derivative not implemented\n";
     
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }

   
 /// Residual of dispersion relation for use in black-box Newton method
 /// which requires global (or static) functions. 
 /// Poisson's ratio is  passed as parameter.
 static void residual_for_dispersion(const Vector<double>& param, 
                                     const Vector<double>& omega,
                                     Vector<double>& residual);

private:

 /// Amplitude of oscillation
 double Ampl;

 /// Period of oscillation
 double T;

 /// Poisson ratio nu
 double Nu;

 /// Eigenfrequency
 double Omega;

}; // end disk_as_geom_object




//==============ic_constructor============================================
/// Constructor: 2 Lagrangian coordinates, 2 Eulerian coords. Pass 
/// amplitude of oscillation, Poisson ratio nu, and pointer to 
/// global timestepper. 
//========================================================================
AxisymOscillatingDisk::AxisymOscillatingDisk(const double& ampl,
                                             const double& nu, 
                                             TimeStepper* time_stepper_pt) : 
 GeomObject(2,2,time_stepper_pt), Ampl(ampl), Nu(nu)
{
 // Parameters for dispersion relation
 Vector<double> param(1);
 param[0]=Nu;
 
 // Initial guess for eigenfrequency
 Vector<double> omega(1);
 omega[0]=2.0; 
 
 // Find eigenfrequency from black box Newton solver
 BlackBoxFDNewtonSolver::black_box_fd_newton_solve(residual_for_dispersion,
                                                   param,omega);

 // Assign eigenfrequency
 Omega=omega[0];
 
 // Assign/doc period of oscillation
 T=2.0*MathematicalConstants::Pi/Omega;

 std::cout << "Period of oscillation: " << T << std::endl;

}



//===============start_position===========================================
/// Position Vector at Lagrangian coordinate xi at present
/// time
//========================================================================
void AxisymOscillatingDisk::position(const Vector<double>& xi,
                                     Vector<double>& r) const
{
 // Parameter values at present time
 double time=Geom_object_time_stepper_pt->time_pt()->time();

 // Radius in Lagrangian coordinates
 double lagr_radius=sqrt( pow(xi[0],2) + pow(xi[1],2) );
 
 if (lagr_radius<1.0e-12)
  {
   // Position Vector
   r[0]=0.0;
   r[1]=0.0;
  }
 else
  {
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
   double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
   CRBond_Bessel::bessjy01a(Omega*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Displacement field 
   double u=Ampl*j1*sin(2.0*MathematicalConstants::Pi*time/T);
   
   // Position Vector
   r[0]=(xi[0]+xi[0]/lagr_radius*u);
   r[1]=(xi[1]+xi[1]/lagr_radius*u);
  }

} //end position




//========================================================================
/// Parametrised velocity on object at current time: 
/// veloc = d r(xi)/dt.
//========================================================================
void AxisymOscillatingDisk::veloc(const Vector<double>& xi,
                                  Vector<double>& veloc)
{
 // Parameter values at present time
 double time=Geom_object_time_stepper_pt->time_pt()->time();
 
 // Radius in Lagrangian coordinates
 double lagr_radius=sqrt(pow(xi[0],2)+pow(xi[1],2));
 
 if (lagr_radius<1.0e-12)
  {
   // Veloc vector
   veloc[0]=0.0;
   veloc[1]=0.0;
  }
 else
  {
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
   double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
   CRBond_Bessel::bessjy01a(Omega*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Deriv of displacement field 
   double udot=2.0*MathematicalConstants::Pi/T*Ampl*j1*
    cos(2.0*MathematicalConstants::Pi*time/T);
   
   // Veloc
   veloc[0]=(xi[0]/lagr_radius*udot);
   veloc[1]=(xi[1]/lagr_radius*udot);
  }
}





//========================================================================
/// Parametrised acceleration on object at current time: 
/// accel = d^2 r(xi)/dt^2.
//========================================================================
void AxisymOscillatingDisk::accel(const Vector<double>& xi,
                                  Vector<double>& accel) 
{
 // Parameter values at present time
 double time=Geom_object_time_stepper_pt->time_pt()->time();
 
 // Radius in Lagrangian coordinates
 double lagr_radius=sqrt(pow(xi[0],2)+pow(xi[1],2));
 
 
 if (lagr_radius<1.0e-12)
  {
   // Veloc vector
   accel[0]=0.0;
   accel[1]=0.0;
  }
 else
  {
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
   double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
   CRBond_Bessel::bessjy01a(Omega*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Deriv of displacement field
   double udotdot=-pow(2.0*MathematicalConstants::Pi/T,2)*Ampl*j1*
    sin(2.0*MathematicalConstants::Pi*time/T);
   
   // Veloc
   accel[0]=(xi[0]/lagr_radius*udotdot);
   accel[1]=(xi[1]/lagr_radius*udotdot);

  }
}
 
  

//======================start_of_dispersion===============================
/// Residual of dispersion relation for use in black box Newton method
/// which requires global (or static) functions. 
/// Poisson's ratio is passed as parameter.
//========================================================================
void AxisymOscillatingDisk::residual_for_dispersion(
 const Vector<double>& param, const Vector<double>& omega,
 Vector<double>& residual)
{
 // Extract parameters
 double nu=param[0];
 
 // Argument of various Bessel functions
 double arg=omega[0];
  
 // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x) and their derivatives
 double j0,j1,y0,y1,j0p,j1p,y0p,y1p;
 CRBond_Bessel::bessjy01a(arg,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
 
 // Residual of dispersion relation
 residual[0]=nu/(1.0-2.0*nu)*(j1+(j0-j1/omega[0])*omega[0])+
  (j0-j1/omega[0])*omega[0];
 
}



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//======================start_mesh================================
/// Elastic quarter circle sector mesh: We "upgrade"
/// the RefineableQuarterCircleSectorMesh to become an
/// SolidMesh and equate the Eulerian and Lagrangian coordinates,
/// thus making the domain represented by the mesh the stress-free 
/// configuration. 
//================================================================
template <class ELEMENT>
class ElasticRefineableQuarterCircleSectorMesh :
 public virtual RefineableQuarterCircleSectorMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
 /// ones so that the initial configuration is the stress-free one.
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(GeomObject* wall_pt,
                                         const double& xi_lo,
                                         const double& fract_mid,
                                         const double& xi_hi,
                                         TimeStepper* time_stepper_pt=
                                         &Mesh::Default_TimeStepper) :
  RefineableQuarterCircleSectorMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                             time_stepper_pt)
  {
#ifdef PARANOID
   /// Check that the element type is derived from the SolidFiniteElement
   SolidFiniteElement* el_pt=dynamic_cast<SolidFiniteElement*>
    (finite_element_pt(0));
   if (el_pt==0)
    {
     throw OomphLibError(
      "Element needs to be derived from SolidFiniteElement\n",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Make the current configuration the undeformed one by
   // setting the nodal Lagrangian coordinates to their current
   // Eulerian ones
   set_lagrangian_nodal_coordinates();
  }

};





/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////



//========start_class=================================================== 
/// Problem class to simulate small-amplitude oscillations of 
/// a circular disk.
//====================================================================== 
template<class ELEMENT>
class DiskOscillationProblem : public Problem
{

public:

 /// Constructor
 DiskOscillationProblem();

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Access function for the solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<ElasticRefineableQuarterCircleSectorMesh<ELEMENT>*>(
    Problem::mesh_pt());
  } 

 /// Run the problem: Pass number of timesteps to be performed.
 void run(const unsigned& nstep);
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Trace file
 ofstream Trace_file;
 
 /// Vector of pointers to nodes whose position we're tracing
 Vector<Node*> Trace_node_pt;

 /// Geometric object that specifies the initial conditions
 AxisymOscillatingDisk* IC_geom_object_pt;

}; // end class





//============start_constructor========================================= 
/// Constructor
//====================================================================== 
template<class ELEMENT>
DiskOscillationProblem<ELEMENT>::DiskOscillationProblem() 
{

 // Allocate the timestepper: The classical Newmark scheme with
 // two history values.
 add_time_stepper_pt(new Newmark<2>);
 
 

 // GeomObject that specifies the curvilinear boundary of the
 // circular disk
 GeomObject* curved_boundary_pt=new Ellipse(1.0,1.0);

 //The start and end intrinsic coordinates on the geometric object
 // that defines the curvilinear boundary of the disk
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 // Fraction along geometric object at which the radial dividing line
 // is placed
 double fract_mid=0.5;

 //Now create the mesh
 Problem::mesh_pt()= new ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(
  curved_boundary_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());


 // Setup trace nodes as the nodes on boundaries 0 (= horizontal symmetry
 // boundary) and 1 (=curved boundary)
 unsigned nnod0=mesh_pt()->nboundary_node(0);
 unsigned nnod1=mesh_pt()->nboundary_node(1);
 Trace_node_pt.resize(nnod0+nnod1);
 for (unsigned j=0;j<nnod0;j++)
  {

   Trace_node_pt[j]=mesh_pt()->boundary_node_pt(0,j);

  }
 for (unsigned j=0;j<nnod1;j++)
  {

   Trace_node_pt[j+nnod0]=mesh_pt()->boundary_node_pt(1,j);

  } //done choosing trace nodes


 // Pin the horizontal boundary in the vertical direction
 unsigned n_hor = mesh_pt()->nboundary_node(0);
 for(unsigned i=0;i<n_hor;i++)
  {
   mesh_pt()->boundary_node_pt(0,i)->pin_position(1);
  }

 // Pin the vertical boundary in the horizontal direction
 unsigned n_vert = mesh_pt()->nboundary_node(2);
 for(unsigned i=0;i<n_vert;i++)
  {

   mesh_pt()->boundary_node_pt(2,i)->pin_position(0);

  } // done bcs


 //Finish build of elements
 unsigned n_element =mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   // Set the timescale ratio
   el_pt->lambda_sq_pt()=&Global_Physical_Variables::Lambda_sq;
  }

 // Refine uniformly
 mesh_pt()->refine_uniformly();

 // Assign equation numbers
 assign_eqn_numbers();

} // end constructor



//=====================start_doc====================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void DiskOscillationProblem<ELEMENT>::doc_solution(
 DocInfo& doc_info)
{

 ofstream some_file, some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output shape of deformed body
 //------------------------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 

 // Write trace file
 //-----------------

 // Get position on IC object (exact solution)
 Vector<double> r_exact(2);
 Vector<double> xi(2);
 xi[0]=1.0;
 xi[1]=0.0;
 IC_geom_object_pt->position(xi,r_exact);
 
 // Exact outer radius for linear elasticity
 double exact_r=r_exact[0]; 
 
 // Add to trace file
 Trace_file << time_pt()->time()  << " " 
            << exact_r << " ";
 
 // Doc radii of control nodes
 unsigned ntrace_node=Trace_node_pt.size();
 for (unsigned j=0;j<ntrace_node;j++)
  {
   Trace_file << sqrt(pow(Trace_node_pt[j]->x(0),2)+
                      pow(Trace_node_pt[j]->x(1),2)) << " ";
  }
 Trace_file << std::endl;


 // Get displacement as a function of the radial coordinate
 //--------------------------------------------------------
 // along boundary 0
 //-----------------
 {
  // Number of elements along boundary 0:
  unsigned nelem=mesh_pt()->nboundary_element(0);

  // Open files
  sprintf(filename,"%s/displ_along_line%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);

  ofstream some_file2;
  sprintf(filename,"%s/exact_displ_along_line%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file2.open(filename);
  
  Vector<double> s(2);
  Vector<double> x(2);
  Vector<double> dxdt(2);
  Vector<double> xi(2);
  Vector<double> r_exact(2);
  Vector<double> v_exact(2);

  for (unsigned e=0;e<nelem;e++)
   {
    some_file << "ZONE " << std::endl;
    some_file2 << "ZONE " << std::endl;

    for (unsigned i=0;i<npts;i++)
     {
      // Move along bottom edge of element
      s[0]=-1.0+2.0*double(i)/double(npts-1);
      s[1]=-1.0;

      // Get pointer to element
      SolidFiniteElement* el_pt=dynamic_cast<SolidFiniteElement*>
       (mesh_pt()->boundary_element_pt(0,e));
      
      // Get Lagrangian coordinate
      el_pt->interpolated_xi(s,xi);

      // Get Eulerian coordinate
      el_pt->interpolated_x(s,x);

      // Get velocity 
      el_pt->interpolated_dxdt(s,1,dxdt);

      // Get exact Eulerian position
      IC_geom_object_pt->position(xi,r_exact);

      // Get exact velocity
      IC_geom_object_pt->veloc(xi,v_exact);
  
      // Plot radial distance and displacement
      some_file << xi[0] << " " << x[0]-xi[0] << " " 
                << dxdt[0] << std::endl;

      some_file2 << xi[0] << " " << r_exact[0]-xi[0] << " " 
                 << v_exact[0] << std::endl;

     }
   }
  some_file.close(); 
  some_file2.close();
 
 } // end line output


 // Get displacement (exact and computed) throughout domain
 //--------------------------------------------------------
 {
  // Number of elements
  unsigned nelem=mesh_pt()->nelement();

  // Open files
  sprintf(filename,"%s/displ%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);

  sprintf(filename,"%s/exact_displ%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file2.open(filename);
  
  Vector<double> s(2);
  Vector<double> x(2);
  Vector<double> dxdt(2);
  Vector<double> xi(2);
  Vector<double> r_exact(2);
  Vector<double> v_exact(2);

  for (unsigned e=0;e<nelem;e++)
   {
    some_file << "ZONE I=" << npts << ", J="<< npts <<  std::endl;
    some_file2  << "ZONE I=" << npts << ", J="<< npts <<  std::endl;

    for (unsigned i=0;i<npts;i++)
     {
      s[0]=-1.0+2.0*double(i)/double(npts-1);
      for (unsigned j=0;j<npts;j++)
       {
        s[1]=-1.0+2.0*double(j)/double(npts-1);

        // Get pointer to element
        SolidFiniteElement* el_pt=dynamic_cast<SolidFiniteElement*>
         (mesh_pt()->element_pt(e));
        
        // Get Lagrangian coordinate
        el_pt->interpolated_xi(s,xi);
        
        // Get Eulerian coordinate
        el_pt->interpolated_x(s,x);
        
        // Get velocity 
        el_pt->interpolated_dxdt(s,1,dxdt);
        
        // Get exact Eulerian position
        IC_geom_object_pt->position(xi,r_exact);
        
        // Get exact velocity
        IC_geom_object_pt->veloc(xi,v_exact);
        
        // Plot Lagrangian position, displacement and velocity
        some_file << xi[0] << " " << xi[1] << " " 
                  << x[0]-xi[0] << " " << x[1]-xi[1] << " " 
                  << dxdt[0] << " " <<  dxdt[1] << std::endl;
        
        
        some_file2 << xi[0] << " " << xi[1] << " " 
                   << r_exact[0]-xi[0] << " " << r_exact[1]-xi[1] << " " 
                   << v_exact[0]  << " " << v_exact[1] << std::endl;
        
       }
     }
   }
  some_file.close(); 
  some_file2.close(); 
 }
 
 
}




//=====================start_run====================================
/// Run the problem: Pass number of timesteps to be performed.
//==================================================================
template<class ELEMENT>
void DiskOscillationProblem<ELEMENT>::run(const unsigned& nstep)
{
 
 // Output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise time
 double time0=1.0;
 time_pt()->time()=time0;

 // Set timestep
 double dt=0.01; 
 time_pt()->initialise_dt(dt);

  // Create geometric object that specifies the initial conditions:
 
 // Amplitude of the oscillation
 double ampl=0.005;

 // Build the GeomObject
 IC_geom_object_pt=new AxisymOscillatingDisk(ampl,
                                             Global_Physical_Variables::Nu,
                                             time_stepper_pt()); 
 
 // Turn into object that specifies the initial conditions:
 SolidInitialCondition* IC_pt = new SolidInitialCondition(IC_geom_object_pt);

 // Assign initial condition 
 SolidMesh::Solid_IC_problem.set_newmark_initial_condition_consistently(
  this,mesh_pt(),time_stepper_pt(),IC_pt,dt,
  Global_Physical_Variables::multiplier); 

 // Doc initial state
 doc_solution(doc_info);
 doc_info.number()++;

 //Timestepping loop
 for(unsigned i=0;i<nstep;i++)
  {
   unsteady_newton_solve(dt);
   doc_solution(doc_info);
   doc_info.number()++;
  }

} // end of run






//========start_main====================================================
/// Driver for disk oscillation problem
//======================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // If there's a command line argument run the validation (i.e. do only 
 // 10 timesteps); otherwise do a few cycles
 unsigned nstep=1000;
 if (CommandLineArgs::Argc!=1)
  {
   nstep=10;
  }

 // Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 //Set up the problem
 DiskOscillationProblem<RefineableQPVDElement<2,3> > problem;
 
 //Run the simulation
 problem.run(nstep);
 
} // end of main









//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include <iostream>
#include <fstream>
#include <cmath>

//My own includes
#include "generic.h"
#include "solid.h"
#include "oomph_crbond_bessel.h"

//Need to instantiate templated mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//================================================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Pointer to strain energy function
 StrainEnergyFunction*Strain_energy_function_pt;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C1=1.3;

 /// Uniform pressure
 double P = 0.00;

 /// Uniform volumetric expansion
 double Uniform_gamma=1.1;

 /// Timescale ratio 
 double Lambda_sq=1.0;

 /// Timescale ration in analytical solution
 double K=1.0;

 /// Constant pressure load
 void constant_pressure(const Vector<double> &xi,const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 }

 /// Growth function: gamma=sqrt(g)=sqrt(g_ij)=const
 void growth_function(const Vector<double>& xi, double& gamma)
 {
  gamma=Uniform_gamma;
 }


 /// \short Multiplier for inertia terms (needed for consistent assignment 
 /// of initial conditions in Newmark scheme)
 double multiplier(const Vector<double>& xi)
 {
  return Global_Physical_Variables::Lambda_sq*
   Global_Physical_Variables::Uniform_gamma;
 }

}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Axisymmetrically oscillating disk
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=========================================================================
/// \short Axisymmetrially oscillating disk:
/// Disk has undeformed radius \f$ a \f$ hierher
//=========================================================================
class AxisymOscillatingDisk : public GeomObject
{

public:

 /// \short Constructor: 2 Lagrangian coordinate, 2 Eulerian coords. Pass 
 /// undeformed radius, amplitude of oscillation, timescale ratio K, 
 /// Poisson ratio nu, and pointer to global timestepper. No geometric 
 /// data. hierher K or K_sq?
 AxisymOscillatingDisk(const double& undef_radius, const double& ampl,
                       const double& k, const double& nu, 
                       TimeStepper* time_stepper_pt);

 /// Destructor (empty)
 ~AxisymOscillatingDisk(){}

 /// Access function for amplitude
 double ampl(){return Ampl;}

 /// Access function for period of oscillation
 double period(){return T;}


 /// \short Residual of dispersion relation for use in black box Newton method
 /// which requires global (or static) functions. Timescale ratio
 /// and Poisson's ratio are passed as parameters (0 and 1, respectively)
 /// since we can't access the member data of any specific object.
 static void residual_for_dispersion(const Vector<double>& param, 
                                     const Vector<double>& omega,
                                     Vector<double>& residual);

 /// \short Position vector at Lagrangian coordinate xi at present
 /// time
 void position(const Vector<double>& xi, Vector<double>& r) const;


 ///\short Parametrised velocity on object at current time: veloc = d r(xi)/dt.
 void veloc(const Vector<double>& xi, Vector<double>& veloc);


 /// \short Parametrised acceleration on object at current time: 
 /// accel = d^2 r(xi)/dt^2.
 void accel(const Vector<double>& xi, Vector<double>& accel);
   

 /// \short Position Vector at Lagrangian coordinate xi at discrete
 /// previous time (t=0: present time; t>0: previous time)   
 void position(const unsigned& t, const Vector<double>& xi, 
               Vector<double>& r) const
  {
   throw OomphLibError("Not implemented",
                       "AxisymmOscillatingDisk::position()",
                       OOMPH_EXCEPTION_LOCATION);
  }


 /// \short Derivative of position Vector w.r.t. to coordinates:
 /// dR_i/dxi_alpha = drdxi[alpha][i]. Evaluated at current time.
 void position(const Vector<double>& xi,
                             Vector<Vector<double> >& drdxi) const
  {
   throw OomphLibError("Not implemented",
                       "AxisymmOscillatingDisk::position()",
                       OOMPH_EXCEPTION_LOCATION);
  }


 /// \short 2nd derivative of position Vector w.r.t. to coordinates:
 /// d^2R_i/dxi_alpha dxi_beta = ddrdxi[alpha][beta][i].
 /// Evaluated at current time.
 void position(const Vector<double>& xi,
                             Vector<Vector<Vector<double> > >& ddrdxi) const
  {
   throw OomphLibError("Not implemented",
                       "AxisymmOscillatingDisk::position()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 /// \short Posn Vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:  dR_i/dxi_alpha = drdxi[alpha][i]
 /// d^2R_i/dxi_alpha dxi_beta = ddrdxi[alpha][beta][i].
 /// Evaluated at current time.
  void position(const Vector<double>& xi, Vector<double>& r,
                Vector<Vector<double> >& drdxi,
                Vector<Vector<Vector<double> > >& ddrdxi) const
  {
   throw OomphLibError("Not implemented",
                       "AxisymmOscillatingDisk::position()",
                       OOMPH_EXCEPTION_LOCATION);
  }


protected:

 /// Undeformed radius
 double A0;

 /// Amplitude of oscillation
 double Ampl;

 /// Period of oscillation
 double T;

 /// Timescale ratio K
 double K;

 /// Poisson ratio nu
 double Nu;

 /// Eigenfrequency
 double Omega;

};




//========================================================================
/// Constructor: 2 Lagrangian coordinates, 2 Eulerian coords. Pass 
/// undeformed radius, amplitude of oscillation, timescale ratio K, 
/// Poisson ratio nu, and pointer to global timestepper. 
/// No geometric data.
//========================================================================
AxisymOscillatingDisk::AxisymOscillatingDisk(const double& undef_radius,
                                             const double& ampl,
                                             const double& k,
                                             const double& nu, 
                                             TimeStepper* time_stepper_pt) : 
 GeomObject(2,2,time_stepper_pt), A0(undef_radius), Ampl(ampl), K(k), Nu(nu)
{
 
 using namespace MathematicalConstants;

 // Get zero of dispersion relation
 //--------------------------------
 
 // Parameters for dispersion relation
 Vector<double> param(2);
 param[0]=K;
 param[1]=Nu;
 
 // Initial guess for eigenfrequency
 Vector<double> omega(1);
 omega[0]=2.125748928;
 
 // Find eigenfrequency from black box Newton solver
 BlackBoxFDNewtonSolver::black_box_fd_newton_solve(residual_for_dispersion,
                                                   param,omega);

 // Assign eigenfrequency
 Omega=omega[0];
 
 // Assign period of oscillation
 T=2.0*Pi/Omega;
 
}



//========================================================================
/// \short Position Vector at Lagrangian coordinate xi at present
/// time
//========================================================================
void AxisymOscillatingDisk::position(const Vector<double>& xi,
                                     Vector<double>& r) const
{
 using namespace MathematicalConstants;

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
   // Get the Bessel functions (all of them) at Lagrangian radius:
   double j0;
   double j1;
   double y0;
   double y1;
   double j0p;
   double j1p;
   double y0p;
   double y1p;
   
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x)
   // and their derivatives
   CRBond_Bessel::bessjy01a(Omega*K*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Displacement field 
   double u=Ampl*j1*sin(2.0*Pi*time/T);
   
   // Position Vector
   r[0]=A0*(xi[0]+xi[0]/lagr_radius*u);
   r[1]=A0*(xi[1]+xi[1]/lagr_radius*u);
  }
}




//========================================================================
///\short Parametrised velocity on object at current time: 
/// veloc = d r(xi)/dt.
//========================================================================
void AxisymOscillatingDisk::veloc(const Vector<double>& xi,
                                  Vector<double>& veloc)
{

 using namespace MathematicalConstants;

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
   // Get the Bessel functions (all of them) at Lagrangian radius:
   double j0;
   double j1;
   double y0;
   double y1;
   double j0p;
   double j1p;
   double y0p;
   double y1p;
   
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x)
   // and their derivatives
   CRBond_Bessel::bessjy01a(Omega*K*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Deriv of displacement field 
   double udot=2.0*Pi/T*Ampl*j1*cos(2.0*Pi*time/T);
   
   // Veloc
   veloc[0]=A0*(xi[0]/lagr_radius*udot);
   veloc[1]=A0*(xi[1]/lagr_radius*udot);
  }
}





//========================================================================
/// Parametrised acceleration on object at current time: 
/// accel = d^2 r(xi)/dt^2.
//========================================================================
void AxisymOscillatingDisk::accel(const Vector<double>& xi,
                                  Vector<double>& accel) 
{

 using namespace MathematicalConstants;

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
   // Get the Bessel functions (all of them) at Lagrangian radius:
   double j0;
   double j1;
   double y0;
   double y1;
   double j0p;
   double j1p;
   double y0p;
   double y1p;
   
   // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x)
   // and their derivatives
   CRBond_Bessel::bessjy01a(Omega*K*lagr_radius,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
   
   // Deriv of displacement field
   double udotdot=-pow(2.0*Pi/T,2)*Ampl*j1*sin(2.0*Pi*time/T);
   
   // Veloc
   accel[0]=A0*(xi[0]/lagr_radius*udotdot);
   accel[1]=A0*(xi[1]/lagr_radius*udotdot);

  }
}
 
  

//========================================================================
/// Residual of dispersion relation for use in black box Newton method
/// which requires global (or static) functions. Timescale ratio
/// and Poisson's ratio are passed as parameters (0 and 1, respectively)
/// since we can't access the member data of any specific object.
//========================================================================
void AxisymOscillatingDisk::residual_for_dispersion(
 const Vector<double>& param, const Vector<double>& omega,
 Vector<double>& residual)
{

 // Extract parameters
 double k=param[0];
 double nu=param[1];
 
 // Argument of various Bessel functions
 double arg=omega[0]*k;
 
 // Get the Bessel functions (all of them)
 double j0;
 double j1;
 double y0;
 double y1;
 double j0p;
 double j1p;
 double y0p;
 double y1p;
 
 // Bessel fcts J_0(x), J_1(x), Y_0(x), Y_1(x)
 // and their derivatives
 CRBond_Bessel::bessjy01a(arg,j0,j1,y0,y1,j0p,j1p,y0p,y1p);
 
 // Residual of dispersion relation
 residual[0]=
  nu/(1.0-2.0*nu)*(j1+(j0-j1/omega[0]/k)*omega[0]*k)+
  (j0-j1/omega[0]/k)*omega[0]*k;
 
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//================================================================
/// \short Mesh of traction elements -- doesn't contain any functionality
/// other than that provided by SolidMesh from which it needs
/// to inherit to ensure that the "elastic equation numbers" are assigned.
/// \n \n
/// Here's a reminder how this works:
/// - \c Mesh::assign_global_eqn_numbers() assigns the global
///   equation numbers for all nodes by calling \c Node::assign_eqn_numbers().
/// - \c SolidMeshes contain SolidNodes which overload this
///   function with \c SolidNode::assign_eqn_numbers(). This function
///   numbers the nodal values \e and the positional dofs, i.e.
///   it also generates the global elastic equation numbers.
/// - Local equation numbers are usually generated by 
///   \c Mesh::assign_local_eqn_numbers() which loops over all 
///   elements and calls
///   \c GeneralisedElement::assign_local_eqn_numbers(). 
/// - SolidMeshes overload this function with
///   \c SolidMesh::assign_local_eqn_numbers() which also
///   assigns local equation numbers for the "elastic equations".
///   by calling \c SolidFiniteElement::assign_solid_local_eqn_numbers().
/// .
/// If we simply inherit from Mesh, the local equation numbers
/// for the positional variables never get assigned. 
//================================================================
class TractionElementMesh : public SolidMesh
{

public:

 /// Empty constructor
 TractionElementMesh(){};

};




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//================================================================
/// Elastic quarter circle sector mesh with functionality to
/// attach traction elements to the curved surface. We "upgrade"
/// the RefineableQuarterCircleSectorMesh to become an
/// SolidMesh and equate the Eulerian and Lagrangian coordinates,
/// thus making the domain represented by the mesh the stress-free 
/// configuration. 
/// \n\n
/// The member function \c make_traction_element_mesh() creates
/// a separate mesh of SolidTractionElements that are attached to the
/// mesh's curved boundary (boundary 1). 
//================================================================
template <class ELEMENT>
class ElasticRefineableQuarterCircleSectorMesh :
 public virtual RefineableQuarterCircleSectorMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// \short Constructor: Build mesh and copy Eulerian coords to Lagrangian
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
      "ElasticRefineableQuarterCircleSectorMesh::Constructor()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Make the current configuration the undeformed one by
   // setting the nodal Lagrangian coordinates to their current
   // Eulerian ones
   set_lagrangian_nodal_coordinates();
  }


 /// Function to create mesh made of traction elements
 void make_traction_element_mesh(TractionElementMesh*& traction_mesh_pt)
  {

   // Make new mesh
   traction_mesh_pt=new TractionElementMesh;

   // Loop over all elements on boundary 1:
   unsigned b=1;
   unsigned n_element = this->nboundary_element(b);
   for (unsigned e=0;e<n_element;e++)
    {
     // The element itself:
     FiniteElement* fe_pt = this->boundary_element_pt(b,e);
     
     // Find the index of the face of element e along boundary b
     int face_index = this->face_index_at_boundary(b,e);
     
     // Create new element
     traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
                                      (fe_pt,face_index));
    }
  }


};





/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//====================================================================== 
/// \short Simulate small-amplitude oscillations of a circular disk
/// in its first eigenmode. Disk is in a state of plane strain.
//====================================================================== 
template<class ELEMENT, class TIMESTEPPER>
class DiskOscillationProblem : public Problem
{

private:


 /// Geometric object that defines the boundary of the undeformed disk
 PseudoBucklingRing* Curved_boundary_pt;

 /// Trace file
 ofstream Trace_file;
 
 /// Vector of pointers to nodes whose position we're tracing
 Vector<Node*> Trace_node_pt;

 /// Pointer to solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointer to mesh of traction elements
 TractionElementMesh*  Traction_mesh_pt;

 /// Geometric object that specifies the initial conditions
 AxisymOscillatingDisk* IC_geom_object_pt;


public:

 /// Constructor:
 DiskOscillationProblem();

 /// \short Run the problem: Pass number of case to allow for labeling of
 /// output and flag to indicate if IC is to be assigned consistently.
 /// Also pass number of timesteps to be performed.
 void run(const unsigned& case_number, const bool& consistent_ic,
          const unsigned& nstep);
 
 /// Access function for the solid mesh
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>*& solid_mesh_pt() 
  {
   return Solid_mesh_pt;  
  } 


 /// Access function for the mesh of surface traction elements
 TractionElementMesh*& traction_mesh_pt() 
  {
   return Traction_mesh_pt;  
  } 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Complete the build of all elements so they are fully functional
 void finish_problem_setup();

 /// \short Dump problem-specific parameters values, then dump
 /// generic problem data.
 void dump_it(ofstream& dump_file);

 /// \short Read problem-specific parameter values, then recover
 /// generic problem data.
 void restart(ifstream& restart_file);

};







//====================================================================== 
/// Constructor
//====================================================================== 
template<class ELEMENT, class TIMESTEPPER>
DiskOscillationProblem<ELEMENT,TIMESTEPPER>::DiskOscillationProblem() 
{

 //Over-ride the default max Newton step
 Max_residuals = 1000.0;

 // Parameters for pseudo-buckling ring that defines the curved domain 
 // boundary  (=quarter circle)
 double eps_buckl=0.0;  
 double ampl_ratio=-0.5;  
 unsigned n_buckl=2;
 double r_0=1.0;
 double T=1.0;

 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER);
 
 // Build wall geometric element
 Curved_boundary_pt=new PseudoBucklingRing(eps_buckl,ampl_ratio,n_buckl,r_0,T,
                                           time_stepper_pt());


 // The curved boundary of the mesh is defined by the geometric object
 // What follows are the start and end coordinates on the geometric object:
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 // Fraction along geometric object at which the radial dividing line
 // is placed
 double fract_mid=0.5;

 //Now create the mesh
 solid_mesh_pt() = new ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(
  Curved_boundary_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());


 // Setup trace nodes as the nodes on boundary 1 (=curved boundary) in 
 // the original mesh (they exist under any refinement!) 
 unsigned nnod0=solid_mesh_pt()->nboundary_node(0);
 unsigned nnod1=solid_mesh_pt()->nboundary_node(1);
 Trace_node_pt.resize(nnod0+nnod1);
 for (unsigned j=0;j<nnod0;j++)
  {
   Trace_node_pt[j]=solid_mesh_pt()->boundary_node_pt(0,j);
  }
 for (unsigned j=0;j<nnod1;j++)
  {
   Trace_node_pt[j+nnod0]=solid_mesh_pt()->boundary_node_pt(1,j);
  }

 // Refine uniformly
 solid_mesh_pt()->refine_uniformly();

 // Mesh has been adapted: Need to setup boundary info again --
 // hierher: shouldn't this be automatic
 solid_mesh_pt()->setup_boundary_element_info();

 // Build traction element mesh
 solid_mesh_pt()->make_traction_element_mesh(traction_mesh_pt());
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Traction mesh is first sub-mesh
 add_sub_mesh(traction_mesh_pt());

 // Build combined "global" mesh
 build_global_mesh();

 // Pin the bottom in the vertical direction
 unsigned n_bottom = solid_mesh_pt()->nboundary_node(0);
 //Loop over the node
 for(unsigned i=0;i<n_bottom;i++)
  {
   solid_mesh_pt()->boundary_node_pt(0,i)->pin_position(1);
  }

 // Pin the left edge in the horizontal direction
 unsigned n_side = solid_mesh_pt()->nboundary_node(2);
 //Loop over the node
 for(unsigned i=0;i<n_side;i++)
  {
   solid_mesh_pt()->boundary_node_pt(2,i)->pin_position(0);
  }

 //Find number of elements in solid mesh
 unsigned n_element =solid_mesh_pt()->nelement();
  
 //Loop over the elements in the main mesh
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   // Set the timescale ratio
   el_pt->lambda_sq_pt()=
    &Global_Physical_Variables::Lambda_sq;

   // Switch on inertia
   el_pt->unsteady()=true;

   // Set the isotropic growth function pointer
   el_pt->isotropic_growth_fct_pt()=Global_Physical_Variables::growth_function;
  }


 //Find number of elements in traction mesh
 n_element=traction_mesh_pt()->nelement();
  
 //Loop over the elements in the traction element mesh
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<SolidTractionElement<ELEMENT>*>
    (traction_mesh_pt()->element_pt(i));

   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 
}



// //==================================================================
// /// Complete the build of all elements so they are fully functional:
// /// Solid elements need to be told about the consitutive law,
// /// the growth function, the timescale ratio, and inertia
// /// needs to be switched on. The traction 
// /// elements need to be told about the traction function. 
// //==================================================================
// template<class ELEMENT, class TIMESTEPPER>
// void DiskOscillationProblem<ELEMENT,TIMESTEPPER>::finish_problem_setup()
// {

//  //Find number of elements in solid mesh
//  unsigned Nelement =solid_mesh_pt()->nelement();
  
//  //Loop over the elements in the main mesh
//  for(unsigned i=0;i<Nelement;i++)
//   {
//    //Cast to a solid element
//    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

//    //Set the constitutive law
//    el_pt->constitutive_law_pt() =
//     Global_Physical_Variables::Constitutive_law_pt;
   
//    // Set the timescale ratio
//    el_pt->lambda_sq_pt()=
//     &Global_Physical_Variables::Lambda_sq;

//    // Switch on inertia
//    el_pt->unsteady()=true;

//    // Set the isotropic growth function pointer
//    el_pt->isotropic_growth_fct_pt()=Global_Physical_Variables::growth_function;
//   }


//  //Find number of elements in traction mesh
//  Nelement=traction_mesh_pt()->nelement();
  
//  //Loop over the elements in the traction element mesh
//  for(unsigned i=0;i<Nelement;i++)
//   {
//    //Cast to a solid traction element
//    SolidTractionElement<ELEMENT> *el_pt = 
//     dynamic_cast<SolidTractionElement<ELEMENT>*>
//     (traction_mesh_pt()->element_pt(i));

//    //Set the traction function
//    el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
//   }

// }


//==================================================================
/// Doc the solution
//==================================================================
template<class ELEMENT, class TIMESTEPPER>
void DiskOscillationProblem<ELEMENT,TIMESTEPPER>::doc_solution(
 DocInfo& doc_info)
{

 ofstream some_file;
 ofstream some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output shape of deformed body
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output traction
 unsigned nel=traction_mesh_pt()->nelement();
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Vector<double> unit_normal(2);
 Vector<double> traction(2);
 Vector<double> x_dummy(2);
 Vector<double> s_dummy(1);
 for (unsigned e=0;e<nel;e++)
  {
   some_file << "ZONE " << std::endl;
   for (unsigned i=0;i<npts;i++)
    {
     s_dummy[0]=-1.0+2.0*double(i)/double(npts-1);
     SolidTractionElement<ELEMENT>* el_pt=
      dynamic_cast<SolidTractionElement<ELEMENT>*>(
       traction_mesh_pt()->finite_element_pt(e));
     el_pt->outer_unit_normal(s_dummy,unit_normal);
     el_pt->traction(s_dummy,traction);
     el_pt->interpolated_x(s_dummy,x_dummy);
     some_file << x_dummy[0] << " " << x_dummy[1] << " " 
               << traction[0] << " " << traction[1] << " "  
               << std::endl;
    }
  }
 some_file.close(); 


 // Get displacement as a function of the radial coordinate
 // along boundary 0
 {

  // Number of elements along boundary 0:
  unsigned nelem=solid_mesh_pt()->nboundary_element(0);

  // Open files
  sprintf(filename,"%s/displ_along_line%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);

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
       (solid_mesh_pt()->boundary_element_pt(0,e));
      
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
 }





 // Get displacement field throughout domain
 {

  // Number of elements along boundary 0:
  unsigned nelem=solid_mesh_pt()->nelement();

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
         (solid_mesh_pt()->element_pt(e));
        
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
 
 
  /// Get position on IC object (exact solution)
 Vector<double> r_exact(2);
 Vector<double> xi(2);
 xi[0]=1.0;
 xi[1]=0.0;
 IC_geom_object_pt->position(xi,r_exact);
  
  // Exact outer radius for linear elasticity
  double exact_r=r_exact[0]; 
  
 // Write trace file
 Trace_file << time_pt()->time()  << " " 
            << Global_Physical_Variables::P  << " " 
            << exact_r << " ";
  
  
 unsigned ntrace_node=Trace_node_pt.size();
 for (unsigned j=0;j<ntrace_node;j++)
  {
   Trace_file << sqrt(pow(Trace_node_pt[j]->x(0),2)+
                      pow(Trace_node_pt[j]->x(1),2)) << " ";
  }
 Trace_file << std::endl;

//  // Write restart file
//  sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  dump_it(some_file);
//  some_file.close();
 





 // Output principal stress vectors at the centre of all elements
 ofstream pos_file;
 ofstream neg_file;
 sprintf(filename,"%s/pos_principal_stress%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 pos_file.open(filename);
 sprintf(filename,"%s/neg_principal_stress%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 neg_file.open(filename);

 // Write dummy data in both so there's at lest one zone in each
 pos_file << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << std::endl;
 neg_file << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << std::endl;


 Vector<double> s(2);
 Vector<double> x(2);
 s[0]=0.0;
 s[1]=0.0;
 unsigned n_solid_element=solid_mesh_pt()->nelement();
 for (unsigned e=0;e<n_solid_element;e++)
  {
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   // Get principal stress
   DenseMatrix<double> principal_stress_vector(2);
   Vector<double> principal_stress(2);
   el_pt->get_principal_stress(s,principal_stress_vector,principal_stress);

   // Get position of centre of element
   el_pt->interpolated_x(s,x);

   // compute vectors at 45 degree for nearly hydrostatic pressure state
   Vector<Vector<double> > rot(2);
   rot[0].resize(2);
   rot[1].resize(2);

   bool hydrostat=false;

   // Max. relative difference between principal stresses
   // required to classify stress state as non-hydrostatic: 1%
   double dev_max=1.0e-2;
   if (principal_stress[0]!=0.0)
    {
     if (abs((principal_stress[0]-principal_stress[1])/
             principal_stress[0])<dev_max)
      {
       hydrostat=true;
       double Cos=cos(0.25*3.14159);
       double Sin=sin(0.25*3.14159);
       rot[0][0]=
        Cos*principal_stress_vector(0,0)-Sin*principal_stress_vector(0,1);
       rot[0][1]=
        Sin*principal_stress_vector(0,0)+Cos*principal_stress_vector(0,1);
       rot[1][0]=
        Cos*principal_stress_vector(1,0)-Sin*principal_stress_vector(1,1);
       rot[1][1]=
        Sin*principal_stress_vector(1,0)+Cos*principal_stress_vector(1,1);
      }
    }

   // Loop over two principal stresses:
   for (unsigned i=0;i<2;i++)
    {
     if (principal_stress[i]>0.0)
      {
       pos_file << x[0] << " " << x[1] << " " 
                << principal_stress_vector(i,0) << " "
                << principal_stress_vector(i,1) << std::endl;
       pos_file << x[0] << " " << x[1] << " " 
                << -principal_stress_vector(i,0) << " "
                << -principal_stress_vector(i,1) << std::endl;
       if (hydrostat)
        {
         pos_file << x[0] << " " << x[1] << " " 
                  << rot[i][0] << " "
                  << rot[i][1] << std::endl;
         pos_file << x[0] << " " << x[1] << " " 
                  << -rot[i][0] << " "
                  << -rot[i][1] << std::endl;
        }
      }
     else
      {
       neg_file << x[0] << " " << x[1] << " " 
                 << principal_stress_vector(i,0) << " "
                 << principal_stress_vector(i,1) << std::endl;
       neg_file << x[0] << " " << x[1] << " " 
                << -principal_stress_vector(i,0) << " "
                << -principal_stress_vector(i,1) << std::endl;
       if (hydrostat)
        {
         neg_file << x[0] << " " << x[1] << " " 
                  << rot[i][0] << " "
                  << rot[i][1] << std::endl;
         neg_file << x[0] << " " << x[1] << " " 
                  << -rot[i][0] << " "
                  << -rot[i][1] << std::endl;
        }
      }
    }
  }
       
 pos_file.close();
 neg_file.close();


 cout << "Doced solution for step " 
      << doc_info.number() 
      << std::endl << std::endl << std::endl;
}



//========================================================================
/// Dump the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void DiskOscillationProblem<ELEMENT,TIMESTEPPER>::dump_it(ofstream& dump_file)
{
  
 // Call generic dump()
 Problem::dump(dump_file);

}



//========================================================================
/// Read solution from disk
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void DiskOscillationProblem<ELEMENT,TIMESTEPPER>::restart(ifstream& restart_file)
{

 // Read generic problem data
 Problem::read(restart_file);

}



//==================================================================
/// Run the problem: Pass number of case to allow for labeling of
/// output and flag to indicate if IC is to be assigned consistently.
/// Also pass number of timesteps to be performed.
//==================================================================
template<class ELEMENT, class TIMESTEPPER>
void DiskOscillationProblem<ELEMENT,TIMESTEPPER>::run(
 const unsigned& case_number, const bool& consistent_ic,
 const unsigned& nstep)
{
 
 // Output
 DocInfo doc_info;

 char dirname[100];
 sprintf(dirname,"RESLT%i",case_number);

 // Output directory
 doc_info.set_directory(dirname);

 // Step number
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

//  // Restart?
//  //---------

//  // Pointer to restart file
//  ifstream* restart_file_pt=0;

//  // No restart
//  //-----------
//  if (CommandLineArgs::Argc==1)
//   {
//    cout << "No restart" << std::endl;
//   }
//  // Restart
//  //--------
//  else if (CommandLineArgs::Argc==2)
//   {
//    // Open restart file
//    restart_file_pt=new ifstream(CommandLineArgs::Argv[1],ios_base::in);
//    if (restart_file_pt!=0)
//     {
//      cout << "Have opened " << CommandLineArgs::Argv[1] << 
//       " for restart. " << std::endl;
//     }
//    else
//     {
//      cout << "ERROR while trying to open " << CommandLineArgs::Argv[1] << 
//       " for restart." << std::endl;
//     }
//    // Do the actual restart
//    pause("need to do the actual restart");
//    //problem.restart(*restart_file_pt);
//   }
//  // More than one restart file specified?
//  else 
//   {
//    cout << "Can only specify one input file " << std::endl;
//    cout << "You specified the following command line arguments: " << std::endl;
//    CommandLineArgs::output();
//    assert(false);
//   }


 // Initial parameter values
 Global_Physical_Variables::P = 0.0; 

 // Initialise time
 double time0=1.0;
 time_pt()->time()=time0;

 // Set initial timestep
 double dt=0.01; 
 time_pt()->initialise_dt(dt);

 //Find out how many timesteppers there are
 unsigned Ntime_steppers = ntime_stepper();

 //Loop over them all and set the weights
 for(unsigned i=0;i<Ntime_steppers;i++)
  {
   time_stepper_pt(i)->set_weights();
  }


 // Geometric object that specifies the initial conditions:
 //--------------------------------------------------------
 
 // Undeformed radius (non-dim on the un-grown, unit radius)
 double a0=sqrt(Global_Physical_Variables::Uniform_gamma);

 // Amplitude of the oscillation
 double ampl=0.005;

 // Timescale ratio
 double K=Global_Physical_Variables::K;

 // Poisson's ratio
 double nu=Global_Physical_Variables::Nu;

 // Build the object
 IC_geom_object_pt=
  new AxisymOscillatingDisk(a0,ampl,K,nu,Problem::time_stepper_pt()); 


 // Create object that specifies the initial conditions:
 SolidInitialCondition* IC_pt = new SolidInitialCondition(IC_geom_object_pt);

 // IC consistent with equations or assigned directly?
 if (consistent_ic)
  {
   // Assign initial condition consistently
   SolidMesh::Solid_IC_problem.set_newmark_initial_condition_consistently(
    this,solid_mesh_pt(), static_cast<TIMESTEPPER*>(time_stepper_pt()),
    IC_pt,dt,&Global_Physical_Variables::multiplier);
  }
 else
  {
   // Assign position, veloc and accel directly
   SolidMesh::Solid_IC_problem.set_newmark_initial_condition_directly(
    this,solid_mesh_pt(), static_cast<TIMESTEPPER*>(time_stepper_pt()),
    IC_pt,dt);
  }

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

}






//======================================================================
/// Driver for simple elastic problem
//======================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);


 // If there's a command line argument run the validation (i.e. do only 
 // 10 timesteps; otherwise do a few cycles
 unsigned nstep=400;
 if (CommandLineArgs::Argc!=1)
  {
   nstep=10;
  }

 //Initialise physical parameters
 Global_Physical_Variables::E = 1.0; // ADJUST 
 Global_Physical_Variables::Nu = 0.3; // ADJUST
 Global_Physical_Variables::C1 = 1.3; // ADJUST

 // Assign timescale ratio for analytical solution hierher play
 Global_Physical_Variables::K = 0.7; // ADJUST

 // Timescale ratio for FE discretisation
 double nu=Global_Physical_Variables::Nu;
 Global_Physical_Variables::Lambda_sq=
  (1.0-nu)/((1.0+nu)*(1.0-2.0*nu))/Global_Physical_Variables::Uniform_gamma
  *pow(Global_Physical_Variables::K,2);


 // Define strain energy function: Generalised Mooney Rivlin
 Global_Physical_Variables::Strain_energy_function_pt = 
  new GeneralisedMooneyRivlin(Global_Physical_Variables::Nu,
                              Global_Physical_Variables::C1,
                              Global_Physical_Variables::E);



 // Case 0: Generalised Mooney Rivlin, consistent assignement of IC
 //----------------------------------------------------------------
 {

  // Define constitutive law (based on strain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
   new IsotropicStrainEnergyFunctionConstitutiveLaw(
    Global_Physical_Variables::Strain_energy_function_pt);
  
  //Set up the problem
  DiskOscillationProblem<RefineableQPVDElement<2,3>, Newmark<3> > problem;
  
  //Run the simulation
  problem.run(0,true,nstep);
  
  // Clean up 
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;

 }


 // Case 1: Generalised Mooney Rivlin, direct assignement of IC
 //----------------------------------------------------------------
 {

  // Define constitutive law (based on strain energy function)
  Global_Physical_Variables::Constitutive_law_pt = 
   new IsotropicStrainEnergyFunctionConstitutiveLaw(
    Global_Physical_Variables::Strain_energy_function_pt);
  
  //Set up the problem
  DiskOscillationProblem<RefineableQPVDElement<2,3>, Newmark<3> > problem;
  
  //Run the simulation
  problem.run(1,false,nstep);
  
  // Clean up 
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;

 }


 // Case 2: Generalised Hooke, consistent assignement of IC
 //--------------------------------------------------------
 {

  // "Big G" Linear constitutive equations:
  Global_Physical_Variables::Constitutive_law_pt = 
   new GeneralisedHookean(Global_Physical_Variables::Nu,
                          Global_Physical_Variables::E);
    
  //Set up the problem
  DiskOscillationProblem<RefineableQPVDElement<2,3>, Newmark<3> > problem;
  
  //Run the simulation
  problem.run(2,true,nstep);
  
  // Clean up 
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;

 }



 // Case 3: Generalised Hooke, direct assignement of IC
 //--------------------------------------------------------
 {

  // "Big G" Linear constitutive equations:
  Global_Physical_Variables::Constitutive_law_pt = 
   new GeneralisedHookean(Global_Physical_Variables::Nu,
                          Global_Physical_Variables::E);
    
  //Set up the problem
  DiskOscillationProblem<RefineableQPVDElement<2,3>, Newmark<3> > problem;
  
  //Run the simulation
  problem.run(3,false,nstep);
  
  // Clean up 
  delete Global_Physical_Variables::Constitutive_law_pt;
  Global_Physical_Variables::Constitutive_law_pt=0;

 }

}









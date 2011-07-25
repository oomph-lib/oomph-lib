//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Driver code for a 2D free-surface hydrostatics problem.
// The system consists of a layer of fluid 
// in a domain of height 1 and half-width 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// circular arcs and the pressure jump across the interface should be
// cos(alpha)/0.5 = 2 cos(alpha)/Ca.

//OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "constitutive.h"
#include "solid.h"

// The mesh
#include "meshes/single_layer_spine_mesh.h"
//Include our special fixed volume interface elements
#include "fix_vol_int_elements.h"
//Include our special fixed volumee elastic elements
#include "fix_vol_int_elastic_elements.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//====start_of_Global_Physical_Variables================
/// Namespace for phyical parameters
//======================================================
namespace Global_Physical_Variables
{

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 /// Pseudo-solid Mooney-Rivlin parameter
 double C1=1.0;

 /// Pseudo-solid Young's modulus
 double E=2.2;

 ///Direction of the wall normal vector
 Vector<double> Wall_normal;

 /// \short Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }
}

//============================================================================
/// Specific mesh for the static cap problem: Essentially a 
/// standard single-layer mesh with an additional element that applies
/// the volume constraint.
//============================================================================
template <class ELEMENT>
class StaticSingleLayerMesh : 
 public SingleLayerSpineMesh<SpineElement<ELEMENT>, 
 FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> > >
{

private:

 /// \short Pointer to the point element that is used to enforce
 /// conservation of mass
 FiniteElement* Point_element_pt;

public:

 /// Constructor: Pass number of elements in axial direction, number
 /// of elements in the layers and half-width.
 StaticSingleLayerMesh(const unsigned &nx, const unsigned &nh,
                       const double & half_width);
 
 /// Return pointer to the volumetric constraint element
 FiniteElement* &point_element_pt() {return Point_element_pt;}
 
};


//======================================================================
/// Constructor: Pass number of elements in horizontal direction, number
/// of elements in the two layers and half-width.
//======================================================================
template<class ELEMENT>
StaticSingleLayerMesh<ELEMENT>::StaticSingleLayerMesh(const unsigned &nx,
                                            const unsigned &nh,
                                            const double & half_width) :
 SingleLayerSpineMesh<SpineElement<ELEMENT>, 
  FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> > >
  (nx,nh,half_width,1.0)
                                               
{
 //Reorder the elements
this->element_reorder();

 // Last interface element:
 FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >* el_pt=
 dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
  SpineElement<ELEMENT> >*>(this->Interface_element_pt[this->Nx-1]);

 //Now make our edge (point)  element
 Point_element_pt  = el_pt->make_bounding_element(1);
 //Add it to the stack
 this->Element_pt.push_back(Point_element_pt);

}



//============================================================================
///A Problem class that solves the Navier--Stokes equations
///in an 2D geometry
//============================================================================
template<class ELEMENT>
class CapProblem : public Problem
{


public:

 //Constructor: Boolean flag indicates if volume constraint is
 //applied by hijacking internal or external pressure
 CapProblem(const bool& hijack_internal);

 /// Peform a parameter study: Solve problem for a range of contact angles
 /// Pass name of output directory as a string
 void parameter_study(const string& dir_name);

 /// Finish full specification of the elements
 void finish_problem_setup();

 /// Overload access function for the mesh
 StaticSingleLayerMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<StaticSingleLayerMesh<ELEMENT>*>(Problem::mesh_pt());}

 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {mesh_pt()->node_update();}

 /// No other actions required after solve step
 void actions_after_newton_solve() {}

 /// No other actions required before after solve step
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 

private:

 /// The Capillary number 
 double Ca;

 /// The volume of the fluid
 double Volume;

 /// The external pressure
 double Pext;

 /// The contact angle
 double Angle;

 /// The normal to the wall
 //Vector<double> Wall_normal;

 /// Trace file
 ofstream Trace_file;

 /// Data object whose single value stores the external pressure
 Data* External_pressure_data_pt;

};



//======================================================================
/// Constructor: Pass boolean flag to indicate if the volume
/// constraint is applied by hijacking an internal pressure
/// or the external pressure
//======================================================================
template<class ELEMENT>
CapProblem<ELEMENT>::CapProblem(const bool& hijack_internal) :
 Ca(2.1),  //Initialise value of Ca to some random value
 Volume(0.5),  //Initialise the value of the volume 
 Pext(1.23),  //Initialise the external pressure to some random value
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = 1.0;
 Global_Physical_Variables::Wall_normal[1] = 0.0;

 // Number of elements in the horizontal direction
 unsigned nx=4;

 // Number of elements in the vertical direction
 unsigned nh=4;

 // Halfwidth of domain
 double half_width=0.5;

 //Construct mesh
 Problem::mesh_pt() = new StaticSingleLayerMesh<ELEMENT>(nx,nh,half_width);

 //Overwrite the linear solver (frontal solver) 
 //Problem::linear_solver_pt() = new HSL_MA42;

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Pext);

 // Create a pointer to the (single value) Data item that
 // will contain the pressure value that we're
 // trading for the volume constraint
 Data* traded_pressure_data_pt;

 // Which pressure are we trading for the volume constraint: We 
 // can either hijack an internal pressure or use the external pressure.
 if (hijack_internal)
  {
   // The external pressure is pinned -- the external pressure
   // sets the pressure throughout the domain -- we do not have
   // the liberty to fix another pressure value!
   External_pressure_data_pt->pin(0);

   //Hijack one of the pressure values in the fluid and use it 
   //as the pressure whose value is determined by the volume constraint.
   //(Its value will affect the residual of that element but it will not
   //be determined by it, i.e. it's hijacked).
   traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
    mesh_pt()->bulk_element_pt(0))->hijack_internal_value(0,0);
  }
 else
  {
   // Regard the external pressure as an unknown and add
   // it to the problem's global data so it gets included
   // in the equation numbering. Note that, at the moment,
   // there's no equation that determines its value!
   add_global_data(External_pressure_data_pt);

   // Declare the external pressure to be the pressure determined
   // by the volume constraint, i.e. the pressure that's "traded":
   traded_pressure_data_pt = External_pressure_data_pt;

   // Since the external pressure is "traded" for the volume constraint,
   // it no longer sets the overall pressure, and we 
   // can add an arbitrary constant to all pressures. To make 
   // the solution unique, we pin a single pressure value in the bulk: 
   // We arbitrarily set the pressure dof 0 in element 0 to zero.
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
  }



 // Loop over the elements on the interface to pass pointer to Ca and
 //-----------------------------------------------------------------
 // the pointers to the Data items that contains the single
 //------------------------------------------------------
 // (pressure) value that is "traded" for the volume constraint 
 //------------------------------------------------------------
 // and the external pressure
 //--------------------------
 unsigned n_interface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<n_interface;e++)
  {
   //Cast to a 1D element
   FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*el_pt=
   dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
    (mesh_pt()->interface_element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;
 
   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(External_pressure_data_pt);

   //Pass the Data item that contains the single (pressure) value
   // that has been "traded" for the volume constraint to the
   // surface elements -- hacky! 
   el_pt->set_traded_pressure_data(traded_pressure_data_pt);
  }

 
 //Finally, pass the Data item that contains the single (pressure)
 //---------------------------------------------------------------
 // value that has been "traded" for the volume constraint
 //-------------------------------------------------------
 // to the volume constraint element.
 //----------------------------------
 {
  SpineVolumeConstraintPointElement<SpineElement<ELEMENT> >* el_pt =
   dynamic_cast<SpineVolumeConstraintPointElement<SpineElement<ELEMENT> >*>
   (mesh_pt()->point_element_pt());
  
  el_pt->volume_pt() = &Volume;
  el_pt->set_traded_pressure_data(traded_pressure_data_pt);
 }    
 


 //Set the boundary conditions
 //---------------------------

 //Pin the velocities on all boundaries apart from the free surface
 //(boundary 2) where all velocities are free, and apart from the symmetry
 //line (boundary 3) where only the horizontal velocity is pinned
 unsigned n_bound=mesh_pt()->nboundary();
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0);
       if (b!=3)
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(1);
        }
      }
    }
  }


 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
 //dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
 //(mesh_pt()->interface_element_pt(n_interface-1))->
 // set_contact_angle_right(&Angle)

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->set_contact_angle(&Angle);

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->ca_pt() = &Ca;

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->wall_unit_normal_fct_pt() 
  =  &Global_Physical_Variables::wall_unit_normal_fct;

 
 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
}









//======================================================================
/// Perform a parameter study. Pass name of output directory as 
/// a string
//======================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::parameter_study(const string& dir_name)
{

 // Create DocInfo object (allows checking if output directory exists)
 DocInfo doc_info;
 doc_info.set_directory(dir_name);
 doc_info.number()=0;


 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 Trace_file << "VARIABLES=\"<greek>a</greek><sub>prescribed</sub>\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>left</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>right</sub>\",";
 Trace_file << "\"p<sub>fluid</sub>-p<sub>ext</sub>\",";
 Trace_file << "\"<greek>D</greek>p<sub>exact</sub>\"";
 Trace_file << std::endl;

 //Gradually increase the contact angle
 for(unsigned i=0;i<6;i++)
  {
   //Solve the problem
   steady_newton_solve();

   //Output result
   doc_solution(doc_info);

   // Bump up counter
   doc_info.number()++;

   //Decrease the contact angle
   Angle -= 5.0*MathematicalConstants::Pi/180.0;
  }

} 




//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 //Output domain
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Number of interface elements
// unsigned ninterface=mesh_pt()->ninterface_element();

 // Number of spines
 unsigned nspine=mesh_pt()->nspine();

 // Doc
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << mesh_pt()->spine_pt(0)->height();
 Trace_file << " "  << mesh_pt()->spine_pt(nspine-1)->height();
 Trace_file << " "  
//             << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                mesh_pt()->interface_element_pt(0))->
//                actual_contact_angle_left()*
//                180.0/MathematicalConstants::Pi << " " ;
//  Trace_file << " "  
//             << dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<
//                SpineElement<ELEMENT> >*>(
//                mesh_pt()->interface_element_pt(ninterface-1))->
//                actual_contact_angle_right()*180.0/MathematicalConstants::Pi 
            << " ";
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(mesh_pt()->bulk_element_pt(0))->p_nst(0)-
               External_pressure_data_pt->value(0);
 Trace_file << " " << -2.0*cos(Angle)/Ca;
 Trace_file << std::endl;

}

 
//======================================================================
///\short We need to upgrade our standard mesh into an Elastic mesh so 
///that it can be used by SolidElements
//======================================================================
template <class ELEMENT, class INTERFACE_ELEMENT>
class ElasticSingleLayerMeshWithInterface :
 public RectangularQuadMesh<ELEMENT>, //If this is virtual then the constructor
             ///does not get called when it should
 public SolidMesh
{
 
public:

 /// \short Constructor: Build mesh and copy Eulerian coords to Lagrangian
 /// ones so that the initial configuration is the stress-free one.
 ElasticSingleLayerMeshWithInterface<ELEMENT,INTERFACE_ELEMENT>(
  const unsigned &nx, 
  const unsigned &nh, 
  const double &half_width) :
  RectangularQuadMesh<ELEMENT>(nx,nh,half_width,1.0)
  {
#ifdef PARANOID
   /// Check that the element type is derived from the SolidFiniteElement
   SolidFiniteElement* tmp_el_pt=dynamic_cast<SolidFiniteElement*>
    (finite_element_pt(0));
   if (tmp_el_pt==0)
    {
     throw OomphLibError(
      "Element needs to be derived from SolidFiniteElement\n",
      "ElasticRefineableQuarterCircleSectorMesh::constructor()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Make the current configuration the undeformed one by
   // setting the nodal Lagrangian coordinates to the current
   // Eulerian values
   set_lagrangian_nodal_coordinates();

   // Add all the standard elements to the bulk_element array
   unsigned n_element = this->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     Bulk_element_pt.push_back(this->finite_element_pt(e));
    }

   // Now add the interface elements
   
   //Loop over the horizontal elements
   for(unsigned i=0;i<nx;i++)
    {
     //Construct the 1D element on the face 2
     FiniteElement *interface_element_pt =
      new INTERFACE_ELEMENT(this->finite_element_pt(nx*(nh-1)+i),2);
     
     //Push it back on the stack
     Element_pt.push_back(interface_element_pt);
     
     //Push it back onto the stack of interface elements
     Interface_element_pt.push_back(interface_element_pt);
    }

   // Last interface element:
   INTERFACE_ELEMENT* el_pt=
    dynamic_cast<INTERFACE_ELEMENT*>(this->Interface_element_pt[nx-1]);

   //Now add the one additional volumetric constraint point
   Point_element_pt = el_pt->make_bounding_element(1);

   this->Element_pt.push_back(Point_element_pt);
  }
  
 /// Return pointer to the volumetric constraint element
 FiniteElement* &point_element_pt() {return Point_element_pt;} 

 /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}
 

 /// Number of elements on interface
 unsigned long ninterface_element() const {return Interface_element_pt.size();}
 
 ///Access functions for pointers to elements in bulk
 FiniteElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}

 ///Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}
 /// Vector of pointers to element in the fluid layer
 Vector <FiniteElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

 /// Pointer to the point element used to enforce conservation of mass
 FiniteElement* Point_element_pt;
};


//==============================================================
/// Wrap up the template classes
//==============================================================
template<class ELEMENT>
class ElasticSingleLayerMesh : public 
ElasticSingleLayerMeshWithInterface<
 ELEMENT, FixedVolumeElasticLineFluidInterfaceElement<ELEMENT> >
{
public:
 ElasticSingleLayerMesh<ELEMENT>(const unsigned &nx, 
                        const unsigned &nh, 
                        const double &half_width) : 
  ElasticSingleLayerMeshWithInterface<ELEMENT,
                                      FixedVolumeElasticLineFluidInterfaceElement<ELEMENT> >(nx,nh,half_width) {}

};

//============================================================================
///A Problem class that solves the Navier--Stokes equations
///in an 2D geometry
//============================================================================
template<class ELEMENT>
class ElasticCapProblem : public Problem
{


public:

 //Constructor: Boolean flag indicates if volume constraint is
 //applied by hijacking internal or external pressure
 ElasticCapProblem(const bool& hijack_internal);

 /// Peform a parameter study: Solve problem for a range of contact angles
 /// Pass name of output directory as a string
 void parameter_study(const string& dir_name);

 /// Overload access function for the mesh
 ElasticSingleLayerMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<ElasticSingleLayerMesh<ELEMENT>*>(Problem::mesh_pt());}

 /// No action required before Newton step
 void actions_before_newton_convergence_check() { }

 void actions_after_newton_step() {}

 /// No other actions required after solve step
 void actions_after_newton_solve() {}

 /// No other actions required before after solve step
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 

private:

 /// The Capillary number 
 double Ca;

 /// The volume of the fluid
 double Volume;

 /// The external pressure
 double Pext;

 /// The contact angle
 double Angle;

 /// The wall normal vector
 Vector<double> Wall_normal;

 //Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Trace file
 ofstream Trace_file;

 /// Data object whose single value stores the external pressure
 Data* External_pressure_data_pt;

};



//======================================================================
/// Constructor: Pass boolean flag to indicate if the volume
/// constraint is applied by hijacking an internal pressure
/// or the external pressure
//======================================================================
template<class ELEMENT>
ElasticCapProblem<ELEMENT>::ElasticCapProblem(const bool& hijack_internal) :
 Ca(2.1),  //Initialise value of Ca to some random value
 Volume(0.5),  //Initialise the value of the volume 
 Pext(1.23),  //Initialise the external pressure to some random value
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 Wall_normal.resize(2);
 Wall_normal[0] = 1.0; Wall_normal[1] = 0.0;

 // Number of elements in the horizontal direction
 unsigned nx=4;

 // Number of elements in the vertical direction
 unsigned nh=4;

 // Halfwidth of domain
 double half_width=0.5;

 //Construct mesh
 Problem::mesh_pt() = new ElasticSingleLayerMesh<ELEMENT>(nx,nh,half_width);

 //Overwrite the linear solver (frontal solver) 
 //Problem::linear_solver_pt() = new HSL_MA42;

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Pext);

 // Create a pointer to the (single value) Data item that
 // will contain the pressure value that we're
 // trading for the volume constraint
 Data* traded_pressure_data_pt;

 // Which pressure are we trading for the volume constraint: We 
 // can either hijack an internal pressure or use the external pressure.
 if (hijack_internal)
  {
   // The external pressure is pinned -- the external pressure
   // sets the pressure throughout the domain -- we do not have
   // the liberty to fix another pressure value!
   External_pressure_data_pt->pin(0);

   //Hijack one of the pressure values in the fluid and use it 
   //as the pressure whose value is determined by the volume constraint.
   //(Its value will affect the residual of that element but it will not
   //be determined by it, i.e. it's hijacked).
   traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
    mesh_pt()->bulk_element_pt(0))->hijack_internal_value(0,0);
  }
 else
  {
   // Regard the external pressure is an unknown and add
   // it to the problem's global data so it gets included
   // in the equation numbering. Note that, at the moment,
   // there's no equation that determines its value!
   add_global_data(External_pressure_data_pt);

   // Declare the external pressure to be the pressure determined
   // by the volume constraint, i.e. the pressure that's "traded":
   traded_pressure_data_pt = External_pressure_data_pt;

   // Since the external pressure is "traded" for the volume constraint,
   // it no longer sets the overall pressure, and we 
   // can add an arbitrary constant to all pressures. To make 
   // the solution unique, we pin a single pressure value in the bulk: 
   // We arbitrarily set the pressure dof 0 in element 0 to zero.
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
  }

 //Set the constituive law
 Constitutive_law_pt = 
  //obsolete: new DeformedMetricTensorLinearElasticConstitutiveLaw(1.0,0.001);
  //new GeneralisedHookean(0.001,1.0);
  //new IsotropicStrainEnergyFunctionConstitutiveLaw(new TongOne(1.0,1.0,0.1));
  new IsotropicStrainEnergyFunctionConstitutiveLaw(
   new GeneralisedMooneyRivlin(&Global_Physical_Variables::Nu,
                               &Global_Physical_Variables::C1,
                               &Global_Physical_Variables::E));
 
//Loop over the elements
 unsigned n_bulk = mesh_pt()->nbulk();
 for(unsigned e=0;e<n_bulk;e++)
  {
   ELEMENT* el_pt = 
   dynamic_cast<ELEMENT*>(mesh_pt()->bulk_element_pt(e));
   
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   // Get Jacobian by FD -- yes for now
   // hierher -- change this when Pseudo-solid elements have 
   // been updated to take availability of analytical solid 
   // Jacobian into account
   el_pt->evaluate_jacobian_by_fd()=true;
  
  }

 // Loop over the elements on the interface to pass pointer to Ca and
 //-----------------------------------------------------------------
 // the pointers to the Data items that contains the single
 //------------------------------------------------------
 // (pressure) value that is "traded" for the volume constraint 
 //------------------------------------------------------------
 // and the external pressure
 //--------------------------
 unsigned n_interface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<n_interface;e++)
  {
   //Cast to a 1D element
   FixedVolumeElasticLineFluidInterfaceElement<ELEMENT>*el_pt=
    dynamic_cast<FixedVolumeElasticLineFluidInterfaceElement<ELEMENT>*>
    (mesh_pt()->interface_element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;
 
   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(External_pressure_data_pt);

   //Pass the Data item that contains the single (pressure) value
   // that has been "traded" for the volume constraint to the
   // surface elements -- hacky! 
   el_pt->set_traded_pressure_data(traded_pressure_data_pt);
  }

 
 //Finally, pass the Data item that contains the single (pressure)
 //---------------------------------------------------------------
 // value that has been "traded" for the volume constraint
 //-------------------------------------------------------
 // to the volume constraint element.
 //----------------------------------
 {
  ElasticVolumeConstraintPointElement<ELEMENT>* el_pt =
   dynamic_cast<ElasticVolumeConstraintPointElement<ELEMENT>*>
   (mesh_pt()->point_element_pt());
  
  el_pt->volume_pt() = &Volume;
  el_pt->set_traded_pressure_data(traded_pressure_data_pt);
 }    
 


 //Set the boundary conditions
 //---------------------------

 //Pin the velocities on all boundaries apart from the free surface
 //(boundary 2) where all velocities are free, and apart from the symmetry
 //line (boundary 3) where only the horizontal velocity is pinned
 unsigned n_bound=mesh_pt()->nboundary();
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0);
       if (b!=3)
        {
         mesh_pt()->boundary_node_pt(b,n)->pin(1);
        }
      }
    }
  }

 //Boundary conditions on the SolidMesh
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       //Pin vertical displacemenet on the bottom
       if(b==0)
        {
         //mesh_pt()->boundary_node_pt(b,n)->pin_position(0);
         mesh_pt()->boundary_node_pt(b,n)->pin_position(1);
        }
       if((b==1) || (b==3))
        {
         //Pin horizontal displacement on the sizes
         mesh_pt()->boundary_node_pt(b,n)->pin_position(0);
        }
      }
    }
  }

 //Constrain all nodes only to move vertically (not horizontally)
 {
  unsigned n_node = mesh_pt()->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->node_pt(n)->pin_position(0);
   }
 }

 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
 //dynamic_cast<FixedVolumeSpineLineFluidInterfaceElement<SpineElement<ELEMENT> >*>
 //(mesh_pt()->interface_element_pt(n_interface-1))->
 // set_contact_angle_right(&Angle)

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))
  ->set_contact_angle(&Angle);
 
 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->ca_pt() = &Ca;

 dynamic_cast<FluidInterfaceBoundingElement*>(
  mesh_pt()->element_pt(mesh_pt()->nelement()-1))->wall_unit_normal_fct_pt() 
  =  &Global_Physical_Variables::wall_unit_normal_fct;

 
 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
  }

//======================================================================
/// Perform a parameter study. Pass name of output directory as 
/// a string
//======================================================================
template<class ELEMENT>
void ElasticCapProblem<ELEMENT>::parameter_study(const string& dir_name)
{
 // Create DocInfo object (allows checking if output directory exists)
 DocInfo doc_info;
 doc_info.set_directory(dir_name);
 doc_info.number()=0;


 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);
 Trace_file << "VARIABLES=\"<greek>a</greek><sub>prescribed</sub>\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>left</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>right</sub>\",";
 Trace_file << "\"p<sub>fluid</sub>-p<sub>ext</sub>\",";
 Trace_file << "\"<greek>D</greek>p<sub>exact</sub>\"";
 Trace_file << std::endl;

 //Gradually increase the contact angle
 for(unsigned i=0;i<6;i++)
  {
   //Solve the problem
   steady_newton_solve();

   //Output result
   doc_solution(doc_info);

   // Bump up counter
   doc_info.number()++;

   //Decrease the contact angle
   Angle -= 5.0*MathematicalConstants::Pi/180.0;
  }

} 




//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ElasticCapProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 //Output domain
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Number of interface elements
 unsigned ninterface=mesh_pt()->ninterface_element();
 //Find number of nodes in the last interface element
 unsigned np = mesh_pt()->interface_element_pt(ninterface-1)->nnode();
 // Doc
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << mesh_pt()->interface_element_pt(0)->node_pt(0)->x(1)
            << " " 
            << mesh_pt()->interface_element_pt(ninterface-1)
               ->node_pt(np-1)->x(1)
            << " "
//             << dynamic_cast<
//   FixedVolumeElasticLineFluidInterfaceElement<ELEMENT>*>(
//                mesh_pt()->interface_element_pt(0))->
//                actual_contact_angle_left()*
//                180.0/MathematicalConstants::Pi
            << " " ;
 // Trace_file << " "  
//             << dynamic_cast<FixedVolumeElasticLineFluidInterfaceElement<
//   ELEMENT>*>(
//                mesh_pt()->interface_element_pt(ninterface-1))->
//                actual_contact_angle_right()*180.0/MathematicalConstants::Pi 
//             << " ";
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(mesh_pt()->bulk_element_pt(0))->p_nst(0)-
               External_pressure_data_pt->value(0);
 Trace_file << " " << -2.0*cos(Angle)/Ca;
 Trace_file << std::endl;

}

 


//======================================================================
///Main driver: Build problem and initiate parameter study
//======================================================================
int main()
{

 // Solve the problem twice, once hijacking an internal, once
 // hijacking the external pressure
 for (unsigned i=0;i<2;i++)
  {

   bool hijack_internal=false;
   if (i==1) hijack_internal=true;
   //Construct the problem
   CapProblem<Hijacked<QCrouzeixRaviartElement<2> >  > problem(hijack_internal);

   string dir_name="RESLT_hijacked_external";
   if (i==1) dir_name="RESLT_hijacked_internal";

   //Do parameter study
   problem.parameter_study(dir_name);
   
  }


  // Solve the elastic problem twice, once hijacking an internal, once
 // hijacking the external pressure
 for (unsigned i=0;i<2;i++)
  {

   bool hijack_internal=false;
   if (i==1) hijack_internal=true;
   //Construct the problem
   ElasticCapProblem<Hijacked<
     PseudoSolidNodeUpdateElement<QCrouzeixRaviartElement<2>,
    QPVDElementWithPressure<2> > > >  problem(hijack_internal);

   string dir_name="RESLT_elastic_hijacked_external";
   if (i==1) dir_name="RESLT_elastic_hijacked_internal";

   //Do parameter study
   problem.parameter_study(dir_name);

  }


}









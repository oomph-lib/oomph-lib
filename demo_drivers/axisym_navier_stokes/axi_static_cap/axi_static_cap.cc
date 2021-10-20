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
// Driver code for an axisymmetric free-surface hydrostatics problem.
// The system consists of a layer of fluid 
// in a domain of height 1 and radius 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// spherical shells and the pressure jump across the interface should be
// 2 cos(alpha)/0.5 = 4 cos(alpha)/Ca.

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"

// The mesh
#include "meshes/single_layer_spine_mesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//====start_of_namespace=================================
/// Namespace for phyical parameters
//=======================================================
namespace Global_Physical_Variables
{

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

 ///Direction of the wall normal vector
 Vector<double> Wall_normal;

 ///  Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }

} //end_of_namespace



//============================================================================
///A Problem class that solves the Navier--Stokes equations + free surface
///in an axisymmetric geometry using a spine-based node update
//============================================================================
template<class ELEMENT>
class CapProblem : public Problem
{

public:

 //Constructor: Boolean flag indicates if volume constraint is
 //applied by hijacking internal or external pressure
 CapProblem(const bool& hijack_internal);

 //Destructor: clean up all allocated memory
 ~CapProblem();

 /// Perform a parameter study: Solve problem for a range of contact angles
 /// Pass name of output directory as a string
 void parameter_study(const string& dir_name);

 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {Bulk_mesh_pt->node_update();}

 /// Create the volume constraint elements
 void create_volume_constraint_elements();

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

 /// The bulk mesh of fluid elements
 SingleLayerSpineMesh<SpineElement<ELEMENT> >* Bulk_mesh_pt;

 /// The mesh for the interface elements
 Mesh* Surface_mesh_pt;

 /// The mesh for the element at the contact point
 Mesh* Point_mesh_pt;

 /// The volume constraint mesh 
 Mesh* Volume_constraint_mesh_pt;

 /// Trace file
 ofstream Trace_file;

 /// Data object whose single value stores the external pressure
 Data* External_pressure_data_pt;

 /// Data that is traded for the volume constraint
 Data* Traded_pressure_data_pt;

};



//======================================================================
/// Constructor: Pass boolean flag to indicate if the volume
/// constraint is applied by hijacking an internal pressure
/// or the external pressure
//======================================================================
template<class ELEMENT>
CapProblem<ELEMENT>::CapProblem(const bool& hijack_internal) :
 Ca(2.1),  //Initialise value of Ca to some random value
 Volume(0.125),  //Initialise the value of the volume:
                 //the physical volume divided by 2pi
 Pext(1.23),  //Initialise the external pressure to some random value
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 //Set the wall normal
 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = 1.0;
 Global_Physical_Variables::Wall_normal[1] = 0.0;

 // Number of elements in the horizontal direction
 unsigned nx=4;

 // Number of elements in the vertical direction
 unsigned nh=4;

 // Halfwidth of domain
 double half_width=0.5;

 //Construct bulk mesh
 Bulk_mesh_pt = 
  new SingleLayerSpineMesh<SpineElement<ELEMENT> >(nx,nh,half_width,1.0);

 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
 //Also create storage for a point mesh that contain the single element
 //responsible for  enforcing the contact angle condition
 Point_mesh_pt = new Mesh;

 //Loop over the horizontal elements adjacent to the upper surface
 //(boundary 2) and create the surface elements
 unsigned n_boundary_element = Bulk_mesh_pt->nboundary_element(2);
 for(unsigned e=0;e<n_boundary_element;e++)
  {
   //Construct a new 1D line element adjacent to boundary 2
   FiniteElement *interface_element_pt =
    new SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >(
     Bulk_mesh_pt->boundary_element_pt(2,e),
     Bulk_mesh_pt->face_index_at_boundary(2,e));
   
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 

   //Find the (single) node that is on the solid boundary (boundary 1)
   unsigned n_node = interface_element_pt->nnode();
   //We only need to check the right-hand nodes (because I know the
   //ordering of the nodes within the element)
   if(interface_element_pt->node_pt(n_node-1)->is_on_boundary(1))
    {
     //Make the point (contact) element from right-hand edge of the element
     FiniteElement* point_element_pt = 
      dynamic_cast<SpineAxisymmetricFluidInterfaceElement<
       SpineElement<ELEMENT> >*>(interface_element_pt)
      ->make_bounding_element(1);

     //Add it to the mesh
     this->Point_mesh_pt->add_element_pt(point_element_pt);
    }
  }

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Pext);

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
   Traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->element_pt(0))->hijack_internal_value(0,0);
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
   Traded_pressure_data_pt = External_pressure_data_pt;

   // Since the external pressure is "traded" for the volume constraint,
   // it no longer sets the overall pressure, and we 
   // can add an arbitrary constant to all pressures. To make 
   // the solution unique, we pin a single pressure value in the bulk: 
   // We arbitrarily set the pressure dof 0 in element 0 to zero.
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->fix_pressure(0,0.0);
  }



 // Loop over the elements on the interface to pass pointer to Ca and
 // the pointers to the Data items that contains the single
 // (pressure) value that is "traded" for the volume constraint 
 // and the external pressure
 unsigned n_interface = Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_interface;e++)
  {
   //Cast to a 1D element
   SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >*el_pt=
    dynamic_cast<SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >*>
    (Surface_mesh_pt->element_pt(e));

   //Set the Capillary number
   el_pt->ca_pt() = &Ca;
   
   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(External_pressure_data_pt);
  }

 
 //Set the boundary conditions

 //Pin the velocities on all boundaries apart from the free surface
 //(boundary 2) where all velocities are free, and apart from the symmetry
 //line (boundary 3) where only the horizontal velocity is pinned
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
       if (b!=3)
        {
         Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
        }
      }
    }
  }

 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
 FluidInterfaceBoundingElement *contact_angle_element_pt
  = dynamic_cast<FluidInterfaceBoundingElement*>(
   Point_mesh_pt->element_pt(0));

 contact_angle_element_pt->set_contact_angle(&Angle);
 contact_angle_element_pt->ca_pt() = &Ca;
 contact_angle_element_pt->wall_unit_normal_fct_pt() 
  =  &Global_Physical_Variables::wall_unit_normal_fct;

 //Now add the volume constraint
 create_volume_constraint_elements();

 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Surface_mesh_pt);
 this->add_sub_mesh(Point_mesh_pt);
 this->add_sub_mesh(Volume_constraint_mesh_pt);

 this->build_global_mesh();

 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
}


//==========================================================================
/// Destructor. Make sure to clean up all allocated memory, so that multiple
/// instances of the problem don't lead to excessive memory usage.
//==========================================================================
template<class ELEMENT>
CapProblem<ELEMENT>::~CapProblem() 
{
 //Loop over the volume mesh and delete the elements
 unsigned n_element = Volume_constraint_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {delete Volume_constraint_mesh_pt->element_pt(e);}
 //Now flush the storage
 Volume_constraint_mesh_pt->flush_element_and_node_storage();
 //Now delete the mesh
 delete Volume_constraint_mesh_pt;

 //Delete the traded pressure if not the same as the external pressure
 if(Traded_pressure_data_pt!=External_pressure_data_pt)
  {delete Traded_pressure_data_pt;}
  //Next delete the external data
 delete External_pressure_data_pt;

 //Loop over the point mesh and delete the elements
 n_element = Point_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++) 
  {delete Point_mesh_pt->element_pt(e);}
 //Now flush the storage
 Point_mesh_pt->flush_element_and_node_storage();
 //Then delete the mesh
 delete Point_mesh_pt;

 //Loop over the surface mesh and delete the elements
 n_element = Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++) 
  {delete Surface_mesh_pt->element_pt(e);}
 //Now flush the storage
 Surface_mesh_pt->flush_element_and_node_storage();
 //Then delete the mesh
 delete Surface_mesh_pt;

 //Then delete the bulk mesh
 delete Bulk_mesh_pt;
}


//=======================================================================
/// Create the volume constraint elements
//========================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::create_volume_constraint_elements()
{
 //The single volume constraint element
 Volume_constraint_mesh_pt = new Mesh;
 VolumeConstraintElement* vol_constraint_element = 
  new VolumeConstraintElement(&Volume,Traded_pressure_data_pt,0);
 Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);

 //Loop over all boundaries (or just 1 and 2 why?)
 for(unsigned b=0;b<4;b++)
  {
   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }
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
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

 //Output domain in paraview format
 sprintf(filename,"%s/soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_paraview(some_file,npts);
 some_file.close();

 // Number of spines
 unsigned nspine=Bulk_mesh_pt->nspine();

 // Doc
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(0)->height();
 Trace_file << " "  << Bulk_mesh_pt->spine_pt(nspine-1)->height();
 Trace_file << " " 
            << dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))
  ->p_axi_nst(0) - External_pressure_data_pt->value(0);
 Trace_file << " " << -4.0*cos(Angle)/Ca;
 Trace_file << std::endl;
 
}


//===========start_of_pseudo_elastic_class====================================
/// A class that solves the Navier--Stokes equations
///to compute the shape of a static interface in a rectangular container
///with imposed contact angle at the boundary.
//============================================================================
template<class ELEMENT>
class PseudoSolidCapProblem : public Problem
{
public:

 ///Constructor: Boolean flag indicates if volume constraint is
 ///applied by hijacking internal or external pressure
 PseudoSolidCapProblem(const bool& hijack_internal);

 /// Destructor: clean up memory allocated by the object
 ~PseudoSolidCapProblem();

 /// Peform a parameter study: Solve problem for a range of contact angles
 /// Pass name of output directory as a string
 void parameter_study(const string& dir_name);

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:
 
 /// Create the free surface elements
 void create_free_surface_elements();

 /// Create the volume constraint elements
 void create_volume_constraint_elements();

 /// Create the contact angle element
 void create_contact_angle_element();

 /// The Capillary number 
 double Ca;

 /// The prescribed volume of the fluid
 double Volume;

 /// The external pressure
 double Pext;

 /// The contact angle
 double Angle;

 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 /// Data object whose single value stores the external pressure
 Data* External_pressure_data_pt;

 // Pointer to the (single valued) Data item that
 // will contain the pressure value that we're
 // trading for the volume constraint
 Data* Traded_pressure_data_pt;

 /// Trace file
 ofstream Trace_file;

 ///Storage for the bulk mesh
 Mesh* Bulk_mesh_pt;

 /// Storage for the free surface mesh
 Mesh* Free_surface_mesh_pt;

 /// Storage for the element bounding the free surface
 Mesh* Free_surface_bounding_mesh_pt;

 /// Storage for the elements that compute the enclosed volume
 Mesh* Volume_computation_mesh_pt;

 /// Storage for the volume constraint
 Mesh* Volume_constraint_mesh_pt;

}; //end_of_pseudo_solid_problem_class



//============start_of_constructor=====================================
/// Constructor: Pass boolean flag to indicate if the volume
/// constraint is applied by hijacking an internal pressure
/// or the external pressure
//======================================================================
template<class ELEMENT>
PseudoSolidCapProblem<ELEMENT>::PseudoSolidCapProblem(const bool& hijack_internal) :
 Ca(2.1),       //Initialise value of Ca to some random value
 Volume(0.125),   //Initialise the value of the volume
                  // the physical volume divided by 2pi
 Pext(1.23),    //Initialise the external pressure to some random value
 Angle(0.5*MathematicalConstants::Pi) //Initialise the contact angle
{
 //Set the wall normal
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
 Bulk_mesh_pt = new ElasticRectangularQuadMesh<ELEMENT>(nx,nh,half_width,1.0);

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Pext);

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
   Traded_pressure_data_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->element_pt(0))->hijack_internal_value(0,0);
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
   Traded_pressure_data_pt = External_pressure_data_pt;

   // Since the external pressure is "traded" for the volume constraint,
   // it no longer sets the overall pressure, and we 
   // can add an arbitrary constant to all pressures. To make 
   // the solution unique, we pin a single pressure value in the bulk: 
   // We arbitrarily set the pressure dof 0 in element 0 to zero.
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->fix_pressure(0,0.0);
  }

 //Set the constituive law
 Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 //Loop over the elements to set the consitutive law and jacobian
 unsigned n_bulk = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_bulk;e++)
  {
   ELEMENT* el_pt = 
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   el_pt->constitutive_law_pt() = Constitutive_law_pt;
  }

 //Set the boundary conditions

 //Fluid velocity conditions
 //Pin the velocities on all boundaries apart from the free surface
 //(boundary 2) where all velocities are free, and apart from the symmetry
 //line (boundary 3) where only the horizontal velocity is pinned
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
       if (b!=3)
        {
         Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
        }
      }
    }
  } //end_of_fluid_boundary_conditions

 //PesudoSolid boundary conditions
 for (unsigned b=0;b<n_bound;b++)
  {
   if (b!=2)
    {
     //Find the number of nodes on the boundary
     unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
     //Loop over the nodes on the boundary
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       //Pin vertical displacement on the bottom
       if(b==0)
        {
         static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n))
          ->pin_position(1);
        }
       if((b==1) || (b==3))
        {
         //Pin horizontal displacement on the sizes
         static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(b,n))
          ->pin_position(0);
        }
      }
    } //end_of_solid_boundary_conditions
  }

 //Constrain all nodes only to move vertically (not horizontally)
 {
  unsigned n_node = Bulk_mesh_pt->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n))->pin_position(0);
   }
 } //end_of_constraint
 
 //Create the free surface elements
 create_free_surface_elements();

 //Create the volume constraint elements
 create_volume_constraint_elements();

 //Need to make the bounding element
 create_contact_angle_element();

 //Now need to add all the meshes
 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Free_surface_mesh_pt);
 this->add_sub_mesh(Volume_computation_mesh_pt);
 this->add_sub_mesh(Volume_constraint_mesh_pt);
 this->add_sub_mesh(Free_surface_bounding_mesh_pt);

 //and build the global mesh
 this->build_global_mesh();
 
 //Setup all the equation numbering and look-up schemes 
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 
 
} //end_of_constructor


//==========================================================================
/// Destructor. Make sure to clean up all allocated memory, so that multiple
/// instances of the problem don't lead to excessive memory usage.
//==========================================================================
template<class ELEMENT>
PseudoSolidCapProblem<ELEMENT>::~PseudoSolidCapProblem() 
{
 //Delete the contact angle element
 delete Free_surface_bounding_mesh_pt->element_pt(0);
 Free_surface_bounding_mesh_pt->flush_element_and_node_storage();
 delete Free_surface_bounding_mesh_pt;
 //Delete the volume constraint mesh
 delete Volume_constraint_mesh_pt;
 //Delete the surface volume computation elements
 unsigned n_element = Volume_computation_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {delete Volume_computation_mesh_pt->element_pt(e);}
 //Now flush the storage
 Volume_computation_mesh_pt->flush_element_and_node_storage();
 //Now delete the mesh
 delete Volume_computation_mesh_pt;
 //Delete the free surface elements
 n_element = Free_surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {delete Free_surface_mesh_pt->element_pt(e);}
 //Now flush the storage
 Free_surface_mesh_pt->flush_element_and_node_storage();
 //Now delete the mesh
 delete Free_surface_mesh_pt;

 //Delete the constitutive law
 delete Constitutive_law_pt;

 //If not the same as the external pressure, delete the traded pressure
 if(Traded_pressure_data_pt!=External_pressure_data_pt)
  {delete Traded_pressure_data_pt;}
 //Next delete the external data
 delete External_pressure_data_pt;
 //Then delete the bulk mesh
 delete Bulk_mesh_pt;
}

//============create_free_surface_element================================
/// Create the free surface elements
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_free_surface_elements()
{
 //Allocate storage for the free surface mesh
 Free_surface_mesh_pt = new Mesh;

 //Loop over boundary 2 and create the interface elements
 unsigned b=2;

 // How many bulk fluid elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
 // Loop over the bulk fluid elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk fluid element that is 
   // adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Create new element
   ElasticAxisymmetricFluidInterfaceElement<ELEMENT>* el_pt =
    new ElasticAxisymmetricFluidInterfaceElement<ELEMENT>(
     bulk_elem_pt,face_index);   
     
   // Add it to the mesh
   Free_surface_mesh_pt->add_element_pt(el_pt);
   
   //Add the appropriate boundary number
   el_pt->set_boundary_number_in_bulk_mesh(b);

   //Add the capillary number
   el_pt->ca_pt() = &Ca;

   //Add the external pressure data
   el_pt->set_external_pressure_data(External_pressure_data_pt);
  } 

}


//============start_of_create_volume_constraint_elements==================
/// Create the volume constraint elements
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_volume_constraint_elements()
{
 //Build the single volume constraint element
 Volume_constraint_mesh_pt = new Mesh;
 VolumeConstraintElement* vol_constraint_element = 
  new VolumeConstraintElement(&Volume,Traded_pressure_data_pt,0);
 Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);
  
  //Now create the volume computation elements
  Volume_computation_mesh_pt = new Mesh;

  //Loop over all boundaries
  for(unsigned b=0;b<4;b++)
   {
    // How many bulk fluid elements are adjacent to boundary b?
    unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
    
    // Loop over the bulk fluid elements adjacent to boundary b?
    for(unsigned e=0;e<n_element;e++)
     {
      // Get pointer to the bulk fluid element that is 
      // adjacent to boundary b
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
       Bulk_mesh_pt->boundary_element_pt(b,e));
      
      //Find the index of the face of element e along boundary b
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
      
      // Create new element
      ElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>* el_pt =
       new ElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(
        bulk_elem_pt,face_index);   

      //Set the "master" volume control element
      el_pt->set_volume_constraint_element(vol_constraint_element);
     
      // Add it to the mesh
      Volume_computation_mesh_pt->add_element_pt(el_pt);     
     }
   }
} //end_of_create_volume_constraint_elements

//==========start_of_create_contact_angle_elements========================
/// Create the contact angle element
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::create_contact_angle_element()
{
 Free_surface_bounding_mesh_pt = new Mesh;

 //Find the element at the end of the free surface
 //The elements are assigned in order of increasing x coordinate
 unsigned n_free_surface = Free_surface_mesh_pt->nelement();
  
 //Make the bounding element for the contact angle constraint
 //which works because the order of elements in the mesh is known
 FluidInterfaceBoundingElement* el_pt = 
  dynamic_cast<ElasticAxisymmetricFluidInterfaceElement<ELEMENT>*> 
  (Free_surface_mesh_pt->element_pt(n_free_surface-1))->
  make_bounding_element(1);
 
 //Set the contact angle (strong imposition)
 el_pt->set_contact_angle(&Angle);
 
 //Set the capillary number
 el_pt->ca_pt() = &Ca;
 
 //Set the wall normal of the external boundary
 el_pt->wall_unit_normal_fct_pt() 
  =  &Global_Physical_Variables::wall_unit_normal_fct;

 //Add the element to the mesh
 Free_surface_bounding_mesh_pt->add_element_pt(el_pt);

} //end_of_create_contact_angle_element



//================start_of_parameter_study===========================
/// Perform a parameter study. Pass name of output directory as 
/// a string
//======================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::parameter_study(const string& dir_name)
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
 Trace_file << "\"p<sub>fluid</sub>-p<sub>ext</sub>\",";
 Trace_file << "\"<greek>D</greek>p<sub>exact</sub>\"";
 Trace_file << std::endl;

 //Solve the problem for six different contact angles
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




//==============start_of_doc_solution=====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PseudoSolidCapProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 //Output stream
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 //Output domain
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();


 // Number of interface elements
 unsigned ninterface=Free_surface_mesh_pt->nelement();
 //Find number of nodes in the last interface element
 unsigned np = Free_surface_mesh_pt->finite_element_pt(ninterface-1)->nnode();
 // Document the contact angle (in degrees), the height of the interface at 
 // the centre of the container, the height at the wall, the computed 
 // pressure drop across the interface and 
 // the analytic prediction of the pressure drop.
 Trace_file << Angle*180.0/MathematicalConstants::Pi;
 Trace_file << " "  << Free_surface_mesh_pt->finite_element_pt(0)
  ->node_pt(0)->x(1)
            << " " 
            << Free_surface_mesh_pt->finite_element_pt(ninterface-1)
  ->node_pt(np-1)->x(1);
 Trace_file << " " 
            << 
  dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->p_axi_nst(0)-
               External_pressure_data_pt->value(0);
 Trace_file << " " << -4.0*cos(Angle)/Ca;
 Trace_file << std::endl;

} //end_of_doc_solution



 
//===start_of_main=======================================================
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
   CapProblem<Hijacked<AxisymmetricQCrouzeixRaviartElement>  > 
    problem(hijack_internal);

   string dir_name="RESLT_hijacked_external";
   if (i==1) dir_name="RESLT_hijacked_internal";

   //Do parameter study
   problem.parameter_study(dir_name);
   
  }


  // Solve the pseudosolid problem twice, once hijacking an internal, once
 // hijacking the external pressure
 for (unsigned i=0;i<2;i++)
  {
   bool hijack_internal=false;
   if (i==1) hijack_internal=true;
   //Construct the problem
   PseudoSolidCapProblem<Hijacked<
     PseudoSolidNodeUpdateElement<AxisymmetricQCrouzeixRaviartElement,
    QPVDElementWithPressure<2> > > >  problem(hijack_internal);

   string dir_name="RESLT_elastic_hijacked_external";
   if (i==1) dir_name="RESLT_elastic_hijacked_internal";

   //Do parameter study
   problem.parameter_study(dir_name);
  }

} //end_of_main

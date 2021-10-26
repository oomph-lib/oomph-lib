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
#ifndef OOMPH_MY_ALGEBRAIC_COLLAPSIBLE_CHANNEL_MESH
#define OOMPH_MY_ALGEBRAIC_COLLAPSIBLE_CHANNEL_MESH

//Include the mesh
#include "meshes/collapsible_channel_mesh.h"





/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////

namespace oomph
{


//===========start_algebraic_mesh==================================
/// Collapsible channel mesh with algebraic node update
//=================================================================
template<class ELEMENT>
class MyAlgebraicCollapsibleChannelMesh : 
  public AlgebraicMesh,
  public virtual CollapsibleChannelMesh<ELEMENT>
{ 

public: 

 /// Constructor: Pass number of elements in upstream/collapsible/
 /// downstream segment and across the channel; lengths of upstream/
 /// collapsible/downstream segments and width of channel, pointer to 
 /// GeomObject that defines the collapsible segment, and pointer to 
 /// TimeStepper (defaults to the default timestepper, Steady). 
 MyAlgebraicCollapsibleChannelMesh(const unsigned& nup, 
                                   const unsigned& ncollapsible, 
                                   const unsigned& ndown, 
                                   const unsigned& ny, 
                                   const double& lup, 
                                   const double& lcollapsible, 
                                   const double& ldown, 
                                   const double& ly,
                                   GeomObject* wall_pt,
                                   TimeStepper* time_stepper_pt=
                                   &Mesh::Default_TimeStepper) :
  CollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                  lup, lcollapsible, ldown, ly,
                                  wall_pt,
                                  time_stepper_pt)
  {
   // Setup algebraic node update operations
   setup_algebraic_node_update();
  }
 


 /// Constructor: Pass number of elements in upstream/collapsible/
 /// downstream segment and across the channel; lengths of upstream/
 /// collapsible/downstream segments and width of channel, pointer to 
 /// GeomObject that defines the collapsible segment, function pointer
 /// to "boundary layer squash function", and pointer to 
 /// TimeStepper (defaults to the default timestepper, Steady). 
 MyAlgebraicCollapsibleChannelMesh(
  const unsigned& nup, 
  const unsigned& ncollapsible, 
  const unsigned& ndown, 
  const unsigned& ny, 
  const double& lup, 
  const double& lcollapsible, 
  const double& ldown, 
  const double& ly,
  GeomObject* wall_pt,
  CollapsibleChannelDomain::BLSquashFctPt bl_squash_function_pt,
  TimeStepper* time_stepper_pt=
  &Mesh::Default_TimeStepper) :
  CollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                  lup, lcollapsible, ldown, ly,
                                  wall_pt,
                                  time_stepper_pt)
  {
   // Set boundary layer squash function
   this->Domain_pt->bl_squash_fct_pt()=bl_squash_function_pt;

   // Do MacroElement-based node update
   CollapsibleChannelMesh<ELEMENT>::node_update();

   // Setup algebraic node update operations
   setup_algebraic_node_update();
  }
 
 
 /// Destructor: empty
 virtual ~MyAlgebraicCollapsibleChannelMesh(){}
 
 /// Function pointer for function that squashes
 /// the mesh near the walls. Default trivial mapping (the identity)
 /// leaves vertical nodal positions unchanged. Mapping is
 /// used in underlying CollapsibleChannelDomain.
 /// This (deliberately broken) function overloads the one
 /// in the CollapsibleChannelMesh base class. 
 CollapsibleChannelDomain::BLSquashFctPt& bl_squash_fct_pt()
  {
   std::ostringstream error_message;
   error_message 
    << "It does not make sense to set the bl_squash_fct_pt \n"
    << "outside the constructor as it's only used to set up the \n"
    << "algebraic remesh data when the algebraic mesh is first built. \n";
    std::string function_name =
    "MyAlgebraicCollapsibleChannelMesh::bl_squash_fct_pt()\n";
   
   throw OomphLibError(error_message.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);

   // Dummy return
   return Dummy_fct_pt;
  }


 /// Update nodal position at time level t (t=0: present; 
 /// t>0: previous)
 void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt);
 
 /// Update the geometric references that are used 
 /// to update node after mesh adaptation.
 /// Empty -- no update of node update required
 void update_node_update(AlgebraicNode*& node_pt) {}

protected:

 /// Function to setup the algebraic node update
 void setup_algebraic_node_update();


 /// Dummy function pointer 
 CollapsibleChannelDomain::BLSquashFctPt Dummy_fct_pt;

};





/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//===========start_refineable_algebraic_collapsible_channel_mesh======
/// Refineable version of the CollapsibleChannel mesh with
/// algebraic node update.
//====================================================================
template<class ELEMENT>
class MyRefineableAlgebraicCollapsibleChannelMesh : 
  public RefineableQuadMesh<ELEMENT>,
  public virtual MyAlgebraicCollapsibleChannelMesh<ELEMENT>
{ 

public: 


 /// Constructor: Pass number of elements in upstream/collapsible/
 /// downstream segment and across the channel; lengths of upstream/
 /// collapsible/downstream segments and width of channel, pointer to 
 /// GeomObject that defines the collapsible segment, function pointer
 /// to "boundary layer squash function", and pointer to 
 /// TimeStepper (defaults to the default timestepper, Steady). 
 MyRefineableAlgebraicCollapsibleChannelMesh(const unsigned& nup, 
                                             const unsigned& ncollapsible, 
                                             const unsigned& ndown, 
                                             const unsigned& ny, 
                                             const double& lup, 
                                             const double& lcollapsible, 
                                             const double& ldown, 
                                             const double& ly,
                                             GeomObject* wall_pt,
                                             TimeStepper* time_stepper_pt=
                                             &Mesh::Default_TimeStepper) :
  CollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                  lup, lcollapsible, ldown, ly,
                                  wall_pt,
                                  time_stepper_pt),
  MyAlgebraicCollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                             lup, lcollapsible, ldown, ly,
                                             wall_pt,
                                             time_stepper_pt)
  {
   // Build quadtree forest
   this->setup_quadtree_forest();
  }
 



 /// Constructor: Pass number of elements in upstream/collapsible/
 /// downstream segment and across the channel; lengths of upstream/
 /// collapsible/downstream segments and width of channel, pointer to 
 /// GeomObject that defines the collapsible segment, function pointer
 /// to "boundary layer squash function",  and pointer to 
 /// TimeStepper (defaults to the default timestepper, Steady). 
 MyRefineableAlgebraicCollapsibleChannelMesh(
  const unsigned& nup, 
  const unsigned& ncollapsible, 
  const unsigned& ndown, 
  const unsigned& ny, 
  const double& lup, 
  const double& lcollapsible, 
  const double& ldown, 
  const double& ly,
  GeomObject* wall_pt,
  CollapsibleChannelDomain::BLSquashFctPt bl_squash_function_pt,
  TimeStepper* time_stepper_pt=
  &Mesh::Default_TimeStepper) :
  CollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                  lup, lcollapsible, ldown, ly,
                                  wall_pt,
                                  time_stepper_pt),
  MyAlgebraicCollapsibleChannelMesh<ELEMENT>(nup, ncollapsible, ndown, ny,
                                             lup, lcollapsible, ldown, ly,
                                             wall_pt,
                                             bl_squash_function_pt,
                                             time_stepper_pt)
  {
   // Build quadtree forest
   this->setup_quadtree_forest();
  }

};








//============start_setup_algebraic_node_update====================
/// Setup algebraic mesh update -- assumes that mesh has
/// initially been set up with the wall in its undeformed position.
//=================================================================
template<class ELEMENT>
void MyAlgebraicCollapsibleChannelMesh<ELEMENT>::setup_algebraic_node_update()
{

 // Extract some reference lengths from the CollapsibleChannelDomain. 
 double l_up=this->domain_pt()->l_up();
 double l_collapsible=this->domain_pt()->l_collapsible();

 // Loop over all nodes in mesh
 unsigned nnod=this->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   // Get pointer to node
   AlgebraicNode* nod_pt=node_pt(j);

   // Get coordinates
   double x=nod_pt->x(0);
   double y=nod_pt->x(1);

   // Check if the node is in the collapsible part:
   if ( (x>=l_up) && (x<=(l_up+l_collapsible)) )
    {

     // Get zeta coordinate on the undeformed wall
     Vector<double> zeta(1);
     zeta[0]=x-l_up;

     // Get position vector to wall:
     Vector<double> r_wall(2);
     this->Wall_pt->position(zeta,r_wall);


     // Sanity check: Confirm that the wall is in its undeformed position
#ifdef PARANOID
     if ((std::abs(r_wall[0]-x)>1.0e-15)&&(std::abs(r_wall[1]-y)>1.0e-15))
      {
       std::ostringstream error_stream;
       error_stream 
        << "Wall must be in its undeformed position when\n"
        << "algebraic node update information is set up!\n "
        << "x-discrepancy: " << std::abs(r_wall[0]-x) << std::endl
        << "y-discrepancy: " << std::abs(r_wall[1]-y) << std::endl;
       
       throw OomphLibError(
        error_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
#endif       

     // Only a single geometric object is involved in the node update operation
     Vector<GeomObject*> geom_object_pt(1);

     // The wall geometric object
     geom_object_pt[0]=this->Wall_pt;
     
     // The update function requires three parameters:
     Vector<double> ref_value(3);
     
     // First reference value: Original x-position on the lower wall 
     ref_value[0]=r_wall[0];
     
     // Second reference value: Fractional position along 
     // straight line from the bottom (at the original x-position)
     // to the point on the wall
     ref_value[1]=y/r_wall[1];
     
     // Third reference value: Zeta coordinate on wall
     ref_value[2]=zeta[0];   

     // Setup algebraic update for node: Pass update information
     // to AlgebraicNode:
     nod_pt->add_node_update_info(
      this,               // mesh
      geom_object_pt,     // vector of geom objects
      ref_value);         // vector of ref. values
    }
   
  }

} //end of setup_algebraic_node_update


//=============start_of_algebraic_node_update======================
/// Perform algebraic mesh update at time level t (t=0: present; 
/// t>0: previous)
//=================================================================
template<class ELEMENT>
void MyAlgebraicCollapsibleChannelMesh<ELEMENT>::algebraic_node_update(
 const unsigned& t, AlgebraicNode*& node_pt)
{

 // Extract reference values for update by copy construction
 Vector<double> ref_value(node_pt->vector_ref_value());

 // Extract geometric objects for update by copy construction
 Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());
 
 // First reference value: Original x-position 
 double x_bottom=ref_value[0];
 
 // Second reference value: Fractional position along 
 // straight line from the bottom (at the original x-position)
 // to the point on the wall
 double fract=ref_value[1];
 
 // Third reference value: Zeta coordinate on wall
 Vector<double> zeta(1);
 zeta[0]=ref_value[2];
     
 // Pointer to wall geom object
 GeomObject* wall_pt=geom_object_pt[0];

 // Get position vector to wall at previous timestep t
 Vector<double> r_wall(2);
 wall_pt->position(t,zeta,r_wall);

 // Assign new nodal coordinates
 node_pt->x(t,0)=x_bottom+fract*(r_wall[0]-x_bottom);
 node_pt->x(t,1)=         fract* r_wall[1];

}




}








#endif

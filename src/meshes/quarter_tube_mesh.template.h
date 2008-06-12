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
#ifndef OOMPH_QUARTER_TUBE_MESH_HEADER
#define OOMPH_QUARTER_TUBE_MESH_HEADER

// Headers
#include "../generic/refineable_brick_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/brick_mesh.h"
#include "../generic/macro_element_node_update_element.h"


//Include the headers file for domain
#include "quarter_tube_domain.h"

namespace oomph
{


//====================================================================
/// \short 3D quarter tube mesh class.
/// The domain is specified by the GeomObject that identifies 
/// boundary 3. Non-refineable base version!
/// \n
/// The mesh boundaries are numbered as follows:
/// - Boundary 0: "Inflow" cross section; located along the 
///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$
///               on the geometric object that specifies the wall.
/// - Boundary 1: Plane x=0
/// - Boundary 2: Plane y=0
/// - Boundary 3: The curved wall
/// - Boundary 4: "Outflow" cross section; located along the 
///               line parametrised by \f$ \xi_0 =  \xi_0^{hi} \f$
///               on the geometric object that specifies the wall.
//====================================================================
template <class ELEMENT>
class QuarterTubeMesh : public virtual BrickMeshBase
{

public:

 /// \short Constructor: Pass pointer to geometric object that
 /// specifies the wall, start and end coordinates on the 
 /// geometric object, and the fraction along
 /// which the dividing line is to be placed, and the timestepper.
 /// Timestepper defaults to Steady dummy timestepper.
  QuarterTubeMesh(GeomObject* wall_pt,
                 const Vector<double>& xi_lo,
                 const double& fract_mid,
                 const Vector<double>& xi_hi,
                 const unsigned& nlayer,
                 TimeStepper* time_stepper_pt=
                 &Mesh::Default_TimeStepper);

 /// \short Destructor: empty
 virtual ~QuarterTubeMesh(){}



 /// \short Function pointer for function that squashes
 /// the outer macro elements towards 
 /// the wall by mapping the input value of the "radial" macro element
 /// coordinate to the return value (defined in the underlying Domain object)
 QuarterTubeDomain::BLSquashFctPt& bl_squash_fct_pt()
  {
   return Domain_pt->bl_squash_fct_pt();
  }


 /// \short Function pointer for function that hierher
 QuarterTubeDomain::AxialSpacingFctPt& axial_spacing_fct_pt()
  {
   return Domain_pt->axial_spacing_fct_pt();
  }

 /// Access function to underlying domain
 QuarterTubeDomain* domain_pt() const {return Domain_pt;}

protected:

 /// Pointer to domain
 QuarterTubeDomain* Domain_pt;

 /// Pointer to the geometric object that represents the curved wall
 GeomObject* Wall_pt;

 /// Lower limits for the coordinates along the wall
 Vector<double> Xi_lo;

 /// Fraction along wall where outer ring is to be divided
 double Fract_mid;

 /// Upper limits for the coordinates along the wall
 Vector<double> Xi_hi;

};




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////





//=============================================================
/// Adaptative version of the QuarterTubeMesh base mesh.
/// The domain is specified by the GeomObject that identifies 
/// boundary 3.
/// \n
/// The mesh boundaries are numbered as follows:
/// - Boundary 0: "Inflow" cross section; located along the 
///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$
///               on the geometric object that specifies the wall.
/// - Boundary 1: Plane x=0
/// - Boundary 2: Plane y=0
/// - Boundary 3: The curved wall
/// - Boundary 4: "Outflow" cross section; located along the 
///               line parametrised by \f$ \xi_0 =  \xi_0^{hi} \f$
///               on the geometric object that specifies the wall.
//=============================================================
template<class ELEMENT> 
class RefineableQuarterTubeMesh : public virtual QuarterTubeMesh<ELEMENT>,
                                public RefineableBrickMesh<ELEMENT>
                                
{

public :

/// \short Constructor for adaptive deformable quarter tube mesh class. 
/// The domain is specified by the GeomObject that 
/// identifies boundary 3. Pass pointer to geometric object that
/// specifies the wall, start and end coordinates on the 
/// geometric object, and the fraction along
/// which the dividing line is to be placed, and the timestepper.
/// Timestepper defaults to Steady dummy timestepper.
 RefineableQuarterTubeMesh(GeomObject* wall_pt,
                           const Vector<double>& xi_lo,
                           const double& fract_mid,
                           const Vector<double>& xi_hi,
                           const unsigned& nlayer,
                           TimeStepper* time_stepper_pt=
                           &Mesh::Default_TimeStepper):
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt)
  {
   // Loop over all elements and set macro element pointer
   for (unsigned ielem=0;ielem<QuarterTubeMesh<ELEMENT>::nelement();ielem++)
    {
     dynamic_cast<RefineableQElement<3>*>(
      QuarterTubeMesh<ELEMENT>::element_pt(ielem))->
      set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
    }
   

   // Setup Octree forest: Turn elements into individual octrees 
   // and plant in forest
   Vector<TreeRoot*> trees_pt;
   for (unsigned iel=0;iel<QuarterTubeMesh<ELEMENT>::nelement();iel++)
    {
     FiniteElement* el_pt=QuarterTubeMesh<ELEMENT>::finite_element_pt(iel);
     ELEMENT* ref_el_pt=dynamic_cast<ELEMENT*>(el_pt);
     OcTreeRoot* octree_root_pt=new OcTreeRoot(ref_el_pt);
     trees_pt.push_back(octree_root_pt);
    }
   this->Forest_pt = new OcTreeForest(trees_pt);

#ifdef PARANOID
   // Run self test
   unsigned success_flag=
    dynamic_cast<OcTreeForest*>(this->Forest_pt)->self_test();
   if (success_flag==0)
    {
     oomph_info << "Successfully built octree forest " << std::endl;
    }
   else
    {
     throw OomphLibError(
      "Trouble in building octree forest ",
      "RefineableQuarterTubeMesh::RefineableQuarterTubeMesh()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
  }
 
 /// \short Destructor: empty
 virtual ~RefineableQuarterTubeMesh(){}
 
}; 





////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// MacroElementNodeUpdate-version of RefineableQuarterTubeMesh
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class MacroElementNodeUpdateNode;

//========================================================================
/// MacroElementNodeUpdate version of RefineableQuarterTubeMesh
//========================================================================
template<class ELEMENT>
class MacroElementNodeUpdateRefineableQuarterTubeMesh : 
public virtual MacroElementNodeUpdateMesh, 
public virtual RefineableQuarterTubeMesh<ELEMENT>
{


public:


 /// \short Constructor: Pass pointer to geometric object, start and
 /// end coordinates on the geometric object and the fraction along
 /// which the dividing line is to be placed when updating the nodal positions,
 /// and timestepper (defaults to (Steady) default timestepper
 /// defined in Mesh). Setup the refineable mesh (by calling the
 /// constructor for the underlying  RefineableQuarterTubeMesh)
 /// and the algebraic node update functions for nodes.
 MacroElementNodeUpdateRefineableQuarterTubeMesh(GeomObject* wall_pt,
                           const Vector<double>& xi_lo,
                           const double& fract_mid,
                           const Vector<double>& xi_hi,
                           const unsigned& nlayer,
                           TimeStepper* time_stepper_pt=
                           &Mesh::Default_TimeStepper) :
  MacroElementNodeUpdateMesh(),
  RefineableQuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                     nlayer,time_stepper_pt),
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt)
  {

#ifdef PARANOID
   ELEMENT* el_pt=new ELEMENT;
   if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt)==0)
    {
     std::ostringstream error_message;
     error_message 
      << "Base class for ELEMENT in "
      << "MacroElementNodeUpdateRefineableQuarterTubeMesh needs" 
      << "to be of type MacroElementNodeUpdateElement!\n";
     error_message << "Whereas it is: typeid(el_pt).name()" 
                   << typeid(el_pt).name() 
                   << std::endl;
     
     std::string function_name =
      "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
     function_name += 
      "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh()";

     throw OomphLibError(error_message.str(),function_name,
                         OOMPH_EXCEPTION_LOCATION);
    }
   delete el_pt;
#endif

   // Setup all the information that's required for MacroElement-based
   // node update: Tell the elements that their geometry depends on the
   // fishback geometric object
   this->setup_macro_element_node_update();
  }

 /// \short Destructor: empty
 virtual ~MacroElementNodeUpdateRefineableQuarterTubeMesh(){}

 /// \short Resolve mesh update: Update current nodal
 /// positions via sparse MacroElement-based update.
 /// [Doesn't make sense to use this mesh with SolidElements anyway,
 /// so we buffer the case if update_all_solid_nodes is set to 
 /// true.]
 void node_update(const bool& update_all_solid_nodes=false)
  {
#ifdef PARANOID
   if (update_all_solid_nodes)
    {
     std::string error_message = 
      "Doesn't make sense to use an MacroElementNodeUpdateMesh with\n";
     error_message += 
      "SolidElements so specifying update_all_solid_nodes=true\n";
     error_message += "doesn't make sense either\n";

     std::string function_name =
      "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
     function_name += "node_update()";

     throw OomphLibError(error_message,function_name,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   MacroElementNodeUpdateMesh::node_update();
  }

  private:

 /// \short Setup all the information that's required for MacroElement-based
 /// node update: Tell the elements that their geometry depends on the
 /// geometric object that parametrises the wall
 void setup_macro_element_node_update()
  {
   unsigned n_element = this->nelement();
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from FiniteElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(this->element_pt(i));
     
#ifdef PARANOID
     // Check if cast is successful
     MacroElementNodeUpdateElementBase* m_el_pt=dynamic_cast<
      MacroElementNodeUpdateElementBase*>(el_pt);
     if (m_el_pt==0)
      {
       std::ostringstream error_message;
       error_message 
        << "Failed to upcast to MacroElementNodeUpdateElementBase\n";
       error_message
        << "Element must be derived from MacroElementNodeUpdateElementBase\n";
       error_message << "but it is of type " << typeid(el_pt).name();
       
       std::string function_name =
        "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
       function_name += "setup_macro_element_node_update()";
       
       throw OomphLibError(error_message.str(),function_name,
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif       
     // There's just one GeomObject
     Vector<GeomObject*> geom_object_pt(1);
     geom_object_pt[0] = this->Wall_pt;
     
     // Tell the element which geom objects its macro-element-based
     // node update depends on     
     el_pt->set_node_update_info(geom_object_pt);
    }
  }
};


}
#endif

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
//Header file for a classes used to represent projectable elements

// TODO:
// -- Force user to specify which field is interpolated isoparametrically
//    and can therefore be abused for coordinate projection
// -- Switch to more efficient solver for projection problem.
   


#ifndef PROJECTION_H
#define PROJECTION_H


#include "mesh.h"
#include "problem.h"
#include "multi_domain.h"
#include "shape.h"
#include "element_with_external_element.h"

namespace oomph
{

 //======================================================================
 /// Helper namespace for objects needed during projection between meshes
 //======================================================================
 namespace ProjectionHelper
 {

  /// \short Struct to backup value and pin status of various degrees of
  /// freedom 
  struct BackupValue
  {
   /// Data object
   Data* Data_pt;

   /// Index of value in Data object
   unsigned Index;

   /// Pin status of value
   bool Pinned;

   /// Original (current) value 
   double Value;
   
  };


  /// \short "Less-than" comparison operator class for BackupValue. Required
  /// so it can be inserted into set.
  class BackupValueLessThan
  {
    public:

   /// Less-than operator
   bool operator() (const BackupValue& lhs, const BackupValue& rhs) const
   {
    // Primary sort on pointer to Data
    if (lhs.Data_pt<rhs.Data_pt)
     {
      return true;
     }
    else if (lhs.Data_pt==rhs.Data_pt)
     {
      // If equality sort on index
      if (lhs.Index<rhs.Index)
       {
        return true;
       }
      // Data is the same and index isn't "less than"
      else
       {
        return false;
       }
     }
    // If we get here comparison based on Data_pt is unique
    else
     {
      return false;
     }
   }
   
  };

 }


//==================================================================
///Template-free Base class for projectable elements
//==================================================================
 class ProjectableElementBase 
{

protected:
 
 /// Field that is currently being projected
 unsigned Projected_field; 

 ///Time level we are projecting  (0=current values; >0: history values)
 unsigned Time_level_for_projection;

 ///\short When projecting the history values of the nodal coordinates,
 /// this is the coordinate we're projecting
 unsigned Projected_coordinate; 

 /// \short Bool to know if we do projection or not. If false (the default)
 /// we solve the element's "real" equations rather than the projection
 /// equations
 bool Do_projection;

 /// \short Bool to indicate if we're projecting the history values of the
 /// nodal coordinates (true) or the values themselves (false)
 bool Project_coordinates; 

 /// \short Store number of "external" interactions that were assigned to
 /// the element before doing the projection.
 unsigned Backup_ninteraction;

 /// \short Remember if the element includes external geometric data 
 /// when used in  non-projection mode (this is temporarily disabled during the
 /// projection)
 bool Backup_external_geometric_data;


 /// \short Remember if the element includes external data when used in 
 /// non-projection mode (this is temporarily disabled during the
 /// projection)
 bool Backup_external_interaction_data;


public:


 /// Constructor: Initialise data so that we don't project but solve
 /// the "normal" equations associated with the element.
 ProjectableElementBase() : Do_projection(false), Project_coordinates(false) {}
 
 ///Virtual destructor
 virtual ~ProjectableElementBase() { }


 /// \short Pure virtual function in which the element writer
 /// must specify the values associated with field fld. 
 /// The information is returned in a vector of pairs which comprise 
 /// the Data object and the value within it, that correspond to field fld. 
 /// E.g. in Taylor Hood elements the fld-th velocities are stored
 /// at the fld-th value of the nodes; the pressures (the DIM-th 
 /// field) are the DIM-th values at the vertex nodes etc. 
 virtual Vector<std::pair<Data*,unsigned> > 
  data_values_of_field(const unsigned& fld)=0;
 
 /// \short Number of fields of the problem, so e.g. for 2D Navier Stokes
 /// this would be 3 (for the two velocities and one pressure)
 virtual unsigned nfields_for_projection()=0;

 ///Number of history values to be stored for fld-th field
 virtual unsigned nhistory_values_for_projection(const unsigned &fld)=0;

 ///\short Number of history values to be stored when projecting 
 /// the history values of the nodal coordinates
 virtual unsigned nhistory_values_for_coordinate_projection()=0;

 /// \short Return number of values (pinned or not) associated with 
 /// field fld within the element. This must correspond to the
 /// number of shape functions returned in jacobian_and_shape_of_field(...).
 virtual unsigned nvalue_of_field(const unsigned &fld)=0;
 
 /// \short Return local equation numbers associated with value ivalue 
 /// of field fld within the element.
 virtual int local_equation(const unsigned &fld,const unsigned &ivalue)=0;

 /// \short Return Jacobian of mapping and the shape functions associated 
 /// with field fld. The number of shape functions must match the
 /// number of values specified in nvalue_of_field(...). For 
 /// Lagrange-type interpolations the shape functinos are simply
 /// the "normal" nodal shape functions; if the element contains
 /// internal Data that is not associated with shape functions,
 /// simply set the corresonding shape function to 1.
 virtual double jacobian_and_shape_of_field
 (const unsigned &fld, const Vector<double> &s, Shape &psi)=0;

 /// \short Return the fld-th field at local coordinates s
 /// at time-level time (time=0: current value; time>0: history values)
 virtual double get_field 
  (const unsigned &time, const unsigned &fld,const Vector<double> &s)=0;

};





//=====================================================================
/// Wrapper class for projectable elements. Adds "projectability"
/// to the underlying ELEMENT.
//=====================================================================
template<class ELEMENT>
class ProjectableElement : public virtual ELEMENT,
                           public virtual ProjectableElementBase,
                           public virtual ElementWithExternalElement
{
 
  protected:
 
 /// Overloaded version of fill_in_contribution_to_residuals
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  //Do projection
  if (Do_projection)
   {
    this->residual_for_projection
     (residuals,GeneralisedElement::Dummy_matrix, 0);
   }
  //solve problem normally
  else
   {
    ELEMENT::fill_in_contribution_to_residuals(residuals);
   }
 }
 
 
 ///Overloaded version of fill_in_contribution_to_jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {
  //Do projection
  if (Do_projection)
   {  
    this->residual_for_projection(residuals,jacobian,1);
   }
  else
   {
    ELEMENT::fill_in_contribution_to_jacobian(residuals,jacobian);
   }
 }
 

  public: 


 /// \short Constructor [this was only required explicitly
 /// from gcc 4.5.2 onwards...]
 ProjectableElement(){}

 /// Residual for the projection step. Flag indicates if we
 /// want the Jacobian (1) or not (0)
 void residual_for_projection(Vector<double> &residuals, 
                              DenseMatrix<double> &jacobian, 
                              const unsigned& flag)
  {
   unsigned n_dim=dim();
   
   //Allocate storage for local coordinates
   Vector<double> s(n_dim);

   //Current field
    unsigned fld=Projected_field;
    
   //Number of dof for current field
   unsigned n_value=nvalue_of_field(fld);  
  
   //Set the value of n_intpt
   const unsigned n_intpt = integral_pt()->nweight();
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {

     // Get the local coordinates of Gauss point
     for(unsigned i=0;i<n_dim;i++) s[i] = integral_pt()->knot(ipt,i);

     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Field shape function
     Shape psi(n_value);

     //Calculate jacobian and shape functions for that field
     double J=jacobian_and_shape_of_field
      (fld, s, psi);
              
     //Premultiply the weights and the Jacobian
     double W = w*J;
       
     //Value of field in current element
     unsigned now=0;
     double interpolated_value_proj = this->get_field(now,fld,s);
                 
     // Find same point in base mesh
     // Using external storage
     FiniteElement* other_el_pt=0;
     other_el_pt=this->external_element_pt(0,ipt);
     Vector<double> other_s(n_dim);
     other_s=this->external_element_local_coord(0,ipt);

     ProjectableElement<ELEMENT>* cast_el_pt = 
      dynamic_cast<ProjectableElement<ELEMENT>*>(other_el_pt);

     //Value of the interpolation of element located in base mesh
     double interpolated_value_bar = 0.0;

     if (Project_coordinates)
      {
       //Interpolation of current dimension
       interpolated_value_bar=
        cast_el_pt->interpolated_x(Time_level_for_projection,
                                   other_s,Projected_coordinate);
      }
     
     else
      {
       //Interpolation of current field
       interpolated_value_bar=
        cast_el_pt->get_field(Time_level_for_projection,fld,other_s);
      }
     
     //Loop over dofs of field fld
     for(unsigned l=0;l<n_value;l++)
      {
       int local_eqn = local_equation(fld,l);      
       if(local_eqn >= 0)
        {
         //calculate residuals
         residuals[local_eqn] += 
          (interpolated_value_proj - interpolated_value_bar)*psi[l]*W;
         
         //Calculate the jacobian
         if(flag==1)
          {
           for(unsigned l2=0;l2<n_value;l2++)
            {
             int local_unknown = local_equation(fld,l2);
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown) 
                += psi[l2]*psi[l]*W;  
              }
            }
          } //end of jacobian
         
        }
      }
    }//End of loop over ipt
   
  }//End of residual_for_projection function


 /// \short Use Eulerian coordinates for matching in locate_zeta
 /// when doing projection
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                   const unsigned &i) const 
 {                                                                        
  if (Do_projection)
   {
    return nodal_position_gen(n,k,i);                             
   }
  else
   {
    return ELEMENT::zeta_nodal(n,k,i);
   }
  }     



 /// \short Backup the element's state and
 /// switch it to projection mode.
 void enable_projection()
  {
   //Backup number of interaction
   Backup_ninteraction= ninteraction();

   //Backup flag for inclusion of geometric data
   if (add_external_geometric_data())
    {
     Backup_external_geometric_data=true;
    }
   else
    {
     Backup_external_geometric_data=false;
    }
   
   //Backup flag for inclusion of interaction data
   if (add_external_interaction_data())
    {
     Backup_external_interaction_data=true;
    }
   else
    {
     Backup_external_interaction_data=false;
    }

   //Actions to enable projection
   Do_projection=true; 
   ignore_external_geometric_data();
   ignore_external_interaction_data();
   set_ninteraction(1);
  }

 ///\short Helper function to restore the element to the state
 /// it was in before we entered the projection mode and switch off
 /// projection mode.
 void disable_projection()
  {
   //Restore number of interaction
   set_ninteraction(Backup_ninteraction);

   //Restore geometric data
   if (Backup_external_geometric_data)
    {
     include_external_geometric_data();
    }
   else
    {
     ignore_external_geometric_data();
    }
   
   //Restore interaction data
   if (Backup_external_interaction_data)
    {
     include_external_interaction_data();
    }
   else
    {
     ignore_external_interaction_data();
    }

   Do_projection=false;
  }


 ///Project (history values of) coordintes (true) or values (false)
 bool& project_coordinates()
  {return Project_coordinates;}


///Field that is currently being projected
 unsigned& projected_field()
  {return Projected_field;}
 
///Which history value are we projecting?
 unsigned& time_level_for_projection()
  {return Time_level_for_projection;}

 /// \short When projecting the history values of the nodal coordinates,
 /// this is the coordinate we're projecting
 unsigned& projected_coordinate()
  {return Projected_coordinate;}


};//End of class





//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
template<class ELEMENT>
class FaceGeometry<ProjectableElement<ELEMENT> > 
: public virtual FaceGeometry<ELEMENT>
 {
   public:
   FaceGeometry() : FaceGeometry<ELEMENT>() {}
 };



//==================================================================
/// Projection problem
//==================================================================
template<class PROJECTABLE_ELEMENT>
class ProjectionProblem : public virtual Problem
{
  public:
  
 ///Default constructor 
 ProjectionProblem() {}


///\short Project from base into the problem's own mesh. hierher: Need to
/// extend this to the case of multiple sub-meshes.
 void project(Mesh* base_mesh_pt)
  {

   //Display stats 
   unsigned n_element = Problem::mesh_pt()->nelement();
   unsigned n_element1=base_mesh_pt->nelement();
   unsigned n_node = Problem::mesh_pt()->nnode();
   unsigned n_node1=base_mesh_pt->nnode();
   oomph_info <<"\n=============================\n";
   oomph_info << "Base mesh has " << n_element1 << " elements\n";
   oomph_info << "Target mesh has " << n_element << " elements\n";
   oomph_info << "Base mesh has " << n_node1 << " nodes\n";
   oomph_info << "Target mesh has " << n_node << " nodes\n";
   oomph_info <<"=============================\n\n";
   
   if (n_element==0)
    {
     oomph_info 
      << "Very odd -- no elements in target mesh; "
      << " not doing anything in ProjectionProblem::project()\n";
     return;
    }

   // How many fields do we have to project?
   unsigned n_fields = dynamic_cast<PROJECTABLE_ELEMENT*>
    (Problem::mesh_pt()->element_pt(0))->nfields_for_projection();

   // Spatial dimension of the problem
   unsigned n_dim = Problem::mesh_pt()->node_pt(0)->ndim();

   // Default number of history values
   unsigned n_history_values=0;

   // Store pinned status and values of current values
   this->store_pinned_status_and_values();

   // Set things up for coordinate projection
   for(unsigned e=0;e<n_element;e++)
    {  
     PROJECTABLE_ELEMENT * el_pt = 
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
     
     // Switch to projection
     el_pt->enable_projection();
     
     //We first project history values of coordinates
     el_pt->project_coordinates()=true;
    }


   // Switch elements in base mesh to projection mode (required
   // to switch to use of Eulerian coordinates when identifying
   // corresponding points in the two meshes)
   for(unsigned e=0;e<n_element1;e++)
    {  
     PROJECTABLE_ELEMENT * el_pt = 
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (base_mesh_pt->element_pt(e));
     
     // Switch to projection
     el_pt->enable_projection();
    }
   
   // Set up multi domain interactions so we can locate the
   // values in the base mesh.
   // Note that it's important to switch elements to projection
   // mode first to ensure that matching is done based on Eulerian
   // rather than Lagrangian coordinates if pseudo-solid elements
   // are used.
   Multi_domain_functions::setup_multi_domain_interaction
    <PROJECTABLE_ELEMENT>(this,Problem::mesh_pt(),base_mesh_pt);

   // Pin all solid position dofs (if required)
   pin_solid_positions();

   // Prepare for projection in value 0: hierher Note that this 
   // assumes that field 0 can actually hold the coordinates
   // i.e. that it's interpolated isoparametrically!
   this->set_current_field_for_projection(0);
   this->unpin_dofs_of_field(0);
   this->pin_dofs_of_other_fields(0);

   //Check number of history values for coordinates
   n_history_values = dynamic_cast<PROJECTABLE_ELEMENT*>
    (Problem::mesh_pt()->element_pt(0))->
    nhistory_values_for_coordinate_projection();
     
   //Projection the coordinates
   for (unsigned i=0;i<n_dim;i++)
    {     
     oomph_info <<"\n\n=============================================\n";
     oomph_info <<    "Projecting history values for coordinate " << i
                << std::endl;
     oomph_info <<"=============================================\n\n";

     // Setup projection for i-th coordinate
     for(unsigned e=0;e<n_element;e++)
      {
       PROJECTABLE_ELEMENT * new_el_pt =
        dynamic_cast<PROJECTABLE_ELEMENT*>
        (Problem::mesh_pt()->element_pt(e));
       
       //Set current field
       new_el_pt->projected_coordinate()=i;
      }
     
     // Loop over number of history values, beginning with the latest one.
     // Don't deal with current time.
     for (unsigned h_tim=n_history_values;h_tim>1;h_tim--)
      {
       unsigned time_level=h_tim-1;
       
       //Set time_level we are dealing with
       this->set_time_level_for_projection(time_level);
       
       //Assign equation number
       oomph_info << "Number of equations for projection of coordinate " 
                  << i << " at time level " << time_level
                  << " : "<< assign_eqn_numbers() <<std::endl << std::endl;
       
       //Projection and interpolation
       Problem::newton_solve();
       
       //Move values back into history value of coordinate
       unsigned n_element = Problem::mesh_pt()->nelement();
       for(unsigned e=0;e<n_element;e++)
        {
         PROJECTABLE_ELEMENT * new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>
          (Problem::mesh_pt()->element_pt(e));
         
         Vector<std::pair<Data*,unsigned> >
          data=new_el_pt->data_values_of_field(0);
         
         unsigned d_size=data.size();
         for(unsigned d=0;d<d_size;d++)
          {
           //Replace as coordinates
           double coord=data[d].first->value(0,0);
           dynamic_cast<Node*>(data[d].first)->x(time_level,i) = coord;
          }
        }
      }
    }
   
   //Restore field 0 to current values
   this->restore_current_values_of_field(0);

   //Disable projection of coordinates
   for(unsigned e=0;e<n_element;e++)
    {  
     PROJECTABLE_ELEMENT * el_pt = 
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
     
     el_pt->project_coordinates()=false;
    }
   
   //Loop over fields 
   for (unsigned fld=0; fld<n_fields ;fld++)
    {
     //Do actions for this field
     this->set_current_field_for_projection(fld);
     this->restore_pin_status_of_field(fld);
     this->pin_dofs_of_other_fields(fld);
      
     //Check number of history values
     n_history_values = dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(0))->nhistory_values_for_projection(fld);
               
     //Loop over number of history values
     //Beginning with the latest one
     for (unsigned h_tim=n_history_values;h_tim>0;h_tim--)
      {
       unsigned time_level=h_tim-1;
       oomph_info <<"\n=========================================\n";
       oomph_info <<   "Projecting field " << fld << " at time level "
                  << time_level<<std::endl;
       oomph_info <<   "========================================\n";
              
       //Set time_level we are dealing with
       this->set_time_level_for_projection(time_level);
       
/*        std::cout<<"NOT APPLYING ANY BOUNDARY CONDITIONS DURING PROJECTION\n"; */
 /*       // Apply boundary conditions for the actual field */
/*        if(time_level==0) */
/*         { */
/*          // Undo the changes to boundary conditions that were */
/*          // applied to deal with the history values (for which */
/*          // we don't apply bcs). Remember that we do these first! */
/*          if (n_history_values>1) */
/*           { */
/*            //Restore pinned status as it has been modified */
/*            //for other history values */
/*            restore_pin_status_of_field(fld); */
/*            pin_dofs_of_other_fields(fld); */
/*           } */
/*          //Restore values */
/*          restore_current_values_of_field(fld); */
/*         } */
/*        // Don't apply any BCs for the history values */
/*        else */
        {
         //Unpin dofs for history values
         this->unpin_dofs_of_field(fld);
        }

       //Assign equation number
       oomph_info << "Number of equations for projection of field " 
                  << fld << " at time level " << time_level
                  << " : "<< assign_eqn_numbers() <<std::endl << std::endl;
       
       //Projection and interpolation
       Problem::newton_solve();
       
       // Move computed values into the required time-level (not needed
       // for  current values which are done last -- they simply
       // stay where they are)
       if(time_level!=0)
        {
         for(unsigned e=0;e<n_element;e++)
          {
           PROJECTABLE_ELEMENT * new_el_pt =
            dynamic_cast<PROJECTABLE_ELEMENT*>
            (Problem::mesh_pt()->element_pt(e));
           
           Vector<std::pair<Data*,unsigned> >
            data=new_el_pt->data_values_of_field(fld);
           
           unsigned d_size=data.size();
           for(unsigned d=0;d<d_size;d++)
            {
             //Move into time level 
             double c_value=data[d].first->value(0,data[d].second);
             data[d].first->set_value(time_level,data[d].second,c_value);
            }
          }
        }
       
      } //End of loop over time levels
     
    } //End of loop over fields
   
   
   //Reset parameters of external storage and interactions
   for(unsigned e=0;e<n_element;e++)
    {
     PROJECTABLE_ELEMENT * new_el_pt =
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
     
     new_el_pt->disable_projection();
    }
   
   for(unsigned e=0;e<n_element1;e++)
    {  
     PROJECTABLE_ELEMENT * el_pt = 
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (base_mesh_pt->element_pt(e));
     
     // Switch to projection
     el_pt->disable_projection();
    }

   //Restore pinned status for all fields
   for (unsigned f=0;f<n_fields;f++)
    {
     this->restore_pin_status_of_field(f);
    }
   restore_solid_pin_status();

   
   // Now cleanup the storage
   Data_backup.clear();
   Solid_pin_backup.clear();

   oomph_info << "Number of unknowns after project: " 
              << this->assign_eqn_numbers() << std::endl;

  } //End of function Projection

  
  private:


 /// \short Helper function to store pinned status and current values
 /// before doing projection
 void store_pinned_status_and_values()
 {

  // No need to do anything if there are no elements (in fact, we
  // probably never get here...)
  if (Problem::mesh_pt()->nelement()==0) return;

  // Get number of fields from first element, and create storage
  PROJECTABLE_ELEMENT * new_el_pt = 
   dynamic_cast<PROJECTABLE_ELEMENT*>
   (Problem::mesh_pt()->element_pt(0));
  unsigned n_fields = new_el_pt->nfields_for_projection();
  Data_backup.resize(n_fields);
  
  // Now extract data from all elements
  unsigned n_element = Problem::mesh_pt()->nelement();
  for(unsigned e=0;e<n_element;e++)
   {  
    PROJECTABLE_ELEMENT * new_el_pt = 
     dynamic_cast<PROJECTABLE_ELEMENT*>
     (Problem::mesh_pt()->element_pt(e));
    
    //Backup pin and values for all fields
    for (unsigned f=0;f<n_fields;f++)
     {
      Vector<std::pair<Data*,unsigned> >
       data=new_el_pt->data_values_of_field(f);

      // Backup
      unsigned d_size=data.size();
      for(unsigned d=0;d<d_size;d++)
       {
        // Package up in Backup item
        ProjectionHelper::BackupValue backup_item;
        backup_item.Data_pt=data[d].first;
        backup_item.Index=data[d].second;
        backup_item.Pinned=backup_item.Data_pt->is_pinned(backup_item.Index);
        backup_item.Value=backup_item.Data_pt->value(0,backup_item.Index);

        // Make backup
        Data_backup[f].insert(backup_item);
       }
     }
   }


  // Deal with positional dofs if (pseudo-)solid element
  SolidFiniteElement* solid_el_pt = dynamic_cast<SolidFiniteElement*>
   (Problem::mesh_pt()->element_pt(0));
  if (solid_el_pt!=0)
   {
    unsigned nnod=this->mesh_pt()->nnode();
    Solid_pin_backup.resize(nnod);
    for (unsigned j=0;j<nnod;j++)
     {      
      SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(
       this->mesh_pt()->node_pt(j));
      unsigned dim=solid_nod_pt->ndim();
      Solid_pin_backup[j].resize(dim);
      for (unsigned i=0;i<dim;i++)
       {
        Solid_pin_backup[j][i]=solid_nod_pt->position_is_pinned(i);
       }
     }
   }
 }
 


 /// \short Pin solid positions (if required)
 void pin_solid_positions()
 {
  // No need to do anything if there are no elements (in fact, we
  // probably never get here...)
  if (Problem::mesh_pt()->nelement()==0) return;

  /// Do we have a solid mesh?
  SolidFiniteElement* solid_el_pt = dynamic_cast<SolidFiniteElement*>
   (Problem::mesh_pt()->element_pt(0));
  if (solid_el_pt!=0)
   {
    unsigned nnod=this->mesh_pt()->nnode();
    for (unsigned j=0;j<nnod;j++)
     {      
      SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(
       this->mesh_pt()->node_pt(j));
      unsigned dim=solid_nod_pt->ndim();
      for (unsigned i=0;i<dim;i++)
       {
        solid_nod_pt->pin_position(i);
       }
     }
   }
 }

 
 /// \short Restore pin status for solid positions (if required)
 void restore_solid_pin_status()
 {
  // No need to do anything if there are no elements (in fact, we
  // probably never get here...)
  if (Problem::mesh_pt()->nelement()==0) return;
  
  /// Do we have a solid mesh?
  SolidFiniteElement* solid_el_pt = dynamic_cast<SolidFiniteElement*>
   (Problem::mesh_pt()->element_pt(0));
  if (solid_el_pt!=0)
   {
    unsigned nnod=this->mesh_pt()->nnode();
    for (unsigned j=0;j<nnod;j++)
     {      
      SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(
       this->mesh_pt()->node_pt(j));
      unsigned dim=solid_nod_pt->ndim();
      for (unsigned i=0;i<dim;i++)
       {
        if (Solid_pin_backup[j][i])
         {
          solid_nod_pt->pin_position(i);
         }
        else
         {
          solid_nod_pt->unpin_position(i);
         }
       }
     }
   }
 }

 /// \short Helper function to restore pinned status of field fld; typically 
 /// called after projection for the field
 void restore_pin_status_of_field(const unsigned &fld)
 {
  std::cout << "Number of unique backed up items for field " << fld << ": "
       << Data_backup[fld].size() << std::endl;
  
  for (std::set<ProjectionHelper::BackupValue>::iterator it=
        Data_backup[fld].begin();it!=Data_backup[fld].end();it++)
   {
    // Recover Backup item
    ProjectionHelper::BackupValue backup=*it;
    
    // Extract Data
    Data* data_pt=backup.Data_pt;
    
    // Extract index
    unsigned index=backup.Index;
    
    // Restore
    if (backup.Pinned)
     {
      data_pt->pin(index);
     }
    else
     {
      data_pt->unpin(index);
     }     
   }
 }


 /// \short Helper function to pin the values of other fields than fld
 void pin_dofs_of_other_fields(const unsigned &fld)
  {
   unsigned n_element = Problem::mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     PROJECTABLE_ELEMENT * new_el_pt =
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));

     //Pin dofs of other fields
     unsigned n_fields = new_el_pt->nfields_for_projection();
     for (unsigned f=0;f<n_fields;f++)
      {
       //Deal with field != fld
       if (f!=fld)
        {
         Vector<std::pair<Data*,unsigned> >
          data=new_el_pt->data_values_of_field(f);         
         unsigned d_size=data.size();         
         for(unsigned d=0;d<d_size;d++)
          {
           data[d].first->pin(data[d].second);
          }
        }
      }
    }
  }



 ///Helper function to unpin dofs of fld-th field
 void unpin_dofs_of_field(const unsigned &fld)
  {
   unsigned n_element = Problem::mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     PROJECTABLE_ELEMENT * new_el_pt =
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
     
     Vector<std::pair<Data*,unsigned> >
      data=new_el_pt->data_values_of_field(fld);
     
     unsigned d_size=data.size();     
     for(unsigned d=0;d<d_size;d++)
      {
       data[d].first->unpin(data[d].second);
      }   
    }
  }
  
 /// Helper function to restore current values of fld-th field
 void restore_current_values_of_field(const unsigned &fld)
 {
  for (std::set<ProjectionHelper::BackupValue>::iterator it=
        Data_backup[fld].begin();it!=Data_backup[fld].end();it++)
   {
    // Recover Backup item
    ProjectionHelper::BackupValue backup=*it;
    
    // Recover Data
    Data* data_pt=backup.Data_pt;
    
    // Recover index
    unsigned index=backup.Index;
    
    // Recover value
    double value=backup.Value;
    
    // Assign
    data_pt->set_value(index,value);
   }
 }




 /// Helper function to set time level for projection
 void set_time_level_for_projection(const unsigned &time_level)
  {
   unsigned n_element = Problem::mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     PROJECTABLE_ELEMENT * el_pt =
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
     
     //Set what time we are dealing with
     el_pt->time_level_for_projection()=time_level;
    }
  }
 

 ///Set current field for projection
 void set_current_field_for_projection(const unsigned &fld)
  {
   unsigned n_element = Problem::mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     PROJECTABLE_ELEMENT * new_el_pt =
      dynamic_cast<PROJECTABLE_ELEMENT*>
      (Problem::mesh_pt()->element_pt(e));
       
     //Set current field
     new_el_pt->projected_field()=fld;
    }
  }

  private:

 /// Backup for current values and pin status
 Vector<std::set<ProjectionHelper::BackupValue,
  ProjectionHelper::BackupValueLessThan> > Data_backup;

 /// Backup for pin status of solid node's position Data
 Vector<std::vector<bool> > Solid_pin_backup;

};


} //End of namespace

#endif

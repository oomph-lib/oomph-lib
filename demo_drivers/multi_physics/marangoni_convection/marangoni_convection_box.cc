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
/// Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations, and also
//implements temperature-dependent surface tension to give
//Marangoni-Bernard convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements and the multi_physics (Boussinesq) elements
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"
#include "multi_physics.h"
#include "fluid_interface.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/single_layer_spine_mesh.h"


//===================================================================
/// Function-type-object to perform comparison of complex data types
/// Needed to sort the complex eigenvalues into order based on the
/// size of the real part.
//==================================================================
template <class T>
class ComplexGreater
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(const std::complex<T> &x, const std::complex<T> &y) const
  {
   return x.real() > y.real();
  }
};


// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

namespace oomph
{
//==========================================================================
/// A \b preliminary class that is used to implement the constraint 
/// that the fluid volume must take a specific value for a Spine-based 
/// formulation. This is required in steady free-surface/interface
/// problems, where, in general, a family of solutions may exist.
/// \n\n
/// Note that as written, the function assumes that the mesh contains
/// three-straight boundaries at (x(0) = 0, x(1) = 0, x[0] = L) and
/// one free surface. In order to implement the constraint in total generality,
/// one would need to introduce surface elements all round the mesh
/// and integrate x.n along them. This is not necessary in the special case
/// described above, because the integrals along the three straight edges
/// are easily evaluated. The first two are zero and the third is
/// included in the residuals calculated in this element.
//=========================================================================
template<class ELEMENT>
class SpineVolumeConstraintPointElement : 
 public SpinePointFluidInterfaceBoundingElement<ELEMENT>
{
  private:
 
 /// Pointer to the desired value of the volume
 double *Volume_pt;

 /// Pointer to the Data item that stores the pressure that has been
 /// "traded" for the volume constraint in its single value.
 Data* Ptraded_data_pt;

 /// The Data that contains the traded pressure is stored
 /// as external Data for the element. Which external Data item is it?
 unsigned External_data_number_of_traded_pressure;

 /// The local eqn number for the traded pressure, which is 
 /// the variable that corresponds to this equation
 inline int ptraded_local_eqn()
  {
   if(Ptraded_data_pt==0) {return -1;}
   else 
    {return 
      this->external_local_eqn(External_data_number_of_traded_pressure,0);}
  }
 
  public:

 /// Constructor, there are no internal values. The pointer to the 
 /// element's (single) spine has to be set manually "from the outside"
 SpineVolumeConstraintPointElement() : 
  SpinePointFluidInterfaceBoundingElement<ELEMENT>()
  {
   // Initialise pointer to prescribed volume of fluid
   Volume_pt=0;
   // Initialise pointer to "traded" pressure Data.
   Ptraded_data_pt=0;
  } 

 /// Access function to the prescribed volume fluid 
 double* &volume_pt() {return Volume_pt;}

 /// Custom overload the additional volume constraint
 void add_additional_residual_contributions_interface_boundary(
   Vector<double> &residuals, 
   DenseMatrix<double> &jacobian,
   const unsigned &flag,
   const Shape &psif,
   const DShape &dpsifds,
   const Vector<double> &interpolated_n, 
   const double &W)
  {
   //If we have an external pressure, add the final term
   //to the volumetric constraint equation
   int local_eqn = ptraded_local_eqn();
   if(local_eqn >= 0)
    {
     //The integral of x.n on the RHS of the boundary is just x(0)*x(1)
     //of the top right-hand corner, we divide by two because we are working
     //in two dimensions and then subtract the desired volume
     residuals[local_eqn] = 
      0.5*this->node_pt(0)->x(0)*this->node_pt(0)->x(1) - *Volume_pt;
    }
  }

 /// Set the Data that contains the single pressure value
 /// that is "traded" for the volume constraint.
 /// The Data item must only contain a single value!
 void set_traded_pressure_data(Data* traded_pressure_data_pt)
  {
#ifdef PARANOID
   if (traded_pressure_data_pt->nvalue()!=1)
    {
     std::ostringstream error_stream;
     error_stream 
      << "The `traded` pressure Data must only contain a single value\n"
      << "This one contains " << traded_pressure_data_pt->nvalue() 
      << std::endl;
     
     std::string function_name = 
      "SpineVolumConstraintPointElement::\n";
     function_name += "set_traded_pressure_data()";
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
     // Store pointer explicitly
     Ptraded_data_pt=traded_pressure_data_pt;
     // Add to the element's external data so it gets included
     // in the black-box local equation numbering scheme
     External_data_number_of_traded_pressure = 
      this->add_external_data(traded_pressure_data_pt);
  }

}; 



//==================================================================
/// New Marangoni surface tension elements that also include a fixed
/// volume constraint for steady problems. 
/// Note that the implementation is rather specific for our meshes used
/// in this demo driver. You may have to alter this class for different
/// problems.
//=================================================================
template<class ELEMENT>
class FixedVolumeSpineLineMarangoniFluidInterfaceElement :
public SpineLineFluidInterfaceElement<ELEMENT>
{
private:

/// Integer to hold the local equation number for the single pressure
 /// value that has been traded for the volume constraint.
 int Ptraded_local_eqn;

 /// The Data that contains the traded pressure is stored
 /// as external Data for the element. Which external Data item is it?
 unsigned External_data_number_of_traded_pressure;

 /// Pointer to the Data item that stores the pressure that has been
 /// "traded" for the volume constraint in its single value.
 Data* Ptraded_data_pt;

 /// The local eqn number for the traded pressure, which is 
 /// the variable that corresponds to this equation
 inline int ptraded_local_eqn()
  {
   if(Ptraded_data_pt==0) {return -1;}
   else 
    {return this->external_local_eqn(
     External_data_number_of_traded_pressure,0);}
  }

 /// Pointer to a Biot number
 double *Bi_pt;

 /// Pointer to a Marangoni number
 double *Ma_pt;

 /// Index at which the temperature is stored at the nodes
 unsigned T_index;

 /// Default value of the physical constants
 static double Default_Physical_Constant_Value;

protected:

 /// Overload the surface tension function
 double sigma(const Vector<double> &s)
  {
   //Find the number of shape functions
   const unsigned n_node = this->nnode();
   //Now get the shape fuctions at the local coordinate
   Shape psi(n_node);
   this->shape(s,psi);
   
   double T = 0.0;
   //Now interpolate the temperature
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,T_index)*psi(l);
    }

   //Get the Marangoni and Capillary numbers
   double Ma = this->ma();
   double Ca = this->ca();
   //Return the variable surface tension
   return (1.0 - Ca*Ma*T);
  }

 /// Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
   {
    //Use finite differences to handle temperature variations
    const unsigned n_node = this->nnode();
    //Find the number of dofs in the element
    const unsigned n_dof = this->ndof();
    //Create newres vector
    Vector<double> newres(n_dof);
    
    //Integer storage for local unknown
    int local_unknown=0;
    
    //Use the default finite difference step
    const double fd_step = this->Default_fd_jacobian_step;
    
    //Loop over the nodes
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned t_index = this->T_index;

      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,t_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(t_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;
       
        //Calculate the new residuals
        this->get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
   }

   //Fill in the external data by finite differences
   this->fill_in_jacobian_from_external_by_fd(jacobian);
   
   //Call the generic routine to handle the spine variables
   this->fill_in_jacobian_from_geometric_data(jacobian);
  }

 
 /// Overload the Helper function to calculate the residuals and 
 /// jacobian entries. This particular function ensures that the
 /// additional entries are calculated inside the integration loop
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif, const DShape &dpsifds,
  const DShape &dpsifdS,
  const DShape &dpsifdS_div,
  const Vector<double> &s,
  const Vector<double> &interpolated_x,
  const Vector<double> &interpolated_n,
  const double &W,
  const double &J)
  {
   //Find the index at which the temperature is stored
   unsigned t_index = this->T_index;
   
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Read out the Bi number
   const double Bi = this->bi();

   //Now calculate the temperature at this point
   //Assuming the same shape functions are used (which they are)
   double T = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     T += this->nodal_value(l,t_index)*psif(l);
    }

   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;

   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,t_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       residuals[local_eqn] -= Bi*T*psif(l)*W*J;

       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,t_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) -= Bi*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over nodes
   
   //Now do the volume constraint
   local_eqn = ptraded_local_eqn(); 
   if(local_eqn >=0)
    {
     //Find the x position
     Vector<double> interpolated_x(2,0.0);
     //Fill in the x position
     for(unsigned l=0;l<n_node;l++)
      {
       const double psi_local = psif(l);
       for(unsigned i=0;i<2;i++)
        {
         interpolated_x[i] += this->nodal_position(l,i)*psi_local;
        }
      }
     
     //Find the dot product
     double dot = 0.0;
     for(unsigned k=0;k<2;k++) {dot += interpolated_x[k]*interpolated_n[k];}
     residuals[local_eqn] += 0.5*dot*W*J;
    }
  } //End of additional jacobian terms
  
 
public:

 /// Constructor that passes the bulk element and the face index
 /// down to the underlying 
 FixedVolumeSpineLineMarangoniFluidInterfaceElement(
 FiniteElement* const &element_pt, const int &face_index) : 
  SpineLineFluidInterfaceElement<ELEMENT>
  (element_pt,face_index)
  {
   //Initialise the values
   Bi_pt = &Default_Physical_Constant_Value;
   Ma_pt = &Default_Physical_Constant_Value;

   Ptraded_data_pt = 0;
   
   //Cast the bulk element 
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   //Now find the index at which the temperature is stored from the 
   //advection-diffusion part of the bulk element
   T_index = cast_element_pt->u_index_adv_diff();
  }

 double bi() {return *Bi_pt;}
 
 double ma() {return *Ma_pt;}
 
 double* &ma_pt() {return Ma_pt;}
 
 double* &bi_pt() {return Bi_pt;}
 

 /// Set the Data that contains the single pressure value
 /// that is "traded" for the volume constraint.
 /// The Data item must only contain a single value!
 void set_traded_pressure_data(Data* traded_pressure_data_pt)
  {
#ifdef PARANOID
   if (traded_pressure_data_pt->nvalue()!=1)
    {
     std::ostringstream error_stream;
     error_stream 
      << "The `traded` pressure Data must only contain a single value\n"
      << "This one contains " << traded_pressure_data_pt->nvalue() 
      << std::endl;
     
     std::string function_name = 
      "FixedVolumeSpineLineFluidInterfaceElement::\n";
     function_name += "set_traded_pressure_data()";
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   // Store pointer explicitly
   Ptraded_data_pt=traded_pressure_data_pt;
   // Add the traded pressure to the external data.
   // If it's a duplicate it will be reported as such
   External_data_number_of_traded_pressure =
    this->add_external_data(traded_pressure_data_pt);
  }
 
 /// / Overload the making of the edge element to create out
 /// volume constraint edge element
 FluidInterfaceBoundingElement* make_bounding_element(const int &face_index)
  {
   //Create a temporary pointer to the appropriate FaceElement
   SpineVolumeConstraintPointElement<ELEMENT> *Temp_pt =
    new SpineVolumeConstraintPointElement<ELEMENT>;

   //Attach the geometrical information to the new element
   this->build_face_element(face_index,Temp_pt);
   
   //Set the index at which the unknowns are stored from the element
   Temp_pt->u_index_interface_boundary() = this->U_index_interface;

   //Set the value of the nbulk_value, the node is not resized
   //in this problem, so it will just be the actual nvalue
   Temp_pt->nbulk_value(0) = Temp_pt->node_pt(0)->nvalue();

   //Set of unique geometric data that is used to update the bulk,
   //but is not used to update the face
   std::set<Data*> unique_additional_geom_data;
   //Get all the geometric data for this (bulk) element
   this->assemble_set_of_all_geometric_data(unique_additional_geom_data);

   //Now assemble the set of geometric data for the face element
   std::set<Data*> unique_face_geom_data_pt;
   Temp_pt->assemble_set_of_all_geometric_data(unique_face_geom_data_pt);
   //Erase the face geometric data from the additional data
   for(std::set<Data*>::iterator it=unique_face_geom_data_pt.begin();
       it!=unique_face_geom_data_pt.end();++it)
    {unique_additional_geom_data.erase(*it);}

   //Finally add all unique additional data as external data
   for(std::set<Data*>::iterator it = unique_additional_geom_data.begin();
       it!= unique_additional_geom_data.end();++it)
    {
     Temp_pt->add_external_data(*it);
     }

   //Return the value of the pointer
   return Temp_pt;
  }

};

//Define the default physical value
template<class ELEMENT>
double FixedVolumeSpineLineMarangoniFluidInterfaceElement<ELEMENT>::
Default_Physical_Constant_Value = 1.0;

} //end of the oomph namespace

//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh = 0.0;//1800.0;


 /// Scaled Bond number (Bo/Ca)
 /// This is set to zero so that there
 /// are no gravitational effects
 double Scaled_Bond = 0.0;

 /// Biot number
 double Biot = 1.0;

 /// Marangoni number
 double Marangoni = 125.0;

 /// Capillary number
 double Capillary = 1.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);

 
 /// Set the value of Pi
 const double Pi=MathematicalConstants::Pi;

 /// The volume of the domain
 double Volume = 1.2;

 /// The contact angle
 double Angle= 0.5*Pi;

 /// The external pressure
 double Pext = 0.0;

 /// Function that specifies the wall unit normal
 void wall_unit_normal_left_fct(const Vector<double> &x, 
                                Vector<double> &normal)
 {
  normal[0]=-1.0;
  normal[1]= 0.0;

 }


 /// Function that specifies the wall unit normal
 void wall_unit_normal_right_fct(const Vector<double> &x, 
                                Vector<double> &normal)
 {
  normal[0]=1.0;
  normal[1]=0.0;
 }



} // end_of_namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT,class INTERFACE_ELEMENT> 
class ConvectionProblem : public Problem
{

public:

 /// Constructor
 ConvectionProblem();

 /// Destructor. Empty
 ~ConvectionProblem() {}

 /// Unpin things for timestepping
 void switch_boundary_conditions()
  {
   //Pin the external pressure
   External_pressure_data_pt->pin(0);
   //Release the internal pressure
   unfix_pressure(0,0);

   //Set the new number of unknowns
   std::cout << "Preparing to timestep: "
             << assign_eqn_numbers() << "\n";
  }


 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve() {}

 /// Remember to update the nodes!
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
   
   // This driver code cannot be allowed to use the analytical form of
   // get_dresidual_dnodal_coordinates(...) that is implemented in the
   // NavierStokesEquations class, since the elemental residuals have
   // contributions from external data which is not taken into account
   // by that routine. We therefore force the bulk elements to use the
   // fully-finite differenced version.
   const unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     ElementWithMovingNodes* el_pt =
      dynamic_cast<ElementWithMovingNodes*>(Bulk_mesh_pt->element_pt(e));
     el_pt->evaluate_shape_derivs_by_direct_fd();
    }
  }

 /// Actions before adapt:(empty)
 void actions_before_adapt(){}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
   set_boundary_conditions(time_pt()->time());
  }

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// UnFix pressure in element e at pressure dof pdof and set to pvalue
 void unfix_pressure(const unsigned &e, const unsigned &pdof)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    unfix_pressure(pdof);
  } // end_of_unfix_pressure


 /// Doc the solution.
 void doc_solution();

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Mesh for the bulk fluid elements
 SingleLayerSpineMesh<ELEMENT> *Bulk_mesh_pt;

 /// The mesh for the interface elements
 Mesh* Surface_mesh_pt;
 
 /// The mesh for the element at the contact points
 Mesh* Point_mesh_pt;
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

 /// Pointer to the external data point
 Data* External_pressure_data_pt;


}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::ConvectionProblem()
{
 //Allocate a timestepper
 add_time_stepper_pt(new BDF<2>);

 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=Global_Physical_Variables::Volume;

 // Domain length in y-direction
 double l_y=1.0;

 Bulk_mesh_pt =
  new SingleLayerSpineMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 for(unsigned i=0;i<n_x;i++)
  {
   //Construct a new 1D line element on the face on which the local
   //coordinate 1 is fixed at its max. value (1) --- This is face 2
   FiniteElement *interface_element_pt =
    new INTERFACE_ELEMENT(Bulk_mesh_pt->finite_element_pt(n_x*(n_y-1)+i),2);
   
   //Push it back onto the stack
   this->Surface_mesh_pt->add_element_pt(interface_element_pt); 
  }
 
 //Create the Point mesh that is responsible for enforcing the contact
 //angle condition
 Point_mesh_pt = new Mesh;
 {
  //Make the point (contact) element from the first surface element
  FiniteElement* point_element_pt = 
   dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(0))
                                    ->make_bounding_element(-1);
 
 //Add it to the mesh
 this->Point_mesh_pt->add_element_pt(point_element_pt);
 
 //Make the point (contact) elemnet from the last surface element
 point_element_pt = 
  dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(n_x-1))
  ->make_bounding_element(1);

//Add it to the mesh
this->Point_mesh_pt->add_element_pt(point_element_pt);
}

 //Create a Data object whose single value stores the
 //external pressure
 External_pressure_data_pt = new Data(1);
 
 // Set external pressure
 External_pressure_data_pt->set_value(0,Global_Physical_Variables::Pext);
 
 // Create a pointer to the (single value) Data item that
 // will contain the pressure value that we're
 // trading for the volume constraint
 Data* traded_pressure_data_pt;
 
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
 fix_pressure(0,0,0.0);
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the minimum index to be pinned (all values by default)
   unsigned val_min=0;
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=3;

   //The left is rigid
   if(ibound==3) {val_max=2;}
   //The right is rigid
   if(ibound==1) {val_max=2;}
   
   //The top wall is free
   if(ibound==2) {val_max=0;}

   //Loop over the number of nodes on the boundary
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=val_min;j<val_max;j++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set Re/Fr (Bo/Ca)
   el_pt->re_invfr_pt() = &Global_Physical_Variables::Scaled_Bond;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
  }


  // Loop over the interface elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   // Upcast from GeneralsedElement to the present element
   INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(
    Surface_mesh_pt->element_pt(i));
   
   // Set the Biot number
   el_pt->bi_pt() = &Global_Physical_Variables::Biot;

   // Set the Marangoni number
   el_pt->ma_pt() =&Global_Physical_Variables::Marangoni;

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Capillary;

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
 // to the volume constraint element. This is the point on the right
 //-----------------------------------------------------------------
 {
  SpineVolumeConstraintPointElement<ELEMENT>* el_pt =
   dynamic_cast<SpineVolumeConstraintPointElement<ELEMENT>*>
   (Point_mesh_pt->element_pt(1));
  
  el_pt->volume_pt() = &Global_Physical_Variables::Volume;
  el_pt->set_traded_pressure_data(traded_pressure_data_pt);
 }    


 // Set the contact angle boundary condition for the leftmost element
 // (pass pointer to double that specifies the contact angle)
{
 FluidInterfaceBoundingElement* left_element_pt = 
  dynamic_cast<FluidInterfaceBoundingElement*>(
   Point_mesh_pt->element_pt(0));
 
 left_element_pt->set_contact_angle(&Global_Physical_Variables::Angle);
// Set the Capillary number
 left_element_pt->ca_pt() = &Global_Physical_Variables::Capillary;
// Set the wall normal function
 left_element_pt->wall_unit_normal_fct_pt() = 
  &Global_Physical_Variables::wall_unit_normal_left_fct; 

 // Set the contact angle boundary condition for the rightmost element
 // (pass pointer to double that specifies the contact angle)
 FluidInterfaceBoundingElement* right_element_pt = 
 dynamic_cast<FluidInterfaceBoundingElement*>(
  Point_mesh_pt->element_pt(1));

 right_element_pt->
  set_contact_angle(&Global_Physical_Variables::Angle);
  // Set the Capillary number
 right_element_pt->ca_pt() = &Global_Physical_Variables::Capillary;
 // Set the wall normal function
 right_element_pt->wall_unit_normal_fct_pt() = 
    &Global_Physical_Variables::wall_unit_normal_right_fct; 
}

 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Surface_mesh_pt);
 this->add_sub_mesh(Point_mesh_pt);
 
 this->build_global_mesh();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //If we are on the top boundary, do not set the velocities
     //(yet)
     if(ibound==2) {vel_max = 0;}

     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the bottom boundary peturb the velocity
     if(ibound==0) 
      {
       //Add small velocity imperfection if desired
       double epsilon = 0.01;

       //Read out the x position
       double x = nod_pt->x(0);

       //Set a sinusoidal perturbation in the vertical velocity
       //This perturbation is mass conserving
       double value = sin(2.0*MathematicalConstants::Pi*x/3.0)*
        epsilon*time*exp(-time);
       nod_pt->set_value(1,value);
      }

            
       //If we are on the bottom boundary, set the temperature
       //to 2 (heated)
     if(ibound==0)
      {
       nod_pt->set_value(2,2.0);
      }
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void ConvectionProblem<ELEMENT,INTERFACE_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 //Calculate the angle (in degrees)
 double angle = 
  Global_Physical_Variables::Angle*180/Global_Physical_Variables::Pi;
 //The filename will have the angle appended here
 sprintf(filename,"%s/soln%i_%g.dat",Doc_info.directory().c_str(),
         Doc_info.number(),angle);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Construct our problem
 ConvectionProblem<SpineElement<BuoyantQCrouzeixRaviartElement<2> >,
  FixedVolumeSpineLineMarangoniFluidInterfaceElement<SpineElement<BuoyantQCrouzeixRaviartElement<2> > > > problem;

 //Set "interesting" Marangoni and Capillary numbers 
 Global_Physical_Variables::Marangoni = 200.0;
 Global_Physical_Variables::Capillary = 2.75e-3;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);
 
 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution();

 //Decrease the contact angle
 for(unsigned i=0;i<3;i++)
  {
   Global_Physical_Variables::Angle -= 5*Global_Physical_Variables::Pi/180;
   problem.steady_newton_solve();
   problem.doc_solution();
  }

 //Now we timestep the system
 //Note that we are not taking any contact line dynamics into
 //account, we are merely evolving between steady solutions
 {
  problem.switch_boundary_conditions();

  //Set the timestep
  double dt = 0.5;
  
  //Initialise the value of the timestep and set initial values 
  //of previous time levels assuming an impulsive start.
  problem.assign_initial_values_impulsive(dt);
  
  
  //Set the number of timesteps to our default value
  unsigned n_steps = 200;
  
  //If we have a command line argument, perform fewer steps 
  //(used for self-test runs)
  if(argc > 1) {n_steps = 2;}
  
  //Perform n_steps timesteps
  for(unsigned i=0;i<n_steps;++i)
   {
    problem.unsteady_newton_solve(dt);
    problem.doc_solution();
   }
 }

} // end of main










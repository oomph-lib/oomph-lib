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
//Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_FLUX_CONTROL_ELEMENTS
#define OOMPH_FLUX_CONTROL_ELEMENTS

// NOTE: This does not appeared to be used. So I have renmaed it from
// flux_control_elements.h to flux_control_elements_bk.h

namespace oomph
{


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//  CHECK WE DETECT OF HANGING NODES SOMEWHERE IN TRACTION MESH

// Forward declaration
template <class ELEMENT> class NavierStokesFluxControlElement;

//======================================================================
/// A class for the elements that applies a net fluid flux across a 
/// boundary by the imposition of an applied traction
/// to the Navier--Stokes equations 
//======================================================================
template <class ELEMENT>
class NetFluxControlElement : public virtual GeneralisedElement 
{
public:

 /// Constructor, which takes a mesh of the face elements which will
 /// impose the pressure to control the flux
 NetFluxControlElement(Mesh* flux_control_mesh_pt,
                       double* prescribed_outflow_value_pt) :
  Flux_control_mesh_pt(flux_control_mesh_pt),
  Prescribed_outflow_value_pt(prescribed_outflow_value_pt)
  {
   // Construct Pressure_data_pt
   Pressure_data_pt = new Data(1);
   
   // Add the new Data to internal Data for this element
   add_internal_data(Pressure_data_pt);

   // pin pressure data 
   // Pressure_data_pt->pin(0);
   
   // Loop over elements in the Flux_control_mesh to add
   // Data from the elements in the flux control
   // mesh to the external data for this element
   unsigned n_el =  Flux_control_mesh_pt->nelement();
   for (unsigned e=0; e<n_el; e++)
    {
     // Get pointer to the element
     FiniteElement * f_el_pt = 
      Flux_control_mesh_pt->finite_element_pt(e);
     
     // Loop over the nodes
     unsigned n_node = f_el_pt->nnode();
     for(unsigned n=0;n<n_node;n++)
      {
       add_external_data(f_el_pt->node_pt(n));
      }
    }
  }
 
 /// Empty Destructor - Data gets deleted automatically
 ~NetFluxControlElement() {}
 
 /// Broken copy constructor
 NetFluxControlElement(const NetFluxControlElement&) 
  { 
   BrokenCopy::broken_copy("NetFluxControlElement");
  } 
 
 /// Function return pointer to the Data object whose
 /// single value is the pressure applied by the elements in 
 /// Flux_control_mesh_pt
 Data* pressure_data_pt() const {return Pressure_data_pt;}


 /// Add the element's contribution to its residual vector:
 /// The flow constraint. [Note: Jacobian is computed 
 /// automatically by finite-differencing]
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  // Initialise volume flux
  double volume_flux  = 0.0;
  
  // Loop over elements in Flux_control_mesh_pt and calculate flux 
  unsigned n_el = Flux_control_mesh_pt->nelement();
  for (unsigned e=0; e<n_el; e++)
   {
    // Get a pointer to the element
    GeneralisedElement* el_pt = Flux_control_mesh_pt->element_pt(e);
    
    // Cast to NavierStokesFluxControlElement
    NavierStokesFluxControlElement<ELEMENT>* flux_control_el_pt = 
     dynamic_cast<NavierStokesFluxControlElement<ELEMENT>* >(el_pt);
   
    // Add the elemental volume flux
    volume_flux += flux_control_el_pt->get_volume_flux();
   }
  
  residuals[0] = volume_flux - *Prescribed_outflow_value_pt;
 }




 /// The number of "blocks" that degrees of freedom in this element
 /// are sub-divided into
 /// 
 /// IMPORTANT:
 /// This is not even the correct function name! Because this appears
 /// to be untested (it will break if tested!), I will comment this out.
 /// Please re-implement if required in the future.
 /// The correct function signature is:
 /// unsigned ndof_types() const
// unsigned nblock_types()
//  {
//   return 2;
//  }

 /// Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "block" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.) Pressure=1 since we expect to add this
 /// unknown to the pressure block
 ///
 /// IMPORTANT:
 /// This is not even the correct function name
 /// I have no self tests for this, because it appears to be not tested,
 /// I will comment this out, the correct function signature is:
 /// void get_dof_numbers_for_unknowns(
 /// std::list<std::pair<unsigned long, unsigned> >& block_lookup_list) const
// void get_block_numbers_for_unknowns(
//  std::list<std::pair<unsigned long, unsigned> >& block_lookup_list)
//  {
//   // pair to store block lookup prior to being added to list
//   std::pair<unsigned,unsigned> block_lookup;
// 
//   block_lookup.first = eqn_number(0);
//   block_lookup.second = 1;
//     
//   // add to list
//   block_lookup_list.push_front(block_lookup);
//  }
 
 
private:
 
 /// Data object whose single value is the pressure 
 /// applied by the elements in the Flux_control_mesh_pt
 Data* Pressure_data_pt;
 
 /// Mesh of elements which impose a pressure which controls
 /// the net flux
 Mesh* Flux_control_mesh_pt;
 
 /// Pointer to the value that stores the prescribed outflow
 double* Prescribed_outflow_value_pt;

};

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////



//======================================================================
/// A class for elements that allows a pressure to be imposed for volume
/// flux control by the the imposition of an applied traction
/// to the Navier--Stokes equations 
/// The geometrical information can be read from the FaceGeometery<ELEMENT> 
/// class and and thus, we can be generic enough without the need to have
/// a separate equations class
//======================================================================
template <class ELEMENT>
class NavierStokesFluxControlElement : 
 public virtual NavierStokesSurfacePowerElement<ELEMENT> 
{
public:
 
 /// Constructor, which takes a "bulk" element and the value of the index
 /// and its limit
 NavierStokesFluxControlElement(FiniteElement* const &element_pt, 
                                const int &face_index) : 
  NavierStokesSurfacePowerElement<ELEMENT>(element_pt, face_index)
  { 
#ifdef PARANOID
   {
    //Check that the element is not a refineable 3d element
    ELEMENT* elem_pt = new ELEMENT;
    //If it's three-d
    if(elem_pt->dim()==3)
     {
      //Is it refineable
      if(dynamic_cast<RefineableElement*>(elem_pt))
       {
        //Issue a warning
        OomphLibWarning(
         "This flux element will not work correctly if nodes are hanging\n",
         "NavierStokesFluxControlElement::Constructor",
         OOMPH_EXCEPTION_LOCATION);
       }
     }
   }
#endif
   
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
   
   //Set the dimension from the dimension of the first node
   Dim = this->node_pt(0)->ndim();
  }
 
 /// Destructor should not delete anything
 ~NavierStokesFluxControlElement() {}

 /// The number of "blocks" that degrees of freedom in this element
 /// are sub-divided into
 unsigned nblock_types()
  {
   return 2;
  }

 /// Create a list of pairs for all unknowns in this element,
 /// Do nothing since this element adds no new dofs
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long, unsigned> >& block_lookup_list)
  {}
 
 /// This function returns just the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function using a dummy matrix argument
   fill_in_generic_residual_contribution_fluid_traction(
    residuals,GeneralisedElement::Dummy_matrix);
  }
 
/*  ///This function returns the residuals and the jacobian */
/*  inline void fill_in_contribution_to_jacobian(Vector<double> &residuals, */
/*                                               DenseMatrix<double> &jacobian) */
/*   { */
/*    //Call the generic routine */
/*    fill_in_generic_residual_contribution_fluid_traction(residuals,jacobian); */
/*   } */
 
 /// Overload the output function
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}
 
 /// Output function: x,y,[z],u,v,[w],p in tecplot format
 void output(std::ostream &outfile, const unsigned &nplot)
  {FiniteElement::output(outfile,nplot);}
 
 /// Function to add to external data the Data object whose
 /// single value is the pressure applied by the element
 void add_pressure_data(Data* pressure_data_pt)
  {
   Pressure_data_id = this->add_external_data(pressure_data_pt);
  }
 
protected:
 
 
 /// Access function that returns the local equation numbers
 /// for velocity components.
 /// u_local_eqn(n,i) = local equation number or < 0 if pinned.
 /// The default is to asssume that n is the local node number
 /// and the i-th velocity component is the i-th unknown stored at the node.
 virtual inline int u_local_eqn(const unsigned &n, const unsigned &i)
  {return this->nodal_local_eqn(n,i);}
 
 /// Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping 
 inline double shape_and_test_at_knot(const unsigned &ipt, 
                                      Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = this->nnode();
   //Calculate the shape functions
   this->shape_at_knot(ipt,psi);
   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}
   //Return the value of the jacobian
   return this->J_eulerian_at_knot(ipt);
  }
 
 
 /// This function returns the residuals for the 
 /// traction function.
 /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
 void fill_in_generic_residual_contribution_fluid_traction(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian)
  {
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions
   Shape psif(n_node), testf(n_node);
   
   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();
   
   //Integers to store local equation numbers
   int local_eqn=0;
   
   // Get the pressure at the outflow
   double pressure = this->external_data_pt(Pressure_data_id)->value(0);
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Find the shape and test functions and return the Jacobian
     //of the mapping
     double J = shape_and_test_at_knot(ipt,psif,testf);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Need to find position to feed into Traction function
     Vector<double> interpolated_x(Dim);
     
     //Initialise to zero
     for(unsigned i=0;i<Dim;i++) {interpolated_x[i] = 0.0;}
     
     //Calculate velocities and derivatives
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over velocity components
       for(unsigned i=0;i<Dim;i++) {interpolated_x[i] += 
                                     this->nodal_position(l,i)*psif[l];}
      }
     
     // Get the outer unit normal
     Vector<double> unit_normal(Dim);
     this->outer_unit_normal(ipt, unit_normal);
     
     // Calculate the traction
     Vector<double> traction(Dim);
     for (unsigned i=0; i<Dim; i++)
      {
       traction[i] = pressure*unit_normal[i];
      }
     
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //Loop over the velocity components
       for(unsigned i=0;i<Dim;i++)
        {
         local_eqn = u_local_eqn(l,i);
         /*IF it's not a boundary condition*/
         if(local_eqn >= 0)
          {
           //Add the user-defined traction terms
           residuals[local_eqn] += traction[i]*testf[l]*W;
           
           //Assuming the the traction DOES NOT depend upon velocities
           //or pressures, the jacobian is always zero, so no jacobian
           //terms are required
           
          }
        } //End of loop over dimension
      } //End of loop over shape functions
   
    }
   
  }
 
private:
 
 /// Id of external Data object whose single value is the 
 /// pressure applied by the elements
 unsigned Pressure_data_id;
 
 /// The highest dimension of the problem 
 unsigned Dim;
 
 
}; 

#endif



}

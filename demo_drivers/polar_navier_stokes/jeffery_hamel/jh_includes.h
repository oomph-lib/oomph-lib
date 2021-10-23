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

//Header file for my FluxConstraint and jh_mesh classes

namespace oomph
{

//=======================================================================
/// Generalised element used to specify the mass flux (=1)
//=======================================================================
class FluxConstraint : public GeneralisedElement
{
 double Flux;

 protected:

 /// Pointer to the Data item that stores the external pressure
 Data* Pext_pt;

 /// The Data that contains the traded pressure is stored
 /// as external Data for the element. Which external Data item is it?
 unsigned External_data_number_of_Pext;

 public:
 /// Constructor there is one bit of internal data, the fixed flux
 FluxConstraint() 
  {
   //Specify flux
   Flux=1.0;

   //Set the external pressure pointer to be zero
   Pext_pt=0;
  }

 ///Function for setting up external pressure
 void set_pressure_data(Data* pext_pt)
  {
   //Set external pressure pointer
   Pext_pt=pext_pt;

   // Add to the element's external data so it gets included
   // in the black-box local equation numbering scheme
   External_data_number_of_Pext = 
    this->add_external_data(Pext_pt);
  }
 
 //Add the contribution to the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Storage for local eqn number
   int pext_local_eqn;

   if(Pext_pt==0)
    {
     pext_local_eqn=-1;
    }
   else
    {
     //If at a non-zero degree of freedom add in the entry
     pext_local_eqn = external_local_eqn(External_data_number_of_Pext,0);
    }

   if(pext_local_eqn >= 0) 
    {
     residuals[pext_local_eqn] = -Flux;
    }
  }
 
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the residuals
   fill_in_contribution_to_residuals(residuals);
   //No Jacobian terms :)
  }

 inline void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
							      DenseMatrix<double> &jacobian,
							      DenseMatrix<double> &mass_matrix)
  {
   //Call the residuals
   fill_in_contribution_to_residuals(residuals);
   //No Jacobian or mass matrix terms :)
  }

};

//======================================================================
///  My fluid mesh (with traction elements)
//======================================================================
template <class ELEMENT>
class jh_mesh : public virtual Refineable_r_mesh<ELEMENT>
{
protected:

 //Storage for pointers to my traction elements
 Vector<PolarNavierStokesTractionElement<ELEMENT>*> Inlet_traction_elt_pt;
 Vector<PolarNavierStokesTractionElement<ELEMENT>*> Outlet_traction_elt_pt;
 /// Generalised element to determine my mass flux (always 1)
 FluxConstraint* Flux_constraint_pt;
 //Storage for pointers to my shear stress integral elements
 Vector<PolarStressIntegralElement<ELEMENT>*> Lower_stress_integral_elt_pt;
 Vector<PolarStressIntegralElement<ELEMENT>*> Upper_stress_integral_elt_pt;

 //Internal storage for traction parameters

public:

 /// Return  pointer to inlet traction element e
 PolarNavierStokesTractionElement<ELEMENT>* inlet_traction_elt_pt(unsigned e) 
   {return Inlet_traction_elt_pt[e];}
 /// Return length of inlet traction element vector
 unsigned inlet_traction_elt_length() {return Inlet_traction_elt_pt.size();}
 /// Return  pointer to outlet traction element e
 PolarNavierStokesTractionElement<ELEMENT>* outlet_traction_elt_pt(unsigned e) 
   {return Outlet_traction_elt_pt[e];}
 /// Return length of outlet traction element vector
 unsigned outlet_traction_elt_length() {return Outlet_traction_elt_pt.size();}

 /// Return  pointer to Flux Cosntraint Element
 FluxConstraint* flux_constraint_elt_pt() 
   {return Flux_constraint_pt;}

 /// Return  pointer to shear integral element
 PolarStressIntegralElement<ELEMENT>* lower_stress_integral_elt_pt(unsigned long e) 
   {return Lower_stress_integral_elt_pt[e];}
 /// Return length of lower stress integral element vector
 unsigned lower_stress_integral_elt_length() {return Lower_stress_integral_elt_pt.size();}
 /// Return  pointer to shear integral element
 PolarStressIntegralElement<ELEMENT>* upper_stress_integral_elt_pt(unsigned long e) 
   {return Upper_stress_integral_elt_pt[e];}
 /// Return length of upper stress integral element vector
 unsigned upper_stress_integral_elt_length() {return Upper_stress_integral_elt_pt.size();}

 /// Constructor, which "builds" the mesh. The arguments are the number
 /// of elements in each direction.
 jh_mesh(const unsigned int &nx,const unsigned int &ny) : 
   Refineable_r_mesh<ELEMENT>(nx,ny)
  {  
   //Now bolt on traction stuff

   //Assign fluid elements to vector
   this->assign_fluid_element_vector();

   //Attach traction elements where needed
   if(Global_Physical_Variables::inlet_traction) make_traction_elements(false);
   if(Global_Physical_Variables::outlet_traction) make_traction_elements(true);
   if(Global_Physical_Variables::inlet_traction && Global_Physical_Variables::outlet_traction) make_flux_element();

   //Attach shear integral elements
   make_shear_elements();
 
  }

 //Function to add the traction boundary elements
 void make_traction_elements(const bool& outlet)
  {
   //Specify inlet/outlet specific quantities
   unsigned ibound; int index; double eta;
   if(outlet) {ibound=1;index=1;eta=Global_Physical_Variables::eta_outlet; }
   else {ibound=3;index=-1;eta=Global_Physical_Variables::eta_inlet; }

   unsigned num_elt = this->nboundary_element(ibound);

   //Loop over the number of elements on the boundary
   for(unsigned ielt=0;ielt<num_elt;ielt++)
    {
     PolarNavierStokesTractionElement<ELEMENT> *surface_element_pt =
       new PolarNavierStokesTractionElement<ELEMENT>
       (this->boundary_element_pt(ibound,ielt),index);
     //Push it back onto the Element_pt Vector
     this->Element_pt.push_back(surface_element_pt);

     if(outlet) { this->Outlet_traction_elt_pt.push_back(surface_element_pt); }
     else { this->Inlet_traction_elt_pt.push_back(surface_element_pt); }

     //Any other information to pass: 
     surface_element_pt->set_eta(eta);
     surface_element_pt->traction_fct_pt() = 
      &Global_Physical_Variables::traction_function;
     surface_element_pt->set_boundary(index);
     surface_element_pt->alpha_pt() = &Global_Physical_Variables::Alpha;
  
    }//End of loop over elements

   std::cout << std::endl << "Traction elements attached to mesh" << std::endl << std::endl;

  }//End of make traction elements

 //Function to add a generalised flux element
 void make_flux_element()
  {
   Flux_constraint_pt = new FluxConstraint();
   //Push it back onto the Element_pt Vector
   this->Element_pt.push_back(Flux_constraint_pt);
  }

 // Function to add shear stress integral elements                      
 void make_shear_elements()
  {
   //loop over boundaries
   for(unsigned ibound=0;ibound<4;ibound+=2)
    {
     int index;
     if(ibound==0) index=-2;
     else index=2;

     unsigned num_elt = this->nboundary_element(ibound);

     //Loop over the number of elements along boundary
     for(unsigned ielt=0;ielt<num_elt;ielt++)
      {
       //Element on lower boundary
       PolarStressIntegralElement<ELEMENT> *surface_element_pt =
       new PolarStressIntegralElement<ELEMENT>
       (this->boundary_element_pt(ibound,ielt),index);

       //Push it back onto the appropriate vector
       if(ibound==0) 
        this->Lower_stress_integral_elt_pt.push_back(surface_element_pt);
       else this->Upper_stress_integral_elt_pt.push_back(surface_element_pt);

      }//End of loop over elements

    }//End of loop over boundaries

   std::cout << std::endl << "Shear elements attached to mesh" << std::endl << std::endl;

  }//End of make shear elements

 //Function to remove the traction boundary elements
 void remove_traction_elements()
  { 
   //Find the number of traction elements
   unsigned Ntraction = this->inlet_traction_elt_length();
   Ntraction += this->outlet_traction_elt_length();

   //If we have traction elements at both ends then we have
   //one extra element to remove (the flux constraint element)
   if(Global_Physical_Variables::inlet_traction && Global_Physical_Variables::outlet_traction) 
    {Ntraction+=1; this->Flux_constraint_pt=0;}

   //The traction elements are ALWAYS? stored at the end
   //So delete and remove them
   for(unsigned e=0;e<Ntraction;e++)
    {
     delete this->Element_pt.back();
     this->Element_pt.pop_back();
    }

   //Now clear the vectors of pointers to traction elements
   Inlet_traction_elt_pt.clear();
   Outlet_traction_elt_pt.clear();

   std::cout << std::endl << "Traction elements removed from mesh" << std::endl << std::endl;

  }//End of remove_traction_elements

 //Function to remove the traction boundary elements
 void remove_shear_elements()
  { 
   //Find the number of shear elements
   unsigned Nshear_lower = this->lower_stress_integral_elt_length();
   unsigned Nshear_upper = this->upper_stress_integral_elt_length();

   //So delete and remove shear elements
   for(unsigned e=0;e<Nshear_lower;e++)
    {
     delete this->Lower_stress_integral_elt_pt.back();
     this->Lower_stress_integral_elt_pt.pop_back();
    }
   //So delete and remove shear elements
   for(unsigned e=0;e<Nshear_upper;e++)
    {
     delete this->Upper_stress_integral_elt_pt.back();
     this->Upper_stress_integral_elt_pt.pop_back();
    }

   //Now clear the vectors of pointers to traction elements
   Lower_stress_integral_elt_pt.clear();
   Upper_stress_integral_elt_pt.clear();

   std::cout << std::endl << "Shear elements removed from mesh" << std::endl << std::endl;

  }//End of remove_shear_elements

  //Function to put current fluid elements into a vector of their own
  void assign_fluid_element_vector()
  {
   unsigned check=inlet_traction_elt_length();
   check+=outlet_traction_elt_length();

   this->Fluid_elt_pt.clear();

   if(check!=0)
    {
     std::cout << "Warning, attempting to assemble fluid element vector ";
     std::cout << "whilst traction elements are attached to the mesh" << std::endl;
    }
   else
    {
     unsigned n_elt = this->Element_pt.size();
     for(unsigned e=0;e<n_elt;e++)
      {
       this->Fluid_elt_pt.push_back(this->Element_pt[e]);
      }
    }
  }//End of assign_fluid_element_vector

}; //End of jh_mesh class


} //End of namespace oomph


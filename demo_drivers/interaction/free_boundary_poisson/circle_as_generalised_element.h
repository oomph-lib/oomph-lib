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
#ifndef CIRCLE_AS_GEN_ELEMENT_HEADER
#define CIRCLE_AS_GEN_ELEMENT_HEADER

#include "circle.h"


namespace oomph
{

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//===========start_of_general_circle=========================================
///  GeneralCircle "upgraded" to a GeneralisedElement: Circular 
/// ring whose position is given by
/// \f[ x = X_c + R \cos(\zeta)  \f]
/// \f[ y = Y_c + R \sin(\zeta)  \f]
/// The ring's vertical position \f$ Y_c \f$ is
/// determined by "pseudo elasticity":
/// \f[
/// 0 = f_{load} - Y_c \ k_{stiff} 
/// \f]
/// This simulates the case where the centre of the ring is mounted on 
/// an elastic spring of stiffness \f$ k_{stiff} \f$ and loaded by 
/// the force \f$ f_{load}. \f$ The "load" is specified by the 
/// Data object \c load_pt(). 
//=========================================================================
class ElasticallySupportedRingElement : public GeneralisedElement, 
                                        public GeneralCircle
{

public:
   
 ///  Constructor: Build  ring from doubles that describe 
 /// the geometry: x and y positions of centre and the radius.
 /// Initialise stiffness to 1.0. By default, no load is set.
 ElasticallySupportedRingElement(const double& x_c, const double& y_c, 
                                 const double& r) : 
  GeneralCircle(x_c,y_c,r), K_stiff(1.0), Load_data_has_been_set(false)
  {
   // The geometric data is internal to the element -- we copy the pointers
   // to the GeomObject's geometric data to the element's internal 
   // data to ensure that any unknown values of geometric data are 
   // given global equation numbers. The add_internal_data(...)
   // function returns the index by which the added Data item
   // is accessible from internal_data_pt(...). 
   Internal_geometric_data_index=add_internal_data(Geom_data_pt[0]);

   // Geometric Data for the GeomObject has been set up (and pinned) in
   // constructor for geometric object. Now free the y-position 
   // of the centre because we want to determine it as an unknown
   internal_data_pt(Internal_geometric_data_index)->unpin(1);
  
   // Change cleanup responsibilities: The GeomData will now be killed
   // by the GeneralisedElement when it wipes its internal Data
   Must_clean_up=false;
  }
 
 
 /// Destructor: 
 virtual ~ElasticallySupportedRingElement()
  {
   // The GeomObject's GeomData is mirrored in the element's
   // Internal Data and therefore gets wiped in the
   // destructor of GeneralisedElement --> No need to kill it here
  }


 ///  Set pointer to Data object that specifies the "load"
 /// on the ElasticallySupportedRingElement
 void set_load_pt(Data* load_pt)
  {
#ifdef PARANOID
   if (load_pt->nvalue()!=1)
    {
     std::ostringstream error_stream;
     error_stream << "The data object that stores the load on the "
                  << "ElasticallySupportedRingElement\n"
                  << "should only contain a single data value\n"
                  << "This one contains " << load_pt->nvalue() << std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Add load to the element's external data and store
   // its index within that storage scheme: Following this assignment, 
   // the load Data is accessible from
   // GeneralisedElement::external_data_pt(External_load_index)
   External_load_index = add_external_data(load_pt);

   // Load has now been set
   Load_data_has_been_set=true;

  } // end of set_load_pt(...)
  

 /// "Load" acting on the ring
 double load()
  {
   // Return the load if it has been set
   if (Load_data_has_been_set)
    {
     return external_data_pt(External_load_index)->value(0);
    }
   // ...otherwise return zero load
   else
    {
     return 0.0;
    }
  } // end of load()
 

 /// Access function for the spring stiffness
 double& k_stiff() {return K_stiff;}


 /// Pin the vertical displacement
 void pin_yc()
  {
   // Vertical position of centre is stored as value 1 in the
   // element's one and only internal Data object.
   internal_data_pt(Internal_geometric_data_index)->pin(1);
  }


 /// Unpin the vertical displacement
 void unpin_yc()
  {
   // Vertical position of centre is stored as value 1 in the
   // element's one and only internal Data object.
   internal_data_pt(Internal_geometric_data_index)->unpin(1);

  } // end of unpin_yc()


 /// Compute element residual vector (wrapper)
 void get_residuals(Vector<double> &residuals)
  {
   //Initialise residuals to zero
   residuals.initialise(0.0);
   //Create a dummy matrix
   DenseMatrix<double> dummy(1);
   //Call the generic residuals function with flag set to 0
   fill_in_generic_residual_contribution(residuals,dummy,0);
   }

  
 /// Compute element residual Vector and element Jacobian matrix (wrapper)
 void get_jacobian(Vector<double> &residuals,
                   DenseMatrix<double> &jacobian)
  {
   //Initialise residuals to zero
   residuals.initialise(0.0);
   //Initialise the jacobian matrix to zero
   jacobian.initialise(0.0);
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,jacobian,1);

  } // end of get_jacobian(...)


 protected:


 ///  Compute element residual Vector (only if flag=0) and also
 /// the element Jacobian matrix (if flag=1)
 void fill_in_generic_residual_contribution(Vector<double> &residuals,
                                            DenseMatrix<double> &jacobian,
                                            unsigned flag)
  { 
   //Find out how may dofs there are in the element
   unsigned n_dof = ndof();
   //If everything is pinned return straight away
   if (n_dof==0) return;
 
   // Pseudo-elastic force balance to determine the position of the
   // ring's centre for a given load.

   // What's the local equation number of the force balance equation
   // [It's the equation that "determines" the value of the internal
   // dof, y_c, which is stored as the second value of the one-and-only
   // internal data object in this element]
   int local_eqn_number_for_yc = 
    internal_local_eqn(Internal_geometric_data_index,1);

   // Add residual to appropriate entry in the element's residual
   // vector:
   residuals[local_eqn_number_for_yc]=load()-K_stiff*y_c();

   // Work out Jacobian: 
   if (flag)
    {
     // Derivative of residual w.r.t. the internal dof, i.e. the vertical
     // position of the ring's centre: d residual[0]/d y_c
     jacobian(local_eqn_number_for_yc,local_eqn_number_for_yc) = -K_stiff;
     
     
     // Derivative with respect to external dof, i.e. the applied 
     // load: d residual[0]/d load -- but only if the load is an unknown
     if (n_dof==2)
      {
       // What's the local equation number of the load parameter?
       // It's stored as the 0th value in the the element's
       // one-and-only external data item:
       int local_eqn_number_for_load = 
        external_local_eqn(External_load_index,0);

#ifdef PARANOID
       if (local_eqn_number_for_load<0)
        {
         throw OomphLibError(
          "Load is pinned and yet n_dof=2?\n This is very fishy!\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
         
       // Add entry into element Jacobian
       jacobian(local_eqn_number_for_yc,local_eqn_number_for_load) = 1.0;
      }
    }
  } // end of get_residuals_generic(...)


private:

 /// Stiffness of the ring's "elastic" support
 double K_stiff;

 ///  Index of the location of the load Data in the element's 
 /// array of external data
 unsigned External_load_index;

 ///  Index of the location of the geometric Data in the element's 
 /// array of internal data
 unsigned Internal_geometric_data_index;

 /// Flag to indicate that load data has been set
 bool Load_data_has_been_set;

};

}

#endif

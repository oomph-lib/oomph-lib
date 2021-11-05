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
#ifndef OOMPH_COMMON_YOUNG_LAPLACE_STUFF_DOC
#define OOMPH_COMMON_YOUNG_LAPLACE_STUFF_DOC
#include <assert.h>

//===== start_of_namespace========================================
/// Namespace for "global" problem parameters
//================================================================
namespace GlobalParameters
{

 // Independent problem parameters:
 //--------------------------------

 /// Use spines (true) or not (false)
 bool Use_spines = true;
 
 /// Use height control (true) or not (false)?
 bool Use_height_control = true;

 /// Enumeration for the possible cases
 enum Cases{Spherical_cap_in_cylinder_pinned,
            All_pinned,
            Barrel_shape,
            T_junction_with_nonzero_contact_angle,
           };
 
 /// What case are we considering: Choose one from the enumeration Cases
 int Case = All_pinned;


 // "Physical parameters"
 //----------------------

 /// Contact angle and its cos (dependent parameter -- is reassigned)
 double Gamma = MathematicalConstants::Pi/4.0;
 double Cos_gamma=cos(Gamma);

 /// Pointer to Data object that stores the prescribed curvature
 Data* Kappa_pt = 0;

 /// Initial value for kappa
 double Kappa_initial = 0.0;

 /// Height control value
 double Controlled_height = 0.0;

 // Resolution parameters
 //----------------------

 /// Increase or decrease the value of the control parameters?
 int Step_sign = 1;

 /// Number of steps
 unsigned Nsteps = 5;

 /// Increment for prescribed curvature
 double Kappa_increment = -0.05;

 /// Increment for height control
 double Controlled_height_increment = 0.1; 

 /// Number of element in bulk mesh at which height control is applied.
 /// Initialise to 0 -- will be overwritte in 
 /// setup_dependent_parameters_and_sanity_check()
 unsigned Control_element = 0;

 // Mesh data 
 // ---------

 /// Length and width of the domain 
 double L_x = 1.0; 
 double L_y = 1.0; 

 /// Number of elements in the mesh
 unsigned N_x = 8;
 unsigned N_y = 8;

 // Spines data 
 // -----------

 /// Min. first spine angle against horizontal plane
 double Alpha_min = MathematicalConstants::Pi/2.0;

 /// Max. first spine angle against horizontal plane
 double Alpha_max = MathematicalConstants::Pi/2.0;

 /// Min. second spine angle against horizontal plane
 double Beta_min = MathematicalConstants::Pi/2.0;

 /// Max. second pine angle against horizontal plane
 double Beta_max = MathematicalConstants::Pi/2.0;
 
 /// Should the spines rotate in the x and y directions (true)?
 bool Rotate_spines_in_both_directions = true;

 // end of parameters

 //-------------------------------------------------------
 /// Setup dependent parameters and perform sanity check
 //-------------------------------------------------------
 void setup_dependent_parameters_and_sanity_check()
 {

  // Reset initial value for kappa
  Kappa_initial=0.0;
 
  // Check that we've got an even number of elements for control element
  if ((N_x%2!=0)||(N_y%2!=0))
   {
    cout << "n_x n_y should even" << endl;
    abort();
   }

  // Find control element
  Control_element=N_y*N_x/2+N_x/2;

  // Set up mesh and spines parameters
  if (Case==Spherical_cap_in_cylinder_pinned)
   {
    // Reset parameters (not realLy used for mesh in this 
    // case but for normalisation of spine rotation)
    L_x=1.0;
    L_y=1.0; 

    // Rotate outwards
    Alpha_min=MathematicalConstants::Pi/2.0;
    Alpha_max=MathematicalConstants::Pi/2.0*0.5;
    Rotate_spines_in_both_directions=true;
   } 
  else if (Case==All_pinned)
   {
    // Spines angles for all pinned boundary conditions
    Alpha_min=MathematicalConstants::Pi/2.0*1.5;
    Alpha_max=MathematicalConstants::Pi/2.0*0.5;
    Rotate_spines_in_both_directions=true;
   }
  else if (Case==Barrel_shape)
   {
    // Spines angles for barrel shaped validation
    Alpha_min=MathematicalConstants::Pi/2.0*1.5;
    Alpha_max=MathematicalConstants::Pi/2.0*0.5;
    Rotate_spines_in_both_directions=false;
   }
  else if (Case==T_junction_with_nonzero_contact_angle)
   {
    // Spines angles for T-junction with non nil contact angle
    Alpha_min=MathematicalConstants::Pi/2.0*1.5;
    Alpha_max=MathematicalConstants::Pi/2.0*0.5;
    Rotate_spines_in_both_directions=false;
   }  
  else
   {
    std::cout << "Never get here: Case = " << Case << std::endl;
    assert(false);
   }

  // Convert angle to cos
  Cos_gamma = cos(Gamma);

 } // end of set up

 // Spine functions
 //----------------

 /// Spine basis: The position vector to the basis of the spine
 /// as a function of the two coordinates x_1 and x_2, and its
 /// derivatives w.r.t. to these coordinates. 
 /// dspine_B[i][j] = d spine_B[j] / dx_i
 /// Spines start in the (x_1,x_2) plane at (x_1,x_2).
 void spine_base_function(const Vector<double>& x, 
                          Vector<double>& spine_B, 
                          Vector< Vector<double> >& dspine_B)
 {

   // Bspines and derivatives 
  spine_B[0]     = x[0];
  spine_B[1]     = x[1];
  spine_B[2]     = 0.0 ;
  dspine_B[0][0] = 1.0 ;
  dspine_B[1][0] = 0.0 ;
  dspine_B[0][1] = 0.0 ; 
  dspine_B[1][1] = 1.0 ;
  dspine_B[0][2] = 0.0 ;
  dspine_B[1][2] = 0.0 ;
  
 } // End of bspine functions

  
 /// Spine: The spine vector field as a function of the two 
 /// coordinates x_1 and x_2, and its derivatives w.r.t. to these coordinates:
 /// dspine[i][j] = d spine[j] / dx_i
 void spine_function(const Vector<double>& xx, 
                     Vector<double>& spine, 
                     Vector< Vector<double> >& dspine)
 {

  // Scale lengths
  Vector<double> x(2,0.0);
  x[0]=xx[0]/L_x;
  x[1]=xx[1]/L_y;

  // Which spine orientation do we have?
  if (!Rotate_spines_in_both_directions)
  {
    /// Spines (and derivatives)  are independent of x[0] and rotate 
    /// in the x[1]-direction
    spine[0]=0.0;     // Sx
    dspine[0][0]=0.0; // dSx/dx[0]
    dspine[1][0]=0.0; // dSx/dx[1]
    
    spine[1]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]);     // Sy
    dspine[0][1]=0.0;                                       // dSy/dx[0]
    dspine[1][1]=-sin(Alpha_min+(Alpha_max-Alpha_min)*x[1])
                    *(Alpha_max-Alpha_min)/L_y;             // dSy/dx[1]

    spine[2]=sin(Alpha_min+(Alpha_max-Alpha_min)*x[1]);     // Sz
    dspine[0][2]=0.0;                                       // dSz/dx[0]
    dspine[1][2]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]) 
                    *(Alpha_max-Alpha_min)/L_y;             // dSz/dx[1]
  }
  else
  {   
    /// Spines are dependent of x[0] AND x[1] and rotate in both directions
   spine[0]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[0]);      // Sx
   dspine[0][0]=-sin(Alpha_min+(Alpha_max-Alpha_min)*x[0])*
                    (Alpha_max-Alpha_min)/L_x;              // dSx/dx[0] 
   dspine[1][0]=0.0;                                        // dSx/dx[1]

   spine[1]=cos(Alpha_min+(Alpha_max-Alpha_min)*x[1]);      // Sy
   dspine[0][1]=0.0;                                        // dSy/dx[0]
   dspine[1][1]=-sin(Alpha_min+(Alpha_max-Alpha_min)*x[1])*
                    (Alpha_max-Alpha_min)/L_y;              // dSy/dx[1]

   spine[2]=1.0;      // Sz
   dspine[0][2]=0.0;  // dSz/dx[0]
   dspine[1][2]=0.0;  // dSz/dx[1]
  }
   

 } // End spine function

 // Exact kappa value
 //------------------
 double get_exact_kappa()
 { 
   if (Use_height_control)
    {
     
     if (Case==Spherical_cap_in_cylinder_pinned)
      {
       // Mean (!) curvature of spherical cap pinned in the 
       // quarter circular mesh (cylindrical tube)
       return 4.0*Controlled_height/
        (Controlled_height*Controlled_height+1.0);
      }
     else if (Case==Barrel_shape)
      {
       // Mean (!) curvature of barrel that goes through
       // the corners of the rectangular domain
       return 2.0*Controlled_height/
        (Controlled_height*Controlled_height+L_y*L_y/4.0);
      }
     else
      {
       std::cout << "No exact solution for this case..." << std::endl;
       return 0.0;
      }
    }
   else
    {
     // Return prescribed kappa no height control
     return 999; //Kappa_pt->value(0);
    }
 }


} // end of namespace

#endif

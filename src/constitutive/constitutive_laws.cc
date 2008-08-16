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
//Non-inline functions for constitutive laws and strain-energy functions

#include "constitutive_laws.h"

namespace oomph
{

//===============================================================
/// This function is used to check whether a matrix is square
//===============================================================
bool ConstitutiveLaw::is_matrix_square(const DenseMatrix<double> &M)
{
 //If the number rows and columns is not equal, the matrix is not square
 if(M.nrow() != M.ncol()) {return false;}
 else {return true;}
}

//========================================================================
/// This function is used to check whether matrices are of equal dimension
//========================================================================
bool ConstitutiveLaw::are_matrices_of_equal_dimensions(
 const DenseMatrix<double> &M1, const DenseMatrix<double> &M2)
{
 //If the numbers of rows and columns are not the same, then the
 //matrices are not of equal dimension
 if((M1.nrow() != M2.nrow()) || (M1.ncol() != M2.ncol())) {return false;}
 else {return true;}
}

//========================================================================
/// This function is used to provide simple error (bounce) checks on the
/// input to any calculate_second_piola_kirchhoff_stress
//=======================================================================
void ConstitutiveLaw::error_checking_in_input(const DenseMatrix<double> &g, 
                                              const DenseMatrix<double> &G,
                                              DenseMatrix<double> &sigma)
{
 //Test whether the undeformed metric tensor is square
 if(!is_matrix_square(g))
  {
   throw OomphLibError(
    "Undeformed metric tensor not square",
    "ConstitutiveLaw::errror_checking_in_input()",
    OOMPH_EXCEPTION_LOCATION);
  }
 
 //If the deformed metric tensor does not have the same dimension as
 //the undeformed tensor, complain
 if(!are_matrices_of_equal_dimensions(g,G))
  {
   std::string error_message =
    "Deformed metric tensor does  \n";
   error_message += 
    "not have same dimensions as the undeformed metric tensor\n";

   throw OomphLibError(error_message,
    "ConstitutiveLaw::errror_checking_in_input()",
    OOMPH_EXCEPTION_LOCATION);
  }

 //If the stress tensor does not have the same dimensions as the others
 //complain.
 if(!are_matrices_of_equal_dimensions(g,sigma))
  {
    std::string error_message =
     "Strain tensor passed to calculate_green_strain() does  \n";
   error_message += 
    "not have same dimensions as the undeformed metric tensor\n";

   throw OomphLibError(error_message,
    "ConstitutiveLaw::errror_checking_in_input()",
    OOMPH_EXCEPTION_LOCATION);
  }
}


//===========================================================================
/// The function to calculate the contravariant tensor from a covariant one
//===========================================================================
double ConstitutiveLaw::
calculate_contravariant( const DenseMatrix<double> &Gdown,
                         DenseMatrix<double> &Gup)
{
 //Initial error checking
#ifdef PARANOID
 //Test that the matrices are of the same dimension 
 if(!are_matrices_of_equal_dimensions(Gdown,Gup))
  {
   throw OomphLibError(
    "Matrices passed to calculate_contravariant() are not of equal dimension",
    "ConstitutiveLaw::calculate_contravariant()",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Find the dimension of the matrix
 unsigned dim = Gdown.ncol();

 //If it's not square, I don't know what to do (yet)
#ifdef PARANOID
 if(!is_matrix_square(Gdown))
  {
   std::string error_message =
    "Matrix passed to calculate_contravariant() is not square\n";
    error_message += 
     "non-square matrix inversion not implemented yet!\n";
    
    throw OomphLibError(error_message,
                        "ConstitutiveLaw::calculate_contravariant()",
                        OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Define the determinant of the matrix
 double det=0.0;

 //Now the inversion depends upon the dimension of the matrix
 switch(dim)
  {
   //Zero dimensions
  case 0:
   throw OomphLibError(
    "Zero dimensional matrix","ConstitutiveLaw::calculate_contravariant()",
    OOMPH_EXCEPTION_LOCATION);
   break;
   
   //One dimension
  case 1:
   //The determinant is just the value of the only entry
   det = Gdown(0,0);
   //The inverse is just the inverse of the value
   Gup(0,0) = 1.0/Gdown(0,0);
   break;

   
   //Two dimensions
  case 2:
   //Calculate the determinant
   det = Gdown(0,0)*Gdown(1,1) - Gdown(0,1)*Gdown(1,0);
   //Calculate entries of the contravariant tensor (inverse)
   Gup(0,0) = Gdown(1,1)/det;
   Gup(0,1) = -Gdown(0,1)/det;
   Gup(1,0) = -Gdown(1,0)/det;
   Gup(1,1) = Gdown(0,0)/det;
   break;

   ///Three dimensions
  case 3:
   //Calculate the determinant of the matrix
   det = Gdown(0,0)*Gdown(1,1)*Gdown(2,2) 
    + Gdown(0,1)*Gdown(1,2)*Gdown(2,0)
    + Gdown(0,2)*Gdown(1,0)*Gdown(2,1) 
    - Gdown(0,0)*Gdown(1,2)*Gdown(2,1)
    - Gdown(0,1)*Gdown(1,0)*Gdown(2,2) 
    - Gdown(0,2)*Gdown(1,1)*Gdown(2,0);

   //Calculate entries of the inverse matrix
   Gup(0,0) =  (Gdown(1,1)*Gdown(2,2) - Gdown(1,2)*Gdown(2,1))/det;
   Gup(0,1) = -(Gdown(0,1)*Gdown(2,2) - Gdown(0,2)*Gdown(2,1))/det;
   Gup(0,2) =  (Gdown(0,1)*Gdown(1,2) - Gdown(0,2)*Gdown(1,1))/det;
   Gup(1,0) = -(Gdown(1,0)*Gdown(2,2) - Gdown(1,2)*Gdown(2,0))/det;
   Gup(1,1) =  (Gdown(0,0)*Gdown(2,2) - Gdown(0,2)*Gdown(2,0))/det;
   Gup(1,2) = -(Gdown(0,0)*Gdown(1,2) - Gdown(0,2)*Gdown(1,0))/det;
   Gup(2,0) =  (Gdown(1,0)*Gdown(2,1) - Gdown(1,1)*Gdown(2,0))/det;
   Gup(2,1) = -(Gdown(0,0)*Gdown(2,1) - Gdown(0,1)*Gdown(2,0))/det;
   Gup(2,2) =  (Gdown(0,0)*Gdown(1,1) - Gdown(0,1)*Gdown(1,0))/det;
   break;
   
  default:
   throw OomphLibError("Dimension of matrix must be 0, 1,2 or 3\n",  
                       "ConstitutiveLaw::calculate_contravariant()",
                       OOMPH_EXCEPTION_LOCATION);
   break;
  }

 //Return the determinant of the matrix
 return(det);
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



//=====================================================================
/// \short Calculate the contravariant 2nd Piola Kirchhoff 
/// stress tensor. Arguments are the 
/// covariant undeformed (stress-free) and deformed metric 
/// tensors, g and G, and the matrix in which to return the stress tensor.
//=====================================================================
void GeneralisedHookean:: 
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g, 
                                        const DenseMatrix<double> &G,
                                        DenseMatrix<double> &sigma)
{
 //Error checking
#ifdef PARANOID
 error_checking_in_input(g,G,sigma);
#endif 

 //Find the dimension of the problem
 unsigned dim = G.nrow();
 
 //Calculate the contravariant Deformed metric tensor 
 DenseMatrix<double> Gup(dim);
 //We don't need the Jacobian so cast the function to void
 (void)calculate_contravariant(G,Gup);
  
 //Premultiply some constants
 double C1 = E/(2.0*(1.0+Nu)), C2 = 2.0*Nu/(1.0-2.0*Nu);
       
//  // Calculate the stiffness tensor; hard-coded as
//  // array because we don't want/have rank four tensor.
//  // Extra dimensions are not used if DIM<3; 2D elasticity
//  // is therefore plain strain version.
//  double Et[3][3][3][3];
   
//  //Loop over dimensions to calculate Et
//  for(unsigned i=0;i<dim;i++)
//   {
//    //Loop backwards over second index, so that we can use symmetry
//    //of first two indicies later
//    for(int j=(dim-1);j>=static_cast<int>(i);j--)
//     {
//      //If k is lower than i, then we can use symmetry of the first and 
//      //second pairs of indicies
//      for(unsigned k=0;k<i;k++)
//       {
//        for(unsigned l=0;l<dim;l++)
//         {
//          Et[i][j][k][l] = Et[k][l][i][j];
//         }
//       }
//      //Loop over cases when k is greater than or equal to i
//      //If k is less than or equal to j then we can again use symmetry
//      //of the first and seocnd indicies
//      for(int k=i;k<=j;k++)   
//       {
//        //If l is greater than j we have this symmetry
//        for(int l=(dim-1);l>j;l--)
//         {
//          Et[i][j][k][l] = Et[k][l][i][j];
//         }
//        //For l less then or equal to j and bigger than k,
//        //we actually have to do some work!
//        //Note that k is less than or equal to j from the above loop
//        for(int l=j;l>=k;l--)
//         {
//          Et[i][j][k][l] = C1*(Gup(i,k)*Gup(j,l) + Gup(i,l)*Gup(j,k) 
//                               + C2*Gup(i,j)*Gup(k,l));
//         }
//        //For l less than k can use symmetry
//        for(int l=(k-1);l>=0;l--)
//         {
//          Et[i][j][k][l] = Et[i][j][l][k];
//         }
//       }
//      //For cases when k is bigger than j can't use inner symmetry
//      for(unsigned k=(j+1);k<dim;k++)
//       {
//        //If l is bigger than k
//        for(int l=(dim-1);l>=static_cast<int>(k);l--)
//         {
//          Et[i][j][k][l] = C1*(Gup(i,k)*Gup(j,l) + Gup(i,l)*Gup(j,k) 
//                               + C2*Gup(i,j)*Gup(k,l));
//         }
//        //For l less than k can use symmetry of last two indicies
//        for(int l=(k-1);l>=0;l--)
//         {
//          Et[i][j][k][l] = Et[i][j][l][k];
//         }
//       }
//     }
   
//    //For j less than i, can use symmetry of the first two indicies
//    for(int j=static_cast<int>(i-1);j>=0;j--)
//     {
//      //Loop over all k and l
//      for(unsigned k=0;k<dim;k++)   
//       {
//        //Only do upper half here
//        for(unsigned l=0;l<dim;l++)
//         {
//          Et[i][j][k][l] = Et[j][i][k][l];
//         }
//       }
//     }
//   }

//  //Now merely calculate the components of the stress sigma
//  for(unsigned i=0;i<dim;i++)
//   {
//    for(unsigned j=0;j<dim;j++)
//     {
//      //Initialise this component of sigma
//      sigma(i,j) = 0.0;
//      //Linear approximation
//      for(unsigned k=0;k<dim;k++)
//       {
//        for(unsigned l=0;l<dim;l++)
//         {
//          sigma(i,j) += Et[i][j][k][l]*0.5*(G(k,l) - g(k,l));
//         }
//       }
//     }
//   }

//---------------------------------------

 // Strain tensor 
 DenseMatrix<double> strain(dim,dim);

 // Upper triangle
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=i;j<dim;j++)
    {
     strain(i,j)=0.5*(G(i,j) - g(i,j));
    }
  }

 // Copy across
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=0;j<i;j++)
    {
     strain(i,j)=strain(j,i);
    }
  }
 

 // hierher further symmetries?

 // Compute upper triangle of stress
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=i;j<dim;j++)
    {
     //Initialise this component of sigma
     sigma(i,j) = 0.0;
     for(unsigned k=0;k<dim;k++)
      {
       for(unsigned l=0;l<dim;l++)
        {
         sigma(i,j) += C1*(Gup(i,k)*Gup(j,l)+Gup(i,l)*Gup(j,k)+
                           C2*Gup(i,j)*Gup(k,l))*strain(k,l);
        }
      }
    }
  }

 // Copy across
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=0;j<i;j++)
    {
     sigma(i,j)=sigma(j,i);
    }
  }

}

//===========================================================================
/// Calculate the deviatoric part of the contravariant 
/// 2nd Piola Kirchoff stress tensor. Also return the contravariant
/// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
/// the inverse of the bulk modulus \f$ \kappa\f$. This form is appropriate
/// for near-incompressible materials for which
/// \f$ \sigma^{ij} = -p G^{ij} + \overline{ \sigma^{ij}}  \f$
/// where the "pressure" \f$ p \f$ is determined from
/// \f$ p / \kappa - d =0 \f$. 
//===========================================================================
void GeneralisedHookean:: 
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g, 
                                        const DenseMatrix<double> &G, 
                                        DenseMatrix<double> &sigma_dev, 
                                        DenseMatrix<double> &Gup,
                                        double &gen_dil, 
                                        double &inv_kappa)
{
 //Find the dimension of the problem
 unsigned dim = G.nrow();

 //Assign memory for the determinant of the deformed metric tensor
 double detG=0.0;

 //Compute deviatoric stress by calling the incompressible
 //version of this function
 calculate_second_piola_kirchhoff_stress(g,G,sigma_dev,Gup,detG);

 //Calculate the inverse of the "bulk" modulus
 inv_kappa = (1.0 - 2.0*Nu)*(1.0 + Nu)/(E*Nu);

 // Finally compute the generalised dilatation (i.e. the term that
 // must be zero if \kappa \to \infty
 gen_dil = 0.0;
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     gen_dil += Gup(i,j)*0.5*(G(i,j) - g(i,j));
    }
  }
 
}

//======================================================================
/// Calculate the deviatoric part 
/// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
/// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
/// Also return the contravariant deformed metric
/// tensor and the determinant of the deformed metric tensor.
/// This form is appropriate 
/// for truly-incompressible materials for which
/// \f$ \sigma^{ij} = - p G^{ij} +\overline{ \sigma^{ij}}  \f$
/// where the "pressure" \f$ p \f$ is determined by
/// \f$ \det G_{ij} - \det g_{ij} = 0 \f$. 
//======================================================================
void GeneralisedHookean:: 
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g, 
                                        const DenseMatrix<double> &G, 
                                        DenseMatrix<double> &sigma_dev, 
                                        DenseMatrix<double> &Gup,
                                        double &detG)
{
 //Error checking
#ifdef PARANOID
 error_checking_in_input(g,G,sigma_dev);
#endif 

 //Find the dimension of the problem
 unsigned dim = G.nrow();
 
 //Calculate the contravariant Deformed metric tensor 
 detG = calculate_contravariant(G,Gup);
      
 //Calculate the stiffness tensor, without the pressure bits
 double Et_NonSingular[3][3][3][3];
 
 //Premultiply the appropriate physical constant
 double C1 = E/(2.0*(1.0+Nu));
     
//  //Loop over dimensions to calculate Et_NonSingular
//  /*for(unsigned i=0;i<dim;i++)
//   {
//    //Loop backwards over second index, so that we can use symmetry
//    //of first two indicies later
//    for(int j=(dim-1);j>=static_cast<int>(i);j--)
//     {
//      //If k is lower than i, then we can use symmetry of the first and 
//      //second pairs of indicies
//      for(unsigned k=0;k<i;k++)
//       {
//        for(unsigned l=0;l<dim;l++)
//         {
//          Et_NonSingular[i][j][k][l] = Et_NonSingular[k][l][i][j];
//         }
//       }
//      //Loop over cases when k is greater than or equal to i
//      //If k is less than or equal to j then we can again use symmetry
//      //of the first and seocnd indicies
//      for(int k=i;k<=j;k++)   
//       {
//        //If l is greater than j we have this symmetry
//        for(int l=(dim-1);l>j;l--)
//         {
//          Et_NonSingular[i][j][k][l] = Et_NonSingular[k][l][i][j];
//         }
//        //For l less then or equal to j and bigger than k,
//        //we actually have to do some work!
//        //Note that k is less than or equal to j from the above loop
//        for(int l=j;l>=k;l--)
//         {
//          Et_NonSingular[i][j][k][l] = 
//           C1*(Gup(i,k)*Gup(j,l) + Gup(i,l)*Gup(j,k));
//         }
//        //For l less than k can use symmetry
//        for(int l=(k-1);l>=0;l--)
//         {
//          Et_NonSingular[i][j][k][l] = Et_NonSingular[i][j][l][k];
//         }
//       }
//      //For cases when k is bigger than j can't use inner symmetry
//      for(unsigned k=(j+1);k<dim;k++)
//       {
//        //If l is bigger than k
//        for(int l=(dim-1);l>=static_cast<int>(k);l--)
//         {
//          Et_NonSingular[i][j][k][l] = 
//           C1*(Gup(i,k)*Gup(j,l) + Gup(i,l)*Gup(j,k));
//         }
//        //For l less than k can use symmetry of last two indicies
//        for(int l=(k-1);l>=0;l--)
//         {
//          Et_NonSingular[i][j][k][l] = Et_NonSingular[i][j][l][k];
//         }
//       }
//     }
   
//    //For j less than i, can use symmetry of the first two indicies
//    for(int j=static_cast<int>(i-1);j>=0;j--)
//     {
//      //Loop over all k and l
//      for(unsigned k=0;k<dim;k++)   
//       {
//        //Only do upper half here
//        for(unsigned l=0;l<dim;l++)
//         {
//          Et_NonSingular[i][j][k][l] = Et_NonSingular[j][i][k][l];
//         }
//       }
//     }
//     }*/

//  for(unsigned i=0;i<dim;i++)
//   {
//    for(unsigned j=0;j<dim;j++)
//     {
//      for(unsigned k=0;k<dim;k++)
//       {
//        for(unsigned l=0;l<dim;l++)
//         {
//          Et_NonSingular[i][j][k][l] = 
//           C1*(Gup(i,k)*Gup(j,l) + Gup(i,l)*Gup(j,k));
//         }
//       }
//     }
//   }

//  //Now merely calculate the components of the stress sigma
//  for(unsigned i=0;i<dim;i++)
//   {
//    for(unsigned j=0;j<dim;j++)
//     {
//      //Initialise this component of sigma
//      sigma_dev(i,j) = 0.0;
//      //Linear approximation
//      for(unsigned k=0;k<dim;k++)
//       {
//        for(unsigned l=0;l<dim;l++)
//         {
//          sigma_dev(i,j) += Et_NonSingular[i][j][k][l]*0.5*(G(k,l) - g(k,l));
//         }
//       }
//     }
//   }

 // Strain tensor 
 DenseMatrix<double> strain(dim,dim);

 // Upper triangle
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=i;j<dim;j++)
    {
     strain(i,j)=0.5*(G(i,j) - g(i,j));
    }
  }

 // Copy across
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=0;j<i;j++)
    {
     strain(i,j)=strain(j,i);
    }
  }
 

 // hierher further symmetries?

 // Compute upper triangle of stress
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=i;j<dim;j++)
    {
     //Initialise this component of sigma
     sigma_dev(i,j) = 0.0;
     for(unsigned k=0;k<dim;k++)
      {
       for(unsigned l=0;l<dim;l++)
        {
         sigma_dev(i,j) += 
          C1*(Gup(i,k)*Gup(j,l)+Gup(i,l)*Gup(j,k))*strain(k,l);
        }
      }
    }
  }

 // Copy across
 for (unsigned i=0;i<dim;i++)
  {
   for (unsigned j=0;j<i;j++)
    {
     sigma_dev(i,j)=sigma_dev(j,i);
    }
  }




}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



//========================================================================
/// Calculate the contravariant 2nd Piola Kirchhoff 
/// stress tensor. Arguments are the 
/// covariant undeformed and deformed metric tensor and the 
/// matrix in which to return the stress tensor.
/// Uses correct 3D invariants for 2D (plane strain) problems.
//=======================================================================
void IsotropicStrainEnergyFunctionConstitutiveLaw::
calculate_second_piola_kirchhoff_stress(
 const DenseMatrix<double> &g, const DenseMatrix<double> &G,
 DenseMatrix<double> &sigma)
{
//Error checking
#ifdef PARANOID
 error_checking_in_input(g,G,sigma);
#endif 

 //Find the dimension of the problem
 unsigned dim = g.nrow();

#ifdef PARANOID
 if (dim==1)
  {
   std::string function_name =     
    "IsotropicStrainEnergyFunctionConstitutiveLaw::";
   function_name += "calculate_second_piola_kirchhoff_stress()";

   throw OomphLibError(
    "Check constitutive equations carefully when dim=1",
    function_name, OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Calculate the contravariant undeformed and deformed metric tensors
 //and get the determinants of the metric tensors 
 DenseMatrix<double> gup(dim), Gup(dim);
 double detg = calculate_contravariant(g,gup);
 double detG = calculate_contravariant(G,Gup);

 //Calculate the strain invariants
 Vector<double> I(3,0.0);
 //The third strain invaraint is the volumetric change
 I[2] = detG/detg;
 //The first and second are a bit more complex --- see G&Z
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     I[0] += gup(i,j)*G(i,j);
     I[1] += g(i,j)*Gup(i,j);
    }
  }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim==2)
  {
   I[0]+=1.0;
   I[1]+=1.0;
  }

 //Second strain invariant is multiplied by the third.
 I[1] *= I[2];

 //Calculate the derivatives of the strain energy function wrt the
 //strain invariants
 Vector<double> dWdI(3,0.0);
 Strain_energy_function_pt->derivatives(I,dWdI);


 //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
 //if the derivative wrt the second strain invariant is non-zero
 // \todo hierher Andrew this is dangerous! If this is to stay, we should
 // at least make the cutoff explicit and allow the user to change it.
 DenseMatrix<double> Bup(dim,dim,0.0);
 if(std::abs(dWdI[1]) > 0.0); //1.0e-10)
  {
   for(unsigned i=0;i<dim;i++)
    {
     for(unsigned j=0;j<dim;j++)
      {
       Bup(i,j) = I[0]*gup(i,j);
       for(unsigned r=0;r<dim;r++)
        {
         for(unsigned s=0;s<dim;s++)
          {
           Bup(i,j) -= gup(i,r)*gup(j,s)*G(r,s);
          }
        }
      }
    }
  }

 //Now set the values of the functions phi, psi and p (Green & Zerna notation)
 double phi = 2.0*dWdI[0]/sqrt(I[2]);
 double psi = 2.0*dWdI[1]/sqrt(I[2]);
 double p = 2.0*dWdI[2]*sqrt(I[2]);

 //Put it all together to get the stress
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     sigma(i,j) = phi*gup(i,j) + psi*Bup(i,j) + p*Gup(i,j);
    }
  }

}

//===========================================================================
/// Calculate the deviatoric part 
/// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
/// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
/// Also return the contravariant deformed metric
/// tensor and the determinant of the deformed metric tensor.
/// Uses correct 3D invariants for 2D (plane strain) problems.
//============================================================================
void IsotropicStrainEnergyFunctionConstitutiveLaw::
calculate_second_piola_kirchhoff_stress(
 const DenseMatrix<double> &g, const DenseMatrix<double> &G,
 DenseMatrix<double> &sigma, DenseMatrix<double> &Gup, double &detG)
{
//Error checking
#ifdef PARANOID
 error_checking_in_input(g,G,sigma);
#endif 

 //Find the dimension of the problem
 unsigned dim = g.nrow();


#ifdef PARANOID
 if (dim==1)
  {
   std::string function_name =     
    "IsotropicStrainEnergyFunctionConstitutiveLaw::";
   function_name += "calculate_second_piola_kirchhoff_stress()";
   
   throw OomphLibError(
    "Check constitutive equations carefully when dim=1",
    function_name, OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Calculate the contravariant undeformed and deformed metric tensors
 DenseMatrix<double> gup(dim);
 //Don't need this determinant
 (void)calculate_contravariant(g,gup);
 //These are passed back
 detG = calculate_contravariant(G,Gup);

 //Calculate the strain invariants
 Vector<double> I(3,0.0);
 //The third strain invaraint must be one (incompressibility)
 I[2] = 1.0;
 //The first and second are a bit more complex
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     I[0] += gup(i,j)*G(i,j);
     I[1] += g(i,j)*Gup(i,j);
    }
  }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim==2)
  {
   I[0]+=1.0;
   I[1]+=1.0;
  }

 //Calculate the derivatives of the strain energy function wrt the
 //strain invariants
 Vector<double> dWdI(3,0.0);
 Strain_energy_function_pt->derivatives(I,dWdI);

 //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
 //if the derivative wrt the second strain invariant is non-zero
 // \todo hierher Andrew this is dangerous! If this is to stay, we should
 // at least make the cutoff explicit and allow the user to change it.
 DenseMatrix<double> Bup(dim,dim,0.0);
 if(std::abs(dWdI[1]) > 0.0) //1.0e-10)
  {
   for(unsigned i=0;i<dim;i++)
    {
     for(unsigned j=0;j<dim;j++)
      {
       Bup(i,j) = I[0]*gup(i,j);
       for(unsigned r=0;r<dim;r++)
        {
         for(unsigned s=0;s<dim;s++)
          {
           Bup(i,j) -= gup(i,r)*gup(j,s)*G(r,s);
          }
        }
      }
    }
  }

 //Now set the values of the functions phi and psi (Green & Zerna notation)
 double phi = 2.0*dWdI[0];
 double psi = 2.0*dWdI[1];

 //Put it all together to get the stress
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     sigma(i,j) = phi*gup(i,j) + psi*Bup(i,j);
    }
  }
}

//===========================================================================
/// Calculate the deviatoric part of the contravariant 
/// 2nd Piola Kirchoff stress tensor. Also return the contravariant
/// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
/// the inverse of the bulk modulus \f$ \kappa\f$. 
/// Uses correct 3D invariants for 2D (plane strain) problems.
//===========================================================================
void IsotropicStrainEnergyFunctionConstitutiveLaw:: 
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g, 
                                        const DenseMatrix<double> &G, 
                                        DenseMatrix<double> &sigma, 
                                        DenseMatrix<double> &Gup,
                                        double &gen_dil, double &inv_kappa)
{

//Error checking
#ifdef PARANOID
 error_checking_in_input(g,G,sigma);
#endif 

 //Find the dimension of the problem
 unsigned dim = g.nrow();

#ifdef PARANOID
 if (dim==1)
  {
   std::string function_name =     
    "IsotropicStrainEnergyFunctionConstitutiveLaw::";
   function_name += "calculate_second_piola_kirchhoff_stress()";
   
   throw OomphLibError(
    "Check constitutive equations carefully when dim=1",
    function_name, OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Calculate the contravariant undeformed and deformed metric tensors
 //and get the determinants of the metric tensors 
 DenseMatrix<double> gup(dim);
 double detg = calculate_contravariant(g,gup);
 double detG = calculate_contravariant(G,Gup);

 //Calculate the strain invariants
 Vector<double> I(3,0.0);
 //The third strain invaraint is the volumetric change
 I[2] = detG/detg;
 //The first and second are a bit more complex --- see G&Z
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     I[0] += gup(i,j)*G(i,j);
     I[1] += g(i,j)*Gup(i,j);
    }
  }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim==2)
  {
   I[0]+=1.0;
   I[1]+=1.0;
  }

 //Second strain invariant is multiplied by the third.
 I[1] *= I[2];

 //Calculate the derivatives of the strain energy function wrt the
 //strain invariants
 Vector<double> dWdI(3,0.0);
 Strain_energy_function_pt->derivatives(I,dWdI);

 //Only bother to calculate the tensor B^{ij} (Green & Zerna notation)
 //if the derivative wrt the second strain invariant is non-zero
 // \todo hierher Andrew this is dangerous! If this is to stay, we should
 // at least make the cutoff explicit and allow the user to change it.
 DenseMatrix<double> Bup(dim,dim,0.0);
 if(std::abs(dWdI[1]) > 0.0) //1.0e-10)
  {
   for(unsigned i=0;i<dim;i++)
    {
     for(unsigned j=0;j<dim;j++)
      {
       Bup(i,j) = I[0]*gup(i,j);
       for(unsigned r=0;r<dim;r++)
        {
         for(unsigned s=0;s<dim;s++)
          {
           Bup(i,j) -= gup(i,r)*gup(j,s)*G(r,s);
          }
        }
      }
    }
  }

 //Now set the values of the functions phi and psi (Green & Zerna notation)
 double phi = 2.0*dWdI[0]/sqrt(I[2]);
 double psi = 2.0*dWdI[1]/sqrt(I[2]);

 //Choose inverse kappa to be one...
 inv_kappa = 1.0;

 //...then the generalised dilation is the same as p  in Green & Zerna's
 // notation
 gen_dil = 2.0*dWdI[2]*sqrt(I[2]);

 //Calculate the non-isotropic part of the stress
 for(unsigned i=0;i<dim;i++)
  {
   for(unsigned j=0;j<dim;j++)
    {
     sigma(i,j) = phi*gup(i,j) + psi*Bup(i,j);
    }
  }
}

}

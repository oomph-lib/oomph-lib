// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

// Oomph-lib includes
#include "generic.h"
#include <random>

using namespace std;
using namespace oomph;


//========================================================================
/// Premultiply mass matrix by this to have shift act on imaginary part.
//========================================================================
namespace SwapReAndIm
{

 // Swap factor for eigenvalues
 std::complex<double> Re_im_swap_factor=std::complex<double>(1.0,0.0);

}


//=====================================================================
/// Namespace for functions to test eigensolvers
//=====================================================================
namespace EigensolverTester
{
  
 /// Check eigensolutions, compute residual norm of complex eigen system,
 /// with (complex) matrices specified as args.
 double check_eigensolutions(
  const ComplexMatrixBase& A,
  const ComplexMatrixBase& M,
  const Vector<std::complex<double>>& eigenvalue,
  const Vector<Vector<std::complex<double>>>& eigenvector)
 {
  
  // Get the dimension of the matrix
  unsigned n = A.nrow();

  // Global residual
  double global_res=0.0;

  // Loop over eigenvalues
  unsigned n_eval=eigenvalue.size();
  for (unsigned eval=0;eval<n_eval;eval++)
   {
    // Skip non-finite ones
    if ((!isfinite(eigenvalue[eval].real()))||
        (!isfinite(eigenvalue[eval].imag())))
     {
      // oomph_info << "Eigenvalue " << eval << " : " << eigenvalue[eval]
      //            << " is not a finite number; skipping" << std::endl;
     }
    else
     {
      // Should really use black box matrix vector product which is likely
      // to be MUCH more efficient for sparse matrices.
      double res_norm=0.0;
      for (unsigned i=0;i<n;i++)
       {
        std::complex<double> res=std::complex<double>(0.0,0.0);
        //Storage for the two matrices multipled by the eigenvectors
        Vector<std::complex<double>> Av, Mv;

        //Multiply both matrices by the eigenvector
        //Need to cast away constness because multiply hasn't been
        //defined as a const function in the base class.
        const_cast<ComplexMatrixBase&>(A).multiply(eigenvector[eval],Av);
        const_cast<ComplexMatrixBase&>(M).multiply(eigenvector[eval],Mv);
        
        //Now calculate the residual of the eigenvalue problem
        //Av - lambda M v
        for (unsigned j=0;j<n;j++)
         {
          //res+=(A(i,j)-eigenvalue[eval]*M(i,j))*eigenvector[eval][j];
          res+= Av[j] - eigenvalue[eval]*Mv[j];
         }
        res_norm+=res.real()*res.real()+res.imag()*res.imag();
       }
      res_norm=sqrt(res_norm/double(n));
      global_res+=res_norm;
     }
   }
  
  global_res/=double(n_eval);
  return global_res;
  
 }


 /// Check eigensolutions, compute residual norm of complex eigen system,
 /// with (complex) matrices specified by problem object (zero shift)
 double check_eigensolutions(
  Problem* const& problem_pt,
  const Vector<std::complex<double>>& eigenvalue,
  const Vector<Vector<std::complex<double>>> & eigenvector)
 {
  // Get the dimension 
  unsigned n = eigenvector[0].size();

  // Out of laziness, assemble all the matrices everywhere so we don't have
  // have to bother about distribution below
  LinearAlgebraDistribution* dist_pt=new LinearAlgebraDistribution(
   problem_pt->communicator_pt(), n, false);
  
  // Allocated row compressed matrices for the mass matrix and main
  // matrix 
  CRDoubleMatrix M(dist_pt);
  CRDoubleMatrix A(dist_pt);
  
  // Assemble the matrices; pass the shift into the assembly (without shift)
  double zero_shift=0.0;
  problem_pt->get_eigenproblem_matrices(M,A,zero_shift);
  
  // Global residual
  double global_res=0.0;
  
  // Loop over eigenvalues
  unsigned n_eval=eigenvalue.size();
  for (unsigned eval=0;eval<n_eval;eval++)
   {
    // Ignore non-finite ones
    if ((!isfinite(eigenvalue[eval].real()))||
        (!isfinite(eigenvalue[eval].imag())))
     {
      // oomph_info << "Eigenvalue " << eval << " : " << eigenvalue[eval]
      //            << " is not a finite number; skipping" << std::endl;
     }
    else
     {
      double res_norm=0.0;
      for (unsigned i=0;i<n;i++)
       {
        // Should really use matrix vector product for this, but too lazy
        // since I'dhave to copy the eigenvectors into separate vectors,
        // deal with real and imag stuff etc.
        std::complex<double> res = std::complex<double>(0.0,0.0);
        for (unsigned j=0;j<n;j++)
         {
          res+=(A(i,j)-eigenvalue[eval]*M(i,j))*eigenvector[eval][j];
         }
        res_norm+=res.real()*res.real()+res.imag()*res.imag();
       }
      res_norm=sqrt(res_norm/double(n));
      global_res+=res_norm;
     }
   }
  global_res/=double(n_eval);

  // Clean up...
  delete dist_pt;

  //...and get the hell out of here
  return global_res;  
 }

}




/// Base eigenproblem element class used to generate the Jacobian and mass
/// matrix which correspond to the eigenproblem (J - lambda M) v = 0
class BaseEigenElement : public GeneralisedElement
{
public:
 BaseEigenElement() : N_value(0), Data_index(0) {}

  void set_size(const unsigned& n)
  {
    N_value = n;

    Data_index = add_internal_data(new Data(N_value));
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix) = 0;

protected:
  unsigned N_value;
  unsigned Data_index;
};

/// RandomAsymmetricEigenElement generates an eigenproblem whose Jacobian and
/// mass matrices have random integer elements between -128 and 127. Seed = 0.
class RandomAsymmetricEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    unsigned seed = 0;
    mt19937 generator(seed);
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += generator() % 256 - 128;
        mass_matrix(local_eqn, local_unknown) += generator() % 256 - 128;
      }
    }
  }
};

/// Fixed RandomAsymmetricEigenElement load an eigenproblem whose Jacobian
/// has random integer elements between -128 and 127 and an identity mass
/// matrix.
class FixedRandomAsymmetricEigenElement : public BaseEigenElement
{
public:
  // Override set_size to ensure the problem is 64x64
  void set_size(const unsigned& n)
  {
   if(n != 64)
    {
     std::ostringstream error_stream;
     error_stream << "Size of problem must be 64."
                  << "It is being set to " << n << std::endl;
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   BaseEigenElement::set_size(n);
  }

  // Implement creation of eigenproblem matrices
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    // Initialise the random number generator to a fixed seed
    ifstream input_stream;
    input_stream.open("random_test_matrix.dat");
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);

        string buffer;
        getline(input_stream, buffer, ',');
        double jac = stod(buffer);

        jacobian(local_eqn, local_unknown) += jac;
        if (local_eqn == local_unknown)
        {
          mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};

/// Eigenproblem class. Creates a mesh with a single eigenelement and assigns
/// the appropriate equation numbers.
template<class ELEMENT>
class Eigenproblem : public Problem
{
 //Storage for the size
 unsigned Size;

public:
 //Choose the default size equal to the stored matrix
 Eigenproblem(const unsigned& size=64) : Size(size)
  {
    this->mesh_pt() = new Mesh;

    ELEMENT* el_pt = new ELEMENT;

    el_pt->set_size(size);

    this->mesh_pt()->add_element_pt(el_pt);

    assign_eqn_numbers();
  }

  ~Eigenproblem()
  {
    delete this->mesh_pt();
  }


void doc_solution(string filename,
                  Vector<std::complex<double>>& eigenvalue,
                  Vector<DoubleVector>& eigenvector_real,
                  Vector<DoubleVector>& eigenvector_imag)
{
  //Open the output stream
  std::ofstream output_stream;
  output_stream.open(filename);
  //output the eigenvalues
  const unsigned N = eigenvalue.size();
  for (unsigned i = 0; i < N; i++)
  {
     output_stream << eigenvalue[i].real() << " , " << eigenvalue[i].imag()
                   << endl;
  }
  output_stream << endl;


  //Set up additional storage for the eigenvectors as complex numbers
  Vector<Vector<std::complex<double>>> eigenvector_size_2n;
  
 // If the problem has been distributed collect all the data by redistribution
 // into non-distributed vectors
  unsigned n_dof = this->ndof();
  LinearAlgebraDistribution* non_distributed_dist_pt=
   new LinearAlgebraDistribution(this->communicator_pt(),
                                 n_dof, false);
  //Read out the number of expected eigenvalues from the size of the problem
  unsigned n_eval = this->Size;
  for (unsigned eval=0;eval<n_eval;eval++)
   {
    eigenvector_real[eval].redistribute(non_distributed_dist_pt);
    eigenvector_imag[eval].redistribute(non_distributed_dist_pt);
   }
  
  // Now copy the eigenvectors across into complex vectors
  eigenvector_size_2n.resize(n_eval);
  for (unsigned eval=0;eval<n_eval;eval++)
   {
    eigenvector_size_2n[eval].resize(n_dof);
    for (unsigned i=0;i<n_dof;i++)
     {
      eigenvector_size_2n[eval][i]=
       std::complex<double>(eigenvector_real[eval][i],
                            eigenvector_imag[eval][i]);
     }
   }
  
  //Rescale the eigenvalue to be the expected size
  eigenvalue.resize(n_eval);
  
  //Find out the average residuals of the eigenvalue problems
  double global_eigenproblem_res = 
   EigensolverTester::check_eigensolutions(this,
                                           eigenvalue,
                                           eigenvector_size_2n);

  output_stream << "Global Eigenproblem Average Error "
                << global_eigenproblem_res << std::endl;

  if(global_eigenproblem_res < 1.0e-13)
   {
    output_stream << "Test Passed" << std::endl;
   }
  else
   {
    output_stream << "Test Failed" << std::endl;
   }
  
  output_stream.close();
}

};


/// Main function. Apply solver tests to each eigensolver.
int main()
{
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  // Matrix dimensions
  const unsigned N = 64;

  // Create eigenproblem
  Eigenproblem<FixedRandomAsymmetricEigenElement> problem(N);

  // Set up additional arguments
  // Output all eigenvalues
  unsigned n_eval = N;

  // Store outputs
  Vector<complex<double>> eigenvalue(N);
  Vector<DoubleVector> eigenvector_real(N);
  Vector<DoubleVector> eigenvector_imag(N);

  /// Test solve_eigenproblem
  problem.solve_eigenproblem(
    n_eval, eigenvalue, eigenvector_real, eigenvector_imag);
  problem.doc_solution(doc_info.directory() + "solve_eigenproblem_test.dat",
               eigenvalue,
               eigenvector_real,
               eigenvector_imag);

  
  return (EXIT_SUCCESS);
}

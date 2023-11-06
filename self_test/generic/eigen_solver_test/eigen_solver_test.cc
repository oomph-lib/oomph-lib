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
/// Base eigenproblem element class used to generate the Jacobian and mass
/// matrix which correspond to the eigenproblem (J - lambda M) v = 0
//========================================================================
class BaseEigenElement : public GeneralisedElement
{
public:
 
 BaseEigenElement() {}

 /// Set problem size
 void set_size(const unsigned& n)
  {
   N_value = n;

   Data_index = add_internal_data(new Data(N_value));
  }

 /// Create the matrices in the eigenproblem, equivalent to the Jacobian and
 /// mass matrix (virtual)
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix) = 0;

protected:
 
 /// Store the problem size and the data index
 unsigned N_value;
 unsigned Data_index;
};

//==========================================================================
/// IdentityEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are equal to the identity matrix
//==========================================================================
class IdentityEigenElement : public BaseEigenElement
{
 
public:
 
 /// Implement the Jacobian and mass matrix construction
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix)
  {
   for (unsigned i = 0; i < N_value; i++)
    {
     unsigned local_eqn = internal_local_eqn(Data_index, i);
     for (unsigned j = 0; j < N_value; j++)
      {
       unsigned local_unknown = internal_local_eqn(Data_index, j);
       if (i == j)
        {
         jacobian(local_eqn, local_unknown) += 1;
         mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};


//==========================================================================
/// AsymmetricEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are simple asymmetric examples such that the eigenvalues should be
/// the integers between 1 and 64 inclusive
//==========================================================================
class AsymmetricEigenElement : public BaseEigenElement
{
public:
 
 /// Implement creation of the matrices
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix)
  {
   for (unsigned i = 0; i < N_value; i++)
    {
     unsigned local_eqn = internal_local_eqn(Data_index, i);
     for (unsigned j = 0; j < N_value; j++)
      {
       unsigned local_unknown = internal_local_eqn(Data_index, j);
       // Set the elements of the Jacobian's diagonal to 1, 2, 3, ..., 64
       // and the mass matrix's diagonal to 1, 1, ..., 1
       if (i == j)
        {
         jacobian(local_eqn, local_unknown) += i + 1.0;
         mass_matrix(local_eqn, local_unknown) += 1.0;
        }
       // Set the upper diagonal elements to one for both matrices
       else if (i > j)
        {
         jacobian(local_eqn, local_unknown) += 1.0;
         mass_matrix(local_eqn, local_unknown) += 1.0;
        }
      }
    }
  }
};



//==========================================================================
/// RandomAsymmetricEigenElement generates an eigenproblem whose Jacobian and
/// mass matrices have random integer elements between -128 and 127. Seed = 0
//==========================================================================
class RandomAsymmetricEigenElement : public BaseEigenElement
{
public:
 
 /// Implement creation of eigenproblem matrices
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix)
  {
   // Initialise the random number generator to a fixed seed
   unsigned seed = 0;
   mt19937 generator(seed);
   for (unsigned i = 0; i < N_value; i++)
    {
     unsigned local_eqn = internal_local_eqn(Data_index, i);
     for (unsigned j = 0; j < N_value; j++)
      {
       unsigned local_unknown = internal_local_eqn(Data_index, j);
       // Create two dense, random, asymmetric matrices
       jacobian(local_eqn, local_unknown) += generator() % 256 - 128;
       mass_matrix(local_eqn, local_unknown) += generator() % 256 - 128;
      }
    }
  }
};



//=========================================================================
/// Fixed RandomAsymmetricEigenElement load an eigenproblem whose Jacobian
/// has random integer elements between -128 and 127 and an identity mass
/// matrix.
//=========================================================================
class FixedRandomAsymmetricEigenElement : public BaseEigenElement
{
public:
 
 /// Override set_size to ensure the problem is 64x64
 void set_size(const unsigned& n)
  {
   /// Override the input argument as the Rosser matrix is fixed at 8x8
   N_value = 64;

   Data_index = add_internal_data(new Data(N_value));
  }

 /// Implement creation of eigenproblem matrices
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


//==========================================================================
/// RosserSymmetricEigenElement generates the classic Rosser eigenproblem
/// matrices. Size 8x8.
//==========================================================================
class RosserSymmetricEigenElement : public BaseEigenElement
{
public:
 
 /// Override set_size to ensure the problem is 8x8
 void set_size(const unsigned& n)
  {
   /// Override the input argument as the Rosser matrix is fixed at 8x8
   N_value = 8;

   Data_index = add_internal_data(new Data(N_value));
  }

 /// Implement creation of Jacobian and mass matrices by setting the Jacobian
 /// equal to the Rosser matrix and the mass matrix to the identity matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix)
  {
   // The Rosser matrix, eigenvalues = 10*sqrt(10405), 1020, 510 +
   // 100*sqrt(26), 1000, 1000, 510 - 100*sqrt(26), 0, -10*sqrt(10405)
   int A[8][8] = {{611, 196, -192, 407, -8, -52, -49, 29},
                  {196, 899, 113, -192, -71, -43, -8, -44},
                  {-192, 113, 899, 196, 61, 49, 8, 52},
                  {407, -192, 196, 611, 8, 44, 59, -23},
                  {-8, -71, 61, 8, 411, -599, 208, 208},
                  {-52, -43, 49, 44, -599, 411, 208, 208},
                  {-49, -8, 8, 59, 208, 208, 99, -911},
                  {29, -44, 52, -23, 208, 208, -911, 99}};

   for (unsigned i = 0; i < N_value; i++)
    {
     unsigned local_eqn = internal_local_eqn(Data_index, i);
     for (unsigned j = 0; j < N_value; j++)
      {
       unsigned local_unknown = internal_local_eqn(Data_index, j);
       jacobian(local_eqn, local_unknown) += double(A[i][j]);
       if (i == j)
        {
         mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};

//=========================================================================
/// Eigenproblem class. Creates a mesh with a single eigenelement and assigns
/// the appropriate equation numbers.
//=========================================================================
template<class ELEMENT>
class Eigenproblem : public Problem
{
public:


 /// Constructor
 Eigenproblem(const unsigned& size)
  {
   // Create mesh
   this->mesh_pt() = new Mesh;
   
   // Create a single eigenproblem element
   ELEMENT* el_pt = new ELEMENT;
   
   // Set problem size
   el_pt->set_size(size);
   
   // Add element to mesh
   this->mesh_pt()->add_element_pt(el_pt);
   
   // Assign eqn numbers
   assign_eqn_numbers();
  }
 
 // Class destructor. Delete mesh
 ~Eigenproblem()
  {
   delete this->mesh_pt();
  }
};



//=========================================================================
// Solve Eigenproblem Test class. Create an eigenproblem on a templated ELEMENT
// and call the solve_eigenproblem function of the passed eigen_solver_pt.
//=========================================================================
template<class ELEMENT>
class SolveEigenProblemTest
{
 
public:
 
 /// Constructor
 SolveEigenProblemTest(EigenSolver* const& eigen_solver_pt,
                       const unsigned& N,
                       const unsigned& n_timing_loops,
                       DocInfo* const& doc_info_pt,
                       bool do_adjoint_problem)
  : Eigen_solver_pt(eigen_solver_pt),
    Problem_pt(0),
    Matrix_size(N),
    N_eval(0),
    Do_adjoint_problem(do_adjoint_problem),
    N_timing_loops(n_timing_loops),
    Doc_info_pt(doc_info_pt)
  {
   // Create Eigenproblem
   Problem_pt = new Eigenproblem<ELEMENT>(Matrix_size);
   
   // Output the first 8 eigenvalues
   N_eval = 8;
   
   // Store outputs
   Vector<complex<double>> eval(N_eval);
   Vector<DoubleVector> eigenvector_real(N_eval);
   Vector<DoubleVector> eigenvector_imag(N_eval);
   
   // Start clock
   clock_t t_start = clock();
   for (unsigned i = 0; i < N_timing_loops; i++)
    {
     // Call solve_eigenproblem
     Eigen_solver_pt->solve_eigenproblem(Problem_pt,
                                         N_eval,
                                         eval,
                                         eigenvector_real,
                                         eigenvector_imag,
                                         Do_adjoint_problem);
    }
   // Stop clock
   clock_t t_end = clock();
   
   
   // Document duration
   double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
   ofstream timing_stream;
   timing_stream.open("timing.dat", ios_base::app);
   timing_stream << "test" << Doc_info_pt->number()
                 << ", time: " << t_length / double(N_timing_loops) << endl;
   timing_stream.close();
   
   // Document solution
   string filename = Doc_info_pt->directory() + "test" +
    to_string(Doc_info_pt->number()) + ".dat";
   ofstream output_stream;
   output_stream.open(filename);
   for (unsigned i = 0; i < N_eval; i++)
    {
     output_stream << eval[i].real() << " " << eval[i].imag() << endl;
    }
   output_stream.close();
   
   // Increment doc info number
   Doc_info_pt->number()++;
  }
 
private:
 
 /// Eigensolver pointer
 EigenSolver* Eigen_solver_pt;
 
 /// Store the eigenproblem 
 Problem* Problem_pt;

 /// Matrix size
 unsigned Matrix_size;
 
 // Number of eigenvalues needed
 unsigned N_eval;

 /// Adjoint or normal problem?
 bool Do_adjoint_problem;
 
 // Store the number of times the function should be call for improved timing
 unsigned N_timing_loops;
 
 // Store a pointer to doc_info
 DocInfo* Doc_info_pt;
};





//============================================================================
// Solve Eigenproblem Test class (Legacy version). Create an eigenproblem on a
// templated ELEMENT and call the solve_eigenproblem function of the passed
// eigen_solver_pt.
//============================================================================
template<class ELEMENT>
class SolveEigenProblemLegacyTest
{
public:

 /// Constructor
 SolveEigenProblemLegacyTest(EigenSolver* const& eigen_solver_pt,
                             const unsigned& N,
                             const unsigned& n_timing_loops,
                             DocInfo* const& doc_info_pt,
                             bool do_adjoint_problem)
  : Eigen_solver_pt(eigen_solver_pt),
    Problem_pt(0),
    Matrix_size(N),
    N_eval(0),
    Do_adjoint_problem(do_adjoint_problem),
    N_timing_loops(n_timing_loops),
    Doc_info_pt(doc_info_pt)
  {
   // Create Eigenproblem
   Problem_pt = new Eigenproblem<ELEMENT>(Matrix_size);
   
   // Output the first 8 eigenvalues
   N_eval = 8;
   
   // Store outputs
   Vector<complex<double>> eval(N_eval);
   Vector<DoubleVector> evec(N_eval);
   
   // Start clock
   clock_t t_start = clock();
   for (unsigned i = 0; i < N_timing_loops; i++)
    {
     // Call solve_eigenproblem
     Eigen_solver_pt->solve_eigenproblem_legacy(
      Problem_pt, N_eval, eval, evec, Do_adjoint_problem);
    }
   // Stop clock
   clock_t t_end = clock();
   
   // Document duration
   double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
   ofstream timing_stream;
   timing_stream.open("timing.dat", ios_base::app);
   timing_stream << "test" << Doc_info_pt->number()
                 << ", time: " << t_length / double(N_timing_loops) << endl;
   timing_stream.close();
   
   // Document solution
   string filename = Doc_info_pt->directory() + "test" +
    to_string(Doc_info_pt->number()) + ".dat";
   
   ofstream output_stream;
   output_stream.open(filename);
   for (unsigned i = 0; i < N_eval; i++)
    {
     output_stream << eval[i].real() << " " << eval[i].imag() << endl;
    }
   output_stream.close();
   
   // Increment doc info number
   Doc_info_pt->number()++;
  }
 
private:

 
 /// Eigen solver pointer
 EigenSolver* Eigen_solver_pt;
 
 /// Eigenproblem
 Problem* Problem_pt;

 /// Matris size
 unsigned Matrix_size;
 
 /// Number of eigenvalues required
 unsigned N_eval;

 /// Do adjoint or normal problem?
 bool Do_adjoint_problem;
 
 /// Store the number of times the function should be call for improved timing
 unsigned N_timing_loops;
 
 /// Pointer to doc_info
 DocInfo* Doc_info_pt;
};



//=======================================================================
// Test the LAPACK_QZ solver against the appropriate problem and methods.
//=======================================================================
void test_lapack_qz(const unsigned N,
                    const unsigned n_timing_loops,
                    DocInfo* doc_info_pt)
{
 // Create a new eigensolver
 EigenSolver* eigen_solver_pt = new LAPACK_QZ;
 
 // Do not test the unimplement adjoint problem
 const bool do_adjoint_problem = false;
 
 // Test the regular solve_eigenproblem
 SolveEigenProblemTest<IdentityEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemTest<RosserSymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemTest<AsymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemTest<FixedRandomAsymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 
 // Test the legacy solve_eigenproblem
 SolveEigenProblemLegacyTest<IdentityEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemLegacyTest<RosserSymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemLegacyTest<AsymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 SolveEigenProblemLegacyTest<FixedRandomAsymmetricEigenElement>(
  eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
 
 // Free the eigen_solver_pt
 delete eigen_solver_pt;
}





//=======================================================================
/// Main function. Call all the testing functions
//=======================================================================
int main(int argc, char** argv)
{
 // Number of times to repeat the operation for better timings
 const unsigned n_timing_loops = 2;
 
 // Matrix dimensions
 const unsigned N = 64;
 
 // Create a DocInfo
 DocInfo* doc_info_pt = new DocInfo;

 // Set directory to lapack
 doc_info_pt->set_directory("RESLT_lapack/");

 // Add a header to the timing data stream
 ofstream timing_stream;
 timing_stream.open("timing.dat", ios_base::app);
 timing_stream << "LAPACK_QZ" << endl;
 timing_stream.close();

 // Call test lapack qz
 test_lapack_qz(N, n_timing_loops, doc_info_pt);

 // Delete doc_info_pt
 delete doc_info_pt;

 return 0;
}

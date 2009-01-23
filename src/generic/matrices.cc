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
//Non-inline member functions for the matrix classes

#include<set> 
#include<map> 

//oomph-lib headers
#include "matrices.h"
#include "linear_solver.h"

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

namespace oomph
{

//=========================================================================
/// Namespaces for (Numerical-Recipes-based) functions for
/// eigensolver, based on Jacobi rotations.
//=========================================================================
namespace JacobiEigenSolver
{

 /// Perform one Jacobi rotation on matrix a
 inline void rot(DenseDoubleMatrix&a, const double s, const double tau, 
                 const unsigned long i, const unsigned long j, 
                 const unsigned long k, const unsigned long l)
 {
  double g,h;
  
  g=a(i,j);
  h=a(k,l);
  a(i,j)=g-s*(h+g*tau);
  a(k,l)=h+s*(g-h*tau);

 }


/// \short Use Jacobi rotations to determine eigenvalues and eigenvectors of 
/// matrix a. d[i]=i-th eigenvalue; v(i,j)=i-th component of j-th eigenvector
/// (note that this is the transpose of what we'd like to have...);
/// nrot=number of rotations used. 
 void jacobi(DenseDoubleMatrix& a, Vector<double>& d, 
             DenseDoubleMatrix& v, unsigned long& nrot)
 {
#ifdef PARANOID
  // Check Matrix a is square
  if (a.ncol()!=a.nrow())
   {
    throw OomphLibError(
     "This matrix is not square, the matrix MUST be square!",
     "JacobiEigenSolver::jacobi()",
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
 
  // If matrix v is wrong size, correct it!
  if (v.ncol()!=a.ncol() || v.nrow()!=a.nrow())
   {
    v.resize(a.nrow(),a.nrow(),0.0);
   }

  unsigned long i,j,ip,iq;
  double tresh,theta,tau,t,sm,s,h,g,c;
 
  unsigned long n=d.size();
  Vector<double> b(n);
  Vector<double> z(n);
  for (ip=0;ip<n;ip++) {
   for (iq=0;iq<n;iq++) v(ip,iq)=0.0;
   v(ip,ip)=1.0;
  }
  for (ip=0;ip<n;ip++) {
   b[ip]=d[ip]=a(ip,ip);
   z[ip]=0.0;
  }
  nrot=0;
  for (i=1;i<=50;i++) {
   sm=0.0;
   for (ip=0;ip<n-1;ip++) {
    for (iq=ip+1;iq<n;iq++)
     sm += std::abs(a(ip,iq));
   }
   if (sm == 0.0)
    return;
   if (i < 4)
    tresh=0.2*sm/(n*n);
   else
    tresh=0.0;
   for (ip=0;ip<n-1;ip++) {
    for (iq=ip+1;iq<n;iq++) {
     g=100.0*std::abs(a(ip,iq));
     if (i > 4 && (std::abs(d[ip])+g) == std::abs(d[ip])
         && (std::abs(d[iq])+g) == std::abs(d[iq]))
      a(ip,iq)=0.0;
     else if (std::abs(a(ip,iq)) > tresh) {
      h=d[iq]-d[ip];
      if ((std::abs(h)+g) == std::abs(h))
       t=(a(ip,iq))/h;
      else {
       theta=0.5*h/(a(ip,iq));
       t=1.0/(std::abs(theta)+std::sqrt(1.0+theta*theta));
       if (theta < 0.0) t = -t;
      }
      c=1.0/std::sqrt(1+t*t);
      s=t*c;
      tau=s/(1.0+c);
      h=t*a(ip,iq);
      z[ip] -= h;
      z[iq] += h;
      d[ip] -= h;
      d[iq] += h;
      a(ip,iq)=0.0;
      for (j=0;j<ip;j++)
       rot(a,s,tau,j,ip,j,iq);
      for (j=ip+1;j<iq;j++)
       rot(a,s,tau,ip,j,j,iq);
      for (j=iq+1;j<n;j++)
       rot(a,s,tau,ip,j,iq,j);
      for (j=0;j<n;j++)
       rot(v,s,tau,j,ip,j,iq);
      ++nrot;
     }
    }
   }
   for (ip=0;ip<n;ip++) {
    b[ip] += z[ip];
    d[ip]=b[ip];
    z[ip]=0.0;
   }
  }
  throw OomphLibError(
   "Too many iterations in routine jacobi",
   "JacobiEigenSolver::jacobi()",
   OOMPH_EXCEPTION_LOCATION);
 }

}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//============================================================================
/// Complete LU solve (overwrites RHS with solution). This is the
/// generic version which should not need to be over-written.
//============================================================================
void DoubleMatrixBase::solve(DoubleVector &rhs)
{
#ifdef PARANOID
 if(Linear_solver_pt==0)
  {
   throw OomphLibError("Linear_solver_pt not set in matrix",
                       "DoubleMatrixBase::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Copy rhs vector into local storage so it doesn't get overwritten
 // if the linear solver decides to initialise the solution vector, say,
 // which it's quite entitled to do!
 DoubleVector actual_rhs(rhs);

 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,actual_rhs,rhs);
}

//============================================================================
/// Complete LU solve (Nothing gets overwritten!). This generic
/// version should never need to be overwritten
//============================================================================
void DoubleMatrixBase::solve(const DoubleVector &rhs, 
                             DoubleVector &soln)
{
#ifdef PARANOID
 if(Linear_solver_pt==0)
  {
   throw OomphLibError("Linear_solver_pt not set in matrix",
                       "DoubleMatrixBase::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,rhs,soln);
}

//============================================================================
/// Complete LU solve (overwrites RHS with solution). This is the
/// generic version which should not need to be over-written.
//============================================================================
void DoubleMatrixBase::solve(Vector<double> &rhs)
{
#ifdef PARANOID
 if(Linear_solver_pt==0)
  {
   throw OomphLibError("Linear_solver_pt not set in matrix",
                       "DoubleMatrixBase::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Copy rhs vector into local storage so it doesn't get overwritten
 // if the linear solver decides to initialise the solution vector, say,
 // which it's quite entitled to do!
 Vector<double> actual_rhs(rhs);

 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,actual_rhs,rhs);
}

//============================================================================
/// Complete LU solve (Nothing gets overwritten!). This generic
/// version should never need to be overwritten
//============================================================================
void DoubleMatrixBase::solve(const Vector<double> &rhs,
                             Vector<double> &soln)
{
#ifdef PARANOID
 if(Linear_solver_pt==0)
  {
   throw OomphLibError("Linear_solver_pt not set in matrix",
                       "DoubleMatrixBase::solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,rhs,soln);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//===============================================================
/// Constructor, set the default linear solver to be the DenseLU 
/// solver
//===============================================================
DenseDoubleMatrix::DenseDoubleMatrix(): DenseMatrix<double>()
{
 Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
}

//==============================================================
/// Constructor to build a square n by n matrix.
/// Set the default linear solver to be DenseLU
//==============================================================
DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n) : 
 DenseMatrix<double>(n)
{
 Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
}
 

//=================================================================
/// Constructor to build a matrix with n rows and m columns.
/// Set the default linear solver to be DenseLU
//=================================================================
 DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n, 
                                      const unsigned long &m) :
  DenseMatrix<double>(n,m)
{
 Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
}

//=====================================================================
/// Constructor to build a matrix with n rows and m columns,
/// with initial value initial_val
/// Set the default linear solver to be DenseLU
//=====================================================================
DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n, 
                                     const unsigned long &m,
                                     const double &initial_val) :
 DenseMatrix<double>(n,m,initial_val) 
{
 Linear_solver_pt = Default_linear_solver_pt = new DenseLU;
}

//=======================================================================
/// Destructor delete the default linear solver
//======================================================================
DenseDoubleMatrix::~DenseDoubleMatrix()
{
 //Delete the default linear solver
 delete Default_linear_solver_pt;
}

//============================================================================
/// LU decompose a matrix, by using the default linear solver
/// (DenseLU)
//============================================================================
void DenseDoubleMatrix::ludecompose()
{
 //Use the default (DenseLU) solver to ludecompose the matrix
 static_cast<DenseLU*>(Default_linear_solver_pt)->factorise(this);
}


//============================================================================
///  Back substitute an LU decomposed matrix.
//============================================================================
void DenseDoubleMatrix::lubksub(DoubleVector &rhs)
{
 //Use the default (DenseLU) solver to perform the backsubstitution
 static_cast<DenseLU*>(Default_linear_solver_pt)->backsub(rhs,rhs);
}

//============================================================================
///  Back substitute an LU decomposed matrix.
//============================================================================
void DenseDoubleMatrix::lubksub(Vector<double> &rhs)
{
 //Use the default (DenseLU) solver to perform the backsubstitution
 static_cast<DenseLU*>(Default_linear_solver_pt)->backsub(rhs,rhs);
}


//============================================================================
///  Determine eigenvalues and eigenvectors, using
/// Jacobi rotations. Only for symmetric matrices. Nothing gets overwritten!
/// - \c eigen_vect(i,j) = j-th component of i-th eigenvector.
/// - \c eigen_val[i] is the i-th eigenvalue; same ordering as in eigenvectors
//============================================================================
void DenseDoubleMatrix::eigenvalues_by_jacobi(Vector<double> & eigen_vals, 
                                              DenseMatrix<double> &eigen_vect)
 const
{
#ifdef PARANOID
 // Check Matrix is square
 if (N!=M)
  {
   throw OomphLibError(
    "This matrix is not square, the matrix MUST be square!",
    "DenseDoubleMatrix::eigenvalues_by_jacobi()",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // Make a copy of the matrix & check that it's symmetric

 // Check that the sizes of eigen_vals and eigen_vect are correct. If not 
 // correct them.
 if (eigen_vals.size()!=N) 
  { 
   eigen_vals.resize(N); 
  }
 if (eigen_vect.ncol()!=N || eigen_vect.nrow()!=N) 
  { 
   eigen_vect.resize(N); 
  }
 
 DenseDoubleMatrix working_matrix(N);
 for (unsigned long i=0;i<N;i++)
  {
   for (unsigned long j=0;j<M;j++)
    {
#ifdef PARANOID
     if (Matrixdata[M*i+j]!=Matrixdata[M*j+i])
      {
       throw OomphLibError(
        "Matrix needs to be symmetric for eigenvalues_by_jacobi()",
        "DenseDoubleMatrix::eigenvalues_by_jacobi()",
        OOMPH_EXCEPTION_LOCATION);
      }
#endif
     working_matrix(i,j)=(*this)(i,j);
    }
  }
 
 DenseDoubleMatrix aux_eigen_vect(N);
 
 // Call Numerical recipies 
 unsigned long nrot;
 JacobiEigenSolver::jacobi(working_matrix, eigen_vals, aux_eigen_vect, 
                           nrot);
 // Copy across (and transpose)
 for (unsigned long i=0;i<N;i++)
  {
   for (unsigned long j=0;j<M;j++)
    {
     eigen_vect(i,j)=aux_eigen_vect(j,i);
    }
  }
}

/*
//============================================================================
///  Find the residual of Ax=b, i.e. r=b-Ax
//============================================================================
void DenseDoubleMatrix::residual(const DoubleVector &x, 
                                 const DoubleVector &rhs, 
                                 DoubleVector &residual)
{
 std::cout << "broken" << std::endl;
 assert(false);

#ifdef PARANOID
 // Check Matrix is square
 if (N!=M)
  {
   throw OomphLibError(
    "This matrix is not square, the matrix MUST be square!",
    "DenseDoubleMatrix::residual()",
    OOMPH_EXCEPTION_LOCATION);
  }
 // Check that size of rhs = nrow() 
 if (rhs.size()!=N)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The rhs vector is not the right size. It is " << rhs.size() 
    << ", it should be " << N << std::endl;
   
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::residual()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Check that the size of x is correct
 if (x.size()!=N)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector is not the right size. It is " << x.size() 
    << ", it should be " << N << std::endl;
   
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::residual()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 // If size of residual is wrong, correct it!
 if (residual.size()!=N) 
  { 
   residual.resize(N);
  }

 // Multiply the matrix by the vector x in residual vector
 for (unsigned long i=0;i<N;i++)
  {
   residual[i]=rhs[i]; 
   for (unsigned long j=0;j<M;j++)
    {
     residual[i] -= Matrixdata[M*i+j]*x[j];
    }
  }
}
*/



//============================================================================
///  Multiply the matrix by the vector x: soln=Ax
//============================================================================
void DenseDoubleMatrix::multiply(const DoubleVector &x, DoubleVector &soln)
{
#ifdef PARANOID
 // Check to see if x.size() = ncol().
 if (x.nrow()!=M)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector is not the right size. It is " << x.nrow() 
    << ", it should be " << M << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that x is not distributed
 if (x.distributed())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector cannot be distributed for DenseDoubleMatrix "
    << "matrix-vector multiply" << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "DenseDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(x.distribution_pt());
  }

 // Multiply the matrix A, by the vector x 
 for (unsigned long i=0;i<N;i++)
  {
   soln[i] = 0.0;
   for (unsigned long j=0;j<M;j++)
    {
     soln[i] += Matrixdata[M*i+j]*x[j];
    }
  }
}


//=================================================================
/// Multiply the transposed matrix by the vector x: soln=A^T x
//=================================================================
void DenseDoubleMatrix::multiply_transpose(const DoubleVector &x, 
                                        DoubleVector &soln)
{
#ifdef PARANOID
 // Check to see if x.size() = ncol().
 if (x.nrow()!=M)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector is not the right size. It is " << x.nrow() 
    << ", it should be " << M << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that x is not distributed
 if (x.distributed())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector cannot be distributed for DenseDoubleMatrix "
    << "matrix-vector multiply" << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "DenseDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "DenseDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif
 
 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(x.distribution_pt());
  }

 // Initialise the solution
 for (unsigned long i=0;i<M;i++)
  {  
   soln[i] = 0.0;
  }

 // Matrix vector product
 for (unsigned long i=0;i<N;i++)
  {  
   for (unsigned long j=0;j<M;j++)
    {
     soln[j] += Matrixdata[N*i+j]*x[i];
    }
  }
}



//=================================================================
/// For every row, find the maximum absolute value of the
/// entries in this row. Set all values that are less than alpha times
/// this maximum to zero and return the resulting matrix in
/// reduced_matrix. Note: Diagonal entries are retained regardless
/// of their size. 
//=================================================================
void DenseDoubleMatrix::matrix_reduction(const double &alpha,
                                       DenseDoubleMatrix &reduced_matrix)
{

 reduced_matrix.resize(N,M,0.0);
 // maximum value in a row
 double max_row;
  
 // Loop over rows
 for(unsigned i=0;i<N;i++)
  { 
   // Initialise max value in row
   max_row=0.0;
   
   //Loop over entries in columns
   for(unsigned long j=0;j<M;j++)
    {
     // Find max. value in row
     if(std::abs( Matrixdata[M*i+j])>max_row)
      {
       max_row=std::abs( Matrixdata[M*i+j]);
      }
    }

   // Decide if we need to retain the entries in the row
   for(unsigned long j=0;j<M;j++)
    {
     // If we're on the diagonal or the value is sufficiently large: retain
     // i.e. copy across.
     if(i==j || std::abs(Matrixdata[M*i+j])>alpha*max_row )
      {
       reduced_matrix(i,j) =Matrixdata[M*i+j];
      }
    }
  }
 
}


//=============================================================================
/// Function to multiply this matrix by the DenseDoubleMatrix  matrix_in.
//=============================================================================
void DenseDoubleMatrix::multiply(const DenseDoubleMatrix &matrix_in,
                                 DenseDoubleMatrix& result)
{
#ifdef PARANOID
 // check matrix dimensions are compatable 
 if ( this->ncol() != matrix_in.nrow()  )
  {
   std::ostringstream error_message;
   error_message
    << "Matrix dimensions incompatable for matrix-matrix multiplication"
    << "ncol() for first matrix:" << this->ncol()
    << "nrow() for second matrix: " << matrix_in.nrow();
   
   throw OomphLibError(error_message.str(),
                       "DenseDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // NB N is number of rows!
 unsigned long n_row = this->nrow();
 unsigned long m_col = matrix_in.ncol();
 
 // resize and intialize result
 result.resize(n_row, m_col, 0.0);
 
 //clock_t clock1 = clock();
 
 // do calculation
 unsigned long n_col=this->ncol();
 for (unsigned long k=0; k<n_col; k++)
  {
   for (unsigned long i=0; i<n_row; i++)
    {
     for (unsigned long j=0; j<m_col; j++)
      {
       result(i,j) += Matrixdata[m_col*i+k] * matrix_in(k,j);
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=======================================================================
/// \short Default constructor, set the default linear solver and 
/// matrix-matrix multiplication method.
//========================================================================
CCDoubleMatrix::CCDoubleMatrix() : CCMatrix<double>()
  {
   Linear_solver_pt = Default_linear_solver_pt = new SuperLU;
   Matrix_matrix_multiply_method = 2;
  }
  
//========================================================================
 /// \short Constructor: Pass vector of values, vector of row indices,
 /// vector of column starts and number of rows (can be suppressed
 /// for square matrices). Number of nonzero entries is read
 /// off from value, so make sure the vector has been shrunk
 /// to its correct length.
//=======================================================================
 CCDoubleMatrix::CCDoubleMatrix(const Vector<double>& value,
                                const Vector<int>& row_index,
                                const Vector<int>& column_start,
                                const unsigned long &n,
                                const unsigned long &m) :
  CCMatrix<double>(value,row_index,column_start,n,m)
  {
   Linear_solver_pt = Default_linear_solver_pt = new SuperLU;
   Matrix_matrix_multiply_method = 2;
  }

 /// Destructor: delete the default linaer solver
 CCDoubleMatrix::~CCDoubleMatrix() {delete Default_linear_solver_pt;}


//===================================================================
/// Perform LU decomposition. Return the sign of the determinant
//===================================================================
void CCDoubleMatrix::ludecompose()
{
 static_cast<SuperLU*>(Default_linear_solver_pt)->factorise(this);
}

//===================================================================
/// Do the backsubstitution
//===================================================================
void CCDoubleMatrix::lubksub(DoubleVector &rhs)
{
 static_cast<SuperLU*>(Default_linear_solver_pt)->backsub(rhs,rhs);
}

//===================================================================
///  Multiply the matrix by the vector x
//===================================================================
void CCDoubleMatrix::multiply(const DoubleVector &x, DoubleVector &soln)
{
#ifdef PARANOID
 // Check to see if x.size() = ncol().
 if (x.nrow()!=M)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector is not the right size. It is " << x.nrow() 
    << ", it should be " << M << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that x is not distributed
 if (x.distributed())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector cannot be distributed for CCDoubleMatrix "
    << "matrix-vector multiply" << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "CCDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(x.distribution_pt());
  }

 for (unsigned i=0;i<N;i++)
  {
   soln[i] = 0.0;
  }
 
 for (unsigned long j=0;j<N;j++)
  {
   for (long k=Column_start[j];k<Column_start[j+1];k++)
    {
     unsigned long i = Row_index[k];
     double a_ij = Value[k];
     soln[i] += a_ij*x[j];
    }
  }
}




//=================================================================
/// Multiply the  transposed matrix by the vector x: soln=A^T x
//=================================================================
void CCDoubleMatrix::multiply_transpose(const DoubleVector &x, 
                                        DoubleVector &soln)
{
#ifdef PARANOID
 // Check to see if x.size() = ncol().
 if (x.nrow()!=M)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector is not the right size. It is " << x.nrow() 
    << ", it should be " << M << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that x is not distributed
 if (x.distributed())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector cannot be distributed for CCDoubleMatrix "
    << "matrix-vector multiply" << std::endl;
   throw OomphLibError(error_message_stream.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "CCDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif
 
 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(x.distribution_pt());
  }

 // Initialise the solution
 for (unsigned long i=0;i<M;i++)
  {  
   soln[i] = 0.0;
  }

 // Matrix vector product
 for (unsigned long i=0;i<N;i++)
  {  
   
   for (long k=Column_start[i];k<Column_start[i+1];k++)
    {
     unsigned long j=Row_index[k];
     double a_ij=Value[k];
     soln[j]+=a_ij*x[i];
    }
  }
}




//===========================================================================
/// Function to multiply this matrix by the CRDoubleMatrix matrix_in
/// The multiplication method used can be selected using the flag
/// Matrix_matrix_multiply_method. By default Method 2 is used.
/// Method 1: First runs through this matrix and matrix_in to find the storage
///           requirements for result - arrays of the correct size are 
///           then allocated before performing the calculation.
///           Minimises memory requirements but more costly.
/// Method 2: Grows storage for values and column indices of result 'on the
///           fly' using an array of maps. Faster but more memory
///           intensive.
/// Method 3: Grows storage for values and column indices of result 'on the
///           fly' using a vector of vectors. Not particularly impressive
///           on the platforms we tried...
//=============================================================================
void CCDoubleMatrix::multiply(const CCDoubleMatrix& matrix_in,
                              CCDoubleMatrix& result)
{
#ifdef PARANOID
 // check matrix dimensions are compatable
 if ( this->ncol() != matrix_in.nrow()  )
  {
   std::ostringstream error_message;
   error_message 
    << "Matrix dimensions incompatable for matrix-matrix multiplication"
    << "ncol() for first matrix:" << this->ncol()
    << "nrow() for second matrix: " << matrix_in.nrow();
   
   throw OomphLibError(error_message.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 // NB N is number of rows!
 unsigned long N = this->nrow();
 unsigned long M = matrix_in.ncol();
 unsigned long Nnz = 0;
 
 // pointers to arrays which store result
 int* Column_start;
 double* Value;
 int* Row_index;

 // get pointers to matrix_in
 const int* matrix_in_col_start = matrix_in.column_start();
 const int* matrix_in_row_index = matrix_in.row_index();
 const double* matrix_in_value = matrix_in.value();

 // get pointers to this matrix
 const double* this_value = this->value();
 const int* this_col_start = this->column_start();
 const int* this_row_index = this->row_index();

 // set method
 unsigned method = Matrix_matrix_multiply_method;

 // clock_t clock1 = clock();

 // METHOD 1
 // --------
 if (method==1)
 {
  // allocate storage for column starts
  Column_start = new int[M+1];
  Column_start[0]=0;

  // a set to store number of non-zero rows in each column of result
  std::set<unsigned> rows;

  // run through columns of this matrix and matrix_in to find number of
  // non-zero entries in each column of result
  for (unsigned long this_col = 0; this_col<M; this_col++)
  {
   // run through non-zeros in this_col of this matrix
   for (int this_ptr = this_col_start[this_col];
        this_ptr < this_col_start[this_col+1];
        this_ptr++)
   {
    // find row index for non-zero
    unsigned matrix_in_col = this_row_index[this_ptr];

    // run through corresponding column in matrix_in
    for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
         matrix_in_ptr < matrix_in_col_start[matrix_in_col+1];
         matrix_in_ptr++)
    {
     // find row index for non-zero in matrix_in and store in rows
     rows.insert(matrix_in_row_index[matrix_in_ptr]);
    }
   }
   // update Column_start
   Column_start[this_col+1] = Column_start[this_col] + rows.size();

   // wipe values in rows
   rows.clear();
  }

  // set Nnz
  Nnz = Column_start[M];

  // allocate arrays for result
  Value = new double[Nnz];
  Row_index = new int[Nnz];

  // set all values of Row_index to -1
  for (unsigned long i=0;i<Nnz;i++)
   Row_index[i] = -1;

  // Calculate values for result - first run through columns of this matrix
  for (unsigned long this_col = 0; this_col<M; this_col++)
   {
    // run through non-zeros in this_column
    for (int this_ptr = this_col_start[this_col];
         this_ptr < this_col_start[this_col+1];
         this_ptr++)
     {
      // find value of non-zero
      double this_val = this_value[this_ptr];
      
      // find row associated with non-zero
      unsigned matrix_in_col = this_row_index[this_ptr];
      
      // run through corresponding column in matrix_in
      for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
           matrix_in_ptr < matrix_in_col_start[matrix_in_col+1];
           matrix_in_ptr++)
       {
        // find row index for non-zero in matrix_in
        int row = matrix_in_row_index[matrix_in_ptr];
        
        // find position in result to insert value
        for(int ptr = Column_start[this_col];
            ptr <= Column_start[this_col+1];
            ptr++)
         {
          if (ptr == Column_start[this_col+1])
           {
            // error - have passed end of column without finding
            // correct row index
            std::ostringstream error_message;
            error_message << "Error inserting value in result";
            
            throw OomphLibError(error_message.str(),
                                "CCDoubleMatrix::multiply()",
                                OOMPH_EXCEPTION_LOCATION);
           }
          else if (Row_index[ptr] == -1 )
           {
            // first entry for this row index
            Row_index[ptr] = row;
            Value[ptr] = this_val * matrix_in_value[matrix_in_ptr];
            break;
           }
          else if ( Row_index[ptr] == row )
           {
            // row index already exists - add value
            Value[ptr] += this_val * matrix_in_value[matrix_in_ptr];
            break;
           }
         }
       }
     }
   }
 }
 
 // METHOD 2
 // --------
 else if (method==2)
  {
   // generate array of maps to store values for result
   std::map<int,double>* result_maps = new std::map<int,double>[M];
   
   // run through columns of this matrix
   for (unsigned long this_col = 0; this_col<M; this_col++)
    {
     // run through non-zeros in this_col
     for (int this_ptr = this_col_start[this_col];
          this_ptr < this_col_start[this_col+1];
          this_ptr++)
      {
       // find value of non-zero
       double this_val = this_value[this_ptr];
       
       // find row index associated with non-zero
       unsigned matrix_in_col = this_row_index[this_ptr];
       
       // run through corresponding column in matrix_in
       for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
            matrix_in_ptr < matrix_in_col_start[matrix_in_col+1];
            matrix_in_ptr++)
        {
         // find row index for non-zero in matrix_in
         int row = matrix_in_row_index[matrix_in_ptr];
         
         // insert value
         result_maps[this_col][row] += 
          this_val * matrix_in_value[matrix_in_ptr];
        }
      }
    }
   
   // allocate Column_start
   Column_start = new int[M+1];
   
   // copy across column starts
   Column_start[0] = 0;
   for (unsigned long col=0; col<M; col++)
    {
     int size = result_maps[col].size();
     Column_start[col+1] = Column_start[col] + size;
    }
   
   // set Nnz
   Nnz = Column_start[M];
   
   // allocate other arrays
   Value = new double[Nnz];
   Row_index = new int[Nnz];
   
   // copy values and row indices
   for (unsigned long col=0; col<M; col++)
    {
     unsigned ptr = Column_start[col];
     for (std::map<int,double>::iterator i = result_maps[col].begin();
          i != result_maps[col].end();
          i ++)
      {
       Row_index[ptr]= i->first;
       Value[ptr] = i->second;
       ptr++;
      }
    }
   
   // tidy up memory
   delete[] result_maps;
  }
 
 // METHOD 3
 // --------
 else if (method==3)
  {
   // vectors of vectors to store results
   std::vector< std::vector<int> > result_rows(N);
   std::vector< std::vector<double> > result_vals(N);
   
   // run through the columns of this matrix
  for (unsigned long this_col = 0; this_col<M; this_col++)
   {
    // run through non-zeros in this_col
    for (int this_ptr = this_col_start[this_col];
         this_ptr < this_col_start[this_col+1];
         this_ptr++)
     {
      // find value of non-zero
      double this_val = this_value[this_ptr];
      
      // find row index associated with non-zero
      unsigned matrix_in_col = this_row_index[this_ptr];
      
      // run through corresponding column in matrix_in
      for (int matrix_in_ptr = matrix_in_col_start[matrix_in_col];
           matrix_in_ptr < matrix_in_col_start[matrix_in_col+1];
           matrix_in_ptr++)
       {
        // find row index for non-zero in matrix_in
        int row = matrix_in_row_index[matrix_in_ptr];
        
        // insert value
        int size = result_rows[this_col].size();
        for (int i = 0; i<=size; i++)
         {
          if (i==size)
           {
            // first entry for this row index
            result_rows[this_col].push_back(row);
            result_vals[this_col].push_back(
             this_val*matrix_in_value[matrix_in_ptr]);
           }
          else if (row==result_rows[this_col][i])
           {
            // row index already exists
            result_vals[this_col][i] += 
             this_val * matrix_in_value[matrix_in_ptr];
            break;
           }
         }
       }
     }
   }
  
  // allocate Column_start
  Column_start = new int[M+1];
  
  // copy across column starts
  Column_start[0] = 0;
  for (unsigned long col=0; col<M; col++)
   {
    int size = result_rows[col].size();
    Column_start[col+1] = Column_start[col] + size;
   }
  
  // set Nnz
  Nnz = Column_start[M];
  
  // allocate other arrays
  Value = new double[Nnz];
  Row_index = new int[Nnz];
  
  // copy across values and row indices
  for (unsigned long col=0; col<N; col++)
   {
    unsigned ptr = Column_start[col];
    unsigned n_rows=result_rows[col].size();
    for (unsigned i = 0; i < n_rows ; i++)
     {
      Row_index[ptr] = result_rows[col][i];
      Value[ptr] = result_vals[col][i];
      ptr++;
     }
   }
  }
 
 // INCORRECT VALUE FOR METHOD
 else
  {
   std::ostringstream error_message;
   error_message << "Incorrect method set in matrix-matrix multiply"
                 << "method=" << method << " not allowed";
   
   throw OomphLibError(error_message.str(),
                       "CCDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 result.build_without_copy(Value, Row_index, Column_start, Nnz, N, M);
}



//=================================================================
/// For every row, find the maximum absolute value of the
/// entries in this row. Set all values that are less than alpha times
/// this maximum to zero and return the resulting matrix in
/// reduced_matrix. Note: Diagonal entries are retained regardless
/// of their size. 
//=================================================================
void CCDoubleMatrix::matrix_reduction(const double &alpha,
                                       CCDoubleMatrix &reduced_matrix)
{
 // number of columns in matrix
 long n_coln=ncol();     

 Vector<double>max_row(nrow(),0.0);

 // Here's the packed format for the new matrix
 Vector<int> B_row_start(1);
 Vector<int> B_column_index;
 Vector<double> B_value;
 

 // k is counter for the number of entries in the reduced matrix
 unsigned k=0;

 // Initialise row start
 B_row_start[0]=0;

 // Loop over columns
 for(long i=0;i<n_coln;i++)
  { 
     
   //Loop over entries in columns
   for(long j=Column_start[i];j<Column_start[i+1];j++)
    {
    
     // Find max. value in row
     if(std::abs(Value[j])>max_row[Row_index[j]])
      {
       max_row[Row_index[j]]=std::abs(Value[j]);
      }
    }

   // Decide if we need to retain the entries in the row
   for(long j=Column_start[i];j<Column_start[i+1];j++)
    {
     // If we're on the diagonal or the value is sufficiently large: retain
     // i.e. copy across.
     if(i==Row_index[j] || std::abs(Value[j])>alpha*max_row[Row_index[j]] )
      {
       B_value.push_back(Value[j]);
       B_column_index.push_back(Row_index[j]);
       k++;
      }
    }
   // This writes the row start for the next row -- equal to 
   // to the number of entries written so far (C++ zero-based indexing!)
   B_row_start.push_back(k);
  }


 // Build the matrix from the compressed format
 dynamic_cast<CCDoubleMatrix&>(reduced_matrix).
  build(B_value,B_column_index,B_row_start,nrow(),ncol());
 }



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/// Default constructor
//=============================================================================
CRDoubleMatrix::CRDoubleMatrix()
  {
    Linear_solver_pt = Default_linear_solver_pt  = 0;
    Built = false;

    // set the serial matrix-matrix multiply method
#ifdef HAVE_TRILINOS
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif
  }

//=============================================================================
/// Constructor: just stores the distribution but does not build the
/// matrix
//=============================================================================
CRDoubleMatrix::CRDoubleMatrix(const LinearAlgebraDistribution* 
                               distribution_pt)
  {
   Distribution_pt->rebuild(distribution_pt);
   Built = false;    

// set the serial matrix-matrix multiply method
#ifdef HAVE_TRILINOS
    Serial_matrix_matrix_multiply_method = 2;
#else
    Serial_matrix_matrix_multiply_method = 2;
#endif

  }

//=============================================================================
/// \short Constructor: Takes the distribution and the number of columns, as 
/// well as the vector of values, vector of column indices,vector of row 
///starts.
//=============================================================================
CRDoubleMatrix::CRDoubleMatrix(const LinearAlgebraDistribution* dist_pt,
                               const unsigned& ncol,
                               const Vector<double>& value, 
                               const Vector<int>& column_index,
                               const Vector<int>& row_start) 
{
 // build the compressed row matrix
 CR_matrix.build(value,column_index,row_start,dist_pt->nrow_local(),ncol);

 // store the Distribution
 Distribution_pt->rebuild(dist_pt);

 // set the linear solver
#ifdef OOMPH_HAS_MPI
 if (Distribution_pt->communicator_pt()->nproc() > 1 && 
     Distribution_pt->distributed() == true)
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU_dist;   
  }
 else
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU;
  }
#else
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU;
#endif

 // set the serial matrix-matrix multiply method
#ifdef HAVE_TRILINOS
 Serial_matrix_matrix_multiply_method = 2;
#else
 Serial_matrix_matrix_multiply_method = 2;
#endif

 // matrix has been built
 Built = true;
}

//=============================================================================
/// rebuild the matrix - assembles an empty matrix with a defined distribution
//=============================================================================
void CRDoubleMatrix::rebuild(const LinearAlgebraDistribution* distribution_pt)
{
 Distribution_pt->rebuild(distribution_pt);
 CR_matrix.clean_up_memory();
 Built = false;
 
 // set the linear solver
#ifdef OOMPH_HAS_MPI
 if (Distribution_pt->communicator_pt()->nproc() > 1 && 
     Distribution_pt->distributed() == true)
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU_dist;   
  }
 else
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU;
  }
#else
 if (Default_linear_solver_pt != 0)
  {
   delete Default_linear_solver_pt;
  }
 Default_linear_solver_pt = new SuperLU;
#endif
}


//=============================================================================
/// clean method
//=============================================================================
void CRDoubleMatrix::clear() 
{
 Distribution_pt->clear();
 CR_matrix.clean_up_memory();
 Built = false;
 delete Default_linear_solver_pt;
 Default_linear_solver_pt = 0;
}

//=============================================================================
/// \short build method: Takes the distribution and the number of columns, as 
/// well as the vector of values, vector of column indices,vector of row 
///starts.
//=============================================================================
void CRDoubleMatrix::rebuild(const LinearAlgebraDistribution* distribution_pt,
                             const unsigned& ncol,
                             const Vector<double>& value, 
                             const Vector<int>& column_index,
                             const Vector<int>& row_start)
{
 // store the Distribution
 Distribution_pt->rebuild(distribution_pt);

 // set the linear solver
#ifdef OOMPH_HAS_MPI
 if (Distribution_pt->communicator_pt()->nproc() > 1 && 
     Distribution_pt->distributed() == true)
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU_dist;   
  }
 else
  {
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU;
  }
#else
   if (Default_linear_solver_pt != 0)
    {
     delete Default_linear_solver_pt;
    }
   Default_linear_solver_pt = new SuperLU;
#endif

 // now build the matrix
 this->rebuild_matrix(ncol,value,column_index,row_start);
}

//=============================================================================
/// \short method to rebuild the matrix, but not the distribution
//=============================================================================
void CRDoubleMatrix::rebuild_matrix(const unsigned& ncol,
                                    const Vector<double>& value,
                                    const Vector<int>& column_index,
                                    const Vector<int>& row_start)
{
 // call the underlying build method
 CR_matrix.clean_up_memory();
 CR_matrix.build(value,column_index,row_start,
                 Distribution_pt->nrow_local(),ncol);

 // matrix has been build
 Built = true;
}

//=============================================================================
/// \short method to rebuild the matrix, but not the distribution
//=============================================================================
void CRDoubleMatrix::rebuild_matrix_without_copy(const unsigned& ncol,
                                                 const unsigned& nnz,
                                                 double* value,
                                                 int* column_index,
                                                 int* row_start)
{
 // call the underlying build method
 CR_matrix.clean_up_memory();
 CR_matrix.build_without_copy(value,column_index,row_start,nnz,
                              Distribution_pt->nrow_local(),ncol);

 // matrix has been build
 Built = true;
}

//=============================================================================
/// Do LU decomposition
//=============================================================================
void CRDoubleMatrix::ludecompose()
{
#ifdef PARANOID
 // check that the this matrix is built
 if (!Built)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix has not been built.";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::ludecompose()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // factorise using superlu or superlu dist if we oomph has mpi
#ifdef OOMPH_HAS_MPI
 if (static_cast<SuperLU* >(Default_linear_solver_pt))
  {
#endif
   static_cast<SuperLU*>(Default_linear_solver_pt)->factorise(this);
#ifdef OOMPH_HAS_MPI
  }
 else
  {
   static_cast<SuperLU_dist*>(Default_linear_solver_pt)->factorise(this);
  }  
#endif
}

//=============================================================================
/// Do back-substitution
//=============================================================================
void CRDoubleMatrix::lubksub(DoubleVector &rhs)
{
#ifdef PARANOID
 // check that the rhs vector is setup
 if (!rhs.distribution_setup())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The vector rhs has not been setup";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::lubksub()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that the rhs vector has the same distribution as this matrix
 if (!(*Distribution_pt == *rhs.distribution_pt()))
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The vector rhs must have the same distribution as the matrix";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::lubksup()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 if (static_cast<SuperLU* >(Default_linear_solver_pt))
  {
   DoubleVector rhs_copy(rhs);
   static_cast<SuperLU*>(Default_linear_solver_pt)->backsub(rhs_copy,rhs);
  }
 else
  {
   DoubleVector rhs_copy(rhs);
   static_cast<SuperLU*>(Default_linear_solver_pt)->backsub(rhs_copy,rhs);
  }  
}

//=============================================================================
///  Multiply the matrix by the vector x
//=============================================================================
void CRDoubleMatrix::multiply(const DoubleVector &x, DoubleVector &soln)
{
#ifdef PARANOID
 // check that this matrix is built
 if (!Built)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that the distribution of x is setup
 if (!x.distribution_setup())
  {
 std::ostringstream error_message_stream;
   error_message_stream 
    << "The distribution of the vector x must be setup";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Check to see if x.size() = ncol().
 if (this->ncol() != x.distribution_pt()->nrow())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The number of rows in the x vector and the number of columns in the "
    << "matrix must be the same";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "CRDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(this->distribution_pt());
  }

#ifdef HAVE_TRILINOS
 // This will only work if we have trilinos on board
 TrilinosHelpers::multiply(*this,x,soln);
#else
 if (this->distributed() && 
     this->distribution_pt()->communicator_pt()->nproc() > 1)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "Matrix-vector product on multiple processors with distributed "
    << "CRDoubleMatrix requires Trilinos.";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 else
  {   
   unsigned n = this->nrow();
   const int* row_start = CR_matrix.row_start();
   const int* column_index = CR_matrix.column_index();
   const double* value = CR_matrix.value();
   for (unsigned long i=0;i<n;i++)
    {  
     soln[i] = 0.0;
     for (long k=row_start[i];k<row_start[i+1];k++)
      {
       unsigned long j=column_index[k];
       double a_ij=value[k];
       soln[i]+=a_ij*x[j];
      }
    }
  }
#endif
}



//=================================================================
/// Multiply the  transposed matrix by the vector x: soln=A^T x
//=================================================================
void CRDoubleMatrix::multiply_transpose(const DoubleVector &x, 
                                        DoubleVector &soln)
{
#ifdef PARANOID
 // check that this matrix is built
 if (!Built)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply_transpose()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // Check to see if x.size() = ncol().
 if (!(*Distribution_pt == *x.distribution_pt()))
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "The x vector and this matrix must have the same distribution.";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply_transpose()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (soln.distribution_setup())
  {
   if (!(*soln.distribution_pt() == *x.distribution_pt()))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The soln vector is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "CRDoubleMatrix::multiply_transpose()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // if soln is not setup then setup the distribution
 if (!soln.distribution_setup())
  {
   // Resize and initialize the solution vector
   soln.rebuild(x.distribution_pt());
  }

 if (this->distributed() && 
     this->distribution_pt()->communicator_pt()->nproc() > 1)
  {
#ifdef HAVE_TRILINOS
   // This will only work if we have trilinos on board
   TrilinosHelpers::multiply(*this,x,soln);
#else
     std::ostringstream error_message_stream;
     error_message_stream 
      << "Matrix-vector product on multiple processors with distributed "
      << "CRDoubleMatrix requires Trilinos.";
     throw OomphLibError(error_message_stream.str(),
                         "CRDoubleMatrix::multiply_transpose()",
                         OOMPH_EXCEPTION_LOCATION);
#endif
  }
 else
  {
   unsigned n = this->nrow();
   const int* row_start = CR_matrix.row_start();
   const int* column_index = CR_matrix.column_index();
   const double* value = CR_matrix.value();
   // Matrix vector product
   for (unsigned long i=0;i<n;i++)
    {  
     for (long k=row_start[i];k<row_start[i+1];k++)
      {
       unsigned long j=column_index[k];
       double a_ij=value[k];
       soln[j]+=a_ij*x[i];
      }
    }
  }
}


//===========================================================================
/// Function to multiply this matrix by the CRDoubleMatrix matrix_in
/// The multiplication method used can be selected using the flag
/// Matrix_matrix_multiply_method. By default Method 2 is used.
/// Method 1: First runs through this matrix and matrix_in to find the storage
///           requirements for result - arrays of the correct size are 
///           then allocated before performing the calculation.
///           Minimises memory requirements but more costly.
/// Method 2: Grows storage for values and column indices of result 'on the
///           fly' using an array of maps. Faster but more memory
///           intensive.
/// Method 3: Grows storage for values and column indices of result 'on the
///           fly' using a vector of vectors. Not particularly impressive
///           on the platforms we tried...
/// Method 4: Trilinos Epetra Matrix Matrix multiply.\n
/// Method 5: Trilinox Epetra Matrix Matrix Mulitply (ml based) \m
//=============================================================================
void CRDoubleMatrix::multiply(CRDoubleMatrix& matrix_in,
                              CRDoubleMatrix& result)
{
#ifdef PARANOID
 // check that this matrix is built
 if (!Built)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // check that this matrix is built
 if (!matrix_in.built())
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "This matrix matrix_in has not been built";
   throw OomphLibError(error_message_stream.str(),
                       "CRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 // if soln is setup then it should have the same distribution as x
 if (result.distribution_setup())
  {
   if (!(*result.distribution_pt() == *Distribution_pt))
    {
     std::ostringstream error_message_stream;
     error_message_stream 
      << "The matrix result is setup and therefore must have the same "
      << "distribution as the vector x";
     throw OomphLibError(error_message_stream.str(),
                         "CRDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

 // if the result has not been setup, then store the distribution
 if (!result.distribution_setup())
  {
   result.rebuild(Distribution_pt);
  }
   
 // short name for Serial_matrix_matrix_multiply_method
 unsigned method = Serial_matrix_matrix_multiply_method;

 // if this matrix is not distributed and matrix in is not distributed
 if (!this->distributed() && !matrix_in.distributed() && 
     ((method == 1) || (method == 2) || (method == 3)))
  {
   // NB N is number of rows!
   unsigned long N = this->nrow();
   unsigned long M = matrix_in.ncol();
   unsigned long Nnz = 0;
   
   // pointers to arrays which store result
   int* Row_start = 0;
   double* Value = 0;
   int* Column_index = 0;
   
   // get pointers to matrix_in
   const int* matrix_in_row_start = matrix_in.row_start();
   const int* matrix_in_column_index = matrix_in.column_index();
   const double* matrix_in_value = matrix_in.value();
   
   // get pointers to this matrix
   const double* this_value = this->value();
   const int* this_row_start = this->row_start();
   const int* this_column_index = this->column_index();
   
   //clock_t clock1 = clock();
   
   // METHOD 1
   // --------
   if (method==1)
    {
     // allocate storage for row starts
     Row_start = new int[N+1];
     Row_start[0]=0;
     
     // a set to store number of non-zero columns in each row of result
     std::set<unsigned> columns;
     
     // run through rows of this matrix and matrix_in to find number of
     // non-zero entries in each row of result
     for (unsigned long this_row = 0; this_row<N; this_row++)
      {
       // run through non-zeros in this_row of this matrix
       for (int this_ptr = this_row_start[this_row];
            this_ptr < this_row_start[this_row+1];
            this_ptr++)
        {
         // find column index for non-zero
         int matrix_in_row = this_column_index[this_ptr];
         
         // run through corresponding row in matrix_in
         for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
              matrix_in_ptr < matrix_in_row_start[matrix_in_row+1];
              matrix_in_ptr++)
          {
           // find column index for non-zero in matrix_in and store in columns
           columns.insert(matrix_in_column_index[matrix_in_ptr]);
          }
        }
       // update Row_start
       Row_start[this_row+1] = Row_start[this_row] + columns.size();
       
       // wipe values in columns
       columns.clear();
      }
     
     // set Nnz
     Nnz = Row_start[N];
     
     // allocate arrays for result
     Value = new double[Nnz];
     Column_index = new int[Nnz];
     
     // set all values of Column_index to -1
     for (unsigned long i=0; i<Nnz; i++)
      {
       Column_index[i] = -1;
      }
     
     // Calculate values for result - first run through rows of this matrix
     for (unsigned long this_row = 0; this_row<N; this_row++)
      {
       // run through non-zeros in this_row
       for (int this_ptr = this_row_start[this_row];
            this_ptr < this_row_start[this_row+1];
            this_ptr++)
        {
         // find value of non-zero
         double this_val = this_value[this_ptr];
         
         // find column associated with non-zero
         int matrix_in_row = this_column_index[this_ptr];
         
         // run through corresponding row in matrix_in
         for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
              matrix_in_ptr < matrix_in_row_start[matrix_in_row+1];
              matrix_in_ptr++)
          {
           // find column index for non-zero in matrix_in
           int col = matrix_in_column_index[matrix_in_ptr];
           
           // find position in result to insert value
           for(int ptr = Row_start[this_row];
               ptr <= Row_start[this_row+1];
               ptr++)
            {
             if (ptr == Row_start[this_row+1])
              {
               // error - have passed end of row without finding
               // correct column
               std::ostringstream error_message;
               error_message << "Error inserting value in result";
               
               throw OomphLibError(error_message.str(),
                           "CRDoubleMatrix::multiply()",
                                   OOMPH_EXCEPTION_LOCATION);
              }
             else if (	Column_index[ptr] == -1 )
              {
               // first entry for this column index
               Column_index[ptr] = col;
               Value[ptr] = this_val * matrix_in_value[matrix_in_ptr];
               break;
              }
             else if ( Column_index[ptr] == col )
              {
               // column index already exists - add value
               Value[ptr] += this_val * matrix_in_value[matrix_in_ptr];
               break;
              }
            }
          }
        }
      }
    }
   
   // METHOD 2
   // --------
   else if (method==2)
    {
     // generate array of maps to store values for result
     std::map<int,double>* result_maps = new std::map<int,double>[N];
     
     // run through rows of this matrix
     for (unsigned long this_row = 0; this_row<N; this_row++)
      {
       // run through non-zeros in this_row
       for (int this_ptr = this_row_start[this_row];
            this_ptr < this_row_start[this_row+1];
            this_ptr++)
        {
         // find value of non-zero
         double this_val = this_value[this_ptr];
         
         // find column index associated with non-zero
         int matrix_in_row = this_column_index[this_ptr];
         
         // run through corresponding row in matrix_in
         for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
              matrix_in_ptr < matrix_in_row_start[matrix_in_row+1];
              matrix_in_ptr++)
          {
           // find column index for non-zero in matrix_in
           int col = matrix_in_column_index[matrix_in_ptr];
           
           // insert value
           result_maps[this_row][col] += this_val * matrix_in_value[matrix_in_ptr];
          }
        }
      }
     
     // allocate Row_start
     Row_start = new int[N+1];
     
     // copy across row starts
     Row_start[0] = 0;
     for (unsigned long row=0; row<N; row++)
      {
       int size = result_maps[row].size();
       Row_start[row+1] = Row_start[row] + size;
      }
     
     // set Nnz
     Nnz = Row_start[N];
     
     // allocate other arrays
     Value = new double[Nnz];
     Column_index = new int[Nnz];
     
     // copy values and column indices
     for (unsigned long row=0; row<N; row++)
      {
       unsigned ptr = Row_start[row];
       for (std::map<int,double>::iterator i = result_maps[row].begin();
            i != result_maps[row].end();
            i ++)
        {
         Column_index[ptr]= i->first;
         Value[ptr] = i->second;
         ptr++;
        }
      }
     
     // tidy up memory
     delete[] result_maps;
    }
 
   // METHOD 3
   // --------
   else if (method==3)
    {
     // vectors of vectors to store results
     std::vector< std::vector<int> > result_cols(N);
     std::vector< std::vector<double> > result_vals(N);
     
     // run through the rows of this matrix
     for (unsigned long this_row = 0; this_row<N; this_row++)
      {
       // run through non-zeros in this_row
       for (int this_ptr = this_row_start[this_row];
            this_ptr < this_row_start[this_row+1];
            this_ptr++)
        {
         // find value of non-zero
         double this_val = this_value[this_ptr];
         
         // find column index associated with non-zero
         int matrix_in_row = this_column_index[this_ptr];
         
         // run through corresponding row in matrix_in
         for (int matrix_in_ptr = matrix_in_row_start[matrix_in_row];
            matrix_in_ptr < matrix_in_row_start[matrix_in_row+1];
              matrix_in_ptr++)
          {
           // find column index for non-zero in matrix_in
           int col = matrix_in_column_index[matrix_in_ptr];
           
           // insert value
           int size = result_cols[this_row].size();
           for (int i = 0; i<=size; i++)
            {
             if (i==size)
              {
               // first entry for this column
               result_cols[this_row].push_back(col);
               result_vals[this_row].push_back(
                this_val*matrix_in_value[matrix_in_ptr]);
              }
             else if (col==result_cols[this_row][i])
              {
               // column already exists
               result_vals[this_row][i] += this_val * 
                matrix_in_value[matrix_in_ptr];
               break;
              }
            }
          }
        }
      }
     
     // allocate Row_start
     Row_start = new int[N+1];
     
     // copy across row starts
     Row_start[0] = 0;
     for (unsigned long row=0; row<N; row++)
      {
       int size = result_cols[row].size();
       Row_start[row+1] = Row_start[row] + size;
      }
     
     // set Nnz
     Nnz = Row_start[N];
     
     // allocate other arrays
     Value = new double[Nnz];
     Column_index = new int[Nnz];
     
     // copy across values and column indices
     for (unsigned long row=0; row<N; row++)
      {
       unsigned ptr = Row_start[row];
       unsigned nnn=result_cols[row].size();
       for (unsigned i = 0; i < nnn; i++) 
        {
         Column_index[ptr] = result_cols[row][i];
         Value[ptr] = result_vals[row][i];
         ptr++;
        }
      }
    }
   
   // build
   result.rebuild_matrix_without_copy(M, Nnz, Value, Column_index, Row_start);
  }
  
 // else we have to use trilinos
 else
  {
#ifdef HAVE_TRILINOS
     bool use_ml = false;
     if (method == 5)
      {
       use_ml = true;
      }
     TrilinosHelpers::multiply(*this,matrix_in,result,use_ml);
#else
     std::ostringstream error_message;
     error_message << "Serial_matrix_matrix_multiply_method = "
                   << Serial_matrix_matrix_multiply_method 
                   << " requires trilinos.";
     throw OomphLibError(error_message.str(),
                         "CRDoubleMatrix::multiply()",
                         OOMPH_EXCEPTION_LOCATION);
#endif
  }
}
  
 

//=================================================================
/// For every row, find the maximum absolute value of the
/// entries in this row. Set all values that are less than alpha times
/// this maximum to zero and return the resulting matrix in
/// reduced_matrix. Note: Diagonal entries are retained regardless
/// of their size. 
//=================================================================
void CRDoubleMatrix::matrix_reduction(const double &alpha,
                                       CRDoubleMatrix &reduced_matrix)
{
 // number of rows in matrix
 long n_row=nrow_local();     
 double max_row;
 
 // Here's the packed format for the new matrix
 Vector<int> B_row_start(1);
 Vector<int> B_column_index;
 Vector<double> B_value;
 
 // get pointers to the underlying data
 const int* row_start = CR_matrix.row_start();
 const int* column_index = CR_matrix.column_index();
 const double* value = CR_matrix.value();
 
 // k is counter for the number of entries in the reduced matrix
 unsigned k=0;

 // Initialise row start
 B_row_start[0]=0;

 // Loop over rows
 for(long i=0;i<n_row;i++)
  { 
   // Initialise max value in row
   max_row=0.0;
   
   //Loop over entries in columns
   for(long j=row_start[i];j<row_start[i+1];j++)
    {
     // Find max. value in row
     if(std::abs(value[j])>max_row)
      {
       max_row=std::abs(value[j]);
      }
    }

   // Decide if we need to retain the entries in the row
   for(long j=row_start[i];j<row_start[i+1];j++)
    {
     // If we're on the diagonal or the value is sufficiently large: retain
     // i.e. copy across.
     if(i==column_index[j] || std::abs(value[j])>alpha*max_row )
      {
       B_value.push_back(value[j]);
       B_column_index.push_back(column_index[j]);
       k++;
      }
    }
   // This writes the row start for the next row -- equal to 
   // to the number of entries written so far (C++ zero-based indexing!)
   B_row_start.push_back(k);
  }
 
 // Build the matrix from the compressed format
 dynamic_cast<CRDoubleMatrix&>(reduced_matrix).
  rebuild_matrix(this->ncol(),B_value,B_column_index,B_row_start);
 }

//=============================================================================
/// if this matrix is distributed then a the equivalent global matrix is built
/// using new and returned. The calling method is responsible for the 
/// destruction of the new matrix.
//=============================================================================
CRDoubleMatrix* CRDoubleMatrix::return_global_matrix()
{
#ifdef OOMPH_HAS_MPI
 // if this matrix is not distributed then this method is redundant
 if (!this->distributed() || 
     this->distribution_pt()->communicator_pt()->nproc() == 1)
  {
   return new CRDoubleMatrix(*this);
  }

 // nnz 
 int nnz = this->nnz();

 // my nrow local
 unsigned nrow_local = this->nrow_local();

 // nrow global
 unsigned nrow = this->nrow();
   
 // cache nproc
 int nproc = Distribution_pt->communicator_pt()->nproc();

 // get the nnzs on the other processors
 int* dist_nnz_pt = new int[nproc];
 MPI_Allgather(&nnz,1,MPI_INT,
               dist_nnz_pt,1,MPI_INT,
               Distribution_pt->communicator_pt()->mpi_comm());
   
 // create a int vector of first rows and nrow local and compute nnz global
 int* dist_first_row = new int[nproc];
 int* dist_nrow_local =  new int[nproc];
 int nnz_global = 0;
 for (int p = 0; p < MPI_Helpers::Nproc; p++)
  {
   nnz_global += dist_nnz_pt[p];
   dist_first_row[p] = this->first_row(p);
   dist_nrow_local[p] = this->nrow_local(p);
  }

 // conpute the offset for the values and column index data
 int* nnz_offset = new int[MPI_Helpers::Nproc];
 nnz_offset[0] = 0;
 for (int p = 1; p < MPI_Helpers::Nproc; p++)
  {
   nnz_offset[p] = nnz_offset[p-1] + dist_nnz_pt[p-1];
  }
   
 // get pointers to the (current) distributed data
 int* dist_row_start = this->row_start();
 int* dist_column_index = this->column_index();
 double* dist_value = this->value();
   
 // space for the global matrix
 int* global_row_start = new int[nrow+1];
 int* global_column_index = new int[nnz_global];
 double* global_value = new double[nnz_global];

 // get the row starts
 MPI_Allgatherv(dist_row_start,nrow_local,MPI_INT,
                global_row_start,dist_nrow_local,dist_first_row,MPI_INT,
                Distribution_pt->communicator_pt()->mpi_comm());
   
 // get the column indexes
 MPI_Allgatherv(dist_column_index,nnz,MPI_INT,
                global_column_index,dist_nnz_pt,nnz_offset,MPI_INT,
                Distribution_pt->communicator_pt()->mpi_comm());
 
 // get the values
 MPI_Allgatherv(dist_value,nnz,MPI_DOUBLE,
                global_value,dist_nnz_pt,nnz_offset,MPI_DOUBLE,
                Distribution_pt->communicator_pt()->mpi_comm());
   
 // finally the last row start
 global_row_start[nrow] = nnz_global;
   
 // update the other row start
 for (int p = 0; p < nproc; p++)
  {
   for (int i = 0; i < dist_nrow_local[p]; i++)
    {
     unsigned j = dist_first_row[p] + i;
     global_row_start[j]+=nnz_offset[p];
    }
  }

 // create the global distribution
 LinearAlgebraDistribution* dist_pt = new 
  LinearAlgebraDistribution(Distribution_pt->communicator_pt(),nrow,false);
 
 // create the matrix
 CRDoubleMatrix* matrix_pt = new CRDoubleMatrix(dist_pt);
 
 // copy of distribution taken so delete
 delete dist_pt;

 // pass data into matrix
 matrix_pt->rebuild_matrix_without_copy(this->ncol(),nnz_global,global_value,
                                        global_column_index,global_row_start);

 // clean up
 delete dist_first_row;
 delete dist_nrow_local;
 delete nnz_offset;
 delete dist_nnz_pt;
 
 // and return
 return matrix_pt;
#else
 return new CRDoubleMatrix(*this);
#endif
}

/*
#ifdef OOMPH_HAS_MPI

///////////////////////////////////////////////////////////////////////////////
// Functions for DistributedCRDoubleMatrix ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


 //============================================================================
 /// Overload round brackets to give access as a(i,j) (const version)
 //============================================================================
 const double& DistributedCRDoubleMatrix::operator()
  (const unsigned long &i, const unsigned long &j) const
 {
#ifdef PARANOID
  if (i >= static_cast<unsigned>(Matrix_distribution.nrow_local()))
   {
    std::ostringstream error_message;
    error_message << "Requested a coefficient in row i = " << i 
                  << " on processor " << MPI_Helpers::My_rank << ".\n"
                  << "This processor contains " 
                  << Matrix_distribution.nrow_local() << " rows.\n";
    throw OomphLibError(error_message.str(),
                        "DistributedCRDoubleMatrix::() const",
                        OOMPH_EXCEPTION_LOCATION);
    }
   else if (j>=ncol())
    {
     std::ostringstream error_message;
     error_message << "Range Error: j=" << j << " is not in the range (0," 
                   << ncol()-1 << ")." << std::endl; 
     throw OomphLibError(error_message.str(),
                         "DistributedCRDoubleMatrix::() const",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   for (long k=Row_start[i];k<Row_start[i+1];k++)
    {
     if (unsigned(Column_index[k])==j)
      {
       return this->Value[k];
      }
    }
   return this->Zero;
  }

//=============================================================================
/// Multiply the matrix by the vector x: soln=Ax.
/// Currently this requires trilinos to be available. Will break
/// otherwise.
/// The vectors x and soln are both distributed vectors
//=============================================================================
void DistributedCRDoubleMatrix::multiply(const DistributedVector<double> &x, 
                                       DistributedVector<double> &soln)
{
#ifdef HAVE_TRILINOS
 // This will only work if we have trilinos on board
 TrilinosHelpers::multiply(*this,x,soln);
#else
   throw OomphLibError("Matrix-vector product requires trilinos",
                       "DistributedCRDoubleMatrix::multiply()",
                       OOMPH_EXCEPTION_LOCATION);
#endif
}


//=============================================================================
///  Function to multiply this matrix by the 
///  DistributedCRDoubleMatrix matrix_in. This requires Trilinos.
//=============================================================================
void DistributedCRDoubleMatrix::multiply(
 DistributedCRDoubleMatrix& matrix_in,
 DistributedCRDoubleMatrix& result)
{
#ifdef HAVE_TRILINOS
 //This will only work if we have trilinos on board
 TrilinosHelpers::multiply(*this,matrix_in,result,
                           Matrix_matrix_multiply_using_ml);
#else
 throw OomphLibError("Matrix-matrix product requires trilinos",
                     "DistributedCRDoubleMatrix::multiply()",
                     OOMPH_EXCEPTION_LOCATION);
#endif
}


//=============================================================================
/// Call the corresponding 
/// constructor of the CRDoubleMatrix, initialise the
/// member data to zero and set the default linear
/// solver to be SuperLU dist
//============================================================================
DistributedCRDoubleMatrix::DistributedCRDoubleMatrix() 
 : CRDoubleMatrix(),
   Distribution_pt(0)
{ 
 //Delete the exisiting default linear solver pt
 delete Default_linear_solver_pt;

 //Replace with SuperLU dist
 Linear_solver_pt = Default_linear_solver_pt = new SuperLU_dist;

 // Set the trilinos matrix matrix multiplication method to not use the
 // ml method as default
 Matrix_matrix_multiply_using_ml=false;

}

//=============================================================================
/// \short Constructor: Pass vector of values, vector of column indices,
/// and vector of row starts. Since the rows of the matrix are 
/// distributed over several processors we also need to pass the 
/// number of the first row and the number of local rows stored here. 
/// n and m are the total number of rows and columns. 
/// \b NOTE: Number of nonzero entries is read
/// off from value, so make sure the vector has been shrunk
/// to its correct length.
//=============================================================================
DistributedCRDoubleMatrix::
DistributedCRDoubleMatrix(const Vector<double>& value, 
                          const Vector<int>& column_index,
                          const Vector<int>& row_start, 
                          const DistributionInfo& distribution,
                          const unsigned& n_col)
 : CRDoubleMatrix(value,column_index,row_start,
                  distribution.nrow_local(),n_col)
{
 // store the distribution
 Distribution_pt = new DistributionInfo(distribution);

 //Delete the exisiting default linear solver pt
 delete Default_linear_solver_pt;

 //Replace with SuperLU dist
 Linear_solver_pt = Default_linear_solver_pt = new SuperLU_dist;

 // Set the trilinos matrix matrix multiplication method to not use the
 // ml method as default
 Matrix_matrix_multiply_using_ml=false;
}


 //============================================================================
 /// \short Function to delete memory associated with the matrix and
 /// set all values to 0
 //============================================================================
 void DistributedCRDoubleMatrix::clean_up_memory()
 {
  CRDoubleMatrix::clean_up_memory();
  delete Distribution_pt;
  Distribution_pt = 0;
 }


 //===========================================================================
 /// \short Build matrix from compressed representation:
 /// Pass vector of values, vector of column indices,
 /// and vector of row starts. Since the rows of the matrix are 
 /// distributed over several processors we also need to pass the 
 /// number of the first row and the number of local rows stored here. 
 /// n and m are the total number of rows and columns. 
 /// \b NOTE: Number of nonzero entries is read
 /// off from value, so make sure the vector has been shrunk
 /// to its correct length.
 //===========================================================================
void DistributedCRDoubleMatrix::build(const Vector<double>& value,
                                      const Vector<int>& column_index,
                                      const Vector<int>& row_start,
                                      const DistributionInfo& distribution,
                                      const unsigned& n_col)
{
 // call build in parent class
 CRMatrix<double>::build(value,column_index, row_start,
                         distribution.nrow_global(),n_col);
 
 // set the distribution
 delete Distribution_pt;
 Distribution_pt = new DistributionInfo(distribution);
}

#endif
*/
}

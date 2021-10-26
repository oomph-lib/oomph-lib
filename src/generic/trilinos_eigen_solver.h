// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_TRILINOS_EIGEN_SOLVER_HEADER
#define OOMPH_TRILINOS_EIGEN_SOLVER_HEADER

// oomph-lib includes
#include "double_multi_vector.h"
#include "problem.h"
#include "linear_solver.h"
#include "eigen_solver.h"

// Trilinos includes
#include "AnasaziOutputManager.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

namespace Anasazi
{
  //========================================================================
  /// Specialize the Anasazi traits class for the oomph-lib DoubleMultiVector.
  /// This provides the interfaces required by the Anasazi eigensolvers.
  //========================================================================
  template<>
  class MultiVecTraits<double, oomph::DoubleMultiVector>
  {
  public:
    /// \brief Creates a new empty DoubleMultiVector containing numvecs columns.
    /// Return a reference-counted pointer to the new multivector
    static Teuchos::RCP<oomph::DoubleMultiVector> Clone(
      const oomph::DoubleMultiVector& mv, const int numvecs)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(numvecs, mv));
    }

    /// \brief Creates a deep copy of the DoubleMultiVector mv
    /// return Reference-counted pointer to the new DoubleMultiVector
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneCopy(
      const oomph::DoubleMultiVector& mv)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv));
    }

    /// \brief Creates a new oomph::DoubleMultiVector and (deep)
    /// copies the contents of
    /// the vectors in index into the new vector
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneCopy(
      const oomph::DoubleMultiVector& mv, const std::vector<int>& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index));
    }

    /// \brief Deep copy of specified columns of oomph::DoubleMultiVector
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneCopy(
      const oomph::DoubleMultiVector& mv, const Teuchos::Range1D& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index));
    }

    /// \brief Creates a new oomph::DoubleMultiVector that contains
    /// shallow copies
    /// of selected entries of the oomph::DoubleMultiVector mv
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneViewNonConst(
      oomph::DoubleMultiVector& mv, const std::vector<int>& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index, false));
    }

    /// \brief Creates a new oomph::DoubleMultiVector that
    /// contains shallow copies
    /// of selected entries of the oomph::DoubleMultiVector mv
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneViewNonConst(
      oomph::DoubleMultiVector& mv, const Teuchos::Range1D& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index, false));
    }

    /// \brief Creates a new oomph::DoubleMultiVector that contains
    /// shallow copies
    /// of selected entries of the oomph::DoubleMultiVector mv
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    /// (const version)
    static Teuchos::RCP<const oomph::DoubleMultiVector> CloneView(
      const oomph::DoubleMultiVector& mv, const std::vector<int>& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index, false));
    }

    /// \brief Creates a new oomph::DoubleMultiVector that contains
    /// shallow copies
    /// of selected entries of the oomph::DoubleMultiVector mv
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    /// (Non-const version for Trilinos 9 interface)
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneView(
      oomph::DoubleMultiVector& mv, const std::vector<int>& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index, false));
    }

    /// \brief Creates a new oomph::DoubleMultiVector that
    /// contains shallow copies
    /// of selected entries of the oomph::DoubleMultiVector mv
    /// return Reference-counted pointer to the new oomph::DoubleMultiVector
    /// (const version)
    static Teuchos::RCP<oomph::DoubleMultiVector> CloneView(
      oomph::DoubleMultiVector& mv, const Teuchos::Range1D& index)
    {
      return Teuchos::rcp(new oomph::DoubleMultiVector(mv, index, false));
    }

    /// Obtain the global length of the vector
    static int GetVecLength(const oomph::DoubleMultiVector& mv)
    {
      return static_cast<int>(mv.nrow());
    }

    /// Obtain the number of vectors in the multivector
    static int GetNumberVecs(const oomph::DoubleMultiVector& mv)
    {
      return static_cast<int>(mv.nvector());
    }

    /// \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
    static void MvTimesMatAddMv(
      const double alpha,
      const oomph::DoubleMultiVector& A,
      const Teuchos::SerialDenseMatrix<int, double>& B,
      const double beta,
      oomph::DoubleMultiVector& mv)
    {
      // For safety let's (deep) copy A
      oomph::DoubleMultiVector C(A);
      // First pre-multiply by the scalars
      mv *= beta;
      C *= alpha;

      // Find the number of vectors in each multivector
      int mv_n_vector = mv.nvector();
      int A_n_vector = A.nvector();
      // Find the number of rows
      int n_row_local = A.nrow_local();
      // Now need to multiply by the matrix and add the result
      // to the appropriate entry in the multivector
      for (int i = 0; i < n_row_local; ++i)
      {
        for (int j = 0; j < mv_n_vector; j++)
        {
          double res = 0.0;
          for (int k = 0; k < A_n_vector; k++)
          {
            res += C(k, i) * B(k, j);
          }
          mv(j, i) += res;
        }
      }
    }

    /// \brief Replace \c mv with \f$\alpha A + \beta B\f$.
    static void MvAddMv(const double alpha,
                        const oomph::DoubleMultiVector& A,
                        const double beta,
                        const oomph::DoubleMultiVector& B,
                        oomph::DoubleMultiVector& mv)
    {
      const unsigned n_vector = A.nvector();
      const unsigned n_row_local = A.nrow_local();
      for (unsigned v = 0; v < n_vector; v++)
      {
        for (unsigned i = 0; i < n_row_local; i++)
        {
          mv(v, i) = alpha * A(v, i) + beta * B(v, i);
        }
      }
    }

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale(oomph::DoubleMultiVector& mv, const double alpha)
    {
      mv *= alpha;
    }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c
     * alpha[i].
     */
    static void MvScale(oomph::DoubleMultiVector& mv,
                        const std::vector<double>& alpha)
    {
      const unsigned n_vector = mv.nvector();
      const unsigned n_row_local = mv.nrow_local();
      for (unsigned v = 0; v < n_vector; v++)
      {
        for (unsigned i = 0; i < n_row_local; i++)
        {
          mv(v, i) *= alpha[v];
        }
      }
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply
     * \f$ \alpha A^Hmv \f$.
     */
    static void MvTransMv(const double alpha,
                          const oomph::DoubleMultiVector& A,
                          const oomph::DoubleMultiVector& mv,
                          Teuchos::SerialDenseMatrix<int, double>& B)
    {
      const unsigned A_nvec = A.nvector();
      const unsigned A_nrow_local = A.nrow_local();
      const unsigned mv_nvec = mv.nvector();
      // const unsigned mv_nrow_local = mv.nrow_local();

      for (unsigned v1 = 0; v1 < A_nvec; v1++)
      {
        for (unsigned v2 = 0; v2 < mv_nvec; v2++)
        {
          B(v1, v2) = 0.0;
          for (unsigned i = 0; i < A_nrow_local; i++)
          {
            B(v1, v2) += A(v1, i) * mv(v2, i);
          }
          B(v1, v2) *= alpha;
        }
      }

      // This will need to be communicated if we're in parallel and distributed
#ifdef OOMPH_HAS_MPI
      if (A.distributed() &&
          A.distribution_pt()->communicator_pt()->nproc() > 1)
      {
        const int n_total_val = A_nvec * mv_nvec;
        // Pack entries into a vector
        double b[n_total_val];
        double b2[n_total_val];
        unsigned count = 0;
        for (unsigned v1 = 0; v1 < A_nvec; v1++)
        {
          for (unsigned v2 = 0; v2 < mv_nvec; v2++)
          {
            b[count] = B(v1, v2);
            b2[count] = b[count];
            ++count;
          }
        }


        MPI_Allreduce(&b,
                      &b2,
                      n_total_val,
                      MPI_DOUBLE,
                      MPI_SUM,
                      A.distribution_pt()->communicator_pt()->mpi_comm());

        count = 0;
        for (unsigned v1 = 0; v1 < A_nvec; v1++)
        {
          for (unsigned v2 = 0; v2 < mv_nvec; v2++)
          {
            B(v1, v2) = b2[count];
            ++count;
          }
        }
      }
#endif
    }


    /*! \brief Compute a vector \c b where the components are the individual
     * dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] =
     * A[i]^Hmv[i]\f$.
     */
    static void MvDot(const oomph::DoubleMultiVector& mv,
                      const oomph::DoubleMultiVector& A,
                      std::vector<double>& b)
    {
      mv.dot(A, b);
    }


    /*! \brief Compute the 2-norm of each individual vector of \c mv.
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c
      i-th column of \c mv.
    */
    static void MvNorm(const oomph::DoubleMultiVector& mv,
                       std::vector<double>& normvec)
    {
      mv.norm(normvec);
    }

    //@}

    //! @name Initialization methods
    //@{
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated
    by the indices given in \c index.

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv
    indicated by the indices given in \c index, i.e.<tt> mv[index[i]] =
    A[i]</tt>.
    */
    static void SetBlock(const oomph::DoubleMultiVector& A,
                         const std::vector<int>& index,
                         oomph::DoubleMultiVector& mv)
    {
      // Check some stuff
      const unsigned n_index = index.size();
      if (n_index == 0)
      {
        return;
      }
      const unsigned n_row_local = mv.nrow_local();
      for (unsigned v = 0; v < n_index; v++)
      {
        for (unsigned i = 0; i < n_row_local; i++)
        {
          mv(index[v], i) = A(v, i);
        }
      }
    }

    /// \brief Deep copy of A into specified columns of mv
    ///
    /// (Deeply) copy the first <tt>index.size()</tt> columns of \c A
    /// into the columns of \c mv specified by the given index range.
    ///
    /// Postcondition: <tt>mv[i] = A[i - index.lbound()]</tt>
    /// for all <tt>i</tt> in <tt>[index.lbound(), index.ubound()]</tt>
    ///
    /// \param A [in] Source multivector
    /// \param index [in] Inclusive index range of columns of mv;
    ///   index set of the target
    /// \param mv [out] Target multivector
    static void SetBlock(const oomph::DoubleMultiVector& A,
                         const Teuchos::Range1D& index,
                         oomph::DoubleMultiVector& mv)
    {
      // Check some stuff
      const unsigned n_index = index.size();
      if (n_index == 0)
      {
        return;
      }
      const unsigned n_row_local = mv.nrow_local();
      unsigned range_index = index.lbound();
      for (unsigned v = 0; v < n_index; v++)
      {
        for (unsigned i = 0; i < n_row_local; i++)
        {
          mv(range_index, i) = A(v, i);
        }
        ++range_index;
      }
    }

    /// \brief mv := A
    ///
    /// Assign (deep copy) A into mv.
    static void Assign(const oomph::DoubleMultiVector& A,
                       oomph::DoubleMultiVector& mv)
    {
      mv = A;
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom(oomph::DoubleMultiVector& mv)
    {
      const unsigned n_vector = mv.nvector();
      const unsigned n_row_local = mv.nrow_local();
      for (unsigned v = 0; v < n_vector; v++)
      {
        for (unsigned i = 0; i < n_row_local; i++)
        {
          mv(v, i) = Teuchos::ScalarTraits<double>::random();
        }
      }
    }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit(
      oomph::DoubleMultiVector& mv,
      const double alpha = Teuchos::ScalarTraits<double>::zero())
    {
      mv.initialise(alpha);
    }

    //@}

    //! @name Print method
    //@{

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint(const oomph::DoubleMultiVector& mv, std::ostream& os)
    {
      mv.output(os);
    }

    //@}
  };


} // namespace Anasazi

namespace oomph
{
  //===================================================================
  /// Base class for Oomph-lib's Vector Operator classes that will be
  /// used with the DoubleMultiVector
  //==================================================================
  class DoubleMultiVectorOperator
  {
  public:
    /// Empty constructor
    DoubleMultiVectorOperator() {}

    /// virtual destructor
    virtual ~DoubleMultiVectorOperator() {}

    /// The apply interface
    virtual void apply(const DoubleMultiVector& x,
                       DoubleMultiVector& y) const = 0;
  };

} // namespace oomph

namespace Anasazi
{
  //======================================================================
  /// Anasazi Traits class that specialises the apply operator based on
  /// oomph-lib's DoubleVector class and double as the primitive type
  //======================================================================
  template<>
  class OperatorTraits<double,
                       oomph::DoubleMultiVector,
                       oomph::DoubleMultiVectorOperator>
  {
  public:
    /// Matrix operator application method
    static void Apply(const oomph::DoubleMultiVectorOperator& Op,
                      const oomph::DoubleMultiVector& x,
                      oomph::DoubleMultiVector& y)
    {
      Op.apply(x, y);
    }

  }; // End of the specialised traits class

} // namespace Anasazi


namespace oomph
{
  //====================================================================
  /// Class for the shift invert operation
  //====================================================================
  class ProblemBasedShiftInvertOperator : public DoubleMultiVectorOperator
  {
  private:
    // Pointer to the problem
    Problem* Problem_pt;

    // Pointer to the linear solver used
    LinearSolver* Linear_solver_pt;

    // Storage for the shift
    double Sigma;

    // Storage for the matrices
    CRDoubleMatrix *M_pt, *AsigmaM_pt;

  public:
    ProblemBasedShiftInvertOperator(Problem* const& problem_pt,
                                    LinearSolver* const& linear_solver_pt,
                                    const double& sigma = 0.0)
      : Problem_pt(problem_pt), Linear_solver_pt(linear_solver_pt), Sigma(sigma)
    {
      // Before we get into the Arnoldi loop solve the shifted matrix problem
      // Allocated Row compressed matrices for the mass matrix and shifted main
      // matrix
      M_pt = new CRDoubleMatrix(problem_pt->dof_distribution_pt());
      AsigmaM_pt = new CRDoubleMatrix(problem_pt->dof_distribution_pt());

      // Assemble the matrices
      problem_pt->get_eigenproblem_matrices(*M_pt, *AsigmaM_pt, Sigma);

      // Do not report the time taken
      Linear_solver_pt->disable_doc_time();
    }


    // Now specify how to apply the operator
    void apply(const DoubleMultiVector& x, DoubleMultiVector& y) const
    {
      const unsigned n_vec = x.nvector();
      const unsigned n_row_local = x.nrow_local();
      if (n_vec > 1)
      {
        Linear_solver_pt->enable_resolve();
      }
      // Solve the first vector

      DoubleVector X(x.distribution_pt());
      // Premultiply by mass matrix
      M_pt->multiply(x.doublevector(0), X);
      // For some reason I need to create a new vector and copy here
      // This is odd and not expected, examine carefully
      DoubleVector Y(x.distribution_pt());
      Linear_solver_pt->solve(AsigmaM_pt, X, Y);
      // Need to synchronise
      //#ifdef OOMPH_HAS_MPI
      //   Problem_pt->synchronise_all_dofs();
      //#endif
      for (unsigned i = 0; i < n_row_local; i++)
      {
        y(0, i) = Y[i];
      }

      // Now loop over the others and resolve
      for (unsigned v = 1; v < n_vec; ++v)
      {
        M_pt->multiply(x.doublevector(v), X);
        Linear_solver_pt->resolve(X, Y);
        //#ifdef OOMPH_HAS_MPI
        //     Problem_pt->synchronise_all_dofs();
        //#endif
        for (unsigned i = 0; i < n_row_local; i++)
        {
          y(v, i) = Y[i];
        }
      }
    }
  };


  //=====================================================================
  /// Class for the Anasazi eigensolver
  //=====================================================================
  class ANASAZI : public EigenSolver
  {
  private:
    typedef double ST;
    typedef Teuchos::ScalarTraits<ST> SCT;
    typedef SCT::magnitudeType MT;
    typedef Anasazi::MultiVec<ST> MV;
    typedef Anasazi::Operator<ST> OP;
    typedef Anasazi::MultiVecTraits<ST, MV> MVT;
    typedef Anasazi::OperatorTraits<ST, MV, OP> OPT;

    /// Pointer to output manager
    Anasazi::OutputManager<ST>* Output_manager_pt;

    /// Pointer to a linear solver
    LinearSolver* Linear_solver_pt;

    /// Pointer to a default linear solver
    LinearSolver* Default_linear_solver_pt;

    /// Integer to set whether the real, imaginary or magnitude is
    /// required
    /// to be small or large.
    int Spectrum;

    /// Number of Arnoldi vectors to compute
    int NArnoldi;

    /// Set the shifted value
    double Sigma;

    /// Boolean to set which part of the spectrum left (default) or right
    /// of the shifted value.
    bool Small;

    /// Boolean to indicate whether or not to compute the eigenvectors
    bool Compute_eigenvectors;


  public:
    /// Constructor
    ANASAZI()
      : Linear_solver_pt(0),
        Default_linear_solver_pt(0),
        Spectrum(0),
        NArnoldi(10),
        Sigma(0.0),
        Compute_eigenvectors(true)
    {
      Output_manager_pt = new Anasazi::BasicOutputManager<ST>();
      // Set verbosity level
      int verbosity =
        Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
      Output_manager_pt->setVerbosity(verbosity);

      // print greeting
      Output_manager_pt->stream(Anasazi::Warnings)
        << Anasazi::Anasazi_Version() << std::endl
        << std::endl;
    }

    /// Empty copy constructor
    ANASAZI(const ANASAZI&) : Sigma(0.0) {}

    /// Destructor, delete the linear solver
    virtual ~ANASAZI() {}

    /// Access function for the number of Arnoldi vectors
    int& narnoldi()
    {
      return NArnoldi;
    }

    /// Access function for the number of Arnoldi vectors (const version)
    const int& narnoldi() const
    {
      return NArnoldi;
    }

    /// Set to enable the computation of the eigenvectors (default)
    void enable_compute_eigenvectors()
    {
      Compute_eigenvectors = true;
    }

    /// Set to disable the computation of the eigenvectors
    void disable_compute_eigenvectors()
    {
      Compute_eigenvectors = false;
    }

    /// Solve the eigen problem
    void solve_eigenproblem(Problem* const& problem_pt,
                            const int& n_eval,
                            Vector<std::complex<double>>& eigenvalue,
                            Vector<DoubleVector>& eigenvector)
    {
      // No access to sigma, so set from sigma real
      Sigma = Sigma_real;

      // Initially be dumb here
      Linear_solver_pt = problem_pt->linear_solver_pt();

      // Let's make the initial vector
      Teuchos::RCP<DoubleMultiVector> initial = Teuchos::rcp(
        new DoubleMultiVector(1, problem_pt->dof_distribution_pt()));
      Anasazi::MultiVecTraits<double, DoubleMultiVector>::MvRandom(*initial);

      // Make the operator
      Teuchos::RCP<DoubleMultiVectorOperator> Op_pt =
        Teuchos::rcp(new ProblemBasedShiftInvertOperator(
          problem_pt, this->linear_solver_pt(), Sigma));

      // Create the basic problem
      Teuchos::RCP<Anasazi::BasicEigenproblem<double,
                                              DoubleMultiVector,
                                              DoubleMultiVectorOperator>>
        anasazi_pt = Teuchos::rcp(
          new Anasazi::BasicEigenproblem<double,
                                         DoubleMultiVector,
                                         DoubleMultiVectorOperator>(Op_pt,
                                                                    initial));

      // Think I have it?
      anasazi_pt->setHermitian(false);

      // set the number of eigenvalues requested
      anasazi_pt->setNEV(n_eval);

      // Inform the eigenproblem that you are done passing it information
      if (anasazi_pt->setProblem() != true)
      {
        oomph_info << "Anasazi::BasicEigenproblem::setProblem() had an error."
                   << std::endl
                   << "End Result: TEST FAILED" << std::endl;
      }

      // Create the solver manager
      // No need to have ncv specificed, Triliinos has a sensible default
      //  int ncv = 10;
      MT tol = 1.0e-10;
      int verbosity =
        Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;

      Teuchos::ParameterList MyPL;
      MyPL.set("Which", "LM");
      MyPL.set("Block Size", 1);
      //  MyPL.set( "Num Blocks", ncv );
      MyPL.set("Maximum Restarts", 500);
      MyPL.set("Orthogonalization", "DGKS");
      MyPL.set("Verbosity", verbosity);
      MyPL.set("Convergence Tolerance", tol);
      Anasazi::BlockKrylovSchurSolMgr<double,
                                      DoubleMultiVector,
                                      DoubleMultiVectorOperator>
        BKS(anasazi_pt, MyPL);

      // Solve the problem to the specified tolerances or length
      Anasazi::ReturnType ret = BKS.solve();
      bool testFailed = false;
      if (ret != Anasazi::Converged)
      {
        testFailed = true;
      }

      if (testFailed)
      {
        oomph_info << "Eigensolver not converged\n";
      }

      // Get the eigenvalues and eigenvectors from the eigenproblem
      Anasazi::Eigensolution<double, DoubleMultiVector> sol =
        anasazi_pt->getSolution();
      std::vector<Anasazi::Value<double>> evals = sol.Evals;
      Teuchos::RCP<DoubleMultiVector> evecs = sol.Evecs;

      eigenvalue.resize(evals.size());
      eigenvector.resize(evals.size());

      for (unsigned i = 0; i < evals.size(); i++)
      {
        // Undo shift and invert
        double a = evals[i].realpart;
        double b = evals[i].imagpart;
        double det = a * a + b * b;
        eigenvalue[i] = std::complex<double>(a / det + Sigma, -b / det);

        // Now set the eigenvectors, I hope
        eigenvector[i].build(evecs->distribution_pt());
        unsigned nrow_local = evecs->nrow_local();
        // Would be faster with pointers, but I'll sort that out later!
        for (unsigned n = 0; n < nrow_local; n++)
        {
          eigenvector[i][n] = (*evecs)(i, n);
        }
      }
    }

    /// Set the desired eigenvalues to be left of the shift
    void get_eigenvalues_left_of_shift()
    {
      Small = true;
    }

    /// Set the desired eigenvalues to be right of the shift
    void get_eigenvalues_right_of_shift()
    {
      Small = false;
    }

    /// Set the real part to be the quantity of interest (default)
    void track_eigenvalue_real_part()
    {
      Spectrum = 1;
    }

    /// Set the imaginary part fo the quantity of interest
    void track_eigenvalue_imaginary_part()
    {
      Spectrum = -1;
    }

    /// Set the magnitude to be the quantity of interest
    void track_eigenvalue_magnitude()
    {
      Spectrum = 0;
    }

    /// Return a pointer to the linear solver object
    LinearSolver*& linear_solver_pt()
    {
      return Linear_solver_pt;
    }

    /// Return a pointer to the linear solver object (const version)
    LinearSolver* const& linear_solver_pt() const
    {
      return Linear_solver_pt;
    }
  };

} // namespace oomph


#endif

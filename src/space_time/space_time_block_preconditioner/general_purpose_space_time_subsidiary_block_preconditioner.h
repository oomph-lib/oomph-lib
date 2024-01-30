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
// Header file for SpaceTimeNavierStokes elements
#ifndef OOMPH_GENERAL_PURPOSE_SPACE_TIME_SUBSIDIARY_BLOCK_PRECONDITIONER_HEADER
#define OOMPH_GENERAL_PURPOSE_SPACE_TIME_SUBSIDIARY_BLOCK_PRECONDITIONER_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "generic/iterative_linear_solver.h"
#include "generic/block_preconditioner.h"
#include "generic/matrices.h"

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //=============================================================================
  /// General purpose block triangular preconditioner. By default this is
  /// Upper triangular. Also, by default ExactPreconditioner
  /// is used to solve the subsidiary systems, but
  /// other preconditioners can be used by setting them using passing a pointer
  /// to a function of type SubsidiaryPreconditionerFctPt to the method
  /// subsidiary_preconditioner_function_pt().
  //=============================================================================
  class SpaceTimeNavierStokesSubsidiaryPreconditioner
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor. (By default this preconditioner is upper triangular).
    SpaceTimeNavierStokesSubsidiaryPreconditioner()
      : BlockPreconditioner<CRDoubleMatrix>()
    {
      // Flag to indicate that the preconditioner has been setup
      // previously -- if setup() is called again, data can
      // be wiped.
      Preconditioner_has_been_setup = false;

      // By default we use SuperLU for both p and f blocks
      Using_default_p_preconditioner = true;
      Using_default_f_preconditioner = true;

      // Set default preconditioners (inexact solvers) -- they are
      // members of this class!
      P_preconditioner_pt = 0;
      F_preconditioner_pt = 0;

      // Null the momentum block matrix-vector product helper
      F_mat_vec_pt = 0;

      // Null the gradient block matrix-vector product helper
      G_mat_vec_pt = 0;
    }

    /// Destructor - delete the preconditioner matrices
    virtual ~SpaceTimeNavierStokesSubsidiaryPreconditioner()
    {
      // Call the auxiliary clean up function
      this->clean_up_memory();
    } // End of ~SpaceTimeNavierStokesSubsidiaryPreconditioner

    /// Clean up the memory
    virtual void clean_up_memory()
    {
      // If we've actually set the preconditioner up
      if (Preconditioner_has_been_setup)
      {
        // Delete matvecs
        delete F_mat_vec_pt;
        F_mat_vec_pt = 0;

        delete G_mat_vec_pt;
        G_mat_vec_pt = 0;

        // Delete the momentum block solver
        delete F_preconditioner_pt;
        F_preconditioner_pt = 0;

        // Delete the Schur complement solver
        delete P_preconditioner_pt;
        P_preconditioner_pt = 0;
      } // if (Preconditioner_has_been_setup)
    } // End of clean_up_memory

    /// Broken copy constructor
    SpaceTimeNavierStokesSubsidiaryPreconditioner(
      const SpaceTimeNavierStokesSubsidiaryPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const SpaceTimeNavierStokesSubsidiaryPreconditioner&) =
      delete;

    /// For some reason we need to remind the compiler that there is
    /// also a function named setup in the base class.
    using Preconditioner::setup;

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

  private:
    /// Pointer to the 'preconditioner' for the F matrix
    Preconditioner* F_preconditioner_pt;

    /// Pointer to the 'preconditioner' for the pressure matrix
    Preconditioner* P_preconditioner_pt;

    /// Flag indicating whether the default F preconditioner is used
    bool Using_default_f_preconditioner;

    /// Flag indicating whether the default P preconditioner is used
    bool Using_default_p_preconditioner;

    /// Control flag is true if the preconditioner has been setup
    /// (used so we can wipe the data when the preconditioner is called again)
    bool Preconditioner_has_been_setup;

    /// MatrixVectorProduct operator for F
    MatrixVectorProduct* F_mat_vec_pt;

    /// MatrixVectorProduct operator for G
    MatrixVectorProduct* G_mat_vec_pt;
  };


  //=============================================================================
  /// The block preconditioner form of GMRES. This version extracts
  /// the blocks from the global systems and assembles the system by
  /// concatenating all the matrices together
  //=============================================================================
  class GMRESBlockPreconditioner
    : public IterativeLinearSolver,
      public virtual BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor (empty)
    GMRESBlockPreconditioner()
      : BlockPreconditioner<CRDoubleMatrix>(),
        Matrix_pt(0),
        Navier_stokes_subsidiary_preconditioner_pt(0),
        Iterations(0),
        Preconditioner_has_been_setup(false),
        Preconditioner_LHS(false)
    {
    }

    /// Destructor
    virtual ~GMRESBlockPreconditioner()
    {
      // Call the auxiliary clean up function
      this->clean_up_memory();
    } // End of ~GMRESBlockPreconditioner

    /// Clean up the memory (empty for now...)
    virtual void clean_up_memory()
    {
      // If the block preconditioner has been set up previously
      if (Preconditioner_has_been_setup)
      {
        // Delete the subsidiary preconditioner
        delete Navier_stokes_subsidiary_preconditioner_pt;

        // Make it a null pointer
        Navier_stokes_subsidiary_preconditioner_pt = 0;

        // Delete the matrix!
        delete Matrix_pt;

        // Make it a null pointer
        Matrix_pt = 0;
      }
    } // End of clean_up_memory

    /// Broken copy constructor
    GMRESBlockPreconditioner(const GMRESBlockPreconditioner&) = delete;

    /// Broken assignment operator
    void operator=(const GMRESBlockPreconditioner&) = delete;

    /// For some reason we need to remind the compiler that there is
    /// also a function named setup in the base class.
    using Preconditioner::setup;

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Solver: Takes pointer to problem and returns the results vector
    /// which contains the solution of the linear system defined by
    /// the problem's fully assembled Jacobian and residual vector.
    void solve(Problem* const& problem_pt, DoubleVector& result)
    {
      // Broken
      throw OomphLibError("Can't use this interface!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // End of solve

    /// Handle to the number of iterations taken
    unsigned iterations() const
    {
      // Return the number of iterations
      return Iterations;
    } // End of iterations

    /// Set left preconditioning (the default)
    void set_preconditioner_LHS()
    {
      Preconditioner_LHS = true;
    }

    /// Enable right preconditioning
    void set_preconditioner_RHS()
    {
      Preconditioner_LHS = false;
    }

    /// Handle to the Navier-Stokes subsidiary block preconditioner
    /// DRAIG: Make sure the desired const-ness is correct later...
    SpaceTimeNavierStokesSubsidiaryPreconditioner* navier_stokes_subsidiary_preconditioner_pt()
      const
    {
      // Return a pointer to the appropriate member data
      return Navier_stokes_subsidiary_preconditioner_pt;
    } // End of navier_stokes_subsidiary_preconditioner_pt

  private:
    /// Helper function to update the result vector using the result,
    /// x=x_0+V_m*y
    void update(const unsigned& k,
                const Vector<Vector<double>>& H,
                const Vector<double>& s,
                const Vector<DoubleVector>& v,
                const DoubleVector& block_x_with_size_of_full_x,
                DoubleVector& x)
    {
      // Make a local copy of s
      Vector<double> y(s);

      // Backsolve:
      for (int i = int(k); i >= 0; i--)
      {
        // Divide the i-th entry of y by the i-th diagonal entry of H
        y[i] /= H[i][i];

        // Loop over the previous values of y and update
        for (int j = i - 1; j >= 0; j--)
        {
          // Update the j-th entry of y
          y[j] -= H[i][j] * y[i];
        }
      } // for (int i=int(k);i>=0;i--)

      // Store the number of rows in the result vector
      unsigned n_x = x.nrow();

      // Build a temporary vector with entries initialised to 0.0
      DoubleVector temp(x.distribution_pt(), 0.0);

      // Build a temporary vector with entries initialised to 0.0
      DoubleVector z(x.distribution_pt(), 0.0);

      // Get access to the underlying values
      double* temp_pt = temp.values_pt();

      // Calculate x=Vy
      for (unsigned j = 0; j <= k; j++)
      {
        // Get access to j-th column of v
        const double* vj_pt = v[j].values_pt();

        // Loop over the entries of the vector, temp
        for (unsigned i = 0; i < n_x; i++)
        {
          temp_pt[i] += vj_pt[i] * y[j];
        }
      } // for (unsigned j=0;j<=k;j++)

      // If we're using LHS preconditioning
      if (Preconditioner_LHS)
      {
        // Since we're using LHS preconditioning the preconditioner is applied
        // to the matrix and RHS vector so we simply update the value of x
        x += temp;
      }
      // If we're using RHS preconditioning
      else
      {
        // This is pretty inefficient but there is no other way to do this with
        // block sub preconditioners as far as I can tell because they demand
        // to have the full r and z vectors...
        DoubleVector block_z_with_size_of_full_z(
          block_x_with_size_of_full_x.distribution_pt());
        DoubleVector temp_with_size_of_full_rhs(
          block_x_with_size_of_full_x.distribution_pt());

        // Reorder the entries of temp and put them into the global-sized vector
        this->return_block_vector(0, temp, temp_with_size_of_full_rhs);

        // Solve on the global-sized vectors
        Navier_stokes_subsidiary_preconditioner_pt->preconditioner_solve(
          temp_with_size_of_full_rhs, block_z_with_size_of_full_z);

        // Now reorder the entries and put them into z
        this->get_block_vector(0, block_z_with_size_of_full_z, z);

        // Since we're using RHS preconditioning the preconditioner is applied
        // to the solution vector
        // preconditioner_pt()->preconditioner_solve(temp,z);

        // Use the update: x_m=x_0+inv(M)Vy [see Saad Y,"Iterative methods for
        // sparse linear systems", p.284]
        x += z;
      }
    } // End of update

    /// Helper function: Generate a plane rotation. This is done by
    /// finding the values of \f$ \cos(\theta) \f$ (i.e. cs) and \sin(\theta)
    /// (i.e. sn) such that:
    /// \f[ \begin{bmatrix} \cos\theta & \sin\theta \newline -\sin\theta & \cos\theta \end{bmatrix} \begin{bmatrix} dx \newline dy \end{bmatrix} = \begin{bmatrix} r \newline 0 \end{bmatrix}, \f]
    /// where \f$ r=\sqrt{pow(dx,2)+pow(dy,2)} \f$. The values of a and b are
    /// given by:
    /// \f[ \cos\theta=\dfrac{dx}{\sqrt{pow(dx,2)+pow(dy,2)}}, \f]
    /// and
    /// \f[ \sin\theta=\dfrac{dy}{\sqrt{pow(dx,2)+pow(dy,2)}}. \f]
    /// Taken from: Saad Y."Iterative methods for sparse linear systems", p.192
    void generate_plane_rotation(double& dx, double& dy, double& cs, double& sn)
    {
      // If dy=0 then we do not need to apply a rotation
      if (dy == 0.0)
      {
        // Using theta=0 gives cos(theta)=1
        cs = 1.0;

        // Using theta=0 gives sin(theta)=0
        sn = 0.0;
      }
      // If dx or dy is large using the normal form of calculting cs and sn
      // is naive since this may overflow or underflow so instead we calculate
      // r=sqrt(pow(dx,2)+pow(dy,2)) by using r=|dy|sqrt(1+pow(dx/dy,2)) if
      // |dy|>|dx| [see <A
      // HREF=https://en.wikipedia.org/wiki/Hypot">Hypot</A>.].
      else if (fabs(dy) > fabs(dx))
      {
        // Since |dy|>|dx| calculate the ratio dx/dy
        double temp = dx / dy;

        // Calculate sin(theta)=dy/sqrt(pow(dx,2)+pow(dy,2))
        sn = 1.0 / sqrt(1.0 + temp * temp);

        // Calculate cos(theta)=dx/sqrt(pow(dx,2)+pow(dy,2))=(dx/dy)*sin(theta)
        cs = temp * sn;
      }
      // Otherwise, we have |dx|>=|dy| so to, again, avoid overflow or underflow
      // calculate the values of cs and sn using the method above
      else
      {
        // Since |dx|>=|dy| calculate the ratio dy/dx
        double temp = dy / dx;

        // Calculate cos(theta)=dx/sqrt(pow(dx,2)+pow(dy,2))
        cs = 1.0 / sqrt(1.0 + temp * temp);

        // Calculate sin(theta)=dy/sqrt(pow(dx,2)+pow(dy,2))=(dy/dx)*cos(theta)
        sn = temp * cs;
      }
    } // End of generate_plane_rotation

    /// Helper function: Apply plane rotation. This is done using the
    /// update:
    /// \f[ \begin{bmatrix} dx \newline dy \end{bmatrix} \leftarrow \begin{bmatrix} \cos\theta & \sin\theta \newline -\sin\theta & \cos\theta \end{bmatrix} \begin{bmatrix} dx \newline dy \end{bmatrix}. \f]
    void apply_plane_rotation(double& dx, double& dy, double& cs, double& sn)
    {
      // Calculate the value of dx but don't update it yet
      double temp = cs * dx + sn * dy;

      // Set the value of dy
      dy = -sn * dx + cs * dy;

      // Set the value of dx using the correct values of dx and dy
      dx = temp;
    } // End of apply_plane_rotation

    /// Pointer to matrix
    CRDoubleMatrix* Matrix_pt;

    /// Pointer to the preconditioner for the block matrix
    SpaceTimeNavierStokesSubsidiaryPreconditioner*
      Navier_stokes_subsidiary_preconditioner_pt;

    /// Number of iterations taken
    unsigned Iterations;

    /// Control flag is true if the preconditioner has been setup (used
    /// so we can wipe the data when the preconditioner is called again)
    bool Preconditioner_has_been_setup;

    /// boolean indicating use of left hand preconditioning (if true)
    /// or right hand preconditioning (if false)
    bool Preconditioner_LHS;
  };
} // End of namespace oomph
#endif

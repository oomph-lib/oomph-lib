// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_FSI_PRECONDITIONERS_HEADER
#define OOMPH_FSI_PRECONDITIONERS_HEADER


#include "../navier_stokes/navier_stokes_preconditioners.h"

namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //============================================================================
  /// FSI preconditioner. This extracts upper/lower triangular
  /// blocks in the 3x3 overall block matrix structure arising from
  /// the monolithic discretisation of FSI problems with algebraic
  /// node updates. Dofs are decomposed into fluid velocity, pressure
  /// and solid unknowns. NavierStokesSchurComplementPreconditioner is used
  /// as the inexact solver for the fluid block; SuperLU (in
  /// its incarnation as an "exact" preconditioner) is used for
  /// the solid block. By default we retain the fluid on solid off
  /// diagonal blocks.
  //=============================================================================
  class FSIPreconditioner : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor: By default use block triangular form with retained
    /// fluid on solid terms. A problem pointer is required for the underlying
    /// NavierStokesSchurComplementPreconditioner.
    FSIPreconditioner(Problem* problem_pt)
    {
      // set the mesh pointers
      this->set_nmesh(2);
      Navier_stokes_mesh_pt = 0;
      Wall_mesh_pt = 0;

      // Initially assume that there are no multiple element types in the
      // meshes.
      Allow_multiple_element_type_in_navier_stokes_mesh = false;
      Allow_multiple_element_type_in_wall_mesh = false;

      // Default setting: Fluid onto solid as it this was shown to be
      // marginally faster than solid onto fluid; see Heil CMAME 193 (2004)
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = true;

      // Create the Navier Stokes Schur complement preconditioner
      Navier_stokes_preconditioner_pt =
        new NavierStokesSchurComplementPreconditioner(problem_pt);

      // Create the Solid preconditioner
      Solid_preconditioner_pt = new SuperLUPreconditioner;

      // Preconditioner hasn't been set up yet.
      Preconditioner_has_been_setup = false;

      // Create the matrix vector product operators
      Matrix_vector_product_0_1_pt = new MatrixVectorProduct;
      Matrix_vector_product_1_0_pt = new MatrixVectorProduct;

      // set Doc_time to false
      Doc_time = false;
    }


    /// Destructor: Clean up.
    ~FSIPreconditioner()
    {
      // Delete the Navier-Stokes preconditioner (inexact solver)
      delete Navier_stokes_preconditioner_pt;

      // Delete the solid preconditioner (inexact solver)
      delete Solid_preconditioner_pt;

      // delete the matrix vector product operators
      delete Matrix_vector_product_0_1_pt;
      delete Matrix_vector_product_1_0_pt;
    }


    /// Broken copy constructor
    FSIPreconditioner(const FSIPreconditioner&) = delete;


    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const FSIPreconditioner&) =
      delete;*/

    /// Set solid preconditioner (deletes existing one)
    void set_solid_preconditioner_pt(Preconditioner* solid_preconditioner_pt)
    {
      // Kill existing one
      if (Solid_preconditioner_pt != 0)
      {
        delete Solid_preconditioner_pt;
      }
      Solid_preconditioner_pt = solid_preconditioner_pt;
    }

    /// Read-only access to solid preconditoner (use set_... to set it)
    Preconditioner* solid_preconditioner_pt() const
    {
      return Solid_preconditioner_pt;
    }


    /// Switch to block-diagonal preconditioner
    void use_block_diagonal_version()
    {
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = false;
    }

    /// Switch to block-triangular preconditioner in which
    /// action of fluid dofs onto solid equations is retained
    void use_block_triangular_version_with_fluid_on_solid()
    {
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = true;
    }

    /// Switch to block-triangular preconditioner in which
    /// action of solid dofs onto fluid equations is retained
    void use_block_triangular_version_with_solid_on_fluid()
    {
      Retain_solid_onto_fluid_terms = true;
      Retain_fluid_onto_solid_terms = false;
    }

    /// Setter function for the mesh containing the
    /// block-preconditionable Navier-Stokes elements. The optional argument
    /// indicates if there are more than one type of elements in same mesh.
    void set_navier_stokes_mesh(
      Mesh* mesh_pt,
      const bool& allow_multiple_element_type_in_navier_stokes_mesh = false)
    {
      // Store the mesh pointer.
      Navier_stokes_mesh_pt = mesh_pt;

      // Are there multiple element types in the Navier-Stokes mesh?
      Allow_multiple_element_type_in_navier_stokes_mesh =
        allow_multiple_element_type_in_navier_stokes_mesh;
    }

    /// Setter function for the mesh containing the
    /// block-preconditionable FSI solid elements. The optional argument
    /// indicates if there are more than one type of elements in the same mesh.
    void set_wall_mesh(
      Mesh* mesh_pt,
      const bool& allow_multiple_element_type_in_wall_mesh = false)
    {
      // Store the mesh pointer
      Wall_mesh_pt = mesh_pt;

      // Are there multiple element types in the wall mesh?
      Allow_multiple_element_type_in_wall_mesh =
        allow_multiple_element_type_in_wall_mesh;
    }

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Access function to the Navier Stokes preconditioner (inexact solver)
    NavierStokesSchurComplementPreconditioner* navier_stokes_preconditioner_pt()
      const
    {
      return Navier_stokes_preconditioner_pt;
    }

    /// Enable documentation of time
    void enable_doc_time()
    {
      Doc_time = true;
    }

    /// Disable documentation of time
    void disable_doc_time()
    {
      Doc_time = false;
    }


  private:
    /// Pointer the Navier Stokes preconditioner (inexact solver)
    NavierStokesSchurComplementPreconditioner* Navier_stokes_preconditioner_pt;

    /// Pointer to the solid preconditioner  (inexact solver)
    Preconditioner* Solid_preconditioner_pt;

    /// Pointer to fluid/solid interaction matrix
    MatrixVectorProduct* Matrix_vector_product_0_1_pt;

    /// Pointer to solid/fluid solid interaction matrix
    MatrixVectorProduct* Matrix_vector_product_1_0_pt;

    /// Boolean indicating the preconditioner has been set up
    bool Preconditioner_has_been_setup;

    /// Boolean flag used to indicate that the solid onto fluid
    /// interaction terms are to be retained
    bool Retain_solid_onto_fluid_terms;

    /// Boolean flag used to indicate that the fluid onto solid
    /// interaction terms are to be retained
    bool Retain_fluid_onto_solid_terms;

    /// Set Doc_time to true for outputting results of timings
    bool Doc_time;

    /// Pointer to the navier stokes mesh
    Mesh* Navier_stokes_mesh_pt;

    /// pointer to the solid mesh
    Mesh* Wall_mesh_pt;

    /// Flag to indicate if there are multiple element types in the
    /// Navier-Stokes mesh.
    bool Allow_multiple_element_type_in_navier_stokes_mesh;

    // Flag to indicate if there are multiple element types in the Wall mesh.
    bool Allow_multiple_element_type_in_wall_mesh;
  };


  /// ///////////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////////
  // FSI preconditioner member functions
  /// ///////////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////////


  //=============================================================================
  /// Setup the preconditioner. Note: Matrix must be a CRDoubleMatrix.
  //=============================================================================
  void FSIPreconditioner::setup()
  {
    // check the meshes have been set
#ifdef PARANOID
    if (Navier_stokes_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Pointer to fluid mesh hasn't been set!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (Wall_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Pointer to solid mesh hasn't been set!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // setup the meshes
    this->set_mesh(0,
                   Navier_stokes_mesh_pt,
                   Allow_multiple_element_type_in_navier_stokes_mesh);
    this->set_mesh(1, Wall_mesh_pt, Allow_multiple_element_type_in_wall_mesh);

    // get the number of fluid dofs from teh first element in the mesh
    unsigned n_fluid_dof = this->ndof_types_in_mesh(0);
    unsigned n_dof = n_fluid_dof + this->ndof_types_in_mesh(1);

    // this fsi preconditioner has two types of DOF fluid dofs and solid dofs
    Vector<unsigned> dof_to_block_map(n_dof, 0);
    for (unsigned i = n_fluid_dof; i < n_dof; i++)
    {
      dof_to_block_map[i] = 1;
    }

    // Call block setup for this preconditioner
    this->block_setup(dof_to_block_map);

    // Block mapping for the subsidiary Navier Stokes preconditioner:
    // blocks 0 and 1 in the FSI preconditioner are also blocks 0 and 1
    // in the subsidiary Navier Stokes one.
    Vector<unsigned> ns_dof_lookup(n_fluid_dof);
    for (unsigned i = 0; i < n_fluid_dof; i++)
    {
      ns_dof_lookup[i] = i;
    }

    // Turn the Navier Stokes Schur complement preconditioner into a
    // subsidiary preconditioner of this preconditioner
    Navier_stokes_preconditioner_pt->turn_into_subsidiary_block_preconditioner(
      this, ns_dof_lookup);

    // Setup the navier stokes preconditioner: Tell it about the
    // Navier Stokes mesh and set it up.
    Navier_stokes_preconditioner_pt->set_navier_stokes_mesh(
      Navier_stokes_mesh_pt);
    Navier_stokes_preconditioner_pt->setup(matrix_pt());

    // Extract the additional blocks we need for FSI:

    // Solid tangent stiffness matrix
    CRDoubleMatrix block_matrix_1_1;
    this->get_block(1, 1, block_matrix_1_1);

    // Setup the solid preconditioner (inexact solver)
    double t_start = TimingHelpers::timer();
    Solid_preconditioner_pt->setup(&block_matrix_1_1);
    double t_end = TimingHelpers::timer();
    block_matrix_1_1.clear();
    double setup_time = t_end - t_start;

    // Solid on fluid terms (if needed)
    if (Retain_solid_onto_fluid_terms)
    {
      CRDoubleMatrix block_matrix_0_1 = get_block(0, 1);
      this->setup_matrix_vector_product(
        Matrix_vector_product_0_1_pt, &block_matrix_0_1, 1);
    }

    // Fluid on solid terms (if needed)
    if (Retain_fluid_onto_solid_terms)
    {
      CRDoubleMatrix block_matrix_1_0 = get_block(1, 0);
      this->setup_matrix_vector_product(
        Matrix_vector_product_1_0_pt, &block_matrix_1_0, 0);
    }

    // Output times
    if (Doc_time)
    {
      oomph_info << "Solid sub-preconditioner setup time [sec]: " << setup_time
                 << "\n";
    }

    // We're done (and we stored some data)
    Preconditioner_has_been_setup = true;
  }


  //======================================================================
  /// Apply preconditioner to Vector r
  //======================================================================
  void FSIPreconditioner::preconditioner_solve(const DoubleVector& r,
                                               DoubleVector& z)
  {
    // if z is not setup then give it the same distribution
    if (!z.built())
    {
      z.build(r.distribution_pt(), 0.0);
    }

    // Make copy of residual vector (to overcome const-ness
    DoubleVector res(r);


    // Retain off-diagonals that represent effect of solid on fluid
    //-------------------------------------------------------------
    if (Retain_solid_onto_fluid_terms)
    {
      // Working vectors
      DoubleVector temp_solid_vec;
      DoubleVector temp_fluid_vec;

      // Copy solid values from residual to temp_vec:
      // Loop over all entries in the global vector (this one
      // includes solid, velocity and pressure dofs in some random fashion)
      get_block_vector(1, res, temp_solid_vec);

      // Solve solid system by back-substitution
      // with LU-decomposed stiffness matrix
      DoubleVector temp_solid_vec2;
      Solid_preconditioner_pt->preconditioner_solve(temp_solid_vec,
                                                    temp_solid_vec2);
      this->return_block_vector(1, temp_solid_vec2, z);

      // NOTE: temp_solid_vec now contains z_s = S^{-1} r_s

      // Multiply C_{us} by z_s
      Matrix_vector_product_0_1_pt->multiply(temp_solid_vec2, temp_fluid_vec);
      temp_solid_vec.clear();

      // Subtract from fluid residual vector for fluid solve
      DoubleVector another_temp_vec;
      this->get_block_vector(0, res, another_temp_vec);
      another_temp_vec -= temp_fluid_vec;
      this->return_block_vector(0, another_temp_vec, res);

      // now apply the navier stokes lsc preconditioner
      Navier_stokes_preconditioner_pt->preconditioner_solve(res, z);
    }


    // Retain off-diagonals that represent effect of fluid on solid
    //-------------------------------------------------------------
    // (or diagonal preconditioner)
    //-----------------------------
    else
    {
      // Call fluid preconditioner for fluid block
      Navier_stokes_preconditioner_pt->preconditioner_solve(res, z);

      // Working vectors
      DoubleVector temp_solid_vec;

      // get the solid vector
      get_block_vector(1, res, temp_solid_vec);

      // Do matrix vector products with fluid onto solid coupling matrices:
      if (Retain_fluid_onto_solid_terms)
      {
        DoubleVector temp_fluid_vec;
        get_block_vector(0, z, temp_fluid_vec);

        // Auxiliary vector to hold the matrix vector product of the
        // fluid-onto-solid coupling matrices with the fluid solutions:
        DoubleVector aux_vec;

        // Multiply C_{su} by z_u
        Matrix_vector_product_1_0_pt->multiply(temp_fluid_vec, aux_vec);

        // ...and subtract from r_s:
        temp_solid_vec -= aux_vec;
      }

      // Solve solid system by back-substitution
      // with LU-decomposed stiffness matrix
      DoubleVector temp_solid_vec2;
      Solid_preconditioner_pt->preconditioner_solve(temp_solid_vec,
                                                    temp_solid_vec2);

      // Now copy result_vec (i.e. z_s) back into the global vector z.
      // Loop over all entries in the global results vector z:
      return_block_vector(1, temp_solid_vec2, z);
    }
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //============================================================================
  /// FSI preconditioner. This extracts upper/lower triangular
  /// blocks in the 3x3 overall block matrix structure arising from
  /// the monolithic discretisation of FSI problems with algebraic
  /// node updates. Dofs are decomposed into fluid velocity, pressure
  /// and solid unknowns. Blocks are then re-assembled into one global
  /// matrix and solved with a direct solver (SuperLU in its incarnation
  /// as an exact preconditioner). By default we retain the fluid on solid off
  /// diagonal blocks.
  //=============================================================================
  template<typename MATRIX>
  class SimpleFSIPreconditioner : public BlockPreconditioner<MATRIX>
  {
  public:
    /// Constructor.
    SimpleFSIPreconditioner() : BlockPreconditioner<MATRIX>()
    {
      // set the mesh pointers
      Navier_stokes_mesh_pt = 0;
      Wall_mesh_pt = 0;
      this->set_nmesh(2);

      // Default setting: Retain fluid on solid
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = true;

      // Null the preconditioner pointer (inexact solver)
      Preconditioner_pt = 0;

      // Initially assume that there are no multiple element types in
      // the same mesh.
      Allow_multiple_element_type_in_navier_stokes_mesh = false;
      Allow_multiple_element_type_in_wall_mesh = false;
    }


    /// Destructor: Clean up
    ~SimpleFSIPreconditioner()
    {
      // Wiping preconditioner (inexact solver) wipes the stored
      // LU decompositon
      if (Preconditioner_pt != 0)
      {
        delete Preconditioner_pt;
        Preconditioner_pt = 0;
      }
    }


    /// Broken copy constructor
    SimpleFSIPreconditioner(const SimpleFSIPreconditioner&) = delete;


    /// Broken assignment operator
    /*void operator=(const SimpleFSIPreconditioner&) = delete;*/

    /// Setter function for the mesh containing the
    /// block-preconditionable Navier-Stokes elements.
    void set_navier_stokes_mesh(
      Mesh* mesh_pt,
      const bool& allow_multiple_element_type_in_navier_stokes_mesh = false)
    {
      // Store the mesh pointer.
      Navier_stokes_mesh_pt = mesh_pt;

      // Are there multiple elements in this mesh?
      Allow_multiple_element_type_in_navier_stokes_mesh =
        allow_multiple_element_type_in_navier_stokes_mesh;
    }

    /// Setter function for the mesh containing the
    /// block-preconditionable FSI solid elements.
    void set_wall_mesh(
      Mesh* mesh_pt,
      const bool& allow_multiple_element_type_in_wall_mesh = false)
    {
      // Store the mesh pointer
      Wall_mesh_pt = mesh_pt;

      // Are the multiple elements in this mesh?
      Allow_multiple_element_type_in_wall_mesh =
        allow_multiple_element_type_in_wall_mesh;
    }

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// Switch to block-diagonal preconditioner
    void use_block_diagonal_version()
    {
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = false;
    }

    /// Switch to block-triangular preconditioner in which
    /// action of fluid dofs onto solid equations is retained
    void use_block_triangular_version_with_fluid_on_solid()
    {
      Retain_solid_onto_fluid_terms = false;
      Retain_fluid_onto_solid_terms = true;
    }

    /// Switch to block-triangular preconditioner in which
    /// action of solid dofs onto fluid equations is retained
    void use_block_triangular_version_with_solid_on_fluid()
    {
      Retain_solid_onto_fluid_terms = true;
      Retain_fluid_onto_solid_terms = false;
    }

  private:
    /// Preconditioner (inexact solver)
    Preconditioner* Preconditioner_pt;

    /// Boolean flag used to indicate that the solid onto fluid
    /// interaction terms are to be retained
    bool Retain_solid_onto_fluid_terms;

    /// Boolean flag used to indicate that the fluid onto solid
    /// interaction terms are to be retained
    bool Retain_fluid_onto_solid_terms;

    /// Identify the required blocks: Here we only need
    /// the momentum, gradient and divergence blocks of the
    /// 2x2 block-structured fluid matrix, the 1x1 solid block
    /// and the selected FSI-off diagonals.
    virtual void identify_required_blocks(DenseMatrix<bool>& required_blocks);

    /// Pointer to the navier stokes mesh
    Mesh* Navier_stokes_mesh_pt;

    /// pointer to the solid mesh
    Mesh* Wall_mesh_pt;

    /// Flag for multiple element types in the Navier-Stokes mesh.
    bool Allow_multiple_element_type_in_navier_stokes_mesh;

    /// Flag for multiple element types in the Wall mesh
    bool Allow_multiple_element_type_in_wall_mesh;
  };


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  // FSI preconditioner member functions
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //===========================================================================
  /// Identify the required blocks: Here we only need
  /// the momentum, gradient and divergence blocks of the
  /// 2x2 block-structured fluid matrix, the 1x1 solid block
  /// and the selected FSI-off diagonals.
  //===========================================================================
  template<typename MATRIX>
  void SimpleFSIPreconditioner<MATRIX>::identify_required_blocks(
    DenseMatrix<bool>& required_blocks)
  {
    // find number of block types
    unsigned n_dof = this->nblock_types();

    // Initialise all blocks to false
    for (unsigned i = 0; i < n_dof; i++)
    {
      for (unsigned j = 0; j < n_dof; j++)
      {
        required_blocks(i, j) = false;
      }
    }

    // Fluid: Only need momentum, gradient and divergence blocks
    required_blocks(0, 0) = true;
    required_blocks(1, 0) = true;
    required_blocks(0, 1) = true;

    // Always retain the solid block
    required_blocks(2, 2) = true;

    // Switch on the required off-diagonals
    if (Retain_solid_onto_fluid_terms)
    {
      required_blocks(0, 2) = true;
      required_blocks(1, 2) = true;
    }
    if (Retain_fluid_onto_solid_terms)
    {
      required_blocks(2, 0) = true;
      required_blocks(2, 1) = true;
      if (Retain_solid_onto_fluid_terms)
      {
        std::ostringstream error_message;
        error_message << "Can't retain all off-diagonal blocks!\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
  }


  //=============================================================================
  /// Setup the preconditioner: Copy the upper/lower triangular
  /// block matrices back into a big matrix (with the entries
  /// re-ordered relative to the original Jacobian matrix).
  //=============================================================================
  template<typename MATRIX>
  void SimpleFSIPreconditioner<MATRIX>::setup()
  {
    // Clean up memory
    if (Preconditioner_pt != 0)
    {
      delete Preconditioner_pt;
      Preconditioner_pt = 0;
    }
#ifdef PARANOID
    if (Navier_stokes_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Pointer to fluid mesh hasn't been set!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (Wall_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Pointer to solid mesh hasn't been set!\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // setup the meshes
    this->set_mesh(0,
                   Navier_stokes_mesh_pt,
                   Allow_multiple_element_type_in_navier_stokes_mesh);
    this->set_mesh(1, Wall_mesh_pt, Allow_multiple_element_type_in_wall_mesh);

    // get the number of fluid dofs from the first element in the mesh
    unsigned n_fluid_dof = this->ndof_types_in_mesh(0);
    unsigned n_dof = n_fluid_dof + this->ndof_types_in_mesh(1);

    // this fsi preconditioner has two types of DOF fluid dofs and solid dofs
    Vector<unsigned> dof_to_block_map(n_dof, 0);
    dof_to_block_map[n_fluid_dof - 1] = 1; // pressure
    for (unsigned i = n_fluid_dof; i < n_dof; i++) // solid
    {
      dof_to_block_map[i] = 2;
    }

    // Set up the blocks look up schemes
    this->block_setup(dof_to_block_map);

    // find number of block types
    n_dof = this->nblock_types();

    // Create matrix that indicates which blocks are required
    DenseMatrix<bool> required_blocks(n_dof, n_dof);

    // Identify required blocks
    identify_required_blocks(required_blocks);

    VectorMatrix<BlockSelector> selected_blocks(n_dof, n_dof);

    for (unsigned dof_i = 0; dof_i < n_dof; dof_i++)
    {
      for (unsigned dof_j = 0; dof_j < n_dof; dof_j++)
      {
        selected_blocks[dof_i][dof_j].select_block(dof_i, dof_j, false, 0);

        if (required_blocks(dof_i, dof_j))
        {
          selected_blocks[dof_i][dof_j].want_block();
        }
      }
    }

    CRDoubleMatrix P_matrix = this->get_concatenated_block(selected_blocks);

    // Setup preconditioner (i.e. inexact solver) -- does the LU decomposition
    Preconditioner_pt = new SuperLUPreconditioner;
    Preconditioner_pt->setup(&P_matrix);
  }


  //======================================================================
  /// Apply preconditioner to Vector r
  //======================================================================
  template<typename MATRIX>
  void SimpleFSIPreconditioner<MATRIX>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // create a temporary vector to hold the result of preconditioning
    DoubleVector temp_vec;

    // get the reordered vector
    this->get_block_ordered_preconditioner_vector(r, temp_vec);

    // apply preconditioner to z and store in r
    Preconditioner_pt->preconditioner_solve(temp_vec, temp_vec);

    // copy the solution back
    this->return_block_ordered_preconditioner_vector(temp_vec, z);
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


} // namespace oomph

#endif

/*----------------------------------------------------------------
   Interface to distributed SuperLU, written by JWB and heavily
   revised by PM by adapting code in the SuperLU distribution. The
   function
              superlu_dist_distributed_matrix(...)
   is based on the file
                      EXAMPLE/pddrive2.c
   and the function
                 superlu_dist_global_matrix(...)
   is based on the file
                   EXAMPLE/pddrive2_ABglobal.c

   To update this driver code for use with later versions of
   SuperLU_DIST look at any changes to the example drivers.

   Adapted from code found in:
   -- Distributed SuperLU routine (version 9.1.0) --
   Lawrence Berkeley National Lab, Univ. of California Berkeley.
   Nov 11, 2024
  ----------------------------------------------------------------
*/
#include <math.h>
#include <superlu_ddefs.h>
#include <superlu_enum_consts.h>


/* ================================================= */
/* Struct for the lu factors  */
/* ================================================= */
typedef struct
{
  gridinfo_t* grid;
  SuperMatrix* A;
  dScalePermstruct_t* ScalePermstruct;
  dLUstruct_t* LUstruct;
  dSOLVEstruct_t* SOLVEstruct;
  superlu_dist_options_t* options;
  int_t rowequ;
  int_t colequ;
  double anorm;
} superlu_dist_data;


/* ================================================= */
/* Can't think of any other way to store the memory  */
/* stats... (PM)                                     */
/* ================================================= */
struct MemoryStatisticsStorage
{
  // Storage for the memory stats
  superlu_dist_mem_usage_t Memory_usage;

  // Boolean
  int Memory_usage_has_been_recorded;
} symbolic_memory_statistics_storage;

/* ========================================================================= */
/* Helper to record memory usage*/
/* ========================================================================= */
double get_lu_factor_memory_usage_in_bytes_dist()
{
  // If the LU decomposition has been stored
  if (symbolic_memory_statistics_storage.Memory_usage_has_been_recorded == 1)
  {
    return symbolic_memory_statistics_storage.Memory_usage.for_lu;
  }
  else
  {
    return 0.0;
  }
} // End of get_lu_factor_memory_usage_in_bytes

/* ========================================================================= */
/* Helper to record memory usage*/
/* ========================================================================= */
double get_total_memory_usage_in_bytes_dist()
{
  // If the LU decomposition has been stored
  if (symbolic_memory_statistics_storage.Memory_usage_has_been_recorded == 1)
  {
    return symbolic_memory_statistics_storage.Memory_usage.total;
  }
  else
  {
    return 0.0;
  }
} // End of get_total_memory_usage_in_bytes

/* ========================================================================= */
/* Helper to record memory usage*/
/* ========================================================================= */
void get_memory_usage_in_bytes_dist(double* lu_factor_memory,
                                    double* total_memory)
{
  (*lu_factor_memory) = symbolic_memory_statistics_storage.Memory_usage.for_lu;
  (*total_memory) = symbolic_memory_statistics_storage.Memory_usage.total;
}

//=============================================================================
// helper method - just calls the superlu method dCompRow_to_CompCol to convert
// the c-style vectors of a cr matrix to a cc matrix
//=============================================================================
void superlu_cr_to_cc(int nrow,
                      int ncol,
                      int nnz,
                      double* cr_values,
                      int* cr_index,
                      int* cr_start,
                      double** cc_values,
                      int** cc_index,
                      int** cc_start)
{
  dCompRow_to_CompCol(nrow,
                      ncol,
                      nnz,
                      cr_values,
                      cr_index,
                      cr_start,
                      cc_values,
                      cc_index,
                      cc_start);
}

/*----------------------------------------------------------------
  Enumeration to select setup, solve, or clean up.
  ----------------------------------------------------------------*/
typedef enum
{
  SETUP_PHASE = 1,
  SOLVE_PHASE = 2,
  CLEAN_UP_PHASE = 3,
} opt_flag_t;

/*----------------------------------------------------------------
   Bridge to distributed SuperLU with distributed memory (version 2.0).
   Requires input of system matrix in compressed row form.

   Parameters:
   op_flag    = int specifies the operation:
                  1, performs LU decomposition for the first time
                  2, performs triangular solve
                  3, free all the storage in the end
   - n = size of system (square matrix)
   - nnz = # of nonzero entries
   - values = 1D C array of nonzero entries
   - row_index =  1D C array of row indices
   - column_start =  1D C array of column start indices
   - b = 1D C array representing the rhs vector, is overwritten
         by solution.
   - nprow = # of rows in process grid
   - npcol = # of columns in process grid
   - doc = 0/1 for doc/no doc
   - data  = pointer to structure to contain LU solver data.
             If *opt_flag == 1, it is an output. Otherwise, it it an input.

   Return value for *info:
           = 0: successful exit
           > 0: if *info = i, and i is
               <= A->ncol: U(i,i) is exactly zero. The factorization has
                  been completed, but the factor U is exactly singular,
                  so the solution could not be computed.
               > A->ncol: number of bytes allocated when memory allocation
                  failure occurred, plus A->ncol.
           < 0: some other error
  ----------------------------------------------------------------
*/
void superlu_dist_distributed_matrix(opt_flag_t opt_flag,
                                     int allow_permutations,
                                     int n,
                                     int nnz_local,
                                     int nrow_local,
                                     int first_row,
                                     double* values,
                                     int* col_index,
                                     int* row_start,
                                     double* b,
                                     int nprow,
                                     int npcol,
                                     int doc,
                                     void** data,
                                     int* info,
                                     MPI_Comm comm)
{
  /* Some SuperLU structures */
  superlu_dist_options_t* options;
  SuperLUStat_t stat;
  SuperMatrix* A;
  dScalePermstruct_t* ScalePermstruct;
  dLUstruct_t* LUstruct;
  dSOLVEstruct_t* SOLVEstruct;
  gridinfo_t* grid;

  /* Structure to hold SuperLU structures and data */
  superlu_dist_data* superlu_data;

  int_t* perm_r; /* row permutations from partial pivoting */
  int_t* perm_c; /* column permutation vector */
  int_t* etree; /* elimination tree */
  int_t* rowptr = NULL;
  int_t* colind; /* Local A in NR*/
  int_t job, rowequ, colequ, iinfo, need_value, i, j, irow, icol;
  int_t m_loc, fst_row, nnz, nnz_loc; /* dist_mem_use; */
  int_t *colptr, *rowind;
  Glu_persist_t* Glu_persist;
  Glu_freeable_t* Glu_freeable = NULL;

  float GA_mem_use = 0.0; /* memory usage by global A */
  float dist_mem_use = 0.0; /* memory usage during distribution */
  superlu_dist_mem_usage_t num_mem_usage, symb_mem_usage;
  int64_t nnzLU;
  int_t nsupers;

  /* Other stuff needed by SuperLU */
  double* berr = NULL;
  double* a = NULL;
  double *X, *b_col;
  double* B = b;
  double *C, *R, *C1, *R1, *x_col; /* *bcol, */
  double amax, t, colcnd, rowcnd;
  double anorm = 0.0;
  char equed[1], norm[1];
  int ldx; /* LDA for matrix X (local). */
  // static superlu_dist_mem_usage_t symb_mem_usage;
  fact_t Fact;
  int_t Equil, factored, notran, permc_spec;

  /* We're only doing single rhs problems
     note: code will need modifying to deal with
     multiple rhs (see function pdgssvx) */

  /* Square matrix */
  int m = n;

  /* Set 'Leading dimension' of rhs vector */
  int ldb = n;

  /* Initialize the statistics variables. */
  PStatInit(&stat);

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     SET UP GRID, FACTORS, ETC
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == SETUP_PHASE)
  {
    /* Allocate data structure to store data between calls to this function */
    superlu_data =
      (superlu_dist_data*)SUPERLU_MALLOC(sizeof(superlu_dist_data));

    /* Initialize the superlu process grid. */
    grid = (gridinfo_t*)SUPERLU_MALLOC(sizeof(gridinfo_t));
    superlu_gridinit(comm, nprow, npcol, grid);
    superlu_data->grid = grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol) return;

    /* Allocate memory for SuperLU_DIST structures */
    options =
      (superlu_dist_options_t*)SUPERLU_MALLOC(sizeof(superlu_dist_options_t));
    A = (SuperMatrix*)SUPERLU_MALLOC(sizeof(SuperMatrix));
    ScalePermstruct =
      (dScalePermstruct_t*)SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
    LUstruct = (dLUstruct_t*)SUPERLU_MALLOC(sizeof(dLUstruct_t));
    SOLVEstruct = (dSOLVEstruct_t*)SUPERLU_MALLOC(sizeof(dSOLVEstruct_t));

    /* Create SuperMatrix from compressed row representation */
    dCreate_CompRowLoc_Matrix_dist(A,
                                   m,
                                   n,
                                   nnz_local,
                                   nrow_local,
                                   first_row,
                                   values,
                                   col_index,
                                   row_start,
                                   SLU_NR_loc,
                                   SLU_D,
                                   SLU_GE);

    /* Set the default options */
    set_default_options_dist(options);

    /* Is the matrix transposed (NOTRANS or TRANS)? */
    options->Trans = NOTRANS;

    /* Row permutations (NATURAL [= do nothing],     */
    /*                   LargeDiag_MC64 [default], ...)?  */
    /*    options->RowPerm=NATURAL; */
    options->RowPerm = LargeDiag_MC64; /* hierher used to be LargeDiag */

    /* note: LargeDiag_HWPM seems to be an opption too */
    /* Column permutations (NATURAL [= do nothing],      */
    /*                      MMD_AT_PLUS_A [default],...) */
    options->ColPerm = MMD_AT_PLUS_A;

    /* Use natural ordering instead? */
    if (allow_permutations == 0)
    {
      options->ColPerm = NATURAL;
      options->RowPerm = NOROWPERM;
    }

    /*    printf("\n\n\nSWITCHING OFF EQUILIBRATION\n\n\n"); */
    /*    options->Equil=NO; */

    /* Iterative refinement (essential as far as I can tell).*/
    /* Can be "NO" or "DOUBLE"*/
    options->IterRefine = SLU_DOUBLE;

    /* Specifies whether to replace the tiny diagonals by sqrt(eps)*||A|| during
     * LU factorization. */
    options->ReplaceTinyPivot = YES;

    /* Print stats during solve? */
    if (doc == 0)
    {
      options->PrintStat = YES;
    }
    else
    {
      options->PrintStat = NO;
    }

    /* Doc output on process 0 if required: */
    if ((!iam) && (doc == 0))
    {
      printf("\nPerforming SuperLU_DIST setup\n");
      printf("Process grid\t%d X %d\n", grid->nprow, grid->npcol);
      print_options_dist(options);
    }

    /* Initialize ScalePermstruct and LUstruct. */
    dScalePermstructInit(m, n, ScalePermstruct);
    dLUstructInit(n, LUstruct);

    /* Call the linear equation solver but only perform the LU factorisation. */
    int nrhs = 0;
    pdgssvx(options,
            A,
            ScalePermstruct,
            b,
            ldb,
            nrhs,
            grid,
            LUstruct,
            SOLVEstruct,
            berr,
            &stat,
            info);

    /* Indicate that A has now been factorised. */
    options->Fact = FACTORED;

    /* Print the statistics.

      PM: A driver can hang if this is only executed when on the root processor
      so make sure not to add "&& (!iam)" to the if condition!
    */
    if (doc == 0)
    {
      printf("\nStats after distributed setup....\n");
      PStatPrint(options, &stat, grid);
    }

    /* ------------------------------------------------------------
       Set up data structure.
       ------------------------------------------------------------*/
    superlu_data->A = A;
    superlu_data->options = options;
    superlu_data->ScalePermstruct = ScalePermstruct;
    superlu_data->LUstruct = LUstruct;
    superlu_data->SOLVEstruct = SOLVEstruct;
    superlu_data->colequ = colequ;
    superlu_data->rowequ = rowequ;
    superlu_data->anorm = anorm;
    *data = superlu_data;
  } /* End of setup */


  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM A SOLVE_PHASE
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == SOLVE_PHASE)
  {
    /* Get pointer to the grid */
    superlu_data = (superlu_dist_data*)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    if ((doc == 0) && (!iam))
    {
      printf("\nPerforming SuperLU_DIST solve\n");
    }

    /* ------------------------------------------------------------
       Set other  pointers to data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    options = superlu_data->options;
    ScalePermstruct = superlu_data->ScalePermstruct;
    LUstruct = superlu_data->LUstruct;
    SOLVEstruct = superlu_data->SOLVEstruct;
    colequ = superlu_data->colequ;
    rowequ = superlu_data->rowequ;
    anorm = superlu_data->anorm;

    /* Solving for a single RHS vector */
    int nrhs = 1;

    if (!(berr = doubleMalloc_dist(nrhs)))
    {
      ABORT("Malloc fails for berr[].");
    }

    /* Call the linear solver to perform the back sub. */
    pdgssvx(options,
            A,
            ScalePermstruct,
            b,
            ldb,
            nrhs,
            grid,
            LUstruct,
            SOLVEstruct,
            berr,
            &stat,
            info);

    SUPERLU_FREE(berr);
  } /* End of solve */

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM CLEAN UP OF MEMORY
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == CLEAN_UP_PHASE)
  {
    /* Get pointer to the process grid */
    superlu_data = (superlu_dist_data*)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      goto out;
    }
    if ((doc == 0) && (!iam))
    {
      printf("\nCleaning up memory allocated for SuperLU_DIST\n");
    }

    /* ------------------------------------------------------------
       Set pointers to the data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    options = superlu_data->options;
    ScalePermstruct = superlu_data->ScalePermstruct;
    LUstruct = superlu_data->LUstruct;
    SOLVEstruct = superlu_data->SOLVEstruct;

    /* -------------------------------
       Set the other pointers required
       -------------------------------*/
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;

    /* Local control paramaters */
    Fact = options->Fact;
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);

    /* Deallocate R and/or C if it was not used. */
    if (Equil && Fact != SamePattern_SameRowPerm)
    {
      switch (ScalePermstruct->DiagScale)
      {
        case NOEQUIL:
          SUPERLU_FREE(R);
          SUPERLU_FREE(C);
          break;
        case ROW:
          SUPERLU_FREE(C);
          break;
        case COL:
          SUPERLU_FREE(R);
          break;
        default:
          break;
      }
    }

    /*  Free storage */
    PStatFree(&stat);
    // Destroy_CompRowLoc_Matrix_dist(&A);
    dScalePermstructFree(ScalePermstruct);
    dDestroy_LU(n, grid, LUstruct);
    dLUstructFree(LUstruct);
    dSolveFinalize(options, SOLVEstruct);

    // Only destroy the store part of the matrix
    Destroy_SuperMatrix_Store_dist(A);

    /* Deallocate memory */
    SUPERLU_FREE(A);
    SUPERLU_FREE(ScalePermstruct);
    SUPERLU_FREE(LUstruct);
    SUPERLU_FREE(SOLVEstruct);
    SUPERLU_FREE(options);

    /*  Release the superlu process grid. */
  out:
    superlu_gridexit(grid);

    SUPERLU_FREE(grid);
    SUPERLU_FREE(superlu_data);
  }
  return;
}


/*----------------------------------------------------------------
   Bridge to distributed SuperLU with distributed memory (version 2.0).
   Requires input of system matrix in compressed row form.

   Parameters:
   op_flag    = int specifies the operation:
                  1, performs LU decomposition for the first time
                  2, performs triangular solve
                  3, free all the storage in the end
   - n = size of system (square matrix)
   - nnz = # of nonzero entries
   - values = 1D C array of nonzero entries
   - row_index =  1D C array of row indices
   - column_start =  1D C array of column start indices
   - b = 1D C array representing the rhs vector, is overwritten
         by solution.
   - nprow = # of rows in process grid
   - npcol = # of columns in process grid
   - doc = 0/1 for doc/no doc
   - data  = pointer to structure to contain LU solver data.
             If *opt_flag == 1, it is an output. Otherwise, it it an input.

   Return value of *info:
           = 0: successful exit
           > 0: if *info = i, and i is
               <= A->ncol: U(i,i) is exactly zero. The factorization has
                  been completed, but the factor U is exactly singular,
                  so the solution could not be computed.
               > A->ncol: number of bytes allocated when memory allocation
                  failure occurred, plus A->ncol.
           < 0: some other error
  ----------------------------------------------------------------
*/
void superlu_dist_global_matrix(opt_flag_t opt_flag,
                                int allow_permutations,
                                int n_in,
                                int nnz_in,
                                double* values,
                                int* row_index,
                                int* col_start,
                                double* b,
                                int nprow,
                                int npcol,
                                int doc,
                                void** data,
                                int* info,
                                MPI_Comm comm)
{
  /* Some SuperLU structures */
  superlu_dist_options_t* options;
  SuperLUStat_t stat;
  SuperMatrix* A;
  dScalePermstruct_t* ScalePermstruct;
  dLUstruct_t* LUstruct;
  gridinfo_t* grid;

  /* Structure to hold SuperLU structures and data */
  superlu_dist_data* superlu_data;

  int_t* perm_r; /* row permutations from partial pivoting */
  int_t* perm_c; /* column permutation vector */
  int_t* etree; /* elimination tree */
  int_t job, rowequ, colequ, iinfo, i, j, irow; /* , need_value */
  int_t m, n, nnz;
  int_t *colptr, *rowind;
  int_t Equil, factored, notran, permc_spec; /*, dist_mem_use; */
  NCformat* Astore;
  NCPformat* ACstore;
  Glu_persist_t* Glu_persist;
  Glu_freeable_t* Glu_freeable = NULL;

  /* Other stuff needed by SuperLU */
  double* berr = NULL;
  double *a, *X, *b_col;
  double* B = b;
  double *C, *R, *C1, *R1, *b_work, *x_col; /* *bcol, */
  double amax, t, colcnd, rowcnd;
  double anorm = 0.0;
  char equed[1], norm[1];
  int ldx; /* LDA for matrix X (local). */
  int iam;
  // static superlu_dist_mem_usage_t  num_mem_usage, symb_mem_usage;
  fact_t Fact;

  /* Square matrix */
  n = n_in;
  m = n_in;
  nnz = nnz_in;

  /* Set 'Leading dimension' of rhs vector */
  int ldb = n;


  /* Initialize the statistics variables. */
  PStatInit(&stat);

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     SET UP GRID, FACTORS, ETC
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == SETUP_PHASE)
  {
    /* Allocate data structure to store data between calls to this function */
    superlu_data =
      (superlu_dist_data*)SUPERLU_MALLOC(sizeof(superlu_dist_data));

    /* Initialize the superlu process grid. */
    grid = (gridinfo_t*)SUPERLU_MALLOC(sizeof(gridinfo_t));
    superlu_gridinit(comm, nprow, npcol, grid);
    superlu_data->grid = grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    /* Allocate memory for SuperLU_DIST structures */
    options =
      (superlu_dist_options_t*)SUPERLU_MALLOC(sizeof(superlu_dist_options_t));
    A = (SuperMatrix*)SUPERLU_MALLOC(sizeof(SuperMatrix));
    ScalePermstruct =
      (dScalePermstruct_t*)SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
    LUstruct = (dLUstruct_t*)SUPERLU_MALLOC(sizeof(dLUstruct_t));

    /* Set the default options */
    set_default_options_dist(options);

    /* Is the matrix transposed (NOTRANS or TRANS)? */
    options->Trans = NOTRANS;

    /* Row permutations (NATURAL [= do nothing],     */
    /*                   LargeDiag_MC64 [default], ...)?  */
    /*    options->RowPerm=NATURAL; */
    options->RowPerm = LargeDiag_MC64;

    /* Column permutations (NATURAL [= do nothing],      */
    /*                      MMD_AT_PLUS_A [default],...) */
    options->ColPerm = MMD_AT_PLUS_A;

    /* Use natural ordering instead? */
    if (allow_permutations == 0)
    {
      options->ColPerm = NATURAL;
      options->RowPerm = NOROWPERM;
    }

    /* Iterative refinement (essential as far as I can tell).*/
    /* Can be "NO" or "DOUBLE"*/
    options->IterRefine = SLU_DOUBLE;

    /* Print stats during solve? */
    if (doc == 0)
    {
      options->PrintStat = YES;
    }
    else
    {
      options->PrintStat = NO;
    }

    /* Doc output on process 0 if required: */
    if ((!iam) && (doc == 0))
    {
      printf("\nPerforming SuperLU_DIST setup\n");
      printf("Process grid\t%d X %d\n", grid->nprow, grid->npcol);
      print_options_dist(options);
    }

    /*  Create SuperMatrix from compressed column representation */
    dCreate_CompCol_Matrix_dist(
      A, m, n, nnz, values, row_index, col_start, SLU_NC, SLU_D, SLU_GE);

    /* Initialize ScalePermstruct and LUstruct. */
    dScalePermstructInit(m, n, ScalePermstruct);
    dLUstructInit(n, LUstruct);

    /* Call the linear equation solver. */
    int nrhs = 0;
    pdgssvx_ABglobal(options,
                     A,
                     ScalePermstruct,
                     b,
                     ldb,
                     nrhs,
                     grid,
                     LUstruct,
                     berr,
                     &stat,
                     info);

    /* Indicate that A has now been factorised. */
    options->Fact = FACTORED;

    if (*info != 0)
    {
      printf("Trouble in  pdgstrf. Info=%i\n", *info);
      if (*info > 0)
      {
        printf(
          "U(%i,%i) is exactly zero. The factorization has\n", *info, *info);
        printf("been completed, but the factor U is exactly singular,\n");
        printf("and division by zero will occur if it is used to solve a\n");
        printf("system of equations.\n");
      }
      else
      {
        printf("The %i-th argument had an illegal value.\n", *info);
      }
    }

    /* Print the statistics.

      PM: A driver can hang if this is only executed when on the root processor
      so make sure not to add "&& (!iam)" to the if condition!
    */
    if (doc == 0)
    {
      printf("\nStats after global setup....\n");
      PStatPrint(options, &stat, grid);
    }

    /* ------------------------------------------------------------
       Set up data structure.
       ------------------------------------------------------------*/
    superlu_data->A = A;
    superlu_data->options = options;
    superlu_data->ScalePermstruct = ScalePermstruct;
    superlu_data->LUstruct = LUstruct;
    superlu_data->colequ = colequ;
    superlu_data->rowequ = rowequ;
    superlu_data->anorm = anorm;
    *data = superlu_data;
  } /* End of setup */


  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM A SOLVE_PHASE
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == SOLVE_PHASE)
  {
    /* Get pointer to the grid */
    superlu_data = (superlu_dist_data*)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    if ((doc == 0) && (!iam))
    {
      printf("\nPerforming SuperLU_DIST solve\n");
    }

    /* ------------------------------------------------------------
       Set other  pointers to data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    options = superlu_data->options;
    ScalePermstruct = superlu_data->ScalePermstruct;
    LUstruct = superlu_data->LUstruct;
    colequ = superlu_data->colequ;
    rowequ = superlu_data->rowequ;
    anorm = superlu_data->anorm;

    /* Solving for a single RHS vector */
    int nrhs = 1;

    if (!(berr = doubleMalloc_dist(nrhs)))
    {
      ABORT("Malloc fails for berr[].");
    }

    /* Call the linear solver to perform the back sub. */
    pdgssvx_ABglobal(options,
                     A,
                     ScalePermstruct,
                     b,
                     ldb,
                     nrhs,
                     grid,
                     LUstruct,
                     berr,
                     &stat,
                     info);

    SUPERLU_FREE(berr);
  } /* End of solve */

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM CLEAN UP OF MEMORY
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag == CLEAN_UP_PHASE)
  {
    /* Get pointer to the process grid */
    superlu_data = (superlu_dist_data*)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      goto out;
    }
    if ((doc == 0) && (!iam))
    {
      printf("\nCleaning up memory allocated for SuperLU_DIST\n");
    }

    /* ------------------------------------------------------------
       Set pointers to the data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    options = superlu_data->options;
    ScalePermstruct = superlu_data->ScalePermstruct;
    LUstruct = superlu_data->LUstruct;

    /* -------------------------------
       Set the other pointers required
       -------------------------------*/
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;

    /* Local control paramaters */
    Fact = options->Fact;
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);

    /* Deallocate storage. */
    if (Equil && Fact != SamePattern_SameRowPerm)
    {
      switch (ScalePermstruct->DiagScale)
      {
        case NOEQUIL:
          SUPERLU_FREE(R);
          SUPERLU_FREE(C);
          break;
        case ROW:
          SUPERLU_FREE(C);
          break;
        case COL:
          SUPERLU_FREE(R);
          break;
        default:
          break;
      }
    }

    /*  Free storage */
    PStatFree(&stat);
    // Destroy_CompCol_Matrix_dist(&A);
    dDestroy_LU(n, grid, LUstruct);
    dScalePermstructFree(ScalePermstruct);
    dLUstructFree(LUstruct);

    //  Only destroy the store part of the matrix
    Destroy_SuperMatrix_Store_dist(A);

    /* Deallocate memory */
    SUPERLU_FREE(A);
    SUPERLU_FREE(ScalePermstruct);
    SUPERLU_FREE(LUstruct);
    SUPERLU_FREE(options);

    /*  Release the superlu process grid. */
  out:
    superlu_gridexit(grid);

    SUPERLU_FREE(grid);
    SUPERLU_FREE(superlu_data);
  }
  return;
}

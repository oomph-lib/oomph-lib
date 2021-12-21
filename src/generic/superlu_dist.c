

/*----------------------------------------------------------------
   Interface to distributed SuperLU, created by JWB by adapting
   code in the superlu distribution, i.e. files /SRC/pdgssvx.c
   (function pdgssvx solves a system of linear equations) and
   /EXAMPLE/pddrive.c demo (a driver program to illustrate how to
   use pdgssvx).
   Essentially the code below performs the same functions as
   pdgssvx in a modified order which allows resolves. Much of the
   code taken from pdgssvx remains essentially unchanged other
   than changing the layout to match the oomph-lib standard. Comments
   from the original code have been left unchanged whenever possible
   to help match this code with that found in pdgssvx.c

   To update this driver code for use with later versions of
   Distributed SuperLU I suggest first looking at changes (if any) to
   the two distributed SuperLU files, and then making the corresponding
   changes to the code below.

   Adapted from code found in:
   -- Distributed SuperLU routine (version 2.0) --
   Lawrence Berkeley National Lab, Univ. of California Berkeley.
   March 15, 2003

  ----------------------------------------------------------------
*/
#include <math.h>
#ifdef USING_OOMPH_SUPERLU_DIST
#include "oomph_superlu_dist_3.0.h"
#else
#include<superlu_defs.h>
#include<superlu_ddefs.h>
#include<Cnames.h>
#include<machines.h>
#include<psymbfact.h>
#include<supermatrix.h>
#include<old_colamd.h>
#include<util_dist.h>
#endif


/* ================================================= */
/* Struct for the lu factors  */
/* ================================================= */
typedef struct
{
  gridinfo_t* grid;
  SuperMatrix* A;
  SuperMatrix* AC;
  ScalePermstruct_t* ScalePermstruct;
  LUstruct_t* LUstruct;
  SOLVEstruct_t* SOLVEstruct;
  superlu_options_t* options;
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
  mem_usage_t Memory_usage;

  // Boolean
  int Memory_usage_has_been_recorded;
} symbolic_memory_statistics_storage;

/* ========================================================================= */
/* Helper to record memory usage*/
/* ========================================================================= */
double get_lu_factor_memory_usage_in_bytes_dist()
{
  // If the LU decomposition has been stored
  if (symbolic_memory_statistics_storage.Memory_usage_has_been_recorded==1)
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
  if (symbolic_memory_statistics_storage.Memory_usage_has_been_recorded==1)
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
  (*lu_factor_memory)=symbolic_memory_statistics_storage.Memory_usage.for_lu;
  (*total_memory)=symbolic_memory_statistics_storage.Memory_usage.total;
}

//=============================================================================
// helper method - just calls the superlu method dCompRow_to_CompCol to convert
// the c-style vectors of a cr matrix to a cc matrix
//=============================================================================
void superlu_cr_to_cc(int nrow, int ncol, int nnz, double* cr_values,
                      int* cr_index, int* cr_start, double** cc_values,
                      int** cc_index, int** cc_start)
{
  dCompRow_to_CompCol(nrow,ncol,nnz,cr_values,cr_index,cr_start,cc_values,
                      cc_index,cc_start);
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
void superlu_dist_distributed_matrix(int opt_flag, int allow_permutations,
                                     int n, int nnz_local,
                                     int nrow_local,int first_row,
                                     double *values, int *col_index,
                                     int *row_start, double *b,
                                     int nprow, int npcol,
                                     int doc, void **data, int *info,
                                     MPI_Comm comm)
{
  /* Some SuperLU structures */
  superlu_options_t *options;
  SuperLUStat_t stat;
  SuperMatrix *A;
  ScalePermstruct_t *ScalePermstruct;
  LUstruct_t *LUstruct;
  SOLVEstruct_t *SOLVEstruct;
  gridinfo_t *grid;

  /* Structure to hold SuperLU structures and data */
  superlu_dist_data *superlu_data;

  int_t *perm_r; /* row permutations from partial pivoting */
  int_t *perm_c; /* column permutation vector */
  int_t *etree;  /* elimination tree */
  int_t *rowptr=NULL; int_t *colind;  /* Local A in NR*/
  int_t job, rowequ, colequ, iinfo, need_value, i, j, irow, icol;
  int_t m_loc, fst_row, nnz, nnz_loc; /* dist_mem_use; */
  int_t *colptr, *rowind;
  NRformat_loc *Astore;
  SuperMatrix GA;      /* Global A in NC format */
  NCformat *GAstore;
  double *a_GA=NULL;
  SuperMatrix GAC;      /* Global A in NCP format (add n end pointers) */
  NCPformat *GACstore;
  Glu_persist_t *Glu_persist;
  Glu_freeable_t *Glu_freeable=NULL;

  /* Other stuff needed by SuperLU */
  double  *berr=NULL;
  double *a=NULL; double *X, *b_col;
  double *B=b;
  double *C, *R, *C1, *R1, *x_col; /* *bcol, */
  double amax, t, colcnd, rowcnd; double anorm=0.0;
  char equed[1], norm[1];
  int ldx;  /* LDA for matrix X (local). */
  //static mem_usage_t symb_mem_usage;
  fact_t Fact;
  int_t Equil, factored, notran, permc_spec;

  /* We're only doing single rhs problems
     note: code will need modifying to deal with
     multiple rhs (see function pdgssvx) */

  int nrhs=1;

  /* Square matrix */
  int m=n;

  /* Set 'Leading dimension' of rhs vector */
  int ldb=n;

  /* Initialize the statistics variables. */
  PStatInit(&stat);

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     SET UP GRID, FACTORS, ETC
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==1)
  {
    /* Allocate data structure to store data between calls to this function */
    superlu_data =
      (superlu_dist_data *) SUPERLU_MALLOC(sizeof(superlu_dist_data));

    /* Initialize the superlu process grid. */
    grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));
    superlu_gridinit(comm, nprow, npcol, grid);
    superlu_data->grid = grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol) return;

    /* Allocate memory for SuperLU_DIST structures */
    options = (superlu_options_t *) SUPERLU_MALLOC(sizeof(superlu_options_t));
    A = (SuperMatrix *) SUPERLU_MALLOC(sizeof(SuperMatrix));
    ScalePermstruct =
      (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
    LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
    SOLVEstruct = (SOLVEstruct_t *) SUPERLU_MALLOC(sizeof(SOLVEstruct_t));

    /* Create SuperMatrix from compressed row representation */
    dCreate_CompRowLoc_Matrix_dist(A,m,n,nnz_local,nrow_local,first_row,
                                   values,col_index,row_start,
                                   SLU_NR_loc,SLU_D,SLU_GE);

    /* Set the default options */
    set_default_options_dist(options);

    /* Is the matrix transposed (NOTRANS or TRANS)? */
    options->Trans = NOTRANS;

    /* Row permutations (NATURAL [= do nothing],     */
    /*                   LargeDiag [default], ...)?  */
    /*    options->RowPerm=NATURAL; */
    options->RowPerm=LargeDiag;

    /* Column permutations (NATURAL [= do nothing],      */
    /*                      MMD_AT_PLUS_A [default],...) */
    options->ColPerm=MMD_AT_PLUS_A;

    /* Use natural ordering instead? */
    if (allow_permutations==0)
    {
      options->ColPerm=NATURAL;
      options->RowPerm=NATURAL;
    }

    /*    printf("\n\n\nSWITCHING OFF EQUILIBRATION\n\n\n"); */
    /*    options->Equil=NO; */

    /* Iterative refinement (essential as far as I can tell).*/
    /* Can be "NO" or "DOUBLE"*/
    options->IterRefine = SLU_DOUBLE;

    /* Print stats during solve? */
    if (doc==0)
    {
      options->PrintStat = YES;
    }
    else
    {
      options->PrintStat = NO;
    }


    /* Doc output on process 0 if required: */
    if ((!iam)&&(doc==0))
    {
      printf("\nPerforming SuperLU_DIST setup\n");
      printf("Process grid\t%d X %d\n", grid->nprow, grid->npcol);
      print_options_dist(options);
    }

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, ScalePermstruct);
    LUstructInit(m, n, LUstruct);

    /* Initialization. */
    Glu_persist = LUstruct->Glu_persist;
    Astore = (NRformat_loc *) A->Store;
    nnz_loc = Astore->nnz_loc;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    a = Astore->nzval;
    rowptr = Astore->rowptr;
    colind = Astore->colind;

    job = 5;
    Fact = options->Fact;
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);
    notran = (options->Trans == NOTRANS);
    rowequ = colequ = FALSE;
    if (factored || (Fact == SamePattern_SameRowPerm && Equil))
    {
      rowequ = (ScalePermstruct->DiagScale == ROW) ||
               (ScalePermstruct->DiagScale == BOTH);
      colequ = (ScalePermstruct->DiagScale == COL) ||
               (ScalePermstruct->DiagScale == BOTH);
    }
    else
    {
      rowequ = colequ = FALSE;
    }

    /* Test the control parameters etc. */
    *info = 0;
    if (Fact < 0 || Fact > FACTORED)
    {
      *info = -1;
    }
    else if (options->RowPerm < 0 || options->RowPerm > MY_PERMR)
    {
      *info = -1;
    }
    else if (options->ColPerm < 0 || options->ColPerm > MY_PERMC)
    {
      *info = -1;
    }
    else if (options->IterRefine < 0 || options->IterRefine > SLU_EXTRA)
    {
      *info = -1;
    }
    else if (options->IterRefine == SLU_EXTRA)
    {
      *info = -1;
      fprintf(stderr, "Extra precise iterative refinement yet to support.\n");
    }
    else if (A->nrow != A->ncol || A->nrow < 0 || A->Stype != SLU_NR_loc
             || A->Dtype != SLU_D || A->Mtype != SLU_GE)
    {
      *info = -2;
    }
    else if (ldb < m_loc)
    {
      *info = -5;
    }
    else if (nrhs < 0)
    {
      *info = -6;
    }
    if (*info)
    {
      printf("Trouble in  pdgstrf. Info=%i\n",*info);
      if (*info==-1)
      {
        printf("Error in options.\n");
      }
      else if (*info==-2)
      {
        printf("Error in matrix.\n");
      }
      else if (*info==-5)
      {
        printf("ldb < m_loc\n");
      }
      else if (*info==-6)
      {
        printf("nrhs < 0\n");
      }
      return;
    }

    /* The following arrays are replicated on all processes. */
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    etree = LUstruct->etree;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;

    /* Allocate storage. */
    if (Equil)
    {
      /* Not factored & ask for equilibration */
      /* Allocate storage if not done so before. */
      switch (ScalePermstruct->DiagScale)
      {
      case NOEQUIL:
        if (!(R = (double *) doubleMalloc_dist(m)))
          ABORT("Malloc fails for R[].");
        if (!(C = (double *) doubleMalloc_dist(n)))
          ABORT("Malloc fails for C[].");
        ScalePermstruct->R = R;
        ScalePermstruct->C = C;
        break;
      case ROW:
        if (!(C = (double *) doubleMalloc_dist(n)))
          ABORT("Malloc fails for C[].");
        ScalePermstruct->C = C;
        break;
      case COL:
        if (!(R = (double *) doubleMalloc_dist(m)))
          ABORT("Malloc fails for R[].");
        ScalePermstruct->R = R;
        break;
      default:
       printf("diagscale: %i %i %i %i\n",ScalePermstruct->DiagScale,NOEQUIL,ROW,COL);
       ABORT("Never get here.");
       break;
      }
    }

    /* ------------------------------------------------------------
       Diagonal scaling to equilibrate the matrix.
       ------------------------------------------------------------*/
    if (Equil)
    {
      t = SuperLU_timer_();

      if (Fact == SamePattern_SameRowPerm)
      {
        /* Reuse R and C. */
        switch (ScalePermstruct->DiagScale)
        {
        case NOEQUIL:
          break;
        case ROW:
          irow = fst_row;
          for (j = 0; j < m_loc; ++j)
          {
            for (i = rowptr[j]; i < rowptr[j+1]; ++i)
            {
              a[i] *= R[irow];       /* Scale rows. */
            }
            ++irow;
          }
          break;
        case COL:
          for (j = 0; j < m_loc; ++j)
          {
            for (i = rowptr[j]; i < rowptr[j+1]; ++i)
            {
              icol = colind[i];
              a[i] *= C[icol];          /* Scale columns. */
            }
          }
          break;
        case BOTH:
          irow = fst_row;
          for (j = 0; j < m_loc; ++j)
          {
            for (i = rowptr[j]; i < rowptr[j+1]; ++i)
            {
              icol = colind[i];
              a[i] *= R[irow] * C[icol]; /* Scale rows and cols. */
            }
            ++irow;
          }
          break;
        }
      }
      else
      {
        /* Compute R & C from scratch */
        /* Compute the row and column scalings. */
        pdgsequ(A, R, C, &rowcnd, &colcnd, &amax, &iinfo, grid);

        /* Equilibrate matrix A if it is badly-scaled. */
        pdlaqgs(A, R, C, rowcnd, colcnd, amax, equed);

        if (lsame_(equed, "R"))
        {
          ScalePermstruct->DiagScale = rowequ = ROW;
        }
        else if (lsame_(equed, "C"))
        {
          ScalePermstruct->DiagScale = colequ = COL;
        }
        else if (lsame_(equed, "B"))
        {
          ScalePermstruct->DiagScale = BOTH;
          rowequ = ROW;
          colequ = COL;
        }
        else
          ScalePermstruct->DiagScale = NOEQUIL;
      } /* if Fact ... */

      stat.utime[EQUIL] = SuperLU_timer_() - t;
    } /* if Equil ... */

    if (!factored)
    {
      /* Skip this if already factored. */
      /*
         Gather A from the distributed compressed row format to
         global A in compressed column format.
         Numerical values are gathered only when a row permutation
         for large diagonal is sought after.
      */
      if (Fact != SamePattern_SameRowPerm)
      {
        need_value = (options->RowPerm == LargeDiag);
        pdCompRow_loc_to_CompCol_global(need_value, A, grid, &GA);
        GAstore = (NCformat *) GA.Store;
        colptr = GAstore->colptr;
        rowind = GAstore->rowind;
        nnz = GAstore->nnz;
        if (need_value) a_GA = GAstore->nzval;
        else assert(GAstore->nzval == NULL);
      }

      /* ------------------------------------------------------------
         Find the row permutation for A.
         ------------------------------------------------------------*/
      if ((int) options->RowPerm != (int) NO)
      {
        t = SuperLU_timer_();
        if (Fact != SamePattern_SameRowPerm)
        {
          if (options->RowPerm == MY_PERMR)
          {
            /* Use user's perm_r. */
            /* Permute the global matrix GA for symbfact() */
            for (i = 0; i < colptr[n]; ++i)
            {
              irow = rowind[i];
              rowind[i] = perm_r[irow];
            }
          }
          else
          {
            /* options->RowPerm == LargeDiag */
            /* Get a new perm_r[] */
            if (job == 5)
            {
              /* Allocate storage for scaling factors. */
              if (!(R1 = doubleMalloc_dist(m)))
              {
                ABORT("SUPERLU_MALLOC fails for R1[]");
              }
              if (!(C1 = doubleMalloc_dist(n)))
              {
                ABORT("SUPERLU_MALLOC fails for C1[]");
              }
            }

            if (!iam)
            {
              /* Process 0 finds a row permutation */
              dldperm(job, m, nnz, colptr, rowind, a_GA, perm_r, R1, C1);

              MPI_Bcast(perm_r, m, mpi_int_t, 0, grid->comm);
              if (job == 5 && Equil)
              {
                MPI_Bcast(R1, m, MPI_DOUBLE, 0, grid->comm);
                MPI_Bcast(C1, n, MPI_DOUBLE, 0, grid->comm);
              }
            }
            else
            {
              MPI_Bcast(perm_r, m, mpi_int_t, 0, grid->comm);
              if (job == 5 && Equil)
              {
                MPI_Bcast(R1, m, MPI_DOUBLE, 0, grid->comm);
                MPI_Bcast(C1, n, MPI_DOUBLE, 0, grid->comm);
              }
            }

            if (job == 5)
            {
              if (Equil)
              {
                for (i = 0; i < n; ++i)
                {
                  R1[i] = exp(R1[i]);
                  C1[i] = exp(C1[i]);
                }

                /* Scale the distributed matrix */
                irow = fst_row;
                for (j = 0; j < m_loc; ++j)
                {
                  for (i = rowptr[j]; i < rowptr[j+1]; ++i)
                  {
                    icol = colind[i];
                    a[i] *= R1[irow] * C1[icol];
                  }
                  ++irow;
                }

                /* Multiply together the scaling factors. */
                if (rowequ)
                {
                  for (i = 0; i < m; ++i)
                  {
                    R[i] *= R1[i];
                  }
                }
                else
                {
                  for (i = 0; i < m; ++i)
                  {
                    R[i] = R1[i];
                  }
                }
                if (colequ)
                {
                  for (i = 0; i < n; ++i)
                  {
                    C[i] *= C1[i];
                  }
                }
                else
                {
                  for (i = 0; i < n; ++i)
                  {
                    C[i] = C1[i];
                  }
                }

                ScalePermstruct->DiagScale = BOTH;
                rowequ = colequ = 1;

              } /* end Equil */

              /* Now permute global A to prepare for symbfact() */
              for (j = 0; j < n; ++j)
              {
                for (i = colptr[j]; i < colptr[j+1]; ++i)
                {
                  irow = rowind[i];
                  rowind[i] = perm_r[irow];
                }
              }
              SUPERLU_FREE(R1);
              SUPERLU_FREE(C1);
            }
            else
            {
              /* job = 2,3,4 */
              for (j = 0; j < n; ++j)
              {
                for (i = colptr[j]; i < colptr[j+1]; ++i)
                {
                  irow = rowind[i];
                  rowind[i] = perm_r[irow];
                } /* end for i ... */
              } /* end for j ... */
            } /* end else job ... */
          } /* end if options->RowPerm ... */

          t = SuperLU_timer_() - t;
          stat.utime[ROWPERM] = t;
        } /* end if Fact ... */
      }
      else
      {
        /* options->RowPerm == NOROWPERM */
        for (i = 0; i <m; ++i) perm_r[i] = i;
      }
    } /* end if (!factored) */

    if (!factored || options->IterRefine)
    {
      /* Compute norm(A), which will be used to adjust small diagonal. */
      if (notran)
        *(unsigned char *)norm = '1';
      else
        *(unsigned char *)norm = 'I';
      anorm = pdlangs(norm, A, grid);
    }


    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    if (!factored)
    {
      t = SuperLU_timer_();
      /*
         Get column permutation vector perm_c[], according to permc_spec:
           permc_spec = NATURAL:  natural ordering
           permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
           permc_spec = MMD_ATA:  minimum degree on structure of A'*A
           permc_spec = COLAMD:   approximate minimum degree column ordering
           permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
      */
      permc_spec = options->ColPerm;
      if (permc_spec != MY_PERMC && Fact == DOFACT)
      {
        get_perm_c_dist(iam, permc_spec, &GA, perm_c);
      }

      stat.utime[COLPERM] = SuperLU_timer_() - t;

      /* Compute the elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'
         (a.k.a. column etree), depending on the choice of ColPerm.
         Adjust perm_c[] to be consistent with a postorder of etree.
         Permute columns of A to form A*Pc'. */
      if (Fact != SamePattern_SameRowPerm)
      {
        int_t *GACcolbeg, *GACcolend, *GACrowind;

        sp_colorder(options, &GA, perm_c, etree, &GAC);

        /* Form Pc*A*Pc' to preserve the diagonal of the matrix GAC. */
        GACstore = GAC.Store;
        GACcolbeg = GACstore->colbeg;
        GACcolend = GACstore->colend;
        GACrowind = GACstore->rowind;
        for (j = 0; j < n; ++j)
        {
          for (i = GACcolbeg[j]; i < GACcolend[j]; ++i)
          {
            irow = GACrowind[i];
            GACrowind[i] = perm_c[irow];
          }
        }

        /* Perform a symbolic factorization on Pc*Pr*A*Pc' and set up the
           nonzero data structures for L & U. */
        t = SuperLU_timer_();
        if (!(Glu_freeable = (Glu_freeable_t *)
                             SUPERLU_MALLOC(sizeof(Glu_freeable_t))))
        {
          ABORT("Malloc fails for Glu_freeable.");
        }

        /* Every process does this. */
        iinfo = symbfact(options, iam, &GAC, perm_c, etree,
                         Glu_persist, Glu_freeable);

        stat.utime[SYMBFAC] = SuperLU_timer_() - t;
        if (iinfo < 0)
        {
          /* Successful return */
          QuerySpace_dist(n, -iinfo, Glu_freeable,
                          &symbolic_memory_statistics_storage.Memory_usage);
        }
        else
        {
          if (!iam)
          {
            fprintf(stderr, "symbfact() error returns %d\n", iinfo);
            exit(-1);
          }
        }
      } /* end if Fact ... */

      /* Apply column permutation to the original distributed A */
      for (j = 0; j < nnz_loc; ++j)
      {
        colind[j] = perm_c[colind[j]];
      }

      /* Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.
         NOTE: the row permutation Pc*Pr is applied internally in the
         distribution routine. */
      t = SuperLU_timer_();
      /* dist_mem_use = */
      pddistribute(Fact, n, A, ScalePermstruct,
                   Glu_freeable, LUstruct, grid);
      stat.utime[DIST] = SuperLU_timer_() - t;

      /* Deallocate storage used in symbolic factorization. */
      if (Fact != SamePattern_SameRowPerm)
      {
        iinfo = symbfact_SubFree(Glu_freeable);
        SUPERLU_FREE(Glu_freeable);
      }

      /* Perform numerical factorization in parallel. */
      t = SuperLU_timer_();
      pdgstrf(options, m, n, anorm, LUstruct, grid, &stat, info);
      stat.utime[FACT] = SuperLU_timer_() - t;

      /* Destroy GA and GAC */
      if (Fact != SamePattern_SameRowPerm)
      {
        Destroy_CompCol_Matrix_dist(&GA);
        Destroy_CompCol_Permuted_dist(&GAC);
      }
    } /* end if (!factored) */

    if (*info!=0)
    {
      printf("Trouble in  pdgstrf. Info=%i\n",*info);
      if (*info>0)
      {
        printf("U(%i,%i) is exactly zero. The factorization has\n",*info,*info);
        printf("been completed, but the factor U is exactly singular,\n");
        printf("and division by zero will occur if it is used to solve a\n");
        printf("system of equations.\n");
      }
      else
      {
        printf("The %i-th argument had an illegal value.\n", *info);
      }
    }

    /* ------------------------------------------------------------
       Initialize the solver.
       ------------------------------------------------------------*/
    if (options->SolveInitialized == NO)
    {
      dSolveInit(options, A, perm_r, perm_c, nrhs, LUstruct, grid,
                 SOLVEstruct);
    }

    if (options->IterRefine)
    {
      if (options->RefineInitialized == NO || Fact == DOFACT)
      {
        /* All these cases need to re-initialize gsmv structure */
        if (options->RefineInitialized)
        {
          pdgsmv_finalize(SOLVEstruct->gsmv_comm);
        }
        pdgsmv_init(A, SOLVEstruct->row_to_proc, grid,
                    SOLVEstruct->gsmv_comm);

        options->RefineInitialized = YES;
      }
    }

    /* Print the statistics. */
    if ((doc==0) && (!iam))
    {
      printf("\nstats after setup....\n");
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
     PERFORM A SOLVE
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==2)
  {
    /* Get pointer to the grid */
    superlu_data = (superlu_dist_data *)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    if ((doc==0)&&(!iam))
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

    /* Initialization. */
    Astore = (NRformat_loc *) A->Store;
    nnz_loc = Astore->nnz_loc;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    colind = Astore->colind;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;

    /* Local control paramaters */
    Fact = options->Fact;
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);
    notran = (options->Trans == NOTRANS);

    /* ------------------------------------------------------------
       Scale the right-hand side if equilibration was performed.
       ------------------------------------------------------------*/
    if (notran)
    {
      if (rowequ)
      {
        b_col = B;
        for (j = 0; j < nrhs; ++j)
        {
          irow = fst_row;
          for (i = 0; i < m_loc; ++i)
          {
            b_col[i] *= R[irow];
            ++irow;
          }
          b_col += ldb;
        }
      }
    }
    else if (colequ)
    {
      b_col = B;
      for (j = 0; j < nrhs; ++j)
      {
        irow = fst_row;
        for (i = 0; i < m_loc; ++i)
        {
          b_col[i] *= C[irow];
          ++irow;
        }
        b_col += ldb;
      }
    }

    /* Save a copy of the right-hand side. */
    ldx = ldb;
    if (!(X = doubleMalloc_dist(((size_t)ldx) * nrhs)))
    {
      ABORT("Malloc fails for X[]");
    }
    x_col = X;
    b_col = B;
    for (j = 0; j < nrhs; ++j)
    {
      for (i = 0; i < m_loc; ++i)
      {
        x_col[i] = b_col[i];
      }
      x_col += ldx;
      b_col += ldb;
    }


    /* ------------------------------------------------------------
       Solve the linear system.
       ------------------------------------------------------------*/
    pdgstrs(n, LUstruct, ScalePermstruct, grid, X, m_loc,
            fst_row, ldb, nrhs, SOLVEstruct, &stat, info);

    if (*info!=0)
    {
      printf("Trouble in pdgstrs. Info=%i\n",*info);
      printf("The %i-th argument had an illegal value.\n", *info);
    }

    /* ------------------------------------------------------------
       Use iterative refinement to improve the computed solution and
       compute error bounds and backward error estimates for it.
       ------------------------------------------------------------*/
    if (options->IterRefine)
    {
      /* Improve the solution by iterative refinement. */
      int_t *it, *colind_gsmv = SOLVEstruct->A_colind_gsmv;
      /*SOLVEstruct_t *SOLVEstruct1;*/  /* Used by refinement. */

      t = SuperLU_timer_();
      if (options->RefineInitialized == NO || Fact == DOFACT)
      {
        /* Save a copy of the transformed local col indices
           in colind_gsmv[]. */
        if (colind_gsmv)
        {
          SUPERLU_FREE(colind_gsmv);
        }
        if (!(it = intMalloc_dist(nnz_loc)))
        {
          ABORT("Malloc fails for colind_gsmv[]");
        }
        colind_gsmv = SOLVEstruct->A_colind_gsmv = it;
        for (i = 0; i < nnz_loc; ++i)
        {
          colind_gsmv[i] = colind[i];
        }
      }
      else if (Fact == SamePattern || Fact == SamePattern_SameRowPerm)
      {
        double at;
        int_t k, jcol, p;
        /* Swap to beginning the part of A corresponding to the
           local part of X, as was done in pdgsmv_init() */
        for (i = 0; i < m_loc; ++i)
        {
          /* Loop through each row */
          k = rowptr[i];
          for (j = rowptr[i]; j < rowptr[i+1]; ++j)
          {
            jcol = colind[j];
            p = SOLVEstruct->row_to_proc[jcol];
            if (p == iam)
            {
              /* Local */
              at = a[k]; a[k] = a[j]; a[j] = at;
              ++k;
            }
          }
        }

        /* Re-use the local col indices of A obtained from the
           previous call to pdgsmv_init() */
        for (i = 0; i < nnz_loc; ++i)
        {
          colind[i] = colind_gsmv[i];
        }
      }

      /* Storage for backward error */
      if (!(berr = doubleMalloc_dist(nrhs)))
      {
        ABORT("Malloc fails for berr[].");
      }

      pdgsrfs(n, A, anorm, LUstruct, ScalePermstruct, grid,
              B, ldb, X, ldx, nrhs, SOLVEstruct, berr, &stat, info);

      stat.utime[REFINE] = SuperLU_timer_() - t;
    }

    if (*info!=0)
    {
      printf("Trouble in pdgsrfs. Info=%i\n",*info);
      printf("The %i-th argument had an illegal value.\n", *info);
    }

    /* Print the statistics. */
    if ((doc==0) && (!iam))
    {
      printf("\nstats after solve....\n");
      PStatPrint(options, &stat, grid);
    }

    /* Permute the solution matrix B <= Pc'*X. */
    pdPermute_Dense_Matrix(fst_row, m_loc, SOLVEstruct->row_to_proc,
                           SOLVEstruct->inv_perm_c,
                           X, ldx, B, ldb, nrhs, grid);

    /* Transform the solution matrix X to a solution of the original
       system before the equilibration. */
    if (notran)
    {
      if (colequ)
      {
        b_col = B;
        for (j = 0; j < nrhs; ++j)
        {
          irow = fst_row;
          for (i = 0; i < m_loc; ++i)
          {
            b_col[i] *= C[irow];
            ++irow;
          }
          b_col += ldb;
        }
      }
    }
    else if (rowequ)
    {
      b_col = B;
      for (j = 0; j < nrhs; ++j)
      {
        irow = fst_row;
        for (i = 0; i < m_loc; ++i)
        {
          b_col[i] *= R[irow];
          ++irow;
        }
        b_col += ldb;
      }
    }

    /* Clean up memory */
    if (options->IterRefine)
    {
      SUPERLU_FREE(berr);
    }
    SUPERLU_FREE(X);

  } /* End of solve */

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM CLEAN UP OF MEMORY
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==3)
  {
    /* Get pointer to the process grid */
    superlu_data = (superlu_dist_data *)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    int iam = grid->iam;
    if (iam >= nprow * npcol) goto out;
    if ((doc==0)&&(!iam))
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
    notran = (options->Trans == NOTRANS);

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
       /* Apparently this one is ok */
       /* printf("diagscale: %i %i %i %i\n",ScalePermstruct->DiagScale,NOEQUIL,ROW,COL); */
       /* ABORT("Never get here. THIS IS THE ONE");*/
       break;
      }
    }

    /*  Free storage */
    ScalePermstructFree(ScalePermstruct);
    Destroy_LU(n, grid, LUstruct);
    LUstructFree(LUstruct);
    dSolveFinalize(options, SOLVEstruct);
    //Destroy_CompRowLoc_Matrix_dist(&A);

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

  /*  Free storage */
  PStatFree(&stat);

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
void superlu_dist_global_matrix(int opt_flag, int allow_permutations,
                                int n_in, int nnz_in,
                                double *values,
                                int *row_index, int *col_start,
                                double *b, int nprow, int npcol,
                                int doc, void **data, int *info,
                                MPI_Comm comm)
{
  /* Some SuperLU structures */
  superlu_options_t *options;
  SuperLUStat_t stat;
  SuperMatrix *A;
  SuperMatrix *AC;
  ScalePermstruct_t *ScalePermstruct;
  LUstruct_t *LUstruct;
  gridinfo_t *grid;

  /* Structure to hold SuperLU structures and data */
  superlu_dist_data *superlu_data;

  int_t *perm_r; /* row permutations from partial pivoting */
  int_t *perm_c; /* column permutation vector */
  int_t *etree;  /* elimination tree */
  int_t job, rowequ, colequ, iinfo, i, j, irow; /* , need_value */
  int_t m, n, nnz;
  int_t *colptr, *rowind;
  int_t Equil, factored, notran, permc_spec; /*, dist_mem_use; */
  NCformat *Astore;
  NCPformat *ACstore;
  Glu_persist_t *Glu_persist;
  Glu_freeable_t *Glu_freeable=NULL;

  /* Other stuff needed by SuperLU */
  double  *berr=NULL;
  double *a, *X, *b_col;
  double *B=b;
  double *C, *R, *C1, *R1, *b_work, *x_col; /* *bcol, */
  double amax, t, colcnd, rowcnd; double anorm=0.0;
  char equed[1], norm[1];
  int ldx;  /* LDA for matrix X (local). */
  int iam;
  //static mem_usage_t  num_mem_usage, symb_mem_usage;
  fact_t Fact;


  /* We're only doing single rhs problems
     note: code will need modifying to deal with
     multiple rhs (see function pdgssvx) */
  int nrhs=1;

  /* Square matrix */
  n=n_in;
  m=n_in;
  nnz=nnz_in;

  /* Set 'Leading dimension' of rhs vector */
  int ldb=n;


  /* Initialize the statistics variables. */
  PStatInit(&stat);

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     SET UP GRID, FACTORS, ETC
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==1)
  {
    /* Allocate data structure to store data between calls to this function */
    superlu_data =
      (superlu_dist_data *) SUPERLU_MALLOC(sizeof(superlu_dist_data));

    /* Initialize the superlu process grid. */
    grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));
    superlu_gridinit(comm, nprow, npcol, grid);
    superlu_data->grid = grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    /* Allocate memory for SuperLU_DIST structures */
    options = (superlu_options_t *) SUPERLU_MALLOC(sizeof(superlu_options_t));
    A = (SuperMatrix *) SUPERLU_MALLOC(sizeof(SuperMatrix));
    AC = (SuperMatrix *) SUPERLU_MALLOC(sizeof(SuperMatrix));
    ScalePermstruct =
      (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
    LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));

    /* Set the default options */
    set_default_options_dist(options);

    /* Is the matrix transposed (NOTRANS or TRANS)? */
    options->Trans = NOTRANS;

    /* Row permutations (NATURAL [= do nothing],     */
    /*                   LargeDiag [default], ...)?  */
    /*    options->RowPerm=NATURAL; */
    options->RowPerm=LargeDiag;

    /* Column permutations (NATURAL [= do nothing],      */
    /*                      MMD_AT_PLUS_A [default],...) */
    options->ColPerm=MMD_AT_PLUS_A;

    /* Use natural ordering instead? */
    if (allow_permutations==0)
    {
      options->ColPerm=NATURAL;
      options->RowPerm=NATURAL;
    }

    /* Iterative refinement (essential as far as I can tell).*/
    /* Can be "NO" or "DOUBLE"*/
    options->IterRefine = SLU_DOUBLE;

    /* Print stats during solve? */
    if (doc==0)
    {
      options->PrintStat = YES;
    }
    else
    {
      options->PrintStat = NO;
    }

    /* Doc output on process 0 if required: */
    if ((!iam)&&(doc==0))
    {
      printf("\nPerforming SuperLU_DIST setup\n");
      printf("Process grid\t%d X %d\n", grid->nprow, grid->npcol);
      print_options_dist(options);
    }


    /*  Create SuperMatrix from compressed column representation */
    dCreate_CompCol_Matrix_dist(A,m,n,nnz,values,row_index,col_start,
                                SLU_NC,SLU_D,SLU_GE);

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, ScalePermstruct);
    LUstructInit(m, n, LUstruct);

    /* Test the control parameters etc. */
    *info = 0;
    Fact = options->Fact;
    if (Fact < 0 || Fact > FACTORED)
    {
      *info = -1;
    }
    else if (options->RowPerm < 0 || options->RowPerm > MY_PERMR)
    {
      *info = -1;
    }
    else if (options->ColPerm < 0 || options->ColPerm > MY_PERMC)
    {
      *info = -1;
    }
    else if (options->IterRefine < 0 || options->IterRefine > SLU_EXTRA)
    {
      *info = -1;
    }
    else if (options->IterRefine == SLU_EXTRA)
    {
      *info = -1;
      fprintf(stderr, "Extra precise iterative refinement yet to support.\n");
    }
    else if (A->nrow != A->ncol || A->nrow < 0 ||
             A->Stype != SLU_NC || A->Dtype != SLU_D || A->Mtype != SLU_GE)
    {
      *info = -2;
    }
    else if (ldb < A->nrow)
    {
      *info = -5;
    }
    else if (nrhs < 0)
    {
      *info = -6;
    }
    if (*info)
    {
      printf("Trouble in  pdgstrf. Info=%i\n",-*info);
      if (*info==-1)
      {
        printf("Error in options.\n");
      }
      else if (*info==-2)
      {
        printf("Error in matrix.\n");
      }
      else if (*info==-5)
      {
        printf("ldb < A->nrow\n");
      }
      else if (*info==-6)
      {
        printf("nrhs < 0\n");
      }
      return;
    }

    /* Initialization. */
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);
    notran = (options->Trans == NOTRANS);
    job = 5;
    Astore = A->Store;
    nnz = Astore->nnz;
    a = Astore->nzval;
    colptr = Astore->colptr;
    rowind = Astore->rowind;
    if (factored || (Fact == SamePattern_SameRowPerm && Equil))
    {
      rowequ = (ScalePermstruct->DiagScale == ROW) ||
               (ScalePermstruct->DiagScale == BOTH);
      colequ = (ScalePermstruct->DiagScale == COL) ||
               (ScalePermstruct->DiagScale == BOTH);
    }
    else
    {
      rowequ = colequ = FALSE;
    }

    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    etree = LUstruct->etree;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;
    Glu_persist = LUstruct->Glu_persist;
    if (Equil)
    {
      /* Allocate storage if not done so before. */
      switch (ScalePermstruct->DiagScale)
      {
      case NOEQUIL:
        if (!(R = (double *) doubleMalloc_dist(m)))
          ABORT("Malloc fails for R[].");
        if (!(C = (double *) doubleMalloc_dist(n)))
          ABORT("Malloc fails for C[].");
        ScalePermstruct->R = R;
        ScalePermstruct->C = C;
        break;
      case ROW:
        if (!(C = (double *) doubleMalloc_dist(n)))
          ABORT("Malloc fails for C[].");
        ScalePermstruct->C = C;
        break;
      case COL:
        if (!(R = (double *) doubleMalloc_dist(m)))
          ABORT("Malloc fails for R[].");
        ScalePermstruct->R = R;
        break;
      default:
       printf("diagscale: %i %i %i %i\n",ScalePermstruct->DiagScale,NOEQUIL,ROW,COL);
       ABORT("Never get here.");
       break;
      }
    }

    /* ------------------------------------------------------------
       Diagonal scaling to equilibrate the matrix.
       ------------------------------------------------------------*/
    if (Equil)
    {
      t = SuperLU_timer_();

      if (Fact == SamePattern_SameRowPerm)
      {
        /* Reuse R and C. */
        switch (ScalePermstruct->DiagScale)
        {
        case NOEQUIL:
          break;
        case ROW:
          for (j = 0; j < n; ++j)
          {
            for (i = colptr[j]; i < colptr[j+1]; ++i)
            {
              irow = rowind[i];
              a[i] *= R[irow];       /* Scale rows. */
            }
          }
          break;
        case COL:
          for (j = 0; j < n; ++j)
          {
            for (i = colptr[j]; i < colptr[j+1]; ++i)
            {
              a[i] *= C[j];          /* Scale columns. */
            }
          }
          break;
        case BOTH:
          for (j = 0; j < n; ++j)
          {
            for (i = colptr[j]; i < colptr[j+1]; ++i)
            {
              irow = rowind[i];
              a[i] *= R[irow] * C[j]; /* Scale rows and columns. */
            }
          }
          break;
        }
      }
      else
      {
        if (!iam)
        {
          /* Compute row and column scalings to equilibrate matrix A. */
          dgsequ_dist(A, R, C, &rowcnd, &colcnd, &amax, &iinfo);

          MPI_Bcast(&iinfo, 1, mpi_int_t, 0, grid->comm);
          if (iinfo == 0)
          {
            MPI_Bcast(R,       m, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(C,       n, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&rowcnd, 1, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&colcnd, 1, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&amax,   1, MPI_DOUBLE, 0, grid->comm);
          }
          else
          {
            if (iinfo > 0)
            {
              if (iinfo <= m)
              {
                fprintf(stderr, "The %d-th row of A is exactly zero\n",
                        iinfo);
              }
              else
              {
                fprintf(stderr, "The %d-th column of A is exactly zero\n",
                        iinfo-n);
              }
              exit(-1);
            }
          }
        }
        else
        {
          MPI_Bcast(&iinfo, 1, mpi_int_t, 0, grid->comm);
          if (iinfo == 0)
          {
            MPI_Bcast(R,       m, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(C,       n, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&rowcnd, 1, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&colcnd, 1, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(&amax,   1, MPI_DOUBLE, 0, grid->comm);
          }
          else
          {
            ABORT("DGSEQU failed\n");
          }
        }

        /* Equilibrate matrix A. */
        dlaqgs_dist(A, R, C, rowcnd, colcnd, amax, equed);
        if (lsame_(equed, "R"))
        {
          ScalePermstruct->DiagScale = rowequ = ROW;
        }
        else if (lsame_(equed, "C"))
        {
          ScalePermstruct->DiagScale = colequ = COL;
        }
        else if (lsame_(equed, "B"))
        {
          ScalePermstruct->DiagScale = BOTH;
          rowequ = ROW;
          colequ = COL;
        }
        else
        {
          ScalePermstruct->DiagScale = NOEQUIL;
        }
      } /* if Fact ... */

      stat.utime[EQUIL] = SuperLU_timer_() - t;
    } /* if Equil ... */


    /* ------------------------------------------------------------
       Permute rows of A.
       ------------------------------------------------------------*/
    if ((int) options->RowPerm != (int) NO)
    {
      t = SuperLU_timer_();

      if (Fact == SamePattern_SameRowPerm  /* Reuse perm_r. */
          || options->RowPerm == MY_PERMR)
      {
        /* Use my perm_r. */
        /*     for (j = 0; j < n; ++j) {
                     for (i = colptr[j]; i < colptr[j+1]; ++i) {*/
        for (i = 0; i < colptr[n]; ++i)
        {
          irow = rowind[i];
          rowind[i] = perm_r[irow];
          /*    }*/
        }
      }
      else if (!factored)
      {
        if (job == 5)
        {
          /* Allocate storage for scaling factors. */
          if (!(R1 = (double *) SUPERLU_MALLOC(m * sizeof(double))))
            ABORT("SUPERLU_MALLOC fails for R1[]");
          if (!(C1 = (double *) SUPERLU_MALLOC(n * sizeof(double))))
            ABORT("SUPERLU_MALLOC fails for C1[]");
        }

        if (!iam)
        {
          /* Process 0 finds a row permutation for large diagonal. */
          dldperm(job, m, nnz, colptr, rowind, a, perm_r, R1, C1);

          MPI_Bcast(perm_r, m, mpi_int_t, 0, grid->comm);
          if (job == 5 && Equil)
          {
            MPI_Bcast(R1, m, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(C1, n, MPI_DOUBLE, 0, grid->comm);
          }
        }
        else
        {
          MPI_Bcast(perm_r, m, mpi_int_t, 0, grid->comm);
          if (job == 5 && Equil)
          {
            MPI_Bcast(R1, m, MPI_DOUBLE, 0, grid->comm);
            MPI_Bcast(C1, n, MPI_DOUBLE, 0, grid->comm);
          }
        }

        if (job == 5)
        {
          if (Equil)
          {
            for (i = 0; i < n; ++i)
            {
              R1[i] = exp(R1[i]);
              C1[i] = exp(C1[i]);
            }
            for (j = 0; j < n; ++j)
            {
              for (i = colptr[j]; i < colptr[j+1]; ++i)
              {
                irow = rowind[i];
                a[i] *= R1[irow] * C1[j]; /* Scale the matrix. */
                rowind[i] = perm_r[irow];
              }
            }

            /* Multiply together the scaling factors. */
            if (rowequ)
            {
              for (i = 0; i < m; ++i)
              {
                R[i] *= R1[i];
              }
            }
            else
            {
              for (i = 0; i < m; ++i)
              {
                R[i] = R1[i];
              }
            }
            if (colequ)
            {
              for (i = 0; i < n; ++i)
              {
                C[i] *= C1[i];
              }
            }
            else
            {
              for (i = 0; i < n; ++i)
              {
                C[i] = C1[i];
              }
            }

            ScalePermstruct->DiagScale = BOTH;
            rowequ = colequ = 1;
          }
          else
          {
            /* No equilibration. */
            /*        for (j = 0; j < n; ++j) {
                                for (i = colptr[j]; i < colptr[j+1]; ++i) {*/
            for (i = colptr[0]; i < colptr[n]; ++i)
            {
              irow = rowind[i];
              rowind[i] = perm_r[irow];
            }
            /*        }*/
          }
          SUPERLU_FREE(R1);
          SUPERLU_FREE(C1);
        }
        else
        {
          /* job = 2,3,4 */
          for (j = 0; j < n; ++j)
          {
            for (i = colptr[j]; i < colptr[j+1]; ++i)
            {
              irow = rowind[i];
              rowind[i] = perm_r[irow];
            }
          }
        }
      } /* else !factored */

      t = SuperLU_timer_() - t;
      stat.utime[ROWPERM] = t;
    } /* if options->RowPerm ... */

    if (!factored || options->IterRefine)
    {
      /* Compute norm(A), which will be used to adjust small diagonal. */
      if (notran)
      {
        *(unsigned char *)norm = '1';
      }
      else
      {
        *(unsigned char *)norm = 'I';
      }
      anorm = dlangs_dist(norm, A);
    }



    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    if (!factored)
    {
      t = SuperLU_timer_();
      /*
         Get column permutation vector perm_c[], according to permc_spec:
           permc_spec = NATURAL:  natural ordering
           permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
           permc_spec = MMD_ATA:  minimum degree on structure of A'*A
           permc_spec = COLAMD:   approximate minimum degree column ordering
           permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
      */
      permc_spec = options->ColPerm;
      if (permc_spec != MY_PERMC && Fact == DOFACT)
      {
        /* Use an ordering provided by SuperLU */
        get_perm_c_dist(iam, permc_spec, A, perm_c);
      }

      /* Compute the elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'
         (a.k.a. column etree), depending on the choice of ColPerm.
         Adjust perm_c[] to be consistent with a postorder of etree.
         Permute columns of A to form A*Pc'. */
      sp_colorder(options, A, perm_c, etree, AC);

      /* Form Pc*A*Pc' to preserve the diagonal of the matrix Pr*A. */
      ACstore = AC->Store;
      for (j = 0; j < n; ++j)
      {
        for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i)
        {
          irow = ACstore->rowind[i];
          ACstore->rowind[i] = perm_c[irow];
        }
      }
      stat.utime[COLPERM] = SuperLU_timer_() - t;

      /* Perform a symbolic factorization on matrix A and set up the
         nonzero data structures which are suitable for supernodal GENP. */
      if (Fact != SamePattern_SameRowPerm)
      {
        t = SuperLU_timer_();
        if (!(Glu_freeable = (Glu_freeable_t *)
                             SUPERLU_MALLOC(sizeof(Glu_freeable_t))))
          ABORT("Malloc fails for Glu_freeable.");

        iinfo = symbfact(options, iam, AC, perm_c, etree,
                         Glu_persist, Glu_freeable);

        stat.utime[SYMBFAC] = SuperLU_timer_() - t;

        if (iinfo < 0)
        {
          QuerySpace_dist(n, -iinfo, Glu_freeable,
                          &symbolic_memory_statistics_storage.Memory_usage);
        }
        else
        {
          if (!iam)
          {
            fprintf(stderr, "symbfact() error returns %d\n", iinfo);
            exit(-1);
          }
        }
      }

      /* Distribute the L and U factors onto the process grid. */
      t = SuperLU_timer_();
      /* dist_mem_use = */
      ddistribute(Fact, n, AC, Glu_freeable, LUstruct, grid);
      stat.utime[DIST] = SuperLU_timer_() - t;

      /* Deallocate storage used in symbolic factor. */
      if (Fact != SamePattern_SameRowPerm)
      {
        iinfo = symbfact_SubFree(Glu_freeable);
        SUPERLU_FREE(Glu_freeable);
      }

      /* Perform numerical factorization in parallel. */
      t = SuperLU_timer_();
      pdgstrf(options, m, n, anorm, LUstruct, grid, &stat, info);
      stat.utime[FACT] = SuperLU_timer_() - t;
    }
    else if (options->IterRefine)
    {
      /* options->Fact==FACTORED */
      /* Permute columns of A to form A*Pc' using the existing perm_c.
         NOTE: rows of A were previously permuted to Pc*A.
      */
      sp_colorder(options, A, perm_c, NULL, AC);
    } /* if !factored ... */

    if (*info!=0)
    {
      printf("Trouble in  pdgstrf. Info=%i\n",*info);
      if (*info>0)
      {
        printf("U(%i,%i) is exactly zero. The factorization has\n",*info,*info);
        printf("been completed, but the factor U is exactly singular,\n");
        printf("and division by zero will occur if it is used to solve a\n");
        printf("system of equations.\n");
      }
      else
      {
        printf("The %i-th argument had an illegal value.\n", *info);
      }
    }

    /* Print the statistics. */
    if ((doc==0) && (!iam))
    {
      printf("\nstats after setup....\n");
      PStatPrint(options, &stat, grid);
    }

    /* ------------------------------------------------------------
       Set up data structure.
       ------------------------------------------------------------*/
    superlu_data->A = A;
    superlu_data->AC = AC;
    superlu_data->options = options;
    superlu_data->ScalePermstruct = ScalePermstruct;
    superlu_data->LUstruct = LUstruct;
    superlu_data->colequ = colequ;
    superlu_data->rowequ = rowequ;
    superlu_data->anorm = anorm;
    *data = superlu_data;

  } /* End of setup */


  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM A SOLVE
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==2)
  {
    /* Get pointer to the grid */
    superlu_data = (superlu_dist_data *)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol)
    {
      return;
    }

    if ((doc==0)&&(!iam))
    {
      printf("\nPerforming SuperLU_DIST solve\n");
    }

    /* ------------------------------------------------------------
       Set other  pointers to data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    AC = superlu_data->AC;
    options = superlu_data->options;
    ScalePermstruct = superlu_data->ScalePermstruct;
    LUstruct = superlu_data->LUstruct;
    colequ = superlu_data->colequ;
    rowequ = superlu_data->rowequ;
    anorm = superlu_data->anorm;

    /* Initialization. */
    Astore = A->Store;
    colptr = Astore->colptr;
    rowind = Astore->rowind;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;

    /* Local control paramaters */
    Fact = options->Fact;
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);
    notran = (options->Trans == NOTRANS);


    /* ------------------------------------------------------------
       Compute the solution matrix X.
       ------------------------------------------------------------*/
    if (!(b_work = doubleMalloc_dist(n)))
    {
      ABORT("Malloc fails for b_work[]");
    }

    /* ------------------------------------------------------------
       Scale the right-hand side if equilibration was performed.
       ------------------------------------------------------------*/
    if (notran)
    {
      if (rowequ)
      {
        b_col = B;
        for (j = 0; j < nrhs; ++j)
        {
          for (i = 0; i < m; ++i)
          {
            b_col[i] *= R[i];
          }
          b_col += ldb;
        }
      }
    }
    else if (colequ)
    {
      b_col = B;
      for (j = 0; j < nrhs; ++j)
      {
        for (i = 0; i < m; ++i)
        {
          b_col[i] *= C[i];
        }
        b_col += ldb;
      }
    }

    /* ------------------------------------------------------------
       Permute the right-hand side to form Pr*B.
       ------------------------------------------------------------*/
    if ((int) options->RowPerm != (int) NO)
    {
      if (notran)
      {
        b_col = B;
        for (j = 0; j < nrhs; ++j)
        {
          for (i = 0; i < m; ++i)
          {
            b_work[perm_r[i]] = b_col[i];
          }
          for (i = 0; i < m; ++i)
          {
            b_col[i] = b_work[i];
          }
          b_col += ldb;
        }
      }
    }


    /* ------------------------------------------------------------
       Permute the right-hand side to form Pc*B.
       ------------------------------------------------------------*/
    if (notran)
    {
      b_col = B;
      for (j = 0; j < nrhs; ++j)
      {
        for (i = 0; i < m; ++i)
        {
          b_work[perm_c[i]] = b_col[i];
        }
        for (i = 0; i < m; ++i)
        {
          b_col[i] = b_work[i];
        }
        b_col += ldb;
      }
    }


    /* Save a copy of the right-hand side. */
    ldx = ldb;
    if (!(X = doubleMalloc_dist(((size_t)ldx) * nrhs)))
    {
      ABORT("Malloc fails for X[]");
    }

    x_col = X;
    b_col = B;
    for (j = 0; j < nrhs; ++j)
    {
      for (i = 0; i < ldb; ++i)
      {
        x_col[i] = b_col[i];
      }
      x_col += ldx;
      b_col += ldb;
    }

    /* ------------------------------------------------------------
       Solve the linear system.
       ------------------------------------------------------------*/
    pdgstrs_Bglobal(n, LUstruct, grid, X, ldb, nrhs, &stat, info);
    if (*info!=0)
    {
      printf("Trouble in pdgstrs_Bglobal. Info=%i\n",*info);
      printf("The %i-th argument had an illegal value.\n", *info);
    }

    /* ------------------------------------------------------------
       Use iterative refinement to improve the computed solution and
       compute error bounds and backward error estimates for it.
       ------------------------------------------------------------*/
    if (options->IterRefine)
    {
      /* Improve the solution by iterative refinement. */
      t = SuperLU_timer_();

      /*  Storage for backward error */
      if (!(berr = doubleMalloc_dist(nrhs)))
      {
        ABORT("Malloc fails for berr[].");
      }

      pdgsrfs_ABXglobal(n, AC, anorm, LUstruct, grid, B, ldb,
                        X, ldx, nrhs, berr, &stat, info);
      if (*info!=0)
      {
        printf("Trouble in pdgsrfs_ABXglobal. Info=%i\n",*info);
        printf("The %i-th argument had an illegal value.\n", *info);
      }
      stat.utime[REFINE] = SuperLU_timer_() - t;
    }

    /* Print the statistics. */
    if ((doc==0) && (!iam))
    {
      printf("\nstats after solve....\n");
      PStatPrint(options, &stat, grid);
    }

    /* Permute the solution matrix X <= Pc'*X. */
    for (j = 0; j < nrhs; j++)
    {
      b_col = &B[j*ldb];
      x_col = &X[j*ldx];
      for (i = 0; i < n; ++i)
      {
        b_col[i] = x_col[perm_c[i]];
      }
    }

    /* Transform the solution matrix X to a solution of the original system
       before the equilibration. */
    if (notran)
    {
      if (colequ)
      {
        b_col = B;
        for (j = 0; j < nrhs; ++j)
        {
          for (i = 0; i < n; ++i)
          {
            b_col[i] *= C[i];
          }
          b_col += ldb;
        }
      }
    }
    else if (rowequ)
    {
      b_col = B;
      for (j = 0; j < nrhs; ++j)
      {
        for (i = 0; i < n; ++i)
        {
          b_col[i] *= R[i];
        }
        b_col += ldb;
      }
    }


    /* Clean up memory */
    if (options->IterRefine)
    {
      SUPERLU_FREE(berr);
    }
    SUPERLU_FREE(b_work);
    SUPERLU_FREE(X);

  } /* End of solve */

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     PERFORM CLEAN UP OF MEMORY
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (opt_flag==3)
  {
    /* Get pointer to the process grid */
    superlu_data = (superlu_dist_data *)*data;
    grid = superlu_data->grid;

    /* Bail out if I do not belong in the grid. */
    iam = grid->iam;
    if (iam >= nprow * npcol) goto out;
    if ((doc==0)&&(!iam))
    {
      printf("\nCleaning up memory allocated for SuperLU_DIST\n");
    }

    /* ------------------------------------------------------------
       Set pointers to the data structure.
       ------------------------------------------------------------*/
    A = superlu_data->A;
    AC = superlu_data->AC;
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
    rowequ = colequ = FALSE;
    if (factored || (Fact == SamePattern_SameRowPerm && Equil))
    {
      rowequ = (ScalePermstruct->DiagScale == ROW) ||
               (ScalePermstruct->DiagScale == BOTH);
      colequ = (ScalePermstruct->DiagScale == COL) ||
               (ScalePermstruct->DiagScale == BOTH);
    }
    else
    {
      rowequ = colequ = FALSE;
    }

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
       /* Apparently this one is ok */
       /* printf("diagscale: %i %i %i %i\n",ScalePermstruct->DiagScale,NOEQUIL,ROW,COL); */
       /* ABORT("Never get here."); */
       break;
      }
    }
    if (!factored || (factored && options->IterRefine))
    {
      Destroy_CompCol_Permuted_dist(AC);
    }

    /*  Free storage */
    ScalePermstructFree(ScalePermstruct);
    Destroy_LU(n, grid, LUstruct);
    LUstructFree(LUstruct);
    //Destroy_CompRowLoc_Matrix_dist(&A);
    // Only destroy the store part of the matrix
    Destroy_SuperMatrix_Store_dist(A);

    /* Deallocate memory */
    SUPERLU_FREE(A);
    SUPERLU_FREE(AC);
    SUPERLU_FREE(ScalePermstruct);
    SUPERLU_FREE(LUstruct);
    SUPERLU_FREE(options);

    /*  Release the superlu process grid. */
out:
    superlu_gridexit(grid);

    SUPERLU_FREE(grid);
    SUPERLU_FREE(superlu_data);
  }

  /*  Free storage */
  PStatFree(&stat);

  return;
}




/*-- OBSOLETE -- CAN BE REMOVED FROM SVN --*/

/*----------------------------------------------------------------
 * Interface to distributed SuperLU written by MH, based on
 * pddrive.c demo code in superlu distribution.
 *
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 *----------------------------------------------------------------
 */
#include <math.h>
#include "oomph_superlu_dist_2.0.h"




/*----------------------------------------------------------------
 * Bridge to distributed SuperLU (version 2.0).
 * Requires input of system matrix in compressed column form.
 *
 * Parameters:
 * - n = size of system (square matrix)
 * - nnz = # of nonzero entries
 * - values = 1D C array of nonzero entries
 * - row_index =  1D C array of row indices
 * - column_start =  1D C array of column start indices
 * - b = 1D C array representing the rhs vector, is overwritten
 *       by solution.
 * - nprow = # of rows in process grid
 * - npcol = # of columns in process grid
 * - doc = 0/1 for doc/no doc
 *
 * Return value:
 *         = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *----------------------------------------------------------------
 */
int superlu_dist_global_matrix_bridge(int n, int nnz, 
                                      double *values, int *row_index,
                                      int *col_start, double *b, 
                                      int nprow, int npcol, int doc)
{

/*  Some SuperLU structures */
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    SOLVEstruct_t SOLVEstruct;
    gridinfo_t grid;

/*  We're only doing single rhs problems */
    int nrhs=1;
    
/*  Square matrix */
    int m=n;
    
/*  Initialize the superlu process grid.  */
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
    
/*  Bail out if I do not belong in the grid. */
    int iam = grid.iam;
    if ( iam >= nprow * npcol )	goto out;
    if ((!iam)&&(doc==0)) 
     {
      printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
     }

/*  Create SuperMatrix from compressed column representation */
    dCreate_CompCol_Matrix_dist(&A,m,n,nnz,values,row_index,col_start,
                                SLU_NC,SLU_D,SLU_GE);

/*  Storage for backward error */
    double  *berr;
    if ( !(berr = doubleMalloc_dist(nrhs)) )
     ABORT("Malloc fails for berr[].");


/*  Set the default input options:  */
    set_default_options_dist(&options); 

/*   Is the matrix transposed (NOTRANS or TRANS)? */
    options.Trans = NOTRANS; 
/*     options.Trans = TRANS; */

/*    Row permutations (NATURAL [= do nothing],  */
/*                      LargeDiag [default], ...)?  */
/*     options.RowPerm=NATURAL; */
    options.RowPerm=LargeDiag;


/*     Column permutations (NATURAL [= do nothing],  */
/*                          MMD_AT_PLUS_A [default],...) */
/*     options.ColPerm=NATURAL; */
    options.ColPerm=MMD_AT_PLUS_A;


/*  Iterative refinement (essential as far as I can tell).*/
/*  Can be "NO" or "DOUBLE"*/
    options.IterRefine = DOUBLE;


/*  Doc options on process 0 if required: */
    if ((!iam) && (doc==0)) print_options_dist(&options);

/*  Set matrix size */
    m = A.nrow;
    n = A.ncol;

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);


/*  Set 'Leading dimension' of rhs vector */
    int ldb=n;

/*  Info flag */
    int info;



/*  Print stats during solve? */
    if (doc==0)
     {
      options.PrintStat = YES;
     }
    else
     {
      options.PrintStat = NO;
     }

    /* Call the linear equation solver. */
/*     pdgssvx(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid, */
/* 	    &LUstruct, &SOLVEstruct, berr, &stat, &info); */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
                     &LUstruct, berr, &stat, &info);

    if (info!=0)
     {
      printf("Trouble in  pdgssvx_ABglobal. Info=%i",info);
      if (info<=A.ncol)
       {
        printf("U(%i,%i) is exactly zero. The factorization has",info,info);
        printf("been completed, but the factor U is exactly singular,");
        printf("so the solution could not be computed.");
       }
      else
       {
        printf("Memory allocation failure occurred");
       }
     }

    /* Print the statistics. */
    if (doc==0) PStatPrint(&options, &stat, &grid);


/*  Free storage */
    PStatFree(&stat);
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(n, &grid, &LUstruct);
    LUstructFree(&LUstruct);
    if ( options.SolveInitialized ) {
        dSolveFinalize(&options, &SOLVEstruct);
    }
    SUPERLU_FREE(berr);

    //Only destroy the store part of the matrix
    //Destroy_CompCol_Matrix_dist(&A);
    Destroy_SuperMatrix_Store_dist(&A);


/*  Release the superlu process grid. */
out:
    superlu_gridexit(&grid);

/*
 * info    (output) int*
 *         = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 */
    return info;

}

/*----------------------------------------------------------------
 * Bridge to distributed SuperLU with distributed memory (version 2.0).
 * Requires input of system matrix in compressed column form.
 *
 * Parameters:
 * - n = size of system (square matrix)
 * - nnz = # of nonzero entries
 * - values = 1D C array of nonzero entries
 * - row_index =  1D C array of row indices
 * - column_start =  1D C array of column start indices
 * - b = 1D C array representing the rhs vector, is overwritten
 *       by solution.
 * - nprow = # of rows in process grid
 * - npcol = # of columns in process grid
 * - doc = 0/1 for doc/no doc
 *
 * Return value:
 *         = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *----------------------------------------------------------------
 */
int superlu_dist_distributed_matrix_bridge(int n, int nnz, int nrow,
                                           int first_row, double *values, 
                                           int *col_index, int *row_start,
                                           double *b, int nprow, int npcol,
                                           int doc)
{
/*  Some SuperLU structures */
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    SOLVEstruct_t SOLVEstruct;
    gridinfo_t grid;

/*  We're only doing single rhs problems */
    int nrhs=1;

/*  Square matrix  */
    int m=n;

/*  Initialize the superlu process grid.  */
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

/*  Bail out if I do not belong in the grid. */
    int iam = grid.iam;
    if ( iam >= nprow * npcol )	goto out;
    if ((!iam)&&(doc==0))
     {
      printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
     }


/*  Create SuperMatrix from compressed row representation */
    dCreate_CompRowLoc_Matrix_dist(&A,m,n,nnz,nrow,first_row,
                                   values,col_index,row_start,
                                   SLU_NR_loc,SLU_D,SLU_GE);

/*  Storage for backward error */
    double  *berr;
    if ( !(berr = doubleMalloc_dist(nrhs)) )
     ABORT("Malloc fails for berr[].");

/*  Set the default input options:  */
    set_default_options_dist(&options);

/*   Is the matrix transposed (NOTRANS or TRANS)? */
    options.Trans = NOTRANS;
/*     options.Trans = TRANS; */

/*    Row permutations (NATURAL [= do nothing],  */
/*                      LargeDiag [default], ...)?  */
/*     options.RowPerm=NATURAL; */
    options.RowPerm=LargeDiag;


/*     Column permutations (NATURAL [= do nothing],  */
/*                          MMD_AT_PLUS_A [default],...) */
/*     options.ColPerm=NATURAL; */
    options.ColPerm=MMD_AT_PLUS_A;


/*  Iterative refinement (essential as far as I can tell).*/
/*  Can be "NO" or "DOUBLE"*/
    options.IterRefine = DOUBLE;


/*  Doc options on process 0 if required: */
    if ((!iam) && (doc==0)) print_options_dist(&options);


/*  Set matrix size */
    m = A.nrow;
    n = A.ncol;


/* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

/* Initialize the statistics variables. */
    PStatInit(&stat);


/*  Set 'Leading dimension' of rhs vector */
    int ldb=n;

/*  Info flag */
    int info;



/*  Print stats during solve? */
    if (doc==0)
     {
      options.PrintStat = YES;
     }
    else
     {
      options.PrintStat = NO;
     }

/* Call the linear equation solver. */
    pdgssvx(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
 	    &LUstruct, &SOLVEstruct, berr, &stat, &info);

    if (info!=0)
     {
      printf("Trouble in  pdgssvx. Info=%i",info);
      if (info<=A.ncol)
       {
        printf("U(%i,%i) is exactly zero. The factorization has",info,info);
        printf("been completed, but the factor U is exactly singular,");
        printf("so the solution could not be computed.");
       }
      else
       {
        printf("Memory allocation failure occurred");
       }
     }

    /* Print the statistics. */
    if (doc==0) PStatPrint(&options, &stat, &grid);


/*  Free storage */
    PStatFree(&stat);
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(n, &grid, &LUstruct);
    LUstructFree(&LUstruct);
    if ( options.SolveInitialized ) {
        dSolveFinalize(&options, &SOLVEstruct);
    }
    SUPERLU_FREE(berr);


/*     hierher */
    //Only destroy the store part of the matrix
    //Destroy_CompRowLoc_Matrix_dist(&A);
    
    Destroy_SuperMatrix_Store_dist(&A);



/*  Release the superlu process grid. */
out:
    superlu_gridexit(&grid);

/*
 * info    (output) int*
 *         = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 */
    return info;

}




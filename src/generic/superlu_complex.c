
/* Wrapper for SuperLU solver. Based on fortran wrapper supplied
 * with SuperLU version 3.0. Given that, it seems appropriate
 * to retain this copyright notice:
 *
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#include "../../external_src/oomph_superlu_4.3/slu_zdefs.h"
#include "math.h"


/* ================================================= */
/* Pointer to the LU factors*/
/* ================================================= */
typedef void* fptr;  


/* ================================================= */
/* Struct for the lu factors  */
/* ================================================= */
typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;



/* =========================================================================
 * Wrapper to superlu solver:
 *
 * op_flag    = int specifies the operation:
 *                1, performs LU decomposition for the first time
 *                2, performs triangular solve
 *                3, free all the storage in the end
 * n          = dimension of matrix
 * nnz        = # of nonzero entries
 * nrhs       = # of RHSs
 * values     = double array containing the nonzero entries 
 * rowind     = int array containing the row indices of the entries 
 * colptr     = int array containing the column start
 * b          = double array containing the rhs vector (gets overwritten
 *              with solution)
 * ldb        = leading dimension of b
 * transpose  =  0/1 if matrix is transposed/not transposed
 * doc        = 0/1 for full doc/no full doc 
 * info       = info flag from superlu
 * f_factors  = pointer to LU factors. (If op_flag == 1, it is an output 
 *              and contains the pointer pointing to the structure of 
 *              the factored matrices. Otherwise, it it an input.
 * Returns the SIGN of the determinant of the matrix 
 * =========================================================================
 */
int superlu_complex(int *op_flag, int *n, int *nnz, int *nrhs,
                    doublecomplex *values, int *rowind, int *colptr,
                    doublecomplex *b, int *ldb, int *transpose, int *doc,
                    fptr *f_factors, int *info)
 
{
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, j, panel_size, permc_spec, relax;
    trans_t  trans;
    //double   drop_tol = 0.0; //No longer needed SuperLU 4.3
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;

    doublecomplex *Lval;
    doublecomplex *diagU, *dblock;
    int_t fsupc, nsupr, nsupc, luptr;
    int_t i2, k2, nsupers; 
    int signature=1;
    int sign = 0;

/*   Do we need to transpose? */
    if (*transpose==0)
     {
      trans = NOTRANS;
     }
    else
     {
      trans = TRANS;
     }

    if ( *op_flag == 1 ) { /* LU decomposition */

        /* Set the default input options. */
        set_default_options(&options);

	/* Initialize the statistics variables. */
	StatInit(&stat);

	zCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
			       SLU_NC, SLU_Z, SLU_GE);
	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	permc_spec = options.ColPerm;
	get_perm_c(permc_spec, &A, perm_c);
	
	sp_preorder(&options, &A, perm_c, etree, &AC);

	panel_size = sp_ienv(1);
	relax = sp_ienv(2);

	zgstrf(&options, &AC, /*drop_tol,*/ relax, panel_size, 
	       etree, NULL, 0, perm_c, perm_r, L, U, &stat, info);

	if ( *info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
            if (*doc!=0)
             {
              printf(" No of nonzeros in factor L = %d\n", Lstore->nnz);
              printf(" No of nonzeros in factor U = %d\n", Ustore->nnz);
              printf(" No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
              zQuerySpace(L, U, &mem_usage);
              printf(" L\\U MB %.3f\ttotal MB needed %.3f\n",
                     //\texpansions %d\n",
                     mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
                     //mem_usage.expansions);
             }
	} else {
	    printf("dgstrf() error returns INFO= %d\n", *info);
	    if ( *info <= *n ) { /* factorization completes */
	       zQuerySpace(L, U, &mem_usage);
	       printf(" L\\U MB %.3f\ttotal MB needed %.3f\n",
                      //\texpansions %d\n\n",
                      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
               //mem_usage.expansions);
	    }
	}
	
	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_c = perm_c;
	LUfactors->perm_r = perm_r;
	*f_factors = (fptr) LUfactors;

        //Work out and print the sign of the determinant
        //This code is hacked from supraLU by  Alex Pletzer
        //and Doug McCune (NTCC) (http://w3.pppl.gov/ntcc/SupraLu/)
        Lstore = L->Store;
        Lval = Lstore->nzval;
        nsupers = Lstore->nsuper + 1;

        //Get the diagonal entries of the U matrix
        //Allocate store for the entries
        if ( !(diagU = SUPERLU_MALLOC( *n * sizeof(SuperMatrix) )) )
         ABORT("Malloc fails for diagU[].");
        //Loop over the number of super diagonal terms(?)
        for(k2=0; k2< nsupers; k2++)
         {
          fsupc = L_FST_SUPC(k2);
          nsupc = L_FST_SUPC(k2+1) - fsupc;
          nsupr = L_SUB_START(fsupc+1) - L_SUB_START(fsupc);
          luptr = L_NZ_START(fsupc);

          dblock = &diagU[fsupc];
          for(i2 = 0; i2 < nsupc; i2++)
           {
            dblock[i2] = Lval[luptr];
            luptr += nsupr + 1;
           }
         }
                  
        //Now multiply all the diagonal terms together to get the determinant
        //Note that we need to use the mantissa, exponent formulation to
        //avoid underflow errors
        double determinant_mantissa=1.0;
        int determinant_exponent = 0, iexp;
        for(i=0; i<*n; i++)
         {
          determinant_mantissa *= frexp(diagU[i].r, &iexp);
          determinant_exponent += iexp;
          /* normalise*/
          determinant_mantissa = frexp(determinant_mantissa,&iexp);
          determinant_exponent += iexp;
                    
          /*Now worry about the permutations 
            (this is done in a stupid, but not too inefficient way)*/
          for(j=i;j<*n;j++)
           {
            if(perm_r[j] < perm_r[i]) {signature *= -1;}
            if(perm_c[j] < perm_c[i]) {signature *= -1;}
           }
         }
        
        //Find the sign of the determinant
        if(determinant_mantissa > 0.0) {sign = 1;}
        if(determinant_mantissa < 0.0) {sign = -1;}
        
        //Multiply the sign by the signature
        sign *= signature;
 
	/* Free un-wanted storage */
        SUPERLU_FREE(diagU);
	SUPERLU_FREE(etree);
	Destroy_SuperMatrix_Store(&A);
	Destroy_CompCol_Permuted(&AC);
	StatFree(&stat);

        //Return the sign of the determinant
        return sign;
        
    } else if ( *op_flag == 2 ) { /* Triangular solve */
	/* Initialize the statistics variables. */
	StatInit(&stat);

	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	L = LUfactors->L;
	U = LUfactors->U;
	perm_c = LUfactors->perm_c;
	perm_r = LUfactors->perm_r;

	zCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_Z, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        zgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

	Destroy_SuperMatrix_Store(&B);
	StatFree(&stat);

        //Return zero
        return 0;

    } else if ( *op_flag == 3 ) { /* Free storage */
	/* Free the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
	Destroy_SuperNode_Matrix(LUfactors->L);
	Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
        return 0;
    } else {
	fprintf(stderr,"Invalid op_flag=%d passed to c_cpp_dgssv()\n",*op_flag);
	exit(-1);
    }
}



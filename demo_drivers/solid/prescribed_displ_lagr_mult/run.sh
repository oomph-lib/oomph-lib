#! /bin/sh


# Check existence of working dir
dir="Runs"
if [ -e `echo $dir`   ] 
then
    echo " "
    echo "ERROR: Please delete directory $dir and try again"
    echo " "
    exit
fi
mkdir $dir


#Build it
make prescribed_displ_lagr_mult_precond

mkdir $dir/CaseSuperLU
mkdir $dir/CaseElasticBlockSuperLU
mkdir $dir/CaseElasticBlockHypreCGLagrange

nel_1d_list="5 10" # "20 40 80"
for nel_1d in `echo $nel_1d_list`; do

    # SuperLU for the lot
    mkdir RESLT
    ./prescribed_displ_lagr_mult_precond --no_adapt --nel_1d $nel_1d > RESLT/OUTPUT
    mv RESLT $dir/CaseSuperLU/RESLT_`echo $nel_1d`
    
    # Block diagonal for elastic block; SuperLU for all solves
    mkdir RESLT
    ./prescribed_displ_lagr_mult_precond --no_adapt --nel_1d $nel_1d --block_upper_for_elastic_block > RESLT/OUTPUT
    mv RESLT $dir/CaseElasticBlockSuperLU/RESLT_`echo $nel_1d`

    
    
    # Block diagonal for elastic block; Hypre for elastic sub-blocks; 
    # Trilinos CG for Lagrange multiplier blocks
    mkdir RESLT
    ./prescribed_displ_lagr_mult_precond --no_adapt --nel_1d $nel_1d --block_upper_for_elastic_block --hypre_for_elastic_blocks --trilinos_cg_for_lagrange_multiplier_blocks > RESLT/OUTPUT
    mv RESLT $dir/CaseElasticBlockHypreCGLagrange/RESLT_`echo $nel_1d`

done
    
ndof_string="<TR> <TD> <CODE>n_dof</CODE> </TD> "
superlu_iter_string="<TR><TD>SuperLU  </TD>"
superlu_block_iter_string="<TR><TD>Upper triangular E</TD>"
full_iter_string="<TR><TD>Upper triangular E; Hypre/CG</TD>"
superlu_solve_string="<TR><TD>SuperLU  </TD>"
superlu_block_solve_string="<TR><TD>Upper triangular E</TD>"
full_solve_string="<TR><TD>Upper triangular E; Hypre/CG</TD>"
for nel_1d in `echo $nel_1d_list`; do

    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    ndofs=`grep "( ndof = " $file  | awk '{print $10}' `
    ndof=`echo $ndofs | awk '{print $1}' `
    ndof_string=`echo -n  $ndof_string " <TD> "$ndof " </TD> "` 


    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_iter_string=`echo $superlu_iter_string " <TD>  " $av " </TD> "`

    file=$dir/CaseElasticBlockSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_iter_string=`echo $superlu_block_iter_string " <TD>  " $av" </TD> "`


    file=$dir/CaseElasticBlockHypreCGLagrange/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_iter_string=`echo $full_iter_string " <TD>  " $av " </TD> "`



    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_solve_string=`echo $superlu_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CaseElasticBlockSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_solve_string=`echo $superlu_block_solve_string "  <TD>  " $av" </TD> "`


    file=$dir/CaseElasticBlockHypreCGLagrange/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_solve_string=`echo $full_solve_string " <TD>   " $av" </TD> "`



done





ndof_string="<TR> <TD> <CODE>n_dof</CODE> </TD> "
superlu_iter_string="<TR><TD>SuperLU  </TD>"
superlu_block_iter_string="<TR><TD>Upper triangular E</TD>"
full_iter_string="<TR><TD>Upper triangular E; Hypre/CG</TD>"
superlu_solve_string="<TR><TD>SuperLU  </TD>"
superlu_block_solve_string="<TR><TD>Upper triangular E</TD>"
full_solve_string="<TR><TD>Upper triangular E; Hypre/CG</TD>"
for nel_1d in `echo $nel_1d_list`; do

    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    ndofs=`grep "( ndof = " $file  | awk '{print $10}' `
    ndof=`echo $ndofs | awk '{print $1}' `
    ndof_string=`echo -n  $ndof_string " <TD> "$ndof " </TD> "` 


    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_iter_string=`echo $superlu_iter_string " <TD>  " $av " </TD> "`

    file=$dir/CaseElasticBlockSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_iter_string=`echo $superlu_block_iter_string " <TD>  " $av" </TD> "`


    file=$dir/CaseElasticBlockHypreCGLagrange/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $7}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_iter_string=`echo $full_iter_string " <TD>  " $av " </TD> "`



    file=$dir/CaseSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_solve_string=`echo $superlu_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CaseElasticBlockSuperLU/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_solve_string=`echo $superlu_block_solve_string "  <TD>  " $av" </TD> "`


    file=$dir/CaseElasticBlockHypreCGLagrange/RESLT_`echo $nel_1d`/OUTPUT
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $10}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_solve_string=`echo $full_solve_string " <TD>   " $av" </TD> "`



done

echo "<CENTER><TABLE>"
echo $ndof_string " </TR>" 
echo $superlu_iter_string" </TR>" 
echo $superlu_block_iter_string" </TR>" 
echo $full_iter_string" </TR>" 
echo "</TABLE></CENTER>"

echo " " 
echo "<CENTER><TABLE>"
echo $ndof_string" </TR>" 
echo $superlu_solve_string" </TR>" 
echo $superlu_block_solve_string" </TR>" 
echo $full_solve_string" </TR>" 
echo "</TABLE></CENTER>"




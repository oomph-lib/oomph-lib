#! /bin/sh

# Check existence of working dir
dir="Runs"

mesh_mult_list="3 6 9 12 " # 12 9 6 3" # "18 12 6 3" # "24"

skip_runs=1

if [ $skip_runs == 0 ]; then
 
if [ -e `echo $dir`   ] 
then
    echo " "
    echo "ERROR: Please delete directory $dir and try again"
    echo " "
    exit
fi
mkdir $dir


#Build it
make fsi_channel_with_leaflet_precond

mkdir $dir/CaseSuperLU
mkdir $dir/CasePrecondSuperLUForBlocks
mkdir $dir/CasePrecondSuperLUForBlocksLSC
mkdir $dir/CasePrecondLSC


for mesh_mult in `echo $mesh_mult_list`; do

    # SuperLU for the lot
    mkdir RESLT
    ./fsi_channel_with_leaflet_precond --use_direct_solver                                      --mesh_multiplier $mesh_mult > RESLT/OUTPUT
    mv RESLT $dir/CaseSuperLU/RESLT_`echo $mesh_mult`
    
    # Block preconditioner; SuperLU for all block solves; no lsc
    mkdir RESLT
    ./fsi_channel_with_leaflet_precond                      --suppress_lsc --superlu_for_blocks --mesh_multiplier $mesh_mult > RESLT/OUTPUT
    mv RESLT $dir/CasePrecondSuperLUForBlocks/RESLT_`echo $mesh_mult`

    # Block preconditioner; SuperLU for all block solves; lsc
    mkdir RESLT
    ./fsi_channel_with_leaflet_precond                                     --superlu_for_blocks --mesh_multiplier $mesh_mult > RESLT/OUTPUT
    mv RESLT $dir/CasePrecondSuperLUForBlocksLSC/RESLT_`echo $mesh_mult`

    # Block preconditioner; optimal block solves
    mkdir RESLT
    ./fsi_channel_with_leaflet_precond                                                          --mesh_multiplier $mesh_mult > RESLT/OUTPUT
    mv RESLT $dir/CasePrecondLSC/RESLT_`echo $mesh_mult`
    
done


fi


echo "" > junk.html
postfix_list="_STEADY _UNSTEADY"
for postfix in `echo $postfix_list`; do

echo $postfix >> junk.html

ndof_string="<TR> <TD> <CODE>n_dof</CODE> </TD> "
superlu_block_iter_string="<TR><TD>GMRES (blocks solved by SuperLU)</TD>"
lsc_superlu_block_iter_string="<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD>"
full_iter_string="<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD>"
superlu_solve_string="<TR><TD>SuperLU</TD>"
superlu_block_solve_string="<TR><TD>GMRES (blocks solved by SuperLU)</TD>"
lsc_superlu_block_solve_string="<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD>"
full_solve_string="<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD>"


for mesh_mult in `echo $mesh_mult_list`; do

    file=$dir/CaseSuperLU/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    ndofs=`grep "( ndof = " $file  | awk '{print $8}' `
    ndof=`echo $ndofs | awk '{print $1}' `
    ndof_string=`echo -n  $ndof_string " <TD> "$ndof " </TD> "` 

    file=$dir/CasePrecondSuperLUForBlocks/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_iter_string=`echo $superlu_block_iter_string " <TD>  " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocksLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    lsc_superlu_block_iter_string=`echo $lsc_superlu_block_iter_string " <TD>  " $av" </TD> "`

    file=$dir/CasePrecondLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_iter_string=`echo $full_iter_string " <TD>  " $av" </TD> "`




    file=$dir/CaseSuperLU/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_solve_string=`echo $superlu_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocks/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_solve_string=`echo $superlu_block_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocksLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    lsc_superlu_block_solve_string=`echo $lsc_superlu_block_solve_string " <TD>   " $av" </TD> "`


    file=$dir/CasePrecondLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_solve_string=`echo $full_solve_string " <TD>   " $av" </TD> "`

done


echo "<CENTER>Average GMRES iteration counts<TABLE BORDER=1>" >> junk.html
echo $ndof_string " </TR>" >> junk.html
#echo $superlu_iter_string" </TR>" >> junk.html
echo $superlu_block_iter_string" </TR>" >> junk.html
echo $lsc_superlu_block_iter_string" </TR>" >> junk.html
echo $full_iter_string" </TR>" >> junk.html
echo "</TABLE></CENTER>">> junk.html
#echo "Average GMRES iteration counts">> junk.html


echo "<CENTER>Average linear solver times (sec)<TABLE BORDER=1>" >> junk.html
echo $ndof_string " </TR>" >> junk.html
echo $superlu_solve_string" </TR>" >> junk.html
echo $superlu_block_solve_string" </TR>" >> junk.html
echo $lsc_superlu_block_solve_string" </TR>" >> junk.html
echo $full_solve_string" </TR>" >> junk.html
echo "</TABLE></CENTER>">> junk.html
#echo "Average linear solver times (sec)">> junk.html


done

exit


###############################

mesh_mult_list="3 4 " #12" # "3 6 9 12" # "20 40 80"


echo "" > junk.html
postfix_list="_STEADY _UNSTEADY"
for postfix in `echo $postfix_list`; do

echo $postfix >> junk.html

ndof_string="<TR> <TD> <CODE>n_dof</CODE> </TD> "
superlu_block_iter_string="<TR><TD>GMRES (blocks solved by SuperLU)</TD>"
lsc_superlu_block_iter_string="<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD>"
full_iter_string="<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD>"
superlu_solve_string="<TR><TD>SuperLU</TD>"
superlu_block_solve_string="<TR><TD>GMRES (blocks solved by SuperLU)</TD>"
lsc_superlu_block_solve_string="<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD>"
full_solve_string="<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD>"


for mesh_mult in `echo $mesh_mult_list`; do

    file=$dir/CaseSuperLU/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    ndofs=`grep "( ndof = " $file  | awk '{print $8}' `
    ndof=`echo $ndofs | awk '{print $1}' `
    ndof_string=`echo -n  $ndof_string " <TD> "$ndof " </TD> "` 

    file=$dir/CasePrecondSuperLUForBlocks/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_iter_string=`echo $superlu_block_iter_string " <TD>  " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocksLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    lsc_superlu_block_iter_string=`echo $lsc_superlu_block_iter_string " <TD>  " $av" </TD> "`

    file=$dir/CasePrecondLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "Linear solver iterations" $file  | awk '{print $5}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_iter_string=`echo $full_iter_string " <TD>  " $av" </TD> "`




    file=$dir/CaseSuperLU/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_solve_string=`echo $superlu_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocks/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    superlu_block_solve_string=`echo $superlu_block_solve_string " <TD>   " $av" </TD> "`

    file=$dir/CasePrecondSuperLUForBlocksLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    lsc_superlu_block_solve_string=`echo $lsc_superlu_block_solve_string " <TD>   " $av" </TD> "`


    file=$dir/CasePrecondLSC/RESLT_`echo $mesh_mult`/OUTPUT`echo $postfix`.0
    echo $file
    iter_list=`grep "time for linear solver"  $file  | awk '{print $8}' `
    av=`echo $iter_list | awk 'BEGIN{sum=0;count=0}{for (j=1;j<=NF;j++){sum+=$j; count++}}END{print sum/count}'`
    full_solve_string=`echo $full_solve_string " <TD>   " $av" </TD> "`

done


echo "<CENTER>Average GMRES iteration counts"<TABLE BORDER=1>" >> junk.html
echo $ndof_string " </TR>" >> junk.html
#echo $superlu_iter_string" </TR>" >> junk.html
echo $superlu_block_iter_string" </TR>" >> junk.html
echo $lsc_superlu_block_iter_string" </TR>" >> junk.html
echo $full_iter_string" </TR>" >> junk.html
echo "</TABLE></CENTER>">> junk.html



echo "<CENTER>Average linear solver times (sec)<TABLE BORDER=1>" >> junk.html
echo $ndof_string " </TR>" >> junk.html
echo $superlu_solve_string" </TR>" >> junk.html
echo $superlu_block_solve_string" </TR>" >> junk.html
echo $lsc_superlu_block_solve_string" </TR>" >> junk.html
echo $full_solve_string" </TR>" >> junk.html
echo "</TABLE></CENTER>">> junk.html



done

exit







########################################################################



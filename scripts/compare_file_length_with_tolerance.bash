#! /bin/sh

#-----------------------------------------------------------------------
# Compare file length with tolerance (Mainly used to compare convergence
# histories of iterative linear solvers against reference data). Such data
# is too sensitive to allow fpdiff-based comparisons. If a solver
# needs one additional iteration, fpdfiff would fail...
#
# Command line arguments: file 1 
#                         file 2 
#                         max. tolerated difference in number of lines
#
#----------------------------------------------------------------------

#Get the number of lines in both files
line_count=`wc -l $1 | awk '{print $1}'`
ref_line_count=`wc -l $2 | awk '{print $1}'`
diff1=`expr $ref_line_count - $line_count`
diff2=`expr $line_count - $ref_line_count`
if [ "$diff1" -ge "$diff2" ]; then
    diff=$diff1 
else
    diff=$diff2
fi

# Tolerance for number of iterations
threshold=`echo $3`
if [ "$diff" -le "$threshold" ]; then
    echo "   [OK] Number of lines in convergence histories agree"\
         "within threshold of "$threshold
else
    echo "   [FAILED] Number of lines in convergence histories differ "\
         "by "$diff" which exceeds the threshold of "$threshold
fi

#! /bin/sh

#---------------------------------------------------------------------------
# Compare two BlockSelector objects with non-null Replacement_block_pt.
# [OK] is outputted if 
#   Row_index is the same
#   Column_index is the same
#   Wanted is the same
#   Replacement_block_pt is NOT null.
#
# Else [FAILED] is outputted.
#
#
# The above definition for equality may seen a bit strange (especially 
# requiring the Replacement_block_pt to not be null). I wanted to test the 
# equality of two BlockSelector objects for which the Replacement_block_pt 
# is not null. Now, I can just ignore the Replacement_block_pt (i.e. it does
# not have to be null), but this is just a bit more prudent.
# 
# Command line arguments: file 1
#                         file 2
# 
# The file formats should be the exact format of the BlockSelector's <<
# operator.
#
#---------------------------------------------------------------------------


FILE1=$1
FILE2=$2

F1Row_index=$(grep "Row_index = " $FILE1 | awk '{print $NF}')
F1Col_index=$(grep "Column_index = " $FILE1 | awk '{print $NF}')
F1Wanted=$(grep "Wanted = " $FILE1 | awk '{print $NF}')
F1Replacement_block_pt=$(grep "Replacement_block_pt = " $FILE1 | awk '{print $NF}')

F2Row_index=$(grep "Row_index = " $FILE2 | awk '{print $NF}')
F2Col_index=$(grep "Column_index = " $FILE2 | awk '{print $NF}')
F2Wanted=$(grep "Wanted = " $FILE2 | awk '{print $NF}')
F2Replacement_block_pt=$(grep "Replacement_block_pt = " $FILE2 | awk '{print $NF}')


if [ "$F1Replacement_block_pt" = "0" ]
then
    echo " [FAILED] BlockSelector 1's Replacement_block_pt is null"
elif [ "$F2Replacement_block_pt" = "0" ]
then
    echo " [FAILED] BlockSelector 2's Replacement_block_pt is null"
elif [ "$F1Row_index" -ne "$F2Row_index" ]
then
    echo " [FAILED] Row index is different"
elif [ "$F1Col_index" -ne "$F2Col_index" ]
then
    echo " [FAILED] Column index is different"
elif [ "$F1Wanted" -ne "$F2Wanted" ]
then
    echo " [FAILED] Wanted is different"
else
    echo " [OK] Good stuff."
fi




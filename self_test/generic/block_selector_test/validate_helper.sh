#! /bin/sh


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

if [ "$F1Row_index" -ne "$F2Row_index" ]
then
    echo " [FAILED] Row index is different"
elif [ "$F1Col_index" -ne "$F2Col_index" ]
then
    echo " [FAILED] Column index is different"
elif [ "$F1Wanted" -ne "$F2Wanted" ]
then
    echo " [FAILED] Wanted is different"
elif [ "$F2Replacement_block_pt" = "0" ]
then
    echo " [FAILED] Replacement_block_pt is null"
else
    echo " [OK] Good stuff."
fi




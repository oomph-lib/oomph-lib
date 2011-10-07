#! /bin/bash


# Remove any existing max_min.dat file
rm -f max_min.dat

# Loop over max and min
min_max_list="min max"
for min_max in `echo  $min_max_list`; do
                
    # Loop over columns in output file
    column_list="1 2 3 4 5" # PATRICKFLAG
    for column in `echo $column_list`; do

        echo -n `awk -v n="$column" -v min_max="$min_max" 'BEGIN{max=-1.0e37; min=1.0e37}{
        if ( ($1 != "ZONE")  && ($1 != "") && ($1 != "TEXT") )
            {
                if (min > $n) min=$n
                    if (max < $n) max=$n
            }
    }
    END{ if ( min_max == "max" ) print max; else print min }' *soln*.dat` " "  >> max_min.dat
  
    done # End of loop over columns in output file
    echo " "  >> max_min.dat

done # End of loop over max and min
    



# NOTE: This script doesn't work if there are any instructions such
# as TEXT to create blue timelines etc. --> add these in to the tecplot
# macro file seperately
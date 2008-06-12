#! /bin/bash

make driven_cavity

rm -rf driven_cavity_runs
mkdir driven_cavity_runs
cp driven_cavity driven_cavity_runs
cd driven_cavity_runs



# Do parameter studies
#---------------------

# Order of arguments: 
# 1: element multiplier
# 2: Use iterative solver (0/1)
# 3: Use hypre for momentum block (0/1)
# 4: Use hypre for pressure block (0/1)
# 5: Do both (0), only CR (1) or only CR (2) elements


# Loop over element types
element_list="1 2"
for element in `echo $element_list`; do

  # Loop over different spatial resolutions
  resolution_list="1 2 3 4 5 6 7"
  for resolution in `echo $resolution_list`; do


    echo "Running $resolution 0 0 0 $element"
    mkdir RESLT
    ./driven_cavity $resolution 0 0 0 $element > OUTPUT_`echo $resolution`_0_0_0_`echo $element` 
    mv RESLT RESLT_`echo $resolution`_0_0_0_`echo $element`
    echo "done"

    echo "Running $resolution 1 0 0 $element"
    mkdir RESLT
    ./driven_cavity $resolution 1 0 0 $element > OUTPUT_`echo $resolution`_1_0_0_`echo $element` 
    mv RESLT RESLT_`echo $resolution`_1_0_0_`echo $element`
    echo "done"

    echo "Running $resolution 1 0 1 $element"
    mkdir RESLT
    ./driven_cavity $resolution 1 0 1 $element > OUTPUT_`echo $resolution`_1_0_1_`echo $element` 
    mv RESLT RESLT_`echo $resolution`_1_0_1_`echo $element`
    echo "done"

    echo "Running $resolution 1 1 0 $element"
    mkdir RESLT
    ./driven_cavity $resolution 1 1 0 $element > OUTPUT_`echo $resolution`_1_1_0_`echo $element` 
    mv RESLT RESLT_`echo $resolution`_1_1_0_`echo $element`
    echo "done"

    echo "Running $resolution 1 1 1 $element"
    mkdir RESLT
    ./driven_cavity $resolution 1 1 1 $element > OUTPUT_`echo $resolution`_1_1_1_`echo $element` 
    mv RESLT RESLT_`echo $resolution`_1_1_1_`echo $element`
    echo "done"

  done

done



# Post-process
#-------------

element_name="Crouzeix Raviart"

echo "<hr><hr>" > stats.html

# Loop over element types
for element in `echo $element_list`; do


    echo "<b><center>$element_name</center></b>" >> stats.html
    echo "<table border=1><th># of dofs</th><th>SuperLU</th><th>GMRES [SuperLU,SuperLU]</th><th>GMRES [SuperLU,AMG]</th><th>GMRES [AMG,SuperLU]</th><th>GMRES [AMG,AMG]</th>" >> stats.html


    # Loop over different spatial resolutions
    for resolution in `echo $resolution_list`; do
        ls -l OUTPUT_`echo $resolution`_*_2 
        grep 'time for Newton' OUTPUT_`echo $resolution`_*_`echo $element` | awk 'BEGIN{printf "<tr><td>"}{printf " </td><td> " $10}END{printf "</td><tr>\n"}' >> stats.html
    done
    echo "</table><hr><hr>" >> stats.html

    element_name="Taylor Hood"
 
done

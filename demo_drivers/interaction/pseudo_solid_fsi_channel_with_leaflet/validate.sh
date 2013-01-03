#! /bin/sh

#Set the number of tests to be checked
NUM_TESTS=6

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation

mkdir RESLT

echo "Running precond pseudo-solid fsi leaflet (SuperLU for blocks)"
../fsi_channel_with_leaflet_precond  --superlu_for_blocks  --mesh_multiplier 3  --validate > OUTPUT_superlu_for_blocks


echo "done"
echo " " >> validation.log
echo "Precond pseudo-solid fsi leaflet (SuperLU for blocks)" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/wall_soln5.dat RESLT/soln5.dat > result_superlu.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_superlu.dat.gz \
    result_superlu.dat  >> validation.log
fi

max_iter=`grep iterations RESLT/OUTPUT_STEADY.0  | awk 'BEGIN{max=0}{if ($NF>max) max=$NF}END{print max}'`
echo " " >> validation.log
echo "Max number of GMRES iterations (steady runs): "$max_iter>> validation.log
echo " " >> validation.log
if [ $max_iter -le 40 ]; then
      echo "   [OK] -- for number of iterations below 40" >> validation.log
else
      echo "   [FAILED] -- number of iterations above 40" >> validation.log
fi


max_iter=`grep iterations RESLT/OUTPUT_UNSTEADY.0  | awk 'BEGIN{max=0}{if ($NF>max) max=$NF}END{print max}'`
echo " " >> validation.log
echo "Max number of GMRES iterations (unsteady runs): "$max_iter>> validation.log
echo " " >> validation.log
if [ $max_iter -le 30 ]; then
      echo "   [OK] -- for number of iterations below 30" >> validation.log
else
      echo "   [FAILED] -- number of iterations above 30" >> validation.log
fi

mv RESLT RESLT_superlu_for_blocks

mkdir RESLT

echo "Running precond pseudo-solid fsi leaflet (optimal solves for blocks)"
../fsi_channel_with_leaflet_precond  --mesh_multiplier 3 --validate > OUTPUT_optimal_solves_for_blocks


echo "done"
echo " " >> validation.log
echo "Precond pseudo-solid fsi leaflet (optimal solves for blocks)" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/wall_soln5.dat RESLT/soln5.dat > result_opt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_opt.dat.gz \
    result_opt.dat  >> validation.log
fi

max_iter=`grep iterations RESLT/OUTPUT_STEADY.0  | awk 'BEGIN{max=0}{if ($NF>max) max=$NF}END{print max}'`
echo " " >> validation.log
echo "Max number of GMRES iterations (steady runs): "$max_iter>> validation.log
echo " " >> validation.log
if [ $max_iter -le 50 ]; then
      echo "   [OK] -- for number of iterations below 50" >> validation.log
else
      echo "   [FAILED] -- number of iterations above 50" >> validation.log
fi


max_iter=`grep iterations RESLT/OUTPUT_UNSTEADY.0  | awk 'BEGIN{max=0}{if ($NF>max) max=$NF}END{print max}'`
echo " " >> validation.log
echo "Max number of GMRES iterations (unsteady runs): "$max_iter>> validation.log
echo " " >> validation.log
if [ $max_iter -le 35 ]; then
      echo "   [OK] -- for number of iterations below 35" >> validation.log
else
      echo "   [FAILED] -- number of iterations above 35" >> validation.log
fi

mv RESLT RESLT_optimal_solve_for_blocks


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo " "
 echo "======================================================================"
 echo " " 
 echo "All tests in" 
 echo " " 
 echo "    `pwd`    "
 echo " "
 echo "passed successfully."
 echo " "
 echo "======================================================================"
 echo " " 
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo " "
   echo "======================================================================"
   echo " " 
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo " " 
   echo "    `pwd`/Validation/validation.log"
   echo " " 
   echo "for details" 
   echo " " 
   echo "======================================================================"
   echo " "
  else 
   echo " "
   echo "======================================================================"
   echo " " 
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo " " 
   echo "    `pwd`/validate.sh"
   echo " "
   echo "======================================================================"
   echo " "
  fi
fi


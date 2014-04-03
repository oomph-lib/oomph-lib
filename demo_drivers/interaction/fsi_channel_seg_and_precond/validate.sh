#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for steady FSI collapsible channel (precond)
#--------------------------------------------------------

echo "Running steady preconditioned solver validation"
./steady_precond.bash
mv STEADY_PRECONDITIONED_COLLAPSIBLE_CHANNEL_RUNS Validation
echo "done"


# Validation for unsteady FSI collapsible channel (precond)
#----------------------------------------------------------

echo "Running unsteady preconditioned solver validation"
./unsteady_precond.bash
mv UNSTEADY_PRECONDITIONED_COLLAPSIBLE_CHANNEL_RUNS Validation
echo "done"

# Validation for FSI collapsible channel (segregated)
#----------------------------------------------------

echo "Running steady segregated solver validation"
./steady_seg.bash
mv STEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS Validation
echo "done"

echo "Running unsteady segregated solver validation"
./unsteady_seg.bash
mv UNSTEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS Validation
echo "done"



cd Validation
mv STEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS/steady.dat .
mv UNSTEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS/unsteady.dat .

cd STEADY_PRECONDITIONED_COLLAPSIBLE_CHANNEL_RUNS
cat RESLT_res2_displ_ctrl1_solver0_subsolver0_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver2_subsolver2_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver1_subsolver0_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver3_subsolver0_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver2_subsolver0_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver3_subsolver1_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver2_subsolver1_Re500_Q1.0e-2/soln6.dat RESLT_res2_displ_ctrl1_solver3_subsolver2_Re500_Q1.0e-2/soln6.dat > ../steady_precond.dat
cd ..


cd UNSTEADY_PRECONDITIONED_COLLAPSIBLE_CHANNEL_RUNS
cat RESLT_steady_res2_solver3_subsolver0_Re500_Q1.0e-2/soln5.dat  RESLT_unsteady_res2_solver3_subsolver0_Re500_Q1.0e-2/soln5.dat RESLT_steady_res2_solver3_subsolver1_Re500_Q1.0e-2/soln5.dat  RESLT_unsteady_res2_solver3_subsolver1_Re500_Q1.0e-2/soln5.dat RESLT_steady_res2_solver3_subsolver2_Re500_Q1.0e-2/soln5.dat  RESLT_unsteady_res2_solver3_subsolver2_Re500_Q1.0e-2/soln5.dat > ../unsteady_precond.dat
cd ..

echo " " >> validation.log
echo "Steady FSI channel preconditioned solver " >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/steady_precond.dat.gz \
    steady_precond.dat 0.1 1.0e-14 >> validation.log
fi


echo " " >> validation.log
echo "Unsteady FSI channel preconditioned solver " >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unst_precond.dat.gz \
    unsteady_precond.dat 0.1 3.0e-14 >> validation.log
fi

echo " " >> validation.log
echo "Steady FSI channel segregated solver " >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/steady.dat.gz \
    steady.dat 0.1 1.0e-14 >> validation.log
fi

echo " " >> validation.log
echo "Unsteady FSI channel segregated solver " >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
#Note that the higher tolerance is to cover the discrepancy between
#different processors with this coarse mesh and the irons_and_tuck
#accelerated case. Going back down to 1% or less doesn't make a significant
#difference in most places
../../../../bin/fpdiff.py ../validata/unsteady.dat.gz \
    unsteady.dat 5.0 5.0e-7 >> validation.log
fi

#Append log to the main validation log
cat validation.log >> ../../../../validation.log

cd ..


#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10

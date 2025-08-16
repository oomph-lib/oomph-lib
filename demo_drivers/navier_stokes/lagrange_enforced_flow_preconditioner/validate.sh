#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Lgr preconditioner
#-------------------------------------
cd Validation

echo "Running 2D Lagrange Enforced Flow preconditioner test "
mkdir RESLT


PROG="two_d_tilted_square"

#GENPARAM="--use_trilinos"
GENPARAM=" "
Prec1=" "
Prec2="--use_lsc"
Visc1="--visc 1"
Ang67="--ang 67"
Re100="--re 100"
Noel8="--noel 8"
Noel16="--noel 16"

../${PROG} ${GENPARAM} ${Prec1} ${Visc1} ${Ang67} ${Re100} ${Noel8} --soln_num 0 > OUTPUT0
../${PROG} ${GENPARAM} ${Prec1} ${Visc1} ${Ang67} ${Re100} ${Noel16} --soln_num 1 > OUTPUT1
../${PROG} ${GENPARAM} ${Prec2} ${Visc1} ${Ang67} ${Re100} ${Noel8} --soln_num 2 > OUTPUT2
../${PROG} ${GENPARAM} ${Prec2} ${Visc1} ${Ang67} ${Re100} ${Noel16} --soln_num 3 > OUTPUT3

echo "done"
echo " " >> validation.log
echo "2D Lagrange Enforced Flow preconditioner test" >> validation.log
echo "-----------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat \
RESLT/iter0.dat \
RESLT/iter1.dat \
RESLT/iter2.dat \
RESLT/iter3.dat \
     > two_d_tilted_square_iter.dat

cat \
RESLT/soln0.dat \
RESLT/soln1.dat \
RESLT/soln2.dat \
RESLT/soln3.dat \
     > two_d_tilted_square_soln.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else

echo "Results: " >> validation.log

../../../../bin/fpdiff.py ../validata/two_d_tilted_square_soln.dat.gz  \
         two_d_tilted_square_soln.dat 0.1 2.0e-12 >> validation.log

echo "Iteration counts (allowing for 20% variation): " >> validation.log
echo " " >> validation.log
echo "[Note: This test fails if  " >> validation.log
echo " " >> validation.log
echo "       (1) the number of linear solver iterations " >> validation.log
echo "           differs by more than 20% " >> validation.log
echo " " >> validation.log
echo "       (2) the linear solver doesn't provide a " >> validation.log
echo "           good enough solution, causing an " >> validation.log
echo "           increase in the number of Newton" >> validation.log
echo "           iterations. " >> validation.log
echo " " >> validation.log
echo "       In the latter case, not even the number of " >> validation.log
echo "       lines in the iter* files will match.] " >> validation.log


../../../../bin/fpdiff.py ../validata/two_d_tilted_square_iter.dat.gz  \
          two_d_tilted_square_iter.dat 20.0 1.0E-10 >> validation.log

fi


# Append log to main validation log
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

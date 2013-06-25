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

cd Validation



# Validation for fig2poly conversion
#-----------------------------------

echo "Checking fig2poly"
cp ../oomph_mesh.fig .
../fig2poly oomph_mesh.fig

echo "done"
echo " " >> validation.log
echo "fig2poly validation" >> validation.log
echo "-------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/oomph_mesh.fig.poly.gz   \
    oomph_mesh.fig.poly  >> validation.log
fi






# Validation for mesh generation from xfig with oomph mesh
#---------------------------------------------------------

echo "Running xfig mesh generation Poisson validation with oomph xfig mesh"
mkdir RESLT
../mesh_from_xfig_poisson ../oomph_mesh.fig.1.node ../oomph_mesh.fig.1.ele ../oomph_mesh.fig.1.poly > OUTPUT_xfig_oomph

echo "done"
echo " " >> validation.log
echo "xfig mesh generation (oomph mesh) Poisson validation" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    RESLT/soln3.dat RESLT/soln4.dat RESLT/soln5.dat \
    > oomph_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/oomph_poisson_results.dat.gz   \
    oomph_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_oomph_xfig_poisson





# Validation for mesh generation from xfig with hole mesh
#---------------------------------------------------------

echo "Running xfig mesh generation Poisson validation with xfig hole mesh"
mkdir RESLT
../mesh_from_xfig_poisson ../hole.fig.1.node ../hole.fig.1.ele ../hole.fig.1.poly > OUTPUT_xfig_hole

echo "done"
echo " " >> validation.log
echo "xfig mesh generation (hole mesh) Poisson validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    RESLT/soln3.dat RESLT/soln4.dat RESLT/soln5.dat \
    > hole_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hole_poisson_results.dat.gz   \
    hole_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_hole_xfig_poisson



# Validation for Navier Stokes on oomph xfig mesh
#------------------------------------------------

echo "Running validation for Navier Stokes on oomph xfig mesh"
mkdir RESLT
../mesh_from_xfig_navier_stokes ../oomph_mesh.fig.1.node ../oomph_mesh.fig.1.ele ../oomph_mesh.fig.1.poly > OUTPUT_xfig_navier_stokes
echo "done"
echo " " >> validation.log
echo "Validation for Navier Stokes on oomph xfig mesh" >> validation.log
echo "-----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1.dat > oomph_xfig_navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ns_results.dat.gz   \
    oomph_xfig_navier_stokes_results.dat  >> validation.log
fi

mv RESLT RESLT_oomph_xfig_navier_stokes



# Append output to global validation log file
#--------------------------------------------
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

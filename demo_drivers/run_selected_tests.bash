backup=`pwd`
#echo $backup

all_dirs=`make list_sub_dirs`
echo " "
echo " "
echo "For your convenience, here's the list of sub-directories"
echo "in the order in which they are visited during the self-tests. "
echo "Paste these into the list variable in this script and edit"
echo "to restart the self-tests following a failure halfway along."
echo " "
echo $all_dirs
echo " "


list=`echo linking optimisation bifurcation_tracking FAQ linear_solvers linear_elasticity womersley biharmonic reaction_diffusion flux_transport young_laplace time_harmonic_fourier_decomposed_linear_elasticity fourier_decomposed_helmholtz foeppl_von_karman axisym_foeppl_von_karman pml_helmholtz axisym_linear_elasticity homogenisation darcy generalised_time_harmonic_linear_elasticity mpi`
 

echo " "
echo "Running self-tests in the following directories: "
echo " "
echo $list
echo " "
echo " "

for dir in `echo $list`; do
    cd $dir
    make check
    cd $backup
done

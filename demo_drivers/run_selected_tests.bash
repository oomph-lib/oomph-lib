backup=`pwd`
#echo $backup

all_dirs=`make list_sub_dirs`
echo " "
echo " "
echo "For your convenience, here's the list of sub-directories"
echo "in the order in which they are visited during the self-tests. "
echo "Past these into the list variable in this script and edit"
echo "to restart the self-tests following a failure halfway along."
echo " "
echo $all_dirs
echo " "



list="./axisym_navier_stokes/unstructured_torus ./interaction/unstructured_fsi ./meshing/mesh_from_inline_triangle ./meshing/mesh_from_triangle "
list=`echo $list "./meshing/mesh_from_xfig_triangle ./solid/unstructured_solid ./solid/unstructured_adaptive_solid "`
list=`echo $list "./navier_stokes/jeffery_orbit ./navier_stokes/unstructured_adaptive_fs "`
list=`echo $list "./navier_stokes/unstructured_adaptive_ALE ./navier_stokes/unstructured_fluid ./navier_stokes/unstructured_adaptive_3d_ALE "`



list=`echo interaction meshing multi_physics linking optimisation bifurcation_tracking FAQ linear_solvers linear_elasticity womersley biharmonic reaction_diffusion flux_transport young_laplace time_harmonic_fourier_decomposed_linear_elasticity fourier_decomposed_helmholtz foeppl_von_karman generalised_helmholtz axisym_linear_elasticity`

echo " "
echo "Running self-tests in the following directories: "
echo " "
echo $list
echo " "
echo " "

for dir in `echo $list`; do
    cd $dir
    make check; cd $backup
done

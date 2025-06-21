#! /bin/bash


oomph-convert -z  soln*.dat; makePvd soln soln.pvd
oomph-convert -z  coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
oomph-convert -z  -p2 solar_radiation*.dat; makePvd solar_radiation solar_radiation.pvd
oomph-convert -z  -p2 atmospheric_radiation*.dat; makePvd atmospheric_radiation atmospheric_radiation.pvd



if [ $# -eq 0 ]; then
    echo "bypassing generation of output files for radiation cones"
    echo "specify step number as command line arg. "
    exit
else
    step=$1
    echo $step
    ls diffuse_radiation_cone`echo $step`_*.dat > .tmp
    ncone=`wc .tmp | awk '{print $2}'`
    rm -f .tmp
    echo "Number of cones for step $step: $ncone"
    for ((j = 0 ; j < $ncone ; j++ )); do
        oomph-convert -p2 -z diffuse_radiation_cone`echo $step`_`echo $j`.dat
        makePvd diffuse_radiation_cone`echo $step`_ diffuse_radiation_cone.pvd 
        
        oomph-convert -p2 -z diffuse_radiation_cone_max_angle`echo $step`_`echo $j`.dat
        makePvd diffuse_radiation_cone_max_angle`echo $step`_ diffuse_radiation_cone_max_angle.pvd 
        
        oomph-convert -p2 -z diffuse_radiation_cone_min_angle`echo $step`_`echo $j`.dat
        makePvd diffuse_radiation_cone_min_angle`echo $step`_ diffuse_radiation_cone_min_angle.pvd 
    done

    oomph-convert soln`echo $step`.dat
    oomph-convert coarse_soln`echo $step`.dat

    rm -f soln_for_cone.vtu
    rm -f coarse_soln_for_cone.vtu
    ln -sf soln`echo $step`.vtu soln_for_cone.vtu
    ln -sf coarse_soln`echo $step`.vtu coarse_soln_for_cone.vtu
fi



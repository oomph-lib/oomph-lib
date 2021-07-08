#! /bin/bash

oomph-convert -z  soln*.dat
makePvd soln soln.pvd
#../replace_step_by_time.bash soln.pvd trace.dat 1 7

oomph-convert -z  coarse_soln*.dat
makePvd coarse_soln coarse_soln.pvd

oomph-convert -z  exact_soln*.dat
makePvd exact_soln exact_soln.pvd

#oomph-convert -z  coarse_exact_soln*.dat
#makePvd coarse_exact_soln coarse_exact_soln.pvd

oomph-convert -p2 -z stefan_boltzmann_rays0_intpt*.dat
makePvd stefan_boltzmann_rays0_intpt stefan_boltzmann_rays.pvd

oomph-convert -p2 sample_points.dat 
oomph-convert populated_bins.dat 

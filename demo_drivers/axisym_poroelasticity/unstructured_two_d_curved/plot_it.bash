oomph-convert -z soln*.dat  
makePvd soln soln.pvd

oomph-convert -z coarse_soln*.dat
makePvd coarse_soln coarse_soln.pvd

oomph-convert -z exact_soln*.dat
makePvd exact_soln exact_soln.pvd

oomph-convert -z -p2 -o bc_elements*.dat
makePvd bc_elements bc_elements.pvd 

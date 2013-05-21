oomph-convert -z -o soln[0-9]*.dat; makePvd soln0 soln.pvd
oomph-convert -z -o coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
oomph-convert -z -p3 newton_iter*.dat; makePvd newton_iter newton_iter.pvd
oomph-convert -z -p2 -o flux_with_melt*.dat;  makePvd flux_with_melt flux_with_melt.pvd

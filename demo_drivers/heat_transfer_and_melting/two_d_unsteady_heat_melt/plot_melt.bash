oomph-convert -z -o soln[0-9]*.dat; makePvd soln0 soln.pvd
oomph-convert -z -o coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
oomph-convert -z -o exact_soln[0-9]*.dat; makePvd exact_soln0 exact_soln.pvd
oomph-convert -z -o -p2 flux_with_melt[0-9]*.dat; makePvd flux_with_melt0 flux_with_melt.pvd
oomph-convert -z -o -p2 exact_height[0-9]*.dat; makePvd exact_height0 exact_height.pvd

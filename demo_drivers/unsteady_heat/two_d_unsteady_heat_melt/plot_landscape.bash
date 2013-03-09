oomph-convert -z -o residual_landscape*.dat; makePvd residual_landscape residual_landscape.pvd
oomph-convert -p2 -z -o soln_landscape*.dat; makePvd soln_landscape soln_landscape.pvd
oomph-convert -z -o -p3 newton_iter*.dat; makePvd newton_iter newton_iter.pvd


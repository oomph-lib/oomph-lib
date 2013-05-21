oomph-convert -z -o landscape*.dat; makePvd landscape landscape.pvd
oomph-convert -z -o -p2 soln*.dat; makePvd soln soln.pvd
oomph-convert -z -o -p3 newton_iter*.dat; makePvd newton_iter newton_iter.pvd

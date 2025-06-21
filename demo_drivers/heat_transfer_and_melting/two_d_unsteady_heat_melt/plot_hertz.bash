#oomph-convert -z -o soln*.dat; makePvd soln soln.pvd
#oomph-convert -z -o coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
#oomph-convert -z -o -p2 penetrator*.dat; makePvd penetrator penetrator.pvd
#oomph-convert -z -o -p2 contact*.dat; makePvd contact contact.pvd


oomph-convert -z  soln*.dat; makePvd soln soln.pvd
oomph-convert -z  coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
oomph-convert -z  -p2 penetrator*.dat; makePvd penetrator penetrator.pvd
oomph-convert -z  -p2 contact?.dat contact??.dat contact???.dat ; makePvd contact contact.pvd
oomph-convert -z -p3 newton_iter*.dat; makePvd newton_iter newton_iter.pvd
oomph-convert -z -p2  hertz*.dat; makePvd hertz hertz.pvd



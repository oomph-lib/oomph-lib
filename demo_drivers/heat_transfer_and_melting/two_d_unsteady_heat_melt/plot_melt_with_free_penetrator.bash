oomph-convert -z  soln*.dat; makePvd soln soln.pvd
oomph-convert -z  coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
#oomph-convert -z  exact_soln*.dat; makePvd exact_soln exact_soln.pvd
#oomph-convert -z  -p2 flux_with_melt*.dat; makePvd flux_with_melt flux_with_melt.pvd
#oomph-convert -z  -p2 exact_height*.dat; makePvd exact_height exact_height.pvd
oomph-convert -z  -p2 penetrator*.dat; makePvd penetrator penetrator.pvd
oomph-convert -z  -p2 contact*.dat; makePvd contact contact.pvd
oomph-convert -z  -p2 the_contact_points*.dat; makePvd the_contact_points the_contact_points.pvd

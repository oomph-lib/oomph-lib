oomph-convert -z -o soln[0-9]*.dat; makePvd soln soln.pvd
oomph-convert -z -o coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd
oomph-convert -z -o exact_soln[0-9]*.dat; makePvd exact_soln0 exact_soln.pvd


oomph-convert -p2 -z stefan_boltzmann_rays_right*.dat
makePvd stefan_boltzmann_rays_right stefan_boltzmann_rays_right.pvd

oomph-convert -p2 -z stefan_boltzmann_rays_left*.dat
makePvd stefan_boltzmann_rays_left stefan_boltzmann_rays_left.pvd

oomph-convert -p2 -z exact_melt_surface*.dat
makePvd exact_melt_surface exact_melt_surface.pvd

oomph-convert -z -p2 sb_radiation*.dat
makePvd sb_radiation sb_radiation.pvd

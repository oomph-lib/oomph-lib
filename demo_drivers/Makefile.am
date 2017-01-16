
# Here are the sub-directories in which we actually want to build/run
#---------------------------------------------------------------------
# demo codes
#-----------

SUBDIRS = \
poisson \
unsteady_heat \
ode \
heat_transfer_and_melting \
advection_diffusion \
helmholtz \
time_harmonic_linear_elasticity \
axisym_advection_diffusion \
spherical_advection_diffusion \
steady_axisym_advection_diffusion \
navier_stokes \
axisym_navier_stokes \
polar_navier_stokes \
spherical_navier_stokes \
solid \
beam \
shell \
linear_wave \
eigenproblems \
interaction \
meshing \
multi_physics \
linking \
optimisation \
bifurcation_tracking \
FAQ \
linear_solvers \
linear_elasticity \
womersley \
biharmonic \
reaction_diffusion \
flux_transport \
young_laplace \
time_harmonic_fourier_decomposed_linear_elasticity \
fourier_decomposed_helmholtz \
foeppl_von_karman \
axisym_foeppl_von_karman \
axisym_displ_based_foeppl_von_karman \
pml_helmholtz \
axisym_linear_elasticity \
pml_fourier_decomposed_helmholtz \
homogenisation \
darcy \
axisym_poroelasticity \
generalised_time_harmonic_linear_elasticity \
linearised_axisym_navier_stokes \
linearised_free_surface_axisym_navier_stokes \
generalised_newtonian_navier_stokes \
generalised_newtonian_axisym_navier_stokes \
multigrid
if WANT_MPI
SUBDIRS += mpi
endif

list_sub_dirs:
	@echo $(SUBDIRS)

EXTRA_DIST=run_tests.bash run_selected_tests.bash









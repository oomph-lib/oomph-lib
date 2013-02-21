#! /bin/bash

#-----------------------------------------------
# Shell script to run partial set of self-tests.
# typically used to resume self-tests after a single 
# test failed halfway through the overall validation 
# procedure
#-----------------------------------------------

#cd poisson; make check; cd ..
#cd unsteady_heat; make check; cd ..
#cd advection_diffusion; make check; cd ..
#cd helmholtz; make check; cd ..
#cd spherical_advection_diffusion; make check; cd ..
#cd steady_axisym_advection_diffusion; make check; cd ..
cd navier_stokes; make check; cd ..
cd axisym_navier_stokes; make check; cd ..
#cd polar_navier_stokes; make check; cd ..
#cd solid; make check; cd ..
#cd beam; make check; cd ..
#cd shell; make check; cd ..
#cd linear_wave; make check; cd ..
#cd eigenproblems; make check; cd ..
#cd interaction; make check; cd ..
#cd meshing; make check; cd ..
cd multi_physics; make check; cd ..
#cd linking; make check; cd ..
#cd optimisation; make check; cd ..
#cd bifurcation_tracking; make check; cd ..
#cd FAQ; make check; cd ..
#cd linear_solvers; make check; cd ..
#cd biharmonic; make check; cd ..
#cd linear_elasticity; make check; cd ..
#cd womersley; make check; cd ..
#cd reaction_diffusion; make check; cd ..
#cd flux_transport; make check; cd ..
#cd young_laplace; make check; cd ..
#cd time_harmonic_fourier_decomposed_linear_elasticity; make check; cd ..
#cd mpi; make check; cd ..

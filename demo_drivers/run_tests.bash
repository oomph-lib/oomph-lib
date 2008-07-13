#! /bin/bash

#-----------------------------------------------
# Shell script to run partial set of self-tests.
# typically used to resume self-tests after a single 
# test failed halfway through the overall validation 
# procedure
#-----------------------------------------------

#cd navier_stokes; make check; cd ..
#cd axisym_navier_stokes; make check; cd ..
#cd solid; make check; cd ..
#cd beam; make check; cd ..
#cd shell; make check; cd ..
#cd linear_wave; make check; cd ..
#cd eigenproblems; make check; cd ..
#cd interaction; make check; cd ..
#cd meshing; make check; cd ..
#cd multi_physics ; pwd; make check; cd ..
#cd linking ; pwd; make check; cd ..
#cd optimisation ; pwd; make check; cd ..
cd bifurcation_tracking ; pwd; make check; cd ..
cd FAQ ; pwd; make check; cd ..
cd linear_solvers ; pwd; make check; cd ..
cd biharmonic ; pwd; make check; cd ..
cd linear_elasticity ; pwd; make check; cd ..
cd womersley ; pwd; make check; cd .. 
cd reaction_diffusion ; pwd; make check; cd .. 

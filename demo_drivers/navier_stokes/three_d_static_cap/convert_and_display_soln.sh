#!/usr/bin/sh

cd $@
ls soln*.dat | while read file
do
    ~/Repositories/oomph-lib-with-mpi/bin/oomph-convert.py $file
done

~/Repositories/oomph-lib-with-mpi/bin/makePvd soln soln.pvd

paraview soln.pvd


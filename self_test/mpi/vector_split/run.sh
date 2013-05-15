#!/bin/bash
# Helper run script
#
# Usage: ./run.sh OPTION NCORE
#
# Option 0 - compile and run the program
#        1 - compile the src directory, compile the program, then run.
#
# NCORE - the number of cores.

CURRENTDIR=`pwd`
FILENAME="vector_split"
COMPILESRC=$1
NPROC=$2

if [ "$COMPILESRC" = "1" ]; then
  cd ../../../src/ \
  && make && make install \
  && cd $CURRENTDIR \
  && make $FILENAME && mpirun -np "$NPROC" ./$FILENAME
elif [ "$COMPILESRC" = "0" ]; then
  make $FILENAME && mpirun -np "$NPROC" ./$FILENAME
else
  echo "I do not recognise argment 0: $COMPILESRC"
fi




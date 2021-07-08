#!/bin/bash

PROG="two_d_tilted_square"

GENPARAM="--use_trilinos"
#GENPARAM=" "

Prec1=" "
Prec2="--use_lsc"
Prec3="--use_lsc --use_amg_for_f --use_amg_for_p"

Visc0="--visc 0"
Visc1="--visc 1"

Ang30="--ang 30"
Ang67="--ang 67"

Re100="--re 100"

Noel8="--noel 8"
Noel16="--noel 16"
Noel32="--noel 32"

function runtest {

mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec1} ${Visc0} ${Ang30} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec1} ${Visc0} ${Ang67} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec1} ${Visc1} ${Ang30} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec1} ${Visc1} ${Ang67} ${Re100} ${Noel16}

mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec2} ${Visc0} ${Ang30} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec2} ${Visc0} ${Ang67} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec2} ${Visc1} ${Ang30} ${Re100} ${Noel16}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec2} ${Visc1} ${Ang67} ${Re100} ${Noel16}

mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec3} ${Visc0} ${Ang30} ${Re100} ${Noel32}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec3} ${Visc0} ${Ang67} ${Re100} ${Noel32}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec3} ${Visc1} ${Ang30} ${Re100} ${Noel32}
mpirun -np 1 ./${PROG} ${GENPARAM} ${Prec3} ${Visc1} ${Ang67} ${Re100} ${Noel32}

}

runtest






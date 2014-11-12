#!/bin/sh

set -o errexit
set -o nounset

dir="Validation/bdf2"
mkdir -p $dir
./basic_ode -ts bdf2 -outdir $dir

dir="Validation/imr"
mkdir -p $dir
./basic_ode -ts imr -outdir $dir

dir="Validation/tr"
mkdir -p $dir
./basic_ode -ts tr -outdir $dir

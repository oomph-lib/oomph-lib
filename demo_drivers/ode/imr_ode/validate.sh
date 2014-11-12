#!/bin/sh

set -o errexit
set -o nounset

dir="Validation/bdf2"
mkdir -p $dir
./imr_ode -ts bdf2 -outdir $dir

dir="Validation/old-imr"
mkdir -p $dir
./imr_ode -ts old-imr -outdir $dir

dir="Validation/imr"
mkdir -p $dir
./imr_ode -ts imr -outdir $dir
 
dir="Validation/tr"
mkdir -p $dir
./imr_ode -ts tr -outdir $dir

#! /bin/bash


# Uncomment this if you just want to compile but not run the
# demo driver codes.
# exit 0

#-----------------------------------------------------------------
# This script is wrapper around the actual validate.sh script
# in the demo_driver directories. It's called via TESTS
# variable in config/makefile_templates/demo_drivers
# Default: Put a timer around the execution of the 
# validate script itself but you can do anything you want here
# including not calling validate.sh at all in which case
# make check will simply check if the demo driver codes compile
#-----------------------------------------------------------------
#(\time ./validate.sh) 2> validate_sh_timing_generated_by_make_check.dat
#(\time -f "\t%E real,\t%U user,\t%S sys" ./validate.sh) 2> validate_sh_timing_generated_by_make_check.dat

(\time -p ./validate.sh) 2> validate_sh_timing_generated_by_make_check.dat

# Process immediately afterwards; no print statements etc. otherwise
# the exit code refers to that!
EXIT_CODE=$?

echo " "
echo "Time for running demo drivers:"
echo " "
cat validate_sh_timing_generated_by_make_check.dat
echo " "


#echo "EXIT_CODE: " $EXIT_CODE
if [ $EXIT_CODE -eq 0 ]; then
    exit 0
else
    exit 1
fi

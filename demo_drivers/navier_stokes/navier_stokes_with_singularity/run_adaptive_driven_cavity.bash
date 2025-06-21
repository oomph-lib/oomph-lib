#! /bin/bash


#-----------------------------
# Run directory
#-----------------------------
run_dir="NEW_RUNS_ADAPTIVE_DRIVEN_CAVITY"
if [ -e `echo $run_dir`   ] 
then
    echo " "
    echo "ERROR: Please delete directory $run_dir and try again"
    echo " "
    exit
fi
mkdir $run_dir


#-----------------------------
# MPI run command
#-----------------------------
mpi_run_command=
#`echo mpirun n2,2,2,2 `

#----------------------------------------
# Declare, make and copy important files
#----------------------------------------
make adaptive_driven_cavity
working_files="adaptive_driven_cavity.cc adaptive_driven_cavity driven_cavity.pvsm"
cp $working_files $run_dir
cd $run_dir


#-----------------------------
# Which executable 
#-----------------------------
executable=adaptive_driven_cavity


#-----------------------------
# Loop over Re 
#-----------------------------
re_list="0 100"
for re in `echo $re_list`; do
    
    command_line_flag=" --re "$re 
    echo "Command line flag: $command_line_flag "    
    case_dir_name="Case_re"$re
    echo "Case dir name    : $case_dir_name "
    mkdir $case_dir_name
    cp $working_files $case_dir_name
    cd $case_dir_name 
    mkdir RESLT
    `echo $mpi_run_command  ./$executable $command_line_flag ` > OUTPUT 
    cd RESLT
    oomph-convert -z soln*.dat; makePvd soln soln.pvd 
    oomph-convert -z coarse_soln*.dat; makePvd coarse_soln coarse_soln.pvd 
    echo " "
    echo " "
    echo "Visualise with: "
    echo " "
    echo "    paraview --state=../driven_cavity.pvsm"
    echo " "
    echo " "
    cd ..
    cd ..

done








#! /bin/bash

set -o errexit
set -o nounset

# Function to test that the given mpi compilation and run commands work
check_mpi_run_command ()
{
    _mpi_run_command="$1"
    _cxx_compile_command="$2"

    # Check that we want to run mpi tests
    if test "$_mpi_run_command" == ""; then
        return 0
    fi

    # Check that we have a complie command
    if test "$_cxx_compile_command" == ""
    then 
        echo " "
        echo "================================================="
        echo " "
        echo "WARNING (ISSUED BY OOMPH-LIB):"
        echo "------------------------------"
        echo " "
        echo "Sorry I got myself confused when parsing " 
        echo "the configure options and could not find"
        echo "the specification of the c++ mpi compiler"
        echo "via the CXX flag. Please make sure this is"
        echo "specified in configure/configure_optinios/current"
        echo "in the form: CXX=mpic++, say, (no quotes, no spaces)."
        echo " " 
        echo "NOTE: This is not necessarily a problem, but will"
        echo "keep me from checking if mpi is up and running"
        echo "before starting the self-tests."
        echo " "
        echo "================================================="
        echo " "
        return 1
    fi

    # Test the compilation of an mpi example
    full_command="$_cxx_compile_command -o bin/minimal_mpi_test bin/minimal_mpi_test.cc"
    printf "Testing mpi compilation command: $full_command\n"
    rm -f bin/minimal_mpi_test

    # Compile
    set +e
    `echo $full_command` 
    set -e

    # run it
    set +e
    `echo $_mpi_run_command bin/minimal_mpi_test ` > bin/minimal_mpi_test.out
    result=`grep 'This worked'  bin/minimal_mpi_test.out | wc | awk '{print $1}'`
    set -e

    echo "RESULT: " $result

    if [ "$result" -eq "2" ]
    then
        printf " [Passed]\n"
        return 0 
    else
        printf " [Failed]\n"
        # Report info about errors
        if [ ! -e bin/minimal_mpi_test ]; then
            echo " " 
            echo "================================================="
            echo " "
            echo "WARNING (ISSUED BY OOMPH-LIB):"
            echo "------------------------------"
            echo " "
            echo "Compilation of bin/minimal_mpi_test.cc failed."
            echo "Are you sure your c++ compiler can compile mpi code?"
            echo " "
            echo "Note: This does not necessarily indicate a problem. "
            echo "      autogen.sh tries to extract the c++ compiler"
            echo "      from the configure options assuming that it is"
            echo "      specified in the form CXX=mpic++, say (no quotes,"
            echo "      no spaces)."
            echo " "
            echo "I will continue regardless but you shouldn't be"
            echo "surprised if large numbers of mpi self-test fail..."
            echo " "
            echo "================================================="
            echo " "
            return 2
        else
            echo " " 
            echo "================================================="
            echo " "
            echo "WARNING (ISSUED BY OOMPH-LIB):"
            echo "------------------------------"
            echo " "
            echo "The mpi test code bin/minimal_mpi_test was not run"
            echo "successfully." 
            echo " " 
            echo "You may want to check the following:" 
            echo "-- Are you sure your mpi demons have been started?"
            echo "   E.g. under lam you have to use the lamboot "
            echo "   command to get mpi up and running...."
            echo "-- Is the mpi run command you specified via "
            echo "   the --with-mpi-self-tests flag in the"
            echo "   configure options valid?"
            echo " "
            echo "I will continue regardless but you shouldn't be"
            echo "surprised if large numbers of mpi self-test fail..."
            echo " "
            echo "================================================="
            echo " "
            return 3
        fi
    fi
}

# Function to test that the given mpi compilation and run commands work
# with variable numbers of processors.
check_mpi_np_run_command()
{
    _mpi_np_run_command="$1"
    _cxx_compile_command="$2"

    # Check if we want variable np mpi tests
    if test "$_mpi_np_run_command" == ""; then
        return 0
    fi

    # Check if OOMPHNP is in the command.
    if echo "$_mpi_np_run_command" | grep -qv "OOMPHNP"; then

        echo "I cannot find OOMPHNP in the run command."
        echo "OOMPHNP will be replaced by the number of processors to run the code on."
        echo "For example, if your mpi run comamnd for two processors is "
        echo " "
        echo "  mpirun -np 2"
        echo " "
        echo "then the mpi run command for self tests on a variable number of"
        echo "processors must be"
        echo " "
        echo "  mpirun -np OOMPHNP"
        echo " "
        echo "I will continue regardless but you shouldn't be"
        echo "surprised if large numbers of mpi self-test fail..."

        return 1
    fi

    # Check we have a compile command
    if test "$_cxx_compile_command" == "" 
    then  
        echo " "
        echo "================================================="
        echo " "
        echo "WARNING (ISSUED BY OOMPH-LIB):"
        echo "------------------------------"
        echo " "
        echo "Sorry I got myself confused when parsing " 
        echo "the configure options and could not find"
        echo "the specification of the c++ mpi compiler"
        echo "via the CXX flag. Please make sure this is"
        echo "specified in configure/configure_optinios/current"
        echo "in the form: CXX=mpic++, say, (no quotes, no spaces)."
        echo " " 
        echo "NOTE: This is not necessarily a problem, but will"
        echo "keep me from checking if mpi is up and running"
        echo "before starting the self-tests."
        echo " "
        echo "================================================="
        echo " "     

        return 2
    fi

    # Check that we can compile
    rm -f bin/minimal_mpi_variablenp_test
    full_command="$_cxx_compile_command -o bin/minimal_mpi_variablenp_test bin/minimal_mpi_variablenp_test.cc"
    eval "$full_command"

    if [ ! -e bin/minimal_mpi_variablenp_test ]; then
        echo " " 
        echo "================================================="
        echo " "
        echo "WARNING (ISSUED BY OOMPH-LIB):"
        echo "------------------------------"
        echo " "
        echo "Compilation of bin/minimal_mpi_test.cc usint the command"
        echo "     $full_command"
        echo "failed. Are you sure your c++ compiler can compile mpi code?"
        echo " "
        echo "Note: This does not necessarily indicate a problem. "
        echo "      autogen.sh tries to extract the c++ compiler"
        echo "      from the configure options assuming that it is"
        echo "      specified in the form CXX=mpic++, say (no quotes,"
        echo "      no spaces)."
        echo " "
        echo "I will continue regardless but you shouldn't be"
        echo "surprised if large numbers of mpi self-test fail..."
        echo " "
        echo "================================================="
        echo " "
        return 3
    fi

    # Check that we can run with variable numbers of cores
    set +e

    MPI_RUN_ON_NP_COMMAND=$(echo $_mpi_np_run_command | sed -e "s/OOMPHNP/1/g")
    $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

    MPI_RUN_ON_NP_COMMAND=$(echo $_mpi_np_run_command | sed -e "s/OOMPHNP/2/g")
    $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

    MPI_RUN_ON_NP_COMMAND=$(echo $_mpi_np_run_command | sed -e "s/OOMPHNP/3/g")
    $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

    MPI_RUN_ON_NP_COMMAND=$(echo $_mpi_np_run_command | sed -e "s/OOMPHNP/4/g")
    $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

    set -e

    rm -f bin/minimal_mpi_variablenp_test
    variablenp_outputs="mpi_seltest_np1rank0 mpi_seltest_np2rank0 mpi_seltest_np2rank1 mpi_seltest_np3rank0 mpi_seltest_np3rank1 mpi_seltest_np3rank2 mpi_seltest_np4rank0 mpi_seltest_np4rank1 mpi_seltest_np4rank2 mpi_seltest_np4rank3"

    for OUTFILE in $variablenp_outputs 
    do
        result=$(grep 'This worked'  $OUTFILE | wc -l)

        if [ "$result" -ne "1" ]
        then
            echo " " 
            echo "================================================="
            echo " "
            echo "WARNING (ISSUED BY OOMPH-LIB):"
            echo "------------------------------"
            echo " "
            echo "The mpi test code bin/minimal_mpi_test was not run"
            echo "successfully." 
            echo " " 
            echo "You may want to check the following:" 
            echo "-- Are you sure your mpi demons have been started?"
            echo "   E.g. under lam you have to use the lamboot "
            echo "   command to get mpi up and running...."
            echo "-- Is the mpi run command you specified via "
            echo "   the --with-mpi-self-tests-variablenp flag in the"
            echo "   configure options valid?"
            echo " "
            echo "I will continue regardless but you shouldn't be"
            echo "surprised if large numbers of mpi self-test fail..."
            echo " "
            echo "================================================="
            echo " "
            return 4
        fi

        rm -f $OUTFILE
    done
}


makefile_location="$1"
get_makefile_variable()
{
    echo "print-var:; @echo \$($1)" | make -f - -f $makefile_location print-var
}

# Extract c++ compilation command: define a new make command which prints
# the variables we want then call it.
cxx_compile_command="$(get_makefile_variable CXX) $(get_makefile_variable CXXFLAGS)"

# Get mpi run command for 2 processors
mpi_run_command="$(get_makefile_variable MPI_RUN_COMMAND)"

# Get mpi run command on variable number of processors.
mpi_np_run_command="$(get_makefile_variable MPI_VARIABLENP_RUN_COMMAND)"


# echo $cxx_compile_command
# echo $mpi_run_command
# echo $mpi_np_run_command

# Finally: run the tests!
check_mpi_run_command "$mpi_run_command" "$cxx_compile_command"

check_mpi_np_run_command "$mpi_np_run_command" "$cxx_compile_command"

# Success!
exit 0

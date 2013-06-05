#!/bin/sh

set -o errexit
set -o nounset

configure_options=$1



  # Extract c++ compilation command
  cxx_compile_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"CXX="); \
    if (start!=0) \
      {  \
        start+=4; rest=substr($0 ,start); \
        end=match(rest," ")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED cxx" \
    } \
  }'` 


  mpi_run_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"--with-mpi-self-tests"); \
    if (start!=0) \
      {  \
        start+=23; rest=substr($0 ,start); \
        end=match(rest,"\"")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

  if test "$mpi_run_command" != "FAILED" 
  then  
    echo " "
    echo "============================================================"
    echo " "
    echo "It appears that you want to excecute the mpi self-tests " 
    echo "by executing the mpi codes with the run command:"
    echo " " 
    echo "       " $mpi_run_command
    echo " "
    echo "I'm now going to check if we can compile and run a basic"
    echo "stand-alone mpi code."
    echo " " 

    if test "$cxx_compile_command" != "FAILED cxx" 
    then  
      echo "OK, let's try to compile the basic mpi test"
      echo "code with the compile command: " 
      echo " "            
      full_command=`echo $cxx_compile_command -o bin/minimal_mpi_test bin/minimal_mpi_test.cc` 
      echo "      "$full_command 
      echo " " 
      rm -f bin/minimal_mpi_test
      `echo $full_command` 
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
      else
        echo "Done: Executable was produced -- that's good!"
        echo " " 
        echo "Now let's run the minimal mpi test to see if mpi is up and running:" 
        echo "I'm going to run the code with the command: "
        full_command=`echo $mpi_run_command ` 
        full_command=$full_command" bin/minimal_mpi_test "
        rm -f bin/minimal_mpi_test.out
        echo " " 
        echo "      " $full_command
        echo " " 
        `echo $full_command`  > bin/minimal_mpi_test.out
        result=`grep 'This worked'  bin/minimal_mpi_test.out | wc | awk '{print $1}'`
        if [ "$result" -ne "2" ]
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
          echo "   the --with-mpi-self-tests flag in the"
          echo "   configure options valid?"
          echo " "
          echo "I will continue regardless but you shouldn't be"
          echo "surprised if large numbers of mpi self-test fail..."
          echo " "
          echo "================================================="
          echo " "
        else
          echo " " 
          echo "Done: mpi test code was executed succesfully. Good stuff."
        fi
      fi
    else
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
    fi
    echo " "
  else
    echo " "
    echo "No mpi self-tests were requested"
    echo " "
  fi  

  # Testing mpi run command on variable number of processors.
  mpi_np_run_command=`echo $configure_options | \
    awk '{ \
    start=match($0,"--with-mpi-self-tests-variablenp"); \
    if (start!=0) \
      {  \
        start+=34; rest=substr($0 ,start); \
        end=match(rest,"\"")-1; \
      } \
    } \
    END{if (start!=0) \
    { \
      print  substr(rest,1,end) \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

  if test "$mpi_np_run_command" != "FAILED" 
  then  
    echo " "
    echo "============================================================"
    echo " "
    echo "It appears that you want to excecute the mpi self-tests " 
    echo " which involves a varying number of processors by executing the "
    echo " mpi codes with the run command:"
    echo " " 
    echo "       " $mpi_np_run_command
    echo " "

    # Check if OOMPHNP is in the mpi run command.
    OOMPHNP_CHECK=`echo $mpi_np_run_command | \
      awk '{ \
      start=match($0,"OOMPHNP"); \
    } \
    END{if (start!=0) \
    { \
      print "PASSED" \
    } \
  else \
    { \
      print "FAILED" \
    } \
  }'` 

    if test "$OOMPHNP_CHECK" != "FAILED"
    then
      echo "I'm now going to check if we can compile and run a basic"
      echo "stand-alone mpi code on 1, 2, 3 and 4 cores."

      if test "$cxx_compile_command" != "FAILED cxx" 
      then  
        echo "OK, let's try to compile the basic mpi test"
        echo "code with the compile command: " 
        echo " "            
        full_command=`echo $cxx_compile_command -o bin/minimal_mpi_variablenp_test bin/minimal_mpi_variablenp_test.cc` 
        echo "      "$full_command 
        echo " " 


        rm -f bin/minimal_mpi_variablenp_test
        `echo $full_command` 
        if [ ! -e bin/minimal_mpi_variablenp_test ]; then
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
        else
          echo "Done: Executable was produced -- that's good!"
          echo " " 
          echo "Now let's run the minimal mpi test to see if mpi is up and running:" 
          echo "I'm going to run the code with the command: "
          full_command=`echo $mpi_np_run_command `
          full_command=$full_command" bin/minimal_mpi_variablenp_test "
          echo " " 
          echo "      " $full_command
          echo "where OOMPHNP = 1, 2, 3, and 4." 
          echo " " 

          # one processor                                                                 
          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/1/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/2/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/3/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          MPI_RUN_ON_NP_COMMAND=`echo $mpi_np_run_command | sed -e "s/OOMPHNP/4/g"`
          $MPI_RUN_ON_NP_COMMAND bin/minimal_mpi_variablenp_test

          rm -f bin/minimal_mpi_variablenp_test

          variablenp_outputs="mpi_seltest_np1rank0
          mpi_seltest_np2rank0 mpi_seltest_np2rank1 \
            mpi_seltest_np3rank0 mpi_seltest_np3rank1 mpi_seltest_np3rank2 \
            mpi_seltest_np4rank0 mpi_seltest_np4rank1 mpi_seltest_np4rank2 mpi_seltest_np4rank3"

          for OUTFILE in $variablenp_outputs 
          do
            result=`grep 'This worked'  $OUTFILE | wc | awk '{print $1}'`

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
            else
              echo " " 
              echo " The file $OUTFILE is correct, good stuff!"
            fi

            rm -f $OUTFILE
          done
        fi
      else
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
      fi

    else
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
    fi # OOMPHNP_CHECK
    echo " "
  else
    echo " "
    echo "No mpi self-tests were requested"
    echo " "
  fi

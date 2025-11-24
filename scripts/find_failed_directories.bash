#! /bin/bash

#=====================================================
# Helper script to list the Validation directories where
# failures are detected. Acts on validation.log
# in current directory.
#=====================================================
full_test_text='tempfile'

awk 'BEGIN{
      next_line_is_dir=0; 
      dir="--"; 
      fail_count=0; 
      total_fail_count=0;
     }
     {
      if ($1=="[FAILED]")
       {
        fail_count++;
        total_fail_count++;
        print dir " : [FAILURE NUMBER " fail_count " ]"
       }
      if ((next_line_is_dir==1)&&($0!=" "))
        {
         dir=$0; next_line_is_dir=0
        }; 
      if ($0=="Validation directory: ")
        {
         next_line_is_dir=1
         fail_count=0
        }
      }
     END{
       print " ";
       print "TOTAL FAIL COUNT: " total_fail_count
       print " "
      }
     ' validation.log > `echo $full_test_text`



failed_count=`grep -c FA validation.log`
if [ $failed_count -gt 0 ]; then 

    echo " "
    echo "Validation directories where failures were detected:"
    echo " "
    cat $full_test_text
    
    
    rerun_script=`mktemp /tmp/rerun_failed_self_tests.XXXXXX.bash`
    awk '{L=length($1); L_without_validation=L-11; if (substr($1,L_without_validation+2)=="Validation"){print "cd "substr($1,1,L_without_validation)" ; make check -k "}}' `echo $full_test_text` | sort -u > `echo $rerun_script`
    
    
    echo " "
    echo "You can re-run only the failed tests by executing the commands in:"
    echo " "
    echo "     "`echo $rerun_script`
    echo " "
    echo "e.g. by doing "
    echo " "
    echo "     source "`echo $rerun_script`
    echo " "

fi




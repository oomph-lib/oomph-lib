#! /bin/sh

# Shell script to wipe the master validation log file




# Two little functions 'borrowed' from the tecplot installation script...
OptionPrompt() 
{ 
  printf "%s " "$1" 
}
OptionRead()
{
  read Opt
  if test "$Opt" = "" ; then
    Opt=$1
  fi
  echo $Opt
}


# Wipe global validation log file?
#---------------------------------
echo " "
echo "================================================================="
echo " "
echo "I'm about to run the self-tests." 
echo " "
echo "Note: The self-test suite involves a large number of separate  "
echo "      test scripts in different directories. We will process all of "
echo "      them and report the results at the end, using the scripts "
echo "      in the directory self_test/analyse_self_test. While the "
echo "      self-tests proceed, autoconf/automake will issue its own "
echo "      success/failure messages of the form "
echo " "
echo "         PASS: validate.sh"
echo "         =================="
echo "         All 1 tests passed"
echo "         =================="
echo " "
echo "     These only indicate that a sub-test has been run, not that"
echo "     it was passed succesfully! "
echo " "
echo "================================================================="
echo " "
if [ -e validation.log ]; then
     echo " I will wipe the existing validation.log file."
     echo " "
     OptionPrompt " Is this OK? [y/n]"
     reply=`OptionRead`
     if test "$reply" != "y" -a "$reply" != "Y" ; then 
         OptionPrompt "Terminating..."
         exit
     fi
     rm -f validation.log
fi
exit 0

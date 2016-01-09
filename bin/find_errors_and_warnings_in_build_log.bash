#! /bin/bash

#======================================================
# Find errors/warnings in specified build logs
#======================================================
if [ $# -eq 0 ]; then
 echo "Usage: find_errors_and_warnings_in_build_log.bash NAMES_OF_BUILD_LOGS"
 exit 1
fi


for FILE in "$@"
do
    echo " " 
    echo "=============================================================================" 
    echo " " 
    echo "Checking for errors/warnings in: "
    echo " " 
    echo "    "$FILE
    echo " " 
    echo "=============================================================================" 
    echo " " 

    # Analyse build script from start of oomph-lib compilation -- second mention of src/generic
    stripped_file=`echo $FILE"_from_start_of_oomph_lib_build"`
    awk 'BEGIN{count_src_generic=0}{if (match($0,"src/generic")!=0){count_src_generic+=1} if (count_src_generic>2){print$0}}' $FILE > $stripped_file
    
    
    echo " "
    echo "Checking for errors in oomph-lib related part of build log:"
    echo " "
    grep -i error $stripped_file | grep -v error_ | grep -v _error | grep -v build_script | grep -v SegregatedSolverError | grep -v ErrorEstimator | grep -v ElementWithZ2ErrorEstimator | grep -v Z2ErrorEstimator | grep -v DummyErrorEstimator | grep -v OomphLibError| grep -v NewtonSolverError | grep -v "If you can't spot any error messages" | grep -v "error.lpk" | grep -v "errors.eps" | grep -v "errors.gif" | grep -v "error.eps" | grep -v "error.gif"
    
 
    echo " "
    echo "Checking for warnings in oomph-lib related part of build log:"
    echo " "
    grep -i warning $stripped_file | grep -v AUTOMATICALLY | grep -v "doxygen could be confused by a macro call without semicolon" | grep -v OomphLibWarning | grep -v  "0 warnings" | grep -v "hyperref Warning" | grep -v "LaTeX Warning" | grep -v "Oomph-lib WARNING" | grep -v "bin/find_errors_and_warnings_in_build_log.bash" | grep -v "an oomph-lib warning whenever a function" | grep -v "gcc-based compilation with debugging (-g) and full warnings"

done

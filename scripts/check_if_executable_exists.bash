#! /bin/bash

#===================================================================
# Script to check if files (executables) specified as command line 
# arguments to this script actually exist; used to find compilation 
# failures
#===================================================================
for executable in `echo $@`; do
    if [ ! -e "$executable" ]; then
        echo "LIKELY COMPILATION ERROR FOR MISSING EXECUTABLE: "`pwd`"/$executable" 
    fi
done

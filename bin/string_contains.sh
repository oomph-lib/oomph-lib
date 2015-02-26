#!/bin/sh

#---------------------------------------------------------------------------
# Taken from http://stackoverflow.com/questions/2829613/how-do-you-tell-if-a-string-contains-another-string-in-unix-shell-scripting
#
# A unix shell string contains function.
# This uses POSIX substring parameter expansion, 
# so it works in bash, dash, ksh...
#
# contains(string, substring)
#
# Returns 0 if the specified string contains the specified substring,
# otherwise returns 1.
#---------------------------------------------------------------------------

contains() {
    string="$1"
    substring="$2"
    if test "${string#*$substring}" != "$string"
    then
        return 0    # $substring is in $string
    else
        return 1    # $substring is not in $string
    fi
}

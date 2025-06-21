#! /bin/sh

#---------------------------------------------------------------------------
# Little helper script to turn include header files into symbolic
# links -- needed because we need to prefix the sources for ln -sf [in 
# its $(LN_S) incarnation] with their full directory.
#
# The script takes the source directory (e.g. ~/oomph-lib-0.0/src/poisson)
# as its argument, reads the list of files that have been copied
# to the include directory from the file include_files.list and
# writes the full filenames into include_files.list.aux
#
# Note: In principle we could perform the entire "replace copies by links"-
#       operation in a shell script but ln -sf is (apparently) not portable
#       so we use $(LN_S) in the Makefile.am to force autoconf/automake
#       to insert the correct version.
#---------------------------------------------------------------------------
rm -f include_files.list.aux
for h in `cat include_files.list`; do 
#old    echo -n " " $1"/$h" >> include_files.list.aux " "
    printf " $1/$h" >> include_files.list.aux " "
done


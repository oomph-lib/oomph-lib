#! /bin/bash

# Get relative address
demo_drivers_dir=`dirname $0`/../demo_drivers
cd $demo_drivers_dir

# Get absolute address
demo_drivers_dir=`pwd`
echo $demo_drivers_dir

# Which directories contain a Makefile.am ?
makefile_list=`find . -name 'Makefile.am'`
for file in `echo $makefile_list`; do
    cd `dirname $file`
    # Now create an index file in that directory and in the validata directory
    # (they are the ones people will want to look at from the documentation)
    # Note that the Validation directory is temporary and will not have been
    # built on the github hosted webpage
    dir_list='. validata'
    for dir in `echo $dir_list`; do
        if [ -e $dir ]; then
            cd $dir
            echo "Creating index.html in: " `pwd`
            echo "<!-- Automatically generated file; don't edit! -->" > index.html
            echo "<!DOCTYPE html>" >> index.html
            echo "<html lang=\"en\">" >> index.html
            echo "  <head>" >> index.html
            echo "    <meta charset=\"utf-8\">" >> index.html
            echo "    <title>directory listing</title>" >> index.html
            echo "    <link rel=\"stylesheet\" href=\"style.css\">" >> index.html
            echo "    <script src=\"script.js\"></script>" >> index.html
            echo "  </head>" >> index.html
            echo "  <body>" >> index.html
            echo "<h3>Directory listing of "`pwd`"</h3>" >> index.html
            echo "    <ul>" >> index.html
            list=`find . -maxdepth 1 ! -name '.*' `
            for file in `echo $list`; do
                echo "<li> <a href=\"$file\">$file</a>" >> index.html
            done
            echo "  </ul> " >> index.html
            echo "  </body>" >> index.html
            echo "</html>" >> index.html
        fi
    done    
    cd $demo_drivers_dir
done

exit

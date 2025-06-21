#! /bin/bash

# Get relative address
demo_drivers_dir=$(dirname $0)/../demo_drivers
cd $demo_drivers_dir

# Which directories contain a CMakeLists.txt ?
cmakefile_list=$(find . -name 'CMakeLists.txt')

if [ -z "${cmakefile_list}" ]; then
    echo "ERROR: Couldn't find any CMakeLists.txt files!"
    exit 1
fi

JUST_LIST_BYPRODUCTS=0
if [ $# -eq 1 ] && ([ $1 == "-l" ] || [ $1 == "--list-byproducts" ]); then
    JUST_LIST_BYPRODUCTS=1
else
    # Get absolute address
    demo_drivers_dir=$(pwd)
    echo $demo_drivers_dir
fi

# Now create an index file in each demo_driver directory and in the validata
# subdirectory (they are the ones people will want to look at from the
# documentation). Note that the Validation directory is temporary and will not
# have been built on the github hosted webpage
for file in $(echo $cmakefile_list); do
    cd $(dirname $file)
    dir_list='. validata'
    for dir in $(echo $dir_list); do
        if [ -e $dir ]; then
            if [ ${JUST_LIST_BYPRODUCTS} -eq 1 ]; then
                echo "${demo_drivers_dir}/$(dirname $file)/${dir}/index.html"
            else
                cd $dir
                echo "Creating index.html in: " $(pwd)
                echo "<!-- Automatically generated file; don't edit! -->" >index.html
                echo "<!DOCTYPE html>" >>index.html
                echo "<html lang=\"en\">" >>index.html
                echo "  <head>" >>index.html
                echo "    <meta charset=\"utf-8\">" >>index.html
                echo "    <title>directory listing</title>" >>index.html
                echo "    <link rel=\"stylesheet\" href=\"style.css\">" >>index.html
                echo "    <script src=\"script.js\"></script>" >>index.html
                echo "  </head>" >>index.html
                echo "  <body>" >>index.html
                echo "<h3>Directory listing of "$(pwd)"</h3>" >>index.html
                echo "    <ul>" >>index.html
                list=$(find . -maxdepth 1 ! -name '.*')
                for file in $(echo $list); do
                    echo "<li> <a href=\"$file\">$file</a>" >>index.html
                done
                echo "  </ul> " >>index.html
                echo "  </body>" >>index.html
                echo "</html>" >>index.html
            fi
        fi
    done
    cd $demo_drivers_dir
done

exit

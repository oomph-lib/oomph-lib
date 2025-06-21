#! /bin/sh

#=======================================================
# Shell script to build js file which contains an array
# of figures, their titles and a link to their
# corresponding documentations. Used in random picture
# show.
#=======================================================

# Variable name, start array literal
echo "var figures = [" > rotate_gifs.js

# Find all figures directories
for dir in `find .. -path '../leftovers' -prune -o -name 'figures' -print `; do

    # Find gifs in figures directories
    for gif in `find $dir -name '*.gif' ` ; do

        # Get corresponing documentation using figure directory
        link=` echo $dir"/../html/index.html " `

        # Get the title by greping for mainpage in txt files,
        # getting the first match, trimming off "\mainpage"
        # and removing control characters and speech marks
        grep -h mainpage ` echo $dir"/../*.txt" ` \
            | head -1 \
            | cut -c 10- \
            | tr -d "[:cntrl:]" \
            | tr -d \" \
            > .tmp.txt
        title=`cat .tmp.txt`
        rm -f .tmp.txt

        # Add this gif to the array
        echo "{" >> rotate_gifs.js
        echo "\"title\": \"$title\"", >> rotate_gifs.js
        echo "\"link\": \"$link\"", >> rotate_gifs.js
        echo "\"gif\": \"$gif\"" >> rotate_gifs.js
        echo "}," >> rotate_gifs.js

    done
done

# End array literal
echo "];" >> rotate_gifs.js

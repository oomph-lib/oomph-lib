#! /bin/sh

#=======================================================
# Shell script to turn build *.php script that
# randomly displays the gif files in the oomph-lib
# documentation and provides a link to the relevant
# tutorial.
# 
# Most of it is straightforard shell-scriptery
# The clever bits are a customisation of the
# script discussed at
#
#    http://www.i-fubar.com/rotation-ad-script.php
# 
#=======================================================


#------------------------------
# WRITE HEADER FOR PHP AND
# INITIALISE COUNTER
#------------------------------
echo "<?php" > rotate_gifs.php
echo "\$bannerCounter= 1;" >> rotate_gifs.php


#------------------------------
# FIND ALL FIGURES DIRECTORIES 
#------------------------------
for dir in `find .. -path '../leftovers' -prune -o -name 'figures' -print `; do 

    #------------------------------
    # FIND GIFS IN FIGURES DIRECTORIES 
    #------------------------------
    for gif in `find $dir -name '*.gif' ` ; do 

        #--------------------------------------
        # THE LINK FROM THE FIGURE SHOULD POINT
        # TO THE CORRESPONDING DOCUMENTATION
        # WHICH IS LOCATED HERE:
        #--------------------------------------
        link=` echo $dir"/../html/index.html " `


        #--------------------------------------
        # FIND THE DOXYGEN SOURCE FILE TO EXTRACT
        # THE TITLE OF THE TUTORIAL; STRIP OUT
        # THE /mainpage AND ANY QUOTATION MARKS
        #--------------------------------------
        grep mainpage ` echo $dir"/../*.txt" ` > .tmp.txt
        sed 's/\\mainpage//g' < .tmp.txt > .tmp2.txt
        sed 's/"//g' < .tmp2.txt > .tmp3.txt
        title=`cat .tmp3.txt`
        rm -f .tmp.txt .tmp2.txt .tmp3.txt

        #--------------------------------------
        # HERE'S THE BANNER CODE ITSELF
        #--------------------------------------
        echo "\$bannerCode[\$bannerCounter] = \"<CENTER><H2>$title</H2></CENTER><A HREF=\\\"$link\\\"><IMG SRC=\\\"$gif\\\" class=\\\"img-responsive centered\\\" border=0></A>\";" >> rotate_gifs.php
        echo "\$bannerCounter++;" >> rotate_gifs.php

    done
done


#------------------------------
# WRITE FOOTER FOR PHP 
#------------------------------
echo "\$bannerAdTotals = \$bannerCounter - 1;" >> rotate_gifs.php
echo "if(\$bannerAdTotals>1)" >> rotate_gifs.php
echo "{" >> rotate_gifs.php
echo "   mt_srand((double)microtime() * 1234567);" >> rotate_gifs.php
echo "   \$bannerPicked = mt_rand(1, \$bannerAdTotals);" >> rotate_gifs.php
echo "}" >> rotate_gifs.php
echo "else" >> rotate_gifs.php
echo "{" >> rotate_gifs.php
echo "   \$bannerPicked = 1;" >> rotate_gifs.php
echo "}" >> rotate_gifs.php
echo "\$bannerAd = \$bannerCode[\$bannerPicked];" >> rotate_gifs.php
echo " " >> rotate_gifs.php

echo "\$bannerAd = str_replace('\"', '\'', \$bannerAd);" >> rotate_gifs.php
echo "   print(\"document.write(\\\"\$bannerAd\\\")\");" >> rotate_gifs.php
echo "?>" >> rotate_gifs.php



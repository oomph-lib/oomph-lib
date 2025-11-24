#! /bin/bash



#======================================================
# Customise latex style file for doxygen pdflatex
# (work-around to replace xtabular by tabular
# in style file). Bug reported and acknowledged
# https://bugzilla.gnome.org/show_bug.cgi?id=732768
# Likely to be fixed at some point at which point this
# becomes unnecessary (xtabluar is better!).
#======================================================
cp doxygen.sty doxygen.sty.orig
sed 's/xtabular/tabular/g' doxygen.sty.orig > doxygen.sty 

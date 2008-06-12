#! /bin/sh


#clear all old pictures from directory
#rm -rf adaptive_driven_cavity_??.png adaptive_driven_cavity_??.eps adaptive_driven_cavity_??.gif


#run tecplot to make png and eps
#tecplot10 ./adaptive_driven_cavity.mcr

#convert pngs to gifs
convert adaptive_driven_cavity_TH.png adaptive_driven_cavity_TH.gif
convert adaptive_driven_cavity_CR.png adaptive_driven_cavity_CR.gif
convert hanging_pressure.png hanging_pressure.gif

mv -f *.gif ../figures
mv -f *.eps ../figures


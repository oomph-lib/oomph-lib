#! /bin/sh


#clear all old pictures from directory
#rm -rf adaptive_driven_cavity_??.png adaptive_driven_cavity_??.eps adaptive_driven_cavity_??.gif


#run tecplot to make png and eps
#tecplot10 ./adaptive_driven_cavity.mcr

rm -rf *.gif

for ((i = 0; i <=108; i++ ))
do
  convert mesh$i.png mesh$i.png.gif
  echo -n "$i "
done

gifmerge -l0 -2 mesh?.png.gif mesh??.png.gif mesh???.png.gif > mesh.gif
rm -rf mesh*.png.gif
mv mesh50.img.eps mesh.eps
rm -rf *.img.eps

#convert pngs to gifs
convert TH.png TH.gif
convert CR.png CR.gif
#convert SpineMesh.png SpineMesh.gif
convert Exercise.png Exercise.gif

mv -f *.gif ../figures
mv -f *.eps ../figures


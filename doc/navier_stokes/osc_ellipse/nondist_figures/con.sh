#! /bin/sh

rm -rf *.gif

for ((i = 0; i <=40; i++ ))
do
  convert wall$i.png wall$i.png.gif
  echo -n "$i "
done

gifmerge -l0 -5 wall?.png.gif wall??.png.gif  > wall.gif

mv -f wall0.img.eps ../figures/wall.eps
mv -f wall.gif ../figures

convert running_CR15.png running_CR15.png.gif
convert running_TH15.png running_TH15.png.gif
mv -f running_CR15.img.eps ../figures/CR.eps
mv -f running_CR15.png.gif ../figures/CR.gif
mv -f running_TH15.img.eps ../figures/TH.eps
mv -f running_TH15.png.gif ../figures/TH.gif

rm -f *.eps *.gif *~

echo " Done "

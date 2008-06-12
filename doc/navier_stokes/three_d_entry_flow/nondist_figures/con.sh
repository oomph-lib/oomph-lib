#! /bin/sh

rm -rf *.gif

for ((i = 1; i <=4; i++ ))
do
    convert axial_veloc$i.png axial_veloc$i.png.gif
    convert full_profiles$i.png full_profiles$i.png.gif
    echo -n "$i "
done

gifmerge -l0 axial_veloc?.png.gif > axial_veloc.gif
gifmerge -l0 full_profiles?.png.gif > full_profiles.gif

mv -f axial_veloc4.img.eps ../figures/axial_veloc.eps
mv -f axial_veloc.gif ../figures/
mv -f full_profiles4.img.eps ../figures/full_profiles.eps
mv -f full_profiles.gif ../figures/

rm -f *.eps *.gif *~

echo " Done "

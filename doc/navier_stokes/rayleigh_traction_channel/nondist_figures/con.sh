#! /bin/sh

rm -rf *.gif

for ((i = 0; i <80; i=i+4 ))
do
    convert velocity_vectors_TH$i.png velocity_vectors_TH$i.png.gif
    convert velocity_vectors_CR$i.png velocity_vectors_CR$i.png.gif
    echo -n "$i "
done

convert error.png error.gif

gifmerge -l0 velocity_vectors_TH?.png.gif velocity_vectors_TH??.png.gif \
    > velocity_vectors_TH.gif

gifmerge -l0 velocity_vectors_CR?.png.gif velocity_vectors_CR??.png.gif \
    > velocity_vectors_CR.gif

mv -f velocity_vectors_TH3.img.eps ../figures/velocity_vectors_TH.eps
mv -f velocity_vectors_TH.gif ../figures
mv -f velocity_vectors_CR3.img.eps ../figures/velocity_vectors_CR.eps
mv -f velocity_vectors_CR.gif ../figures
mv -f error.gif ../figures
mv -f error.eps ../figures

rm -f *.eps *.gif *~

echo " Done "

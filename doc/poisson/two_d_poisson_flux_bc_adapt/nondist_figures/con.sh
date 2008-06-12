#! /bin/sh

rm -rf *.gif

for ((i = 0; i <=11; i++ ))
do
convert rotate$i.png rotate$i.gif
done

gifmerge -l0 rotate?.gif rotate??.gif > rotate.gif

cp rotate0.img.eps rotate.eps

for ((i = 0; i <=11; i++ ))
do
rm -rf rotate$i.img.eps rotate$i.gif
done






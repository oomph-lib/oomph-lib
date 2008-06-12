#! /bin/sh

rm -rf *.gif

convert wave_sol.png wave_sol.gif
convert wave_error.png wave_error.gif

mv -f ./wave_sol.gif ../figures
mv -f ./wave_sol.eps ../figures

mv -f ./wave_error.gif ../figures
mv -f ./wave_error.eps ../figures

rm -f *.eps *.gif *~

echo " Done "

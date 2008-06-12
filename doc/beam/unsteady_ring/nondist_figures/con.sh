#! /bin/sh

rm -rf *.gif

convert unsteady_ring.png unsteady_ring.gif
convert trace_file.png trace_file.gif

mv -f unsteady_ring.eps ../figures
mv -f unsteady_ring.gif ../figures
mv -f trace_file.gif ../figures
mv -f trace_file.eps ../figures


rm -f *.eps *.gif *~

echo " Done "

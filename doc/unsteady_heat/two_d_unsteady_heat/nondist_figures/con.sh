rm -rf *.gif

for ((i = 0; i <=2; i++ ))
do
    convert un_heat_unforced_soln$i.png un_heat_unforced_soln$i.gif
    echo -n "$i "
done

cp un_heat_unforced_soln1.img.eps un_heat_unforced_soln.eps
cp un_heat_unforced_soln1.gif un_heat_unforced_soln.gif


for ((i = 0; i <=2; i++ ))
do
    rm -f un_heat_unforced_soln$i.img.eps un_heat_unforced_soln$i.gif
    echo -n "$i "
done

echo " Done "

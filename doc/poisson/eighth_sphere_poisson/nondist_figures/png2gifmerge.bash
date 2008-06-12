#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 poisson?.png.gif  > eight_sphere_poisson.gif



#gifmerge -100 -l0 refine?.png.gif refine??.png.gif > all.gif


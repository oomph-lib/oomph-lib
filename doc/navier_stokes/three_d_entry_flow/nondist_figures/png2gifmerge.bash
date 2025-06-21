#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 full_profiles?.png.gif  > full_profiles.gif
gifmerge -100 -l0 axial_veloc?.png.gif  > axial_veloc.gif



#gifmerge -100 -l0 refine?.png.gif refine??.png.gif > all.gif


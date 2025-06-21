#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -20 -l0 displ_adapt?.png.gif displ_adapt??.png.gif  displ_adapt???.png.gif   > all.gif

#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 refine?.png.gif  > all.gif



#gifmerge -100 -l0 refine?.png.gif refine??.png.gif > all.gif


#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 string?.png.gif string??.png.gif > all.gif


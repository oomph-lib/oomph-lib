#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 rotate?.png.gif rotate??.png.gif > all.gif


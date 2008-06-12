#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 elastic_fish?.png.gif    > all.gif

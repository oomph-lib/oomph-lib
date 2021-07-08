#! /bin/bash

mcom "convert # #.gif" *png
gifmerge -100 -l0 circle_as_element?.png.gif    > all.gif

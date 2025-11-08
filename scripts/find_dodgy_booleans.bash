#!/bin/bash

#Script to find the now forbidden pattern of member functions returning
#references to boolean flags;
# as well as any points in driver codes in which such functions are used.
find . -name '*.h' -exec grep -H '^ *bool *&[ A-Za-z_]*[(]' \{\} \;
find . \( -name '*.h' -o -name '*.cc' \) -exec grep -H '() *= *true' \{\} \;
find . \( -name '*.h' -o -name '*.cc' \) -exec grep -H '() *= *false' \{\} \;

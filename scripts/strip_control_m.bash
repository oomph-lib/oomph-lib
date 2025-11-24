#!/bin/bash

#Little script to remove any dangling carriage returns ^M's that
#get introduced by DOS

for FILE in \
`find . \( -name '*.cc' -o -name '*.h' \) -exec grep -l $'\r' \{\} \;`
 do
 echo "Stripping ^M's from $FILE"
 cat $FILE | tr -d '\015' > $FILE.tmp
 mv $FILE.tmp $FILE
 done

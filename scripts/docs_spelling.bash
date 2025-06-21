#!/bin/bash
# Check spelling in all the .txt files in the doc directory, using custom dictionary to avoid Doxygen markup and common oomph words
find ../doc -name \*.txt | while read line; do printf "\n\n\n$line\n"; cat $line | aspell --mode=html -W 2 -p ./docs_spelling_dict.aspell.en.pws list | sort | uniq | LC_COLLATE=C sort; done

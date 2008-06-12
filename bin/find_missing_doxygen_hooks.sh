#! /bin/sh

# run this after documentation has been built to find
# empty boxes that doxygen puts where it can't find
# source code that's supposed to have been included:

find . -name 'index.html' -exec grep -H '<pre class="fragment"></pre>' {} \;


# find dangling "<a" most likely from emacs-induced linebreak after
# a long url. doxgyen doesn't build proper links from those.
find . -name '*.txt' -exec grep -l '<[Aa] *$' {} \;

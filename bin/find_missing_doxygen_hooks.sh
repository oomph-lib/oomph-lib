#! /bin/sh

# run this after documentation has been built to find
# empty boxes that doxygen puts where it can't find
# source code that's supposed to have been included:

echo " "
echo "Checking for broken doxygen includes in"
echo " "
echo "    "`pwd`
echo " "
echo "Any files listed below contain html code"
echo "that looks suspiciously as if it's been generated from"
echo "a doxygen code include that went wrong -- usually because "
echo "the hook in the c++ source code got changed and doxygen "
echo "doesn't know what to include any more, resulting in an empty box. "
echo "Output contains one previous line of text (from the generated html file)"
echo "to facilitate tracking down where the problem occurs."
echo " "
echo "Pattern 1: "
echo " "
find . -name 'index.html' -exec grep -H -B 1 '<pre class="fragment"></pre>' {} \;
echo " "
echo "Pattern 2: "
echo " "
find . -name 'index.html' -exec grep -H -B 1 '<div class="fragment"></div><!-- fragment -->' {} \;
echo " "
echo "Done"

# find dangling "<a" most likely from emacs-induced linebreak after
# a long url. doxgyen doesn't build proper links from those.
echo " "
echo "Checking for broken hrefs:"
echo "--------------------------"
echo " "
echo "Any files listed below contain a dangling \"<a\" or \"<A\""
echo "This is most likely to arise from an emacs-induced linebreak"
echo "before a long url -- check these out:"
echo " "
find . -name '*.txt' -exec grep -l '<[Aa] *$' {} \;
echo " "
echo "Done"

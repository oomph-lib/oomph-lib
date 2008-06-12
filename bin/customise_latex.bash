#! /bin/bash

#======================================================
# Customise latex.
#
# LINE NUMBERS ARE HARD CODED FOR Doxygen 1.4.4
#
#======================================================
awk '
BEGIN{ref_count_down=3; count_down=-2; print_it=1; ref_print_it=1; found=0}
{

# Copy everything until title page
if ($1=="\\begin{titlepage}")
{
  # print $0 " is title page"
  print_it=0
}


# Recommence printing
if ($1=="\\pagenumbering{arabic}")
{
  print_it=1
}

# Stop after first bit of input
if (found==0)
{
if (substr($1,0,8)=="\\chapter")
{
  found=1
  print_it=ref_print_it
  count_down=ref_count_down
#  print "contdown printit: " count_down " " print_it
}
}


count_down=count_down-1

if (count_down==0)
{
  print_it=0
  ref_print_it=0
  ref_count_down=-2
  count_down=ref_count_down
}

if (print_it==1)
{
    print $0
}



}
END{print "\\end{document}"}
' $1

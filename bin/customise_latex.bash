#! /bin/bash



#======================================================
# Customise latex.
#
# LINE NUMBERS ARE HARD CODED FOR Doxygen 1.4.4
#
#======================================================
awk '
BEGIN{print_it=1}
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

if (print_it==1)
{
    print $0
}

}
' $1

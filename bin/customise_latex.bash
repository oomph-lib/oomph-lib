#! /bin/bash



#======================================================
# Customise latex.
#======================================================



#------------------------------------------------------
# Stage 1: Get rid of title page stuff
#------------------------------------------------------
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
' $1 > tmp.tex

#------------------------------------------------------
# Stage 2: Get rid of pdfencoding info since it causes
# problems with some latex installations.
# This suppresses the line containing "pdfencoding-unicode"
# and removes the trailing comma (actually the last characer)
# in the previous line.
#------------------------------------------------------
awk '{if (NR==1){prev_line=$0}else{if ($1=="pdfencoding=unicode"){prev_line=substr(prev_line,1,length(prev_line)-1)}else{print prev_line; prev_line=$0}}}' tmp.tex

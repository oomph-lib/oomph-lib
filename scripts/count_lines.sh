#! /bin/sh

#--------------------------------------------------------------
# Count the number of lines in library (only our stuff --
# exclude external sources)
#--------------------------------------------------------------
echo " " 
echo "----------------------------------------------------------------- " 
echo " " 
echo "[This is slow -- wait for script to complete!]" 
echo " " 
echo "----------------------------------------------------------------- " 
echo " " 
echo  "Total number of source lines [*.cc and *.h        "
echo  "in ./src and ./demo_drivers                  ] : "
wc `find src demo_drivers  \( -name '*.cc' -o -name '*.h' \) | awk '{printf "%s ", $1}' ` | awk 'END{print $1}'
echo " "

echo  "Total number of lines in doc [*.txt in doc,  "
echo  "ignoring .svn subdirectories]                  : "
wc `find doc  -name '*.txt' -o -path '*svn*' -prune | awk '{printf "%s ", $1}' ` >& .oomph_lib_doc_line_counter_aux.txt ;  awk 'END{print $1}' .oomph_lib_doc_line_counter_aux.txt
#rm .oomph_lib_doc_line_counter_aux.txt 
echo " "
echo " " 
echo "----------------------------------------------------------------- " 
echo " " 

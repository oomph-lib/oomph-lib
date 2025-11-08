#! /bin/bash

#-----------------------------------------------------------------------------
# Script to extract time-trace of free (and total) memory from
# continuously recorded top session, e.g. with continuous runs
# of top -b -n 2  (say) via oomph-lib's memory monitoring tools
# Needs to be adjusted if output from top differs. The script
# below operates on this type of output from top. Adjust
# accordingly for your machine
#-----------------------------------------------------------------------------
#
#top - 10:58:34 up 9 days, 16:24,  2 users,  load average: 0.00, 0.01, 0.00
#Tasks: 289 total,   1 running, 288 sleeping,   0 stopped,   0 zombie
#Cpu(s):  1.5%us,  0.5%sy,  0.0%ni, 96.8%id,  1.1%wa,  0.0%hi,  0.0%si,  0.0%st
#Mem:     16148M total,     2906M used,    13242M free,        0M buffers
#Swap:        0M total,        0M used,        0M free,     1364M cached
#
#  PID USER      PR  NI  VIRT  RES  SHR S %CPU %MEM    TIME+  COMMAND 
#29975 mheil     20   0  8912 1212  756 R    4  0.0   0:00.02 top  
#20595 mheil     20   0  144m  22m  11m S    0  0.1   0:02.22 emacs-gtk
#28711 mheil     20   0 95240 4164 1168 S    0  0.0   0:00.91 sshd
#28712 mheil     20   0 14008 3124 1612 S    0  0.0   0:01.17 bash 
#
#-----------------------------------------------------------------------------
echo "Processing: " $@
for file in `ls $@`; do 

    awk '{if ($1=="top"){
          # Hours, minutes and seconds from third entry in 
          # in first line of top output
          h=substr($3,1,2)
          m=substr($3,4,2)
          s=substr($3,7,2)
          # Combine to total seconds since midnight (obviously 
          # will not work if run goes over midnight...)
          total_s=3600*h+60*m+s
          #printf $3" "h" "m" "s" "total_s" " } 
          printf total_s " " } 
          # Free and total memory (with "M" stripped out)
          if ($1=="Mem:")
            {
             printf substr($6, 0, length($6)-1) " "
             printf substr($2, 0, length($2)-1) "\n"
            }
         }
         END{print" "}h' $file > `echo $file`_extracted

done


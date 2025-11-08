#===================================================
# Awk script to remove licence blocks
# (line starting with //LIC//)
#===================================================
{
   # Does the line start with comment '//LIC//' ?
   lic_flag=index($0,"//LIC//")
   if (lic_flag==0)
    {
     print $0
    }
}


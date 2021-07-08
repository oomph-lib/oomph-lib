#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=circular_disk

# Sources for executable
circular_disk_SOURCES = circular_disk.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_disk_LDADD = -L@libdir@ -lfoeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


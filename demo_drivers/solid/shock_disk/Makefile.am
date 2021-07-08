#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=shock_disk

# Sources for executable
shock_disk_SOURCES = shock_disk.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
shock_disk_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                   $(EXTERNAL_LIBS) $(FLIBS)

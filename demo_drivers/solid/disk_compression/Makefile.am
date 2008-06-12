#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=disk_compression

# Sources for executable
disk_compression_SOURCES = disk_compression.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
disk_compression_LDADD = -L@libdir@ -lsolid -lconstitutive \
                         -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

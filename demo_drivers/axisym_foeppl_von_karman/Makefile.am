#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=axisym_fvk

#----------------------------------------------------------------------

# Sources for executable
axisym_fvk_SOURCES = axisym_fvk.cc 

# Required libraries: 
axisym_fvk_LDADD = -L@libdir@  -laxisym_foeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

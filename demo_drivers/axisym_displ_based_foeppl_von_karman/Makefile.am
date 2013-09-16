#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=axisym_displ_based_fvk

#----------------------------------------------------------------------

# Sources for executable
axisym_displ_based_fvk_SOURCES = axisym_displ_based_fvk.cc 

# Required libraries: 
axisym_displ_based_fvk_LDADD = -L@libdir@  -laxisym_displ_based_foeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

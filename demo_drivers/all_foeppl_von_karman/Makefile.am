#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=circular_disk displacement_based_circular_disk \
axisym_fvk axisym_displ_based_fvk

#----------------------------------------------------------------------

# Sources for executable
circular_disk_SOURCES = circular_disk.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_disk_LDADD = -L@libdir@ -lfoeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
displacement_based_circular_disk_SOURCES = displacement_based_circular_disk.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
displacement_based_circular_disk_LDADD = -L@libdir@ -lfoeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
axisym_displ_based_fvk_SOURCES = axisym_displ_based_fvk.cc 

# Required libraries: 
axisym_displ_based_fvk_LDADD = -L@libdir@  -laxisym_displ_based_foeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#----------------------------------------------------------------------

# Sources for executable
axisym_fvk_SOURCES = axisym_fvk.cc 

# Required libraries: 
axisym_fvk_LDADD = -L@libdir@  -laxisym_foeppl_von_karman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


EXTRA_DIST = run_all_versions.bash compare_u_r.lay  compare_w.lay


#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
turek_flag 



#---------------------------------------------------------------------

# Sources for first executable
turek_flag_SOURCES = turek_flag.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
turek_flag_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------


EXTRA_DIST += animate_cylinder_with_flag.mcr fsi1.mcr fsi2.mcr
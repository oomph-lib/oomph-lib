#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
elastic_bretherton



#---------------------------------------------------------------------

# Sources for first executable
elastic_bretherton_SOURCES = elastic_bretherton.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
elastic_bretherton_LDADD = \
 -L@libdir@ -lbeam -lnavier_stokes -lfluid_interface -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

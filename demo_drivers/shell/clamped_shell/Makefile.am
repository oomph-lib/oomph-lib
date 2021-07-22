#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=clamped_shell clamped_shell_with_arclength_cont

#-------------------------------------------------------------------

# Sources for executable
clamped_shell_SOURCES = clamped_shell.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
clamped_shell_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

# Sources for executable
clamped_shell_with_arclength_cont_SOURCES = \
 clamped_shell_with_arclength_cont.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
clamped_shell_with_arclength_cont_LDADD = \
-L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

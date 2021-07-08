#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=oscillating_shell 

#-------------------------------------------------------------------

# Sources for executable
oscillating_shell_SOURCES = oscillating_shell.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
oscillating_shell_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------------------------


EXTRA_DIST += shell.mcr

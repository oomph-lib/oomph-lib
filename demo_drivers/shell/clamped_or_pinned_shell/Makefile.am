#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=clamped_or_pinned_shell 

#-------------------------------------------------------------------

# Sources for executable
clamped_or_pinned_shell_SOURCES = clamped_or_pinned_shell.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
clamped_or_pinned_shell_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------------------------

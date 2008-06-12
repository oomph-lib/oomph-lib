#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=fold hopf pitchfork

#----------------------------------------------------------------------

# Sources for executable
fold_SOURCES = fold.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fold_LDADD = -L@libdir@  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
hopf_SOURCES = hopf.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hopf_LDADD = -L@libdir@  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
pitchfork_SOURCES = pitchfork.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
pitchfork_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

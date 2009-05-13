#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=curved_pipe helical_pipe

# Sources for executable
curved_pipe_SOURCES = curved_pipe.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
curved_pipe_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
helical_pipe_SOURCES = helical_pipe.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
helical_pipe_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

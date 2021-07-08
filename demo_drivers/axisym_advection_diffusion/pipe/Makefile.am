#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS= pipe pipe_variable_diff

# Sources for first executable
pipe_SOURCES = pipe.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pipe_LDADD = -L@libdir@ -laxisym_advection_diffusion -lgeneric \
                             $(EXTERNAL_LIBS) $(FLIBS)

# Sources for first executable
pipe_variable_diff_SOURCES = pipe_variable_diff.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pipe_variable_diff_LDADD = -L@libdir@ -laxisym_advection_diffusion -lgeneric \
                             $(EXTERNAL_LIBS) $(FLIBS)

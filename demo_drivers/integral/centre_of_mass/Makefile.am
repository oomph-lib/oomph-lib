#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
	backward_step_quad \
	semi_circle_tri

# Sources for executable
backward_step_quad_SOURCES = backward_step_quad.cc
semi_circle_tri_SOURCES = semi_circle_tri.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
backward_step_quad_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
semi_circle_tri_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

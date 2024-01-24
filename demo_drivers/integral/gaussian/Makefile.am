#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
	string_quad \
	backward_step_quad \
	simple_cubic_quad \
	semi_circle_tri

# Sources for executable
string_quad_SOURCES = string_quad.cc
backward_step_quad_SOURCES = backward_step_quad.cc
simple_cubic_quad_SOURCES = simple_cubic_quad.cc
semi_circle_tri_SOURCES = semi_circle_tri.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
string_quad_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
backward_step_quad_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
simple_cubic_quad_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
semi_circle_tri_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = faceted_surface_vtu_output_tester 

# Sources for executable
faceted_surface_vtu_output_tester_SOURCES = faceted_surface_vtu_output_tester.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
faceted_surface_vtu_output_tester_LDADD = -L@libdir@ -lgeneric -lpoisson $(EXTERNAL_LIBS) $(FLIBS)
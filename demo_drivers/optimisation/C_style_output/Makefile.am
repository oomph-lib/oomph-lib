#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=c_style_output

# Sources for executable
c_style_output_SOURCES = c_style_output.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
c_style_output_LDADD = -L@libdir@  -lunsteady_heat -lnavier_stokes \
                                   -lpoisson -lbeam -lsolid -lconstitutive \
	        	           -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)




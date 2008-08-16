#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=static_fish

#-------------------------------------------------------------------

# Sources for executable
static_fish_SOURCES = static_fish.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
static_fish_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

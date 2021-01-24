#---------------------------------------------------------------------------
# Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Name of executable
check_PROGRAMS = flow_past_oscillating_cylinder

#---------------------------------------------------------------------------

flow_past_oscillating_cylinder_SOURCES = flow_past_oscillating_cylinder.cc

flow_past_oscillating_cylinder_LDADD = -L@libdir@ \
-lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

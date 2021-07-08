#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=simple_spine_channel spine_channel spine_channel2



# Sources for executable
simple_spine_channel_SOURCES = simple_spine_channel.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
simple_spine_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
spine_channel_SOURCES = spine_channel.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spine_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
spine_channel2_SOURCES = spine_channel2.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spine_channel2_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=adaptive_bubble_in_channel adaptive_drop_in_channel

# Sources for executable
adaptive_bubble_in_channel_SOURCES = adaptive_bubble_in_channel.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_bubble_in_channel_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lfluid_interface  -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
adaptive_drop_in_channel_SOURCES = adaptive_drop_in_channel.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_drop_in_channel_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lfluid_interface -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST+= bubble.mcr transfer.mcr

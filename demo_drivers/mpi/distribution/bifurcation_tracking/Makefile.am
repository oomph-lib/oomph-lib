#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=pitchfork track_pitch

#-------------------------------------------------------------------

# Sources for executable
pitchfork_SOURCES = pitchfork.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pitchfork_LDADD = \
-L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------


# Sources for executable
track_pitch_SOURCES = track_pitch.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
track_pitch_LDADD = -L@libdir@ -laxisym_navier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
#---------------------------------------------------------------------

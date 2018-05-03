#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=fold hopf pitchfork track_pitch adaptive_pitchfork \
	adaptive_hopf periodic_orbit adaptive_hopf_with_separate_meshes

#----------------------------------------------------------------------

# Sources for executable
fold_SOURCES = fold.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fold_LDADD = -L@libdir@  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
hopf_SOURCES = hopf.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hopf_LDADD = -L@libdir@  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
pitchfork_SOURCES = pitchfork.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
pitchfork_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
track_pitch_SOURCES = track_pitch.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
track_pitch_LDADD = -L@libdir@ -laxisym_navier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
#---------------------------------------------------------------------

# Sources for executable
adaptive_pitchfork_SOURCES = adaptive_pitchfork.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
adaptive_pitchfork_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
adaptive_hopf_SOURCES = adaptive_hopf.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
adaptive_hopf_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
adaptive_hopf_with_separate_meshes_SOURCES = adaptive_hopf_with_separate_meshes.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
adaptive_hopf_with_separate_meshes_LDADD = -L@libdir@ -llinearised_navier_stokes -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
periodic_orbit_SOURCES = periodic_orbit.cc 
# Required libraries
# $(FLIBS) is included in case the solver involves fortran sources
periodic_orbit_LDADD = -L@libdir@ -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += adapt_hopf_eigen.dat			   

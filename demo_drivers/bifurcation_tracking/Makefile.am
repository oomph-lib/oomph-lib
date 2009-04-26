#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=fold hopf pitchfork adaptive_pitchfork adaptive_hopf

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


EXTRA_DIST += adapt_hopf_eigen.dat			   

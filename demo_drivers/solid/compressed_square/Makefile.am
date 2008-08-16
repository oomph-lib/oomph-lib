#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=compressed_square compressed_square_adapt

#-------------------------------------------------------------------

# Sources for executable
compressed_square_SOURCES = compressed_square.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
compressed_square_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)
#-------------------------------------------------------------------

# Sources for executable
compressed_square_adapt_SOURCES = compressed_square.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
compressed_square_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)


compressed_square_adapt_CXXFLAGS = 	-DREFINE

#-------------------------------------------------------------------

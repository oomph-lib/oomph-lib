#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=spin_up spin_up_cyl refineable_spin_up refineable_spin_up_cyl 

# Sources for executable
spin_up_SOURCES = spin_up.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spin_up_LDADD = -L@libdir@ -lspherical_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
spin_up_cyl_SOURCES = spin_up_cyl.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spin_up_cyl_LDADD = -L@libdir@ -laxisym_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
refineable_spin_up_SOURCES = refineable_spin_up.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_spin_up_LDADD = -L@libdir@ -lspherical_navier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
refineable_spin_up_cyl_SOURCES = refineable_spin_up_cyl.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_spin_up_cyl_LDADD = -L@libdir@ -laxisym_navier_stokes\
                               -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)

#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=static_two_layer static_single_layer

# Sources for first executable
static_two_layer_SOURCES = static_two_layer.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
static_two_layer_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
                         -lsolid -lconstitutive \
                         -lgeneric \
                         $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
static_single_layer_SOURCES = static_single_layer.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
static_single_layer_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
                            -lsolid -lconstitutive \
                            -lgeneric \
                            $(EXTERNAL_LIBS) $(FLIBS)

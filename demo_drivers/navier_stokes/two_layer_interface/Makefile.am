#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = spine_two_layer_interface elastic_two_layer_interface

# Sources for first executable
spine_two_layer_interface_SOURCES = spine_two_layer_interface.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
spine_two_layer_interface_LDADD = -L@libdir@ -lnavier_stokes \
                                  -lfluid_interface -lgeneric \
                                  $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
elastic_two_layer_interface_SOURCES = elastic_two_layer_interface.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
elastic_two_layer_interface_LDADD = -L@libdir@ -lnavier_stokes \
                                    -lfluid_interface -lsolid -lconstitutive \
                                    -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

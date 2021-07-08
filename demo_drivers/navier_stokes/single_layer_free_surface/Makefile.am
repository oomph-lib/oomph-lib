#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = spine_single_layer elastic_single_layer

# Sources for first executable
spine_single_layer_SOURCES = spine_single_layer.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
spine_single_layer_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
			   -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
elastic_single_layer_SOURCES = elastic_single_layer.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
elastic_single_layer_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
	                     -lsolid -lconstitutive -lgeneric \
			     $(EXTERNAL_LIBS) $(FLIBS)
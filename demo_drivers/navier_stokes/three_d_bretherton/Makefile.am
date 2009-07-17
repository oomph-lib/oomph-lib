#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=three_d_breth

# Sources for executable
three_d_breth_SOURCES = three_d_breth.cc \
canyon_spine_mesh.h comb_tip_spine_mesh.h  merge_meshes.h  \
tip_spine_mesh.h comb_can_spine_mesh.h  extra_elements.h st_mesh.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_breth_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

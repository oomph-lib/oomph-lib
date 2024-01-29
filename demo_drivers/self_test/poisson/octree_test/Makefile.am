#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS=octree_test

# THE EXECUTABLE:
#----------------
# Sources the executable depends on:
octree_test_SOURCES = octree_test.cc

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
octree_test_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST+=octree_neighbours.mcr octree_edge_neighbours.mcr \
            no_edge_neighbours.mcr orientation.mcr



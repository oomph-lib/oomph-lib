#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
tree_2d tree_3d


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
tree_2d_SOURCES = tree_2d.cc

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "poisson" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
tree_2d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)



# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
tree_3d_SOURCES = tree_3d.cc

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "poisson" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
tree_3d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


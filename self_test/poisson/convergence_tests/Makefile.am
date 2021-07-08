#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
q_convergence_2d \
q_convergence_3d \
t_convergence_2d \
t_convergence_3d 


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
q_convergence_3d_SOURCES = q_convergence_3d.cc

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
q_convergence_3d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
q_convergence_2d_SOURCES = q_convergence_2d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
q_convergence_2d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
t_convergence_2d_SOURCES = t_convergence_2d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
t_convergence_2d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)




# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
t_convergence_3d_SOURCES = t_convergence_3d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
t_convergence_3d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)



EXTRA_DIST+=convergence_2d.lay  convergence_3d.lay
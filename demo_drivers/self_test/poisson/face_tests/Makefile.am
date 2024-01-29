#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
q_faces_2d \
q_faces_3d \
t_faces_2d \
t_faces_3d 

# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
q_faces_3d_SOURCES = q_faces_3d.cc

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
q_faces_3d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
q_faces_2d_SOURCES = q_faces_2d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
q_faces_2d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
t_faces_2d_SOURCES = t_faces_2d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
t_faces_2d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)




# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
t_faces_3d_SOURCES = t_faces_3d.cc 

# Note: The following only works if the libraries have been installed!

# Required libraries: Only the "generic" and "navier_stokes" libraries,
# which are accessible via the general library directory which
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
t_faces_3d_LDADD = -L@libdir@ -lpoisson  \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------
# The stuff below was introduced to speed up the compilation
# but breaks when the library is compiled with paranoid
# because the paranoid flag is not used in the compilation of
# the drivers (which include headers) but it is used
# when the library is built. FaceElements contain some
# additional variables that are only used in paranoid mode
# so the memory layout assumed by driver and library
# clashes causing a nasty seg fault.
#----------------------------------------------------------------

# #Automake doesn't have an easy way to make without optimisation
# #This is the automake generated option without CXX and CPP flags
# q_faces_2d.o : q_faces_2d.cc
# 	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CXXFLAGS) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
# 	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

# q_faces_3d.o : q_faces_3d.cc
# 	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CXXFLAGS) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
# 	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

# t_faces_2d.o : t_faces_2d.cc
# 	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CXXFLAGS) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
# 	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

# t_faces_3d.o : t_faces_3d.cc
# 	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CXXFLAGS) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
# 	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

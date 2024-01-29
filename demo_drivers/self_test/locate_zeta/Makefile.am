#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=locate_zeta_tester  locate_zeta_tester_3d locate_zeta_tester_triangle locate_zeta_tester_tetgen

#----------------


# Sources for executable
locate_zeta_tester_SOURCES = locate_zeta_tester.cc locate_zeta_tester.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
locate_zeta_tester_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------

# Sources for executable
locate_zeta_tester_3d_SOURCES = locate_zeta_tester_3d.cc locate_zeta_tester.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
locate_zeta_tester_3d_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------

# Sources the executable depends on:
locate_zeta_tester_triangle_SOURCES = locate_zeta_tester_triangle.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
locate_zeta_tester_triangle_LDADD = -L@libdir@ -lpoisson \
 -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#----------------

# Sources the executable depends on:
locate_zeta_tester_tetgen_SOURCES = locate_zeta_tester_tetgen.cc 

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "poisson" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
locate_zeta_tester_tetgen_LDADD = -L@libdir@ -lpoisson \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST += cube_hole.1.node cube_hole.1.ele cube_hole.poly cube_hole.1.face \
cube_hole.1.edge distance_2d.lay distance_3d.lay \
analyse_and_plot_spiraling.bash \
spiraling_2D.pvsm spiraling_3D.pvsm

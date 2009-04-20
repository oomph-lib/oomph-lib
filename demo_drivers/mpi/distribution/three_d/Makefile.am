#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
three_d_mesh_dist \
three_d_entry_flow \
three_d_cantilever_adapt

#----------------------------------------------------------------------

# Sources for executable
three_d_mesh_dist_SOURCES = three_d_mesh_dist.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_mesh_dist_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
three_d_entry_flow_SOURCES = three_d_entry_flow.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_entry_flow_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
three_d_cantilever_adapt_SOURCES = three_d_cantilever.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_cantilever_adapt_LDADD = \
                -L@libdir@ -lsolid -lconstitutive  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
three_d_cantilever_adapt_CXXFLAGS = -DREFINE

#----------------------------------------------------------------------



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  


#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= unstructured_three_d_fsi unstructured_three_d_fsi_work \
create_poly_for_bifurcation_from_tecplot \
potential_for_bifurcation_mesh unstructured_three_d_fsi_for_doc


#---------------------------------------------------------------------

# Sources for first executable
unstructured_three_d_fsi_for_doc_SOURCES = unstructured_three_d_fsi_for_doc.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_fsi_for_doc_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for first executable
unstructured_three_d_fsi_SOURCES = unstructured_three_d_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_fsi_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------


# Sources for executable
create_poly_for_bifurcation_from_tecplot_SOURCES = \
 create_poly_for_bifurcation_from_tecplot.cc


#---------------------------------------------------------------------


# Sources for executable
potential_for_bifurcation_mesh_SOURCES = potential_for_bifurcation_mesh.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
potential_for_bifurcation_mesh_LDADD = -L@libdir@ -lgeneric \
                                        $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for first executable
unstructured_three_d_fsi_work_SOURCES = unstructured_three_d_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_fsi_work_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------


EXTRA_DIST += fsi_fluid.1.ele  fsi_fluid.1.face  fsi_fluid.1.node \
 fsi_fluid.poly  fsi_solid.1.ele  fsi_solid.1.face \
 fsi_solid.1.node  fsi_solid.poly \
 three_d_unstructured_fluid.mcr	three_d_unstructured_solid.mcr \
fsi_bifurcation_fluid.1.ele   fsi_bifurcation_fluid.1.node  \
fsi_bifurcation_solid.1.face  fsi_bifurcation_solid.poly \
fsi_bifurcation_fluid.1.face  fsi_bifurcation_fluid.poly  \
fsi_bifurcation_solid.1.ele  fsi_bifurcation_solid.1.node \
create_poly_for_bifurcation_from_potential.bash \
create_bifurcation_from_potential.mcr \
fsi_potential_bifurcation_fluid.poly \
fsi_potential_bifurcation_solid.poly only_solid.mcr






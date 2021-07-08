#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
algebraic_free_boundary_poisson \
elastic_poisson \
geom_object_element circle \
free_boundary_poisson_with_elastic_mesh_update \
macro_element_free_boundary_poisson \
old_macro_element_free_boundary_poisson_for_doc \
doc_sparse_macro_node_update \
macro_element_free_boundary_poisson_non_ref


##################################################################

# Sources 
doc_sparse_macro_node_update_SOURCES = doc_sparse_macro_node_update.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
doc_sparse_macro_node_update_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

##################################################################

# Sources for circle code
circle_SOURCES = circle.cc circle.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circle_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

# Sources for circle as geometric object
geom_object_element_SOURCES = geom_object_element.cc\
                              circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
geom_object_element_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

# Sources for free-boundary Poisson with algebraic mesh update
algebraic_free_boundary_poisson_SOURCES = algebraic_free_boundary_poisson.cc\
                                    circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
algebraic_free_boundary_poisson_LDADD = -L@libdir@ -lpoisson \
                                        -lconstitutive -lgeneric \
                                        $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

# Sources for free-boundary Poisson with macro-element mesh update
macro_element_free_boundary_poisson_SOURCES = \
                            macro_element_free_boundary_poisson.cc\
                            circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
macro_element_free_boundary_poisson_LDADD = -L@libdir@ -lpoisson \
                                      -lgeneric \
                                      $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

# Sources for old free-boundary Poisson with macro-element mesh update
old_macro_element_free_boundary_poisson_for_doc_SOURCES = \
                       old_for_doc.cc \
                       circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
old_macro_element_free_boundary_poisson_for_doc_LDADD = \
                                      -L@libdir@ -lpoisson \
                                      -lgeneric \
                                      $(EXTERNAL_LIBS) $(FLIBS)

##################################################################

# Sources for elastic poisson
elastic_poisson_SOURCES = elastic_poisson.cc\
                          circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
elastic_poisson_LDADD = -L@libdir@ -lpoisson -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

# Sources for free-surface Poisson with elastic mesh update
free_boundary_poisson_with_elastic_mesh_update_SOURCES = \
elastic_mesh_update.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
free_boundary_poisson_with_elastic_mesh_update_LDADD = \
                                         -L@libdir@ -lpoisson -lsolid \
                                         -lconstitutive -lgeneric \
                                         $(EXTERNAL_LIBS) $(FLIBS)


##################################################################

##################################################################

# Sources for free-boundary Poisson with macro-element mesh update
macro_element_free_boundary_poisson_non_ref_SOURCES = \
                            macro_element_free_boundary_poisson_non_ref.cc\
                            circle_as_generalised_element.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
macro_element_free_boundary_poisson_non_ref_LDADD = -L@libdir@ -lpoisson \
                                      -lgeneric \
                                      $(EXTERNAL_LIBS) $(FLIBS)
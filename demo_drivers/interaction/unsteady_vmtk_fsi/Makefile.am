#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

SUBDIRS = \
vmtk_files

# Name of executables
check_PROGRAMS=\
unsteady_vmtk_fsi  

#------------------------------------------------------------------------

# Sources for executable
unsteady_vmtk_fsi_SOURCES = unsteady_vmtk_fsi.cc 

#Required libraries
unsteady_vmtk_fsi_LDADD = -L@libdir@ -lsolid -lconstitutive \
                          -lnavier_stokes -lpoisson -llinear_elasticity \
                          -lmeshes -lgeneric   $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

EXTRA_DIST+=fluid.poly solid.poly fluid.1.ele  fluid.1.face  \
            fluid.1.node fluid_quadratic_nodes.dat \
            solid.1.ele  solid.1.face  solid.1.node \
            solid_quadratic_nodes.dat \
            boundary_enumeration.dat \
            three_d_unstructured_fluid.mcr \
            three_d_unstructured_solid.mcr \
            show_snap.mcr .vtp_to_tetgen_replies.txt vmtk_files

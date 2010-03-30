#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers



# Name of executables
check_PROGRAMS=\
unsteady_vmtk_fsi  \
xda_to_poly_fsi

#------------------------------------------------------------------------


# Sources for executable
xda_to_poly_fsi_SOURCES = \
xda_to_poly_fsi.cc

#Required libraries
xda_to_poly_fsi_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------

# Sources for executable
unsteady_vmtk_fsi_SOURCES = unsteady_vmtk_fsi.cc 

#Required libraries
unsteady_vmtk_fsi_LDADD = -L@libdir@ -lsolid -lconstitutive \
                          -lnavier_stokes -lpoisson -llinear_elasticity \
                          -lmeshes -lgeneric   $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------

vmtk_clean:
	rm -f a.lua boundary_enumeration.dat \
              fluid_iliac.1.ele  fluid_iliac.1.face  fluid_iliac.1.node \
              solid_iliac.1.ele  solid_iliac.1.face  solid_iliac.1.node \
              iliac.vtu quadratic_outer_solid_boundary.dat \
              quadratic_fsi_boundary.dat \
              solid_iliac.poly fluid_iliac.poly

EXTRA_DIST += fluid_iliac_short_fine.1.ele   fluid_iliac_short_fine.1.node \
              short_fine_iliac_quadratic_fsi_boundary.dat \
              solid_iliac_short_fine.1.ele   solid_iliac_short_fine.1.node \
              fluid_iliac_short_fine.1.face  \
              short_fine_iliac_boundary_enumeration.dat  \
              short_fine_iliac_quadratic_outer_solid_boundary.dat \
              solid_iliac_short_fine.1.face \
              fluid_iliac_short_coarse.1.ele   fluid_iliac_short_coarse.1.node \
              short_coarse_iliac_quadratic_fsi_boundary.dat \
              solid_iliac_short_coarse.1.ele   solid_iliac_short_coarse.1.node \
              fluid_iliac_short_coarse.1.face  \
              short_coarse_iliac_boundary_enumeration.dat  \
              short_coarse_iliac_quadratic_outer_solid_boundary.dat \
              solid_iliac_short_coarse.1.face \
              fluid_iliac_long_coarse.1.ele   fluid_iliac_long_coarse.1.node \
              long_coarse_iliac_quadratic_fsi_boundary.dat \
              solid_iliac_long_coarse.1.ele   solid_iliac_long_coarse.1.node \
              fluid_iliac_long_coarse.1.face  \
              long_coarse_iliac_boundary_enumeration.dat  \
              long_coarse_iliac_quadratic_outer_solid_boundary.dat \
              solid_iliac_long_coarse.1.face \
              show_snap.mcr three_d_unstructured_fluid.mcr \
              three_d_unstructured_solid.mcr get_files.bash



SUBDIRS = vmtk_files
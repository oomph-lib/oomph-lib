#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=\
xda_to_poly_fsi \
snap_mesh


#------------------------------------------------------------------------


# Sources for executable
xda_to_poly_fsi_SOURCES = \
xda_to_poly_fsi.cc

#Required libraries
xda_to_poly_fsi_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------

# Sources for executable
snap_mesh_SOURCES = snap_mesh.cc 

#Required libraries
snap_mesh_LDADD = -L@libdir@ -lmesh_smoothing -lsolid -lconstitutive \
                          -lnavier_stokes -lpoisson -llinear_elasticity \
                          -lmeshes -lgeneric   $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

#clean-local: 
#	rm -rf vmtk_working_dir 

EXTRA_DIST += README vtp_to_tetgen.bash expand.mcr show_snap.mcr \
              short_coarse_mesh_files vmtk_files


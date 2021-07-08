#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Name of executable
check_PROGRAMS=turek_flag_non_fsi_alg  turek_flag_non_fsi_macro

#-----------------------------------------------------------------------

# Sources for executable
turek_flag_non_fsi_alg_SOURCES = turek_flag_non_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
turek_flag_non_fsi_alg_LDADD = -L@libdir@ -lnavier_stokes -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-----------------------------------------------------------------------


# Sources for executable
turek_flag_non_fsi_macro_SOURCES = turek_flag_non_fsi.cc

turek_flag_non_fsi_macro_CXXFLAGS = -DUSE_MACRO_ELEMENTS

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
turek_flag_non_fsi_macro_LDADD = -L@libdir@ -lnavier_stokes -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-----------------------------------------------------------------------




EXTRA_DIST += animate_cylinder_with_flag.mcr

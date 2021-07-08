#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
fsi_driven_cavity_driver \
fsi_driven_cavity_segregated_driver 


#---------------------------------------------------------------------

# Sources for first executable
fsi_driven_cavity_driver_SOURCES = fsi_driven_cavity_driver.cc \
 fsi_driven_cavity_problem.h 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_driven_cavity_driver_LDADD = \
 -L@libdir@ -lbeam -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for first executable
fsi_driven_cavity_segregated_driver_SOURCES = \
fsi_driven_cavity_segregated_driver.cc \
fsi_driven_cavity_problem.h 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_driven_cavity_segregated_driver_LDADD = \
 -L@libdir@ -lmulti_physics -lbeam -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

EXTRA_DIST += anim_flow.mcr
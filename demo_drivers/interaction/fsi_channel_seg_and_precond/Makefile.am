#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable: The prefix "noinst" means the executable doesn't
# have to be installed anywhere; check programs are only compiled/run
# during make check
check_PROGRAMS= \
fsi_collapsible_channel_segregated_driver simple_segregated_driver \
fsi_collapsible_channel_precond_driver 

# THE EXECUTABLE:
#-----------------------

# Sources for  executable
fsi_collapsible_channel_precond_driver_SOURCES = \
fsi_chan_precond_driver.cc \
fsi_chan_problem.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_precond_driver_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS)  $(FLIBS)

# THE EXECUTABLE:
#----------------

# Sources for  executable
fsi_collapsible_channel_segregated_driver_SOURCES = \
  fsi_chan_seg_driver.cc \
  fsi_chan_problem.h


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_segregated_driver_LDADD = -L@libdir@ \
                     -lmulti_physics -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)


# Sources for  executable
simple_segregated_driver_SOURCES = \
  simple_segregated_driver.cc \
  fsi_chan_problem.h


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
simple_segregated_driver_LDADD = -L@libdir@ \
                     -lmulti_physics -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST +=  \
steady_seg.bash \
steady_seg_paper.bash \
unsteady_seg.bash \
unsteady_seg_paper.bash \
steady_precond.bash \
unsteady_precond.bash

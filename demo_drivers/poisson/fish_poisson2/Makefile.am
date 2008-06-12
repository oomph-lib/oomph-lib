#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
fish_poisson_simple_adapt \
fish_poisson_adapt \
fish_poisson_no_adapt \
fish_domain \
fish_poisson_node_update

# Sources for first executable
fish_poisson_simple_adapt_SOURCES = fish_poisson_simple_adapt.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_simple_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric \
                                  $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
fish_poisson_adapt_SOURCES = fish_poisson_adapt.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for third executable
fish_domain_SOURCES = fish_domain.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
fish_domain_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for fourth executable
fish_poisson_no_adapt_SOURCES = fish_poisson_no_adapt.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_no_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric \
                              $(EXTERNAL_LIBS) $(FLIBS)

# Sources for fifth executable
fish_poisson_node_update_SOURCES = fish_poisson_node_update.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_node_update_LDADD = -L@libdir@ -lpoisson -lgeneric \
                                 $(EXTERNAL_LIBS) $(FLIBS)

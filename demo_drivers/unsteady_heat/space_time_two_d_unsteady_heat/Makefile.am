#---------------------------------------------------------------------------
# Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Name of executable
check_PROGRAMS=\
test_mixed_order_galerkin_petrov\
test_equal_order_galerkin_petrov\
test_equal_order_galerkin

#---------------------------------------------------------------------------

# Sources for executable
test_mixed_order_galerkin_petrov_SOURCES = \
test_mixed_order_galerkin_petrov.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
test_mixed_order_galerkin_petrov_LDADD = -L@libdir@ \
-lspace_time_unsteady_heat_mixed_order_galerkin_petrov \
-lspace_time_block_preconditioner \
-lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
test_equal_order_galerkin_petrov_SOURCES = \
test_equal_order_galerkin_petrov.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
test_equal_order_galerkin_petrov_LDADD = -L@libdir@ \
-lspace_time_unsteady_heat_equal_order_galerkin_petrov \
-lspace_time_block_preconditioner \
-lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
test_equal_order_galerkin_SOURCES = \
test_equal_order_galerkin.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
test_equal_order_galerkin_LDADD = -L@libdir@ \
-lspace_time_unsteady_heat_equal_order_galerkin \
-lspace_time_block_preconditioner \
-lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

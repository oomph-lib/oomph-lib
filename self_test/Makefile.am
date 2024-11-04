

# Do we want to run the (possibly rather expensive) comprehensive self tests
#---------------------------------------------------------------------------
# for all available libraries? Add additional libraries here.
#------------------------------------------------------------
if WANT_EXTENDED_SELF_TESTS
  full_self_test = \
    navier_stokes poisson matrix_matrix_multiply line_visualiser deliberately_broken_code_for_self_test_test locate_zeta faceted_surface_vtu_output

if WANT_MPI
  full_self_test += mpi
endif
else
  full_self_test = 
endif


# Append the generic folder
other_self_test = \
generic

# Here are the sub-directories in which we actually want to run self-tests
#-------------------------------------------------------------------------
# This is followed by the final analysis which scans validation.log
#------------------------------------------------------------------
# for any failures
#-----------------

SUBDIRS = \
$(full_self_test) \
$(other_self_test) \
analyse_self_tests




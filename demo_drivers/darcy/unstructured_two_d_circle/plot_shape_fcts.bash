oomph-convert -z q_shape_fct?.dat
makePvd q_shape_fct q_shape_fct.pvd

oomph-convert -z p_shape_fct*.dat
makePvd p_shape_fct p_shape_fct.pvd

oomph-convert -z coarse_soln*.dat
makePvd coarse_soln coarse_soln.pvd

oomph-convert -z -p2 outer_unit_normal*.dat
makePvd outer_unit_normal outer_unit_normal.pvd

oomph-convert -z q_shape_fct_with_project*.dat
makePvd q_shape_fct_with_project q_shape_fct_with_project.pvd

oomph-convert -z -p2 flux_interpolation_point*.dat
makePvd flux_interpolation_point flux_interpolation_point.pvd

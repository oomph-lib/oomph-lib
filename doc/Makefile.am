SUBDIRS = \
people \
index \
picture_show \
optimisation \
contact \
publications \
coding_conventions \
copyright \
intro \
FAQ \
meshes \
tar_files \
the_distribution \
example_code_list \
the_data_structure \
quick_guide \
poisson \
unsteady_heat \
linear_wave \
helmholtz \
navier_stokes \
axisym_navier_stokes \
advection_diffusion \
beam \
shell \
solid \
interaction \
multi_physics \
eigenproblems \
leftovers \
order_of_action_functions \
hijacking \
change_log \
bugs \
linear_solvers \
paraview \
preconditioners \
to_be_written \
young_laplace \
mpi \
creating_doc \
linear_elasticity \
time_harmonic_fourier_decomposed_linear_elasticity \
fourier_decomposed_helmholtz \
time_harmonic_linear_elasticity \
pml_fourier_decomposed_helmholtz \
pml_helmholtz \
axisym_linear_elasticity \
acoustic_fsi \
fourier_decomposed_acoustic_fsi \
pml_time_harmonic_linear_elasticity \
deliberately_broken_doc_for_checking \
acknowledgements \
search_results

include $(top_srcdir)/config/makefile_templates/doc

docfile = doc

all:
        # Create index.html in demo_drivers directory (needed for github
        # github-hosted webpages
	../bin/make_index_html_in_demo_driver_dirs.bash


EXTRA_DIST += oomph-lib_header.html.template oomph-lib_footer.html.template \
css js\
extra_latex_style_files navbar_toggle.js raw_doxygen


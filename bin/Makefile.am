
# Scripts do not need to be installed
noinst_SCRIPTS = txt2h.sh fig2poly 

# Run this when checking

# Instructions how  to build txt2h.sh from txt2h.sh.in
txt2h.sh: txt2h.sh.in
	rm -f txt2h.sh.aux
	echo "#! /bin/sh" > txt2h.sh.aux
	echo "#---------------------------------------------------" \
	 >> txt2h.sh.aux
	echo "#Automatically generated shell script -- do not edit" \
	 >> txt2h.sh.aux
	echo "#---------------------------------------------------" \
	 >> txt2h.sh.aux
	cat txt2h.sh.in >> txt2h.sh.aux
	echo "s/AWK/"$(AWK)"/g" > update_awk.sed
	sed  -f update_awk.sed txt2h.sh.aux > txt2h.sh
	chmod ugo+x txt2h.sh
	rm -f txt2h.sh.aux
	rm -f update_awk.sed

EXTRA_DIST=build_mesh_makefile.sh change_headers_to_links.sh fpdiff.py \
           dummy_validate.sh headers.awk autogen_helpers.sh \
           customise_dist.sh txt2h.sh.in apply_licence.sh \
	   remove_licence.sh remove_licence.awk \
           cc_and_h_licence_block.txt \
	   build_oomph_html_header.sh count_lines.sh \
           compare_file_length_with_tolerance.bash customise_latex.bash \
           oomph-convert  oomph-convert.py makePvd \
           find_missing_doxygen_hooks.sh big_diff.bash \
           make_all_distributions.bash check_mpi_command.sh \
           destruct_test_distribution.bash \
           minimal_mpi_test.cc check_oomph-lib_links.bash strip_down.bash \
           validate_ok_count minimal_mpi_variablenp_test.cc \
           make_index.py check_links_in_example_list.py \
           find_dodgy_booleans.bash \
           extract_free_mem_from_top_trace.bash strip_control_m.bash \
           parallel_self_test.py tweak_doxygen_latex_style_file.bash \
           check-makefile-subdirs.py fpsmall.py string_contains.sh \
           meld_changed_files.bash find_failed_directories.bash \
           run_svnversion_if_it_exists.bash \
           find_errors_and_warnings_in_build_log.bash \
           find_compilation_failures_for_demo_drivers.bash \
           check_if_executable_exists.bash bump_doc_txt.bash \
           run_triangle_tests_only.bash make_combined_pvd.bash \
           run_selected_self_tests_based_on_string.bash \
           wrapper_for_validate.sh meld_against_files_in_different_dir.bash \
           docs_spelling.bash  docs_spelling_dict.aspell.en.pws \
           check_for_dynamic_linking_flag.bash

fig2poly:
	(cd ../demo_drivers/meshing/mesh_from_xfig_triangle; make install)

#Cleanup -- Note: - fig2poly comes from 
#                   demo_drivers/meshing/mesh_from_xfig_triangle
#                 - create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
#                   comes from demo_drivers/meshing/mesh_from_vmtk/

clean-local:
	rm -f *~ txt2h.sh fig2poly create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh minimal_mpi_test minimal_mpi_test.out
distclean-local:
	rm -f *~ txt2h.sh fig2poly create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh minimal_mpi_test minimal_mpi_test.out

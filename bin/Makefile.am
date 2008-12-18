

# Scripts do not need to be installed
noinst_SCRIPTS = txt2h.sh

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

EXTRA_DIST=build_mesh_makefile.sh change_headers_to_links.sh fpdiff.py\
           regenerate_config_files.sh dummy_validate.sh headers.awk\
           customise_dist.sh txt2h.sh.in apply_licence.sh\
	   remove_licence.sh remove_licence.awk \
           cc_and_h_licence_block.txt \
	   build_oomph_html_header.sh count_lines.sh \
           compare_file_length_with_tolerance.bash customise_latex.bash \
           oomph-convert  oomph-convert.py makePvd \
           find_missing_doxygen_hooks.sh


#Cleanup
clean-local:
	rm -f *~ txt2h.sh
distclean-local:
	rm -f *~ txt2h.sh

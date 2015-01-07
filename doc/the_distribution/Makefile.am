the_distribution.txt: linked_make_file.am joe_cool_makefile.am

linked_make_file.am:
	rm -f linked_make_file.am
	$(LN_S) ../../user_src/jack_cool/Makefile.am linked_make_file.am

joe_cool_makefile.am:
	rm -f joe_cool_makefile.am
	$(LN_S) ../../user_drivers/joe_cool/Makefile.am joe_cool_makefile.am

include $(top_srcdir)/config/makefile_templates/doc



docfile = the_distribution


EXTRA_DIST+= makefile.sample.biowulf

# Additional files to clean (there's already a clean-local target in
# the generic makefile included above).
CLEANFILES=linked_make_file.am

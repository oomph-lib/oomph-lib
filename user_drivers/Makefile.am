# These are the src directories that will be built by make when it is run
# from the oomph-lib root directory. If you want your src directories to be
# included by make mydist or make doc you need to add them to this list. On
# the other hand you may want to build your src directories independently,
# in which case do not add them to this list.
SUBDIRS = joe_cool jack_cool

# Extra actions for clean and distclean
distclean-local:
	rm -f *~
clean-local:
	rm -f *~

doc_dirs:
	echo $(SUBDIRS)

# These are the src directories we're actually going to build:
# They are defined by configure.ac which scans 
# config/configure.ac_scripts/user_src.dir_list
# to set the USER_SRC_SUBDIRECTORIES variable
SUBDIRS =  $(USER_SRC_SUBDIRECTORIES)

# Extra actions for clean and distclean
distclean-local:
	rm -f *~
clean-local:
	rm -f *~
doc_dirs:
	echo $(USER_SRC_SUBDIRECTORIES)
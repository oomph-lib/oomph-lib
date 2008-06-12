# These are the driver directories we're actually going to build:
# They are defined by configure.ac which scans 
# config/configure.ac_scripts/user_drivers.dir_list
# to set the USER_DRIVER_SUBDIRECTORIES variable
SUBDIRS =  $(USER_DRIVER_SUBDIRECTORIES)

# Extra actions for clean and distclean
distclean-local:
	rm -f *~
clean-local:
	rm -f *~

doc_dirs:
	echo $(USER_DRIVER_SUBDIRECTORIES)
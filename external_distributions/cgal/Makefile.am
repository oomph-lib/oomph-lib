


all: 
	./build_script.bash ${GMP_INCLUDE} ${GMP_LIB} ${MPFR_INCLUDE} ${MPFR_LIB} ${BOOST_INCLUDE} ${CGAL_PERMANENT_INSTALLATION_DIR}

install:
	echo " " 
	echo "================================================================================" 
	echo "Doing nothing during install; cgal already installs itself during make process."
	echo "================================================================================" 
	echo " " 

clean: 
	rm -rf CGAL-4.11

clean-local:
	rm -rf CGAL-4.11

distclean-local:
	rm -rf CGAL-4.11
	rm -rf cgal_default_installation

EXTRA_DIST = \
CGAL-4.11.tar.xz \
build_script.bash

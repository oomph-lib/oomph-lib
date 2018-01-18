


all: 
	./build_script.bash ${GMP_LIB} ${GMP_INCLUDE} ${GMP_DIR} `echo " export CXX=${CXX};" "export CXXFLAGS=\" ${CXXFLAGS} \";" "export CC=${CC};" "export CFLAGS=\" ${CFLAGS} \";"`


install:
	echo " " 
	echo "================================================================================" 
	echo "Doing nothing during install; gmp already installs itself during make process."
	echo "================================================================================" 
	echo " " 

clean:
	rm -rf gmp-6.1.2
	rm -rf gmp_default_installation

clean-local:
	rm -rf gmp-6.1.2
	rm -rf gmp_default_installation

distclean-local:
	rm -rf gmp-6.1.2
	rm -rf gmp_default_installation

EXTRA_DIST = \
gmp-6.1.2.tar.xz \
build_script.bash

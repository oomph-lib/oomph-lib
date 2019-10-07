


all:
	./build_script.bash `echo  ${MPFR_LIB} ${MPFR_INCLUDE} ${MPFR_DIR} ${GMP_INCLUDE} ${GMP_LIB} "export CC=${CC};" "export CFLAGS=\" ${CFLAGS} \";"`


install:
	echo " " 
	echo "================================================================================" 
	echo "Doing nothing during install; mpfr already installs itself during make process."
	echo "================================================================================" 
	echo " " 

clean:
	rm -rf mpfr-3.1.6

clean-local:
	rm -rf mpfr-3.1.6

distclean-local:
	rm -rf mpfr-3.1.6
	rm -rf mpfr_default_installation

EXTRA_DIST = \
mpfr-3.1.6.tar.gz \
build_script.bash

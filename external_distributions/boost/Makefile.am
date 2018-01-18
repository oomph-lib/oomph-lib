


all: 
	./build_script.bash ${BOOST_LIB} ${BOOST_INCLUDE} ${BOOST_DIR}

install:
	echo " " 
	echo "================================================================================" 
	echo "Doing nothing during install; boost already installs itself during make process."
	echo "================================================================================" 
	echo " " 


clean:
	rm -rf boost_1_65_1
	rm -rf boost_default_installation

clean-local:
	rm -rf boost_1_65_1
	rm -rf boost_default_installation

distclean-local:
	rm -rf boost_1_65_1
	rm -rf boost_default_installation


EXTRA_DIST = \
boost_1_65_1.tar.gz \
build_script.bash



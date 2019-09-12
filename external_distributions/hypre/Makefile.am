EXTRA_DIST = \
build_scripts

#makefile variable to hold the configure command (part 1)
run_command1 = " ./configure "
if WANT_MPI
run_command1 += " --with-MPI  \\" 
else
run_command1 += " --without-MPI  \\"  
endif

#makefile variable to hold the configure command (part 2)
run_command2 = " "
if ALREADY_HAVE_BLAS
run_command2 += " --with-blas=$(BLAS_LIB)  "
else
run_command2 += "  --with-blas=dummy  "
endif
if ALREADY_HAVE_LAPACK
run_command2 += " --with-lapack=$(LAPACK_LIB) "
else
run_command2 += " --with-lapack=dummy" 
endif
run_command2 += " F77=\"$(F77)\"  CC=\"$(CC)\" CXX=\"$(CXX)\" FFLAGS=\"$(FFLAGS)\"  CFLAGS=\"$(CFLAGS)\" CXXFLAGS=\"$(CXXFLAGS)\" LDFLAGS=\"$(LDFLAGS)\" " 

all:
	echo "#! /bin/bash" > build_script.bash
	echo " " >> build_script.bash
	echo "# WARNING: THIS IS AN AUTOMATICALLY GENERATED FILE" >> build_script.bash
	echo "#          DO NOT EDIT!" >> build_script.bash
	echo " " >> build_script.bash
	cat $(srcdir)/build_scripts/build_script.part1 >> build_script.bash
	echo " " >> build_script.bash
	echo $(run_command1) >> build_script.bash
	cat $(srcdir)/build_scripts/build_script.part2 >> build_script.bash
	echo $(run_command2) >> build_script.bash
	echo " " >> build_script.bash 
	echo "# Build and install" >> build_script.bash 
	echo "make">> build_script.bash 
	echo "make install" >> build_script.bash 
	echo " ">> build_script.bash 
	echo "#Cleanup">> build_script.bash 
	echo "cd ../..">> build_script.bash 
	echo "# rm -rf hypre-2.0.0">> build_script.bash 

	chmod a+x build_script.bash

	./build_script.bash 

clean-local:
	rm -f build_script.bash
	rm -rf hypre-2.0.0


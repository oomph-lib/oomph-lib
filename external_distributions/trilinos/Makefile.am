EXTRA_DIST = \
build_scripts

# makefile variable to hold absolute address of top build directory
oomph_abs_top_srcdir=`cd $(top_builddir); pwd`

#makefile variable to hold the configure command (part 1)
run_command1 = " ./configure "
if WANT_MPI
run_command1 += " --enable-mpi  \\" 
else
run_command1 += " --disable-mpi  \\"  
endif

#makefile variable to hold the configure command (part 2)
run_command2 = " "
if ALREADY_HAVE_BLAS
run_command2 += " --with-blas=$(BLAS_LIB) "
else
run_command2 += "  --with-blas=$(oomph_abs_top_srcdir)/external_src/oomph_blas/.libs/liboomph_blas.a "
endif
if ALREADY_HAVE_LAPACK
run_command2 += " --with-lapack=$(LAPACK_LIB) "
else
run_command2 += " --with-lapack=$(oomph_abs_top_srcdir)/external_src/oomph_flapack/.libs/liboomph_flapack.a" 
endif
run_command2 += "--with-gnumake --disable-default-packages --enable-teuchos --enable-epetra --enable-epetraext  --enable-triutils --enable-ifpack --enable-aztecoo --enable-amesos --enable-ml --enable-anasazi"

run_command2 += " F77=$(F77)  CC=$(CC) CXX=$(CXX) FFLAGS=\"$(FFLAGS)\"  CFLAGS=\"$(CFLAGS)\" CXXFLAGS=\"$(CXXFLAGS)\" " 

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
	cp build_script.bash junk.bash
	cat junk.bash $(srcdir)/build_scripts/build_script.part3 > build_script.bash
	rm -f junk.bash

	chmod a+x build_script.bash

	./build_script.bash 

clean-local:
	rm -f build_script.bash
	rm -rf trilinos-9.0.2


EXTRA_DIST = \
build_scripts

# makefile variable to hold absolute address of top build directory
oomph_abs_top_srcdir=`cd $(top_builddir); pwd`

#============================TRILINOS 9=================================

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

run_command2 += " F77=\"$(F77)\"  CC=\"$(CC)\" CXX=\"$(CXX)\" FFLAGS=\"$(FFLAGS)\"  CFLAGS=\"$(CFLAGS)\" CXXFLAGS=\"$(CXXFLAGS)\" LDFLAGS=\"$(LDFLAGS)\" " 

#===============================TRILINOS 10================================

#makefile variable to hold the cmake command (part 1)
run_command3 = " cmake "
run_command3 += "-D CMAKE_BUILD_TYPE:STRING=NONE "
run_command3 += "-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF "
run_command3 += "-D Trilinos_ENABLE_TESTS:BOOL=OFF "
run_command3 += "-D Trilinos_ENABLE_EXAMPLES:BOOL=OFF "
run_command3 += "-D BUILD_SHARED_LIBS=ON "
if WANT_MPI
run_command3 += "-D TPL_ENABLE_MPI:BOOL=ON \\" 
else
run_command3 += "-D TPL_ENABLE_MPI:BOOL=OFF \\"  
endif

#makefile variable to hold the configure command (part 2)
run_command4 = " "
if ALREADY_HAVE_BLAS
run_command4 += " -D TPL_ENABLE_BLAS:BOOL=ON -D TPL_BLAS_LIBRARIES:STRING=\"$(BLAS_LIB)\" "
else
run_command4 += " -D TPL_ENABLE_BLAS:BOOL=ON  -D TPL_BLAS_LIBRARIES:STRING=\"$(oomph_abs_top_srcdir)/external_src/oomph_blas/.libs/liboomph_blas.a\" "
endif
if ALREADY_HAVE_LAPACK
run_command4 += " -D TPL_ENABLE_LAPACK:BOOL=ON  -D TPL_LAPACK_LIBRARIES:STRING=\"$(LAPACK_LIB)\"  "
else
run_command4 += " -D TPL_ENABLE_LAPACK:BOOL=ON -D TPL_LAPACK_LIBRARIES:STRING=\"$(oomph_abs_top_srcdir)/external_src/oomph_flapack/.libs/liboomph_flapack.a\" " 
endif
run_command4 += " -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF -D Trilinos_ENABLE_Teuchos:BOOL=ON -D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_Triutils:BOOL=ON -D Trilinos_ENABLE_Ifpack:BOOL=ON -D Trilinos_ENABLE_AztecOO:BOOL=ON -D Trilinos_ENABLE_Amesos:BOOL=ON -D Trilinos_ENABLE_ML:BOOL=ON -D Trilinos_ENABLE_Anasazi:BOOL=ON "
run_command4 += " -D CMAKE_Fortran_COMPILER=\"$(F77)\"  -D CMAKE_C_COMPILER=\"$(CC)\" -D CMAKE_CXX_COMPILER=\"$(CXX)\" -D CMAKE_Fortran_FLAGS:STRING=\"$(FFLAGS)\"  -D CMAKE_C_FLAGS:STRING=\"$(CFLAGS)\" -D CMAKE_CXX_FLAGS:STRING=\"$(CXXFLAGS)\" \\"


all:
	echo "#! /bin/bash" > build_script.bash
	echo " " >> build_script.bash
	echo "# WARNING: THIS IS AN AUTOMATICALLY GENERATED FILE" >> build_script.bash
	echo "#          DO NOT EDIT!" >> build_script.bash
	echo " " >> build_script.bash
if HAVE_CMAKE
	echo "CMAKE_INSTALLED=1" >> build_script.bash
else
	echo "CMAKE_INSTALLED=0" >> build_script.bash
endif
	cat $(srcdir)/build_scripts/build_script.part1 >> build_script.bash
	echo " " >> build_script.bash
	echo $(run_command1) >> build_script.bash
	cat $(srcdir)/build_scripts/build_script.part2 >> build_script.bash
	echo $(run_command2) >> build_script.bash
	cp build_script.bash junk.bash
	cat junk.bash $(srcdir)/build_scripts/build_script.part3 > build_script.bash
	rm -f junk.bash
	sed 's/--with-blas=auto//g'  build_script.bash > junk.bash
	sed 's/--with-lapack=auto//g'  junk.bash > build_script.bash
	rm -f junk.bash
	echo $(run_command3) >> build_script.bash
	cat $(srcdir)/build_scripts/build_script.part4 >> build_script.bash
	echo $(run_command4) >> build_script.bash
	cp build_script.bash junk.bash
	cat junk.bash $(srcdir)/build_scripts/build_script.part5 > build_script.bash
	rm -f junk.bash
	sed 's/-D TPL_BLAS_LIBRARIES:STRING="auto"//g' build_script.bash > junk.bash
	sed 's/-D TPL_LAPACK_LIBRARIES:STRING="auto"//g' junk.bash > build_script.bash
	rm -f junk.bash
	chmod a+x build_script.bash

	./build_script.bash 

clean-local:
	rm -f build_script.bash
	rm -rf trilinos-9.0.2
	rm -rf trilinos-10.8.5-Source
	rm -rf trilinos-11.4.3-Source
	rm -rf trilinos-11.8.1-Source

distclean-local:
	rm -f build_script.bash
	rm -rf trilinos-9.0.2
	rm -rf trilinos-10.8.5-Source
	rm -rf trilinos-11.4.3-Source
	rm -rf trilinos-11.8.1-Source
	rm -rf trilinos_default_installation

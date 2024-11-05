EXTRA_DIST = \
build_scripts

# makefile variable to hold absolute address of top build directory
oomph_abs_top_srcdir=`cd $(top_builddir); pwd`

# Define mumps tar file
mumps_tar_file=MUMPS_4.10.0.tar.gz

# Define gzipped scalapack installer
scalapack_tar_file=scalapack_installer.tgz

# Mumps/scalapack need blas
if ALREADY_HAVE_BLAS
blas_lib="$(BLAS_LIB)"
else
blas_lib="$(oomph_abs_top_srcdir)/external_src/oomph_blas/.libs/liboomph_blas.a"
endif

# Mumps/scalapack need lapack
if ALREADY_HAVE_LAPACK
lapack_lib="$(LAPACK_LIB)"
else
lapack_lib="$(oomph_abs_top_srcdir)/external_src/oomph_flapack/.libs/liboomph_flapack.a"
endif

# Default for install directory
INSTALL_DIRECTORY=`pwd`/mumps_and_scalapack_default_installation

# Scalapack build needs directory of mpi include files
mpi_include=$(MPI_INCLUDE_DIRECTORY)

all:
	echo "#! /bin/bash" > build_script.bash
	echo " " >> build_script.bash
	echo "# WARNING: THIS IS AN AUTOMATICALLY GENERATED FILE" >> build_script.bash
	echo "#          DO NOT EDIT!" >> build_script.bash
	echo " " >> build_script.bash
if WILL_BUILD_MUMPS_AND_SCALAPACK
	echo "printf \"About to build mumps etc. from within build script.\n\"" >> build_script.bash
else
	echo "printf \"Bailing out of build from script because we already have mumps\n\"" >> build_script.bash
	echo "printf \"and associated libraries installed elsewhere.\n\"" >> build_script.bash
	echo "exit" >> build_script.bash
endif
	echo " "
	echo "if [ -e ${INSTALL_DIRECTORY}/include ]; then" >> build_script.bash
	echo " printf \"Scalapack and MUMPS already appear to have been installed\n\" " >> build_script.bash
	echo " printf \"because the directory \n\n\" " >> build_script.bash
	echo " printf \"  ${INSTALL_DIRECTORY}/include \n\n\"" >> build_script.bash
	echo " printf \"already exits. I\'m not installing them again\n\n\" " >> build_script.bash
	echo " exit " >> build_script.bash
	echo " else " >> build_script.bash
	echo " printf \"Installing Scalapack and MUMPS in \n\n${INSTALL_DIRECTORY} \n\" " >> build_script.bash
	echo " fi" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "#Check settings" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ ! -e "$(lapack_lib)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"ERROR: lapack library  "$(lapack_lib)" doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo \"Well done! Found the lapack library "$(lapack_lib)" \"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(blas_lib)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"ERROR: blas library "$(blas_lib)" doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo \"Well done! Found the blas library  "$(blas_lib)" \"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ x$(mpi_include) = x ]; then"  >> build_script.bash
	echo "    echo \"ERROR: mpi include directory not set \"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo " else" >> build_script.bash
	echo "if [ ! -e "$(mpi_include)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"ERROR: mpi include directory  "$(mpi_include)" doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo \"Well done! Found the mpi include directory "$(mpi_include)"\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(mpi_include)"/mpi.h ] " >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"ERROR: mpi header file "$(mpi_include)"/mpi.h doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo  \"Well done! Found the mpi header file "$(mpi_include)"/mpi.h\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(scalapack_tar_file)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"ERROR: scalapack installer tar file "$(scalapack_tar_file)" doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo  \"Well done!  scalapack installer tar file "$(scalapack_tar_file)" \"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(mumps_tar_file)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo  \"ERROR: mumps tar file "$(mumps_tar_file)" doesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo \"Well done! mumps tar file "$(mumps_tar_file)" \"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "# Install scalapack" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ -d scalapack ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo  \"Info:   \"`pwd`\"/scalapack already exists; overwriting\"" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    mkdir scalapack" >> build_script.bash
	echo "    echo  \"Made scalapack directory for installation \" " >> build_script.bash
	echo "fi" >> build_script.bash
	echo "cp  "$(scalapack_tar_file)" scalapack" >> build_script.bash
	echo "cd scalapack" >> build_script.bash
	echo "tar xfz "$(scalapack_tar_file) >> build_script.bash
	echo "cd  scalapack_installer_0.94" >> build_script.bash
	echo "" >> build_script.bash
	echo "#Install with blacs (downloaded on the fly)" >> build_script.bash
	echo "echo \"Bypassing test of executables (because it breaks on hector, or any\"">> build_script.bash
	echo "echo \"other machine where we cross compile)\"">> build_script.bash
#OLD (PYTHON2)	echo "cp script/framework.py script/framework.py.orig" >> build_script.bash
	echo "cp ../../scalapack_python3_scripts/*py script" >> build_script.bash
#OLD (PYTHON2)	echo "sed 's/.\/tmpf/date/g' script/framework.py.orig | sed 's/.\/tmpc/date/g' > script/framework.py ">> build_script.bash
#OLD (PYTHON2)	echo "patch script/blacs.py < ../../build_scripts/blacs.py.patch  ">> build_script.bash
	echo "echo "Building... [this can take a while!]"" >> build_script.bash
	echo "python3 setup.py --notesting --downblacs --lapacklib=$(lapack_lib) --blaslib=$(blas_lib) --mpiincdir=$(mpi_include) --mpicc=$(CC) --mpif77=$(F77) --ccflags=\"$(CFLAGS)\" --fcflags=\"$(FFLAGS) -fallow-argument-mismatch\"  --ldflags_c=\"$(CFLAGS)\" --ldflags_f77=\"$(FFLAGS)\" --noopt=\"$(FFLAGS_NO_OPT)\" > scalapack.log" >> build_script.bash
	echo "" >> build_script.bash
	echo "echo \"...done!\"" >> build_script.bash
	echo "echo \"here's the scalapack.log:\"" >> build_script.bash
	echo "cat scalapack.log" >> build_script.bash
	echo "#Test success" >> build_script.bash
	echo "success=\`grep 'ScaLAPACK installation completed' scalapack.log | wc -l\`" >> build_script.bash
	echo  -e "if [ \044success -eq 1 ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo  \"Scalapack installation appears to have been successful.\"" >> build_script.bash
	echo >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo  \"ERROR(?): Scalapack installation does not seem to have been successful.\"" >> build_script.bash
	echo >> build_script.bash
	echo "fi" >> build_script.bash
	echo "echo  \"Check   \"`pwd`\"scalapack/scalapack_installer_0.94/scalapack.log for details.\"" >> build_script.bash
	echo "" >> build_script.bash
	echo "# Get handle to scalapack" >> build_script.bash
	echo "scalapack1=\`grep -A1 'Your ScaLAPACK library is:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "" >> build_script.bash
	echo -e "echo  \"Scalapack library: \044scalapack1 \" "  >> build_script.bash
	echo "" >> build_script.bash
	echo "# Handle to blasc" >> build_script.bash
	echo "blacs1=\`grep -A1 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "blacs2=\`grep -A2 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\`" >> build_script.bash
	echo "blacs3=\`grep -A3 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "" >> build_script.bash
	echo -e "echo  \"Blacs libraries: \044blacs1  \044blacs2  \044blacs3 \" "  >> build_script.bash
	echo " " >> build_script.bash
	echo "cd ../.." >> build_script.bash
	echo "" >> build_script.bash
	echo "# Read out the Fortran mangling from scalapack" >> build_script.bash
	printf "fortran_mangle=\`grep 'mangling' scalapack/scalapack_installer_*/scalapack.log | awk '{print \044NF;}'\`\n" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "# Make/install mumps" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ -d mumps ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo \"Info:   \"`pwd`\"/mumps already exists; overwriting\"" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    mkdir mumps" >> build_script.bash
	echo "    echo  \"Made mumps directory for installation \" " >> build_script.bash
	echo "fi" >> build_script.bash
	echo "cp "$(mumps_tar_file)" mumps" >> build_script.bash
	echo "cp build_scripts/Makefile.inc.template mumps" >> build_script.bash
	echo "cd mumps" >> build_script.bash
	echo "tar xfz "$(mumps_tar_file) >> build_script.bash
	echo "" >> build_script.bash
	echo "" >> build_script.bash
	echo "# Update Makefile" >> build_script.bash
	echo "escaped_blas_lib=\`echo $(blas_lib) | awk '{gsub(\"/\",\"\\\\\\/\"); print}' \`" >> build_script.bash
	echo -e "sed \"s/@mpicc@/"$(CC)"/g\" Makefile.inc.template  | sed \"s/@mpif90@/"$(F77)"/g\"  | sed \"s/@blas_lib@/\044escaped_blas_lib/g\" | sed \"s/@fopts@/"$(FFLAGS) -fallow-argument-mismatch"/g\" | sed \"s/@lopts@/"$(FFLAGS)"/g\" | sed \"s/@copts@/"$(CFLAGS)"/g\" | sed \"s/@mangle@/\044fortran_mangle/g\" > MUMPS_4.10.0/Makefile.inc" >> build_script.bash
	echo "" >> build_script.bash
	echo "cd MUMPS_4.10.0" >> build_script.bash
	echo "echo \"Building...\""
	echo "make d" >> build_script.bash
	echo "echo \"...done\" " >> build_script.bash
	echo "mumps1=\`pwd\`/lib/libdmumps.a" >> build_script.bash
	echo "mumps2=\`pwd\`/lib/libmumps_common.a" >> build_script.bash
	echo "mumps3=\`pwd\`/lib/libpord.a">> build_script.bash
	echo -e "echo  \"Mumps/pord libraries: \044mumps1  \044mumps2  \044mumps3 \" "  >> build_script.bash
	echo "cd .." >> build_script.bash
	echo "cd .." >> build_script.bash
	echo "rm -rf mumps_and_scalapack_default_installation"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation/lib"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation/include"  >> build_script.bash
	echo "cp mumps/MUMPS_4.10.0/lib/libdmumps.a mumps/MUMPS_4.10.0/lib/libmumps_common.a mumps/MUMPS_4.10.0/lib/libpord.a mumps_and_scalapack_default_installation/lib"  >> build_script.bash
	echo "cp scalapack/scalapack_installer_0.94/lib/libscalapack.a scalapack/scalapack_installer_0.94/lib/blacs.a scalapack/scalapack_installer_0.94/lib/blacsC.a scalapack/scalapack_installer_0.94/lib/blacsF77.a mumps_and_scalapack_default_installation/lib "  >> build_script.bash
	echo "cp mumps_and_scalapack_default_installation/lib/blacs.a mumps_and_scalapack_default_installation/lib/blacs_copy.a "  >> build_script.bash
	echo "cp mumps/MUMPS_4.10.0/lib/libdmumps.a mumps/MUMPS_4.10.0/include/*  mumps_and_scalapack_default_installation/include"  >> build_script.bash
	chmod a+x build_script.bash
	./build_script.bash

clean-local:
	rm -f build_script.bash
	rm -rf mumps
	rm -rf scalapack

distclean-local:
	rm -f build_script.bash
	rm -rf mumps
	rm -rf scalapack
	rm -rf mumps_and_scalapack_default_installation

junk:
	echo $(FFLAGS)

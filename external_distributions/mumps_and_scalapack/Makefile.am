EXTRA_DIST = \
build_scripts

# makefile variable to hold absolute address of top build directory
oomph_abs_top_srcdir=`cd $(top_builddir); pwd`

# Define mumps tar file
mumps_tar_file=MUMPS_4.8.4.tar.gz

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
	echo "#------------------------" >> build_script.bash
	echo "#Check settings" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ ! -e "$(lapack_lib)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: lapack library \n\n   "$(lapack_lib)" \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done! Found the lapack library \n\n   "$(lapack_lib)" \n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(blas_lib)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: blas library \n\n   "$(blas_lib)" \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done! Found the blas library \n\n   "$(blas_lib)" \n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(mpi_include)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: mpi include directory \n\n   "$(mpi_include)" \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done! Found the mpi include directory \n\n   "$(mpi_include)"\n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(mpi_include)"/mpi.h ] " >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: mpi header file \n\n   "$(mpi_include)"/mpi.h \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done! Found the mpi header file \n\n   "$(mpi_include)"/mpi.h\n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(scalapack_tar_file)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: scalapack installer tar file \n\n  "$(scalapack_tar_file)" \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done!  scalapack installer tar file \n\n  "$(scalapack_tar_file)" \n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "if [ ! -e "$(mumps_tar_file)" ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"ERROR: mumps tar file \n\n  "$(mumps_tar_file)" \n\ndoesn't exist\"" >> build_script.bash
	echo "    exit 1" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"Well done! mumps tar file \n\n  "$(mumps_tar_file)" \n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "# Install scalapack" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ -d scalapack ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"Info: \n\n   \"`pwd`\"/scalapack\n\nalready exists; overwriting\n\n\"" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    mkdir scalapack" >> build_script.bash
	echo "    echo -e \"Made scalapack directory for installation \n\n\" " >> build_script.bash
	echo "fi" >> build_script.bash
	echo "cp  "$(scalapack_tar_file)" scalapack" >> build_script.bash
	echo "cd scalapack" >> build_script.bash
	echo "tar xvfz "$(scalapack_tar_file) >> build_script.bash
	echo "cd  scalapack_installer_0.94" >> build_script.bash
	echo "" >> build_script.bash
	echo "#Install with blacs (downloaded on the fly)" >> build_script.bash
	echo "echo "Building..."" >> build_script.bash
	echo "./setup.py --downblacs --lapacklib=$(lapack_lib) --blaslib=$(blas_lib) --mpiincdir=$(mpi_include) > scalapack.log" >> build_script.bash
	echo "" >> build_script.bash
	echo "echo "...done!"" >> build_script.bash
	echo "#Test success" >> build_script.bash
	echo "success=\`grep 'ScaLAPACK installation completed' scalapack.log | wc -l\`" >> build_script.bash
	echo -e "if [ \044success -eq 1 ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"\n\nScalapack installation appears to have been successful.\n\n\"" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    echo -e \"\n\nERROR(?): Scalapack installation does not seem to have been successful.\n\n\"" >> build_script.bash
	echo "fi" >> build_script.bash
	echo "echo -e \"Check \n\n  \"`pwd`\"/scalapack.log \n\nfor details.\n\n\"" >> build_script.bash
	echo "" >> build_script.bash
	echo "# Get handle to scalapack" >> build_script.bash
	echo "scalapack1=\`grep -A1 'Your ScaLAPACK library is:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "" >> build_script.bash
	echo -e "echo \"Scalapack library: \044scalapack1 \" "  >> build_script.bash
	echo "" >> build_script.bash
	echo "# Handle to blasc" >> build_script.bash
	echo "blacs1=\`grep -A1 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "blacs2=\`grep -A2 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\`" >> build_script.bash
	echo "blacs3=\`grep -A3 'Your BLACS libraries are:' scalapack.log  | awk 'END{print $0}'\` " >> build_script.bash
	echo "" >> build_script.bash
	echo -e "echo -e \"Blacs libraries: \044blacs1  \044blacs2  \044blacs3 \" "  >> build_script.bash
	echo " " >> build_script.bash
	echo "cd ../.." >> build_script.bash
	echo "" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "# Make/install mumps" >> build_script.bash
	echo "#------------------------" >> build_script.bash
	echo "if [ -d mumps ]" >> build_script.bash
	echo "then" >> build_script.bash
	echo "    echo -e \"Info: \n\n   \"`pwd`\"/mumps\n\nalready exists; overwriting\n\n\"" >> build_script.bash
	echo "else" >> build_script.bash
	echo "    mkdir mumps" >> build_script.bash
	echo "    echo -e \"Made mumps directory for installation \n\n\" " >> build_script.bash
	echo "fi" >> build_script.bash
	echo -e "cp "$(mumps_tar_file)" mumps" >> build_script.bash
	echo "cp build_scripts/Makefile.inc.template mumps" >> build_script.bash
	echo "cd mumps" >> build_script.bash
	echo -e "tar xvfz "$(mumps_tar_file) >> build_script.bash
	echo "" >> build_script.bash
	echo "" >> build_script.bash
	echo "# Update Makefile" >> build_script.bash
	echo -e "escaped_blas_lib=\`echo \044blas_lib | awk '{gsub(\"/\",\"\0134\0134\0134/\"); print}' \`" >> build_script.bash
	echo -e "sed \"s/@mpicc@/"$(CC)"/g\" Makefile.inc.template  | sed \"s/@mpif90@/"$(F77)"/g\"  | sed \"s/@blas_lib@/\044escaped_blas_lib/g\" > MUMPS_4.8.4/Makefile.inc" >> build_script.bash
	echo "" >> build_script.bash
	echo "cd MUMPS_4.8.4" >> build_script.bash
	echo "echo \"Building...\""
	echo "make d" >> build_script.bash
	echo "echo \"...done\" " >> build_script.bash
	echo "mumps1=\`pwd\`/lib/libdmumps.a" >> build_script.bash
	echo "mumps2=\`pwd\`/lib/libmumps_common.a" >> build_script.bash
	echo "mumps3=\`pwd\`/lib/libpord.a">> build_script.bash
	echo -e "echo -e \"Mumps/pord libraries: \044mumps1  \044mumps2  \044mumps3 \" "  >> build_script.bash
	echo "cd .." >> build_script.bash
	echo "cd .." >> build_script.bash
	echo "rm -rf mumps_and_scalapack_default_installation"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation/lib"  >> build_script.bash
	echo "mkdir mumps_and_scalapack_default_installation/include"  >> build_script.bash
	echo "cp mumps/MUMPS_4.8.4/lib/libdmumps.a mumps/MUMPS_4.8.4/lib/libmumps_common.a mumps/MUMPS_4.8.4/lib/libpord.a mumps_and_scalapack_default_installation/lib"  >> build_script.bash
	echo "cp scalapack/scalapack_installer_0.94/lib/libscalapack.a scalapack/scalapack_installer_0.94/lib/blacs.a scalapack/scalapack_installer_0.94/lib/blacsC.a scalapack/scalapack_installer_0.94/lib/blacsF77.a mumps_and_scalapack_default_installation/lib "  >> build_script.bash
	echo "cp mumps_and_scalapack_default_installation/lib/blacs.a  mumps_and_scalapack_default_installation/lib/blacs_copy.a"  >> build_script.bash
	echo "cp mumps/MUMPS_4.8.4/lib/libdmumps.a mumps/MUMPS_4.8.4/include/*  mumps_and_scalapack_default_installation/include"  >> build_script.bash
	./build_script.bash 

clean-local:
	rm -f build_script.bash
	rm -rf mumps
	rm -rf scalapack


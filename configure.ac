#Process this file with autoconf to produce a configure script.

# Initialise autoconfigure: Pass package name, version number and contact email
# The M4 magic to automatically use the svn revision number is from
# http://caulfield.info/emmet/2008/05/get-the-subversion-revision-nu.html
# Slightly wrapped by script which checks if svnversion exists locally
# and, if it does, if we're working on revision controlled version.
# Useful because make mydist is often used to backup entire distribution
# so this should even work for non-developers. If we can't get a
# meaninful revision number for either reason we use 0000 instead.
define([svnversion], esyscmd([sh -c "bin/run_svnversion_if_it_exists.bash|tr -d '\n'"]))dnl
AC_INIT([oomph-lib],[1.0.svnversion],[oomph-lib@maths.manchester.ac.uk])

# config header file
AC_CONFIG_HEADERS(config.h)

# Intialise automake (AM_INIT_AUTOMAKE) using the command specified in this
# file. This is necessary so that we can select different options for
# different versions of automake (newer automake versions are not
# backwards compatible!).
m4_include([config/configure.ac_scripts/automake_init_command_file])

# Checks for programs:
#---------------------
# C++ compiler
AC_PROG_CXX

# C compiler
AC_PROG_CC

# C preprocessor
AC_PROG_CPP

# F77 compiler
AC_PROG_F77
AC_PROG_FC

# Install
AC_PROG_INSTALL

# Create symbolic links
AC_PROG_LN_S

# awk
AC_PROG_AWK

# Sets make macro
AC_PROG_MAKE_SET

# LIBTOOL
AC_PROG_LIBTOOL

# Determine linker flags for fortran
AC_F77_LIBRARY_LDFLAGS
if test -z $ac_cv_prog_f77_v; then
 echo "==================================================================="
 echo "Empty verbose linking string."
 echo "It is likely that the compilation has failed."
 echo "Either the Fortran compiler has not been found or it doesn't work"
 echo "If your Fortran compiler is in a non-standard location"
 echo "please specify by setting the variable F77=/path/to/fortran_compiler"
 echo "in your configuration file or as an environment variable."
 echo "The value of F77 found used here was $F77"
 echo "===================================================================="
 exit
fi

# Defined the no-optimisation environment flag for use in oomph_flapack
AC_ARG_VAR([FFLAGS_NO_OPT],[Fortran compiler flags without optimisation])

# Set additional flags for C++ compilation:
#------------------------------------------
# cfortran stuff:
case "$G77" in
  yes)
  #Test the version number by writing compiling and running
  #a small test C++ program
  #If the version number is different for the c-compiler and the
  #fortran compiler then we're in real trouble!

  cat > version_test.cc << EOF
#include<iostream>
int main() {if(__GNUC__ > 3) {std::cout << "yes";}else{std::cout << "no";}}
EOF

  $CXX version_test.cc
  gcc_version4=`./a.out`
  rm ./a.out version_test.cc
  if test "$gcc_version4" = yes; then
   echo "GNU Compilers > v 4.0, linking with -DgFortran"
   AM_CXXFLAGS='-DgFortran'
  else
   echo "GNU Compilers < v 4.0, linking with -Df2cFortran"
   AM_CXXFLAGS='-Df2cFortran'
  fi
  ;;

 no)
 echo "========================================================"
 echo "Not using GNU compiler suite"
 echo "You must specify a flag for cfortran.h"
 echo "in the CXXFLAGS of your config file."
 echo "e.g. -DINTEL_COMPILER for the Intel suite of compilers"
 echo "See cfortran.h documentation for further details"
 echo "========================================================"
 ;;
esac

#case "$F77 " in
# *gfortran*)
#  AM_CXXFLAGS='-DgFortran'
#  AM_FFLAGS=' '
#  ;;
# *ifort*)
#  AM_CXXFLAGS='-DINTEL_COMPILER'
#  AM_FFLAGS=' '
#  ;;
# *[[A-Za-z0-9]]*)
#  AM_CXXFLAGS='-Df2cFortran'
#  AM_FFLAGS=' '
#  ;;
# *)
#  echo "==================================================";
#  echo "Fortran compiler not found. ";
#  echo "If your Fortran compiler is in a non-standard location"
#  echo "please specify by setting the variable F77=/path/to/fortran_compiler"
#  echo "in your configuration file or as an environment variable."
#  echo "==================================================";
#  exit
#  ;;
#esac


AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_FFLAGS)


# Prepare variable that can be used to drag along C++ compilation
# flags
accumulated_cpp_flags=''



#...........................................................................
# Start of commands suggested by autoscan.
#
# Andrew: Will need to re-run autoscan on the
# full set of sources *and* decide on what to do if headers are not found.
# This is where the real work starts on making the code portable, I suppose.
# At the moment we simply bail out...


# Checks for header files.
#
# Note: old version was
#
#   AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h malloc.h],
#
# but malloc.h doesn't exist on BSD Linux. Since malloc.h was only suggested
# by autoscan, I've experimentally taken it out to see if this breaks
# anything under GNU< Linux (MH 12/01/2007)
#
AC_HEADER_STDC

#AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h],
#                  [echo "Congratulations -- required headers were found";],
#                  [echo "==================================================";
#                   echo "Trouble -- missing header[s]! ";
#                   echo "==================================================";
#                   exit])


#Check for the execinfo header
AC_CHECK_HEADER(execinfo.h,[have_execinfo=true;
                            echo "Can use stacktrace in debugging";],
                           [have_execinfo=false;
                            echo "Can't use stacktrace in debugging";])

# Check to see if we're on a Mac (the only system which has
# the OS name "Darwin")
AM_CONDITIONAL(HAVE_DARWIN, test $(uname) == "Darwin")

# Pass result of test to automake (makefiles can now check for
# status of HAVE_STACKTRACE in any Makefile.am:
AM_CONDITIONAL(HAVE_STACKTRACE, test x$have_execinfo = xtrue)


# Add flag (don't use the obvious HAVE_STACKTRACE as this is used by
# other packages too and leads to clashes!)
if test x$have_execinfo = xtrue; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_STACKTRACE"`; \
fi;


#Check for the unistd.h header (needed for getpid() in memory monitoring)
AC_CHECK_HEADER(unistd.h,[have_unistdh=true;
                            echo "Found unistd.h -- can do memory monitoring";],
                           [have_unistdh=false;
                            echo "Didn't find unistd.h -- cannot do memory monitoring";])

# Pass result of test to automake (makefiles can now check for
# status of HAVE_UNISTDH in any Makefile.am:
AM_CONDITIONAL(HAVE_UNISTDH, test x$have_unistdh = xtrue)

# Add flag (don't use the obvious HAVE_UNISTDH as this is used by
# other packages too and leads to clashes!)
if test x$have_unistdh = xtrue; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_UNISTDH"`; \
fi;




# Check for the fpu_control.h. Needed by triangle. If we don't have it
# use an empty replacement (included via symbolic link during build process
# and hope for the best...)
AC_CHECK_HEADER(fpu_control.h,[have_fpucontrolh=true;
                            echo "Found fpu_control.h -- can do triangle properly!";],
                           [have_fpucontrolh=false;
                            echo "Didn't find fpu_control.h -- using an empty replacement for triangle which may struggle... Keep an eye on it!";])

# Pass result of test to automake (makefiles can now check for
# status of HAVE_FPUCONTROLH in any Makefile.am:
AM_CONDITIONAL(HAVE_FPUCONTROLH, test x$have_fpucontrolh = xtrue)

# Add flag (don't use the obvious HAVE_FPUCONTROLH as this may be used by
# other packages too and may lead to clashes!)
if test x$have_fpucontrolh = xtrue; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_FPUCONTROLH"`; \
fi;





# Check for the malloc.h. Needed by metis. If we don't have it
# use an empty replacement
AC_CHECK_HEADER(malloc.h,[have_malloch=true;
                            echo "Found malloc.h!";],
                           [have_malloch=false;
                            echo "Didn't find malloc.h -- using an empty replacement metis... Keep an eye on it!";])


# Pass result of test to automake (makefiles can now check for
# status of HAVE_MALLOCH in any Makefile.am:
AM_CONDITIONAL(HAVE_MALLOCH, test x$have_malloch = xtrue)

# Add flag (don't use the obvious HAVE_MALLOCH as this may be used by
# other packages too and may lead to clashes!)
if test x$have_malloch = xtrue; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_MALLOCH"`; \
fi;





# Does C compiler support const?
AC_C_CONST

# Something to do with a stdbool.h header file. Don't quite understand
# the significance; see
# http://web.umr.edu/~gnudoc/single/autoconf-2.5.7/autoconf.html
# for full(er) description and suggested action.
AC_HEADER_STDBOOL

# Does the C compiler support inlining?
AC_C_INLINE

# If size_t is not defined, define size_t to be unsigned
AC_TYPE_SIZE_T

#If pid_t is not defined, define pid_t to be int.
AC_TYPE_PID_T

# Hmm, including this produces warnings during ./configure because
# RANLIB is made obsolete by LIBTOOL. You can either keep in here (in which
# case autoscan doesn't complain) or comment it out (in which case
# you should ignore the warnings in autoscan).
#AC_PROG_RANLIB

# Only needed by matlab part of superlu -- ignore when
# autoscan flags this up...
#AC_FUNC_ERROR_AT_LINE

# Something to do with the way Gnu C handles malloc -- see
# http://web.umr.edu/~gnudoc/single/autoconf-2.5.7/autoconf.html
# Only an issue in external C code (Superlu)
AC_FUNC_MALLOC

# autoscan suggests I should check for these functions. autoconf doesn't
# seem to be able to find them but this causes no problem --> ?
#AC_CHECK_FUNCS([floor pow gethrtime sqrt memchr memset],
#               [echo "Congratulations -- required function was found";],
#               [echo "==================================================";
#                echo "Trouble -- missing function! ";
#                echo "==================================================";
#                exit])



# END OF COMMANDS SUGGESTED BY autoscan
#...........................................................................


# Cray hack: Static C++ libraries
#--------------------------------
AC_ARG_ENABLE(static-cxx-libs,
              [  --enable-static-cxx-libs Pass -all-static to CXX when used as linker (needed when using gnu complier on Crays)],
        [want_static_cxx_libs=true],
        [want_static_cxx_libs=false])

if test x$want_static_cxx_libs = xtrue; then
        AC_SUBST([CXXLD], ['$(CXX) -all-static'])
        echo "Using -all-static flag for CXX linker"
else
        AC_SUBST([CXXLD], ['$(CXX)'])
        echo "Not using -all-static flag for CXX linker"
fi


# Do we have a sufficiently recent version of libtool installed?
#---------------------------------------------------------------
# Check to see if we're on a Mac (the only system which has
# the OS name "Darwin")
if test $(uname) == "Darwin"
then
    LIBTOOL_COMMAND=glibtool;
    # MacOS has problems compiling some external libraries because
    # of implicit function definitions
    accumulated_cpp_flags=`echo $accumulated_cpp_flags `" -Wno-implicit-function-declaration"
else
    LIBTOOL_COMMAND=libtool;
fi

# Convert libtool version number into single int -- strip dots and
# accomodate up to three digits for each part version id
libtool_version=`$LIBTOOL_COMMAND --version | \
   head -n 1 | $AWK  '{ print $NF}' | $AWK 'BEGIN { FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`

# Need at least version 242 --> 2004002
if test "$libtool_version" -ge 2004002; then
  echo "Congratulations!"
  echo "I have located a (sufficiently) up-to-date version of libtool ( "
  echo $libtool_version" ) and can build cgal and associated libraries."
  have_up_to_date_libtool=true
else
  echo " "; \
  echo "==============================================================="; \
  echo "Trouble: Your version of libtool is: " `$LIBTOOL_COMMAND --version | head -n 1 | $AWK  '{ print $NF}'`
  echo "Please upgrade to version 2.4.2 or higher."
  echo "I cannot/will not build cgal and associated libraries."
  echo "==============================================================="; \
  echo " "; \
  have_up_to_date_libtool=false
fi
AM_CONDITIONAL(HAVE_UP_TO_DATE_LIBTOOL, test x$have_up_to_date_libtool = xtrue)


# Do we have libtool at all?
AC_CHECK_PROG(have_libtool,$LIBTOOL_COMMAND,true,false)
AM_CONDITIONAL(HAVE_LIBTOOL, test x$have_libtool = xtrue)


# Do we have cmake installed?
#-----------------------------
AC_CHECK_PROG(have_cmake,cmake,true,false)
if test x$have_cmake = xfalse; then \
  echo " "; \
  echo "================================================================="; \
  echo "Cannot find cmake on your system"; \
  echo "---> Will not be able to build Trilinos > 9 and cgal etc."; \
  echo "If you want to build Trilinos 10 or higher "; \
  echo "please get your system administrator to install cmake for you. "; \
  echo " "; \
  echo "=================================================================="; \
  echo " "; \
else \
  echo "Congratulations! I have located cmake."; \
fi
AM_CONDITIONAL(HAVE_CMAKE, test x$have_cmake = xtrue)


# Do we have a sufficiently recent version of cmake installed?
#---------------------------------------------------------------
have_up_to_date_cmake=false
if test x$have_cmake = xtrue; then \

# Convert libtool version number into single int -- strip dots and
# accomodate up to three digits for each part version id
cmake_version=`cmake --version | \
   head -n 1 | $AWK  '{ print $NF}' | $AWK 'BEGIN { FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`
# Need at least version 2.8.11 --> 2008011
if test "$cmake_version" -ge 2008011; then
  echo "Congratulations!"
  echo "I have located a (sufficiently) up-to-date version of cmake ( "
  echo $cmake_version" ) and can build cgal and associated libraries."
  have_up_to_date_cmake=true
else
  echo " "; \
  echo "==============================================================="; \
  echo "Trouble: Your version of cmake is: " `cmake --version | head -n 1 | $AWK  '{ print $NF}'`
  echo "Please upgrade to version 2.8.1 or higher."
  echo "I cannot/will not build cgal and associated libraries."
  echo "==============================================================="; \
  echo " "; \
  have_up_to_date_cmake=false
fi


fi
AM_CONDITIONAL(HAVE_UP_TO_DATE_CMAKE, test x$have_up_to_date_cmake = xtrue)




# Do we have (a sufficiently recent version of) doxygen installed?
#-----------------------------------------------------------------
AC_CHECK_PROG(have_doxygen,doxygen,true,false)
if test x$have_doxygen = xfalse; then \
  echo " "; \
  echo "===================================================="; \
  echo "Cannot find doxygen on your system"; \
  echo "---> documentation cannot be generated. "; \
  echo " "; \
  echo "You can either obtain the (built) oomph-lib "; \
  echo "documentation online at "; \
  echo " "; \
  echo "          http://www.oomph-lib.org/    "; \
  echo " "; \
  echo "or download doxygen from  "; \
  echo " "; \
  echo "          http://www.doxygen.org/ "; \
  echo " "; \
  echo "===================================================="; \
  echo " "; \
else \
  # Convert doxygen version number into single int -- strip dots and
  # accomodate up to three digits for each part version id
  doxygen_version=`doxygen --version | \
    $AWK 'BEGIN { FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`
  # Need at least version 1.9.2 --> 1009002
  if test "$doxygen_version" -ge 1009002; then
    echo "Congratulations!"
    echo "I have located a (sufficiently) up-to-date version of doxygen:" \
         `doxygen --version`
    echo "and can build the documentation."
  else
    echo " "; \
    echo "==============================================================="; \
    echo "Trouble: Your version of doxygen is: " `doxygen --version`
    echo "Please upgrade to version 1.9.3 or higher."
    echo " "; \
    echo "          http://www.doxygen.org/ "; \
    echo " "; \
    echo "I cannot build the documentation..."
    echo " "; \
    echo "Meanwhile you can obtain the built oomph-lib documentation from"; \
    echo " "; \
    echo "          http://www.oomph-lib.org/    "; \
    echo " "; \
    echo "==============================================================="; \
    echo " "; \
    have_doxygen=false
  fi
fi
AM_CONDITIONAL(HAVE_DOXYGEN, test x$have_doxygen = xtrue)




# Do we have triangle installed online?
#---------------------------------------
AC_CHECK_PROG(have_online_triangle,triangle,true,false)
if test x$have_online_triangle = xfalse; then \
  echo " "; \
  echo "===================================================="; \
  echo "Cannot find triangle on your system"; \
  echo "---> cannot run fig2poly test. "; \
  echo " "; \
  echo "The source code for triangle may be obtained from "; \
  echo " "; \
  echo "     http://www.cs.cmu.edu/~quake/triangle.html    "; \
  echo " "; \
  echo "Please ask your systems administrator to install it "; \
  echo "on your machine. "; \
  echo " "; \
  echo "NOTE: Maybe you do have triangle installed in an  "; \
  echo "unusual place -- installing a symbolic link from"; \
  echo "a directory in the normal search path to your actual "; \
  echo "installation should then all that's required. "; \
  echo "===================================================="; \
  echo " "; \
else \
  echo "Congratulations! I have located triangle and can test the "; \
  echo "triangle-based mesh generation procedures."; \
fi
AM_CONDITIONAL(HAVE_ONLINE_TRIANGLE, test x$have_online_triangle = xtrue)


# Do we have tetgen installed?
#-------------------------------
AC_CHECK_PROG(have_tetgen,tetgen,true,false)
if test x$have_tetgen = xfalse; then \
  echo " "; \
  echo "===================================================="; \
  echo "Cannot find tetgen on your system"; \
  echo " "; \
  echo "The source code for tetgen may be obtained from "; \
  echo " "; \
  echo "     http://tetgen.berlios.de/    "; \
  echo " or "; \
  echo "     http://tetgen.org/    "; \
  echo " "; \
  echo "Please ask your systems administrator to install it "; \
  echo "on your machine. "; \
  echo " "; \
  echo "NOTE: Maybe you do have tetgen installed in an  "; \
  echo "unusual place -- installing a symbolic link from"; \
  echo "a directory in the normal search path to your actual "; \
  echo "installation should then all that's required. "; \
  echo "===================================================="; \
  echo " "; \
else \
  echo "Congratulations! I have located tetgen and can test the "; \
  echo "tetgen-based mesh generation procedures."; \
fi
AM_CONDITIONAL(HAVE_TETGEN, test x$have_tetgen = xtrue)



# Suppress build of pdf-based documentation?
# If yes, run configure as './configure --enable-suppress-pdf-doc'
AC_ARG_ENABLE(suppress-pdf-doc,
              [  --enable-suppress-pdf-doc Suppress building of oomph-lib's doxygen-based pdf documentation],
              [want_pdf_doc=false],
              [want_pdf_doc=true])

# Do we have pdflatex installed?
#-------------------------------
AC_CHECK_PROG(have_pdf_latex,pdflatex,true,false)
if test x$have_pdf_latex = xfalse; then \
  echo " "; \
  echo "===================================================="; \
  echo "Cannot find pdflatex on your system"; \
  echo "---> latex documentation cannot be generated. "; \
  echo " "; \
  echo "You can either lookup the documentation online at "; \
  echo " "; \
  echo "          http://www.oomph-lib.org/    "; \
  echo " "; \
  echo "or ask your systems administrator to install pdflatex "; \
  echo "on your machine. "; \
  echo " "; \
  echo "NOTE: Maybe you do have pdflatex installed in an  "; \
  echo "unusual place -- installing a symbolic link from"; \
  echo "a directory in the normal search path to your actual "; \
  echo "installation should then all that's required. "; \
  echo "===================================================="; \
  echo " "; \
else \
  echo "Congratulations! I have located pdflatex and can build the "; \
  echo "latex documentation in the doc directory."; \
  if test x$want_pdf_doc = xfalse; then \
     echo "BUILDING OF PDF DOC SUPPRESSED VIA CONFIGURE FLAG!"; \
     have_pdf_latex=false; \
  fi \
fi

echo "have_pdf_latex $have_pdf_latex"
AM_CONDITIONAL(HAVE_PDF_LATEX, test x$have_pdf_latex = xtrue)

# Do we have python installed?
#-----------------------------
AC_CHECK_PROG(have_python,python3,true,false)
if test x$have_python = xfalse; then \
  echo " "; \
  echo "================================================================="; \
  echo "Cannot find python on your system"; \
  echo "---> Will not be able to check the accuracy of "; \
  echo "     the self tests with bin/fpdiff.py :( "; \
  echo " "; \
  echo "Please get your system administrator to install python for you. "; \
  echo " "; \
  echo "=================================================================="; \
  echo " "; \
else \
  echo "Congratulations! I have located python."; \
fi
AM_CONDITIONAL(HAVE_PYTHON, test x$have_python = xtrue)



# Do we have validata?
#---------------------
AC_CHECK_FILES(demo_drivers/poisson/two_d_poisson_flux_bc_adapt/validata/results.dat.gz,
               echo "Found validata"; have_validata=true,
               echo "Not found validata"; have_validata=false)

if test x$have_validata = xfalse; then \
  echo " "; \
  echo "================================================================="; \
  echo "This distribution does not seem to have validata"; \
  echo "(based on existence of demo_drivers/poisson/two_d_poisson_flux_bc_adapt/validata/results.dat.gz)"; \
  echo "---> Will not be able to check the accuracy of "; \
  echo "     the self tests with bin/fpdiff.py :( "; \
  echo " "; \
  echo "=================================================================="; \
  echo " "; \
else \
  echo "This distribution does seem to have validata."; \
  echo "(based on existence of demo_drivers/poisson/two_d_poisson_flux_bc_adapt/validata/results.dat.gz)"; \
fi
AM_CONDITIONAL(HAVE_VALIDATA, test x$have_validata = xtrue)



# Can we compare against validata?
#---------------------------------
compare_against_validata=false;
if test x$have_python = xtrue; then \
if test x$have_validata = xtrue; then \
  compare_against_validata=true; \
fi \
fi

if test x$compare_against_validata = xtrue; then \
  echo "Will compare against validata "; \
else \
  echo "Will not compare against validata "; \
fi

AM_CONDITIONAL(COMPARE_AGAINST_VALIDATA, test x$compare_against_validata = xtrue)




# Do we have the sources for the HSL frontal solver in external_src?
#-------------------------------------------------------------------
AC_CHECK_FILES(external_src/oomph_hsl/frontal.f,
               echo "Found external_src/oomph_hsl/frontal.f"; have_hsl_sources=true,
               echo "Not found external_src/oomph_hsl/frontal.f"; have_hsl_sources=false)


#If we don't have them check if they are private/external_src and make sym link
#------------------------------------------------------------------------------
if test x$have_hsl_sources = xfalse; then \

  AC_CHECK_FILES(private/external_src/hsl/frontal.f,
                 echo "Found private/external_src/hsl/frontal.f"; have_hsl_sources_in_private_external_src=true,
                 echo "Not found private/external_src/hsl/frontal.f"; have_hsl_sources_in_private_external_src=false)

  if test x$have_hsl_sources_in_private_external_src = xtrue; then \
    echo "Making symlink"; \
    ln -sf ../../private/external_src/hsl/frontal.f external_src/oomph_hsl/frontal.f;\
    have_hsl_sources=true;\
  else \
    echo "Not making symlink";\
  fi; \

fi;

# Now do the final check if we have the sources for the HSL frontal solver
#-------------------------------------------------------------------------
AM_CONDITIONAL(HAVE_HSL_SOURCES, test x$have_hsl_sources = xtrue)


# Build flags that can be set during configure stage
#---------------------------------------------------

# Dummies to produce better output for ./configure --help
AC_ARG_ENABLE(dummy1,[  ])
AC_ARG_ENABLE(dummy2,[Oomph-lib options:  ])







# Is MPI build required?
# If yes, run configure as './configure --enable-MPI'
AC_ARG_ENABLE(MPI,
              [  --enable-MPI Build oomph-lib MPI routines],
              [want_mpi=true],
              [want_mpi=false])
# Pass result of test to automake (makefiles can now check for
# status of WANT_MPI in any Makefile.am, like this:
#
#   # WANT_MPI flag is set by autoconf
#   #---------------------------------
#   if WANT_MPI
#     SUBDIRS = $(non_mpi_sub_dirs) $(mpi_sub_dirs)
#   else
#     SUBDIRS = $(non_mpi_sub_dirs)
#   endif
#
AM_CONDITIONAL(WANT_MPI, test x$want_mpi = xtrue)


# Add flag (don't use the obvious HAVE_MPI as this is used by
# other packages too and leads to clashes!)
if test x$want_mpi = xtrue; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_MPI"`; \
fi;




# Do we want to run the gmsh tests?
AC_ARG_WITH(gmsh-self-tests,
[  --with-gmsh-self-tests[=GMSH_RUN_COMMAND] run gmsh self-tests.  Requires
                          specification of the string required to run gmsh
                          from the command line. GMSH_RUN_COMMAND. E.g. if
                          gmsh on your machine needs to be run with

                             /usr/bin/gmsh

                          then set GMSH_RUN_COMMAND to "/usr/bin/gmsh"],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    run_gmsh_self_tests=true
    GMSH_RUN_COMMAND=$withval
  fi
])
AC_SUBST(GMSH_RUN_COMMAND)




# Do we want to run the mpi tests?
AC_ARG_WITH(mpi-self-tests,
[  --with-mpi-self-tests[=MPI_RUN_COMMAND] run mpi self-tests on two processors
                          using MPI_RUN_COMMAND. E.g. if a parallel executable
                          on your machine needs to be run with

                             mpirun -np 2 a.out

                          then set MPI_RUN_COMMAND to "mpirun -np 2"],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    run_mpi_self_tests=true
    MPI_RUN_COMMAND=$withval
  fi
])
AC_SUBST(MPI_RUN_COMMAND)


# Do we want to run the mpi self tests which involves more than two cores?
AC_ARG_WITH(mpi-self-tests-variablenp,
[  --with-mpi-self-tests-variablenp[=MPI_VARIABLENP_RUN_COMMAND] run mpi self-tests on multiple processors
                          using MPI_VARIABLENP_RUN_COMMAND. E.g. if a parallel executable
                          on your machine needs to be run with

                             mpirun -np 2 a.out

                          then set MPI_VARIABLENP_RUN_COMMAND to "mpirun -np OOMPHNP"],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    run_mpi_self_tests_variablenp=true
    MPI_VARIABLENP_RUN_COMMAND=$withval
  fi
],[run_mpi_self_tests_variablenp=false])

# Perform the substitution
AC_SUBST(MPI_VARIABLENP_RUN_COMMAND)

#   # WANT_MPI_VARIABLENP flag is set by autoconf
#   #---------------------------------
#   if WANT_MPI_VARIABLENP
#     SUBDIRS = $(non_mpi_variablenp_sub_dirs) $(mpi_variablenp_sub_dirs)
#   else
#     SUBDIRS = $(non_mpi_variablenp_sub_dirs)
#   endif
AM_CONDITIONAL(WANT_MPI_VARIABLENP, test x$run_mpi_self_tests_variablenp = xtrue)


# Do we want to build Trilinos dynamically?
 AC_ARG_ENABLE(dynamic_linking_for_trilinos,
               [  --enable-dynamic-linking-for-trilinos Build trilinos as a
 		   					shared library instead
							of a static library],
               [dynamic_linking_for_trilinos_requested=true],
               [dynamic_linking_for_trilinos_requested=false])
  

# Pass result of test to automake (makefiles can now check for
# status of WANT_DYNAMIC_BUILD_OF_TRILINOS in any Makefile.am)
AM_CONDITIONAL(WANT_DYNAMIC_BUILD_OF_TRILINOS, test x$dynamic_linking_for_trilinos_requested = xtrue)


# Do we need multiple teuchos sub-libraries
AC_ARG_ENABLE(multiple_teuchos_libraries,
              [  --enable-multiple_teuchos_libraries Assume that the
                                                     version of trilinos
                                                     provides multiple
                                                     teuchos (sub-)libraries],
              [multiple_teuchos_libs_requested=true],
              [multiple_teuchos_libs_requested=false])


# Suppress build of extensive documentation?
# If yes, run configure as './configure --enable-suppress-doc'
AC_ARG_ENABLE(suppress-doc,
              [  --enable-suppress-doc Suppress building of oomph-lib's doxygen-based documentation],
              [want_doc=false],
              [want_doc=true])

#If the doc directory does not exist, then you can't have it
AC_CHECK_FILE(doc,[],[want_doc=false])

# Pass result of test to automake (makefiles can now check for
# status of WANT_DOC in any Makefile.am
AM_CONDITIONAL(WANT_DOC, test x$want_doc = xtrue)


# Suppress extended self-tests?
# If yes, run configure as './configure --enable-suppress-extended-self-tests'
AC_ARG_ENABLE(suppress-extended-self_tests,
              [  --enable-suppress-extended-self-tests Suppress build/running of oomph-lib's extended self tests for all libraries],
              [want_extended_self_tests=false],
              [want_extended_self_tests=true])
# Pass result of test to automake (makefiles can now check for
# status of WANT_DOC in any Makefile.am
AM_CONDITIONAL(WANT_EXTENDED_SELF_TESTS, test x$want_extended_self_tests = xtrue)




# Suppress build of extensive demo codes?
# If yes, run configure as './configure --enable-suppress-demo'
AC_ARG_ENABLE(suppress-demo,
              [  --enable-suppress-demo Suppress build/check of oomph-lib's demo codes],
              [want_demo=false],
              [want_demo=true])
# Pass result of test to automake (makefiles can now check for
# status of WANT_DEMO in any Makefile.am
AM_CONDITIONAL(WANT_DEMO, test x$want_demo = xtrue)



# Suppress build of user sources/libraries?
# If yes, run configure as './configure --enable-suppress-user-src'
AC_ARG_ENABLE(suppress-user-src,
              [  --enable-suppress-user-src Suppress build/check of user sources/libraries],
              [want_user_src=false],
              [want_user_src=true])
# Pass result of test to automake (makefiles can now check for
# status of WANT_USER_SRC in any Makefile.am
AM_CONDITIONAL(WANT_USER_SRC, test x$want_user_src = xtrue)



# Replace include files by links?
# If yes, run configure as './configure --enable-symbolic_links_for_headers'
AC_ARG_ENABLE(symbolic-links-for-headers,
              [  --enable-symbolic-links-for-headers Use symbolic links for header files; useful during development],
              [want_sym_link=true],
              [want_sym_link=false])
# Pass result of test to automake (makefiles can now check for
# status of SYMBOLIC_LINKS_FOR_HEADERS in any Makefile.am
AM_CONDITIONAL(SYMBOLIC_LINKS_FOR_HEADERS, test x$want_sym_link = xtrue)



# Find user_src subdirectories
USER_SRC_SUBDIRECTORIES=`cd user_src/ && echo */`
AC_SUBST(USER_SRC_SUBDIRECTORIES)

# Find user_drivers subdirectories
USER_DRIVERS_SUBDIRECTORIES=`cd user_drivers/ && echo */`
AC_SUBST(USER_DRIVERS_SUBDIRECTORIES)


# Link gfortran libraries when building oomph_flapack library
# This is a workaround to a problem with the latest versions of
# gfortran (ubuntu13.10 -- seems like a bug to use so we're not
# hardcoding it here..
# Run configure with './configure --enable_use_gfortran_libraries_for_oomph_flapack'
AC_ARG_ENABLE(use-gfortran-libraries-for-oomph-flapack,
              [  --enable-use-gfortran-libraries-for-oomph-flapack Link gfortran libraries when building oomph_flapack library (workaround for ubuntu 13.10 feature(bug?)],
              [use_gfortran_libraries_for_oomph_flapack=true],
              [use_gfortran_libraries_for_oomph_flapack=false])

# Pass result of test to automake (makefiles can now check for
# status of ADD_LGFORTRAN_TO_OOMPH_LAPACK in any Makefile.am
AM_CONDITIONAL(ADD_LGFORTRAN_TO_OOMPH_LAPACK, test x$use_gfortran_libraries_for_oomph_flapack = xtrue)




#----------------------
# CGAL stuff
#----------------------

# Suppress build of cgal and underlying third party distributions?
# If yes, run configure as './configure --enable-suppress-cgal-build'
AC_ARG_ENABLE(suppress-cgal-build,
              [  --enable-suppress-cgal-build Suppress building of cgal and underlying third party distributions],
              [want_cgal_build=false],
              [want_cgal_build=true])

if  test x$have_up_to_date_libtool = xfalse; then
    want_cgal_build=false
fi
if  test x$have_up_to_date_cmake = xfalse; then
    want_cgal_build=false
fi

# Pass result of test to automake (makefiles can now check for
# status of WANT_CGAL in any Makefile.am)
AM_CONDITIONAL(WANT_CGAL_BUILD, test x$want_cgal_build = xtrue)


# Do we want to install cgal and associated third party distributions
# in permament location?
AC_ARG_WITH(cgal-permanent-installation-dir,
[  --with-cgal-permanent-installation-dir[=ABSOLUTE_PATH_TO_PERMAMENT_INSTALL_DIRECTORY] Install cgal and associated third party distributions permamently
                          in specified directory.],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    install_cgal_in_permanent_dir=true
    CGAL_PERMANENT_INSTALLATION_DIR=$withval
  fi
])
AC_SUBST(CGAL_PERMANENT_INSTALLATION_DIR)

# Pass result of test to automake (makefiles can now check for
# status of CGAL_PERMANENT_INSTALLATION_DIR_SPECIFIED in any Makefile.am)
AM_CONDITIONAL(CGAL_PERMANENT_INSTALLATION_DIR_SPECIFIED, test x$install_cgal_in_permanent_dir = xtrue)

if  test x$have_up_to_date_libtool = xtrue; then
if  test x$have_up_to_date_cmake = xtrue; then
if test x$install_cgal_in_permanent_dir = xtrue; then
if test -e $CGAL_PERMANENT_INSTALLATION_DIR; then
   echo "Permanent install directory for CGAL and associated third party"
   echo "distributions"
   echo " "
   echo "   "$CGAL_PERMANENT_INSTALLATION_DIR
   echo " "
   echo "already exists."
else
   mkdir $CGAL_PERMANENT_INSTALLATION_DIR
   echo " "
   echo "Created permanent install directory for CGAL and associated"
   echo "third party distributions:"
   echo " "
   echo "   "$CGAL_PERMANENT_INSTALLATION_DIR
   echo " "
fi
fi
fi
fi

#----------------------
# End CGAL stuff
#----------------------






# Auto-detect directories in private/*
# -----------------------------------------------------------------------

# This information is used in "Makefile.am"s within private/ to include
# directories from private/ when running `make dist`. In particular it has
# nothing to do with the creation of Makefile s by configure.


# Do we have private directory?
AC_CHECK_FILES(private,
               echo "Found private"; have_private_directory=true,
               have_private_directory=false)
AM_CONDITIONAL(HAVE_PRIVATE_DIRECTORY, test x$have_private_directory = xtrue)


# Define subdirectories in private directory
if test x$have_private_directory = xtrue; then  \

  PRIVATE_DIR="private"; \

  # A hack to exclude the private/external_src directory from DIST_SUBDIRS
  # (it's handled by EXTRA_DIST)
  PRIVATE_SUBDIRECTORIES=`cd private/ && echo */ | tr ' ' '\n' | grep -v 'external_src' | tr '\n' ' '`; \
  echo "Do have private subdirectories" $PRIVATE_SUBDIRECTORIES; \

else \

  PRIVATE_DIR="";\
  PRIVATE_SUBDIRECTORIES=""; \
  echo "Don't have private subdirectories: " $PRIVATE_SUBDIRECTORIES; \

fi;

AC_SUBST(PRIVATE_DIR)
AC_SUBST(PRIVATE_SUBDIRECTORIES)

# Do we have private user sources?
AC_CHECK_FILES(private/user_src,
               echo "Found private/user_src"; have_private_user_sources=true,
               have_private_user_sources=false)
AM_CONDITIONAL(HAVE_PRIVATE_USER_SRC, test x$have_private_user_sources = xtrue)


# Define subdirectories in private user src directory
if test x$have_private_user_sources = xtrue; then  \

  PRIVATE_USER_SRC_SUBDIRECTORIES=`cd private/user_src/ && echo */`
  echo "Do have private/user_src: " $PRIVATE_USER_SRC_SUBDIRECTORIES; \

else \

  PRIVATE_USER_SRC_SUBDIRECTORIES=""; \
  echo "Don't have private/user_src: " $PRIVATE_USER_SRC_SUBDIRECTORIES; \

fi;

AC_SUBST(PRIVATE_USER_SRC_SUBDIRECTORIES)

# Do we have private user drivers?
AC_CHECK_FILES(private/user_drivers,
               echo "Found private/user_drivers"; have_private_user_drivers=true,
               have_private_user_drivers=false)
AM_CONDITIONAL(HAVE_PRIVATE_USER_DRIVERS, test x$have_private_user_drivers = xtrue)


# Define subdirectories in private user drivers directory
if test x$have_private_user_drivers = xtrue; then  \

  PRIVATE_USER_DRIVER_SUBDIRECTORIES=`cd private/user_drivers/ && echo */`
  echo "Do have private/user_drivers: " $PRIVATE_USER_DRIVER_SUBDIRECTORIES; \

else \

  PRIVATE_USER_DRIVER_SUBDIRECTORIES=""; \
  echo "Don't have private/user_drivers: " $PRIVATE_USER_DRIVER_SUBDIRECTORIES; \

fi;

AC_SUBST(PRIVATE_USER_DRIVER_SUBDIRECTORIES)



# Initialiase flags that contains all external libraries that come with
#----------------------------------------------------------------------
# their own distribution and may be installed anywhere.
#------------------------------------------------------
# These are distributed to all Makefile.am-s via the AM_CPPFLAGS
#---------------------------------------------------------------
# and AM_LDFLAGS variables so *everything* links against them
#------------------------------------------------------------
external_dist_libraries=''
external_dist_libraries_lib=''
external_dist_libraries_include=''




#................................................................
#................................................................



#---------------------------------------------------
# Potential location for the blas libraries.
#
# If the complete blas library is available on
# your machine, specify the absolute path to the
# library with --with-blas
#----------------------------------------------------
BLAS_LIB=
check_for_blas=false

AC_ARG_WITH(blas,
[  --with-blas[=ABSOLUTE_PATH] use blas library specified via ABSOLUTE PATH
                               e.g. /usr/local/lib/blas.a
                               or specify "auto" (without quotation marks!)
                               if the library is specified automatically
                               by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_blas=true
    BLAS_LIB=$withval
  fi
])



# Do they exist?
have_blas=false
if test x$check_for_blas = xtrue; then \
if test -e $BLAS_LIB; then \
    have_blas=true; \
    echo "Found the blas library "; \
    echo " "; \
    echo "   "  $BLAS_LIB ; \
    echo " "; \
    # Make the variable available to all Makefile.am-s
    AC_SUBST(BLAS_LIB)
else \
   if test "$BLAS_LIB" == "auto"; then
       have_blas=true; \
       echo "Will assume that blas library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
       have_blas=false; \
       echo "Did NOT find the blas library"; \
       echo " "; \
       echo "   "  $BLAS_LIB ; \
       echo " "; \
    fi; \
fi; \
else \
    echo "No blas library specified"; \
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_BLAS, test x$have_blas = xtrue)




#................................................................
#................................................................



#---------------------------------------------------
# Potential location for the lapack libraries.
#
# If the complete lapack library is available on
# your machine, specify the absolute path to the
# library with --with-lapack
#----------------------------------------------------
LAPACK_LIB=
check_for_lapack=false

AC_ARG_WITH(lapack,
[  --with-lapack[=ABSOLUTE_PATH] use lapack library specified via ABSOLUTE PATH
                               e.g. /usr/local/lib/lapack.a
                               or specify "auto" (without quotation marks!)
                               if the library is specified automatically
                               by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_lapack=true
    LAPACK_LIB=$withval
  fi
])



# Do they exist?
have_lapack=false
if test x$check_for_lapack = xtrue; then \
if test -e $LAPACK_LIB; then \
    have_lapack=true; \
    echo "Found the lapack library "; \
    echo " "; \
    echo "   "  $LAPACK_LIB ; \
    echo " "; \
    # Make the variable available to all Makefile.am-s
    AC_SUBST(LAPACK_LIB)
else \
   if test "$LAPACK_LIB" == "auto"; then
       have_lapack=true; \
       echo "Will assume that lapack library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
       have_lapack=false; \
       echo "Did NOT find the lapack library"; \
       echo " "; \
       echo "   "  $LAPACK_LIB ; \
       echo " "; \
    fi \
fi; \
else \
    echo "No lapack library specified"; \
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_LAPACK, test x$have_lapack = xtrue)



#................................................................
#................................................................



#----------------------------------------------------------------
# Allow suppression of oomph-lib's triangle build (there are some subtle
# comiler problems on some machines; if they exist we'd rather bin
# the third-party code rather than endangering the entire oomph-lib
# build
#----------------------------------------------------------------
have_triangle_lib=false
AC_ARG_ENABLE(suppress-triangle,
              [  --enable-suppress-triangle Supress build/use of
                                 oomph-lib's own copy of the triangle library. ],
        [suppress_triangle=true],
        [suppress_triangle=false])

if test x$suppress_triangle = xtrue; then
        echo "Suppressing build of oomph-lib's own copy of the triangle library."
        echo "This disables a lot of useful functionality but may be"
        echo "necessary where the compiler struggles to compile the"
        echo "triangle sources. In that case, triangle should be "
        echo "installed separately (by whatever method it takes) "
        echo "and included in the linker path (similar to the way we handle "
        echo "blas and lapack, say), but this hasn't been implemented"
        echo "in our build machinery. Please contact the developers if you"
        echo "need this..."
else
        echo "Building/using oomph-lib's own copy of the triangle library"
fi

#Tell everybody about it...
AM_CONDITIONAL(SUPPRESS_TRIANGLE_LIB, test x$suppress_triangle = xtrue)


# Add flag
if test x$suppress_triangle = xfalse; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_TRIANGLE_LIB"`; \
fi;


#................................................................



#----------------------------------------------------------------
# Allow suppression of oomph-lib's tetgen build (there may be subtle
# comiler problems on some machines; if they exist we'd rather bin
# the third-party code rather than endangering the entire oomph-lib
# build
#----------------------------------------------------------------
have_tetgen_lib=false
AC_ARG_ENABLE(suppress-tetgen,
              [  --enable-suppress-tetgen Supress build/use of
                                 oomph-lib's own copy of the tetgen library. ],
        [suppress_tetgen=true],
        [suppress_tetgen=false])

if test x$suppress_tetgen = xtrue; then
        echo "Suppressing build of oomph-lib's own copy of the tetgen library."
        echo "This disables a lot of useful functionality but may be"
        echo "necessary where the compiler struggles to compile the"
        echo "tetgen sources. In that case, tetgen should be "
        echo "installed separately (by whatever method it takes) "
        echo "and included in the linker path (similar to the way we handle "
        echo "blas and lapack, say), but this hasn't been implemented"
        echo "in our build machinery. Please contact the developers if you"
        echo "need this..."
else
        echo "Building/using oomph-lib's own copy of the tetgen library"
fi

#Tell everybody about it...
AM_CONDITIONAL(SUPPRESS_TETGEN_LIB, test x$suppress_tetgen = xtrue)


# Add flag
if test x$suppress_tetgen = xfalse; then \
    accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DOOMPH_HAS_TETGEN_LIB"`; \
fi;



#----------------------------------------------------------------
# Allow automatic specification of metis by compiler/linker
# (e.g. for Hector)
#----------------------------------------------------------------
have_metis=false
AC_ARG_ENABLE(suppress-oomph-metis,
              [  --enable-suppress-oomph-metis Supress building of oomph-lib's
                                 own copy of metis. The library
                                 is then assumed to be provided automatically
                                 by the compiler/linker.],
        [have_metis=true],
        [have_metis=false])

if test x$have_metis = xtrue; then
        echo "Assuming that metis is provided automatically by compiler/linker"
else
        echo "Building/using oomph-lib's own copy of metis"
fi

#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_METIS, test x$have_metis = xtrue)



#----------------------------------------------------------------
# Allow automatic specification of parmetis by compiler/linker
# (e.g. for Hector)
#----------------------------------------------------------------
have_parmetis=false
AC_ARG_ENABLE(suppress-oomph-parmetis,
              [  --enable-suppress-oomph-parmetis Supress building of oomph-lib's
                                 own copy of parmetis. The library
                                 is then assumed to be provided automatically
                                 by the compiler/linker.],
        [have_parmetis=true],
        [have_parmetis=false])

if test x$have_parmetis = xtrue; then
        echo "Assuming that parmetis is provided automatically by compiler/linker"
else
        echo "Building/using oomph-lib's own copy of parmetis"
fi

#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_PARMETIS, test x$have_parmetis = xtrue)


#----------------------------------------------------------------
# Allow automatic specification of superlu by compiler/linker
# (e.g. for Hector)
#----------------------------------------------------------------
have_superlu=false
AC_ARG_ENABLE(suppress-oomph-superlu,
              [  --enable-suppress-oomph-superlu Supress building of oomph-lib's
                                 own copy of serial SuperLU. The library
                                 is then assumed to be provided automatically
                                 by the compiler/linker.],
        [have_superlu=true],
        [have_superlu=false])

if test x$have_superlu = xtrue; then
        echo "Assuming that SuperLU is provided automatically by compiler/linker"
else
        echo "Building/using oomph-lib's own copy of SuperLU"
        accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DUSING_OOMPH_SUPERLU"`; \
fi

#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_SUPERLU, test x$have_superlu = xtrue)

#----------------------------------------------------------------
# Allow automatic specification of superlu_dist by compiler/linker
# (e.g. for Hector)
#----------------------------------------------------------------
have_superlu_dist=false
AC_ARG_ENABLE(suppress-oomph-superlu_dist,
              [  --enable-suppress-oomph-superlu_dist Supress building of oomph-lib's
                                 own copy of serial SuperLU_dist. The library
                                 is then assumed to be provided automatically
                                 by the compiler/linker.],
        [have_superlu_dist=true],
        [have_superlu_dist=false])

if test x$have_superlu_dist = xtrue; then
        echo "Assuming that SuperLU_dist is provided automatically by compiler/linker"
else
        echo "Building/using oomph-lib's own copy of SuperLU_dist"
        accumulated_cpp_flags=`echo $accumulated_cpp_flags " -DUSING_OOMPH_SUPERLU_DIST"`; \
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_SUPERLU_DIST, test x$have_superlu_dist = xtrue)




#................................................................
#................................................................



#---------------------------------------------------
# Potential location of HYPRE libraries.
#
# If the HYPRE is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-hypre
#----------------------------------------------------
HYPRE_DIR=
HYPRE_LIB=
HYPRE_INCLUDE=
check_for_hypre=false

AC_ARG_WITH(hypre,
[  --with-hypre[=DIR] use hypre library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_hypre=true
    HYPRE_DIR=$withval
    HYPRE_LIB="$withval/lib"
    HYPRE_INCLUDE="$withval/include"
  fi
])

# Has it already been built?
hypre_has_already_been_built=false;

# Do they exist?
have_hypre=false
if test x$check_for_hypre = xtrue; then \
if test -d $HYPRE_LIB; then \
    have_hypre=true; \
    echo "Found the hypre libraries in"; \
    echo " "; \
    echo "   "  $HYPRE_LIB ; \
    echo " "; \
else \
   if test "$HYPRE_DIR" == "auto"; then \
       have_hypre=true; \
       echo "Will assume that hypre library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_hypre=false; \
       echo "Did NOT find the hypre libraries in"; \
       echo " "; \
       echo "   "  $HYPRE_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $HYPRE_INCLUDE; then \
    have_hypre=true; \
    echo "Found the hypre include files in"; \
    echo " "; \
    echo "   "  $HYPRE_INCLUDE ; \
    echo " "; \
else \
   if test "$HYPRE_DIR" == "auto"; then \
       have_hypre=true; \
       echo "Will assume that hypre library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
      have_hypre=false; \
      echo "Did NOT find the hypre include files in"; \
      echo " "; \
      echo "   "  $HYPRE_INCLUDE ; \
      echo " "; \
    fi; \
fi; \
hypre_has_already_been_built=$have_hypre; \
else \
    echo "No hypre library directory specified"; \

    # Will hypre library be built?
    #-----------------------------
    AC_CHECK_FILES(external_distributions/hypre/hypre-2.0.0.tar.gz,
                   echo "Found hypre-2.0.0.tar.gz --> will build hypre library"; will_have_hypre=true,
                   echo "Not found hypre-2.0.0.tar.gz --> will not build hypre library"; will_have_hypre=false)

   if test x$will_have_hypre = xtrue; then \
    have_hypre=true; \
    HYPRE_DIR="`pwd`/external_distributions/hypre/hypre_default_installation"
    HYPRE_LIB="`pwd`/external_distributions/hypre/hypre_default_installation/lib"
    HYPRE_INCLUDE="`pwd`/external_distributions/hypre/hypre_default_installation/include"
   fi;
fi



# Do we have the hypre libraries? If so add them to the
#------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_hypre = xfalse; then
  HYPRE_LIB='' ; \
  HYPRE_INCLUDE='' ;\
else \
  if test "$HYPRE_DIR" == "auto"; then \
    HYPRE_LIB='' ; \
    HYPRE_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_HYPRE"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -lHYPRE"` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$HYPRE_LIB -Wl,-rpath,$HYPRE_LIB` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$HYPRE_INCLUDE -DOOMPH_HAS_HYPRE` ; \
  fi
fi

# Indicate that Trilinos has already been built (so we don't bother trying# to install        
# it again if there's a tar file available!)                                                  
AM_CONDITIONAL(ALREADY_HAVE_BUILT_HYPRE, test x$hypre_has_already_been_built = xtrue)  

#Tell everybody about it...
AM_CONDITIONAL(OOMPH_HAS_HYPRE, test x$have_hypre = xtrue)

# Make thes variables available to all Makefile.am-s
AC_SUBST(HYPRE_LIB)
AC_SUBST(HYPRE_INCLUDE)


#................................................................




#---------------------------------------------------
# Potential location of the TRILINOS libraries.
#
# If the TRILINOS is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-trilinos
#----------------------------------------------------
TRILINOS_DIR=
TRILINOS_LIB=
TRILINOS_INCLUDE=
check_for_trilinos=false

AC_ARG_WITH(trilinos,
[  --with-trilinos[=DIR] use trilinos library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_trilinos=true
    TRILINOS_DIR=$withval
    TRILINOS_LIB="$withval/lib"
    TRILINOS_INCLUDE="$withval/include"
  fi
])

# Has Trilinos already been built? Assume it hasn't...
trilinos_has_already_been_built=false;

# Do they exist?
have_trilinos=false
if test x$check_for_trilinos = xtrue; then \
if test -d $TRILINOS_LIB; then \
    have_trilinos=true; \
    trilinos_has_already_been_built=true; \
    echo "Found the trilinos libraries in"; \
    echo " "; \
    echo "   "  $TRILINOS_LIB ; \
    echo " "; \
else \
   if test "$TRILINOS_DIR" == "auto"; then \
       have_trilinos=true; \
       trilinos_has_already_been_built=true; \
       echo "Will assume that trilinos library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_trilinos=false; \
       trilinos_has_already_been_built=false; \
       echo "Did NOT find the trilinos libraries in"; \
       echo " "; \
       echo "   "  $TRILINOS_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $TRILINOS_INCLUDE; then \
    have_trilinos=true; \
    trilinos_has_already_been_built=true; \
    echo "Found the trilinos include files in"; \
    echo " "; \
    echo "   "  $TRILINOS_INCLUDE ; \
    echo " "; \
else \
   if test "$TRILINOS_DIR" == "auto"; then \
       have_trilinos=true; \
       trilinos_has_already_been_built=true; \
       echo "Will assume that trilinos library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
       have_trilinos=false; \
       trilinos_has_already_been_built=false; \
       echo "Did NOT find the trilinos include files in"; \
       echo " "; \
       echo "   "  $TRILINOS_INCLUDE ; \
       echo " "; \
   fi;\
fi; \
echo ""; \
echo ""; \
echo "========================================================================"; \
echo ""; \
echo "NOTE: We're going to link against an existing trilinos installation."; \
echo "      Unfortunately different versions of that library require"; \
echo "      explicit linking against different sets of (sub-)libraries."; \
echo "      If your version of Trilinos is greater than 10.8.5 then we're"; \
echo "      going to assume that you have multiple teuchos libraries, not just"; \
echo "      a single teuchos library, i.e. -lteuchos. This is the earliest";
echo "      version that we know of that has multiple teuchos libraries. If "; \
echo "      your version of Trilinos is older than that but still has multiple"; \
echo "      teuchos then you may still have linking problems that indicate"; \
echo "      missing teuchos functions. In this case, try to re-install oomph-lib"; \
echo "      with the"; \
echo ""; \
echo "                   --enable_multiple_teuchos_libraries "; \
echo ""; \
echo "      configuration flag."; \
echo ""; \
echo "========================================================================"; \
echo ""; \
echo ""; \
echo ""; \
else \
    echo "No trilinos library directory specified"; \

    TRILINOS_DEFAULT_DIR="`pwd`/external_distributions/trilinos/trilinos_default_installation"
    TRILINOS_DEFAULT_LIB="`pwd`/external_distributions/trilinos/trilinos_default_installation/lib"
    TRILINOS_DEFAULT_INCLUDE="`pwd`/external_distributions/trilinos/trilinos_default_installation/include"

    # Does the default library folder exist?
    #---------------------------------------
    if test -d $TRILINOS_DEFAULT_LIB; then \
       have_trilinos=true; \
       trilinos_has_already_been_built=true; \
       echo "Found the trilinos libraries in default location: "; \
       echo " "; \
       echo "   "  $TRILINOS_DEFAULT_LIB ; \
       echo " "; \
       TRILINOS_DIR=$TRILINOS_DEFAULT_DIR; \
       TRILINOS_LIB=$TRILINOS_DEFAULT_LIB; \
       TRILINOS_INCLUDE=$TRILINOS_DEFAULT_INCLUDE; \

       if test -d $TRILINOS_INCLUDE; then \
          have_trilinos=true; \
	  trilinos_has_already_been_built=true; \
 	  echo "Found the trilinos include files in"; \
 	  echo " "; \
 	  echo "   "  $TRILINOS_INCLUDE ; \
 	  echo " "; \
       else \
          have_trilinos=false; \
	  trilinos_has_already_been_built=true; \
          echo "Did NOT find the trilinos include files in"; \
          echo " "; \
          echo "   "  $TRILINOS_INCLUDE ; \
          echo " "; \
       fi; \

    else \

       # Will trilinos library be built?
       #--------------------------------
       will_have_trilinos=false;
       AC_CHECK_FILES([external_distributions/trilinos/trilinos-9.0.2.tar.gz  external_distributions/trilinos/trilinos-10.8.5-Source.tar.gz external_distributions/trilinos/trilinos-11.4.3-Source.tar.gz external_distributions/trilinos/trilinos-11.8.1-Source.tar.gz], will_have_trilinos=true, )

       if test x$will_have_trilinos = xtrue; then \
           echo "Trilinos distribution .tar.gz file found ... will build trilinos"; \
           have_trilinos=true; \
           TRILINOS_DIR="`pwd`/external_distributions/trilinos/trilinos_default_installation"
           TRILINOS_LIB="`pwd`/external_distributions/trilinos/trilinos_default_installation/lib"
           TRILINOS_INCLUDE="`pwd`/external_distributions/trilinos/trilinos_default_installation/include"

           tar_file_needs_single_teuchos_library=false;
           AC_CHECK_FILES([external_distributions/trilinos/trilinos-9.0.2.tar.gz  external_distributions/trilinos/trilinos-10.8.5-Source.tar.gz], tar_file_needs_single_teuchos_library=true, )

           if test x$tar_file_needs_single_teuchos_library = xtrue; then \
                echo "We believe that this version of the trilinos library needs single teuchos library"; \
           else \
                echo "We believe that this version of the trilinos library needs multiple teuchos libraries"; \
           fi; \
       else \
           echo "Trilinos distribution .tar.gz file not found ... can't build trilinos"; \
       fi; \
    fi; \
fi


# Using multiple techos libraries?
#---------------------------------
use_single_teuchos_library=true;

# PM: Do we have a single Teuchos library? Placed here because if we
# find that the library has *already been installed* in the default
# location then we need to check what version it is.
version_file=$TRILINOS_INCLUDE/Trilinos_version.h;
if test x$have_trilinos = xtrue; then \
  if test -e $TRILINOS_INCLUDE; then \
       trilinos_version_number=`grep TRILINOS_MAJOR_MINOR_VERSION $version_file | awk -F' ' '{print $3}'`; \
    else \
       # If we're here then we don't have Trilinos unpacked yet however I can sneak into the tar
       # file without unpacking it (using zgrep) to get the version number...
       echo "I don't know what your Trilinos tar file is explicitly but I'm going to assume you only have one of them!"; \
       trilinos_directory_check="`pwd`/external_distributions/trilinos/"; \
       trilinos_version_number=`zgrep -a Trilinos_MAJOR_MINOR_VERSION "$trilinos_directory_check"trilinos*.tar.gz | sed 's/)//g' | cut -d' ' -f 2`; \
       echo "I managed to search through your tarred file and I found the version number."; \
    fi; \
    echo "Library version number is: "$trilinos_version_number; \
    if test "$trilinos_version_number" -gt 100805; then \
        multiple_teuchos_libs_requested=true; \
	echo "The version number is greater than 10.8.5. Therefore requesting multiple teuchos libraries!"; \
    fi; \
fi;

# If the user asked to make sure we assume multiple teuches libraries
if test x$multiple_teuchos_libs_requested = xtrue; then \
    echo "Requested multiple teuchos libraries via config flag"; \
    use_single_teuchos_library=false; \
    if test x$will_have_trilinos = xtrue; then \
        if test x$tar_file_needs_single_teuchos_library = xtrue; then \
           echo "Note that this changes the (sensible!) default inferred from version of tar file!"; \
        fi; \
    fi; \
else \
    echo "Not requested multiple teuchos library via config flag"; \
    if test x$will_have_trilinos = xtrue; then \
      if test x$tar_file_needs_single_teuchos_library = xtrue; then \
         echo "Building with single teuchos library"; \
         use_single_teuchos_library=true; \
      else
         echo "Building with multi teuchos libraries"; \
         use_single_teuchos_library=false; \
      fi; \
    fi; \
fi;


# Do we have the trilinos libraries? If so add them to the
#---------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_trilinos = xfalse; then
   TRILINOS_LIB='' ; \
   TRILINOS_INCLUDE='' ;\
else \
   if test "$TRILINOS_DIR" == "auto"; then \
     TRILINOS_LIB='' ; \
     TRILINOS_INCLUDE='' ;\
     external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_TRILINOS"` ; \
   else \
     if test x$use_single_teuchos_library = xtrue; then \
        echo "Actually building with single teuchos library"; \
        external_dist_libraries=`echo $external_dist_libraries " -lml -lifpack -lamesos -lanasazi -laztecoo -lepetraext -ltriutils -lepetra -lteuchos"` ; \
     else \
        echo "Actually building with multiple teuchos libraries"; \
        external_dist_libraries=`echo $external_dist_libraries " -lml -lifpack -lamesos -lanasazi -laztecoo -lepetraext -ltriutils -lepetra -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore "` ; \
     fi
     external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$TRILINOS_LIB` ; \
     if test x$dynamic_linking_for_trilinos_requested = xtrue; then \
        external_dist_libraries_lib=`echo $external_dist_libraries_lib -Wl,-rpath,$TRILINOS_LIB` ; \
     fi; \
     external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$TRILINOS_INCLUDE -DOOMPH_HAS_TRILINOS` ; \
  fi
fi

# Indicate that Trilinos has already been built (so we don't bother trying# to install 
# it again if there's a tar file available!)
AM_CONDITIONAL(ALREADY_HAVE_BUILT_TRILINOS, test x$trilinos_has_already_been_built = xtrue)

# Indicate that this build of oomph-lib will have Trilinos installed (which means it
# already has access to a built version or it will get built through the autogen process)
AM_CONDITIONAL(OOMPH_HAS_TRILINOS, test x$have_trilinos = xtrue)

# Make these variables available to all Makefile.am-s
AC_SUBST(TRILINOS_LIB)
AC_SUBST(TRILINOS_INCLUDE)


#................................................................



#---------------------------------------------------
# Potential location of the HLIB libraries.
#
# If the HLIB is available on your machine, specify the library directory
# i.e. the one above the lib and include directories with --with-hlib
#----------------------------------------------------
HLIB_LIB=
HLIB_INCLUDE=
check_for_hlib=false

AC_ARG_WITH(hlib,
[AS_HELP_STRING([--with-hlib], [Use hlib library. Assumed to be installed in the default install
               directorys: headers in /usr/local/include/HLib and library in
               /usr/local/lib.])]
      ,[
    check_for_hlib=true
    HLIB_LIB="/usr/local/lib"
    HLIB_INCLUDE="/usr/local/include/HLib"
])


# Do they exist?
have_hlib=false
if test x$check_for_hlib = xtrue; then \
if test -d $HLIB_LIB; then \
    have_hlib=true; \
    echo "Found the hlib libraries in"; \
    echo " "; \
    echo "   "  $HLIB_LIB ; \
    echo " "; \
else \
    have_hlib=false; \
    echo "Did NOT find the hlib libraries in"; \
    echo " "; \
    echo "   "  $HLIB_LIB ; \
    echo " "; \
fi; \
if test -d $HLIB_INCLUDE; then \
    have_hlib=true; \
          echo "Found the hlib include files in"; \
    echo " "; \
    echo "   "  $HLIB_INCLUDE ; \
    echo " "; \
else \
     have_hlib=false; \
     echo "Did NOT find the hlib include files in"; \
     echo " "; \
     echo "   "  $HLIB_INCLUDE ; \
     echo " "; \
fi; \
else \
    echo "No hlib library directory specified"; \

dnl dnl Not building hlib for people at the moment, probably not worth it
dnl     # Will hlib library be built?
dnl     #--------------------------------
dnl will_have_hlib=false;
dnl     AC_CHECK_FILES([external_distributions/hlib/HLib-1.3p20.tar.gz], will_have_hlib=true, )

dnl    if test x$will_have_hlib = xtrue; then \
dnl     echo "Hlib distribution .tar.gz file found ... will build hlib"; \
dnl     have_hlib=true; \
dnl     HLIB_DIR="`pwd`/external_distributions/hlib/hlib_default_installation"
dnl     HLIB_LIB="`pwd`/external_distributions/hlib/hlib_default_installation/lib"
dnl     HLIB_INCLUDE="`pwd`/external_distributions/hlib/hlib_default_installation/include"
dnl    else
dnl     echo "Hlib distribution .tar.gz file not found ... can't build hlib";
dnl    fi;

fi



# Do we have the hlib libraries? If so add them to the external
# distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_hlib = xfalse; then
  HLIB_LIB='' ; \
  HLIB_INCLUDE='' ;\
else \
  external_dist_libraries=`echo $external_dist_libraries " -lhmatrix"` ; \
  external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$HLIB_LIB` ; \
  external_dist_libraries_include=`echo $external_dist_libraries_include "-I"$HLIB_INCLUDE -DOOMPH_HAS_HLIB` ; \
fi

#Tell everybody about it...
AM_CONDITIONAL(OOMPH_HAS_HLIB, test x$have_hlib = xtrue)

# Make these variables available to all Makefile.am-s
AC_SUBST(HLIB_LIB)
AC_SUBST(HLIB_INCLUDE)


#................................................................
#................................................................

# Will mumps and scalapack be built?
#-----------------------------------

will_build_mumps_and_scalapack=false;

# Initialise test if all the specified libraries are available;
# gets switched to true when the first one (for scalapack) passes;
# then reset to false if/when any of the subsequent ones fail
all_specifiable_mumps_etc_libraries_available=false;


mpi_include_dir=

# Specify the mpi include directory that contains mpi.h
AC_ARG_WITH(mpi-include-directory,
[  --with-mpi-include-directory[=MPI_INCLUDE_DIRECTORY] specify mpi include
                          directory that contains mpi.h. Needed for scalapack
                          build.],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    mpi_include_dir=$withval
    have_mpi_include_directory=true
    MPI_INCLUDE_DIRECTORY=$withval
  fi
])
AC_SUBST(MPI_INCLUDE_DIRECTORY)


if test x$have_mpi_include_directory = xfalse; then \
        echo "Nobody has specified mpi include directory (the directory"; \
        echo "that contains mpi.h) with --with-mpi-include-directory"; \
        echo "therefore I cannot do mumps/scalapack build from source"; \
else \
        echo "mpi include directory (the directory that contains mpi.h) "; \
        echo "has been specified with --with-mpi-include-directory"; \
        echo "therefore I can do mumps/scalapack build from source"; \
fi;


AC_CHECK_FILES($mpi_include_dir/mpi.h,
                echo "Found mpi.h for scalapack build";
                found_mpi_h_file=true,
                echo "Not found mpi.h for scalapack build";
                found_mpi_h_file=false)


if test x$found_mpi_h_file = xfalse; then \
        echo "Did not find mpi.h in $mpi_include_dir for scalapack build"; \
fi;


AC_CHECK_FILES(external_distributions/mumps_and_scalapack/MUMPS_4.10.0.tar.gz,
                echo "Found MUMPS_4.10.0.tar.gz";
                found_mumps_tar_file=true,
                echo "Not found MUMPS_4.10.0.tar.gz";
                found_mumps_tar_file=false)

AC_CHECK_FILES(external_distributions/mumps_and_scalapack/scalapack_installer.tgz,
                echo "Found scalapack_installer.tgz";
                found_scalapack_tar_file=true,
                echo "Not found scalapack_installer.tgz";
                found_scalapack_tar_file=false)


if test x$found_mumps_tar_file = xtrue; then \
if test x$found_scalapack_tar_file = xfalse; then \
        echo "Need external_distributions/mumps_and_scalapack/MUMPS_4.10.0.tar.gz *and* "; \
        echo "external_distributions/mumps_and_scalapack/scalapack_installer.tgz "; \
        echo "to build mumps"; \
fi \
fi;

if test x$found_mumps_tar_file = xfalse; then \
if test x$found_scalapack_tar_file = xtrue; then \
        echo "Need external_distributions/mumps_and_scalapack/MUMPS_4.10.0.tar.gz *and* "; \
        echo "external_distributions/mumps_and_scalapack/scalapack_installer.tgz "; \
        echo "to build mumps"; \
fi \
fi;

if test x$found_mpi_h_file = xfalse; then \
        echo "Not doing mumps/scalapack build because location of mpi.h file not specified"; \
        will_build_mumps_and_scalapack=false; \
else \
        if test x$found_mumps_tar_file = xtrue; then \
        if test x$found_scalapack_tar_file = xtrue; then \
        if test x$want_mpi = xfalse; then \
             echo "Currently I can only build mumps/scalapack with mpi support"; \
             will_build_mumps_and_scalapack=false; \
        else \
             echo "Will build mumps/scalapack with mpi support"; \
             will_build_mumps_and_scalapack=true; \
             MUMPS_AND_SCALAPACK_DIR="`pwd`/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation"
             MUMPS_AND_SCALAPACK_LIB="`pwd`/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/lib"
             MUMPS_AND_SCALAPACK_INCLUDE="`pwd`/external_distributions/mumps_and_scalapack/mumps_and_scalapack_default_installation/include"
        fi \
        fi \
        fi \
fi;

# Make thes variables available to all Makefile.am-s
AC_SUBST(MUMPS_AND_SCALAPACK_LIB)
AC_SUBST(MUMPS_AND_SCALAPACK_INCLUDE)

#................................................................
#................................................................



#---------------------------------------------------
# Potential location for the scalapack libraries.
#
# If the complete scalapack library is available on
# your machine, specify the absolute path to the
# library with --with-scalapack
#----------------------------------------------------
SCALAPACK_LIB=
check_for_scalapack=false

AC_ARG_WITH(scalapack,
[  --with-scalapack[=ABSOLUTE_PATH] use scalapack library specified via ABSOLUTE PATH
                               e.g. /usr/local/lib/scalapack.a],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_scalapack=true
    SCALAPACK_LIB=$withval
  fi
])



# Do they exist?
have_scalapack=false
if test x$check_for_scalapack = xtrue; then \
if test -e $SCALAPACK_LIB; then \
    have_scalapack=true; \
    echo "Found the scalapack library "; \
    echo " "; \
    echo "   "  $SCALAPACK_LIB ; \
    echo " "; \
    # Make the variable available to all Makefile.am-s
    AC_SUBST(SCALAPACK_LIB)
    # This switches it to true because this test goes first!
    all_specifiable_mumps_etc_libraries_available=true;
else \
    have_scalapack=false; \
    echo "Did NOT find the scalapack library"; \
    echo " "; \
    echo "   "  $SCALAPACK_LIB ; \
    echo " "; \
    all_specifiable_mumps_etc_libraries_available=false; \
fi; \
else \
    echo "No scalapack library specified"; \
    if test x$will_build_mumps_and_scalapack = xtrue; then \
        echo "...but building it from source"; \
        SCALAPACK_LIB=$MUMPS_AND_SCALAPACK_LIB/libscalapack.a ; \
        echo "after which it will be: " $SCALAPACK_LIB ; \
        AC_SUBST(SCALAPACK_LIB) \
        have_scalapack=true; \
    fi; \
    all_specifiable_mumps_etc_libraries_available=false; \     
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_SCALAPACK, test x$have_scalapack = xtrue)





#................................................................
#................................................................



#---------------------------------------------------
# Potential location for the blacs libraries
#
# If the complete blacs library is available on
# your machine, specify the absolute path to the
# library directory (!) with --with-blacs
#----------------------------------------------------
BLACS_LIB=
check_for_blacs=false

AC_ARG_WITH(blacs,
[  --with-blacs[=ABSOLUTE_PATH] ABSOLUTE PATH to blacs libraries
                               e.g. /usr/local/lib/ [with trailing forward slash!]
                               This specified directory must contain
                               the files blacs.a blacsF77.a and blacs_copy.a
                               the latter being a copy of blacs.a.],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_blacs=true
    BLACS_LIB=$withval
  fi
])


# Do they exist?
have_blacs=false
if test x$check_for_blacs = xtrue; then \
if test -e $BLACS_LIB; then \
    have_blacs=true; \
    echo "Found the blacs library "; \
    echo " "; \
    echo "   "  $BLACS_LIB ; \
    echo " "; \
    # Make the variable available to all Makefile.am-s
    AC_SUBST(BLACS_LIB)
else \
    have_blacs=false; \
    echo "Did NOT find the blacs library"; \
    echo " "; \
    echo "   "  $BLACS_LIB ; \
    echo " "; \
    all_specifiable_mumps_etc_libraries_available=false; \
fi; \
else \
    echo "No blacs library specified"; \
    if test x$will_build_mumps_and_scalapack = xtrue; then \
        echo "...but building it from source"; \
        BLACS_LIB=$MUMPS_AND_SCALAPACK_LIB ; \
        echo "after which it will be: " $BLACS_LIB ; \
        AC_SUBST(BLACS_LIB) \
        have_blacs=true; \
    fi; \
    all_specifiable_mumps_etc_libraries_available=false; \      
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_BLACS, test x$have_blacs = xtrue)



#................................................................
#................................................................



#---------------------------------------------------
# Potential location for the pord library
#
# If the complete pord library is available on
# your machine, specify the absolute path to the
# library with --with-pord
#----------------------------------------------------
PORD_LIB=
check_for_pord=false

AC_ARG_WITH(pord,
[  --with-pord[=ABSOLUTE_PATH] ABSOLUTE PATH to pord libraries
                               e.g. /usr/local/lib/libpord.a ],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_pord=true
    PORD_LIB=$withval
  fi
])


# Do they exist?
have_pord=false
if test x$check_for_pord = xtrue; then \
if test -e $PORD_LIB; then \
    have_pord=true; \
    echo "Found the pord library "; \
    echo " "; \
    echo "   "  $PORD_LIB ; \
    echo " "; \
    # Make the variable available to all Makefile.am-s
    AC_SUBST(PORD_LIB)
else \
    have_pord=false; \
    echo "Did NOT find the pord library"; \
    echo " "; \
    echo "   "  $PORD_LIB ; \
    echo " "; \
    all_specifiable_mumps_etc_libraries_available=false; \
fi; \
else \
    echo "No pord library specified"; \
    if test x$will_build_mumps_and_scalapack = xtrue; then \
        echo "...but building it from source"; \
        PORD_LIB=$MUMPS_AND_SCALAPACK_LIB/libpord.a ; \
        echo "after which it will be: " $PORD_LIB ; \
        AC_SUBST(PORD_LIB) \
        have_pord=true; \
    fi; \
    all_specifiable_mumps_etc_libraries_available=false; \     
fi


#Tell everybody about it...
AM_CONDITIONAL(ALREADY_HAVE_PORD, test x$have_pord = xtrue)


#................................................................
#................................................................



#---------------------------------------------------
# Potential location of the MUMPS libraries.
#
# If the MUMPS is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-mumps
#----------------------------------------------------
MUMPS_DIR=
MUMPS_LIB=
MUMPS_INCLUDE=
check_for_mumps=false

AC_ARG_WITH(mumps,
[  --with-mumps[=DIR] use mumps library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_mumps=true
    MUMPS_DIR=$withval
    MUMPS_LIB="$withval/lib"
    MUMPS_INCLUDE="$withval/include"
  fi
])


# Do they exist?
have_mumps=false
if test x$check_for_mumps = xtrue; then \
if test -d $MUMPS_LIB; then \
    have_mumps=true; \
    echo "Found the mumps libraries in"; \
    echo " "; \
    echo "   "  $MUMPS_LIB ; \
    echo " "; \
else \
   if test "$MUMPS_DIR" == "auto"; then \
       have_mumps=true; \
       echo "Will assume that mumps library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
      have_mumps=false; \
      echo "Did NOT find the mumps libraries in"; \
      echo " "; \
      echo "   "  $MUMPS_LIB ; \
      echo " "; \
      all_specifiable_mumps_etc_libraries_available=false; \
   fi; \
fi; \
if test -d $MUMPS_INCLUDE; then \
    have_mumps=true; \
    echo "Found the mumps include files in"; \
    echo " "; \
    echo "   "  $MUMPS_INCLUDE ; \
    echo " "; \
else \
   if test "$MUMPS_DIR" == "auto"; then \
       have_mumps=true; \
       echo "Will assume that mumps library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
       have_mumps=false; \
       echo "Did NOT find the mumps include files in"; \
       echo " "; \
       echo "   "  $MUMPS_INCLUDE ; \
       echo " "; \
      all_specifiable_mumps_etc_libraries_available=false; \
    fi; \
fi; \
else \
    echo "No mumps library directory specified"; \
    if test x$will_build_mumps_and_scalapack = xtrue; then \
        echo "...but building it from source"; \
        MUMPS_DIR=$MUMPS_AND_SCALAPACK_DIR ; \
        MUMPS_LIB=$MUMPS_AND_SCALAPACK_LIB ; \
        MUMPS_INCLUDE=$MUMPS_AND_SCALAPACK_INCLUDE ; \
        echo "after which it will be: " $MUMPS_DIR ; \
        have_mumps=true; \
    fi; \
    all_specifiable_mumps_etc_libraries_available=false; \  
fi




if test x$all_specifiable_mumps_etc_libraries_available = xfalse; then \
   if test x$will_build_mumps_and_scalapack = xfalse; then \
      echo " "; \
      echo "====================================================================";\
      echo "Trouble: I don't have the tar files needed to build mumps";\
      echo "         and the associated libraries, but you haven't (correctly? ";\
      echo "         specified all the required library locations with the";\
      echo " " ;\
      echo "            --with-scalapack";\
      echo "            --with-blacs";\
      echo "            --with-pord";\
      echo "            --with-mumps";\
      echo " ";\
      echo "         flags in the configure options file.";\
      echo "         To stop stuff from breaking, I'll tell oomph-lib that";\
      echo "         we don't have mumps at all (though the existing";\
      echo "         will still be listed (but this should be harmless).";\
      echo "=====================================================================";\
      have_mumps=false;
   fi; \ 
fi


# Do we have the mumps libraries? If so add them to the
#---------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
# Note: for some reason we can't use -isystem (rather than -I) to include
# MUMPS libraries without showing warnings.
if test x$have_mumps = xfalse; then
  MUMPS_LIB='' ; \
  MUMPS_INCLUDE='' ;\
else \
  if test "$MUMPS_DIR" == "auto"; then \
    MUMPS_LIB='' ; \
    MUMPS_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_MUMPS"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -ldmumps -lmumps_common "` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$MUMPS_LIB` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-I"$MUMPS_INCLUDE -DOOMPH_HAS_MUMPS` ; \
  fi
fi


#Tell everybody about it...
AM_CONDITIONAL(OOMPH_HAS_MUMPS, test x$have_mumps = xtrue)

# Make these variables available to all Makefile.am-s
AC_SUBST(MUMPS_LIB)
AC_SUBST(MUMPS_INCLUDE)



# If we have mumps and associated libraries available (pre-installed; don't rebuild even
# if we have the tar files!
echo "all_specifiable_mumps_etc_libraries_available " $all_specifiable_mumps_etc_libraries_available
if test x$all_specifiable_mumps_etc_libraries_available = xfalse; then \
   echo "Possibly rebuilding mumps etc. libraries (if there's a tar file)"; \
else \   
   echo "Already have mumps etc libraries built elsewhere so not doing it again"; \
   will_build_mumps_and_scalapack=false; \
fi

AM_CONDITIONAL(WILL_BUILD_MUMPS_AND_SCALAPACK, test x$will_build_mumps_and_scalapack = xtrue)

# Make variable available to all Makefile.am-s
AC_SUBST(WILL_BUILD_MUMPS_AND_SCALAPACK)


#................................................................


#................................................................
# Start build cgal and associated third party distributions
#................................................................



#---------------------------------------------------
# Potential location for the boost libraries.
#
# If the BOOST is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-boost
#----------------------------------------------------
BOOST_DIR=
BOOST_LIB=
BOOST_LIB=
check_for_boost=false
dont_build_boost=false
AC_ARG_WITH(boost,
[  --with-boost[=DIR]  use boost library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_boost=true
    dont_build_boost=true
    BOOST_DIR=$withval
    BOOST_LIB="$withval/lib"
    BOOST_INCLUDE="$withval/include"
  fi
])


# Do they exist?
have_boost=false
if test x$check_for_boost = xtrue; then \
if test -d $BOOST_LIB; then \
    have_boost=true; \
    echo "Found the boost libraries in"; \
    echo " "; \
    echo "   "  $BOOST_LIB ; \
    echo " "; \
else \
   if test "$BOOST_DIR" == "auto"; then \
       have_boost=true; \
       echo "Will assume that boost library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_boost=false; \
       echo "Did NOT find the boost libraries in"; \
       echo " "; \
       echo "   "  $BOOST_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $BOOST_INCLUDE; then \
    have_boost=true; \
    echo "Found the boost include files in"; \
    echo " "; \
    echo "   "  $BOOST_INCLUDE ; \
    echo " "; \
else \
   if test "$BOOST_DIR" == "auto"; then \
       have_boost=true; \
       echo "Will assume that boost library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
      have_boost=false; \
      echo "Did NOT find the boost include files in"; \
      echo " "; \
      echo "   "  $BOOST_INCLUDE ; \
      echo " "; \
    fi; \
fi; \
else \
   echo "No boost library directory specified"; \

   # Will boost library be built?
   #-----------------------------
   if test x$want_cgal_build = xtrue; then  \
      echo "Building boost within oomph-lib as part of cgal build process"; \
      will_have_boost=true; \
   else \
      echo "Not building boost within oomph-lib as part of cgal build process"; \
      will_have_boost=false; \
   fi; \
   if test x$will_have_boost = xtrue; then \
    have_boost=true; \
    if test x$install_cgal_in_permanent_dir = xtrue; then  \
      BOOST_DIR="$CGAL_PERMANENT_INSTALLATION_DIR/boost_default_installation";
      BOOST_LIB="$CGAL_PERMANENT_INSTALLATION_DIR/boost_default_installation/lib";
      BOOST_INCLUDE="$CGAL_PERMANENT_INSTALLATION_DIR/boost_default_installation/include";
    else \

      BOOST_DIR="`pwd`/external_distributions/boost/boost_default_installation";
      BOOST_LIB="`pwd`/external_distributions/boost/boost_default_installation/lib";
      BOOST_INCLUDE="`pwd`/external_distributions/boost/boost_default_installation/include";
    fi;
    echo "BOOST_DIR    : " $BOOST_DIR
    echo "BOOST_LIB    : " $BOOST_LIB
    echo "BOOST_INCLUDE: " $BOOST_INCLUDE
   fi;
fi



# Do we have the boost libraries? If so add them to the
#------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_boost = xfalse; then
  BOOST_LIB='' ; \
  BOOST_INCLUDE='' ;\
else \
  if test "$BOOST_DIR" == "auto"; then \
    BOOST_LIB='' ; \
    BOOST_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_BOOST"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -lpthread -lboost_thread -lboost_system "` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$BOOST_LIB -Wl,-rpath,$BOOST_LIB ` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$BOOST_INCLUDE -DOOMPH_HAS_BOOST` ; \
  fi
fi



#Tell everybody about it...
#AM_CONDITIONAL(OOMPH_HAS_BOOST, test x$have_boost = xtrue)

#Tell everybody about it...
AM_CONDITIONAL(DONT_BUILD_BOOST, test x$dont_build_boost = xtrue)

# Make thes variables available to all Makefile.am-s
AC_SUBST(BOOST_LIB)
AC_SUBST(BOOST_INCLUDE)
AC_SUBST(BOOST_DIR)

#................................................................

#................................................................


#---------------------------------------------------
# Potential location for the gmp libraries.
#
# If the GMP is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-gmp
#----------------------------------------------------
GMP_DIR=
GMP_LIB=
GMP_LIB=
check_for_gmp=false
dont_build_gmp=false
AC_ARG_WITH(gmp,
[  --with-gmp[=DIR]  use gmp library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_gmp=true
    dont_build_gmp=true
    GMP_DIR=$withval
    GMP_LIB="$withval/lib"
    GMP_INCLUDE="$withval/include"
  fi
])


# Do they exist?
have_gmp=false
if test x$check_for_gmp = xtrue; then \
if test -d $GMP_LIB; then \
    have_gmp=true; \
    echo "Found the gmp libraries in"; \
    echo " "; \
    echo "   "  $GMP_LIB ; \
    echo " "; \
else \
   if test "$GMP_DIR" == "auto"; then \
       have_gmp=true; \
       echo "Will assume that gmp library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_gmp=false; \
       echo "Did NOT find the gmp libraries in"; \
       echo " "; \
       echo "   "  $GMP_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $GMP_INCLUDE; then \
    have_gmp=true; \
    echo "Found the gmp include files in"; \
    echo " "; \
    echo "   "  $GMP_INCLUDE ; \
    echo " "; \
else \
   if test "$GMP_DIR" == "auto"; then \
       have_gmp=true; \
       echo "Will assume that gmp library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
      have_gmp=false; \
      echo "Did NOT find the gmp include files in"; \
      echo " "; \
      echo "   "  $GMP_INCLUDE ; \
      echo " "; \
    fi; \
fi; \
else \
   echo "No gmp library directory specified"; \

   # Will gmp library be built?
   #---------------------------
   if test x$want_cgal_build = xtrue; then  \
      echo "Building gmp within oomph-lib as part of cgal build process"; \
      will_have_gmp=true; \
   else \
      echo "Not building gmp within oomph-lib as part of cgal build process"; \
      will_have_gmp=false; \
   fi; \
   if test x$will_have_gmp = xtrue; then \
    have_gmp=true; \
    if test x$install_cgal_in_permanent_dir = xtrue; then  \
      GMP_DIR="$CGAL_PERMANENT_INSTALLATION_DIR/gmp_default_installation";
      GMP_LIB="$CGAL_PERMANENT_INSTALLATION_DIR/gmp_default_installation/lib";
      GMP_INCLUDE="$CGAL_PERMANENT_INSTALLATION_DIR/gmp_default_installation/include";
    else \
      GMP_DIR="`pwd`/external_distributions/gmp/gmp_default_installation"
      GMP_LIB="`pwd`/external_distributions/gmp/gmp_default_installation/lib"
      GMP_INCLUDE="`pwd`/external_distributions/gmp/gmp_default_installation/include"
    fi;
    echo "GMP_DIR    : " $GMP_DIR
    echo "GMP_LIB    : " $GMP_LIB
    echo "GMP_INCLUDE: " $GMP_INCLUDE
   fi;
fi



# Do we have the gmp libraries? If so add them to the
#------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_gmp = xfalse; then
  GMP_LIB='' ; \
  GMP_INCLUDE='' ;\
else \
  if test "$GMP_DIR" == "auto"; then \
    GMP_LIB='' ; \
    GMP_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_GMP"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -lgmp"` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$GMP_LIB` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$GMP_INCLUDE -DOOMPH_HAS_GMP` ; \
  fi
fi

#Tell everybody about it...
#AM_CONDITIONAL(OOMPH_HAS_GMP, test x$have_gmp = xtrue)

#Tell everybody about it...
AM_CONDITIONAL(DONT_BUILD_GMP, test x$dont_build_gmp = xtrue)

# Make thes variables available to all Makefile.am-s
AC_SUBST(GMP_LIB)
AC_SUBST(GMP_INCLUDE)
AC_SUBST(GMP_DIR)

#................................................................


#---------------------------------------------------
# Potential location for the mpfr libraries.
#
# If the MPFR is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-mpfr
#----------------------------------------------------
MPFR_DIR=
MPFR_LIB=
MPFR_LIB=
check_for_mpfr=false
dont_build_mpfr=false
AC_ARG_WITH(mpfr,
[  --with-mpfr[=DIR]  use mpfr library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_mpfr=true
    dont_build_mpfr=true
    MPFR_DIR=$withval
    MPFR_LIB="$withval/lib"
    MPFR_INCLUDE="$withval/include"
  fi
])


# Do they exist?
have_mpfr=false
if test x$check_for_mpfr = xtrue; then \
if test -d $MPFR_LIB; then \
    have_mpfr=true; \
    echo "Found the mpfr libraries in"; \
    echo " "; \
    echo "   "  $MPFR_LIB ; \
    echo " "; \
else \
   if test "$MPFR_DIR" == "auto"; then \
       have_mpfr=true; \
       echo "Will assume that mpfr library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_mpfr=false; \
       echo "Did NOT find the mpfr libraries in"; \
       echo " "; \
       echo "   "  $MPFR_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $MPFR_INCLUDE; then \
    have_mpfr=true; \
    echo "Found the mpfr include files in"; \
    echo " "; \
    echo "   "  $MPFR_INCLUDE ; \
    echo " "; \
else \
   if test "$MPFR_DIR" == "auto"; then \
       have_mpfr=true; \
       echo "Will assume that mpfr library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
      have_mpfr=false; \
      echo "Did NOT find the mpfr include files in"; \
      echo " "; \
      echo "   "  $MPFR_INCLUDE ; \
      echo " "; \
    fi; \
fi; \
else \
   echo "No mpfr library directory specified"; \

   # Will mpfr library be built?
   #----------------------------
   if test x$want_cgal_build = xtrue; then  \
      echo "Building mpfr within oomph-lib as part of cgal build process"; \
      will_have_mpfr=true; \
   else \
      echo "Not building mpfr within oomph-lib as part of cgal build process"; \
      will_have_mpfr=false; \
   fi; \
   if test x$will_have_mpfr = xtrue; then \
    have_mpfr=true; \
    if test x$install_cgal_in_permanent_dir = xtrue; then  \
      MPFR_DIR="$CGAL_PERMANENT_INSTALLATION_DIR/mpfr_default_installation";
      MPFR_LIB="$CGAL_PERMANENT_INSTALLATION_DIR/mpfr_default_installation/lib";
      MPFR_INCLUDE="$CGAL_PERMANENT_INSTALLATION_DIR/mpfr_default_installation/include";
    else \
      MPFR_DIR="`pwd`/external_distributions/mpfr/mpfr_default_installation"
      MPFR_LIB="`pwd`/external_distributions/mpfr/mpfr_default_installation/lib"
      MPFR_INCLUDE="`pwd`/external_distributions/mpfr/mpfr_default_installation/include"
    fi;
    echo "MPFR_DIR    : " $MPFR_DIR
    echo "MPFR_LIB    : " $MPFR_LIB
    echo "MPFR_INCLUDE: " $MPFR_INCLUDE
   fi;
fi



# Do we have the mpfr libraries? If so add them to the
#------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_mpfr = xfalse; then
  MPFR_LIB='' ; \
  MPFR_INCLUDE='' ;\
else \
  if test "$MPFR_DIR" == "auto"; then \
    MPFR_LIB='' ; \
    MPFR_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_MPFR"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -lmpfr"` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$MPFR_LIB` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$MPFR_INCLUDE -DOOMPH_HAS_MPFR` ; \
  fi
fi



#Tell everybody about it...
#AM_CONDITIONAL(OOMPH_HAS_MPFR, test x$have_mpfr = xtrue)

#Tell everybody about it...
AM_CONDITIONAL(DONT_BUILD_MPFR, test x$dont_build_mpfr = xtrue)

# Make thes variables available to all Makefile.am-s
AC_SUBST(MPFR_LIB)
AC_SUBST(MPFR_INCLUDE)
AC_SUBST(MPFR_DIR)

#................................................................


#---------------------------------------------------
# Potential location for the cgal libraries.
#
# If the CGAL is available on
# your machine, specify the library directory
# i.e. the one above the lib and include directories
# with --with-cgal
#----------------------------------------------------
CGAL_DIR=
CGAL_LIB=
CGAL_LIB=
check_for_cgal=false
dont_build_cgal=false
AC_ARG_WITH(cgal,
[  --with-cgal[=DIR]  use cgal library in directory DIR
                          i.e. library is in DIR/lib
                          and include files are in DIR/include
                          or specify "auto" (without quotation marks!)
                          if the library is specified automatically
                          by your compiler/linker],[
  if test "$withval" != "no" -a "$withval" != "yes"; then
    check_for_cgal=true
    dont_build_cgal=true
    CGAL_DIR=$withval
    CGAL_LIB="$withval/lib"
    CGAL_INCLUDE="$withval/include"
  fi
])


# Do they exist?
have_cgal=false
if test x$check_for_cgal = xtrue; then \
if test -d $CGAL_LIB; then \
    have_cgal=true; \
    echo "Found the cgal libraries in"; \
    echo " "; \
    echo "   "  $CGAL_LIB ; \
    echo " "; \
else \
   if test "$CGAL_DIR" == "auto"; then \
       have_cgal=true; \
       echo "Will assume that cgal library is specified automatically by"; \
       echo "compiler/linker."; \
   else \
       have_cgal=false; \
       echo "Did NOT find the cgal libraries in"; \
       echo " "; \
       echo "   "  $CGAL_LIB ; \
       echo " "; \
   fi;\
fi; \
if test -d $CGAL_INCLUDE; then \
    have_cgal=true; \
    echo "Found the cgal include files in"; \
    echo " "; \
    echo "   "  $CGAL_INCLUDE ; \
    echo " "; \
else \
   if test "$CGAL_DIR" == "auto"; then \
       have_cgal=true; \
       echo "Will assume that cgal library is specified automatically by"; \
       echo "compiler/linker."; \
    else \
      have_cgal=false; \
      echo "Did NOT find the cgal include files in"; \
      echo " "; \
      echo "   "  $CGAL_INCLUDE ; \
      echo " "; \
    fi; \
fi; \
else \
   echo "No cgal library directory specified"; \

   # Will cgal library be built?
   #----------------------------
   if test x$want_cgal_build = xtrue; then  \
      echo "Building cgal within oomph-lib as part of cgal build process"; \
      will_have_cgal=true; \
   else \
      echo "Not building cgal within oomph-lib as part of cgal build process"; \
      will_have_cgal=false; \
   fi; \
   if test x$will_have_cgal = xtrue; then \
    have_cgal=true; \
    if test x$install_cgal_in_permanent_dir = xtrue; then  \
      CGAL_DIR="$CGAL_PERMANENT_INSTALLATION_DIR/cgal_default_installation";
      CGAL_LIB="$CGAL_PERMANENT_INSTALLATION_DIR/cgal_default_installation/lib";
      CGAL_INCLUDE="$CGAL_PERMANENT_INSTALLATION_DIR/cgal_default_installation/include";
    else \
      CGAL_DIR="`pwd`/external_distributions/cgal/cgal_default_installation"
      CGAL_LIB="`pwd`/external_distributions/cgal/cgal_default_installation/lib"
      CGAL_INCLUDE="`pwd`/external_distributions/cgal/cgal_default_installation/include"
    fi;
    echo "CGAL_DIR    : " $CGAL_DIR
    echo "CGAL_LIB    : " $CGAL_LIB
    echo "CGAL_INCLUDE: " $CGAL_INCLUDE
   fi;
fi



# Do we have the cgal libraries? If so add them to the
#------------------------------------------------------
# external distributions; if not wipe the variables
#--------------------------------------------------
if test x$have_cgal = xfalse; then
  CGAL_LIB='' ; \
  CGAL_INCLUDE='' ;\
else \
  if test "$CGAL_DIR" == "auto"; then \
    CGAL_LIB='' ; \
    CGAL_INCLUDE='' ;\
    external_dist_libraries_include=`echo $external_dist_libraries_include " -DOOMPH_HAS_CGAL"` ; \
  else \
    external_dist_libraries=`echo $external_dist_libraries " -lCGAL_Core -lCGAL "` ; \
    # external_dist_libraries=`echo $external_dist_libraries " -lCGAL_Core -lCGAL_ImageIO -lCGAL "` ; \
    external_dist_libraries_lib=`echo $external_dist_libraries_lib "-L"$CGAL_LIB -Wl,-rpath,$CGAL_LIB` ; \
    external_dist_libraries_include=`echo $external_dist_libraries_include "-isystem "$CGAL_INCLUDE -DOOMPH_HAS_CGAL` ; \
  fi
fi



#Tell everybody about it...
AM_CONDITIONAL(OOMPH_HAS_CGAL, test x$have_cgal = xtrue)

#Tell everybody about it...
AM_CONDITIONAL(DONT_BUILD_CGAL, test x$dont_build_cgal = xtrue)

# Make thes variables available to all Makefile.am-s
AC_SUBST(CGAL_LIB)
AC_SUBST(CGAL_INCLUDE)
AC_SUBST(CGAL_DIR)


#... doc ...

# Do we want to build cgal and associated third party distributions ?
if test x$want_cgal_build = xtrue; then  \

    if test x$check_for_cgal = xfalse; then \
      echo "We're building cgal and associated third party distributions" ; \
    fi

else \

  echo "Not building cgal and associated third party distributions "; \
  echo " "; \

  # Install cgal and associated third party distributions in permament
  # directory?
  if test x$install_cgal_in_permanent_dir = xtrue; then  \
     echo " "; \
     echo "================================================================ "; \
     echo "WARNING: Makes no sense to specify cgal install directory"; \
     echo "         with --with-cgal-permanent-installation-dir when"; \
     echo "         suppressing cgal build with --enable-suppress-cgal-build"; \
     echo "         We won't build it (and therefore not install it either)."; \
     echo "================================================================ "; \
     echo " "; \
  fi; \
  echo " "; \

fi;

# Install cgal and associated third party distributions in permament
# directory?
if test x$want_cgal_build = xtrue; then  \
if test x$install_cgal_in_permanent_dir = xtrue; then  \

  echo "Installing cgal and associated third party distributions in:" $CGAL_PERMANENT_INSTALLATION_DIR; \

else \

  echo "Installing cgal and associated third party distributions within oomph-lib's external_distributions directory. "; \

fi;
fi;

#................................................................
# end cgal build
#................................................................



# Make external libraries with their own distributions available to
#------------------------------------------------------------------
# all Makfile.am-s. Note that by adding to AM_LDFLAGS and AM_CPPFLAGS
#--------------------------------------------------------------------
# everything automatically links against these libraries.
#--------------------------------------------------------
# Also add any other accumulated flags for these variables.
#----------------------------------------------------------
AM_LDFLAGS=`echo $external_dist_libraries_lib`
AC_SUBST(AM_LDFLAGS)

AM_CPPFLAGS=`echo $external_dist_libraries_include  $accumulated_cpp_flags`
AC_SUBST(AM_CPPFLAGS)


# Define "external" libraries
#----------------------------
# These are the third-party libraries distributed with oomph-lib.
# They are only used if the user hasn't specified his/her own
# local version of these.
external_libs='-loomph_hsl -loomph_crbond_bessel'

# Triangle?
if test x$suppress_triangle = xfalse; then \
  external_libs=`echo $external_libs " -loomph_triangle "`
else
  echo "Omitting specification of triangle library as we haven't built"
  echo "it locally. Could/should provide option to link against external"
  echo "version..."
fi

# Tetgen?
if test x$suppress_tetgen = xfalse; then \
  external_libs=`echo $external_libs " -loomph_tetgen "`
else
  echo "Omitting specification of tetgen library as we haven't built"
  echo "it locally. Could/should provide option to link against external"
  echo "version..."
fi


# SuperLU?
if test x$have_superlu = xfalse; then
  external_libs=`echo $external_libs " -loomph_superlu_4.3"`
else
  echo "Omitting specification of superlu library -- assumed to be"
  echo "provided automatically by compiler/linker..."
fi



# Include superlu_dist and parallel metis if compiled in parallel
if test x$want_mpi = xtrue; then

  # parmetis?
  if test x$have_parmetis = xfalse; then
    external_libs=`echo $external_libs " -loomph_parmetis_3.1.1"`
  else
    echo "Omitting specification of parmetis library -- assumed to be"
    echo "provided automatically by compiler/linker..."
  fi

  # SuperLU_dist?
  if test x$have_superlu_dist = xfalse; then
    external_libs=`echo $external_libs " -loomph_superlu_dist_3.0"`
  else
    echo "Omitting specification of superlu_dist library -- assumed to be"
    echo "provided automatically by compiler/linker..."
  fi
fi


#Note: serial metis has to go last
if test x$have_metis = xfalse; then
  external_libs=`echo $external_libs " -loomph_metis_from_parmetis_3.1.1"`
else
  echo "Omitting specification of metis library -- assumed to be"
  echo "provided automatically by compiler/linker..."
fi


# Pord?
if test x$have_pord = xtrue; then
  external_libs=`echo " $PORD_LIB " $external_libs `
fi



# Blacs?
if test x$have_blacs = xtrue; then
  external_libs=`echo  " $BLACS_LIB/blacs.a $BLACS_LIB/blacsF77.a  $BLACS_LIB/blacs_copy.a" $external_libs`
fi


# Scalapack?
if test x$have_scalapack = xtrue; then
  external_libs=`echo " $SCALAPACK_LIB " $external_libs `
fi


# Lapack?
if test x$have_lapack = xfalse; then
  external_libs=`echo $external_libs " -loomph_lapack -loomph_flapack"`
else
  if test "$LAPACK_LIB" != "auto"; then
      external_libs=`echo $external_libs " $LAPACK_LIB"`
  else
      echo "Omitting specification of lapack library -- assumed to be"
      echo "provided automatically by compiler/linker..."
  fi
fi


# Blas?
if test x$have_blas = xfalse; then
  external_libs=`echo $external_libs " -loomph_blas"`
else
  if test "$BLAS_LIB" != "auto"; then
      external_libs=`echo $external_libs " $BLAS_LIB"`
  else
      echo "Omitting specification of blas library -- assumed to be"
      echo "provided automatically by compiler/linker..."
  fi
fi



echo " "
echo "/////////////////////////////////////////////////////////////////////////////////"
echo " "
echo "external_dist_libraries        : $external_dist_libraries"
echo "external_dist_libraries_lib    : $external_dist_libraries_lib"
echo "external_dist_libraries_include: $external_dist_libraries_include"
echo " "
echo "/////////////////////////////////////////////////////////////////////////////////"
echo " "

# Now add the third-party external libraries that are shared by all
#------------------------------------------------------------------
# Makefile.am-s
#--------------
external_libs=`echo $external_dist_libraries " "  $external_libs`
EXTERNAL_LIBS=`echo $external_libs`
AC_SUBST(EXTERNAL_LIBS)

# Include the file containing the constructed AC_CONFIG_FILES([....]) command
m4_include([config/configure.ac_scripts/makefile_list])


# Create the output
AC_OUTPUT()



#-----------------------------------------------------------------
# This macro adjusts config.h to oomph-lib-config.h
# in which all macros are package specific to avoid clashes
# with those declared by other third-party libraries.
#
# Downloaded and adjusted from
#
#  http://autoconf-archive.cryp.to/ac_create_prefix_config_h.html
#
# Changes:
#  (1) ac_prefix_conf_INP wasn't set correctly. Commented the
#      the relevant bit of code out and hard-coded the
#      name as config.h
#  (2) In the AC_MSG_ERROR macro rather than referring to the variable
#      $ac_prefix_conf_INP the non-existent (and therefore empty)
#      variable $ac_prefix_conf_IN was referred to. Therefore no
#      meaningful error message was produced.
#-----------------------------------------------------------------
AC_DEFUN([AC_CREATE_PREFIX_CONFIG_H],
[changequote({, })dnl
ac_prefix_conf_OUT=`echo ifelse($1, , $PACKAGE-config.h, $1)`
ac_prefix_conf_DEF=`echo _$ac_prefix_conf_OUT | sed -e 'y:abcdefghijklmnopqrstuvwxyz./,-:ABCDEFGHIJKLMNOPQRSTUVWXYZ____:'`
ac_prefix_conf_PKG=`echo ifelse($2, , $PACKAGE, $2)`
ac_prefix_conf_LOW=`echo _$ac_prefix_conf_PKG | sed -e 'y:ABCDEFGHIJKLMNOPQRSTUVWXYZ-:abcdefghijklmnopqrstuvwxyz_:'`
ac_prefix_conf_UPP=`echo $ac_prefix_conf_PKG | sed -e 'y:abcdefghijklmnopqrstuvwxyz-:ABCDEFGHIJKLMNOPQRSTUVWXYZ_:'  -e '/^[0-9]/s/^/_/'`
#ac_prefix_conf_INP=`echo ifelse($3, , _, $3)`
#if test "$ac_prefix_conf_INP" = "_"; then
#   case $ac_prefix_conf_OUT in
#      */*) ac_prefix_conf_INP=`basename $ac_prefix_conf_OUT`
#      ;;
#      *-*) ac_prefix_conf_INP=`echo $ac_prefix_conf_OUT | sed -e 's/[a-zA-Z0-9_]*-//'`
#      ;;
#      *) ac_prefix_conf_INP=config.h
#      ;;
#   esac
#fi
ac_prefix_conf_INP=config.h
changequote([, ])dnl
if test -z "$ac_prefix_conf_PKG" ; then
   AC_MSG_ERROR([no prefix for _PREFIX_PKG_CONFIG_H])
else
  AC_MSG_RESULT(creating $ac_prefix_conf_OUT - prefix $ac_prefix_conf_UPP for $ac_prefix_conf_INP defines)
  if test -f $ac_prefix_conf_INP ; then
#    AC_ECHO_MKFILE([/* automatically generated */], $ac_prefix_conf_OUT)
changequote({, })dnl
    echo '#ifndef '$ac_prefix_conf_DEF >>$ac_prefix_conf_OUT
    echo '#define '$ac_prefix_conf_DEF' 1' >>$ac_prefix_conf_OUT
    echo ' ' >>$ac_prefix_conf_OUT
    echo /'*' $ac_prefix_conf_OUT. Generated automatically at end of configure. '*'/ >>$ac_prefix_conf_OUT

    echo 's/#undef  *\([A-Z_]\)/#undef '$ac_prefix_conf_UPP'_\1/' >conftest.sed
    echo 's/#undef  *\([a-z]\)/#undef '$ac_prefix_conf_LOW'_\1/' >>conftest.sed
    echo 's/#define  *\([A-Z_][A-Za-z0-9_]*\)\(.*\)/#ifndef '$ac_prefix_conf_UPP"_\\1 \\" >>conftest.sed
    echo '#define '$ac_prefix_conf_UPP"_\\1 \\2 \\" >>conftest.sed
    echo '#endif/' >>conftest.sed
    echo 's/#define  *\([a-z][A-Za-z0-9_]*\)\(.*\)/#ifndef '$ac_prefix_conf_LOW"_\\1 \\" >>conftest.sed
    echo '#define '$ac_prefix_conf_LOW"_\\1 \\2 \\" >>conftest.sed
    echo '#endif/' >>conftest.sed
    sed -f conftest.sed $ac_prefix_conf_INP >>$ac_prefix_conf_OUT
    echo ' ' >>$ac_prefix_conf_OUT
    echo '/*' $ac_prefix_conf_DEF '*/' >>$ac_prefix_conf_OUT
    echo '#endif' >>$ac_prefix_conf_OUT
    touch -r $ac_prefix_conf_INP $ac_prefix_conf_OUT
changequote([, ])dnl
  else
    AC_MSG_ERROR([input file $ac_prefix_conf_INP does not exist, dnl
    skip generating $ac_prefix_conf_OUT])
  fi
  rm -f conftest.*
fi])



AC_CREATE_PREFIX_CONFIG_H

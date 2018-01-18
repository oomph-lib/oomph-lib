SUBDIRS = hypre trilinos 

if WANT_CGAL_BUILD

if ! DONT_BUILD_GMP
SUBDIRS += gmp
endif

# Order matters! mpfr needs gmp
if ! DONT_BUILD_MPFR
SUBDIRS += mpfr
endif

# doesn't actually seem to be needed...
#if ! DONT_BUILD_ZLIB
#SUBDIRS += zlib
#endif

if ! DONT_BUILD_BOOST
SUBDIRS += boost
endif

if ! DONT_BUILD_CGAL
SUBDIRS += cgal
endif

all:
	rm -f cgal_configure_flags.txt
if CGAL_PERMANENT_INSTALLATION_DIR_SPECIFIED
	echo "--with-boost="${BOOST_DIR} >  cgal_configure_flags.txt
	echo "--with-gmp="${GMP_DIR}     >> cgal_configure_flags.txt
	echo "--with-mpfr="${MPFR_DIR}   >> cgal_configure_flags.txt
	echo "--with-cgal="${CGAL_DIR}  >> cgal_configure_flags.txt
endif

endif

if WANT_MPI
if WILL_BUILD_MUMPS_AND_SCALAPACK
SUBDIRS += mumps_and_scalapack
endif
endif

clean-local:
	rm -f cgal_configure_flags.txt

distclean-local:
	rm -f cgal_configure_flags.txt


EXTRA_DIST=hlib

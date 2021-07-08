#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Name of executables: 
check_PROGRAMS= two_d_poisson  two_d_poisson_without_singularity_subtracted_off

#---------------------------------------------------------------------------

# Local sources that the Poisson code depends on:
two_d_poisson_without_singularity_subtracted_off_SOURCES = two_d_poisson.cc poisson_elements_with_singularity.h

# Required libraries: Only the "generic" and "poisson" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. The generic sources also require the "external" libraries
# that are shipped with oomph-lib. The Fortran libraries, $(FLIBS), get 
# included just in case we decide to use a solver that involves Fortran 
# sources. 
# NOTE: The order in which libraries are specified is important! 
#       The most specific ones (that might involve references to other
#       libraries) need to be listed first, the more basic ones 
#       later. In this example we have (from right to left, i.e. from
#       general to specific):
#       -- The fortran libraries: They  are compiler specific and
#          obviously can't depend on any code that we (or others)
#          have written. $(FLIBS) is a variable that automake will translate
#          to the actual fortran libraries.
#       -- The external (third party) libraries: They cannot depend on 
#          any of our code. The variable $(EXTERNAL_LIBS) is defined
#          in the machine-generated file configure.ac in oomph-lib's 
#          home directory. 
#       -- Oomph-lib's generic library contains oomph-lib's fundamental 
#          objects which do not depend on any specific system of PDEs 
#          or element types.
#       -- Finally, oomph-lib's poisson library contains oomph-lib's 
#          Poisson elements which refer to objects from the generic
#          library. 
two_d_poisson_without_singularity_subtracted_off_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

two_d_poisson_without_singularity_subtracted_off_CXXFLAGS = -DDONT_USE_SINGULARITY

#---------------------------------------------------------------------------

# Local sources that the Poisson code depends on:
two_d_poisson_SOURCES = two_d_poisson.cc poisson_elements_with_singularity.h

# Required libraries: Only the "generic" and "poisson" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. The generic sources also require the "external" libraries
# that are shipped with oomph-lib. The Fortran libraries, $(FLIBS), get 
# included just in case we decide to use a solver that involves Fortran 
# sources. 
# NOTE: The order in which libraries are specified is important! 
#       The most specific ones (that might involve references to other
#       libraries) need to be listed first, the more basic ones 
#       later. In this example we have (from right to left, i.e. from
#       general to specific):
#       -- The fortran libraries: They  are compiler specific and
#          obviously can't depend on any code that we (or others)
#          have written. $(FLIBS) is a variable that automake will translate
#          to the actual fortran libraries.
#       -- The external (third party) libraries: They cannot depend on 
#          any of our code. The variable $(EXTERNAL_LIBS) is defined
#          in the machine-generated file configure.ac in oomph-lib's 
#          home directory. 
#       -- Oomph-lib's generic library contains oomph-lib's fundamental 
#          objects which do not depend on any specific system of PDEs 
#          or element types.
#       -- Finally, oomph-lib's poisson library contains oomph-lib's 
#          Poisson elements which refer to objects from the generic
#          library. 
two_d_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


EXTRA_DIST+=run_convergence_test.bash conv.lay 

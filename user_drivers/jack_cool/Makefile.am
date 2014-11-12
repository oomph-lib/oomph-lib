# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= jacks_poisson_code jacks_own_code

#---------------------------------------------------------------------------

# Local sources that Jack's Poisson code depends on:
jacks_poisson_code_SOURCES = jacks_poisson_code.cc

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
#          any of our code.
#       -- Oomph-lib's generic library contains oomph-lib's fundamental 
#          objects which do not depend on any specific system of PDEs 
#          or element types.
#       -- Finally, oomph-lib's poisson library contains oomph-lib's 
#          Poisson elements which refer to objects from the generic
#          library. 
jacks_poisson_code_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  

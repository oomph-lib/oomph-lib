
# A few file definitions
#-----------------------

# Define the sources
sources = segregated_fsi_solver.cc pseudo_elastic_preconditioner.cc \
pseudo_elastic_fsi_preconditioner.cc

# Include files which shouldn't be compiled
incl_cc_files =

# Define the headers
headers =  \
fsi_preconditioners.h segregated_fsi_solver.h pseudo_elastic_preconditioner.h \
pseudo_elastic_fsi_preconditioner.h multi_domain_boussinesq_elements.h \
boussinesq_elements.h helmholtz_time_harmonic_linear_elasticity_interaction.h \
fourier_decomposed_helmholtz_time_harmonic_linear_elasticity_interaction.h \
pml_helmholtz_time_harmonic_linear_elasticity_interaction.h
# old_multi_domain_boussinesq_elements.h

# Template only files. These should be included in include directory
# but should not be listed in combined header file
templates =

# Define name of library
libname = multi_physics

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libmulti_physics.la

# Sources that the library depends on:
#-------------------------------------
libmulti_physics_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/multi_physics

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files) $(templates)


# Required libraries
#-------------------
# Include path for library headers -- need to refer to include files
# in their respective source directoriesas they will not have been
# installed yet!
AM_CPPFLAGS += -I$(top_srcdir)/src/generic -I$(top_srcdir)/src/navier_stokes


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers) $(incl_cc_files) $(templates) )
	(echo "$(headers) $(incl_cc_files) $(templates)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list )
else
install-data-hook:
endif


# Cleanup: Kill automatically generated combined header file
clean-local:
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~

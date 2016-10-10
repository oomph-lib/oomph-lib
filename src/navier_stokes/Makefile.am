
# A few file definitions
#-----------------------

# Define the sources
sources =  \
navier_stokes_elements.cc \
refineable_navier_stokes_elements.cc \
Tnavier_stokes_elements.cc \
navier_stokes_preconditioners.cc \
lagrange_enforced_flow_preconditioner.h

# fluid_traction_elements.cc 

# Include files which shouldn't be compiled
incl_cc_files =

# Define the headers
headers =  \
navier_stokes_elements.h  \
fluid_traction_elements.h \
refineable_navier_stokes_elements.h  \
Tnavier_stokes_elements.h \
navier_stokes_preconditioners.h \
navier_stokes_surface_power_elements.h \
navier_stokes_surface_drag_torque_elements.h \
navier_stokes_flux_control_elements.h \
impose_parallel_outflow_element.h \
impose_impenetrability_element.h \
lagrange_enforced_flow_preconditioner.h

# Template only files. These should be included in include directory
# but should not be listed in combined header file
templates = 

#fluid_traction_elements.template.cc

# Define name of library
libname = navier_stokes

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libnavier_stokes.la

# Sources that the library depends on: 
#-------------------------------------
libnavier_stokes_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h 

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/navier_stokes

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files) $(templates)


# Required libraries 
#-------------------
# Include path for library headers -- need to refer to include files
# in src as they will not have been installed yet! 
AM_CPPFLAGS += -I$(top_srcdir)/src/generic 

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
	rm -f $(libname).h  *~
distclean-local:
	rm -f $(libname).h  *~



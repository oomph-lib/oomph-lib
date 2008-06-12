#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=3d_static_cap	 

# Sources for executable
3d_static_cap_SOURCES = 3d_static_cap.cc fix_vol_int_elements.h \
                        fix_vol_int_elastic_elements.h

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
3d_static_cap_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
                           -lsolid -lconstitutive -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

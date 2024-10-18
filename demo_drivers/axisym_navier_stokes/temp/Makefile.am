#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=axi_static_cap

# Sources for second executable
axi_static_cap_SOURCES = axi_static_cap.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
axi_static_cap_LDADD = -L@libdir@ -laxisym_navier_stokes -lfluid_interface \
                            -lsolid -lconstitutive \
                            -lgeneric \
                            $(EXTERNAL_LIBS) $(FLIBS)

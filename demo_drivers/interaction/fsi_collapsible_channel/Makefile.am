#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
fsi_collapsible_channel_macro_CR \
fsi_collapsible_channel_macro_adapt_CR \
fsi_collapsible_channel_algebraic_CR \
fsi_collapsible_channel_algebraic_adapt_CR \
fsi_collapsible_channel_macro_TH \
fsi_collapsible_channel_macro_adapt_TH \
fsi_collapsible_channel_algebraic_TH \
fsi_collapsible_channel_algebraic_adapt_TH  \
fsi_pseudo_solid_collapsible_channel_adapt_TH \
fsi_pseudo_solid_collapsible_channel_adapt_CR \
fsi_pseudo_solid_collapsible_channel_TH \
fsi_pseudo_solid_collapsible_channel_CR 

#---------------------------------------------------------------------

# Sources for first executable
fsi_pseudo_solid_collapsible_channel_adapt_TH_SOURCES = fsi_pseudo_solid_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_pseudo_solid_collapsible_channel_adapt_TH_LDADD = \
 -L@libdir@ -lbeam -lconstitutive -lsolid -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)
fsi_pseudo_solid_collapsible_channel_adapt_TH_CXXFLAGS = -DTAYLOR_HOOD

#---------------------------------------------------------------------

# Sources for first executable
fsi_pseudo_solid_collapsible_channel_adapt_CR_SOURCES = fsi_pseudo_solid_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_pseudo_solid_collapsible_channel_adapt_CR_LDADD = \
 -L@libdir@ -lbeam -lconstitutive -lsolid -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for first executable
fsi_pseudo_solid_collapsible_channel_TH_SOURCES = fsi_pseudo_solid_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_pseudo_solid_collapsible_channel_TH_LDADD = \
 -L@libdir@ -lbeam -lconstitutive -lsolid -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)
fsi_pseudo_solid_collapsible_channel_TH_CXXFLAGS = -DTAYLOR_HOOD

#---------------------------------------------------------------------

# Sources for first executable
fsi_pseudo_solid_collapsible_channel_CR_SOURCES = fsi_pseudo_solid_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_pseudo_solid_collapsible_channel_CR_LDADD = \
 -L@libdir@ -lbeam -lconstitutive -lsolid -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_macro_CR_SOURCES = fsi_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_CR_LDADD = \
 -L@libdir@ -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_macro_CR_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE

#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_macro_adapt_CR_SOURCES = fsi_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_adapt_CR_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_macro_adapt_CR_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE

#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_algebraic_CR_SOURCES = fsi_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_algebraic_CR_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_algebraic_adapt_CR_SOURCES = \
fsi_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_algebraic_adapt_CR_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------


#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_macro_TH_SOURCES = fsi_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_TH_LDADD = \
 -L@libdir@ -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_macro_TH_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE \
                                            -DTAYLOR_HOOD

#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_macro_adapt_TH_SOURCES = fsi_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_adapt_TH_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_macro_adapt_TH_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE \
                                                  -DTAYLOR_HOOD

#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_algebraic_TH_SOURCES = fsi_collapsible_channel.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_algebraic_TH_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_algebraic_TH_CXXFLAGS = -DTAYLOR_HOOD

#---------------------------------------------------------------------

# Sources for first executable
fsi_collapsible_channel_algebraic_adapt_TH_SOURCES = \
fsi_collapsible_channel_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_algebraic_adapt_TH_LDADD = -L@libdir@ -lbeam \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

fsi_collapsible_channel_algebraic_adapt_TH_CXXFLAGS = -DTAYLOR_HOOD

#---------------------------------------------------------------------

EXTRA_DIST+=anim.mcr 

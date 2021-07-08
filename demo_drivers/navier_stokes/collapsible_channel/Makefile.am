#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=collapsible_channel collapsible_channel_bl_squash \
collapsible_channel_algebraic collapsible_channel_algebraic_bl_squash \
collapsible_channel_adaptive_algebraic collapsible_channel_adaptive_algebraic_bl_squash

#---------------------------------------------------------------------------

# Sources for executable
collapsible_channel_adaptive_algebraic_bl_squash_SOURCES = \
collapsible_channel_algebraic.cc my_alg_channel_mesh.h

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_adaptive_algebraic_bl_squash_LDADD = \
                               -L@libdir@ -lnavier_stokes -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)


# Comile bl squash adaptive version with approriate flag
collapsible_channel_adaptive_algebraic_bl_squash_CXXFLAGS = \
-DUSE_BL_SQUASH_FCT \
-DADAPTIVE


#---------------------------------------------------------------------------



# Sources for executable
collapsible_channel_adaptive_algebraic_SOURCES = \
collapsible_channel_algebraic.cc my_alg_channel_mesh.h

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_adaptive_algebraic_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                                       $(EXTERNAL_LIBS) $(FLIBS)

# Comile  adaptive version with approriate flag
collapsible_channel_adaptive_algebraic_CXXFLAGS = -DADAPTIVE


#---------------------------------------------------------------------------

# Sources for executable
collapsible_channel_algebraic_bl_squash_SOURCES = \
collapsible_channel_algebraic.cc my_alg_channel_mesh.h

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_algebraic_bl_squash_LDADD = \
                               -L@libdir@ -lnavier_stokes -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)


# Comile bl squash version with approriate flag
collapsible_channel_algebraic_bl_squash_CXXFLAGS = -DUSE_BL_SQUASH_FCT


#---------------------------------------------------------------------------



# Sources for executable
collapsible_channel_algebraic_SOURCES = \
collapsible_channel_algebraic.cc my_alg_channel_mesh.h 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_algebraic_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                                       $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------



# Sources for executable
collapsible_channel_SOURCES = collapsible_channel.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                            $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Sources for executable
collapsible_channel_bl_squash_SOURCES = collapsible_channel.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
collapsible_channel_bl_squash_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                                      $(EXTERNAL_LIBS) $(FLIBS)

# Comile bl squash version with approriate flag
collapsible_channel_bl_squash_CXXFLAGS = -DUSE_BL_SQUASH_FCT


#---------------------------------------------------------------------------

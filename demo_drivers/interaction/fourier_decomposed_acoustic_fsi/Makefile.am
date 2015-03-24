#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
fourier_decomposed_acoustic_fsi \
unstructured_fourier_decomposed_acoustic_fsi \
adaptive_unstructured_fourier_decomposed_acoustic_fsi 


#---------------------------------------------------------------------

# Sources for first executable
fourier_decomposed_acoustic_fsi_SOURCES = fourier_decomposed_acoustic_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fourier_decomposed_acoustic_fsi_LDADD =  -L@libdir@ \
	-lfourier_decomposed_helmholtz \
        -ltime_harmonic_fourier_decomposed_linear_elasticity  \
        -lnavier_stokes \
        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for first executable
unstructured_fourier_decomposed_acoustic_fsi_SOURCES = unstructured_fourier_decomposed_acoustic_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_fourier_decomposed_acoustic_fsi_LDADD =  -L@libdir@ \
	-lfourier_decomposed_helmholtz \
        -ltime_harmonic_fourier_decomposed_linear_elasticity  \
        -lnavier_stokes \
        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#---------------------------------------------------------------------

# Sources for first executable
adaptive_unstructured_fourier_decomposed_acoustic_fsi_SOURCES = unstructured_fourier_decomposed_acoustic_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_fourier_decomposed_acoustic_fsi_LDADD =  -L@libdir@ \
	-lfourier_decomposed_helmholtz \
        -ltime_harmonic_fourier_decomposed_linear_elasticity  \
        -lnavier_stokes \
        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_unstructured_fourier_decomposed_acoustic_fsi_CXXFLAGS = -DADAPTIVE

#---------------------------------------------------------------------

EXTRA_DIST += exact.map animate_displacement_with_mesh_fd.mcr


SUBDIRS = \
driven_cavity \
adaptive_driven_cavity \
hp_adaptive_driven_cavity \
circular_driven_cavity \
osc_ring \
single_layer_free_surface \
two_layer_interface \
three_d_free_surface \
adaptive_interface \
inclined_plane \
bretherton \
three_d_bretherton \
static_cap \
three_d_entry_flow \
curved_pipe \
rayleigh_channel \
rayleigh_traction_channel \
osc_ellipse \
spine_channel \
collapsible_channel \
flow_past_cylinder \
turek_flag_non_fsi \
channel_with_leaflet \
free_surface_rotation \
three_d_static_cap \
falling_jet \
vmtk_fluid \
unstructured_fluid \
unstructured_three_d_fluid \
flux_control \
schur_complement_preconditioner \
lagrange_enforced_flow_preconditioner

if SUPPRESS_TRIANGLE_LIB
else
  SUBDIRS += unstructured_adaptive_ALE unstructured_adaptive_fs jeffery_orbit
endif

if SUPPRESS_TETGEN_LIB
else
  SUBDIRS += unstructured_adaptive_3d_ALE
endif

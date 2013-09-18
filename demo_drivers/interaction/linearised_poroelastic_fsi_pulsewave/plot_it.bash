
#oomph-convert -z -p2 poro_traction*.dat; makePvd poro_traction poro_traction.pvd

#oomph-convert -z soln-poro*.dat; 
makePvd soln-poro soln-poro.pvd

#oomph-convert -z -p2 fluid_traction*.dat; makePvd fluid_traction fluid_traction.pvd

#oomph-convert -z -p2 fluid_fpsi_traction*.dat; makePvd fluid_fpsi_traction fluid_fpsi_traction.pvd

#oomph-convert -z soln-fluid*.dat; 
makePvd soln-fluid soln-fluid.pvd

oomph-convert -z -p2 regular_fluid*.dat; makePvd regular_fluid regular_fluid.pvd
oomph-convert -z -p2 regular_poro*.dat; makePvd regular_poro regular_poro.pvd

#oomph-convert -z pulse_wave*.dat; makePvd pulse_wave pulse_wave.pvd

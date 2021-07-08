
oomph-convert -z -p2 solid_traction*.dat; makePvd solid_traction solid_traction.pvd

oomph-convert -z soln-solid*.dat; makePvd soln-solid soln-solid.pvd

oomph-convert -z -p2 fluid_traction*.dat; makePvd fluid_traction fluid_traction.pvd

oomph-convert -z -p2 fluid_fsi_traction*.dat; makePvd fluid_fsi_traction fluid_fsi_traction.pvd

oomph-convert -z soln-fluid*.dat; makePvd soln-fluid soln-fluid.pvd

oomph-convert -z -p2 regular_fluid*.dat; makePvd regular_fluid regular_fluid.pvd

 
oomph-convert -z pulse_wave*.dat; makePvd pulse_wave pulse_wave.pvd

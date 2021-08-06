# Contributors

The `oomph-lib` "architects" (in alphabetical order) are
[**Andrew Hazel**](https://github.com/alhazel) (**co-founder**) &
[**Matthias Heil**](https://github.com/MatthiasHeilManchester) (**co-founder**) assisted by
former/current project/MSc/PhD students and collaborators who made (or are still
making) significant contributions to the development of the library (listed in
reverse chronological order):


Person                              | Position        | Contribution
------------------------------------|-----------------|--------------------------------------------------------------------------------------------
[**Puneet Matharu**](https://github.com/PuneetMatharu) | **Maintainer**  | Worked on the implementation of a space-time finite-element modelling of the flow past an oscillating cylinder, geometric multigrid solvers, particularly for the Helmholtz equations, and various other small additions.
[**Jonathan Deakin**](https://github.com/jondea) | **Maintainer**  | Worked on the numerical solution of acoustic fluid-structure interaction problems as PhD student, overhauled the `oomph-lib` website, and worked on a glaciology-related melt problem.
[**Christian Vaquero-Stainer**](https://github.com/orgs/oomph-lib/people/wobblygomboc) | **Contributor** | _**PM to MH: Fill this in.**_
**Thierry Gonon**                   | **Contributor** | Implemented the methodology to subtract singular (or non-singular) functions off solutions to the Poisson and Navier-Stokes equations.
**Chris Johnson**                   | **Contributor** | Has provided many bug fixes.
**Chihebeddine Hammami**            | **Contributor** | Worked on the implementation of Yulii Shikhmurzaev's interface formation theory.
**Narjes Akriche**                  | **Contributor** | Worked on pseudo-resonances in acoustic fluid-structure interaction problems.
**Aman Rajvardhan**                 | **Contributor** | Worked on implementing surfactant transport equations in two- and three-dimensional geometries.
**Jordan Rosso**                    | **Contributor** | Worked on topological fluid mechanics of the Karman vortex street.
**Florian Molinier**                | **Contributor** | Did some early work on the coupled solution of the axisymmetric free-surface Navier-Stokes equations and the axisymmetric Foeppl von Karman equations.
**Draga Pihler-Puzovic**            | **Contributor** | Worked on the the coupled solution of the Foeppl von Karman equations and the Reynolds lubrication equation to model wrinkling/fingering in elastic-walled Hele-Shaw cells.
**Joris Ferrand**                   | **Contributor** | Worked on the solution of the Foeppl von Karman equations.
**Harsh Ranjan**                    | **Contributor** | Worked on multiple solutions of Navier-Stokes flows in curved tubes.
**Anton Martinsson**                | **Contributor** | Implemented the machinery required to output `oomph-lib` data in Paraview format, bypassing the need for running the time-consuming Tecplot to Paraview conversion scripts. He also implemented the displacement-based axisymmetric Foeppl von Karman equations.
**Andr&eacute; Von Borries**        | **Contributor** | Is working on free-surface Navier-Stokes and lubrication theory problems.
**Matthew Walker**                  | **Contributor** | Implemented PML methods for the azimuthally Fourier-decomposed Helmholtz equations.
**Joris Ferrand**                   | **Contributor** | Implemented the axisymmetric Foeppl von Karman equations.
**Philippe Mesnard**                | **Contributor** | Worked on acoustic FSI problems and introduced many improvements to `oomph-lib`'s machinery for handling such problems.
**Florian Molinier**                | **Contributor** | Worked on the coupling of the free surface Navier-Stokes equations and the axisymmetric Foeppl von Karman equations (in the context of simulating flows in elastic-walled Hele-Shaw problems).
**David Nigro**                     | **Contributor** | Developed and implemented much of the machinery for acoustic fluid-structure interaction problems.
**Matthew Russell**                 | **Contributor** | Implemented the Foeppl-von-Karman equations; he now continues to work on poro-elastic FSI problems.
**Raphael Perillat**                | **Contributor** | Worked on the simulation of flows in elastic-walled Hele-Shaw cells.
**Robert Harter**                   | **Contributor** | Works on acoustic fluid-structure interaction problems.
**Radu Cimpeanu**                   | **Contributor** | Implemented the PML boundary conditions for the Helmholtz equations and the time-harmonic equations of linear elasticity.
**Julio Perez Sansalvador**         | **Contributor** | Works on parallel unstructured mesh adaptation.
**David Shepherd**                  | **Contributor** | Works on the numerical solution of micromagnetic problems.
**Ray White**                       | **Contributor** | Is working on block preconditioners.
**Nico Bergemann**                  | **Contributor** | Made (and continues to make) significant contributions to the adaptive unstructured mesh (re-)generation capabilities for free-surface problems.
**Ben Saxby**                       | **Contributor** | Works on hp adaptivity and XFEM.
**Michael Crabb**                   | **Contributor** | Worked on Discontinuous Galerkin (DG) methods.
**Peter Ashcroft**                  | **Contributor** | Worked on eigenvalue problems.
**Jeremy van Chu**                  | **Contributor** | Contributed to the completion the Tecplot to Paraview conversion scripts and significantly extended the [Paraview tutorial](doc/Paraview/html/index.html). He also developed the `LineVisualiser` machinery (which allows the extraction of computational data along lines in a higher-dimensional domain) and wrote the domain-based tube mesh.
**Guilherme Rocha**                 | **Contributor** | Developed elements to simulate Hele-Shaw problems (by solving the free-surface Reynolds lubrication equations).
**Ahmed Wassfi**                    | **Contributor** | Extended **Tarak Kharrat's** work on the Helmholtz equation and implemented the Fourier-decomposed version of this equation.
**Alexandre Raczynski**             | **Contributor** | Keeps providing bug fixes and contributed to the completion the Tecplot to Paraview conversion scripts discussed in the [Paraview tutorial](doc/Paraview/html/index.html).
**David Rutter**                    | **Contributor** | Wrote the [tutorial for the linear elasticity equations](doc/linear_elasticity/periodic_load/html/index.html).
**Tarak Kharrat**                   | **Contributor** | Implemented the Helmholtz elements and the methodology to apply the Sommerfeld radiation condition.
**Luigi Collucci**                  | **Contributor** | Continued Benjamin Metz's work and developed the interface from `oomph-lib` to `Triangle`.
**Francisco Jose Blanco Rodriguez** | **Contributor** | Worked on free-surface problems and wrote the driver code that simulates the Rayleigh instability of an axisymmetric jet.
**Wassamon Phusakulkajorn**         | **Contributor** | Worked on C1-continuous triangular finite elements for shell, beam and biharmonic problems.
**Benjamin Metz**                   | **Contributor** | Worked on adaptivity and solution transfer for unstructured meshes.
**Amine Massit**                    | **Contributor** | Worked on outflow boundary conditions for Navier-Stokes problems and physiological FSI problems based on meshes generated by [vmtk](http://www.vmtk.org/).
**Patrick Hurley**                  | **Contributor** | Works on free surface Navier-Stokes problems.
**Andy Gait**                       | **Contributor** | worked on parallelisation, in particular the problem distribution and the subsequent distributed mesh adaptation.
**Angelo Simone**                   | **Contributor** | Wrote Python scripts that convert `oomph-lib`'s output to the VTU format that can be read by [Paraview](http://www.Paraview.org); see the [Paraview tutorial](doc/Paraview/html/index.html) for details.
**Sophie Kershaw**                  | **Contributor** | Worked on the Navier-Stokes equations in spherical coordinates.
**Floraine Cordier**                | **Contributor** | Developed the driver codes and tutorials for the [flow past the elastic leaflet](doc/interaction/fsi_channel_with_leaflet/html/index.html) and [Turek & Hron's FSI benchmark](doc/interaction/turek_flag/html/index.html). In the process she significantly extended `oomph-lib`'s FSI capabilities.
**Stefan Kollmannsberger**          | **Contributor** | Stefan and his students Iason Papaioannou and Orkun Oezkan Doenmez developed early versions of the code for [Turek & Hron's FSI benchmark](doc/interaction/turek_flag/html/index.html) and its [non-FSI counterpart](doc/navier_stokes/turek_flag_non_fsi/html/index.html).
**Cedric Ody**                      | **Contributor** | Developed the `YoungLaplace` elements and their refineable counterparts to study capillary statics problems.
**Alice Gaertig**                   | **Contributor** | Developed interfaces to the third-party mesh generators [`Triangle`](http://www.cs.cmu.edu/~quake/triangle.html), [`TetGen`](http://wias-berlin.de/software/tetgen//), [`Geompack++`](http://members.shaw.ca/bjoe/) [`CQMesh`](http://www.dimap.ufrn.br/~mfsiqueira/Marcelo_Siqueiras_Web_Spot/cqmesh.html).
**Claire Blancon**                  | **Contributor** | Developed the demo drivers for the collapsible channel problem (with and without fluid-structure interaction).
**Nick Chapman**                    | **Contributor** | Worked on the implementation of triangular and tet-elements.
**Chris Gold**                      | **Contributor** | Implemented explicit timestepping schemes.
**Phil Haines**                     | **Contributor** | Worked on bifurcation detection and tracking for the Navier-Stokes equations and developed the formulation of the equations in plane polar coordinates.
**Richard Muddle**                  | **Contributor** | Worked on the block preconditioning techniques for the biharmonic (and many other) equations, and parallel solvers.
**Glyn Rees**                       | **Contributor** | Worked on iterative linear solvers and multigrid
**Alberto de Lozar**                | **Contributor** | Worked on 3D free-surface Navier-Stokes problems.
**Jonathan Boyle**                  | **Contributor** | Developed the initial interfaces to third-party iterative solvers and is now involved in the further parallelisation of the code and the implementation and application of block-preconditioning techniques for Navier-Stokes and fluid-structure interaction problems.
**Renaud Schleck**                  | **Contributor** | Completed the octree-based mesh refinement procedures and wrote the MPI-based parallel assembly routines and the interfaces to [`SuperLU_DIST`](https://portal.nersc.gov/project/sparse/superlu/#superlu_dist).
**Sharaf Al-Sharif**                | **Contributor** | Provided the initial implementation of nodal spectral elements.
**Daniel Meyer**                    | **Contributor** | Used oomph-lib to study a variety of axisymmetric Navier-Stokes problems, with and without free surfaces, and developed drafts for many of our tutorials.
**Alexandre Klimowicz**             | **Contributor** | Worked on block-preconditioning methods.
**Jean-Michel Lenoir**              | **Contributor** | Implemented the first part of the octree-based 3D mesh refinement procedures.
**Gemma Barson**                    | **Contributor** | Provided the initial implementation for the 2D Delaunay mesh generation procedures.

We're always looking for more help! Get in touch if you're interested
in joining the team.

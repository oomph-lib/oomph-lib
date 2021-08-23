/**

\mainpage People

The \c oomph-lib "architects" are (in no particular order)
 
\image html architects.gif " " 
\image latex architects.eps " " width=0.75\textwidth

...assisted by former/current project/MSc/PhD students and
collaborators who made 
(or are still making) significant contributions to the development of
the library (listed in reverse chronological order):
- \b Patrick \b Keuchel worked on forced oscillations and resonances
  in axisymmetric fluid-conveying tubes.
- \b Ben \b Gavan worked on mesh adaptation procedures using gmsh.
- \b Rupinder \b Matharu worked on the simulation of creep processes
  during the annealing of HDPE materials.
- \b Aidan \b Retallick works on the simulation of graphene-based
  pressure transducers.
- \b Christian \b Vaquero-Stainer did his Masters project on modelling
  fingering in Hele-Shaw cells with time-dependent gap-widths, and then
  moved on to this PhD on the subtraction of singularities in the flow
  past arbitrarily-shaped disk-like objects in Stokes flow.
- \b David \b Robinson developed and implemented C1 continuous
  triangular elements for the solution of fourth-order PDES. He used
  them for the implementation of the Foeppl-von-Karman and Koiter-Steigman
  plate theories.
- <b>Thomas Brion, Simon Finney and Hannah Chamberlain</b> all worked
  on Foeppl-von-Karman based models of graphene-based microphones.
- \b Louis \b Calot-Plaetevoet improved the methodology used to
  transfer solutions between different meshes.
- \b Thierry \b Gonon implemented the methodology to subtract
  singular (or non-singular) functions off solutions to the Poisson
  and Navier-Stokes equations. 
- \b Chris \b Johnson has provided many bug fixes. 
- \b Puneet \b Matharu works on the implementation of geometric
  multigrid solvers, particularly for Helmholtz equations. He then did
  his PhD on exotic wakes in the flow past circular cylinders and has
  become one of the maintainers of oomph-lib,
- \b Chihebeddine \b Hammami worked on the implementation of 
  Yulii Shikhmurzaev's interface formation theory. 
- \b Narjes \b Akriche worked on pseudo-resonances in acoustic 
  fluid-structure interaction problems. 
- \b Aman \b Rajvardhan worked on implementing surfactant transport
  equations in two- and three-dimensional geometries.
- \b Jordan \b Rosso worked on topological fluid mechanics of the
  the Karman vortex street.
- \b Florian \b Molinier did some early work on the coupled solution
  of the axisymmetric free-surface Navier-Stokes equations and 
  the axisymmetric Foeppl von Karman equations.
- \b Jonathan \b Deakin worked on a glaciology-related melt problem
  (and has since returned as PhD student to work on the numerical
  solution of acoustic fluid-structure interaction problems and
  optimal PML methods. He has since become one of the maintainers of oomph-lib.
- \b Draga \b Pihler-Puzovic worked on the the coupled solution of the
   Foeppl von Karman equations and the Reynolds lubrication equation
   to model wrinkling/fingering in elastic-walled Hele-Shaw cells.
- \b Joris \b Ferrand worked on the solution of the Foeppl von Karman
  equations.
- \b Harsh \b Ranjan worked on multiple solutions of Navier--Stokes flows in curved tubes.
- \b Anton \b Martinsson implemented the machinery required to 
  output \c oomph-lib data in paraview format, bypassing the need
  for running the time-consuming tecplot to paraview conversion scripts.
  He also implemented the displacement-based axisymmetric Foeppl von
  Karman equations. 
- \b Andr&eacute; \b Von \b Borries is working on free-surface Navier--Stokes and
  lubrication theory problems.
- \b Matthew \b Walker implemented PML methods for the azimuthally
  Fourier-decomposed Helmholtz equations.
- \b Joris \b Ferrand implemented the axisymmetric Foeppl von Karman
  equations.
- \b Philippe \b Mesnard worked on acoustic FSI problems and introduced many
  improvements to \c oomph-lib's machinery for handling such problems.
- \b Florian \b Molinier worked on the coupling of the free surface
  Navier-Stokes equations and the axisymmetric Foeppl von Karman 
  equations (in the context of simulating flows in elastic-walled
  Hele-Shaw problems).
- \b David \b Nigro developed and implemented much of the machinery
  for acoustic fluid-structure interaction problems.
- \b Matthew \b Russell implemented the Foeppl-von-Karman equations;
  he now continues to work on poro-elastic FSI problems.
- \b Raphael \b Perillat worked on the simulation of flows in
  elastic-walled Hele-Shaw cells.
- \b Robert \b Harter works on acoustic fluid-structure interaction
  problems.
- \b Radu \b Cimpeanu implemented the PML boundary conditions for the
  Helmholtz equations and the time-harmonic equations of linear
  elasticity.
- \b Julio \b Perez \b Sansalvador works on parallel unstructured
  mesh adaptation.
- \b David \b Shepherd works on the numerical solution of
  micromagnetic problems.
- \b Ray \b White is working on block preconditioners.
- \b Nico \b Bergemann made significant
  contributions to the adaptive unstructured mesh (re-)generation
  capabilities for free-surface problems. He then did his PhD on
  the simulation of viscous and visco-plastic free surface flows. 
- \b Ben \b Saxby works on hp adaptivity and XFEM.
- \b Michael \b Crabb worked on Discontinuous Galerkin (DG) methods.
- \b Peter \b Ashcroft worked on eigenvalue problems. 
- \b Jeremy \b van \b Chu contributed to the completion the tecplot to paraview
  conversion scripts and significantly extended the 
  <a href="../../paraview/html/index.html">the paraview tutorial.</a>
  He also developed the \c LineVisualiser machinery (which allows
  the extraction of computational data along lines in a
  higher-dimensional domain) and wrote the domain-based tube mesh.
- \b Guilherme \b Rocha developed elements to simulate Hele-Shaw 
  problems (by solving the free-surface Reynolds lubrication
  equations).
- \b Ahmed \b Wassfi extended \b Tarak \b Kharrat's work on the Helmholtz
  equation and implemented the Fourier-decomposed version of this
  equation.
- \b Alexandre \b Raczynski keeps providing bug fixes and contributed 
   to the completion the tecplot to paraview conversion scripts
   discussed in the <a href="../../paraview/html/index.html">the 
   paraview tutorial.</a>
- \b David \b Rutter wrote the 
   <a href="../../linear_elasticity/periodic_load/html/index.html">
   tutorial for the linear elasticity equations.</a>
- \b Tarak \b Kharrat implemented the Helmholtz elements and the
  methodology to apply the Sommerfeld radiation condition.
- \b Luigi \b Collucci continued Benjamin Metz's work and developed the
  interface from \c oomph-lib to \c Triangle. 
- \b Francisco \b Jose \b Blanco \b Rodriguez worked on free-surface
   problems and wrote the driver code that simulates the Rayleigh 
   instability of an axisymmetric jet. 
- \b Wassamon \b Phusakulkajorn worked on C1-continuous triangular
   finite elements for shell, beam and biharmonic problems.
- \b Benjamin \b Metz worked on adaptivity and solution transfer for 
  unstructured meshes. 
- \b Amine \b Massit worked on outflow boundary conditions for
  Navier-Stokes problems and physiological FSI problems based
  on meshes generated by <a href="http://www.vmtk.org/">vmtk</a>.
- \b Patrick \b Hurley works on free surface Navier-Stokes problems.
- <B>Andy Gait</B>  worked on parallelisation, in particular
  the problem distribution and the subsequent distributed mesh adaptation.
- <B>Angelo 
  Simone</B> wrote python scripts that convert \c oomph-lib's
  output to the vtu format that can be read by 
  <a href="http://www.paraview.org">paraview</a>;
  see <a href="../../paraview/html/index.html">the paraview tutorial</a>
  for details. 
- \b Sophie \b Kershaw worked on the Navier-Stokes equations
  in spherical coordinates.
- \b Floraine \b Cordier developed the driver codes and tutorials for the 
  <a href="../../interaction/fsi_channel_with_leaflet/html/index.html"> 
  flow past the elastic leaflet</a> and
  <a href="../../interaction/turek_flag/html/index.html"> 
  Turek & Hron's FSI benchmark</a>. In the process she significantly
  extended \c oomph-lib's FSI capabilities.
- <B> Stefan Kollmannsberger</B> and his students Iason Papaioannou and  
  Orkun Oezkan Doenmez developed early versions of the code
  for <a href="../../interaction/turek_flag/html/index.html"> 
  Turek & Hron's FSI benchmark</a> and its 
  <a href="../../navier_stokes/turek_flag_non_fsi/html/index.html">
   non-FSI counterpart.</a>
- <B>Cedric Ody</B> developed the \c YoungLaplace elements and their refineable
  counterparts to study capillary statics problems.
- \b Alice \b Gaertig developed interfaces to the third-party mesh
  generators <A HREF="http://www.cs.cmu.edu/~quake/triangle.html"> 
  \c Triangle, </A>  <A HREF="http://wias-berlin.de/software/tetgen//"> 
  \c TetGen, </A>, <A HREF="http://members.shaw.ca/bjoe/">\c
  Geompack++, </A> and <A HREF="http://www.dimap.ufrn.br/~mfsiqueira/Marcelo_Siqueiras_Web_Spot/cqmesh.html">
  \c CQMesh. </A>
- \b Claire \b Blancon developed the demo drivers for the collapsible   
  channel problem (with and without fluid-structure interaction). 
- \b Nick \b Chapman worked on the implementation of triangular
  and tet-elements. 
- \b Chris \b Gold implemented explicit timestepping
  schemes.
- \b Phil \b Haines worked on bifurcation detection and tracking for
   the Navier-Stokes equations and developed the formulation of
   the equations in plane polar coordinates.
- \b Richard \b Muddle worked on the block preconditioning techniques for
  the biharmonic (and many other) equations, and parallel solvers. 
- \b Glyn \b Rees worked on iterative linear solvers and multigrid
- \b Alberto \b de \b Lozar worked on 3D free-surface Navier-Stokes problems.
- \b Jonathan \b Boyle developed the initial interfaces to third-party 
  iterative solvers and is now involved in the further parallelisation of the
  code and the implementation and application of block-preconditioning 
  techniques for Navier-Stokes and fluid-structure interaction problems.
- \b Renaud \b Schleck completed the octree-based mesh refinement procedures
  and wrote the MPI-based parallel assembly routines and the
  interfaces to SuperLU_dist. 
- \b Sharaf \b Al-Sharif provided the initial implementation of nodal
  spectral elements.
- \b Daniel \b Meyer used oomph-lib to study a variety of axisymmetric 
  Navier-Stokes problems, with and without free surfaces, and developed
  drafts for many of our tutorials.
- \b Alexandre \b Klimowicz worked on block-preconditioning methods.
- \b Jean-Michel \b Lenoir implemented the first part of the
  octree-based 3D mesh refinement procedures.
- \b Gemma \b Barson provided the initial implementation for the 2D Delaunay 
  mesh generation procedures.

We're always looking for more help! Get in touch if you're interested
in joining the team.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Example codes and Tutorials

This document provides a complete list of the example codes that are 
distributed with the \c oomph-lib library. For each
code we give a brief description of the problem solved and
provide a link to the detailed documentation. The
codes are listed in order of increasing complexity. The bullet-point
list in the right column lists the new \c oomph-lib features
that are introduced in the example. You may either work
through the examples one-by-one, treating the example codes and their
documentation as chapters in a self-study course, or use the list of 
topics in the right column as a quick reference to example codes
that provide an introduction to a specific feature.

You may also wish to consult the following documents:
- The <A HREF="../../intro/html/index.html">Introduction </A> provides
  a "low tech" review of the theory of finite elements and presents a 
  "top-down" discussion of the method's implementation within \c oomph-lib.
- The document <A HREF="../../the_data_structure/html/index.html">
  The Data Structure </A> provides a "bottom-up" discussion of 
  \c oomph-lib's data structure.
- The <A HREF="../../quick_guide/html/index.html">(Not-So-)Quick
  Guide</A> provides a "quick" introduction on how to create new
  instances of \c oomph-lib's fundamental objects: \c Problems,
  \c Meshes, and \c Elements. 
- The <A HREF="../../mpi/general_mpi/html/index.html">
  tutorial discussing <code>oomph-lib</code>'s parallel processing
  capabilities.</a>.
.

<CENTER>
<TABLE class="oomph-table">
<TR>
<TD class="panel panel-default">
<CENTER class="panel-heading">
<B>Work in progress</B> 
</CENTER>

<div class="panel-body">
We're still working on the detailed documentation for many
of the demo problems listed below. The fully-documented 
demo problems are accessible via the links. If you
are particularly interested in a specific problem for which the
detailed documentation is incomplete, let us know -- we might be able
to give it a slightly higher priority. We are happy to let you have 
driver codes before the documentation is complete. Such codes
usually need a bit of tidying to make them acceptable for 
"general release", but they are fully functional. In fact, they
are run on a regular basis as part of <code> oomph-lib's </code> self-test
routines (activated by typing <code> make check </code> in the top-level
directory). 
</div>
</TD>
</TR>
</TABLE>
</CENTER>

<H1 NAME="overview">Overview:</H1>
-# <A HREF="#problems">Example codes for specific problem/equations</A>
  - Single-physics problems:
    - <A HREF="#poisson">Poisson problems</A>
    - <A HREF="#poisson_adapt">Mesh adaptation illustrated for 
               Poisson problems</A>
    - <A HREF="#adv_diff">The advection-diffusion equation</A>
    - <A HREF="#unsteady_heat">The unsteady heat equation; with an
               introduction to timestepping</A>
    - <A HREF="#wave">The linear wave equation</A>
    - <A HREF="#helmholtz">The Helmholtz equation</A>
    - <A HREF="#fourier_decomposed_helmholtz"> The azimuthally 
               Fourier-decomposed 3D Helmholtz equation </A>
    - <A HREF="#pml_helmholtz">The Helmholtz equation and
               perfectly matched layers (PMLs)</A>
    - <A HREF="#pml_fourier_decomposed_helmholtz">The azimuthally 
               Fourier-decomposed 3D Helmholtz equation and
               perfectly matched layers (PMLs)</A>
    - <A HREF="#young_laplace">The Young-Laplace equations</A>
    - <A HREF="#nst">The Navier-Stokes equations</A>
    - <A HREF="#axisym_nst">The axisymmetric Navier-Stokes equations</A>
    - <A HREF="#free_surface_nst">The free-surface Navier-Stokes equations</A>
    - <A HREF="#axi_free_surface_nst">The axisymmetric free-surface Navier-Stokes equations</A>
    - <A HREF="#solid">Nonlinear solid mechanics problems</A>
    - <A HREF="#linear_elasticity">Linear elasticity</A>
    - <A HREF="#axisym_linear_elasticity">Axisymmetric linear
               elasticity</A>
    - <A HREF="#time_harmonic_lin_elas">Time-harmonic linear elasticity</A>
    - <A HREF="#gen_time_harmonic_lin_elas">Generalised time-harmonic
               linear elasticity and perfectly matched layers (PMLs)</A>
    - <A HREF="#fourier_decomp_lin_elas"> Azimuthally Fourier-decomposed,
                3D time-harmonic linear elasticity </A>
    - <A HREF="#beam">Beam structures</A>
    - <A HREF="#shell">Shell structures</A>
    .
  - Multi-physics problems:
    - <A HREF="#fsi">Large-displacement fluid-structure interaction problems</A>
    - <A HREF="#acoustic_fsi">Acoustic fluid-structure interaction problems</A>
    - <A HREF="#fd_acoustic_fsi">Azimuthally Fourier-decomposed 3D
               acoustic fluid-structure interaction problems</A>
    - <A HREF="#multi"> Multi-physics problems</A>
    .
  - Eigenproblems:
   - <A HREF="#eigen"> How to formulate and solve an eigenproblem</A>
   .
  .
-# <A HREF="#meshes">Mesh generation</A>
  - <A HREF="#available_meshes">Structured meshes</A>
  - <A HREF="#third_party_meshes">Unstructured meshes generated via input from
             third-party mesh generators</A>
  .
-# <A HREF="#solvers">Linear solvers and preconditioners</A>
  - <A HREF="#linear_solvers">Direct and iterative 
    linear solvers and general-purpose 
    preconditioners</A>
  - <A HREF="#lin_alg">(Distributed) linear algebra</A>
  - <A HREF="#specific_preconditioners">Problem-specific preconditioners</A>
  .
-# <A HREF="#visualisation">Visualisation</A>
  - <A HREF="#paraview">Visualising oomph-lib's output files with Paraview</A>
  .
-# <A HREF="#parallel">Parallel driver codes</A>
  - <A HREF="#distributed">Distributed problems</A>
  .
.

\latexonly
\begin{center}
\bf
[Sorry -- the table listing the various tutorials cannot be rendered properly i#n
latex; please consult the html-based documentation.]
\end{center}
\endlatexonly

<A NAME="problems"></A>


  \htmlonly
<div class="panel panel-default">
  <div class="panel-heading">
    <center><b>How to use example code list</b></center>
  </div>
  <div class="panel-body" >
\endhtmlonly
<TABLE class="example-table">
<TR>
<TD>
<A HREF="#"><B>Problem solved by example code</B></A>

Short description of problem.
</TD>
<TD>
- \c oomph-lib features/conventions illustrated 
by the example code.
</TD>
</TR>
</TABLE>
\htmlonly
</div>
</div>
\endhtmlonly
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="poisson"><H2 class="header">Poisson problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../poisson/one_d_poisson/html/index.html"><B>
The 1D Poisson equation
</B></A>

We (re-)solve the problem considered in the 
<A HREF="../../quick_guide/html/index.html">
Quick Guide</A>, this time using existing \c oomph-lib objects:
The \c OneDMesh and finite elements from the \c QPoisson family.
</TD>
<TD>
- General post-processing routines. 
- General conventions:
  - Use of public typedefs to specify function pointers.
  - Element constructors should not have any arguments.
  .
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../poisson/two_d_poisson/html/index.html"><B>
The 2D Poisson equation.
</B></A>

We solve a 2D Poisson problem with Dirichlet boundary conditions
and compare the results against an exact solution.
</TD>
<TD>
- How to apply Dirichlet boundary conditions in complex meshes.
- How to use the \c DocInfo object to label output files.
- How does one change the linear solver in the Newton method?
- General conventions:
  - \c oomph-lib mesh objects are templated by the element type.
    How does one pre-compile mesh objects?
</TD>


</TR>
<TR>
<TD>
<A HREF="../../poisson/two_d_poisson_flux_bc/html/index.html"><B>
The 2D Poisson equation with flux boundary conditions (I)
</B></A>

Another 2D Poisson problem -- this time with Dirichlet and Neumann
boundary conditions.
</TD>
<TD>
- How to apply non-Dirichlet (flux) boundary conditions with \c
  FaceElements.
- General conventions:
  - What are broken virtual functions and why/when/where are they used?
</TD>

</TR>
<TR>
<TD>
<A HREF="../../poisson/two_d_poisson_flux_bc2/html/index.html"><B>
The 2D Poisson equation with flux boundary conditions (II)
</B></A>

An alternative solution for the previous problem, using multiple meshes.
</TD>
<TD>
- How to use multiple meshes. 
- General conventions:
  - In problems with multiple sub-meshes, \c Nodes retain 
    the boundary numbers of the mesh in which they were created.


- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/two_d_poisson_flux_bc_adapt/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="poisson_adapt"><H2 class="header">Poisson problems with adaptivity</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../poisson/fish_poisson/html/index.html"><B>
Adaptive solution of Poisson's equation in a fish-shaped domain
</B></A> 

We solve a 2D Poisson equation in a nontrivial, fish-shaped domain and 
demonstrate \c oomph-lib's fully-automatic mesh adaptation 
routines.
</TD>
<TD>
- How to perform automatic mesh adaptation. 
- \c oomph-lib's "black-box" adaptive Newton solver.
- \c The functions \c Problem::actions_before_adapt() and 
  \c Problem::actions_after_adapt().
</TD>

</TR>
<TR>
<TD>
<A HREF="../../poisson/two_d_poisson_adapt/html/index.html"><B>
The 2D Poisson equation revisited -- how to create a refineable mesh
</B></A> 

We revisit an <A HREF="../../poisson/two_d_poisson/html/index.html">
earlier example</A> and demonstrate how easy it is to "upgrade" an
existing mesh to a mesh that can be used with \c oomph-lib's
automatic mesh adaptation routines.
</TD>
<TD>
- "Upgrading" meshes to make them refineable.
- General conventions:
  - Hanging nodes -- the functions \c Node::position(...)
    and \c Node::value(...).
</TD>

</TR>
<TR>
<TD>
<A HREF="../../poisson/fish_poisson2/html/index.html">
<B>
Poisson's equation in a fish-shaped domain revisited -- mesh
adaptation in deformable domains with curvilinear and/or moving boundaries.
</B>
</A>

We revisit an <A HREF="../../poisson/fish_poisson/html/index.html">
earlier example</A> and demonstrate how to create refineable
meshes for problems with curvilinear and/or moving domain boundaries.
</TD>
<TD>
- How to create refineable meshes for problems with
  curvilinear and/or moving domain boundaries.
- General conventions:
  - The \c GeomObject, \c Domain and \c MacroElement objects.
  - The \c Mesh::node_update() function. 
  - It is good practice to store boundary coordinates for 
    (\c Boundary) \c Nodes that are located on curvilinear domain boundaries. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
<B>
Adaptive solution of Poisson's equation with flux boundary
conditions. 
</B>
</A>

We revisit an <A HREF="../../poisson/two_d_poisson_flux_bc2/html/index.html">
earlier example</A> and demonstrate how to apply flux boundary
conditions in problems with spatial adaptivity. 
</TD>
<TD>
- How to apply flux boundary conditions in problems with spatial 
  adaptivity.
- General conventions:
  - The \c Mesh::flush_element_and_node_storage() function: "Emptying"
    a mesh without deleting its constituent nodes and elements (they
    might be shared with other meshes!).


- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/two_d_poisson_flux_bc_adapt/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 

</TR>
<TR>
<TD>
<A HREF="../../poisson/eighth_sphere_poisson/html/index.html">
<B>
Adaptive solution of a 3D Poisson equations in a spherical domain
</B>
</A>

We demonstrate \c oomph-lib's octree-based 3D mesh adaptation
routines.
</TD>
<TD>
- Setting up and solving 3D problems isn't any harder than 
  doing it in 2D. 
- General conventions:
  - The namespace \c CommandLineArgs provides storage for the
    command line arguments to make them accessible throughout the
    code.
  - Plotting mesh boundaries.
</TD>

</TR>

</TABLE>
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="adv_diff"><H2 class="header">The advection-diffusion equation</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../advection_diffusion/two_d_adv_diff_adapt/html/index.html">
<B>The 2D advection diffusion equation with spatial adaptivity
</B>
</A>

We solve a 2D advection-diffusion equation and illustrate the
characteristic features of solutions at large Peclet number.
</TD>
<TD>
- The adaptive discretisation of the advection diffusion equation.
- How to specify the "wind" and the Peclet number for the 
  advection diffusion equation.
- How to document the progress of \c oomph-lib's "black box" adaptive
  Newton solver.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../advection_diffusion/two_d_adv_diff_flux_bc/html/index.html">
<B>2D advection diffusion equation with Neumann (flux) boundary conditions.
</B>
</A>

We solve a 2D advection-diffusion equation with flux boundary conditions.
</TD>
<TD>
- How to specify Neumann (flux) boundary conditions for the
  advection diffusion equation.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../advection_diffusion/two_d_adv_diff_SUPG/html/index.html">
<B>The 2D advection diffusion equation revisited: Petrov-Galerkin
methods and SUPG stabilisation.
</B>
</A>

We demonstrate how to implement a stabilised Petrov-Galerkin discretisation
of the advection diffusion equation. 
</TD>
<TD>
- Petrov-Galerkin discretisation of the advection-diffusion equation.
- General conventions:
  - The role of shape, basis and test functions.
  - The element building blocks: Geometric elements, equation classes
    and specific elements.
- <span style="color:red">Note: the driver code is unannotated.</span>
</TD>
</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="unsteady_heat"><H2 class="header">The unsteady heat equation and 
an introduction to time-stepping</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat/html/index.html">
<B>
The 2D unsteady heat equation
</B>
</A>

We solve the 2D unsteady heat equation and 
demonstrate \c oomph-lib's time-stepping procedures for 
parabolic problems.
</TD>
<TD>
- Solving time-dependent problems with \c Problem::unsteady_newton_solve(...)
- The functions \c Problem::actions_before_implicit_timestep() and
  \c Problem::actions_before_implicit_timestep().
- The BDF timesteppers and how to set up initial conditions 
  for parabolic problems.
- Initialising the "previous" nodal positions for elements that are
  based on an ALE formulation. 
- General conventions:
  - Providing a default \c Steady timestepper for all \c Mesh
    constructors. 
  - Steady and unsteady versions of functions -- position and
    interpretation of the (discrete) time index.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat2/html/index.html">
<B>
The 2D unsteady heat equation with restarts
</B>
</A>

We demonstrate \c oomph-lib's dump/restart capabilities which
allow time-dependent simulations to be restarted.
</TD>
<TD>
- The functions \c Problem::dump(...) and
  \c Problem::read(...).
- How to customise the generic dump/restart functions.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat_t_adapt/html/index.html">
<B>
The 2D unsteady heat equation with adaptive timestepping
</B>
</A>

We demonstrate \c oomph-lib's adaptive timestepping capabilities.
</TD>
<TD>
- How to enable temporal adaptivity.
- The function \c Problem::adaptive_unsteady_newton_solve().
- The function \c Problem::global_temporal_error_norm().
- How to choose the target for the global temporal error norm.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">
<B>
Spatially adaptive solution of the 2D unsteady heat equation with
Neumann (flux) boundary conditions.
</B>
</A>

We solve a 2D unsteady heat equation in a non-trivial domain
with flux boundary conditions and compare the computed results
against the exact solution.
</TD>
<TD>
- Spatial adaptivity in time-dependent problems.
- Choosing the maximum number of spatial adaptations per timestep.
- Using \c Problem::set_initial_condition() to assign initial
  conditions ensures that the initial conditions are re-assigned when
  mesh adaptations are performed while the first timestep is
  computed.
- The functions \c Problem::dump(...) and
  \c Problem::read(...) can handle adaptive meshes. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat_ALE/html/index.html">
<B>
Spatially adaptive solution of the 2D unsteady heat equation in a
moving domain with Neumann (flux) boundary conditions.
</B>
</A>

We demonstrate the  spatially adaptive solution of a 2D
unsteady heat equation in a nontrivial moving domain.
</TD>
<TD>
- The ALE form of the unsteady heat equation and its implementation in
  \c oomph-lib's unsteady heat elements.
- The role of the positional \c TimeStepper and the importance of
  assigning history values for the nodal positions. 
- General conventions:
  - The function \c Mesh::node_update() automatically updates 
    the nodal positions in response to the deformation/motion of time-dependent
    \c GeomObjects that define the \c Domain and \c Mesh boundaries.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../unsteady_heat/two_d_unsteady_heat_2adapt/html/index.html">
<B>
Spatially and temporally adaptive solution of the 2D unsteady heat
equation in a moving domain with flux boundary conditions.
</B>
</A>

We demonstrate the use of combined spatial and temporal adaptivity for the
solution of a 2D unsteady heat equation in a nontrivial moving domain.
</TD>
<TD>
- Adaptive timestepping combined with spatial adaptivity.
</TD>

</TR>

</TABLE>
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="wave"><H2 class="header">The linear wave equation</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../linear_wave/two_d_linear_wave/html/index.html">
<B>The 2D linear wave equation.
</B>
</A>

We solve a 2D wave equation and demonstrate \c oomph-lib's 
time-stepping procedures for hyperbolic problems.
</TD>
<TD>
- Timestepping for hyperbolic problems: The Newmark scheme.
- How to set up initial conditions for hyperbolic problems.
- Default settings for the linear wave equation elements. 
- How to apply Neumann (flux) boundary conditions for the
  linear wave equation. 
</TD>

</TR>

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="helmholtz"><H2 class="header">The Helmholtz equation</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../helmholtz/scattering/html/index.html">
<B>The Helmholtz equation.
</B>
</A>

We solve a 2D Helmholtz problem, simulating scattering of a planar
wave from a circular cylinder. 
</TD>
<TD>
- The Helmholtz equation and its discretisation.
- The Sommerfeld radiation condition and its approximation by
  approximate/absorbing boundary conditions and Dirichlet-to-Neumann mappings.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../helmholtz/unstructured_scattering/html/index.html">
<B>Solving the Helmholtz equation on an unstructured mesh.
</B>
</A>

We solve a 2D Helmholtz problem, simulating scattering of a planar
wave from a circular cylinder -- this time using an unstructured mesh.
</TD>
<TD>
- How to solve the Helmholtz equation on an unstructured mesh.
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="fourier_decomposed_helmholtz"><H2 class="header">The
 azimuthally Fourier-decomposed 3D Helmholtz equation</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../fourier_decomposed_helmholtz/sphere_scattering/html/index.html">
<B>The azimuthally Fourier-decomposed 3D Helmholtz equation.
</B>
</A>

We solve the 3D Helmholtz equation in cylindrical polar coordinates,
  using a Fourier-decomposition in the azimuthal direction.
</TD>
<TD>
- The Helmholtz equation and its azimuthal Fourier decomposition.
- Dirichlet to Neumann mapping.
- Validation with an exact solution formed by a superposition of
  outgoing waves from a sphere.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../fourier_decomposed_helmholtz/adaptive_sphere_scattering/html/index.html">
<B>The spatially adaptive solution of the azimuthally Fourier-decomposed 3D Helmholtz equation on unstructured meshes.
</B>
</A>

We re-visit the 3D Helmholtz equation in cylindrical polar coordinates,
using a Fourier-decomposition in the azimuthal direction -- this
time using spatial adaptivity and an unstructured mesh
</TD>
<TD>
- Spatial adaptivity on unstructured meshes for the azimuthally 
  Fourier-decomposed 3D Helmholtz equation.
</TD>

</TR>

</TABLE>
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="pml_helmholtz"><H2 class="header">The Helmholtz equation and perfectly 
matched layers</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../pml_helmholtz/scattering/html/index.html">
<B>The Helmholtz equation and perfectly matched layers (PMLs).
</B>
</A>
We demonstrate the imposition of the Sommerfeld radiation condition
by perfectly matched layers (PMLs) using the example of a radiating cylinder.
</TD>
<TD>
- The Helmholtz equation.
- The Sommerfeld radiation condition and its approximation by
  perfectly matched layers (PMLs).
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="pml_fourier_decomposed_helmholtz">
<H2 class="header">The azimuthally 
Fourier-decomposed 3D Helmholtz equation and
perfectly matched layers (PMLs)</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../pml_fourier_decomposed_helmholtz/oscillating_sphere/html/index.html">
<B>The azimuthally 
Fourier-decomposed 3D Helmholtz equation and perfectly matched layers (PMLs).
</B>
</A>

We consider the azimuthally 
Fourier-decomposed 3D Helmholtz equation and
demonstrate the imposition of the Sommerfeld radiation condition
by perfectly matched layers (PMLs).
</TD>
<TD>
- The azimuthally 
  Fourier-decomposed 3D Helmholtz equation.
- The Sommerfeld radiation condition and its approximation by
  perfectly matched layers (PMLs).
</TD>

</TR>



<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="young_laplace"><H2 class="header">The Young Laplace equation</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../young_laplace/young_laplace/html/index.html">
<B>The solution of the Young-Laplace equation.
</B>
</A>

We solve the Young Laplace equation that governs the shape of
static air-liquid interfaces. 
</TD>
<TD>
- Theory and implementation.
- The use of spines for the representation of complicated
  air-liquid interfaces.
- Natural boundary conditions along free contact lines.
- How to use displacement control to compute (past) limit
  points in the load-displacement characteristics.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../young_laplace/contact_angle/html/index.html">
<B>
Contact-angle boundary conditions for the Young-Laplace equation
</B>
</A>

We demonstrate how to apply contact angle-boundary conditions for the
Young-Laplace equation. 
</TD>
<TD>
- Theory and implementation of contact angle-boundary conditions.
- How to generate initial guesses for the solution
- Limitations of the current approach and suggestions for improvement
- An inherent difficulty in problems with zero contact angles.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="nst"><H2 class="header">The Navier-Stokes equations</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/driven_cavity/html/index.html">
<B>The 2D Navier-Stokes equations: Driven cavity flow
</B>
</A>

Probably the most-solved problem in computational fluid dynamics:
Steady driven cavity flow. We illustrate the problem's discretisation
with Taylor-Hood and Crouzeix-Raviart elements.
</TD>
<TD>
- Discretising the steady Navier-Stokes equations: The governing
  equations and their implementation in the stress-divergence and
  simplified forms.
- Non-dimensional parameters and their default values.
- The pressure representation in Taylor-Hood and Crouzeix-Raviart elements.
- Pinning a pressure value in problems with Dirichlet boundary
  conditions for the velocity on all boundaries.      

</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/adaptive_driven_cavity/html/index.html">
<B>The 2D Navier-Stokes equations: Adaptive solution of the 2D driven cavity problem
</B>
</A>

We employ \c oomph-lib's mesh adaptation routines 
to refine the mesh in the neighbourhood of the 
pressure singularities. 
</TD>
<TD>
- Treatment of pressure degrees of freedom in Navier-Stokes 
  simulations with adaptive mesh refinement -- pinning "redundant"
  pressure degrees of freedom.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/adaptive_driven_cavity/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span>
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/circular_driven_cavity/html/index.html">
<B>The 2D Navier-Stokes equations: Driven cavity flow in a
quarter-circle domain with mesh adaptation
</B>
</A>

We re-solve the driven-cavity problem in a different domain,
demonstrate how to apply body forces and show how to switch
between the stress-divergence and simplified forms of the
Navier-Stokes equations.
</TD>
<TD>
- Adapting the driven cavity problem to different domains.
- How to apply body forces in the Navier-Stokes equations.
- How to switch between the stress-divergence and the simplified
  form of the incompressible Navier-Stokes equations.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/three_d_entry_flow/html/index.html">
<B>Adaptive simulation of 3D finite Reynolds number entry flow into a circular pipe
</B>
</A>

We solve the classical problem of entry flow into a 3D tube.
</TD>
<TD>
- Adaptivity for 3D Navier-Stokes problems
- How to determine the numbering scheme for mesh boundaries.
- How to adjust parameters that control the behaviour of
  \c oomph-lib's Newton solver.
- The natural (traction-free) boundary conditions for the
  Navier-Stokes equations.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/rayleigh_channel/html/index.html">
<B>A variant of Rayleigh's oscillating plate problem: The unsteady 
2D Navier-Stokes equations  with periodic boundary conditions
</B>
</A>

We solve a variant of the classical Rayleigh plate problem to demonstrate
the use of periodic boundary conditions and time-stepping for the
Navier-Stokes equations.
</TD>
<TD>
- Timestepping for the Navier-Stokes equations. 
- How to apply periodic boundary conditions.
- General conventions:
  - Periodic boundary conditions should be applied in the \c Mesh
    constructor.
  .
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/rayleigh_traction_channel/html/index.html">
<B>Another variant of Rayleigh's oscillating plate problem: The unsteady 
2D Navier-Stokes equations with periodic boundary conditions, driven
by an applied traction.
</B>
</A>

We demonstrate how to apply traction boundary conditions for the
Navier-Stokes equations. 
</TD>
<TD>
- How to apply traction boundary conditions for the Navier-Stokes
  equations.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/osc_ellipse/html/index.html">
<B> 2D finite-Reynolds-number-flow driven by an oscillating ellipse
</B>
</A>

We study the 2D finite-Reynolds number flow contained inside an 
oscillating elliptical ring and compare the computed results against
an exact solution (an unsteady stagnation point flow).
</TD>
<TD>
- Solving the Navier-Stokes equations in a moving domain. 
- How to apply no-slip boundary conditions on moving walls, using
  the function
  \c FSI_functions::apply_no_slip_on_moving_wall(...)
</TD>

</TR>

</TABLE>
<TABLE class="example-table">
<TR>
<TD>
<A HREF="../../navier_stokes/collapsible_channel/html/index.html">
<B> 2D finite-Reynolds-number-flow in a 2D channel with a moving wall
</B>
</A>

This is a "warm-up" problem for the classical fluid-structure
interaction problem of 
<A HREF="../../interaction/fsi_collapsible_channel/html/index.html"> 
flow in a 2D collapsible channel</A>. Here we
compute the flow through a 2D channel in which part of one wall is 
replaced by a moving "membrane" whose motion is prescribed. 
</TD>
<TD>
- The adaptive solution of the Navier-Stokes equations in a moving 
  domain with traction boundary conditions. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/algebraic_collapsible_channel/html/index.html">
<B> 2D finite-Reynolds-number-flow in a 2D channel with a moving wall
revisited: Algebraic Node updates.
</B>
</A>

We re-visit the problem studied in the previous example and
demonstrate an alternative node-update procedure, based on
\c oomph-lib's \c AlgebraicNode, \c AlgebraicElement and
\c AlgebraicMesh classes. Algebraic node updates will turn out to be essential
for the efficient implementation of fluid-structure interaction 
problems.
</TD>
<TD>
- How to customise the node-update, using \c oomph-lib's \c AlgebraicNode, 
\c AlgebraicElement and \c AlgebraicMesh classes. 
- Existing \c AlgebraicMeshes are easy to use: Simply "upgrade" the
  required element (of type \c ELEMENT, say) in the templated wrapper class
  \c AlgebraicElement<ELEMENT>. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/spine_channel/html/index.html">
<B> Steady 2D finite-Reynolds-number-flow in a 2D channel of
  non-uniform width: An Introduction to Spine meshes
</B>
</A>

We consider a variation of the problem studied in the previous example and
demonstrate an alternative node-update procedure, based on
\c oomph-lib's \c SpineNode, \c SpineElement and
\c SpineMesh classes. Spine node updates are an alternative to
  algebraic nodes updates for the efficient implementation of
  fluid problems in deforming domains.
</TD>
<TD>
- Introduction to the method of spines.
- How to create a custom \c SpineMesh.

</TD>

</TR>

<TR>
<TD>
<A HREF="../../navier_stokes/channel_with_leaflet/html/index.html">
<B> 2D finite-Reynolds-number-flow in a 2D channel that is partially
obstructed by an oscillating leaflet
</B>
</A>

This is a "warm-up" problem for 
<A HREF="../../interaction/fsi_channel_with_leaflet/html/index.html"> 
the corresponding fluid-structure interaction problem</A> where
the leaflet deforms in response to the fluid traction. Here we
consider the case where the motion of the leaflet is prescribed.
</TD>
<TD>
- Another example illustrating the use of algebraic and 
  \c MacroElement/Domain-based node update techniques.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/turek_flag_non_fsi/html/index.html">
<B> Flow past a cylinder with a waving flag
</B>
</A>

This is a "warm-up" problem for 
<A HREF="../../interaction/turek_flag/html/index.html"> 
Turek & Hron's FSI benchmark problem</A> where
the flag deforms in response to the fluid traction. Here we
consider the case where the motion of the flag is prescribed.
</TD>
<TD>
- Another example illustrating the use of algebraic and 
  \c MacroElement/Domain-based node update techniques.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">
<TR>
<TD>
<A HREF="../../navier_stokes/unstructured_fluid/html/index.html">
<B> Unstructured meshes for fluids problems
</B>
</A>

This is a "warm-up" problem for 
<A HREF="../../interaction/unstructured_fsi/html/index.html"> 
another tutorial</A> in which we demonstrate the use of
unstructured meshes for FSI problems.
</TD>
<TD>
- How to use xfig/triangle-generated, unstructured meshes
  for flow problems.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/unstructured_three_d_fluid/html/index.html">
<B> Unstructured meshes for 3D fluids problems
</B>
</A>

This is a "warm-up" problem for 
<A HREF="../../interaction/unstructured_three_d_fsi/html/index.html"> 
another tutorial</A> in which we demonstrate the use of
unstructured 3D meshes for FSI problems.
</TD>
<TD>
- How to use tetgen-generated, unstructured meshes
  for 3D flow problems.
- Avoiding "locking" with the \c split_corner_elements flag.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/vmtk_fluid/html/index.html">
<B> Steady finite-Reynolds-number flow through an iliac bifurcation
</B>
</A>

We show how to simulate physiological flow problems, using
the  <a href="http://www.vmtk.org">Vascular 
Modeling Toolkit (VMTK).</a>
This is a "warm-up" problem for 
<A HREF="../../interaction/vmtk_fsi/html/index.html"> 
another tutorial</A> in which we consider the corresponding FSI 
problems in which the vessel wall is elastic.
</TD>
<TD>
- How to use \c ImposeParallelOutflowElements to enforce parallel 
  in- and outflow from cross-sections that are not aligned with 
  any coordinate planes.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/jeffery_orbit/html/index.html">
<B>Motion of elliptical particles in shear flow: unstructured adaptivity</B>
</A>

We solve the classical problem of shear flow past a immersed ellipse
</TD>
<TD>
- An example of using inline mesh generation for adaptivity in
  unstructured meshes
- \c ImmersedRigidBodyElements to describe the interaction of fluids with rigid bodies
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/curved_pipe/html/index.html">
<B>Adaptive simulation of flow at finite Reynolds number in a curved circular pipe
</B>
</A>

We solve the classical problem of flow into a 3D curved tube.
</TD>
<TD>
- Another example of adaptivity for 3D Navier-Stokes problems
- <span style="color:red">Note: the documentation for this problem is incomplete.</span>
</TD>
</TR>


<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="axisym_nst"><H2 class="header">The axisymmetric Navier-Stokes equations</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../axisym_navier_stokes/spin_up/html/index.html">
<B>Spin-up of a viscous fluid -- the spatially adaptive solution of the unsteady axisymmetric Navier-Stokes equations.</B>
</A>

A classical fluid mechanics problem: Spin-up of a viscous fluid.
A key feature of the flow is the development of thin Ekman (boundary)
layers during the early stages of the spin-up. We demonstrate how
the use of spatial adaptivity helps to resolve these layers. At 
large times, the flow field approaches a rigid-body rotation -- 
this poses a subtle problem for the spatial adaptivity as its default 
behaviour would cause strong spatially uniform refinement.
</TD>
<TD>
- The axisymmetric Navier-Stokes equations
- How to prescribe a constant reference value for the normalisation of
  the error in spatially-adaptive computations in which the solution
  approaches a "trivial" solution.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="free_surface_nst"><H2 class="header">The free-surface Navier-Stokes equations</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/surface_theory/html/index.html">
<B>Interfaces, Free Surfaces and Surface Transport: Theory and Implementation.</B>
</A>
In this document, we introduce the basic theory of moving surfaces,
surface calculus and surface transport equations.
In addition, we describe how \c oomph-lib's free-surface and
surface-transport capabilities are implemented.
</TD>
<TD>
- Theory
 - Geometry of Surfaces
 - Differential Operators On A Surface
 - Free Surface and Interface Boundary Conditions
 - Surface Transport Equations
 .
- Implementation
 - The FluidInterfaceElement class
 - Spine and Elastic formulations of free surface elements
</TD>
</TR>


<TR>
<TD>
<A HREF="../../navier_stokes/single_layer_free_surface/html/index.html">
<B>Free-surface relaxation oscillations of a viscous fluid layer.</B>
</A>

We study the oscillations of a perturbed fluid layer and compare the
results to the
analytic dispersion relation based on linearised disturbances.
</TD>
<TD>
- Boundary conditions at a free surface
- How to solve free-boundary problems using a pseudo-solid node update approach.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/two_layer_interface/html/index.html">
<B>Relaxation oscillations of an interface between two viscous fluids.</B>
</A>

We study the oscillations of a two-layer fluid system and compare the
results to the analytic dispersion relation based on linearised disturbances.
</TD>
<TD>
- Boundary conditions at an interface between two immiscible fluids.
- Adaptivity in free-surface/interfacial problems with a 
  pseudo-solid node update approach.
</TD>

</TR>


<TR>
<TD>
<A HREF="../../navier_stokes/static_single_layer/html/index.html">
<B>A static free surface bounding a layer of viscous fluid.</B>
</A>

A hydrostatics problem: Compute the static free surface that
bounds a layer of viscous fluid -- harder than you might think!
</TD>
<TD>
 - Imposing a volume constraint in steady free-surface problems.
 - Imposing static contact angle constraints when a free-surface meets
 a solid boundary.
 - Hijacking (overwriting) data from within elements.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../navier_stokes/static_two_layer/html/index.html">
<B>A static interface between two viscous fluids.</B>
</A>

A hydrostatics problem: Compute the static interface between two
viscous fluids -- harder than you might think!
</TD>
<TD>
- Adapting an existing mesh to include an interface.
- Imposing a static contact angle when an interface meets a solid boundary.
</TD>

</TR>

<TR>
<TD>
<A HREF="../../navier_stokes/inclined_plane/html/index.html">
<B>Flow of a fluid film down an inclined plane</B>
</A>

A classical fluid mechanics problem: We study the flow of a film of fluid 
down an inclined plane and compare the results to the exact solution
of Nusselt. The stability is assessed by simulating the time-evolution
of a perturbation to the free surface.
</TD>
<TD>
- Comparison between spine-based and pseudo-solid
node update strategies for free surface problems.
- Applying time-dependent perturbations to steady solutions.
- Open flow boundary conditions.
</TD>

</TR>

<!-- <TR>
<TD>
<A HREF="../../navier_stokes/adaptive_interface/html/index.html">
<B>A steadily-rotating cylinder below a free surface in a finite box.</B>
</A>

Compute the free surface position and fluid velocity and pressure
fields about a fixed, steadily-rotating cylinder immersed in a viscous
fluid in a finite box. 
Uses a pseudo-elastic remesh strategy and spatial adaptivity in a 
non-trivial free surface problem.
</TD>
<TD>
- ***
- General conventions:
  - ***
<span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR> -->


</TABLE>
<TABLE class="example-table">

<TR>
<TD>
<A HREF="../../navier_stokes/bretherton/html/index.html">
<B>The Bretherton problem: An air finger propagates into a 2D
fluid-filled channel.</B>
</A>

A classical fluid mechanics problem: We study the propagation
of an inviscid (air) finger into a 2D fluid-filled channel and
compare our results against those from Bretherton's theoretical
analysis.
</TD>
<TD>
- Another example for the solution of a free-surface problem using
  the Method of Spines.
- Use of external \c Data to implement non-local interactions
  between elements.
</TD>

</TR>


<TR>
<TD>
<A HREF="../../navier_stokes/adaptive_bubble_in_channel/html/index.html">
<B>A finite bubble propagates in a 2D fluid-filled channel.</B>
</A>

Uses a pseudo-elastic remesh strategy and unstructured spatial adaptivity in a 
non-trivial free surface problem.
</TD>
<TD>
- Mesh adaptation for closed free boundaries -- what happens "under
  the hood" and how to customise the default behaviour, specifically:
  - Selecting the refinement and unrefinement tolerances on
    the free boundaries.
  - Redistributing segments between polylines in cases when the
    nodes on the free surface all get convected into one region.
  - Why you shouldn't use \c TriangleMeshCurviLines to describe
    the initial shape of (initially curvilinear) free boundaries.
  .
- How to impose a volume constraint onto an enclosed "bubble" in 
  a viscous fluid.
.
</TD>

</TR>

<TR>
<TD>
<A HREF="../../navier_stokes/adaptive_droplet_in_channel/html/index.html">
<B>A finite droplet propagates in a 2D fluid-filled channel.</B>
</A>

Uses a pseudo-elastic remesh strategy and unstructured spatial
adaptivity in 
a non-trivial interfacial (two-fluid) problem.
</TD>
<TD>
- How to use regions in an unstructured  two-fluid problem.
- How to impose and then remove a volume constraint onto an enclosed
"droplet" in a viscous fluid.
.
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="axi_free_surface_nst"><H2 class="header">The axisymmetric free-surface Navier-Stokes equations</H2></A>
</TD>
</TR>

<TR>
<TD>
<A HREF="../../axisym_navier_stokes/two_layer_interface_axisym/html/index.html">
<B>Relaxation oscillations of an interface between two viscous fluids in an axisymmetric domain.</B>
</A>

We study the oscillations of a two-layer fluid system in an
axisymmetric domain and compare the
results to the analytic dispersion relation based on linearised disturbances.
</TD>
<TD>
- Free-surface/interfacial problems in an axisymmetric domain with a 
  pseudo-solid node update approach and spatial adaptivity.
</TD>

</TR>

<TR>
<TD>
<A HREF="../../axisym_navier_stokes/axi_static_cap/html/index.html">
<B>An axisymmetric  static free surface bounding a layer of viscous fluid.</B>
</A>

A hydrostatics problem: Compute the static free surface that
bounds a layer of viscous fluid in an axisymmetric geometry.
</TD>
<TD>
 - Trivial differences between two-dimensional and axisymmetric static problems.
</TD>

</TR>

</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="solid"><H2 class="header">Nonlinear solid mechanics problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../solid/solid_theory/html/index.html">
<B>Nonlinear solid mechanics: Theory and implementation
</B>
</A>

In this document we discuss the theoretical background and the
practical implementation of \c oomph-lib's nonlinear solid mechanics
capabilities. 
</TD>
<TD>
- Theory:
 - Nonlinear solid mechanics problems -- Lagrangian coordinates
 - The geometry
 - Equilibrium and the Principle of Virtual Displacements
 - Constitutive Equations for Purely Elastic Behaviour
 - Non-dimensionalisation
 - 2D problems: Plane strain.
 - Isotropic growth.
 - Specialisation to a cartesian basis and finite element discretisation
 . 
- Implementation:
 - The \c SolidNode class
 - The \c SolidFiniteElement class
 - The \c SolidMesh class
 - The \c SolidTractionElement class
 .
- Timestepping and the generation of initial conditions for solid mechanics problems
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/airy_cantilever/html/index.html">
<B>Bending of a cantilever beam
</B>
</A>

We study a classical solid mechanics problem:  the bending
of a cantilever beam subject to a uniform pressure loading on its
upper face and/or gravity. We compare the results for zero-gravity
against the (approximate) analytical St. Venant solution for
the stress field. 
</TD>
<TD>
- How to formulate solid mechanics problems.
- How to choose a constitutive equation.
- How to apply traction boundary conditions, using
  \c SolidTractionElements.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/disk_compression/html/index.html">
<B>Axisymmetric compression of a circular disk
</B>
</A>

We study the axisymmetric compression of a circular, elastic disk,
loaded by an external traction. The results are compared against 
the predictions from small-displacement elasticity. 
</TD>
<TD>
- How to upgrade an existing mesh to \c SolidMesh.
- Why it is necessary to use "undeformed MacroElements" to ensure 
  that the numerical results converge to the correct solution under 
  mesh refinement if the domain has curvilinear boundaries.
- How to switch between different constitutive equations
- how to incorporate isotropic growth into the model.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/compressed_square/html/index.html">
<B>Compressible and incompressible behaviour
</B>
</A>

We discuss various issues related to (in)compressible
material behaviour and illustrate various solution
techniques in a simple test-problem: The compression of
a square block of (compressible or incompressible) material by 
a gravitational body force. The results are compared against 
the predictions from small-displacement elasticity. 
</TD>
<TD>
- Different formulations of the constitutive laws for
  compressible, near-incompressible and incompressible behaviour.
- How to combine the various constitutive laws with the displacement and
  pressure/displacement formulations of the principle of virtual
  displacements. 
- The default setting: Incompressibility is not enforced automatically.
  It must be requested explicitly.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/disk_oscillation/html/index.html">
<B>Axisymmetric oscillations of a circular disk
</B>
</A>

We study the free axisymmetric oscillations of a circular, 
elastic disk and compare the eigenfrequencies and modes
against the predictions from small-displacement elasticity. 
</TD>
<TD>
- How to assign initial conditions for unsteady solid mechanics problems.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/prescribed_displ_lagr_mult/html/index.html">
<B>Deformation of a solid by a prescribed boundary motion
</B>
</A>

We study the large deformations of a 2D elastic domain,
driven by the prescribed deformation of its boundary. 
The boundary motion is imposed by Lagrange multipliers. 
This technique is important for the solution of fluid-structure
interaction problems in which the deformation of the fluid mesh
is controlled by (pseudo-)elasticity.
</TD>
<TD>
- How to impose displacement boundary conditions in solid mechanics
  problems by Lagrange multipliers.
</TD>

</TR>

</TABLE>
<TABLE class="example-table">
<TR>
<TD>
<A HREF="../../solid/three_d_cantilever/html/index.html">
<B>Large-amplitude bending of an asymmetric 3D cantilever beam
made of incompressible material.
</B>
</A>

We study the deformation of an asymmetric 3D cantilever beam
made of incompressible material.
</TD>
<TD>
- How to enforce incompressible behaviour. 
- How to solve 3D solid mechanics problems with spatial adaptivity.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/unstructured_solid/html/index.html">
<B>Unstructured meshes for 2D and 3D solid mechanics problems.
</B>
</A>

We demonstrate how to use unstructured meshes to solve 2D and
3D solid mechanics problems. This tutorial acts as a "warm-up"
problem for the solution of unstructured FSI problems. 
</TD>
<TD>
- How to solve 2D and 3D solid mechanics problems on unstructured 
  meshes. 
- How to identify domain boundaries in xfig-generated unstructured
  meshes. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/unstructured_three_d_solid/html/index.html">
<B>Unstructured meshes for 3D solid mechanics problems.
</B>
</A>

We demonstrate how to use unstructured meshes to solve 
3D solid mechanics problems. This is a "warm-up" problem for 
<A HREF="../../interaction/unstructured_three_d_fsi/html/index.html"> 
another tutorial</A> in which we demonstrate the use of
unstructured 3D meshes for FSI problems.
</TD>
<TD>
- How to solve 3D solid mechanics problems on unstructured 
  meshes. 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/vmtk_solid/html/index.html">
<B>Inflation of a blood vessel
</B>
</A>

We show how to simulate physiological solid mechanics problems, using
the  <a href="http://www.vmtk.org">Vascular 
Modeling Toolkit (VMTK).</a>
This is a "warm-up" problem for 
<A HREF="../../interaction/vmtk_fsi/html/index.html"> 
another tutorial</A> in which we consider the corresponding FSI 
problems in which the vessel conveys (and is loaded by) a viscous
fluid.
</TD>
<TD>
- How to solve physiological solid mechanics problems.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/shock_disk/html/index.html">
<B>Large-amplitude shock waves in a circular disk
</B>
</A>

We study the propagation of shock waves in an elastic 2D circular disk.
</TD>
<TD>
- How to employ spatial adaptivity in time-dependent solid mechanics problems.
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../solid/unstructured_adaptive_solid/html/index.html">
<B>Solid Mechanics using unstructured meshes with adaptivity
</B>
</A>

We study the deflection of a 2D rectangular solid under a lateral
pressure load.
</TD>
<TD>
- How to employ spatial adaptivity in solid mechanics problems using
unstructured meshes.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../solid/simple_shear/html/index.html">
<B>Large shearing deformations of a hyper-elastic, incompressible
block of material 
</B>
</A>

We solve a classical problem in large-displacement elasticity and
compare against Green and Zerna's exact solution. 
</TD>
<TD>
- A validation problem.
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="linear_elasticity"><H2 class="header">Linear elasticity</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../linear_elasticity/periodic_load/html/index.html">
<B>Linear Elasticity: Theory and implementation
</B>
</A>

In this document we discuss the theoretical background and the
practical implementation of \c oomph-lib's linear elasticity
elements and demonstrate their use in a 2D model problem.
</TD>
<TD>
- The equations of linear elasticity and their non-dimensionalisation.
- Implementation of the equations, based on a displacement formulation.
- The solution of a 2D model problem: The deformation of a linearly
  elastic strip, loaded by a spatially periodic surface traction.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../linear_elasticity/refineable_periodic_load/html/index.html">
<B>Spatially-adaptive simulation of the deformation of a linearly
elastic strip, loaded by a spatially periodic surface traction.
</B>
</A>

We demonstrate how to compute the deformation of a linearly
elastic strip, loaded by a spatially periodic surface traction, using
spatial adaptivity.
</TD>
<TD>
- How to apply periodic boundary conditions in spatially adaptive
  computations.
</TD>

</TR>





<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="axisym_linear_elasticity"><H2 class="header">Axisymmetric linear elasticity</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../axisym_linear_elasticity/cylinder/html/index.html">
<B>Axisymmetric linear elasticity: Theory, implementation
and a time-dependent demo problem.
</B>
</A>

In this document we discuss the theoretical background and the
practical implementation of \c oomph-lib's axisymmetric linear elasticity
elements and demonstrate their use in time-dependent test problem.
</TD>
<TD>
- The equations of axisymmetric linear elasticity and their 
  non-dimensionalisation.
- Implementation of the equations, based on a displacement formulation.
- The solution of a 2D model problem: The time-dependent deformation
  of an annular region and the validation against a manufactured
  exact solution.
.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="time_harmonic_lin_elas"><H2 class="header">
Time-harmonic linear elasticity</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../time_harmonic_linear_elasticity/elastic_annulus/html/index.html">
<B>The equations of time-harmonic linear
elasticity: Theory and implementation
</B>
</A>

In this document we discuss the theoretical background and the
practical implementation of equations describing forced,
time-harmonic oscillations of elastic bodies.
</TD>
<TD>
- The time-harmonic equations of linear elasticity and their
  non-dimensionalisation.
- Implementation of the equations, based on a displacement formulation.
- The solution of test problem: The forced oscillations of an
  annular elastic region, computed with and without spatial adaptation.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../time_harmonic_linear_elasticity/unstructured_elastic_annulus/html/index.html">
<B>Solving the equations of time-harmonic linear
elasticity on unstructured meshes
</B>
</A>

We re-visit the solution of the equations of time-harmonic linear
elasticity -- this time using an unstructured mesh.
</TD>
<TD>
- Solving the equations of time-harmonic linear
  elasticity on unstructured meshes, with and without adaptation.
- How to assign different material properties to different parts of
  the domain.
.
</TD>

</TR>

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="gen_time_harmonic_lin_elas"><H2 class="header">
Generalised time-harmonic linear elasticity and perfectly matched 
layers (PMLs)</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../generalised_time_harmonic_linear_elasticity/html/index.html">
<B>The generalised equations of time-harmonic linear
elasticity and perfectly matched layers (PMLs).
</B>
</A>

In this document we discuss a generalisation of the 
equations of time-harmonic linear elasticity that allows the
implementation of far field boundary condition by perfectly
matched layers
</TD>
<TD>
- The generalised time-harmonic equations of linear elasticity and their
  non-dimensionalisation.
- Perfectly matched layers
- The solution of test problem: Oscillations of an infinite 2D medium
  forced by the prescribed oscillations of an embedded circular disk. 
.
</TD>

</TR>

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="fourier_decomp_lin_elas"><H2 class="header">
Azimuthally Fourier-decomposed 3D time-harmonic linear elasticity</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../time_harmonic_fourier_decomposed_linear_elasticity/cylinder/html/index.html">
<B>The azimuthally Fourier-decomposed equations of 3D time-harmonic linear
elasticity: Theory and implementation
</B>
</A>

In this document we discuss the theoretical background and the
practical implementation of equations describing  forced,
time-harmonic, non-axisymmetric oscillations of axisymmetric elastic
bodies.
</TD>
<TD>
- The time-harmonic equations of linear elasticity and their
non-dimensionalisation.
- Fourier decomposition of general three-dimensional disturbances.
- Implementation of the equations, based on a displacement formulation.
- The solution of a "manufactured" test problem.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../time_harmonic_fourier_decomposed_linear_elasticity/adaptive_pressure_loaded_cylinder/html/index.html">
<B>The spatially-adaptive solution of the azimuthally Fourier-decomposed equations of 3D time-harmonic linear
elasticity on unstructured meshes.
</B>
</A>

We simulate the forced, time-harmonic oscillations of a
hollow cylinder loaded by a spatially-constant pressure load
on its inner surface.
</TD>
<TD>
- The use of spatial adaptivity in the solution of the azimuthally 
  Fourier-decomposed equations of 3D time-harmonic linear 
  elasticity on unstructured meshes.
.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="beam"><H2 class="header">Beam structures</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../beam/tensioned_string/html/index.html">
<B>The deformation of a pre-stressed elastic beam, loaded by a pressure load
</B>
</A>

We study the lateral deformation of a pre-stressed elastic beam,
using \c oomph-lib's geometrically non-linear Kirchhoff-Love-type
\c HermiteBeamElement and compare the results against an
(approximate) analytical solution. 
</TD>
<TD>
- How to specify the undeformed reference shape for 
  the Kirchhoff-Love-type beam elements.
- How to apply boundary conditions and loads for the 
  \c HermiteBeamElement.
- General conventions:
  - How to change control parameters for the Newton solver.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../beam/steady_ring/html/index.html">
<B>Large-displacement post-buckling of a pressure-loaded,
thin-walled elastic ring
</B>
</A>

We compute the post-buckling deformation of a thin-walled elastic
ring, subjected to a pressure load and compare the
results against results from the literature.
</TD>
<TD>
- How to use \c oomph-lib's \c DisplacementControlElement to apply
  displacement control in solid mechanics problems.
- General conventions:
  - What should be stored in a \c GeneralisedElement's "external" \c Data?
  - What should be stored in a \c Problem's "global" \c Data?
.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD>
<A HREF="../../beam/unsteady_ring/html/index.html">
<B>Large-amplitude oscillations of a thin-walled elastic ring.
</B>
</A>

We compute the free, large-amplitude oscillations of a thin-walled elastic
ring and demonstrate that Newmark's method is energy conserving.
</TD>
<TD>
- How to assign initial conditions for beam structures.
- How to use the dump/restart function for \c HermiteBeamElements.
- Demonstrate that \c Newmark timesteppers can be used with
  variable timesteps. 
- How to retrieve solutions at previous timesteps in computations 
  with \c Newmark timesteppers.
- Changing the default non-dimensionalisation for time. 
- The non-dimensionalisation of the kinetic and potential (strain)
  energies of \c HermiteBeamElements.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../beam/lin_unsteady_ring/html/index.html">
<B>Small-amplitude oscillations of a thin-walled elastic ring.
</B>
</A>

We compute the free, small-amplitude oscillations of a thin-walled elastic
ring, demonstrate that Newmark's method is energy conserving,
and compare the oscillation frequencies and mode shapes against
analytical predictions.
</TD>
<TD>
- How to assign initial conditions for beam structures.
- General conventions:
  - ***
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="shell"><H2 class="header">Shell structures</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../shell/clamped_shell/html/index.html">
<B>Large-displacement post-buckling of a clamped, circular cylindrical
shell.
</B>
</A>

We simulate the post-buckling deformation of a pressure-loaded,
clamped, thin-walled elastic shell. 
</TD>
<TD>
- How to specify the undeformed reference configuration with
  \c HermiteShellElement structures.
- Using displacement control for  \c HermiteShellElements.
- General conventions:
  - ***
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="fsi"><H2 class="header">Large-displacement fluid-structure interaction problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../interaction/circle_as_element/html/index.html">
<B>
Warm-up problem for free-boundary problems: How to parametrise unknown
boundaries.
</B>
</A>

We demonstrate how to "upgrade" a \c GeomObject to a \c
GeneralisedElement so that it can be used to parametrise an unknown 
domain boundary.
</TD>
<TD>
- How to use multiple inheritance to combine \c GeomObjects and \c
  GeneralisedElements. 
- How to "upgrade" a \c GeomObject to a \c
  GeneralisedElement so that it can be used to parametrise an unknown 
  domain boundary.
- What is a \c GeomObject's geometric \c Data?
- What is a \c GeneralisedElement's external and internal \c Data?
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/macro_element_free_boundary_poisson/html/index.html">
<B>
A toy interaction problem: The solution of a 2D Poisson equation coupled to the deformation of the domain boundary
</B>
</A>

We show how to use \c MacroElementNodeUpdateElements and
\c MacroElementNodeUpdateMeshes
to implement sparse node update operations in free-boundary problems. 
We demonstrate their
use in a simple free-boundary problem: The solution of Poisson's
equation, coupled to an equation that determines the position of the
domain boundary.
</TD>
<TD>
- How to use \c MacroElementNodeUpdateElements and
  \c MacroElementNodeUpdateMeshes to implement
  sparse node update operations in free-boundary problems.
- Basic free-boundary problems: Making the domain boundary dependent on
  the solution in the domain.
- General conventions:
  - \c MacroElementNodeUpdateElements and \c MacroElementNodeUpdateMeshes.
</TD>

</TR>

<TR>
<TD>
<A HREF="../../interaction/fsi_collapsible_channel/html/index.html">
<B>A classical fluid-structure interaction problem: Finite Reynolds
number flow in a 2D channel with an elastic wall.
</B>
</A>

We demonstrate the solution of this classical fluid-structure interaction
problem and demonstrate how easy it is to combine the two
single-physics problems (the 
<A HREF="../../beam/tensioned_string/html/index.html">
deformation of an elastic beam under pressure loading</A> and 
the <A HREF="../../navier_stokes/collapsible_channel/html/index.html">
flow in a 2D channel with a moving wall</A>) to a fully-coupled 
fluid-structure interaction problem.
</TD>
<TD>
- The \c FSIFluidElements and \c FSIWallElement base classes.
- Representing a discretised beam/shell structure as a "compound"
  \c GeomObject: The \c MeshAsGeomObject class.
- Using the function 
  \c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
  to set up the fluid-structure interaction.
- The pros (convenient!) and cons (slow!)  of the 
  \c MacroElement/Domain - based node-update procedures
  in fluid-structure interaction problems.
.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/fsi_channel_with_leaflet/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/fsi_collapsible_channel_algebraic/html/index.html">
<B>Finite Reynolds number flow in a 2D channel with an elastic wall revisited:
Sparse algebraic node updates
</B>
</A>

We revisit the 
<A HREF="../../interaction/fsi_collapsible_channel/html/index.html">
problem of flow in a collapsible channel</A> to 
demonstrate that the sparse algebraic node update procedures
first discussed in an 
<A HREF="../../navier_stokes/algebraic_collapsible_channel/html/index.html">
earlier non-FSI example</A> lead to a much more efficient code.
</TD>
<TD>
- How to "sparsify" the node update with algebraic node update
  procedures.
- The \c GeomObject::locate_zeta(...) function and its default
  implementation.
.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/fsi_channel_with_leaflet/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/fsi_collapsible_channel_adapt/html/index.html">
<B>Finite Reynolds number flow in a 2D channel with an elastic wall
revisited again: Spatial adaptivity in fluid-structure interaction problems.
</B>
</A>

We revisit the 
<A HREF="../../interaction/fsi_collapsible_channel/html/index.html">
problem of flow in a collapsible channel</A> yet again to 
demonstrate the use of spatial adaptivity in fluid-structure
interaction problems.
</TD>
<TD>
- How to use spatial adaptivity in fluid-structure interaction problems.
- The \c Steady<NSTEPS> timestepper: How to assign positional
  history values for newly created nodes.
- Updating the node-update data in refineable \c AlgebraicMeshes.
.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/fsi_channel_with_leaflet/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 
</TD>

</TR>


</TABLE>
<TABLE class="example-table">
<TR>
<TD>
<A HREF="../../interaction/fsi_channel_segregated_solver/html/index.html">
<B>
Segregated solvers for fluid-structure-interaction problems: 
Revisiting the flow in a 2D collapsible channel
</B>
</A>

We revisit the 
<a HREF="../../interaction/fsi_collapsible_channel/html/index.html">
problem of flow in a collapsible channel</a> once more to 
demonstrate the use of segregated solvers in fluid-structure
interaction problems.
</TD>
<TD>
- The base SegregatableFSIProblem class
- How to construct and solve a segregated problem from a 
(standard) "monolithic" problem 
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../preconditioners/fsi/html/index.html">
<B>
Preconditioning monolithic solvers for fluid-structure-interaction problems: 
Revisiting the flow in a 2D collapsible channel yet again
</B>
</A>

We revisit the 
<a HREF="../../interaction/fsi_collapsible_channel/html/index.html">
problem of flow in a collapsible channel</a> yet again to 
demonstrate the use of \c oomph-lib's FSI preconditioner
for the monolithic solution of fluid-structure
interaction problems in which the node update in the fluid mesh
is performed by algebraic node updates.
</TD>
<TD>
- How to use \c oomph-lib's \c FSIPreconditioner
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/fsi_channel_with_leaflet/html/index.html">
<B>Flow past a flexible leaflet
</B>
</A>

We study the flow in a 2D channel that is partially obstructed 
by an elastic leaflet.
</TD>
<TD>
- FSI problems with fully immersed beam structures (i.e. beams
  that are subjected to the fluid traction on both faces).
- Another application of \c oomph-lib's \c FSIPreconditioner. 

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/fsi_channel_with_leaflet/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span> 
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/turek_flag/html/index.html">
<B>Turek & Hron's FSI benchmark: Flow past an elastic flag
attached to a cylinder
</B>
</A>

We demonstrate how to discretise and solve this benchmark problem with
\c oomph-lib.
</TD>
<TD>
- FSI problems with "proper" 2D solids (rather than beam or shell structures)
- FSI problems with wall inertia.
- Another application of \c oomph-lib's \c FSIPreconditioner. 


- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/turek_flag/html/index.html">
  separate tutorial</a> that discusses the parallelisation of this code
  and shows to modify it to allow spatial adaptivity in the fluid 
  and solid meshes.</span>
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/unstructured_fsi/html/index.html">
<B>Using unstructured meshes for FSI problems.
</B>
</A>

We demonstrate how to use xfig/triangle-generated unstructured
meshes (together with a pseudo-solid node update strategy for the fluid mesh) 
in fluid-structure interaction problems.
</TD>
<TD>
- How to use pseudo-elasticity to update fluid meshes in FSI problems.
- How to use xfig/triangle-generated unstructured
  meshes for fluid-structure interaction problems.
- The automatic generation of boundary coordinates for
  xfig/triangle-generated, unstructured meshes. How it works
  and what can go wrong...
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/unstructured_three_d_fsi/html/index.html">
<B>Using unstructured meshes for 3D FSI problems.
</B>
</A>

We demonstrate how to use tetgen-generated unstructured
meshes for 3D fluid-structure interaction problems.
</TD>
<TD>
- How to use tetgen-generated unstructured
  meshes for 3D fluid-structure interaction problems.
- The automatic generation of boundary coordinates for
  tetgen-generated, unstructured meshes. How it works
  and what can go wrong...
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD>
<A HREF="../../interaction/unstructured_adaptive_fsi/html/index.html">
<B>2D FSI on unstructured meshes  with adaptivity
</B>
</A>

We demonstrate how to use spatial adaptivity on unstructured
meshes for 2D fluid-structure interaction problems.
</TD>
<TD>
- How to use inline mesh generation to provide spatial adaptivity for
  fluid-structure interaction problems
</TD>

</TR>

<TR>
<TD>
<A HREF="../../interaction/vmtk_fsi/html/index.html">
<B>Finite-Reynolds-number flow through an elastic iliac bifurcation
</B>
</A>

We show how to simulate physiological fluid-structure interaction 
problems, using the  <a href="http://www.vmtk.org">Vascular 
Modeling Toolkit (VMTK).</a>
</TD>
<TD>
- How to attach multiple \c FaceElements to the same node --
  distinguishing different Lagrange multipliers.
- How to solve physiological fluid-structure interaction problems.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../preconditioners/pseudo_solid_fsi/html/index.html"><B>
<code>oomph-lib</code>'s 
preconditioner for the solution of fluid-structure interaction problems with
pseudo-solid node updates for the fluid mesh</B></A> 

We discuss <code>oomph-lib</code>'s preconditioner for the
solution of FSI problems in which the fluid node update is performed
by pseudo-elasticity.
</TD>
<TD>
- How to use <code>oomph-lib</code>'s  for the solution of FSI 
  problems in which the fluid node update is performed
  by pseudo-elasticity.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../interaction/osc_ring_macro/html/index.html">
<B>A simple fluid-structure interaction problem: Finite Reynolds
number flow, driven by an oscillating ring.
</B>
</A>

This is a very simple fluid-structure interaction problem: We study
the finite-Reynolds number internal flow generated by an 
oscillating ring. The wall motion only has a single degree of freedom:
The ring's average radius, which needs to be adjusted to 
conserve mass. The nodal positions in the fluid domain is updated by \c
MacroElements. [This is a warm-up problem for the full
fluid structure interaction problem discussed in the next example].
We compare the predictions for the flow field against
asymptotic results. 
</TD>
<TD>
- Fluid-structure interaction.
- General conventions:
  - ***
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../interaction/osc_ring_algebraic/html/index.html">
<B>A simple fluid-structure interaction problem re-visited: Finite Reynolds
number flow, driven by an oscillating ring -- this time with
algebraic updates for the nodal positions.
</B>
</A>

We re-visit the simple fluid-structure interaction problem
considered in the <A HREF="../../interaction/osc_ring_macro/html/index.html">
earlier example</A>.This time we perform the update of the
nodal positions with \c AlgebraicElements.
</TD>
<TD>
- Fluid-structure interaction.
- General conventions:
  - ***
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../interaction/fsi_osc_ring/html/index.html">
<B>A real fluid-structure interaction problem: Finite Reynolds
number flow in an oscillating elastic ring.
</B>
</A>

Our first "real" fluid-structure interaction problem: We study
the finite-Reynolds number internal flow generated by the motion 
of an oscillating elastic ring and compare the results against
asymptotic predictions. 
</TD>
<TD>
- Fluid-structure interaction.
- General conventions:
  - ***
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="acoustic_fsi"><H2 class="header">Acoustic fluid-structure interaction problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../acoustic_fsi/acoustic_fsi_annulus/html/index.html">
<B>
A time-harmonic acoustic fluid-structure interaction problem: Sound radiation
from an immersed oscillating cylinder that is coated with an elastic
layer.
</B>
</A>

We provide an overview of \c oomph-lib's methodology for the
solution of time-harmonic acoustic fluid-structure interaction problems which
are based on a coupled solution of the 
<a href="../../time_harmonic_linear_elasticity/elastic_annulus/html/index.html">time-harmonic equations of linear elasticity</a>
and
<a href="../../helmholtz/scattering/html/index.html">the Helmholtz
equation.</a>
</TD>
<TD>
- Theory: The formulation and non-dimensionalisation 
  of time-harmonic acoustic fluid-structure 
  interaction problems. 
- Coupling 
  <a href="../../time_harmonic_linear_elasticity/elastic_annulus/html/index.html">time-harmonic equations of linear elasticity</a>
  and <a href="../../helmholtz/scattering/html/index.html">the Helmholtz
  equation,</a> using the 
  \c
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh(...)
  helper function
- Validation of the methodology via a comparison against an analytical solution.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../acoustic_fsi/unstructured_acoustic_fsi_annulus/html/index.html">
<B>
A time-harmonic acoustic fluid-structure interaction problem: Sound radiation
from an immersed oscillating cylinder that is coated with an elastic
layer -- this time solved on an unstructured mesh.
</B>
</A>

A brief extension of the 
<A HREF="../../acoustic_fsi/acoustic_fsi_annulus/html/index.html">
previous tutorial</a> illustrating how to solve the problem 
with unstructured meshes.
</TD>
<TD>
- The solution of time-harmonic acoustic fluid-structure interaction
  problems on unstructured meshes.
- How to assign different material properties to different regions
  of the elastic body.
</TD>

</TR>

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="fd_acoustic_fsi"><H2 class="header">Azimuthally Fourier-decomposed 
3D acoustic fluid-structure interaction problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../fourier_decomposed_acoustic_fsi/sphere/html/index.html">
<B>
An azimuthally Fourier-decomposed time-harmonic acoustic 
fluid-structure interaction problem: Sound radiation
from an immersed oscillating sphere that is coated with an elastic
layer.
</B>
</A>

We provide an overview of \c oomph-lib's methodology for the
solution of azimuthally Fourier-decomposed 
time-harmonic acoustic fluid-structure interaction problems which
are based on a coupled solution of the 
<a href="../../time_harmonic_fourier_decomposed_linear_elasticity/cylinder/html/index.html">time-harmonic equations of linear elasticity</a>
and
<a href="../../fourier_decomposed_helmholtz/sphere_scattering/html/index.html">the Helmholtz
equation.</a>
</TD>
<TD>
- Theory: The formulation and non-dimensionalisation 
  of time-harmonic acoustic fluid-structure 
  interaction problems in cylindrical polar coordinates, using a 
  Fourier-decomposition of all fields in the azimuthal direction.
- Coupling 
  <a href="../../time_harmonic_fourier_decomposed_linear_elasticity/cylinder/html/index.html">time-harmonic equations of linear elasticity</a>
  and <a href="../../fourier_decomposed_helmholtz/sphere_scattering/html/index.html">the Helmholtz
  equation,</a> using the 
  \c
  Multi_domain_functions::setup_bulk_elements_adjacent_to_face_mesh(...)
  helper function
- Validation of the methodology via a comparison against an analytical solution.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../fourier_decomposed_acoustic_fsi/unstructured_sphere/html/index.html">
<B>
A time-harmonic acoustic fluid-structure interaction problem: Sound radiation
from an immersed oscillating sphere that is coated with an elastic
layer -- this time solved on an unstructured mesh with spatial adaptivity.
</B>
</A>

A brief extension of the 
<A HREF="../../fourier_decomposed_acoustic_fsi/sphere/html/index.html">
previous tutorial</a> illustrating how to solve the problem 
on adaptive, unstructured meshes.
</TD>
<TD>
- The spatially-adaptive solution of Fourier-decomposed 
  time-harmonic acoustic fluid-structure interaction
  problems on unstructured meshes.
- How to assign different material properties to different regions
  of the elastic body.
</TD>

</TR>


</TABLE>
<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="multi"><H2 class="header">Multi-physics problems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A NAME="bous">
<A HREF="../../multi_physics/b_convection/html/index.html">
<B>
Simple multi-physics problem: How to combine existing single-physics 
elements into new multi-physics elements.</B>
</A></A>

We demonstrate how to "combine" a \c CrouzeixRaviartElements and 
\c QAdvectionDiffusionElements into a single \c BuoyantCrouzeixRaviartElement
that solves the Navier--Stokes equations under the Boussinesq approximation
coupled to an energy equation.
</TD>
<TD>
- How to use multiple inheritance to combine two single-physics elements.
- How to write single-physics elements that can be combined into multi-physics
  elements.
- How to use the \c Problem::steady_newton_solve(...) function to find
  steady solutions of unsteady problems.


- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/boussinesq_convection/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span>
</TD>

</TR>
<TR>
<TD>
<A HREF="../../multi_physics/refine_b_convect/html/index.html">
<B>
Refineable multi-physics problem: How to combine existing refineable 
single-physics elements into new refineable multi-physics elements.</B>
</A>

We demonstrate how to "combine" a \c RefineableCrouzeixRaviartElements and 
\c RefineableQAdvectionDiffusionElements 
into a single \c RefineableBuoyantCrouzeixRaviartElement
that solves the Navier--Stokes equations under the Boussinesq approximation
coupled to an energy equation.
</TD>
<TD>
- How to use multiple inheritance to combine two refineable
  single-physics elements.
- How to choose the "Z2 flux" for multi-physics elements that are both
  derived from the \c ElementWithZ2ErrorEstimator class.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/boussinesq_convection/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span>
</TD>

</TR>
<TR>
<TD>
<A HREF="../../multi_physics/multi_domain_ref_b_convect/html/index.html">
<B>
Solving multi-field problems with multi-domain discretisations.</B>
</A>

We demonstrate an alternative approach to the solution of multi-field
problems, in which the governing PDEs are discretised in separate
meshes and interact via "external elements".
</TD>
<TD>
- How to discretise multi-field problems with multi-domain approaches.
- The \c Multi_domain_functions namespace.

- <span style="color:red"><b>Note:</b> There is a 
  <A HREF="../../mpi/boussinesq_convection/html/index.html">
  separate tutorial</a>
  that discusses the parallelisation of this code.</span>
</TD>

</TR>
<TR>
<TD>
<A NAME="thermo">
<A HREF="../../multi_physics/thermo/html/index.html">
<B>
Thermoelasticity: How to combine single-physics elements with solid  
mechanics elements.</B>
</A></A>

We demonstrate how to "combine" a \c QUnsteadyHeatElement and 
\c QPVDElement into a single \c QThermalPVDElement
that solves the equations governing elastic deformations coupled to
 uniform thermal expansion. The geometric coupling back to the heat
equation is completely hidden.
</TD>
<TD>
- How to use multiple inheritance to combine a single-physics element
and a solid element
- <span style="color:red">Note: this driver code is currently undocumented.</span>
</TD>
</TR>
<TR>
<TD>
<A NAME="surfactant">
<A HREF="../../multi_physics/rayleigh_instability_surfactant/html/index.html">
<B>
Surfactant Transport: How to add surface transport equations to 
free surface elements.</B>
</A></A>

We demonstrate how to add general surface transport equations to the
\c FluidInterfaceElements and apply them to determine the effects of 
insoluble surfactant
on the Rayleigh--Plateau instability.
</TD>
<TD>
- How to use the member function \c 
add_additional_residual_contributions_interface()
to include additional terms in free surface problems.
</TD>
</TR>


</TABLE>
<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="eigen"><H2 class="header">Eigenproblems</H2></A>
</TD>
</TR>
<TR>
<TD>
<A NAME="harmonic">
<A HREF="../../eigenproblems/harmonic/html/index.html">
<B>
How to formulate and solve an eigenproblem.</B>
</A></A>

We demonstrate how to write \c QHarmonicElements that solve the
  eigenvalues and eigenfunctions of the one-dimensional Laplace operator.
</TD>
<TD>
- A brief introduction to linear stability theory.
- How to use \c oomph-lib's interfaces to eigensolvers.
- How to write elements that can be used in eigenproblems.
- How to use the \c Problem::solve_eigenproblem(...) function to find
  solutions of the eigenproblems.
- Post-processing of eigenvalue problems.
</TD>

</TR>

<TR>
<TD>
<A NAME="complex harmonic">
<A HREF="../../eigenproblems/complex_harmonic/html/index.html">
<B>
How to formulate and solve an eigenproblem involving complex eigenvalues.</B>
</A></A>

We demonstrate how to write \c QComplexHarmonicElements that determine the
  eigenvalues and eigenfunctions of the shifted one-dimensional 
  Laplace operator as a quadratic eigenvalue problem.
</TD>
<TD>
- How to write eigenproblems that involve more than one field.
- How to write a quadratic eigenproblem as two linear eigenproblems.
- How to output complex eigenvalues and eigenvectors.
</TD>

</TR>


</TABLE>  

<A NAME="meshes"><H1 class="section-header">Mesh generation</H1></A>


<TABLE class="example-table">

<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="available_meshes"><H2 class="header">Structured meshes</H2></A>
</TD>
</TR>
<TR>
<TD >
<B>
<A HREF="../../meshes/mesh_list/html/index.html">
Structured meshes
</A>
</B>

We list \c oomph-lib's existing structured
meshes and provide a quick overview of their common features.
</TD>
<TD>
- Reminder of \c oomph-lib's design features that facilitate
  the re-use of meshes.
- What structured meshes are available?
.        
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="third_party_meshes"><H2 class="header">Unstructured meshes using
input from third-party mesh generators</H2></A>
</TD>
</TR>
<TR>
<TD>
<B>
<A NAME="third_party_meshes">
<A HREF="../../meshes/third_party_meshes/html/index.html">
Unstructured meshes generated via input from
           third-party mesh generators
</A>
</A>
</B>

We describe \c oomph-lib's wrappers to third-party 
(unstructured) mesh generators. 
</TD>
<TD>
- <A HREF="../../meshes/mesh_from_triangle/html/index.html"><B>
  TriangleMesh<ELEMENT>:</B></A>  A Mesh based on the output from
  <A HREF="http://www.cs.cmu.edu/~jrs">J.R.Shewchuk's</A>
  Delaunay mesh generator  
  <A HREF="http://www.cs.cmu.edu/~quake/triangle.html">
  <TT>Triangle</TT></A>
- <A HREF="../../meshes/mesh_from_tetgen/html/index.html"><B>
  TetgenMesh<ELEMENT>:</B></A> A Mesh based on the output from
  <A HREF="http://www.wias-berlin.de/~si">Hang Si's</A> 
  unstructured tetrahedral mesh generator
  <A HREF="http://wias-berlin.de/software/tetgen//index.html"><TT>TetGen.</TT></A>
- We provide a <a href="../../meshes/mesh_from_vmtk/html/index.html">separate 
  tutorial</a> that shows how to generate \c oomph-lib meshes
  from medical images, using the 
  <a href="http://www.vmtk.org">Vascular Modeling
  Toolkit (VMTK).</a>
- <A HREF="../../meshes/mesh_from_geompack/html/index.html"><B>
  GeompackQuadMesh<ELEMENT>:</B></A> A Mesh based on the output from 
  Barry Joe's mesh generator 
  <A HREF="http://members.shaw.ca/bjoe/">\c Geompack++, </A>
  available as freeware at 
  <A HREF="http://members.shaw.ca/bjoe/">
  http://members.shaw.ca/bjoe/.</A>
</TD>

</TR>
<TR>
<TD>
<B>
<A HREF="../../meshes/mesh_from_inline_triangle/html/index.html">
Inline unstructured mesh generation </A>
</B>

We describe how to generate unstructured 2D meshes from 
 within an \c oomph-lib driver code.
</TD>
<TD>
- How to generated inline unstructured 2D meshes using <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle.</a>
- How to create unstructured meshes with curvilinear boundaries
- How to perform spatial adaptivity on unstructured meshes via complete remeshing
</TD>

</TR>

<TR>
<TD>
<B>
<A HREF="../../meshes/mesh_from_inline_triangle_internal_boundaries/html/index.html">
Inline unstructured mesh generation including internal boundaries
</A>
</B>

We describe how to generate unstructured 2D meshes that contain internal
boundaries, delineating different regions of space, from within 
an \c oomph-lib driver code.
</TD>
<TD>
- How to generated more complicated inline unstructured
2D meshes using <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle.</a>
- How to create additional regions in space
</TD>

</TR>


<TR>
<TD>
<B>
<A NAME="xfig_mesh">
<A HREF="../../meshes/mesh_from_xfig/html/index.html">
Mesh generation with <code>xfig</code>
</a>
</B>

\c oomph-lib's one-and-only GUI: Generating unstructured
triangular meshes using <a href="http://en.wikipedia.org/wiki/Xfig">xfig</a>
and <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle</a>
</TD>
<TD>
- How to generated unstructured
triangular meshes using <a href="http://en.wikipedia.org/wiki/Xfig">xfig</a>
and <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle.</a>
</TD>

</TR>
</TABLE>
  



<A NAME="solvers"><H1 class="section-header">Linear solvers and preconditioners</H1></A>


<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="linear_solvers"><H2 class="header">Direct and iterative linear solvers 
and general-purpose preconditioners</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../linear_solvers/html/index.html"><B>
Overview</B></A> 

We provide an overview of <code>oomph-lib's</code> direct and iterative
linear solvers and preconditioners.
</TD>
<TD>
- How to change the linear solver for \c oomph-lib's Newton
  solver.
- How to use \c oomph-lib's \c IterativeLinearSolvers and
  \c Preconditioners.
- How to use \c oomph-lib's wrappers to the third-party 
  iterative linear solvers/preconditioners from the
  \c Hypre and \c Trilinos libraries.
.
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="lin_alg"><H2 class="header">(Distributed) linear algebra and
oomph-lib's block preconditioning framework</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../mpi/distributed_linear_algebra_infrastructure/html/index.html"><B>
(Distributed) Linear Algebra Infrastructure</B></A>
We provide an overview of <code>oomph-lib's</code> (distributed)
linear algebra infrastructure.
</TD>
<TD>
- The \c OomphCommunicator.
- The \c LinearAlgebraDistribution and the \c
  DistributedLinearAlgebraObject base class.
- The \c CRDoubleMatrix and the \c DoubleVector.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/block_preconditioners/html/index.html"><B>
(Distributed) Block preconditioners</B></A> 
We provide an overview of <code>oomph-lib's</code> (distributed)
block preconditioning framework and demonstrate how to write
a new block preconditioner.
</TD>
<TD>
- Theory
- Block preconditionable elements: Classifying the block and dof
  types.
- Master and subsidiary preconditioners.
- An example: A simple implementation of an FSI preconditioner.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/distributed_general_purpose_block_preconditioners/html/index.html"><B>
(Distributed) General-purpose block preconditioners</B></A>
We provide an overview of <code>oomph-lib's</code> (distributed)
general purpose block preconditioners
</TD>
<TD>
- Theory
- Block diagonal and block triangular preconditioners.
- Two-level parallelisation for block diagonal preconditioning.
- How to specify subsidiary preconditioners for the (approximate) 
  solution of the linear systems involving the diagonal blocks.
- How to specify the dof types.
.
</TD>

</TR>
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="specific_preconditioners"><H2 class="header">Problem-specific 
preconditioners</H2></A>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../preconditioners/lsc_navier_stokes/html/index.html"><B>
<code>oomph-lib</code>'s Least-Squares-Commutator (LSC) 
Navier-Stokes preconditioner</B></A>

We discuss <code>oomph-lib</code>'s implementation of 
Elman, Silvester & Wathen's Least-Squares-Commutator (LSC) 
Navier-Stokes preconditioner.
</TD>
<TD>
- How to use <code>oomph-lib</code>'s Least-Squares-Commutator (LSC) 
Navier-Stokes preconditioner.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../preconditioners/fsi/html/index.html"><B>
<code>oomph-lib</code>'s fluid-structure interaction 
preconditioner</B></A>

We discuss <code>oomph-lib</code>'s preconditioner for the
solution of monolithically-discretised fluid-structure interaction 
problems with algebraic node updates.
</TD>
<TD>
- How to use <code>oomph-lib</code>'s FSI preconditioner for problems
  with algebraic node updates.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../preconditioners/prescribed_displ_lagr_mult/html/index.html"><B>
<code>oomph-lib</code>'s 
preconditioner for the solution of solid mechanics problems with
prescribed boundary displacements</B></A> 

We discuss <code>oomph-lib</code>'s preconditioner for the
solution of solid mechanics problems in which the displacement
of a boundary is prescribed and imposed by Lagrange multipliers.
This preconditioner is an important building block for the
solution of FSI problems in which the fluid node update is performed
by pseudo-elasticity.
</TD>
<TD>
- How to use <code>oomph-lib</code>'s  for the solution of solid 
  mechanics problems with prescribed boundary displacements.
- How to specify subsidiary preconditioners (inexact solvers) 
  with any of <code>oomph-lib</code>'s existing block preconditioners.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../preconditioners/pseudo_solid_fsi/html/index.html"><B>
<code>oomph-lib</code>'s 
preconditioner for the solution of fluid-structure interaction problems with
pseudo-solid node updates for the fluid mesh</B></A> 

We discuss <code>oomph-lib</code>'s preconditioner for the
solution of FSI problems in which the fluid node update is performed
by pseudo-elasticity.
</TD>
<TD>
- How to use <code>oomph-lib</code>'s  for the solution of FSI 
  problems in which the fluid node update is performed
  by pseudo-elasticity.
.
</TD>

</TR>
</TABLE>
  






<A NAME="visualisation"><H1 class="section-header">Visualisation of the results</H1></A>

<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<A NAME="paraview"><H2 class="header">Paraview</H2></A> 
</TD>
</TR>
<TR>
<TD>

<A HREF="../../paraview/html/index.html"><B>
Displaying results with paraview</B></A>

We demonstrate how to use 
Angelo Simone's conversion scripts that allow the \c oomph-lib
results to be displayed by 
<a href="http://www.paraview.org">paraview.</a>
</TD>
<TD>
- How to display \c oomph-lib's results with
  <a href="http://www.paraview.org">paraview.</a>
.
</TD>

</TR>
</TABLE>
  



<A NAME="parallel"><H1 class="section-header">Parallel driver codes</H1></A>

Please consult the 
<A HREF="../../mpi/general_mpi/html/index.html">
<b>general tutorial on <code>oomph-lib</code>'s parallel processing
capabilities.</b></a>


<TABLE class="example-table">
<TR>
<TD COLSPAN=2 class="header-cell">
<H2 class="header">Distributed problems</H2>
</TD>
</TR>
<TR>
<TD>
<B>Example code</B>
</TD>
<TD>
<B>\c oomph-lib features/conventions illustrated 
by the example code</B>
</TD>

</TR>
<TR>
<TD>
<A NAME="distributed">
<A HREF="../../mpi/adaptive_driven_cavity/html/index.html"><B>
Parallel solution of the adaptive driven cavity problem</B></A></a> 

We demonstrate how to distribute a straightforward single-physics
problem.
</TD>
<TD>
- Initialising and finalising MPI
- Distributing the problem
- Specifying a pre-determined partition of the problem
- Modifying the output filename
- Pinning values in specific elements
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/adaptive_driven_cavity_load_balance/html/index.html"><B>
Parallel solution of the adaptive driven cavity problem with load
balancing</B></A></a>

We demonstrate the modifications required to perform a load balancing
step within the distributed adaptive driven cavity problem
</TD>
<TD>
- Load balancing a problem
- The build_mesh function
- The use of a default partition when performing a test run
- The actions before and after load balance functions
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/two_d_poisson_flux_bc_adapt/html/index.html"><B>
Parallel solution of the 2D Poisson problem with flux boundary
conditions</B></A></a>

We demonstrate the modifications required to distribute a problem
involving \c FaceElements.
</TD>
<TD>
- The actions before and after distribute functions and their use to 
  strip off and re-attach \c FaceElements that are used to enforce
  Neumann/flux boundary conditions.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/boussinesq_convection/html/index.html"><B>
Parallel solution of the Boussinesq convection problem</B></A></a> 

We demonstrate how to distribute a straightforward multi-physics
problem where two domains interact.
</TD>
<TD>
- More details on the \c Multi_domain_functions helper functions and 
  their parallel implementation. 
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/fsi_channel_with_leaflet/html/index.html"><B>
Parallel solution of an FSI problem: Channel with an elastic
leaflet</B></A></a> 

We demonstrate how to distribute FSI problems that use algebraic
update methods. 
</TD>
<TD>
- How to distribute fluid-structure interaction problems in which 
  algebraic node update methods are used to deform the fluid mesh 
  in response to changes in the shape of the domain boundary
- How to retain all elements in a mesh on all processors when the
  problem is distributed.
.
</TD>

</TR>
<TR>
<TD>
<A HREF="../../mpi/turek_flag/html/index.html"><B>
Parallel solution of Turek and Hron's FSI benchmark
problem</B></A></a> 

We demonstrate how to distribute a problem involving refineable 2D
solid and fluid meshes that interact along interface boundaries.
</TD>
<TD>
- How to "upgrade" 
  <A HREF="../../interaction/turek_flag/html/index.html">
  Turek & Hron's FSI benchmark problem</a> to allow
  spatial adaptivity in the fluid and solid meshes.
- How to retain selected elements on all processors when the
  problem is distributed.
- How to implement load-balancing of a distributed multi-domain
  problem, with particular reference to the actions before and after
  load balance functions.
.
</TD>

</TR>
</TABLE>
  



**/


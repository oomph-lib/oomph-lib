/**

\mainpage Example problem: The 2D Driven Cavity Problem

This is our first Navier-Stokes example problem. We discuss the
non-dimensionalisation of the equations and their implementation in \c
oomph-lib, and demonstrate the solution of the 2D driven cavity
problem. 


<HR>
<HR>

\section equation The Navier-Stokes equations
In dimensional form the 2D [3D] Navier-Stokes equations (in cartesian
coordinates \f$ x_i^*; \ i=1,2 [,3] \f$) are given by the momentum equations
<CENTER>
\f[
\rho\left(\frac{\partial u_i^*}{\partial t^*} + 
u_j^*\frac{\partial u_i^*}{\partial x_j^*}\right) 
=
- \frac{\partial p^*}{\partial x_i^*} + 
B_i^*(x_j^*,t^*) + \rho \, G^*_i +
\frac{\partial }{\partial x_j^*} \left[ \mu\left(
\frac{\partial u_i^*}{\partial x_j^*} +  
\frac{\partial u_j^*}{\partial x_i^*} \right) \right],
\f]
</CENTER>
and the continuity equation
<CENTER>
\f[
\frac{\partial u_i^*}{\partial x_i^*} = Q^*,
\f]
</CENTER>
where we have used index notation and the summation convention. 

Here, the velocity components are denoted by \f$ u_i^* \f$,
the pressure by \f$ p^* \f$, and time by \f$ t^* \f$, and we
have split the body force into two components: A constant
vector \f$ \rho \, G_i^* \f$ which typically represents gravitational
forces; and a variable body force, \f$ B_i^*(x^*_j,t^*) \f$.
\f$ Q^*(x^*_j,t^*) \f$ is a volumetric source term for the continuity equation 
and is typically equal to zero.


We non-dimensionalise the equations, using problem-specific reference
quantities for the velocity, \f$ {\cal U},\f$ length, \f$
{\cal L}, \f$ and time, \f${\cal T},\f$ and 
scale the constant body force vector on the gravitational
acceleration, \f$g, \f$ so that
<CENTER>
\f[
u_i^* = {\cal U} \, u_i, \qquad        
x_i^* = {\cal L}\, x_i, \qquad        
t^* = {\cal T}\, t,     \qquad
G_i^* = g \, G_i,
\f]
\f[
p^* = \frac{\mu_{ref} {\cal U}}{{\cal L}} \, p, \qquad
B_i^* = \frac{{\cal U}\mu_{ref}}{{\cal L}^2} \, B_i, \qquad 
 Q^* = \frac{{\cal U}}{{\cal L}}\, Q, \qquad
\f]
</CENTER>
where we note that the pressure and the variable body force have 
been non-dimensionalised on the viscous scale.
\f$ \mu_{ref} \f$ and \f$ \rho_{ref} \f$ (used below) are reference values 
for the fluid viscosity and density, respectively. In single-fluid
problems, they are identical to the viscosity \f$ \mu \f$ and 
density \f$ \rho \f$ of the (one and only) fluid in the problem.


The non-dimensional form of the Navier-Stokes equations
is then given by
<CENTER>
\f[
R_\rho Re \left( St \frac{\partial u_i}{\partial t} + 
u_j \frac{\partial u_i}{\partial x_j} \right) =
- \frac{\partial p}{\partial x_i}  + B_i(x_j,t) +
R_\rho \frac{Re}{Fr} G_i +
\frac{\partial }{\partial x_j} \left[  R_\mu \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right) \right],
\f]
</CENTER>
and
<CENTER>
\f[
\frac{\partial u_i}{\partial x_i} = Q,
\f]
</CENTER>
where the dimensionless parameters
<CENTER>
\f[
Re = \frac{{\cal U}{\cal L}\rho_{ref}}{\mu_{ref}}, \qquad 
St = \frac{{\cal L}}{{\cal U}{\cal T}}, \qquad 
Fr = \frac{{\cal U}^2}{g{\cal L}},
\f]
</CENTER>
are the Reynolds number, Strouhal number and Froude number
respectively. \f$ R_\rho=\rho/\rho_{ref} \f$ and 
\f$ R_\mu =\mu/\mu_{ref}\f$ represent the ratios
of the fluid's density and its dynamic viscosity, relative to the
density and viscosity values used to form the non-dimensional
parameters (By default, \f$ R_\rho  = R_\mu = 1 \f$; other values
tend to be used in problems involving multiple fluids). 

The above equations are typically augmented by Dirichlet boundary conditions
for (some of) the velocity components. On boundaries where 
no velocity boundary conditions are applied, the flow satisfies the
"traction free" natural boundary condition
\f$ t_i = -p \, n_i + R_\mu (\partial u_i/\partial x_j + 
\partial u_j/\partial x_i ) \, n_j = 0. \f$ 
(We refer to <A HREF="../../rayleigh_traction_channel/html/index.html"> another
example </A> for an illustration of how to apply traction  boundary
conditions for  the Navier-Stokes equations.)

If the velocity is prescribed along the entire domain boundary, 
the fluid pressure \f$ p \f$ is only determined up to an arbitrary
constant. This indeterminacy may be overcome by prescribing the 
value of the pressure at a single point in the domain (see 
<A HREF="#exercises"> exercise</A>).
<HR>
<HR>


\section element_types Implementation 

<HR>

\subsection el The elements
\c oomph-lib provides two LBB-stable isoparametric Navier-Stokes 
elements that are based on the \c QElement<DIM,3> family of 
geometric finite elements. They are nine-node quadrilateral
(for \c DIM=2), and 27-node brick (for \c DIM=3) elements  in which
the mapping between local and global (Eulerian) coordinates
is given by
<CENTER>
\f[
x_i = \sum_{j=1}^{N^{(E)}} X^{(E)}_{ij} \, \psi_j, \qquad
i=1,2\quad [\mbox{and }3].
\f]
</CENTER>
Here \f$ N^{(E)} = 3^{DIM}\f$ is the number of nodes in the element, \f$ X^{(E)}_{ij} \f$ is the \f$ i \f$-th global (Eulerian) coordinate
of the \f$ j \f$-th \c Node in the element, and the \f$ \psi_j \f$ 
are the element's geometric shape functions, defined in the \c QElement<DIM,3>
class.

In both elements the velocity components \f$ u_1, \f$ \f$ u_2, \f$ 
[and \f$ u_3 \f$] are stored as nodal values and the geometric
shape functions are used to interpolate the velocities inside the
element, 
<CENTER>
\f[
u_i = \sum_{j=1}^{N^{(E)}} U^{(E)}_{ij} \, \psi_j, \qquad 
i=1,2\quad [\mbox{and }3],
\f]
</CENTER>
where \f$  U^{(E)}_{ij} \f$  is the \f$ i \f$-th velocity component 
at \f$ j \f$-th \c Node in the element. Nodal values of the velocity 
components are accessible via the access function

\code
NavierStokesEquations<DIM>::u(i,j)
\endcode

which returns the \c i-th velocity component stored at the element's 
\c j-th \c Node.

The two elements differ in the way in which the pressure is
represented:

<HR>

\subsubsection CR Crouzeix-Raviart elements
In \c oomph-lib's \c QCrouzeixRaviartElements the pressure
is represented by,
<CENTER>
\f[
p = P_0^{(E)} + P_1^{(E)} \, s_1 + P_2^{(E)} \, s_2 
\quad \left[+ P_3^{(E)}  \, s_3\right],
\f]
</CENTER>
where the \f$ s_i \in [-1,1]\f$ are the element's local coordinates.
This provides a discontinuous, piecewise bi-[tri-]linear
representation of the pressure in terms of 3 [4] pressure degrees
of freedom per element.  Crouzeix-Raviart elements ensure that the 
continuity equation is satisfied within each element.

The pressure degrees of freedom are local to the element and are 
stored in the element's internal \c Data. They are accessible
via the member function

\code
QCrouzeixRaviartElement<DIM>::p(j)
\endcode

which returns the value of the \c j-th pressure degree of freedom 
in this element.

Each \c Node in a 2D [3D] Crouzeix-Raviart element stores 2 [3] nodal
values, representing the two [three] velocity components at that \c Node.

<HR>

\subsubsection TH Taylor-Hood elements
In \c oomph-lib's \c QTaylorHoodElements the pressure
is represented by a globally-continuous, piecewise bi-[tri-]linear
interpolation between the pressure values \f$ P_j^{(E)} \f$ 
that are stored at the elements' \f$ N_p^{(E)} = 2^{DIM}\f$
corner/vertex nodes,
<CENTER>
\f[
p = \sum_{j=1}^{N_p^{(E)}} P_j^{(E)} \, \psi_j^{(p)},
\f]
</CENTER>
where the \f$\psi_j^{(p)}\f$ are the bi-[tri-]linear pressure shape functions.

The first 2 [3] values of each \c Node in a 2D [3D] Taylor-Hood element store
the two [three] velocity components at that \c Node. The corner
[vertex] nodes store an additional value which represents the
pressure at that \c Node. The access function 
\code
QTaylorHoodElement<DIM>::p(j)
\endcode
returns the nodal pressure value at the element's \c j-th corner
[vertex] \c Node.


In sufficiently fine meshes, Taylor-Hood elements generate a 
much smaller number of pressure degrees of freedom than the 
corresponding Crouzeix-Raviart elements. However, Taylor-Hood elements 
do not conserve mass locally. 

<HR>

\section params Non-dimensional parameters and their default values

The Reynolds number, Strouhal number, inverse-Froude number, density
ratio and viscosity ratio are assumed to be constant within
each element. Their values are accessed via pointers which are
accessible via the member functions \c
re_pt() for \f$ Re, \f$ \c re_st_pt() for \f$ Re \ St, \f$ \c
re_invfr_pt() for \f$ Re/Fr, \f$ \c density_ratio_pt() for \f$ R_\rho \f$
and \c viscosity_ratio_pt() for \f$ R_\mu \f$.


By default the pointers point to default values (implemented as
static member data in the \c NavierStokesEquations class), therefore
they only need to be over-written if the default values
are not appropriate. The default values are:
- Default Reynolds number:
\f[
Re = 0 \mbox{\ \ \ \ (Stokes flow)}
\f] 
- Default Womersley number (product of Reynolds and Strouhal number):
\f[
Re St = 0 \mbox{\ \ \ \ (steady flow)}
\f]
- Default product of Reynolds and inverse Froude number (a measure of
gravity on the viscous scale)
\f[
Re/Fr = 0 \mbox{\ \ \ \ (no gravity)}
\f]
- Default viscosity ratio:
\f[
R_\mu = 1 \mbox{\ \ \ \ (viscosity = viscosity used in the definition
of $Re$)}
\f]
- Default density ratio:
\f[
R_\rho = 1 \mbox{\ \ \ \ (density = density used in the definition
of $Re$)}
\f]
.

We use the same approach for the specification of the body force
vector \f$ G_i \f$, accessible via the function \c g_pt (),
the variable body force \f$ B_i \f$, accessible via the 
function pointer \c body_force_fct_pt(), and the volumetric
source function \f$ Q \f$, accessible via the function pointer \c
source_fct_pt(time,x). By default the (function) pointers are set
such that 
- Default gravity vector:
\f[
G_i = 0 \mbox{\ \ \ \ (no gravity)}
\f] 
- Default body force function:
\f[
B_i(x_j,t) = 0 \mbox{\ \ \ \ (no body force)}
\f] 
- Default volumetric source function
\f[
Q(x_j,t) = 0 \mbox{\ \ \ \ (flow is divergence free)}
\f] 
.
<HR>
<HR>

\section example The example problem
We will illustrate the solution of the steady 2D Navier-Stokes
equations using the well-known example of the driven cavity.

<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The 2D steady driven-cavity problem in a square domain.</B>
</CENTER> 
Solve
\f[
Re\phantom{i}u_j\frac{\partial u_i}{\partial x_j} =
- \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
in the square domain \f$ D = \{x_i \in [0,1]; i=1,2 \} \f$,
subject to the Dirichlet boundary conditions:
\f[
\left. \mathbf{u}\right|_{\partial D}=(0,0),
\ \ \ \ \ \ \ \ \ \ (2)
\f]
on the right, top and left boundaries and 
\f[
\left. \mathbf{u}\right|_{\partial D}=(1,0),
\ \ \ \ \ \ \ \ \ \ (3)
\f]
on the bottom boundary, \f$ x_2=0.\f$
</TD>
</TR>
</TABLE>  
</CENTER>

<HR>
<HR>

\section results Results
\subsection results_CR Crouzeix-Raviart elements
The figure below shows "carpet plots" of the velocity and pressure
fields as well as a contour plot of the pressure distribution with
superimposed streamlines. The velocity vanishes along the entire
domain boundary, apart from the bottom boundary \f$ (x_2=0) \f$
where the moving "lid" imposes a unit tangential velocity which drives a
large vortex, centred at \f$ (x_1,x_2) \approx (0.62,0.26)\f$. The
discontinuity in the velocity boundary conditions creates pressure
singularities at \f$ (x_1,x_2) = (0,0)\f$ and \f$ (x_1,x_2) =
(1,0)\f$. The rapidly varying pressure in the vicinity of these points
clearly shows the discontinuous pressure interpolation employed by the
Crouzeix-Raviart elements.

\image html CrouzeixRaviart.gif "Plot of the velocity and pressure fields for Re=100 computed with Crouzeix-Raviart (Q2P-1) elements. " 
\image latex CrouzeixRaviart.eps "Plot of the velocity and pressure fields for Re=100 computed with Crouzeix-Raviart (Q2P-1) elements. " width=0.75\textwidth

\subsection results_TH Taylor-Hood elements
The next figure shows the corresponding results obtained from a 
computation with \c QTaylorHoodElements. The pressure plot illustrates
how the interpolation between the corner nodes creates
a globally continuous representation of the pressure. 

\image html TaylorHood.gif "Plot of the velocity and pressure fields for Re=100 computed with Taylor-Hood (Q2Q1) elements. " 
\image latex TaylorHood.eps "Plot of the velocity and pressure fields for Re=100 computed with Taylor-Hood (Q2Q1) elements. " width=0.75\textwidth

Note that in both simulations, the flow field is clearly
under-resolved near the ends of the "lid". In 
<A HREF="../../adaptive_driven_cavity/html/index.html">
another example</A> we will demonstrate the use of spatial
adaptivity to obtain much better solutions for this problem.

<HR>
<HR>

\section namespace Global parameters and functions
The Reynolds number is the only non-dimensional parameter
needed in this problem. As usual, we define it in a namespace:

\dontinclude driven_cavity.cc
\skipline start_of_namespace
\until end_of_namespace

<HR>
<HR>

\section main The driver code

We start by creating a \c DocInfo object to store the output directory and
the label for the output files.

\skipline start_of_main
\until doc_info.number()=0

We build the problem using \c QCrouzeixRaviartElements, solve
using the \c Problem::newton_solve() function, and document
the result before incrementing the label for the output files.

\skipline Doing QCrouzeixRaviartElements
\until end of QCrouzeixRaviartElements

Finally, we repeat the process with \c QTaylorHoodElements.

\skipline Doing QTaylorHoodElements
\until end_of_main

<HR>
<HR>

\section problem The problem class 
The \c Problem class for our steady Navier-Stokes problem is very similar to
those used for the steady scalar problems (Poisson and advection-diffusion)
that we considered in previous examples. We provide a helper function
\c fix_pressure(...) which pins a pressure value in a specified element
and assigns a specific value. 

\dontinclude driven_cavity.cc
\skipline start_of_problem_class
\until end of fix_pressure

No actions are performed after the solution is found, since the
solution is documented in \c main. However, before solving, the
boundary conditions must be set.

\skipline empty
\until end_of_actions_before_newton_solve

Finally, we provide an access function to the specific mesh and define the
post-processing function \c doc_solution(...).

\skipline Access function for the specific mesh
\until end_of_problem_class

<HR>
<HR>

\section constructor The problem constructor
Since this is a steady problem, the constructor is quite simple. We
begin by building the mesh and pin the velocities on the boundaries.

\skipline start_of_constructor
\until end loop over boundaries

Next we pass a pointer to the Reynolds number (stored in
\c Global_Physical_Variables::Re) to all elements.

\skipline Complete the build
\until end loop over elements

Since Dirichlet conditions are applied to both velocity components
on all boundaries, the pressure is only determined up to an 
arbitrary constant. We use the \c fix_pressure(...) function to
pin the first pressure value in the first element and set its value
to zero.

\skipline Now set the
\until fix_pressure

Finally, the equation numbering scheme is set up, using the function \c
assign_eqn_numbers().

\skipline Setup equation
\until end_of_constructor


<HR>
<HR>

\section doc Post-processing
As expected, this member function documents the computed solution.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>

\section comments Comments and Exercises

\subsection stress_divergence The stress-divergence form

As discussed in the introduction, by default \c oomph-lib's Navier-Stokes
elements use the stress-divergence form of the momentum equations,
<CENTER>
\f[
R_\rho Re \left( St \frac{\partial u_i}{\partial t} + 
u_j \frac{\partial u_i}{\partial x_j} \right) =
- \frac{\partial p}{\partial x_i}  +  B_i(x_j) +
R_\rho \frac{Re}{Fr} G_i +
\frac{\partial }{\partial x_j} \left[  R_\mu \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right) \right],
\f]
</CENTER>
as this form is required in problems with free surfaces
or problems in which traction boundary conditions are applied.

If the flow is divergence free (\f$ Q = 0 \Longrightarrow \partial
u_i/\partial x_i =0\f$), the viscous term may be simplified to
<CENTER>
\f[
R_\rho Re \left( St \frac{\partial u_i}{\partial t} + 
u_j \frac{\partial u_i}{\partial x_j} \right) =
- \frac{\partial p}{\partial x_i}  +  B_i(x_j) +
R_\rho \frac{Re}{Fr} G_i +  R_\mu
\frac{\partial^2 u_i }{\partial x_j^2},
\f]
</CENTER>
assuming that the viscosity ratio \f$ R_\mu \f$ remains constant

This simpler form of the equations can be used to solve problems that
do not incorporate traction boundaries or free surfaces. We illustrate
the use of these equations in 
<A HREF="../../circular_driven_cavity/html/index.html">another example.</A>

\subsection exercises Exercises

-# Compare the pressure distributions obtained with Taylor-Hood
   elements to that computed with 
   Crouzeix-Raviart elements. Why do they differ? [Hint: Consider how
   the pressure is represented in the two elements.]
-# Confirm that the velocities stored at boundary nodes must be
   pinned. Investigate what happens if you do not apply \e any 
   velocity boundary conditions [Hint: You should still be able to    
   compute a solution -- what does this solution represent?]
-# Investigate what happens when no pressure value is fixed. 

   


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/driven_cavity/">
demo_drivers/navier_stokes/driven_cavity/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/driven_cavity/driven_cavity.cc">
demo_drivers/navier_stokes/driven_cavity/driven_cavity.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


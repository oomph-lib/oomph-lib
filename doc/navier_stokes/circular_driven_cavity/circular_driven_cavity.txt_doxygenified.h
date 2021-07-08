/**

\mainpage Example problem: 2D driven cavity flow in a quarter-circle domain with spatial adaptation.

In this example we shall demonstrate
- how easy it is to adapt the code for the solution of the driven 
  cavity problem in a square domain, discussed in a
  <A HREF="../../adaptive_driven_cavity/html/index.html"> previous
  example </A>, to a different domain shape,
- how to apply body forces (e.g. gravity) in a Navier-Stokes problem,
- how to switch between the stress-divergence and the simplified forms
  of the incompressible Navier-Stokes equations.
.
<HR>
<HR>

\section example The example problem
In this example we shall illustrate the solution of the steady 2D Navier-Stokes
equations in a modified driven cavity problem: The fluid is contained
in a quarter-circle domain and is subject to gravity which acts 
in the vertical direction. We solve the problem in two different
formulations, using the stress-divergence and the simplified form of
the Navier-Stokes equations, respectively, and by applying the 
gravitational body  force via the gravity vector, \f$ {\bf G} \f$, and via the
body force  function, \f$ {\bf B}\f$, respectively.

<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Problem 1:</B>

<B>The 2D driven cavity problem in a quarter circle domain with
gravity, <BR>
using the stress-divergence form of the Navier-Stokes equations</B>

</CENTER> 
Solve
\f[
Re\phantom{i}u_j\frac{\partial u_i}{\partial x_j} =
- \frac{\partial p}{\partial x_i} +
\frac{Re}{Fr} \, G_i +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
in the quarter-circle domain \f$ D = \{x_1 \geq 0 \f$, \f$ x_2 \geq 0
\f$ and \f$ x_1^2+x_2^2 \leq 1 \} \f$,
subject to the Dirichlet boundary conditions
\f[
\left. \mathbf{u}\right|_{\partial D}=(0,0),
\ \ \ \ \ \ \ \ \ \ (2)
\f]
on the curved and left boundaries; and 
\f[
\left. \mathbf{u}\right|_{\partial D}=(1,0),
\ \ \ \ \ \ \ \ \ \ (3)
\f]
on the bottom boundary, \f$ x_2 = 0 \f$. 
Gravity acts vertically downwards so that \f$ (G_1, G_2) = (0,-1) \f$.

</TD>
</TR>
</TABLE>  
</CENTER>


When discussing the implementation of the Navier-Stokes equations in
an <A HREF="../../driven_cavity/html/index.html#stress_divergence"> 
earlier example </A>, we mentioned that \c oomph-lib allows the incompressible
Navier-Stokes equations to be solved in the simplified, rather than
the (default) stress-divergence form. We will demonstrate the use of 
this feature by solving the following problem:
<CENTER>
<TABLE>
<TR> 
<TD>

<CENTER>
<B>Problem 2:</B>

<B>The 2D driven cavity problem in a quarter circle domain with
gravity, <BR>
using the simplified form of the Navier-Stokes equations</B>

</CENTER> 
Solve
\f[
Re\phantom{i}u_j\frac{\partial u_i}{\partial x_j} =
- \frac{\partial p}{\partial x_i} +
B_i + 
\frac{\partial^2 u_i}{\partial x_j^2},
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
in the quarter-circle domain \f$ D = \{x_1 \geq 0 \f$, \f$ x_2 \geq 0
\f$ and \f$ x_1^2+x_2^2 \leq 1 \} \f$,
subject to the Dirichlet boundary conditions
\f[
\left. \mathbf{u}\right|_{\partial D}=(0,0),
\ \ \ \ \ \ \ \ \ \ (2)
\f]
on the curved and left boundaries; and 
\f[
\left. \mathbf{u}\right|_{\partial D}=(1,0),
\ \ \ \ \ \ \ \ \ \ (3)
\f]
on the bottom boundary, \f$ x_2 = 0 \f$. To make this consistent
with Problem 1, we define the body force function as
\f$ (B_1, B_2) = (0,-Re/Fr) \f$.
</TD>
</TR>
</TABLE>  
</CENTER>

Note that in Problem 2, the gravitational body force is 
represented by the body force rather than the gravity vector.

<HR>

\subsection stress_divergence Switching between the stress-divergence and the simplified forms of the Navier-Stokes equations 

The two forms of the Navier-Stokes equations differ in the
implementation of the viscous terms, which may be represented as
\f[
\frac{\partial^2 u_i}{\partial x_j^2} \mbox{ \ \ \ or \ \ \ }
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right).
\f]
For an incompressible flow, \f$ \partial u_i/\partial x_i = 0\f$, both
forms are mathematically equivalent but the stress-divergence form
is required for 
<A HREF="../../single_layer_free_surface/html/index.html">problems with free surfaces
</A>, or for
<A HREF="../../rayleigh_traction_channel//html/index.html">problems in which 
traction boundary conditions are to be applied.</A>

In order to be able do deal with both cases, \c oomph-lib's Navier-Stokes 
elements actually implement the viscous term as
\f[
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +
\Gamma_i \, \frac{\partial u_j}{\partial x_i} \right).
\f]

By default the components of the vector \f$ \Gamma_i \f$,
are set to 1.0, so that the stress-divergence form is used.
The components \f$ \Gamma_i \f$ are stored in the 
static data member
\code
static Vector<double> NavierStokesEquations<DIM>::Gamma
\endcode
of the \c NavierStokesEquations<DIM> class which forms the basis
for all Navier-Stokes elements in \c oomph-lib.
Its entries are initialised to 1.0. The user may over-write these
assignments and thus re-define the values of \f$ \Gamma \f$ being used
for a specific problem. [In principle, it is possible to use
stress-divergence form for the first component of the momentum
equations, and the simplified form for the second one, say. However,
we do not believe that this is  a particularly useful/desirable option
and have certainly never used such (slightly bizarre) assignments in
any of our own computations.] 

<HR>
<HR>

\subsection sol1 Solution to problem 1
The figure below shows "carpet plots" of the velocity and pressure
fields as well as a contour plot of the pressure distribution with
superimposed streamlines for Problem 1 at a Reynolds number of 
\f$ Re=100\f$
and a ratio of Reynolds and Froude numbers (a measure of gravity on 
the viscous scale) of \f$ Re/Fr = 100\f$. The velocity vanishes 
along the entire
domain boundary, apart from the bottom boundary \f$ (x_2 = 0) \f$
where the moving "lid" imposes a unit tangential velocity which drives
a large vortex, centred at \f$ (x_1,x_2) \approx (0.59,0.22) \f$. The
pressure singularities created by the velocity discontinuities at 
\f$ (x_1,x_2)=(0,0) \f$ and \f$ (x_1, x_2)=(1,0) \f$ are well
resolved. The pressure plot shows that away from the singularities,
the pressure decreases linearly with \f$ x_2 \f$, reflecting the
effect of the gravitational body forces which acts in the
negative \f$ x_2-\f$ direction.

\image html sol1.gif "Plot of the velocity and pressure fields for problem 1 with Re=100 and Re/Fr=100, computed with adaptive Taylor-Hood elements. " 
\image latex sol1.eps "Plot of the velocity and pressure fields for problem 1 with Re=100 and Re/Fr=100, computed with adaptive Taylor-Hood elements. " width=0.75\textwidth

<HR>

\subsection sol2 Solution to problem 2
The next figure shows the computational results for Problem 2, obtained from
a computation with adaptive Crouzeix-Raviart elements. 

\image html sol2.gif "Plot of the velocity and pressure fields for problem 2 with Re=100 and Re/Fr=100, computed with adaptive Crouzeix-Raviart elements. " 
\image latex sol2.eps "Plot of the velocity and pressure fields for problem 2 with Re=100 and Re/Fr=100, computed with adaptive Crouzeix-Raviart elements. " width=0.75\textwidth

<HR>
<HR>

\section code The code
We use a namespace \c Global_Physical_Variables to define the
various parameters: The Reynolds number,

\dontinclude circular_driven_cavity.cc
\skipline start_of_namespace
\until Re=100
 
the gravity vector \f$ {\bf G}\f$, and the ratio of Reynolds and Froude number,
\f$ Re/Fr \f$, which represents the ratio of gravitational and viscous forces,

\until Vector

In Problem 2, gravity is introduced via the body force function
\f$ {\bf B}\f$ which we define such that Problems 1 and 2 are
equivalent. (We use the gravity vector \f$ {\bf G} = (0,-1)\f$ 
to specify the direction of gravity, while indicating it magnitude
by \f$ Re/Fr.\f$)

\until }

Finally we define a body force function, which returns zero values,
for use when solving Problem 1.

\until end_of_namespace

<HR>
<HR>

\section main The driver code

First we create a \c DocInfo object to control the output, and set the maximum
number of spatial adaptations to three.

\skipline start_of_main
\until max_adapt=3

To solve problem 1 we define the direction of gravity, \f$ {\bf G} = 
(0,-1)\f$, and set the entries in the \c NavierStokesEquations<2>::Gamma
vector to \f$ (1,1)\f$, so that the stress-divergence form of
the equation is used [In fact, this step is not strictly necessary
as it simply re-assigns the default values.] 

\until Gamma[1]

Next we build problem 1 using Taylor-Hood elements and passing a
function pointer to the \c zero_body_force(...) function (defined 
in the namespace \c Global_Physical_Variables) as the argument.

\until problem(&

Now problem 1 can be solved as in 
<A HREF="../../adaptive_driven_cavity/html/index.html#main">the previous
example</A>.

\until end of problem 1

<HR>

To solve problem 2 we set the entries in the \c
NavierStokesEquations<2>::Gamma  vector to zero (thus choosing the
simplified version of the Navier-Stokes equations), define \f$ {\bf G} = 
(0,0)\f$, and pass a function pointer to the \c body_force(...)
function to the problem constructor.

\until problem(&

Problem 2 may then be solved as before.

\until end_of_main

<HR>
<HR>

\section problem The problem class 
The problem class is very similar to that used in the 
<A HREF="../../adaptive_driven_cavity/html/index.html#problem"> previous
example </A>, with two exceptions:

- We pass a function pointer to the body force function \f$ {\bf B}\f$
  to the constructor and
- store the function pointer to the body force function in the
  problem's private member data.

\dontinclude circular_driven_cavity.cc
\skipline start_of_problem_class
\until end_of_problem_class 

<HR>
<HR>

\section constructor The problem constructor
We store the function pointer to the body force function in the
private data member \c Body_force_fct_pt.

\skipline start_of_constructor
\until Body_force_fct_pt

As usual the first task is to create the mesh. We now use the \c
RefineableQuarterCircleSectorMesh<ELEMENT>, which requires the
creation of a \c GeomObject to describe geometry of the curved wall:
We choose an ellipse with unit half axes (i.e. a unit circle).

\until Wall_pt,xi_lo,

Next the error estimator is set, the boundary nodes are pinned and the
Reynolds number is assigned, as 
<A HREF="../../adaptive_driven_cavity/html/index.html#constructor">
before </A>.

\until re_pt()

Within this loop we also pass the pointers to \f$ Re/Fr \f$, the
gravity vector and the body-force function to the elements.

\until end loop over elements

The \c RefineableQuarterCircleSectorMesh<ELEMENT> contains only 
three elements and therefore provides a very coarse discretisation of
the domain. We refine the mesh uniformly twice
before pinning the redundant pressure degrees of freedom, pinning
a single pressure degree of freedom, and assigning the equation
numbers, as <A HREF="../../adaptive_driven_cavity/html/index.html#doc">
before</A>.

\until end_of_constructor

<HR>
<HR>

\section doc Post processing
The post processing function remains the same as in the 
<A HREF="../../driven_cavity/html/index.html#doc"> previous
examples </A>.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>

\section comments Comments and Exercises
-# Try making the curved boundary the driving wall [Hint: this
   requires a change in the wall velocities prescribed in \c
   Problem::actions_before_newton_solve(). The figure below shows what you
   should expect.]

\image html curved_driving_wall.gif "Plot of the velocity and pressure distribution for a circular driven cavity in which the flow is driven by the tangential motion of the curvilinear boundary. " 
\image latex curved_driving_wall.eps "Plot of the velocity and pressure distribution for a circular driven cavity in which the flow is driven by the tangential motion of the curvilinear boundary. " width=0.75\textwidth


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/circular_driven_cavity/">
demo_drivers/navier_stokes/circular_driven_cavity/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/circular_driven_cavity/circular_driven_cavity.cc">
demo_drivers/navier_stokes/circular_driven_cavity/circular_driven_cavity.cc
</A>
</CENTER>
.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


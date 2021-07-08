/**

\mainpage Example problem: Finite-Reynolds-number flow inside an oscillating ellipse

In this example we consider our first moving-boundary Navier-Stokes
problem: The flow of a viscous fluid contained in an elliptical
ring whose walls perform periodic oscillations.  
 
\c oomph-lib's Navier-Stokes elements are based on the Arbitrary 
Lagrangian Eulerian (ALE) form of the Navier-Stokes equations and can 
therefore be used in moving domain problems. In this example we 
illustrate their use in \c Domain - based meshes (first discussed in 
the example demonstrating the solution of the 
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_ALE/html/index.html">
unsteady heat equation in a moving domain</A>) in which \c
MacroElements are used to update the nodal positions in response 
to changes in the domain boundary. In subsequent examples, we 
will discuss alternative, sparse mesh update techniques that 
are useful in problems with free boundaries and in fluid-structure 
interaction problems.


<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>Finite-Reynolds-number-flow driven by an oscillating ellipse</B>
</CENTER> 

We consider the unsteady 2D flow of a Newtonian fluid that is contained in an 
oscillating elliptical ring whose wall shape is parametrised
by the Lagrangian coordinate \f$ \xi \f$ as
\f[
{\bf R}_w(\xi) = \left( 
\begin{array}{l}
a(t) \cos(\xi) \\
a^{-1}(t) \sin(\xi)
\end{array}
\right)
\f]
where 
\f[
a(t) = A + \widehat{A} \cos\left(\frac{2\pi t}{T}\right).
\f]
\f$ A \f$ represents the average half-axis of the elliptical 
ring in the \f$ x_1 \f$-direction, and \f$ \widehat{A} \f$ is the amplitude 
of its periodic variation. The ring has constant cross-sectional area
-- consistent with the incompressibility of the fluid whose motion
is governed by the ALE form of the Navier-Stokes equations, 
\f[
Re\left(St\frac{\partial u_i}{\partial t} + 
\left(u_j-u_j^M\right)\frac{\partial u_i}{\partial x_j}\right) =
- \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and the continuity equation
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
where \f$ u_j^M \f$ is the mesh velocity. We exploit the symmetry of
the problem and solve the equations in the quarter domain
\f[
D = \left\{(x_1,x_2) \, \bigg| \, 
x_1\geq0, \, x_2\geq0,\,  \left(\frac{x_1}{a(t)}\right)^2 + 
\big(x_2 \, a(t)\big)^2 \leq 1\right\},
\f]
shown in this sketch (for \f$ A=1, \ \widehat{A}=1 \f$ and \f$T=1\f$), 
\image html wall.gif "Sketch of the computational domain. " 
\image latex wall.eps "Sketch of the computational domain. " width=0.75\textwidth

The fluid is subject to no-slip boundary conditions on the curved wall, 
\f[
\left. {\bf u} \right|_{\partial D_{ellipse}} = \frac{\partial {\bf
R}_w(\xi, t)}{\partial t} 
\f]
and symmetry conditions on the symmetry boundaries, 
\f[
\left.u_1\right|_{x_2=0}=0, \qquad \left.u_2\right|_{x_1=0}=0.
\f]
The initial conditions for the velocity are given by
\f[
{\bf u}(x_1,x_2,t=0) = {\bf u}_{IC}(x_1,x_2),
\f]
where \f$  {\bf u}_{IC}(x_1,x_2) \f$ is a given, divergence-free
velocity field that satisfies the velocity boundary conditions 
at \f$ t=0 \f$. No initial conditions are required for the pressure.
</TD>
</TR>
</TABLE>  
</CENTER>

<HR>

\section solution An exact solution
It is easy to show (by inspection) that the unsteady stagnation point flow
\f[
u_1(x_1,x_2,t) = \frac{1}{a}\frac{da}{dt}\, x_1 = 
-\frac{2\pi\widehat{A}\sin\left(\frac{2\pi t}{T}\right)}
      {T\left(A+\widehat{A}\cos\left(\frac{2\pi
      t}{T}\right)\right)} \, x_1
\mbox{\ \ \ \ and \ \ \ \ \ }
u_2(x_1,x_2,t) = -\frac{1}{a}\frac{da}{dt} \, x_2 = 
-\frac{2\pi\widehat{A}\sin\left(\frac{2\pi t}{T}\right)}
      {T\left(A+\widehat{A}\cos\left(\frac{2\pi
      t}{T}\right)\right)} \, x_2, 
\f]
is an exact solution of the above problem as it satisfies the
Navier-Stokes equations and the velocity boundary conditions. 
The pressure is given by
\f[
p =
\frac{2\pi\widehat{A}\, Re\, \left(x_1^2\, St\, 
      \cos\left(\frac{2\pi t}{T}\right)A + x_1^2\, St\, \widehat{A} -
      x_1^2\widehat{A} + 
      x_1^2\widehat{A}\cos^2\left(\frac{2\pi t}{T}\right) - 
      x_2^2\, St\, \cos\left(\frac{2\pi t}{T}\right) - 
      x_2^2\, St\, \widehat{A} - 
      x_2^2\widehat{A} +
      x_2^2\widehat{A}\cos^2\left(\frac{2\pi t}{T}\right)\right)}
     {T^2\left(A^2 + 2A\widehat{A}\cos\left(\frac{2\pi t}{T}\right) +
      \widehat{A}^2\cos^2\left(\frac{2\pi t}{T}\right)\right)}.
\f]

<HR>

\section results Results

The two figures below show two snapshots of the solution for \f$ Re=Re \
St=100 \f$ , extracted from an animations of the results computed
with <A HREF="../figures/TH.avi"> Taylor-Hood</A>
and <A HREF="../figures/CR.avi"> Crouzeix-Raviart
elements</A>. In both cases, the exact solution was used as the 
initial condition for the velocities. The figures 
show "carpet plots" of the two velocity components and the 
pressure, respectively, and a contour plot of the pressure,
superimposed on the moving mesh.
The carpet plot of the velocities clearly shows that the flow is of
stagnation-point type as the horizontal velocity, \f$u_1\f$ , is a 
linear function of \f$x_1\f$ while the vertical velocity, \f$u_2\f$ , is a 
linear function of \f$ - x_2\f$. 

\image html CR.gif "Plot of the velocity and pressure fields computed with 2D Crouzeix-Raviart elements, with Re=100 and St=1. " 
\image latex CR.eps "Plot of the velocity and pressure fields computed with 2D Crouzeix-Raviart elements, with Re=100 and St=1. " width=0.75\textwidth

\image html TH.gif "Plot of the velocity and pressure fields computed with 2D Taylor-Hood elements, with Re=100 and St=1. " 
\image latex TH.eps "Plot of the velocity and pressure fields computed with 2D Taylor-Hood elements, with Re=100 and St=1. " width=0.75\textwidth

<HR>
<HR>

\section wall The moving wall
As usual, we represent the moving wall as a \c GeomObject and define
its shape by implementing the pure virtual function 
\c GeomObject::position(...). The arguments to the constructor
specify the mean half-axis of the ellipse, \f$ A \f$, the amplitude
of its variations, \f$ \widehat{A} \f$, and the period of the oscillation,
\f$ T \f$ . We also pass the pointer to a \c Time object to the
constructor and store it in a private data member, to 
allow the \c position(...) functions to access the current 
value of the continuous time. 

\dontinclude osc_quarter_ellipse.cc
\skipline start_of_MyEllipse
\until  end of MyEllipse

<HR>
<HR>

\section namespace The global parameters

As in most previous examples, we use a namespace to define and initialise
global problem parameters such as the Reynolds and Strouhal numbers:

\dontinclude osc_quarter_ellipse.cc
\skipline start_of_namespace
\until ReSt=

We also define and initialise the parameters that specify the
motion of the domain boundary and specify the exact solution.

\until end of namespace

<HR>
<HR>

\section main The driver code
As in most previous unsteady demo codes, we allow the code to be run
in a validation mode (in which we use a coarser mesh and execute
fewer timesteps). This mode is selected by specifying an (arbitrary)
command line argument that we store in the namespace 
\c CommandLineArgs.


\skipline start_of_main
\until CommandLineArgs

We create a \c DocInfo object to specify the output
directory, build the problem with adaptive Crouzeix-Raviart 
elements and the \c BDF<2> timestepper and perform the unsteady
simulation.

\until }

Then we repeat this process for adaptive Taylor-Hood elements.

\until end of main

<HR>
<HR>

\section problem The problem class 
Most of the problem class is a straightforward combination of the 
problem classes employed in the simulation of 
<A HREF="../../driven_cavity/html/index.html"> the adaptive driven cavity
</A> and 
<A HREF="../../rayleigh_channel/html/index.html"> Rayleigh channel
</A> problems, in that the problem combines unsteadiness with spatial 
adaptivity
(though in the present problem the adaptivity is only used to 
uniformly refine the very coarse base mesh; we refer to ***another
example*** for the use of full spatial adaptivity in a moving-domain
Navier-Stokes problem). 


\dontinclude osc_quarter_ellipse.cc
\skipline start_of_problem_class
\until end of actions_after_adapt

 
The key new feature in the current problem
is the presence of the moving domain which requires updates of
-# all nodal positions
-# the prescribed velocities on the 
   moving wall via the no-slip condition.
.
before every timestep. Since the nodal positions of the 
\c QuarterCircleSectorMesh are determined via its \c MacroElement
 / \c Domain representation (which updates the nodal position in
response to changes in the geometry of the \c GeomObjects that define
its boundaries), the former task may be accomplished by
executing the \c Mesh::node_update() function; the update of
the no-slip condition may be performed by calling the function
\c FSI_functions::apply_no_slip_on_moving_wall(Node* node_pt),
a helper function, defined in the namespace \c FSI_functions,
which updates the velocity components \f$ u_1, u_2 [, u_3] \f$
according to the no-slip boundary condition
\f[
{\bf u}_{Node} = \frac{\partial {\bf x}_{Node}}{\partial t}
\f] 
where the time-derivative of the nodal positions is evaluated by
the \c Node's positional timestepper. [\b Note: The function 
\c FSI_functions::apply_no_slip_on_moving_wall(...) 
assumes that the velocity components are stored in the \c Node's
first 2 [3] values. This is consistent with the storage of
the velocity component in all existing Navier-Stokes elements.
If you develop your own Navier-Stokes elements and use a different 
storage scheme you use this function at your own risk.]

Here is the implementation of these tasks:

\until actions_after_implicit_timestep()

The remaining functions are similar to those used in 
our previous Navier-Stokes examples and require no further 
explanation.

\until end of problem_class 

<HR>
<HR> 

\section constructor The problem constructor
We start by creating a timestepper of the type specified by 
the \c Problem's template parameter and add (a pointer to) it
to the \c Problem's collection of \c Timesteppers. Recall that 
this function also creates the \c Problem's \c Time object.

\until add_time_stepper_pt

Next we create the \c GeomObject that defines the curvilinear
domain boundary and pass it to the Mesh constructor.
(Since we will only use adaptivity to refine the mesh uniformly, 
it is not necessary to define an error estimator.)

\until patial_error_estimator_pt()=error_estimator_pt

Both velocity components on the curvilinear mesh boundary are
determined by the no-slip condition and must therefore be pinned,

\until end boundary 1

whereas on the symmetry boundaries only one of the two velocity
components is set to zero:

\skipline Bottom boundary: 
\until end boundary 2

Finally, we pass the pointers to \f$ Re \f$, \f$ Re \ St \f$
and the global \c Time object (automatically created by the \c Problem when
the timestepper was passed to it at the beginning of the constructor)
to the elements, pin the redundant nodal pressure degrees of freedom
(see the discussion of the <A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">adaptive driven-cavity problem</A> for more
details), pin one pressure degree of freedom,
and set up the equation numbering scheme.

\until end of constructor

<HR>
<HR>

\section IC Assigning the initial conditions
This function assigns "history values" for the velocities and the
nodal positions from the exact solution. It is implemented in exactly the
same way as in the solution of the
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_ALE/html/index.html#IC">
unsteady heat equation in a moving domain</A>. Note that because
the domain is moving, the nodal positions must be updated
(according to the position of the domain boundary at the relevant previous
timestep), before evaluating the exact solution at the nodal position.

\skipline start_of_set_initial_condition
\until end of set_initial_condition

<HR>
<HR>

\section doc Post processing
The function \c doc_solution(...) is similar to that in the 
<A HREF="../../../unsteady_heat/two_d_unsteady_heat/html/index.html#doc">
unsteady heat examples </A> and the previous Navier-Stokes examples. 
We add dummy zones and tecplot geometries to facilitate the
post-processing of the results with tecplot.
 
\skipline start_of_doc_solution
\until end of doc_solution 

<HR>
<HR>

\section unsteady_run The timestepping loop
The timestepping loop is extremely straightforward: We choose a
timestep and the overall length of the simulation, initialise the
timestepper(s) by calling \c Problem::initialise_dt(...) and
assign the initial condition. 

\skipline start_of_unsteady_run
\until set_

Next we set the number of timesteps for a normal run.

\until nstep =

We over-write this number and perform a single uniform mesh refinement
if the code is run in self-test mode (indicated by a non-zero number 
of command line arguments),

\until }

otherwise we refine the mesh three times and output the initial conditions
 
\until doc_solution

Finally we execute the proper timestepping loop and document
the solution after every timestep

\until end of unsteady_run

<HR>
<HR>

\section comments Comments and Exercises
-# Compare the results of the numerical simulation in which
   \f$ {\bf u}_{IC} \f$ is given by the exact solution (an unsteady
   stagnation point flow) to that obtained from an "impulsive start"
   where \f$ {\bf u}_{IC} = {\bf 0}. \f$  (This is most easily 
   implemented by replacing the call to \c set_initial_condition() 
   with a call to \c Problem::assign_initial_values_impulsive().
   \n \n
   Why do we obtain the same velocity with both initial conditions
   and why does the pressure take a few timesteps (How many exactly?
   Compare simulations with \c BDF<4> and \c BDF<2> timesteppers.)
   to "catch up" with the exact solution? [Hint: The 
   unsteady stagnation point flow is a potential flow, therefore 
   the viscous terms in the Navier-Stokes equations disappear. See 
   also chapter 3.19 in Volume 2 of Gresho & Sani's wonderful book 
   "Incompressible Flow and the Finite Element Method".]




<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/osc_ellipse/">
demo_drivers/navier_stokes/osc_ellipse/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/osc_ellipse/osc_quarter_ellipse.cc">
demo_drivers/navier_stokes/osc_ellipse/osc_quarter_ellipse.cc
</A>
</CENTER>
.











<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Example problem: Spatially adaptive solution of the 2D unsteady heat equation with flux boundary conditions.


This is a slightly more advanced example in which we demonstrate
the use of spatial adaptivity in time-dependent problems.
We discuss the implementation of the spatially adaptive version of
\c oomph-lib's unsteady Newton solver, \c
Problem::unsteady_newton_solve(...), and explain why the assignment
of initial conditions should be performed by overloading the 
\c Problem::set_initial_condition() function. We also
discuss briefly how \c oomph-lib's generic dump and restart
functions deal with adaptive meshes.

For this purpose consider the following problem:


<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The two-dimensional unsteady heat equation with flux boundary
conditions in a quarter circle domain</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2\frac{\partial^2 u}{\partial x_i^2} 
= \frac{\partial u}{\partial t} + f\left(x_1,x_2,t\right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the quarter-circle domain \f$ D \f$, bounded by the coordinate axes
and the unit circle, subject to Neumann boundary conditions,
\f[
\left. \frac{\partial u}{\partial n}\right|_{\partial D_{Neumann}}=
- \left. \frac{\partial u}{\partial x_2}\right|_{\partial D_{Neumann}}=
g_0, 
\ \ \ \ \ \ \ \ \ \ (2)
\f]
along the horizontal domain boundary \f$ \partial D_{Neumann} = 
\{ (x_1,x_2) | x_1 \in [0,1], x_2=0 \} \f$,
and to Dirichlet boundary conditions,
\f[
\left. u\right|_{\partial D_{Dirichlet}}=h_0, 
\ \ \ \ \ \ \ \ \ \ (3)
\f]
elsewhere.
\image html domain.gif "Sketch of the domain and the boundary conditions. " 
\image latex domain.eps "Sketch of the domain and the boundary conditions. " width=0.75\textwidth
The initial conditions are given by
\f[
u(x_1,x_2,t=0)=k_0(x_1,x_2),
\ \ \ \ \ \ \ \ \ \ (4)
\f]
where the functions \f$ f, g_0, \ h_0\f$ and  \f$ k_0\f$ are given. 
</TD>
</TR>
</TABLE>  
</CENTER>

We choose the functions \f$ f, g_0, \ h_0\f$ and  \f$ k_0\f$ 
so that
\f[
u_0(x_1,x_2,t) = \tanh\bigg[1-\alpha\bigg(\tan\Phi
\big(x_1-\beta\tanh[ \gamma\cos\left(2\pi
t\right)]\big)- x_2\bigg)\bigg]
 \ \ \ \ \ \ \ \ \ \ (5)
\f]
is the exact solution.


The solution represents the "usual" tanh profile, whose steepness
is controlled by the parameter \f$ \alpha \f$ so that for \f$ \alpha
\gg 1 \f$ the solution approaches a step. The step is oriented
at an angle \f$ \Phi \f$ against the \f$ x_1-\f$ axis and its position varies
periodically. The parameter \f$ \beta \f$ controls the amplitude
of the step's lateral displacement, while \f$ \gamma \f$ determines the
rate at which its position changes. For \f$ \gamma \gg 1 \f$ ,
the step remains stationary for most of the period and then
translates rapidly parallel to the \f$ x_1- \f$ axis, making this a very 
challenging problem.





The figure below shows a snapshot of the <A HREF="../figures/step_soln.avi">
animated solution</A>, obtained from the spatially adaptive simulation
discussed below, for the parameter values \f$ \alpha=10, \ 
\Phi=45^o, \ \beta=0.3, \ \gamma=5. \f$

\image html step_soln.gif "Snapshot of the solution. " 
\image latex step_soln.eps "Snapshot of the solution. " width=0.75\textwidth


  The mesh adaptation in response to the translation of the step can be
seen more clearly in this contour plot, taken from
<A HREF="../figures/unsteady_heat_contour.avi">
another animation of the solution</A>.

\image html unsteady_heat_contour.gif "Contour plot of the solution. " 
\image latex unsteady_heat_contour.eps "Contour plot of the solution. " width=0.75\textwidth

<HR>
<HR>



\section spatial_adapt Background: Spatial adaptivity in time-dependent problems

Enabling spatial adaptivity in time-dependent problems involves
essentially the same steps as for steady problems:
- The domain must be discretised with a mesh that is derived
  from the \c RefineableMesh base class.
- An \c ErrorEstimator object must be created and passed to the
  mesh.
- The empty virtual functions \c Problem::actions_before_adapt() and 
  \c Problem::actions_after_adapt() may be overloaded to 
  perform any actions that are required before or after the mesh
  adaptation, such as the deletion or recreation of any \c
  FaceElements that apply flux boundary conditions.
.
Once these steps have been performed, a spatially adaptive
solution can be computed with a three-argument version of \c oomph-lib's
unsteady Newton solver \c Problem::unsteady_newton_solve(...):

\code 
Problem::unsteady_newton_solve(dt,max_adapt,first)
\endcode
  
The arguments to this function are as follows:
- The \c double \c dt specifies the (fixed) timestep.
- The \c unsigned \c max_adapt specifies the maximum number of 
  spatial adaptations allowed.
- The \c bool \c first indicates if the first timestep is performed.
  This argument is required to allow the automatic re-assignment of
  the initial conditions following any mesh adaptations during the
  computation of the first timestep.
. 
Given these arguments, the unsteady Newton solver solves the
non-linear system of spatially and temporally discretised equations 
to advance the solution from time  \f$ t \f$ to \f$ t + dt \f$ . 
Once the solution  at time \f$ t + dt \f$ has been obtained, 
error estimates are computed for all elements. If any 
elemental error estimates are outside the target range, the solution
is rejected  and the mesh is adapted. In the course of mesh adaptation
the existing solution (the nodal values \e and the history
values) at time  \f$ t \f$ are interpolated onto the new mesh
 before recomputing the solution.
This process is repeated until the error estimates are within
the target range, or until the maximum number of adaptations, specified
by the parameter \c max_adapt, is exceeded, just as in the steady case.

Here is an illustration of the procedure for a 1D problem:

\image html adapted.gif "Sketch of the mesh adaptation for time-dependent problems. " 
\image latex adapted.eps "Sketch of the mesh adaptation for time-dependent problems. " width=0.75\textwidth

 This procedure is the obvious generalisation of the procedure for steady
problems. However, in time-dependent problems two additional 
issues arise:
-# In a steady problem the interpolation of the solution onto the
   adapted mesh (step 4 in the above sketch) merely serves to 
   provide an initial guess for the solution on the refined 
   mesh. It is irrelevant if the interpolation from the coarse mesh
   provides a poor approximation of the actual solution as 
   the solution is completely recomputed anyway. 
   \n\n
   In an unsteady problem, we also have to interpolate the history 
   values (the solution at previous timesteps in a BDF scheme) onto 
   the adapted mesh. Their values are \e not
   changed when the solution is advanced from time \f$ t \f$ to \f$ t+ dt. \f$
   In time-dependent problems, the benefit of repeated mesh
   adaptations (i.e. \c max_adapt > 1) is therefore 
   limited by the fact that mesh refinement cannot improve their 
   accuracy -- the history values are always given by
   the (possibly poor) approximations obtained by interpolation 
   from the coarser mesh employed at the previous timestep. We
   therefore recommend limiting the number of spatial adaptations to
   \c max_adapt = 1. We stress that, in practice, this is not a serious
   restriction because the time-integration procedure will only provide
   (temporally) accurate results if the timestep \c dt is so small that the
   solution at time \f$ t \f$ only differs slightly from that
   at time  \f$ t + dt \f$ . One level of mesh adaptation per timestep 
   should therefore be sufficient to adapt the mesh
   in response to these changes.
-# The only exception to this recommendation arises during the computation of
   the first timestep, illustrated in the following sketch:
\image html adapted_ic1.gif "Sketch of the mesh adaptation during the computation of the first timestep. " 
\image latex adapted_ic1.eps "Sketch of the mesh adaptation during the computation of the first timestep. " width=0.75\textwidth
   When computing the first timestep, the solution on the initial mesh
   will have been created by assigning the nodal values according to the
   analytical initial condition
   (4). If the initial mesh is very coarse
   (as it should be), the finite-element representation of the
   initial condition is likely to be very poor, as shown in the above
   sketch. Clearly, the interpolation
   from the coarse onto the fine mesh cannot recover any small-scale
   features in the initial condition that were missed by its 
   representation on the coarse mesh. It is therefore better to
   re-assign the initial condition (the values \e and the history values!) 
   on the adapted mesh, as shown in this sketch:
\image html adapted_ic2.gif "Sketch of the modified mesh adaptation during the computation of the first timestep. " 
\image latex adapted_ic2.eps "Sketch of the modified mesh adaptation during the computation of the first timestep. " width=0.75\textwidth
   With this procedure, repeated mesh adaptations will
   improve the accuracy of the solution, therefore much larger 
   values of \c max_adapt can (and should!) be specified
   when the first timestep is computed.
   The unsteady Newton solver \c Problem::unsteady_newton_solve(...)
   performs the revised procedure if the boolean argument
   \c first is set \c true. In that case,
   the values and history values on the adapted mesh are 
   (re-)assigned by calling the function
   \code
   Problem::set_initial_condition()
   \endcode
   which is defined as an 
   empty virtual function in the \c Problem base class. You should
   overload it in your derived \c Problem to ensure that your specific
   initial conditions are assigned
   by the mesh adaptation procedures. [In fact, the function  \c
   Problem::set_initial_condition() is not quite empty -- not 
   re-setting the initial condition when performing mesh adaptations 
   during the first timestep of a time-dependent simulation seems "so 
   wrong" that the function issues a warning message. Although the
   overloading of this function is not strictly necessary if the initial
   conditions can be represented exactly by the interpolation from the
   coarse mesh onto the fine mesh, we consider it good practice to
   do so, for reasons discussed in <A HREF="../../../axisym_navier_stokes/spin_up/html/index.html#good_practice_ics">another tutorial</A>.]
.
 

<HR>
<HR>

\section overview Overview of the driver code 

Equipped with this background information, the driver code for our
example problem is easy to understand, if somewhat 
lengthy. [Using an example with Dirichlet boundary conditions 
along the entire domain boundary would have shortened
the code significantly but we deliberately chose an example with
Neumann boundary conditions to demonstrate that the functions 
\c Problem::actions_before_adapt()
and \c Problem::actions_after_adapt() may be used exactly as in
the steady computations.] We will not discuss the
methodology for applying flux-type boundary conditions in problems
with spatial adaptivity in detail, but refer to the discussion provided
in the <A HREF="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
 earlier steady example.</A>

Overall, the code is a straightforward combination of the driver code for
the <A HREF="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
steady Poisson problem with flux boundary conditions and spatial
adaptivity</a> and the driver code for the 
<A HREF="../../two_d_unsteady_heat/html/index.html"> unsteady heat equation
without spatial adaptivity. </A> 

<HR>
<HR>

\section namespace Global parameters and functions

As usual, we store the problem parameters in a namespace,
\c TanhSolnForUnsteadyHeat, 
in which we also specify the source function, the prescribed flux 
along the Neumann boundary and the exact solution.

<HR>
<HR>
 
\section ellipse Representing the curvilinear domain boundary by a  GeomObject

As discussed <A HREF="../../../poisson/fish_poisson2/html/index.html">
elsewhere,</A> \c oomph-lib's mesh adaptation procedures require
curvilinear domain boundaries to be represented by
\c GeomObjects which describe the object's shape via their member 
function \c GeomObject::position(...). This function exists
in two versions: 
- The two-argument version,  \c GeomObject::position(xi,r)
  computes the position vector, \c r, to the point on/in the \c GeomObject,
  parametrised by the vector of intrinsic coordinates, \c xi.
- The three-argument version \c GeomObject::position(t,xi,r), where
  \c t is an \c unsigned, computes the position vector at 
  the \c t - th previous timestep. 
.
In the current problem, the domain boundary
is stationary, therefore the steady and unsteady versions of the
function are identical. Here is the complete source code for
the \c MyUnitCircle object which we will use to represent the 
curvilinear domain boundary:

\dontinclude two_d_unsteady_heat_adapt.cc
\skipline start_of_MyUnitCircle
\until end of MyUnitCircle


<HR>
<HR> 


\section main The main function

As before, we use command line arguments to (optionally) specify 
a restart file. We store the command line arguments in
the namespace \c CommandLineArgs and build the \c Problem object,
passing the pointer to the source function. Next we specify the
time-interval for the simulation and set the error targets for 
the spatial adaptation.

\dontinclude two_d_unsteady_heat_adapt.cc
\skipline start_of_main
\until min_permitted_error

We create and initialise the boolean flag that indicates if the first timestep
is computed, and choose a large initial value for the number of permitted mesh 
adaptations. We then assign the initial conditions on the coarse
initial mesh and retrieve the timestep (chosen when the
initial conditions are assigned in \c set_initial_condition() ) 
from the problem's \c Time object.

\skipline Set IC
\until ->dt

If the simulation has been restarted, the first timestep is not the step
at which the initial condition has to be assigned, therefore
we reset the \c first and \c max_adapt parameters to their
appropriate values. If the run is not restarted, the 
problem will have been built with a very coarse initial mesh
(comprising just three elements). We don't need an error estimator
to tell us that this is too coarse to represent the solution 
accurately and apply two levels of uniform refinement 
before solving the problem. Note that we refine the entire 
problem, not just the mesh to ensure that 
\c Problem::actions_before_adapt() and \c
Problem::actions_after_adapt() are executed and the
equation numbering scheme is re-generated. \c Problem::refine_uniformly()
also interpolates the solution from the coarse initial mesh onto the
refined mesh but, as discussed above, this will lead to a very poor
representation of the initial condition. Therefore we re-assign
the initial condition on the refined mesh and document the
finite-element representation of the initial condition.
 
\skipline If restart:
\until doc_solution

The time-stepping loop itself is very similar to that used in the
<A HREF="../../two_d_unsteady_heat/html/index.html">
example without spatial adaptivity</A>. Here we call the
three-argument version of the unsteady Newton solver
\c Problem::unsteady_newton_solve(...) and re-set the 
parameters \c max_adapt and \c first to their appropriate values
once the first step has been performed.

\dontinclude two_d_unsteady_heat_adapt.cc
\skipline  Find number of steps
\until end of main


<HR>
<HR>

\section problem The problem class
As discussed above, the problem class mainly contains verbatim copies of
the member functions in the corresponding 
 <A HREF="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
steady </a> and 
<A HREF="../../two_d_unsteady_heat/html/index.html">unsteady</A> 
problems:

\dontinclude two_d_unsteady_heat_adapt.cc
\skipline start_of_problem_class
\until // end of problem_class


<HR>
<HR>

\section constructor The problem constructor
The problem constructor combines the constructors of the 
 <A HREF="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
steady </a> and 
<A HREF="../../two_d_unsteady_heat/html/index.html">unsteady</A> 
problems. We start by creating a \c DocInfo object to control
the output, set the parameters for the exact solution and create
the \c TimeStepper:

\skipline start_of_constructor
\until add_time_stepper


We create the \c GeomObject that describes the curvilinear
domain boundary and pass it to the mesh constructor:

\skipline Setup mesh
\until Boundary_pt,xi_lo,fract_mid

Next, we create the surface mesh that contains the prescribed flux
elements and combine the two submeshes to the \c Problem's global
mesh. We create an instance of the \c Z2ErrorEstimator
and pass it to the bulk mesh. 


\skipline Create the surface
\until Bulk_mesh_pt->spatial_error_estimator_pt()

We pin the nodal values on the Dirichlet boundaries and select the
central node in the unrefined three-element mesh as the control node
at which the solution is documented in the trace file.


\skipline Set the boundary conditions for 
\until Doc_node_pt=el0_pt->node_pt

Finally, we complete the build of all elements by passing the relevant
function pointers to the elements, and assign the equation numbers.

\skipline Complete the build of all elements
\until end of constructor

<HR>
<HR>

\section old_functions Other member functions
The remaining member functions
- \c  actions_after_newton_solve()
- \c  actions_before_newton_solve()
- \c  actions_after_implicit_timestep()
- \c  actions_before_implicit_timestep()
- \c  actions_before_adapt()
- \c  actions_after_adapt()
- \c  set_initial_condition()
- \c  create_flux_elements(...)
- \c  delete_flux_elements(...)
- \c  doc_solution()
- \c  dump_it(...)
- \c  restart(...)
.
are identical (or at least extremely similar) to those in previous 
examples, so we do not list them here. You can examine the functions in detail
in the source code
<A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat_adapt/two_d_unsteady_heat_adapt.cc">two_d_unsteady_heat_adapt.cc</A>.


<HR>
<HR>


\section restart Dump/restart with spatial adaptivity

It is worth examining the
dump and restart functions, however, as they demonstrate that
the generic versions defined in the \c Problem base class
can also deal with adaptive problems -- a non-trivial task!

\skipline start_of_dump_it
\until end of restart

Since details of their implementation are hidden from the user, we
briefly comment on the various tasks performed by these functions.
The main task of the \c Problem::read(...) function is to 
read values (and history values) of all \c Data objects from a file and to
assign these values to the appropriate \c Data (and \c Node) objects in the
\c Problem. This assumes that the \c Problem's 
constituent \c Meshes, elements, \c Nodes and \c Data objects have been 
created, and that the \c Problem's various pointer-based lookup schemes
access them in the order they were in when the \c Problem was dumped
to the restart file. In a non-adaptive computation, the number of 
elements and the number of \c Data objects remain constant throughout the 
simulation and the \c Problem::read(...) function can be called 
as soon as the \c Problem has been built -- usually by its
constructor. (The \c Problem constructor always builds and enumerates 
its constituent objects in the same order.)

In a simulation with spatial adaptivity the number of elements,  \c
Nodes and \c Data objects varies throughout the computation. 
It is therefore necessary to re-generate the \c Problem's refinement pattern 
before the \c Data values can be read from the restart file. This is achieved
(internally) by calling the function \c
RefineableMesh::dump_refinement(...) for all refineable meshes
before the \c Data is dumped. This function writes the \c Mesh's
refinement pattern to the restart file, using a format that
can be read by the corresponding member function 
\c RefineableMesh::refine(...) which adapts an unrefined mesh so that
its topology and the order of its \c Nodes and elements is recreated. 


<HR>
<HR>


\section comments Comments and exercises

The plots below show the time history of various parameters.
- The upper graph  compares the solution at the control node 
  (red line) against the exact solution (green line). 
- The middle graph shows the position of the step by plotting its
  intercept with the \f$ x_1- \f$ axis as a function of time, and
  the error of the solution. 
- The lower graph illustrates the evolution of the adaptive spatial 
  refinement process: The green line illustrates the total number of 
  elements; the blue and red lines show the number of elements that 
  are refined and unrefined at each timestep.

\image html trace.gif "Time history of the solution. " 
\image latex trace.eps "Time history of the solution. " width=0.75\textwidth


The plots illustrate clearly how the mesh is adapted as the step
moves through the domain -- the peaks in the number of refined/unrefined 
elements per timestep coincide with the periods during which the step moves 
very rapidly. The increase in the error during these phases is mainly 
due to the temporal error -- the <A HREF="../figures/step_soln.avi">
animation</A> shows that the computed solution lags behind
the exact one. We will address this by adding adaptive time-stepping
in <A HREF="../../two_d_unsteady_heat_2adapt/html/index.html">
another example.</A>


\subsection ex Exercises
-# Confirm that the error during the periods of rapid change in the
   solution is due to the temporal error by repeating the 
   simulation with a smaller/larger timestep and/or a time-stepping scheme
   with higher/lower order (e.g. BDF<1> or BDF<4>).
-# Assess the importance of re-assigning the initial conditions when
   spatial adaptations are performed during the computation of the
   first timestep.
   -# Compare the finite-element representation of the initial
      condition(contained in the file \c
      RESLT/soln0.dat) against that obtained when the 
      re-assignment of the initial conditions 
      after the two calls to \c problem.refine_uniformly()
      in the \c main function is suppressed.
   -# Comment out the calls to \c problem.refine_uniformly()
      and set \c first=false throughout the \c main function and
      compare the computed results against those obtained with  
      the correct procedure.

 



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="
../../../../
demo_drivers/unsteady_heat/two_d_unsteady_heat_adapt/
">
demo_drivers/unsteady_heat/two_d_unsteady_heat_adapt/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="
../../../../
demo_drivers/unsteady_heat/two_d_unsteady_heat_adapt/two_d_unsteady_heat_adapt.cc
">
demo_drivers/unsteady_heat/two_d_unsteady_heat_adapt/two_d_unsteady_heat_adapt.cc
</A>
</CENTER>
.











<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


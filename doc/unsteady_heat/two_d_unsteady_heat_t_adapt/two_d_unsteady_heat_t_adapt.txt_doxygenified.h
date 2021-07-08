/**

\mainpage Example problem: Solution of the 2D unsteady heat equation with temporal adaptivity

This example illustrates \c oomph-lib's adaptive timestepping
capabilities. We consider, yet again, the 2D unsteady 
heat equation
<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The two-dimensional unsteady heat equation in a square domain.</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2\frac{\partial^2 u}{\partial x_i^2} 
= \frac{\partial u}{\partial t} + f\left(x_1,x_2,t\right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the square domain \f$ D = \{x_i \in [0,1]; i=1,2 \} \f$,
subject to the Dirichlet boundary conditions
\f[
\left. u\right|_{\partial D}=g_0
\ \ \ \ \ \ \ \ \ \ (2)
\f]
and initial conditions
\f[
u(x_1,x_2,t=0)=h_0(x_1,x_2),
\ \ \ \ \ \ \ \ \ \ (3)
\f]
where the functions \f$ g_0\f$ and  \f$ h_0\f$ are given.
</TD>
</TR>
</TABLE>  
</CENTER>

Here we choose the forcing function and the boundary
and initial conditions so that
\f[
u_0(x_1,x_2,t) = \frac{1}{2} \bigg( 1 + \tanh\big(\gamma \cos(2\pi
t)\big)\bigg) \sin\left(K \left( x_1 \cos \Phi + x_2 \sin
\Phi\right)\right),
 \ \ \ \ \ \ \ \ \ \ (4)
\f]
is the exact solution. The solution represents a 1D sin 
profile with wavenumber \f$ K \f$, rotated against 
the \f$x_1\f$-axis by an angle \f$ \Phi \f$, and modulated 
by a time-periodically varying amplitude. The parameter \f$ \gamma \f$
controls the rate at which the amplitude changes. For large  
values of \f$ \gamma \f$, the amplitude remains constant for
most of the period but changes rapidly at 
\f$ t = \left( 1/4 + i/2 \right), i=0,1,2,... \f$.
To resolve these rapid changes accurately, very small timesteps
are required. Conversely, relatively large timesteps could be employed
during the phases when the solution remains approximately constant.
The problem therefore presents an ideal test case for 
an adaptive timestepping scheme.
 
  The figure below shows a snapshot from the
<A HREF="../figures/unsteady_heat_soln.avi">animation of
the exact and computed
solutions,</A> obtained from a simulation with \c oomph-lib's 
adaptive \c BDF<2> timestepper. Since the time interval between 
subsequent frames in this animation varies, each frame shows a blue line (in 
the top left corner) whose length is proportional to the elapsed time.
The different rate at which the length of the line increases 
reflects the fact that smaller timesteps are taken when the 
solution varies rapidly.

\image html unsteady_heat_soln.gif "Snapshot from the animation of the exact and computed solutions. " 
\image latex unsteady_heat_soln.eps "Snapshot from the animation of the exact and computed solutions. " width=0.75\textwidth


  The figure below shows the time trace of the solution and 
documents the timesteps chosen by the adaptive timestepping scheme.
Note how smaller timesteps are chosen when the solution undergoes
rapid changes.

\image html trace.gif "Upper plot: Time evolution of the computed and exact solutions at a control node. Lower plot: The timestep dt and the norm of the error. " 
\image latex trace.eps "Upper plot: Time evolution of the computed and exact solutions at a control node. Lower plot: The timestep dt and the norm of the error. " width=0.75\textwidth
 

Most of the driver code for this example is identical to that 
discussed in the 
<A HREF="../../two_d_unsteady_heat/html/index.html">previous example,</A>
therefore we only discuss the modifications required to enable
temporal adaptivity:
- We pass a boolean flag
  to the constructor of the \c BDF<2> timestepper.
- We define a global error norm for the adaptive timestepper by overloading
  the (empty) virtual function \c Problem::global_temporal_error_norm().
- We replace the call to \c Problem::unsteady_newton_solve(...) by
  a call to \c Problem::adaptive_unsteady_newton_solve(...) and specify
  a target for the temporal error norm.  
.
<HR>
<HR> 

\section namespace Global parameters and functions

We store the problem parameters and define the source function 
and the exact solution in the usual namespace. 

\dontinclude two_d_unsteady_heat_t_adapt.cc
\skipline start_of_ExactSolnForUnsteadyHeat
\until end of ExactSolnForUnsteadyHeat


<HR>
<HR>

\section main The driver code

Temporal adaptivity only requires a few straightforward changes
to the time-stepping loop. Since the number of timesteps required to
reach the end of the simulation is not known a priori, we
replace the \c for - loop over the fixed number of timesteps, employed
in the <A HREF="../../two_d_unsteady_heat/html/index.html">
non-adaptive version of the code,</A> by a \c while - loop that continues the
time-integration until \f$ t \ge  t_{max}\f$. 
 
  The adaptive timestepper, \c Problem::adaptive_unsteady_newton_solve(...)
takes two arguments. The first one is a suggestion for the
size of the next timestep; the second specifies the target error
The adaptive timestepper automatically adjusts the timestep
until the error estimate computed by \c Problem::global_temporal_error_norm()
is less than that target. If the error estimate for the
solution computed with the suggested value of \c dt is too large,
\c dt is reduced by a factor of 2 and the solution is recomputed.
This process is repeated until
- the estimated error has become sufficiently small 
.
or 
- \c dt has been reduced below a threshold. 
.
The threshold is stored
in the private member data \c Problem::Minimum_dt and is initialised
to \f$ 10^{-12} \f$. This default can be changed with the
access function \c Problem::minimum_dt(). It is also possible to
specify a maximum value for the timestep by overwriting the default
value \f$ 10^{12} \f$ for the corresponding data member,
 \c Problem::Maximum_dt, using the function  \c Problem::maximum_dt().

 The adaptive unsteady Newton solver returns a suggestion for 
the size of the next timestep.
 
  Here is the relevant code fragment from the otherwise unchanged
\c main function:

\dontinclude two_d_unsteady_heat_t_adapt.cc
\skipline  Target error for adaptive timestepping
\until // end of timestepping loop

The rest of the main function is identical
to that in the <A HREF="../../two_d_unsteady_heat/html/index.html">
previous, non-adaptive example.</A>

<HR>
<HR>

\section problem The problem class 
 
The problem class contains a single additional member function

\dontinclude two_d_unsteady_heat_t_adapt.cc
\skipline  Global error norm
\until global_temporal_error_norm();

which we will discuss below.

<HR>
<HR>

\section constructor The problem constructor

The  problem constructor is almost identical to that in the
<A HREF="../../two_d_unsteady_heat/html/index.html">previous
example.</A> The only difference is that the boolean flag \c true is
passed to the constructor of the \c BDF<2> timestepper, which means that
a predictor step is computed for each timestep, see <a href="#back">
background </a> for more details.

<HR>
<HR>

\section destructor The problem destructor

The  problem destructor is identical to that in the 
<A HREF="../../two_d_unsteady_heat/html/index.html">previous example.</A>


<HR>
<HR>

\section before_timestep Actions before timestep
This function is identical to that in the 
<A HREF="../../two_d_unsteady_heat/html/index.html">previous example.</A>


<HR>
<HR>

\section IC Set initial condition
This function
is identical to that in the
<A HREF="../../two_d_unsteady_heat2/html/index.html">previous example.</A>

<HR>
<HR>

\section doc Post-processing

The Problem member function \c doc_solution(...) is virtually identical to that
in the  <A HREF="../../two_d_unsteady_heat2/html/index.html">previous 
example.</A> We merely add the timestep \c dt to the trace file.

<HR>
<HR>

\section dump Dumping the solution

This function is identical to that in the 
<A HREF="../../two_d_unsteady_heat2/html/index.html">previous example,</A>
indicating that the generic \c Problem::dump() function can
deal with time-dependent simulations.

 
<HR>
<HR>

\section read Reading a solution from disk

This function is identical to that in the 
<A HREF="../../two_d_unsteady_heat2/html/index.html">previous
example, </A>indicating that the generic \c Problem::read() function can
deal with time-dependent simulations.


<HR>
<HR>

\section error Defining the global error norm for the adaptive timestepper

\subsection back Background

\c oomph-lib's adaptive timesteppers employ a predictor-corrector
scheme to control the size of the timestep. In these schemes a 
low-order explicit timestepper
is used to predict the solution at the next timestep. 
This prediction is compared to the solution computed with the 
actual (usually implicit) timestepper itself. The difference
between the two predictions is then used to derive an estimate of 
the error (exploiting the different truncation errors of the two
timestepping schemes). 


Interfaces to the functions that compute the temporal error estimates 
are defined as broken virtual function in the \c TimeStepper base class.
Specific \c TimeSteppers that allow adaptive timestepping  
therefore overload the broken virtual function


\code
TimeStepper::temporal_error_in_value(data_pt,i)
\endcode

which computes an estimate of the error of \c i - th value stored
in the \c Data object pointed to by \c data_pt. In free-boundary 
problems in which the position of the nodes
is an unknown, the corresponding function 

\code 
TimeStepper::temporal_error_in_position(node_pt,i)
\endcode

may be used to obtain an estimate of the error in the \c i - th nodal 
coordinate of the  \c Node pointed to by \c node_pt.


These individual error estimates must be combined into 
a problem-specific, scalar error norm, \f${\cal E}\f$ say,
which forms the basis of the adaptive adjustment of the timestep in 
\c Problem::adaptive_unsteady_newton_solve(...).

\subsection impl Implementation

In the present problem, we choose the RMS of the errors
 \f$ e_j \ (j=1,...,N_{node})\f$ at the nodes as the 
global error norm
\f[
{\cal E} = \sqrt{ \frac{1}{N_{node}} \sum_{j=1}^{N_{node}} e_j^2 }.
\f]
This may be implemented in a few lines of code:

\dontinclude two_d_unsteady_heat_t_adapt.cc
\skipline start_of_global_temporal_error_norm
\until end of global_temporal_error_norm



<HR>
<HR>

\section comments Comments and Exercises
As demonstrated above, enabling temporal adaptivity
for a given problem is extremely straightforward as it only
requires the implementation of the problem-specific 
function \c Problem::global_temporal_error_norm().
In most cases, the RMS of the nodal errors (or some suitable
generalisation for vector-valued problems) is an obvious
choice. In free-boundary problems, the RMS of the error estimate for 
the nodal positions is often a useful error measure. 

<HR>

\subsection norm How to choose the target for the temporal error norm
 Having decided on an error norm, how does one choose the target
for the error norm? The answer is the same as in a simulation with
a fixed timestep: Trial and error, followed by careful convergence 
tests. We usually employ the following strategy: 
- Implement the function \c Problem::global_temporal_error_norm()
  and perform an initial computation with a fixed timestep, choosing
  its size heuristically, exploiting prior knowledge that we (usually!)
  have about our problem. For instance, if we expect a periodic solution 
  with an approximate period \f$ T \f$, we may start with a fixed timestep 
  \f$ dt = T/20,\f$ say. The computed results are likely to be very
  inaccurate but the time-trace will usually reveal the
  characteristic features of the solution and thus identify phases 
  during which the solution varies rapidly.  
- Now repeat the simulation with a smaller time-step, \f$ dt =
  T/40,\f$ say, and check if any new features develop in the
  time-trace. If the time-trace appears to be robust, monitor the 
  temporal error norm, e.g. by including the output of 
  \c Problem::global_temporal_error_norm() into the trace file.
- Use the maximum of the temporal error norm observed during the
  simulation with the fixed timestep as the target in a first
  simulation with temporal
  adaptivity. The timestepper should now increase the size of the
  timestep in regions where the error estimate was small.
- Now repeat the simulation with smaller and smaller target errors
  until further reductions do not lead to further changes in the
  computed results. 

\subsubsection ex Exercise:
Employ the above procedure to determine the target error \f$ \hat{\cal
E}\f$ required for the computed solution to be graphically
indistinguishable from that obtained with a target error of
\f$ \hat{\cal E}/2\f$.

<HR>

\subsection restarts Restarting from a simulation with temporal adaptivity
We mentioned above that the data recorded/read by the generic 
\c Problem::dump(...) and \c Problem::read(...) functions is sufficient to 
restart a temporally adaptive simulations as the functions record the
history values and the history of previous timesteps. 
Here is an illustration of a simulation that was started from the
restart file produced at \f$ t=0.175\f$ in the original run. The 
solution computed in the restarted run follows that obtained in 
the original simulation but it does not employ exactly the same timesteps.
  
\image html restart_trace.gif "Time-trace of the solution and the timestep chosen by the adaptive timestepper in the original (green) and restarted (red) simulation. " 
\image latex restart_trace.eps "Time-trace of the solution and the timestep chosen by the adaptive timestepper in the original (green) and restarted (red) simulation. " width=0.75\textwidth


\subsubsection ex Exercise:
Explain why the restarted simulation uses slightly different timesteps
and modify the \c dump_it(...) and \c restart(...) functions 
to solve this problem. [\b Hints: <B> (i) </B> Recall that the
adaptive timestepper
returns a suggestion for the next timestep. This is not recorded
in the restart data! <B> (ii) </B> If you can't solve the problem, 
have a look at the 
<A HREF="../../two_d_unsteady_heat_2adapt/html/index.html">
discussion of \c oomph-lib's doubly-adaptive unsteady Newton 
solver</A> where an improved dump/restart procedure is implemented.]





<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat_t_adapt/">
demo_drivers/unsteady_heat/two_d_unsteady_heat_t_adapt/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat_t_adapt/two_d_unsteady_heat_t_adapt.cc">
demo_drivers/unsteady_heat/two_d_unsteady_heat_t_adapt/two_d_unsteady_heat_t_adapt.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Demo problem: Large-amplitude non-axisymmetric oscillations of a thin-walled elastic ring

In this document we discuss the solution of a time-dependent
beam problem: The large-amplitude oscillations of a thin-walled
elastic ring. Specifically we shall
- demonstrate how to specify initial conditions for time-dependent
  simulations with \c oomph-lib's \c KirchhoffLoveBeam elements,
- demonstrate how to use the \c dump/restart functions for \c KirchhoffLoveBeam
  elements,
- show that timesteppers from the \c Newmark family
  - can be used with variable timesteps,
  - conserve energy,
  - can allocate and maintain storage for the solution at previous timesteps. 
.

<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Large-amplitude oscillations of a thin-walled elastic beam</B>
</CENTER> 

We wish to compute the large-amplitude oscillations of a
linearly-elastic, circular 
ring of undeformed radius \f$ R_0 \f$ and wall thickness \f$ h^* \f$ ,
subject to a transient pressure load
\f[
p_{ext}\left(\xi ,t \right) = 
\left\{ 
\begin{array}{ll}
-P_{cos}\cos\left(2\xi\right) & \mbox{for $t<T_{kick}$} \\
0  & \mbox{for $t>T_{kick}$} 
\end{array}
\right.
\mbox{\ \ \ \ \ \ \ \ \ \ (1) }
\f]
which initiates an oscillation in which the ring deforms into a
non-axisymmetric mode, as indicated in the sketch below:

\image html ring_sketch.gif "Sketch of the buckling ring. " 
\image latex ring_sketch.eps "Sketch of the buckling ring. " width=0.75\textwidth

We choose the same non-dimensionalisation as in 
<A HREF="../../steady_ring/html/index.html"> the previous steady 
example</A>  and parametrise the position vector 
to the ring's undeformed centreline as
\f[
{\bf r}_w(\xi) = {\bf R}_w(\xi,t=0) = \left( 
\begin{array}{c}
\cos(\xi) \\
\sin(\xi) 
\end{array}
\right),
\f]
where the non-dimensional arclength \f$ \xi \in [0,2\pi] \f$ 
along the ring's undeformed centreline acts as the Lagrangian
coordinate. Assuming that the ring is at rest for \f$ t\le 0 \f$ ,
we wish to compute the position vector to the deformed ring's 
centreline, \f$ {\bf R}_w(\xi,t) \f$ , for \f$ t > 0 \f$ .


</TD> 
</TR>
</TABLE>  
</CENTER>  

<HR>
<HR>

\section results Results


The figure below shows a snapshot, taken from 
<A HREF="../figures/unsteady_ring.avi">an animation of the ring's
computed deformation </A> for a wall thickness of \f$ h^*/R_0 =
1/20\f$ , a pressure load of \f$ p_{cos} = 1.0 \times 10^{-4}\f$ , and
\f$ T_{kick} = 15 \f$ . The arrows represent the instantaneous
velocity of the ring and show that at this point in time
the ring is still collapsing inwards. 

\image html unsteady_ring.gif "Shape and velocity of the ring. " 
\image latex unsteady_ring.eps "Shape and velocity of the ring. " width=0.75\textwidth


The red line with square markers in the figure below shows the 
time-history of the ring's control displacement: At \f$ t=0 \f$
the ring is in its initial undeformed configuration and we have 
\f$ R_{ctrl} = 1 \f$ . Application of the cosinusoidal pressure 
load during the interval \f$ 0 < t < T_{kick} = 15 \f$ causes 
the ring to deform non-axisymmetrically in the mode shape shown 
in the plot above. Because of inertia, the ring continues to 
collapse inwards even after the pressure load has been 
"switched off". The ring reaches its most strongly collapsed configuration,
in which the radius of the control point is reduced to just 
above 20% of its original value, at \f$ t \approx 50 \f$ .
Subsequently, the elastic restoring forces cause the ring to "reopen"
to its axisymmetric configuration (where \f$ R_{ctrl}=1 \f$ again),
which it traverses with finite velocity at \f$ t \approx 100 \f$ .
The ring overshoots the axisymmetric state and deforms in the
"opposite" direction to the deformation during the initial stages of collapse. 
It reverses its motion again when it reaches a second non-axisymmetric 
extreme at \f$ t \approx 150 \f$ . This is seen most clearly in the 
<A HREF="../figures/unsteady_ring.avi">animation of the ring's motion. </A>
 


\image html trace_file.gif "Plot of the control radius, the kinetic and potential (=strain) energy and their sum. " 
\image latex trace_file.eps "Plot of the control radius, the kinetic and potential (=strain) energy and their sum. " width=0.75\textwidth


The remaining lines in the plot show the ring's kinetic and potential (i.e. the
strain) energy and their sum. Up to \f$ t = T_{kick} \f$ , the
external load does work on the ring and increases its kinetic and
potential (strain) energy. Once the external load is "switched off"
the total energy stored in the system should (and indeed does) remain 
constant. See the section \ref energies for further details.


<HR>
<HR>

\section namespace Global parameters and functions
As usual, we employ a namespace to define the problem's physical
parameters and the load that acts on the ring. For the pressure loading
defined in equation (1), the load is given by
\f[
\mathbf{f} = -P_{cos}\cos\left(2\xi\right)\mathbf{N},
\f]
where \f$ \mathbf{N} \f$ is the outer unit normal to the ring's
deformed centreline. 


\dontinclude unsteady_ring.cc
\skipline start_of_namespace
\until end of load

We also define the non-dimensional wall thickness \f$ h \f$ and the timescale
ratio \f$ \Lambda^2 \f$. These are multiplied by powers 
of a scaling factor whose role will become apparent in the \ref 
exercises. (By default, the scaling factor is set to 1.0 and does 
not play any role.)

\until end of namespace

<HR>
<HR>

\section main The driver code
The main function is very simple. We store the (up to) two optional 
command line arguments which (if present) specify (i) a flag that indicates 
if the code is run in validation mode, and (ii) the name of a 
restart file.
 
\skipline start_of_main
\until CommandLineArgs::

Next, we build the problem with thirteen \c HermiteBeamElements 
and a \c Newmark<3> timestepper (
<A HREF="../../../linear_wave/two_d_linear_wave/html/index.html#IC">recall </A>
that a \c Newmark<NSTEPS> timestepper allocates and manages storage
for the solution at \c NSTEPS previous timesteps; we shall illustrate
this capability in the section \ref prev ), before executing
the timestepping loop.


\until end of main

<HR>
<HR>

\section problem The problem class 
The problem class is very similar to that used in the 
<A HREF="../../steady_ring/html/index.html"> previous, steady example</A> ,
but includes a few (obvious) additional functions that
specify the initial conditions (\ref IC) and perform the timestepping 
(\ref unsteady_run). We also provide two functions that allow us
to dump the solution to disk (\ref dump) and to restart the
time-dependent simulation (\ref restart).

\dontinclude unsteady_ring.cc
\skipline start_of_problem_class
\until restart(

The private member data includes an output stream that we shall use
to write a trace file. The two boolean flags indicate if the
code is run in  validation mode, and if the simulation has been
restarted, respectively. 

\until end of problem class

<HR>
<HR>

\section constructor The problem constructor
The constructor assigns default values for the two control flags
corresponding to a non-validation run without restart. We create
a timestepper of the type specified by the template parameter and
add it to the \c Problem's collection of timesteppers. 

\skipline start_of_constructor
\until add_time_stepper_pt

Next we build the geometric object that defines the shape of the
ring's undeformed centreline (an ellipse with unit half axes, i.e.
a unit circle) and use it to build the mesh.
As in <A HREF="../../steady_ring/html/index.html"> the previous steady 
example</A> we exploit the symmetry of the deformation and only
discretise a quarter of the domain.

\until ::time_stepper_pt()

The boundary conditions are identical to those in 
<A HREF="../../steady_ring/html/index.html"> the steady 
example</A>. 


\until pin_position(1,1); 

Finally, we pass the pointers to the load function and the pointer to the
geometric object that specifies the ring's initial shape to the elements
and assign the equation numbers. 


\until end of constructor

<HR>
<HR>

\section doc Post-processing
We compute the total kinetic and potential (=strain) energies stored in the 
(quarter-)ring and document them, their sum, and the control radius 
in the trace file.

\skipline start_of_doc_solution
\until // end of output to trace file

Next we use the default output function to document the 
ring shape and add a few tecplot commands to facilitate the
animation of the results.

\until some_file.close();

\subsection prev How to retrieve the solution at previous timesteps

The next few lines illustrate how to retrieve (and document) the ring shape at
previous timesteps. 
<A HREF="../../../linear_wave/two_d_linear_wave/html/index.html#IC">Recall
</A> that \c Newmark timesteppers are implicit, single-step
timestepping schemes that compute approximations for the 
time-derivatives, based on the solution at the current time level, 
and on "history values" that represent quantities at a single previous 
timestep. In some applications (particularly in fluid-structure 
interaction problems) it is necessary to keep track of the solution
at additional previous timesteps. Storage for such additional history values
is allocated (and managed) by the generalised \c Newmark<NSTEPS>
timesteppers if \c NSTEPS > 1. We stress that these
additional history values are not involved in the approximation 
of the time-derivatives; they are simply stored and updated by the
timestepper when the solution is advanced to the next timestep.

  
<A HREF="../../../unsteady_heat/two_d_unsteady_heat/html/index.html#IC">Recall
</A> also that the member function \c TimeStepper::nprev_values() may
be used to determine how many of the history values that are stored in an 
associated \c Data object represent the solution at previous
timesteps. Finally, 
<A HREF="../../../unsteady_heat/two_d_unsteady_heat/html/index.html#IC">recall
</A>  that history values that represent the solution at previous
timesteps are always stored before those that represent "generalised" history
values (such as approximations of the first time-derivative at the
previous timestep, etc).  It is therefore always possible to determine
how many previous solutions are stored in a \c Data object, and where
they are stored.

To document the shape of a \c HermiteBeamElement at a previous timestep,
the \c HermiteBeamElement provides an additional three-argument 
output function that may be called as follows:

\until end of output of previous solutions

At timestep \c t, these statements create the output files
\c RESLT/ring t\c -0\c .dat,   \c RESLT/ring t\c -1\c .dat,   \c
RESLT/ring t\c -2\c .dat,  \c RESLT/ring t\c -3\c .dat,  which contain the
shape of the ring at the \c t -th, ( \c t -1)th,  ( \c t -2)th and  
( \c t -3)th, timestep, respectively. 

Finally, we write a restart file that will allow us to restart the
simulation. 

\until end of doc solution

<HR>
<HR>

\section dump Writing a restart file
Writing the restart file for the present problem is as easy as in
the previous examples, as the generic \c Problem data may again be 
written with the \c Problem::dump(...) function. We customise
the restart file slightly by adding the value \f$ p_{cos} \f$
and the flag that indicates if the code is run in validation mode.

\skipline start_of_dump_it
\until end of dump it

<HR>
<HR>

\section restart Restarting from a file
The restart operation reverses the steps performed in the dump
function: We recover the two problem-specific parameters and then
read the generic \c Problem data with the 
\c Problem::read(...) function.

\skipline start_of_restart
\until end of restart

<HR>
<HR>

\section IC Setting the initial condition
The assignment of initial conditions depends on whether or not
a restart from a previous computation is performed. If no restart is
performed, we specify the initial timestep, \c dt,
and assign history values that are consistent with an impulsive
start from the ring's initial shape.

\skipline start_of_set_ic 
\until }

If the computation is restarted, the name of the restart file
will have been specified on the command line. We try to open 
the restart file

\until }

and display an error message and terminate the program execution
if the file cannot be opened.

\until }

If the file can be opened we call the restart function which returns the
the \c Problem into the state it was in when the restart file 
was written. No further steps are required. 

\until end of set ic

<HR>
<HR>

\section unsteady_run The timestepping loop
We start by converting the (optional) command line arguments
into the flags that determine what mode the code is run in:
Without command line arguments, we use the default assignments,
as specified in \ref constructor.

\skipline start_of_unsteady_run
\until }

A single command line argument is interpreted as the "validation run"
flag (1 for true, 0 for false) which will be used to limit the number of
timesteps. 

\until }

The presence of two command line arguments indicates 
that a restart is performed. In this case the second argument 
specifies the name of the restart file.

\until }

We print an error message if the code is run with any other number of
command line arguments.

\until }

We create a \c DocInfo object to specify the name of the output
directory, and open the trace file.

\until << std::endl;

Next, we set the problem parameters and the number of timesteps
to be performed, before assigning the initial conditions.

\until set_initial

The timestepping loop itself is practically identical to
that used in 
<A HREF="../../../unsteady_heat/two_d_unsteady_heat/html/index.html#main"> 
driver codes for other unsteady problems</A>. To demonstrate
that \c Newmark timesteppers can deal with variable timesteps,
we reduce the timestep slightly after every step. 

\until end of unsteady run

<HR>
<HR> 

\section comments Comments and Exercises
\subsection timescale The default non-dimensionalisation of time
The non-dimensionalisation of the principle of virtual displacements
that forms the basis of \c oomph-lib's beam elements, was 
discussed in detail in <A HREF="../../tensioned_string/html/index.html">
an earlier example</A>. However, since this is the first 
time-dependent beam problem it is worth re-iterating that, by default, time is
non-dimensionalised on the timescale for extensional deformations
i.e. on the natural timescale 
\f[ 
{\cal T}_{natural} = {\cal L} \sqrt{\frac{\rho}{E_{eff}}}
\f] 
of oscillations in which the beam 
is stretched/compressed along its centreline. The relation between
the dimensional time \f$ t^* \f$ and its non-dimensional equivalent \f$ t
\f$ is given by
\f[
t^* = \underbrace{{\cal L} 
\sqrt{\frac{\rho}{E_{eff}}}}_{{\cal T}_{natural}} \ t,
\f]
where \f$ {\cal L} \f$ is the lengthscale used to non-dimensionalise all
lengths (in the present example \f$ {\cal L} = R_0 \f$ , 
the undeformed radius of the ring), 
\f$ \rho \f$ is the density of the ring, and \f$ E_{eff} = E/(1-\nu^2)\f$ 
is the "effective" 1D Young's modulus of the beam, formed with
its 3D Young's modulus \f$ E \f$ , and its Poisson ratio \f$ \nu \f$ .

This non-dimensionalisation of time is consistent with the 
non-dimensionalisation  of all stresses/loads on  \f$ E_{eff}\f$. 
It implies that if the beam deforms in a mode in which its
deformation is dominated by bending effects, the numerical 
values for the non-dimensional load are relatively 
small (indicating that the loads required to induce a deformation 
of a given size are much smaller if the ring deforms in a
bending-dominated mode, than in a mode in which it is dominated
by extensional deformations), while the non-dimensional 
period of the oscillation is  relatively large (indicating that bending 
oscillations occur at a much smaller frequency than oscillations in
which the beam's deformation is dominated by extensional deformations).  

<HR>


\subsection energies The default non-dimensionalisation for the kinetic and potential (strain) energies

With the default non-dimensionalisation discussed above, the non-dimensional
kinetic and potential (strain) energies are given by
\f[
\Pi_{strain} = \frac{\Pi_{strain}^*}{{\cal L} \, h^* E_{eff} } =
\frac{1}{2}\int \left( \gamma^2 + 
\frac{1}{12}\left(\frac{h^*}{{\cal L}}\right)^2
\kappa^2 \right) \ d \xi
\f]
and 
\f[
\Pi_{kin} = \frac{\Pi_{kin}^*}{{\cal L} \, h^* E_{eff}}
= \frac{1}{2} 
\int 
\frac{\partial {\bf R}_w}{\partial t} \cdot
\frac{\partial {\bf R}_w}{\partial t} \ d \xi,
\f]
respectively. Conservation of energy implies that
\f[
\Pi_{kin} + \Pi_{strain} = const.
\f]
if there is no external forcing. The plot of the energies
shown at the beginning of this document shows that the
time-integration with the Newmark method is energy-conserving.

<HR>
\subsection change_T Changing the timescale used to non-dimensionalise the equations

It is possible to non-dimensionalise the governing equations on 
a different timescale, \f$ {\cal T} \f$ , so that
\f[
t^* = {\cal T} \ t.
\f]
This is achieved by overwriting the default assignment 
for the ratio
\f[
\Lambda = \frac{{\cal T}_{natural}}{{\cal T}} = 
\frac{{\cal L}}{{\cal T}} \sqrt{\frac{\rho}{E_{eff}}},
\f]
of the natural timescale \f$ {\cal T}_{natural}\f$ and
the time \f$ {\cal T} \f$ used to non-dimensionalise the equations.
 
By default, we have \f$ \Lambda = 1 \f$ but the member function
\code
KirchhoffLoveBeamEquations::lambda_sq_pt()
\endcode
may be used to assign a different value for the \e square of the
timescale ratio which may also be interpreted as the
non-dimensional density. The case \f$ \Lambda^2=0 \f$  therefore
corresponds to the case of zero wall inertia. 
(We store \f$ \Lambda^2 \f$ rather than \f$ \Lambda
\f$ itself because the governing equations contain only the square of the
timescale ratio).  As with most other physical parameters,
\f$ \Lambda^2 \f$ must be defined as a global variable, preferably
by adding it to the namespace that contains the problem parameters, 
e.g.

\code 
//====start_of_namespace============================
/// Namespace for global parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Square of timescale ratio, i.e. the non-dimensional density
 double Lambda_sq=4.0;

 [...]

} // end of namespace
\endcode

...and passing a pointer to \f$ \Lambda^2 \f$ to the elements.
The statement can be added to the loop over the elements that
passes the pointer to the pressure load to the elements.

\code

 //Loop over the elements to set physical parameters etc.
 for(unsigned i=0;i<n_element;i++)
  {
   // Cast to proper element type
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   // Pass pointer to square of timescale ratio (non-dimensional density)
   elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   [...]

  }

\endcode

The definition of the non-dimensional kinetic energy
used in \c KirchhoffLoveBeamEquations::get_energy(...) incorporates
the timescale ratio by computing the non-dimensional kinetic energy as
\f[
\Pi_{kin} = 
\frac{\Pi_{kin}^*}{{\cal L} \, h^* E_{eff}} = 
\frac{1}{2} \Lambda^2
\int 
\frac{\partial {\bf R}_w}{\partial t} \cdot
\frac{\partial {\bf R}_w}{\partial t} \ d \xi
\f]
which, for the default assignment \f$ \Lambda^2 = 1 \f$ , reduces
to the definition given above.
<HR>

\subsection exercises Exercises
-# As discussed in \ref timescale , the period of the oscillation
   ( \f$ \approx 170 \f$ ) computed in the above example has a 
   large numerical value because
   time was non-dimensionalised on a timescale that is representative
   for oscillations in which the ring's motion is dominated by
   extensional (rather than bending) deformations. Furthermore,
   because the load is non-dimensionalised on the ring's extensional
   (rather than its bending) stiffness, a very small load 
   \f$ p_{cos} = 10^{-4} \f$ was sufficient to induce large-amplitude
   oscillations.
   \n \n
   Change the load function to a spatially-constant external pressure
   (which will deform the ring axisymmetrically so that its
   motion is dominated by extensional deformations) to confirm 
   that in this mode 
   -# \f$ {\cal O}(1) \f$ pressures are required to induce
      \f$ {\cal O}(1) \f$ axisymmetric displacements, and
   -# the period of the axisymmetric oscillations is
      \f$ {\cal O}(1) \f$. 
   . 
   [\b Hint: You will have to reduce the timestep to capture the
   much faster axisymmetric oscillations.] 
-# Use the principle of virtual displacements (without prestress)
   \f[
   \int_0^{L}   \left[ 
   \gamma \ \delta  \gamma +  
   \frac{1}{12} h^2 \kappa 
   \ \delta \kappa   - 
   \left(\frac{1}{h} \sqrt{\frac{A}{a}} \   {\bf f} - \Lambda^2 
   \frac{\partial^2 {\bf R}_w}{\partial t^2} \right) \cdot  
   \delta {\bf R}_w 
   \right] \ \sqrt{a} \ d\xi  = 0,
   \ \ \ \ \ \ \ \ \  (2)
   \f]
   to show (analytically) that for bending-dominated deformations (i.e. 
   deformations for which \f$ |\gamma| \ll |\kappa| \f$)
   -# an increase in the wall thickness \f$ h \f$ by a factor 2, say, will 
      -# reduce the amplitude of the ring's static deformation 
         by a factor of 8 (thicker rings are stiffer)
      -# reduces the period of the unforced oscillations 
         (i.e. the oscillations it performs when \f$ {\bf f} = 
         {\bf 0} \f$ ) by factor of  4 (thicker 
         rings oscillate more rapidly). 
      .
   -# an increase in the ring's density \f$ \rho \sim \Lambda^2\f$ 
      by a factor of 2, say, increases the period of its
      bending-dominated oscillations by a factor of 4 (heavier rings 
      oscillate more slowly).
   .
   Combine these results to show that
   a ring of wall thickness \f$ h \f$ and a density \f$ \rho \f$, 
   subject to a forcing of magnitude \f$ f \f$ will deform 
   (approximately) as a ring of wall thickness \f$ \alpha h \f$ with 
   a density  \f$ \alpha^2 \rho \f$, subject to a forcing of magnitude
   \f$ \alpha^3 f\f$.
   Confirm these theoretical predictions computationally: Change
   the initial assignment for the scaling factor 
   \c Global_Physical_Variables::Alpha and repeat the computation. 
-# Confirm that the restart procedure works correctly by plotting
   the time-trace obtained from a restarted simulation on top
   of the original time-trace. 
-# Compare the results obtained from a simulation with a variable
   timestep against the results obtained from a computation with a 
   fixed timestep (comment out the line that reduces the timestep
   after every solve). 
-# Compare the various output files generated in \c doc_solution()
   to confirm that the \c Newmark<3> timestepper correctly maintains
   the history of the solution at three previous timesteps. 
   E.g. confirm that the output file \c RESLT/ring3-0.dat which contains the
   ring shape at timestep 3 is identical to \c RESLT/ring5-2.dat which
   contains the ring shape computed two timesteps before timestep 5. 
   Explore \c oomph-lib's internal use of the history values by
   analysing the functions 
   \code 
   HermiteBeamElement::output(...)
   \endcode
   \code
   FiniteElement::interpolated_x(...)
   \endcode
   and
   \code
   FiniteElement::x_gen(...)
   \endcode
   


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/beam/unsteady_ring/">
demo_drivers/beam/unsteady_ring/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/beam/unsteady_ring/unsteady_ring.cc">
demo_drivers/beam/unsteady_ring/unsteady_ring.cc
</A>
</CENTER>
.






<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


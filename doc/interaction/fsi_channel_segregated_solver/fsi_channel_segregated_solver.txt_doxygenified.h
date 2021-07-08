/**

\mainpage Segregated solvers for fluid-structure-interaction problems: Revisiting the flow in a 2D collapsible channel
 
In this document we discuss the implementation of segregated solution
strategies for multi-physics problems, in particular fluid-structure
interaction, within \c oomph-lib. The method illustrated by
revisiting the fluid-structure interaction problem of
<A HREF="../../fsi_collapsible_channel/html/index.html">
 finite-Reynolds-number flow in a 2D collapsible channel;</A> an
 example discussed in detail in <CENTER>
<a href="http://www.springerlink.com/content/m3r6318701g338g4/">Heil, M., Hazel, A.L. & Boyle, J. (2008): Solvers for large-displacement fluid-structure interaction problems: Segregated vs. monolithic approaches. Computational Mechanics.</a>
</CENTER>

where we compare the relative performance of segregated and monolithic
solvers. Since the paper comes to the conclusion that, despite
various claims in the literature, segregated solvers are not
necessarily more efficient than fully-coupled monolithic schemes
(of the type employed in \c oomph-lib) you should also consult
the <a href="../../../preconditioners/fsi/html/index.html">related
tutorial on the monolithic solution of the problem
with \c oomph-lib's FSI preconditioner.</A>
 
<HR>
<HR> 


\section problem The problem 


<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>Flow in a 2D collapsible channel </B>
</CENTER> 

The figure below shows a sketch of the problem: 
Flow is driven by a 
prescribed Poiseuille flow \f$U^{*}_{p}\f$ through a 2D channel of width 
\f$ H^* \f$ and total length
\f$ L^*_{total} = L^*_{up} + L^*_{collapsible} + L^*_{down}. \f$
The upstream and downstream lengths of the channel are rigid, whereas
the upper wall in the central section is an elastic membrane whose
shape is parametrised by a Lagrangian coordinate, \f$ \xi^* \f$ , so 
that the position vector to the moving wall is given by 
\f$ {\bf R}_w^*(\xi^*,t^*) \f$ . The wall is loaded by the
external pressure \f$ p_{ext}^* \f$ and by the traction 
that the viscous fluid exerts on it.

\image html collapsible_channel_sketch.gif "Sketch of the problem. " 
\image latex collapsible_channel_sketch.eps "Sketch of the problem. " width=0.75\textwidth

The non-dimensionalisation and governing equations have already been
discussed in the 
<A HREF="../../fsi_collapsible_channel/html/index.html"> previous
(monolithic) example. </A>
The problem is not quite the same, however, because the upstream
boundary condition is now one of prescribed flow, rather than prescribed
pressure:
- Prescribed inflow, 
  \f[
  {\bf u}(x_1,x_2) = {\bf u}_{p}(x_1,x_2) =   
  6 \ x_2 \  (1-x_2) \ {\bf e}_1.
  \ \ \ \ \ \ \ \ \ \ (1)
  \f] at \f$ x_{1}=0\f$.
.
All other boundary conditions remain the same.

</TD>
</TR>
</TABLE>  
</CENTER>
 
<HR>
<HR>



 
\section reslt Results
The behaviour of the system under the prescribed-inflow boundary
conditions is somewhat different to its behaviour when the pressure
drop is prescribed. In the first instance, we consider steady states,
in which all time-derivatives are neglected. The figure below shows 
steady flows at a Reynolds number of \f$ Re=500 \f$ and two values
of the fluid-structure-interaction parameter, \f$ Q =
10^{-4} \f$ (upper) and \f$ Q = 10^{-2}\f$ (lower). For low values of
\f$Q\f$, corresponding to weak fluid-structure interaction, the
deformation of the wall is approximately symmetric, being dominated by
the external pressure. As \f$Q\f$
increases, the influence of fluid traction can been seen in the
asymmetric deformation of the elastic wall. The viscous pressure
drop along the tube leads to higher pressure upstream (causing an
outward deformation) and lower pressures downstream (causing an inward
deflection).

\image html steady_flows.gif "Steady flows at Re=500 and Q=10e-4 (upper), Q=10e-2 (lower). " 
\image latex steady_flows.eps "Steady flows at Re=500 and Q=10e-4 (upper), Q=10e-2 (lower). " width=0.8\textwidth

The overall behaviour of the system can be characterised by steady
load-displacement curves in which the vertical position of a 
control point on the elastic
section of the channel wall is plotted as a function of the external
pressure. 

\image html steady_trace.gif "Load-displacement curve: the vertical position of a control point on the elastic wall (located at 50, 50, 60 and 70 percent of its length for Q = 0, 10e-4, 10e-3 and 10e-2, respectively) as a function of the external pressure. " 
\image latex steady_trace.eps "Load-displacement curve: the vertical position of a control point on the elastic wall (located at 50, 50, 60 and 70 percent of its length for Q = 0, 10e-4, 10e-3 and 10e-2, respectively) as a function of the external pressure. " width=0.5\textwidth

 At low \f$Q\f$, the displacement is directly proportional to
the external pressure. As \f$Q\f$ increases the curves shift to the
right because a large external pressure is required to keep the wall
in its undeformed position;  
a consequence of the increased viscous pressure drop \em and the 
boundary condition that \f$p=0\f$ at the outlet. A second consequence
of increasing \f$Q\f$ is that (at finite Reynolds number) a 
smaller increase in external pressure is required to achieve 
a given degree of collapse. This is because the Bernoulli
effect reduces the fluid pressure in the region that is
most strongly collapsed and therefore increases the compressive
load on the wall. For \f$Q = 10^{-2}\f$ two limit
points develop on the load-displacement curve, indicating that the wall
"snaps through" into a collapsed buckled configuration when
\f$p_{ext}\f$ becomes sufficiently large. The appearance of the 
limit points means that it is no longer possible to perform the
steady parameter study by slowly increasing \f$p_{ext}: \f$ 
At sufficiently large values of \f$Q\f$ the displacement of the control 
point is not a single-valued function of the external pressure 
\f$ p_{ext}\f$. However, the application of "displacement control",
described in the tutorial discussing the
<A HREF="../../../beam/steady_ring/html/index.html"> 
large-displacement post-buckling of an elastic ring</A> is sufficient
to circumvent this difficulty: We treat the external pressure as an unknown
and control the channel's collapse by prescribing the vertical position
of the control point, \f$ x_2^{[ctrl]}.\f$ This resolves the problem
because the curve \f$ p_{ext}(x_2^{[ctrl]})\f$ is single-valued,
allowing us to perform the parameter study by slowly increasing
the wall collapse by reducing \f$ x_2^{[ctrl]}\f$, computing the
pressure required to achieve this deformation as part of the solution.

<HR>
<HR>

\section overview Overview: Segregated solution strategies with oomph-lib

The general methodology for setting up fluid-structure-interaction
problems is discussed in 
<A HREF="../../fsi_collapsible_channel/html/index.html#overview">another 
tutorial; </A>
and we shall assume that the standard monolithic problem has already
been written. In the present example, the monolithic problem class \c
FSICollapsibleChannelProblem is
specified in the header file 
<A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/fsi_chan_problem.h">
fsi_chan_problem.h</A> 

Having specified the monolithic (fully-coupled) discretisation, our segregated
solution strategy proceeds by alternating between fluid and solid
solves: Initially, the degrees of freedom
associated with the (pure) solid mechanics problem are "pinned" and the
global assembly procedure is modified to omit the corresponding 
solid elements. The Newton solver will, therefore, solve the fluid
equations with a "frozen" wall shape. Next, the degrees of freedom
associated with the (pure) fluid  mechanics problem are pinned and the
original boundary conditions for the solid mechanics problem are
re-assigned.  The assembly procedure
is now modified so that only solid elements contribute to the 
global system. The Newton solver will then solve for a new wall shape
corresponding to the tractions exerted by the given flow
field. At this point we allow for under-relaxation, i.e. we provide
the option to increment the solid mechanics degrees of freedom
by a fraction of the change computed by the Newton solver.
These two steps are repeated in a fixed-point iteration which
continues until a given convergence criterion is satisfied, or a
maximum number of iterations is exceeded. 
We note that different linear solvers/preconditioners 
may be specified for solution of the linear systems 
arising during the Newton iteration for the isolated "fluid" and 
"solid" problems, allowing the re-use of optimal solution methods 
for individual sub-problems. This is generally perceived to be one
of the key advantages of segregated solvers.

<HR>
<HR>

\section how_its_done Brief discussion of the implementation

\subsection seg_problem The SegregatableFSIProblem

<b>a. Overall structure</b>

The 
<A HREF="../../../the_data_structure/html/classoomph_1_1SegregatableFSIProblem.html">
 \c SegregatableFSIProblem </A> class is used to implement our
segregated solution strategy within \c oomph-lib. The most important
problem-specific task is to divide all the problem data 
into distinct fluid and solid degrees of freedom
and to partition the monolithic mesh into a mesh of fluid
elements and a mesh of solid elements. The problem-specific
partitioning should be implemented in the (pure) virtual member
function
\code 
 /// Identify fluid and solid data
 virtual void identify_fluid_and_solid_dofs(Vector<Data*>& fluid_data_pt,
                                            Vector<Data*>& solid_data_pt,
                                            Mesh*& fluid_mesh_pt,
                                            Mesh*& solid_mesh_pt)=0;
\endcode 
 which returns vectors of fluid and solid data and the meshes of fluid
 and solid elements. This virtual function is called within the 
 member function 
\code
  /// \short Set up segregated solver. The optional boolean argument
  /// defaults to true and causes the identify_fluid_and_solid_dofs(...)
  /// to be called again. This is required, e.g. if any of the
  /// meshes were adapted since the previous call to the segregated solver.
  void setup_segregated_solver(const bool &full_setup_of_fluid_and_solid_dofs=true)
\endcode
 which \em must be
 called immediately before every segregated solve. The optional boolean
 flag may be set to \c false if the solid and fluid meshes 
 have not changed between solves (i.e. if no spatial adaptation was
 performed since the last call to the segregated solver). 
 The \c setup_segregated_solver(...)
 function must still be called, however, in order that 
 data associated with convergence acceleration techniques is reset to
 its default values.

<b>b. The segregated solvers</b>

 The class inherits from the standard \c Problem
 class, which provides the standard (monolithic) \c newton_solve() and
 related functions. Thus any \c SegregatedFSIProblem can be solved
 "monolithically" as normal and, moreover, it is straightforward to ensure that
 exactly the same system is being solved when comparing monolithic and
 segregated solutions. The segregated solution strategy is implemented in
 the analogous member functions:
 - The equivalent of the monolithic \c Problem::newton_solve() is
   \code 
   SegregatedFSIProblem::segregated_solve();
   \endcode
 - The equivalent of the monolithic \c Problem::steady_newton_solve() is
   \code
   SegregatedFSIProblem::steady_segregated_solve();
   \endcode
 - Finally, the equivalent of  \c Problem::unsteady_newton_solve(dt) is
   \code
   PicardConvergenceData unsteady_segregated_solve(const double &dt);
   \endcode
 .
All three functions return an instance of a \c PicardConvergenceData
object which stores the convergence statistics of the segregated solve.

 In addition, the virtual member functions
 \code
 SegregatedFSIProblem::actions_before_segregated_solve() 
 \endcode
 \code
 SegregatedFSIProblem::actions_after_segregated_solve() 
 \endcode
 \code SegregatedFSIProblem::actions_before_segregated_convergence_check() 
 \endcode
 are provided to allow the user to specify any actions, such as 
initialisation of
 counters, mesh updates, output, etc, that should be performed before or after
 each complete segregated solve. Note that the \c Problem member functions
 \code
 Problem::actions_before_newton_solve() 
\endcode
\code
 Problem::actions_after_newton_solve() 
\endcode
\code
 Problem::actions_before_newton_convergence_check() 
 \endcode
 are called as usual during the Newton solve of each
 sub-problem and may be used for fine-grained operations that should
 be performed before or after each fluid or solid solve. 
 For this purpose, the \c SegregatedFSIProblem provides a flag, 
\c int \c SegregatedFSIProblem::Solve_type that indicates which
(sub-)solve is currently being performed. The flag can take
the (enumerated) values \c SegregatedFSIProblem::Full_solve, 
\c SegregatedFSIProblem::Fluid_solve and 
\c SegregatedFSIProblem::Solid_solve, allowing the user to perform
specific actions during the distinct sub-solves.

<b>c. Choosing the convergence criterion</b>


 Other public member functions provided
 by the \c SegregatedFSIProblem class are used to specify the
 convergence criterion for the global fixed-point iteration: 
 \code
 /// Base convergence based on max. global residual
 void assess_convergence_based_on_max_global_residual(const double &tol)
 \endcode
 \code
 /// Base convergence on maximum absolute change of solid degrees of freedom
 void assess_convergence_based_on_absolute_solid_change(const double &tol)
 \endcode
 \code
 /// Base convergence on maximum relative change of solid degrees of freedom
 void assess_convergence_based_on_relative_solid_change(const double &tol)
 \endcode
 If a tolerance is not specified the default \c
 Problem::Newton_solver_tolerance is used.

<b>d. Under-relaxation</b>

 Finally, there are several member functions that are used to specify
 the convergence-acceleration techniques:
 -# <b>Static under-relaxation:</b>
    \n\n
    \code
    //Use under-relaxation for solid degrees of freedom and specify 
    //the optional under-relaxation parameter. The default of 1.0 
    //corresponds to no under-relaxation.
    void use_under_relaxation (const double &omega=1.0)
    \endcode
    If this function is called, under-relaxation is performed after the
    solid sub-solve, i.e. each solid degree of freedom, \f$ s \f$, say
    is updated via
    \f[
    s = s_{new} + (1-\omega) (s_{old} - s_{new})
    \f]
    where \f$ s_{new} \f$ is the new value computed by the Newton solver
    and \f$ s_{old} \f$ is its previous value.
 -# <b>Adaptive under-relaxation:</b>
    \n\n
    \code
    //Boolean flag that controls whether Irons & Tuck extrapolation
    //is used to dynamically modify the under-relaxation parameter for
    //the under-relaxation of the solid degrees of freedom.
    void enable_irons_and_tuck_extrapolation ()
    \endcode
    If this function is called (and if under-relaxation is enabled)
    the under-relaxation parameter \f$ \omega \f$ is adjusted
    throughout the fixed-point iteration, using Irons & Tucks
    convergence acceleration procedure; see Irons, B.M. & Tuck, R.C. 
    "A version of the Aitken accelerator for computer iteration".
    International Journal of Numerical Methods in Engineering <b>1</b>,
    275-277 (1969).
 -# <b>Pointwise Aitken-acceleration:</b>
    \n\n
    \code
    //Set a boolean flag that controls whether pointwise Aitken
    //extrapolation is used. The optional argument specifies the Picard
    //Iteration after which the extrapolation is to be used for the first
    //time. The default value is zero.
    void enable_pointwise_aitken (const unsigned &pointwise_aitken_start)
    \endcode
    If this function is called,  the classical Aitken extrapolation is
    used to accelerate the convergence of (individual) solid degrees of freedom
    after every three iterations.
 .
<hr>
<hr>

\section channel_seg_problem The SegregatedFSICollapsibleChannelProblem

We shall now briefly discuss the application of the segregated
solver for the collapsible channel problem. 
The \c SegregatedFSICollapsibleChannelProblem is defined in the driver
code <code><A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc">simple_segregated_driver.cc</a> 
</code> and inherits from the 
"monolithic" \c FSICollapsibleChannelProblem and also from the 
<A HREF="../../../the_data_structure/html/classoomph_1_1SegregatableFSIProblem.html">
 \c SegregatableFSIProblem </A> class. The code <code>
 <A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc">simple_segregated_driver.cc</a> </code>
is specifically designed for ease of
 exposition and does not contain any timing statements or
 documentation of convergence histories. The alternative driver code 
 <code>
 <A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/fsi_chan_seg_driver.cc">fsi_chan_seg_driver.cc</a></code> contains complete
 timing and documentation statements and is the code that was used by
 <a href="http://www.springerlink.com/content/m3r6318701g338g4/"> Heil,
 Hazel \& Boyle (2008). </a>


 The simplified \c SegregatedFSICollapsibleChannelProblem  class 
contains six member functions
- The constructor
- The destructor
- \c void \c identify_fluid_and_solid_dofs(...)
- \c void \c actions_before_newton_convergence_check()
- \c void \c actions_before_segregated_convergence_check()
- \c void \c steady_run()
- \c void \c doc_solution(DocInfo& doc_info)
.
The \c doc_solution(...) function simply writes the bulk (fluid)
elements and wall (solid) elements to two separate files and the 
destructor is empty. We discuss the other four member functions below.

<hr>


\subsection constructor The constructor
The constructor calls the constructor of the underlying "monolithic"
problem and then selects the convergence criterion and
convergence-acceleration technique based on the values of control
flags defined in the namespace \c Flags.

\dontinclude simple_segregated_driver.cc
\skipline start_of_constructor
\until end_of

<hr>

\subsection identify Identifying the fluid and solid degrees of freedom

The underlying monolithic problem provides pointers to the 
fluid and solid (sub-)meshes via the member data
\code
 AlgebraicCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;
\endcode
\code
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 
\endcode
which are accessible via the member functions  
\c SegregatedFSICollapsibleChannelProblem::bulk_mesh_pt()
and \c SegregatedFSICollapsibleChannelProblem::wall_mesh_pt(),
and so the identification of fluid and solid degrees of freedom is
reasonably straightforward. The only complication arises because we
may, or may not, be using displacement control which introduces
a further element into the global mesh. Displacement control
affects the solid problem suggesting that the (variable) external 
pressure should be regarded as a solid degrees of freedom and the \c
DisplacementControlElement should be included in the solid mesh.

\dontinclude simple_segregated_driver.cc
\skipline start_of_identify_fluid_and_solid
\until end_of_identify_fluid_and_solid

<hr>

\subsection actions Actions before convergence checks

During a monolithic solve the function 
\c actions_before_newton_convergence_check() must 
update the nodal positions in the bulk (fluid) mesh. 
In principle, it should remain empty during a
segregated solve, but we found it beneficial to update the bulk 
mesh, and hence the fluid load on the wall, during the 
solution of the solid problem.

The function \c actions_before_segregated_convergence_check() 
contains an update of the nodal positions in the bulk mesh 
in order that the segregated solution is self-consistent.

\dontinclude simple_segregated_driver.cc
\skipline Update nodal positions in the fluid mesh in
\until end_of_convergence_checks

<hr>

\subsection steady Solving a steady problem

The function \c steady_run() conducts a simple parameter study in
which the external pressure (or prescribed displacement) is 
varied. After specification of the initial conditions, parameter
increments and output directories, the parameter study is
straightforward
\dontinclude simple_segregated_driver.cc
\skipline Parameter study
\until End of parameter

<hr>
<hr>

\section main The driver code

 Having written our \c SegregatedFSICollapsibleChannelProblem, the
 driver code is extremely simple. We specify number of elements and
 dimensions of our computational domain, construct the problem and
 perform a steady parameter study.

\dontinclude simple_segregated_driver.cc
\skipline start_of_main
\until  end of main

<hr>
<hr> 


\section comments_and_ex Comments and Exercises

\subsection comments Comments


- <b>(In-)efficiency of \c setup_segregated_solver()</b> \n\n
  In our simple example code, we did not employ spatial adaptivity. It is
  not necessary, therefore, to (re-)identify the fluid and solid degrees of
  freedom before each solve, the default (safe) behaviour of \c
  setup_segregated_solver(). Nonetheless, data associated with
  the techniques used to accelerate the convergence of the Picard
  iterations must be reset before each segregated solve. In the 
  more complex driver code, a boolean
  flag \c bool \c full_setup is used as an argument to \c
  setup_segregated_solver() which modifies the behaviour, as
  indicated below.
  \code 
  // Boolean flag used to specify whether a full setup of solid and fluid dofs
  // is required
  bool full_setup = true;

  // Parameter study
  for (unsigned istep=0;istep<Flags::Nsteps;istep++)
   {
    // Setup segregated solver
    setup_segregated_solver(full_setup);
 
     [...]
   
    steady_segregated_solve()
  
     [...] 
  
    //We no longer need a full setup of the dofs
    full_setup = false;
   }
  \endcode  
.  

<hr>
<hr>

\subsection ex Exercises
-# Modify the control flags in <code><A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc">simple_segregated_driver.cc</a></code> to
   verify that the monolithic solution is the same (to within finite
   precision) as the segregated solution.
-# Modify the control flags in <code><A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc">simple_segregated_driver.cc</a></code> to
   investigate the influence of the convergence acceleration
   techniques and convergence criterion on the segregated solution. Which
   combination of parameters gives convergence in the fewest Picard
   iterations? 
-# Investigate the behaviour of the system if the fluid (bulk) mesh is \em not
   updated after each Newton step in the solution of the solid
   problem. Can you obtain converged solutions?
-# Write your own \c
   SegregatedFSICollapsibleChannelFlow::unsteady_run()
   member function that computes the time evolution of the system after
   a perturbation to the external pressure. Compare your answer with
   the equivalent member function in the much more comprehensive
   driver code
   <code><A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/fsi_chan_seg_driver.cc">fsi_chan_seg_driver.cc</a></code> that was used in <a href="http://www.springerlink.com/content/m3r6318701g338g4/"> Heil, Hazel \& Boyle (2008). </a>
.   

<hr>
<hr>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/">
demo_drivers/interaction/fsi_channel_seg_and_precond/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc">
demo_drivers/interaction/fsi_channel_seg_and_precond/simple_segregated_driver.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Example problem: Solution of the 2D unsteady heat equation.

This is our first time-dependent example problem. We will demonstrate
that, compared to the solution of steady problems, the solution of 
time-dependent problems only requires a few additional steps:
- The creation of a suitable \c TimeStepper object and its addition 
  to the \c Problem's list of timesteppers. (\c Problems may 
  employ multiple \c TimeSteppers -- a key requirement for the simulation
  of multiphysics problems.)
- The initialisation of the timestep, \c dt.
- Setting the initial conditions by assigning suitable values for
  the \c Data objects' "history values" and the \c Nodes' "positional
  history values".
- Optionally: (Re-)implementing the empty virtual functions
  \c Problem::actions_before_implicit_timestep() and 
  \c Problem::actions_after_implicit_timestep(), e.g. to update time-dependent
  boundary conditions before each timestep.
.

Once these steps have been performed, \c oomph-lib's unsteady Newton
solver, \c Problem::unsteady_newton_solve(...), may be called to advance 
the solution from its state at time \f$ t \f$ to its new state at
 \f$ t \ + \f$ \c dt.

<HR>
<HR>

\section example The example problem
We will illustrate the basic timestepping procedures by considering
the solution of the 2D unsteady heat equation in a square domain:

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
where the functions \f$ g_0 \f$ and  \f$ h_0 \f$ are given.
</TD>
</TR>
</TABLE>  
</CENTER>

Here we consider the unforced case, \f$ f=0 \f$, and choose boundary
and initial conditions that are consistent with the exact solution
\f[
u_0(x_1,x_2,t) = e^{-Kt}\sin\left( \sqrt{K} \left( x_1 \cos \Phi + x_2 \sin
\Phi\right)\right),
 \ \ \ \ \ \ \ \ \ \ (4)
\f]
where \f$ K \f$ and \f$ \Phi \f$ are constants, controlling the
decay rate of the solution and its spatial orientation, respectively.

The figure below shows a plot of computed and exact solutions at
time \f$ t = 0.01 \f$ (for an animation  
<A HREF="../figures/unsteady_heat_soln.avi"> click
here</A>). 

\image html unsteady_heat_soln.gif "Plot of the exact and computed solutions. " 
\image latex unsteady_heat_soln.eps "Plot of the exact and computed solutions. " width=0.75\textwidth

\image html trace.gif "Time evolution of the computed and exact solutions at a control node, the global error norm and the norm of the solution. " 
\image latex trace.eps "Time evolution of the computed and exact solutions at a control node, the global error norm and the norm of the solution. " width=0.75\textwidth
 
<HR>
<HR> 

\section namespace Global parameters and functions

As usual, we store the problem parameters in a namespace. 

\dontinclude two_d_unsteady_heat.cc
\skipline start_of_ExactSolnForUnsteadyHeat
\until end of ExactSolnForUnsteadyHeat

<HR>
<HR>

\section main The driver code
We start by building the \c Problem object, create a \c  DocInfo
object to label the output files, and open a trace file in which we
will record the time evolution of the solution and the error. 
We choose the length \c t_max of the simulation and 
the (constant) timestep, \c dt.

\skipline start_of_main
\until dt=

Before using any of \c oomph-lib's timestepping functions, 
the timestep \c dt \b must be passed to the \c Problem's
timestepping routines by calling the function
\c Problem::initialise_dt(...) which sets the weights for 
all timesteppers in the problem.

Next we assign the initial conditions by calling 
\c set_initial_condition(), to be discussed in more detail below,
and document the initial conditions. 

\skipline Initialise
\until number()++

Finally, we execute the timestepping loop and document the computed
solutions. 


\skipline Find number
\until end of main

In this loop, each call to the unsteady Newton solver, 
\c Problem::unsteady_newton_solve(...), advances the solution
from its current state, at time \f$ t = \f$ \c Problem::time_pt()->time(),
to \f$ t \ + \f$ \c dt. The unsteady Newton solver automatically
"shifts" the history values "backwards" and advances the
value of the continuous time. 

<HR>
<HR>

\section problem The problem class 

The \c Problem classes for time-dependent problems are very similar 
to those for steady problems. The most important 
additional member functions are \c Problem::actions_before_implicit_timestep()
and  \c Problem::actions_after_implicit_timestep() (both defined as empty
virtual functions in the \c Problem base class) and 
\c set_initial_condition(). The functions 
\c Problem::actions_before_implicit_timestep() and  \c Problem::actions_after_implicit_timestep()
are called automatically by \c oomph-lib's unsteady Newton solver
\c Problem::unsteady_newton_solve() and may be used to update any
time-dependent boundary conditions before the Newton solve, or to perform
any postprocessing steps after a timestep has been completed. 
Here we only use the first of these two functions. 

We note that the (self-explanatory) function \c
set_initial_condition() overwrites an empty virtual function 
in the \c Problem base class. While the assignment of initial 
conditions could, in principle, be performed by any other function, 
e.g. the \c Problem constructor, we strongly recommend using
this function to facilitate the extension to spatial 
adaptivity. (In spatially adaptive computations of time-dependent
problems, a standard interface for the re-assignment of 
initial conditions following mesh adaptations is required; we will
discuss this aspect in  
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">
another example</A>).

Finally, the private member data \c Control_node_pt provides
storage for a pointer to a control \c Node  at which we shall 
document the evolution of the solution. 

\dontinclude two_d_unsteady_heat.cc
\skipline start_of_problem_class
\until end of problem class

<HR>
<HR>

\section constructor The problem constructor

We start by constructing the \c TimeStepper, the second-order accurate
\c BDF<2> timestepper from the \c BDF family, and pass a pointer to it
to the \c Problem, using the member function
\c Problem::add_time_stepper_pt(...).
As the name of this function indicates, \c oomph-lib can operate with
multiple timesteppers -- an essential feature in multi-physics
problems. (For instance, in fluid-structure interaction-problems
timestepping for the solid equations might be performed with a
timestepper from the \c Newmark family, while a \c BDF timestepper
might be used for the Navier--Stokes equations.)
When called for the first time, the function  
\c Problem::add_time_stepper_pt(...) creates the
\c Problem's \c Time object (accessible via \c Problem::time_pt())
with sufficient storage for the history of previous 
timesteps. This is required if the timestep is adjusted during the 
simulation, e.g. when an adaptive timestepper is used. 
(If further \c TimeSteppers which require more storage are added 
subsequently, \c Problem::add_time_stepper(...) updates the amount of 
storage in the \c Problem's \c Time object accordingly).


\dontinclude two_d_unsteady_heat.cc
\skipline start_of_constructor
\until add_time_stepper_pt


Next we set the problem parameters and build the mesh, passing the
pointer to the \c TimeStepper as the \e last argument to the mesh
constructor.

\skipline Setup parameters for exact solution
\until mesh_pt()

The position of the pointer to the timestepper in the list of
arguments to the mesh constructor reflects another \c oomph-lib
convention:

<TABLE>
<TR>
<TD bgcolor="cornsilk">
\anchor conv_mesh_const
<CENTER><B>A general convention</B></CENTER> 

Recall that all \c Data objects store a pointer to a \c TimeStepper
that translates their "history values" into approximations of
the values' time derivatives. The required number of "history values" 
depends on the specific timestepper. For instance, 
a BDF<1> timestepper (the backward Euler scheme) requires storage
of the solution at the previous timestep; the BDF<2> timestepper computes
an approximation of the time-derivative, based on the solution at two
previous timesteps, etc. 

\c Nodes (which \e are \c Data!) are typically built by the \c Mesh
constructor using the \c FiniteElement's member function
\c FiniteElement::construct_node(...). This function takes
a pointer to the timestepper from which the required amount of storage
for the "history values" is extracted. To facilitate the (re-)use of
meshes in steady and time-dependent problems, we adopt the convention 
that

<CENTER>
<TABLE>
<TR>
<TD WIDTH=500 bgcolor="pink">
<CENTER>
The final argument of all mesh
constructors should be the pointer to a \c TimeStepper which 
defaults to \c &Mesh::Default_TimeStepper -- a (static) instantiation
of \c oomph-lib's dummy steady timestepper, \c Steady. 
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>

This convention allows the use of meshes in steady problems
without having to (artificially) create a timestepper. 
The following code fragment illustrates the implementation
of this approach in a mesh constructor:

\code

//=======================================================================
/// Some mesh class
//=======================================================================
template <class ELEMENT>
class SomeMesh : public virtual Mesh
{

public:

  /// \short Constructor: Pass number of elements and pointer to timestepper.
  /// Note that the timestepper defaults to the Steady default timestepper.
  SomeMesh(const unsigned& n_element, TimeStepper* time_stepper_pt=
           &Mesh::Default_TimeStepper)
   {
   
     [...]


     // Allocate storage for all n_element elements in the mesh
     Element_pt.resize(n_element); 

     // Create first element and store it (in its incarnation as        
     // a GeneralisedElement) in the Mesh's Element_pt[] array
     Element_pt[0] = new ELEMENT;
 
     // Create the element's first node and store it in the
     // Mesh's Node_pt[] array. [The member function
     // Mesh::finite_element_pt(...) recasts the pointer to the
     // GeneralisedElement to a pointer to a FiniteElement -- only      
     // FiniteElements have a member function construct_node(...)]
     Node_pt[0] = finite_element_pt(0)->construct_node(0,time_stepper_pt);

     [...]

   }

[...]

};
\endcode

</TD>
</TR>
</TABLE>


Next, we apply the boundary conditions, pinning the values
at all boundary nodes.
\skipline Set the boundary
\until end of set boundary conditions 

Finally, we loop over the elements and pass the pointer to the source
function.


\skipline Complete the build of all elements so they are fully
\until end of constructor


<HR>
<HR>

\section before_timestep Actions before (implicit) timestep
We overload the (empty) virtual function 
\c Problem::actions_before_implicit_timestep() to update the time-dependent
boundary conditions (2), using the current 
value of the continuous time from the \c Problem's \c Time object.

\skipline start of actions_before_implicit_timestep
\until end of actions_before_implicit_timestep

<HR>
<HR>

\section IC Set initial condition
Before starting a time-dependent simulation, 
the current and history values of all \c Data objects must be 
initialised. In a \c BDF timestepping scheme, the history values
represent the solution at previous discrete timesteps. In the
present problem (where the exact solution is known -- admittedly,
a somewhat artificial situation) we can therefore assign the history
values by looping over the previous timesteps
and setting the history values with \c Data::set_value(...).

  \b Important: \c oomph-lib's \c UnsteadyHeatEquations are based on the
Arbitrary-Lagrangian-Eulerian (ALE) formulation of the unsteady 
heat equation to permit computations in moving domains; we will 
illustrate this capability in  
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_ALE/html/index.html">
another example. </a>
In such problems, the nodal positions may 
vary as a function of time. In the present problem, the computation 
is performed in a fixed domain, therefore we initialise the previous 
nodal positions with their current values, accessed via
the member function \c Node::x(t,i) which returns (a reference to) 
the \c i -th nodal coordinate at previous timestep \c t.

<CENTER>
<TABLE>
<TR>
<TD bgcolor="cornsilk">
\anchor conv_t_arg
<CENTER><B>A general convention</B></CENTER> 
Many functions in \c oomph-lib have steady and time-dependent
versions which usually differ in their first argument. Typically,
the first argument of the time-dependent function is an 
additional (unsigned) integer, \c t, 
say. As a general convention, the time-dependent versions return
values (or perform actions) that are appropriate for the current 
time if called with \c t=0, and return values (or perform actions) that are 
appropriate for previous time levels if the argument is set to \c
t>0. Note that we refer
to previous time levels, rather than previous timesteps, because  history
values do not necessarily represent values at previous timesteps,
as they do for \c BDF schemes. For instance, in \c Newmark
timestepping schemes, the history values include approximations of
the first and second time-derivatives at the previous timestep. 

<HR>
\anchor conv_hist_vals
<CENTER><B>Another general convention</B></CENTER> 
While, in general, not all "history values" represent
the solution at previous timesteps, the "history values"
that do, should be (and, for any existing \c TimeSteppers, are) 
listed before those that represent other quantities.
The number of history values required/used by a \c TimeStepper
may be obtained from its member function

\code
TimeStepper::nprev_values()
\endcode

As an example, \c BDF<4>::nprev_values() returns 4. 
The total number of history values (including the current value!)
is returned by

\code
TimeStepper::ntstorage()
\endcode

As an example, \c BDF<4>::ntstorage() returns 5.
</TD>
</TR>
</TABLE>
</CENTER>


Here is the source code for the \c  set_initial_condition() function:
\skipline start_of_set_initial_condition
\until end of set_initial_condition

<HR> 
<HR>

\section doc Post-processing
As in many previous examples, this member function outputs 
the computed solution, the exact solution and the error. 
We augment the solution data by tecplot text and geometries
to facilitate the visualisation and record the time evolution of
the solution and the error in the trace file.

\skipline start_of_doc_solution
\until end of doc_solution

<HR>
<HR>

\section comments Comments and Exercises
The current example only illustrates the most basic timestepping procedures.
In subsequent examples we will demonstrate 
\c oomph-lib's  
<A HREF="../../two_d_unsteady_heat2/html/index.html">
dump and restart functions,</A> 
<A HREF="../../two_d_unsteady_heat_t_adapt/html/index.html">the use of adaptive
timestepping,</A> 
<A HREF="../../two_d_unsteady_heat_adapt/html/index.html">the use of spatial
adaptivity,</A> 
<A HREF="../../two_d_unsteady_heat_ALE/html/index.html">
computations in moving domains,</A> and 
<A HREF="../../two_d_unsteady_heat_2adapt/html/index.html">
the combination of temporal and spatial adaptivity</A>.

We stress that setting the initial conditions in a "real" 
problem often presents a delicate step, especially if 
higher-order timesteppers from the BDF family
are used. This is because in the absence of a known exact solution,
the initial condition (3) only provides enough 
information to determine a single "history value" at each node -- the 
value at the initial time. \c oomph-lib's timestepping procedures 
provide several functions that allow the simulation to be initiated
with an "impulsive start", corresponding to a past history in which
the boundary condition (3) describes the system's
state for all \f$ t \le 0\f$ rather than only \e at \f$ t =
0\f$ . For instance, the top-level function \c
Problem::assign_initial_values_impulsive() sets the  "history values"
of all \c Data objects and the \c Nodes' "positional history values" 
to values that are appropriate for an impulsive start from the
currently assigned nodal values and positions. 
The following exercises aim to explore this functionality.

\subsection exercises Exercises
-# Replace the call to \c problem.set_initial_condition() 
   in the main function by a call to 
   \c Problem::assign_initial_values_impulsive() and analyse
   the results. [Hint: When \c Data objects are created, their 
   values are initialised to zero.]
-# Confirm that the loop over the coordinate directions
   \dontinclude two_d_unsteady_heat.cc 
   \skipline Loop over coordinate directions
   \until }
   in \c set_initial_condition(), can be replaced by
   \code
   time_stepper_pt()->assign_initial_positions_impulsive(
                        mesh_pt()->node_pt(n));
   \endcode
   [This statement could then be moved outside the loop 
   over the previous time levels.] 
-# Confirm that the initialisation of the previous nodal positions
   is essential by commenting out this step -- see the  
   <A HREF="../../two_d_unsteady_heat_ALE/html/index.html">ALE 
   example</A> for further details on computations in moving domains.
-# Overwrite the correct assignment of the "history values" in
   \c set_initial_condition() by adding the statement
   \c Problem::assign_initial_values_impulsive() at the end
   of this function (rather than bypassing the assignment
   completely as in the first exercise). Repeat this with the \c BDF<1>
   and \c BDF<4> timesteppers and explain the different behaviour. 
-# Examine the accuracy of the various \c BDF timesteppers by
   re-running the simulations with various timesteppers and with
   different timesteps.
        


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat/">
demo_drivers/unsteady_heat/two_d_unsteady_heat/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat/two_d_unsteady_heat.cc">
demo_drivers/unsteady_heat/two_d_unsteady_heat/two_d_unsteady_heat.cc
</A>
</CENTER>
.



















<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


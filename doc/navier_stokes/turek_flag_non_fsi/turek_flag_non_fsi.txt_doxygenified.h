/**

\mainpage Demo problem: Flow past a cylinder with a waving flag


In this example we consider the flow in a 2D channel that contains a cylinder
with a waving "flag". This is a "warm-up problem"
for the solution of Turek & Hron's FSI benchmark problem
discussed in <A HREF="../../../interaction/turek_flag/html/index.html">
another tutorial.</A>

<HR> 
<HR>

\section the_problem The Problem

The figure below shows a sketch of the problem: A 2D channel of
height \f$ H^*\f$  and length \f$ L^* \f$ contains 
a cylinder of diameter \f$ d^* \f$, centred at 
\f$ (X^*_c, Y^*_c) \f$ to which a "flag" of thickness 
\f$ H^*_{flag}\f$  and length \f$ L^*_{flag} \f$ is attached.
We assume that the flag performs time-periodic oscillations with 
period \f$ T^*. \f$  Steady Poiseuille flow with average velocity \f$ U^* \f$
is imposed at the left end of the channel while we assume the 
outflow to be parallel and axially traction-free.

\image html turek_flag_non_fsi_dim.gif "Sketch of the problem in dimensional variables. " 
\image latex turek_flag_non_fsi_dim.eps "Sketch of the problem in dimensional variables. " width=0.75\textwidth
   
 

We non-dimensionalise all length and coordinates 
on the diameter of the cylinder, \f$ d^* \f$, time on the natural timescale of the flow, \f$ d^*/U^* \f$, the
velocities on the mean velocity,  \f$ U^* \f$, and the pressure
on the viscous scale. The problem is then
governed by the non-dimensional Navier-Stokes equations
\f[
Re \left( St \frac{\partial u_i}{\partial t} + 
u_j \frac{\partial u_i}{\partial x_j} \right) =
- \frac{\partial p}{\partial x_i}  +
\frac{\partial }{\partial x_j} \left[ \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right) \right],
\f] 
where \f$ Re = \rho U^* H_0^* / \mu \f$ and \f$ St = 1 \f$,
and
<CENTER>
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
</CENTER>
subject to parabolic inflow
<CENTER>
\f[
{\bf u} = 6 x_2 (1-x_2) {\bf e}_1
\f] 
</CENTER>
at the inflow cross-section; parallel, axially-traction-free outflow
at the outlet; and no-slip on the stationary channel walls
and the surface of the cylinder, \f$ {\bf u} = {\bf 0} \f$. 
The no-slip condition on the moving flag is
<CENTER>
\f[
{\bf u} = \frac{\partial {\bf R}_w(\xi_{[top,tip,bottom]},t)}{\partial t}
\f] 
</CENTER>
where \f$ \xi_{[top,tip,bottom]} \f$ are Lagrangian coordinates 
parametrising the three faces of the flag.
The flag performs oscillations with non-dimensional period
\f$ T = T^* U^* / H^*_{tot} \f$. Here is a sketch of 
the non-dimensional version of the problem:

\image html turek_flag_non_fsi.gif "Sketch of the problem in dimensionless variables, showing the Lagrangian coordinates that parametrise the three faces of the flag. " 
\image latex turek_flag_non_fsi.eps "Sketch of the problem in dimensionless variables, showing the Lagrangian coordinates that parametrise the three faces of the flag. " width=0.75\textwidth
 


<HR> 
<HR>

\section results Results

The figure below shows a snapshot of the flow field  
(pressure contours and instantaneous streamlines)  for a Reynolds
number of \f$ Re=100 \f$ and an oscillation period of \f$ T= 10 \f$, as 
well as the corresponding fluid mesh. Note how \c oomph-lib's automatic 
mesh adaptation has refined the mesh  in the high-shear
regions near the front of the cylinder and at the trailing edge
of the flag.

\image html flow_cropped.gif "Mesh (top) and flow field (bottom). " 
\image latex flow_cropped.eps "Mesh (top) and flow field (bottom). " width=0.75\textwidth

The corresponding <A HREF="../figures/turek_flag_non_fsi.avi"> 
animation</A> illustrates the algebraic node update strategy 
(implemented with an \c AlgebraicMesh, discussed in more detail
in  <A HREF="../../algebraic_collapsible_channel/html/index.html">
another tutorial</A>) and the evolution of
the flow field. Note that the instantaneous streamlines intersect
the (impermeable) flag because the flag is not stationary. 



<HR> 
<HR>

\section parameters The global parameters

As usual we use a namespace to define the (single) global parameter,
the Reynolds number.

\dontinclude turek_flag_non_fsi.cc
\skipline start_of_global_parameters
\until }

<HR> 
<HR>


\section leaflet Specification of the flag geometry

We specify the flag geometry and its time-dependent motion by
representing its three faces  by \c GeomObjects, each parametrised
by its own Lagrangian coordinate, as shown in the sketch above.
The geometry of the cylinder is represented by one of \c oomph-lib's
generic \c GeomObjects, the \c Circle.

We enclose the relevant data in its own namespace and 
start by defining the period of the oscillation, and the geometric
parameters controlling the flag's initial shape and its subsequent motion.

\dontinclude turek_flag_non_fsi.cc
\skipline start_of_flag_definition
\until Time_pt

We choose the motion of the flag such that it vaguely resembles 
that expected in the corresponding FSI problem: We subject the flag's 
upper and lower faces to purely vertical displacements that deform 
them into a fraction of a sine wave, while keeping the face 
at the tip of the flag straight and vertical. We implement 
this by prescribing the time-dependent motion 
of the two vertices at the tip of the flag with two functions:

\until end of bottom tip

This information is then used in three custom-written \c GeomObject that
define the time-dependent shape of the three faces. Here is the
one describing the shape of the upper face:

\until };

[We omit the definition of the other two \c GeomObjects, 
\c BottomOfFlag and \c TipOfFlag in in the interest of brevity.
Their definitions can be found in the <A HREF="
../../../../demo_drivers/navier_stokes/turek_flag_non_fsi/turek_flag_non_fsi.cc
">source code</A>.] We provide storage for the pointers to the 
four \c GeomObjects required for the representation of the flag,

\skipline Pointer to GeomObject
\until Cylinder_pt

and provide a setup function that generates the \c GeomObjects 
and stores the pointer to the global \c Time object that 
will be created by the \c Problem:

\until end of namespace

<HR>
<HR>



\section mesh The mesh

The figure below shows a sketch of the (unrefined) mesh and the
enumeration of its boundaries. Various different implementations 
of the mesh exist, allowing the user to perform the node-update in response
to the motion of the flag by a \c Domain/MacroElement - based 
procedure, or by using an algebraic node update.
In either case, only the nodes in the elements that are shaded
in yellow (or refined elements that are generated from these)
participate in the node-update. 

\image html cylinder_with_flag_mesh.gif "The (unrefined) mesh. Only nodes in the yellow regions participate in the node-update. " 
\image latex cylinder_with_flag_mesh.eps "The (unrefined) mesh. Only nodes in the yellow regions participate in the node-update. " width=0.75\textwidth

The node update strategy is illustrated in the 
<A HREF="../figures/turek_flag_non_fsi.avi"> animation</A> 
of the flow field and the mesh motion.
 
<HR>
<HR>

\section main The driver code

The driver code has the usual structure for a time-dependent problem:
We store the command line arguments, create a \c DocInfo object,
and assign the parameters that specify the dimensions of the channel.

\dontinclude turek_flag_non_fsi.cc
\skipline start_of_main
\until height=4.1;

We build the problem using a mesh in which
the node update is performed by the \c AlgebraicMesh class. (The 
driver code also provides the option to use a \c Domain/MacroElement -
based node update -- this option is chosen via suitable \#ifdefs;
see the <A HREF="
../../../../demo_drivers/navier_stokes/turek_flag_non_fsi/turek_flag_non_fsi.cc
">source code</A> for details).
 

\skipline Create Problem with AlgebraicMesh-based node update
\until height);

Next we set up the time-stepping (as usual, fewer timesteps are 
performed during a validation run which is identified by a non-zero 
number of command line arguments):

\skipline Number of timesteps per period
\until initialise_dt(dt);


We start the simulation with a steady solve, allowing up to
three levels of adaptive refinement (fewer if we are performing 
a validation run):

\until doc_info.number()++;

Finally, we enter the proper timestepping loop, allowing one
spatial adaptation per timestep and suppressing the re-assignment
of initial conditions following an adaptation by setting the
parameter \c first to false. (See the discussion of timestepping
with automatic mesh adaptation in <a href="../../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">another tutorial</a>.)
 
\until }//end of main


<HR> 
<HR>



\section problem The Problem class

The problem class provides an access function to the mesh (note
the use of \#ifdefs to choose the correct one), and defines the
interfaces for the usual member functions, either empty
or explained in more detail below.

\dontinclude turek_flag_non_fsi.cc
\skipline start_of_problem_class
\until end_of_problem_class


<HR>
<HR>

 

\section constructor The problem constructor

The initial assignment of zero velocity and pressure provides a
very poor initial guess for the preliminary steady Newton solve.
We therefore increase the maximum residuals and iteration counts
allowed in the Newton iteration before creating the timestepper
and setting up the \c GeomObjects required to parametrise the
flag:


\skipline start_of_constructor
\until Flag_definition

Next we build the mesh, passing the pointers to the various \c GeomObjects
and the geometric parameters to its constructor.

\skipline // Build mesh (with AlgebraicMesh-based node update)
\until );

We perform two rounds of uniform refinement (the Newton solver
does not converge on coarser meshes) before creating an
error estimator for the subsequent automatic mesh adaptation.

\skipline Refine uniformly
\until error_estimator_pt;

Both velocity components are imposed by Dirichlet conditions
(either via a no-slip condition or via the imposed inflow profile) on all
boundaries, apart from the outflow cross-section (mesh boundary 1), 
where the axial velocity is unknown.

\until done bc

Next we complete the build of the elements, passing the relevant
pointers to physical parameters,

\until }

and impose the steady Poiseuille profile at the inlet (mesh boundary 3).

\until }

Finally, we pin the redundant pressure degrees of
freedom (see 
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details), and assign the equations numbers.


\until end of constructor





 
<HR>
<HR>

\section after_adapt Actions before adapt

After each adaptation, we unpin and re-pin all redundant
pressures degrees of freedom (see 
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details). Since the inflow profile is parabolic,
it is interpolated correctly from "father" to "son" elements during
mesh refinement so no further action is required.

\skipline actions_after_adapt
\until end_of_actions_after_adapt


<HR>
<HR>

\section timestep Actions before the timestep

Before each timestep we update the nodal positions in the mesh and
re-apply the no-slip condition at the nodes on mesh boundaries 5 - 7.

\skipline start_of_actions_before_implicit_timestep
\until end_of_actions_before_implicit_timestep


<HR>
<HR>

\section doc Post Processing

The function doc_solution(...) documents the results.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>

\section ex_com Comments and Exercises

- Compare the results obtained when the node update is
  performed with the algebraic or the \c MacroElement/Domain-based
  node update. (\c oomph-lib's build machinery will automatically
  generate both versions of the code, using the \c -DUSE_MACRO_ELEMENTS
  compilation flag). 
. 


<HR>
<HR>

\section ackn Acknowledgements


- This code was originally developed by 
  Stefan Kollmannsberger and his students Iason Papaioannou and  
  Orkun Oezkan Doenmez. It was completed by Floraine Cordier.

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/navier_stokes/turek_flag_non_fsi/
">
demo_drivers/navier_stokes/turek_flag_non_fsi/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/navier_stokes/turek_flag_non_fsi/turek_flag_non_fsi.cc
">
demo_drivers/navier_stokes/turek_flag_non_fsi/turek_flag_non_fsi.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


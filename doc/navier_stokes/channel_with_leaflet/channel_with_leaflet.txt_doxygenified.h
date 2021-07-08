/**

\mainpage Demo problem: Flow in a 2D channel with an oscillating leaflet
 
In this example we consider the flow in a 2D channel which is partially
obstructed by an oscillating leaflet. We consider the case where the 
motion of the leaflet is prescribed -- this is a "warm-up exercise" for 
the <A HREF="../../../interaction/fsi_channel_with_leaflet/html/index.html">
corresponding FSI problem </A> in which the leaflet is
an elastic structure.
  
<HR> 
<HR>

\section the_problem The Problem

The figure below shows a sketch of the problem: A 2D channel of
height  \f$ H^*_{tot}\f$  and length \f$ L^*_{left} + L^*_{right} \f$ 
is partially occluded by a (zero-thickness) leaflet of height 
\f$ H^*_{leaflet} \f$.  The leaflet is
parametrised by a Lagrangian coordinate \f$ \xi^* \f$ so that the
position vector to a material point on the leaflet is given by
\f$ {\bf R}^*_w(\xi^*,t^*) \f$, and we assume that the leaflet
performs time-periodic oscillations with period \f$ T^*. \f$ 
Steady Poiseuille flow with average velocity \f$ U^* \f$
is imposed at the left end of the channel while we assume that the 
outflow is parallel and axially traction-free.

\image html channel_with_leaflet_dim.gif "Sketch of the problem in dimensional " 
\image latex channel_with_leaflet_dim.eps "Sketch of the problem in dimensional " width=0.75\textwidth
  
 

We non-dimensionalise all length and coordinates 
on the channel width, \f$ H^*_{tot}
\f$ , time on the natural timescale of the flow, \f$ H^*_{tot}/U^* \f$, the
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
at the outlet; and no-slip on the stationary channel walls, 
\f$ {\bf u} = {\bf 0} \f$. The no-slip condition on the
leaflet is
<CENTER>
\f[
{\bf u} = \frac{\partial {\bf R}_w(\xi,t)}{\partial t}
\f] 
</CENTER>
and the leaflet performs oscillations with non-dimensional period
\f$ T = T^* U^* / H^*_{tot} \f$. Here is a sketch of 
the non-dimensional version of the problem:

\image html channel_with_leaflet.gif "Sketch of the problem in dimensionless variables. " 
\image latex channel_with_leaflet.eps "Sketch of the problem in dimensionless variables. " width=0.75\textwidth
 
An interesting feature of this problem is that even though the leaflet
is assumed to have negligible thickness its presence will generate a pressure
jump between its two faces. (The velocities are continuous
because the no-slip condition imposes the same velocity
on both faces.) When discretising the problem with
\c QTaylorHood elements (for which the pressure varies continuously 
across the element boundaries), the mesh must therefore be
"opened up" with a cut along the position of the leaflet. 
This is done in the constructor of the 
<A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
<CODE>ChannelWithLeafletMesh</CODE></a>
which forms the basis of \c RefineableAlgebraicChannelWithLeafletMesh
used to discretise this problem.  (See \ref comm_and_ex for a more
detailed discussion of the mesh.)

<HR> 
<HR>

\section results Results

The figure below shows a snapshot of the flow field  
(pressure contours and instantaneous streamlines)  for a Reynolds
number of \f$ Re=20 \f$ and an oscillation period of \f$ T=20 \f$, as 
well as the corresponding fluid mesh. Note how \c oomph-lib's automatic 
mesh adaptation has refined the mesh near the tip of the leaflet 
where the pressure has a singularity.
 
\image html channel_with_leaflet_flow_cropped.gif "Mesh (top) and flow field (bottom). " 
\image latex channel_with_leaflet_flow_cropped.eps "Mesh (top) and flow field (bottom). " width=0.75\textwidth

The corresponding <A HREF="../figures/channel_with_leaflet_flow.avi">
animation</A> illustrates the algebraic node update strategy 
(implemented with an \c AlgebraicMesh, discussed in more detail
in  <A HREF="../../algebraic_collapsible_channel/html/index.html">
another tutorial</A>) and the evolution of
the flow field. Note that the instantaneous streamlines intersect
the (impermeable) leaflet because the leaflet is not stationary. 

<HR> 
<HR>

\section parameters The global parameters

As usual we use a namespace to define the (single) global parameter,
the Reynolds number.

\dontinclude channel_with_leaflet.cc
\skipline start_of_global_parameters
\until end_of_namespace

<HR> 
<HR>


\section leaflet Specification of the leaflet geometry

We specify the leaflet geometry and its time-dependent motion by
representing it as  a \c GeomObject. The \c GeomObject has one Lagrangian
and two Eulerian coordinates, and its geometry is characterised
by its length, the x-coordinate of its origin, \f$ X_0 \f$, 
and the period and  amplitude of the horizontal and vertical tip deflection. 

\dontinclude channel_with_leaflet.cc
\skipline start_of_leaflet_class
\until end_of_the_GeomObject


<HR>
<HR>
\section main The driver code

We store the command line arguments, create a \c DocInfo object,
and assign the parameters that specify the domain and the
leaflet geometry:
 

\dontinclude channel_with_leaflet.cc
\skipline start_of_main
\until  unsigned ny2=2; 

Next we build the problem and assign the time-stepping parameters
(as usual, fewer timesteps are used during a validation run which is
identified by a non-zero number of command line arguments):

\until initialise_dt(dt);


We start the simulation with a steady solve, allowing up to
five levels of adaptive refinement (fewer if we are performing 
a validation run):

\until doc_info.number()++;

Finally, we enter the proper timestepping loop, allowing one
spatial adaptation per timestep and suppressing the re-assignment
of initial conditions following an adaptation by setting the
parameter \c first to false (see the discussion of timestepping
with automatic mesh adaptation in <a href="../../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">another tutorial</a>.)

\until }//end of main


<HR>
<HR>

\section problem The Problem class

The problem class has the usual member functions to perform
actions after the mesh adaptation and before every implicit timestep:

\dontinclude channel_with_leaflet.cc
\skipline start_of_problem_class
\until };

<HR>
<HR>



\section constructor The problem constructor


We construct the timestepper and the \c GeomObject that represents
the leaflet, and pass it pointers to them to
the constructor of the 
\c RefineableAlgebraicChannelWithLeafletMesh (discussed in more
detail in \ref comm_and_ex).

\until  time_stepper_pt()); 

Next we create the spatial error estimator and 
loop over the  elements to set the pointers
to the relevant physical parameters.

\until // end loop over elements

The velocity is prescribed everywhere apart from the outflow boundary
(boundary 1; see the sketch in \ref comm_and_ex for the enumeration
of the mesh boundaries). Along the inflow (boundary 3) we apply a parabolic
velocity profile with unit flux:

\until  }// end of setup boundary condition

Finally, we pin the redundant pressure degrees of
freedom (see 
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details), and assign the equations numbers.

\until }//end of constructor


<HR>
<HR>



\section timestep Actions before the timestep

Before each timestep we update the nodal positions in the mesh and
re-apply the no-slip condition on the nodes of the moving leaflet
(boundaries 4 and 5; see the sketch in \ref comm_and_ex for the enumeration
of the mesh boundaries).

\skipline start_of_actions_before_implicit_timestep
\until end_of_actions_before_implicit_timestep


<HR>
<HR>




\section adapt Actions after the mesh adaptation

Once the mesh has been adapted, we free all pressure degrees of
freedom and then (re-)pin any redundant ones (see 
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details):

\skipline start_of_actions_after_adaptation
\until end_of_actions_after_adapt

Note that the default interpolation of the (quadratic!) inflow velocity 
profile from father to son elements during the mesh adaptation 
already ensures that the inflow profile remains quadratic, therefore
no further action is required.


<HR>
<HR>

\section doc Post-processing

The function \c doc_solution(...) documents the results.

\skipline start_of_doc_solution
\until end_of_doc_solution



<HR>
<HR>

\section comm_and_ex Further comments and exercises

\subsection comm Further comments: The algebraic node update procedure

The figure below illustrates the algebraic node update procedure
employed in the \c RefineableAlgebraicChannelWithLeafletMesh. 
The mesh employs four different node update functions,
depending on which region a node is located in: Nodes in 
region I (or II) are located on straight lines that connect the
upstream (or downstream) boundary with the leaflet; nodes in
region III (or IV) are located on straight lines that connect the
upstream and (or downstream) boundary with the straight line
from the tip of the leaflet to upper channel wall. 

\image html channel_with_leaflet_mesh.gif "Sketch illustrating the algebraic node update procedure. " 
\image latex channel_with_leaflet_mesh.eps "Sketch illustrating the algebraic node update procedure. " width=0.75\textwidth

The implementation of the node update functions is straightforward
and can be found in the source files
<CENTER>
<A HREF="
../../../../src/meshes/channel_with_leaflet_mesh.template.h">
src/meshes/channel_with_leaflet_mesh.template.h
</A>
</CENTER>\n
and 
<CENTER>
<A HREF="
../../../../src/meshes/channel_with_leaflet_mesh.template.cc">
src/meshes/channel_with_leaflet_mesh.template.cc
</A>
</CENTER>\n
which also illustrate how the mesh is constructed by inheritance
from the \c SimpleRectangularQuadMesh (the main task being the
creation of additional nodes in the interior to "cut open" the mesh
along the position of the leaflet). The source files also
contain other versions of the mesh in which the node update is
performed with Domain/MacroElements, using the technique described in 
<A HREF="../../../poisson/fish_poisson2/html/index.html"> another
tutorial</A>.

\subsection ex Exercise

With the node update strategy illustrated above, the position of 
\e all nodes in the fluid mesh has to be updated when the leaflet 
moves. This is not a particular problem in the current application
where the node-update is only performed once per timestep. However,
in the <A HREF="../../../interaction/fsi_channel_with_leaflet/html/index.html">
corresponding FSI problem </A>, the approach is costly because
of the of large number of shape derivatives to be computed. 

As an exercise, we suggest to make the node-update procedure 
more efficient by sub-dividing the
regions upstream and downstream of the leaflet into a central
section in which the nodes move in response to the motion of the leaflet
(the old regions I-IV) and two additional regions (regions V and VI)
in which they remain stationary. This is easy because, as explained <A HREF="../../algebraic_collapsible_channel/html/index.html#comm">elsewhere,</A>
all \c AlgebraicNodes already have a default node update function
that leaves them stationary.


   
\image html channel_with_leaflet_better_mesh.gif "A better node update strategy. " 
\image latex channel_with_leaflet_better_mesh.eps "A better node update strategy. " width=0.75\textwidth



<HR>
<HR>

\section ackn Acknowledgements

- This code was originally developed by Floraine Cordier.

<HR>
<HR> 


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/navier_stokes/channel_with_leaflet/
">
demo_drivers/navier_stokes/channel_with_leaflet/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/navier_stokes/channel_with_leaflet/channel_with_leaflet.cc
">
demo_drivers/navier_stokes/channel_with_leaflet/channel_with_leaflet.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


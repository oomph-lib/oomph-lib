/**

\mainpage Demo problem: Flow in a 2D channel with an elastic leaflet 



In this example we consider the flow in a 2D channel which is partially
obstructed by an elastic leaflet -- the  FSI generalisation of the
<A HREF="../../../navier_stokes/channel_with_leaflet/html/index.html">
Navier-Stokes problem in which the motion of the leaflet is
prescribed.</A> A particular feature of this problem is that the
leaflet (modelled as a thin-walled beam structure) is totally immersed
in the fluid and is therefore exposed to the fluid traction from
both sides. 

The problem presented here was used as one of the test cases
for \c oomph-lib's FSI preconditioner; see \n\n

<CENTER>
<A href="http://www.maths.man.ac.uk/~mheil/oomph_lib_additional_material/HeilHazelBoyleCompMech/HeilHazelBoyleCompMech.pdf">Heil, M., Hazel, A.L. & Boyle, J. (2008): Solvers for large-displacement fluid-structure interaction problems: Segregated vs. monolithic approaches. Computational Mechanics.</A>
</CENTER>

In this tutorial we concentrate on the problem formulation. 
The application of the preconditioner is discussed 
<a href="../../../preconditioners/fsi/html/index.html">elsewhere</a> --
the required source code is contained in the
<A href="../../../../demo_drivers/interaction/fsi_channel_with_leaflet/fsi_channel_with_leaflet.cc">driver code.</A>

<HR> 
<HR>

\section the_problem The Problem

The figure below shows a sketch of the problem:  A 2D channel of
height  \f$ H^*_{tot}\f$ and length \f$ L^*_{left} + L^*_{right} \f$ 
is partially occluded by a thin-walled elastic leaflet of height 
\f$ H^*_{leaflet} \f$. (As usual, we use asterisks to distinguish
dimensional quantities from their non-dimensional equivalents to be
introduced later.) We assume that the leaflet is
sufficiently thin so that, as far as the boundary conditions
for the fluid are concerned, the leaflet can be assumed to 
be infinitely thin. This allows us to parametrise its shape 
by a single Lagrangian coordinate \f$ \xi^* \f$. Hence we write the
position vector to a material point on the leaflet as
\f$ {\bf R}^*_w(\xi^*,t^*) \f$. 
A pulsatile Poiseuille flow whose mean velocity fluctuates between
\f$ U^* \f$ and \f$ 2 U^* \f$ is imposed at the upstream end of the
channel. The outflow is assumed to be parallel and axially 
traction-free.

\image html fsi_channel_with_leaflet_dim.gif "Sketch of the problem. " 
\image latex fsi_channel_with_leaflet_dim.eps "Sketch of the problem. " width=0.75\textwidth

We non-dimensionalise all length and coordinates on the channel 
width, \f$ H^*_{tot} \f$ , time on the natural timescale of the flow,
\f$ H^*_{tot}/U^* \f$, the velocities on the (minimum) mean velocity,  
\f$ U^* \f$, and the pressure on the viscous scale. 

The fluid flow is then governed by the non-dimensional Navier-Stokes equations
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
subject to a time-periodic, parabolic inflow with non-dimensional
period \f$ T \f$,
<CENTER>
\f[
{\bf u} = 6 x_2 (1-x_2) \left( 1 + \frac{1}{2} (1-\cos(2 \pi t/T)) \right) 
{\bf e}_1
\f]
</CENTER>
in the inflow cross-section; parallel, axially-traction-free outflow
at the outlet; and no-slip on the stationary channel walls, 
\f$ {\bf u} = {\bf 0} \f$. The no-slip condition on the
leaflet is
<CENTER>
\f[
{\bf u} = \frac{\partial {\bf R}_w(\xi,t)}{\partial t}.
\f] 
</CENTER>

We model the leaflet as a thin-walled, massless, elastic Kirchhoff-Love
beam of wall thickness \f$ h^* \f$. The beam's effective (1D) elastic modulus
is given by \f$ E_{eff} = E/(1-\nu^2), \f$ where \f$ E \f$ and \f$ \nu \f$
are its 3D Young's modulus and Poisson's ratio, respectively.
The beam's deformation is governed by the principle of virtual
displacements, discussed in detail 
<A HREF="../../../beam/tensioned_string/html/index.html">in another 
tutorial.</A> 
As in the Navier-Stokes equations, we scale all lengths in
the beam problem on the channel's width, \f$ H^*_{tot}. \f$ The
non-dimensional position vector
\f$ {\bf r}_w(\xi) \f$ to the undeformed wall is then given by
\f[
{\bf r}_w(\xi) = 
{\bf R}_w(\xi,t=0) = 
\left(
\begin{array}{c}
X_0  \\
\xi
\end{array}
\right) \ \ \ \ \ \ \mbox{where $\xi \in [0,H_{leaflet}].$}
\f]
Our non-dimensionalisation of the principle of virtual displacements
requires all stresses and tractions to be non-dimensionalised 
on the beam's (effective 1D) elastic modulus, \f$ E_{eff} \f$. 
The non-dimensional load vector \f$ {\bf f} = {\bf f}^{*}/E_{eff} \f$ 
that acts on the leaflet (combining the fluid tractions acting
on its front and back faces) is then given by
\f[
f_i = Q \bigg\{ 
\left( p\big|_{front} \ N_i^{[front]} -
\left(
\frac{\partial u_i}{\partial x_j} +
\frac{\partial u_j}{\partial x_i} 
\right)\bigg|_{front} N_j^{[front]} \right)
+
\f]
\f[
. \hspace{3cm} + \left( p\big|_{back} \ N_i^{[back]} -
\left(
\frac{\partial u_i}{\partial x_j} +
\frac{\partial u_j}{\partial x_i} 
\right)\bigg|_{back} N_j^{[back]} \right) 
\bigg\}
\ \ \ \ \ \ \mbox{for $i=1,2,$}
\f]
where 
\f[
Q = \frac{\mu U}{E_{eff} H_{tot}^*}
\f]
is the ratio of the fluid pressure scale, \f$ \mu U/H_{tot}^* \f$ ,
used to non-dimensionalise the Navier-Stokes equations, to the
beam's effective elastic modulus, \f$ E_{eff} \f$. The parameter \f$ Q
\f$ therefore indicates the strength of the fluid-structure
interaction. In particular, if \f$ Q=0 \f$ the leaflet does not
"feel" the fluid traction. \f$ {\bf N}_{front} \f$ and 
\f$ {\bf N}_{back} \f$ are the outer unit normals on the "front" and "back"
faces of the deformed leaflet, as shown in this sketch of 
the non-dimensional version of the problem: 

\image html fsi_channel_with_leaflet.gif "Sketch of the problem in dimensionless variables. " 
\image latex fsi_channel_with_leaflet.eps "Sketch of the problem in dimensionless variables. " width=0.75\textwidth
 
<HR> 
<HR>

\section results Results

The figure below shows (a) a snapshot of the flow field  
(pressure contours and instantaneous streamlines) and (b) the time-trace
of the horizontal position of the leaflet's tip for a Reynolds number
of \f$ Re = 200 \f$, a non-dimensional period of \f$ T=2 \f$,
and an interaction parameter of \f$ Q = 10^{-6}
\f$. Following the decay of initial transients the leaflet performs
periodic large-amplitude oscillations with the period of the pulsating
inflow.  Large velocity gradients develop at the front of the leaflet and 
in the shear layer that emanates from its tip and separates the 
recirculating flow region behind the leaflet from the main flow.  Fig. (c) 
illustrates the non-uniform mesh refinement and shows
the improved resolution in the high-shear regions, particularly near 
the leaflet's tip where the pressure is singular. The mesh was 
continuously adapted throughout the simulation and contained an average of
about 32,000 degrees of freedom. This is a fraction of the 1,324,343 degrees 
of freedom that would be required to achieve the same local resolution via 
uniform mesh refinement.
 

 
\image html flow_from_paper.gif "Computational results. " 
\image latex flow_from_paper.eps "Computational results. " width=0.75\textwidth

The corresponding <A HREF="../figures/fsi_channel_with_leaflet_flow.avi">
animation</A> illustrates the algebraic node update strategy 
(implemented with an \c AlgebraicMesh, discussed in more detail
in  <A HREF="../../fsi_collapsible_channel_algebraic/html/index.html">
another tutorial</A>) and the evolution of
the flow field. Note that the instantaneous streamlines intersect
the (impermeable) leaflet because the leaflet is not stationary. 
The animation shows how the mesh adapts itself to changes in the flow
field -- using much smaller elements in high-shear regions,
including the artificial outflow boundary layer that is created
by the imposition of parallel outflow in a region where the 
flow is far from fully-developed.

<HR> 
<HR>


\section parameters The global parameters

As usual we use a namespace to define the problem parameters:
The Reynolds number and its product with the Strouhal number 
(both initialised to 50),

\dontinclude fsi_channel_with_leaflet.cc
\skipline start_of_global_parameters
\until ReSt=50.0;

the leaflet's non-dimensional wall thickness, \f$ h = h^*/H^*_{tot}
\f$, and the FSI parameter,

\until Q=1.0e-6; 

and the parameters that define the magnitude of the pulsating inflow.

\until end_of_namespace


<HR>  
<HR>

\section leaflet The undeformed leaflet

We use a \c GeomObject to describe the initial, stress-free shape
of the elastic leaflet: a vertical straight line. The member function
\c d2position(...) provides the first and second derivatives of
the position vector, as required by the variational principle
that governs the beam's deformation;
see <A HREF="../../../beam/tensioned_string/html/index.html">
the beam theory tutorial for details. </A>

\skipline start_of_undeformed_leaflet
\until end_of_undeformed_wall

<HR>
<HR>


\section main The driver code

As with most time-dependent demo codes, the specification of 
a non-zero number of command line arguments will be interpreted 
as the code being run as part  of \c oomph-lib's self-test in which
case only a few timesteps will be performed. Therefore we start by storing
the command line arguments in the namespace \c CommandLineArgs
to make them accessible throughout the code.

\skipline start_of_main
\until ::setup(

Next we specify the geometry and mesh parameters, and build 
the \c Problem object, using the \c AlgebraicElement version of
the two-dimensional \c RefineableQTaylorHoodElements.

\until ,ny2,x_0);

We prepare a \c DocInfo object and open a trace file to document
the system's evolution.

\until trace.open

Next we assign the timestepping parameters (using fewer timesteps if the
code is run during a self-test) and document the system's initial 
state which provides the initial guess for the Newton iteration in 
the subsequent steady solve.

\until doc_info.number()++; 

We wish to start the time-dependent simulation from an initial
condition that corresponds to the steady solution of the problem 
with unit influx (i.e. the steady solution obtained if the 
time-dependent boundary condition that is applied at the inlet is held
fixed at its value for \f$ t= 0 \f$). For the required
Reynolds number of \f$ Re=200 \f$, the currently assigned
initial guess (zero velocity and pressure, with the wall  
in its undeformed position) is "too far" from the actual solution
and the Newton method would diverge (try it!). Therefore we
generate the initial steady solution in a preliminary sequence
of steady computations in which we first compute the solution for
a Reynolds number of  \f$ Re=50 \f$ (for this value the Newton 
method does converge).
This solution is then used as the initial guess for the computation
at  \f$ Re=100 \f$, etc. We allow for up to 3 mesh adaptations 
for each Reynolds number to allow the mesh to adapt itself to the flow
conditions before starting the time-dependent run.

\until reached final Reynolds number


Finally, we start the proper timestepping procedure, allowing for
one mesh adaptation per timestep and suppressing the re-assignment
of the initial condition after the mesh adaptation by setting
the \c first flag to false.
 
\until }//end of main

<HR>
<HR>


\section problem The Problem class

The \c Problem class contains the usual member functions, most of which
are either empty or explained in more detail below. 

\dontinclude fsi_channel_with_leaflet.cc
\skipline start_of_problem_class
\until doc_solution


Two member functions are implemented here: The function 
\c actions_before_implicit_timestep() updates the time-dependent
velocity profile at the upstream boundary (boundary 3; see the
enumeration of the mesh boundaries in the
<A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
<CODE>ChannelWithLeafletMesh</CODE></A>):
\until   } // end of actions_before_implicit_timestep

Since the node-update is performed by an algebraic node-update
procedure the nodal positions in the fluid mesh must be 
updated whenever any of the (solid) degrees-of-freedom change.
This is done automatically during the computation of the shape
derivatives (implemented in the \c AlgebraicElement wrapper class
described <A HREF="../../fsi_collapsible_channel_algebraic/html/index.html">
another tutorial</A>), but an additional node update must be performed
when the unknowns are updated by the Newton solver. This
is achieved by performing a further node-update in the function 
\c actions_before_newton_convergence_check():

\until }

The private member data contains the usual pointers to the \c
Problem's two sub-meshes,
the pointer to \c GeomObject that represents the leaflet, and
the total height of the channel.

\until };

<HR>
<HR>

\section constructor The problem constructor

We start the problem construction by creating two timesteppers:
A \c BFD<2> timestepper for the fluid and the fake timestepper
\c Steady<2> for the (massless) solid. 
(<A HREF="../../fsi_collapsible_channel_adapt/html/index.html">Recall 
</A> that timesteppers
from the \c Steady family return zero time-derivatives but keep 
track of the past histories. These are needed during the adaptive
refinement of the fluid mesh: when assigning the history of the
previous nodal positions for newly created fluid nodes, we must
evaluate the position of the leaflet at previous timesteps; this is
discussed in more detail in 
<A HREF="../../fsi_collapsible_channel_adapt/html/index.html">
another tutorial</A>.)

\skipline start_of_constructor
\until add_time_stepper_pt(wall_time_stepper_pt);

Next we create the mesh of 1D Hermite beam elements that represents
the leaflet, passing a pointer to an instance of the \c
UndeformedLeaflet to its constructor.

\until ,undeformed_wall_pt,wall_time_stepper_pt);

The wall mesh defines the geometry of the deformed leaflet, therefore
we create a \c GeomObject representation of that mesh, using the
\c MeshAsGeomObject class, and pass it to the constructor of the fluid mesh:

\until  fluid_time_stepper_pt);

We specify an error estimator for the fluid mesh
and add both meshes to the \c Problem before combining them
into a global mesh.

\until  build_global_mesh();

Dirichlet conditions (prescribed velocity) are applied on virtually
all boundaries of the fluid domain (prescribed inflow at the inlet;
zero velocity on the rigid channel walls; fluid velocity prescribed
by the wall motion on the leaflet), apart from the outflow
(boundary 1; see the enumeration of the mesh boundaries in the
<A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
<CODE>ChannelWithLeafletMesh</CODE></A>) where the axial velocity 
is unknown (and determined indirectly
by the "axially-traction-free" condition) while the vertical
velocity has to be set to zero to ensure parallel outflow.

\until end loop over boundaries

Next we assign the parabolic velocity profile at the inlet
(recall that all values are initialised to zero so no further action
is required on any of the other Dirichlet boundaries). 

\until  }// end of setup boundary condition

The leaflet is clamped at its lower end so we pin its \f$x_1\f$- and
\f$x_2\f$-positions, and impose a vertical slope by setting \f$ dx_1/ds = 0 \f$
(where s is the local coordinate along the element; see the
discussion of the boundary conditions for beam elements in 
<A HREF="../../../beam/tensioned_string/html/index.html">another
tutorial</A> for details).

\until pin_position(1,0)

Next, we complete the build of the fluid elements by passing 
pointers to the relevant physical parameters to the elements
and pinning any redundant pressure degrees of freedom; see
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details.

\until pin_redundant

 When completing the build of the wall elements we use the function
\c enable_fluid_loading_on_both_sides() to indicate that
the leaflet is completely immersed in the fluid so that
fluid tractions act on both of its faces. When setting up the
fluid-structure interaction below, one of the two faces
will have to be identified as the "front" (face 0) and the other
one as the "back" (face 1). The function \c normal_points_into_fluid()
allows us to indicate if the outer unit normal on the leaflet 
(as computed by \c FSIHermiteBeamElement::get_normal(...) points
into the fluid, when viewed from the "front" face. Here it does not --
see \ref comm for more details on this slightly subtle point). 


\until // end of loop over elements


We can now set up the fluid structure interaction. The motion of the
leaflet determines the fluid velocity of all nodes on boundaries 4 and 5
via the no-slip condition (see the enumeration of the mesh boundaries in the
<A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
<CODE>ChannelWithLeafletMesh</CODE></A>). The fluid velocity at these nodes can be updated
automatically whenever a fluid node is moved (by the fluid mesh's
node-update function) by specifying an auxiliary node update function.

\until aux node update fct has been set

Finally, we have to determine which fluid elements are adjacent
to the two faces of the \c FSIHermiteBeamElements 
to allow them to compute the fluid traction they are exposed to.
This is done separately for the "front" and "back" faces. The
"front" of the leaflet (face 0) is assumed to coincide with the fluid mesh
boundary 4; the "back" (face 1) is assumed to coincide with the fluid mesh
boundary 5.

\until   (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 


Following the setup of the fluid-structure interaction we
assign the equation numbers:

\skipline   // Setup equation numbering scheme
\until cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
\skipline end of constructor

<HR>
<HR>

\section adapt Actions after the mesh adaptation

Once the fluid mesh has been adapted we free all pressure
degrees of freedom and then (re-)pin any redundant ones;
see the discussion in 
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</a> for details.

\skipline start_of_actions_after_adapt
\until pin_redundant

Any newly created fluid nodes on the leaflet (i.e. on the
fluid mesh boundaries 4 and 5) must be subjected to the automatic
application of the no-slip condition. For simplicity 
we (re-)specify the auxiliary node update function pointer
for all of the nodes on those boundaries.

\until  aux node update fct has been (re-)set

Finally, the adaptation of the fluid mesh may change the fluid 
elements that are adjacent to the wall elements so we re-generate the
corresponding FSI lookup schemes.

\until end_of_actions_after_adapt


<HR>
<HR>

\section doc Post Processing

The function \c doc_solution(...) documents the results. We output
the fluid and solid meshes and document selected additional quantities
in the trace file.

\skipline start_of_doc_solution
\until << std::endl;

The remaining output is created to determine which elements
are located next to the fluid mesh boundaries 4 and 5
(Yes, <a href="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
the documentation for the mesh</A> illustrates the enumeration of the mesh
boundaries but we generally prefer to be paranoid; see \ref comm ), 
and to establish the direction of the unit normal on the beam. 


\until end_of_doc_solution


<HR>
<HR>

\section comm_ex Further comments and exercises

\subsection comm Further comments

-# <B>How does one identify the "front" and "back" of an immersed beam
   structure?</B> \n\n
   When setting up the fluid-structure interaction in the \ref constructor
   we had to associate the fluid mesh boundaries 4 and 5 (the
   boundaries adjacent to the leaflet (see the enumeration of 
   the mesh boundaries in the
   <A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
   <CODE>ChannelWithLeafletMesh</CODE></A>) with the "front" and "back"
   of the fully immersed beam structure. Furthermore, since the
   computation of the fluid traction acting on the leaflet depends on the
   direction of the normal, and the normal on the "front" of the leaflet
   points in the opposite direction to that on its "back" it is
   important to assess which normal is used in the automatic 
   computation of the fluid traction.
   There are two ways of establishing this: \n\n 
   -# <B>Quick and dirty:</B> \n Set 
      \c FSIHermiteBeamElement::normal_points_into_fluid() to \c true
      and perform a steady computation. If the leaflet is sucked
      towards the high-pressure region you got it wrong. \n\n
   -# <B>Do it properly:</B> \n Use the output generated in 
      \c doc_solution(...) to identify the elements adjacent
      to the fluid mesh boundaries 4 and 5 (i.e. the fluid elements next to the
      leaflet's "front" and "back" faces), respectively, and the normal vector
      returned by \c FSIHermiteBeamElement::get_normal(...) -- this
      is the normal that is used in the computation of the fluid
      traction. The figure below shows a plot of these
      within an adaptively refined fluid mesh. \n\n
\image html front_and_back.gif "Elements near the `front' (red) and `back' (green) of the leaflet and the outer unit normal on the leaflet. " 
\image latex front_and_back.eps "Elements near the `front' (red) and `back' (green) of the leaflet and the outer unit normal on the leaflet. " width=0.75\textwidth
      The red and green elements are the fluid elements adjacent to
      the fluid mesh boundaries 4 and 5 (i.e. the "front" and "back" 
      of the leaflet), respectively. The normal
      vectors (shown in blue) point towards the upstream end of the channel so 
      they point \e away \e from (rather than \e into) the fluid
      that is adjacent to the "front" (the red region). This is why we set
      \c FSIHermiteBeamElement::normal_points_into_fluid()=false
      (Admittedly, we used the former method to get this right and only
      documented the proper way to do it when writing this tutorial...)
      \n\n
-# <B>Use of faster solvers: \c oomph-lib's FSI preconditioner</B> \n\n
   The problem presented in this tutorial was used as one of the test cases
   for \c oomph-lib's FSI preconditioner; see \n\n
   <CENTER>
   <a href="http://www.maths.man.ac.uk/~mheil/oomph_lib_additional_material/HeilHazelBoyleCompMech/HeilHazelBoyleCompMech.pdf">
   Heil, M., Hazel, A.L. & Boyle, J. (2008): Solvers for
   large-displacement fluid-structure interaction problems: 
   Segregated vs. monolithic approaches. Computational Mechanics.</a>
   </CENTER>
   \n\n
   The use of this preconditioner (which leads to much faster solution
   times) is described in another tutorial. However, the required source
   code is already part of the driver code and has simply been ignored
   in this tutorial. Feel free to experiment with the preconditioner 
   by modifying
   <A HREF="../../../../demo_drivers/interaction/fsi_channel_with_leaflet/fsi_channel_with_leaflet.cc">the source code.</A>
   
\subsection ex Exercises

-# Confirm that choosing the wrong value for the flag  
   \c FSIHermiteBeamElement::normal_points_into_fluid() makes
   the leaflet move in the wrong direction as it perceives positive
   fluid tractions as negative ones and vice versa. Here is 
   what you should get:
\image html normal_points_in_wrong_direction.gif "Deformation of the leaflet with the wrong choice of FSIHermiteBeamElement::normal_points_into_fluid(): The leaflet gets sucked towards the high-pressure region. " 
\image latex normal_points_in_wrong_direction.eps "Deformation of the leaflet with the wrong choice of FSIHermiteBeamElement::normal_points_into_fluid(): The leaflet gets sucked towards the high-pressure region. " width=0.75\textwidth
   Clearly, the normal points in the wrong direction and the
   tractions have the wrong sign. This can be made "consistent",
   however, by also swapping the association between the mesh boundaries
   and the leaflet's faces. I.e. if we associate mesh boundary 5
   with the "front" via  \n\n
   \code
   // Front of leaflet: face 0
   face=0; 
   FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
    (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 
   \endcode
   \n and mesh boundary 4 with the "back" via \n\n
    \code
   // Back of leaflet: face 1
   face=1; 
   FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
    (this,4,Fluid_mesh_pt,Wall_mesh_pt,face); 
   \endcode
   \n the code gives the right results again (make sure you change 
   the code in the problem constructor and in the actions
   after adaptation!) -- sometimes two wrongs do
   cancel each other out! 
   \n\n 
-# <A HREF="../figures/fsi_channel_with_leaflet_flow.avi">
   The animation</A> of the flow field shows that the
   outflow boundary condition is applied "too close"
   to the leaflet: in the relatively short domain chosen here,
   the outflow is far from fully developed and
   as a result the imposition of a parallel flow leads to the
   development of a artificial outflow boundary layer. \n\n
\image html outflow_bl.gif "Artificial outflow boundary layer. " 
\image latex outflow_bl.eps "Artificial outflow boundary layer. " width=0.75\textwidth
   Confirm 
   that an increase in the downstream length of the domain (or a reduction
   in the Reynolds number) improves the situation. 
   \n\n 
-# The algebraic node-update technique employed in the 
   <A HREF="../../../meshes/mesh_list/html/index.html#channel_with_leaflet">
   <CODE>ChannelWithLeafletMesh</CODE></A> forces \e all fluid nodes
   to move when the leaflet deforms. An increase in the downstream
   length of the fluid domain therefore leads to a significant
   increase in computational cost. Use the improved node-update
   technique suggested in the 
   <A HREF="../../../navier_stokes/channel_with_leaflet/html/index.html#ex">
   corresponding Navier-Stokes example</A> (only move the nodes
   in the vicinity of the leaflet) to improve the
   efficiency of the code. 
.
 
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
demo_drivers/interaction/fsi_channel_with_leaflet/
">
demo_drivers/interaction/fsi_channel_with_leaflet/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/interaction/fsi_channel_with_leaflet/fsi_channel_with_leaflet.cc
">
demo_drivers/interaction/fsi_channel_with_leaflet/fsi_channel_with_leaflet.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


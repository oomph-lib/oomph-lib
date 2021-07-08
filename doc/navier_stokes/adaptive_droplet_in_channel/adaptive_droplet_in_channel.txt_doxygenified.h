/**

\mainpage Propagation of a droplet in a channel - mesh generation and adaptation for free surface problems

In this tutorial we demonstrate another adaptive solution of free surface
problems on unstructured meshes, using the example of a droplet
propagating along a straight channel. The problem is extremely similar
to <a href="../../adaptive_bubble_in_channel/html/index.html"> the
propagation of a bubble in a channel tutorial </a>. Thus, we shall
only discuss the differences from that tutorial. The key physical
difference is that instead of the uniform pressure state in an
inviscid bubble, the droplet
consists of a viscous fluid that can support internal stress variations.


<HR> <HR> 

\section example The example problem 
We illustrate the solution
of the unsteady 2D Navier-Stokes equations by considering the
propagation of a single droplet along a straight channel as shown in
the sketch below. The non-dimensionalisation is
the same as in the 
<a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
tutorial, </a> and we choose the viscosity and density of the
surrounding liquid to be the reference values.

\image html problem.gif "The problem setup. " 
\image latex problem.eps "The problem setup. " width=0.9\textwidth

The governing equations and boundary conditions are 
the same as those in the  
<a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
tutorial. </a>
The only difference is that the Navier--Stokes equations must also be
solved within the droplet. In fact, this is a two-fluid problem, a
class of problems that is first introduced in 
<A HREF="../../two_layer_interface/html/index.html"> another tutorial. </A>

<HR>

The constraint that the droplet
volume remains constant must be enforced explicitly in the static
case, as in the bubble problem. In the time-simulations, however, a
constant drop volume is implicitly enforced by the continuity
equation. For the bubble, the continuity equation is
not solved within the interior, which is why the volume constraint
must always be explicitly enforced in that case.

<HR> <HR>

\section implementation Implementation

We use the same method as in the 
<a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
problem, </a> an ALE-based finite-element method with
<a href="../../single_layer_free_surface/html/index.html#kinematic_condition_implementation">a pseudo-elastic
node-update procedure</a>.  In this case,
there is a pressure jump across the interface between the
fluids, which means that we use triangular Crouzeix--Raviart 
elements rather than the continuous-pressure Taylor--Hood
elements. Again, we impose the kinematic and dynamic boundary conditions
with \c FaceElements. The volume constraint for the static problem
is also imposed in a similar
way: we attach \c LineVolumeConstraintBoundingSolidElements to the 
droplet surface and create an additional \c VolumeConstraintElement. 
Here, a pressure degree of
freedom within the droplet is hijacked, see 
<A HREF="../../static_two_layer/html/index.html"> another tutorial, </A> 
to be used as the unknown associated with the volume constraint. Once
the static initial problem has been solved, the
"volume constraint" elements are deleted and the pressure degree of
freedom is unhijacked.

<HR> <HR>

\section results Results
We perform the simulation in a two-stage procedure. We start by 
performing a steady solve with the inflow switched off. This
deforms the droplet into its steady state (approximately) circular
configuration with the required volume. The actual time-dependent 
simulation is then performed with an impulsive
start from this configuration. 

The figure below shows the location of the droplet and mesh (upper figure) and 
a contour plot of the pressure
distribution with overlaid velocity vectors of the difference between
the background Poiseuille flow and the velocity field (lower
figure). The figure is a snapshot for the parameters \f$ Re = ReSt =
0.0 \f$, \f$ Ca = 10.0 \f$ and a droplet that is ten times as viscous
as the surrounding liquid.

\image html drop_snapshot.gif "Snapshot of the flow field (velocity vectors and pressure contours) for a propagating droplet. " 
\image latex drop_snapshot.eps "Snapshot of the flow field (velocity vectors and pressure contours) for a propagating droplet. " width=0.7\textwidth



<HR> <HR>


\section parameters Global parameters

The namespace containing the dimensionless parameters contains an
additional viscosity ratio parameter, compared to that in the 
 <a href="../../adaptive_bubble_in_channel/html/index.html">
bubble problem </a>.

\dontinclude adaptive_drop_in_channel.cc
\skipline Viscosity ratio
\until 10.0


<HR> <HR>

\section main The driver code 

The first difference from the bubble problem is that the steady solver
does not converge when \f$ Ca = 10 \f$. Instead, we start with \f$
Ca = 1 \f$, solve the steady problem, set \f$ Ca=10\f$ and then resolve.

\skipline problem.steady_newton_solve(1)
\until problem.steady_newton_solve()

After documenting the solution, we remove the (explicit) volume
constraint and then the remainder of the code is identical to that in
the <a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
tutorial. </a>

\skipline Switch off
\until volume_constraint()
 
<HR> <HR>

\section problem The problem class

Other than trivial name changes from "bubble" to "drop", there are
a few significant changes between this problem and that in 
the  
<a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
tutorial. </a> There is an additional boolean \c Use_volume_constraint
and an additional function \c
remove_volume_constraint(), which are used to manage the switch from
the explicit enforcement of the volume constraint in the static case
to the implicit enforcement in the time simulations. The other key
difference is that we use \c Triangle's 
<a href="../../../meshes/mesh_from_inline_triangle_internal_boundaries/html/index.html#def_extra_regions"> region attributes </a> 
to distinguish the elements inside
the droplet from those outside. The default behaviour is that all
elements are in region 0, but we label those element within the drop
as region 1.

<HR> <HR>

\section constructor The problem constructor 

 The construction of the mesh proceeds in exactly the same way as in
 <a href="../../adaptive_bubble_in_channel/html/index.html"> bubble
tutorial, </a> except that we add a region tag "1" to label the
  elements within the droplet (so we specify a coordinate within the drop)
  and we must tell \c Triangle to use the assigned attributes.
\dontinclude adaptive_drop_in_channel.cc
\skipline Define the region
\until drop_center)
 The remainder of the constructor is the same as the other tutorial.

<HR> <HR>

\section problem_setup Problem setup 

 When the bulk elements are made fully functional, we add the pointer to
 the viscosity ratio to all elements in the drop (region 1).

\skipline within the droplet 
\until } 

<HR> <HR>

\section face_elements Generation of face elements 

As usual we impose
the kinematic and dynamic boundary condition at the interface by
attaching \c FaceElements to the relevant boundaries of the bulk
elements. However, we must be careful to add only a \b single layer of
elements. If we use the standard "boundary element" functions then we
will be creating face elements on both sides of the internal
boundary. Instead, we use the elements adjacent to the boundary
within region 0, which ensures that only a single layer of interface
elements are added.
 
\dontinclude adaptive_drop_in_channel.cc
\skipline start_of_create_free_surface_elements 
\until end of create_free_surface_elements
 
The volume constraint elements are only created if the boolean flag \c
Use_volume_constraint is true (the default on construction of the problem).
We hijack the pressure degree of freedom associated with the first
element in region 1 and then the construction of the elements again
uses the regions to ensure that a single layer of elements is created.
 
\skipline start_of_create_volume_constraint_elements
\until end of create_volume_constraint_elements

<HR> <HR>

\section remove Removal of the volume constraint

The function \c remove_volume_constraint(), resets the boolean flag to
false, clears the hijacked data, deletes the volume constraint
elements and mesh and then removes the volume constraint mesh from the
problem's list of sub meshes, before reassigning the equation numbers.

\dontinclude adaptive_drop_in_channel.cc
\skipline Change the boundary conditions
\until }


<HR> <HR>

\section comments Comments and Exercises


 The computation of the initial static solution means that we must
 still include the calls to the \c create_volume_constraint_elements() and
 \c delete_volume_constraint_elements() in \c actions_before_adapt()
 and \c actions_after_adapt(). We have chosen to have the functions
 return immediately if the volume constraint is not being enforced,
 rather than using \c if blocks within the \c
 actions_before/after_adapt() functions.
 
<HR>

\subsection exercises Exercises

-# Explore what happens as the viscosity ratio is varied. Are the
   results as you expect? Does the solution tend to a steadily
   propagating state? Does the solution approach the bubble solution
   as this viscosity ratio tends to zero? (Is this a sensible limit to take?)
   \n\n
-# How could you modify the code to compute the steadily-propagating
   solutions directly, rather than using time simulation?
   \n\n
.

<HR> <HR>


\section sources Source files for this tutorial 
- The source files for
  this tutorial are located in the directory:\n\n 
  <CENTER> 
  <A HREF="../../../../demo_drivers/navier_stokes/unstructured_adaptive_fs/">
  demo_drivers/navier_stokes/unstructured_adaptive_fs/ </A> 
  </CENTER>\n
- The driver code is: \n\n 
  <CENTER> 
  <A HREF="../../../../demo_drivers/navier_stokes/unstructured_adaptive_fs/adaptive_drop_in_channel.cc">
   demo_drivers/navier_stokes/unstructured_adaptive_fs/adaptive_drop_in_channel.cc</A> 
  </CENTER> .


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


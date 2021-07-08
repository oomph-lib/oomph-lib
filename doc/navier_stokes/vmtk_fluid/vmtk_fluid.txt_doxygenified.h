/**

\mainpage Demo problem: Steady finite-Reynolds-number flow through an iliac bifurcation

The purpose of this tutorial is to demonstrate the simulation of
cardiovascular fluid mechanics problems with \c oomph-lib. We showed in
<a href="../../../meshes/mesh_from_vmtk/html/index.html">another tutorial</a> 
how to use <a href="http://www.vmtk.org">VMTK</a>
together with \c oomph-lib's conversion code \n\n
<CENTER>
<a href="../../../../demo_drivers/meshing/mesh_from_vmtk/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc">demo_drivers/meshing/mesh_from_vmtk/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc</a>
</CENTER>
\n\n
to generate unstructured fluid and solid meshes for the simulation 
of physiological fluid-structure interaction problems based on data from 
medical images. Here we show how to simulate steady finite-Reynolds-number flow
through a (rigid) iliac bifurcation. A particular feature of this
problem is that, unlike the problem considered in
<a href="../../unstructured_three_d_fluid/html/index.html">
another tutorial,</a> the in- and outflow boundaries are not
aligned with any coordinate planes. Parallel in- and outflow
is therefore enforced by a Lagrange multiplier method, implemented
using \c oomph-lib's \c FaceElements.


The problem studied here also serves as a "warm-up problem" for the 
<a href="../../../interaction/vmtk_fsi/html/index.html">
corresponding fluid-structure interaction problem </a> in which the vessel
wall is elastic and deforms in response to the traction 
that the fluid exerts onto it.
 
We stress that the tutorial focuses on the implementation 
aspects, not the actual physics. Since the driver code discussed
here is also used in the library's self-tests we
deliberately use a very coarse mesh and restrict ourselves to
steady flows.  The results shown below are therefore unlikely 
to bear much resemblance to the actual flows that arise
\e in \e vivo. The section \ref realism at the end of this
tutorial provides several suggestions on how to make the simulation
more realistic.

<HR>
<HR>

\section problem The problem (and results)
The two figures below show the geometry of the blood vessel 
(obtained from a scan of an iliac bifurcation, using the procedure 
discussed in \c oomph-lib's 
<a href="../../../meshes/mesh_from_vmtk/html/index.html">VMTK tutorial</a>) 
and the flow field (velocity vectors and pressure contours) for a 
nominal Reynolds number of \f$ Re = 10. \f$  (See \ref re for a more
detailed discussion of the Reynolds number.) The flow is driven by an
applied pressure drop between the in- and outflow boundaries,
as in the <a href="../../unstructured_three_d_fluid/html/index.html">
previous example.</a>

\image html vmtk_fluid.gif "Velocity field and pressure contours. " 
\image latex vmtk_fluid.eps "Velocity field and pressure contours. " width=0.75\textwidth

The figure below shows more clearly that the in- and outflow from the 
upstream and downstream cross-sections is parallel, even though
the cross-sections are not aligned with any the coordinate planes. 

\image html vmtk_fluid2.gif "Velocity field and pressure contours. Note the parallel in- and outflow. " 
\image latex vmtk_fluid2.eps "Velocity field and pressure contours. Note the parallel in- and outflow. " width=0.75\textwidth


<HR>
<HR> 

\section parallel_flow Imposing parallel in- and outflow
In most of the Navier-Stokes problems we have considered so far, 
the geometry of the fluid domain was such that the in- and outflow 
cross-sections were aligned with the Cartesian coordinate planes. 
In such geometries the imposition of parallel in- and outflow is 
straightforward as it only requires pinning of the transverse 
velocity components.  A pressure drop between 
upstream and downstream ends can be applied by attaching 
\c NavierStokesTractionElements to the bulk elements that are 
adjacent to the relevant domain boundaries.

In the current problem, the in- and outflow cross-sections have 
some arbitrary orientation, implying that when the in- or outflow is 
parallel (or, more accurately, orthogonal to the in- or outflow cross
sections), neither of the three velocity components vanishes. 
The easiest way to enforce parallel outflow in such situations is
to employ Lagrange multipliers to enforce the two constraints
\f[     
{\bf u} \cdot {\bf t}_\alpha = 0 \mbox{\ \ \ \ for $\alpha = 1,2$,}
\f] 
where \f$ {\bf t}_\alpha \f$ (for \f$\alpha = 1,2\f$) are the
two tangent vectors spanning the in- or outflow cross-sections. 
Physically, the Lagrange multipliers act as tangential tractions 
that enforce the parallel flow. The Lagrange multipliers 
introduce additional degrees of freedom into the problem and their 
implementation as \c ImposeParallelOutflowElements employs a 
technique similar to that used to enforce prescribed boundary 
displacements in solid mechanics problems. (This is discussed in 
<a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
another tutorial</a>.) The \c ImposeParallelOutflowElements also
allow the specification of a pressure that acts on the fluid
(in the direction opposite to the outer unit normal on the bulk
fluid element). 


<HR>
<HR>

\section driver_code The driver code
The driver code is very similar to that discussed in 
<a href="../../unstructured_three_d_fluid/html/index.html">
another example</a> where we imposed parallel outflow
in coordinate-aligned in- and outflow cross-sections.
We will therefore only discuss the changes required 
to impose parallel flow in cross-sections with arbitrary
orientation.

<HR>


\subsection namespace Problem parameters
The namespace specifying the problem parameters contains the Reynolds 
number and the in- and outflow pressures (rather than the complete traction
vector):

\dontinclude vmtk_fluid.cc
\skipline start_namespace
\until end namespace

<HR>


\subsection main The main function
The main function remains identical to that in the 
<a href="../../unstructured_three_d_fluid/html/index.html">
problem with axis-aligned outflow.</a>

<HR>


\section class The Problem class
The \c Problem class is practically identical to that in the 
<a href="../../unstructured_three_d_fluid/html/index.html">
problem with axis-aligned outflow</a>, apart from some trivial changes
such as the provision of storage for meshes of 
\c ImposeParallelOutflowElements rather than
\c NavierStokesTractionElements.

<HR>

\section constructor The Problem constructor
The problem constructor is also very similar. We start by building 
the fluid mesh, using the files created by 
<A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A>:

\skipline start_constructor
\until split_corner_elements);

 
Next, we set up a boundary lookup scheme that records which elements
are located next to which domain boundaries, and specify the 
IDs of the mesh boundaries that coincide with the in- and outflow
cross-sections. Note that this information reflects the
specification of the boundary IDs in the \c tetgen <code>*.poly</code>
file. [The conversion code 
<a href="../../../../demo_drivers/meshing/mesh_from_vmtk/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc">create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc</a>
lists the relation between the original boundary IDs and the
new ones (obtained by giving each surface facet a separate 
boundary ID) at the end of the <code>*.poly</code> file.]

\until // done outflow boundaries

We create the meshes containing the Lagrange multiplier elements
and add all sub-meshes to the Problem's global mesh.

\until build_global_mesh();


Next we apply the boundary conditions. We start by identifying
the IDs of the boundaries that are subject to no-slip boundary 
conditions.

\until done identification of boundaries where velocities are pinned

Next we loop over all boundaries, visit their nodes and pin all three velocity
components if the boundary is subject to a no-slip
condition:

\until pin(2);

We now check if the node in question is also located on the
in- and outflow boundaries...

\until now we know


...and if it is, we pin the 
Lagrange multipliers. They are stored after
the values allocated by the "bulk" elements and we obtain
the index of the first value associated with the Lagrange
multipliers from the \c BoundaryNodeBase::index_of_first_value_assigned_by_face_element() 
function.

\until end of BC

The rest of the constructor is unchanged. We pass the pointer to 
the Reynolds number to the elements and assign the equation numbers:

\until end constructor

<HR>

\subsection lagr Creating the Lagrange multiplier elements

The helper function \c create_parallel_outflow_lagrange_elements()
loops over the bulk elements that are adjacent to the in- and outflow
cross-sections and attaches \c ImposeParallelOutflowElements
to the relevant faces. We store pointers to the newly-created
elements in the appropriate meshes, and pass pointers to the 
doubles that specify the imposed pressure to the elements.


\skipline start_of_lagrange_multiplier_elements
\until done

 
<HR>
<HR>


\subsection doc Post-processing

The post-processing routine is unchanged.

<HR>

\section comm_ex Comments and Exercises

\subsection re What does the Reynolds number mean in this problem?
\c oomph-lib's implementation of the Navier-Stokes equations is 
based on <a href="../../driven_cavity/html/index.html#equation">their 
non-dimensional form</a> so we typically expect
the geometry of the problem to have been non-dimensionalised on
a representative lengthscale, \f$ {\cal L} \f$. When dealing with
geometries that are obtained from medical images, the vessel coordinates are
typically provided as dimensional quantities, e.g. in millimetres. 
Discarding the unit of the coordinates (i.e. using
a non-dimensional coordinate \f$ (x,y,z)=(1,2,3)\f$ to represent 
the point located at \f$
(x^*,y^*,z^*)=(1\mbox{mm},2\mbox{mm},3\mbox{mm}),\f$ say) 
is therefore equivalent to non-dimensionalising
all lengths on a reference length of \f$ {\cal L} = 1\mbox{mm}. \f$
<a href="../../driven_cavity/html/index.html#equation">Recall</a> 
that the Reynolds number is defined as 
\f[
Re = \frac{\rho {\cal U} {\cal L}}{\mu}
\f]
where \f$ \rho \f$ and \f$ \mu \f$ are the fluid density and viscosity,
respectively. The lengthscale \f$ {\cal L} = 1\mbox{mm} \f$ does obviously
not provide a measure of the dimension of our blood vessel whose
diameter (at the inlet) is about \f$ D \approx 16\mbox{mm}. \f$ The 
nominal Reynolds number of \f$ Re=10 \f$ used in the computations 
is therefore equivalent to an actual Reynolds number (formed with 
the vessel diameter)
of 
\f[
Re_D = \frac{\rho {\cal U} D }{\mu} = Re \ \frac{D}{\cal L} = 160,
\f]
where the velocity scale \f$ {\cal U} \f$ is formed with the 
(dimensional) applied pressure drop \f$ \Delta P^* \f$
between the in- and outflow cross-sections, as in the 
<a href="../../unstructured_three_d_fluid/html/index.html">problem
with axis-aligned in- and outflow cross-sections,</a>
\f[
{\cal U} = \frac{\Delta P^*\  {\cal L} }{\mu}.
\f] 
<HR> 
 

\subsection realism How to make the simulation more realistic
The simulation presented above is obviously very crude and serves
primarily as a proof of concept. However, it is straightforward
to address most of the shortcomings and we encourage you to explore
the following improvements:
- Generate a finer fluid mesh using the instructions 
  in \c oomph-lib's 
  <a href="../../../meshes/mesh_from_vmtk/html/index.html">VMTK
  tutorial </a> We are hoping to make the image files
  used in this tutorial available soon. Please 
  <a href="../../../contact/html/index.html">contact us</a>
  if you can't wait. Remember that you will have to update the 
  enumeration of the domain boundaries if you change the mesh.
  \n\n
- Attach "flow extensions" to the in- and outflow
  cross-sections, using the technique described in \c oomph-lib's
  <a href="../../../meshes/mesh_from_vmtk/html/index.html#add_extensions">
  VMTK tutorial.</a>
  \n\n
- Make the problem time-dependent and apply a period driving pressure
  drop.
  \n\n
.
<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/vmtk_fluid/">
demo_drivers/navier_stokes/vmtk_fluid/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/vmtk_fluid/vmtk_fluid.cc">
demo_drivers/navier_stokes/vmtk_fluid/vmtk_fluid.cc
</A>
</CENTER>
.











<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


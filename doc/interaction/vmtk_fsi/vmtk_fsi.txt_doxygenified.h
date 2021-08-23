/**

\mainpage Demo problem: Steady finite-Reynolds-number flow through an elastic iliac bifurcation   
This tutorial demonstrates how to simulate physiological
fluid-structure interaction problems, based on medical imaging
data, processed with the <a href="http://www.vmtk.org">Vascular 
Modeling Toolkit.</a> 

 We combine two single-physics problems, namely
- <a href="../../../solid/vmtk_solid/html/index.html">
  The inflation of a blood vessel</a>
- <a href="../../../navier_stokes/vmtk_fluid/html/index.html">
  Finite-Reynolds number flow through a rigid iliac bifurcation</a>
.
to study the steady flow through an elastic iliac bifurcation.
(Please refer to \c oomph-lib's 
<a href="../../../meshes/mesh_from_vmtk/html/index.html">VMTK
tutorial</a> to find out how to generate \c oomph-lib meshes
from medical images.)

The tutorial also provides a detailed discussion of the mechanism by which
\c FaceElements introduce additional unknowns into a problem. 
This is important in the problem considered here because we
employ two different types of \c FaceElements, both of which
introduce additional unknowns. When dealing with nodes that
are attached to both types of \c FaceElements we must therefore be able 
to determine which nodal values are associated with which \c FaceElement.
The section \ref face below shows how this is done.


<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
</CENTER>
This tutorial and the associated driver codes were developed jointly
with Amine Massit (ENSTA, Paris).
</TD>
</TR>
</TABLE>
</CENTER>


<hr>
<hr>
  

\section prob The problem (and some results)
The animations below show the deformation of an iliac bifurcation
that conveys viscous fluid and deforms in response to the traction 
that the fluid exerts onto it. As in the previous examples
where we considered the constituent single-physics 
<a href="../../../navier_stokes/vmtk_fluid/html/index.html">
fluid</a> and <a href="../../../solid/vmtk_solid/html/index.html">
solid</a> problems, the meshes are obviously very coarse and the
boundary conditions are far from "physiologically realistic":
We keep the upstream and downstream "ends" of the vessel
wall in a fixed position and drive the (steady!) flow by applying
a constant pressure drop between the in- and outflow cross sections.
The different frames shown in the animation show results
for different wall stiffnesses, using the same setup as in
<a href="../../unstructured_three_d_fsi/html/index.html">another
tutorial.</a> 

We stress that the computations presented here serve as
quick-to-compute proof of concept and refer to the section 
\ref comm_ex  at the end of this tutorial for a discussion 
on how to make the simulation more "realistic". 

\image html elastic_iliac_solid.gif "The flow-induced deformation of an iliac bifurcation. The animation shows the effect of variations in the wall stiffness. The arrows show the magnitude of the fluid traction on the wall. " 
\image latex elastic_iliac_solid.eps "The flow-induced deformation of an iliac bifurcation. The animation shows the effect of variations in the wall stiffness. The arrows show the magnitude of the fluid traction on the wall. " width=0.75\textwidth

\image html elastic_iliac_fluid.gif "Steady finite-Reynolds number flow through an elastic iliac bifurcation (velocity vectors and pressure contours). The animation shows the effect of variations in the wall stiffness. " 
\image latex elastic_iliac_fluid.eps "Steady finite-Reynolds number flow through an elastic iliac bifurcation (velocity vectors and pressure contours). The animation shows the effect of variations in the wall stiffness. " width=0.75\textwidth


The driver code used for this computation is very similar to the one
discussed in <a href="../../unstructured_three_d_fsi/html/index.html">another
tutorial</a> where we used a much simpler geometry in which the
three branches of the bifurcation had rectangular cross-sections. 
Before discussing the changes to the driver code in section
\ref driver_code below, we provide a more detailed discussion
of the way in which \c FaceElements introduce additional unknowns into 
a problem. Feel free to skip the next section if you're not keen on theory.

 

<HR>
<HR>

\section face How FaceElements introduce additional unknowns into a problem

\c FaceElements are used widely throughout \c oomph-lib to apply 
Neumann/flux/traction-type boundary conditions on the faces of
higher-dimensional "bulk" elements. Examples include:
- <a href="../../../poisson/two_d_poisson_flux_bc2/html/index.html">
  the application of a flux boundary condition in a Poisson problem;</a> 
- the application of traction boundary conditions in 
  <a href="../../../navier_stokes/rayleigh_traction_channel/html/index.html">
  fluid</a> and 
  <a href="../../../solid/airy_cantilever/html/index.html">solid</a>
  mechanics problems;
- <a href="../../../young_laplace/contact_angle/html/index.html">the
  application of a contact-angle boundary condition for the
  Young-Laplace equation.</a>
.
In all the examples listed above, the boundary conditions simply
add a contribution to the elements' residuals but they do not
introduce any additional unknowns into the problem.

\c FaceElements may also be used to apply boundary conditions
via Lagrange multipliers. An example is given in the tutorial
that discusses
- <a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
   the application of displacement boundary conditions for solid
   mechanics problems</a>
.
In such problems, the Lagrange multipliers must be determined
as part of the solution, and storage for the associated discrete
unknowns is created at the nodes of the \c FaceElements.

 To explain the relevant details of the implementation we consider a 
simple 2D Navier-Stokes problem discretised using nine-node Taylor-Hood
elements (in these elements each vertex node stores two discrete velocities and
one pressure; the other nodes store only two velocity degrees of
freedom). We assume that boundaries 0 and 1 are subject to 
boundary conditions imposed via \c FaceElements, and that
each boundary condition introduces its own Lagrange multipliers field<b>s</b>. 
[Yes, the plural is correct. As an example, consider the case of 
imposing displacement constraints in a 2D solid mechanics problem via 
Lagrange multipliers. In this approach 
the imposition of the boundary condition requires \e two Lagrange 
multipliers along each constrained boundary. Physically, the Lagrange
multipliers represent the two components of the surface traction required to 
deform the boundary into the required shape; see 
<a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
the relevant solid mechanics problem</a> for details.]

The sketch below shows the discretisation of the domain, with the black 
circles representing the nodes. The enlargement of the top right
corner also shows the discrete unknowns (nodal velocities and
pressures) stored at each node after the creation of the 
"bulk" Navier-Stokes elements.



\image html multiple_face_elements1.gif "Sketch of a problem that is subject to flux-type boundary conditions along boundaries 0 and 1. The sketch illustrates the degrees of freedom at each node before any FaceElements are attached. " 
\image latex multiple_face_elements1.eps "Sketch of a problem that is subject to flux-type boundary conditions along boundaries 0 and 1. The sketch illustrates the degrees of freedom at each node before any FaceElements are attached. " width=0.75\textwidth


The next figure shows the nodal degrees of freedom after the \c
FaceElements on boundary 0 (shown in red) have been attached. The
\c FaceElements share the existing nodes of the underlying
"bulk" elements and automatically create storage for any additional
nodal unknowns. Here we provide storage for two discrete Lagrange multipliers, 
\f$ \Lambda_x \f$ and \f$ \Lambda_y. \f$ Provided that a single
\c FaceElement is attached to a node, the function
\code
unsigned FaceElement::nbulk_value(const unsigned& j)
\endcode
can be used to determine the number of nodal values at the \c
FaceElement's \c j -th node created by the underlying 
"bulk" element \e before the \c FaceElement was attached. 
It is then easy to identify the additional nodal values associated 
with the \c FaceElement in order to apply the boundary conditions
for the Lagrange multipliers, say.
The methodology is illustrated in the 
<a href="../../../solid/prescribed_displ_lagr_mult/html/index.html#bcs_for_lagrange_multipliers">
the solid mechanics problem referred to earlier.</a>
 



\image html multiple_face_elements2.gif "Sketch illustrating the degrees of freedom at each node after FaceElements have been attached to boundary 0. " 
\image latex multiple_face_elements2.eps "Sketch illustrating the degrees of freedom at each node after FaceElements have been attached to boundary 0. " width=0.75\textwidth



The next figure shows the degrees of freedom after the \c
FaceElements on boundary 1 (shown in green) have also been attached.
These \c FaceElements must create storage for their own
two Lagrange multipliers, \f$ \lambda_x \f$ and 
\f$ \lambda_y \f$. Thus, the corner node (which is attached
to both types of \c FaceElements) has four additional degrees of
freedom after all the \c FaceElements have been created.



\image html multiple_face_elements3.gif "Sketch illustrating the degrees of freedom at each node after all FaceElements have been attached. " 
\image latex multiple_face_elements3.eps "Sketch illustrating the degrees of freedom at each node after all FaceElements have been attached. " width=0.75\textwidth

 
The identification of the additional degrees of freedom via a simple
offset from the degrees of freedom created by the "bulk" element
is now no longer possible. We therefore provide an alternative
mechanism to access the relevant information from the nodes
themselves via the function
\code
unsigned BoundaryNodeBase::index_of_first_value_assigned_by_face_element()
\endcode
which does exactly what it says. If only a single type of \c
FaceElement is attached to a (boundary) node, the unsigned that is returned
by this function is exactly the same as the unsigned that is returned
by the corresponding call to \c FaceElement::nbulk_value(...). 
To cater for the case where multiple \c FaceElements are 
attached to the same node, the above function can take an ID
(which defaults to zero) that identifies which type of
\c FaceElement we are dealing with, so the full interface is, in fact,
\code
unsigned BoundaryNodeBase::index_of_first_value_assigned_by_face_element(const unsigned& id=0)
\endcode
The ID must be established by the user, typically when the 
constructor of the specific \c FaceElement is called. It can then 
be passed on to the \c Nodes when the number of values at the nodes is 
adjusted to accommodate the additional values required by the
\c FaceElement. 


To illustrate this, the code extract shown below provides a 
(partial) listing of the constructor of 
the \c ImposeDisplacementByLagrangeMultiplierElement
that was used in 
<a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
the solid mechanics problem referred to earlier</a>. The constructor
has the usual two arguments that specify the pointer to the "bulk"
element, and the index of the face that the \c FaceElement is to 
be attached to. The final (optional) argument allows the specification
of the ID referred to above. We store the ID in a private 
the element's private member data.
 
\dontinclude solid_traction_elements.h
\skipline Constructor takes a "bulk"
\until FaceElement()
\skipline {
\until Id =


[We omit a few lines of code that are irrelevant for the
present discussion]. Next we specify the number of additional
values required at each of the element's nodes and store them
in a vector. For a \c dim -dimensional bulk element, we need \c dim 
additional values at each node to store the Lagrange multipliers. 

\skipline Dimension of the bulk
\until n_additional_values

Finally, we pass this vector, together with ID that identifies the
type of the \c FaceElement to the function
\c  FaceElement::add_additional_values(...):

\until }

This function creates the additional storage at each node and
updates the node's internal lookup scheme that provides access to the
first value associated with the specified ID.


The fact that the ID is specified as an optional argument means that
the user does not have to provide a dummy ID in cases where
none is required, i.e. in problems involving just a single type 
of \c FaceElement, as in
<a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
the solid mechanics problem referred to earlier.</a>
If a problem does involve multiple \c FaceElements, the 
user will be aware of this when writing the driver code and
can then provide distinct IDs as and when required; see the 
section \ref face_elements_code below.



<HR>
<HR>

\section driver_code The driver code for the FSI problem

 The driver code for the flow through the elastic iliac is almost 
identical to that in the 
<a href="../../unstructured_three_d_fsi/html/index.html">
tutorial considering the same problem in a more simplistic
geometry.</a> Comparing the two driver codes using

\code
sdiff demo_drivers/interaction/vmtk_fsi/vmtk_fsi.cc \
      demo_drivers/interaction/unstructured_three_d_fsi/unstructured_three_d_fsi.cc
\endcode 

shows that the only differences relate to the facts that:
- although the meshes for the two problems are topologically 
  equivalent, the boundary numbers assigned by
  \c Tetgen are different.
- we impose parallel flow at the vessel's in- and outflow
  cross-sections, none of which are not aligned with any of the Cartesian
  coordinate planes. The parallel flow and the imposed pressure
  drop are therefore enforced by attaching 
  \c ImposeParallelOutflowElements, rather than
  \c NavierStokesTractionElements to the in- and outflow boundaries.
  We refer to the corresponding 
  <a href="../../../navier_stokes/vmtk_fluid/html/index.html#parallel_flow">
  single-physics fluids problem</a> for more details on this aspect.
- The problem involves two different types of \c FaceElements: 
  - \c ImposeParallelOutflowElements are used
    to impose parallel flow at the in- and outflow cross-sections.
  - \c ImposeDisplacementByLagrangeMultiplerElements are used to
    deform the boundary of the pseudo-solid fluid mesh 
    to reflect the changes in the geometry of the 
    vessel wall.
  .
  Since both \c FaceElements introduce additional unknowns 
  into the problem, we use the methodology described
  <a href="#face">above</a> to distinguish between the two
  types of Lagrange multipliers.
.


Once again, we shall discuss only those aspects of the code
that are changed from the driver code
discussed in the <a href="../../unstructured_three_d_fsi/html/index.html">
previous tutorial.</a>

<hr>

\subsection namespace The namespace for global parameters

As usual, global parameters are specified in a namespace, which now
includes an \c enum that specifies the IDs for the two different
\c FaceElements.

\dontinclude vmtk_fsi.cc
\skipline start_of_namespace
\until end_of_namespace

<hr>

\subsection constructor The Problem constructor

The general structure of the problem constructor remains
unchanged. There are a few trivial changes in the import of fluid and solid
meshes to reflect the fact that the tetgen boundary numbers are different.

When imposing the boundary conditions for the fluid mesh 
we do <b> not </b> pin the transverse velocities at the in- and outlets because
the parallel flow is now imposed by Lagrange multiplier elements.
Hence the only boundary condition to be applied at the in- and
outflow cross-sections is the pinning of the nodal positions.

The only significant change arises in the application of the boundary
conditions for the Lagrange multipliers. Recall that
\c ImposeDisplacementByLagrangeMultiplierElements are used
to deform the FSI boundary of the fluid mesh so that it
stays in contact with the vessel wall. This constraint 
must be applied along the entire FSI boundary, apart from 
the lines along which it meets the in- and outflow boundaries
where the position of the fluid nodes is already pinned. 
Consequently, we pin the Lagrange multipliers in this part of
the mesh. 

   Similarly, the \c ImposeParallelOutflowElements ensure that the
velocity in the in- and outflow cross-sections is orthogonal to 
these cross-sections. This constraint  must be applied along the 
entirety of the in- and outflow boundaries, apart from the 
lines along which they meet the FSI interface where the fluid 
velocity is already determined by the no-slip condition. 
Consequently, we also pin the second set of Lagrange multipliers 
along this part of the fluid mesh boundary. 

 We loop over all the fluid nodes on the FSI boundary:

\skipline // Loop over nodes on the FSI boundary in the fluid mes
\until boundary_node_pt(b,inod);

For each node we apply the no-slip condition on the wall by 
pinning all three velocity components,

\until pin(2);

and determine whether the node also happens to be located on the
in- or outflow cross-sections:

\until now we know if the node is 

If it is, we pin the Lagrange multipliers associated
with the \c ImposeParallelFlowElements, using
the \c BoundaryNodeBase::index_of_first_value_assigned_by_face_element(...)
function referred to <a href="#face">above</a>, and specifying 
the appropriate ID.
 
\until }

We repeat the same procedure for the Lagrange multipliers associated
with the \c ImposeDisplacementByLagrangeMultiplierElements:

\until }

As usual we document the position of the nodes at which we pinned
the Lagrange multipliers in a file to allow for an external sanity check.

\until end of BC for fluid mesh

The rest of the problem constructor is unchanged.

<hr> 

\subsection face_elements_code Creating the Lagrange multiplier elements
The creation of the \c ImposeDisplacementByLagrangeMultiplierElements
is virtually identical to that in the 
<a href="../../unstructured_three_d_fsi/html/index.html">
previous tutorial.</a>

\skipline start_of_create_lagrange_multiplier_elements
\until face_index_at_boundary(b,e);

The only difference is that we pass the ID that identifies
the type of the \c FaceElement to the constructor of the 
\c ImposeDisplacementByLagrangeMultiplierElements:

\until multiplier_id);

The rest of the function is unchanged:

\until end of create

An equivalent procedure is adopted in the function 
\c create_parallel_flow_lagrange_elements() which follows exactly 
the same steps as in the corresponding 
<a href="../../../navier_stokes/vmtk_fluid/html/index.html#parallel_flow">
single-physics fluids problem</a>, apart from the fact that
we pass the other enumerated ID to the constructor 
of the \c ImposeParallelFlowElements. We therefore omit the listing
of the function and refer to 
<A HREF="../../../../demo_drivers/interaction/vmtk_fsi/vmtk_fsi.cc">
the source code.</A>



<HR>
<HR>

\section comm_ex Comments and Exercises

\subsection realism Making the simulation more "realistic"
The simulation shown at the beginning of this tutorial is 
obviously very crude and suffers from (at least) the sum of the 
shortcomings that we identified in the tutorials for the constituent 
single-physics 
<a href="../../../navier_stokes/vmtk_fluid/html/index.html#realism">fluid</a>
and <a href="../../../solid/vmtk_solid/html/index.html#finer"> solid</a>
problems. You should consider repeating the computation using finer
meshes (consult \c oomph-lib's 
<a href="../../../meshes/mesh_from_vmtk/html/index.html">VMTK
tutorial</a> for details) and explore the use of "flow extensions" 
which allow the (inevitably artificial) boundary conditions to be applied 
further from the region of interest. Adding time-dependence to the
problem, e.g. by subjecting the flow to a periodic fluctuation
in the applied pressure drop would be another interesting
exercise.

 
<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/interaction/vmtk_fsi/">
demo_drivers/interaction/vmtk_fsi/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/interaction/vmtk_fsi/vmtk_fsi.cc">
demo_drivers/interaction/vmtk_fsi/vmtk_fsi.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


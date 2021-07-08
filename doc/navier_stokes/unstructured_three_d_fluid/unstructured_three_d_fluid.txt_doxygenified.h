/**

\mainpage Demo problem: Fluid Mechanics on unstructured 3D meshes

This tutorial provides another demonstration of how to use
3D unstructured meshes for the solution of fluid flow problems.
(The main <a href="../../../meshes/mesh_from_tetgen/html/index.html">
tetgen tutorial</a> already contains a 3D unstructured fluid example.)

The specific problem considered here also serves as a "warm-up problem" for the
<a href="../../../interaction/unstructured_three_d_fsi/html/index.html">
corresponding fluid-structure interaction problem</a> in which the
domain boundary is replaced by an elastic vessel.
 

<HR>
<HR>

\section problem The problem 
Here is a sketch of the problem: Flow is driven through a 3D rigid vessel
made of three rectangular tubes that meet at a common junction. 
The flow is driven by a prescribed pressure
drop between the upstream and the two downstream ends, \f$ \Delta P^* = 
P^*_{in} - P^*_{out}, \f$ and we assume/impose parallel in- and outflow
in the inlet and outlet cross-sections, all of which are parallel to 
\f$ x-y \f$ coordinate plane. 

\image html problem_sketch.gif "Sketch of the domain with boundary conditions. " 
\image latex problem_sketch.eps "Sketch of the domain with boundary conditions. " width=0.75\textwidth

We non-dimensionalise all lengths on the half-width, \f$ W \f$, of the
square main vessel and use the overall pressure drop,
\f$ \Delta P^* \f$ to define the (viscous) 
velocity scale
\f[
{\cal U} = \frac{\Delta P^*\  W}{\mu}.
\f] 
With this choice the Reynolds number becomes
\f[
Re = \frac{\rho {\cal U} W }{\mu} = \frac{ \Delta P^* \rho W^2 }{\mu^2},
\f]
and we choose to drive the flow with  a dimensionless pressure drop of
\f$ \Delta P = 1. \f$ An increase in Reynolds number may 
therefore be interpreted as in increase in the applied 
(dimensional) pressure drop along the vessel.

<HR>
<HR> 


\section mesh 3D unstructured mesh generation
We use 
<A HREF="http://www.wias-berlin.de/~si ">Hang Si's</A>
open-source mesh generator 
<A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A> to generate the 
unstructured tetrahedral mesh "offline". We then
process the output files produced by  
<A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A>
to generate an unstructured \c oomph-lib mesh. 

<A HREF="http://wias-berlin.de/software/tetgen//"> \c Tetgen </A>
requires the specification of the domain boundaries
via so-called facets -- planar surface patches that are bounded
by closed polygonal line segments. In our simple geometry
each of the three tube segments has four external
faces. Together with the three in- and outflow sections
this results in a total of 15 facets.


The 15 facets are defined in a <code>*.poly</code> file that
specifies the position of the vertices, and identifies the 
facets via a "face list" that establishes their bounding vertices. 
The well-annotated <code>*.poly</code> file for this problem 
is located at:

<center>
<a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.poly">demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.poly</a>
</center>

We refer to the <A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen webpages </A>
and \c oomph-lib's own 
<a href="../../../meshes/mesh_from_tetgen/html/index.html">
tetgen tutorial</a> for further details on how to 
create <code>*.poly</code> files.


Here is a plot of the domain specified by 
<a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.poly">fsi_bifurcation_fluid.poly</a>.
The plot was created using \c tetview which is distributed with
<A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A>.

\image html tetgen_boundaries.gif "The domain and its bounding facets. " 
\image latex tetgen_boundaries.eps "The domain and its bounding facets. " width=0.75\textwidth


Note that we have deliberately assigned a different boundary ID to 
each facet. This will make the assignment of the boundary
condition somewhat tedious as the domain boundaries of interest 
tend to be represented by multiple, separate mesh boundaries. However, 
the assignment of distinct boundary IDs for the different facets 
is essential for the automatic generation of boundary coordinates in the
<a href="../../../interaction/unstructured_three_d_fsi/html/index.html">
corresponding fluid-structure interaction problem </a> and is
therefore <b>strongly recommended</b>.


<A HREF="http://wias-berlin.de/software/tetgen//"> \c Tetgen </A> generates
an unstructured volumetric mesh from the information contained
in the <code>*.poly</code> file and outputs the mesh's nodes, 
elements and faces in the files
- <a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.node">demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.node</a>
- <a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.ele">demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.ele</a>
- <a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.face">demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.1.face</a>
.
These files can be used as input to \c oomph-lib's \c TetgenMesh
class, using the procedure discussed in 
<a href="../../../meshes/mesh_from_tetgen/html/index.html">
another tutorial.</a>
   
  The figure below shows a \c tetview plot of the mesh, created with a volume
constraint of 0.2 (i.e. the maximum volume of each tetrahedron is
guaranteed to be less than 0.2 units), using the command
\code
tetgen -a0.2 fsi_bifurcation_fluid.poly 
\endcode
 
\image html mesh_boundaries.gif "Plot of the mesh, generated by tetgen. " 
\image latex mesh_boundaries.eps "Plot of the mesh, generated by tetgen. " width=0.75\textwidth
 
Note how <A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A>
has subdivided each of the 15 original facets specified in the 
<code>*.poly</code> file into a surface triangulation. The 
nodes and tetrahedral elements that are located on (or adjacent to) the
15 original facets inherit their boundary IDs. This is important
when we assign the boundary conditions for the actual computation.

<HR>
<HR>



\section results Results
The plot shown below illustrates the flow field (streamribbons 
coloured by pressure contours) for a Reynolds number of \f$ Re = 100. \f$
The transparent faces show the boundaries of the fluid elements and
illustrate that the mesh is very coarse. As a result, the flow is 
clearly under-resolved, particularly near the two outflow 
cross-sections where the imposition of parallel outflow forces
the fluid velocity to re-adjust rapidly as it approaches the outlet.
(See \ref finer_mesh in \ref comm_ex for a more detailed discussion of this aspect.)

\image html flow_with_stream_ribbons.gif "Flowfield (streamribbons, coloured by the pressure contours) and element boundaries. " 
\image latex flow_with_stream_ribbons.eps "Flowfield (streamribbons, coloured by the pressure contours) and element boundaries. " width=0.75\textwidth


<HR>
<HR>

\section namespace Problem parameters
As usual we define the various problem parameters in a 
global namespace. We define the Reynolds number and 
specify the tractions to be applied at the in- and outflow 
cross-sections:

\dontinclude unstructured_three_d_fluid.cc
\skipline start_namespace
\until end namespace

<HR>
<HR>

\section main The driver code

We specify an output directory, create the \c Problem object 
using ten-node tetrahedral Taylor-Hood elements, and output the
initial guess for the flow field:


\skipline start_main
\until doc_info.number()++;

Next we perform a parameter study in which we increase the Reynolds
number of the flow -- corresponding to an increase in the
applied pressure drop. (As usual we perform a smaller number
of steps in a validation run -- performed when the code is run
with nonzero number of command-line arguments.) 


\until end_of_main


<HR>
<HR>

\section class The Problem class
The \c Problem class has the usual member functions and provides
explicit storage for the fluid mesh and the meshes containing
the \c FaceElements that apply the traction conditions at the
in- and outflow boundaries. We also provide storage for the
IDs of the mesh boundaries that constitute the in- and outflow
boundaries to facilitate the application of the boundary conditions.


\dontinclude unstructured_three_d_fluid.cc
\skipline start_problem_class
\until };



<HR>
<HR>

\section constructor The Problem constructor

We start by building the fluid mesh, using the files created
by <A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen </A>:

\skipline start_constructor
\until split_corner_elements);

(We refer to the subsection \ref split in the 
section \ref comm_ex for a discussion of the \c split_corner_elements flag). 

 
Next, we set up the boundary lookup schemes that determine which elements
are located next to which domain boundaries, and specify the 
IDs of the mesh boundaries that coincide with the in- and outflow
cross-sections. Note that this information reflects the
specification of the boundary IDs in the \c tetgen <code>*.poly</code>
file. 

\until Outflow_boundary_id[1]=2;

Next we apply the boundary conditions. We impose parallel in- and
outflow by pinning the transverse velocities at all nodes
that are located on the in- and outflow boundaries, identifying the
nodes via the boundary IDs just set up. We use the boolean map \c done
to indicate which boundaries we have visited already.

\until done in and outflow

The nodes on all other boundaries (i.e. the ones for which \c done[b] is
still \c false) are subjected to no-slip conditions by pinning
all three velocity components. (This approach facilitates
the "extension" of the mesh discussed in section \ref finer_mesh, but
we note that, in general, keeping track of the boundary IDs associated
with each physical boundary must be done "by hand". )

\until done no slip elsewhere

We complete the build of the Navier-Stokes elements by specifying
the pointer to the Reynolds number,

\until }

and attach the \c FaceElements that apply the imposed in- and outflow
tractions to the appropriate faces of the elements on the in- and
outflow boundaries:

\until create_fluid_traction_elements();

Finally, we combine the various sub-meshes to a combined global mesh
and assign the equation numbers.

\until end constructor



<HR>
<HR>

\section traction Creating the fluid traction elements

The helper function \c create_fluid_traction_elements() 
loops over the bulk elements that are adjacent to the in- and outflow
cross-sections and attaches \c NavierStokesTractionElements
to the relevant faces. We store pointers to the newly-created
elements in the appropriate meshes, and pass pointers to the functions
that specify the imposed traction to the elements.


\skipline start_of_fluid_traction_elements
\until end of create


<HR>
<HR>


\section doc Post-processing

The post-processing routine outputs the flow field.
\until }


<HR>
<HR>

\section comm_ex Comments and Exercises



<HR> 
 

\subsection split Splitting corner elements in unstructured meshes to avoid locking
Meshes generated by <A HREF="http://wias-berlin.de/software/tetgen//"> \c tetgen
</A> tend to be of very high quality and can usually be used
without further modification. However, in Navier-Stokes computations
(or other problems in which the elements have to satisfy an LBB-type 
stability constraint) problems can arise if too many of an 
element's nodes are constrained by boundary conditions, causing
the discretisation to "lock". This tends to 
happen when three of the element's four faces are located on
domain boundaries, as in the case of the top-left and bottom-left
elements in the inflow cross-section shown below.


\image html mesh_without_split.gif "Plot of the mesh generated by tetgen. " 
\image latex mesh_without_split.eps "Plot of the mesh generated by tetgen. " width=0.75\textwidth

Ten-noded tetrahedral elements have nodes at their vertices
and on their edges only. Therefore all the nodes in these particular
two elements are located on
domain boundaries. Furthermore, only one of these (the node 
halfway along the edge that traverses the inflow cross-section) is
unconstrained; the imposition of parallel flow at
the inflow face only  constrains the node's transverse
velocities. Thus, both elements 
have a single velocity degree of freedom (the axial 
velocity at the node on the inflow face) but they retain their four pressure 
degrees of freedom (the pressures at their vertices). This does not
\e necessarily over-constrain the problem (and in the present problem
it does not) -- "locking", due to the presence of too many
pressure degrees of freedom, arises at the level of the global, fully-assembled
problem, not on an element-by-element basis. However, the presence
of such elements makes the occurrence of locking more likely and we
have occasionally come across examples in which locking does occur. 
The remedy would be to pin exactly the required number of
superfluous pressure degrees of freedom (no more and no fewer!) 
to "unlock" the problem. However, the diagnosis of the problem 
is as difficult as its practical resolution: How do we determine 
in advance if a problem will lock and, if so, which pressures 
should be pinned, etc.?

To avoid the problem altogether, \c oomph-lib's \c TetgenMesh 
constructor provides
an optional boolean flag, \c split_corner_elements. If this flag is 
set to \c true, the mesh constructor identifies all elements in which at least 
three faces are located on mesh boundaries. These elements are then
split into four smaller ones, using a localised refinement as 
shown in this plot:



\image html mesh_with_split.gif "Plot of the mesh after splitting corner elements. " 
\image latex mesh_with_split.eps "Plot of the mesh after splitting corner elements. " width=0.75\textwidth

Note that the localised refinement leads to a small deterioration 
in the element quality but this is usually acceptable. If you are concerned 
about this aspect, allow the splitting of corner elements only if 
locking actually occurs (see below), or try to generate a different 
\c tetgen mesh that does not suffer from this problem, e.g. by imposing a 
different volume constraint. 

<center>
<b>How to spot "locking":</b>
</center>
- If a solution can be computed (using a direct solver) the occurrence of
  "locking" is easy to spot: Typically the pressure in over-constrained 
  elements becomes extremely large while the rest of the flow field
  looks fairly normal. Iterative solvers may fail to converge. If you
  suspect that locking may have occurred, try re-computing the
  solution on a mesh with the \c split_corner_elements flag set to 
  \c true.
.

<hr>
 
\subsection finer_mesh Creating a finer mesh and applying more appropriate outflow boundary conditions


The flow field shown in the \ref results section is clearly 
under-resolved and a much
finer mesh would have to be used to fully resolve all flow features.
The problem is particularly bad because we have imposed parallel outflow
(i.e. flow in the \f$ z \f$-direction) at the ends of tubes whose axes are not
aligned with the \f$ z \f$-axis. This creates thin
outflow boundary layers within which the flow
is forced to change direction as it exits the tubes. These
observations motivate the following exercises that 
allow you to explore unstructured mesh generation.

- <b>Exercise 1: Refining the mesh</b> \n\n
  Create a finer mesh by specifying a smaller volume constraint for
  tetgen. For instance, using 
  \n \n
  \code
  tetgen -a0.02 fsi_bifurcation_fluid.poly
  \endcode
  \n
  will generate a much finer mesh, containing 3345 tets. Note that the
  driver code can remain completely unchanged.
  \n\n
- <b>Exercise 2: Add straight outflow vessels</b> \n\n
  The imposition of parallel flow at the outlet boundaries would be
  less problematic if at least the final part of the two outflow tubes was
  aligned with the \f$ z \f$-axis. Modify the file 
  <code><a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid.poly">fsi_bifurcation_fluid.poly</a>
  </code> so that two straight vessel 
  segments, parallel to the \f$ z \f$-axis, are added to the
  mesh. The modification to the <code>*.poly</code> file should be 
  straightforward. You have to add eight additional vertices and faces,
  and re-define the two outflow faces. 
  \n\n
  Here is a sample plot of a modified mesh in which two additional 
  straight segments of different lengths have been attached to the
  downstream tubes. Note that the enumeration of the in- 
  and outflow boundaries was retained, allowing the mesh to be used with 
  the same driver code.
  \n\n
\image html mesh_boundaries_extended.gif "Plot of the mesh, with two straight segments added to the outflow branches. " 
\image latex mesh_boundaries_extended.eps "Plot of the mesh, with two straight segments added to the outflow branches. " width=0.75\textwidth
  \n
  The modified mesh was created with the file
  <code><a href="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/fsi_bifurcation_fluid_with_extended_tubes.poly">fsi_bifurcation_fluid_with_extended_tubes.poly</a>
  </code> which you may wish to consult. However, we do encourage you 
  to do this exercise by yourself first, in order to familiarise 
  yourself with \c tetgen -- you will notice that \c tetgen is not very 
  forgiving (or verbose) when it encounters errors in the <code>*.poly</code>
  file.
  \n\n
- <b>Exercise 3: Use Lagrange multipliers to allow parallel outflow in
  cross-sections that are not aligned with the coordinate axes.</b>
  \n\n
  The main (only?) reason why we tend to impose parallel in- or outflow
  in cross-sections that are aligned with the Cartesian coordinate
  planes is, of course, that such boundary conditions are easiest to apply in a 
  discretisation that is based on the Cartesian form of the
  Navier-Stokes equations.  What if we wished to impose parallel
  outflow from cross-sections that are not aligned with such planes?
  The proper way to do this is to use Lagrange multipliers 
  to enforce the two constraints
  \f[
  {\bf u} \cdot {\bf t}_\alpha = 0 \mbox{\ \ \ \ for $\alpha = 1,2$,}
  \f] 
  where \f$ {\bf t}_\alpha \f$ (for \f$\alpha = 1,2\f$) are the
  two tangent vectors spanning the in- or outflow cross-sections. 
  The imposition of these constraints requires two Lagrange multiplier 
  fields, defined in the in-/outflow cross-sections. Physically, these 
  Lagrange multipliers act as tangential tractions that enforce the 
  parallel in- or outflow.
  The Lagrange multipliers introduce additional degrees of freedom
  into the problem and an implementation can be based 
  on an extension of the existing \c NavierStokesTractionElements,
  using a technique similar to that used to enforce prescribed
  boundary displacements in solid mechanics problems, discussed
  in   <a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
  another tutorial.</a>
  Note that, compared to the previous two exercises, this one is not entirely 
  trivial. An implementation is provided in \c
  ImposeParallelOutflowElement, but, as always, it is more instructive to
  try to write it yourself!
.



<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/">
demo_drivers/navier_stokes/unstructured_three_d_fluid/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/unstructured_three_d_fluid/unstructured_three_d_fluid.cc">
demo_drivers/navier_stokes/unstructured_three_d_fluid/unstructured_three_d_fluid.cc
</A>
</CENTER>
.











<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


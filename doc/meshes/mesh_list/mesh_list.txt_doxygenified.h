/**

\mainpage Existing (structured) meshes

This document provides a brief description of the various 
structured meshes that are distributed with the library.
Most of these meshes were developed for specific example codes but we expect
them to be useful in other problems too. Many of the meshes 
exist in many different variants, usually constructed 
by multiple inheritance. When meshes are relatively trivial
variations of each other, e.g. a basic mesh and its refineable equivalent,
we only list the mesh once. The detailed documentation for the mesh
(obtained by following the link) contains the full inheritance
diagram, showing the mesh's own base classes and any meshes that are derived
from it. For each mesh, we provide a link to a fully-documented
example problem that illustrates its use.


 We stress that the <B>key feature of any given mesh</B> is
<B>its topology</B>, rather than the specific shape for which it
was originally developed. For instance, \c oomph-lib does not
provide a mesh for the discretisation of an annular
domain with quadrilateral elements. However, such a mesh is
trivial to construct by deriving it from, e.g., the
\c SimpleRectangularQuadMesh (a mesh that discretises a
rectangular domain), and then adjusting its nodal positions.
Consult the example in the 
<a HREF="../../../quick_guide/html/index.html#distorted_mesh">
(Not-So-)Quick-Guide</a> for details.

<HR>
<HR>

\section conventions Reminder: General conventions to facilitate the re-use of meshes
Since mesh generation tends to be most tedious part of any numerical
simulation, \c oomph-lib's overall data structure, described in
detail <a HREF="../../../the_data_structure/html/index.html">elsewhere,
</a> aims to facilitate the re-use of meshes in many different
applications. For this purpose most specific \c FiniteElements in \c oomph-lib
are derived by multiple inheritance, combining a "geometric" \c FiniteElement
(e.g. a line/quadrilateral/brick-shaped element from the 
\c QElement<DIM,NNODE_1D> family) with an equations class 
(such as \c PoissonEquations&lt;DIM&gt;) that implements the
weak form of a specific PDE.  For instance, \c oomph-lib's
quadrilateral nine-node Poisson element, \c QPoissonElement<2,3>, 
is derived from the \c PoissonEquations<2> equations class and the
\c QElement<2,3> geometric \c FiniteElement.

The mesh generation process
is mainly concerned with the geometric properties of the mesh's
constituent \c FiniteElements (their topology, number of
nodes, etc.) which are defined by the geometric
\c FiniteElement. This makes it possible to use a mesh that 
was originally developed for the solution of
a Poisson equation with a \c QPoissonElement<2,NNODE_1D>, say, 
for the solution of an advection-diffusion problem
with a \c QAdvctionDiffusionElement<2,NNODE_1D> since both
elements are derived from the same geometric \c FiniteElement. 
Only two aspects of the mesh generation process require
information that is not provided by the geometric \c FiniteElement:
- <B>The number of values to be stored at the various nodes in the mesh:</B>
  For instance, if a mesh that is designed for quadrilateral elements
  from the \c QElement<2,NNODE_1D> family of geometric elements is used
  to solve a (scalar) Poisson equation with nine-noded elements, 
  each \c Node has to store a single value. However, if the same mesh is 
  used for the solution of the 2D Navier-Stokes equations with a
  nine-node quadrilateral elements of type \c QTaylorHoodElement<2>, 
  nodes located
  at the elements' four vertices have to store three values (two
  velocity components and one pressure) whereas nodes located
  at the elements' interior and on their edges only have to store 
  two values (the two velocity components).  \n\n
- <B>The number of "history values" required by the timestepping
  procedure (if any):</B> For instance, if the mesh is used for the 
  solution of a Poisson problem, no "history" values are required. If the same
  mesh is used for the solution of an unsteady heat problem, the
  number of history values required is determined by the \c
  TimeStepper used to approximate the time-derivative of the
  nodal values.
.
The meshes used in our example codes (and all the meshes listed below)
have a common structure that allows the required information
to be become available to the mesh constructor:
-# All meshes are templated by the element type, \c ELEMENT.\n\n
-# The final argument of all mesh constructors is a pointer to a 
   \c TimeStepper. We provide a default for this argument --  a 
   pointer to the \c Steady<0> timestepper, defined (as static
   member data) in the \c Mesh base class. 
. 
The availability of the template parameter allows the mesh generator
to build elements of the required type. \c Nodes are generally built by the 
\c FiniteElement::construct_node(...) function whose arguments are
the \c Node's local node number within the current element,
and a pointer to the timestepper. These arguments provide all the 
information that is required 
to build \c Nodes with the right number of values (as required by 
the element) and history values (as required by the \c TimeStepper).
When the function \c FiniteElement::construct_node(...) is called,
it creates the new \c Node, stores a pointer to the newly created \c Node
in the \c FiniteElement's own lookup scheme, and returns that pointer.
This allows the pointer to the newly created \c Node to be stored in 
the \c Mesh's own lookup scheme for its constituent \c Nodes.
The <a HREF="../../../quick_guide/html/index.html#mesh">
(Not-So-)Quick-Guide</a> contains a section that explains
how to write simple meshes.
<HR>
<HR>

\section faq Mesh FAQ
When using a mesh that was originally developed for a different
application, it is sometimes necessary establish the node/element/boundary
numbering scheme employed by the mesh writer. While we generally
assume that the mesh writer will have carefully documented his/her
code, here is what to do if he/she hasn't: \n\n

- <B>How do I find out how the elements are numbered?</B> \n\n
  The function \c Mesh::output(...) outputs the elements in the order in
  which they are stored internally. If you prefer a different 
  element numbering scheme you can re-number the elements;
  see e.g. the member function \c element_reorder() in the
  <a HREF="classoomph_1_1RectangularQuadMesh.html"><B>
  RectangularQuadMesh</B></a> class.
  \n\n         
- <B>How do I find out how the nodes are numbered?</B> \n\n 
  The function \c Mesh::node_pt(j) provides pointer-based access to the
  \c j - th \c Node in the mesh. To plot a node's position, you can
  determine its coordinates from the \c Node::x(...) function.
  \n\n
- <B>How do I find out how the mesh boundaries are numbered?</B> \n\n
  The function \c Mesh::output_boundaries(...) outputs the nodes located
  on the mesh boundaries in a 
  <a HREF="http://www.tecplot.com">tecplot</a>-able format. Nodes 
  that are located on separate mesh boundaries are contained in separate 
  <a HREF="http://www.tecplot.com">tecplot</a> zones. 
  \n\n
.

With this information it should be straightforward to use any
of the meshes listed below in one of your own problems. 
The example code listed next to each mesh illustrates its use
in an actual driver code.
If you develop a new mesh, let us know! If it is written
according to \c oomph-lib's
<a HREF="../../../coding_conventions/html/index.html">coding
standards</a>,  we'll be delighted to
include it into the library. 

<HR>
<HR>

\section mesh_list Mesh list

<CENTER>
<TABLE>
<tr>
<td>
<B>Mesh</B>
</td>
<td>
<B>Representative Mesh plot</B>
</td>
</tr>
<tr>
<td>
<a HREF="classoomph_1_1OneDMesh.html"><B>
<CODE>OneDMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<1,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes.
- A refineable version of this mesh exists.
. 
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../poisson/one_d_poisson/html/index.html">
  the solution of the 1D Poisson equation. 
  </a>
.
</td>
<td>
\image html one_d_mesh.gif " " 
\image latex one_d_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>
<td>
<a HREF="classoomph_1_1SimpleRectangularQuadMesh.html"><B>
<CODE>SimpleRectangularQuadMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes.
.
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../poisson/two_d_poisson/html/index.html">
  the solution of the 2D Poisson equation. 
  </a>
.
</td>
<td>
\image html simple_rectangular_quadmesh.gif " " 
\image latex simple_rectangular_quadmesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>
<td>
<a HREF="classoomph_1_1RectangularQuadMesh.html"><B>
<CODE>RectangularQuadMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This is a slightly more sophisticated version of the
  \c SimpleRectangularQuadMesh discussed above; it allows for
  non-uniform spacing of the nodes, and periodicity in the x and
  y-directions. 
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes.
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable variant of this mesh is used in the driver code for 
  <a HREF="../../../advection_diffusion/two_d_adv_diff_adapt/html/index.html">
  the solution of the 2D Advection Diffusion equation. 
  </a>
.
</td>
<td>
\image html simple_rectangular_quadmesh.gif " " 
\image latex simple_rectangular_quadmesh.eps " " width=0.75\textwidth
</td> 
</tr>
<tr>
<td>
<a HREF="classoomph_1_1TwoDAnnularMesh.html"><B>
<CODE>TwoDAnnularMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This is a "wrapped around" version of 
  \c RectangularQuadMesh discussed above. It can either be used
  as a complete annulus (in which case periodicity is enforced)
  or a gap can appear in the annulus.
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../helmholtz/scattering/html/index.html">
  the 2D Helmholtz equation</a >and the 
  <a HREF="../../../time_harmonic_linear_elasticity/elastic_annulus/html/index.html">time-harmonic 
  linear elasticity equations</a>.
  </a>
.
</td>
<td>
\image html annular_meshes.gif " " 
\image latex annular_meshes.eps " " width=0.75\textwidth
</td> 
</tr>
<tr>
<td>
<a NAME="channel_with_leaflet"></a>
<a HREF="classoomph_1_1ChannelWithLeafletMesh.html"><B>
<CODE>ChannelWithLeafletMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- Mesh for the simulation of flow in a 2D channel that is
  partially occluded by a moving leaflet. 
- Leaflet must be represented by a \c GeomObject
- Nodes along the leaflet are duplicated to allow for a pressure
  jump across the leaflet even for discretisations that impose 
  continuous pressure distributions across element boundaries. 
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes.
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable variant of this mesh with an algebraic node update
  is used in the driver code for 
  <a HREF="../../../navier_stokes/channel_with_leaflet/html/index.html">
  the simulation of flow in a 2D channel that is partially occluded
  by a moving leaflet</a> and also in 
  <a HREF="../../../interaction/fsi_channel_with_leaflet/html/index.html">
  corresponding FSI problem. </a>
.
</td> 
<td>
\image html channel_with_leaflet_mesh.gif " " 
\image latex channel_with_leaflet_mesh.eps " " width=0.75\textwidth
</td> 
</tr>
<tr>
<td>
<a HREF="classoomph_1_1SimpleRectangularTriMesh.html"><B>
<CODE>SimpleRectangularTriMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This is a simple structured mesh made of triangular elements.
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c TElement<2,NNODE_1D>.
.
<B>Example driver code:</B> \n
- This mesh is used in the self-test code
  <a HREF="../../../../demo_drivers/self_test/poisson/convergence_tests/t_convergence_2d.cc">
  <CODE>t_convergence_2d.cc</CODE>
  </a>
.
</td>
<td>
\image html simple_rectangular_trimesh.gif " " 
\image latex simple_rectangular_trimesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>
<td> 
<a HREF="classoomph_1_1FishMesh.html"><B>
<CODE>FishMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes.
- The curvilinear boundaries are represented by a \c GeomObject
  and the mesh has a \c Domain representation, allowing a
  \c MacroElement - based node update. 
- There is also a version of the mesh that performs the node update
  in response to changes in the domain boundary by an algebraic 
  node update. 
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable variant of this mesh is used in the driver code for
  <a HREF="../../../poisson/fish_poisson/html/index.html">the adaptive solution of the 2D Poisson equation in a fish-shaped domain.</a>
.
</td>
<td>
\image html fish_mesh.gif " " 
\image latex fish_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a NAME="collapsible_channel">
<a HREF="classoomph_1_1CollapsibleChannelMesh.html"><B>
<CODE>CollapsibleChannelMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes
- The curvilinear boundary is represented by a \c GeomObject
  and the mesh has a \c Domain representation, allowing a \c
  MacroElement - based node update. 
- There is also a version of the mesh that performs the node update
  in response to changes in the domain boundary by an algebraic 
  node update.  
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../navier_stokes/collapsible_channel/html/index.html">
  the solution of the 2D unsteady Navier-Stokes equations in 2D
  channel with an oscillating wall.
  </a>
.
</td>
<td>
\image html collapsible_channel_mesh.gif " " 
\image latex collapsible_channel_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>
<td>
<a HREF="classoomph_1_1CylinderWithFlagMesh.html"><B>
<CODE>CylinderWithFlagMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This \c Mesh was mainly developed for the solution 
  of Turek & Hron's  FSI benchmark problems. The curvilinear
  boundaries of the cylinder and the "flag" are represented
  by \c GeomObjects. 
- A refineable version of the mesh exists.
- The node-update in response to changes in the shape of the "flag"
  can be performed by a version based on an \c AlgebraicMesh or
  using a \c Domain/MacroElement - based node update. 
- The bulk elements have to be derived from the
  geometric finite element \c QElement<2,NNODE_1D>.
.
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../interaction/turek_flag/html/index.html">
  Turek & Hron's FSI benchmark problems
  </a> and their 
  <a HREF="../../../navier_stokes/turek_flag_non_fsi/html/index.html"> 
  non-FSI equivalents.</a>
. 
</td>
<td>
\image html cylinder_with_flag_mesh.gif " " 
\image latex cylinder_with_flag_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>
<td>
<a HREF="classoomph_1_1BrethertonSpineMesh.html"><B>
<CODE>BrethertonSpineMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This \c SpineMesh was mainly developed for the simulation
  of the Bretherton problem but it can, of course, also be used 
  in other problems. The mesh topology would be suitable
  for the simulation of flows in a bifurcating channel, say. 
- The bulk elements have to be derived from the
  geometric finite element \c QElement<2,NNODE_1D>.
.
<B>Example driver code:</B> \n
- This mesh is used in the driver code for 
  <a HREF="../../../navier_stokes/bretherton/html/index.html">
  the simulation of the Bretherton problem.
  </a>
.
</td>
<td>
\image html bretherton_spine_mesh.gif " " 
\image latex bretherton_spine_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1QuarterCircleSectorMesh.html"><B>
<CODE>QuarterCircleSectorMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes
- The curvilinear boundary is represented by a \c GeomObject
  and the mesh has a \c Domain representation, allowing a \c
  MacroElement - based node update. 
- There is also a version of the mesh that performs the node update
  in response to changes in the domain boundary by an algebraic 
  node update.  
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable version of this mesh is used in the driver code for 
  <a HREF="../../../navier_stokes/osc_ellipse/html/index.html">
  the simulation of flow inside a oscillating ellipse.
  </a>
.
</td>
<td>
\image html quarter_circle_sector_mesh.gif " " 
\image latex quarter_circle_sector_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1SimpleCubicMesh.html"><B>
<CODE>SimpleCubicMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<3,NNODE_1D>.
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- This mesh is used in the self-test code
  <a HREF="../../../../demo_drivers/self_test/poisson/convergence_tests/q_convergence_3d.cc">
  <CODE>q_convergence_3d.cc</CODE>
  </a>
.
</td>
<td>
\image html simple_cubic_mesh.gif " " 
\image latex simple_cubic_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1SimpleCubicTetMesh.html"><B>
<CODE>SimpleCubicTetMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This is a simple structured mesh for tet elements.
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c TElement<3,NNODE_1D>.
.
<B>Example driver code:</B> \n
- This mesh is used in the self-test code
  <a HREF="../../../../demo_drivers/self_test/poisson/convergence_tests/t_convergence_3d.cc">
  <CODE>t_convergence_3d.cc</CODE>
  </a>
.
</td>
<td>
\image html simple_cubic_tetmesh.gif " " 
\image latex simple_cubic_tetmesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1QuarterTubeMesh.html"><B>
<CODE>QuarterTubeMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<3,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes
- The curvilinear boundary is represented by a \c GeomObject
  and the mesh has a \c Domain representation, allowing a \c
  MacroElement - based node update. 
- There is also a version of the mesh that performs the node update
  in response to changes in the domain boundary by an algebraic 
  node update.  
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable version of this mesh is used in the driver code for 
  <a HREF="../../../navier_stokes/three_d_entry_flow/html/index.html">
  the simulation of 3D entry flow into a cylindrical tube.
  </a>
.
</td>
<td>
\image html quarter_tube_mesh.gif " " 
\image latex quarter_tube_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1TubeMesh.html"><B>
<CODE>TubeMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<3,NNODE_1D>.
- This mesh forms the basis for numerous derived meshes that describe
  topologically-tube-shaped domains.
- The entire domain is represented by a \c GeomObject
  and the mesh has a \c Domain representation, allowing a \c
  MacroElement - based node update. 
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable version of this mesh is used in the driver code for 
  <a HREF="../../../navier_stokes/curved_pipe/html/index.html">
  the simulation of 3D flow in a curved cylindrical pipe.
  </a>
.
</td>
<td>
\image html tube_mesh.gif " " 
\image latex tube_mesh.eps " " width=0.75\textwidth
</td>
</tr>
<tr>   
<td>
<a HREF="classoomph_1_1EighthSphereMesh.html"><B>
<CODE>EighthSphereMesh&lt;ELEMENT&gt;</CODE>
</B></a> \n\n
- This mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<3,NNODE_1D>.
- A refineable version of this mesh exists.
.
<B>Example driver code:</B> \n
- The refineable version of this mesh is used in the driver code for 
  <a HREF="../../../poisson/eighth_sphere_poisson/html/index.html">
  the adaptive solution of the 3D Poisson equation.
  </a>
.
</td>
<td>
\image html eighth_sphere_mesh.gif " " 
\image latex eighth_sphere_mesh.eps " " width=0.75\textwidth
</td>
</tr>
</TABLE>
</CENTER>   

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
\*/


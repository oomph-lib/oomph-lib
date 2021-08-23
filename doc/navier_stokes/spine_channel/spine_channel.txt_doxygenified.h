/**

\mainpage Example problem: Steady 2D finite-Reynolds-number flow in a channel of non-uniform width -- An introduction to Spine meshes.



Many previous examples demonstrated \c oomph-lib's ability to solve
problems on domains with moving, curvilinear boundaries. These
examples had the following common features: 
- The motion of the curvilinear domain boundaries was prescribed.
- The domain was discretised by \c Domain / \c MacroElement - based
  meshes. Recall that in such meshes
  the function \c Mesh::node_update() updates the position of
  \e all of its constituent nodes in response to changes in the
  shape/position of the geometric objects that define its curvilinear
  boundaries. The update of the nodal positions is performed on an 
  element-by-element basis and each element determines the new 
  positions of its nodes by referring to the \c MacroElement 
  representation of the domain. 
- The governing equations were implemented in their Arbitrary 
  Eulerian Lagrangian (ALE) form, in which the mesh velocity is
  determined from the "history values" of the nodal positions. 
.

We will now consider problems in which the position of the domain
boundary is unknown and has to be determined as part of the overall
solution. This situation arises, e.g., in free-surface fluid flow
problems and in fluid-structure interaction problems. We shall explain
why \c Domain / \c MacroElement - based node update strategies 
are unlikely to be efficient for such problems and then introduce
the "Method of Spines" as one of a number of sparse (and therefore more
efficient) node-update strategies available in \c oomph-lib.



<HR>

\section why_sparse Why we need sparse node updates

The sketch below shows a free-surface fluids problem in which
the "height" of the fluid domain (parametrised by the scalar function
\f$ x_2 = h(x_1) \f$) is unknown. The lower half of the sketch shows
a body-fitted finite-element mesh (the nodes and elements are shown in
dark blue) that discretises the fluid domain.

\image html film_sketch.gif "Sketch of a free-surface fluid problem. " 
\image latex film_sketch.eps "Sketch of a free-surface fluid problem. " width=0.75\textwidth


 Assume now that we have some discrete representation of the unknown
free surface so that \f$ h(x_1) \f$ is approximated by
a function that involves a finite number of discrete unknowns 
\f$ H_i \  (i=1,...,N_H)\f$ . In principle, this allows us to 
represent the unknown boundary by a \c GeomObject in which the 
unknowns \f$ H_i \  (i=1,...,N_H)\f$ play the role of "geometric
\c Data", i.e. \c Data whose values determine the shape of the
geometric object. Once the curvilinear
boundary is represented by a \c GeomObject, the update of the
nodal positions in the "bulk" mesh could be performed 
by the  \c Domain / \c MacroElement - based methods referred to above.


How exactly the unknowns \f$ H_i \  (i=1,...,N_H)\f$ are 
determined is irrelevant for the purpose of this discussion (
in free surface flow problems, the relevant equation is the
kinematic free-surface condition discussed in 
<a href="../../../navier_stokes/single_layer_free_surface/html/index.html#kinematic_condition_theory">another tutorial</a>) -- we simply assume that there \e are some equations
that determine their values. The feature we wish to focus on here is
that the solution of the problem by Newton's method requires the
computation of the derivatives of \e all discrete residuals with 
respect to \e all unknowns in the problem. \c oomph-lib's
Navier-Stokes elements compute the element residual vectors 
(the residuals of the discretised momentum and continuity equations, 
evaluated for the current values of the unknowns), and the derivatives of these
residuals with respect to the elements' velocity and pressure
degrees of freedom. Clearly, the entries in the element's residual
vector also depend on the position of the element's constituent
nodes, which, in a free-boundary problem, are determined (via the
node update function) by the unknowns \f$ H_i \  (i=1,...,N_H)\f$
that discretise the position of the free surface. 





  

 The main purpose of this example is to demonstrate the use (and the
 creation) of so-called "spine meshes". Such meshes
are similar to the \c MacroElement / \c Domain - based meshes employed
in many previous examples, in that they allow the nodal positions to be updated
in response to changes in the shape of their (curvilinear) domain boundaries. 
The key feature of "spine meshes" is that the node update can be
performed on a node-by-node basis -- this an important requirement for the
efficient solution of free-boundary and fluid-structure
interaction problems in which the position of the nodes in the "bulk
mesh" is determined by the (unknown) position of the domain boundary. 
The efficient evaluation of the so-called "shape derivatives"
(the derivatives of the residuals of the equations discretised by
the elements in the "bulk mesh" with respect to the unknowns
that determine the position of the domain boundary) 


The idea behind spine-based node updates is illustrated in the sketch
below. Assume that the position of the domain boundary is parametrised
by a scalar function, so that, for instance, \f$ x_2 = h(x_1)\f$, 
where \f$ h(x_1)\f$ may have to be determined as part of the solution 
(e.g. in free-surface fluids problems -- the origin of the 
"Method of Spines").

\image html spine_sketch.gif "Sketch of the Method of Spines. " 
\image latex spine_sketch.eps "Sketch of the Method of Spines. " width=0.75\textwidth

Further, assume that the mesh topology is such that the
mesh's \f$ N_{node} \f$ nodes are distributed along \f$ N_{spine} \f$ 
lines that are (topologically) orthogonal to the free boundary. 
We refer to these lines as the "spines" and denote
the "height" of the domain, measured along spine \f$ s \f$ by
\f$ H_s \ (s=1,...,N_{spine}) \f$. We associate each node with 
a particular spine (so that node \f$ j \f$ is located on spine
\f$ s_j\f$)  and locate it along a fixed fraction 
\f$ \omega_j \f$ along "its" spine. The position of node \f$ j \f$
may therefore be written as
\f[
{\bf x}_j =  {\bf B}_{s_j} + \omega_j \, H_{s_j}  \, 
{\bf S}_{s_j} \ \ \ \ \ \ \ \ \ \ \ \ \  (1) 
\f]
where \f$ {\bf B}_s \f$ is the vector to the "base" of spine \f$ s \f$,
and \f$ {\bf S}_s \f$ the unit vector along that spine. 

A key feature of this method is that 

Determining the nodal positions via the "Method of Spines" 
equation (1)

 

Spine-based node updates 
This document has two main parts:
- In \ref part1 we demonstrate how to use a \c SpineMesh
- In \ref part2 we explain the general "philosophy" behind
  spine-based node-updates and demonstrate their implementation.
.

<HR>
<HR>

\section part1 Part 1: Flow through a channel of non-uniform width


\subsection example The example problem
We shall illustrate the use of \c SpineMeshes by considering 
the problem of steady 2D flow through a channel of non-uniform width.

<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The steady 2D Navier-Stokes equations in a channel of non-uniform width.</B>
</CENTER> 
Solve
\f[
Re \ u_j\frac{\partial u_i}{\partial x_j} 
= - \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right),
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
in the region \f$ D = \left\{(x_1,x_2) \ \bigg| \ x_1 \in [0,L], 
x_2 \in [0,h(x_1)] \right\} \f$,
where
\f[
h(x_1) = \left\{
\begin{array}{cl}
H & 0 \leq x_1 \leq L_1 \\
H + A\sin\left(\frac{x_1-L_1}{L_2-L_1}\right) & L_1 \leq x_1 \leq L_2 \\
H & L_2 \leq x_1 \leq L
\end{array} \right.
\f]
shown in this sketch
\image html channel_sketch.gif "Sketch of the problem. " 
\image latex channel_sketch.eps "Sketch of the problem. " width=0.75\textwidth

subject to the no-slip Dirichlet boundary conditions on the top and bottom rigid walls
\f[
\left. \mathbf{u}\right|_{\partial D_{wall}}=(0,0),
\f]
parallel, parabolic inflow on the left inflow boundary, 
\f$ \partial D_{inflow} = \{ (x_1,x_2) \ |  \ x_1=0 \}\f$ ,
\f[
\left. \mathbf{u}\right|_{\partial D_{inflow}}
=\left(x_2\left(H-x_2\right),0\right),
\f]
and axially traction-free,
parallel outflow on the outflow boundary, 
\f$ \partial D_{outflow} = \{ (x_1,x_2) \ | \ x_1=L \}\f$ ,
\f[
\left. u_2\right|_{\partial D_{outflow}}=0.
\f]
</TD>
</TR>
</TABLE>  
</CENTER>

<HR>
<HR>

\section results Results

The figures below show the results (carpet plots of the two 
velocity components and the pressure, and a contour plot of the
pressure distribution with superimposed streamlines), 
obtained from computations with Taylor-Hood and Crouzeix-Raviart 
elements for a channel of length \f$ L = 2.7\f$ , height \f$ H=1.0 \f$ , with 
deflection amplitude \f$ A = 0.4 \f$, and a Reynolds number of 
\f$ Re=100\f$.

\image html TH.gif "Plot of results computed with 2D Taylor-Hood elements. " 
\image latex TH.eps "Plot of results computed with 2D Taylor-Hood elements. " width=0.75\textwidth
\image html CR.gif "Plot of results computed with 2D Crouzeix-Raviart elements. " 
\image latex CR.eps "Plot of results computed with 2D Crouzeix-Raviart elements. " width=0.75\textwidth


<HR>
<HR>

\section namespace Global parameters
The Reynolds number is the only parameter in 
this problem. As usual, we define and initialise it in a namespace:

\dontinclude spine_channel.cc
\skipline start_of_namespace
\until end_of_namespace

<HR>
<HR>

\section main The driver code

We start by creating a \c DocInfo object to define the output 
directory.

\skipline start_of_main
\until doc_info.number()=0

When using spines, we must use elements augmented by the \c
SpineElement<ELEMENT> class. This class, adds the functionality to be
updated using the method of spines (i.e storing a vector of pointers
to the spines and allocating equations numbers associated with the
spines degrees of freedom). We build the problem using \c
SpineElement<TayloorHoodElement<2>>.

We now build and solve the problem with Spine-Taylor-Hood elements,
then repeat for Spine-Crouzeix-Raviart elements.

\skipline Solve problem with Taylor Hood elements
\until end_of_main

<HR>
<HR>

\section problem The problem class 
The problem class for this example is very similar to our previous
<A HREF="../../driven_cavity/html/index.html#problem"> steady
Navier-Stokes examples </A>. We store the
height of the channel as private data, this is because we need it to
set the inflow boundary condition.

\dontinclude spine_channel.cc
\skipline start_of_problem_class
\until end_of_problem_class

Note that the absence of boundary conditions on the right boundary
(1), causes a zero traction condition to be applied there. This
implies that we should not fix a pressure degree of freedom.

<HR>
<HR>

\section constructor The problem constructor
We begin by setting all the mesh parameters, and building it.

\skipline start_of_constructor
\until UpperWall);

We then pin the velocities on the left, top and bottom boundaries (3,2
and 0).

\until end loop over boundaries

Finally, we pass a pointer to the Reynolds number to each element and
assign the equation numbers.

\until end_of_constructor

<HR>
<HR>


\section part2 Part 2: How to create a SpineMesh


Spine-based meshes have their origin 
in free-surface fluid-mechanics problems where they were first (?)
introduced by Kistler & Scriven in their paper 

<center>
 Kistler, S.F. & Scriven, L.E. ``Coating Flows.'' In:
``Computational Analysis of Polymer Processing,'' Pearson, J.R.A. &
Richardson, S.M. (eds.); Applied Science Publishers, London (1983).
</center>

\c oomph-lib's
 SpineMeshes provide a generalisation of their node-update
techniques.

In the past when solving a problem in a domain with curved boundaries,
we have made a specific Mesh for the domain, and use a geometric object to
define its curved wall(s).

\image html mesh.gif "Diagram of a ChannelSpineMesh, going through the process of updating its central nodes, where the heights have been changed. " 
\image latex mesh.eps "Diagram of a ChannelSpineMesh, going through the process of updating its central nodes, where the heights have been changed. " width=0.75\textwidth

To generalise this approach \c oomph-lib makes use of spines. 
- A \c Spine is most easily visualised a line of a certain height, in
  one coordinate direction.
- A \c SpineNode is a \c Node located at a certain fraction along a \c
  Spine. 
- A \c SpineMesh is a \c Mesh which will update using spines.
- A \c SpineElement<ELEMENT> takes a "normal" element and adds the
  functionality to work with spines.
.

The creation of a \c SpineMesh is discussed in detail below \ref
spine_mesh. While the picture above demonstrates the ability of
SpineNode to individually update using the function \c
spine_node_update(spine_node_pt).


\subsection spine_mesh Making a SpineMesh 

In this example we use a \c SpineMesh to model the domain. This mesh
constitutes three regions: left (0), centre (1) and right (2), the
left and right regions have a constant height, while in the centre
region the height varies. These heights are defined by two
geometric objects, a 
<A HREF="../../../the_data_structure/html/classoomph_1_1StraightLine.html"> 
straight
line </A> and a deflected line respectively. 
We will discuss the necessary steps taken to create
this mesh (the complete documentation for this specific mesh can be
found <A HREF="../../../the_data_structure/html/classoomph_1_1ChannelSpineMesh.html">
here </A>).

To begin, we create a new class templated by a element and inheriting
from \c RectangularQuadMesh<ELEMENT> and \c SpineMesh. The latter adds
the functionality needed for using Spines.

\dontinclude channel_spine_mesh.template.h
\skipline template<
\until public SpineMesh

All SpineMeshs must include a function \c spine_node_update(SpineNode*
spine_node_pt), this will describe the operations performed when
updating every SpineNode in the mesh. First we find the SpineNodes
fraction along the spine.

\skipline virtual void spine_node_update
\until fraction()

We then get the local coordinate on the geometric object that defines the
upper wall.

\until geom_parameter(0)

Finally we use the use the local coordinate to get the position of the
geometric object and set the first coordinate value of the node.

\until }

We store the number of elements in the x direction in each region, the
number of elements in the y direction, the lengths of each region and
the height of the uniform boundary, as well as pointers to the two
geometric objects.

\skipline Number of elements in the left region
\until Straight_wall_pt;

All these details are passed to the constructor, except for the
pointer to the geometric object for the uniform wall. The constructor calls
the empty constructor for \c RectangularQuadMesh<ELEMENT>, copies
these values to their storage in the mesh.

\dontinclude channel_spine_mesh.template.cc
\skipline template<class ELEMENT>
\until Wall_pt

We then assign all the parameters for the \c
RectangularQuadMesh<ELEMENT>, create the geometric object for the uniform
wall and call the function \c build_channel_spine_mesh(...)

\until }

When we call the function \c build_channel_spine_mesh(...), it calls
its counterpart in the \c RectangularQuadMesh<ELEMENT>, then store the
numbers of elements in each direction in each region (and all at once).

\dontinclude channel_spine_mesh.template.cc
\skipline ::build_channel_spine_mesh
\until n_x2

We then allocate memory for the elements and spines in each region in
each region.

\until end Allocating memory

Now we allocate storage for the parameters used to build the spines.

\until n_prev_elements = 0;

Now we create the first Spine with unit length, pin the height (since
it is not a degree of freedom in this mesh) and push the spine back
onto the \c Spine_pt.

\until push_back(new_spine_pt)
 
We then set the \c spine_pt() of the first node to the Spine we just
created, assign the nodes \c fraction() to zero and define the node as
part of this mesh.

\until nod_pt->spine_mesh_pt() = this;

We then mark the node as part of the left (0) region.

\until node_update_fct_id() = 0;

When we built the Spine, we set its height to 1.0. We now need to
assign its height from the \c Straight_wall_pt and assign all the
information needed to update the mesh to the Spine.

First we set the value of \f$\zeta\f$ and get the geometric object and the
local coordinate.

\until locate_zeta

Then we store these geometric parameters in the Spine.

\until nod_pt

We then set the height of the Spine according to the geometric object.

\until height()

Finally we set the Spines' pointer to the geometric object.

\until }

Now we loop vertically along the spine, adding a pointer to the spine
to each element, and assigning the fraction for each node on this
spine, define each node as part of this mesh, and mark it as part of
the left region.

\until end loop over elements

We then loop over the remaining spines in the left region repeating
this process, except that the first spine in each element (except the
first) is copied from the last spine in the previous element.

We then repeat this process for the centre and right regions, using
the correct geometric objects to define the upper wall, which can be
examined in detail 
<A HREF="../../../the_data_structure/html/classoomph_1_1ChannelSpineMesh.html">
here </A>.

<HR>
<HR>

\subsection exercises Exercises

-# Investigate what happens when a pressure degree of freedom is fixed. 
-# Try creating a new geometric object, which creates a triangular
   indentation in the central region of the upper wall as shown.

\image html Exercise.gif "Plot of the solution to the problem specified in the above exercise computed with 3x3 Taylor-Hood elements and Re=100. " 
\image latex Exercise.eps "Plot of the solution to the problem specified in the above exercise computed with 3x3 Taylor-Hood elements and Re=100. " width=0.75\textwidth


   




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


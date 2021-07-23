/**

\mainpage Demo problem: How to create refineable meshes in domains with curvilinear and/or moving boundaries

 

In an 
<A HREF="../../../../doc/poisson/two_d_poisson_adapt/html/index.html">
earlier example</A> we demonstrated how easy it is to "upgrade" 
an existing quad mesh to a \c RefineableMesh that can be used
with \c oomph-lib's mesh adaptation routines. The "upgrade" 
was achieved by multiple inheritance: We combined the basic 
(non-refineable) mesh object with \c oomph-lib's \c RefineableQuadMesh
-- a class that implements the required mesh adaptation procedures,
using  \c QuadTree - based refinement techniques for meshes
that contain quadrilateral elements. During the
refinement process, selected elements are split into four "son"
elements and the nodal values and coordinates of any newly created 
nodes are determined by interpolation from the "father" element.
This procedure is perfectly adequate for problems with polygonal
domain boundaries in which the initial coarse mesh provides a
perfect representation of the domain. The situation is more
complicated in problems with curvilinear domain boundaries since
we must ensure that successive mesh refinements lead to an
increasingly accurate representation of the domain boundary.


To illustrate these issues we (re-)consider the 2D Poisson problem
<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Two-dimensional model Poisson problem in a non-trivial domain</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2 \frac{\partial^2u}{\partial x_i^2} = -1,
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the fish-shaped domain \f$D_{fish} \f$ with homogeneous
Dirichlet boundary conditions
\f[
\left. u\right|_{\partial D_{fish}}=0.
\ \ \ \ \ \ \ \ \ \ (2)
\f]
</TD>
</TR>
</TABLE>  
</CENTER>  

In Part 1 of this document we shall explain how \c oomph-lib's mesh 
adaptation procedures employ the \c Domain and \c MacroElement objects to adapt
meshes in domains with curvilinear  boundaries. In Part 2, 
we demonstrate how to create new \c Domain objects.


<HR>

\section adapt Part 1: Mesh adaptation in domains with curvilinear boundaries, using Domain and MacroElement objects



The plot below shows the domain  \f$D_{fish} \f$, represented by the
multi-coloured, shaded region and its (extremely coarse) 
discretisation with four four-node quad elements. The elements' 
edges and nodes are shown in black.

\image html coarse_fish_mesh_and_domain.gif "The fish-shaped domain and its discretisation with four four-node quad elements " 
\image latex coarse_fish_mesh_and_domain.eps "The fish-shaped domain and its discretisation with four four-node quad elements " width=0.75\textwidth


Obviously, the curvilinear boundaries of the fish-shaped domain
(arcs of circles) are very poorly resolved by the elements' straight
edges. Simple mesh adaptation, based on the techniques
described in the 
<A HREF="../../../../doc/poisson/two_d_poisson_adapt/html/index.html">
earlier example</A> will not result in convergence to the
exact solution since the refined mesh never approaches the exact
domain geometry:

\image html fish_poisson_soln_no_macro.gif "Plot of the mesh adaptation without MacroElements " 
\image latex fish_poisson_soln_no_macro.eps "Plot of the mesh adaptation without MacroElements " width=0.75\textwidth


To overcome this problem, the mesh adaptation routines must be given access to
an exact, analytical representation of the actual domain. This is
the purpose of \c oomph-lib's \c Domain object. A \c Domain
provides an analytical description of a mathematical
domain, by decomposing it into a number of so-called \c MacroElements.
Each \c MacroElement provides a mapping between a set of local
and global coordinates \f$ {\bf r}_{macro}({\bf s})\f$ -- similar to the 
mapping between the local and global coordinates in a finite element. 
The key difference between the two types of element is that the
\c MacroElement mapping resolves curvilinear domain boundaries
exactly, whereas the finite element mapping interpolates the global coordinates
between the coordinates of its nodes. The topology of \c MacroElements 
mirrors that of the associated (geometric) finite elements: 
For instance, the \c QMacroElement family is the counterpart 
of the \c QElement family of geometric finite elements. Both are 
templated by the spatial dimension, and the local coordinates 
(in their right-handed local coordinate systems) are in the range 
between -1 and +1. 


  The different-coloured, shaded regions in the above sketch
represent the four two-dimensional \c QMacroElements 
by which the \c FishDomain represents the fish-shaped
domain \f$ D_{fish} \f$. For instance, \c MacroElement 0
(shown in orange) represents the lower half of the fish's body; 
within this \c MacroElement, the curved "belly" is represented by the line
\f$ {\bf r}_{macro}(s_0,s_1=-1)\f$ for \f$s_0\in[-1,1]\f$; the lower
"jaw" is represented by \f$ {\bf r}_{macro}(s_0=-1,s_1)\f$ 
for \f$s_1\in[-1,1]\f$; etc.

To illustrate the use of \c MacroElements / \c Domains, 
the following code fragment (from  
<A HREF="../../../../src/meshes/fish_mesh.template.cc">fish_mesh.template.cc
</A>) demonstrates how the constructor of the original, 
non-refineable \c FishMesh assigns the nodal positions. 
Each of the \c FishDomain's four \c
QMacroElements is associated with one of the four finite elements in the mesh.
Since both types of elements are parametrised by the same
local coordinate systems, we determine the position of
the node that is located at \f$(s_0, s_1)\f$
(in the finite element's local coordinate system) from the 
corresponding \c MacroElement mapping, \f${\bf r}_{macro}(s_0, s_1)\f$:

\dontinclude fish_mesh.template.cc
\skipline Create elements
\until end of loop 

This technique ensures that the mesh's boundary nodes
are placed on the exact domain boundary when the mesh is created.
  
  To retain this functionality during the mesh adaptation,
each \c FiniteElement provides storage for a pointer to an 
associated \c MacroElement. By default,  the \c MacroElement pointer 
is set to \c NULL, indicating that the element is not associated 
with a \c MacroElement. In that case, the coordinates of newly created 
nodes are determined by interpolation from the father element, as discussed
above. If the \c MacroElement pointer is non-\c NULL, the refinement process
refers to the element's \c MacroElement representation to determine
the new nodal positions.

  To enable the mesh adaptation process to respect the domain's
curvilinear boundaries, each element in the coarse base mesh must 
therefore be given a pointer to its associated \c MacroElement,
e.g. by using the following loop:
\dontinclude fish_mesh.template.cc
\skipline Loop over all elements and set macro element pointer
\until }

Once the mesh is aware of the curvilinear boundaries, each level of 
mesh refinement produces a better representation of the 
curvilinear domain, ensuring the convergence to the exact solution:


\image html fish_poisson_adapt_soln.gif "Plot of the mesh adaptation with MacroElements " 
\image latex fish_poisson_adapt_soln.eps "Plot of the mesh adaptation with MacroElements " width=0.75\textwidth

The  results shown in this animation were computed with the demo code 
<A HREF="../../../../demo_drivers/poisson/fish_poisson2/fish_poisson_adapt.cc">
fish_poisson_adapt.cc</A> -- a simple modification
of the code  
<A HREF="../../../../demo_drivers/poisson/fish_poisson/fish_poisson.cc">
fish_poisson.cc</A> that we used in the
<A HREF="../../../../doc/poisson/fish_poisson/html/index.html">
earlier example</A>. The only difference between the two codes
is that in the present example, the \c FishDomain is discretised with 
four-node rather than nine-node \c RefineableQPoissonElements to 
highlight the inadequacy of the basic mesh refinement process.
Note that, as a result of lower accuracy of the four-node elements,
we require a much finer discretisation in the interior of the domain. 



<HR>

\section domain Part 2: How to represent domains with curvilinear boundaries by Domain and MacroElement objects

The above example demonstrated that "upgrading" existing 
meshes to \c RefineableMeshes that can be used with \c oomph-lib's
mesh adaptation procedures, can be achieved in two trivial steps:
-# Associate each \c RefineableQElement with a \c QuadTree --
   this can done completely automatically by calling the function
   \c RefineableQuadMesh::setup_quadtree_forest().
-# If the problem's domain has curvilinear boundaries, associate 
   each \c RefineableQElement with a \c MacroElement -- 
   defined in the \c Domain object that provides an     
   analytical representation of the domain. 
. 
While this looks (and indeed is) impressively simple, we still have to 
explain how to create \c Domain objects. We start by introducing
yet another useful \c oomph-lib class, the \c GeomObject.

\subsection geom_object The geometric object, GeomObject

As the name suggests, \c GeomObjects are \c oomph-lib objects that 
provide an analytical description/parametrisation of geometric
objects. Mathematically, \c GeomObjects define a mapping from 
a set of "Lagrangian" (intrinsic) coordinates to the global 
"Eulerian" coordinates of the object. The number of Lagrangian
and Eulerian coordinates can differ. For instance, the unit circle,
centred at the origin may be parametrised by a single coordinate, \f$ \xi
\f$ (representing the polar angle), as
\f[
{\bf r}_{circle} = 
\left( 
\begin{array}{c}
\cos \xi \\
\sin \xi
\end{array}
\right),
\f] 
while a 2D disk may be parametrised by two coordinates 
 \f$ \xi_1 \f$ and \f$ \xi_2 \f$ (representing the radius and the
polar angle, respectively) as
\f[
{\bf r}_{disk} = 
\xi_1 \left(
\begin{array}{c}
\cos \xi_2 \\
\sin \xi_2
\end{array}
\right).
\f] 
All specific \c GeomObjects must implement the pure virtual function
\c GeomObject::position(...) which computes the Eulerian
position vector \f$ {\bf r}\f$ as a function of the 
(vector of) Lagrangian coordinates \f$ {\bf \xi}\f$.
(The \c GeomObject base class also provides interfaces for
a multitude of other functions, such as functions that compute
the spatial and temporal derivatives of the position vector. 
These functions are implemented as "broken" virtual functions and 
their implementation is optional; see the 
<A HREF="../../../../doc/poisson/two_d_poisson_flux_bc/html/index.html">
earlier example</A> for a discussion of "broken" virtual functions.)

  Here is a complete example of a specific \c GeomObject:

\dontinclude fish_domain.cc
\skipline start_of_unit_circle
\until end of unit circle


[The dummy time-dependent version of the \c position(...) function is
required to stop the compiler from complaining about "only partially 
overridden" virtual functions]. 

\subsection domain_objects Domains

\c GeomObjects provide a natural way of representing a \c Domain's
curvilinear boundaries. For instance, the fish's body in  \f$
D_{fish}\f$ is bounded by two circular arcs. These may be
represented by \c GeomObjects of type \c Circle -- a slight generalisation
of the \c UnitCircle class shown above. The \c FishDomain constructor
therefore takes a pointer to a 2D \c GeomObject and the "start" 
and "end" values of the Lagrangian coordinate along this object. 
The \c GeomObject represents the curvilinear boundary 
of the fish's  (upper) body and the two coordinates represent the 
Lagrangian coordinates of the "nose" and the "tail" 
on this \c GeomObject, as shown in this sketch:


\image html fish_domain.gif "The fish-shaped domain and its MacroElement-based representation by the FishDomain object. The arrows show the orientation of the MacroElements' local coordinate systems. " 
\image latex fish_domain.eps "The fish-shaped domain and its MacroElement-based representation by the FishDomain object. The arrows show the orientation of the MacroElements' local coordinate systems. " width=0.75\textwidth


To construct a \c FishDomain whose curvilinear boundaries
are arcs of unit circles, centred at \f$ (x_0,x_1) = (1/2 , 0) \f$
we create a \c GeomObject of type \c Circle, passing the
appropriate parameters to its constructor:

\code

 // Fish back is a circle of radius 1, centred at (0.5,0.0)
 double x_c=0.5;
 double y_c=0.0;
 double r_back=1.0;
 GeomObject* back_pt=new Circle(x_c,y_c,r_back);
 
\endcode

Next, we pass the (pointer to the)  \c Circle object
to the constructor of the \c FishDomain, locating
the "nose end" of the fish's back at \f$\xi=2.4\f$ 
and  its "tail end" at \f$\xi=0.4\f$:

\code

 double xi_nose=2.6; 
 double xi_tail=0.4; 
 Domain* domain_pt=new FishDomain(back_pt,xi_nose,xi_tail);
 
\endcode

 
 To see how this works internally, let us have a look at the
\c FishDomain constructor. The constructor stores the 
pointer to the fish's "back", and the start and end values of the
Lagrangian coordinates in the private data members \c Back_pt, 
\c Xi_nose and \c Xi_tail. Next we set some additional parameters,
that define the geometry (the mouth is located at the origin; 
the fin is a vertical line at \f$ x=1.7\f$, ranging from 
\f$ y=-0.9\f$ to  \f$ y=+0.9\f$). Finally, we allocate storage for 
the four \c MacroElements and build them. Note that the constructor 
of the \c MacroElement takes a pointer to the Domain, and the
\c MacroElement's number within that \c Domain:

 
\dontinclude fish_domain.h
\skipline Constructor: Pass pointer to GeomObject that 
\until end of constructor

Most of the remaining public member functions are equally
straightforward. We provide various access functions to the geometric parameters
such as \c X_mouth, etc -- we will not list these explicitly.
All the "real work" is done in the implementation of the pure 
virtual function \c Domain::macro_element_boundary(...). Given 
- the number of the \c MacroElement in its \c Domain
- the direction of its boundary (N[orth], S[outh], E[ast], W[est],
  enumerated in the namespace \c QuadTreeNames)
.
this function must compute the vector \f$ {\bf r}(\zeta)\f$ to the 
\c MacroElement's boundary.
Here \f$ \zeta \in [-1,1]\f$ is the 1D coordinate along the element
boundary, aligned with the direction of the \c 
MacroElement's 2D coordinates \f$(s_0,s_1) \f$ as indicated in this sketch:
 
\image html macro_element_sketch.gif "Sketch illustrating the parametrisation of the MacroElement's four boundaries. " 
\image latex macro_element_sketch.eps "Sketch illustrating the parametrisation of the MacroElement's four boundaries. " width=0.75\textwidth


Since the shape of the domain can evolve in time, the full interface 
for the function includes an additional parameter, \c t, 
which indicates the (discrete) time level at which the domain shape
is to be evaluated. If \c t=0 the function computes the domain shape
at the current time; if \c t>0 it computes the shape at the \c t -th 
previous timestep. 
(<A HREF="../../../unsteady_heat/two_d_unsteady_heat_ALE/html/index.html">
Another example</A> in which we solve the unsteady
heat equation in a moving domain, provides a a more detailed discussion of 
this aspect.) Here is the full
interface for the \c FishDomain::macro_element_boundary(...) function:

\skipline Vector repre
\until r);

The implementation of this function is the only tedious task that
needs to be performed by the "mesh writer". Once 
\c Domain::macro_element_boundary(...) 
is implemented, the \c Domain's constituent \c MacroElements
can refer to this function to establish the positions of their
boundaries (recall that we passed the pointer to the \c Domain and
the \c MacroElement's number in the \c Domain to the \c MacroElement
constructor). The \c MacroElement::macro_map(...) functions
interpolate the position of the \c MacroElement's boundaries into their
interior. 


To illustrate the general procedure, here is the complete listing of 
the \c FishDomain::macro_element_boundary(...) function. 
The function employs switch statements to identify the private member functions
that provide the parametrisation of individual \c MacroElement
boundaries. Some of these functions are listed below.

\dontinclude fish_domain.h
\skipline start_of_macro_element_boundary
\until end of macro_element_boundary


Here are a few of the private member functions that define individual 
\c MacroElement boundaries:

- The N[orthern] boundary
  of macro element 2 (which represents the upper body) 
  coincides with the domain boundary that is parametrised by the
  geometric object pointed to by \c Back_pt. The function 
  translates the coordinate \f$ \zeta \in [-1,1] \f$ to the
  Lagrangian coordinate \f$ \xi \in [\xi_{nose}, \xi_{tail}] \f$
  along the geometric object. We use this Lagrangian coordinate
  to obtain the position vector to the domain boundary via a call to the
  \c GeomObject::position(...) function of the geometric object 
  pointed to by \c Back_pt:
\dontinclude fish_domain.h
\skipline start_of_r_upper_body_N
\until  end of r_upper_body_N 


- The E[astern] boundary of macro element 2 is a straight vertical
  line from the "tail end" of the curved fish back to the
  x-axis: 
\dontinclude fish_domain.h
\skipline start_of_r_upper_body_E
\until  end of r_upper_body_E

- The S[outhern] boundary of macro element 2 is a straight horizontal
  line from the "mouth" to the end of the
  body: 
\dontinclude fish_domain.h
\skipline start_of_r_upper_body_S
\until  end of r_upper_body_S


- The W[estern] boundary of macro element 2 is a straight
  line from the "mouth" to the "mouth" end of the
  curved upper boundary of the body:
\dontinclude fish_domain.h
\skipline start_of_r_upper_body_W
\until  end of r_upper_body_W

- The S[outhern] boundary of macro element 0 (which represents the
  lower body) is simply a reflection of the N[orthern] boundary
  of macro element 2:
\dontinclude fish_domain.h
\skipline  Southern boundary of lower body macro element
\until  }

- etc.

Tedious? Yes! Rocket Science? No!
 
<HR>
<HR>

\section com Further comments

\subsection nod_update Node updates in response to changes in the Domain shape.
You may have noticed that, even though we introduced
\c MacroElements in the context of adaptive mesh refinement, 
the pointer to a refineable element's \c MacroElement is 
stored in the \c FiniteElement, rather than the (derived) \c
RefineableElement class, suggesting that \c MacroElements
have additional uses outside the context of mesh adaptation. 
Indeed, the code fragment that illustrated the use of \c Domains 
and \c MacroElements during mesh generation, was taken from 
the constructor of the \e non-refineable \c FishMesh,
rather than its adaptive counterpart. During the mesh generation
process, the \c FiniteElement's \c MacroElement representation 
was used to determine the position of its \c Nodes within the \c Domain. 
The same procedure can be employed to \e update the nodal positions in 
response to changes in the domain shape. This is implemented, generically,
in the function

\code
Mesh::node_update()
\endcode

This function loops over all elements in a \c Mesh and updates
their nodal positions in response to changes in the domain boundary.
(If the \c Mesh's constituent elements's are not associated with 
\c MacroElements and if the \c Mesh does not implement the node update
by other means, this function does not change the mesh.)

The following code fragment illustrates the trivial modifications
to the driver code required to compute the solution of Poisson's 
equation in fish-shaped  domain of various widths. We simply
change the position of the \c GeomObject that specifies the
curvilinear boundary (by changing the position of the circle's
centre), call the \c Mesh::node_update() function, and recompute
the solution.


\dontinclude fish_poisson_node_update.cc
\skipline start_of_main
\until end of main

The 
 <A HREF="../../../../demo_drivers/poisson/fish_poisson2/fish_poisson_node_update.cc">rest
 of the code</A> remains unchanged. Here is a plot of the solution for 
various widths of the domain (computed with nine-node elements).

\image html fish_poisson_node_update.gif "Adaptive solution of Poisson's equation in fish-shaped domains of varying width. " 
\image latex fish_poisson_node_update.eps "Adaptive solution of Poisson's equation in fish-shaped domains of varying width. " width=0.75\textwidth

<HR>
<HR>

\subsection boundary_coords Good practice: Storing boundary coordinates

The above example demonstrated how the representation of curvilinear 
domain boundaries by \c GeomObjects
allows \c oomph-lib's mesh generation and adaptation procedures to 
place nodes on these boundaries. We note that the 
Lagrangian coordinate(s) that parametrise(s) the relevant \c
GeomObjects also provide a parametrisation of the
corresponding domain boundaries. In certain applications (such as 
free-boundary or fluid-structure interaction problems) it is useful
to have direct access to these boundary coordinates. For this purpose
the \c Node class provides the function
\code
Node::set_coordinates_on_boundary(const unsigned& b, 
                                  const Vector<double>& xi);
\endcode
which allows the mesh writer to store the (vector of) boundary
coordinates that a given (\c Boundary)\c Node is located at. The argument \c b 
specifies the number of the mesh boundary, reflecting the fact that nodes may
be located on multiple domain boundaries, each of which is likely to have
a different set of surface coordinates. [\b Note: The function is
implemented as a broken virtual function in the \c Node base class.
The actual functionality to store boundary coordinates is only
provided in (and required by) the derived \c BoundaryNode class.]

Since the storage of boundary coordinates is optional, 
the \c Mesh base class provides a protected vector of bools,
\code
std::vector<bool> Mesh::Boundary_coordinate_exists;
\endcode
that indicates if the boundary coordinates have been stored for
all \c Nodes on a specific mesh boundary. This vector
is resized and its entries are initialised to \c false, when
the number of mesh boundaries is declared with a call
to \c Mesh::set_nboundary(...). If, during mesh refinement, 
a new  \c BoundaryNode is created on the mesh's boundary \c b, 
its boundary coordinates are computed by interpolation from the
corresponding values at the nodes in the father element, if 
\c Mesh::Boundary_coordinate_exists[b] has been set to \c true.

We regard it as good practice to set boundary coordinates for
all \c BoundaryNodes that are located on curvlinear mesh boundaries. 
The source code
<A HREF="../../../../src/meshes/fish_mesh.template.cc">fish_mesh.template.cc</A>
for the refineable FishMesh illustrates the methodology.


<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/fish_poisson2/">
demo_drivers/poisson/fish_poisson2/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/fish_poisson2/fish_poisson_adapt.cc">
demo_drivers/poisson/fish_poisson2/fish_poisson_adapt.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Demo problem: Solution of a "free-boundary" Poisson problem -- a simple model for "fluid"-structure interaction.

In this example we shall consider our first (toy!) interaction
problem. The problem combines two single-physics problems, studied
in earlier examples, and combines them into a coupled free-boundary problem.

- In <A HREF="../../../poisson/fish_poisson/html/index.html">
  one of our first examples</A> we demonstrated the solution of
  Poisson's equation in a fish-shaped domain, \f$ D_{fish} \f$, in 
  which the curvilinear
  upper and lower boundaries of the fish's body were given by
  circular arcs which we represented by \c GeomObjects.
  Given the position of the two circular arcs, whose centres are 
  located at \f$ (X_c, \pm Y_c) \f$, the single-physics Poisson code
  computes the solution of Poisson's equation in the corresponding
  domain.  
  <A HREF="../../../poisson/fish_poisson2/html/index.html#com">We 
  have already demonstrated</A> how \c oomph-lib's
  \c MacroElement / \c Domain representation of  \f$ D_{fish} \f$
  allows an update of the nodal positions in response to changes
  in the domain boundary by a simple call to \c Mesh::node_update().
  The example code therefore allowed us to compute the solution at a control 
  node, \f$ u_{ctrl}\f$, as a function of the "height" of the domain,
  characterised by \f$ Y_c \f$. 
  \n\n
- In <A HREF="../../../interaction/circle_as_element/html/index.html">
  another example</A>, we demonstrated how to solve a (trivial)
  solid mechanics problem: The vertical displacement of a circular ring
  (represented by a \c GeomObject) that is mounted on an elastic 
  foundation of  spring stiffness \f$ k \f$. The example code allowed 
  us to compute the displacement of the ring, characterised by
  \f$ Y_c \f$, as a function of the load \f$ f \f$
  acting on it.   
.
We will now consider the coupled problem obtained by using the
solution of Poisson's equation at the control node, \f$ u_{ctrl} \f$, 
as the "load", \f$ f \f$, that acts on the two circular arcs that define the
curvilinear boundaries of \f$ D_{fish} \f$. The resulting coupled
problem is sketched in the figure below.  While this problem  is 
obviously somewhat artificial, it has many of the key features that
arise in genuine fluid-structure interaction problems. In particular, the 
displacement of the domain boundary is driven by the solution 
of the "bulk" (Poisson) equations, just as the deformation of an elastic
structure in an FSI problem is driven by the fluid pressure and shear
stresses, i.e. quantities that are derived from the solution of the
Navier-Stokes equations in the "bulk" domain.


\image html fish_fsi_sketch.gif "Sketch of the two individual single-physics problems (top) and the coupled problem (bottom). " 
\image latex fish_fsi_sketch.eps "Sketch of the two individual single-physics problems (top) and the coupled problem (bottom). " width=0.75\textwidth
 

The two single-physics problems involve two uncoupled sets of 
equations and unknowns:
- The residual vector of the Poisson elements in the single-physics
  Poisson problem depends on the nodal values in the "bulk" mesh. 
  These nodal values are the only unknowns in the problem since the
  position of the domain boundary, and hence the position of the nodes 
  are fixed. An update of the nodal positions in response to any
  changes in the domain boundary is a mere pre-processing step, to be
  performed just once, before computing the solution.
- The residual vector of the \c ElasticallySupportedRingElement 
  depends on the position of ring's centre, \f$ Y_c \f$, which
  is the only unknown in the problem as the load on the ring
  is fixed. 
.
The coupling between the two single-physics problem introduces
additional dependencies:
- The residuals of the Poisson elements also depend on the nodal
  positions which in turn depend (via the  \c MacroElement/Domain - based
  node-update function) on the the position of the domain boundary. 
  The boundary position is controlled by the 
  \c ElasticallySupportedRingElement's geometric \c Data, which 
  stores the value of \f$ Y_c \f$. 
- The residual vector of the \c ElasticallySupportedRingElement also
  depends on the load, which is now given by the unknown nodal value
  at a control node in the "bulk" mesh.
.

We note that most of the methodology required to solve this coupled
problem is already available: 
- The \c MacroElement/Domain representation of the Mesh makes it
  possible to update the nodal positions in the bulk mesh 
  in response to changes in the shape/position of the curvilinear 
  domain boundary.
- Multiple inheritance allows the \c ElasticallySupportedRingElement
  to act as a \c GeomObject (a role in which it can be used to 
  parametrise the unknown curvilinear domain boundary) and 
  as a \c GeneralisedElement (a role in which its unknown geometric \c
  Data value, \f$ Y_c \f$, can be determined as part of the overall
  solution).
- The load  \f$ f \f$ on the \c ElasticallySupportedRingElement
  is stored in the element's external \c Data, and derivatives of the
  element's residual vector with respect to \f$ f \f$ are 
  automatically taken into account when the element's 
  Jacobian matrix is computed.
.
The only interaction that still has to be incorporated into the problem
formulation is the dependence of the Poisson element's 
residual vectors on the geometric \c Data in the 
\c ElasticallySupportedRingElement. This interaction arises
through the \c MacroElement/Domain - based node-update function
which translates changes in the \c GeomObject's geometric \c Data 
into changes in the nodal positions. Such dependencies may 
be added to \e any existing element by "wrapping" the element
into the templated wrapper class \c MacroElementNodeUpdateElement
which has the following inheritance structure:

\code
template<class ELEMENT>
class MacroElementNodeUpdateElement : public virtual ELEMENT, 
                                      public virtual MacroElementNodeUpdateElementBase
\endcode

An element of type 
\c MacroElementNodeUpdateElement<ELEMENT> is an element of type \c 
ELEMENT, and inherits the additional functionality provided by the
\c  MacroElementNodeUpdateElementBase base class. The most important
additional functionality provided by this class is the ability to 
add the values stored in 
the geometric \c Data of associated \c GeomObjects to the element's 
list of unknowns. Once added, the derivatives of the element's residual 
vector with respect to these additional unknowns are automatically 
included into the element's Jacobian matrix. This is achieved by 
overloading the \c ELEMENT::get_jacobian(...) function
and evaluating the additional derivatives by finite differencing. 
See \ref comments for details on the implementation. 


The solution of the coupled problem therefore only requires a 
few trivial changes to the single-physics (Poisson) code:
-# The element type used for the solution of the "bulk" equations 
   must be changed to its "wrapped" counterpart,
   as discussed above. For instance, if the
   single-physics code used a nine-node refineable Poisson element
   of type \c RefineableQPoissonElement<2,3>, the coupled problem 
   must be discretised by elements of type
   \c MacroElementNodeUpdateElement<RefineableQPoissonElement<2,3> > 
   (Yes, it's a bit of a mouthful...).
-# The "bulk" mesh must be "upgraded" (again via multiple
   inheritance) to a Mesh that is derived from the
   \c MacroElementNodeUpdateMesh base class. 
-# A vector of pointers to those \c GeomObjects that are involved
   in an element's \c MacroElement/Domain - based node update
   operation must be passed to the elements. (This is done most easily 
   in the constructor of the "upgraded" mesh.) The geometric \c Data
   contained in these \c GeomObjects is then automatically included
   in the elements' list of unknowns.
-# The Mesh's \c node_update() function must be executed whenever the
   Newton method has changed the values of the unknowns: This is
   because changing a value that is stored in a \c GeomObject's 
   geometric \c Data does not automatically update the positions 
   of any dependent nodes. This is done most easily be including the
   \c node_update() function into the 
   \c Problem::actions_before_newton_convergence_check() function; we refer to
   <A HREF="../../../order_of_action_functions/html/index.html">
   another document</A> for a more detailed discussion of the order in 
   which the various "action" functions are called by \c oomph-lib's 
   Newton solver. 
.

<HR>
<HR>


\section results Results

The animation below shows the results of a spatially-adaptive solution
of Poisson's equations in the fish-shaped domain, for a variety of
domain "heights". This animation was produced with the 
<A HREF="../../../poisson/fish_poisson/html/index.html">single-physics
Poisson solver discussed in an earlier example</A>. 


\image html elastic_fish.gif "Spatially adaptive solution of Poisson's equation in a fish-shaped domain for various `widths' of the domain. " 
\image latex elastic_fish.eps "Spatially adaptive solution of Poisson's equation in a fish-shaped domain for various `widths' of the domain. " width=0.75\textwidth
 

An increase in the height of the domain increases the amplitude of the 
solution. This is reflected by the red line in the figure below 
which shows a plot of \f$ u_{ctrl} \f$ as a function of \f$ Y_c \f$. 
The green marker shows the solution of the coupled problem for a spring 
stiffness of \f$ k=1 \f$. For this value of the spring stiffness, 
the solution of the coupled problem should be (and indeed is) located
at the intersection of the  curve \f$ u_{ctrl}(Y_c) \f$ with the diagonal,
\f$ u_{ctrl} = Y_c \f$, shown by the dashed blue line.

\image html trace.gif "Solution of Poisson's equation at a control node as a function of the `height' of the domain. " 
\image latex trace.eps "Solution of Poisson's equation at a control node as a function of the `height' of the domain. " width=0.75\textwidth
 
<HR>   
<HR>

\section impl Implementation in oomph-lib
The sections below provide the usual annotated listing of the
driver code. We stress that only a few trivial changes are required
to incorporate the presence of the free boundary into the
existing single-physics code:
- \ref mesh : Upgrading the \c RefineableFishMesh via multiple inheritance. 
- \ref main : Changing the element type for the solution of the
  Poisson equation.        
- \ref constructor : Storing the element that represents the free boundary
  in a (sub-)mesh.
- \ref problem : Implementing the function 
  \c Problem::actions_before_newton_convergence_check() to update the 
  nodal positions after each Newton step.
.

<HR>   
<HR>


\section global Global parameters and functions
The namespace \c ConstSourceForPoisson defines 
the constant source function, exactly as in the 
<A HREF="../../../poisson/fish_poisson/html/index.html#global">corresponding 
single-physics code.</A>


\dontinclude macro_element_free_boundary_poisson.cc
\skipline start_of_namespace
\until end of namespace 


<HR>
<HR>

\section mesh The Mesh
Meshes that are to be used with \c MacroElementNodeUpdateElements 
should be derived (typically by multiple inheritance) from the  
\c MacroElementNodeUpdateMesh class. This class overloads the 
generic \c Mesh::node_update() function and ensures that the node update
is performed by calling the \c node_update() function of the
\c Mesh's constituent nodes, rather than simply updating their
positions, using the \c FiniteElement::get_x(...) function. 
The overloaded version is not only more efficient but also ensures
that any auxiliary node update functions (e.g. functions that update the
no-slip condition on a moving fluid node on a solid boundary)
are performed too. 

 In our driver code we add the additional functionality
provided by the \c MacroElementNodeUpdateMesh class to 
the \c RefineableFishMesh class used in the 
<A HREF="../../../poisson/fish_poisson/html/index.html">
single-physics Poisson problem considered earlier</A>. 


\skipline start_of_mesh
\until {

The constructor calls the constructors of the underlying
\c RefineableFishMesh. [Note the explicit call to the \c FishMesh 
constructor prior to calling the constructor of the 
\c RefineableFishMesh. Without this
call, only the default (argument-free) constructor of the
\c FishMesh would be called! Consult your favourite C++ book
to check on constructors for derived classes if you don't understand
this. We recommend <A HREF="http://www.math.wayne.edu/~yang/">Daoqi Yang's</A> 
brilliant book 
<A HREF="http://www.springeronline.com/sgw/cda/frontpage/0,11855,4-40007-22-2105335-0,00.html?changeHeader=true"> C++ and Object-Oriented Numeric 
  Computing for Scientists and Engineers.</A>)


\until {

To activate the \c MacroElementNodeUpdateElement's ability to
automatically compute the derivatives of the residual vectors
with respect to the geometric \c Data that determines
its nodal positions, we must pass the pointers to the
\c GeomObjects that are involved in the element's \c MacroElement - based 
node-update to the elements. In general, an element's node-update will
be affected by multiple \c GeomObjects therefore the 
\c set_node_update_info(...) function expects a vector of
pointers to \c GeomObjects. In the present example, only a single
\c GeomObject (the \c GeomObject that represents the fish's curved "back")
determines the nodal position of all elements:
 
\until end of constructor
 
The destructor can remain empty but we provide a final overload
for the \c Mesh's \c node_update() function to avoid any ambiguities
as to which one is to be used.
\until // end of mesh class

<HR>
<HR>

\section main The driver code
The driver code is very simple: We build the problem with the
"wrapped" version of the refineable quadrilateral nine-node 
Poisson element. Since the initial mesh is very coarse we 
perform two uniform mesh refinements before solving the problem 
with automatic spatial adaptivity, allowing for up to two further
mesh adaptations. 

\dontinclude macro_element_free_boundary_poisson.cc
\skipline start_of_main
\until end of main 

<HR> 
<HR>
 

\section problem The problem class
Apart from a few trivial additions, 
the problem class is virtually identical to that used in the 
<A HREF="../../../poisson/fish_poisson/html/index.html">single-physics 
Poisson problem</A>.
The most important addition to the single-physics problem class is 
the function \c Problem::actions_before_newton_convergence_check() which
updates the nodal positions in the "bulk" Poisson mesh following an update
of the geometric \c Data that controls the position of the
curvilinear domain boundary; we refer to
<A HREF="../../../order_of_action_functions/html/index.html">
another document</A> for a more detailed discussion of the order in 
which the various "action" functions are called by \c oomph-lib's 
Newton solver. 


\dontinclude macro_element_free_boundary_poisson.cc
\skipline start_of_problem_class
\until end of problem class
 
<HR>
<HR>

\section constructor The Problem constructor
We start by creating the \c GeomObject/GeneralisedElement that
will represent the unknown curvilinear domain boundary and 
pass it (in its role as a \c GeomObject) to the constructor
of the bulk mesh. We then add the pointer to the bulk mesh to the 
\c Problem's collection of submeshes and create an error estimator
for the adaptive solution of the Poisson equation.


\skipline start_of_constructor
\until Z2


Next we store the pointer to the \c ElasticallySupportedRingElement
in its own Mesh and add it to the \c Problem's collection of submeshes 
before building the \c Problem's global \c Mesh from its two submeshes:

\skipline Build mesh
\until build_global_mesh()

We choose the central node in the Poisson mesh as the control node
and use it (in its role as \c Data) as the "load" for the
\c ElasticallySupportedRingElement.

\skipline Choose
\until set_load_pt

Finally,  we pin the nodal values on all boundaries, apply the 
homogeneous Dirichlet boundary conditions,
pass the pointer to the source function to the elements,
and set up the equation
numbering scheme. 

\skipline Set the 
\until end of constructor

<HR>
<HR>

\section doc Post-processing
The post-processing routine writes the computed result to 
an output file.

\skipline start_of_doc
\until end of doc



<HR>
<HR>

\section comments Comments 

A more detailed description of the theory and the implementation 
can be found in the paper
- Heil, M. & Hazel, A. L. "<TT>oomph-lib</TT> -- An 
  <EM>O</EM>bject-<EM>O</EM>riented <EM>M</EM>ulti-<EM>Ph</EM>ysics
  Finite-Element <EM>Lib</EM>rary". In: <EM>Fluid-Structure
  Interaction</EM>, Editors: M. Schafer and H.-J. Bungartz. 
  Springer (Lecture Notes on Computational Science and Engineering 53, 2006), 
  (32 pages) 
  <A HREF="http://www.maths.man.ac.uk/~mheil/MATTHIAS/ABSTRACTS/HeilHazelOomph2006.html">(abstract)</A> 
  <A HREF="http://www.maths.man.ac.uk/~mheil/MATTHIAS/PDF/oomph_for_www.pdf">
  (pdf preprint)</A>.
. 
and in this talk:
- Heil, M. & Hazel, A. L.  "An object-oriented approach to the 
  evaluation of the `shape derivatives' in monolithic fluid-structure 
  interaction solvers". 7th World Congress on Computational Mechanics,
  LA, July 2006.
  <A HREF="http://www.maths.man.ac.uk/~mheil/oomph_lib_additional_material/LA_talk_2006/LA_talk.pdf">(pdf)</A>.
.
The following subsections provide a brief description of the
main features.

<HR>

\subsection sparse_node_updates Sparse node updates

The key feature of our implementation which allows the efficient
computation of the "shape derivatives" is 
the ability of \c MacroElementNodeUpdateNodes (discussed in more detail below) 
to "update their own position" in response to changes in
shape/position of the domain boundary.
This capability is demonstrated in the following simple example code.

We start by building the Mesh as before

\dontinclude doc_sparse_macro_node_update.cc
\skipline start_of_main
\until MacroElementNodeUpdateRefineableFishMesh<ELEMENT>(Fish_back_pt);

and document the mesh (i.e. the shape of its constituent finite elements
and the nodal positions):

\skipline Number of plot
\until count++;

Next, we "manually" increment \f$ Y_c \f$, i.e. the y-coordinate of
the centre of the circular arc that defines the upper curvilinear
boundary of the fish mesh. 

\skipline Increment
\until y_c()

This step mimics the incrementation of one
of the \c Problems's unknowns (recall that in the free-boundary
problem considered above,\f$ Y_c \f$ has to be determined as part 
of the solution!) during the finite-difference based computation of
the shape derivatives. 

For meshes that are not derived from the
\c MacroElementNodeUpdateMeshBase  class, the only way to update the
nodal positions in response to a change in the boundary position,
is to call the \c Mesh::node_update() function. This updates the
position of \e all nodes in the mesh -- a very costly operation. 

Meshes that are derived from the \c MacroElementNodeUpdateMeshBase
class contain \c MacroElementNodeUpdateNodes which can update their
own position, as shown here:

\until end of main

We note that the \c Node::node_update() function is defined as an 
empty virtual function in the \c Node base class, 
indicating that "normal" \c Nodes cannot "update their own position". 
The function is overloaded in the \c MacroElementNodeUpdateNode class,
details of which are given below. Overloaded versions of this function 
also exist in various other derived \c Node classes (such as as the 
\c AlgebraicNodes and the \c SpineNodes) for which algebraic node update
operations are defined. 

Here is an animation that illustrates how the successive update of the
individual nodal positions in response to the change in the boundary
position gradually updates the entire mesh.

\image html sparse_node_update.gif "Illustration of the sparse node-update procedure. " 
\image latex sparse_node_update.eps "Illustration of the sparse node-update procedure. " width=0.75\textwidth
 

<HR>

\subsection how_it_works How it works
The implementation employs three key components:
- <B><TT>MacroElementNodeUpdateNodes</TT></B> are derived from the 
  \c Node base class. Their main purpose is to provide the 
  \c MacroElementNodeUpdateNodes::node_update() function which updates 
  the nodal position in response to changes in the domain
  boundary. This capability was demonstrated above and is achieved
  by allowing the \c MacroElementNodeUpdateNodes to store
  a pointer to the \c MacroElementNodeUpdateElement that determines
  its position (using its own \c MacroElement - based representation) 
  and its local coordinates in that element.
  \n\n \c MacroElementNodeUpdateNodes also store a function pointer to an
  auxiliary node update function that allows additional tasks
  to be performed whenever a node update is performed. This is 
  useful, e.g. in unsteady fluid-structure interaction problems 
  in which a change in the position of nodes that are located on 
  a no-slip boundary also requires an update of the
  fluid velocities at that node. By default, the function pointer
  is initialised to NULL, indicating that no auxiliary node update
  functions have to be executed.
  \n\n Finally, the \c MacroElementNodeUpdateNodes store pointers to
  the \c GeomObjects that affect their node update.
  While this information is not required by the node update
  function itself, it must be available to correctly set up the
  equation numbering scheme in the presence of hanging nodes.
  (Details are too messy to explain here but it's true!).
  \n\n
- The <B><TT>MacroElementNodeUpdateElement<ELEMENT></TT></B> class was already
  discussed in the main part of this document. These elements
  "wrap around" the element specified by the template argument,
  \c ELEMENT, overload some of its member functions and add 
  some new ones.
  \n \n Overloaded functions include: \n\n
  - The \c FiniteElement::construct_node(...) functions
    create an element's local \c Nodes. This is overloaded by a version
    that creates \c MacroElementNodeUpdateNodes instead. \n\n
  - The functions \c GeneralisedElement::get_jacobian(...) and 
    \c GeneralisedElement::fill_in_contribution_to_jacobian(...) 
    are overloaded by versions that add the shape derivatives
    to the Jacobian matrices computed by the underlying ELEMENT. \n\n
  - Similarly, the function 
    \c FiniteElement::assign_all_generic_local_eqn_numbers()
    is overloaded to add the unknowns associated with the node update
    functions into the element's equation numbering scheme. \n\n
  .
  Additional member functions are provided to specify (and access) the
  \c GeomObjects that affect an element's \c MacroElement - based node 
  update. 
  Full details may be found in the 
  <A HREF="../../../the_data_structure/html/index.html">"bottom
  up"</A> discussion of \c oomph-lib's data structure.
  \n\n
- Finally, the <B><TT>MacroElementNodeUpdateMeshBase</TT></B> class
  overloads the \c Mesh::node_update() function to ensure that
  node updates are performed node-by-node, using the
  \c MacroElementNodeUpdateNode::node_update() function. This ensures
  that the node update not only updates the nodal positions but also
  executes any auxiliary update functions. 
.

<HR>

\subsection fsi The method also works for non-"toy" problems!

The above example demonstrated how easy it is to "upgrade" a driver
code for the solution of a single-physics problem to a
fluid-structure-interaction-like free-boundary problem. It is
important to stress that the methodology employed in our 
"toy" free-boundary problem can also be used for genuine 
fluid-structure interaction problems. For instance, the 
driver code for the simulation of 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
2D unsteady finite-Reynolds number flow in a channel with an
oscillating wall whose motion is prescribed</A> can easily be
extended to a driver code for the corresponding 
<A HREF="../../../interaction/fsi_collapsible_channel/html/index.html">fluid-structure interaction problem in which the
wall is replaced by a flexible membrane that is loaded by the fluid
traction.</A>

 

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/free_boundary_poisson/">
demo_drivers/interaction/free_boundary_poisson/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/free_boundary_poisson/macro_element_free_boundary_poisson.cc">
demo_drivers/interaction/free_boundary_poisson/macro_element_free_boundary_poisson.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


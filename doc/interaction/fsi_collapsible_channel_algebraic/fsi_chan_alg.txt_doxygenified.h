/**

\mainpage Flow in a 2D collapsible channel revisited -- sparse algebraic node updates

In an <A HREF="../../fsi_collapsible_channel/html/index.html">earlier
example,</A> we demonstrated how the \c MacroElement/Domain - based
node-update procedure that we originally developed for problems 
with moving, curvilinear domain boundaries may also be used in
fluid-structure interaction problems in which the position of the 
domain boundary has to be determined as part of the solution. 
We demonstrated that the driver code for the coupled multi-physics
problem was a straightforward combination of the driver codes
for the two constituent single-physics problems. The two key 
steps required to couple the two single-physics codes were:
-# Recast the wall mesh to a \c GeomObject, using the \c
   MeshAsGeomObject class. This class turns an existing solid
   mechanics mesh
   into a "compound" \c GeomObject in which material points on the      
   wall are identified by their Lagrangian coordinate, \f$ \xi \f$ ,
   which doubles as the \c GeomObject's intrinsic coordinate,
   \f$ \zeta \f$ . \n\n
-# Use the "compound" \c GeomObject 
   to represent the moving boundary of the fluid mesh. "Upgrade" the
   fluid elements (of type \c FLUID_ELEMENT, say) to the "wrapped"
   version \c MacroElementNodeUpdateElement<FLUID_ELEMENT> to allow the
   the node-update to be performed node-by-node, and to automatically
   evaluate the "shape derivatives" -- the derivatives of the fluid
   equations with respect to the (solid mechanics) degrees of freedom
   that determine their nodal positions. \n\n
.

While the implementation of these steps is very straightforward, we
pointed out that the resulting code was not particularly efficient 
as the fluid-node update is not as sparse as it could (should!)
be: Since it is impossible to distinguish between the various sub-objects in 
the "compound" \c GeomObject, we can do no better than to assume the 
worst-case scenario, namely that the positional degrees 
of freedom of all \c SolidNodes in the wall mesh potentially affect 
the nodal position in all fluid elements. This dramatically increases 
the size of the elemental Jacobian matrices, and creates many
nonzero entries in the off-diagonal blocks in the global Jacobian 
matrix.

\section alg_sparse Sparse algebraic node updates in FSI problems

To avoid this problem we need a node-update strategy in which the
position of each fluid node is determined by only a small number
of solid mechanics degrees of freedom. The algebraic node-update
strategy discussed in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
non-FSI version of the collapsible channel problem,</A> 
provides an ideal framework for this, as 
it allows each node to update its own position, using a node-specific
update function. Recall that in the \c AlgebraicMesh - version of the
\c CollapsibleChannelMesh, each \c AlgebraicNode stored a pointer
to the (single) \c GeomObject that represented the moving curvilinear
boundary, and the Lagrangian coordinate of a reference point on this
\c GeomObject. The node's node-update function then placed the node 
at a fixed vertical fraction on the line connecting the reference
point on the "elastic" wall to a second reference point on the fixed
lower channel wall.
Furthermore, the "wrapped" element, \c AlgebraicElement<FLUID_ELEMENT>, 
automatically computes the "shape derivatives" by finite-differencing the 
fluid residuals with respect to the degrees of freedom stored
in the \c GeomObject's geometric \c Data, just as in the
case of the \c MacroElement - based node-update procedure. 

If used in the form discussed in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
earlier example,</A> this methodology does not (yet!) improve
the sparsity of the node update: The geometric \c Data of the
"compound" \c GeomObject that represents the wall still contains
the positional degrees of freedom of \e all of the mesh's constituent
\c FSIHermiteBeamElements. This is wasteful because the position of a material
point on the (discretised) wall depends only on the positional degrees 
of freedom of the element that this point is located in. The 
(costly-to-compute) derivatives with respect to all other solid mechanics 
degree of freedom are zero. We will therefore modify the node-update 
procedure as follows: Each \c AlgebraicNode stores a pointer to the \c
FSIHermiteBeamElement that its reference point is located in.
This is possible because \c FSIHermiteBeamElements are derived from the
\c FiniteElement class which, in turn, is derived from
the \c GeomObject class. In other words, the sub-objects of
the compound \c MeshAsGeomObject are \c GeomObjects themselves. 
Their shape is parametrised by the \c FSIHermiteBeamElement's
local coordinate, \f$ s \f$ , which acts as the (sub-)\c GeomObject's 
intrinsic coordinate, \f$ \zeta \f$ . 


Given a pointer to a compound \c GeomObject, \c geom_obj_pt, say, and 
the intrinsic coordinate \f$ \zeta_{compound} \f$ = \c zeta_compound 
of a point in that \c GeomObject, the function \c GeomObject::locate_zeta(...)
may be used to determine a pointer, \c sub_obj_pt, to the sub-object
that this point is located in, and the vector of intrinsic 
coordinates \c zeta_sub_obj of the point in this sub-object.
This procedure is illustrated in this code fragment:

\code

 [...] 

// Vector containing the (1D) intrinsic coordinate in the 
// compound GeomObject:
Vector<double> zeta_compound(1);
zeta_compound[0]=0.3;

// Pointer to the sub-GeomObject:
GeomObject* sub_geom_obj_pt;

// Vector containing the (1D) intrinsic coordinate in the 
// sub-GeomObject:
Vector<double> zeta_sub_obj(1);

// Get it...
geom_obj_pt->locate_zeta(zeta_compound, sub_geom_obj_pt, zeta_sub_obj);


// Check the result: 

// Position vector to the point when viewed from
// the compound GeomObject
Vector<double> r_compound(2);
geom_obj_pt->position(zeta_compound,r_compound);

// Position vector to the point when viewed from
// the sub-GeomObject
Vector<double> r_sub(2);
sub_geom_obj_pt->position(zeta_sub,r_sub);

// With a bit of luck we should now have r_sub == r_compound...

[...]

\endcode

Here is an illustration of the relation between the various
coordinates and \c GeomObjects:


\image html locate_zeta.gif "Sketch of the various coordinates and GeomObjects. The (continuous) beam is parametrised by its Lagrangian coordinate xi which doubles as the intrinsic coordinate zeta for its role as a GeomObject. The (discretised) beam is a compound GeomObject, parametrised by the Lagrangian coordinate xi; its constituent FSI beam elements are sub-GeomObjects that are parametrised by their local coordinates, s. " 
\image latex locate_zeta.eps "Sketch of the various coordinates and GeomObjects. The (continuous) beam is parametrised by its Lagrangian coordinate xi which doubles as the intrinsic coordinate zeta for its role as a GeomObject. The (discretised) beam is a compound GeomObject, parametrised by the Lagrangian coordinate xi; its constituent FSI beam elements are sub-GeomObjects that are parametrised by their local coordinates, s. " width=0.75\textwidth

We note that the \c GeomObject base class provides a default
implementation for the \c GeomObject::locate_zeta(...) function
as a virtual member function which returns the \c GeomObject's 
<CODE>"this"</CODE> pointer and sets \f$ \zeta_{compound} =
\zeta_{sub} \f$ . Unless the function is overloaded in a
specific derived class, the \c GeomObject therefore acts as its own 
sub-object. This is a sensible default as it ensures that
(\c geom_obj_pt, \f$ \zeta_{compound}\f$ ) and (\c sub_geom_obj_pt, \c
\f$ \zeta_{sub}\f$ ) always identify the same point, regardless of 
whether nor not the \c GeomObject pointed to by \c geom_obj_pt is a 
"compound" \c GeomObject.

<HR>
<HR>

\section impl The implementation

The implementation of the sparse node-update strategy requires only
a few minor modifications to the \c
AlgebraicCollapsibleChannelMesh, first discussed in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
non-FSI example.</A>

<HR>


\subsection mesh The AlgebraicCollapsibleChannelMesh
We construct the mesh by multiple inheritance, combining the
already existing \c CollapsibleChannelMesh with the \c AlgebraicMesh
base class:

\dontinclude collapsible_channel_mesh.template.h
\skipline start_of_algebraic_collapsible_channel_mesh
\until {

The constructor calls the constructor of the underlying
\c CollapsibleChannelMesh and then calls the private member function
\c setup_algebraic_node_update() to initialise the data for the 
algebraic node update procedures.  (The initialisation is implemented in 
a separate function so it can be called from additional
mesh constructors that are not discussed here.) The destructor
remains empty.

\until virtual ~AlgebraicCollapsibleChannelMesh()

The function \c algebraic_node_update(...) is defined as a pure virtual
function in the \c AlgebraicMesh base class and therefore must be 
implemented, whereas the virtual function \c update_node_update(...)
is only required for refineable meshes and can remain
empty.

\skipline Update nodal position at time level t
\until };


The setup of the algebraic node update is very similar to that
used in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
non-FSI example discussed earlier.</A> The main difference between the
two versions of the mesh is that we use the function 
\c GeomObject::locate_zeta(...) to determine
the sub-\c GeomObject within which the reference point on the wall
is located. As discussed above, the default implementation of this
function in the \c GeomObject base class ensures that the mesh
can be used with compound and non-compound \c GeomObjects.

We start by determining the x and y-coordinates of the nodes and
decide if they are located in the collapsible part of
the mesh. (The positions of nodes that are located in the rigid upstream 
and downstream channel segments do not have to be updated; for such
nodes we skip the assignment of the node-update data. See the
discussion in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
non-FSI example</A> for details.) 

\dontinclude collapsible_channel_mesh.template.cc
\skipline start_setup
\until x <= (l_up + l_collapsible)

Assuming that the wall is in its undeformed position (we'll check
this in a second...), we determine the intrinsic coordinate 
of the reference point on the upper wall (taking the offset 
between \f$ x\f$ and \f$ \zeta\f$ into account: The left end of
the elastic wall is located at\f$ \zeta=0\f$ and at \f$ x = L_{up}
\f$ ), and identify the sub - \c GeomObject within which the reference 
point is located. 

\until locate_zeta

Just to be on the safe side, we double check that the wall
is still in its undeformed position:

\until #endif

Now we can create the node update data for the present \c
AlgebraicNode. The node update function involves a single \c
GeomObject: The (sub-)\c GeomObject within which the reference
point on the upper wall is located.

\until geom_object_pt[0] =

As in the mesh used in the 
<A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
non-FSI example</A> we store the x-coordinate of the reference
point on the lower wall, the fractional height of the node, and its intrinsic
coordinate in the (sub-)\c GeomObject on the upper wall.
We also store the intrinsic coordinate of the reference point 
in the compound \c GeomObject (i.e. the Lagrangian coordinate of the
reference point in the continuous beam). This will turn out to be useful in 
the refineable version of this mesh, to be discussed in 
<A HREF="../../fsi_collapsible_channel_adapt/html/index.html">the next
example.</A>

\until ref_value[3]

Finally, we create the node update information by passing the pointer 
to the mesh, the pointer to the \c GeomObject and the reference
values to the \c AlgebraicNode.

\until  end of setup_algebraic_node_update 


<HR>

\section driver The driver code

Since \c oomph-lib's various node update procedures use the same
interfaces, changing the node update strategy from the 
\c Domain/MacroElement - based procedure, discussed in the 
<A HREF="../../fsi_collapsible_channel/html/index.html">previous example,</A>
to the procedure implemented in the \c AlgebraicCollapsibleChannelMesh, 
only requires minimal changes to the driver code. In fact, the
changes are so trivial, that both versions are implemented in the
same driver code, 
<A HREF="../../../../demo_drivers/interaction/fsi_collapsible_channel/fsi_collapsible_channel.cc">fsi_collapsible_channel.cc</A>,
using compiler flags to switch from one version to the other. 
If the code is compiled with the flag \c -DMACRO_ELEMENT_NODE_UPDATE
the \c Domain/MacroElement - based node-update procedure, implemented in the
\c  MacroElementNodeUpdateCollapsibleChannelMesh is used, otherwise
the code uses the  \c AlgebraicCollapsibleChannelMesh, discussed
above. Here is one of the few portions of the code where 
the distinction between the two versions is required: The access
function to the "bulk" (fluid) mesh in the problem class.
 
\dontinclude fsi_collapsible_channel.cc
\skipline #ifdef MACRO_ELEMENT
\until #endif

Incidentally, the driver code also uses compiler flags to switch
between Crouzeix-Raviart and Taylor-Hood elements for the
discretisation of the Navier-Stokes equations. By default,
Crouzeix-Raviart elements are used; Taylor-Hood elements are used
if the code is compiled with with the flag \c -DTAYLOR_HOOD.



<HR>
<HR>

\section results Results

The animations shown below illustrate the interaction between fluid and 
solid mechanics degrees of freedom in the computations with the
algebraic node update. Comparison with the corresponding animations
for the \c Domain/MacroElement - based procedures, shown in the
<A HREF="../../fsi_collapsible_channel/html/index.html#comments">earlier
example</A> demonstrates the greatly improved sparsity of the
node update. With the algebraic node-update procedures, the residuals 
of the \c FSIHermiteBeamElements now only depend on the fluid degrees
of freedom in the adjacent fluid elements and on the solid mechanics 
degree of freedom in the  \c FSIHermiteBeamElements that affect the 
nodal position in these fluid elements.


\image html cr_alg.gif "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Crouzeix-Raviart elements.) " 
\image latex cr_alg.eps "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Crouzeix-Raviart elements.) " width=0.75\textwidth


Here is the corresponding animation for a discretisation with
2D Taylor-Hood elements. These elements have no internal \c Data but
the pressure degrees of freedom are stored at the fluid element's
corner nodes:

\image html th_alg.gif "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Taylor-Hood elements.) " 
\image latex th_alg.eps "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Taylor-Hood elements.) " width=0.75\textwidth


Finally, here is an animation that shows the (solid mechanics) degrees
of freedom that affect the node-update of a given fluid node. 
The red square marker shows the fluid node; the green numbers
show the number of the degrees of freedom at the \c SolidNodes
that are involved that fluid node's node update. With the
algebraic node update, the position of each fluid node is
only affected by the solid mechanics degree of freedom 
in the \c FSIHermiteBeamElement that contains its reference point.


\image html fsi_fluid_nodes.gif "Animation of the Data values that affect the node update of the fluid nodes. " 
\image latex fsi_fluid_nodes.eps "Animation of the Data values that affect the node update of the fluid nodes. " width=0.75\textwidth

The improved sparsity leads to a very significant speedup compared
to the \c MacroElement/Domain - based node update procedure. 

<HR>
<HR>

\section ex Exercises
-# Demonstrate that the dramatically improved execution speed
   achieved with the \c AlgebraicCollapsibleChannelMesh is mainly due to the
   improved sparsity of the node update, achieved by using the 
   \c GeomObject::locate_zeta(...) function. \n\n \b Hint: You can either
   copy the basic \c MyAlgebraicCollapsibleChannelMesh in the file 
   <A HREF="../../../../demo_drivers/navier_stokes/collapsible_channel/my_alg_channel_mesh.h">
   my_algebraic_collapsible_channel_mesh.h</A>,
   developed for the 
   <A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
   non-FSI version of the collapsible channel problem,</A> into 
   the FSI driver code
   <A HREF="../../../../demo_drivers/interaction/fsi_collapsible_channel/fsi_collapsible_channel.cc">fsi_collapsible_channel.cc</A> and use that mesh
   instead of the \c  AlgebraicCollapsibleChannelMesh, or
   replace the line
   \n\n
   \code
   this->Wall_pt->locate_zeta(zeta,geom_obj_pt,s);
   \endcode
   \n\n
   in the function \c 
   AlgebraicCollapsibleChannelMesh<ELEMENT>::setup_algebraic_node_update()
   in <A HREF="../../../../src/meshes/collapsible_channel_mesh.template.cc">
   collapsible_channel_mesh.template.cc</A> by \n\n
   \code
   this->Wall_pt->GeomObject::locate_zeta(zeta,geom_obj_pt,s);
   \endcode
   \n\n
   thus bypassing the "sparsification".
   \n\n
-# Explore how the speedup achievable with the algebraic node
   update procedure depends on the mesh resolution. A speedup by a
   factor of ten is typical for computations on the coarse mesh
   used for the validation runs; much more dramatic speedups 
   tend to be obtained on finer meshes. \n\n        
.


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/fsi_collapsible_channel/">
demo_drivers/interaction/fsi_collapsible_channel/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/fsi_collapsible_channel/fsi_collapsible_channel.cc">
demo_drivers/interaction/fsi_collapsible_channel/fsi_collapsible_channel.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage A fluid-structure-interaction problem: Flow in a 2D collapsible channel


 
In this document we discuss the solution of a standard 
fluid-structure interaction problem -- finite-Reynolds-number flow 
in a 2D collapsible channel. We shall demonstrate that the driver code
for this multi-physics problem is a straightforward combination of the driver
codes for the two corresponding single-physics problems
discussed earlier:
- <B><A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
  Flow in a 2D channel with prescribed wall motion:</A></B>\n\n
  In this single-physics problem, we represented the moving wall by
  a \c GeomObject, and created a \c Domain object to provide     
  an analytical representation of the domain boundary.
  We discretised the Navier-Stokes equations with 2D
  Crouzeix-Raviart elements and updated their nodal positions
  (in response to the prescribed changes in the wall position) by their
  \c MacroElement representation.
.
and
- <B><A HREF="../../../beam/tensioned_string/html/index.html">
  The deformation of a pressure-loaded elastic beam:</A></B> \n\n
  In this single-physics problem, we discretised the elastic
  beam with \c HermiteBeamElements and computed its deformation
  (determined by the positional degrees of freedom, stored
  at the \c HermiteBeamElement's \c SolidNodes) in response to the prescribed
  pressure load.
.


<HR>
<HR> 


\section problem The problem 




<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>Flow in a 2D collapsible channel </B>
</CENTER> 

The figure below shows a sketch of the problem: Flow is driven by a 
prescribed pressure drop through a 2D channel of width 
\f$ H^* \f$ and total length
\f$ L^*_{total} = L^*_{up} + L^*_{collapsible} + L^*_{down}. \f$
The upstream and downstream lengths of the channel are rigid, whereas
the upper wall in the central section is an elastic membrane whose
shape is parametrised by a Lagrangian coordinate, \f$ \xi^* \f$ , so 
that the position vector to the moving wall is given by 
\f$ {\bf R}_w^*(\xi^*,t^*) \f$ . The wall is loaded by the
external pressure \f$ p_{ext}^* \f$ and by the traction 
that the viscous fluid exerts on it. The components of the load 
vector \f$ {\bf f}^* \f$ that acts on the wall are 
therefore given by 
\f[  
f^*_{i} = (p^* - p^*_{ext}) N_i -
\mu \left(
\frac{\partial u^*_i}{\partial x^*_j} +
\frac{\partial u^*_j}{\partial x^*_i} 
\right)
N_j \ \ \ \ \ \ \mbox{for $i=1,2$,}
\f]
where \f$ N_i \ (\mbox{\ for \ } i=1,2) \f$ are 
the components of the outer unit normal on the fluid domain. 
 

\image html collapsible_channel_sketch.gif "Sketch of the problem. " 
\image latex collapsible_channel_sketch.eps "Sketch of the problem. " width=0.75\textwidth

We scale all lengths on the channel width, \f$ H^* \f$ , use the
average velocity through the undeformed 
channel, \f$ U =P^*_{up} H^{*2}/(12 \mu L^*_{total}) \f$ , to scale 
the velocities, and use \f$ H^{*}/U \f$ to non-dimensionalise
time. Finally, the fluid pressure is non-dimensionalised on the viscous scale 
\f$ p^{*} = p \mu U/H^{*}\f$.
(As usual, asterisks distinguish dimensional 
parameters from their non-dimensional equivalents.) 

With this non-dimensionalisation, the Navier-Stokes equations have
the same form as in the <A HREF="../../../navier_stokes/collapsible_channel/html/index.html">earlier
example with prescribed wall motion:</A>
\f[
Re\left(St\frac{\partial u_i}{\partial t} +  
u_j\frac{\partial u_i}{\partial x_j}\right) =
- \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right)
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
 \ \ \ \ \ \ \ \ \ \ (2)
\f]
with \f$ St=1 \f$. As 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
before</A>, the flow is subject to the following boundary and initial
conditions:
- Initial condition: Poiseuille flow, i.e. 
  \f[
  {\bf u}(x_1,x_2,t=0) = {\bf u}_{Poiseuille}(x_1,x_2) =   
  6 \ x_2 \  (1-x_2) \ {\bf e}_1.
  \ \ \ \ \ \ \ \ \ \ (3)
  \f]

- Parallel inflow, \f$ {\bf u} \cdot {\bf e}_2 = {\bf 0}, \f$
  and an applied axial traction of \f$ {\bf t} \cdot {\bf e}_1 = 
  p_{up} = 12 \ L_{total} \f$ at the upstream end, \f$ x_1=0\f$.

- Parallel, axially traction-free outflow at the downstream end, i.e.
  \f$ {\bf u} \cdot {\bf e}_2 = {\bf 0} \f$ and  
  \f$ {\bf t} \cdot {\bf e}_1 = p_{down} = 0  \f$ at  \f$
  x_1=L_{total}. \f$ 

- No slip on all channel walls, i.e.  \f$ {\bf u} = {\bf 0} \f$ on 
  the rigid walls and 
  \f[
  {\bf u} = \frac{\partial {\bf R}_w}{\partial t}
  \mbox{\ \ \ \ on the moving wall,}
  \ \ \ \ \ \ \ \ \ \ (4)
  \f]
These boundary conditions are identical to those in the problem with
prescribed wall motion, apart from the fact that in the present
problem the wall motion, described by \f$ {\bf R}_w(\xi, t), \f$
has to be determined as part of the solution.

We model the elastic membrane as a thin-walled elastic Kirchhoff-Love
beam of wall thickness \f$ h^* \f$ , subject to an axial (2nd Piola-Kirchhoff)
pre-stress \f$ \sigma^*_0. \f$ The beam' effective (1D) elastic modulus
is given by \f$ E_{eff} = E/(1-\nu^2), \f$ where \f$ E \f$ and \f$ \nu \f$
are its 3D Young's modulus and Poisson's ratio, respectively.
The beam's deformation is governed by the principle of virtual
displacements, discussed in detail 
<A HREF="../../../beam/tensioned_string/html/index.html">elsewhere.</A> 
As in the Navier-Stokes equations, we scale all lengths in
the beam problem on the channel's width, \f$ H^*. \f$ The
non-dimensional position vector
\f$ {\bf r}_w(\xi) \f$ to the undeformed wall is then given by
\f[
{\bf r}_w(\xi) = 
\left(
\begin{array}{c}
L_{up} + \xi  \\
1 
\end{array}
\right) \ \ \ \ \ \ \mbox{where $\xi \in [0,L_{collapsible}].$}
 \ \ \ \ \ \ \ \ \ \ (5)
\f]
Our non-dimensionalisation of the principle of virtual displacements
requires all stresses and tractions to be non-dimensionalised 
on the beam's (effective 1D) elastic modulus, \f$ E_{eff} \f$ ,
therefore we define the non-dimensional external pressure
as \f$ p_{ext} = p_{ext}^{*}/E_{eff}. \f$ The non-dimensional 
load vector \f$ {\bf f} = {\bf f}^{*}/E_{eff} \f$ 
that acts on the wall is then given by
\f[
f_i = - p_{ext} N_i  + Q \left( p N_i -
\left(
\frac{\partial u_i}{\partial x_j} +
\frac{\partial u_j}{\partial x_i} 
\right) N_j \right) 
\ \ \ \ \ \ \mbox{for $i=1,2,$}
\f]
where 
\f[
Q = \frac{\mu U^*}{E_{eff} H^*_{tot}}
\f]
is the ratio of the fluid pressure scale, \f$ \mu U/H^* \f$ ,
used to non-dimensionalise the Navier-Stokes equations, to the
beam's effective elastic modulus, \f$ E_{eff} \f$. The parameter \f$ Q
\f$ therefore indicates the strength of the fluid-structure
interaction. In particular, if \f$ Q=0 \f$ the wall deformation is 
not affected by the fluid flow. 
</TD>
</TR>
</TABLE>  
</CENTER>
 
<HR>
<HR>



 
\section reslt Results
The figure below shows a snapshot of the flow field, taken from 
<A HREF="../figures/fsi_taylor_hood_flow8.avi">the animation of the
computational results</A>. The first four figures show (from top left to bottom
right) "carpet plots" of the axial and transverse velocities,
the axial component of the perturbation velocity 
\f$ {\bf u} - {\bf u}_{Poiseuille} \f$ , and the pressure
distribution. The 2D contour plot at the bottom
of the figure shows a contour plot of the pressure and a few
instantaneous streamlines.

The overall structure of the flow field is very similar to that observed
in the corresponding  
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
problem with prescribed wall motion:</A> The wall oscillation generates
a large-amplitude sloshing flow that is superimposed on the pressure-driven
Poiseuille flow. At the instant shown in this figure, the wall is moving
inwards. Consequently, the sloshing flow generated in the upstream 
(downstream)  rigid sections is directed against (in the same)  
direction as that of the pressure-driven mean flow.

\image html fsi_taylor_hood_flow.gif "Snapshot from the animation of the flow field. " 
\image latex fsi_taylor_hood_flow.eps "Snapshot from the animation of the flow field. " width=0.8\textwidth


In the present problem, the wall motion is, of course, not prescribed
but determined as part of the overall solution. 
<A HREF="http://www.maths.man.ac.uk/~mheil/MATTHIAS/PDF/JensenHeil2003.pdf">
Jensen \& Heil's (2003)</A> asymptotic analysis of the problem shows that
the period of the oscillations is determined by the balance between fluid
inertia and the elastic restoring forces.
In certain parameter regimes (at sufficiently large Reynolds number), 
the wall can extract energy from the pressure-driven mean flow, 
causing the oscillations to grow in amplitude. In the present
example, the Reynolds number is too small for this to happen and
viscous dissipation causes the oscillations to decay, as shown in this
plot.


\image html trace_fsi_taylor_hood.gif "Time-trace of the axial velocities at two control points in the upstream and downstream cross-sections, and the vertical position of a control point on the wall. " 
\image latex trace_fsi_taylor_hood.eps "Time-trace of the axial velocities at two control points in the upstream and downstream cross-sections, and the vertical position of a control point on the wall. " width=0.20\textwidth


<HR>
<HR>

\section overview Overview: How to solve fluid-structure interaction problems with oomph-lib
Before attempting to solve a fluid-structure interaction problem,
we generally recommend to first study the constituent
single-physics problems in isolation, i.e.
- Discretise the solid mechanics problem with a suitable 
  \c SolidFiniteElement (e.g.  a \c HermiteBeamElement) and determine the
  wall deformation in response to a prescribed external load. If possible,
  choose a load that is vaguely "representative" of the loads expected
  in the actual fluid-structure interaction problem. 
- Discretise the fluid mechanics problem in a domain in which 
  the deformation of the "elastic boundary" is described by 
  a \c GeomObject that performs a prescribed motion. 
  Choose a suitable node-update strategy to adjust the position
  of the fluid nodes in response to the (prescribed) motion of the 
  "elastic boundary". (As demonstrated in many previous examples,
  a \c Domain/MacroElement - based node-update method is
  very easy to implement.) If possible, try to impose a wall 
  motion that is "representative" of the type of wall motion expected in 
  the actual fluid-structure interaction problem.
.


Once the behaviour of the two isolated single-physics problems is 
sufficiently well understood (e.g. what spatial discretisations are
required, etc.), combine the two single-physics problems to a
fully-coupled problem. The coupling introduces two types of interactions 
that must be incorporated into the computational framework:
-# The position of the nodes in the fluid mesh depends on the wall
   shape. Since the wall shape is now determined by the degrees of
   freedom in the \c SolidFiniteElements that we use to discretise the
   wall, we must ensure that: 
   -# The node-update operations that we developed for problems
      with prescribed wall motion (in which the position of the
      curvilinear domain boundaries is determined by \c GeomObjects)
      also work for problems in which the curvilinear domain
      boundaries are represented by \c SolidFiniteElements.
   -# The dependence of the residuals of the fluid elements
      on the solid-mechanics degrees of freedom that affect the
      positions of the fluid nodes (via the node-update 
      operation) are taken into account when computing 
      the fluid elements' Jacobian matrices. (See the discussion of
      the "shape derivatives" in the 
      <A HREF="../../macro_element_free_boundary_poisson/html/index.html"> 
      "toy" free-boundary Poisson problem</A> for details.)  
-# The \c SolidFiniteElements are not only loaded by the external load
   but also by the traction that the fluid exerts on
   them. When computing the residuals of the
   \c SolidFiniteElements we have to evaluate the combined load 
   vector at the \c SolidFiniteElement's Gauss points. We therefore 
   have to provide a lookup scheme that specifies: 
   -# which fluid element is adjacent to 
      a given Gauss point in the \c SolidFiniteElements, and
   -# which (fluid) degrees of freedom affect the fluid traction at that 
      point. 
   .
   The dependence of the residuals of the \c SolidFiniteElements on these
   (fluid) degrees of freedom must be taken into account when 
   computing the \c SolidFiniteElement's Jacobian matrix.
. 
   


\c oomph-lib provides a number of high-level functions that
allow the required lookup schemes to be generated completely
automatically. In the following sections we shall provide a brief 
discussion of the methodology but we stress that the details
are not particularly important for the "user". If you just
want to "use" \c oomph-lib's fluid-structure interaction 
capabilities and don't care too much about the technical
details, you may wish to skip the next few sections and continue with the \ref
driver_code, where we demonstrate that, apart from a few trivial 
modifications,  the driver code for the
fully-coupled fluid-structure interaction problem is
a straightforward combination of the two single-physics driver codes.


<HR>
<HR>

\section how_its_done Brief discussion of the implementation


\subsection mesh MacroElement-based (fluid-)node updates in FSI problems 
\subsubsection shape_deriv The shape derivatives

When discussing our 
<A HREF="../../macro_element_free_boundary_poisson/html/index.html"> "toy"
free-boundary Poisson problem</A> we demonstrated how
a "bulk" element's \c MacroElement - representation allows the efficient 
automatic evaluation of the "shape derivatives" -- the derivatives of the 
"bulk" equations (here the Navier-Stokes equations) with 
respect to the \c Data values (here the nodal positions in the beam elements) 
that determine the shape of the domain boundary. The methodology 
discussed in  the context of the "toy" problem may also be used 
for genuine fluid-structure interaction problems, such as the problem 
considered here, provided
-# A \c MacroElement/Domain - based Mesh is used to discretise the
   fluid domain.
-# The moving boundary is represented by a \c GeomObject whose 
   member function \c GeomObject::geom_data_pt(...) provides
   access to its "geometric" \c Data, i.e. the \c Data that affects 
   its shape.
.
Condition 1 is satisfied as the <A HREF="../../../meshes/mesh_list/html/index.html#collapsible_channel">
<CODE>CollapsibleChannelMesh</CODE></A> used for the simulation of the
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
single-physics fluids problem with prescribed wall motion</A> employs the \c
CollapsibleChannelDomain to perform the node update in response to 
changes in the domain boundary. Multiple inheritance may therefore
be used to upgrade the existing <A HREF="../../../meshes/mesh_list/html/index.html#collapsible_channel">
<CODE>CollapsibleChannelMesh</CODE></A> to a mesh that is derived
from the \c MacroElementNodeUpdateMesh base class. The dependence
of the residuals of the "bulk" (fluid) elements on the geometric \c
Data that affects its nodal positions is automatically taken into 
account if the "bulk" (fluid) elements (of type \c BULK_ELEMENT, say)
are "upgraded" to the "wrapped" class 
\c MacroElementNodeUpdateElement<BULK_ELEMENT>. We refer to the discussion
of the <A HREF="../../macro_element_free_boundary_poisson/html/index.html"> 
"toy" free-boundary Poisson problem</A> for details of this methodology. 

\subsubsection wall_geom_object Representing the wall mesh as a GeomObject

What remains to be done is to represent the wall (discretised
by the \c SolidFiniteElements contained in the wall mesh) as a \c GeomObject.
For this purpose, \c oomph-lib provides the class 

\code 
class MeshAsGeomObject : public GeomObject
\endcode 

whose constructor takes a pointer to the \c Mesh that is to be
represented as a \c GeomObject. The Lagrangian and Eulerian dimensions 
of the mesh are taken to be the dimension of elements and the dimension
of the nodes, respectively. The conversion from \c Mesh into \c
GeomObject only makes sense if the element is derived from a 
\c FiniteElement, itself a \c GeomObject, which has
associated Lagrangian (local) and Eulerian coordinate systems.
In our collapsible channel problem, the Lagrangian dimension is one 
because the wall 
is parametrised by a single intrinsic coordinate (the shape of the
mesh's constituent elements is parametrised
by a single local coordinate), and the Eulerian dimension is two 
because the nodes
in the wall mesh have two Eulerian coordinates. Assuming
that \c Wall_mesh_pt stores a pointer to the wall mesh, the
code

\dontinclude fsi_collapsible_channel.cc
\skipline Build a geometric object
\until new MeshAsGeomObject
 
creates a \c GeomObject representation of the wall mesh that
allows us to obtain the position of a material point on the wall
(parametrised by its Lagrangian coordinate \f$\xi \f$) from
the \c GeomObject::position(...) function. For instance,
the following code computes the position vector \f$ {\bf r} \f$ to the
material point on the deformed wall, located at  \f$ \xi=0.5. \f$
\code
Vector<double> xi(1);
xi[0]=0.5;
Vector<double> r(2);
wall_geom_object->position(xi,r);
\endcode

Here is a graphical illustration of the various representations of
the domain boundary:

<B>1. A "normal" GeomObject</B>

In single-physics problems with prescribed boundary motion
the domain boundary may be represented by a \c GeomObject. 
The \c GeomObject provides a parametrisation of its shape
in terms of an intrinsic coordinate, \f$ \zeta \f$ , as shown in this
sketch:

\image html wall_mesh_geometry_sketch1.gif "A geometric object parametrised by an intrinsic coordinate. " 
\image latex wall_mesh_geometry_sketch1.eps "A geometric object parametrised by an intrinsic coordinate. " width=0.75\textwidth

The function \c GeomObject::position(zeta,r)
computes the position vector \c r to a point on the \c GeomObject,
as identified by its intrinsic coordinate \c zeta.

<B>2. A beam/shell structure</B>

In beam/shell problems, the shape of the deformed structure
is parametrised by its Lagrangian coordinate , \f$\xi. \f$


\image html wall_mesh_geometry_sketch2.gif "A (continuous) beam structure, parametrised by a Lagrangian coordinate. " 
\image latex wall_mesh_geometry_sketch2.eps "A (continuous) beam structure, parametrised by a Lagrangian coordinate. " width=0.75\textwidth

Beam/shell structures may therefore act as \c GeomObjects if
we interpret their Lagrangian coordinate, \f$ \xi \f$ , as the
\c GeomObject's intrinsic coordinate, \f$ \zeta \f$. 
 
 
<B>3. A discretised beam/shell structure</B>

In an \c oomph-lib computation, the beam/shell structure will,
of course, have been discretised by a number of \c SolidFiniteElements. The \c
MeshAsGeomObject discussed above, is therefore a "compound" 
\c GeomObject that contains a number of sub-objects -- the 
mesh's constituent \c SolidFiniteElements. Within the "compound"
\c GeomObject each sub-object acts as a \c GeomObject in its own 
right -- the shape of a \c SolidFiniteElement is parametrised 
by its local coordinate \c s. 

\image html wall_mesh_geometry_sketch3.gif "A discretised beam structure. " 
\image latex wall_mesh_geometry_sketch3.eps "A discretised beam structure. " width=0.75\textwidth

The \c MeshAsGeomObject::position(...) function therefore determines
the position vector to the point labelled by the Lagrangian coordinate
\f$ \xi \f$ (i.e. the \c GeomObject's intrinsic coordinate, \f$ \zeta \f$) in
a two-stage process: First it determines which of its constituent
\c SolidFiniteElements "contains" the relevant material point (This is possible
because the function \c SolidFiniteElement::interpolated_xi(...)
provides access to the Lagrangian coordinate inside the element)
and then uses the \c SolidFiniteElement::interpolated_x(...)
function to determine the Eulerian coordinates of that point.

<HR>

\subsection fsi_wall_elements Applying the fluid-traction to the wall elements: FSIWallElements

Next, we shall discuss how the fluid traction is added to the load terms
in the wall equations. As mentioned above, the computation of the
residuals of the wall equations requires the evaluation of the
combined load vector at the Gauss points in the wall elements. Furthermore,
the dependence of the residuals on those (fluid) degrees of freedom
that affect the traction must be taken into account when computing
the wall element's Jacobian matrix. Storage for the various lookup schemes 
required for such computations is provided in the virtual base class
\c FSIWallElement whose inheritance structure is as follows:

\code
class FSIWallElement : public virtual SolidFiniteElement, 
                       public virtual ElementWithExternalElement
\endcode

The \c FiniteElement class inherits from \c GeomObject and by
default the \c GeomObject::position(...) function calls the function \c
FiniteElement::interpolate_x(..); in other words,
the element's local coordinate is
regarded as the intrinsic
coordinate that parametrises its shape. Thus we already have a standard 
interface through which an \c FSIWallElement can be used to parametrise 
the shape of (part of) the domain boundary. 

By inheriting from the \c SolidFiniteElement class,
we establish that the shape of the \c FSIWallElement is determined 
by the positional \c Data stored at its constituent \c SolidNodes. Recall that 
this information is required during the computation of the shape 
derivatives of the fluid equations.

The \c ElementWithExternalElement class provides the generic storage and helper
functions required to keep track of the external elements that are
adjacent to Gauss points in any \c FiniteElement, see also the
tutorial on 
<a href="../../../multi_physics/multi_domain_ref_b_convect/html/index.html">
Boussinesq convection with a multi-domain approach. </a>
In the \c
FSIWallElement, it is the fluid elements that load the structure which
are adjacent to the Gauss points.
The various lookup schemes required to determine these fluid elements
may be generated completely automatically by the helper function

\code
template<class FLUID_ELEMENT, unsigned DIM_FLUID>
void FSI_functions::setup_fluid_load_info_for_solid_elements(
     Problem* problem_pt,
     const unsigned &boundary_in_fluid_mesh,
     Mesh* const& fluid_mesh_pt,
     Mesh* const& solid_mesh_pt);
\endcode

which is defined in the namespace \c FSI_functions. The template
parameters \c FLUID_ELEMENT and \c DIM_FLUID specify the type of 
the fluid element and its spatial (Eulerian) dimension, respectively. 
The arguments are the pointer to the problem, the number of the
boundary in the fluid mesh adjacent to the elastic wall, and the 
pointers to the fluid and wall meshes. The function assumes that 
boundary coordinates have been set for the fluid nodes on the mesh 
boundary specified by the argument \c boundary_in_fluid_mesh, and 
that these boundary coordinates are consistent with the
parametrisation of the wall mesh by its Lagrangian coordinate. (See 
<A HREF="../../../poisson/fish_poisson2/html/index.html#boundary_coords">the
discussion of the mesh generation procedures for domains with
curvilinear boundaries</A> for a more detailed discussion 
of boundary coordinates for nodes.)

The \c FSIWallElement provides the protected member function
\c FSIWallElement::fluid_load_vector(...) which may be used
in a specific \c FSIWallElement (such as the \c FSIHermiteBeamElement)
to compute the fluid traction (on the solid mechanics stress scale), 
and to add it to any external load that may already be acting on 
the element. The conversion from the fluid to the solid 
non-dimensionalisation of the traction is performed automatically 
by multiplying the traction vector (on the fluid stress scale) 
obtained from the "adjacent" \c FSIFluidElement by the stress 
ratio \f$ Q \f$ which has a default value of \f$ Q=1. \f$ This 
default assignment may be overwritten by setting a pointer to a
variable that specifies \f$ Q \f$ by using the function 
\c FSIWallElement::q_pt().


The class also overloads the 
\c SolidFiniteElement::fill_in_contribution_to_jacobian(...) function
so that the derivatives of the residuals with 
respect to those unknowns that affect the fluid traction on the wall
are included when computing the element's Jacobian matrix. 

In the driver code, discussed below, we will discretise the wall
with \c FSHermiteBeamElements, a class that is composed
(by multiple inheritance) from the single-physics 
\c HermiteBeamElement and the \c FSIWallElement base class:

\code
class FSIHermiteBeamElement : public virtual HermiteBeamElement, 
                              public virtual FSIWallElement
\endcode

\subsubsection fsi_fluid_elements Obtaining the fluid traction from "adjacent" fluid elements: FSIFluidElements

Having provided a function that allows us to determine which 
fluid elements are located next to a given \c FSIWallElement, we have 
to define standard interfaces through which we can obtain the traction 
that the "adjacent" fluid element exerts onto the wall. Furthermore,
we have to determine which unknowns affect the fluid traction to
enable us to evaluate the derivatives of the wall residuals with
respect to these unknowns. Interfaces for the relevant functions 
are provided in the base class \c FSIFluidElement, whose most 
important member function (for the purpose of the present discussion) is 
\c FSIFluidElement::get_load(...). The purpose of this function is
to compute the traction exerted by the \c FSIFluidElement 
onto the adjacent \c FSIWallElement, given the outer unit
vector onto the \c FSIFluidElement. The \c FSIFluidElement class
has two further pure virtual member functions whose role is to 
determine the unknowns (e.g. velocity and pressure values) that 
affect the fluid traction. For newly developed fluid elements
these functions have to implemented on a case-by-case basis.
However, all existing fluid elements in \c oomph-lib are already derived
from the \c FSIFluidElement class and therefore provide
a suitable implementation of these functions. It is therefore
not necessary to explicitly "upgrade" \c oomph-lib's fluid elements 
before using them in FSI computations.




<HR>
<HR>


\section driver_code Overview of the driver code

The driver code for the fully-coupled
fluid-structure interaction problem is a straightforward
combination of the two single-physics codes with a few
trivial additional steps. The main steps in the problem setup are:
-# Create the wall mesh, using \c FSIHermiteBeamElements instead of
   \c HermiteBeamElements. 
-# <B>NEW:</B> Create a \c GeomObject representation of the wall mesh, using
   the \c MeshAsGeomObject class. 
-# Create the fluid mesh, using elements of type
   \c MacroElementNodeUpdateElement<QCrouzeixRaviartElement<2>> instead of
   \c QCrouzeixRaviartElement<2>. Use the \c MeshAsGeomObject
   representation of the elastic wall, created in the previous step
   to represent the curvilinear domain boundary.
-# Build a mesh of traction elements that apply the
   prescribed-traction boundary condition at the inflow. 
   Add all three meshes (fluid, solid and traction mesh)
   to the \c Problem's collection of sub-meshes and build the global
   mesh.
-# Pass the relevant function pointers (for Reynolds number, 
   external loads, etc) to the various elements and apply the boundary 
   conditions. 
-# <B>NEW:</B> Call the function 
   \c FSI_functions::setup_fluid_load_info_for_solid_elements(...) to
   set up the lookup scheme that establishes which fluid
   elements affect the traction on the wall. 
-# Set up the equation numbering scheme. 
-# Done! . 

Most steps in this sequence are either identical to those 
in the corresponding single-physics codes, or require only trivial
modifications. Each of the two new steps (2 and 6) can be implemented
with a single line of code. Consequently, most of the driver
code, discussed in detail below, contains verbatim copies of code
segments from the respective single-physics codes. 

<HR>
<HR>



\section variables  Namespace for the global physical variables
The namespace for the "global" physical parameters contains
the same fluid parameters as in the <A HREF="../../../navier_stokes/collapsible_channel/html/index.html">collapsible
channel problem with prescribed wall motion</A>: We define the Reynolds and
Womersley numbers and the fluid pressure at the upstream end, and
provide a function that specifies the applied traction at the inflow:

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_physical_parameters
\until end traction

Next we define the wall parameters (wall thickness, prestress and
external pressure) and assign default values. Note that the function
that specifies the load on the wall contains only the load due to 
the external pressure -- the additional load due to the fluid traction
will be added automatically by the \c FSIWallElements, using the
lookup scheme set up by the function 
\c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
discussed earlier.

\until end of load

Finally, we define the interaction parameter \f$ Q \f$ and give it 
a default value. 

\until end of namespace 


<HR>
<HR> 



\section undeformed The undeformed wall
We represent the undeformed geometry of the elastic wall, defined by 
equation (5), as a \c GeomObject, specifying
the \f$ x \f$-coordinate of its left end and its (constant) 
\f$ y \f$-coordinate as arguments to the constructor:

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_underformed_wall
\until  }

The two versions of the \c position(...) function are straightforward:

\until end of position

Since the \c GeomObject is used to specify the undeformed shape of a
\c HermiteBeamElement, the function \c
GeomObject::d2position(...) must be implemented to define
the beam's curvature in the undeformed configuration. 

\until end of d2position

The private member data contains the two geometric parameters.

\until end_of_undeformed_wall

<HR>
<HR>

\section main The driver code

 As with most previous time-dependent codes, we use command line
arguments to indicate if the code is run during \c oomph-lib's
self-test procedures. If command line arguments are specified, we use a
coarser discretisation and perform fewer timesteps.
After storing the command line arguments, we choose the number 
of elements in the mesh, and set the lengths of the domain.

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_main
\until  double ly=1.0;

We assign values for the external pressure (on the
wall stiffness scale) and for the upstream fluid pressure (on the
fluid pressure scale). The latter is again chosen so that
in the absence of any wall deformation, the applied pressure 
difference would drive steady Poiseuille flow through the channel.

\until P_up

We build the problem with 2D quadrilateral Crouzeix-Raviart elements, 
"upgraded" to \c MacroElementNodeUpdateElements. This ensures that the "shape
derivatives" (the derivatives of the fluid residuals with respect
to the solid mechanics degrees of freedom that affect the nodal 
positions in the fluid elements), are incorporated into the
fluid elements' Jacobian matrices.

\skipline Build the problem with QCrouzeixRaviartElements
\until ly);

We choose the timestepping parameters before assigning the
initial conditions. (Preliminary computations showed that
the system performs oscillations with approximately  
unit period, so the chosen value for the timestep \c dt
corresponds to a time-integration with about 40 timesteps per
period.)


\skipline Timestep.
\until  set_initial_condition()


Next we specify the output directory, open a trace file 
and document the initial conditions

\until doc_info.number()++;

The timestepping loop is identical to that in 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
the problem with prescribed wall motion:</A>


\until end of main

<HR>
<HR> 

\section problemclass The problem class

The problem class is very similar to that used for the 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">problem
with prescribed wall motion</A>. We specify the type of the fluid
element as a template parameter and pass the number of
elements and the lengths of the domain to the constructor:

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_problem_class
\until  }

We provide access functions to the (pointers to) the fluid mesh,
\skipline Access function
\until }

and the wall mesh:

\skipline Access function for the wall mesh
\until end of access to wall mesh

 
Unlike the <A HREF="../../../navier_stokes/collapsible_channel/html/index.html">problem
with prescribed wall motion,</A> the FSI problem does not
have any time-dependent boundary conditions, therefore the
pure virtual functions \c Problem::action_before_solve() and 
\c Problem::action_after_solve() can remain empty, and the
function \c Problem::actions_before_implicit_timestep() is not needed.

\until actions_after_newton_solve

However, since the wall displacement (which is determined as part of 
the solution!) affects the nodal positions in the fluid mesh
via the \c MacroElement/\c Domain - based node-update,
the position of the fluid nodes must be updated whenever
the Newton solver updates the unknowns. This is precisely
what the function \c Problem::actions_before_newton_convergence_check()
is for; see <A HREF="../../../order_of_action_functions/html/index.html">
the discussion of \c oomph-lib's various "action" functions</A> 
for more details.


\until }

The functions \c doc_solution(...) and \c set_initial_condition()
do what they always do. 

\until set_initial_condition()

The private member function \c create_traction_elements(...) 
is used to attach the applied traction elements to the upstream end of
the channel, exactly as in the <A HREF="../../../navier_stokes/collapsible_channel/html/index.html">problem
with prescribed wall motion.</A>

\until &traction_mesh_pt

The private member data stores the problem parameters,

\until  double Ly;

the pointers to the fluid mesh,

\skipline Pointer to the "bulk"
\until MacroElementNode

the surface mesh that contains the applied-traction elements,

\skipline Pointer to the "surface
\until Applied_fluid

and the wall mesh,

\until OneDLagrangian

as well as pointers to various control nodes


\until end of problem class

<HR>
<HR>

\section problemcontr The problem constructor

We copy the various mesh parameters to the Problem's private data

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_constructor
\until  Ly=ly;

and increase the maximum value of the residual that is permitted 
during the Newton iteration to accommodate possible poor initial
guesses for the solution.
\until Max_residual


We construct a \c BDF<2> - timestepper for the time-integration 
of the fluid equations and add it to the Problem's collection of
timesteppers:

\until  add_time_stepper_pt(new BDF<2>);

The wall mesh is built as in the corresponding 
<A HREF="../../../beam/tensioned_string/html/index.html">
single-physics beam problem:</A> We create the \c GeomObject that 
describes the undeformed wall shape and construct the wall mesh,
this time with \c FSIHermiteBeamElements:

\until (Ncollapsible

We note that even though the same number of fluid and beam elements
are used to discretise the common boundary, the discretisation of
the two domains is \b not matching as the shape of the domain boundary is 
represented by piecewise cubic Hermite polynomials within the beam
elements, and by piecewise quadratic Lagrange polynomials within the
fluid elements. Mathematically, this does not cause any problems as 
both representations converge to the same boundary shape
as the meshes are refined further and further. We stress that \c oomph-lib
does not even require the number of elements along the interface
to match; see \ref ex .


The \c MacroElement/Domain - based
fluid-mesh update requires the wall shape (which is now parametrised
by the wall mesh's constituent elements) to be represented
by single \c GeomObject. We create the required "compound" \c
GeomObject using the \c MeshAsGeomObject class:

\until new MeshAsGeomObject
 
This "compound" \c GeomObject can now be used to build the 
fluid mesh, exactly as in  
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
the corresponding single-physics fluids problem.</A> (The
\c MacroElementNodeUpdateCollapsibleChannelMesh is a trivial
extension of the \c CollapsibleChannelMesh from which it is 
derived; see the discussion of the 
<A HREF="../../macro_element_free_boundary_poisson/html/index.html#mesh"> 
"toy" free-boundary Poisson problem</A> for details.)

\skipline Build bulk (fluid)
\until time_stepper_pt()

As in <A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
the fluids problem with prescribed wall motion </A> we use a
"boundary-layer squash function" to distribute the fluid elements
non-uniformly across the channel so that more elements are located
inside the thin Stokes layers that are likely to develop near the wall.

\until node_update();

We create the sub-mesh that stores the applied traction elements
and attach the elements to the upstream end of the channel. The
three sub-meshes are then combined into the \c Problem's single 
global mesh.

\skipline Create "surface mesh"
\until build_global_mesh()

We complete the build process for the fluid elements
by passing the pointers to the relevant problem parameters to the
elements,

\until   } // end loop over elements

then we apply the boundary conditions for the fluid velocities:
- both axial and transverse velocities are pinned along the bottom
  and the top boundaries (boundaries 0, 2, 3 and 4)
- the transverse velocity is pinned along the in- and outflow
  (boundaries 1 and 5).
.

\until pin_velocity


The applied traction elements require a pointer to the prescribed
traction function:

\until }

The wall elements require pointers to the various problem parameters
and the pointer to the \c GeomObject that defines the beam's
undeformed shape. Depending on the relative orientation of the fluid
and solid meshes, the normal vector used by the \c FSIWallElement
to determine the traction exerted by the adjacent \c FSIFluidElement
may either point into or out of the fluid domain. 
By default, it is assumed that the normal points into the
fluid -- in this case the traction computed by the
adjacent \c FSIFluidElement is added to the 
external load that acts on the \c FSIWallElement. 
Evaluating the direction of the normal (e.g by plotting the vector
obtained from \c FSIHermiteBeamElement::get_normal(...) ) in the 
present problem shows that the normal to the wall elements actually 
points out of the fluid domain, therefore the fluid traction acts in
the opposite direction to that assumed in the original
formulation. This may be rectified
by setting the boolean flag 
\c FSIHermiteBeamElement::normal_points_into_fluid() 
to \c false.


\until   } // end of loop over elements

Both ends of the beam are pinned:

\until }


We choose two fluid control nodes in the middle of the 
inflow and outflow cross-sections to document the velocities,
and choose a central node in the wall mesh to document its
displacement.

\until node_pt(Ncollapsible/2);

Finally, we set up the remaining fluid-structure interaction: 
The fluid nodes that are located on the moving
wall remain attached to material particles on the wall. The no-slip
condition (4) therefore implies that
the fluid velocity at each of these nodes must be equal to the 
nodes' velocity. Hence the fluid velocity
must be updated whenever a node update function
changes the nodal position. This is done most easily by means of
the auxiliary node update function -- a function that is executed
automatically whenever a  node's \c node_update() function
is called. To achieve this we pass a function pointer to the
\c FSI_functions::apply_no_slip_on_moving_wall() function
to the nodes on the fluid mesh's boundary 3: 

\until }

Next, the \c FSIHermiteBeamElements have to be "told" which fluid
elements are located next to their Gauss points to allow them
to work out the fluid traction. The required lookup tables are
created by the function 
\c FSI_functions::setup_fluid_load_info_for_solid_elements(...):

\until    (this,3,Bulk_mesh_pt,Wall_mesh_pt);

Finally, we set up the equation numbering scheme.

\until end of constructor


<HR> 
<HR>  

\section doc Post processing
The function \c doc_solution(...) outputs the velocity
and wall displacement fields and  records the time-trace of 
the axial velocity at the control nodes and the position of the wall's
midpoint. The function \c FSI_functions::doc_fsi(...)
is a helper function that can be used to document/validate the
various FSI lookup schemes; see \ref comments_and_ex for an
illustration of their output.

\dontinclude fsi_collapsible_channel.cc
\skipline start_of_doc_solution
\until {
\skipline Doc fsi
\skipline FSI_functions

\skipline ofstream
\until end_of_doc_solution

<HR>
<HR>

 \section tractioncre Creation of the traction elements
This function is the same as the one used in 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">the
problem with prescribed wall motion</A>. 

<HR>
<HR>
 

 \section IC Applying the initial conditions
This function is the same as the one used in 
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">the
problem with prescribed wall motion</A>. 

<HR>
<HR>
\section comments_and_ex Comments and Exercises

\subsection comments Comments


- <B>Variables that affect the fluid traction are not just velocities
  and pressures!</B> \n\n
  In the section \ref fsi_wall_elements we briefly discussed
  how the function 
  \c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
  determines the variables that affect the fluid traction 
  onto the \c FSIWallElements. For a Newtonian fluid, 
  the components of the fluid traction vector (on the viscous scale) 
  onto the wall is given by
  \f[
  t_i = -p N_i + \left(\frac{\partial u_i}{\partial x_j} + 
                       \frac{\partial u_j}{\partial x_i} \right) N_j
  \mbox{ \ \ \ for $i=1,2[,3]$} 
  \f]
  where \f$ N_i \ (\mbox{for } i=1,2[,3]) \f$ are the components of the
  normal to the fluid domain, pointing into the fluid. This equation
  shows that the fluid traction is primarily affected by the
  velocity and pressure degrees of freedom in the 
  fluid elements that are "adjacent" to a given \c FSIWallElement.
  \n\n
  However, since the traction involves derivatives of the
  velocity, the traction is also affected by changes to the  
  geometry of the fluid element. Therefore, the list of
  variables that affect the fluid traction must (and indeed does) 
  include all those \c Data objects that are involved in the adjacent 
  fluid elements' node update operations.
  \n\n
  Here is an animation that illustrates these dependencies for
  a relatively coarse discretisation in which the collapsible section
  of the fluid mesh is discretised with 10 "vertical columns" of 
  \c QCrouzeixRaviartElement<2> elements, while the wall is
  discretised with 22 \c FSIHermiteBeamElements, each of which
  contains 3 Gauss points. The animation shows the region of the fluid
  mesh close to the (strongly deformed) elastic wall.
  Each different frame illustrates the FSI lookup schemes
  for a different wall element. 
  - The position of the wall Gauss points are displayed by "gradient" 
    markers while the corresponding points in the adjacent fluid
    elements (i.e. the points at which the fluid traction is computed) 
    are displayed by 
    "delta" markers. Since the fluid and solid discretisations are 
    non-matching the points do not coincide exactly though they 
    will continue to approach each other under further mesh
    refinement. 
  - The coloured numbers indicate the number of values that affect
    the fluid traction on this \c FSIWallElement: 
    - The red numbers represent the number of nodal values at the
      nodes of the adjacent fluid element(s) that affect the traction: For a
      2D Crouzeix-Raviart element, each fluid node stores two velocity
      degrees of freedom, both of which affect the traction.  
    - The blue numbers represent the number of internal \c Data values
      stored in an adjacent fluid element that affect the traction:
      In a 2D Crouzeix-Raviart element, each element stores three
      pressure values in its internal \c Data and all three
      pressure values affect the traction. 
    - Finally, the green numbers indicate the number of \c Data values 
      that are (potentially) involved in the node-update operation for
      the nodes in the adjacent fluid elements. Recall that during the
      \c MacroElement - based node-update we only refer to the wall via
      its representation as a \c WallAsGeomObject, i.e. as a
      "compound" \c GeomObject. The geometric \c Data of a compound
      \c GeomObject is given by the geometric \c Data of \e all
      its sub-objects. Therefore, the geometric \c Data 
      of the \c WallAsGeomObject includes the positional
      \c Data of all \c SolidNodes stored in this mesh.
      In an \c FSIHermiteBeamElement, each \c SolidNode stores
      four values, representing the node's x- and y-positions
      and their derivatives with respect to the element's local coordinate. 
    .
  .
\image html crozier_raviart_fsi_macro.gif "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Crouzeix-Raviart elements.) " 
\image latex crozier_raviart_fsi_macro.eps "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Crouzeix-Raviart elements.) " width=0.75\textwidth

   Here is the corresponding animation for a discretisation with
   2D Taylor-Hood elements. These elements have no internal \c Data but
   the pressure degrees of freedom are stored at the fluid element's
   corner nodes:
\image html taylor_hood_fsi_macro.gif "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Taylor-Hood elements.) " 
\image latex taylor_hood_fsi_macro.eps "Animation of the Data values that affect the fluid traction that the adjacent fluid elements exert onto the various FSIHermiteBeamElements in the wall mesh. (The fluid elements are 2D Taylor-Hood elements.) " width=0.75\textwidth


Finally, here is an animation that shows the (wall) degrees of freedom
that affect the node-update of a given fluid node. 
The red square marker shows the fluid node; the green numbers
show the number of the degrees of freedom at the \c SolidNodes
that are involved that fluid node's node update. Again it is 
clear that the \c MacroElement/Domain-based node update 
procedure in which the wall mesh is represented by a compound
\c GeomObject does not result in a sparse node-update procedure:
Each \c SolidNode in the wall mesh is assumed to affect the
position of all fluid nodes.

\image html fsi_nodes.gif "Animation of the Data values that affect the node update of the fluid nodes. " 
\image latex fsi_nodes.eps "Animation of the Data values that affect the node update of the fluid nodes. " width=0.75\textwidth


<HR> 

- <B>(In-)efficiency of the MacroElement-based node-update</B> \n\n
  The animations shown above illustrate very graphically
  that the implementation of the fluid-structure interaction
  via \c MacroElement/Domain - based node updates does not lead to
  a particularly efficient algorithm. The current approach suffers
  from two main problems: 
  -# The fluid-node update is not sparse: Since we cannot distinguish between
     the various sub-objects in the "compound" \c GeomObject, we can do
     no better than assuming the worst-case scenario, namely 
     that \e all 
     positional degrees of freedom of \e all \c SolidNodes in the wall 
      mesh potentially affect the nodal position in the fluid elements
     adjacent to the wall. Consequently, each \c FSIHermiteBeamElement
     in the wall mesh depends on all solid mechanics degrees of 
     freedom in the wall mesh. As a result, the wall discretisation 
     completely loses its sparsity!
  -# When updating the position of the fluid nodes via the fluid element's
     \c Domain/MacroElement representation, we obtain the wall
     shape from the \c GeomObject::position(...) function of the
     "compound" \c GeomObject, using the wall's Lagrangian coordinate
     \f$ \xi \f$ as the "compound" \c GeomObject's intrinsic
     coordinate. As discussed in the section \ref wall_geom_object, this
     is a very costly operation, since we first have to determine
     which of the constituent \c FSIHermiteBeamElements "contains"
     the required Lagrangian coordinate, and then evaluate the
     Eulerian position of the relevant point in the element.
  .
  In <A HREF="../../../interaction/fsi_collapsible_channel_algebraic/html/index.html">
  the next  example</A> we will demonstrate how the use of the algebraic node update
  procedure, described in 
  <A HREF="../../../navier_stokes/algebraic_collapsible_channel/html/index.html">
  an earlier example</A>, allows us to avoid both problems, resulting
  in a much more efficient code.
    
<HR>
<HR>

\subsection ex Exercises
-# Use the function \c FSIHermiteBeamElement::get_normal(...) to plot 
   the unit normal vector to the wall and thus confirm that the value
   for \c FSIHermiteBeamElement::normal_points_into_fluid() is correct.
   [\b Hint: Another way to sanity-check that the correct value for this flag
   has been set is to change the velocity boundary conditions at the upstream 
   end to a pure Dirichlet condition by prescribing the axial velocity
   profile. With Dirichlet conditions everywhere, one fluid pressure degree 
   of freedom, \f$ p_{fix} \f$ , say, can (indeed must!) then be assigned 
   arbitrarily. Now set the inflow velocity to zero and increase the value of 
   \f$ p_{fix} \f$ from zero, say. If the wall collapses inwards
   as \f$ p_{fix} \f$ is increased the 
   direction of normal was chosen wrongly!] 
-# In section \ref fsi_wall_elements we mentioned that the function
   \c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
   assumes that the fluid nodes on the FSI boundary store 
   the boundary coordinate. Investigate what happens if this 
   step is omitted, e.g. by commenting out the assignment of boundary
   coordinates with \c Node::set_coordinates_on_boundary(...) in 
   <A HREF="../../../../src/meshes/collapsible_channel_mesh.template.cc">
   collapsible_channel_mesh.template.cc</A>. 
-# In section \ref undeformed we mentioned that the function 
   \c GeomObject::d2position(...) must be implemented for all \c
   GeomObjects that specify the undeformed shape of a beam element.
   Check what happens if this function is not implemented, e.g. by      
   commenting out its definition in the \c UndeformedWall class. 
-# In section \ref problemcontr we commented that \c oomph-lib
   does not require the discretisations of the fluid and solid
   meshes to match along the common boundary. Confirm this, e.g., by 
   increasing the number of elements in the wall mesh. 
.   



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/interaction/fsi_collapsible_channel/">
demo_drivers/interaction/fsi_collapsible_channel/
</A>
</CENTER>
- The driver code is: 
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


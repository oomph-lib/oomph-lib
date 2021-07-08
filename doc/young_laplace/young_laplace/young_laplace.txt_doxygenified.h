/**

\mainpage Example problem: The Young Laplace equation


This document discusses the finite-element-based solution of the 
<A HREF="http://en.wikipedia.org/wiki/Young-Laplace_equation">
Young Laplace equation</A>, a nonlinear PDE that determines the 
static equilibrium shapes of droplets or bubbles. We start by 
reviewing the relevant theory and then present the solution of a 
simple model problem.


<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
</CENTER>
This tutorial and the associated driver codes were developed jointly
with Cedric Ody (Ecole Polytechnique, Paris; now Rennes).
\n\n
</TD>
</TR>
</TABLE>
</CENTER>


\section theory Theory

The figure below illustrates a representative problem: A small 
droplet is extruded slowly from the outlet of a cylindrical tube. The 
air-liquid interface is pinned at the end of the tube. We assume that the 
size of the droplet and the rate of extrusion are so small that
gravitational, viscous and inertial forces may be neglected compared
to the interfacial forces acting at the air-liquid interface.

\image html spherical_cap.gif "A typical problem: The extrusion of a small droplet from a cylindrical tube. " 
\image latex spherical_cap.eps "A typical problem: The extrusion of a small droplet from a cylindrical tube. " width=0.6\textwidth

The shape of the air-liquid interface (the meniscus) is then determined
by  Laplace's law which expresses the balance between the (spatially
constant) pressure 
drop across the meniscus, \f$ \Delta p^* \f$, and the surface 
tension forces acting at the air-liquid interface,
\f[
\Delta p^* = \sigma \kappa^*,
\f]
where  \f$ \kappa^* \f$ is the mean curvature
and \f$ \sigma \f$ the surface tension.

Non-dimensionalising all lengths on some problem-specific lengthscale
\f$ {\cal L} \f$ (e.g. the radius of the cylindrical tube) and 
scaling the pressure on the associated capillary scale, 
\f$ \Delta p^* = \sigma/{\cal L} \ \Delta p \f$, yields 
the non-dimensional form of the Young-Laplace equation 
\f[
\Delta p =  \kappa.
\f]
This must be augmented by suitable boundary conditions, to be applied
at the contact line -- the line along which the meniscus
meets the solid surface. We can either "pin" the position of the 
contact line (as in the above example), or prescribe the  
<A HREF="http://en.wikipedia.org/wiki/Contact_angle">contact
angle</A> at which the air-liquid interface meets the solid surface.


The mathematical problem may be interpreted as follows: 
Given the prescribed pressure drop \f$ \Delta p \f$ (and thus 
\f$ \kappa \f$),  we wish to determine an interface shape of 
the required mean curvature \f$ \kappa \f$ that satisfies the 
boundary conditions along the contact line.

\subsection cartesian Cartesian PDE-based formulation

Expressing the (non-dimensional) height 
of the interface above the \f$ (x_1, x_2) \f$-plane as 
\f$ x_3 = u(x_1,x_2) \f$ yields the cartesian form of the 
Young-Laplace equation:
\f[
\Delta p =  \kappa = - \nabla \cdot \frac{\nabla u}{\sqrt{1+|\nabla u|^2}},
\f]
where \f$ \nabla = (\partial/\partial x_1,\partial/\partial x_2)^T \f$
Given the imposed pressure drop across the interface, this equation 
must be solved for \f$ u(x_1, x_2) \f$. Note that this is only possible
if the interface can be projected onto the \f$ (x_1, x_2) \f$-plane. 


\subsection pvd The principle of virtual displacements

The interface shape can also be determined
from the <A HREF="http://en.wikipedia.org/wiki/Virtual_displacement"> 
principle of virtual displacements</A>
\f[ 
\delta \left( \int_S  \mathit{\, d} s \right)= \int_S \Delta p {\bf N}
\cdot \delta {\bf N} \mathit{\, d} s + \oint_{L} \ \ {\bf T}_n \cdot
\delta {\bf R } \mathit{\, d} l \ \ \ \ \ \ \ \ \ \ \ \ (1)
\f]  
which may be derived from energetic considerations. Here the 
symbol \f$ \delta \f$ denotes a variation, \f$ {\bf R } \f$ is the 
position vector to the interface, and the vector \f$ {\bf N}
\f$ is the unit normal to the meniscus. The left
hand side of this equation represents the variation of the interfacial energy 
during a virtual displacement, and \f$ ds \f$ 
is an infinitesimal element of the meniscus surface, \f$S\f$. 
The terms on the right hand side represent the
virtual work done by the pressure and the virtual work done by the
surface tension forces acting at the free contact line, respectively. 
\f$ dl \f$ is the length of an element of the contact line
\f$L\f$, and \f$ {\bf T}_n \f$ is the vector tangent to the interface
and normal to the contact line. Note that, if the contact line is 
pinned, the variation of its position is zero, \f$ \delta {\bf R} = 
{\bf 0} \f$, and the line integral vanishes.


The variational formulation of the problem is particularly 
convenient for a finite-element-based solution; see, e.g. the
<a href="../../../solid/solid_theory/html/index.html#CartesianLagrangian">
Solid Mechanics Theory document</a> for a discussion of how
to discretise variational principles with finite elements. The
variational and PDE-based formulations are, of course, related to each
other: The Young-Laplace equation is the Euler-Lagrange equation of the
variational principle.


\subsection intrinsic Parametric representation 

To deal with cases in which the interface cannot be projected 
onto the  \f$ (x_1, x_2) \f$-plane, we parametrise the meniscus 
by two intrinsic coordinates as 
\f$ {\bf R}(\zeta_1,\zeta_2) \in { R}^3\f$, where 
\f$(\zeta_1,\zeta_2) \in D \in { R}^2\f$.  Furthermore, 
we parametrise the domain boundary, \f$\partial D\f$, by a 
scalar coordinate \f$\xi\f$ so that, 
\f[
{\partial D} = \bigg\{ (\zeta_1,\zeta_2) \ \bigg| \ 
 (\zeta_1,\zeta_2) = \left( \zeta_1^{[\partial D]}(\xi), \ 
\zeta_2^{[\partial D]}(\xi) \right) \bigg\}. 
\f] 
The normal to the meniscus is then given by 
\f[ 
{\bf N} = \frac{{\bf R}_{,1} \times {\bf R}_{,2} }
{\mathcal{A}^{1/2}}, 
\f] 
where the commas denote partial differentiation with respect to the
intrinsic coordinates, and 
\f$ \mathcal{A}^{1/2}= |{\bf R}_{,1} \times {\bf R}_{,2}| \f$ 
is the square root of the surface metric tensor.

The area and length differentials required in variational principle
are
\f[ 
d s=\mathcal{A}^{1/2} d \zeta_1 d \zeta_2 
\f] 
and 
\f[ 
d l=\left| \frac{d {\bf R}\left( \zeta_1^{[\partial D]}(\xi), \ 
\zeta_2^{[\partial D]}(\xi) \right)}{d \xi} \right| d \xi,
\f] 
allowing us to write the principle of virtual displacements as
\f[ 
\int_S \left( \frac{\delta \mathcal{A}^{1/2}}{\mathcal{A}^{1/2}} -
\kappa \ {\bf N} \cdot \delta {\bf R} \right) \mathcal{A}^{1/2} d \zeta_1 d \zeta_2 =
\oint_{L} {\bf T}_n \cdot \delta {\bf R} \ 
\left| \frac{d {\bf R}\left( \zeta_1^{[\partial D]}(\xi), \ 
\zeta_2^{[\partial D]}(\xi) \right)}{d \xi} \right| d \xi.
\f]

\subsection spines The method of spines
In the current form, the variational principle cannot yield
a unique solution because there are infinitely many 
vector fields  \f$ {\bf R}(\zeta_1,\zeta_2) \f$ that parametrise
the same interface shape. To remove this ambiguity, and to 
allow for interface shapes that cannot be projected onto the
\f$ (x_1,x_2) \f$-plane, we employ the so-called "Method of Spines".
This method was originally proposed by Kistler and Scriven
for the computation of free surface flows. We decompose 
the vector \f$ {\bf R} \f$ into two parts by writing it as 
\f[ 
{\bf R}(\zeta_1,\zeta_2) = {\bf B}(\zeta_1,\zeta_2) + 
\mathit{u}(\zeta_1,\zeta_2) \ {\bf S}(\zeta_1,\zeta_2). 
\f] 
Here the "spine basis" \f${\bf B}(\zeta_1,\zeta_2)\f$
and the "spines" \f${\bf S}(\zeta_1,\zeta_2)\f$ are pre-determined vector
fields that must be chosen by the user. Using this
decomposition, the meniscus shape is determined by the scalar function
\f$\mathit{u}(\zeta_1,\zeta_2)\f$ which represents the meniscus' displacement 
along the spines \f${\bf S}\f$.

The idea is illustrated in the simple 2D sketch below: the spine
basis vectors, \f$ {\bf B}(\zeta), \f$  parametrise the straight line 
corresponding to a flat meniscus. Positive values of \f$ u \f$
displace the meniscus along the spines, \f$ {\bf S}(\zeta) \f$ (the red 
vectors), whose orientation allows the representation  of meniscus 
shapes that cannot be represented in cartesian form as 
\f$ x_2 = u(x_1). \f$


\image html spine_sketch.gif "Sketch illustrating the parametrisation of the meniscus by the Method of Spines. " 
\image latex spine_sketch.eps "Sketch illustrating the parametrisation of the meniscus by the Method of Spines. " width=0.35\textwidth


The spine basis and the spines themselves must be chosen such that the 
mapping from \f$(\zeta_1,\zeta_2)\f$ to \f${\bf R}(\zeta_1,\zeta_2)\f$ is
one-to-one, at least for the meniscus shapes
of interest. Pinned boundary conditions of the form
\f$ \left. {\bf R}\right|_{\partial D} = {\bf R}_{\rm pinned} \f$
are most easily imposed by choosing the spine basis such that \f$ \left. {\bf
B}\right|_{\partial D} = {\bf R}_{\rm pinned}, \f$ implying
that \f$ \left. u\right|_{\partial D} = 0. \f$ 

The simplest possible choice for the spines and spine basis
is one that returns the problem to its original cartesian formulation.
This is achieved by setting \f$ {\bf B}(\zeta_1,\zeta_2) = 
(\zeta_1,\zeta_2,0)^T \f$ and \f$ {\bf S}(\zeta_1,\zeta_2) = (0,0,1)^T. \f$


\subsection displ Displacement Control
The Young-Laplace equation is a highly nonlinear PDE. We therefore 
require a good initial guess for the solution in order to ensure the
convergence of the Newton iteration. In many cases good initial guesses can be
provided by a simple, physically motivated continuation method.
For instance in the model problem shown above, the computation was
started by computing the solution for \f$ \Delta p = \kappa = 0 \f$ --
a flat interface. This solution was then used as the initial guess for 
the solution at a small positive value of \f$ \Delta p \f$. This process
was continued, allowing us to compute strongly deformed
meniscus shapes. This method works well, provided small increments in 
the control parameter \f$ \Delta p \f$ (or equivalently, \f$ \kappa \f$)
create small changes in the interface shape. This is not always
the case, however, as we shall see in the example below. In 
such cases it is often possible to re-formulate the problem, using the 
displacement control method discussed in the
<a href="../../../beam/steady_ring/html/index.html">solid
mechanics tutorials</a>. Rather than prescribing the pressure 
drop \f$ \Delta p \f$ we prescribe the displacement of a control 
point on the meniscus and
regard the pressure drop required to achieve this displacement
as an unknown. Since the implementation of the method is very similar
to that used for solid mechanics problems, we shall not
discuss it in detail here but refer to the appropriate 
<a href="../../../beam/steady_ring/html/index.html">solid mechanics 
tutorial.</a>


<HR>
<HR>

\section barrel An example problem: A barrel-shaped meniscus
As an example, we consider the following problem: Fluid is extruded
from an infinitely long, parallel-sided slot of width \f$ 2a \f$.
If we assume that the air-liquid interface is pinned at the edges of
the slot, as shown in the sketch below, the problem has an 
obvious exact solution. The air-liquid interface 
must have constant mean curvature, so, assuming that its shape does 
not vary along the slot, the meniscus must be a circular
cylinder. If we characterise the meniscus' shape by its vertical
displacement along the centreline, \f$ H \f$, the (dimensional) 
curvature of the air-liquid interface is given by
\f[
\kappa^*  = \frac{2/a}{H/a + a/H}.
\f] 
The plot of this function, shown in the right half of the figure
below, may be interpreted as a "load-displacement
diagram" as it shows the deflection of the meniscus as a function 
of the imposed non-dimensional pressure drop across the air-liquid interface, 
\f$ \Delta p = a \kappa^*. \f$
 
\image html meniscus_sketch_and_limit_point.gif "Sketch of the meniscus above a slot of width 2a. The figure on the right illustrates the relation between the meniscus curvature and its vertical displacement along the centreline. " 
\image latex meniscus_sketch_and_limit_point.eps "Sketch of the meniscus above a slot of width 2a. The figure on the right illustrates the relation between the meniscus curvature and its vertical displacement along the centreline. " width=0.75\textwidth

The plot shows that the load-displacement curve is 
not single-valued. Furthermore, the presence of a limit point 
at \f$ a \kappa^* = 1 \f$ indicates that the maximum curvature of 
the meniscus is given by \f$ \kappa_{max}^* = 1/a. \f$ This implies 
that the maximum (dimensional) pressure that the meniscus can
withstand is given by \f$ \Delta p^*_{max} = \sigma \kappa_{max}^* = 
\sigma/a \f$.

The presence of the limit point makes it impossible to compute the entire 
solution curve by simply increasing \f$ \kappa \f$ in small increments. 
However, use of a displacement control approach, by prescribing 
\f$ H \f$ while regarding \f$ \Delta p \f$ (and 
thus \f$\kappa\f$) as an unknown, yields a single-valued function 
that can be computed by a straightforward continuation method in 
which \f$ H \f$ is increased in small increments. 
 
This approach was employed to compute the meniscus shapes shown 
in the plot below. The initial, zero-curvature configuration
of the meniscus (corresponding to a vanishing pressure drop across the 
air-liquid interface) is the unit square. The meniscus is pinned along 
the lines \f$ y=0 \f$ and \f$ y=1 \f$, and symmetry (natural) 
boundary conditions were applied along the lines \f$ x=0 \f$ and \f$
x=1; \f$ see \ref comm_ex for further details on the natural 
boundary conditions.
As the pressure drop increases, the meniscus is deflected upwards
until it reaches the configuration of maximum curvature when 
\f$ H/a=1 \f$. Beyond this point, the pressure drop has to be reduced
to compute the "bulging" meniscus shapes obtained for \f$ H/a > 1. \f$

\image html barrel.gif "Deformation of a meniscus that is pinned along the lines y=0 and y=1. " 
\image latex barrel.eps "Deformation of a meniscus that is pinned along the lines y=0 and y=1. " width=0.75\textwidth

The comparison of the computed "load-displacement curve" 
against the exact solution, shown below, indicates that the two agree
to within plotting accuracy. 


\image html trace.gif "Load-displacement diagram for a meniscus that is pinned along the lines y=0 and y=1. " 
\image latex trace.eps "Load-displacement diagram for a meniscus that is pinned along the lines y=0 and y=1. " width=0.75\textwidth
 
<HR> 
<HR> 

\section impl Implementation

We shall now discuss the solution of the above problem with
\c oomph-lib's  Young Laplace elements.

\subsection namespace The global namespace

As usual we define the problem parameters in a global namespace.
The key parameter is the value of the prescribed control displacement
which we initialise to zero.  The function \c get_exact_kappa()
returns the exact solution for the mean curvature of the meniscus.
We will use this function for the validation of our results.

\dontinclude barrel.cc
\skipline start_of_namespace
\until end exact kappa 

Next we define the orientation of the spines. Anticipating the shape
of the meniscus, we set \f$ {\bf B}(\zeta_1,\zeta_2)=
(\zeta_1,\zeta_2,0) \f$ so that \f$ u=0 \f$ corresponds to a flat
meniscus in the \f$ (x,y) \f$-plane,

\until bspine functions

and rotate the spines in the \f$ y \f$ direction by setting
\f[ 
{\bf S}(\zeta_1,\zeta_2) = 
\left(
\begin{array}{c}
0 \\
-\sin\big(\alpha(\zeta_2)\big) \\
\ \ \cos\big(\alpha(\zeta_2)\big) \\
\end{array}
\right)
\f]
where \f$ \alpha(\zeta_2) = \alpha_{min} + (\alpha_{max}-\alpha_{min})
\zeta_2 \f$. With this choice, the spines are unit vectors that
form an angle \f$ \alpha \f$ (varying between \f$ \alpha_{min} \f$ and
\f$ \alpha_{max} \f$) with the y-axis.

\until spine function

 
<HR>
<HR>

\subsection main The driver code
We start by preparing an output directory and open a trace file to 
record the control displacement, and the computed and
exact values of the interface curvature \f$ \kappa \f$. 

\dontinclude barrel.cc
\skipline start_of_main 
\until ZONE

Next we build the problem object and document the initial
configuration: a flat meniscus. 

\until doc_info.number()++;

Finally, we perform a parameter study by increasing the control
displacement in small increments and re-computing the 
meniscus shapes and the associated interface curvatures.

\until end of main


 
<HR>
<HR>

\subsection class The problem class

The problem class has the usual member functions. (Ignore the lines
in \c actions_before_newton_solve() as they are irrelevant in the
current context. They are discussed in one of the exercises in 
\ref comm_ex.) The problem's private member
data include a pointer to the node at which the meniscus displacement
is controlled by the displacement control element, and a pointer to 
the \c Data object whose one-and-only value contains the unknown
interface curvature, \f$ \kappa. \f$

\dontinclude barrel.cc
\skipline start_of_problem_class
\until end of problem

 
<HR>
<HR>

\subsection constr The problem constructor

We start by building the mesh, discretising the two-dimensional
parameter space \f$ (\zeta_1,\zeta_2) \in [0,1] \times [0,1]\f$ 
with 8x8 elements.

\until mesh_pt()

Next, we choose the central node in the mesh as the node whose
displacement is imposed by the displacement control method. 

\until Control_node_pt->x(1)
 
We pass the pointer to that node and the pointer to the double
that specifies the imposed displacement to the constructor
of the displacement control element. The constructor automatically creates
a \c Data object whose one-and-only value stores the unknown
curvature, \f$ \kappa. \f$ We store the pointer to this \c Data
object in the private member data to facilitate its output.

\until Kappa_pt=

The meniscus is pinned along mesh boundaries 0 and 2:

\skipline Boundary conditions
\until end bc

We complete the build of the Young Laplace elements by passing the
pointer to the spine functions, and the prescribed curvature. 

\until }

Finally, we add the displacement control element to the mesh
and assign the equation numbers.

\until end of constructor


<HR>
<HR>

\subsection post Postprocessing
We document the exact and computed meniscus curvatures in the
trace file and output the meniscus shape.

\skipline start_of_doc
\until end of doc


<HR>
<HR>

\section comm_ex Comments and Exercises

-# <b>Choice of spines:</b> \n\n
   We discussed earlier that the spine basis and the spines themselves 
   must be chosen such that the mapping from \f$(\zeta_1,\zeta_2)\f$
   to \f${\bf R}(\zeta_1,\zeta_2)\f$ is one-to-one, at least for 
   the meniscus shapes of interest. This requires some prior knowledge
   of the expected interface shapes.  \n \n
   The spine basis and the spines must be defined via function pointers
   that are passed to the Young Laplace elements. If the function
   pointers are not specified, the Young-Laplace elements revert to 
   a cartesian formulation. \n\n 
   Experiment with different spine orientations and explore the
   interface shapes that are obtained if no spines are specified (by commenting
   out the lines in the constructor that pass the relevant function
   pointers to the Young Laplace elements).  
   \n\n
-# <b>Natural boundary conditions:</b> \n\n
   As usual in any finite-element computation, we only enforced the 
   essential boundary conditions by pinning the meniscus displacement
   along the "pinned contact lines" at \f$ y=0 \f$ and  \f$ y=1. \f$
   No constraints were applied along the two other domain boundaries
   (at \f$ x=0 \f$ and  \f$ x=1 \f$), indicating that these boundaries
   are controlled by implied, "natural" boundary conditions. The
   variational principle (1) shows what these are: since
   we neglected the boundary integral on the right hand side of
   equation (1), 
   the meniscus shape must satisfy 
   \f$ {\bf T}_n \cdot \delta {\bf R} = \delta u \ {\bf T}_n \cdot
   {\bf S} = 0, \f$
   implying that outer unit normal to the meniscus boundary, 
   \f$ {\bf T}_n, \f$ must be orthogonal to the direction of the spines.
   Since the spines do not have an
   \f$ x \f$ - component, the meniscus must therefore have
   zero slope in that direction -- just what we need for our problem. 
   \n\n
   To convince yourself that the this argument is correct, rotate 
   the spines in the  \f$ x \f$ -direction, e.g. by changing
   their definition to
   \f[ 
   {\bf S}(\zeta_1,\zeta_2) = 
   \left(
   \begin{array}{c}
   1/2 \\
   -\sin\big(\alpha(\zeta_2)\big) \\
   \ \ \cos\big(\alpha(\zeta_2)\big) \\
   \end{array}
   \right).
   \f]
   The natural boundary condition will still force the meniscus to be 
   normal to the spines along the "free" contact line, resulting in 
   interface shapes similar to the one shown in the figure below.
\image html funky_meniscus.gif "Meniscus shape created by the natural boundary conditions when the spines (shown as vectors) are rotated in the x-direction. " 
\image latex funky_meniscus.eps "Meniscus shape created by the natural boundary conditions when the spines (shown as vectors) are rotated in the x-direction. " width=0.6\textwidth
   This demonstrates yet again that the orientation of spines must 
   reflect the relevant features of the problem. 
   We refer to <a href="../../contact_angle/html/index.html">another
   tutorial</a> for a more detailed discussion 
   of the boundary condition and its relation to contact angles.
   \n\n
-# <b>Displacement control:</b> \n\n
   Explore what happens if you disable displacement control and
   prescribe the pressure drop (i.e. \f$ \kappa \f$) directly. 
   The relevant code is already contained in the driver code.
   You'll have to comment out the lines in the problem constructor 
   that create the displacement control element and the line that adds 
   it to the problem's mesh.  Replace them by the lines
   \n\n
   \dontinclude barrel.cc
   \skipline Comment out the pre
   \until pin(0);
   \n\n
   which create the \c Data object that stores the prescribed curvature. 
   (Note that the value of \f$ \kappa \f$ is already incremented in \c
   actions_before_newton_step(). With displacement control this 
   step has no real effect as the Newton method will overwrite this
   assignment). Check what happens if the prescribed
   curvature exceeds the maximum possible curvature of the meniscus. 
   \n\n
-# <b>Inefficient implementation:</b> \n\n
   Note that the current implementation of the Young Laplace elements   
   is inefficient as the elemental Jacobian matrices are computed by
   finite-differencing. You are invited to implement the analytical
   computation of the Jacobian matrix as an exercise.
   \n\n
-# <b>Other problems:</b> \n\n
   We provide a number of additional demo driver codes that
   demonstrate the solution of other, related problems.
   \n\n 
   - The code \n\n
     <CENTER>
     <A HREF="../../../../demo_drivers/young_laplace/young_laplace.cc">
     demo_drivers/young_laplace/young_laplace.cc
     </A>
     </CENTER>
     \n\n
     and its adaptive counterpart
     \n\n
     <CENTER>
     <A HREF="../../../../demo_drivers/young_laplace/refineable_young_laplace.cc">
     demo_drivers/young_laplace/refineable_young_laplace.cc
     </A>
     </CENTER>
     \n\n
     demonstrate the solution of three problems: (i) the barrel-shaped meniscus
     problem already discussed above; (ii) the deformation of a meniscus that 
     is pinned at all four edges of a square tube; and (iii) the solution of
     a problem with contact-angle boundary conditions. The latter one
     is discussed in a <a href="../../contact_angle/html/index.html">
     separate tutorial.</a>     
     \n\n
   - The spherical meniscus that emanates from a circular tube, shown 
     in the animation at the beginning of current document, was computed with:
      \n\n
     <CENTER>
     <A HREF="../../../../demo_drivers/young_laplace/spherical_cap_in_cylinder.cc">
     demo_drivers/young_laplace/spherical_cap_in_cylinder.cc
     </A>
     </CENTER>
     \n\n
   .
.

 
<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
  <CENTER>
  <A HREF="../../../../demo_drivers/young_laplace/">
  demo_drivers/young_laplace/
  </A>
  </CENTER>\n
- The driver code is: \n\n
  <CENTER>
  <A HREF="../../../../demo_drivers/young_laplace/barrel.cc">
  demo_drivers/young_laplace/barrel.cc
  </A>
  </CENTER>
.


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Demo problem: Deformation of a string under tension, using Kirchhoff-Love Beam elements.

In this document, we discuss the solution of a simple one-dimensional
solid mechanics problem, using \c oomph-lib's Kirchhoff-Love beam elements:

<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>One-dimensional beam problem</B>
</CENTER> 
Determine the deformation of a pre-stressed, pressure-loaded elastic 
beam of length \f$ L^*
\f$ and wall thickness \f$ h^* \f$ . In its undeformed configuration,
the beam is subject to an axial (2nd Piola-Kirchhoff) pre-stress of 
size \f$ \sigma^*_0 \f$ . We wish to compute its deformation when it
is loaded by a spatially uniform pressure of size \f$ p^*_{ext} \f$ .

\image html string_sketch.gif "A pre-stressed elastic beam under constant pressure loading. " 
\image latex string_sketch.eps "A pre-stressed elastic beam under constant pressure loading. " width=0.75\textwidth
</TD> 
</TR>
</TABLE>  
</CENTER>  

 
<HR> 
  
\section non_dim Theory and non-dimensionalisation

\c oomph-lib's beam elements are based on geometrically 
nonlinear Kirchhoff-Love beam theory with incrementally linear 
constitutive equations (Young's modulus \f$ E \f$ and Poisson's ratio 
\f$ \nu \f$). The equations are implemented in non-dimensional form,
obtained by non-dimensionalising all length on some lengthscale 
\f$ {\cal L} \f$, by scaling the stresses and the applied traction on the
beam's effective Young's modulus \f$ E_{eff} = E/(1-\nu^2) \f$
and by non-dimensionalising time on some timescale \f$ {\cal T} \f$ so that
the dimensional (identified by an asterisk) and non-dimensional 
variables are related by
\f[
L^* = {\cal L} \ L,
\f]
\f[
h^* = {\cal L} \ h,
\f]
\f[
t^* = {\cal T} \ t,
\f]
\f[
p_{ext}^* = E_{eff} \ p_{ext},
\f]
and 
\f[
\sigma_0^* = E_{eff} \ \sigma_0.
\f]
The beam's undeformed shape is parametrised by a non-dimensional Lagrangian 
coordinate \f$ \xi = \xi^*/{\cal L} \f$ so that the non-dimensional 
position vector to a material particle on the beam's centreline in the 
undeformed 
configuration is given by \f$ {\bf r}_w(\xi) \f$ . We denote the
unit normal to the beam's undeformed centreline by \f$ {\bf n} \f$. 
The applied traction \f$ {\bf f} = {\bf f}^*/E_{eff} \f$ 
(a force per unit deformed length of the beam) deforms the beam, 
causing its material particles to be displaced to their 
new positions \f$ {\bf R}_w(\xi) \f$; the unit
normal to the beam's deformed centreline is \f$ {\bf N} \f$.


The non-dimensional form of the principle of virtual displacements
that governs the beams deformation is then given by
\f[
\int_0^{L}   \left[ 
(\sigma_0 + \gamma) \ \delta  
\gamma +  
\frac{1}{12} h^2 \kappa 
 \ \delta \kappa   - 
\left(\frac{1}{h} \sqrt{\frac{A}{a}} \   {\bf f} - \Lambda^2 
\frac{\partial^2 {\bf R}_w}{\partial t^2} \right) \cdot  
\delta {\bf R}_w 
 \right] \ \sqrt{a} \ d\xi  = 0,
\ \ \ \ \ \ \ \ \  (1)
\f]
where 
\f[
a =   \frac{\partial {\bf r}_w}{\partial \xi}
\cdot \frac{\partial {\bf r}_w}{\partial \xi},
\f]
and
\f[
A =   \frac{\partial {\bf R}_w}{\partial \xi}
\cdot \frac{\partial {\bf R}_w}{\partial \xi},
\f]
represent the squares of the lengths of infinitesimal material line elements
in the undeformed and the deformed configurations, respectively.
\f[
ds = \sqrt{a} \ d\xi \mbox{\ \ \ \ and \ \ \ } dS = \sqrt{A} \ d\xi.
\f]
\f$ A \f$ and  \f$ a \f$  may be interpreted as the "1x1 metric
tensors" of the beam's centreline, in the deformed and undeformed
configurations, respectively. The ratio \f$ \sqrt{A/a} \f$ represents
the "extension ratio" or the "stretch" of the beam's centreline. 

 We represent the curvature of the beam's centreline before and 
after the deformation by
\f[
b = {\bf n} \cdot \frac{\partial^2 {\bf
r}_w}{d\xi^2}
\f]
and
\f[
B = {\bf N} \cdot \frac{\partial^2 {\bf
R}_w}{d\xi^2}
\f]
respectively. The ("1x1") strain and and bending "tensors" \f$ \gamma \f$ and
\f$ \kappa \f$ are then given by
\f[
\gamma = \frac{1}{2}\left(A-a\right)
\mbox{\ \ \  and  \ \ \ }
\kappa = - \left( B-b  \right).
\f]
Finally,
\f[
\Lambda = \frac{{\cal L}}{{\cal T}} \sqrt{\frac{\rho}{E_{eff}}}
\f] 
is the ratio of the natural timescale of the beam's in-plane
extensional oscillations,
\f[
{\cal T}_{natural} = {\cal L} \sqrt{\frac{\rho}{E_{eff}}},
\f] 
to the timescale \f$ {\cal T}\f$ used in the non-dimensionalisation of the
equations. \f$ \Lambda^2 \f$ may be interpreted as
the non-dimensional wall density, therefore \f$ \Lambda=0 \f$
corresponds to the case without wall inertia. 


\anchor hermite
\c oomph-lib's \c HermiteBeamElement provides a discretisation of the
variational principle  (1) with one-dimensional, isoparametric, 
two-node Hermite solid mechanics elements. In these elements, the Eulerian
positions of the nodes, accessible via \c Node::x(...), are regarded 
as unknowns, and Hermite interpolation is used to interpolate 
the position between the nodes, so that the Eulerian position of material
points within an element is given by
\f[
x_i(s) =  \sum_{j=1}^2  \sum_{k=1}^2 X_{ijk} \ \psi_{jk}(s) 
\mbox{\ \ \ for $i=1,2$}, 
\ \ \ \ \ \ \ \ (2)
\f]
where \f$ s \in [-1,1] \f$ is the element's 1D local coordinate.
The functions \f$ \psi_{jk}(s) \f$ are the one-dimensional
Hermite shape functions
\f[ \psi_{11}(s) =  \frac{1}{4}\left(s^3 - 3s + 2\right), \f]
\f[ \psi_{12}(s) =  \frac{1}{4}\left(s^3 - s^2 - s +1 \right), \f]
\f[ \psi_{21}(s) =  \frac{1}{4}\left(2 + 3 s - s^3\right), \f]
\f[ \psi_{22}(s) =  \frac{1}{4}\left(s^3 + s^2 -s -1 \right). \f]
They have the property that 
- \f$ \psi_{j1} = 1 \f$ at node \f$ j \f$ and  \f$ \psi_{j1} = 0 \f$ 
  at the other node. Furthermore,  \f$ d\psi_{j1}/ds = 0 \f$ at 
  both nodes.
- \f$ d\psi_{j2}/ds = 1 \f$ at node \f$ j \f$ and  \f$ d\psi_{j2}/ds = 0 \f$ 
  at the other node. Furthermore,  \f$ \psi_{j2} = 0 \f$ at 
  both nodes.
.
The mapping (2) therefore provide an
independent interpolation for the position and its derivative
with respect to the local coordinate \f$ s \f$. As a result
we have two types of (generalised) nodal coordinates:
-  \f$ X_{ij1} \f$ represents the \f$ i \f$-th coordinate of the element's
   local node \f$ j \f$. 
-  \f$ X_{ij2} \f$ represents the derivative of \f$ i \f$-th coordinate 
   with respect to the local coordinate \f$ s \f$, evaluated
   at the element's local node \f$ j \f$. 
.
This representation ensures the
\f$ C^1 \f$ continuity of the wall shape, required by
the variational principle (1) which contains second derivatives
of \f$ {\bf R}_w \f$. Physically, the inter-element continuity 
of the slope is enforced by the beam's nonzero bending stiffness. 

\anchor solid_bound
The two "types" of "generalised" positional degrees of freedom
are accessible via the function \c Node::x_gen(k,i), where
\c Node::x_gen(0,i) \f$ \equiv \f$ \c Node::x(i). 
The two types of positional degrees of freedom correspond to 
the two types of boundary conditions (pinned and clamped) that can be applied
at the ends of the beam. Mathematically, the number of boundary
conditions reflects the fact that the Euler-Lagrange equations 
of the variational principle are of fourth-order.


The nodes of solid mechanics elements are \c SolidNodes,
a generalisation of the basic \c Node class that allows
the nodal positions to be treated as unknowns. Its member 
function \c SolidNode::pin_position(...) allows the application
of boundary conditions for the (generalised) nodal positions. The
table below lists several common boundary conditions and their application:
<TABLE BORDER=1>
<TR>
<TD>
Boundary condition
</TD>
<TD>
Mathematical condition
</TD>
<TD>
Implementation
</TD>
</TR>
<TR>
<TD>
<B>Pinned:</B>
<img src="../figures/pinned.gif">
</TD>
<TD>
\f$ {\bf R}_w \cdot {\bf e}_x = const. \f$ \n
\f$ {\bf R}_w \cdot {\bf e}_y = const. \f$ 
</TD>
<TD>
\c SolidNode::pin_position(0);\n 
\c SolidNode::pin_position(1); \n \n
or, equivalently, \n \n
\c SolidNode::pin_position(0,0);\n 
\c SolidNode::pin_position(0,1); \n
</TD>
</TR>
<TR>
<TD>
<B>Pinned, sliding in the x-direction:</B>
<img src="../figures/pinned_sliding_x.gif">
</TD>
<TD>
\f$ {\bf R}_w \cdot {\bf e}_y = const. \f$ 
</TD>
<TD>
\c SolidNode::pin_position(1); \n \n
or, equivalently, \n \n
\c SolidNode::pin_position(0,1); \n
</TD>
</TR>
<TR>
<TD>
<B>Pinned, sliding in the y-direction:</B>
<img src="../figures/pinned_sliding_y.gif">
</TD>
<TD>
\f$ {\bf R}_w \cdot {\bf e}_x = const. \f$ 
</TD>
<TD>
\c SolidNode::pin_position(0); \n \n
or, equivalently, \n \n
\c SolidNode::pin_position(0,0); \n
</TD>
</TR>
<TR>
<TD>
<B>Clamped:</B>
<img src="../figures/clamped.gif">
</TD>
<TD>
\f$ {\bf R}_w \cdot {\bf e}_x = const. \f$\n
\f$ {\bf R}_w \cdot {\bf e}_y = const. \f$ \n
\f$ d({\bf R}_w \cdot {\bf e}_y)/d\xi = 0. \f$ 
</TD>
<TD>
\c SolidNode::pin_position(0); \n 
\c SolidNode::pin_position(1); \n 
\c SolidNode::pin_position(1,1); \n \n
or, equivalently, \n \n
\c SolidNode::pin_position(0,0); \n 
\c SolidNode::pin_position(0,1); \n 
\c SolidNode::pin_position(1,1); \n
</TD>
</TR>
<TR>
<TD>
<B>Clamped, sliding in the y-direction\n
(symmetry boundary condition!)</B>
<img src="../figures/clamped_sliding_y.gif">
</TD>
<TD>
\f$ {\bf R}_w \cdot {\bf e}_x = const. \f$\n
\f$ d({\bf R}_w \cdot {\bf e}_y)/d\xi = 0. \f$ 
</TD>
<TD>
\c SolidNode::pin_position(0); \n 
\c SolidNode::pin_position(1,1); \n \n
or, equivalently, \n \n
\c SolidNode::pin_position(0,0); \n 
\c SolidNode::pin_position(1,1); \n
</TD>
</TR>
</TABLE> 


The \c HermiteBeamElement provides default values for all
non-dimensional physical parameters:
- the non-dimensional 2nd Piola Kirchhoff pre-stress \f$ \sigma_0 \f$ is zero.
- the non-dimensional beam thickness \f$ h \f$ is 1/20.
- the timescale ratio \f$ \Lambda^2 \f$ is 1.
- the non-dimensional traction vector \f$ {\bf f} \f$ evaluates to zero.
.
These values can be over-written via suitable access functions. 
[Time-dependent computations also require the specification of
a timestepper for the elements. This is demonstrated in 
<a href="../../../beam/unsteady_ring/html/index.html"> another example.</a>]
The "user" must specify:
- the undeformed wall shape \f$ {\bf r}_w(\xi) \f$ as a \c GeomObject
  -- see the <a href="../../../poisson/fish_poisson2/html/index.html">
  earlier example</a> for a discussion of \c oomph-lib's geometric
  objects.
.

<HR>

\section results Results
The animation shown below illustrates the deformation of the beam
under increasing pressure. An increase in pressure initially deflects 
the beam vertically downwards. The pressure acts as a "follower
load" since it always acts in the direction normal to the deformed
beam. This causes the beam to deform into an approximately circular shape 
whose radius increases rapidly. 

\image html string_animation.gif "Deformation of a pre-stressed elastic beam under constant pressure loading. " 
\image latex string_animation.eps "Deformation of a pre-stressed elastic beam under constant pressure loading. " width=0.75\textwidth


<HR>

\section approx An approximate analytical solution

Before discussing the details of the numerical solution, we will derive an
approximate analytical solution of the problem. The analytical solution
will be useful to validate the computational results, and the derivation
will allow us to discuss certain aspects of the theory in more detail.

We start by parametrising the beam's undeformed shape as
\f[
{\bf r}_w(\xi) = 
\left(
\begin{array}{c}
\xi \\
0
\end{array}
\right)
\mbox{ \ \ \ \ where $\xi  \in [0,L].$}
\f]
The 1x1 "metric tensor" associated with this parametrisation is given by
\f[
a = \frac{{\bf r}_w(\xi)}{d\xi} \cdot \frac{{\bf r}_w(\xi)}{d\xi} = 1,
\f]
consistent with the fact that the Lagrangian coordinate \f$ \xi \f$ 
is the arclength along the beam's undeformed centreline.


If the beam is thin (so that \f$ h \ll 1 \f$) bending effects
will be confined to thin boundary layers near its
ends. The beam will therefore behave (approximately)
like a "string under tension" and its deformed shape will
be an arc of a circle. All material line elements will be stretched
by the same amount so that the tension is spatially uniform. The position 
vector to the material particles on the beam's deformed centreline is 
therefore given by
\f[
{\bf R}_w(\xi) = 
\left(
\begin{array}{c}
\frac{1}{2}\left(1+\frac{\sin(-\alpha/2+\alpha\xi/L)}
{\sin(\alpha/2)}\right)  \\
\frac{\cos(-\alpha/2+\alpha\xi/L)-\cos(\alpha/2)}
{2\sin(\alpha/2)}
\end{array}
\right)
\mbox{ \ \ \ \ where $\xi  \in [0,L].$}
\ \ \ \ \ \ \ \ \  (3)
\f]
(See the \ref comments for a more detailed discussion of this
analytical solution.) Here \f$ \alpha \f$ is the opening angle 
of the circular arc as shown in this sketch:
 

\image html string.gif "Sketch illustrating the approximate analytical solution for the problem. " 
\image latex string.eps "Sketch illustrating the approximate analytical solution for the problem. " width=0.75\textwidth

This deformation generates a uniform stretch of
\f[
\frac{dS}{ds} = \sqrt{\frac{A}{a}} = \frac{1}{2} \frac{\alpha}{\sin(\alpha/2)},
\f]
corresponding to a uniform strain 
\f[
\gamma = \frac{1}{2}\left(\frac{1}{4}\frac{\alpha^2}{\sin^2(\alpha/2)}
- 1 \right).
\f]
The incremental Hooke's law predicts a linear relation between
the (2nd Piola Kirchhoff) stress, \f$ \sigma \f$, and the
(Green) strain, \f$ \gamma \f$, so that
\f[
\sigma = \sigma_0 + \gamma,
\f]  
where \f$ \sigma_0 \f$ is the axial pre-stress that acts in the undeformed
configuration.

An elementary force balance shows that the pressure 
\f$ p_{ext} \f$ required to deform the beam into the
shape specified by (3), is given by "Laplace's law"  
\f[ 
p_{ext} =\frac{{\tt T}}{R} 
\f]
where 
\f[ 
R = \frac{L}{2 \sin(\alpha/2)}
\f] 
is the radius of the circular arc,
and the axial tension \f$ {\tt T} \f$ is given by
\f[
{\tt T} = h \ (\sigma_0 + \gamma) \ \left|\frac{d{\bf R}_w}{d\xi}\right| 
= h \ (\sigma_0 + \gamma) \ \sqrt{A}.
\f]
In this expression we have used the fact that the 2nd Piola Kirchhoff 
stress \f$ \sigma \f$ decomposes the (physical) stress vector into the
vector \f$ d{\bf R}_w/d\xi \f$. This vector is tangent to the 
deformed centreline but not necessarily a unit vector.
 
The pressure required to deform the beam into a circular arc 
with opening angle \f$ \alpha \f$ is therefore given by
\f[ 
p_{ext} =
 \frac{\alpha h}{L}
 \ \left(\sigma_0 + 
\frac{1}{2}\left(\frac{1}{4}\frac{\alpha^2}{\sin^2(\alpha/2)} - 1
\right) 
\right).
\f]


To facilitate comparisons with the numerical solutions, we determine
the opening angle \f$ \alpha \f$ as a function of the vertical 
displacement \f$ d \f$ of the beam's midpoint by determining the
auxiliary angle \f$ \beta \f$ from 
 \f[
\tan\beta = \frac{2d}{L}.
\f]
The opening angle then follows from
\f[
\tan\left(\frac{\alpha}{2}\right) = \tan(\pi -
 2\beta) = - \tan 2\beta = \frac{2\tan\beta}{\tan^{2}\beta -1}.
\f]

Here is a comparison between the computed and analytically predicted
values for the pressure \f$ p_{ext} \f$, required to deflect
the beam's midpoint by the specified displacement, \f$ d \f$.

 
\image html trace.gif "Comparison between the computed and analytical solutions. " 
\image latex trace.eps "Comparison between the computed and analytical solutions. " width=0.75\textwidth

<TABLE>
<TR>
<TD bgcolor="cornsilk">
<CENTER><B>Comment: Geometric and constitutive nonlinearities</B></CENTER>
The comparison between the computational results and the analytical 
predictions is very satisfying but is important to realise that
the agreement only validates the numerical solution, not the
physical model. \c oomph-lib's \c HermiteBeamElement is based on a \b geometrically 
\b  nonlinear theory, implying that the kinematics of the deformation
are captured exactly for arbitrarily large displacements and
rotations. However, the use of \b incrementally \b linear \b constitutive 
\b equations in the variational principle
(1) can only be justified if the strain is small. This
is clearly not the case for the deformations shown above. The
combination of geometric nonlinearity with linear constitutive
equations can, however, be justified in applications in which the beam 
undergoes large displacements with little extension of its centreline.
This typically occurs in stability problems, such as the
<a href="../../steady_ring/html/index.html">
buckling of a circular ring under external pressure</a>, considered
in <a href="../../steady_ring/html/index.html">another example</a>.





</TD> 
</TR>
</TABLE>




<HR>   
<HR>

\section global Global parameters and functions
The namespace \c Global_Physical_Variables contains the dimensionless beam
thickness, \f$ h \f$, the dimensionless pre-stress,
\f$ \sigma_{0} \f$, and the pressure load, \f$ p_{ext} \f$, as well as
the function \c Global_Physical_Variables::load() 
which computes the load vector in the form required by 
the \c HermiteBeamElements. (The function's arguments allow the load vector
to be a function of the Lagrangian and Eulerian coordinates, and
the unit normal to the deformed beam, \f$ {\bf N} \f$.)

 

\dontinclude tensioned_string.cc
\skipline start_of_namespace
\until end of namespace

<HR>
<HR>

\section main The driver code

The main code is very short. The physical parameters
\f$ h \f$, \f$ \sigma_{0} \f$ and \f$ L \f$ are initialised and the
problem is constructed using 10 elements. Following the usual
self-test, we call the function ElasticBeamProblem::parameter_study()
to compute the deformation of the beam for a range of external
pressures.


\dontinclude tensioned_string.cc
\skipline start_of_main
\until end of main

<HR>
<HR>

\section problem The problem class

The problem class has five member functions, only two of which are
non-trivial:
 - the problem constructor, \c ElasticBeamProblem(...), whose
   arguments specify the number of elements and the beam's
   undeformed length.
 - the function \c parameter_study(), which computes the beam's
   deformation for a range of external pressures.
.

 In the present problem, the functions 
\c Problem::actions_before_newton_solve() and \c
Problem::actions_after_newton_solve() are not required, so remain
empty. The function 
\c ElasticBeamProblem::mesh_pt() overloads the
(virtual) function \c Problem::mesh_pt() to return a pointer
to the specific mesh used in this problem. This avoids explicit re-casts
when member functions of the specific mesh need to be accessed.

 The class also includes three private data members which store
a pointer to a node at which the displacement is documented,
the length of the domain, and a pointer to the geometric object
that specifies the beam's undeformed shape.

\dontinclude tensioned_string.cc
\skipline start_of_problem_class
\until end of problem class
 
<HR>
<HR>

\section constructor The Problem constructor
We start by creating the undeformed centreline of the beam as a
\c StraightLine, one of \c oomph-lib's standard geometric objects. 

\skipline start_of_constructor
\until Undef_beam_pt

 We then construct the a one-dimensional Lagrangian mesh in
 two-dimensional space, using the previously-constructed geometric
 object to set the initial positions of the nodes. 

 \skipline Create
 \until OneDLagrangianMesh

The \c OneDLagrangianMesh is a \c SolidMesh whose constituent nodes are
\c SolidNodes. These nodes store not only their (variable) 2D Eulerian 
position, accessible via \c SolidNode::x(...), but also their (fixed) 
1D Lagrangian coordinates, accessible via \c SolidNode::xi(...).
The  \c OneDLagrangianMesh constructor assigns the nodes' Lagrangian
coordinate, \f$ \xi \f$, by spacing them evenly in the
range \f$ \xi \in [0,L] \f$. The \c GeomObject pointed to by
\c Undef_beam_pt, provides a parametrisation of the beam's
 undeformed shape in the form \f$ {\bf r}(\xi) \f$, and this is used
to determine the nodes's initial Eulerian position.



 Next we pin the nodal positions on both boundaries

\skipline Set the 
\until }

 We then loop over the elements and set the pointers to the
physical parameters (the pre-stress and the thickness), the 
function pointer to the load vector, and the pointer to the 
geometric object that specifies the undeformed beam shape.
 \skipline Loop
 \until end of loop over elements

 We choose a node near the centre of the beam to monitor the
displacements. (If the total number of nodes is even, the control
node will not be located at the beam's exact centre; its vertical 
displacement will therefore differ from the analytical solution
that we output in \c doc_solution(...) -- in this case we issue a 
suitable warning.)
 \skipline Choose
 \until n_nod+1

 Finally, we assign the equation numbers
 \skipline Assign
 \until end of


<HR>
<HR>

\section param The Parameter Study

 The function ElasticBeamProblem::parameter_study() is used to perform a
parameter study, computing the beam's deformation for a range
of external pressures. During the solution of this
particular problem, the maximum residual in the Newton iteration
can be greater than the default maximum value of 10.0. 
We increase the default value by assigning a (much) larger value to
\c Problem::Max_residuals.

Next, we choose the increment in the control parameter (the external pressure),
set its initial value and open an output file that will contain
the value of the external pressure, the mid-point displacement and
external pressure computed from the analytical solution. 
We also create an output stream and a string 
that will be used to write the complete solution for each value of the 
external pressure.

 In the loop, we increment the external pressure \f$ P_{ext} \f$, 
solve the problem, calculate the analytical prediction for the
pressure that is required to achieve the computed deformation, plot 
the solution and write the pressure, the displacement and exact 
pressure to the trace file.

 \skipline start_of_parameter
 \until end of parameter
  

<HR>
<HR>


\section comments Exercises and Comments

-# Modify the code so that one end of the beam is no longer fixed in
   space. What happens? Why? 
-# Increase the bending effects by increasing the beam's
   thickness to \f$ h=0.1 \f$, say, and by "clamping" its ends,
   so that \f$ d({\bf R}_w(\xi)\cdot {\bf e}_y)/d\xi = 0 \f$ at \f$ \xi=0 \f$
   and \f$ \xi=L \f$. This condition can be enforced by pinning the
   "type 1" (slope) positional degree of freedom in the vertical
   (1) direction at both ends; this requires the insertion of the statement
   \code
   mesh_pt()->boundary_node_pt(b,0)->pin_position(1,1); 
   \endcode
   in the loop in the \c Problem constructor.
   You will have to adjust the number of elements
   to fully resolve the bending boundary layers. 
-# \c oomph-lib's default \anchor newton_cust nonlinear solver, \c
   Problem::newton_solve(...), provides an implementation of
   the Newton method, which has the attractive feature
   that it converges quadratically -- provided 
   a good initial guess for the solution is available. Good initial guesses 
   can often (usually?) be generated by computing the solution via a 
   sequence of substeps, as in the above example where we started with
   a known solution (the undeformed beam -- the exact solution 
   for \f$ p_{ext} = p_{ext}^{(0)} =0 \f$ ) and used it as the
   initial guess for the solution at a small external
   pressure, \f$ p_{ext}^{(1)} \f$. When the Newton method converged,
   we used the computed
   solution as the initial guess for the solution at a
   slightly larger pressure, \f$ p_{ext}^{(2)} \f$, etc. 
   If the increase in the load (or some other control parameter) is
   too large, the Newton method will diverge. To avoid unnecessary
   computations, the Newton iteration is terminated if:
   - the number of iterations exceeds \c
     Problem::Max_newton_iterations (which has a default value of 10)
   - the residual exceeds \c Problem::Max_residuals (which has a 
     default value of 10.0)
   .
   The Newton method continues until the maximum residual
   has been reduced to  \c Problem::Newton_solver_tolerance, which has
   a default value of \f$ 10^{-8} \f$.
   \n\n
   All three values are protected data members of the \c Problem base class
   and can therefore be changed in any specific \c Problem. For
   instance, in the problem considered above, the undeformed beam
   provides a poor approximation of its equilibrium shape at the 
   first pressure value.
   The Newton method still converges (very slowly initially, then
   quadratically as it approaches the exact solution), even though
   the initial maximum residual has a relatively large value of 
   19.6. Here are some exercises that explore the convergence
   characteristics of the Newton method:
   -# Experiment with the Newton solver and find the largest
      value for the load increment, \c pext_increment, for which
      the Newton method still converges.
   -# Explain why the Newton method converges 
      very slowly for small values of \f$ p_{ext} \f$ and much more rapidly
      at larger values, even though the load increment, \c
      pext_increment, remains constant.
   -# Compare the solutions for different values of 
      \c Problem::Newton_solver_tolerance. Is the default value  
      \f$ 10^{-8} \f$ adequate? Reduce it to  \f$ 10^{-12} \f$ and
      \f$ 10^{-18} \f$. What do you observe?  



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="
../../../../
demo_drivers/beam/tensioned_string/
">
demo_drivers/beam/tensioned_string/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="
../../../../
demo_drivers/beam/tensioned_string/tensioned_string.cc
">
demo_drivers/beam/tensioned_string/tensioned_string.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Demo problem: A preconditioner for the solution of Navier-Stokes equations with weakly imposed boundary conditions via Lagrange multipliers


The purpose of this tutorial is to show how to use \c oomph-lib's  Lagrange
Enforced Flow Navier-Stokes preconditioner. Similarly to the problem considered
in the <a href="../../../navier_stokes/vmtk_fluid/html/index.html"> Steady
finite-Reynolds-number flow through an iliac bifurcation</a> tutorial, the
outflow boundary of the demo problem (discussed below) are not aligned with any
coordinate planes. Parallel outflow is therefore enforced by a Lagrange
multiplier method, implemented using \c oomph-lib's \c FaceElement framework.

<HR>
<HR>

\section model The model problem, theory and preconditioner

We will demonstrate the development and application of the preconditioner using
the Poiseuille flow through a unit square domain \f$\Omega^{[\alpha]} \in
\mathbb{R}^2\f$ rotated by an arbitrary angle \f$\alpha\f$ (see the figure
below). The domain \f$\Omega^{[\alpha]}\f$ is obtained by rotating the discrete
points \f$(x_1,x_2)\f$ in the unit square \f$\Omega = [0,1]^2\f$ by the
following transformation
\f[
R(\alpha)=
\left( 
\begin{array}{cc}
{\cos(\alpha)} & {-\sin(\alpha)} \\ {\sin(\alpha)} & {\cos(\alpha)} 
\end{array} 
\right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
where \f$\alpha\f$ is the angle of rotation.  The figure below show the flow
field (velocity vectors and pressure contours) for a unit square domain rotated
by an angle of \f$\alpha=30^\circ \f$ and a Reynolds number of \f$ Re = 100
\f$. The flow is driven by a prescribed parabolic boundary condition.  

\image html SimAng30Rey100Noel4.gif "Velocity field and pressure " 
\image latex SimAng30Rey100Noel4.eps "Velocity field and pressure " width=0.4\textwidth

For convenience, we present the boundary conditions for the non-rotated unit
square \f$\alpha = 0^{\circ}\f$. In order to obtain the boundary conditions for
\f$\alpha \ne 0^\circ\f$, we only have to apply the rotation
(1).  The flow is driven by imposing a parabolic velocity
profile along the inflow boundary \f$\Omega_I\f$. Along the characteristic
boundary, \f$\Omega_C\f$, the no-slip condition \f$u_i = 0\f$, \f$i = 1,2\f$,
is prescribed.  We impose `parallel outflow' along the outlet \f$\Omega_O\f$ by
insisting that
\f[     
{\bf u}\cdot {\bf t}_{} = 0 \mbox{\ \ \ \ on $\Omega_{O}$,}
 \ \ \ \ \ \ \ \ \ \ (2)
\f]
where \f$ {\bf t}\f$ is the tangent vector at each discrete point on the
boundary \f$ \Omega_{O}\f$.  We weakly enforce the flow constraint by
augmenting the Navier-Stokes momentum residual equation (introduced in
the
<a href="../../../navier_stokes/rayleigh_traction_channel/html/index.html">
Unsteady flow in a 2D channel, driven by an applied traction</a> tutorial) with
a Lagrange multiplier term so that it becomes
\f[
r_{il}^{u}
=
\int_\Omega
\left[
  Re
    \left(
      St\frac{\partial u_i}{\partial t}
      + u_j\frac{\partial u_i}{\partial x_j}
    \right) \psi_l
  + \tau_{ij}\frac{\partial \psi_l}{\partial x_j}
\right]
d\Omega
-
\int_{\partial\Omega}
\tau_{ij} n_j\psi_l dS
+ \delta \Pi_{constraint} = 0,
 \ \ \ \ \ \ \ \ \ \ (3)
\f]
where
\f[
\Pi_{constraint}
=
\int_{\partial \Omega} \lambda u_i t_i dS,
 \ \ \ \ \ \ \ \ \ \ (4)
\f]
and \f$\lambda\f$ is the Lagrange multiplier. Upon taking the first variation
of the constraint with respect to the unknown velocity and the Lagrange
multiplier, the residual form of the constrained momentum equation is
\f[
r_{il}^{u}
=
\int_\Omega
\left[
  Re
    \left(
      St\frac{\partial u_i}{\partial t}
      + u_j\frac{\partial u_i}{\partial x_j}
    \right) \psi_l
  + \tau_{ij}\frac{\partial \psi_l}{\partial x_j}
\right]
d\Omega
-
\int_{\partial\Omega}
\tau_{ij} n_j\psi_l dS
+ 
\int_{\partial\Omega}
\lambda\psi_l t_i = 0.
 \ \ \ \ \ \ \ \ \ \ (5)
\f]
The weak formulation of (2) is simply
\f[
r_{l}^{\lambda}
=
\int_\Omega
u_i t_i \psi^\lambda
dS
=
0,
 \ \ \ \ \ \ \ \ \ \ (6)
\f]
where \f$\psi^\lambda\f$ is a suitable basis function.  Equation
(5) reveals that the Lagrange multipliers act as
the (negative) tangential traction \f$(\lambda =
-{\mathbf{n}}^T{\mathbf{\tau}}{\mathbf{t}})\f$ that enforce the parallel flow
across the boundary \f$\partial\Omega_{O}\f$.  We discretise this constraint by
attaching \c ImposeParallelOutflowElements to the boundaries of the "bulk"
Navier-Stokes elements that are adjacent to \f$ \partial \Omega_O \f$ as shown
in the <a href="../../../navier_stokes/vmtk_fluid/html/index.html"> Steady
finite-Reynolds-number flow through an iliac bifurcation</a> tutorial, also see
the <a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
Deformation of a solid by a prescribed boundary motion </a> tutorial which
employs a similar technique used to enforce prescribed boundary displacements
in solid mechanics problems. We discretise the Navier-Stokes equations using \c
oomph-lib's \c QTaylorHoodElements, see the <a href="../../../navier_stokes/driven_cavity/html/index.html"> 2D Driven Cavity
Problem</a> tutorial for more information.


The discretised problem therefore contains the following types of discrete 
unknowns:
- The fluid degrees of freedom (velocity and pressure).\n\n
- The nodal values representing the components of the (vector-valued)
  Lagrange multipliers. These only exist for the nodes on 
  \f$ \partial \Omega_O\f$. (The nodes are re-sized to
  accommodate the additional unknowns when the \c ImposeParallelOutflowElements are
  attached to the bulk elements.)
. 


The preconditioner requires a further sub-division of these degrees of
freedom into the following categories:
- the unconstrained velocity in the x-direction
- the unconstrained velocity in the y-direction
- [the unconstrained velocity in the z-direction (only in 3D)]
- the constrained velocity in the x-direction
- the constrained velocity in the y-direction
- [the constrained velocity in the z-direction (only in 3D)]
- the Lagrange multiplier at the constrained nodes
- [the other Lagrange multiplier at the constrained nodes
  (only in 3D)].
.
For a 2D problem, the linear system to be solved in the course of the Newton
iteration can then be (formally) re-ordered into the following block structure:

\f[
\left[
\begin{array}{ccccc|c}
{ F_{\rm xx}}&{ F_{\rm x\bar{\rm
      x}}}&{ F_{\rm xy}}&{ F_{\rm x\bar{\rm y}}}&{ B_{\rm x}^{T}}&\\
{ F_{\bar{\rm x}\rm x}}&{ F_{\bar{\rm x}\bar{\rm x}}}&{ F_{\bar{\rm x}\rm y}}&{ F_{\bar{\rm x}\bar{\rm y}}}&{ B_{\bar{\rm x}}^{T}}&{ M_{\rm x}}\\
{ F_{\rm yx}}&{ F_{\rm y\bar{\rm
      x}}}&{ F_{\rm yy}}&{ F_{\rm y\bar{\rm y}}}&{ B_{\rm y}^{T}}&\\
{ F_{\bar{\rm y}\rm x}}&{ F_{\bar{\rm
      y}\bar{\rm x}}}&{ F_{\bar{\rm y}\rm y}}&{ 
  F_{\bar{\rm y}\bar{\rm y}}}&{ B_{\bar{\rm y}}^{T}}&{ M_{\rm y}}\\
{ B_{\rm x}}&{ B_{\bar{\rm x}}}&{ B_{\rm y}}&{ B_{\bar{\rm y}}}&\\
\hline
&{ M_{\rm x}}&&{ M_{\rm y}}&&
\end{array}
\right]
\left[
\begin{array}{c}
\Delta \mathbf{U}_{\rm x}\\
\Delta \mathbf{\overline{U}}_{\rm x}\\
\Delta \mathbf{U}_{\rm y}\\
\Delta \mathbf{\overline{U}}_{\rm y}\\
\Delta \mathbf{P}\\
\Delta \mathbf{\Lambda}
\end{array}
\right]
 =
-
\left[
\begin{array}{c}
\mathbf{r}_{\rm x}\\
\mathbf{r}_{\bar{\rm x}}\\
\mathbf{r}_{\rm y}\\
\mathbf{r}_{\bar{\rm y}}\\
\mathbf{r}_{\rm p}\\
\mathbf{r}_{\rm \Lambda}
\end{array}
\right].
\ \ \ \ \ \   (7)
\f]
Here the vectors \f${\bf U}_{\rm x}\f$, \f${\bf U}_{\rm y}\f$, \f${\bf P}\f$
and \f${\bf \Lambda}\f$ contain the \f$x\f$ and \f$y\f$ components of the
velocity unknowns, the pressure unknowns and Lagrange multipliers unknowns,
respectively. The overbars identify the unknown nodal positions that are
constrained by the Lagrange multiplier. The matrices \f$M_{\rm x}\f$ and
\f$M_{\rm y}\f$ are mass-like matrices whose entries are formed from products
of the basis functions multiplied by a component of the tangent vector at each
discrete point on \f$\partial\Omega_{O}\f$, for example,
\f$[M_{{\rm x}}]_{ij}
=
\int_{\partial \Omega_{O}} t_x {\psi_i} \ {\psi_j} \ dS \f$. Denote
\f[
J_{\rm NS}
=
\left[
\begin{array}{ccccc}
{ F_{\rm xx}}&{ F_{\rm x\bar{\rm
      x}}}&{ F_{\rm xy}}&{ F_{\rm x\bar{\rm y}}}&{ B_{\rm x}^{T}}\\
{ F_{\bar{\rm x}\rm x}}&{ F_{\bar{\rm x}\bar{\rm x}}}&{ F_{\bar{\rm x}\rm y}}&{ F_{\bar{\rm x}\bar{\rm y}}}&{ B_{\bar{\rm x}}^{T}}\\
{ F_{\rm yx}}&{ F_{\rm y\bar{\rm
      x}}}&{ F_{\rm yy}}&{ F_{\rm y\bar{\rm y}}}&{ B_{\rm y}^{T}}\\
{ F_{\bar{\rm y}\rm x}}&{ F_{\bar{\rm
      y}\bar{\rm x}}}&{ F_{\bar{\rm y}\rm y}}&{ 
  F_{\bar{\rm y}\bar{\rm y}}}&{ B_{\bar{\rm y}}^{T}}\\
{ B_{\rm x}}&{ B_{\bar{\rm x}}}&{ B_{\rm y}}&{ B_{\bar{\rm y}}}
\end{array}
\right],
 \ \ \ 
L = 
\left[
\begin{array}{ccccc}
 &{M_{\rm x}}& &{ M_{\rm x}} &
\end{array}
\right],
 \ \ \ 
\Delta \mathbf{X}_{\rm NS}
=
\left[
\begin{array}{c}
\Delta \mathbf{U}_{\rm x}\\
\Delta \mathbf{\overline{U}}_{\rm x}\\
\Delta \mathbf{U}_{\rm y}\\
\Delta \mathbf{\overline{U}}_{\rm y}\\
\Delta \mathbf{P}
\end{array}
\right],\mbox{\ \ \ \ and }
 \ \ \ 
\mathbf{r}_{\rm NS}
=
\left[
\begin{array}{c}
\mathbf{r}_{\rm x}\\
\mathbf{r}_{\bar{\rm x}}\\
\mathbf{r}_{\rm y}\\
\mathbf{r}_{\bar{\rm y}}\\
\mathbf{r}_{\rm p}
\end{array}
\right].
\f]
Then we can re-write (7) as
\f[
\left[
\begin{array}{cc}
{ J_{\rm NS}}&L^{T}\\
L& 
\end{array}
\right]
\left[
\begin{array}{c}
\Delta \mathbf{X}_{\rm NS}\\
\Delta \mathbf{\Lambda}
\end{array}
\right]
 =
-
\left[
\begin{array}{c}
\mathbf{r}_{\rm NS}\\
\mathbf{r}_{\rm \Lambda}
\end{array}
\right].
\ \ \ \ \ \   (8)
\f]
We have shown that 
\f[
P=
\left[
\begin{array}{cc}
{ J_{\rm NS}} + L^{T}W^{-1} L& \\
 & W
\end{array}
\right],
\ \ \ \ \ \   (9)
\f]
where \f$W = \frac{1}{\sigma}LL^{T}\f$ is an optimal preconditioner for the
linear system (8) if we set \f$ \sigma=\|F\|_{\infty}
\f$ where \f$F\f$ is the compound \f$4\times 4\f$ top-left block 
\f[
F
=
\left[
\begin{array}{ccccc}
{ F_{\rm xx}}&{ F_{\rm x\bar{\rm
      x}}}&{ F_{\rm xy}}&{ F_{\rm x\bar{\rm y}}}\\
{ F_{\bar{\rm x}\rm x}}&{ F_{\bar{\rm x}\bar{\rm x}}}&{ F_{\bar{\rm x}\rm y}}&{ F_{\bar{\rm x}\bar{\rm y}}}\\
{ F_{\rm yx}}&{ F_{\rm y\bar{\rm
      x}}}&{ F_{\rm yy}}&{ F_{\rm y\bar{\rm y}}}\\
{ F_{\bar{\rm y}\rm x}}&{ F_{\bar{\rm
      y}\bar{\rm x}}}&{ F_{\bar{\rm y}\rm y}}&{ 
  F_{\bar{\rm y}\bar{\rm y}}}
\end{array}
\right]
\f]
in the Jacobian matrix.  Application of the preconditioner \f$P\f$ requires the
repeated solution of linear systems involving the diagonal blocks \f${J_{\rm
NS}} + L^{T}W^{-1} L\f$ and \f$W\f$.  The matrix \f$W^{-1} = \sigma
(LL^{T})^{-1} = \sigma ({M_{\rm x}}^{2}+{M_{\rm y}^{2}})^{-1}\f$ is dense and
will cause the addition of dense sub-matrices to the Jacobian matrix:
\f[
\sigma L^{T}(LL^{T})^{-1}L
=
\sigma
\left[
\begin{array}{ccccc}
\cal O &\cal O &\cal O &\cal O &\cal O \\
 \cal O&{M_{\rm x}}({M_{\rm x}}^{2}+{M_{\rm y}^{2}})^{-1}{M_{\rm x}}&\cal O &{M_{\rm x}}({M_{\rm x}}^{2}+{M_{\rm y}^{2}})^{-1}{M_{\rm y}}&\cal O \\
 \cal O&\cal O &\cal O &\cal O &\cal O \\
 \cal O&{M_{\rm y}}({M_{\rm x}}^{2}+{M_{\rm y}^{2}})^{-1}{M_{\rm x}}&\cal O &{M_{\rm y}}({M_{\rm x}}^{2}+{M_{\rm y}^{2}})^{-1}{M_{\rm y}}& \cal O\\
 \cal O&\cal O &\cal O &\cal O &\cal O
\end{array}
\right].
\f]
Numerical experiments show that an efficient implementation can be obtained by
replacing \f$W\f$ by its diagonal approximation \f$\widehat{W} =
\mbox{diag}(M_{\rm x}^{2}+M_{\rm y}^{2})\f$.  Then the inversion of \f$W\f$ is
straight forward and the addition of \f$L^{T}{\widehat{W}}^{-1}L\f$ to the
Jacobian does not significantly increase the number of non-zero entries in the
matrix \f$J_{\rm NS}\f$.  Denote the efficient implementation by
\f[
\tilde{P}
=
\left[
\begin{array}{cc}
{ \tilde{J}_{\rm NS}}& \\
 & \widehat{W}
\end{array}
\right],
\ \ \ \ \ \   (10)
\f]
where \f$\tilde{J}_{\rm NS} =J_{\rm NS}+ L^{T}{\widehat{W}}^{-1}L\f$ is the
augmented Navier-Stokes Jacobian matrix.  In our implementation of the
preconditioner, the linear system involving \f$\tilde{J}_{\rm NS}\f$ can either
be solved "exactly", using \c SuperLU (in its incarnation as an exact
preconditioner; this is the default) or by any other \c Preconditioner
(interpreted as an "inexact solver") specified via the access function
\code
LagrangeEnforcedFlowPreconditioner::set_navier_stokes_preconditioner(...)
\endcode
Numerical experiments show that a nearly optimal preconditioner is obtained by
replacing the solution of the linear system involving the augmented
Navier-Stokes Jacobian \f$ \tilde{J}_{\rm NS} \f$ by an application of Elman,
Silvester and Wathen's <a href="../../lsc_navier_stokes/html/index.html">
Least-Squares Commutator (LSC) preconditioner</a>, and by replacing the
remaining block-solves within these preconditioners by a small number of AMG
cycles.

With these approximations, the computational cost of one application of
\f$\tilde P\f$ is linear in the number of unknowns.  The optimality of the
preconditioner can therefore be assessed by demonstrating that the number of
GMRES iterations remains constant under mesh refinement.

<HR>
<HR>

\section example Demo driver and use of the preconditioner



To demonstrate how to use the preconditioner, here are the relevant extracts
from the driver code <A HREF="../../../../demo_drivers/navier_stokes/lagrange_enforced_flow_preconditioner/two_d_tilted_square.cc">
two_d_tilted_square.cc</A> which solves the model problem described above.  As
explained in the <a href="../../../linear_solvers/html/index.html">Linear
Solvers Tutorial</A> switching to an iterative linear solver is typically
performed in the \c Problem constructor and involves a few straightforward
steps:

-# <b>Create an instance of the IterativeLinearSolver and pass it to the
   Problem</b> \n\n
   In our problem, we choose right preconditioned \c GMRES as the linear 
   solver: \n\n
   \dontinclude two_d_tilted_square.cc
   \skipline Create oomph-lib
   \until solver_pt;
   \n\n
-# <b>Create an instance of the Preconditioner and give it access to the 
   meshes</b> \n\n
   The \c LagrangeEnforceFlowPreconditioner takes a pointer of meshes. 
   It is important that the bulk mesh is in position \c 0 :
   \skipline Create the preconditioner
   \until set_meshes(mesh_pt);
   By default, \c SuperLUPreconditioner is used for all subsidiary block
   solves. To use the LSC preconditioner to approximately solve the
   sub-block system involving the momentum block, we invoke the following:
   \skipline Create the NS LSC preconditioner
   \until set_navier_stokes_preconditioner(lsc_prec_pt);
   The LSC preconditioner is discussed in
   <a href="../../lsc_navier_stokes/html/index.html">
   another tutorial</A>.
   \n\n
-# <b>Pass the preconditioner to the solver, and the solver to the problem
   </b> \n\n
   \skipline Pass the preconditioner to the solver.
   \until linear_solver_pt() = Solver_pt;
   \n\n

.

<HR>
<HR> 






\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/lagrange_enforced_flow_preconditioner">
demo_drivers/navier_stokes/lagrange_enforced_flow_preconditioner
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/lagrange_enforced_flow_preconditioner/two_d_tilted_square.cc">
demo_drivers/navier_stokes/lagrange_enforced_flow_preconditioner/two_d_tilted_square.cc
</A>
</CENTER>







<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


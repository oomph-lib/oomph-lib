/**

\mainpage Demo problem: A preconditioner for the solution of fluid-structure interaction problems with (pseudo-)solid fluid mesh updates

In this tutorial we re-visit the solution of fluid-structure 
interaction problems with (pseudo-)solid fluid mesh updates,
focusing on the efficient solution of the governing equations
by GMRES, using the problem-specific block-preconditioner developed
in 
<center>
Muddle, R.L., Mihajlovic, M. & Heil, M. (2012)
An efficient preconditioner for monolithically-coupled 
large-displacement fluid-structure interaction problems 
with pseudo-solid mesh updates. <em>Journal of Computational Physics</em>
<b>231</b>, 7315-7334. DOI: 
<a href="http://dx.doi.org/10.1016/j.jcp.2012.07.001">
10.1016/j.jcp.2012.07.001</a>
</center>
where full details of the analysis and further results from a 
variety of other test problems can be found.

 
<HR> 
<HR>

\section model Background: Problem formulation, the resulting linear system and the preconditioner

We will demonstrate the development and application 
of the preconditioner using the
problem of flow in a 2D channel with an elastic leaflet discussed in
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
another tutorial</a> and illustrated in the sketch below:

\image html fsi_channel_with_leaflet.gif "Sketch of the problem. " 
\image latex fsi_channel_with_leaflet.eps "Sketch of the problem. " width=0.75\textwidth

Flow is driven (via an imposed parabolic inflow profile) through a 2D channel 
which is partially obstructed by a thin-walled elastic leaflet.
The fluid leaves the domain at the far downstream end of the channel
where we assume parallel, axially traction-free outflow. 

We use the same discretisation (2D quadrilateral
Taylor-Hood elements for the fluid; 1D Hermite beam elements for the
leaflet) as in the 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
original tutorial</a>. However, here we perform the fluid mesh update 
in response to the deformation of the elastic leaflet using a pseudo-elastic
approach. This is done by "wrapping" the underlying  \c
QTaylorHoodElement into a \c PseudoSolidNodeUpdateElement
and employing \c ImposeDisplacementByLagrangeMultiplierElements 
to impose the deformation of the elastic leaflet onto the (pseudo-solid)
fluid mesh. The following tutorials discuss the
relevant methodologies:
- <a href="../../../solid/prescribed_displ_lagr_mult/html/index.html">
  The use of Lagrange multipliers to impose displacement constraints
  onto (pseudo-)solids.</a>\n\n
- <a href="../../../interaction/unstructured_fsi/html/index.html">
  The use of pseudo-elasticity to deform the fluid mesh in FSI
  problems.</a> \n\n
.

The discretised problem contains the 
following types of discrete unknowns:
- The fluid degrees of freedom (velocities and pressures).\n\n
- The nodal positions (real and generalised; the latter representing
  derivative degrees of freedom in the Hermite elements used to
  discretise the leaflet; see the 
  <a href="../../../beam/tensioned_string/html/index.html">
  beam theory tutorial</a> for details) of the elastic leaflet.\n\n
- The nodal positions in the fluid mesh. \n\n
- The nodal values representing the components of the (vector-valued)
  Lagrange multipliers which impose the motion of the elastic leaflet
  on the pseudo-elastic fluid mesh. These unknowns only exist for the
  fluid nodes on the FSI boundary.
.
Using this classification of the unknowns,
the linear system to be solved in the course 
of the Newton iteration can be re-ordered 
into the following block structure:
\f[
J_{\rm FSI} \ \Delta {\bf x} = 
 \left[
\begin{array}{c|c|cc}
F& &C_{\rm fx}&\\ \hline C_{\rm sf}&S&C_{\rm sx}&\\ \hline
&&E&C_{\rm xl}\\ &C_{\rm ls}&C_{\rm lx}&\\
\end{array}
\right] \left[
\begin{array}{c}
 \Delta {\bf F}\\ 
\hline \Delta {\bf S}\\ 
\hline \Delta {\bf X}\\ 
\Delta {\bf L}\\
\end{array}
\right] =- \left[
\begin{array}{c}
{\bf r}_{\rm f}\\ 
\hline {\bf r}_{\rm s}\\ 
\hline {\bf r}_{\rm x}\\ 
{\bf r}_{\rm l}\\
\end{array}
\right] 
\ \ \ \ \ \ \ \ (1)
\f]
where the vector \f${\bf F}\f$ contains the Navier-Stokes (fluid) unknowns 
(velocities and pressure), \f${\bf S}\f$ represents unknowns describing the
deformation of the fluid-loaded solid, \f${\bf X}\f$ represents 
the nodal positions in the fluid mesh, and \f${\bf L}\f$ contains 
the discrete Lagrange multipliers which impose
the deformation of the FSI boundary in the fluid mesh.
The three diagonal blocks in the Jacobian matrix, \f$J_{\rm FSI}\f$, are 
the two ``single-physics'' 
Jacobian matrices (\f$F\f$, the Navier-Stokes Jacobian; \f$S\f$, the 
tangent stiffness matrix of the fluid-loaded solid) and 
the Jacobian, 
\f[
J_{\rm PS} = 
\left[
\begin{array}{cc}
E&C_{\rm xl}\\ 
C_{\rm lx}&\\
\end{array}
\right]
\f]
associated with the Lagrange-multiplier-constrained pseudo-elasticity
problem governing the fluid mesh update, discussed in more detail
in <a href="../../prescribed_displ_lagr_mult/html/index.html">
another tutorial</a>.
The non-zero off-diagonal blocks arise through the interactions between 
the various single-physics
problems: \f$C_{\rm fx}\f$ represents the effect of the pseudo-solid unknowns 
(the nodal positions in the fluid mesh) on the discretised 
Navier--Stokes equations -- this incorporates the so-called shape-derivatives
and the time-discretised version of the no-slip condition on the FSI
boundary which expresses the fluid velocities in 
terms of the nodal velocities of the fluid mesh. 
\f$C_{\rm sf}\f$ captures the effect of the fluid traction (pressure and shear 
stresses) on the fluid-loaded solid. Since the shear stresses depend
on gradients of the fluid velocity, the traction 
is also affected by the nodal positions in the fluid mesh -- this
dependency gives rise to the matrix \f$C_{\rm sx}\f$. Finally, \f$C_{\rm
ls}\f$ arises because the prescribed boundary displacement of the fluid 
mesh depends
on the displacement field in the fluid-loaded solid. The zero \f$ "fs"
\f$ block
reflects the fact that the discretised Navier-Stokes equations do not
depend directly on the displacement field of the actual solid. The
zero \f$ "xf" \f$ block indicates that the 
pseudo-solid equations are not affected by the fluid unknowns. The 
\f$ "xs"\f$ block is zero
because the real solid affects the pseudo-solid only indirectly
via the Lagrange multiplier constraint (which gives rise to  \f$C_{\rm
ls}\f$). Finally, the non-zero blocks in the last block column
reflect the fact that the Lagrange multiplier is only used to 
enforce the displacement constraint 
in a one-way coupling in which the real solid affects the pseudo-solid 
but not vice versa.
See <a href="http://dx.doi.org/10.1016/j.jcp.2012.07.001">
Muddle, Mihajlovic & Heil (2012) </a> for more detail.

To construct a preconditioner for the linear system
(1) we replace the
bottom-right \f$2 \times 2\f$ block in the Jacobian with the pseudo-solid
preconditioner, 
\f[
{\cal P}_{\rm PS} = \left[
\begin{array}{cc}
E_{\rm aug}&\\
&W
\end{array}
\right]
\ \ \ \ \ \ \ \ \ \ \ \ \ \ (2)
\f] 
which we discussed in
<a href="../../prescribed_displ_lagr_mult/html/index.html">
another tutorial</a>. This yields the pseudo-solid FSI preconditioner
\f[
{\cal P}_{\rm FSI}=\left[
\begin{array}{c|c|cc}
F& &C_{\rm fx}&\\ \hline
C_{\rm sf}&S&C_{\rm sx}&\\ \hline
&&E_{\rm aug}&\\
&C_{\rm ls}&&W
\end{array}
\right].
\ \ \ \ \ \ (3)
\f]

The preconditioner can be shown to have a block triangular structure.
Its  application therefore
requires the solution  of four subsidiary linear systems involving the
matrices
\f$F\f$, \f$S\f$, \f$E_{\rm aug}\f$ and \f$W\f$ and 
four sparse matrix-vector products with the interaction 
matrices \f$C_{\rm fx}\f$, \f$C_{\rm sx}\f$, \f$C_{\rm sf}\f$ and
\f$C_{ls}\f$. 

Numerical experiments in <a href="http://dx.doi.org/10.1016/j.jcp.2012.07.001">
Muddle, Mihajlovic & Heil (2012) </a> show that an efficient 
implementation of the preconditioner is obtained by
- replacing the solution of the linear system involving 
  \f${\cal P}_{\rm PS}\f$ by the approximation \f$\widetilde{\cal P}_{\rm PS}\f$
  discussed in <a href="../../prescribed_displ_lagr_mult/html/index.html">
  another tutorial,</a> \n\n
- replacing the solution of the linear system involving the
  Navier-Stokes Jacobian \f$ F \f$ by an application of 
  Elman, Silvester and Wathen's
  <a href="../../lsc_navier_stokes/html/index.html">
  Least-Squares Commutator (LSC)
  preconditioner</a>,
.
and by replacing the remaining block-solves within these preconditioners
by a small number of AMG cycles or CG iterations. We retain
SuperLU as the (exact) direct solver for the linear system involving
the real solid's tangent stiffness matrix \f$ S \f$.

With these approximations, the computational
cost of one application of the preconditioner is linear
in the number of unknowns (see \ref comm_and_ex
for a more detailed discussion of this issue).
The optimality of the preconditioner
can therefore be assessed by demonstrating that 
the number of GMRES iterations remains constant under 
mesh refinement.


<HR>
<HR>

\section implementation Implementation and use of the preconditioner

The preconditioner described above is implemented within \c oomph-lib's
(parallel) block preconditioning framework which is described in 
<a href="../../../mpi/block_preconditioners/html/index.html">
another tutorial</a>. For the purpose of the implementation,
we decompose the preconditioning matrix into the 3x3
main blocks indicated by the vertical and horizontal lines
in (1) and (3).

The degrees of freedom within these sub-blocks are classified
into "dof-types" on an element-by-element basis, using the functions
- \c GeneralisedElement::ndof_types(...) which returns the number of
  dof types classified by that element.
- \c GeneralisedElement::get_block_numbers_for_unknowns(...) which
  associates each degree of freedom (identified by its global equation
  number) with the dof type within its block. 
  \n\n
.
The default implementation of these functions within the
Navier-Stokes elements (which differentiate between the 
fluid velocity components and the pressure)
and the beam elements (which do not distinguish between the different types
of degree of freedom) is appropriate for the use with the
(subsidiary) preconditioners employed here. The degrees of freedom
representing the nodal positions in the (pseudo-solid) fluid mesh
need to be sub-divided into those that are and are not constrained
by the Lagrange multipliers that impose the deformation of the
FSI boundary. This is most conveniently done by overloading the relevant
functions in a templated wrapper element, \c PseudoElasticBulkElement<...>,  
as discussed in <a href="../../prescribed_displ_lagr_mult/html/index.html">
another tutorial.</a>




<HR>
<HR>

 
\section results Results
We examine the performance of the preconditioner in steady and
unsteady test problems. First we perform a sequence of steady solves,
incrementing the Reynolds number in steps of 25. The final steady
solution is then used as the initial condition for an unsteady
simulation in which the inflow is subjected to a time-periodic 
oscillation. The tables below show the GMRES iteration counts (averaged 
over all linear solves performed in the course of all Newton
iterations) as a function of the mesh refinement (represented by the
total number of unknowns) for different implementations of the
preconditioner.


<CENTER>Average GMRES iteration counts (steady runs)<TABLE BORDER=1>
<TR> <TD> <CODE>n_dof</CODE> </TD> <TD> 9570 </TD> <TD> 38724 </TD> <TD> 87462 </TD> <TD> 155784 </TD>  </TR>
<TR><TD>GMRES (blocks solved by SuperLU)</TD> <TD> 5.88889 </TD> <TD> 5.88889 </TD> <TD> 5.88889 </TD> <TD> 6.11111 </TD>  </TR>
<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD> <TD> 27.8889 </TD> <TD> 36.1111 </TD> <TD> 42.3333 </TD> <TD> 47.2222 </TD>  </TR>
<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD> <TD> 37.3333 </TD> <TD> 43.2222 </TD> <TD> 49.3333 </TD> <TD> 55.3333 </TD>  </TR>
</TABLE></CENTER>



<CENTER><BR><BR>Average GMRES iteration counts (unsteady runs)<TABLE BORDER=1>
<TR> <TD> <CODE>n_dof</CODE> </TD> <TD> 9570 </TD> <TD> 38724 </TD> <TD> 87462 </TD> <TD> 155784 </TD>  </TR>
<TR><TD>GMRES (blocks solved by SuperLU)</TD> <TD> 6.5124 </TD> <TD> 6.95868 </TD> <TD> 7.15323 </TD> <TD> 7.2619 </TD>  </TR>
<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD> <TD> 15 </TD> <TD> 15.9256 </TD> <TD> 16.2097 </TD> <TD> 16.6429 </TD>  </TR>
<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD> <TD> 22.2066 </TD> <TD> 25.3967 </TD> <TD> 27.0484 </TD> <TD> 27.9921 </TD>  </TR>
</TABLE></CENTER>




As expected from the theory (see 
<a href="http://dx.doi.org/10.1016/j.jcp.2012.07.001">
  Muddle, Mihajlovic & Heil (2012)</a>), the GMRES iteration counts 
are small and virtually mesh independent for the
exact implementation of the preconditioner. Replacing the (costly) exact 
solves by (faster) approximate solves leads to a modest
increase in the (absolute) number of GMRES iterations. The most
significant deterioration of the iteration counts results from the
use of the LSC Navier Stokes preconditioner as the inexact solver
(subsidiary preconditioner) for the Navier-Stokes block. 
The behaviour observed here (slight mesh dependence for 
steady solves; virtual mesh independence (once the mesh is
sufficiently fine) for unsteady solves) mirrors that observed
in single-physics Navier-Stokes problems.

The benefit of switching to the approximate solves becomes apparent
in the next two tables which shows the average cpu times required 
for the solution of the linear systems by GMRES.
(For comparison the times labelled as <b>"SuperLU"</b> show the times 
when the linear systems (1)
are solved directly by SuperLU rather than by GMRES; they illustrate
how quickly direct solvers become uncompetitive.) For sufficiently 
fine discretisations the larger number of GMRES iterations for the
inexact solves is more than compensated for by the much lower cost
of the preconditioning operations. The final
implementation yields cpu times that (at least for unsteady problems)
 are proportional to the number of unknowns -- the hallmark of an optimal solver.

<CENTER>Average linear solver times (steady runs; sec)<TABLE BORDER=1>
<TR> <TD> <CODE>n_dof</CODE> </TD> <TD> 9570 </TD> <TD> 38724 </TD> <TD> 87462 </TD> <TD> 155784 </TD>  </TR>
<TR><TD>SuperLU</TD> <TD> 2.71863 </TD> <TD> 33.7939 </TD> <TD> 129.248 </TD> <TD> 325.384 </TD>  </TR>
<TR><TD>GMRES (blocks solved by SuperLU)</TD> <TD> 2.31567 </TD> <TD> 15.1042 </TD> <TD> 58.5491 </TD> <TD> 147.531 </TD>  </TR>
<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD> <TD> 2.8051 </TD> <TD> 18.3197 </TD> <TD> 58.7335 </TD> <TD> 138.558 </TD>  </TR>
<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD> <TD> 3.47733 </TD> <TD> 15.1812 </TD> <TD> 37.6856 </TD> <TD> 73.0526 </TD>  </TR>
</TABLE></CENTER>

<CENTER><BR><BR>Average linear solver times (unsteady runs; sec)<TABLE BORDER=1>
<TR> <TD> <CODE>n_dof</CODE> </TD> <TD> 9570 </TD> <TD> 38724 </TD> <TD> 87462 </TD> <TD> 155784 </TD>  </TR>
<TR><TD>SuperLU</TD> <TD> 2.06807 </TD> <TD> 25.7064 </TD> <TD> 99.7166 </TD> <TD> 248.082 </TD>  </TR>
<TR><TD>GMRES (blocks solved by SuperLU)</TD> <TD> 1.73668 </TD> <TD> 11.1679 </TD> <TD> 42.2051 </TD> <TD> 105.97 </TD>  </TR>
<TR><TD>GMRES (LSC; blocks solved by SuperLU)</TD> <TD> 1.87403 </TD> <TD> 11.0113 </TD> <TD> 34.2697 </TD> <TD> 81.4382 </TD>  </TR>
<TR><TD>GMRES (LSC; pseudo-solid; blocks solved by Hypre/CG)</TD> <TD> 2.26537 </TD> <TD> 9.40194 </TD> <TD> 22.5644 </TD> <TD> 41.8251 </TD>  </TR>
</TABLE></CENTER>

<HR>
<HR>



\section modifications Modifications to the driver code
Most of the driver code is unchanged from the implementation
discussed in the  <a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">original
tutorial.</a> We therefore only discuss the modifications required
- to perform the fluid mesh update using (pseudo-)elasticity,
.
and
- to employ preconditioned GMRES as the linear solver in the Newton method.
.
<HR>

\subsection pseudo-solid Fluid-mesh udate by pseudo-elasticity
To perform the node update in the fluid mesh in response to the
deformation of the leaflet using the equations of
large-displacement elasticity we "wrap" the Navier-Stokes
element in the doubly-templated 
\c PseudoSolidNodeUpdateElement<WRAPPED_ELEMENT, SOLID_ELEMENT> class.
The first template argument of this class, \c WRAPPED_ELEMENT,
specifies the type of the "wrapped"
element -- here the nine-node, 2D quadrilateral \c
QTaylorHoodElement<2>. The second template argument, \c SOLID_ELEMENT,
specifies the \c SolidElement used to perform the node update -- here 
the nine-node, 2D quadrilateral \c QPVDElement<2,3>. Since the
preconditioner requires the (re-)classification of the pseudo-solid
degrees of freedom into constrained and unconstrained nodal positions,
we wrap this element too, using the
\c PseudoElasticBulkElement already discussed in
<a href="../../prescribed_displ_lagr_mult/html/index.html">
another tutorial</a>. The full specification of the 
fluid element is therefore given by
\code
PseudoSolidNodeUpdateElement<QTaylorHoodElement<2>,
                             PseudoElasticBulkElement<QPVDElement<2,3> > > 
\endcode


<HR> 


\subsection driver Modifications to the main code
The main code only requires minor changes, all associated with
the specification of different solver options. Since we provide the option to 
use Hypre and Trilinos solvers, we need to activate MPI
if \c oomph-lib has been compiled with MPI support (even 
if the code is run in serial):


\dontinclude fsi_channel_with_leaflet_precond.cc
\skipline start_of_main
\until #endif

Next we define and process the possible command line flags 
which are used to select solver options and to specify the spatial 
resolution 

\skipline Store command
\until validate

We then assign the command line flags that were specified
and document the ones that were recognised, before building the
problem, specifying the heavily wrapped version of the \c
Navier-Stokes element as the template parameter.

\skipline Parse command line
\until (mesh_multiplier);

The rest of the main code is very similar to the original version 
discussed in 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
another tutorial</a> and is therefore omitted here.

<HR>

\subsection it_solve Choosing the iterative solver and preconditioner
We provide a new function, \c set_iterative_solver() to switch the
solver to GMRES and to instantiate the FSI preconditioner.
We start by creating the GMRES solver, using Trilinos'
implementation if it is available:

\dontinclude fsi_channel_with_leaflet_precond.cc
\skipline start_iterative_solver
\until linear_solver_pt() = solver_pt;
 
Next, we create an instance of the FSI preconditioner for a 2D problem
and pass it to the solver:

\until  solver_pt->preconditioner_pt() = prec_pt;

The classification of the degrees of freedom within the block
preconditioning framework requires the specification of the
(pseudo-solid) fluid mesh, the "real" solid mesh and the
mesh containing the Lagrange multiplier elements that impose
the deformation of the FSI boundary in the fluid mesh:

\skipline Specify
\until Lagrange_multiplier_mesh_pt);

We provide the option to bypass the use of the LSC Schur complement
preconditioner as the subsidiary preconditioner (inexact solver)
for the Navier-Stokes block:

\skipline Use oomph-lib's
\until done disable lsc

By default the preconditioner performs the
preconditioning operations as described above, using SuperLU
for the block solves in (3). It is possible to specify 
alternative (approximate) solvers for these. For instance, the behaviour
of the LSC preconditioner can be modified by following the pointer
to the Schur complement preconditioner

\until navier_stokes_schur_complement_preconditioner_pt();

whose behaviour can be modified as discussed in 
<a href="../../../preconditioners/lsc_navier_stokes/html/index.html">
LSC preconditioner tutorial</a>. For instance, to use a block-triangular
approximation for the Navier-Stokes momentum block we specify

\until set_f_preconditioner(f_prec_pt);

If Hypre is available, we can then (approximately) solve the diagonal blocks 
using a small number of AMG V-cycles (with settings specified in the
helper function \c LSC_Preconditioner_Helper::set_hypre_preconditioner):

\until (LSC_Preconditioner_Helper::set_hypre_preconditioner);


Similarly, we can use Hypre (with settings appropriate for a 2D
Poisson problem) to approximately solve the pressure Poisson systems
in the computations of the Schur complement approximation:
 

\until #endif

Approximate solvers for the diagonal blocks in the preconditioner 
(2) for the the Lagrange-multiplier
constrained pseudo-elastic equations can be
chosen as discussed in 
<a href="../../prescribed_displ_lagr_mult/html/index.html">
another tutorial</a>. 
For instance, we can employ a block-upper triangular approximation
for the augmented pseudo-elastic block \f$E_{\rm aug}\f$:

\until PseudoElasticPreconditioner::Block_upper_triangular_preconditioner;

Its diagonal blocks can then again be solved by Hypre,

\until #endif

Similarly, the linear systems involving the mass matrices 
in \f$ W \f$ can be solved using diagonally preconditioned CG:

\until end set_iterative_solver

 
<HR>
<HR>

\section comm_and_ex Comments and Exercises

\subsection comm Comments
- Further details of the theory, such as the proof of the optimality of
  the preconditioner, can be found in
  <a href="http://dx.doi.org/10.1016/j.jcp.2012.07.001">
  Muddle, Mihajlovic & Heil (2012).</a>\n\n
- Note that we retained SuperLU as the solver for the linear system
  involving the tangent stiffness matrix \f$ S \f$ of the actual,
  fluid-load solid. In the present example this still leads to an
  optimal solver since the leaflet is modelled as a 1D beam structure
  whose number of unknowns increases much more slowly than the number
  of unknowns in the fluid mesh when the meshes are refined uniformly.
  If the solid is modelled as a "proper" solid, the use of a direct
  solver for the solution of the linear systems involving \f$ S \f$
  becomes sub-optimal. The member function 
  \code
  PseudoElasticFSIPreconditioner::set_solid_preconditioner(...)
  \endcode
  can then be used to specify a subsidiary preconditioner
  (inexact solver) for the linear systems involving \f$ S \f$.
  \n\n
.
      
\subsection ex Exercises

-# Experiment with different preconditioner settings and explore
   the performance of the preconditioner at different Reynolds
   numbers. The script 
   <A HREF="../../../../demo_drivers/interaction/pseudo_solid_fsi_channel_with_leaflet/run.sh"><CODE>run.sh</CODE></a>
   used to generate the data presented above, may be helpful.
.
<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/pseudo_solid_fsi_channel_with_leaflet">
demo_drivers/interaction/pseudo_solid_fsi_channel_with_leaflet
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/interaction/pseudo_solid_fsi_channel_with_leaflet/fsi_channel_with_leaflet_precond.cc">
demo_drivers/interaction/pseudo_solid_fsi_channel_with_leaflet/fsi_channel_with_leaflet_precond.cc</A>
</CENTER>

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


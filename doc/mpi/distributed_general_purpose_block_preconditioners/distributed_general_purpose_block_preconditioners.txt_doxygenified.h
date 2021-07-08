/**

\mainpage General-Purpose Block Preconditioners 

In this document we will demonstrate how to use the general-purpose block
preconditioners implemented in \c oomph-lib. This tutorial follows from the 
 <a href="../../../mpi/block_preconditioners/html/index.html">
Block Preconditioners</A> tutorial, which provides
 an overview of \c oomph-lib's generic block preconditioning framework. 


We use the \c Problem described in the
 <a href="../../../solid/airy_cantilever/html/index.html">
Bending of a Cantilever Beam</A>  tutorial to illustrate the key concepts.

\section introduction Introduction

In this section we define the four (distributed) general purpose block
preconditioning methodologies. To recap, all \c oomph-lib problems are
solved in a Newton iteration which requires the repeated solution of
linear systems of the form
\f[
J\delta {\bf x}=-{\bf r}
\f]
where \f$J\f$ is the Jacobian matrix, \f$\bf r\f$ is the vector of residuals and
\f$\bf \delta x\f$ is the Newton correction. We divide the DOFs in the two-dimensional cantilever problem into two
subsets corresponding to the \f$x\f$ and \f$y\f$ nodal positions.

\f[
\left[
\begin{array}{cc}
J_{xx}&J_{xy}\\
J_{yx}&J_{yy}
\end{array}
\right].
\left[
\begin{array}{c}
\bf \delta x \\
\bf \delta y \\
\end{array}
\right]
=
-\left[
\begin{array}{c}
\bf r_x \\
\bf r_y \\
\end{array}
\right]
\f]

Utilising this partitioning we will describe four (distributed)
general purpose block preconditioning methodologies. (Left) 
preconditioning represents a transformation of the original
linear system to

\f[
P^{-1}J\;{\bf \delta x}=-P^{-1}{\bf r}
\f]

with the aim of accelerating the convergence of Krylov subspace iterative methods
such as GMRES or CG. The application of the preconditioner requires the solution of

\f[
P{\bf z}={\bf w}
\f]

for \f$\bf z\f$ at each Krylov iteration.

\subsection block_diagonal Block Diagonal Preconditioning

We drop the off-diagonal blocks to form the block diagonal preconditioner

\f[
P_{BD}=\left[
\begin{array}{cc}
J_{xx}&\\
&J_{yy}
\end{array}
\right].
\f]

the application of this preconditioner requires the solution of the
subsidiary systems \f$J_{xx}\f$ and \f$J_{yy}\f$.

\subsection two_level_block_diagonal Block Diagonal Preconditioning with Two-Level Parallelisation

The two-subsidiary systems in the block diagonal
preconditioner (involving \f$J_{xx}\f$ and \f$J_{yy}\f$) can be solved
in any order. In a parallel computation we can either solve the two
systems one after the other using the full set of processes for the
solution of each linear system. An alternative is to solve all the 
subsidiary systems simultaneously, using only a
subset of processes for each system. We refer to this technique as
two-level parallelisation and note that this approach is particularly useful if the
linear solvers do not have good parallel scaling properties.

\subsection upper_block_triangular Upper Block Triangular Preconditioning

An alternative to block diagonal preconditioning is block triangular
preconditioning in which only off diagonal blocks on one side of the
diagonal are dropped. For example, in the block-upper triangular preconditioner

\f[
P_{BUT}=\left[
\begin{array}{cc}
J_{xx}&J_{xy}\\
&J_{yy}
\end{array}
\right]
\f]

the block below the diagonal (\f$J_{yx}\f$) has been dropped. In addition to the two subsidiary solves for \f$J_{xx}\f$  and
\f$J_{yy}\f$ this preconditioner requires a matrix-vector product involving \f$J_{xy}\f$.

\subsection lower_block_triangular Lower Block Triangular Preconditioning

Similarly we can define a lower triangular block preconditioner

\f[
P_{BLT}=\left[
\begin{array}{cc}
J_{xx}&\\
J_{yx}&J_{yy}
\end{array}
\right].
\f]


\section application Application


In this section we demonstrate the use of \c oomph-lib's general-purpose
block preconditioners. All general-purpose block preconditioners
are derived from the base class \c GeneralPurposeBlockPreconditioner (which
is itself derived from the \c BlockPreconditioner class).

By default all general purpose block preconditioners use \c SuperLUPreconditioner as
the preconditioner for the subsidiary systems (\f$J_{xx}\f$ and
\f$J_{yy}\f$ in the \ref introduction). \c
SuperLUPreconditioner is a 
 wrapper to both the 
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/#superlu"><code>SuperLU</code></a>
  direct solver and the 
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/#superlu_dist">
   <code>SuperLU Dist</code></a> distributed direct
  solver. Often we seek to replace this direct solver preconditioning
  with an inexact solver to make the preconditioner more efficient. To use an alternative subsidiary preconditioner we must
  define a function to return new instances of the chosen type of
  preconditioner (inexact solver). For example
\dontinclude airy_cantilever.cc
\skipline hypre_helper
\until end_of_hypre_helper

would return instances of \c HyprePreconditioner, a wrapper to the
  distributed <a href="https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html">
   <code>Hypre BoomerAMG</code></a> implementation of classical
  AMG. Later we will pass a pointer to this function to the block
  preconditioner to enable the use of \c HyprePreconditioner as a
  subsidiary preconditioner.
  Note that the function only creates the subsidiary preconditioner --
  it will be deleted automatically by the master preconditioner
  when it is no longer required.

The rest of the section is concerned with the main function, and in
particular setting up the preconditioner for use.

\dontinclude airy_cantilever.cc
\skipline start_of_main
\until {

Given an instance of the problem,

\dontinclude airy_cantilever.cc
\skipline Set up the problem
\until CantileverProblem<MySolidElement<RefineableQPVDElement<2,3> > > problem;

we specify GMRES as the linear solver. If available, we use the \c
TrilinosAztecOOSolver wrapper to the 
<a href="http://trilinos.sandia.gov/packages/aztecoo/"><code>Trilinos
AztecOO</code></a> implementation of GMRES. (This is the only
distributed implementation of GMRES in \c oomph-lib.)

\dontinclude airy_cantilever.cc
\skipline use trilinos gmres if available
\until #endif

\c GeneralPurposeBlockPreconditioner is the base class for all general
purpose block preconditioners. 

\dontinclude airy_cantilever.cc
\skipline Pointer to general purpose block preconditioner base class
\until GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* prec_pt = 0;

We introduced four general purpose block preconditioning methodologies in the
\ref introduction. The next step is to construct one of these preconditioners.
- \b Block \b Diagonal \b Preconditioning. This is implemented in the class \c BlockDiagonalPreconditioner.
\dontinclude airy_cantilever.cc
\skipline Standard Block Diagonal
\until prec_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
- \b Enabling \b Two-\b Level \b Block \b Diagonal \b Preconditioning. By default
  two-level preconditioning is disabled and hence \c
  enable_two_level_parallelisation() must have been called.
  Once this is done, each subsidiary system will be solved on an (as near to) equal size
  subset of processes.
\dontinclude airy_cantilever.cc
\skipline Two Level Block Diagonal
\until (prec_pt)->enable_two_level_parallelisation();
- \b Block \b Upper \b Triangular \b Preconditioning. Both block
  triangular preconditioners are implemented in the class \c
  BlockTriangularPreconditioner. By default this employs the
  upper-triangular version of the preconditioner.
\dontinclude airy_cantilever.cc
\skipline Block Upper Triangular
\until prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;
- \b Block \b Lower \b Triangular \b Preconditioning. The lower
triangular version of the preconditioner can be selected with a call
to the method \c lower_triangular().
\dontinclude airy_cantilever.cc
\skipline Block Lower Triangular
\until (prec_pt)->lower_triangular();


Having chosen a preconditioner structure, the next stage is to choose
the preconditioner for the subsidiary systems (\f$J_{xx}\f$ and
\f$J_{yy}\f$ in the \ref introduction ). By default this is \c
SuperLUPreconditioner, but we wish to use \c HyprePreconditioner so we
pass the previously specified function \c
Hypre_Subsidiary_Preconditioner_Helper::get_hypre_preconditioner() to
the preconditioner.

\dontinclude airy_cantilever.cc
\skipline  Specify Hypre as the subsidiary block preconditioner
\until (Hypre_Subsidiary_Preconditioner_Helper::get_hypre_preconditioner);

The same subsidiary preconditioner is used for all subsidiary
  systems in a general purpose block preconditioner.

As discussed in the 
<a href="../../../mpi/block_preconditioners/html/index.html">
Block Preconditioners</A>
 tutorial, the classification of the DOFs
is implemented at an elemental level so we pass a pointer to the mesh
containing the elements to the preconditioner. (Note that this problem 
contains two meshes, one containing the bulk elements and one
containing the FaceElements that apply the traction boundary condition.
Since the latter do not introduce any new DOFs, all the DOFs are 
classified by the bulk elements. Therefore, we do not need to pass the 
traction element mesh to the block preconditioner.)

\dontinclude airy_cantilever.cc
\skipline The preconditioner only requires the bulk mesh since its
\until prec_pt->add_mesh(problem.solid_mesh_pt());
 
Finally, we pass the preconditioner to the solver

\dontinclude airy_cantilever.cc
\skipline pass the preconditioner to the solver
\until solver_pt->preconditioner_pt() = prec_pt;

and solve the problem:

\dontinclude airy_cantilever.cc
\skipline solve the problem
\until problem.newton_solve();

\section Parallelisation

Given that \c BlockPreconditioner, \c
TrilinosAztecOOSolver, \c SuperLUPreconditioner, \c
HyprePreconditioner and \c MatrixVectorProduct are all automatically 
distributed, all that is
required for a distributed solution is to run the executable under
MPI with multiple processes.


<HR>
<HR>

\section sources Source files for this tutorial
- The source files for the driver code are in
<CENTER>
<A HREF="../../../../demo_drivers/mpi/solvers/">
demo_drivers/mpi/solvers/
</A>
</CENTER>
\n
- The driver code is
<CENTER>
<A HREF="../../../../demo_drivers/mpi/solvers/airy_cantilever.cc">
demo_drivers/mpi/solvers/airy_cantilever.cc
</A>
</CENTER>
\n
.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


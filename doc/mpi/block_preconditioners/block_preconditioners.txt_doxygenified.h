/**

\mainpage oomph-lib's Block Preconditioning Framework
 
   
\c oomph-lib's block preconditioning framework
provides an environment for the development of block 
preconditioners for the iterative solution of linear systems
by Krylov subspace methods. The framework is based on the
classification of the problem's unknowns (degrees of freedom; here
abbreviated as dofs)
into different "dof types" which, in a multi-physics context, 
typically represent different physical quantities. A key feature of
the framework is that
it allows existing block preconditioners (which were developed for a particular
single-physics application, say) to be re-used, in a hierarchical 
fashion, in block preconditioners for related
multi-physics problems. This means that existing Navier-Stokes
and solid mechanics preconditioners can be used to create
preconditioners for fluid-structure interaction problems, say. 

Following a brief overview of the underlying
ideas and their implementation in \c oomph-lib this tutorial discusses 
a sequence of increasingly complex block preconditioners that illustrate the
framework's capabilities in the context of a (rather artificial) model
problem. The final example illustrates a simple implementation of
a block preconditioner for an FSI problem. We conclude with a few 
comments on the use of block preconditioners in parallel. 
Other tutorials discuss how the methodology is used in "real"
preconditioners. See, for instance, the tutorials discussing
 
- \c oomph-lib's 
  <a href="../../distributed_general_purpose_block_preconditioners/html/index.html">"general purpose" 
  block preconditioners.</a> \n\n
- The <a href="../../../preconditioners/lsc_navier_stokes/html/index.html">
  NavierStokesSchurComplementPreconditioner for Navier-Stokes problems
  </a> \n\n
- The <a href="../../../preconditioners/fsi/html/index.html"> FSIPreconditioner
  for monolithically-discretised fluid-structure interaction problems.</a>
  \n\n
- The 
  <a href="../../../preconditioners/prescribed_displ_lagr_mult/html/index.html">
  preconditioner for large-displacement solid mechanics problems 
  in which boundary displacements are prescribed.</a>\n\n
- The previous preconditioner is mainly used as a subsidiary
  block preconditioner for the 
  <a href="../../../preconditioners/pseudo_solid_fsi/html/index.html">
  solution of fluid-structure interaction problems with (pseudo-)solid 
  fluid mesh updates.</a> \n\n
.

\section theory Theoretical background

In \c oomph-lib, all problems are solved by Newton's method,
which requires the repeated solution of linear systems of the form

\f[
{\bf J}\;{\bf \delta x}=-{\bf r}
\f]

for the Newton correction \f$\bf \delta x\f$ where \f${\bf J}\f$ is the
Jacobian matrix and \f$\bf r\f$ is the vector of residuals. (Left) 
preconditioning represents a transformation of the original linear system to

\f[
{\bf P}^{-1}{\bf J}\;{\bf \delta x}=-{\bf P}^{-1}{\bf r},
\f]

introduced with the aim of accelerating the convergence of Krylov subspace 
solvers such as GMRES or CG. The application of the preconditioner 
requires the solution of

\f[
{\bf P}{\bf z}={\bf y}
\f]

for \f$\bf z\f$ at each Krylov iteration.

Block preconditioners are based (at least formally) on a
reordering of the linear system such that related unknowns
(e.g. dofs representing the same physical quantity)
are grouped together and enumerated consecutively. 

For instance, in linear elasticity problems (discussed in 
<a href="../../../linear_elasticity/periodic_load/html/index.html">another
tutorial</A>) where we compute the displacement field of an elastic
body in response to an applied traction, the (discrete) unknowns
can be sub-divided according to which component of the displacement
vector they represent. Using this classification of the dofs, 
the re-ordered linear system for a two-dimensional problem then has the form

\f[
\left[ 
\begin{array}{cc}
{\bf J}_{xx}&{\bf J}_{xy}\\
{\bf J}_{xy}&{\bf J}_{yy}
\end{array}
\right]
\left[ 
\begin{array}{c}
\bf \delta x_x\\
\bf \delta x_y\\
\end{array}
\right]
=
-
\left[
\begin{array}{c}
\bf r_x\\
\bf r_y\\
\end{array}
\right].
\f]

A simple (and, in fact, quite effective) block preconditioner for this
linear system can be formed by retaining only the diagonal blocks
of the system matrix, leading to the block diagonal preconditioner 

\f[
{\bf P}_{diag}=
\left[
\begin{array}{cc}
{\bf J}_{xx}& \\
      &{\bf J}_{yy}
\end{array}
\right].
\f]


The application of this preconditioner requires the solution of 
the linear system

\f[
\left[
\begin{array}{cc}
{\bf J}_{xx}& \\
      &{\bf J}_{yy}
\end{array}
\right]
\left[
\begin{array}{c}
{\bf z_x}\\
\bf z_y\\
\end{array}
\right]
=
\left[
\begin{array}{c}
\bf y_x\\
\bf y_y
\end{array}
\right],
\f]

which requires the (exact or approximate) solution of the two 
smaller linear systems \f$ {\bf J}_{xx} \ {\bf z_x} = {\bf y_x}\f$
and \f$ {\bf J}_{yy} \ {\bf z_y} = {\bf y_y}\f$.


\section generic_implementation Overview

The above example shows that the application of block preconditioners
typically require several generic steps:

-# The classification of the dofs.
-# The application of subsidiary preconditioning operations such as 
   the solution of (smaller) linear systems or the evaluation of matrix-vector 
   products with some of the blocks that are extracted from the
   original linear system.
   
The following subsections describe how these tasks are performed
within \c oomph-lib's block preconditioning framework.

\subsection block_preconditionable_elements The classification of dof types via block preconditionable elements

The classification of dofs is specified by the elements since they are
the only objects within \c oomph-lib's data structure 
that "know" what role a specific dof plays in "their" equations. 
During the setup phase, the block preconditioner loops over "all
elements" (specified via one or more \c Meshes -- here simply used as 
containers for elements; see below for further details) to establish the 
"dof type" for each global unknown.
 
To achieve this, the class \c GeneralisedElement 
contains two broken virtual methods that must be re-implemented/overloaded to
label each of the element's dofs with its type. These methods are:

- \c GeneralisedElement::ndof_types() must return the number of dof types
  associated with an element.
- \c GeneralisedElement::get_dof_numbers_for_unknowns(...) must return
  a list of pairs comprising a map from global equation number to 
  dof type for all unknowns in the element.

These are already implemented for many elements. If not, the functions
are easy to write. For instance, \c oomph-lib's <c>DIM</c>-dimensional 
linear elasticity elements from the \c QLinearElasticityElement family
can be made block-preconditionable by using the following
wrapper class:

\dontinclude two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc
\skipline  start_of_mylinearelasticityelement
\until };

Thus, in the two-dimensional \c MyLinearElasticityElement<2> we 
have two types of dofs, corresponding to the displacements in the 
\f$x\f$ and \f$ y \f$ directions, respectively. They are enumerated 
as dof types 0 and 1, respectively.

\subsection dof_types_and_block_types dof types, blocks, compound blocks and meshes

In the block diagonal preconditioner for the two-dimensional linear
elasticity problem, discussed above, we have dof types that 
correspond directly to
the blocks in the (re-ordered) Jacobian matrix. However, as we will
demonstrate <a href="#compound">below</a>, it is also 
possible to combine the blocks associated with multiple dofs into a
single (compound) block in which case the number of blocks is smaller than the
number of dof types. The relationship between dof types, block types, 
the elemental dof type classification and meshes are as follows

- \b Elemental \b dof \b type \b classification: Each element
  classifies its own dof types in the function 
  \c get_dof_numbers_for_unknowns(...). In the case of the two-dimensional
  \c MyLinearElasticityElement<2> elements, the dof types are classified
  as \c 0 and \c 1; for two-dimensional \c QTaylorHoodElement<2> 
  Navier-Stokes 
  elements, the dof types
  are classified as \c 0 and \c 1 for the \f$ x \f$ and \f$y
  \f$-velocities, and \c 2 for the pressure \f$ p \f$; etc.
  \n\n
- \b Role \b of \b meshes: When classifying the degrees of freedom into
  dof types, the block preconditioning framework
  visits all elements that make contributions to the Jacobian matrix
  and associates the global equation number of each dof with the dof type
  specified by the element. The block preconditioning framework is
  given access to the elements via (possibly multiple) meshes (here
  simply interpreted as containers for elements), each of 
  which is assumed to contain elements of a single type. The total
  number of dof types in the block preconditioner is the
  sum of the dof types of the elements in the meshes. For
  instance, in a 2D fluid-structure interaction problem we have two
  different element types, the solid elements (which contain the
  \f$ x \f$ and \f$ y\f$ solid displacements, \f$ u_x \f$ and \f$
  u_y \f$, respectively, assumed to be 
  enumerated as dof types
  0 and 1 by these elements) and the fluid elements (which contain the
  \f$ x \f$- and \f$ y \f$ - fluid velocities, \f$ v_x \f$ and \f$ v_y
  \f$ ,  and the pressure, \f$ p \f$, assumed to 
  be enumerated as dof types
  0, 1 and 2 by these elements). Assuming the mesh of solid elements
  is specified as mesh 0 and the mesh of fluid elements is mesh 1, the
  block preconditioner has a total of five dof types which represent,
  in order, \f$ u_x, u_y, v_x, v_y, p\f$. Note that if certain degrees of
  freedom are classified by multiple elements, the most recent
  assignment of the dof type over-writes previous assignments.
  The order in which meshes are specified therefore matters. 
  \n
  A corollary to this is that a block preconditioner does not
  need to "know" about elements that do not introduce any new 
  unknowns. For instance, \c FaceElements that apply Neumann/flux 
  boundary conditions operate on dofs that are already
  contained in (and therefore classified by) the elements in the 
  "bulk" mesh. Conversely, if a \c FaceElement
  imposes a boundary condition via Lagrange-multipliers, the dofs
  that represent these Lagrange multipliers must be
  classified by the \c FaceElements since the "bulk elements" are not
  aware of them. 
  \n
  If \c oomph-lib is compiled with the \c PARANOID flag, 
  an error is thrown if any of the global unknowns are not associated
  with a dof type.\n\n
- \b Blocks: The blocks are the sub-blocks of the system matrix
  (usually the Jacobian matrix from the Newton method) that the block
  preconditioner works with. By default, each block is associated with
  exactly one dof type. However, it is possible create "compound blocks" that 
  are associated with more than one dof type. For example, in the 
  Navier-Stokes LSC preconditioner (in 2D) we have three dof 
  types (the \f$x\f$ and
  \f$y\f$-velocities and the pressure), but the preconditioner works with 
  just two block types (forming the velocity and pressure blocks). The setup
  of the block types is handled by the function \c block_setup(...) discussed
  below.
.

\section multi_poisson Simple preconditioner examples

We will now illustrate the capabilities of the block preconditioning
framework by considering the system of \f$ N \f$ coupled PDEs 
\f[ 
\left( \frac{\partial^2 u_i}{\partial x_j^2} 
+ \beta \sum_{k=1}^{N} u_k \right)  = f_i(x_j) \ \ \ \ i=1,...,N
\ \ \ \ \ \ \ \ \ \ \ \ (1)
\f]
for the  \f$ N \f$ fields  \f$ u_i(x_j) \f$. If  \f$ \beta=0 \f$,
the system represents  \f$ N \f$ (uncoupled) Poisson equations,
each with their own source function  \f$ f_i(x_j). \f$ If  \f$ \beta
\ne 0 \f$ the PDE for \f$ u_i(x_j) \f$ is affected by all other
fields via the Helmholtz-like second term on the left-hand-side.


The \c MultiPoissonElements discretise the equations 
with standard Galerkin-type finite elements in which each field
is treated as its own dof type. If  \f$ N = 5, \f$ the linear system 
to be solved in the course of the Newton method,
\f[
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
{\bf J} \ \delta {\bf x} = -{\bf r},
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
(2)
\f]
has a \f$5 \times 5\f$  block structure implying that, following a
formal re-numbering
of the unknowns, the matrix and the vectors can be written as
\f[
\hspace{3cm}
{\bf J} = 
\left(
\begin{array}{ccccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
J_{21} & J_{22} & J_{23} & J_{24} & J_{25} \\
J_{31} & J_{32} & J_{33} & J_{34} & J_{35} \\
J_{41} & J_{42} & J_{43} & J_{44} & J_{45} \\
J_{51} & J_{52} & J_{53} & J_{54} & J_{55} \\
\end{array}
\right),
 \ \
\delta {\bf x} = 
\left(
\begin{array}{c}
\delta x_{1}  \\
\delta x_{2}  \\
\delta x_{3}  \\
\delta x_{4}  \\
\delta x_{5}  \\
\end{array}
\right)
 \ \mbox{ and } \ 
{\bf r} = 
\left(
\begin{array}{c}
r_{1}  \\
r_{2}  \\
r_{3}  \\
r_{4}  \\
r_{5}  \\
\end{array}
\right).
\hspace{3cm}
(3)
\f]
We wish to solve this linear system by preconditioned Krylov subspace methods,
using a block preconditioner \f${\bf P}\f$ formed (formally) 
from the blocks of the system matrix \f${\bf }J\f$. As
discussed above, the
application of the preconditioner (typically once per iteration of the
Krylov solver) then requires the solution of linear systems of the form
\f$
{\bf P} {\bf y} = {\bf z},
\f$ 
for \f${\bf y}\f$. The preconditioning operation can also be written as
\f$
{\bf y} = {\bf P}^{-1} {\bf z}
\f$ 
where the operator \f${\bf P}^{-1}\f$ represents 
the application of the preconditioner to a vector \f${\bf z}\f$. 
Formally, the operator \f${\bf P}^{-1}\f$ represents the inverse of 
the matrix \f${\bf P}\f$ but its application may, of course, be performed  
approximately by another ``subsidiary'' preconditioner/inexact solver
e.g. by performing a small number of multigrid cycles, say.
(Note that we say ``formally'' because the preconditioner does not
actually have to be associated with a specific matrix -- it simply has
to act as a linear operator that ``turns \f${\bf z}\f$ into \f${\bf y}\f$'').

A specific block preconditioner must be derived from the
\c BlockPreconditioner base class and must implement two 
pure virtual member functions of the underlying \c Preconditioner
class:
- \c void \c Preconditioner::setup(): This function is
called once during the solution of a given linear system
by any of \c oomph-lib 's Krylov subspace solvers. It typically 
extracts a  certain number of blocks from the matrix \f${\bf J}\f$, possibly 
manipulates its local copies of these blocks, and performs any 
preliminary computations required 
to allow the rapid subsequent application of \f${\bf P}^{-1}\f$.
- \c void \c Preconditioner::preconditioner_solve(\f$\bf{z}\f$,\f$\bf{y}\f$):
This function applies \f${\bf P}^{-1}\f$ to the input argument
\f${\bf z}\f$ and computes \f${\bf y}\f$, typically using some data that has
been pre-computed in the \c setup() function.

To allow a block preconditioner to classify all dofs, the
preconditioner must be given access to all elements that contribute 
to the linear system to be solved. This is done by specifying pointers to these
elements via one or more \c Mesh objects (which simply act as
containers for the element pointers), using the functions \c set_nmesh(...)
(which specifies how many meshes the preconditioner works with) and
\c set_mesh(...).

We will discuss the implementation of the required functions (and associated
capabilities of the block preconditioning framework)
in a number of increasingly complex block preconditioners
for the solution of the \f$5 \times 5\f$ linear system 
defined by equations (2) and
(3). We stress that the purpose
of this exercise is not the development of particularly clever
preconditioners but simply an excuse to demonstrate the use of the
available ``machinery''. Specifically we will demonstrate how to:

- extract selected blocks from the system matrix (usually the Jacobian
  matrix assembled by the Newton solver).

- perform matrix vector products with selected off-diagonal blocks.

- solve linear systems associated with selected blocks, using either a
  direct solver and/or subsidiary preconditioners (inexact solvers), 
  including cases where the
  subsidiary preconditioners are block preconditioners themselves.

- replace and modify selected blocks and how to make such
  modified blocks available to subsidiary block preconditioners.

- concatenate and coarsen blocks.




<HR>
<HR>




\subsection diagonal A block diagonal preconditioner

NEW FEATURES: How to extract matrix blocks and corresponding block vectors
from their full-size counterparts

\subsubsection diag_theory Theory

The simplest possible block preconditioner is a block-diagonal 
preconditioner, formed by retaining only the diagonal blocks of
\f${\bf J}\f$, so that
\f[
\hspace{3cm}
{\bf P} = 
\left(
\begin{array}{ccccc}
J_{11} & & & & \\
& J_{22} &  &  & \\
 & & J_{33} & & \\
 & & & J_{44} & \\
 &  & &  & J_{55} \\
\end{array}
\right).
\hspace{3cm}
(4)
\f]
The application of this preconditioner (i.e. the solution of the
linear system \f${\bf P} {\bf
  y} ={\bf z}\f$ for \f${\bf y}\f$) requires the solution of the 
five much smaller linear systems
\f[
\hspace{3cm}
\begin{array}{c}
{\bf J}_{11} \ {\bf y}_1 = {\bf z}_1, \\
{\bf J}_{22} \ {\bf y}_2 = {\bf z}_2, \\
{\bf J}_{33} \ {\bf y}_3 = {\bf z}_3, \\
{\bf J}_{44} \ {\bf y}_4 = {\bf z}_4, \\
{\bf J}_{55} \ {\bf y}_5 = {\bf z}_5, \\
\end{array}
\hspace{3cm}
(5)
\f]
where we have assumed that the two vectors \f${\bf y}\f$ and \f${\bf z}\f$
are re-ordered into ``block vectors'' in the same way as 
the vectors \f$\delta {\bf x}\f$ and \f${\bf r}\f$ in 
"the original linear system" (2) are re-ordered into
the ``block vectors'' in (3).

The implementation of the preconditioning operations in 
(5)  can naturally be subdivided into two
distinct \c setup() and \c preconditioner\_solve(...) phases.
Assuming that the linear systems in (5) are solved exactly 
by a direct solver (an ``exact preconditioner'') that 
can pre-compute and store the LU decomposition of the diagonal matrix blocks,
the \c setup() phase involves the following operations
[text in square brackets refers to their \c oomph-lib -specific 
implementation]:
- Set up any data structures/lookup tables that are required to 
extract matrix blocks from the original matrix
\f${\bf J}\f$ [by calling
the \c BlockPreconditioner::block\_setup() function]. 
- Extract the five diagonal blocks \f${\bf J}_{ii}\f$ 
(for \f$i=1,...,5\f$) [using the 
\c BlockPreconditioner::get\_block(...) function].
- Compute and store the LU decomposition of the diagonal blocks
to allow the rapid solution of the systems \f${\bf J}_{ii} \ {\bf y}_i = 
{\bf z}_i\f$ (for \f$i=1,...,5\f$) during the \c preconditioner\_solve(...)
phase by back-substitution. [This is done by calling the \c setup(...)
function of the subsidiary preconditioner/inexact solver. Following this,
the diagonal matrix blocks are longer required and can be deleted.]

Once the \c setup() phase has been completed, the solution of the linear
system \f${\bf P} {\bf y} = {\bf z}\f$ by the \c preconditioner\_solve(...) 
function involves the following steps:
- Extract the five ``block vectors'' \f${\bf z}_i\f$ (for \f$i=1,...,5\f$)
from the vector \f${\bf z}\f$ [using the 
\c BlockPreconditioner::get\_block\_vectors(...) function].
- Solve the linear systems \f${\bf J}_{ii} \ {\bf y}_i = 
{\bf z}_i\f$ for the vectors \f${\bf y}_i\f$ (for \f$i=1,...,5\f$)  using 
the precomputed LU decomposition of the diagonal blocks \f${\bf J}_{ii}\f$ 
(for \f$i=1,...,5\f$) created during the \c setup() phase.
- Combine the five ``block vectors'' \f${\bf y}_i\f$ (for \f$i=1,...,5\f$)
to the full-length vector \f${\bf y}\f$ [using the 
\c BlockPreconditioner::return\_block\_vectors(...) function].

\subsubsection diag_implementation Implementation as a BlockPreconditioner
Here is a sample implementation of the diagonal block preconditioner
as a class \c Diagonal, derived from the \c BlockPreconditioner base 
class. The class provides storage for the subsidiary preconditioners
that solve the linear systems associated with the diagonal blocks, 
implements the \c setup() and \c preconditioner\_solve(...)
functions, and provides a helper function \c clean\_up\_my\_memory()
which does what it says. We also provide an access function
which allows the user to specify the pointer to the \c Mesh that
contains the \c MultiPoissonElements which classify the dofs.

\dontinclude multi_poisson_block_preconditioners.h
\skipline start_of_diagonal_class
\until };


\subsubsection diag_setup The setup() function

As mentioned above, a \c Preconditioner's \c setup() function is
called at the beginning of the \c IterativeLinearSolver's 
\c solve(...) function. In time-dependent and/or nonlinear
problems many (different) linear systems have to be solved by the same
linear solver (and the associated preconditioner) throughout the code 
execution. To avoid memory leaks it is therefore 
important to free up any memory that may have been allocated in any
previous use of the preconditioner. The \c setup() function of 
all block preconditioner should therefore always 
start by freeing up such memory. This is best done by using a helper
function that can also be called from the destructor. 

\skip start_of_setup_for_simple
\until this->clean_up_my_memory();

Next we set the pointer to the preconditioner's one-and-only mesh,
and call the \c block_setup() function to set up the internal data
structures and lookup tables required to extract blocks from the
system matrix.

\until this->block_setup();

We create five subsidiary preconditioners (all exact solvers -- SuperLU
in its incarnation as an "exact" preconditioner) for the solution of the
linear systems involving the diagonal blocks:

\skip //
\until }


Next we set up the subsidiary preconditioner by
extracting the diagonal blocks from the system matrix
and passing them to the subsidiary preconditioners. 

Note that each preconditioner is expected to retain a copy
of whatever data it needs to subsequently perform its 
\c preconditioner_solve(...)
function. The deep copy of the block that is returned by the
\c get\_block(...) function can therefore be deleted (here simply
go out of scope) once the subsidiary preconditioner has been set up.
(In the specific case of the \c SuperLUPreconditioner, the \c setup(...)
function computes and stores the LU decomposition of the matrix;
the matrix itself is then no longer required). 

\skip //
\until }


\subsubsection diag_solve The preconditioner_solve() function

To apply the preconditioner to a given vector, r, we first
extract the five block-vectors whose sizes (and permutations) match
that of the diagonal matrix blocks, using the
\c get_block_vectors(...) function.


\skip //======
\until this->get_block

We then provide storage for the five solution vectors and compute them
by applying the subsidiary preconditioners' \c
preconditioner_solve(...) function:

\skip //
\until }

Finally the solutions in \c block_z are returned into the 
full-length solution vector \c z via a call to 
\c return_block_vectors(...). 

\skip //
\until }

\subsubsection diag_clean The clean_up_my_memory() function
This function (which is called by the \c setup()
function and the destructor) frees the memory that is
allocated when a new linear system is solved -- here the
subsidiary preconditioners and their associated data (the
LU decompositions of the diagonal blocks).


\skip //
\until // End of clean_up_my_memory function.

\subsubsection diag_comments Comments and Exercises

- \c The function \c get_block_vectors(r,block_r) extracts the five
  (or, in general, \c nblock_type() ) block vectors \c block_r from the
  full-length vector \c r. The sizes of the block vectors (and the 
  permutation of their entries relative
  to their order in the full length vector \c r) match that of the
  matrix blocks. There is an alternative function \c get_block_vector(...)
  (note the missing s) which extracts a single block vector.
  An equivalent version exists for the \c return_block_vector[s]
  functions. 
- In the example above we used an "exact preconditioner" (direct solver)
  to solve the five linear systems associated with the diagonal blocks.
  However, the (approximate) solution of these linear systems can 
  be performed by any other matrix-based preconditioner, such as 
  \c oomph-lib's diagonal preconditioner, \c
  MatrixBasedDiagPreconditioner, discussed in 
  <a href="../../../linear_solvers/html/index.html">another tutorial</a>.
  The setup and application of this
  preconditioner is obviously much faster than for the \c
  SuperLUPreconditioner. Its setup merely requires the extraction of the 
  diagonal  entries and storage of their inverses (rather than the
  computation of the LU decomposition), while the application 
  simply requires the multiplication of the input vector by the 
  pre-computed inverses of the diagonal entries (rather than a
  back-substitution). However, the preconditioner is clearly not
  "as good" and therefore results in a larger number of iterations
  in the iterative linear solver. In fact, using the diagonal preconditioner
  for the approximate solution of the five linear systems involving the
  diagonal blocks is mathematically equivalent to using the
  diagonal preconditioner on the entire matrix. Try it out! 



<HR>
<HR>




\subsection upper_triangular A block upper triangular preconditioner

NEW FEATURES: How to set up matrix vector products with off-diagonal blocks

\subsubsection upper_triangular_theory Theory

Next we consider the implementation of an upper triangular
preconditioner, formed by retaining only the blocks in the upper right
hand part of \f${\bf J}\f$, including the diagonals.
\f[
\hspace{3cm}
{\bf P} = 
\left(
\begin{array}{ccccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
& J_{22} &  J_{23} &  J_{24} & J_{25} \\
 & & J_{33} & J_{34} & J_{35} \\
 & & & J_{44} & J_{45} \\
 &  & &  & J_{55} \\
\end{array}
\right).
\hspace{3cm}
(6)
\f]
The application of this preconditioner (i.e. the solution of the
linear system \f${\bf P} {\bf
  y} ={\bf z}\f$ for \f${\bf y}\f$) again requires the solution of 
five much smaller linear systems
\f[
\hspace{3cm}
\begin{array}{l}
{\bf J}_{11} \ {\bf y}_1 = \widetilde{\bf z}_1 = {\bf z}_1 - {\bf J}_{15} \ {\bf y}_5 - {\bf J}_{14} \ {\bf y}_4 - {\bf J}_{13} \ {\bf y}_3 - {\bf J}_{12} \ {\bf y}_2, \\
{\bf J}_{22} \ {\bf y}_2 = \widetilde{\bf z}_2 = {\bf z}_2 - {\bf J}_{25} \ {\bf y}_5 - {\bf J}_{24} \ {\bf y}_4 - {\bf J}_{23} \ {\bf y}_3, \\
{\bf J}_{33} \ {\bf y}_3 = \widetilde{\bf z}_3 ={\bf z}_3 - {\bf J}_{35} \ {\bf y}_5 - {\bf J}_{34} \ {\bf y}_4, \\
{\bf J}_{44} \ {\bf y}_4 = \widetilde{\bf z}_4 ={\bf z}_4 - {\bf J}_{45} \ {\bf y}_5, \\
{\bf J}_{55} \ {\bf y}_5 = \widetilde{\bf z}_5 = {\bf z}_5, \\
\end{array}
\hspace{3cm}
(7)
\f]
where we have again assumed that the two vectors 
\f${\bf y}\f$ and \f${\bf z}\f$
are re-ordered into ``block vectors'' in the same way as 
the vectors \f$\delta {\bf x}\f$ and \f${\bf r}\f$ in 
"the original linear system" (2) are re-ordered into
the ``block vectors'' in (3).

The main difference to the block diagonal preconditioner considered
before is that the right hand sides of the linear systems have to be
modified. We start by solving the final equation for \f$ {\bf y}_5
\f$. We then multiply this vector by the off-diagonal block \f$ {\bf J}_{45}
\f$, subtract the result from \f$ {\bf z}_4 \f$ and use the result of
this operation as the right-hand-side for the linear system that determines
\f$ {\bf y}_4\f$,  etc.

The implementation of the preconditioning operations in 
(7)  can again be subdivided into two
distinct \c setup() and \c preconditioner\_solve(...) phases.
Assuming that the linear systems in (7)
are solved exactly by a direct solver (an ``exact preconditioner'')
that can pre-compute and store the LU decomposition of the diagonal
matrix blocks, the \c setup() phase involves the following operations
[text in square brackets refers to their \c oomph-lib specific 
implementation]:
- Set up any data structures/lookup tables that are required to 
extract matrix blocks from the original matrix
\f${\bf J}\f$ [by calling
the \c BlockPreconditioner::block\_setup() function]. 
- Extract the five diagonal blocks \f${\bf J}_{ii}\f$ 
(for \f$i=1,...,5\f$) [using the 
\c BlockPreconditioner::get\_block(...) function].
- Compute and store the LU decomposition of the diagonal blocks
to allow the rapid solution of the systems \f${\bf J}_{ii} \ {\bf y}_i = 
\widetilde{\bf z}_i\f$ (for \f$i=1,...,5\f$) during the \c preconditioner\_solve(...)
phase by back-substitution. [This is done by calling the \c setup(...)
function of the subsidiary preconditioner/inexact solver. Following this,
the diagonal matrix blocks are longer required and can be deleted.]
- Extract the relevant off-diagonal blocks from \f${\bf J}\f$ and create 
\c MatrixVectorProduct operators. [The matrix vector products are 
set up using the 
\c setup_matrix_vector_product(...) function. As with the subsidiary
preconditioners, the \c MatrixVectorProduct operators retain their own
copy of any required data, so the off-diagonal matrix blocks can be
deleted (or be allowed to go out of scope) following the setup.]

\subsubsection upper_triangular_implementation Implementation as a BlockPreconditioner
Here is a sample implementation of the upper triangular block preconditioner
as a class \c UpperTriangular, derived from the \c BlockPreconditioner base 
class. The class provides storage for the subsidiary preconditioners
that solve the linear systems associated with the diagonal blocks, and
the \c MatrixVectorProduct operators. We also implement the \c
setup() and \c preconditioner\_solve(...) functions, and provide a
helper function \c clean\_up\_my\_memory() which does what it says.
As before we also provide an access function
which allows the user to specify the pointer to the \c Mesh that
contains the \c MultiPoissonElements which classify the dofs.

\dontinclude multi_poisson_block_preconditioners.h
\skipline start_of_upper_triangular_class
\until };

\subsubsection upper_triangular_setup The setup() function
As before, we start by cleaning up the memory, set the pointer to the 
mesh, and set up
the generic block preconditioner functionality by calling 
\c block_setup().

\skip start_of_setup_for_upper_triangular
\until this->block_setup();

We provide storage for the (pointers to the) matrix vector products
and the subsidiary preconditioners.

\skip //
\until Block_preconditioner_pt.resize(nblock_types);

Next we create the subsidiary preconditioners which we will use to 
solve the linear systems involving the diagonal blocks.

\until end of brace

We then extract the relevant off-diagonal blocks (those above the
diagonal) from the full matrix, create a \c MatrixVectorProduct 
operator for each and 
use the \c BlockPreconditioner::setup_matrix_vector_product(...) function
to make them fully functional. Note that the final argument to this
function (the column index of the off-diagonal block in its block enumeration
within the current preconditioner) is required to set up additional lookup
tables that are required to ensure the correct operation of this object
in cases when the preconditioner operates in parallel. The details are
messy and not worth explaining here -- just do it!

\skip //
\until // End setup

\subsubsection upper_triangular_solve The preconditioner_solve() function
As in the block diagonal preconditioner, we start by extracting the
block vectors from the full-length vector, \c r.

\skip //======
\until this->get_block

Next we provide storage for the solution vectors and 
work backwards through the (block)-rows of the (block-)linear
system (7). Following each linear solve
we update the right-hand-side of the next linear system, as discussed
above.

\skip //
\until // End for over i

Finally, the solutions in \c block_z are combined via \c
return_block_vectors(...) which places the results back
into the full-length vector \c z that is returned by this function.

\skip //
\until }

\subsubsection upper_triangular_clean The clean_up_my_memory() function
This  function again deletes any data that was
allocated in the setup function -- here the subsidiary preconditioners
(and their LU decompositions) and the matrix-vector product operators.

\skip //
\until // End of clean_up_my_memory function.






<HR>
<HR>



<a name="compound"></a>
\subsection two_plus_three Combining multiple dof types into compound blocks. Part 1


NEW FEATURES: How to combine multiple dof types into compound blocks


\subsubsection two_plus_three_theory Theory
So far we have illustrated how to implement block preconditioners 
for cases where the dof types (as identified by the elements) 
correspond directly to 
the block types. This is appropriate for our model PDE system 
(1) in which the five fields (and the governing
equations) are all of the same type. In many applications,
particularly in multi-physics problems, it
may be desirable to combine similar/related dof types into single
blocks. For instance, in a 2D fluid-structure interaction problem,
we may wish to distinguish between the two solid (x and y solid
displacements) and three fluid (x and y fluid velocities and the pressure)
dofs and employ subsidiary preconditioners that act directly on 
the two distinct solid and fluid blocks. A basic block diagonal 
preconditioner for such a problem that ignores the coupling between 
fluid and solid dofs has the following structure 
\f[
\hspace{3cm}
{\bf P} = 
\left(
\begin{array}{cc|ccc}
J_{11} & J_{12} & & & \\
J_{21} & J_{22} &  &  & \\
\hline 
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{43} & J_{44} & J_{45} \\
 & & J_{53} & J_{54}  & J_{55} \\
\end{array}
\right)
=
\left(
\begin{array}{c|c}
B_{11}  & \\
\hline
        & B_{22} \\
\end{array}
\right)
\hspace{3cm}
\f]
where \f$ B_{11} \f$ and \f$ B_{22} \f$ are the blocks
formed from the corresponding "dof blocks" (the \f$
J_{ij} \f$ matrices).
The application of this preconditioner (i.e. the solution of the
linear system \f${\bf P} {\bf
  y} ={\bf z}\f$ for \f${\bf y}\f$) requires the solution of the 
two smaller linear systems
\f[
\hspace{3cm}
\left(
\begin{array}{cc}
J_{11} & J_{12} \\
J_{21} & J_{22} \\
\end{array}
\right) \left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right)
\mbox{\ \ \ \ \ \ \ or \ \ \ \ \ \ }
B_{11} \ Y_1 = Z_1
\hspace{3cm}
(8)
\f]
and
\f[
\hspace{3cm}
\left(
\begin{array}{ccc}
J_{33} & J_{34} & J_{35} \\
J_{43} & J_{44} & J_{45} \\
J_{53} & J_{54} & J_{55} \\
\end{array}
\right) \left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right)
\mbox{\ \ \ \ \ \ \ or \ \ \ \ \ \ }
B_{22} \ Y_2 = Z_2.
\hspace{3cm}
(9)
\f]

A key feature of the block preconditioning framework is the ability
to combine dof types in this manner so that the preconditioner can 
operate directly with blocks \f$ B_{11} \f$ and \f$ B_{22} \f$ and
the corresponding block vectors  \f$ Y_1, Y_2, Z_1 \f$ and \f$ Z_2 \f$.

Assuming again that the linear systems in
(8) and (9) are
solved exactly by a direct solver (an ``exact preconditioner'') that 
can pre-compute and store the LU decomposition of the diagonal matrix blocks,
\f$ B_{11} \f$ and \f$ B_{22} \f$, 
the \c setup() phase involves the following operations
[text in square brackets refers to their \c oomph-lib specific 
implementation]:
- Set up any data structures/lookup tables that are required to 
extract the matrix blocks \f$ B_{11} \f$ and \f$ B_{22} \f$ and
the associated block vectors [by calling
the \c BlockPreconditioner::block\_setup(...) function -- this time 
with arguments that specify the mapping between "dof types" and 
"block types"]. 
- Extract the two diagonal blocks, \f$ B_{11} \f$ and \f$ B_{22}\f$  [using the 
\c BlockPreconditioner::get\_block(...) function].
- Compute and store the LU decomposition of the diagonal blocks
to allow the rapid solution of the systems during the
 \c preconditioner\_solve(...)
phase by back-substitution. [This is done by calling the \c setup(...)
function of the subsidiary preconditioner/inexact solver. Following this,
the diagonal matrix blocks are longer required and can be deleted.]

Once the \c setup() phase has been completed, the solution of the linear
system \f${\bf P} {\bf y} = {\bf z}\f$ by the \c preconditioner\_solve(...) 
function involves the following steps:
- Extract the two ``block vectors'' \f${\bf Z}_i\f$ (for \f$i=1,2\f$)
from the vector \f${\bf z}\f$ [using the 
\c BlockPreconditioner::get\_block\_vectors(...) function].
- Solve the linear systems \f${\bf B}_{ii} \ {\bf Y}_i = 
{\bf Z}_i\f$ for the vectors \f${\bf Y}_i\f$ (for \f$i=1,2\f$)  using 
the precomputed LU decomposition of the diagonal blocks \f${\bf B}_{ii}\f$ 
(for \f$i=1,2\f$) created during the \c setup() phase.
- Combine the two ``block vectors'' \f${\bf Y}_i\f$ (for \f$i=1,...,2\f$)
to the full-length vector \f${\bf y}\f$ [using the 
\c BlockPreconditioner::return\_block\_vectors(...) function].

\subsubsection two_plus_three_implementation Implementation as a BlockPreconditioner

The implementation of the preconditioner closely follows that of the 
block diagonal preconditioner discussed above, the main difference
being that the current preconditioner only ever operates 
with exactly two blocks. Therefore we store pointers to the two 
subsidiary preconditioners (rather than a vector of pointers 
that can store an arbitrary number of these).

\skip start_of_two_plus_three_class
\until };

\subsubsection two_plus_three_setup The setup() function

As usual, we start by freeing up any previously allocated memory,
and set the pointer to the mesh:

\skip start_of_setup_for_two_plus_three
\until set_mesh


Since this preconditioner assumes explicitly that the problem involves
five dof types we check that this is actually the case.

\until #endif

To indicate that several dof types are to be combined into 
single blocks, we specify the mapping between dof types and block types
as an argument to the \c block_setup(...) function
This is done by creating vector of length \c ndof_type()
in which each entry indicates the block that the corresponding
dof is supposed to end up in:

\skip //
\until this->block_setup(dof_to_block_map);

To show that this actually worked, we output the number of blocks
(which should be -- and indeed is -- equal to two).

\skip Show that it 
\until types

Next we create the two subsidiary preconditioners and call their
\c setup(...) functions, passing the two diagonal blocks
 \f${\bf B}_{11} \f$  and \f${\bf B}_{22} \f$  to them.

\skip // Create the subsidiary preconditioners
\until // End of setup

\subsubsection two_plus_three_solve The preconditioner_solve() function

The \c preconditioner_solve(...) function is equivalent to that
in the \c Diagonal preconditioner discussed above, though here it simply
acts on a 2x2 block system. 



\skip //======
\until }

\subsubsection two_plus_three_clean The clean_up_my_memory() function

This function again deletes the allocated storage --
here the subsidiary preconditioners.

\skip //
\until // End of clean_up_my_memory function.





<HR>
<HR>




\subsection two_plus_three_upper_triangular Combining multiple dof types into compound blocks. Part 2: How to deal with off-diagonal blocks



NEW FEATURES: How to set up matrix vector products when multiple dof types
have been combined into compound blocks

\subsubsection two_plus_three_upper_triangular_theory Theory
The extension of the preconditioner introduced in the
previous section to block-triangular form is
straightforward: We use the same dof-to-block mapping as before but 
now retain the off-diagonal block \f$ B_{12} \f$ so that the 
preconditioner has the structure:
\f[
\hspace{3cm}
{\bf P} = 
\left(
\begin{array}{cc|ccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
J_{21} & J_{22} & J_{13} & J_{14} & J_{15} \\
\hline
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{43} & J_{44} & J_{45} \\
 & & J_{53} & J_{54} & J_{55} \\
\end{array}
\right)
=
\left(
\begin{array}{c|c}
B_{11}  & B_{12}\\
\hline
        & B_{22} \\
\end{array}
\right).
\hspace{3cm}
(10)
\f]

In the FSI context where \f$ B_{11} \f$  and \f$ B_{22} \f$ 
represent the solid and fluid sub-blocks, respectively, the inclusion
of the off-diagonal block \f$ B_{12} \f$ incorporates the effect of
fluid dofs (via pressure and shear stress) onto the solid equations.
Since this captures "more of the physics" the preconditioner can be
expected to be better than its block diagonal counterpart.

The application of the preconditioner (i.e. the solution of the
linear system \f${\bf P} {\bf
  y} ={\bf z}\f$ for \f${\bf y}\f$) requires the solution of the 
two smaller linear systems
\f[
\ \ \ \ \
\left(
\begin{array}{cc}
J_{11} & J_{12} \\
J_{21} & J_{22} \\
\end{array}
\right) \left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right) - \left(
\begin{array}{ccc}
J_{13} & J_{14} & J_{15} \\
J_{23} & J_{24} & J_{25}  \\
\end{array}
\right) \left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right)
\mbox{\ \ \ \ \ \ \ \ \ or \ \ \ \ \ \  \ \ }
B_{11} Y_1 =  Z_1 - B_{12} Y_2
\hspace{3cm}
(11)
\f]
and 
\f[
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 
\left(
\begin{array}{ccc}
J_{33} & J_{34} & J_{35} \\
J_{43} & J_{44} & J_{45} \\
J_{53} & J_{54} & J_{55} \\
\end{array}
\right) \left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right)
\mbox{\ \ \ \ \ \ \ \ \ or \ \ \ \ \ \  \ \ }
B_{22} Y_2 =  Z_2.
\hspace{3cm}
(12)
\f]

\subsubsection two_plus_three_upper_triangular_implementation Implementation as a BlockPreconditioner

The implementation is very similar to that in the previous example --
we simply provide additional storage for the (single) matrix vector
product operator required for the multiplication with \f$ B_{12} \f$
when updating the right-hand-side in equation
(11).

\skip start_of_two_plus_three_upper_triangular_class
\until }; 

\subsubsection two_plus_three_upper_triangular_setup The setup() function

As before, we start by freeing up any previously allocated memory
and set the pointer to the mesh,

\skip start_of_setup_for_two_plus_three_upper_triangular
\until set_mesh

and check that the number of dof types is correct.
\until #endif

The block setup is again performed with a dof-to-block mapping
that results in a block preconditioner with 2x2 blocks.

\skip //
\until this->block_setup(dof_to_block_map);

We create the two subsidiary preconditioners and pass the 
two diagonal blocks \f$ B_{11} \f$ and  \f$ B_{22} \f$
to their \c setup() functions. As before, the deep copies
of these matrices are then allowed to go out of scope, freeing
up the memory, since the subsidiary preconditioners retain whatever
information they require.

\skip // Create the subsidiary preconditioners
\until end setup of last

Finally we create and set up the off-diagonal vector product. Note
that the block column index refers to the block enumeration, 
so the block column index of \f$ B_{12} \f$ is 1 (in a C++ zero-based
enumeration!). 

\skip //
\until }


\subsubsection two_plus_three_upper_triangular_solve The preconditioner_solve() function

The application of the preconditioner is the exact equivalent of that
of the general-purpose block triangular preconditioner discussed
above, restricted to a 2x2 system:

\skip //======
\until }

\subsubsection two_plus_three_upper_triangular_clean The clean_up_my_memory() function

As before, this function frees up any memory that has been allocated
in the \c setup() function.

\skip //
\until // End of clean_up_my_memory function.


<HR>
<HR>


\subsection two_plus_three_upper_triangular_with_sub Using subsidiary block preconditioners

NEW FEATURES: How to use subsidiary block preconditioners to
(approximately) solve linear systems constructed from subsets of
dof-blocks.
 
\subsubsection two_plus_three_upper_triangular_with_sub_theory Theory

The two previous examples were motivated by the observation that in 
multi-physics problems (such as fluid-structure interaction) it is
natural to combine "related" dof blocks into compound block
matrices. We showed that the block preconditioning framework makes
it easy to extract such matrices from the original
system matrix and demonstrated how to solve linear systems involving these 
matrices with separate subsidiary preconditioners.  
One problem with this approach is that, once a compound matrix
has been created (by the \c get_block(...) function), all 
information about its dof types is lost, making it impossible to 
employ block preconditioners as subsidiary preconditioners. 

We will now revisit the 2x2 block triangular preconditioner 
described in the previous example and demonstrate how to employ
subsidiary block preconditioners to (approximately) solve linear
systems involving matrices formed (formally) by compound matrices
that are constructed from multiple dof-level blocks. From a
mathematical point of view, the structure of
the preconditioner therefore remains unchanged and is given by
\f[
\hspace{3cm}
{\bf P} = 
\left(
\begin{array}{cc|ccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
J_{21} & J_{22} & J_{13} & J_{14} & J_{15} \\
\hline
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{43} & J_{44} & J_{45} \\
 & & J_{53} & J_{54} & J_{55} \\
\end{array}
\right)=
\left(
\begin{array}{c|c}
B_{11}  & B_{12}\\
\hline
        & B_{22} \\
\end{array}
\right).
\hspace{3cm}
(13)
\f]
We will continue to use a dof-to-block mapping to view this as
the 2x2 block matrix shown on the right. This makes it easy to 
extract the compound off-diagonal block  \f$ B_{12} \f$ from the
system matrix when setting up the matrix-vector product (as before).
The setup of the subsidiary block preconditioners used to (approximately)
solve the linear systems involving \f$ B_{11} \f$ and \f$ B_{22}\f$ 
is handled differently:
- When calling the subsidiary block preconditioner's 
  \c setup(...) function we pass a pointer to the 
  entire system matrix, i.e. the matrix containing, formally, 
  all the dof-level blocks in equation (3).
.
- We then turn the preconditioner into a subsidiary block
  preconditioner, using its member function 
  \c turn_into_subsidiary_block_preconditioner(...) whose arguments 
  specify which of the dof-level blocks in the current (master) 
  preconditioner are to be used by the subsidiary block preconditioner.
.
The subsidiary block preconditioner is thus given access to all
the information required to extract the relevant data directly from the
original system matrix (and any associated full-length vectors).
It is in fact a key design principle of the block preconditioning
framework that <b>subsidiary block preconditioners
are given access to the "full size" matrices and vectors, but only
operate on the subset of data that they are "in charge of".</b>

When employing subsidiary block preconditioners for the approximate 
solution of the two smaller linear systems
\f[
\begin{array}{c}
\underbrace{
\left(
\begin{array}{cc}
J_{11} & J_{12} \\
J_{21} & J_{22} \\
\end{array}
\right) 
}_{B_{11}}
\underbrace{
\left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right)
}_{Y_{1}}
 = 
\underbrace{
\left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right)
}_{{\bf Z}_1} - 
\underbrace{
\left(
\begin{array}{ccc}
J_{13} & J_{14} & J_{15} \\
J_{23} & J_{24} & J_{25}  \\
\end{array}
\right)
}_{B_{12}}
\underbrace{
\left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right)
}_{{\bf Y}_2}
\\
\hspace{3cm}
(14)
\\
\underbrace{
\left(
\begin{array}{ccc}
J_{33} & J_{34} & J_{35} \\
J_{43} & J_{44} & J_{45} \\
J_{53} & J_{54} & J_{55} \\
\end{array}
\right)
}_{B_{22}}
\underbrace{
\left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right)
}_{{\bf Y}_2}
 = 
\underbrace{
\left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right)
}_{{\bf Z}_2},
\end{array}
\f]
the subsidiary preconditioners that operate on 
the linear systems involving \f$ B_{11} \f$ and \f$ B_{22} \f$ 
therefore retain access to the 
relevant dof-level blocks. Hence, if we employ the block triangular 
preconditioner discussed above to (approximately) solve the
two linear systems in equation 
(14),
the complete preconditioning operation is described by the
following equations:
\f[
\begin{array}{c}
\left(
\begin{array}{cc}
J_{11} & J_{12} \\
 & J_{22} \\
\end{array}
\right) 
\left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right) 
= 
\underbrace{
\underbrace{
\left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right)
}_{{\bf Z}_1}
 - 
\underbrace{
\left(
\begin{array}{ccc}
J_{13} & J_{14} & J_{15} \\
J_{23} & J_{24} & J_{25}  \\
\end{array}
\right)
}_{B_{12}}
\underbrace{
\left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right)
}_{{\bf Y}_2}
}_{\widehat{\bf Z}_1}
\\
\left(
\begin{array}{ccc}
J_{33} & J_{34} & J_{35} \\
 & J_{44} & J_{45} \\
 &  & J_{55} \\
\end{array}
\right) \left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right).
\end{array}
\f]
Note that when we wrote the block triangular preconditioner 
we did not have to be aware of the fact that it may subsequently be used
as a subsidiary block preconditioner. The
internal data structures implemented in the \c BlockPreconditioner
base class ensure that when we call \c get_block(0,0,block_matrix)
in the subsidiary block preconditioner acting on \f$ B_{22} \f$, 
\c block_matrix will 
receive a deep copy of \f$ J_{33} \f$, extracted from the full system matrix. 
Similarly, a call to \c get_block_vectors(r,block_r) will
extract the three block vectors \f$ {\bf r}_3, {\bf r}_4 \f$
and \f${\bf r}_5 \f$ from the full-length vector \f$ {\bf r} \f$,
while \c return_block_vectors(block_z,z) will return the
three solution vectors \f$ {\bf z_3}, {\bf z_4} \f$ and 
\f$ {\bf z_5} \f$ to the appropriate entries in the full-length
vector \f$ {\bf z}. \f$


The implementation of the preconditioning operations
can again be subdivided into two
distinct \c setup() and \c preconditioner\_solve(...) phases.
- Set up the data structures/lookup tables that map dof types
  0 and 1 to block 0 and dof types 2, 3 and 4 to block 1 [by calling
  the \c BlockPreconditioner::block\_setup(...) function
  with arguments that specify the mapping between "dof types" and 
  "block types" as before]. 
- Create two instances of the block triangular preconditioner
  (or any other block preconditioner) and turn them into the 
  subsidiary preconditioners for the current (master) preconditioner,
  specifying which dof types in the master preconditioner the
  subsidiary block preconditioners are to work with.
- Extract the compound off-diagonal block \f$ B_{12}\f$  and create 
  a \c MatrixVectorProduct operator. 

Once the \c setup() phase has been completed, the solution of the linear
system \f${\bf P} {\bf y} = {\bf z}\f$ by the \c preconditioner\_solve(...) 
function involves the following steps:
- Solve the linear systems \f$ {\bf B}_{22} \ {\bf Y}_2 = 
{\bf Z}_2\f$ using the subsidiary block preconditioner that works
with  \f$ {\bf B}_{22} \f$. [The subsidiary block preconditioner's
\c preconditioner_solve(...) function is given access to the 
full-size vectors \f$ {\bf z} \f$ and \f$ {\bf y} \f$ and
extracts/returns \f$ {\bf Z}_2 \f$
and \f$ {\bf Y}_2\f$ directly from/into these.]
- Extract the solution vector \f$ {\bf Y}_2 \f$ from the just undated full-length
vector \f$ {\bf y} \f$, perform the matrix vector product with
\f$ B_{12} \f$ and store the result in a temporary vector \f$ {\bf t} \f$.
- Extract the block vector \f$ {\bf Z}_1\f$ from the full-length
vector \f$ {\bf z} \f$, subtract \f$ {\bf t} \f$
from it, and return the result, \f$ \widehat{\bf Z}_1 = 
{\bf Z}_1 - B_{12} {\bf Z}_2 \f$
into the appropriate entries into the
full-length vector \f$ {\bf z} \f$.
- Solve the linear systems \f$ {\bf B}_{11} \ {\bf Y}_1 = 
\widehat{\bf Z}_1 = {\bf Z}_1 - B_{12} {\bf Z}_2 \f$ using the subsidiary block 
preconditioner that works
with  \f$ {\bf B}_{11} \f$. [The subsidiary block preconditioner's
\c preconditioner_solve(...) function is given access to the 
full-size vectors \f$ {\bf z} \f$ and \f$ {\bf y} \f$ and
extracts/returns \f$ \widehat{\bf Z}_1 \f$
and \f$ {\bf Y}_1\f$ directly from/into these; recall that the 
relevant entries in \f$ {\bf z}\f$ have been over-written 
in the previous step so that \f$ \widehat{\bf Z}_1 \f$ contains 
the updated right hand side.]


\subsubsection two_plus_three_upper_triangular_with_sub_implementation Implementation as a BlockPreconditioner

The implementation of the preconditioner is completely equivalent to the
corresponding block triangular preconditioner considered in the
previous example:

\skip start_of_two_plus_three_upper_triangular_with_sub_class
\until }; 



\subsubsection two_plus_three_upper_triangular_with_sub_setup The setup() function
As usual we free up any memory and set the pointer to the mesh:

\skip start_of_setup_for_two_plus_three_upper_triangular
\until set_mesh

We check that the number of dof types is appropriate for this preconditioner:
\until #endif

Next we define the block structure of the preconditioner, using
a dof-to-block mapping to combine dofs 0 and 1 into block 0,
and dofs 2, 3 and 4 into block 1:

\skip //
\until this->block_setup(dof_to_block_map);
 
Next we create the block triangular preconditioner used
to (approximately) solve linear systems involving 
the compound "top left" 2x2 block:
\skip // Create the subsidiary block preconditioners
\until First_sub

Next we specify the pointer to the mesh that contains the elements
that classify the degrees of freedom. We note, that, strictly speaking
this is not necessary since the preconditioner will only be used
as a subsidiary preconditioner --  the enumeration of the dof types
is always handled by the top-most master preconditioner. One (or more)
mesh pointers must be set for the master preconditioner, and, if
compiled in \c PARANOID mode, \c oomph-lib will throw an error 
if this is not done. Some (but not all!) \c oomph-lib developers
regard it as "good practice" to set the mesh 
pointer anyway, so one is less likely to forget...

\skip Set mesh
\until set_multi_poisson_mesh


We turn this preconditioner into a subsidiary block preconditioner, specifying
the pointer to the current (master) preconditioner and the
mapping between dof types in the present and the subsidiary
block preconditioners (here the identity):
\until turn_into_subsid
 
When calling the subsidiary block preconditioners \c setup(...)
function we pass a pointer to the full matrix:
\until }


The second subsidiary block preconditioner (for the 3x3 "bottom right"
compound matrix) is created similarly,  though the mapping between
dof-types is now no longer the identity but maps dof types 2, 3 and 4
in the current (master) preconditioner to dof types 0, 1 and 2
in the subsidiary block preconditioner:
\until }

The setup of the matrix-vector product with the off-diagonal matrix 
is unchanged from the previous example:

\skip //
\until }
\until }

\subsubsection two_plus_three_upper_triangular_with_sub_solve The preconditioner_solve() function

As discussed in the theory section, we start by (approximately)
solving the system \f$ B_{22} Y_2 = Z_2 \f$, using the second
subsidiary block preconditioner which automatically 
extracts \f$ Z_2 \f$ from the full length vector \c z
and returns the result \f$ Y_2 \f$ into the appropriate entries
of the full length vector \c y.

\skip //
\until Second_subsidiary_preconditioner_pt->preconditioner_solve(z,y);

We now extract the block vector \f$ Y_2 \f$ from the full-length
vector \c y,
\until get_block_vector

multiply it by \f$ B_{12} \f$, using the \c MatrixVectorProduct
operator,
\until Off_diagonal

and subtract the result from \f$ Z_1 \f$ (which we extract from the
full length vector \c z ):

\skip //
\until block_z -= temp; 

\c block_z now contains the updated right hand side, \f$ \widehat{{\bf
Z}_1}\f$, for the linear system to be (approximately) solved by the
first subsidiary block preconditioner. We therefore return 
\f$ \widehat{{\bf Z}_1}\f$ to the appropriate entries into a full
length vector of the same size as right hand side vector \c z:

\until return_block_vector

We then pass this vector to first subsidiary preconditioner which 
updates the appropriate entries in the full-length solution vector \c
y which can therefore be returned directly by this function:

\skip //
\until }

\subsubsection two_plus_three_upper_triangular_with_sub_clean The clean_up_my_memory() function

As usual, we use this helper function to free up any memory 
allocated in the \c setup() function to avoid memory leaks.

\skip //
\until // End of clean_up_my_memory function.




<HR>
<HR>




\subsection two_plus_three_upper_triangular_with_replace Replacing/modifying blocks


NEW FEATURES: How to replace/modify matrix blocks


\subsubsection two_plus_three_upper_triangular_with_replace_theory Theory

So far, we have demonstrated how to extract matrix blocks from the
full-sized system matrix (typically the Jacobian matrix used in 
Newton's method) and how to apply a preconditioner via operations
involving these blocks. Many preconditioners do not operate directly
with the matrix blocks themselves, but on matrices that are derived from them.
For instance, <c>oomph-lib</c>'s  <a href="../../../preconditioners/lsc_navier_stokes/html/index.html">Schur complement Navier-Stokes preconditioner</a>
operates on an (approximate) Schur complement;
augmentation preconditioners involve
operations on matrices that are obtained by the addition of a diagonal
matrix to some of the matrix blocks; etc. Within a given
preconditioner such derived matrices are typically pre-computed 
by the preconditioner's \c setup() function and then stored as 
private member data which makes them available to the \c
preconditioner_solve() function. Unfortunately, this approach does not 
work if the modified block is to be used in a subsidiary block preconditioner
because, as discussed in the previous example, 
by default the subsidiary block preconditioner will extract its block
matrices directly from the full-size system matrix and will therefore
ignore any (local) modifications made by its master preconditioner(s). 
What is therefore required is a method that indicates to the block 
preconditioning framework that a given sub-block is not to be 
extracted from the full system matrix but to be represented by suitable 
replacement matrix.

We demonstrate this methodology by re-visiting the preconditioner
considered in the previous example, namely
\f[
\hspace{3cm}
{\bf P}_{\rm previous} = 
\left(
\begin{array}{cc|ccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
J_{21} & J_{22} & J_{13} & J_{14} & J_{15} \\
\hline
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{43} & J_{44} & J_{45} \\
 & & J_{53} & J_{54} & J_{55} \\
\end{array}
\right).
\hspace{3cm}
(15)
\f]
However, here we want to modify the off-diagonal blocks by
"replacing" each block \f$J_{ij}\f$ (for \f$ i \ne j \f$) by 
a "replacement matrix" \f$R_{ij}\f$ so that the preconditioner becomes
\f[
{\bf P} = 
\left(
\begin{array}{cc|ccc}
J_{11} & R_{12} & R_{13} & R_{14} & R_{15} \\
R_{21} & J_{22} & R_{23} & R_{24} & R_{25} \\
\hline
 & & J_{33} & R_{34} & R_{35} \\
 & & R_{43} & J_{44} & R_{45} \\
 & & R_{53} & R_{54} & J_{55} \\
\end{array}
\right) =
\left(
\begin{array}{c|c}
B_{11}  & B_{12}\\
\hline
        & B_{22} \\
\end{array}
\right).
\hspace{3cm}
(16)
\f]
The application of this preconditioner (i.e. the solution of the
linear system \f${\bf P} {\bf
  y} ={\bf z}\f$ for \f${\bf y}\f$) still requires the solution of the 
two smaller linear systems
\f[
\begin{array}{c}
\underbrace{
\left(
\begin{array}{cc}
J_{11} & R_{12} \\
R_{21} & J_{22} \\
\end{array}
\right)
}_{B_{11}}
\left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right) - 
\underbrace{
\left(
\begin{array}{ccc}
R_{13} & R_{14} & R_{15} \\
R_{23} & R_{24} & R_{25} \\
\end{array}
\right)
}_{B_{12}}
\left(
\begin{array}{c}
z_3 \\
z_4 \\
z_5 \\
\end{array}
\right)
\\
\hspace{3cm}
(17)
\\
\underbrace{
\left(
\begin{array}{ccc}
J_{33} & R_{34} & R_{35} \\
R_{43} & J_{44} & R_{45} \\
R_{53} & R_{54} & J_{55} \\
\end{array}
\right) 
}_{B_{22}}
\left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right)
\end{array}
\f]
where we have again assumed that the two vectors \f${\bf y}\f$ and \f${\bf z}\f$
are re-ordered into ``block vectors'' in the same way as 
the vectors \f$\delta {\bf x}\f$ and \f${\bf r}\f$ in 
"the original linear system" (3) are re-ordered into
the ``block vectors'' in
(17).
We wish to continue to solve the linear systems involving
the compound matrices \f$ B_{11}  \f$ and \f$ B_{22} \f$ 
(which involve "replaced" blocks) by two subsidiary block 
preconditioners (which operate on 3x3 and 2x2 dof blocks, 
respectively).


In the specific example below we replace all of the diagonal 
matrices by suitably sized zero matrices, so that the actual 
preconditioning operation is defined by the following
linear systems
\f[
\begin{array}{c}
\underbrace{
\left(
\begin{array}{cc}
J_{11} &  \\
& J_{22} \\
\end{array}
\right) 
}_{B_{11}}
\left(
\begin{array}{c}
y_1 \\
y_2 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_1 \\ 
z_2 \\
\end{array}
\right) - 
\underbrace{
\left(
\begin{array}{ccc}
& &  \\
 & &  \\
\end{array}
\right) 
}_{B_{12}}
\left(
\begin{array}{c}
z_3 \\
z_4 \\
z_5 \\
\end{array}
\right)
\\
\underbrace{
\left(
\begin{array}{ccc}
J_{33} & &  \\
& J_{44} &  \\
 &  & J_{55} \\
\end{array}
\right) 
}_{B_{22}}
\left(
\begin{array}{c}
y_3 \\
y_4 \\
y_5 \\
\end{array}
\right) = \left(
\begin{array}{c}
z_3 \\ 
z_4 \\
z_5 \\
\end{array}
\right)
\end{array}
\f]
which, in effect, turns the preconditioner into the 
block-diagonal preconditioner considered at the very
beginning of this tutorial.

\subsubsection two_plus_three_upper_triangular_with_replace_implementation Implementation as a BlockPreconditioner

The implementation of the preconditioner is completely equivalent to the
preconditioner considered in the previous example. The only additional
feature is the provision a matrix of pointers to the
replacement matrices, \c Replacement_matrix_pt.

\skip start_of_two_plus_three_upper_triangular_with_replace_class
\until }; 



\subsubsection two_plus_three_upper_triangular_with_replace_setup The setup() function

As usual, we start by cleaning up any memory using a call to a
\c clean_up_my_memory() function, and set the pointer to the mesh

\skip start_of_setup_for_two_plus_three_upper_triangular_with_replace
\until set_mesh

Next we check that the number of dof types is 5, as the
preconditioner is designed to only work for that number. 

\skip How many dof types do we have
\until #endif

The block setup follows exactly the same pattern as in the previous
example: Dof types 0 and 1 are combined into compound block 0, while
dof types 2, 3 and 4 are combined into compound block 1.  
On return from the block setup function we should therefore have
two block types:
\skip Call block setup
\until #endif

Now we perform the replacement of the off-diagonal dof blocks.
(Note that there are still five of these. Dof-blocks and compound
blocks are not the same -- if you get them confused you will get into
trouble!). We allocate storage for the pointers to the replacement 
matrices and loop over the off-diagonal blocks:

\skip //
\until (i!=j)
\until {

  
Given that the replacement matrices are zero matrices, we could
simply create them without ever looking at the original blocks.
Sadly the creation of zero matrices turns out to be slightly
more painful than one would wish because they have
to be created as a (possibly distributed) \c CRDoubleMatrix. The 
relevant code is contained in the source code but we won't discuss 
it here since the more common situation is one where we actually
want to modify the already existing entries of an already 
existing block matrix. 
Therefore we simply extract the matrix and set its initially nonzero entries
to zero (admittedly a bit silly -- we now have a sparse matrix full of
zeroes, but it's just a demonstration!):

\until this is just an example!

We then pass the pointer to the replacement dof block to the
block preconditioner

\skip Replace (i,j)
\until loop of i

The rest of the setup works exactly as in the previous example, only
this time, the subsidiary preconditioners and the matrix vector 
products will work with the replacement dof blocks that we've just
defined.

We create and set up the first subsidiary block preconditioner 
which operates on our dof types 0 and 1 (and treats them as its own
dof  types 0 and 1):
\skip First subsidiary precond
\until }

The second  subsidiary block preconditioner 
which operates on our dof types 2, 3 and 4 (and treats them as its own
dof  types 0, 1 and 2):

\skip Second subsidiary precond 
\until }

Finally, we create the matrix vector product operator:
\skip Next setup the off diagonal mat vec operators:
\until }
\until }



\subsubsection two_plus_three_upper_triangular_with_replace_solve The preconditioner_solve() function

The \c preconditioner_solve() function is completely identical to the
one used in the previous preconditioner, so we omit the code
listing -- the subsidiary preconditioners and the matrix vector 
product operator work in the same way but now simply operate on 
the replacement dof blocks where they have been set.

\subsubsection two_plus_three_upper_triangular_with_replace_clean The clean_up_my_memory() function

Memory is cleaned up as before, so we omit the code listing.

<HR>
<HR>


\subsection coarse_two_plus_two_plus_one Coarsening/combining dof types

NEW FEATURES: How to coarsen/combine dof types for use by subsidiary
block preconditioners.

\subsubsection coarse_two_plus_two_plus_one_theory Theory

In the examples presented so far we have demonstrated how
to combine various dof-blocks into compound blocks in order to 
facilitate the application of certain preconditioning operations.
For instance, in many of the previous examples we performed a matrix
vector product using the compound matrix \f$ B_{12} \f$ that was
(formally) formed by the concatenation of the 2x3 "top right"
off-diagonal dof blocks in the full-sized system.

We also showed how subsidiary block preconditioners which operate on
a specific number of dof blocks can be instructed to operate on 
selected dof types from the full-sized system. Our standard
example for this was a 2D Navier-Stokes preconditioner which operates on
three dof types (two fluid velocities and one pressure) and 
is used as a subsidiary block preconditioner in an FSI problem that 
also involves additional dofs associated with the solid mechanics
(e.g. the two solid displacement components). This was done
by informing the subsidiary preconditioner which of the dof types
in the full-sized system to regard as "its own" when calling
its \c turn_into_subsidiary_block_preconditioner(...) function.
This implies that the subsidiary block preconditioner remains unaware of any 
 compound blocks that may have been formed in its master 
 preconditioner. The functionality presented so far only allows us
to associate dof-blocks in the master preconditioner with dof
blocks in the subsidiary block preconditioner.
It is therefore not possible (without further functionality which we explain
in this example) to use a subsidiary block preconditioner if the dof-types
in the master preconditioner are "too fine-grained". This arises,
for instance, in Navier-Stokes problems where the master preconditioner
sub-divides the two components of the fluid velocity into degrees of freedom
on the domain boundary and those in the interior. It is then necessary
to make the subsidiary preconditioner act on the combined 
dof types, a process that we describe as "coarsening".



We illustrate the procedure by returning, yet again, to our 5x5 block
linear system that we wish to precondition with
\f[
\hspace{3cm}
{\bf P}_{\rm previous} = 
\left(
\begin{array}{cc|ccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
J_{21} & J_{22} & J_{13} & J_{14} & J_{15} \\
\hline
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{43} & J_{44} & J_{45} \\
 & & J_{53} & J_{54} & J_{55} \\
\end{array}
\right) =
\left(
\begin{array}{c|c}
B_{11}  & B_{12}\\
\hline
        & B_{22} \\
\end{array}
\right).
\hspace{3cm}
(18)
\f]
However, now we wish to solve the two linear systems involving the
compound matrices \f$ B_{11} \f$ and \f$ B_{22} \f$ with a 
2x2 upper triangular subsidiary block
preconditioner. To make this possible, we "coarsen" the dof types 
such that the subsidiary  block preconditioner acting on 
\f$ B_{22}\f$ treats the global dof types 3 and 4 as a single dof type
so that the block structure can be viewed as
\f[
{\bf P} = 
\left(
\begin{array}{c||c}
\begin{array}{cc}
J_{11} & J_{12} \\
J_{21} & J_{22}  \\
\end{array} &
\begin{array}{ccc}
J_{13} & J_{14} & J_{15} \\
J_{23} & J_{24} & J_{25} \\
\end{array} \\
\hline
\hline
&
\begin{array}{cc|c}
J_{33} & J_{34} & J_{35} \\
J_{43} & J_{44} & J_{45} \\
\hline
J_{53} & J_{54} & J_{55} \\
\end{array}
\end{array}
\right) =
\left(
\begin{array}{c|c}
B_{11}  & B_{12}\\
\hline
        & B_{22} \\
\end{array}
\right).
\hspace{3cm}
(19)
\f]

If we now use a 2x2 upper triangular block preconditioner 
to (approximately) solve the linear systems involving the
diagonal blocks \f$ B_{11} \f$ and \f$ B_{22} \f$ 
the preconditioner is given (mathematically) by
\f[
{\bf P} = 
\left(
\begin{array}{ccccc}
J_{11} & J_{12} & J_{13} & J_{14} & J_{15} \\
       & J_{22} & J_{23} & J_{24} & J_{25} \\
 & & J_{33} & J_{34} & J_{35} \\
 & & J_{34} & J_{44} & J_{45} \\
 & & & & J_{55} \\
\end{array}
\right).
\f]

[Note that, In the actual implementation discussed below, we also set
the off diagonal dof-blocks to zero, using the replacement methodology
discussed in the previous example. The preconditioner therefore
becomes mathematically equivalent to the \f$ 5x5  \f$ block 
diagonal preconditioner discussed at the very beginning of this tutorial.]

\subsubsection coarse_two_plus_two_plus_one_implementation Implementation as a BlockPreconditioner

The implementation of the preconditioner is completely equivalent to the
preconditioner considered in the previous examples:

\skip start_of_coarse_two_plus_two_plus_one_class
\until }; 

\subsubsection coarse_two_plus_two_plus_one_setup The setup() function

As usual we clean up any previously allocated memory and set the pointer to the
mesh:

\skip start_of_setup_for_coarse_two_plus_two_plus_one
\until set_mesh

Next we check that the number of degrees of freedom is 5, as the
preconditioner is designed to only work for that number.


\skip //
\until #endif

The block setup follows exactly the same pattern as in the previous 
examples: Dof types 0 and 1 are combined into compound block 0, 
while dof types 2, 3 and 4 are combined into compound block 1.


\skip //
\until this->block_setup(dof_to_block_map);

[We omit the code listing the replacement of the off-diagonal dof
blocks with zero matrices since it is identical to what we already
discussed in the previous example.]

Next we create the two subsidiary preconditioners that (approximately)
solve the linear systems involving the diagonal blocks \f$ B_{11} \f$
and \f$ B_{22} \f$. The first subsidiary preconditioner is a
standard upper triangular block preconditioner which acts on the
compound block formed by dof types 0 and 1:

\skip // Create the subsidiary
\until }

The second subsidiary preconditioner is more interesting. It's a 
block preconditioner that only operates on a 2x2 block system, yet
we want to use it to solve the linear system involving the compound
block formed the three dof types 2, 3 and 4. To do this we wish
to combine the dof blocks associated with dof types 2 and 3 into 
a single block. We start by setting the mesh pointer and by
setting up the usual mapping that identifies the dof types (in the 
current preconditioner) that we wish the 
subsidiary preconditioner to act on.

\skip // Second 
\until  dof_map[2]

To combine/coarsen dof types 2 and 3 (in the current preconditioner) 
into a single dof type for the subsidiary preconditioner we create 
a vector of vectors, <code> doftype_coarsening </code> whose
entries are to be interpreted as

<code> doftype_coarsening[coarsened_dof_type][i]=dof_type</code>

where \c i ranges from 0 to the number of dof types (minus one, because
of the zero-based indexing...) in the enumeration of the
subsidiary preconditioner that are to be
combined/coarsened into dof type \c coarsened_dof_type:

\skip The subsidiary block preconditioner
\until doftype_coarsening[1][0]=2

We pass both lookup schemes to the function that turns the preconditioner
into a subsidiary block preconditioner and then call its own setup
function, as usual.

\until }

Finally, we set up the of diagonal matrix-vector product which acts on the
compound (0,1) block (formed from dof types {0,1}x{2,3,4}) in the current
preconditioner. 

\until End of setup


\subsubsection coarse_two_plus_two_plus_one__solve The preconditioner_solve() function

The \c preconditioner_solve() function is completely identical to the
one used in the previous example, so we omit the code
listing.

\subsubsection coarse_two_plus_two_plus_one_clean The clean_up_my_memory() function

Memory is cleaned up as before, so we omit the code listing.


<HR>
<HR>


\subsection fsi_multiple_meshes Using multiple meshes -- explained for a genuine fluid-structure interaction problem

NEW FEATURES: How to use multiple meshes
 
\subsubsection fsi_multiple_meshes_theory Theory
Finally, we demonstrate the use of multiple meshes by discussing a simple 
implementation of the FSI preconditioner described in the 
<a href="../../../preconditioners/fsi/html/index.html">FSI Preconditioner
Tutorial</A>. We refer to the tutorial discussing the 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">FSI
channel with leaflet problem</A> for the overall problem setup.


FSI problems involve fluid (velocities and pressures from the 
Navier-Stokes equations) and solid (the nodal positions in the solid
domain) degrees of freedom (dofs). We begin by reordering the linear 
system to group together the two types of dof 

\f[
\left[ 
\begin{array}{cc}
F&C_{fs}\\
C_{sf}&S
\end{array}
\right]
\left[ 
\begin{array}{c}
\bf \delta f\\
\bf \delta s
\end{array}
\right]
=
-
\left[
\begin{array}{c}
\bf r_f\\
\bf r_s
\end{array}
\right] ,
\f]


where \f$\bf f\f$ and \f$\bf s\f$ denote the fluid and solid dofs, \f$F\f$
is the Navier-Stokes Jacobian (representing the derivatives of the 
discretised fluid equations with respect to the fluid dofs), \f$S\f$ 
is the solid Jacobian, and the
blocks \f$C_{fs}\f$ and \f$C_{sf}\f$ arise from the interaction
between fluid and solid equations.

The Navier Stokes Jacobian \f$F\f$ has its own block
structure. Decomposing the fluid dofs into velocity and pressure dofs so that

\f[
{\bf f}=
\left[
\begin{array}{c}
\bf u\\
\bf p
\end{array}
\right],
\f]

we obtain the well known saddle-point structure of \f$F\f$

\f[
F=
\left[
\begin{array}{cc}
A&B^T\\
B& 
\end{array}
\right],
\f]

where \f$A\f$ is the momentum block, \f$B^T\f$ the discrete gradient
operator, and \f$B\f$ the discrete divergence operator (see 
<a href="../../../preconditioners/lsc_navier_stokes/html/index.html">Navier
Stokes Preconditioner Tutorial</A>). 

This FSI preconditioner 
takes the form of a block triangular preconditioner. Here we
only consider the lower block triangular version 
\f[
P_{FSI}=
\left[
\begin{array}{cc}
F& \\
C_{sf}&S
\end{array}
\right]
\f]
obtained by omitting the \f$C_{fs}\f$ block from the Jacobian.


The application of the preconditioner requires the solution of the
linear 
system

\f[
\left[
\begin{array}{cc}
F& \\
C_{sf}&S
\end{array}
\right]
\left[
\begin{array}{c}
\bf z_f\\
\bf z_s
\end{array}
\right]
=
\left[
\begin{array}{c}
\bf y_f\\
\bf y_s
\end{array}
\right].
\f]

However, for preconditioning purposes this system does not have to be
solved exactly. We therefore replace the solution of the linear
systems involving the diagonal blocks (representing the single-physics
fluid and solid Jacobians \f$F\f$ and \f$S\f$) by existing
preconditioners (interpreted as inexact solvers). Formally, we write
this as

\f[
\left[
\begin{array}{cc}
\tilde F& \\
C_{sf}&\tilde S
\end{array}
\right]
\left[
\begin{array}{c}
\bf z_f\\
\bf z_s
\end{array}
\right]
=
\left[
\begin{array}{c}
\bf y_f\\
\bf y_s
\end{array}
\right] . \ \ \ \ \ \ \ \ \ \ \ \  (20)
\f]

where  \f$ \tilde F \f$ is the fluid preconditioner and \f$ \tilde{S} \f$
the solid preconditioner, both used as 
subsidiary preconditioners. 

The application of the preconditioner can be accomplished in four
distinct steps:

-# Apply the fluid preconditioner \f$ \tilde F \f$ to the fluid dofs of
   the RHS vector \f$ \bf y_f \f$ and store the result in the fluid
   solution \f$ {\bf z_f}=\tilde F^{-1}{\bf y_f} \f$ .
-# Multiply the fluid-solid coupling matrix \f$ C_{sf} \f$ with the
   fluid solution \f$ \bf z_f \f$ and store the result in the temporary
   vector \f$ {\bf w}=C_{sf}\bf z_f \f$ .
-# Subtract \f$ \bf w \f$ from the solid dofs of the RHS vector 
   \f$ \bf y_s \f$ and store the result in the temporary \f$ \bf w \f$ to
   complete the action of the \f$ C_{sf} \f$ matrix vector product, 
   \f$ {\bf w}={\bf y_s}-{\bf w} \f$ .
-# Apply the solid preconditioner \f$\tilde S\f$ to the temporary 
   \f$ \bf w \f$ to compute the solid solution 
   \f$ {\bf z_s}=\tilde S^{-1}{\bf w} \f$ .
.
This is, of course, extremely similar to the methodology explained in 
the section \ref two_plus_three_upper_triangular_with_sub, the main
difference being that the fluid and solid dofs are classified
by two different elements. In the two-dimensional 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
FSI channel with leaflet problem</A> these are:
- The fluid elements are of type \c
   RefineableQTaylorHoodElement<2>. These  elements have 
  three types of dof; \f$x\f$-velocity dofs are labelled \c 0, 
  \f$y\f$-velocity dofs are labelled \c 1 and
  the pressure dofs are labelled \c 2.\n\n
- The solid elements are of type \c FSIHermiteBeamElement.  
  They have  one type of dof (the nodal position) labelled \c 0.
.
When classifying the dofs we specify the elements via two separate
meshes, the first one containing the pointers to the fluid elements,
the second one the pointers to the solid elements. This means that
in the global enumeration of the dof types the fluid dofs appear
before the solid dofs.

  
\subsubsection implementing_the_fsi_preconditioner The Implementation of the FSI Preconditioner
We implement the FSI preconditioner in the
class \c SimpleFSIPreconditioner. This class inherits
from the base class \c BlockPreconditioner which provides the 
generic functionality required for common block preconditioning operations.

The overall structure of the class is similar to that of the
preconditioners considered before, the main difference being that
we now store pointers to two meshes.

\dontinclude fsi_channel_with_leaflet.cc
\skipline start_of_simple_fsi_preconditioner
\until };

\subsubsection setup Preconditioner Setup

We start by setting up the meshes,
choosing the fluid mesh to be
mesh \c 0 and the solid mesh to be mesh \c 1. The preconditioner therefore
has four dof types enumerated in mesh order:

- \c 0 fluid \f$x\f$ velocity (dof type 0 in mesh 0)
- \c 1 fluid \f$y\f$ velocity (dof type 1 in mesh 0)
- \c 2 fluid pressure (dof type 2 in mesh 0)
- \c 3 solid (dof type 0 in mesh 1)

\dontinclude fsi_channel_with_leaflet.cc
\skipline start_of_setup
\until n_fluid_dof_type + this->ndof_types_in_mesh(1);


Next we define the mapping from dof number to block
number. The preconditioner has two block types -- fluid and solid --
therefore we group the fluid dofs into block type \c 0 and the solid dofs
into block type \c 1. We define a map from dof type to block
type in a vector (the vector indices denote the dof type and the
vector elements denote the block type) and pass it to \c
block_setup(...) to complete the setup of the \c BlockPreconditioner 
infrastructure.

\dontinclude fsi_channel_with_leaflet.cc
\skipline This fsi preconditioner has two types of block
\until this->block_setup(dof_to_block_map);

Next we set up the subsidiary
operators required by the preconditioner. 
We start with the solid subsidiary preconditioner (\f$\tilde S\f$). We  extract
the solid matrix block \f$S\f$ from the Jacobian using the \c
BlockPreconditioner method \c get_block(...) and then set up the solid
subsidiary preconditioner:

\dontinclude fsi_channel_with_leaflet.cc
\skipline First the solid preconditioner
\until delete solid_matrix_pt;

Note that, compared to the previous examples, we have used an 
alternative, pointer-based version of the \c get_block(...) function.
However, as before, the block matrix can be deleted once the 
subsidiary preconditioner has been set up since the latter
retains whatever data it requires. 

The fluid subsidiary preconditioner (\f$\tilde F\f$) a block preconditioner
itself. Its setup is therefore performed in two steps:

-# First we turn the \c NavierStokesSchurComplementPreconditioner
   into a subsidiary block preconditioner. We assemble a list a fluid
   dof types in the current (master) preconditioner, and pass this list to the
   Navier-Stokes preconditioner to indicate that dof type \c i in the
   master FSI preconditioner is dof type \c i in the subsidiary
   fluid preconditioner (for \c i \c = \c 0, \c 1, \c 2) (Note that
   the fact that this mapping is the identity mapping is a result of
   choosing the fluid mesh to be mesh \c 0; in general the index of \c
   ns_dof_list corresponds to the dof type number in the Navier Stokes
   subsidiary preconditioner and the value corresponds to the index in
   this master preconditioner).
\dontinclude fsi_channel_with_leaflet.cc
\skipline Next the fluid preconditioner
\until turn_into_subsidiary_block_preconditioner(this,ns_dof_list);
-# Next we set up the \c NavierStokesSchurComplementPreconditioner.  We
   pass the Navier-Stokes mesh to the the subsidiary preconditioner
   and set up the preconditioner. Note that the pointer to the full FSI
   Jacobian is passed to the subsidiary block preconditioner. This
   allows the subsidiary preconditioner to extract the relevant
   sub-blocks, using the lookup schemes established by the call to \c
   turn_into_subsidiary_block_preconditioner(...).
\dontinclude fsi_channel_with_leaflet.cc
\skipline Set up the NavierStokesSchurComplement preconditioner.
\until Navier_stokes_preconditioner_pt->setup(this->matrix_pt());

Finally, we set up is the matrix-vector product. This mirrors
the set up of the solid subsidiary preconditioner. First the subsidiary
matrix is extracted from the Jacobian and then the operator is set up:

\dontinclude fsi_channel_with_leaflet.cc
\skipline Finally the fluid onto solid matrix vector product operator
\until delete fluid_onto_solid_matrix_pt;

Again, the extracted block can be deleted since the matrix vector
product operator retains the relevant data.
The FSI preconditioner is now ready to be used. 

\subsection solve Preconditioner Solve


The \c preconditioner_solve(...) method applies  
the preconditioner to the input vector \f$\bf y\f$ and returns the 
result in \f$\bf z\f$.


We start by applying the Navier-Stokes preconditioner \f$\tilde
F\f$ to the fluid elements \c y_f of \c y. Since \f$\tilde F\f$ is a 
subsidiary block preconditioner  we apply it to the full-length
 \c y and \c z vectors which contain both the fluid and solid 
unknowns. The block preconditioning infrastructure utilised 
within the \c NavierStokesSchurComplementPreconditioner will ensure that the 
preconditioner only operates on fluid dofs.

\dontinclude fsi_channel_with_leaflet.cc
\skipline start_of_preconditioner_solve
\until Navier_stokes_preconditioner_pt->preconditioner_solve(y,z);

The fluid elements \c z_f of the vector \c z will now have been updated to
contain the action of the SchurComplement preconditioner on the fluid 
elements \c y_f of the vector \c y.

To apply the fluid-solid coupling matrix vector product \f$C_{sf}\f$, 
we copy the fluid elements from \c z into another vector \c z_f. We then apply
the matrix-vector product operator to  \c z_f and store the result in a
vector \c w. Finally, we subtract \c w from the solid residuals 
\c y_s and store the result in \c w to complete the application of 
the matrix-vector product.

\dontinclude fsi_channel_with_leaflet.cc
\skipline Fluid Onto Solid Matrix Vector Product Operator
\until w = y_s;


Finally, we apply the solid subsidiary preconditioner \f$\tilde S\f$ to \c
w and return the result to \c z. We note that because the solid subsidiary
preconditioner is not a block preconditioner, the preconditioner solve
method must be called with the solid block vectors. The result is then
copied to the full-length vector \c z which contains the fluid and solid dofs.

\dontinclude fsi_channel_with_leaflet.cc
\skipline Solid Subsidiary Preconditioner
\until }




<HR>
<HR>


\section para Parallelisation

We note that the above discussion did not address the parallelisation 
of the preconditioners. This is because all the required parallel features
are "hidden" within the block preconditioning framework which relies heavily on the
library's 
<a href="../../../mpi/distributed_linear_algebra_infrastructure/html/index.html">
distributed linear algebra infrastructure.</a> Any of the
preconditioners discussed in this tutorial can therefore be used 
without change when \c oomph-lib is compiled with MPI support and if the 
the executable is run on multiple processes.


<HR>
<HR>

\section sources Source files for this tutorial
- The source file for the simple block diagonal preconditioner
  for the linear elasticity problem is
<CENTER>
<A HREF="../../../../demo_drivers/linear_solvers/simple_block_preconditioners.h">
demo_drivers/linear_solvers/simple_block_preconditioners.h
</A>
</CENTER>
\n
- The driver code demonstrating the use of the simple block 
  diagonal preconditioner for the linear elasticity problem is
<CENTER>
<A HREF="../../../../demo_drivers/linear_solvers/two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc">
demo_drivers/linear_solvers/two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc
</A>
</CENTER>
\n
- The source files for the "multi-poisson" preconditioners and the
serial driver codes are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/linear_solvers/">
demo_drivers/linear_solvers/
</A>
</CENTER>
\n
- The serial "multi-poisson" driver code (which demonstrates the use
  of the various "multi-poisson" preconditioners discussed above) is:
<CENTER>
<A HREF="../../../../demo_drivers/linear_solvers/two_d_multi_poisson.cc">
demo_drivers/poisson/two_d_multi_poisson.cc
</A>
</CENTER>
\n
- The parallel counterpart is here (note that, as claimed, this code
  uses exactly the same preconditioners as the serial version):
<CENTER>
<A HREF="../../../../demo_drivers/mpi/solvers/two_d_multi_poisson.cc">
demo_drivers/mpi/solvers/two_d_multi_poisson.cc
</A>
</CENTER>
\n
- The (parallel) driver code which demonstrates the implementation and 
  use of the simple FSI preconditioner (for the "channel with leaflet" problem)
  is here:
<CENTER>
<A HREF="../../../../demo_drivers/mpi/solvers/fsi_channel_with_leaflet.cc">
demo_drivers/mpi/solvers/fsi_channel_with_leaflet.cc
</A>
</CENTER>
\n
.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


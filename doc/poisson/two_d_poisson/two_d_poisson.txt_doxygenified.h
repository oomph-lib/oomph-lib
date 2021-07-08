/**

\mainpage Demo problem: A two-dimensional Poisson problem
In this document, we demonstrate how to solve a 2D Poisson problem
using existing objects from the \c oomph-lib library: 
<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Two-dimensional model Poisson problem</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2 \frac{\partial^2u}{\partial x_i^2} = f(x_1,x_2),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the rectangular domain \f$D =\left\{ (x_1,x_2) \in 
[0,1] \times [0,2]\right\}\f$, 
with Dirichlet boundary conditions
\f[
\left. u\right|_{\partial D}=u_0 
\ \ \ \ \ \ \ \ \ \ (2)
\f]
where the function \f$ u_0 \f$ is given.
</TD>
</TR>
</TABLE>  
</CENTER>


We provide a detailed discussion of the driver code 
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson.cc">
two_d_poisson.cc</A> which solves the problem for
\f[
u_0(x_1,x_2) = \tanh(1-\alpha(x_1 \tan\Phi - x_2)),
\ \ \ \ \ \ \ \ \  (3)
\f]
and
\f[
f(x_1,x_2) = \sum_{i=1}^2 \frac{\partial^2 u_0}{\partial x_i^2},
\ \ \ \ \ \ \ \ \  (4)
\f]
so that \f$u_0(x_1,x_2) \f$ represents the exact solution of the
problem. For large values of \f$ \alpha \f$ the solution approaches
a step function 
\f[
u_{step}(x_1,x_2) = \left\{
\begin{array}{rl}
-1 & \mbox {for $x_2 < x_1 \ \tan\Phi$} \\
1 & \mbox {for $x_2 > x_1 \ \tan\Phi$} 
\end{array}
 \right.
\f]
which presents a serious challenge for any numerical method.
The figure below compares the numerical and exact solutions for
\f$ \alpha = 1 \f$ and \f$\Phi = 45^o\f$.

\image html two_d_poisson_soln.gif "Plot of the solution " 
\image latex two_d_poisson_soln.eps "Plot of the solution " width=0.75\textwidth

<HR>   
<HR>


\section global Global parameters and functions
Following our usual practice, we use a namespace, 
\c TanhSolnForPoisson, to define the source function (4)
and the exact solution (3). Both functions permit
arbitrary values of the tangent \f$ \tan \Phi \f$ and the 
steepness parameter  \f$ \alpha \f$, which are stored in
TanhSolnForPoisson::TanPhi and TanhSolnForPoisson::Alpha, respectively,
so that the "user" can set their values from the driver code.

\dontinclude two_d_poisson.cc
\skipline start_of_namespace
\until end of namespace


<HR>
<HR>

\section main The driver code
In order to solve the 2D Poisson problem using \c oomph-lib, we represent the 
mathematical problem defined by equations (1) 
and (2) in a specific \c Problem object, \c
PoissonProblem. \c oomph-lib provides a variety 
of 2D Poisson elements (e.g. 2D quadrilateral elements with bi-linear,
bi-quadratic and bi-cubic representations for the unknown function)
and we pass the specific element type as a template parameter to the 
\c Problem. In the driver code, listed below, we use the
\c QPoissonElement<2,3>, a nine-node (bi-quadratic) 2D Poisson element. 

The next few lines of the \c main() function
create a \c DocInfo object -- an \c oomph-lib
object that collates various items of data that can be used
to label output files: Here we specify that the output files
are to be written to the directory "RESLT", and that the first batch
of output files should be labelled with the identifier "0". 
See the discussion of the postprocessing routine \c doc_solution(...)
for details.
[\b Note: While the ability to specify an output directory from the driver
code is useful, it does rely on the "user" having created the
directory before the code is executed. We could use the
C++ \c system(...) function to issue a system command which creates
the directory if it does not exist. Since this would make
the code non-portable, we only issue a warning suggesting the 
likely cause of the problem if the output file cannot be opened.
If you want to make absolutely sure that the output directory
does exist and can be written to, you can change this forgiving 
behaviour with the function \c DocInfo::directory_must_exist(). 
This function provides access to a boolean flag which is set to 
\c false by default. If set to \c true, the code execution terminates 
with \c assert(false) if the directory specified with
\c DocInfo::set_directory(...) cannot be written to.] 




Next we execute the \c Problem::self_test()
function to check whether the \c Problem has been correctly initialised.
If this test is passed, we proceed to the solution. 
We choose the angle of the "step" as 45 degrees (corresponding to 
\f$ \tan \Phi = 1 \f$) and then solve the problem for a number of 
values of the steepness parameter \f$ \alpha. \f$ We document
each solution with the post-processing routine \c doc_solution(...)
which accesses the step number and the output directory
via the \c DocInfo object. 

\dontinclude two_d_poisson.cc
\skipline start_of_main
\until end of main


<HR>
<HR>

\section problem The problem class
The \c PoissonProblem is derived from \c oomph-lib's generic
\c Problem class and the specific element
type is specified as a template parameter to make it easy for 
the "user" to change the element type from the driver code.

\dontinclude two_d_poisson.cc
\skipline start_of_problem_class
\until public Problem

The problem class has five member functions, only three of which
are non-trivial:
- the constructor \c PoissonProblem(...)
- the function \c actions_before_newton_solve()
- the function \c doc_solution(...)
.
The function \c Problem::actions_after_newton_solve() is a pure
virtual member function of the \c Problem base class and must
be provided. However, it is not required in the present problem
and we leave it empty. Similarly, the problem destructor can remain
empty as all memory de-allocation is handled in the destructor of the 
\c Problem base class. The Problem only stores one private data
member, the pointer to the source function. 


\skipline public:
\until end of problem class

[See the discussion of the 
<A HREF="../../../poisson/one_d_poisson/html/index.html">
1D Poisson problem</A> for a more detailed discussion of the
function type PoissonEquations<2>::PoissonSourceFctPt.]

<HR>
<HR>

\section constructor The Problem constructor
In the \c Problem constructor, we start by discretising the rectangular domain,
using \c oomph-lib's \c SimpleRectangularQuadMesh object. The arguments of
this object's constructor are the number of elements (whose type is specified
by the template parameter), and the domain lengths in the \f$x_1 \f$ 
and \f$x_2 \f$ directions, respectively. 


The subsequent lines of code pin the nodal values along the 
entire domain boundary. In the 
<A HREF="../../../poisson/one_d_poisson/html/index.html"> 1D example</A>
considered earlier, the identification of the nodes
on the domain boundaries was trivial. In higher-dimensional problems,
this task can become rather involved. \c oomph-lib's \c Mesh base class
provides the helper function \c Mesh::boundary_node_pt(...)
giving (pointer-based) access to nodes on specified 
mesh boundaries. [The total number of boundaries can be obtained from
\c Mesh::nboundary(), while the number of nodes on a specific boundary
is available from \c Mesh::nboundary_node(...).]
The nested loops over the mesh boundaries and 
the nodes on these boundaries therefore provide a convenient and
completely generic method of accessing all boundary nodes. 


Finally we loop over all elements to assign the source function
pointer, and then call the generic \c
Problem::assign_eqn_numbers() routine to set up the equation
numbers.


\skipline start_of_constructor
\until end of constructor



<HR>
<HR>

\section actions_before "Actions before solve" 
We use \c Problem::actions_before_newton_solve() to update the boundary 
conditions in response to possible changes in the problem parameters.
We use the exact solution, specified in \c
TanhSolnForPoisson::get_exact_u(...),
to determine the boundary values that are appropriate for the current
values of \f$ \alpha \f$ and \f$ \tan \Phi \f$.

\skipline start_of_actions_before_newton_solve
\until end of actions before solve


[See the discussion of the  
<A HREF="../../../poisson//one_d_poisson/html/index.html">
1D Poisson problem</A> for a more detailed discussion of the
pure virtual functions  \c Problem::actions_before_newton_solve() and 
 \c Problem::actions_after_newton_solve().]


<HR>
<HR>


\section doc Post-processing
The function \c doc_solution(...) writes the FE solution and the corresponding 
exact solution, defined in \c TanhSolnForPoisson::get_exact_u(...)
to disk. The \c DocInfo object specifies the output directory
and the label for the file names. [See the discussion of the  
<A HREF="../../../poisson/one_d_poisson/html/index.html">
1D Poisson problem</A> for a more detailed discussion of the
generic \c Mesh member functions \c Mesh::output(...),
 \c Mesh::output_fct(...) and \c Mesh::compute_error(...)].

\skipline start_of_doc
\until end of doc


<HR>
<HR>


\section exercises Comments and exercises
-# In its current form, the number of elements and the 
   dimensions of the domain are hard-coded in the \c Problem
   constructor. Change the \c Problem constructor so that
   these quantities become input parameters that can be set from the
   \c main() function. 
-# Note how the accuracy of the FE solution decreases as
   the steepness of the "step" is increased:
\image html alpha.gif "Plot of the solution for different values of the steepness parameter " 
\image latex alpha.eps "Plot of the solution for different values of the steepness parameter " width=0.75\textwidth
   How many elements are required to resolve the solution 
   with \f$ \alpha = 10 \f$ as accurately as in the case
   of  \f$ \alpha = 1\f$? [\b Note: Since the solution
   has steep gradients only in a very narrow region, uniform mesh refinement
   is an extremely wasteful method of improving the accuracy of the
   computed solution. \c oomph-lib provides powerful mesh
   adaptation routines which perform fully-automatic mesh refinement
   and unrefinement, based on a posteriori error estimates of
   the solution. We will demonstrate these in 
   <A HREF="../../../poisson/fish_poisson/html/index.html">
   another example.</A>]
-# Repeat the numerical experiments with different element types.
   Replace the nine-node Poisson element, \c QPoissonElement<2,3>,
   by its lower- and higher-order counterparts \c QPoissonElement<2,2> and
   \c QPoissonElement<2,4>, respectively. Compare the total number of
   degrees of freedom, the errors, and the run-times for the 
   different discretisations.

\subsection precompile_mesh Header files and precompiled meshes
We have repeatedly stressed that \c oomph-lib \c Mesh objects
are (and any user-written ones should be) templated by the element
type, so that meshes can be used with all finite elements that are derived from
the same geometric element (2D quad elements from the \c QElement
family, say). Typically, the element type is specified in the
driver code. Consequently, the compiler must instantiate 
the \c Mesh object for a particular element type when the driver
code is compiled -- there is no point in trying
to "pre-compile" a \c Mesh object for "all possible element types".
The source code for \c Mesh objects is therefore usually contained
in a single (header) file which must be included in the 
driver code. The first few lines of the driver code
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson.cc">
two_d_poisson.cc</A> illustrate the technique:

 
\dontinclude two_d_poisson.cc
\skipline Driver for
\until meshes/simple

The code uses objects from the \c generic and \c poisson libraries
whose function prototypes are contained in the header files
\c generic.h and \c poisson.h, located in the \c oomph-lib \c include
directory. All objects in these libraries are fully instantiated
and no re-compilation is required -- we simply link against the
libraries which are located in \c oomph-lib's \c lib directory. 
The mesh header files (which include the entire source code for
each mesh) are located in the include (sub-)directory 
\c include/meshes, and are included into the driver code with
a C++ include directive

  While this strategy greatly facilitates code reuse, it can
incur significant compile-time overheads as the (possibly very
lengthy) mesh sources must be recompiled whenever the driver code is changed.
During code development, this overhead can become unacceptable.
To avoid the constant re-compilation of the mesh sources,
all \c oomph-lib mesh objects are contained in two
separate source files. In the case of the \c SimpleRectangularMesh,
the class definition and function prototypes are contained in 
the small auxiliary header file \c simple_rectangular_quadmesh.template.h,
while the actual function definitions are contained in 
\c simple_rectangular_quadmesh.template.cc. These are the only
sources that the mesh-writer has to provide. The header file
\c simple_rectangular_quadmesh.h is generated (automatically)
by concatenating the two *.template.* files and all three files
are contained in the mesh include directory. This allows the "user"
to pre-compile the mesh for a specific element type (or for 
a range of specific elements) to produce a separate object
file that can be linked against when the driver code is built.

 The procedure is illustrated in the alternative source code  
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson2.cc">
two_d_poisson2.cc</A> and the associated mesh file, 
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson2_mesh.cc">
two_d_poisson2_mesh.cc</A>. In the original version of the
code, 
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson.cc">
two_d_poisson.cc</A>, the mesh was instantiated with
the element type \c QPoissonElement<2,3> and we will assume
that this is the only element type required in the driver code.
We force the instantiation of the \c SimpleRectangularQuadMesh 
for this element type by employing the C++ "template" statement in
the mesh file,
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson2_mesh.cc">
two_d_poisson2_mesh.cc</A>, which is listed in its entirety here:

\include two_d_poisson2_mesh.cc

This source file can be pre-compiled into an object file,
\c two_d_poisson2_mesh.o, say. 

The driver code only needs to include the templated header
file (which contains the class definition and the
function prototypes) so that the first few lines of the
modified driver code look like this:

\dontinclude two_d_poisson2.cc
\skipline Driver for
\until meshes/simple

The driver code can now be compiled separately (without having to
recompile the mesh sources every time) and the correctly instantiated
version of the \c SimpleRectangularQuadMesh can be made available
by including \c two_d_poisson2_mesh.o during the linking phase.


\subsection different_solvers How to choose the linear solver for the Newton method

\c oomph-lib treats all problems as nonlinear problems and provides
steady (and unsteady) Newton solvers to solve the system of nonlinear algebraic
equations that arise from the spatial (and temporal) discretisation
of the governing equations. Typically, the repeated assembly of the
Jacobian matrix and the solution of the linear systems during 
the Newton iteration provides the major part of the computational work. 
Within this framework linear problems are simply special cases
of nonlinear problems for which the Newton method converges in one
iteration. The assembly of the Jacobian matrix and the
solution of the linear system is performed by \c oomph-lib's 
\c LinearSolver objects. These typically provide interfaces
to general purpose linear solvers such as \c SuperLUSolver (our default
solver). The list of solvers includes:
- \c SuperLUSolver: An interface to  Demmel, Eistenstat, Gilbert, 
  Li & Liu's serial SuperLU solver. See 
  <A HREF="http://crd.lbl.gov/~xiaoye/SuperLU/">
           http://crd.lbl.gov/~xiaoye/SuperLU/</A> for details.
- \c HSL_MA42: An interface to the MA42 frontal solver from the HSL
  library. See <A HREF="http://www.hsl.rl.ac.uk/">
  http://www.hsl.rl.ac.uk</A>
  for details.
- \c FD_LU: An extremely inefficient solver which computes the
  Jacobian matrix by finite differencing and stores it in a dense
  matrix. This solver is mainly provided to facilitate sanity checks 
  during code development. (The residuals are easier to compute than the 
  the Jacobian matrix!)
- ...and many others. See the 
  <a href="../../../linear_solvers/html/index.html">Linear Solvers 
  Tutorial</a> for a more detailed discussion of \c oomph-lib's
  various direct and iterative solvers.
.

By default the \c oomph-lib Newton solvers use \c SuperLU with
compressed row storage for the Jacobian matrix as the linear solver.
To change the linear solver to another type you can over-write
the \c Problem's pointer to its linear solver. For instance,
to change the linear solver to \c HSL_MA42, add the following
lines to the \c Problem constructor:

\code

 /// Build a linear solver: Use HSL's MA42 frontal solver
 Problem::linear_solver_pt() = new HSL_MA42;

 /// Switch on full doc for frontal solver
 static_cast<HSL_MA42*>(Problem::linear_solver_pt())->doc_stats()=true;

\endcode

\c HSL_MA42 can document various statistics such 
the memory usage etc. This is enabled with 
the second command. Other solvers have similar member functions.
See the <A HREF="../../../the_data_structure/html/index.html">
full documentation of all \c oomph-lib classes</A> for details.

We provide an example code 
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson_compare_solvers.cc">
two_d_poisson_compare_solvers.cc</A> which can be used to 
explore the performance of various linear solvers for the
2D Poisson problem considered above.



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/">
demo_drivers/poisson/two_d_poisson/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/poisson/two_d_poisson/two_d_poisson.cc">
demo_drivers/poisson/two_d_poisson/two_d_poisson.cc
</A>
</CENTER>
.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


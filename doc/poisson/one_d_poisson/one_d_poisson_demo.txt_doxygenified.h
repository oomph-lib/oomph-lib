/**

\mainpage Demo problem: A one-dimensional Poisson problem
In this document, we demonstrate how a solve the 1D Poisson problem
using existing objects from the \c oomph-lib library:
<CENTER>
<TABLE>  
<TR> 
<TD>  
<CENTER>
<B>One-dimensional model Poisson problem</B>
</CENTER> 
Solve
\f[
\frac{\mbox{d}^2u}{\mbox{d}x^2} = f(x),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in a one-dimensional domain \f$ x\in[0,1] \f$, with boundary conditions 
\f[
u(0)=u_0 \ \mbox{\ \ \ and \ \ \ } \ \ u(1)=u_1, 
\ \ \ \ \ \ \ \ \ \ (2)
\f]
where \f$ u_0 \f$ and \f$ u_1 \f$  are given.
</TD>
</TR>
</TABLE>  
</CENTER>

We provide a detailed discussion of the driver code 
<A HREF="../../../../demo_drivers/poisson/one_d_poisson/one_d_poisson.cc">
one_d_poisson.cc</A> which solves the 
problem for the specific source function
\f[
f(x) = \pm 30 \sin(\sqrt{30} x),  \ \ \ \  (3)
\f]
and boundary conditions
\f[
u_0= 0; \ \ \  u_1 = \mp 1.
\f]
In this case, the problem has the (fish-shaped) exact solution
\f[
u(x) = \pm \left[ \left(\sin(\sqrt{30})-1\right) x - \sin(\sqrt{30} x)\right].
\ \ \ \ \  (4)\f]

\image html one_d_fish.gif "The fish-shaped solution u(x). " 
\image latex one_d_fish.eps "The fish-shaped solution u(x). " width=0.75\textwidth

<HR>   
<HR>


\section global Global parameters and functions
Source functions and problem parameters generally need
to be accessible to a variety of \c oomph-lib objects and we 
tend to define such functions/quantities in C++ namespaces.
Here we use the namespace \c FishSolnOneDPoisson to define
the source function  (3) and the exact
solution (4). Both functions use the integer value
FishSolnOneDPoisson::Sign, which can be set by the "user" in 
the driver code.

\dontinclude one_d_poisson.cc
\skipline start_of_namespace
\until end of namespace


<HR>
<HR>

\section main The driver code
In order to solve the Poisson problem with \c oomph-lib, 
we represent the mathematical
problem defined by equations (1) and (2)
in a specific \c Problem object, \c OneDPoissonProblem. 
\c oomph-lib provides a variety 
of 1D Poisson elements (1D elements with linear, quadratic and cubic
representations for the unknown function) and we pass the
specific element type as a template parameter to the \c Problem.
In the driver code, listed below, we use the \c QPoissonElement<1,4>, 
a four-node (cubic) 1D Poisson element. We pass the number of elements in the
mesh as the first argument to the \c Problem constructor and the pointer
to the source function, defined in the namespace \c FishSolnOneDPoisson, as the second.

Once the problem has been created, we execute \c Problem::self_test(),
which checks that the \c Problem has been properly set up.
If this test is passed, we proceed to solve the problem:
Initially, we set the sign of the source function (defined in the
variable \c FishSolnOneDPoisson::Sign)  to -1 and solve, 
using \c oomph-lib's Newton solver. We write the
solution to output files, using the \c OneDPoissonProblem's 
member function \c doc_solution(), discussed below. We then repeat 
the process, using a positive sign for the source function.

\dontinclude one_d_poisson.cc
\skipline start_of_main
\until end of main


<HR>
<HR>

\section problem The problem class
\c oomph-lib driver codes tend to be 
very compact and "high-level" because all the "hard work" is done 
in the \c Problem specification. For our simple Poisson 
problem, this step is completely straightforward:

The \c OneDPoissonProblem is derived from \c oomph-lib's generic
\c Problem class and, as discussed above, the specific element
type is specified as a template parameter to make it easy for 
the "user" to change the element type in the driver code.

\dontinclude one_d_poisson.cc
\skipline start_of_problem_class
\until public Problem

The \c OneDPoissonProblemclass has five member functions, only three of which
are non-trivial:
- the constructor \c OneDPoissonProblem(...)
- the function \c actions_before_newton_solve()
- the function \c doc_solution(...)
.
The function \c Problem::actions_after_newton_solve() is a pure
virtual member function of the \c Problem base class and must 
be provided. However, it is not required in the present problem so
we leave it empty. Similarly, the problem destructor can remain
empty as all memory de-allocation is handled in the destructor of the 
\c Problem base class. The Problem only stores one private data
member, the pointer to the source function.

\skipline public:
\until end of problem class
 

 <TABLE>
 <TR>
 <TD bgcolor="cornsilk">
<CENTER><B>\anchor gen_conv A general convention</B></CENTER>
The type \c PoissonEquations<1>::PoissonSourceFctPt, used to 
define the type of the source function pointer,
is a public typedef, defined in \c oomph-lib's Poisson equation
class, as follows:

\code

[...]

public:

 /// \short Function pointer to source function has the form fct(x,f(x)). 
 /// Note that x is a Vector! 
 typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);

[...]

\endcode

This reflects a \b general \b convention in \c oomph-lib: The function types 
of source functions etc. that are required by specific elements,
are always declared as public types in the element classes. 
The logic behind this is that only the element writer/maintainer knows 
what type of function (i.e. the type of its arguments and its return 
value) a specific element requires.  For instance, the source function 
for the Poisson equation requires the spatial coordinate \c x as 
input and computes the (scalar) value of the source function.
Since the syntax for C++ function pointers is somewhat "non-obvious" 
(that's according to Bjarne Stroustrup, the designer of C++,
himself!), we use typedefs to give the function pointers more
intuitive names. Since the typedefs are public, they can be used
anywhere in the "user's" code.

 While we're at it, here's \b another
\b convention: If an \c oomph-lib function has input and output
arguments in its argument list, the input arguments appear first
(and are usually passed as constant references), while the output
arguments appear last (and are passed as references).
</TD>
</TR>
</TABLE>

<HR>
<HR>

\section constructor The Problem constructor
In the \c Problem constructor, we define the domain length and 
build a \c  Mesh, using \c oomph-lib's \c OneDMesh object
which is templated by the element type. The required number of elements and 
the domain length are passed as arguments to the \c OneDMesh
constructors. The subsequent lines of code 
pin the nodal values at the two boundary nodes. Next we pass the 
source function pointer to the elements. Finally, we call the generic \c
Problem::assign_eqn_numbers() routine which does precisely 
what it says...

\skipline start_of_constructor
\until end of constructor

The cast in the loop over the elements
is required because \c Mesh::element_pt(...) returns 
a pointer to the element base class \c GeneralisedElement, which does not 
have an access function to the Poisson element's source function pointer.

 <TABLE>
 <TR>
 <TD bgcolor="cornsilk">
<CENTER><B>\anchor const_conv A general convention</B></CENTER> 
 It might seem more natural to pass
essential parameters (such as the Poisson element's source function 
pointer) to an element when the element is created. If we 
made the source function pointer an argument of the element constructor,
we would not have to execute this loop in the constructor.

However, \c oomph-lib employs a general \b convention 
that element constructors should not have any arguments.
This is because  elements are usually created
by the \c Mesh constructor. To 
allow \c Meshes that were originally developed for one particular 
element type (e.g. a quadrilateral Poisson element) to be used with
other elements of same geometry (e.g. a quadrilateral Navier-Stokes 
element), \c Mesh objects are usually templated by the 
element type. The \c Mesh constructor creates elements 
by calling the element's default (argument-free!) 
constructor. If we were to set any element-specific arguments via arguments
to the element constructor it would be impossible to re-use the
\c Mesh with other element types -- Navier-Stokes elements, for instance, do
not have a source function pointer. The actions performed in
this loop are therefore
fairly typical as most non-trivial elements need to be passed
some additional information to become fully functional.

</TD> 
</TR>
</TABLE>

<HR>
<HR>

\section actions_before "Actions before solve" 
The pure virtual function \c Problem::actions_before_newton_solve() must
be implemented for all specific \c Problems and, as the name
suggests, should perform any actions that need to be performed
before the system of equations is solved. In the current problem,
we use \c Problem::actions_before_newton_solve() to update the boundary 
conditions in response to possible changes in the sign of the source function.
We use the exact solution (specified in the namespace \c FishSolnOneDPoisson)
to determine the boundary values that are appropriate for the sign 
specified in \c FishSolnOneDPoisson::Sign.

\skipline start_of_actions_before_newton_solve
\until end of actions before solve



<HR>
<HR>


\section doc Post-processing
The function \c doc_solution(...) writes the FE solution and the corresponding 
exact solution, defined in \c FishSolnOneDPoisson::get_exact_u(...)
to disk. The argument \c label is used to add identifiers 
to the output file names. Note that all output functions
are implemented in the generic \c Mesh class: 
- The function \c Mesh::output(...) executes the \c 
  FiniteElement::output(...) function for each
  element in a mesh. For 1D Poisson elements, this function
  writes the values of \f$ x \f$ and  \f$ u(x) \f$ at \c npts uniformly spaced
  points in the element to the specified file. 
- The function \c Mesh::output_fct(...) plots the function specified
  by the function pointer in its last argument at the specified
  number of points in each of the constituent elements. This allows
  point-by-point comparisons between exact and FE solutions. Here
  we plot the exact solution at a larger number of points to ensure
  that the exact solution looks smooth even if only a small number of
  elements are used for the discretisation of the ODE.
.
Finally, we call the function \c Mesh::compute_error(...)
which determines the square of the L2 error, based on the
difference between the exact solution (specified by a function
pointer) and the FE solution. We also plot the pointwise error
in the specified output file.

\skipline start_of_doc
\until end of doc


<HR>
<HR>


\section exercises Exercises
-# Run the code with different numbers of elements. 
   How does the error between the exact and the analytical solution 
   change? 
-# Compare the error obtained with different element types --
   replace the four-node Poisson element, \c QPoissonElement<1,4>,
   by its lower-order counterparts   \c QPoissonElement<1,3> and
   \c QPoissonElement<1,2>.
-# The fish-shaped exact solution (4) is fairly smooth.
   Postulate a more rapidly varying "exact" solution, such as
   \f[
   u(x) = \tanh\left[\alpha \left(x-\frac{1}{2}\right)\right]
   \f]
   which produces a "step" at \f$ x=1/2 \f$ when \f$ \alpha \f$ becomes
   sufficiently large. Calculate the source function required for this
   function to be an exact solution and implement both functions in another
   namespace, \c TanhSolnOneDPoisson, say. Replace the reference
   to \c FishSolnOneDPoisson by \c TanhSolnOneDPoisson and repeat 
   the above exercises.
-# Remove the Dirichlet boundary condition at the left end of the
   domain. What do you observe? [We shall return to this question
   in <A HREF="../../../poisson/two_d_poisson_flux_bc/html/index.html">
   another example</A>
   where we discuss the 
   <A HREF="../../../poisson/two_d_poisson_flux_bc/html/index.html">
   application of Neumann-type boundary conditions</A>.]


<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/one_d_poisson/">
demo_drivers/poisson/one_d_poisson/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/poisson/one_d_poisson/one_d_poisson.cc">
demo_drivers/poisson/one_d_poisson/one_d_poisson.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


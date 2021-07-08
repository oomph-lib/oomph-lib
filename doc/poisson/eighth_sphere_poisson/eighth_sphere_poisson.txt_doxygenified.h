/**

\mainpage Example problem: Adaptive solution of the 3D Poisson equation in a spherical domain

Following the numerous 2D problems discussed in earlier examples we now
demonstrate that the solution of 3D problems is just as easy. 
For this purpose we discuss the adaptive solution of the 3D Poisson problem
<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Three-dimensional model Poisson problem</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^3 \frac{\partial^2u}{\partial x_i^2} = f(x_1,x_2,x_3),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the "eighth-of-a-sphere" domain \f$D\f$, 
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

We choose a source function and boundary conditions for which
\f[
u_0(x_1,x_2,x_3) = 
\tanh\left(1-\alpha(({\bf x}-{\bf x_0}) \cdot {\bf N})\right),
\ \ \ \ \ \ \ \ \  (3)
\f]
is the exact solution. Here where \f${\bf x} = (x_1,x_2,x_3)\f$  
is the vector of the spatial coordinates, and the vectors 
\f${\bf x}_0 = \left(x_1^{(0)},x_2^{(0)},x_3^{(0)}\right)\f$  and 
\f${\bf N}_0 = (N_1,N_2,N_3)\f$ are constants. For large values of the
constant \f$ \alpha \f$ the solution varies rapidly across the plane through
\f${\bf x}_0\f$ whose normal is given by  \f${\bf N}\f$.


Here are some plots of the exact and computed solutions for 
\f${\bf x}_0 = (0,0,0)\f$, \f${\bf N}_0 = 1/\sqrt{3} \ (-1,-1,1)\f$,
and \f$ \alpha=50\f$ at various levels of mesh refinement.
Note that the plot of the exact solution was produced by
setting the nodal values to the exact solution, obtained by evaluating 
(3) at the nodal positions. The elements' basis functions 
were then used to interpolate between the nodal values. On the coarse
meshes, the interpolation between the "exact" nodal values is clearly
inadequate to resolve the rapid variation of the solution. 


\image html eight_sphere_poisson.gif "Plot of the solution " 
\image latex eight_sphere_poisson.eps "Plot of the solution " width=0.75\textwidth

<HR>
<HR>


\section global Global parameters and functions
Following our usual practice, we use a namespace, 
\c TanhSolnForPoisson, to define the source function, the exact 
solution and various problem parameters.

\dontinclude eighth_sphere_poisson.cc
\skipline start_of_namespace
\until end of namespace



<HR>
<HR>

\section main The driver code

The driver code solves the 3D Poisson problem with full spatial
adaptivity -- a fairly time-consuming process. To minimise the
run-times when the code is executed during \c oomph-lib's self-tests,
we use command line arguments to optionally limit the number of 
adaptive refinements. If the code is run with a(ny) command line
arguments, only a single adaptive refinement is performed;
otherwise up to four levels of refinement are permitted. 
\c oomph-lib provides storage for the command line arguments
in the namespace \c CommandLineArgs to make them accessible to 
other parts of the code.

  Otherwise the driver code is very similar to that used in the
<A HREF="../../../poisson/fish_poisson/html/index.html">corresponding 2D Poisson problems</A>: We construct the problem, passing the
pointer to the source function. Next, we create a \c DocInfo object to 
specify the output directory, and execute the global self-test to 
assert that the problem has been set up correctly. Next we
solve the problem on the coarse initial mesh (comprising
four 27-node brick elements) and then adapt the problem
based on the elemental error estimates, until the maximum number of
adaptations has been reached or until the adaptation ceases to
changes the mesh.

\dontinclude eighth_sphere_poisson.cc
\skipline start_of_main
\until end of main


<HR>
<HR>

\section problem The problem class

The problem class has the usual structure -- the only difference to
the corresponding 2D codes is that the assignment of the boundary 
conditions in \c actions_before_newton_solve() now involves thee 
nodal coordinates rather than two.

\dontinclude eighth_sphere_poisson.cc
\skipline start_of_class_definition
\until  end of class definition

[See the discussion of the 
<A HREF="../../../poisson/one_d_poisson/html/index.html">
1D Poisson problem</A> for a more detailed discussion of the
function type PoissonEquations<3>::PoissonSourceFctPt.]

<HR>
<HR>

\section constructor The Problem constructor
In the \c Problem constructor, we set the "steepness parameter" 
\f$ \alpha \f$ to a large value and create the mesh for a 
a sphere of radius 5. Next, we create the error estimator and
pass it to the adaptive mesh. 

\skipline start_of_constructor
\until error_estimator_pt;

We adjust the targets for the mesh adaptation so that 
the single mesh adaptation performed during a validation run
produces a non-uniform refinement pattern. (The error targets
for this case were determined by trial and error.) The tighter
error tolerances specified otherwise are appropriate to 
properly resolve the solution, as shown in the animated gif files
at the beginning of this document.
\skipline  Adjust error targets for adaptive refinement
\until end adjustment

 Next, we assign the boundary conditions. In the present problem all boundaries
are Dirichlet boundaries, therefore we loop over all nodes 
on all boundaries and pin their values. If only a subset of the mesh 
boundaries were of Dirichlet type, only the nodes on those boundaries
would have to be pinned. "Usually" the numbering of the mesh
boundaries is (or at least should be!) documented in the mesh 
constructor but it can also be obtained from the function
\c Mesh::output_boundaries(...) whose use is illustrated here.

\skipline Doc the mesh
\until end of 


Finally we loop over all elements to assign the source function
pointer, and then call the generic \c
Problem::assign_eqn_numbers() routine to set up the equation
numbers.


\skipline Find number
\until end of constructor
 


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


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="
../../../../
demo_drivers/poisson/eighth_sphere_poisson/
">
demo_drivers/poisson/eighth_sphere_poisson/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="
../../../../
demo_drivers/poisson/eighth_sphere_poisson/eighth_sphere_poisson.cc
">
demo_drivers/poisson/eighth_sphere_poisson/eighth_sphere_poisson.cc
</A>
</CENTER>
.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


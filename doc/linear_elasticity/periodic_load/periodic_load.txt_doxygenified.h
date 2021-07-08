/**

\mainpage Example problem: The deformation of an elastic strip by a periodic traction

This is our first linear elasticity example problem. We discuss the 
non-dimensionalisation of the governing equations and their
implementation in \c oomph-lib, and then demonstrate the solution 
of a 2D problem: the deformation of an elastic strip by a periodic traction.


<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
This tutorial and the associated driver code was developed jointly
with David Rutter.
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>


<HR>
<HR>  


\section equation The governing equations

The figure below shows a sketch of a general elasticity 
problem. A linearly elastic solid body occupies the domain 
\f$ V \f$ and is loaded by a body
force \f$ {\bf F}^* \f$ and by a surface traction
\f$ {\bf t}^* \f$ which is applied along part of its boundary, 
\f$ A_{\rm tract} \f$. The displacement is prescribed along the 
remainder of the boundary, \f$ A_{\rm displ} \f$, where 
\f$ \partial V = A_{\rm tract} \cup  A_{\rm displ} \f$. 

\image html Solid_boundary_conds.gif "Sketch of a general elasticity problem: A linearly elastic body is loaded by a body force and is exposed to a prescribed traction along part of its boundary while the displacement is prescribed along the remainder of the boundary. " 
\image latex Solid_boundary_conds.eps "Sketch of a general elasticity problem: A linearly elastic body is loaded by a body force and is exposed to a prescribed traction along part of its boundary while the displacement is prescribed along the remainder of the boundary. " width=0.5\textwidth

We adopt an Eulerian approach and describe the 
deformation in terms of the displacement field 
\f$ {\bf u}^*\left( x_i^*,t^* \right) \f$ where 
\f$ x_i^* \f$ and \f$ t^* \f$ are the spatial coordinates and time, 
respectively. Throughout this document we will use index notation 
and the summation convention, and use asterisks to distinguish 
dimensional quantities from their non-dimensional counterparts.  
Denoting the density of the body by \f$ \rho \f$, the deformation
is governed by the Cauchy equations,
<CENTER>
\f[
\frac{\partial \tau_{ij}^*}{\partial x_j^*} + \rho F_i^*
=
\rho \frac{\partial^2 u_i^*}{\partial t^{2*}},
\f]
</CENTER>
where \f$ \tau_{ij}^* \f$ is the Cauchy stress tensor which, for a 
linearly elastic solid, is given by
<CENTER>
\f[
\tau_{ij}^* = E_{ijkl}^* \ e_{kl},
\f]
</CENTER>
where \f$ e_{kl} \f$ is the strain tensor,
<CENTER>
\f[
e_{ij} = \frac{1}{2} 
\left( \frac{\partial u_i^*}{\partial x_j^*}+
       \frac{\partial u_j^*}{\partial x_i^*}\right).
\f]
</CENTER>
\f$ E_{ijkl} \f$ is the 4th order elasticity tensor, which for 
a homogeneous and isotropic solid is
<CENTER>
\f[
E_{ijkl}^* = \frac{E}{1+\nu} \left(\frac{\nu}{1-2\nu} \delta_{ij}
 \delta_{kl} + \delta_{ik} \delta_{jl} \right),
\f]
</CENTER>
where \f$ E \f$ is Young's modulus, \f$ \nu \f$ is the 
Poisson ratio and \f$ \delta_{ij} \f$ is the Kronecker delta. 
Thus the Cauchy stress is given in terms of the displacement
derivatives 
by
<CENTER>
\f[
\tau_{ij}^* = \frac{E}{1+\nu} \left(\frac{\nu}{1-2\nu} \ 
\delta_{ij} \ \frac{\partial u_k^*}{\partial x_k^*} + 
\frac{1}{2}\left(\frac{\partial u_i^*}{\partial x_j^*}+
\frac{\partial u_j^*}{\partial x_i^*} \right) \right).
\f]
</CENTER>


We non-dimensionalise the equations, using a problem specific reference
length, \f$ {\cal L} \f$, and a timescale \f$ {\cal T} \f$, and  use
Young's modulus to non-dimensionalise the body force and the stress,
<CENTER>
\f[
\tau_{ij}^* = E \, \tau_{ij}, \qquad
x_i^* = {\cal L}\, x_i, \qquad
u_i^* = {\cal L}\, u_i,
\f]
\f[
t^* = {\cal T}\, t, \qquad
F_i^* = \frac{E}{\rho \cal L} \, F_i. \qquad
\f]
</CENTER>

The non-dimensional form of the Cauchy equations is then given by
<CENTER>
\f[
\frac{\partial \tau_{ij}}{\partial x_j} + F_i
=
\Lambda^2  \frac{\partial^2 u_i}{\partial t^2},
\ \ \ \ \ \ \ \ \ \ (1)
\f]
</CENTER>
where
<CENTER>
\f[
\tau_{ij} = \frac{1}{1+\nu} \left(\frac{\nu}{1-2\nu} \ 
\delta_{ij} \ \frac{\partial u_k}{\partial x_k} + 
\frac{1}{2} \left(\frac{\partial u_i}{\partial x_j}+
\frac{\partial u_j}{\partial x_i} \right) \right).
\ \ \ \ \ \ \ \ \ \ (2)
\f]
</CENTER>
The parameter
<CENTER>
\f[
\Lambda = \frac{{\cal L} \sqrt{\frac{\rho}{E}} }{\cal T},
\f]
</CENTER>
is the ratio of the elastic body's intrinsic timescale,
 \f$ {\cal L} \sqrt{\frac{\rho}{E}} \f$, to the problem-specific 
timescale, \f$ {\cal T} \f$, that we used to non-dimensionalise time.

The displacement constraints provide
a Dirichlet condition for the displacements,
\f[
u_i = u_i^{[\rm prescribed]} \ \ \ \ \ \ \mbox{on $A_{\mathrm{displ}}$},

\f]
while the traction boundary conditions require that
\f[
t_i = \tau_{ij} n_j \ \ \ \ \ \ \mbox{on $A_{\mathrm{tract}},$}
\f]
where the \f$ n_j \f$ are the components of the outer unit normal
to the boundary.

In this tutorial we only consider steady problems for which 
the equations reduce to
<CENTER>
\f[
\frac{\partial \tau_{ij}}{\partial x_j} + F_i = 0.
\f]
</CENTER>

<HR>
<HR>


\section element_types Implementation 

<HR>

\subsection el The elements
Within \c oomph-lib, the non-dimensional version of the 
<CODE>DIM</CODE>-dimensional Cauchy equations (1) with 
the constitutive equations (2) are implemented 
in the \c LinearElasticityEquations<DIM> equations class. Following
our usual approach, discussed in the
<A HREF="../../../quick_guide/html/index.html"> (Not-So-)Quick
Guide,</A> this equation class is then combined with a geometric finite
element to form a fully-functional finite element. 
For instance, the combination of the \c LinearElasticityEquations<2>
class with the geometric finite element \c QElement<2,3> yields a nine-node 
quadrilateral linear elasticity element. As usual, the mapping 
between local and global (Eulerian) coordinates within an element is given by,
<CENTER>
\f[
x_i = \sum_{j=1}^{N^{(E)}} X^{(E)}_{ij} \, \psi_j, \qquad
i=1,2\quad [\mbox{and }3],
\f]
</CENTER>
where \f$ N^{(E)} \f$ is the number of nodes in the element, 
\f$ X^{(E)}_{ij} \f$ is the \f$ i \f$-th global (Eulerian) coordinate
of the \f$ j \f$-th \c Node in the element, and the \f$ \psi_j \f$ are 
the element's shape functions, defined in the geometric
finite element. 

The cartesian displacement components \f$ u_1, \f$ \f$ u_2, \f$  
[and \f$ u_3 \f$] are stored as nodal values, and the 
shape functions are used to interpolate the displacements as
<CENTER>
\f[
u_i = \sum_{j=1}^{N^{(E)}} U^{(E)}_{ij} \, \psi_j, \qquad 
i=1,2\quad [\mbox{and }3],
\f]
</CENTER>
where \f$  U^{(E)}_{ij} \f$  is the \f$ i \f$-th displacement component 
at the \f$ j \f$-th \c Node in the element. Nodal values of the displacement 
components are accessible via the access function

\code
LinearElasticityEquations<DIM>::u(i,j)
\endcode

which returns the \f$ i \f$-th displacement component stored at the element's 
\f$ j \f$-th \c Node.



<HR>
<HR>

\section example The example problem
To illustrate the solution of the steady equations of linear
elasticity, we consider the 2D problem shown in the sketch below.

<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The problem.</B>
</CENTER> 
\image html attempt3.gif "Infinitely long strip loaded by a periodic traction. " 
\image latex attempt3.eps "Infinitely long strip loaded by a periodic traction. " width=0.5\textwidth
Solve
\f[
\frac{\partial \tau_{ij}}{\partial x_j} + F_i=0,
 \ \ \ \ \ \ \ \ \ \ (3)
\f]
in the domain \f$ D = \{x_1 \in [-\infty,+\infty], x_2 \in [0,L_y]\} \f$,
subject to the Dirichlet boundary conditions
\f[
\left. \mathbf{u}\right|_{x_2 = 0}=(0,0),
\ \ \ \ \ \ \ \ \ \ (4)
\f]
on the bottom boundary, the Neumann (traction) boundary conditions
\f[
\left. \mathbf{t}\right|_{x_2 = L_y}=\left(-A \cos{\left(\frac{2 \pi x_1}{L_x}\right)}, -A \sin{\left(\frac{2 \pi x_1}{L_x}\right)}\right),
\ \ \ \ \ \ \ \ \ \ (5)
\f]
on the top boundary, and symmetry conditions at \f$ x_1 = 0 \f$ and \f$ x_1 = L_x \f$,
\f[
\left. {\bf u}\right|_{x_1=0} = \left.{\bf u}\right|_{x_1=L_x}.
\f]
</TD>
</TR>
</TABLE>  
</CENTER>


We note that for \f$ L_y \to \infty \f$ the problem converges to the
analytical solution.
\f[
u_1^{[exact]} = -\frac{A(1+\nu)}{2 \pi} 
\cos{\left(\frac{2 \pi x_1}{L_x}\right)} \exp{\left(2 \pi (x_2-L_y)\right)},
\f]
\f[
u_2^{[exact]} = -\frac{A(1+\nu)}{2 \pi} 
\sin{\left(\frac{2 \pi x_1}{L_x}\right)} \exp{\left(2 \pi (x_2-L_y)\right)},
\f]

<HR>
<HR>

\section results Results
The figure below shows a vector plot of the displacement field
near the upper domain boundary for \f$ L_x = 1 \f$ and \f$ L_y=2 \f$.  
Note that we only discretised the infinite strip over one period of
the applied, spatially-periodic surface traction, and imposed symmetry
conditions on the left and right mesh boundaries.

The plot shows that the displacements decay rapidly with
distance from the loaded surface -- as suggested by the analytical
solution for the infinite depth case. This suggests that the computation
could greatly benefit from the use of spatial adaptivity. This is
indeed the case and is explored in 
<a href="../../refineable_periodic_load/html/index.html">
another tutorial.</a>

\image html displ.gif "Plot of the displacement field. " 
\image latex displ.eps "Plot of the displacement field. " width=0.75\textwidth


<HR>
<HR>

\section namespace Global parameters and functions
As usual, we define all non-dimensional parameters in a namespace.

\dontinclude periodic_load.cc
\skipline start_of_namespace
\until end_of_namespace

<HR>
<HR>

\section main The driver code

We start by setting the number of elements in each of the two
coordinate directions before creating a  \c DocInfo object to 
store the output directory.

\skipline start_of_main
\until doc_info.set_directory("RESLT")

We build the problem using two-dimensional \c QLinearElasticityElements, solve 
using the \c Problem::newton_solve() function, and document the 
results.

\skipline // Set up problem
\until end_of_main


<HR>
<HR>

\section problem The problem class 
The \c Problem class is very simple. As in other problems with Neumann 
boundary conditions, we provide separate meshes for the "bulk" 
elements and the face elements that apply the traction 
boundary conditions.  The latter are attached to the relevant 
faces of the bulk elements by the function \c assign_traction_elements().

\dontinclude periodic_load.cc
\skipline start_of_problem_class
\until end_of_problem_class



<HR>
<HR>

\section constructor The problem constructor
Since this is a steady problem, the constructor is quite simple. We
begin by building the meshes and pin the displacements on 
the appropriate boundaries. We then assign the boundary values
for the displacements along the bottom boundary. We either 
set the displacements to zero or assign their values 
from the exact solution for the infinite depth case.



\skipline start_of_constructor
\until end_loop_over_boundary_nodes

Next we pass a pointer to the elasticity tensor (stored in
\c Global_Physical_Variables::E) to all elements.

\skipline Complete the problem
\until end loop over elements

We loop over the traction elements and specify the applied traction.

\skipline Loop over the traction elements
\until end loop over traction elements

The two submeshes are now added to the problem and a global mesh is constructed before the equation numbering scheme is set up, using the function \c assign_eqn_numbers().

\skipline Add the submeshes to the problem
\until end of constructor
 

<HR>
<HR>

\section traction_elements The traction elements
In anticipation of the extension of this code to its 
<a href="../../refineable_periodic_load/html/index.html">adaptive
counterpart</a>, we create the face elements that apply the traction to
the upper boundary in a separate function. 

\skipline start_of_traction
\until end of assign_traction_elements


<HR>
<HR>

\section doc Post-processing
As expected, this member function documents the computed solution.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>


\section comments Comments and Exercises

\subsection nondim Comments

As discussed in the introduction, the non-dimensional version of the
steady Cauchy equations only contains a single non-dimensional
parameter, the Poisson ratio \f$ \nu \f$ which is passed to the
constructor of the \c IsotropicElasticityTensor. If you inspect
the relevant source code 
<a href="../../../../src/linear_elasticity/elasticity_tensor.h">
src/linear_elasticity/elasticity_tensor.h</a> you will find that this 
constructor has a second
argument which defaults to one. This argument plays the role of
Young's modulus and is best interpreted as the ratio of the material's
actual Young's modulus to the (nominal) Young's modulus used
in the non-dimensionalisation of the equations. The ability to provide
this ratio is important if different regions of the body contain
materials with different material properties.

\subsection exercises Exercises

-# Fix the size of the domain and set the displacements along the 
   bottom boundary to the exact solution for the infinite depth case,
   i.e. \f$ {\bf u} = {\bf u}^{[\mathrm{exact}]} \f$, using the \c
   Global_Parameters::Finite flag. Then investigate how the solution 
   converges to the exact solution for increasing numbers of elements.
-# Try varying the depth of the domain by changing 
   \c Global_Parameters::Ly while maintaining a constant spatial
   resolution (i.e. increasing the number of elements -- this is
   already done in the driver code where we compute \c ny in terms of \c
   Global_Parameters::Ly ) and compare how the solution converges to the 
   exact solution of the infinite depth case.
   

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/linear_elasticity/periodic_load/">
demo_drivers/linear_elasticity/periodic_load/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="../../../../demo_drivers/linear_elasticity/periodic_load/periodic_load.cc">
demo_drivers/linear_elasticity/periodic_load/periodic_load.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


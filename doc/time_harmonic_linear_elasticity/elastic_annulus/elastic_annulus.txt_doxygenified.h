/**

\mainpage The equations of time-harmonic linear elasticity 

The aim of this tutorial is to demonstrate the solution of the
time-harmonic equations of linear elasticity in cartesian
coordinates. These equations are useful to describe forced, time-harmonic
oscillations of elastic bodies.



<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
\n\n
The implementation of the equations and the documentation 
were developed jointly with David Nigro.
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>

<HR>
<HR>

\section theory Theory
Consider a linearly elastic body (of density \f$ \rho \f$, Young's
modulus \f$ E \f$ and Poisson's ratio \f$ \nu \f$), occupying 
the region \f$ D \f$
whose boundary is \f$ \partial D \f$.  Assuming that the body performs
time-harmonic oscillations of frequency of \f$ \omega \f$ 
its motion is governed by the equations of time-harmonic 
linear elasticity
<CENTER>
\f[
\pmb{\nabla}^*\cdot\pmb{\tau}^*+
\rho \mathbf{F}^*=-\rho\omega^2 \mathbf{u}^*,
\f]
</CENTER> 
where the \f$ x_{i}^* \f$ are the cartesian coordinates, and 
the time-periodic stresses, body force and displacements are given by 
\f$ {\rm Re}\{\pmb{\tau^*}(x_{i}^*){\rm e}^{-{\rm i}\omega
t^*}\} \f$, 
\f$ {\rm Re}\{\mathbf{F}^*(x_{i}^*){\rm e}^{-{\rm i}\omega
t^*}\} \f$ and \f$ {\rm Re}\{\mathbf{u}^*(x_{i}^*){\rm e}^{-{\rm i}\omega
t^*}\} \f$ respectively.  Note that, as usual, the superscript
asterisk notation is used to distinguish dimensional quantities
from their non-dimensional counterparts where required.

The body is subject to imposed time-harmonic displacements
\f$  {\rm Re}\{\mathbf{\hat{u}}^* {\rm e}^{-{\rm i}\omega
t^*}\} \f$
along  \f$ \partial D_d \f$, and to an imposed time-harmonic traction
\f$  {\rm Re}\{ \pmb{\hat{\tau}}^* {\rm e}^{-{\rm i}\omega
t^*}\} \f$ along \f$ \partial D_n \f$
where \f$ \partial D=\partial D_d\cup\partial D_n \f$. This requires that
<CENTER>
\f[
\mathbf{u}^*=\mathbf{\hat{u}}^*\,\,\textrm{on $
\partial D_d $ },\quad
\pmb{\tau}^* \cdot 
{\mathbf{n}}=\pmb{\hat{\tau}}^*\,\,\textrm{on $
\partial D_n $ }
\f]
</CENTER>
where \f$ {\bf n} \f$ is the outer unit normal on the boundary.

The stresses and displacements are related by the constitutive equations
<CENTER>
\f[
\pmb{\tau}^*=\frac{E}{1+\nu}\left(
\frac{\nu}{1-2\nu}(\pmb{\nabla}^*\cdot\mathbf{u}^*)\textbf{I}+
\frac{1}{2}(\pmb{\nabla}^*\mathbf{u}^*+\pmb{\nabla}^*\mathbf{u}^{*{\rm T}})\right),
\f]
</CENTER>
where \f$ \pmb{\nabla}^*\mathbf{u}^{*{\rm T}}
\f$ represents the transpose of \f$ \pmb{\nabla}^*\mathbf{u}^*
\f$.

We non-dimensionalise the equations, using a problem specific reference
length, \f$ {\cal L} \f$, and a timescale \f$ \displaystyle {\cal T}=\frac{1}{\omega} \f$, and  use
Young's modulus to non-dimensionalise the body force and the stress tensor:
<CENTER>
\f[
\pmb{\tau}^* = E \, \pmb{\tau}, \qquad
x_{i}^* = {\cal L}\, x_{i}
\f]
\f[
\mathbf{u}^* = {\cal L}\, \mathbf{u}, \qquad
\mathbf{F}^* = \frac{E}{\rho \cal L} \, \mathbf{F},
\qquad
t^* = {\cal T}\, t.
\f]
</CENTER>

The non-dimensional form of the equations is then given by
<CENTER>
\f[
\pmb{\nabla}\cdot\pmb{\tau}+\mathbf{F}=-\Omega^2\mathbf{u},
\ \ \ \ \ \ \ \ \ \ (1)
\f]
</CENTER>
with the non-dimensional constitutive relation,
<CENTER>
\f[
\pmb{\tau}=\frac{1}{1+\nu}\left(
\frac{\nu}{1-2\nu}(\pmb{\nabla}\cdot\mathbf{u})\textbf{I}+
\frac{1}{2}(\pmb{\nabla}\mathbf{u}+\pmb{\nabla}\mathbf{u}^{{\rm T}})\right).
\ \ \ \ \ \ \ \ \ \ (2)
\f]
</CENTER>
The non-dimensional parameter
<CENTER>
\f[
\Omega = {\cal L}\omega \sqrt{\frac{\rho}{E}}
\f]
</CENTER>
is the ratio of the elastic body's intrinsic timescale,
 \f$\displaystyle {\cal L} \sqrt{\frac{\rho}{E}}\f$, to the problem-specific 
timescale,\f$ \displaystyle {\cal T}=\frac{1}{\omega} \f$, that we 
used to non-dimensionalise time. \f$ \Omega \f$ can be interpreted as a
non-dimensional version of the excitation frequency; alternatively/equivalently 
\f$ \Omega^2 \f$ may be interpreted as a non-dimensional density.
The boundary conditions are
<CENTER>
\f[
\mathbf{u}=\mathbf{\hat{u}}\,\,\textrm{on $\partial D_d $ },\quad
\pmb{\tau} \cdot 
{\mathbf{n}}=\pmb{\hat{\tau}}\,\,\textrm{on $\partial D_n $ }.
\f]
</CENTER>



<HR>
<HR>



\section element_types Implementation 

Within \c oomph-lib, the non-dimensional version of equations 
(1) with 
the constitutive equations (2) are implemented 
in the \c TimeHarmonicLinearElasticityEquations<DIM> equations 
class, where the template parameter \c DIM indicates the spatial dimension. 
Following our usual approach, discussed in the
<A HREF="../../../quick_guide/html/index.html"> (Not-So-)Quick
Guide,</A> this equation class is then combined with a geometric finite
element to form a fully-functional finite element. 
For instance, the combination of the 
\c TimeHarmonicLinearElasticityEquations<2>
class with the geometric finite element \c QElement<2,3> yields a nine-node 
quadrilateral element. As usual, the mapping 
between local and global (Eulerian) coordinates within an element is given by,
<CENTER>
\f[
x_i = \sum_{j=1}^{N^{(E)}} X^{(E)}_{ij} \, \psi_j, \qquad
i=1,2,
\f]
</CENTER> 
where \f$ N^{(E)} \f$ is the number of nodes in the element, \f$ X^{(E)}_{ij}
\f$ is the \f$ i \f$-th global (Eulerian) coordinate of 
the \f$ j \f$-th \c Node 
in the element, and the \f$ \psi_j \f$ are 
the element's shape functions, defined in the geometric
finite element. 

All the constitutive parameters are real. The two
components of the displacement field have a real and imaginary part.
We store the four real-valued nodal unknowns in the order
\f$ {\rm Re}\{u_x\}, {\rm Re}\{u_y\}, {\rm Im}\{u_x\}, {\rm Im}\{u_y\}
\f$ and use the shape functions to interpolate the displacements as
<CENTER>
\f[
u_i^{(n)} = \sum_{j=1}^{N^{(E)}} U^{(E)}_{ij} \, \psi_j, \qquad 
i=1,...4,
\f]
</CENTER>
where \f$  U^{(E)}_{ij} \f$  is the \f$ i \f$-th displacement component 
(enumerated as described above) at the \f$ j \f$-th \c Node in the element.


<HR>
<HR>

\section test A test problem: Oscillations of an elastic annulus

We consider the time-harmonic axisymmetric deformation of
a 2D annular elastic body that occupies the region \f$ r_{\rm min}\leq
r\leq r_{\rm max}, 0\leq\theta\leq 2\pi \f$, shown in the left half
of the sketch below. We impose a constant-amplitude
axisymmetric displacement \f$ {\mathbf{u}}(r_{\rm min},\theta)=
u_0 (1-{\rm i}) {\mathbf{e}}_r\f$ on the inner boundary 
and a constant-amplitude pressure load 
\f$ \pmb{\hat{\tau}} = -P_0 {\mathbf{e}}_r \f$  on the outer boundary.
(\f$ {\mathbf{e}}_r \f$ is the unit vector in the radial
direction).


\image html sketch.gif "Sketch of the test problems. " 
\image latex sketch.eps "Sketch of the test problems. " width=0.8\textwidth


It is easy to find an analytical solution of this problem by working
in polar coordinates and exploiting the axisymmetry of the solution by
writing the displacement as
\f$ \displaystyle {\mathbf{u}} = U(r) {\mathbf{e}}_r \f$. 
The radial displacement  \f$ U(r) \f$ is then governed by 
<CENTER>
\f[
\frac{d}{d r}\left(\frac{U}{r}+\frac{dU}{dr}\right) + k^2 U=0,
\f]
</CENTER>
where \f$ \displaystyle k^2= \frac{\Omega^2}{\lambda+2\mu}\f$ and
<CENTER>
\f[
\lambda = \frac{\nu}{(1+\nu)(1-2\nu)} \qquad \textrm{and} 
\qquad \mu = \frac{1}{2(1+\nu)}
\f]
</CENTER>
are the non-dimensional Lame parameters.
The solution of this equation is given by:
<CENTER>
\f[
U(r) = aJ_1(kr)+bY_1(kr).
\f]
</CENTER>
where \f$ J_1 \f$ and \f$ Y_1 \f$ are Bessel functions
of the first and second kind, respectively.
The amplitudes \f$ a \f$ and \f$ b \f$ can be found using the 
boundary conditions on \f$ r_{\rm min} \f$ and \f$ r_{\rm max} \f$. 
In the driver code discussed below, the (lengthy) expressions 
for \f$ a \f$ and  \f$ b \f$ in terms
of the problem parameters can be found in 
the <code> GlobalParameters::exact_u() </code> function.

We note that even though a relatively simple analytical solution
(in polar coordinates!) exists for this problem, it is a non-trivial
test case for our code which solves the governing equations in cartesian
coordinates. However, to show that we can also compute non-trivial
solutions, we also consider the case where the annular region
has a "gap" and therefore occupies only a fraction (90%) of the circumference.
This creates two additional boundaries (the radial lines bounding
the "gap" and we subject these to the same pressure
that acts on the outer boundary, as shown in the right half of the
sketch above.

<HR>
<HR>


\section results Results
The figures below show "carpet plots" of the real and imaginary 
parts of the exact (green) and computed (red) horizontal displacement 
(\f$ u_1 \f$ -- the plots for \f$ u_2 \f$ obviously look very similar)
for the continuous coating and \f$ \Omega^2=10 \f$, \f$ \nu=0.3 \f$, 
\f$ P_0 = 0.3 \f$, \f$ r_{\rm min}=1\f$ and \f$ r_{\rm max}=1.5 \f$.

\image html re_ux.gif "Real part of the horizontal displacement. Green: exact; red: computed. " 
\image latex re_ux.eps "Real part of the horizontal displacement. Green: exact; red: computed. " width=0.8\textwidth

\image html im_ux.gif "Imaginary part of the horizontal displacement. Green: exact; red: computed. " 
\image latex im_ux.eps "Imaginary part of the horizontal displacement. Green: exact; red: computed. " width=0.8\textwidth


To demonstrate that the resulting displacement field is indeed
axisymmetric, here is a plot of the real part of the radial
displacement, \f$ ( Re(u_1)^2 + Re(u_2)^2 )^{1/2} \f$.

\image html re_ur.gif "Real part of the radial displacement. Green: exact; red: computed. " 
\image latex re_ur.eps "Real part of the radial displacement. Green: exact; red: computed. " width=0.8\textwidth

Finally, here is a plot of the real part of the horizontal
displacement for the case when there is a 10% "gap" in the annular
region. The presence of the gap clearly breaks the axisymmetry of the solution 
and creates waves that propagate in all directions:

\image html re_ux_hole.gif "Real part of the horizontal displacement for an incomplete annular region. " 
\image latex re_ux_hole.eps "Real part of the horizontal displacement for an incomplete annular region. " width=0.8\textwidth

Unsurprisingly, we are not aware of an analytical solution for this
problem. 

<HR> 
<HR>
 

\section namespace Global parameters and functions
As usual, we define all non-dimensional parameters in a namespace
where we also define the displacement to be applied on 
the inner boundary, the traction (pressure) to be applied on the 
outer boundary, and the exact solution (which we don't list here
because it is very lengthy).

\dontinclude time_harmonic_elastic_annulus.cc
\skipline start_namespace
\until end_of_pressure_load

We also define the output directory and the number of
elements in the mesh.

\skipline Output
\until Nr

 
<HR>
<HR> 

\section main The driver code

We start by defining command line arguments which specify the
number of elements in the mesh and indicate the presence or absence of the
"gap" in the coating. 

\dontinclude time_harmonic_elastic_annulus.cc
\skipline start_of_main
\until have_gap

The code performs a parameter study in which we compute the solution
for a range of pressures. We specify the pressure increment and
the number of steps to be performed, parse the command line arguments 
and document them.

\skipline P increment
\until doc_specified 

Next, we create the problem (discretising the domain with nine-noded
quadrilateral elements), 
\skipline Set up the problem
\until AnnularDisk

and perform the parameter study:

\skipline Initial values
\until {

\skipline Solve the problem using Newton's method
\until newton_solve();

\skipline Doc solution
\until end of main
 
<HR>
<HR>

\section problem The problem class 
The \c Problem class is very simple. As in other problems with Neumann 
boundary conditions, we provide separate meshes for the "bulk" 
elements and the face elements that apply the traction 
boundary conditions.  The latter are attached to the relevant 
faces of the bulk elements by the function \c create_traction_elements().

\dontinclude time_harmonic_elastic_annulus.cc
\skipline begin_problem
\until actions_before_newton_solve()

\skipline Doc the solution
\until delete_traction

\skipline Pointer to solid mesh
\until Solid_mesh_pt

\skipline Pointer to mesh of traction
\until };



<HR>
<HR>

\section constructor The problem constructor
We begin by building the "bulk" solid mesh, specifying the 
presence of the gap (and its width) if necessary. If there
is no gap, the mesh is periodic; see \ref comments for a more
detailed discussion of how the mesh for this problem is constructed.

\dontinclude time_harmonic_elastic_annulus.cc
\skipline start_of_constructor
\until }

\skipline Build solid mesh
\until H_annulus

It is always a good idea to have a look at the mesh
and its boundary enumeration:

\skipline Let's have
\until solid_mesh_boundary

We loop over the elements and specify the relevant physical parameters
via the pointer to the tensor of constitutive parameters and 
the (square of the) non-dimensional frequency,
\until }

Next we create the mesh that contains the FaceElements that apply 
the traction boundary conditions, and combine all meshes into a
global mesh:

\until build_global_mesh()

We apply the displacement boundary conditions by pinning
the real and imaginary part of the two displacement components
at all nodes on boundary 0. We then impose a purely
radial displacement on the real and imaginary parts of the
displacement field, using the function defined in the \c Global_Parameters
namespace:

\skipline Solid boundary conditions
\until }

Finally we assign the equation numbers and specify the output directory:

\until end_of_constructor
 

<HR>
<HR>


\section traction_elements The traction elements
We create the face elements that apply the traction on the
outer boundary (boundary 2). If there is a gap in the
annular region we also apply the pressure loading on boundaries 
1 and 3.

\dontinclude time_harmonic_elastic_annulus.cc
\skipline start_of_create_traction_elements
\until end_of_create_traction_elements


<HR>
<HR>

\section doc Post-processing
As expected, this member function documents the computed solution.

\dontinclude time_harmonic_elastic_annulus.cc
\skipline start_doc
\until end doc

<HR>
<HR>
 

\section conclusion Comments and Exercises

\subsection comments Comments

- We did not discuss the mesh generation in detail here: 
  The mesh is created by straightforward overloading
  of \c oomph-lib's existing rectangular quad mesh -- the constructor
  simply adjusts the nodal positions (exactly as suggested in the
  <A HREF="../../../quick_guide/html/index.html#distorted_mesh">
  "Not-so-quick" guide</a>.) A little bit of extra work is required
  to enforce periodicity on the mesh for the case without a gap
  in the annular region because two of the boundaries in 
  the original mesh then overlap. How this is dealt with is discussed in 
  <A  HREF="../../../navier_stokes/rayleigh_channel/html/index.html#periodic">
  another tutorial.</a>
  \n\n
- If you inspect the <A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/time_harmonic_elastic_annulus.cc">driver code</a>
  you will notice that it also contains relevant code to perform 
  spatially adaptive simulations of the problem -- the adaptive
  version of the code is selected with \c \#ifdefs. Dealing with the
  periodic boundary conditions for spatially adaptive meshes
  requires a few additional steps, but they are also discussed 
  <A HREF="../../../linear_elasticity/refineable_periodic_load/html/index.html#constructor">elsewhere,</a>
  so we won't discuss them here.
.


\subsection exercises Exercises
- Change the parameter study performed by the driver code such that
  the loop varies the frequency parameter \f$ \Omega^2 \f$. Assess the effect
  of an increase in \f$ \Omega^2 \f$ on the accuracy of the 
  solution by comparing the computed results against 
  the exact solution at fixed spatial resolution.
  \n\n
- Explore the use of spatial adaptivity in the problem
  \n\n
  - for the same parameter study suggested in the previous exercise
    (increase in \f$ \Omega^2 \f$ ) 
  .
  and 
  - for the problem with the "gap" in the annular region.
  .
  \n\n
- Modify the code to exploit at least some of the problem's symmetry,
  e.g. by solving the problem for the continuous annulus in a
  quarter of the original domain, \f$ x_1, x_2 \ge 0 \f$, say, using 
  appropriate symmetry boundary conditions along the coordinate
  axes.
.
<HR>
<HR> 


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/time_harmonic_elastic_annulus.cc">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/time_harmonic_elastic_annulus.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


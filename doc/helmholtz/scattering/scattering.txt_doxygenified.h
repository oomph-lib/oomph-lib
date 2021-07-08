/**

\mainpage Example problem: The Helmholtz equation -- scattering problems


In this document we discuss the finite-element-based solution of the 
Helmholtz equation, an elliptic PDE that describes time-harmonic wave propagation 
problems. We start by  reviewing the relevant theory and then 
present the solution of a simple model problem -- the scattering of a
planar wave from a circular cylinder.


<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
This tutorial and the associated driver codes were developed jointly
with Tarak Kharrat (EnstaParisTech, Paris).
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>

<HR>
<HR>

\section theory Theory: The Helmholtz equation for time-harmonic scattering problems

The Helmholtz equation governs time-harmonic solutions of problems
governed by the linear wave equation
\f[
\nabla^2 U(x,y,t) = \frac{1}{c^2} \frac{\partial^2 U(x,y,t)}{\partial t^2},
 \ \ \ \ \ \ \ \ \ \ \ \ (1)
\f]
where \f$ c \f$ is the wavespeed. Assuming that \f$ U(x,y,t) \f$
is time-harmonic, with frequency \f$ \omega \f$, we write the real
function \f$ U(x,y,t) \f$ as
\f[
U(x,y,t) =Re (u(x,y) \ e^{-i \omega t})
\f]
where \f$ u(x,y) \f$ is complex-valued. This transforms 
(1) into the Helmholtz equation
\f[
\nabla^2 u(x,y) + k^2 u(x,y) = 0
 \ \ \ \ \ \ \ \ \ \ \ \ (2)
\f]
where 
\f[
k = \frac{\omega}{c}
 \ \ \ \ \ \ \ \ \ \ \ \ (3)
\f]
is the wave number. Like other elliptic PDEs the Helmholtz equation
admits Dirichlet, Neumann (flux) and Robin boundary conditions.

If the equation is solved in an infinite domain (e.g. in scattering
problems) the solution must satisfy the so-called
<a href="http://en.wikipedia.org/wiki/Sommerfeld_radiation_condition">
Sommerfeld radiation condition</a> which in 2D has the form
\f[
\lim_{r\to \infty} \sqrt{r} \left(\frac{\partial u}{\partial r} - iku
\right) =0. 
\f]
Mathematically, this conditions is required to ensure the uniqueness
of the solution (and hence the well-posedness of the problem).
In a physical context, such as a scattering problem, the condition 
ensures that scattering of an incoming wave only produces outgoing not 
incoming waves from infinity.

<HR>
<HR>

\section discr Discretisation by finite elements
The discretisation of the Helmholtz equation itself only requires a trivial
modification of \c oomph-lib's Poisson elements -- we simply add
the term  \f$ k^2 u \f$ to the residual. Since most practical applications
of the Helmholtz equation involve complex-valued solutions, we provide 
separate storage for the real and imaginary parts of the solution -- 
each \c Node therefore stores two unknowns values. By default,
the real and imaginary parts are stored as values 0 and 1,
respectively; see the section \ref numbering for details.

The application of Dirichlet and Neumann boundary conditions is 
straightforward and follows the pattern employed for the solution
of the Poisson equation: 
- Dirichlet conditions are imposed by pinning the relevant nodal
  values and setting them to the appropriate prescribed values. 
- Neumann (flux) boundary conditions are imposed via 
  \c FaceElements (here the \c HelmholtzFluxElements). 
   <a href="../../../poisson/two_d_poisson_flux_bc/html/index.html">
  As usual</a> we attach these to the faces of the "bulk" elements
  that are subject to the Neumann boundary conditions.
.

The imposition of the Sommerfeld radiation condition for problems in
infinite domains is slightly more complicated. In the following 
discussion we will restrict ourselves to two dimensions and assume 
that the infinite domain is truncated at a circular artificial
boundary \f$ \Gamma \f$ of radius \f$ R. \f$ [This assumption is also 
made in the implementation of \c oomph-lib's \c FaceElements
that allow the (approximate) imposition of the Sommerfeld radiation
condition. The methodology can easily be modified to deal with other 
geometries but this has not been done yet -- any volunteers?]
All methods exploit the fact that the relevant solution of the 
Helmholtz equation can be written in polar coordinates as
\f[
u(r,\varphi) = \sum_{n=-\infty}^{+\infty} A_n \ H_n^{(1)}(kr)  
\ e^{i n \varphi}, \ \ \ \ \ \ \ (4)
\f]
where the \f$A_n \f$ are suitable coefficients and
\f$ H_n^{(1)}(r) \f$ is the \f$ n \f$-th-order 
Hankel function of the first kind.


<HR>

\subsection ABCs Approximate/absorbing boundary conditions (ABCs)
It is possible to derive approximate versions of the Sommerfeld
radiation condition in which the normal derivative of the
solution on the artificial boundary is related to its
value and possibly its tangential derivatives. Such boundary
conditions (sometimes referred to as approximate or absorbing 
boundary conditions -- ABCs) are typically derived from asymptotic expansions
of the solution at large distances from the origin and become
more accurate the larger the radius \f$ R \f$ of the artificial boundary
\f$ \Gamma \f$ is. Higher accuracy can therefore only be achieved 
by increasing the size of the computational domain, with an 
associated increase in computational cost.

 \c oomph-lib provides an implementation of the following three 
boundary conditions (all taken from
J. J. Shirron & I. Babuska's paper "A comparison of approximate boundary 
conditions and infinite element methods for exterior Helmholtz
problems", Computer Methods in Applied Mechanics and Engineering
<b>164</b> 121-139 (1998), in which the authors
compare the accuracy of these and many other approximate boundary
conditions).

- <b>Feng's first order ABC:</b>
\f[
 \frac {\partial u}{\partial n}-\bigg( ik - \frac{1}{2R}\bigg) u
 =0   \mbox{\ \ \  \ \  on\ }\Gamma  
\f] 
(This is identical to the first-order Bayliss and Turkel boundary condition).
- <b>Feng's second order ABC:</b>
\f[
 \frac{\partial u}{\partial n}-\bigg[ik-\frac
 {1}{2R}+\frac{i}{8kR^{2}}\bigg(1+4 \frac {\partial^{2}}{\partial
 \varphi^{2}}\bigg)\bigg]u
=0  \mbox{\ \ \  \ \  on\ }\Gamma  
\f]
- <b>Feng's third order ABC:</b>
\f[
 \frac {\partial u}{\partial n}-\bigg[ik-\frac
 {1}{2R}+\frac{1}{8k^{2}R^{2}}\bigg(ik+ \frac {1}{R}\bigg)
 \bigg(1+4 \frac {\partial^{2}}{\partial \varphi^{2}}\bigg)\bigg]u
=0   \mbox{\ \ \  \ \  on\ }\Gamma  
\f]
.
All three boundary conditions are implemented in the class
\c HelmholtzAbsorbingBCElement. The order of the approximation 
can be set via the member function
\c HelmholtzAbsorbingBCElement::abc_order(). All three boundary
conditions are local (relating the function to its normal derivative)
and do therefore not change the sparsity of the resulting finite
element equations.

<HR>

\subsection DtN The Dirichlet-to-Neumann mapping (DtN)
Using (4), it is easy to show (see, e.g., 
J. Jin "The Finite Element Method in Electromagnetics (second
edition)", Wiley (2002) p. 501ff -- but note that Jin assumes that
the potential varies like \f$ \exp(i\omega t)\f$ rather than
\f$ \exp(-i\omega t)\f$ as assumed here) that the 
normal (radial) derivative, \f$ \partial u / \partial n =   
\partial u / \partial r, \f$ on the artificial boundary 
\f$ \Gamma \f$ is given by
\f[
 \frac {\partial u}{\partial r}\bigg|_{r=R} =
 \frac {\partial u}{\partial n}\bigg|_{r=R} =
 \gamma (u) \ \ \ \ \ \ \ (5)   
\f]
where 
\f[
\gamma (u) = \frac {k}{2 \pi} \sum_{n=-\infty}^{+\infty} 
\frac {H_n^{(1)^{'}}(kR)}{H_n^{(1)}(kR)} \quad 
\int_0^{2\pi}u(R,\varphi^{'}) \ e^{in(\varphi-\varphi^{'})}
\,d\varphi^{'}.  
 \ \ \ \ \ \ \ (6)   
\f]
Equation (5) again provides a condition on the
normal derivative of the solution along the artificial boundary 
and is implemented in the \c HelmholtzDtNBoundaryElement class.
Since \f$ \gamma \f$ depends on the solution everywhere 
along the artificial boundary (see (6)), the application of the 
boundary condition (5) introduces a non-local coupling
between all the degrees of freedom located on that boundary.
This is handled by classifying the unknowns that affect \f$ \gamma \f$ 
but are not associated with the element's own nodes as external \c Data. 

To facilitate the setup of the interaction between the 
\c HelmholtzDtNBoundaryElements, \c oomph-lib provides the
class \c HelmholtzDtNMesh which provides storage for 
(the pointers to) the  \c HelmholtzDtNBoundaryElements 
that discretise the artificial boundary. The 
member function \c HelmholtzDtNMesh::setup_gamma() pre-computes the 
\f$ \gamma \f$ values required for the imposition of equation (5).
The radius \f$ R \f$ of the artificial boundary and the (finite) number of
(Fourier) terms used in the sum in (6) are specified
as arguments to the constructor of the \c HelmholtzDtNMesh.

<b>NOTE:</b> Since \f$ \gamma \f$  depends on the solution, 
it must be recomputed whenever the unknowns are updated during
the Newton iteration. This is best done by adding a call to
\c HelmholtzDtNMesh::setup_gamma() to 
\c Problem::actions_before_newton_convergence_check().
[If Helmholtz's equation is solved in isolation (or within
a coupled, but linear problem), Newton's method will converge 
in one iteration. In such cases the unnecessary recomputation of 
\f$ \gamma \f$ after the one-and-only Newton iteration can be 
suppressed by setting \c Problem::Problem_is_nonlinear to 
\c false.]

<HR>
<HR>

\section scattering A specific example: Scattering of an acoustic wave from a sound-hard obstacle
We will now demonstrate the methodology for a specific example:
the scattering of sound waves in an acoustic medium 
of density \f$ \rho  \f$ and bulk modulus \f$ B \f$.
Assuming that an incoming sound wave impacts a rigid, impermeable
obstacle as shown in this sketch,

\image html scattering.gif "Scattering of an incoming wave from a sound-hard obstacle -- the scatterer. " 
\image latex scattering.eps "Scattering of an incoming wave from a sound-hard obstacle -- the scatterer. " width=0.6\textwidth

we wish to find the wave field that is scattered from the body.

For this purpose we denote the time-dependent displacement
of the fluid particle in the acoustic medium by 
\f$ {\bf u}^*(x^*,y^*,t^*) \f$ and introduce a displacement
potential \f$ \Phi^*(x^*,y^*,t^*) \f$ such that
\f[
{\bf u}^* = \nabla^* \Phi^*.
\f]
(As usual we employ asterisks to distinguish dimensional quantities
from their non-dimensional equivalents, to be introduced below.) 
It is easy to show that \f$ \Phi^* \f$ satisfies the linear wave
equation (1) with wave speed \f$ c = \sqrt{B/\rho}\f$.
 
Since the surface \f$ \partial D_{bound}\f$ of the scatterer is
impenetrable, the normal displacement of the fluid has to vanish
on  \f$ \partial D_{bound}\f$  and the boundary condition for the 
displacement potential becomes
\f[
\left. \frac{\partial \Phi^*}{\partial n^*}\right|_{\partial D_{bound}} = 0.
\ \ \ \ \ \ \ \ \ \ \ \ \  (7)
\f]



We non-dimensionalise all lengths and displacements on some problem-dependent
lengthscale \f$ {\cal L}\f$ (e.g. the radius of the scatterer),
non-dimensionalise the potential as \f$ \Phi^* =  a^2 \Phi \f$ and
scale time on the period of the oscillation,
\f$ t^* = \frac{2\pi}{\omega} t.\f$ The governing equation
then becomes
\f[
\nabla^2 \Phi + k^2 \Phi = 0, \ \ \ \ \ \ \ \ \ \ \ (8)
\f]
where the square of the wavenumber  is given by
\f[
k^2  = \frac{\rho (a\omega)^2}{B}.
\f]

Assuming that the incoming wave (already satisfying (8))
is described by a (known) non-dimensional displacement potential of the form
\f[
\Phi_{inc}(x,y,t) = \phi_{inc}(x,y) \  e^{-i 2\pi t},
\f]
we write the total potential as
\f[
\Phi(x,y,t) = \bigg( \phi_{inc}(x,y) + u(x,y) \bigg) \ e^{-i 2\pi t},
\f]
where \f$ u(x,y) \  e^{-i 2\pi t} \f$ represents the displacement potential 
associated with the scattered field which must satisfy (2). 
The boundary condition (7) then becomes a Neumann 
(flux) boundary condition for the scattered field,
\f[
\left. \frac{\partial u}{\partial n}\right|_{\partial D_{bound}} =
- \left. \frac{\partial \phi_{inc}}{\partial n}\right|_{\partial D_{bound}}.
\ \ \ \ \ \ \ \ \ \ \ \ \  (9)
\f]

For the special case of the incoming wave being a planar wave,
propagating along the x-axis, the incoming field can be written 
in polar coordinates as
\f[
  \phi_{inc}(r, \varphi)  
= \sum_{n=-\infty}^{+\infty} i^n  J_n(kr) e^{in\varphi}, 
\f]
where \f$ J_n \f$ is the Bessel function of the first kind
of order \f$ n \f$. The exact solution for the scattering
of such a wave from a circular disk is given by the series
\f[
  u_{ex}(r,\varphi) 
= -\sum_{n=-\infty}^{+\infty} i^n \frac
{H^{'}_{n}(k)}{J^{'}_{n}(k)}  H_n(kr) e^{in\varphi},
\ \ \ \ \ \ \ \ \ \ \ \ \  (10)
\f]
where we have chosen the disk's radius, \f$ a \f$,  as the
lengthscale by setting \f$ {\cal L} = a\f$. In the above expression,
\f$ H_n \f$ denotes the Hankel function of the first kind
of order \f$ n \f$ and the prime denotes differentiation with respect 
to the function's argument.

A quantity that is of particular interest in wave propagation
problems is the time-average of the power radiated by the scatterer,
\f[
\overline{\cal P}^* = 
\frac{\omega}{2\pi} \int_{0}^{2\pi/\omega} {\cal P}^*(t) \ dt^*.
\f]
In the context of an acoustic wave, the total instantaneous power,
\f$ {\cal P}^*(t), \f$ radiated over a closed boundary is
\f[
{\cal P}^*(t) = \oint \frac{\partial {\bf u^*} }{\partial t^*} 
\cdot p^* {\bf n} \ dS^*,
\f]
where the pressure is related to the displacement potential via
\f[
p^* = \rho \omega^2 \Phi^*.
\f]
The non-dimensional time-averaged radiated power can be expressed in 
terms of the complex potential \f$ \phi \f$ as
\f[
\overline{\cal P} = \frac{\overline{\cal P}^*}{\rho \omega^3 {\cal L}^4} =
\frac{1}{2}
\oint\bigg[Im\bigg(\frac{\partial \phi}{\partial n}\bigg) \ Re(\phi) -
           Re\bigg(\frac{\partial \phi}{\partial n}\bigg) \ Im(\phi)
           \bigg]
 \ dS.
\f]
<HR>
<HR>

\section results Results
The figure below shows an animation of the displacement potential
\f$ Re(u(x,y,t)) \f$ for scattering from a circular disk 
for a non-dimensional wavenumber of \f$ k=1 \f$ over one 
period of the oscillation. The simulation was performed in an 
annular computational domain, bounded by the outer surface the 
(unit) disk and an artificial outer boundary of non-dimensional 
radius \f$ R=1.5. \f$ The Sommerfeld radiation
condition was imposed using the DtN mapping and the simulation
was performed with spatial adaptivity (note the non-uniform 
refinement).

The "carpet plot" compares the exact (green) and computed (red)
solutions for the displacement potential. The colours in the contour 
plot at the bottom of the figure provide an alternative
visualisation of the magnitude of the scattered field. 


\image html scattering_animation.gif "The displacement potential associated with the scattered wave, animated over one period of the oscillation. " 
\image latex scattering_animation.eps "The displacement potential associated with the scattered wave, animated over one period of the oscillation. " width=0.6\textwidth

<HR>
<HR>

\section num_soln The numerical solution

\subsection namespace The global namespace

As usual, we define the problem parameters in a global namespace.
The main physical parameter is the (square of the) wave number,
\f$ k^2 \f$. \c N_fourier is the number of
(Fourier) terms  to be used in evaluation of the 
series in equations (6) and (10).
The remaining parameters determine how the Sommerfeld radiation
condition is applied. 

\dontinclude scattering.cc
\skipline start_of_namespace
\until std::complex<double> I(0.0,1.0);  


The function \c get_exact_u
returns the exact solution for the scattering problem.
We will use this function for the validation of our results.

\dontinclude scattering.cc
\skipline  Exact solution for scattered
\until  }// end of get_exact_u
  
Next we provide a function that computes the prescribed flux 
(normal derivative) of the solution, \f$ \partial u/\partial n =
-\partial \phi_{inc}/\partial n \f$, evaluated on the surface of the
unit disk.

\dontinclude scattering.cc
\skipline Flux (normal derivative)
\until // end of namespace

<HR>
<HR>

\subsection main The driver code
The driver code is very straightforward. We parse the command line
to determine which boundary condition to use and set the 
flags in the global namespace accordingly.

\dontinclude scattering.cc
\skipline start_of_main
\until }

Next we build the problem, either with or without enabling spatial
adaptivity and define the output directory.

\until RESLT

Finally, we solve the problem and document the results.

\until } //end of main

 
<HR>
<HR>

\subsection class The problem class

The problem class is very similar to that employed for 
the <a href="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
adaptive solution of the 2D Poisson equation with flux boundary
conditions.</a> The only difference is that we provide 
two separate meshes of \c FaceElements: one for the 
inner boundary where the \c HelmholtzFluxElements apply the
Neumann condition (9), and one for 
the outer boundary where we apply the (approximate) Sommerfeld 
radiation condition. As discussed in section \ref DtN , we use the function
\c actions_before_newton_convergence_check() to recompute the
\f$ \gamma \f$ integral whenever the unknowns are updated
during the Newton iteration.

\dontinclude scattering.cc
\skipline start_of_problem_class
\until }; // end of problem class


<HR>
<HR>

\subsection constr The problem constructor

We start by building the bulk mesh, using the refineable or
non-refineable version of the \c TwoDAnnularMesh, depending on the
macro \c ADAPTIVE. (The error tolerances for the adaptive version
are chosen such that the mesh is refined non-uniformly -- with the
default tolerances, \c oomph-lib's automatic mesh adaptation
procedure refine the mesh uniformly.)


\dontinclude scattering.cc
\skipline start_of_constructo
\until #endif

Next we create the two (empty) meshes for the \c FaceElements,

\dontinclude scattering.cc
\skipline  // Pointer to mesh containing the Helmholtz outer boundary condition
\until  new Mesh; 

and populate them using the functions \c create_flux_elements(...) and \c
create_outer_bc_elements(...).

\dontinclude scattering.cc
\skipline // Create prescribed 
\until outer_boundary_mesh_pt);

We add the various (sub-)meshes to the problem and build the global
mesh

\dontinclude scattering.cc
\skipline // Add the several  sub
\until build_global_mesh();

Finally, we complete the build of the various elements by
by passing pointers to the relevant quantities to them,
and assign the equation numbers. 
\dontinclude scattering.cc
\skipline  // Complete the build
\until } // end of constructor

The problem is now ready to be solved.

<HR>
<HR>

\subsection before_adapt Actions before adapt
The mesh adaptation is driven by the error estimates for the bulk
elements. The various \c FaceElements must therefore be removed from
the global mesh before the adaptation takes place. We do this by 
calling the function \c delete_flux_elements(...) for the two
face meshes, before rebuilding the Problem's global mesh.

\dontinclude scattering.cc
\skipline start_of_actions_before_adapt
\until }// end of actions_before_adapt


<HR>
<HR>

\subsection after_adapt Actions after adapt
After the (bulk-)mesh has been adapted, the flux elements must 
be re-attached. This is done by calling the functions
\c create_flux_elements(...) and \c create_outer_bc_elements,
followed by rebuilding the Problem's global mesh. Finally, 
we complete the build of the \c FaceElements by calling
the functions \c setup_outer_boundary() and
\c set_prescribed_incoming_flux_pt().

\dontinclude scattering.cc
\skipline start_of_actions_after_adapt
\until }// end of actions_after_adapt


<HR>
<HR>

\subsection delete Delete flux elements
The helper function \c delete_face_elements() is used
to delete all \c FaceElements in a given surface mesh
before the mesh adaptation.

\dontinclude scattering.cc
\skipline start_of_delete
\until } // end of delete_outer_face_elements

<HR>
<HR>


\subsection create_flux Creating the face elements
The functions \c create_flux_elements(...) and 
\c create_outer_bc_elements(...) create the \c FaceElements
required to apply the boundary conditions on the inner and
outer boundaries of the annular computational domain.
They both loop over the bulk elements that are adjacent to the appropriate
mesh boundary and attach the required \c FaceElements to their
faces. The newly created \c FaceElements are then added
to the appropriate mesh.

\dontinclude scattering.cc
\skipline start_of_create_outer
\until } // end of create_outer

(We omit the listing of the function \c create_flux_elements(...) because
it is very similar.  Feel free to inspect in the 
<A HREF="../../../../demo_drivers/helmholtz/scattering/scattering.cc">
source code.</A>)


<HR>
<HR>


\subsection doc Post-processing
The post-processing function \c doc_solution(...) computes and outputs 
the total radiated power, and plots the computed and exact solutions
(real and complex parts). 

\dontinclude scattering.cc
\skipline start_of_doc
\until Norm of solution

Finally, we create the data required to produce an animation of 
the actual (real) potential at 40 instants during a period of the
oscillation.

\until } // end of doc


<HR>
<HR>

\section comm_ex Comments and Exercises
\subsection numbering The enumeration of the unknowns
As discussed in the introduction, most practically relevant
solutions of the Helmholtz equation are complex valued. Since \c oomph-lib's 
solvers only deal with real (double precision) unknowns, the equations
are separated into their real and imaginary parts.
In the implementation of the Helmholtz elements, we store the real
and imaginary parts of the solution as two separate values at each 
node. By default, the real and imaginary parts are accessible via 
\c Node::value(0) and \c Node::value(1). However, to facilitate the use of the elements
in multi-physics problems we avoid accessing the unknowns
directly in this manner but provide the virtual function
\code
std::complex<unsigned> HelmholtzEquations<DIM>::u_index_helmholtz()
\endcode
which returns a complex number made of the two unsigneds that indicate
which nodal value represents the real and imaginary parts of the solution.
This function may be overloaded in combined multi-physics elements
in which a Helmholtz element is combined (by multiple inheritance) 
with another element, using the strategy described in 
<a href="../../../multi_physics/b_convection/html/index.html">
the Boussinesq convection tutorial</a>.

<HR>

\subsection ex Exercises
\subsubsection lin Exploiting linearity
Confirm that the (costly) re-computation of the 
\f$ \gamma \f$ integral in 
\c actions_before_newton_convergence_check() after the
first (and only) linear solve in the Newton iteration can be 
avoided by declaring the problem to be linear.

\subsubsection acc The accuracy of the boundary condition elements
Explore the accuracy (and computational cost) of the various
\c FaceElements that apply the Sommmerfeld radiation condition.
In particular, confirm that the accuracy of the DtN boundary 
condition is (nearly) independent of the radius of the artificial outer
boundary, whereas the accuracy of the ABC boundary condition can
only be improved by increasing the size of the computational domain.

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
  <CENTER>
  <A HREF="../../../../demo_drivers/helmholtz/scattering">
  demo_drivers/helmholtz/scattering/
  </A>
  </CENTER>
- The driver code is: 
  <CENTER>
  <A HREF="../../../../demo_drivers/helmholtz/scattering/scattering.cc">
  demo_drivers/helmholtz/scattering/scattering.cc
  </A>
  </CENTER>
.


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


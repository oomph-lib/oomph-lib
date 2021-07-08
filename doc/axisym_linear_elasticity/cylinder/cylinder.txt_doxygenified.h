/**

\mainpage The axisymmetric equations of linear elasticity 

The aim of this tutorial is to demonstrate the solution of the
axisymmetric equations of linear elasticity in cylindrical polar
coordinates.

<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
\n\n
This implementation of the equations and the documentation were developed
jointly with Matthew Russell with financial support to Chris Bertram 
from the Chiari and Syringomyelia Foundation.
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>

<HR>
<HR>

\section theory Theory
Consider a three-dimensional, axisymmetric body (of density \f$ \rho
\f$, Young's modulus \f$ E \f$, and Poisson's ratio \f$ \nu \f$),
occupying the region \f$ D \f$
whose boundary is \f$ \partial D \f$. Using cylindrical coordinates \f$
(r^*,\theta,z^*) \f$, the equations of linear elasticity can be written as
\f[
\nabla^* \cdot \bm{\tau}^* + \rho \mathbf{F}^* = \rho\frac{\partial^2\mathbf{u}^*}
{\partial {t^*}^2},
\f]
where \f$\nabla^* = \left(\frac{\partial}{\partial r^*},
\frac{1}{r^*}\frac{\partial}
{\partial \theta},\frac{\partial}{\partial z^*}\right)\f$,
\f$ \bm{\tau}^*(r^*,z^*,t^*) \f$ is the stress tensor,
\f$\mathbf{F}^*(r^*,z^*,t^*)\f$ is the body force and
\f$ \mathbf{u^*}(r^*,z^*,t^*) \f$ is the displacement field.

Note that, despite the fact that none of the above physical quantities depend
on the azimuthal angle \f$ \theta \f$, each can have a non-zero \f$ \theta \f$
component. Also note that variables written with a superscript asterisk are
dimensional, and their non-dimensional counterparts will be written without an
asterisk. (The coordinate \f$ \theta \f$ is, by definition, non-dimensional, so
it will always be written without an asterisk.)

A boundary traction \f$ \hat{\bm{\tau}}^* \f$ and boundary displacement \f$
\hat{\mathbf{u}}^* \f$ are imposed along the boundaries \f$ \partial
D_\mathrm{n} \f$ and \f$ \partial D_\mathrm{d} \f$ respectively, where \f$
\partial D = \partial D_\mathrm{d} \cup \partial D_\mathrm{n} \f$ so that
\f[
\mathbf{u}^* = \hat{\mathbf{u}}^*\text{ on }\partial D_\mathrm{d},
\qquad
\bm{\tau}^*\cdot\mathbf{n} = \hat{\bm{\tau}}^*\text{ on }\partial
D_\mathrm{n},
\f]
where \f$ \mathbf{n} \f$ is the outer unit normal vector.

The constitutive equations relating the stresses to the displacements are
\f[
\bm{\tau}^* =
\frac{E}{1+\nu}\left(\frac{\nu}{1-2\nu}\left(\nabla^*\cdot\mathbf{u}^*\right)
\mathbf{I} + \frac{1}{2}\left(\nabla^*\mathbf{u}^* +
\left(\nabla^*\mathbf{u}^*\right)^\mathrm{T}\right)\right),
\f]
where 
\f$ \mathbf{I} \f$ is the identity tensor and superscript \f$ \mathrm{T} \f$
denotes the transpose. In cylindrical coordinates, the matrix representation of
the tensor \f$ \nabla^*\mathbf{u}^* \f$ is
\f[
\nabla^*\mathbf{u}^* =
\begin{pmatrix}
\dfrac{\strut \partial u_r^*}{\strut \partial r^*} &
-\dfrac{\strut u_\theta^*}{\strut r^*} & \dfrac{\strut \partial
u_r^*}{\strut \partial z^*}\\
\dfrac{\strut \partial u_\theta^*}{\strut \partial r^*} &
\dfrac{\strut u_r^*}{\strut r^*} & \dfrac{\strut \partial
u_\theta^*}{\strut \partial z^*}\\
\dfrac{\strut \partial u_z^*}{\strut \partial r^*} & 0
& \dfrac{\strut \partial u_z^*}{\strut \partial z^*}
\end{pmatrix}
\f]
and \f$ \nabla^*\cdot \mathbf{u}^* \f$ is equal to the trace
of this matrix.

We non-dimensionalise the equations, using a problem specific reference
length, \f$ \mathcal{L} \f$, and a timescale \f$ \mathcal{T} \f$, and use
Young's modulus to non-dimensionalise the body force and the stress tensor:
\f[
\bm{\tau}^* = E \bm{\tau}, \qquad
r^* = \mathcal{L}r, \qquad
z^* = \mathcal{L}z
\f]
\f[
\mathbf{u}^* = \mathcal{L} \mathbf{u} \qquad
\mathbf{F}^* = \frac{E}{\rho \mathcal{L}}\mathbf{F},
\qquad
t^* = \mathcal{T}t.
\f]

The non-dimensional form of the axisymmetric linear elasticity equations is
then given by
\f[
\nabla\cdot \bm{\tau} + \mathbf{F} = \Lambda^2\frac{\partial^2\mathbf{u}}
{\partial t^2},\qquad\qquad (1)
\f]
where \f$ \nabla = \left(\frac{\partial}{\partial r},
\frac{1}{r}\frac{\partial}
{\partial \theta},\frac{\partial}{\partial z}\right) \f$,
\f[
\bm{\tau} =
\frac{1}{1+\nu}\left(\frac{\nu}{1-2\nu}\left(\nabla\cdot\mathbf{u}\right)
\mathbf{I} + \frac{1}{2}\left(\nabla\mathbf{u} +
\left(\nabla\mathbf{u}\right)^\mathrm{T}\right)\right),\qquad\qquad
(2)
\f]
and the non-dimensional parameter
\f[
\Lambda = \frac{\mathcal{L}}{\mathcal{T}}\sqrt{\frac{\rho}{E}}
\f]
is the ratio of the elastic body's intrinsic timescale,
\f$ \mathcal{L} \sqrt{\frac{\rho}{E}} \f$, to the problem-specific
timescale, \f$ \mathcal{T} \f$, that we use for time-dependent problems.
The boundary conditions are
\f[
\mathbf{u} = \hat{\mathbf{u}}\text{ on } \partial D_\mathrm{d}\qquad
\bm{\tau}\cdot\mathbf{n} = \hat{\bm{\tau}}\text{ on } \partial D_\mathrm{n}.
\f]
We must also specify initial conditions:
\f[
\left.\mathbf{u}\right|_{t=t_0} = \mathbf{u}^0,\qquad
\left.\frac{\partial\mathbf{u}}{\partial t}\right|_{t=t_0} =
\mathbf{v}^0.\qquad\qquad (3)
\f]

<HR>
<HR>

\section element_types Implementation 

Within \c oomph-lib, the non-dimensional version of the axisymmetric linear
elasticity equations (1) combined with the constitutive
equations (2) are implemented in the \c
AxisymmetricLinearElasticityEquations class. This class implements the
equations in a way which is general with respect to the specific element
geometry. To obtain a fully functioning element class, we must combine the
equations class with a specific geometric element class, as discussed in the
<A HREF="../../../quick_guide/html/index.html"> (Not-So-)Quick Guide</A>.
For example, we will combine the \c AxisymmetricLinearElasticityEquations class
with \c QElement<2,3> elements, which are 9-node quadrilateral elements, in our
example problem. As usual, the mapping between local and global (Eulerian)
coordinates within an element is given by
\f[
x_i = \sum_{j=1}^{N^{(E)}} X_{ij}^{(E)}\psi_j,\qquad i = 1,2,
\f]
where \f$ x_1 = r, x_2 = z \f$; \f$ N^{(E)} \f$ is the number of nodes in the
element. \f$ X_{ij}^{(E)} \f$ is the \f$ i \f$-th global (Eulerian) coordinate
of the \f$ j \f$-th node in the element and the \f$ \psi_j \f$ are the
element's shape functions, which are specific to each type of geometric
element.

We store the three components of the displacement vector as nodal data in the
order \f$ u_r, u_z, u_\theta \f$ and use the shape functions to interpolate the
displacements as
\f[
u_i = \sum_{j=1}^{N^{(E)}}U_{ij}^{(E)}\psi_{ij},\qquad i = 1,\dotsc,3,
\f]
where \f$ U_{ij} \f$ is the \f$ i \f$-th displacement component at the \f$ j
\f$-th node in the element, i.e., \f$ u_1=u_r, u_2=u_z,u_3=u_\theta \f$.

The solution of time dependent problems requires the specification of a \c
TimeStepper that is capable of approximating second time derivatives. In the
example problem below we use the Newmark timestepper.

<HR>
<HR>

\section test The test problem

As a test problem we consider forced oscillations of the circular cylinder shown
in the sketch below:
\image html figures/geom.gif "Azimuthal cross-section of the geometry. " 
\image latex figures/geom.eps "Azimuthal cross-section of the geometry. " width=0.35\textwidth

It is difficult to find non-trivial exact solutions of the governing equations
(1), (2), so we manufacture a
time-harmonic solution:

\f[
\mathbf{u} =
\begin{pmatrix}
u_r\\
u_\theta\\
u_z
\end{pmatrix} =
\cos(t)
\begin{pmatrix}
r^3\cos(z)\\
r^3z^3\\
r^3\sin(z)
\end{pmatrix}.\qquad\qquad (4)
\f]
This is an exact solution if we set the body force to
\f[
\mathbf{F} = 
\cos(t)
\begin{pmatrix}
-r\cos(z)\left\{\left(8+3r\right)\lambda +
\left(16 - r\left(r-3\right)\right)\mu + r^2\Lambda^2\right\}\\
-r\left\{8z^3\mu + r^2\left(z^3\Lambda^2 + 6\mu z\right)\right\}\\
r\sin(z)\left\{- 9\mu + 4r\left(\lambda + \mu\right) +
r^2\left(\lambda + 2\mu -
\Lambda^2\right)\right\}
\end{pmatrix},\qquad\qquad (5)
\f]
where \f$ \lambda = \nu/((1+\nu)(1-2\nu)) \f$ and \f$ \mu = 1/(2(1+\nu)) \f$ are the nondimensional Lam&eacute; parameters.
We impose the displacement along the boundaries \f$ r=r_\mathrm{max},
z=z_\mathrm{min}, z=z_\mathrm{max} \f$ according to (4),
and impose the traction
\f[
\hat{\bm{\tau}}_3 = \bm{\tau}(r_\mathrm{min},z,t) =
\cos(\omega t)
\begin{pmatrix}
-\cos(z)r_\mathrm{min}^2\left(6\mu +
\lambda\left(4+r_\mathrm{min}\right)\right)\\
-2\mu r_\mathrm{min}^2z^3\\
-\mu r_\mathrm{min}^2\sin(z)\left(3-r_\mathrm{min}\right)
\end{pmatrix},\qquad\qquad (6)
\f]
along the boundary \f$ r = r_\mathrm{min} \f$.

The problem we are solving then consists of equations
(1), (2) along with the body
force, (5) and boundary traction
(6). The initial conditions for the problem
are the exact displacement, velocity (and acceleration; see below) according to
the solution
(4).

\section results Results

The animation below shows the time dependent deformation of the cylinder in the
\f$ r \f$-\f$ z \f$ plane, while the colour contours indicate the azimuthal
displacement component.
The animation is for \f$ t \in [0,2\pi] \f$, since the time scale is
nondimensionalised on the reciprocal of the angular frequency of the
oscillations.

\image html figures/a.gif "Animation (HTML only) of the resulting displacement field. " 
\image latex figures/a.eps "Animation (HTML only) of the resulting displacement field. " width=0.75\textwidth
The next three figures show plots of the radial, axial and azimuthal
displacements as functions of \f$ (r,z) \f$ at \f$ t = 2.64 \f$. Note the
excellent agreement between the numerical and exact solutions.
\image html figures/u_r.gif "Comparison between exact and FE solutions for the r-component of displacement at t = 2.64 " 
\image latex figures/u_r.eps "Comparison between exact and FE solutions for the r-component of displacement at t = 2.64 " width=0.5\textwidth
\image html figures/u_z.gif "Comparison between exact and FE solutions for the z-component of displacement at t = 2.64 " 
\image latex figures/u_z.eps "Comparison between exact and FE solutions for the z-component of displacement at t = 2.64 " width=0.5\textwidth
\image html figures/u_theta.gif "Comparison between exact and FE solutions for the theta-component of displacement at t = 2.64 " 
\image latex figures/u_theta.eps "Comparison between exact and FE solutions for the theta-component of displacement at t = 2.64 " width=0.5\textwidth

<HR>
<HR>

\section global_namespace Global parameters and functions
As usual, we define all non-dimensional parameters in a namespace. In this
namespace, we also define the body force, the traction to be applied on the 
boundary \f$ r=r_\mathrm{min} \f$, and the exact solution.
Note that, consistent with the enumeration of the unknowns, discussed above, the
order of the components in the functions that specify the body force and the
surface traction is \f$ (r,z,\theta) \f$.
 
\dontinclude cylinder.cc
\skipline start_of_Global_Parameters_namespace
\until end of body force

In addition, the namespace includes the necessary machinery for providing the
time dependent equations with their initial data from the exact solution. There
are 9 functions, one for each of the components of displacement, velocity and
acceleration, and a helper function. For brevity, we list only one of these
functions; the others are similar.

\skipline Calculate the time dependent
\until end_of_u_r

<HR>
<HR>

\section main The driver code

We start by creating a \c DocInfo object that will be used to output the
solution, and then build the problem.

\dontinclude cylinder.cc
\skipline start_of_main
\until problem;

Next we set up a timestepper and assign the initial conditions.

\skipline Set the initial time
\until doc_info.number()

We calculate the number of timesteps to perform - if we are validating, just
do small number of timesteps; otherwise do a full period the time-harmonic
oscillation.

\skipline Find the number of timesteps to perform
\until end_of_calculate_number_of_timesteps

Finally we perform a time dependent simulation.

\skipline Do the timestepping
\until end_of_main

<HR>
<HR>

\section problem The problem class 

The problem class is very simple, and similarly to other problems with Neumann
conditions, there are separate meshes for the "bulk" elements and the "face"
elements that apply the traction boundary conditions. The function \c
assign_traction_elements() attaches the traction elements to the appropriate
bulk elements.

\dontinclude cylinder.cc
\skipline start_of_problem_class
\until end_of_problem_class

<HR>
<HR>

\section constructor The problem constructor

The problem constructor creates the mesh objects (which in turn create the
elements), pins the appropriate boundary nodes and assigns the boundary
conditions according to the functions defined in the \c Global_Parameters
namespace.

\skipline start_of_constructor
\until set_boundary_conditions()

Then the physical parameters are set for each element in the bulk mesh.

\until end_loop_over_elements

We then loop over the traction elements and set the applied traction.

\until end_loop_over_traction_elements

Finally, we add the two meshes as sub-meshes, build a global mesh from these and
assign the equation numbers.

\skipline Add
\until assign_eqn_numbers()
\skip end_of_constructor
\until end_of_constructor

<HR>
<HR>

\section traction_elements The traction elements
We create the face elements that apply the traction to
the boundary \f$ r=r_\mathrm{min} \f$.

\skipline start_of_traction
\until end of assign_traction_elements

<HR>
<HR>

\section initial_data Initial data
The time integration in this problem is performed using the Newmark scheme
which, in addition to the standard initial conditions
(3), requires an initial value for the acceleration.
Since we will be solving a test case in which the exact solution is known, we
can use the exact solution to provide the complete set of initial data required.
For the details of the Newmark scheme, see the tutorial on the
<A HREF="../../../linear_wave/two_d_linear_wave/html/index.html">linear wave equation</A>.

If we're doing an impulsive start, set the displacement, velocity and
acceleration to zero, and fill in the time history to be consistent with this.

\dontinclude cylinder.cc
\skipline start_of_set_initial_conditions
\until end_of_impulsive_start

If we are not doing an impulsive start, we must provide the timestepper with
time history values for the displacement, velocity and acceleration. Each
component of the these vectors is represented by a function pointer, and in this
case, the function pointers return values based on the exact solution.

\skipline Smooth start
\until d2_u_theta_dt2

Then we loop over all nodes in the bulk mesh and set the initial data values
from the exact solution.

\skipline Number of nodes
\until end_of_loop_over_nodes

<HR>
<HR>

\section doc Post-processing
This member function documents the computed solution to file and calculates the
error between the computed solution and the exact solution.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>

\section comments Comments
- Given that we non-dimensionalised all stresses on Young's modulus it seems odd
that we provide the option to specify a non-dimensional Young's modulus via
the member function \c AxisymmetricLinearElasticity::youngs_modulus_pt().
The explanation for this is that this function specifies the ratio of the
material's actual Young's modulus to the Young's modulus used in the
non-dimensionalisation of the equations. The capability to specify such ratios
is important in problems where the elastic body is made of multiple materials
with different constitutive properties. If the body is made of a single,
homogeneous material, the specification of the non-dimensional Young's modulus
is not required -- it defaults to 1.0.

\section exercises Exercises
- Try setting the boolean flag \c impulsive_start to \c true in the \c
AxisymmetricLinearElasticityProblem::set_initial_conditions function and
compare the system's evolution to that obtained when a "smooth" start from the
exact solution is performed.
- Omit the specification of the Young's modulus and verify that the default
value gives the same solution.
- Confirm that the assignment of the history values for the Newmark timestepper
  in \c AxisymmetricLinearElasticityProblem::set_initial_conditions sets the
  correct initial values for the displacement, velocity and acceleration. (Hint:
  the relevant code is already contained in the driver code, but was omitted in
  the code listings shown above.)

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/axisym_linear_elasticity/cylinder/">
demo_drivers/axisym_linear_elasticity/cylinder/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/axisym_linear_elasticity/cylinder/cylinder.cc">
demo_drivers/axisym_linear_elasticity/cylinder/cylinder.cc
</A>
</CENTER>
.


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


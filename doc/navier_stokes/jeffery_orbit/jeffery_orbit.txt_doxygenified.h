/**

\mainpage Motion of an elliptical particle in shear flow: Jeffery orbits

 We study the problem of the motion of a rigid elliptical particle
freely suspended in a shear flow as described by Jeffery (1922) 
[The motion of elliptical particles immersed in viscous fluid, <I>
 Proc. Roy. Soc. A </I> <B> 102 </B> 161-179]. 
 The problem is solved using \c oomph-lib's 
<a href="../../../../doc/meshes/mesh_from_inline_triangle/html/index.html">
inline unstructured mesh generation </a> procedures to modify 
the fluid mesh in response to changes in orientation of the ellipse.
 
\section overview Overview of the problem

\image html jeff_sketch.gif "A rigid ellipse immersed in a shear flow " 
\image latex jeff_sketch.eps "A rigid ellipse immersed in a shear flow " width=0.75\textwidth

We consider an ellipse with centre of mass fixed at the origin of a
Cartesian coordinate system \f$ (x_{1}^{*}, x_{2}^{*}) \f$, but
immersed in a viscous fluid undergoing a linear shear flow 
with shear rate \f$ G \f$. In the absence of the ellipse 
the fluid velocity field would be 
 \f$ u^{*}_{1} = Gx^{*}_{2} \f$, \f$ u^{*}_{2} = 0 \f$, where
 \f$ u^{*}_{i} \f$ is the dimensional velocity component in the
 \f$ x^{*}_{i} \f$ coordinate direction.

 The configuration of a rigid body in two dimensions is determined entirely
 by the position of its centre of mass, \f$(X^{*}_{1},X^{*}_{2})\f$ 
 and an angle, \f$\chi\f$ that specifies its orientation to a fixed
 axis. The equations governing the motion of the particle are then simply
 conservation of linear and angular momentum:
<center>
\f[ M \frac{\mbox{d}^{2} X^{*}_{i}}{\mbox{d} t^{*2}} = F^{*}_{i} \quad
 \mbox{and} \quad I \frac{\mbox{d}^{2} \chi}{\mbox{d} t^{*2}} = T^{*}, \f]
</center>
 where \f$M\f$ is the mass of the body, \f$I\f$ is its moment of
 inertia about the centre of mass, \f$(F^{*}_{1},F^{*}_{2})\f$ is the
 resultant force on the body and \f$T^{*}\f$ is the resultant torque about
 the centre of mass.

 In the present context, the force and torque on the body are entirely
 due to the viscous fluid loading on its surface in which case
<center>
\f[ F^{*}_{i} = \int \tau^{*}_{ij}\, n_{j}\, \mbox{d} s \quad\mbox{and}\quad
    T^{*} = \oint  (x^{*}_{1} - X^{*}_{1}) \tau^{*}_{2j}\,n_{j} - 
               (x^{*}_{2} - X^{*}_{2}) \tau^{*}_{1j}\,n_{j}\, \mbox{d}
 s, \f]
</center>
 where the integral is around the perimeter of the ellipse,
 \f$\tau^{*}_{ij}\f$ is the fluid stress tensor and
 \f$(n_{1},n_{2})\f$ is the unit normal to the ellipse surface, directed
 away from the solid body. 

 We non-dimensionalise the rigid-body equations, using the same 
problem-specific reference quantities as used in the
non-dimensionalisation of the Navier--Stokes equations, described 
<a href="../../../../doc/navier_stokes/driven_cavity/html/index.html">
in another tutorial. </a> Thus, \f$ {\cal U}\f$ is a typical fluid velocity, 
\f$ {\cal L} \f$ is the length scale, \f${\cal T}\f$ is the time scale
and the fluid pressure is non-dimensionalised on the viscous scale, \f$
\mu_{ref} {\cal U}/{\cal L} \f$, where \f$ \mu_{ref} \f$ is the
reference fluid viscosity. Hence,
<center>
\f[
\tau^{*}_{ij} = \frac{\mu_{ref} {\cal U}}{{\cal L}} \, \tau_{ij}, \qquad        
x_i^{*} = {\cal L}\, x_i, \qquad
X_i^* = {\cal L}\, X_i, \qquad        
t^* = {\cal T}\, t.
\f]
</center>
The external forces and torques are non-dimensionalised on the
viscous scales per unit length, \f$ F^{*}_{i} = \mu_{ref} {\cal U}
F_{i} \f$ and \f$ T^{*} = \mu_{ref} {\cal U} {\cal L} T \f$.
The dimensionless rigid-body equations are then
<center>
\f[ Re St^{2} \frac{\rho_{s}}{\rho_{f}}\, \bar{M}\, 
\frac{\mbox{d}^{2} X_{i}}{\mbox{d} t^{2}} = F_{i} 
\quad
 \mbox{and} \quad Re St^{2} \frac{\rho_{s}}{\rho_{f}}\, \bar{I}\, 
 \frac{\mbox{d}^{2} \chi}{\mbox{d} t^{2}} = T, \f]
</center>
 where the dimensionless parameters
<center>
\f[ Re = \frac{\rho_{f} {\cal U}{\cal L}}{\mu_{ref}}, \quad
    St = \frac{\cal L}{{\cal U}{\cal T}},\quad
    \frac{\rho_{s}}{\rho_{f}}, \quad
    \bar{M} = \frac{M}{\rho_{s} {\cal L}^{2}}, \quad
    \bar{I} = \frac{I}{\rho_{s} {\cal L}^{4}}, \f]
</center>
are the Reynolds number, the Strouhal number, the density ratio, and the
dimensionless mass and moment of inertia, respectively. In the
above \f$ \rho_{f} \f$ is the fluid density and \f$ \rho_{s} \f$ is
the solid density.

 In the specific problem considered here, the centre of mass is
 fixed, and the only possible motion of the particle is free
 rotation. The particle motion is therefore reduced to the solution of
 a single equation for the unknown angle. We choose \f$ {\cal L} = 2b
 \f$, the major axis of the ellipse, \f$ {\cal U} = 2Gb \f$ and \f$
 {\cal T} = G^{-1} \f$, so that \f$ St = 1 \f$ and the governing
 equations for the fluid and solid become
 <center>
 \f[  Re \left (\frac{\partial u_i}{\partial t} + 
u_j \frac{\partial u_i}{\partial x_j} \right) =
- \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left[
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right]
\quad\mbox{and}\quad
\frac{\partial u_i}{\partial x_i} = 0.
\f]
</center>
 and
<center>
\f[ Re \left(\frac{\rho_{s}}{\rho_{f}}\right) \bar{I} 
=    \oint  \left[x_{1}\tau_{2j}\,n_{j} - 
                x_{2} \tau_{1j}\,n_{j}\, \right] \mbox{d} s, \f].
</center>

We perform the computations in the domain \f$ x_{1} \in [-L,
 L]\f$ and \f$ x_{2} \in [-H, H] \f$, and apply the boundary conditions
<center>
 \f[ u_{1} = H f(t), u_{2} = 0, \quad\mbox{on}\quad x_{2} = -H; \f]
 \f[ u_{1} = -H f(t), u_{2} = 0, \quad\mbox{on}\quad x_{2} = H; \f]
 \f[ \tau_{11} = 0, u_{2} = 0, \quad\mbox{on}\quad x_{1} = -L, L; \f]
 where \f$f(t)\f$ is a smooth ramp function such that \f$f(0) = \dot{f}(0) =
 \ddot{f}(0) = 0\f$ and \f$f(t) \to 1\f$ as \f$ t \to \infty\f$.
</center>

\section results Results

 Jeffery (1922) showed that for a two-dimensional ellipse in Stokes
 flow (\f$ Re = 0 \f$), the exact solution for the angle as a function
 of time is
 \f[ \chi = \tan^{-1}\left(\frac{b}{a}\tan \frac{abGt}{a^{2} +
 b^{2}}\right),\quad\mbox{and}\quad \dot{\chi} = \frac{G}{a^{2} +
 b^{2}}\left(b^{2}\cos^{2}\chi + a^{2}\sin^{2}\chi\right). \f]
 Thus, the ellipse performs periodic orbits but with a non-uniform
 velocity. For sufficiently small \f$ Re \f$, Ding \& Aidun (2000)
 [The dynamics and scaling law for particles suspended in shear flow
 with inertia, <I> J. Fluid Mech. </I> <B> 423 </B> 317-344]
 showed that the system approximates the Jeffery orbits but with an
 increased period. Typical solutions for \f$ Re = 1\f$ are shown
 below.
 
\image html orbit_mesh.gif "An ellipse performing approximate Jeffery orbits at Re = 1 " 
\image latex orbit_mesh.eps "An ellipse performing approximate Jeffery orbits at Re = 1 " width=0.75\textwidth


\image html jeff_graph.gif "Angle and angular velocity as functions of time for Re = 1 " 
\image latex jeff_graph.eps "Angle and angular velocity as functions of time for Re = 1 " width=0.75\textwidth

\section fsi Implementation as a fluid-structure interaction problem

 The problem is a fluid-structure interaction problem in
which the structural dynamics are particularly simple, depending only
on a single degree of freedom. 
Nonetheless, the two types of physical coupling between the fluid 
and the solid remain:
-# The position of the free boundary depends on the position
   of the rigid body.
-# The rigid body is loaded by the fluid traction.
. 

 As in other 
<a href="../../../../doc/interaction/unstructured_fsi/html/index.html">
 2D unstructured FSI </a> problems, we treat the fluid mesh as a
 pseudo-solid body and determine the position of the boundary nodes 
 on the fluid-solid interface using  
 \c ImposeDisplacementByLagrangeMultiplierElements. 

 The rigid body mechanics is handled by using a \c GeomObject that
represents the perimeter of the rigid body to create
an \c ImmersedRigidBodyElement
that solves the three equations of motion for
a rigid body; and the load is applied to the rigid body
using \c NavierStokesSurfaceDragTorqueElements.

<HR>
<HR>

\section global Problem Parameters

 We use a namespace to define the parameters used in the problem

\dontinclude jeffery_orbit.cc
\skipline start_of_namespace
\until end_of_namespace

<HR>
<HR>

\section ellipse Defining the ellipse as a GeomObject

 We create a basic \c GeomObject to represent the ellipse whose
boundary we parametrise by the polar angle, measured from its centre
of mass.

\dontinclude jeffery_orbit.cc 
\skipline start_of_general_ellipse
\until end_of_general_ellipse

<HR>
<HR>

\section driver The driver code

 After parsing the command-line arguments, which are used to modify
 certain parameters for validation runs, a single instance of the
 \c UnstructuredImmersedEllipseProblem (described below) is
 constructed using Taylor Hood elements.
 
\dontinclude jeffery_orbit.cc
\skipline start_of_main
\until problem;

 After construction the \c Nodes on the boundary of the ellipse will
 have been directly mapped onto the curvilinear 
 surface using a strong (collocation)
 condition, \f$ x_{i} = R_{i}\f$, where \f$R_{i}\f$ is the
 corresponding boundary of the ellipse.
 In the full problem 
 the displacement boundary condition is enforced weakly 
 via Lagrange multipliers 
 \f$ \oint \left\{x_{i} - R_{i}\right\}\psi \mbox{d} s =
 0\f$. In order to ensure consistency, we initially solve the problem
 in which the rigid body is pinned so that the boundary nodes are
 adjusted to be consistent with the weak form of the boundary
 condition. We note that for sufficiently fine initial meshes 
 the difference is minimal.

 \until solve_for_consistent

 Now that we have a consistent initial condition, we initialise the
 timestepper and set conditions consistent with a impulsive start from
 rest.

 \until doc_solution();

 We then take a fixed number of timesteps on the initial mesh,
 documenting the solution after each solve.

 \until }

 Finally, we loop over a number of ``cycles'' in which we adapt the problem
 and then solve for a fixed number of time steps on the each mesh.

 \skipline for
 \until //end of main

<HR>
<HR>
\section problem The Problem class

 The Problem class follows the usual pattern. The time-dependent
 boundary conditions are applied using the \c
 actions_before_implicit_timestep() function and the no-slip boundary
 condition is applied in \c actions_before_newton_convergence_check()
 via an auxiliary node update function.

 <a href="../../../../doc/meshes/mesh_from_inline_triangle/html/index.html#adapt">
 Recall</a> that when adapting an unstructured mesh, its constituent
 elements are completely re-generated. Physical
 parameters and boundary conditions must therefore be reassigned, which is the
 task of the \c complete_problem_setup() function, called in \c
 actions_after_adapt().
 Helper functions are also provided to solve the initial problem to
 move the boundary nodes [\c solve_for_consistent_nodal_positions()]; 
 to apply the boundary conditions [\c set_boundary_velocity()];
 and to construct and delete the surface
 elements that impose the Lagrange multiplier constraints and compute the load
 on the rigid body [\c create_lagrange_multiplier_elements(), \c
 delete_lagrange_multiplier_elements(), \c create_drag_elements(), 
 \c delete_drag_elements()].

 The class also provides storage for the meshes, the rigid body 
and file handles for documentation.

\dontinclude jeffery_orbit.cc
\skipline start_of_problem_class
\until end_of_problem_class

<HR>
<HR>
\section constructor The Problem Constructor

 We begin by opening the output files and allocating two time
 steppers, one for the fluid problem and one for the rigid body
 problem.

\skipline start_constructor
\until Newmark

 We then define the geometry that defines the outer boundary of the
 unstructured mesh by constructing a \c TriangleMeshPolygon that 
 consists of four straight-line boundaries.

\until new TriangleMeshPolygon
 
 Next we build the single \c ImmersedRigidBodyElement from an
 instantiation of a \c GeneralEllipse
 geometric object.

\until this->time_stepper_pt(1)

 The \c ImmersedRigidBodyElement is used to define a 
\c TriangleMeshCurvilinearClosedCurve in exactly the same way as if it
were simply a (passive) \c GeomObject, as discussed in 
<a href="../../../../doc/meshes/mesh_from_inline_triangle/html/index.html">
another tutorial. </a>

\until linear_boundary_pt,hole_

 We then build the unstructured fluid mesh using the boundary
 information, set a spatial error estimator and complete the setup of
 the problem

\until complete_problem_setup()

 The \c ImmersedRigidBodyElement is not deleted during the adaptation
 process and so its physical parameters can be set once in the
 constructor. We set the initial position of the centre of mass,
 as well as the non-dimensional mass and moment of inertia shape, the 
 Reynolds and Strouhal numbers, and the density ratio, which appear in the
 governing equations above. For this problem, we also fix the location
 of the centre of mass. (The section \ref comments contains an
 exercise that asks you to explore what happens when you omit this step).

\until pin_centre_of_mass_coordinate(1)

 For later reference, we store the single \c ImmersedRigidBodyElement in a mesh
 
\until add_element_pt

 We then create the elements that apply the load on the rigid body and
 pass the entire mesh of elements to the rigid body. This is the
 equivalent of the function 
 \c FSI_functions::setup_fluid_load_info_for_solid_elements(..), but
 here the procedure is very simple because <b> all </b> the surface
 elements affect the single rigid body. 

\until rigid_element1_pt

 We next create the mesh of Lagrange-multiplier elements that drive
the deformation of the fluid mesh in response to the motion of the
ellipse

 \until create_lagrange
 
 and then construct the global mesh and assign equation numbers.

 \until end_of_constructor

 Note that the \c Drag_mesh_pt does not need to be added as a sub-mesh
 because its elements do not contribute <i> directly </i> to the
 residuals and Jacobian. 
 
<HR>
<HR>

\section complete Completing the problem setup

 The helper function \c complete_problem_setup() starts by
 (re-)applying the boundary conditions by pinning the fluid velocity
 in the \f$ x_{2}\f$-direction on all boundaries and that in the \f$
 x_{1} \f$-direction on the top and bottom (boundaries 1 and 3).

\dontinclude jeffery_orbit.cc
\skipline start_complete_problem_setup
\until ->pin(1)

 The boundary conditions for the solid degrees of freedom that
 describe the mesh deformation are assigned next by pinning the
 nodal positions on the fixed domain boundaries (boundaries 0, 1, 2,
 3):

 \until }

 The nodes on the boundary of the rigid body should be free to move,
 so they are unpinned and an auxiliary node update function is set
 to apply the no-slip boundary condition from the \c Node's positional
 history values. 

\until // End loop over boundaries

 We then loop over the fluid elements and set pointers to the physical
 parameters and the apply the velocity boundary conditions

\until end_of_complete_problem_setup

<HR>
<HR>

\section surface Creating and destroying the surface elements
 
 The general
 procedure for creating, attaching and deleting\c FaceElements is
 exactly the same as described 
<a href="../../../../doc/poisson/two_d_poisson_flux_bc/html/index.html">
 in another tutorial, </a> so is not described in detail here.

 The functions \c create_lagrange_multiplier_elements() and \c
 create_drag_elements() construct surface 
\c ImposeDisplacementByLagrangeMulitiplierElements and
\c NavierStokesSurfaceDragTorqueElements, respectively, around the
 boundary of the rigid body, setting any required member data. For
 example, the \c NavierStokesSurfaceDragTorqueElements require the
 location of the centre of mass in order to compute the torque. The
 elements are added to the internal storage containers 
\c Lagrange_multiplier_mesh_pt and \c Drag_mesh_pt. 
 The corresponding functions \c delete_lagrange_multiplier_elements()
 and \c delete_drag_elements() are used to delete and remove 
 the elements before adaptation.


<HR>
<HR>

\section bc Setting the boundary velocity

 The function \c set_boundary_velocity()  is used to apply the 
 time-dependent boundary conditions to the external boundaries of the
 fluid domain. The only subtlety is that after a remesh the history
 values for the boundary nodes must also be (re-)applied; and, for
 simplicity, the history values are always reset.

<HR>
<HR>

\section initial Solving for consistent initial nodal positions

 The initial nodal positions are made consistent with the
 weakly-imposed displacement boundary condition by pinning the rigid
 body degrees of freedom, performing a steady Newton solve and then
 releasing the rigid body degrees of freedom.

 \skipline start_solve_for_consistent_nodal_positions
 \until end_solve_for_consistent_nodal_positions

<HR>
<HR>

\section comments Comments and Exercises

\subsection exercises Exercises

-# Confirm that for sufficiently large times the solution agrees 
   with Jeffery's analytic solution when you set \f$ Re=0 \f$. 
   Explain why you expect there to be a discrepancy at early times.\n\n
-# What happens when the centre of mass is not fixed?
   Can you explain the observed behaviour?\n\n
-# What happens if you don't call the function 
   \c solve_for_consistent_nodal_positions()? 
   Can you explain the observed behaviour?\n\n
-# Investigate the behaviour of the system with increasing
   \f$Re\f$. What happens to the oscillations for \f$ Re > 30 \f$?\n\n
.

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/jeffery_orbit/">
demo_drivers/navier_stokes/jeffery_orbit/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/jeffery_orbit/jeffery_orbit.cc">
demo_drivers/navier_stokes/jeffery_orbit/jeffery_orbit.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


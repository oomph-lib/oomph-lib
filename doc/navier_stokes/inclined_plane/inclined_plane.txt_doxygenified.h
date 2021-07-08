/**

\mainpage Demo problem: Flow of a fluid film down an inclined plane

 The two-dimensional flow of a free surface down an inclined plane is
 a simple exact solution of the Navier--Stokes equations. We describe
 two different ways of solving the problem using either spines or a
 pseudo-elastic method to define the bulk mesh motion. Reassuringly, 
 the results are the same irrespective of the method chosen.

<HR>
<HR>

\section form Problem formulation

\image html sketch.gif "A film of incompressible viscous fluid of a given thickness flows down a plane inclined at a prescribed angle to the gravitational field. " 
\image latex sketch.eps "A film of incompressible viscous fluid of a given thickness flows down a plane inclined at a prescribed angle to the gravitational field. " width=0.75\textwidth

 Formulating the problem in coordinates tangential (\f$x^{*}\f$) and
 normal (\f$y^{*}\f$) to the plane and assuming that the flow is steady
 and only in the tangential direction, but independent of the tangential
 coordinate, reduces the momentum equations to 
  
<CENTER>
\f[ \frac{\partial p^*}{\partial x^*}  = \rho \,g \sin\alpha +
\mu\,\frac{\partial^{2} u^*}{\partial y^{*2}}, \f]
\f[ \frac{\partial p^*}{\partial y^*} = -\rho \,g \cos\alpha, \f]
</CENTER>
where \f$u^*\f$ is the velocity component tangential to the plane and
\f$p^*\f$ is the fluid pressure. Note that 
the continuity equation is automatically satisfied. 

We non-dimensionalise using the only length-scale in the
problem \f$H\f$, choosing the viscous scale for the pressure and
choosing a reference velocity scale \f$U\f$:
\f[ x^{*} = H x,\qquad y^{*} = H y,\qquad u^* = U u,\qquad p^* = \mu
U/H,\f] 
and the governing equations become
<CENTER>
\f[ \frac{\partial p}{\partial x}  = \frac{\rho \,g H^{2}}{\mu U} \sin\alpha +
    \frac{\partial^{2} u}{\partial y^{2}} = \frac{Re}{Fr} \sin \alpha + 
    \frac{\partial^{2} u}{\partial y^{2}}, \f]
\f[ \frac{\partial p}{\partial y} = -\frac{\rho \,g H^{2}}{\mu U} 
    \cos\alpha = - \frac{Re}{Fr} \cos\alpha. \f]
</CENTER>
The dimensionless grouping \f$ \rho g H^{2} / (\mu U) \f$ represents
the ratio of gravitational forces to viscous forces and we choose to
identify it as a Reynolds number \f$ Re = \rho U H / \mu\f$ divided by
a Froude number \f$ Fr = U^{2} / (g H)\f$. 

We proceed by assuming that the flow is driven entirely by the
gravitational body force and that there is no additional tangential 
pressure gradient. Then, integrating the tangential momentum balance
twice and using the boundary conditions of no-slip at the plane (\f$y=0\f$)and
that the free surface (\f$y=1\f$) is tangentially stress-free gives
\f[ u =  \frac{1}{2}\frac{Re}{Fr}\sin\alpha \left(2 y - y^{2}\right). \f]
Integrating the normal momentum balance and setting the reference
external pressure to be zero at the free surface gives
\f[ p = \frac{Re}{Fr}\cos\alpha (1 - y).\f]

Finally, we specify a "natural" velocity scale by setting \f$Re/Fr =
2\f$, corresponding to the velocity of the free-surface for a vertical
film (\f$\alpha = \pi/2\f$).

We shall assess the stability of the flat-film solution by applying a
small, short-duration perturbation to the wall velocity 
and evolving the system in time. If the interface is stable, 
the perturbation should decay, if not it should grow. A linear
stability  analysis for this problem was performed by Benjamin (1957)
and Yih (1963), who both found that for long waves in the absence of
surface tension, the interface was unstable when
\f[ Re > \frac{5}{4} \sin\alpha. \f]
(If you read the papers you will see that 
the Reynolds number was defined such that the
average fluid downslope velocity was one; to convert to our Reynolds
number, we must multiply by 3/2.)

The figure below shows the time evolution of the interface on a slope of
\f$ \pi/4 \f$ for
Reynolds numbers of zero (red line) and \f$ 4 \sin\alpha \f$ (green
line). The perturbation wavenumber is \f$ K = 0.1 \f$ and the interface rapidly develops waves that grow as they are
convected downslope for the higher Reynolds number, but decay when \f$
Re = 0 \f$.

\image html animated_surface.gif "Time evolution (or static snapshot at t = 7.5) of the interface shape for Re = 0 (Red) and 4 sin(alpha) (Green). " 
\image latex animated_surface.eps "Time evolution (or static snapshot at t = 7.5) of the interface shape for Re = 0 (Red) and 4 sin(alpha) (Green). " width=0.75\textwidth

The decay rate of the interfacial perturbation at \f$ Re = 0 \f$ is
slow, but can be seen in the next figure, which shows the height of
the interface at the downstream end of the domain plotted against
time. The domain is chosen so that it will contain three waves and the
decay or growth of successive crests and troughs can be seen.

\image html time_trace.gif "Time history of the interface position at the downstream end of the computational domain. " 
\image latex time_trace.eps "Time history of the interface position at the downstream end of the computational domain. " width=0.75\textwidth

\subsection bound A note on the boundary conditions

 Resolving the above analytic solution in a finite computational
domain requires some thought about boundary conditions. We are only
ever free to set one pressure value and setting the external pressure
to zero fixes the pressure within the fluid. The boundary conditions
at the plane are those of no-slip and at the free-surface 
the usual dynamic and kinematic conditions apply.
Nonetheless, we have a number of
possibilities for the boundary conditions at the "artificial" upstream
and downstream computational boundaries.
 - Prescribe periodic boundary conditions.
 - Prescribe the velocity profile as a Dirichlet condition at both
 ends.
 - Prescribe the appropriate hydrostatic pressure gradient and zero
 normal velocity.
.
 
We have chosen the last option, in which case the hydrostatic
 pressure gradient must be consistent with the external pressure. In
 other words, the pressure must be zero at the free surface
 (\f$y=0\f$). Changing the external pressure would correspond to
 changing the film thickness, so the external pressure is directly
 responsible for enforcing a specific volume constraint, unless 
 \f$\alpha = \pi/2\f$. When \f$\alpha = \pi/2\f$ there is no variation in 
  hydrostatic pressure through the film and its thickness is
  not specified by the external pressure.

We must also worry about the boundary conditions on the free surface
itself and we choose to impose a contact angle condition of \f$ \pi/2 \f$
at the upstream end, which ensures that the film remains flat. 
At the downstream end, we add a line tension term that arises
from use of the surface divergence theorem to integrate the
contribution of the dynamic boundary condition. This term can be used to
enforce contact angle conditions in a weak formulation, but here we
simply add the term using the angle calculated from the current
position of the free surface. 

<HR>
<HR>

\section global Global parameters and functions
 The global parameters are the Reynolds number, 
 the dimensionless grouping \f$Re/Fr\f$, the angle of inclination of
 the slope \f$\alpha\f$, the direction of the gravity vector \f$G\f$
 and the capillary number \f$Ca\f$, which only influences the dynamics.

\dontinclude inclined_plane.cc
\skipline Global_Physical_Variables
\until Ca=

 The hydrostatic pressure field is specified as an applied
 traction. At the outlet (inlet), the outer unit normal is in the 
 positive (negative)
 \f$x\f$ direction and so the required traction is given by \f$-p\f$
 (\f$p\f$),
 \f[ \mbox{\boldmath$t$}|_{\mbox{outlet}} = (- (Re/Fr)\cos\alpha (1 -
 y), 0), \qquad 
 \mbox{\boldmath$t$}|_{\mbox{inlet}} = ((Re/Fr)\cos\alpha (1 -
 y), 0).\f] 
 These tractions are specified by the two different functions
 \skipline hydrostatic
 \until end of traction
 Note that \c G [ \c 1 ] is the component of the gravitational body force in
 the vertical direction, so \f$ G[1] = - \cos\alpha\f$.  

 We must also specify the direction of the normals (directed out of the
 fluid) to the notional walls that form the inlet and outlet 
 and a contact angle of \f$ \pi /2 \f$ that will be
 used as a boundary condition on the free surface at the upstream end
 of the domain. In this case the normal to the inlet is in the
 negative x-direction and the normal to the outlet is in the positive
 x-direction. The actual value of the \c Wall_normal vector is set in
 \c main()

\dontinclude inclined_plane.cc
\skipline  wall normal
\until Inlet_Angle

<HR>
<HR>
\section main The driver code
 We start by specifying the constitutive law used to define the mesh
 motion when pseudo-elastic deformation is used.
\skipline start of main
\until ::Nu
 Next, the type of fluid element is 
 chosen according to specified compiler flags
\skipline CR_ELEMENT
\until #endif
 We then initialise the physical parameters, the Reynolds number and
 the direction of the gravitational body force, both based on the angle
 of inclination \f$\alpha\f$.
\skipline Initialise physical
\until G[1] 
 We also set the direction of the notional wall normal vector.
\until Wall_normal[1]

 We now create the spine version of the problem, solve the steady
 problem, assign initial conditions by assuming that the problem has
 been at the steady state for all previous times, 
 and then evolve the system in time.
 \skipline Spine problem
 \until End of spine problem
 Finally, exactly the same procedure is performed for the elastic problem
 \skipline Elastic problem
 \until End of elastic problem

<HR>
<HR>
\section mesh The mesh classes

 The base mesh class is the \c SimpleRectangularQuadMesh: boundary 0
 will be the wall; boundary  2 will be the free surface; and the
 remaining boundaries will be the inlet (3) and outlet (1). Below we shall
 demonstrate how to convert an existing mesh into a \c SpineMesh and
 \c ElasticMesh suitable for free-surface problems.

\subsection spine_mesh Creating the spine mesh
The \c SpineInclinedPlaneMesh  inherits from the generic \c
SimpleRectangularQuadMesh and adds vertical spines to the Nodes
within the mesh in the constructor. 
Note that the resulting mesh is essentially the same
as the \c SingleLayerSpineMesh, but has a somewhat simpler interface.

\dontinclude inclined_plane.cc
\skipline Create a spine mesh
\until end of constructor
In addition, a \c spine_node_update() function must be provided that
determines how the \c Nodes move as functions of the \c Spines.
\skipline General node
\until }

\subsection elastic_mesh Creating the ElasticMesh

The \c ElasticInclinedPlaneMesh inherits from the \c
SimpleRectangularQuadMesh and the undeformed (reference)
configuration is set to be the current position of the \c Nodes.
\skipline Create an Elastic mesh
\until };
Note that the specification of the ElasticMesh is much simpler than
that of a SpineMesh because no decision needs to be taken about how to
describe the motion using Spines.

<HR>
<HR>

\section problem The problem classes
\subsection generic_prob The generic problem
 For ease of exposition, all generic functionality is included
 in the \c InclinedPlaneProblem class, which is templated by the bulk
 \c ELEMENT and the \c INTERFACE_ELEMENT. The class includes storage
 for the different sub-meshes: Bulk, the Traction elements associated
 with the inlet and outlet, the (free) Surface elements and the point
 elements associated with the ends of the interface. In addition, a
 string \c Output_prefix is used to distinguish
 between the output files from different formulations.
 \dontinclude inclined_plane.cc
 \skipline Generic problem class
 \until Output_prefix

 The time-dependent perturbation is introduced in the function \c
 actions_before_implicit_timestep(), which sets the vertical
 velocity on the wall (boundary 0)
 \f[ v = \epsilon \sin(K x) t \mbox{e}^{-t} \f]
 \skipline actions_before_implicit_timestep()
\until end_of_actions_before_implicit
 
 The function \c make_traction_elements() creates
 \c NavierStokesTractionElement s adjacent to the mesh boundaries 3
 (the inlet) and 1 (the inlet). These elements are added to the \c
 Mesh \c Traction_mesh_pt, which is itself constructed 
 in the function and pointers to the 
 appropriate traction functions are assigned.
 \skipline traction boundary elements
 \until end of make_traction_elements

 The function \c make_free_surface_elements() creates the appropriate 
\c INTERFACE_ELEMENTs adjacent to the free surface (boundary 2), sets 
the capillary number and also creates free-surface boundary elements
at the left- and right-hand ends of the interface. If these "point"
elements are not included then the surface tension is not
applied correctly at the edges of the domain. The contact angle is set
to be the value \c Inlet_Angle at the left-hand edge of the domain.
\skipline Make the free
\until end of make_free_surface_elements


 The function \c
 complete_build() assigns physical parameters to the fluid elements,
 sets the boundary conditions and assigns equation numbers.
 \dontinclude inclined_plane.cc
 \skipline complete_build(
 \until end of complete_build
 Note that boundary conditions for the nodal positions in the
 pseudo-elastic formulation are specified by testing whether the \c
 Nodes are \c SolidNodes. In this case, the \c Nodes on the
 inlet and outlet boundaries are constrained to remain at the same
 horizontal position and the \c Nodes on the plane wall are fixed. 

 The function \c solve_steady() initialises the velocity of at
 all \c Nodes to the flat-film solution, solves the steady equations
 and writes the solution to a file.
 \skipline :solve_steady()
 \until end of solve_steady

 Finally, the function \c timestep() takes a number of fixed timesteps
 writing vertical positions and the time to a trace file and writing
 the complete flow field to disk after a given number of timesteps.
\skipline timestep(
\until end of timestep

\subsection spine_prob The spine-based formulation

The class \c SpineInclinedPlaneProblem  inherits
from the generic \c InclinedPlaneProblem  class
and requires only minor modification. The constructor sets the string
\c Output_prefix, builds a timestepper, builds the specific \c
SpineMesh, creates the appropriate \c FaceElements, adds all
sub-meshes to the \c Problem, builds the global mesh 
 and then calls \c InclinedPlaneProblem::complete_build(). 
\skipline SpineInclinedPlaneProblem(
\until }

In a spine-based formulation, the nodal positions must be updated
after every Newton step, which is achieved by overloading the function
\c Problem::actions_before_newton_convergence_check()
\skipline Spine heights/
\until }

We also specify a destructor to clean up memory allocated by the
class.

\subsection elastic_probe The pseudo-solid-based formulation
The class \c ElasticInclinedPlaneProblem inherits from the generic \c
InclinedPlaneProblem class and also requires only minor
modification. The constructor sets the string \c Output_prefix, builds
a timestepper, builds the specific \c SolidMesh, sets the constitutive
law for the bulk elements,  creates the appropriate \c FaceElements, adds all
sub-meshes to the \c Problem, builds the global mesh and
then calls \c
InclinedPlaneProblem::complete_build()
\skipline ElasticInclinedPlaneProblem(
\until end of constructor
 
In a pseudo-solid formulation, it is advantageous to reset the
undeformed configuration after every timestep (an updated Lagrangian
formulation). Hence, the \c Problem::actions_after_implicit_timestep()
function is overloaded
\skipline actions_after_implicit
\until  end of actions_

We also specify a destructor to clean up memory allocated by the
class.

<HR>
<HR>

\section exercises Exercises
-# Confirm that the steady solution agrees with the exact solution.
-# Investigate what happens when the angle is varied. What happens
when the angle is set to zero? What happens when the angle is set to
\f$\pi/2\f$?
-# What happens if the hydrostatic pressure boundary conditions are
not applied?
-# How does the stability of the system to the perturbation change
with angle, \f$ Ca \f$ and \f$ K \f$? Are the results in agreement
with the theoretical predictions?
-# Are the results independent of the length of the domain? 
-# Compare the spine-based and pseudo-elastic-based formulations? What
is the same and what is different? Which method do you prefer?
<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/inclined_plane/">
demo_drivers/navier_stokes/inclined_plane/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/inclined_plane/inclined_plane.cc">
demo_drivers/navier_stokes/inclined_plane/inclined_plane.cc
</A>
</CENTER>

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


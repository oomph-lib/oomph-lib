/**

\mainpage Example problem: Adaptive simulation of finite-Reynolds number entry flow into a 3D tube

In this example we shall demonstrate the spatially adaptive solution of the 
steady 3D Navier-Stokes equations using the problem of developing pipe flow.

<HR>
<HR>

\section example The example problem

<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>The 3D developing pipe flow in a quarter-tube domain.</B>
</CENTER> 
Solve the steady Navier-Stokes equations:
\f[
Re\phantom{i}u_j\frac{\partial u_i}{\partial x_j} =
- \frac{\partial p}{\partial x_i} +
\frac{\partial }{\partial x_j} \left(
\frac{\partial u_i}{\partial x_j} +  
\frac{\partial u_j}{\partial x_i} \right),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
and
\f[
\frac{\partial u_i}{\partial x_i} = 0,
\f]
in the quarter-tube domain \f$
D = \left\{x_1 \geq 0, x_2 \geq 0, x_1^2+x_2^2 \leq 1, x_3 \in [0,L]
\right\}\f$,
subject to the Dirichlet boundary conditions:
\f[
\left. \mathbf{u}\right|_{\partial D_{wall}}=(0,0,0),
\ \ \ \ \ \ \ \ \ \ (2)
\f]
on the curved wall 
\f$ \partial D_{wall} = \{(x_1,x_2,x_3) \ | \ x_1^2+x_2^2=1\}, \f$
\f[
\left. \mathbf{u}\cdot\mathbf{n}\right|_{\partial D_{sym[1,2]}}=0,
\ \ \ \ \ \ \ \ \ \ (3)
\f]
(where \f$ \mathbf{n} \f$ is the outer unit normal vector) on the
symmetry boundaries 
\f$ \partial D_{sym[1]} = \{(x_1,x_2,x_3) \ | \ x_1=0\} \f$
and
\f$ \partial D_{sym[2]} = \{(x_1,x_2,x_3) \ | \ x_2=0\}, \f$
\f[
\left. \mathbf{u}\right|_{\partial D}=(0,0,1-(x_1^2+x^2)^\alpha),
\ \ \ \ \ \ \ \ \ \ (4)
\f]
on the inflow boundary,
\f$ \partial D_{inflow} = \{(x_1,x_2,x_3) \ | \ x_3=0\}, \f$ and finally
\f[
\left. u_1\right|_{\partial D_{outflow}}= 
\left. u_2\right|_{\partial D_{outflow}}=0,
\ \ \ \ \ \ \ \ \ \ (5)
\f]
(parallel flow) on the outflow boundary \f$ \partial D_{outflow} = 
\{(x_1,x_2,x_3) \ | \ x_3=L\}.\f$
Note that the axial velocity component,
\f$u_3\f$, is not constrained at the outflow. Implicitly, we are therefore 
setting the axial component of traction on the fluid to zero,
\f[ 
t_3|_{\partial D_{outflow}} = \left. \left(-p +  2 \, 
\frac{\partial u_3}{\partial x_3}\right)
\right|_{\partial D_{outflow}} = 0.
\f]  
Since \f$ \partial u_1 / \partial x_1 = \partial u_2
/\partial x_2 = 0\f$ in the outflow cross-section [see 
(5], and the flow is incompressible, 
this is equivalent to (weakly) setting the pressure at the outflow to zero,
\f[
p|_{\partial D_{outflow}}=0.
\f]
</TD>
</TR>
</TABLE>  
</CENTER>

<HR>
<HR>

\subsection results Results for Taylor-Hood elements

The figure below shows the results computed with \c oomph-lib's 3x3x3-node 
3D adaptive Taylor-Hood elements for the
parameters \f$ \alpha=20\f$, \f$ L=7\f$ and \f$Re=100 \f$. 
The large exponent \f$ \alpha=20\f$ imposes a very blunt inflow
profile, which creates a thin boundary layer near the wall.
Diffusion of vorticity into the centre of the tube smooths the 
velocity profile which ultimately approaches a parabolic Poiseuille
profile. If you are viewing these results online you will be able
to see how successive mesh adaptations refine the mesh -- predominantly
near the entry region where the large velocity gradients in the
boundary layer require a fine spatial discretisation. Note also that
on the coarsest mesh, even the (imposed) inflow profile is 
represented very poorly.

\image html axial_veloc.gif "Contour plot of the axial velocity distribution for Re=100. Flow is from right to left. " 
\image latex axial_veloc.eps "Contour plot of the axial velocity distribution for Re=100. Flow is from right to left. " width=0.75\textwidth
\image html full_profiles.gif "Axial velocity profiles in equally-spaced cross-sections along the tube for Re=100. Flow is from left to right. " 
\image latex full_profiles.eps "Axial velocity profiles in equally-spaced cross-sections along the tube for Re=100. Flow is from left to right. " width=0.75\textwidth

<HR>
<HR>

\section namespace Global parameters and functions
The problem only contains one global parameter, the Reynolds number,
which we define in a namespace, as usual.

\dontinclude three_d_entry_flow.cc
\skipline start_of_namespace
\until end_of_namespace

<HR>
<HR>

\section main The driver code
Since the 3D computations can take a long time, and since all demo codes
are executed during \c oomph-lib's self-test procedures,
we allow the code to operate in two modes:
- By default, we specify error targets for which the code
  refines the mesh near the inflow region,  and allow up to 
  five successive mesh adaptations. The code is executed
  in this mode if the executable is run without any command
  line arguments.
- If the code is run during the self-test procedure (indicated by specifying 
  some random command line argument), we only perform one level of
  adaptation to speed up the self-test. However, because the original mesh
  is very coarse, the first mesh adaptation refines \e all elements
  in the mesh (cf. the animation of the adaptive mesh refinement
  shown above), so that no hanging nodes are created -- not a good
  test-case for a validation run! Therefore adjust the error
  targets so that the first (and only) mesh adaption only refines 
  a few elements and therefore creates a few hanging nodes.
.
The main code therefore starts by storing the command line arguments
and setting the adaptation targets accordingly:

\skipline start_of_main
\until end max_adapt setup

We then create a \c DocInfo object to specify the labels for the
output files, and solve the problem, first with Taylor-Hood and then
with Crouzeix-Raviart elements, writing the results from the
two discretisations to different directories:

\until end_of_main

<HR>
<HR>

\section problem The problem class
The problem class is very similar to the ones used in the 
<A HREF="../../circular_driven_cavity/html/index.html#problem">
2D examples</A>. We pass the \c DocInfo object and the target errors
to the \c Problem constructor.

\dontinclude three_d_entry_flow.cc
\skipline start_of_problem_class
\until ~

The function \c Problem::actions_after_newton_solve() is used
to document the solutions computed at various levels of
mesh refinement:

\until }

The function \c Problem::actions_before_newton_solve() is discussed 
<A HREF="#before_solve">below</A>, and, as in all adaptive
Navier-Stokes computations, we use the function 
\c Problem::actions_after_adapt() to pin any redundant pressure
degrees of freedom; see
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details.

\until }

Finally, we have the usual \c doc_solution() function and 
include an access function to the mesh. The private member data
\c Alpha determines the bluntness of the inflow profile.

\until end_of_problem_class

<HR>
<HR>

\section constructor The constructor

We start by building the adaptive mesh for the quarter tube domain.
As for most meshes with curvilinear boundaries, 
the \c RefineablequarterTubeMesh expects the curved boundary to be 
represented by a \c GeomObject. We therefore create an \c EllipticalTube with
unit half axes, i.e. a unit cylinder and a pass a pointer to the \c
GeomObject to the mesh constructor. The "ends" of 
the curvilinear boundary (in terms of the maximum and minimum values
of the Lagrangian coordinates that parametrise the shape of the 
\c GeomObject) are such that it represents a quarter of a cylindrical tube 
of length \f$L=7\f$.

\skipline start_of_constructor
\until RefineableQuarterTubeMesh

Next, we build an error estimator and specify the target errors
for the mesh adaptation:

\until min_permitted_error

Now we have to apply boundary conditions on the various mesh
boundaries. [\b Reminder: If the numbering of the
mesh boundaries is not apparent from its documentation (as it should
be!), you can use the function \c Mesh::output_boundaries(...)
to output them in a tecplot-readable form. 

\until some_file.close();

If the mesh has
N boundaries, the output file will contain N different zones, each
containing the \f$(x,y,z)\f$ coordinates of the nodes on the
boundary.]


For the \c RefineableQuarterTubeMesh, the boundaries are numbered
as follows:
- Boundary 0: "Inflow" cross section; located along the line
  parametrised by \f$ x_3=\xi_0 = \xi_0^{lo} \f$ on the 
  \c GeomObject that specifies the wall.
- Boundary 1: Plane \f$ x_1=0 \f$
- Boundary 2: Plane \f$ x_2=0 \f$
- Boundary 3: The curved wall
- Boundary 4: "Outflow" cross section; located along the line
  parametrised by \f$ x_3=\xi_0 = \xi_0^{hi} \f$ on the 
  \c GeomObject that specifies the wall.
.

\image html boundaries.gif "Plot of the mesh boundaries. " 
\image latex boundaries.eps "Plot of the mesh boundaries. " width=0.75\textwidth


We apply the following boundary conditions:
- Boundary 0: (\f$ x_3=0\f$) pin all three velocities.
- Boundary 1: (\f$ x_1=0\f$) pin \f$ u_1\f$.
- Boundary 2: (\f$ x_2=0\f$) pin \f$ u_2\f$.
- Boundary 3: (\f$ x_1^2+x_2^2=1\f$) pin all three velocities.
- Boundary 4: (\f$ x_3=L\f$) pin \f$ u_1 \f$ and \f$ u_2 \f$.
.

\until end loop over boundaries

Now we assign the \c re_pt() for each element and pin the redundant
nodal pressures (see
<A HREF="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
another tutorial</A> for details).

\until pin_redundant_nodal_pressures

We provide an initial guess for the velocity field by
initialising all velocity components with their Poiseuille flow values.

\until }

Finally, we set the value of \c Alpha, the exponent that specifies
the bluntness of the inflow profile and assign
the equation numbers.

\until end_of_constructor

<HR>
<HR>

\section before_solve Actions before solve
We use the function \c Problem::actions_before_newton_solve() to re-assign
the inflow boundary conditions before every solve. In the present
problem this is an \e essential step because the blunt inflow profile
(4)
cannot be represented accurately on the initial coarse mesh (see the
animation of the axial velocity profiles at the beginning of this
document). As discussed in the 
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">
example that illustrates the use of spatial
adaptivity for time-dependent problems</A>,
\c oomph-lib automatically <B>(i)</B> applies the correct boundary conditions
for newly created nodes that are located on the mesh boundaries
and <B>(ii)</B> assigns the nodal values at such nodes by interpolation from
the previously computed solution. This procedure is adequate on
boundaries where homogeneous boundary conditions are applied,
e.g. on the curved wall, the symmetry and the outflow boundaries.
However, on the inflow boundary, the interpolation from the FE
representation of the blunt velocity profile (imposed on 
the coarse initial mesh) onto the refined mesh, does not yield a 
more accurate representation of the prescribed inflow profile. It is
therefore necessary to re-assign the nodal values on this boundary 
after every adaptation, i.e. before every solve.

\skipline start_of_actions_before_newton_solve
\until }


<HR>
<HR>

\section doc Post processing
This function remains exactly the same as in the 
<A HREF="../../driven_cavity/html/index.html"> 2D examples</A>.

\skipline start_of_doc_solution
\until end_of_doc_solution

<HR>
<HR>

\section comments Comments and exercises
-# Suppress the reassignment of the prescribed inflow profile in
   \c Problem::actions_before_newton_solve() to confirm that this step
   is essential if the computation is to converge to the exact
   solution.
-# Suppress the specification of the parabolic (Poiseuille) velocity
   profile as the initial guess for the velocity field in the problem 
   constructor to confirm that the assignment of a "good" initial 
   guess for the solution is essential for the convergence of 
   the Newton method. [Hint: You can simply comment out the initialisation
   of the velocities -- they then retain their default initial values of
   0.0. When you re-run the code, the Newton iteration will "die" 
   immediately with an error message stating that the maximum residual
   exceeds the default threshold of 10.0, stored in the protected data member
   \c Problem::Max_residuals. Try increasing the value of this
   threshold in the Problem constructor. Is this sufficient to make the
   Newton method converge?]
   
   


<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/three_d_entry_flow/">
demo_drivers/navier_stokes/three_d_entry_flow/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/three_d_entry_flow/three_d_entry_flow.cc">
demo_drivers/navier_stokes/three_d_entry_flow/three_d_entry_flow.cc
</A>
</CENTER>
.
















<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


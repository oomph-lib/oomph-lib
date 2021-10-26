/**

\mainpage Example problem: Adaptive solution of the 2D Poisson equation with flux boundary conditions

In this document we discuss the adaptive solution of a 2D Poisson
problem with Neumann boundary conditions.
<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B>Two-dimensional model Poisson problem with Neumann boundary conditions</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2 \frac{\partial^2u}{\partial x_i^2} = f(x_1,x_2),
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the rectangular domain \f$D = \left\{ (x_1,x_2) \in 
[0,1] \times [0,2]\right\} \f$. The domain boundary
\f$ \partial D = \partial D_{Neumann} \cup \partial D_{Dirichlet} \f$,
where \f$ \partial D_{Neumann} 
= \left\{ (x_1,x_2) | x_1=1, \ x_2\in [0,2] \right\} \f$.
On \f$ \partial D_{Dirichlet}\f$ we apply the Dirichlet boundary conditions
\f[
\left. u\right|_{\partial D_{Dirichlet}}=u_0,
\ \ \ \ \ \ \ \ \ \ (2)
\f]
where the function \f$ u_0 \f$ is given. On  \f$ \partial D_{Neumann}\f$
we apply the Neumann conditions 
\f[
\left. \frac{\partial u}{\partial n}\right|_{\partial D_{Neumann}} = 
\left. \frac{\partial u}{\partial x_1}\right|_{\partial D_{Neumann}}
=g_0,
\ \ \ \ \ \ \ \ \ \ (3)
\f]
where the function \f$ g_0 \f$ is given.
</TD>
</TR>
</TABLE>  
</CENTER>  

In two previous examples we demonstrated two approaches for 
the \e non-adaptive solution of this problem.
In both cases we created a "bulk" mesh of \c QPoissonElements 
and applied the flux boundary conditions by "attaching" 
\c PoissonFluxElements to the appropriate boundaries 
of the "bulk" Poisson elements. In the 
<A HREF="../../two_d_poisson_flux_bc/html/index.html"> first
implementation</A> we simply added the pointers to the flux
elements to the bulk mesh; in the 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html"> second
implementation</A> we stored the surface and bulk elements in 
separate meshes and combined them to a single global mesh. We 
will now demonstrate that the second approach greatly facilitates the 
automatic problem adaptation. We use the mesh adaptation procedures 
of the \c RefineableQuadMesh class to adapt the bulk mesh, 
driven by the spatial error estimates in that mesh. 
The flux elements are neither involved in nor adapted by 
the refinement process of the bulk mesh. We therefore use
the function \c
Problem::actions_before_adapt() to delete the flux elements
before the adaptation, and \c Problem::actions_after_adapt() 
to re-attach them once the bulk mesh has been adapted.

As in the 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html">
previous example</A> we choose a source function and
boundary conditions for which the function 
\f[
u_0(x_1,x_2) = \tanh(1-\alpha(x_1 \tan\Phi - x_2)),
\ \ \ \ \ \ \ \ \  (4)
\f]
is the exact solution of the problem.

\image html rotate.gif "Plot of the solution obtained with automatic mesh adaptation " 
\image latex rotate.eps "Plot of the solution obtained with automatic mesh adaptation " width=0.75\textwidth

Since many functions in the driver code are identical to that in
the non-adaptive version, discussed in the
<A HREF="../../two_d_poisson_flux_bc2/html/index.html">
previous example</A>, we only list those functions that
differ. Please consult the source code 
<A HREF="../../../../demo_drivers/poisson/two_d_poisson_flux_bc_adapt/two_d_poisson_flux_bc_adapt.cc">two_d_poisson_flux_bc_adapt.cc
</A> for full details of the implementation.

<HR>   
<HR>


\section global Global parameters and functions
The specification of the source function and the exact solution
in the namespace \c TanhSolnForPoisson is identical to that
in the single-mesh version discussed in the
<A HREF="../../two_d_poisson_flux_bc/html/index.html">
previous example</A>. 
 
<HR>
<HR>

\section main The driver code
The main code is virtually identical to that in the 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html"> previous 
non-adaptive example </A>. The only change is the provision of
an argument to the Newton solver

\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline Solve the problem
\until (3);

which indicates that the problem should be adapted up to three times.

<HR>
<HR>

\section problem The problem class
The problem class is very similar to that in the 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html#problem">
non-adaptive implementation</A>: The only difference is the
provision of the functions \c actions_before_adapt(), \c
actions_after_adapt(),  \c set_prescribed_flux_pt(), 
and \c delete_flux_elements(...) which we discuss in more detail
below.

 
\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline start_of_problem_class
\until end of problem class

[See the discussion of the 
<A HREF="../../one_d_poisson/html/index.html">
1D Poisson problem</A> for a more detailed discussion of the
function type PoissonEquations<2>::PoissonSourceFctPt.]
 
<HR>
<HR>

\section constructor The Problem constructor
We create the bulk mesh and surface mesh as 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html#constructor">
before</A>. Next we create the spatial error estimator and
pass it to the bulk mesh.

\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline start_of_constructor
\until new Z2ErrorEstimator

Apart from this, the problem is constructed as in the 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html#constructor">
non-adaptive previous example.</A>

\until Surface_mesh
\until end of constructor

<HR>
<HR>



\section before_adapt Actions before adaptation
The mesh adaptation is driven by the error estimates in the bulk
elements and only performed for that mesh. The flux elements must 
therefore be removed before adaptation. We do this by calling 
the function \c delete_flux_elements(...), and
then rebuilding the \c Problem's global mesh.

\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline start_of_actions_before_adapt
\until end of actions_before_adapt

<HR>
<HR>

\section after_adapt Actions after adapt
After the (bulk-)mesh has been adapted, the flux elements must be re-attached.
This is done by calling the function \c
create_flux_elements(...), followed by a rebuild of the \c Problem's global
mesh. Finally, we set the function pointer to the prescribed flux
function for each flux element.
 
\skipline start_of_actions_after_adapt
\until end of actions_after_adapt

<HR>
<HR>


\section delete_flux Delete flux elements
This function  loops over all the flux elements (i.e. those in the
surface mesh) and deletes them and their storage.

\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline start_of_delete_flux_elements
\until end of delete_flux_elements

\b IMPORTANT: Note how the elements are first deleted "manually" and then "flushed"
from the surface mesh, using the function 
\c Mesh::flush_element_and_node_storage(). This is necessary
because deleting the surface mesh directly (by \c delete \c
surface_mesh_pt;) would  also delete its constituent \c Nodes
which are shared with the bulk mesh and must therefore be 
retained!




<HR>
<HR>



\section actions_before Actions before solve
This remains as 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html#actions_before">
before. </A> 

<HR>
<HR>

\section doc Post-processing
This remains as 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html#doc">
before. </A> 

<HR>
<HR>

\section create_flux Create flux elements
This remains as 
<A HREF="../../two_d_poisson_flux_bc2/html/index.html">
before. </A> 

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/two_d_poisson_flux_bc_adapt/">
demo_drivers/poisson/two_d_poisson_flux_bc_adapt/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/poisson/two_d_poisson_flux_bc_adapt/two_d_poisson_flux_bc_adapt.cc">
demo_drivers/poisson/two_d_poisson_flux_bc_adapt/two_d_poisson_flux_bc_adapt.cc
</A>
</CENTER>
.


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


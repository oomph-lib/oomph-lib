/**

\mainpage Parallel solution of Turek & Hron's FSI benchmark problem with spatial adaptivity for the fluid and solid meshes




This document provides an overview of 
- how to change the serial driver code for 
  <a href="../../../interaction/turek_flag/html/index.html">Turek &
  Hron's FSI benchmark problem</a> so that both the fluid and solid meshes
  can be adapted, 
. 
- how to distribute the problem across multiple processors,
.
and
- how to enable load balancing of the problem once it is distributed.
.
The document is part of a
<a href="../../../example_code_list/html/index.html#distributed">
series of tutorials</a>
that discuss how to modify existing serial driver codes so that the
\c Problem object can be distributed across multiple processors.



<HR>
<HR>

\section double_adapt Enabling spatial adaptivity for the fluid and solid meshes
In the original (serial) driver code for 
<a href="../../../interaction/turek_flag/html/index.html">Turek &
Hron's FSI benchmark problem</a> we only adapted the fluid mesh. 
Before discussing how to modify the code to refine the 
fluid and solid meshes simultaneously, we provide a
brief reminder of the procedure used to discretise fluid-structure
interaction problems that involve fluid and solid domains of equal spatial
dimension (e.g. a 2D fluid domain interacting with a 2D solid domain)
when using algebraic node update methods to adjust the position of
the nodes in the fluid mesh.  
We refer to <a href="../../fsi_channel_with_leaflet/html/index.html">
another tutorial</a> for a discussion of FSI
problems involving the interaction of fluids with (lower-dimensional)
shell and beam structures.

\subsection orig General methodology
  
The figure below shows a sketch of  a simple(r)
fluid-structure interaction problem involving fluid and solid domains 
that meet along a single mesh boundary. We assume that the fluid mesh 
uses an algebraic node 
update function to adjust the position of its nodes 
in response to changes in the domain 
boundary, represented by the \c GeomObject shown in 
magenta. (You may wish to consult
<a href="../../../interaction/fsi_collapsible_channel_algebraic/html/index.html">
another tutorial</a> for a reminder of how \c oomph-lib's algebraic
node update methods work). 


\image html fsi_with_adaptive_wall_basic.gif "Basic setup for FSI problems involving algebraic node updates for the fluid mesh. " 
\image latex fsi_with_adaptive_wall_basic.eps "Basic setup for FSI problems involving algebraic node updates for the fluid mesh. " width=0.75\textwidth
 
In an FSI problem, the fluid mesh's free
boundary is a boundary of the solid mesh, <em> i.e. </em> the
boundary along which the fluid exerts a traction onto the solid.
Within \c oomph-lib, the fluid traction is applied to the solid domain
by attaching \c FSISolidTractionElements to the faces of the 
"bulk" solid elements adjacent to the FSI boundary. 
(In the above sketch the \c FSISolidTractionElements are shown
in blue.)
The deformation of the fluid and solid meshes 
is coupled by using the \c MeshAsGeomObject formed from 
the \c FSISolidTractionElements as the \c GeomObject that defines
the moving boundary of the fluid mesh. (In sketch above, this is 
indicated by the magenta arrow.)


\subsection modification1 Modifications to allow adaptivity of the fluid and solid meshes
  

If the solid mesh is not adapted, the adaptation for the fluid
mesh  is straightforward and proceeds fully automatically
as described 
<a href="../../../interaction/fsi_collapsible_channel_algebraic/html/index.html">elsewhere</a>.
In particular, the node update data for
newly-created fluid nodes is created automatically by a call to 
the \c AlgebraicMesh::update_node_update(...) function during the adaptation. 
This function obtains the required
information about the boundary by using the \c MeshAsGeomObject 
built from the \c FSISolidTractionElements.

If the solid mesh is also adapted, then the existing 
\c FSISolidTractionElements must (at some point) be deleted and new ones
must be attached to the adapted "bulk" solid mesh. In all other problems, this
is done by
deleting the \c FSISolidTractionElements in \c
Problem::actions_before_adapt() and attaching new ones in 
\c Problem::actions_after_adapt(); see, e.g. the 
<a href="../../../poisson/two_d_poisson_flux_bc/html/index.html">tutorial
on the solution of a Poisson problem with flux boundary 
conditions.</a> However, in the present problem
this is not possible because, once the \c FSISolidTractionElements 
have been deleted, the \c MeshAsGeomObject can no longer be used to 
represent the shape and position of the FSI boundary, which would cause
the adaptation of the fluid mesh to fail. 

To avoid this problem, we adopt the following strategy:
-# When adding the various meshes to the \c Problem's
   collection of sub-meshes, we add the fluid mesh
   \b before the solid mesh. (This happens to be what was done
   already in the original driver code.) Usually, the order in which 
   sub-meshes are added to the \c Problem is irrelevant. Here the order
   \b does matter because we will exploit the fact that
   the sub-meshes are adapted individually, in the 
   order in which they were added to the \c Problem.
   \n\n
-# The \c FSISolidTractionElements are not deleted in 
   \c Problem::actions_before_adapt() and 
   remain attached to the "bulk" solid elements throughout 
   the "bulk" mesh adaptation procedure. 
   When the fluid mesh is adapted, the appropriate \c MeshAsGeomObject
   is, therefore, still fully-functional
   (and refers to the boundary as represented by the solid domain
   \b before the "bulk" solid mesh is adapted). 
   \n\n
   Here is a sketch of problem after adaptation of the fluid mesh:
   \n\n
\image html fsi_with_adaptive_wall_adapted_fluid.gif "Sketch of the problem following the adaptation of the fluid mesh. The solid mesh has not yet been refined. " 
\image latex fsi_with_adaptive_wall_adapted_fluid.eps "Sketch of the problem following the adaptation of the fluid mesh. The solid mesh has not yet been refined. " width=0.75\textwidth
   \n\n
-# The subsequent adaptation of the "bulk" solid mesh is likely 
   to turn some of the \c FSISolidTractionElements into "dangling"
   elements. (This occurs whenever a \c FSISolidTractionElements
   is attached to a "bulk" solid elements that disappears during the
   adaptation, <em> e.g. </em> by being refined.) 
   \n\n
   Here is a plot of the problem following the adaptation of the
   solid mesh :
   \n\n
\image html fsi_with_adaptive_wall_adapted_solid.gif "Sketch of the problem following the adaptation of the solid mesh -- the `dangling' FSISolidTractionElements are represented by dotted lines. " 
\image latex fsi_with_adaptive_wall_adapted_solid.eps "Sketch of the problem following the adaptation of the solid mesh -- the `dangling' FSISolidTractionElements are represented by dotted lines. " width=0.75\textwidth
   \n\n
-# Hence, in \c Problem::actions_after_adapt() we 
   delete the existing  \c FSISolidTractionElements and immediately 
   (re-)attach new ones. Now, the \c MeshAsGeomObject that represents the
   FSI boundary is broken because it still refers to the just deleted 
   \c FSISolidTractionElements. 
   \n\n
\image html fsi_with_adaptive_wall_new_face_elements.gif "Sketch of the problem following the creation of new FSISolidTractionElements. The fact that the MeshAsGeomObject is broken is indicated by the dashed lines. " 
\image latex fsi_with_adaptive_wall_new_face_elements.eps "Sketch of the problem following the creation of new FSISolidTractionElements. The fact that the MeshAsGeomObject is broken is indicated by the dashed lines. " width=0.75\textwidth
   \n\n
-# Thus, we rebuild the \c MeshAsGeomObject from the 
   newly-created \c FSISolidTractionElements, and update the
   fluid mesh's pointer to this new \c GeomObject that describes the
   boundary shape. 
   \n\n
\image html fsi_with_adaptive_wall_done.gif "Sketch of the problem with re-built MeshAsGeomObject. " 
\image latex fsi_with_adaptive_wall_done.eps "Sketch of the problem with re-built MeshAsGeomObject. " width=0.75\textwidth
   \n\n
-# Finally, we execute the \c AlgebraicMesh::update_node_update(...)
   function for all nodes in the fluid mesh to ensure that their
   node update data refers to the new \c FSISolidTractionElements.
   \n\n
-# The remaining tasks (such as the renewed setup of the fluid load on the
   \c FSISolidTractionElements via a call to 
   \c FSI_functions::setup_fluid_load_info_for_solid_elements(...), 
   etc.) remain the same as in the previous version of the code.
.

<HR> 
<HR>

\section distr Distributing the Problem
In the present example, there are two "bulk"
meshes corresponding to the fluid and solid domains
 and three "surface" meshes of traction
elements. The traction elements are \c FaceElements 
created from the "bulk" fluid elements and should be deleted before the
problem is distributed, see the tutorial on applying 
<a href="../../two_d_poisson_flux_bc_adapt/html/index.html"> flux
boundary conditions in a Poisson problem </a> for more details.
In the  <a href="../../fsi_channel_with_leaflet/html/index.html">
previous example</a> 
involving the interaction of a 2D fluid domain with a 1D beam
structure there were only two meshes: a "bulk" fluid mesh
and a "surface" solid mesh. In that problem
\b all elements in the
1D mesh of \c FSIHermiteBeamElements were retained on all processors 
as halo elements by using the function 
\c Mesh::keep_all_elements_as_halos(). The 
same methodology could be used here, but it would be extremely wasteful
to retain all the solid elements in the "bulk" solid mesh because only
the elements next to the FSI boundary are required. 
Instead, we use a more fine-grained
method of retaining elements via the function \c
GeneralisedElement::must_be_kept_as_halo(). 

<HR>
<HR>


\section impl Implementation
Most of the driver code is identical to the original serial version 
discussed in <a href="../../fsi_channel_with_leaflet/html/index.html">
another tutorial</a>. We therefore only discuss those parts of the
code that have to be changed to allow (i) the simultaneous adaptation of
the fluid and solid meshes, and (ii) the problem distribution. 

<HR>

\subsection main_body The main function

As usual in a parallel driver code, the only addition to the \c main()
function is the inclusion of calls to \c MPI_Helpers::init(), \c
MPI_Helpers::finalize(), and the \c Problem::distribute()
functions. 

<HR>

\subsection problem_class The problem class

The only additions to the serial version of the problem class are the
functions \c actions_before_distribute() and \c
actions_after_distribute(), and the helper function \c
delete_fsi_traction_elements(), discussed below.


<HR>

\subsection delete_fsi Deleting the FSISolidTractionElements

To facilitate the deletion and re-creation of the 
\c FSISolidTractionElements before and after the adaptation
(and distribution) we provide a new helper function 
\c delete_fsi_traction_elements() which complements the 
already-existing \c create_fsi_traction_elements() function:

\dontinclude turek_flag.cc
\skipline start_of_delete_traction_elements
\until end of delete traction elements


<HR>


\subsection actions_before_distribute Actions before distribute

 As discussed above, we must ensure that the "bulk" solid
elements adjacent to the FSI boundary are retained on all processors.
Hence, the \c actions_before_distribute()
function starts with a loop over the 
\c FSISolidTractionElements within which we use the 
function \c GeneralisedElement::must_be_kept_as_halo() 
to indicate that the associated bulk elements must be retained.

\dontinclude turek_flag.cc
\skipline start_of_actions_before_distribute
\until end of loop over meshes of fsi traction elements

Next, we flush all the meshes from the problem's
collection of sub-meshes and add only the "bulk" fluid and solid
meshes (in that order!). The \c FaceElements do not need to be
distributed, because they will be re-created in \c actions_after_distribute().

\until end of actions before distribute


<HR>


\subsection actions_after_distribute Actions after distribute

Following the problem distribution, we delete the old 
\c FSISolidTractionElements and then (re-)attach new ones, which will
be created as halo elements where necessary.

\dontinclude turek_flag.cc
\skipline start_of_actions_after_distribute
\until create_fsi_traction_elements()

We complete the build of the \c
FSISolidTractionElements by passing the FSI parameter and the
boundary number in the bulk mesh. The relevant code is identical to 
the serial version and we omit its listing here.

Next, we create new \c MeshAsGeomObjects from the 
newly-created \c FSISolidTractionElements and pass them to the 
(algebraic) fluid mesh:

\skipline Turn the three meshes of FSI traction elements
\until set_tip_flag_pt

The \c MeshAsGeomObjects have changed, so we
must call the \c update_node_update() function again for each
node in the fluid mesh:

\until }

Now we add the FSI traction meshes back to the problem and rebuild the
global mesh.

\until rebuild_global_mesh

Finally, we re-set the fluid load on the solid elements by
calling \c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
before re-assigning the auxiliary node update function that imposes the 
no-slip condition for all fluid nodes on the FSI boundaries.
[<a href="../../fsi_channel_with_leaflet/html/index.html#action_after_adapt">
Recall</a>
that the (re-)assignment of the auxiliary node-update function
must be performed \b after the call to  \c
FSI_functions::setup_fluid_load_info_for_solid_elements(...).]

\until end of (re-)assignment of the auxiliary node update fct

The remainder of the function identifies which processors
contain the fluid control node whose velocities 
we document in the trace file.

<HR>


\subsection actions_after_adapt Actions after adapt

The \c actions_after_adapt() function is very similar to 
\c actions_after_distribute() function, so we omit is listing here.
 The only significant differences are that (i) the redundant fluid and
 solid pressures are (re)-pinned; (ii) the identification of the 
fluid control node does not need to
 be setup;  and (iii) the traction meshes were never removed from the
 problem, so do not need to be added back in. 

<HR>

\subsection doc_solution The doc_solution() function

As with the other parallel driver codes, the main modification to the
post-processing function is the
addition of the processor number to all output files. 
Furthermore, we only write the trace file on the processors that
contain the fluid control node. In the interest of brevity 
we omit the listing of the modified function.

<HR>
<HR>

\section result Results

The figure below illustrates the distribution of the problem across four
processors, represented by the four colours, with the fluid elements
outlined in black and the solid elements outlined in white.

\image html turek_partition.gif "Distribution of the Turek \& Hron benchmark problem over four processors. " 
\image latex turek_partition.eps "Distribution of the Turek \& Hron benchmark problem over four processors. " width=0.75\textwidth

Zooming in near the "flag" shows how both fluid and solid meshes
are refined and distributed independently:

\image html turek_partition_zoom.gif "Distribution of the Turek \& Hron benchmark problem over four processors; zoomed in view near the `flag'. " 
\image latex turek_partition_zoom.eps "Distribution of the Turek \& Hron benchmark problem over four processors; zoomed in view near the `flag'. " width=0.75\textwidth

<HR>
<HR>

\section load_balance Load balancing

When employing load balancing in this problem,
we modify the time-stepping loop to perform the procedure
after each timestep:

\dontinclude turek_flag_load_balance.cc
\skipline Start of timestepping loop
\until end of timestepping loop

\subsection build_mesh The build_mesh() function

The function \c Problem::build_mesh() must be supplied by the user if
they wish to use the load balancing capability.  Thus, in this driver
code, we move all the required code to build the entire global mesh
into this function, and call it from within the problem constructor:

\dontinclude turek_flag_load_balance.cc
\skipline start_of_constructor
\until end_of_constructor

The \c build_mesh() function itself contains all the relevant code
from within the previous parallel driver code's problem constructor.

\subsection actions_functions Actions before and after load balancing 

In this example, all that is required for the \c
actions_after_load_balance() function is the
addition of the unpin-repin procedure from the \c actions_after_adapt()
function to the appropriate part of the \c actions_after_distribute()
function, since all the other functionality is already identical.  The
\c actions_before_load_balance() function is identical to the \c
actions_before_distribute() function.

<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/multi_domain/turek_flag/">
demo_drivers/mpi/multi_domain/turek_flag/
</A>
</CENTER>\n
- The main driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/multi_domain/turek_flag/turek_flag.cc">
demo_drivers/mpi/multi_domain/turek_flag/turek_flag.cc
</A>
</CENTER>
- The driver code for the load balancing example is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/multi_domain/turek_flag/turek_flag_load_balance.cc">
demo_drivers/mpi/multi_domain/turek_flag/turek_flag_load_balance.cc
</A>
</CENTER>
.
 


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/


/**

\mainpage Parallel adaptive driven cavity problem with load balancing

In this tutorial we demonstrate how to perform load balancing
to (re-)distribute elements between processors in parallel, distributed
computations. Following a brief discussion of the underlying
methodology, we illustrate the application in the
<a href="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
adaptive driven cavity problem</a> where the spatially non-uniform
refinement of the mesh leads to a significant load imbalance.

Most of the driver code is identical to the codes discussed in the
tutorials explaining the 
<a href="../../../navier_stokes/adaptive_driven_cavity/html/index.html">
serial</a> and 
<a href="../../adaptive_driven_cavity/html/index.html">distributed, parallel
</a> solution of the problem. Therefore we only discuss the changes 
required to perform load balancing on the problem.

<HR>

\section load_balance Load balancing

The initial distribution of a problem via a call to 
\code
Problem::distribute(...)
\endcode
attempts to distribute the elements in the Problem's mesh over the available
processors such that (i) each processor stores approximately the same
number of elements and (ii) the anticipated volume of inter-processor 
communication required to synchronise the solution across processor 
boundaries is minimised. Typically, this procedure works very well in
the sense that the assembly times for the Jacobian matrix in a 
distributed computation scale extremely well with the number of processors.

A re-distribution of elements may be required because
-# strongly non-uniform mesh adaptation may lead to a drastic
   increase in the number of elements on some processors;
   \n\n
-# in multi-physics computations the cpu times required to compute
   the elements' contributions to the global Jacobian matrix
   may differ significantly between different element types. In such 
   cases, the mere equidistribution of elements between processors
   will not achieve good parallel scaling;
   \n\n
-# when restarting a computation on a larger number of processors,
   elements need to be re-distributed to populate the (otherwise
   empty) additional processors.
.

As with other methods within \c oomph-lib, load balancing is
implemented such that only minimal user intervention is required.
The function
\code
Problem::load_balance() 
\endcode
may be called at any point following the distribution of the problem. The only
change required to an existing driver code is the provision (via overloading
of a broken virtual function in the Problem base class) of the
function
\code
Problem::build_mesh()
\endcode
This function must
-# Build the mesh (and in the case of multiple sub-meshes, build these
   and combine them to a global mesh using the 
   \c Problem::build_global_mesh() function).
  \n\n
-# Complete the build of the elements by setting pointers to 
   physical parameters (such as Reynolds numbers) or functions
   (such as source function pointers), etc.
   \n\n
-# Apply the boundary conditions.
   \n\n
-# Perform any uniform mesh refinement that was applied before
   calling \c Problem::distribute().
.
Typically this requires no more than a straightforward cut-and-paste 
of code from the problem constructor into the \c Problem::build_mesh() 
function; see also the discussion \ref What for more details. 

The load balancing routines then perform the following steps:
-# Take the current (distributed) global mesh, and calculate a new
   partition for each of the current elements, taking into account
   the cpu time that each element spent on the most recent computation of its 
   contribution to the problem's Jacobian matrix.
   \n\n
-# Build a new (global) mesh using the \c Problem::build_mesh()
   function.
   \n\n
-# Distribute the new (global) mesh according to the new
   partitioning.
   \n\n
-# Once distributed, refine the new (global) mesh on each processor
   to achieve the same refinement pattern as in the original problem.
   \n\n
-# Copy the \c Data values from the old to the new (global) meshes.
   \n\n
. 
We note that for "structured" meshes (i.e. meshes whose refinement 
pattern is represented by 
<a href="../../../the_data_structure/html/index.html#RefineableQuadMesh_setup_section">
"tree forests"</a>) only complete trees can be
moved between processors. If the mesh was refined uniformly after
being distributed, a more fine-grained tree-forest (which may allow
better load balancing) can be built
by calling \c Problem::prune_halo_elements_and_nodes() 
as discussed in <a href="../../../mpi/general_mpi/html/index.html#how_it_works">
another tutorial.</a>
<HR>

\section example An example: Revisiting the adaptive driven cavity problem 

In this section we outline the required changes to the
<a href="../../adaptive_driven_cavity/html/index.html">parallel version
of the adaptive driven cavity problem</a> so that the problem can use
the load balancing method described above.

The figure below demonstrates the advantages of load balancing
in that problem: The left hand panel shows the distribution of the mesh across
four processors (indicated by the colours) after three 
spatially adaptive solves. Note how the singularities in the 
bottom corners result in 
a strongly non-uniform spatial refinement which leads to a significant
load imbalance because the "pink" and "grey" processors contain far
more elements than then "green" and "cyan" ones. The right hand panel 
shows the distribution of the mesh across four processors when 
load balancing is performed in between the second and third mesh adaptation.

\image html load_balance_partition.gif "Plot illustrating the distribution of the mesh for the adaptive driven cavity problem, both without (left image) and with (right image) the use of load balancing. " 
\image latex load_balance_partition.eps "Plot illustrating the distribution of the mesh for the adaptive driven cavity problem, both without (left image) and with (right image) the use of load balancing. " width=0.75\textwidth

\subsection build_mesh The build_mesh() function

As discussed above, the function \c Problem::build_mesh()
is created most easily by moving the code that (i) creates the mesh, (ii)
applies the relevant boundary conditions, and (iii) completes the
build of all the elements in the problem from
the problem constructor.

\dontinclude adaptive_driven_cavity_load_balance.cc
\skipline start_of_build_mesh
\until end_of_build_mesh

\subsection problem_constructor Changes to the problem constructor

The problem constructor becomes very short since the bulk of the
code has been moved into the \c Problem::build_mesh() function.

\dontinclude adaptive_driven_cavity_load_balance.cc
\skipline start_of_constructor
\until end_of_constructor

<HR>

\section customising_load_balance Customising the load balancing

The <A HREF="../../../../demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity_load_balance.cc">driver code</a> demonstrates some of the different options available 
for the \c Problem::load_balance(...) function. 
- Passing a boolean argument to \c
  Problem::load_balance(...) enables the output of extended
  statistics.
  \n\n
- By default, the mesh partitioning for the initial problem
  distribution and any subsequent load balancing operations is
  determined by METIS.
  Unfortunately, METIS is not completely deterministic. This makes
  it difficult to compute reference data that can be used to assert
  the correct execution of the code during self-tests. When the driver code
  is run with the \c --validate command line flag, we therefore 
  call the \c Problem::distribute(...) function with a pre-determined
  (and deterministic but non-optimal) distribution of the elements. 
  The use of METIS
  during the load balancing operations can be bypassed by calling
  \code
  Problem::set_default_partition_in_load_balance()
  \endcode


<HR>
<HR>

\section comm_and_ex Comments and Exercises

\subsection comm Comments
\subsubsection What What goes into the build_mesh() function?
The main rule regarding what should (and should not) be moved from the
problem constructor into the \c build_mesh() function is that
following the return from the \c build_mesh() function, all meshes
should be re-generated (and refined to the same degree as when
\c Problem::distribute() was called), their constituent elements made fully
functional, and all boundary conditions applied.

  It is \b not necessary (and would, in fact, be undesirable) to 
re-create \c Timesteppers and/or \c GeomObjects that are used to 
define curvilinear mesh boundaries. We recommend 
generating such objects once (in the Problem constructor) and making 
them available to the \c build_mesh() function by storing pointers to
them in the problem's private member data.

It is not necessary to re-generate the error estimator, though
the pointer to it (and any non-default target errors) need to be
passed to all (newly re-generated) adaptive meshes.

\subsubsection pointers Updating pointers
Since load balancing re-generates the meshes -- and thus 
their constituent elements and nodes -- pointers to such objects
must be re-assigned on return from \c Problem::load_balance().
In our experience such pointers tend to be used predominantly in
- Post-processing functions (e.g. to document the solution at a fixed
  node)
  \n\n
- Multi-physics/block preconditioners which tend to use 
  pointers to meshes to classify the degrees of freedom. 
  \n\n
.
Failure to re-assign any dangling pointers will cause segmentation faults --
recompile \c oomph-lib with debugging enabled and use 
<a href="http://www.gnu.org/software/ddd">ddd</a> to see
where the code crashes. 

      
\subsection ex Exercises

-# Explore what happens if you "forget" to implement the
   \c Problem::build_mesh() function (e.g. by renaming it
   \c my_build_mesh(), say, so that it no longer overloads 
   the function the \c Problem base class.
   \n\n
-# The \c Problem::build_mesh() shown above
   (accidentally) illustrates a common problem with a mere
   cut-and-paste approach -- it creates a memory leak!
   Where is it and how would you fix it?
   \n\n
-# The \c demo_drivers directory contains a few additional
   driver codes that employ load balancing. These codes 
   exist mainly for self-test purposes and do not have separate 
   tutorials. It may be instructive to compare the
   different versions of these codes to further clarify
   the modifications required to enable load balancing.
   We suggest you compare \n\n
   - The original version of the distributed, adaptive driven cavity code \n\n
     <center>
     <A HREF="../../../../demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity.cc">
     demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity.cc
     </A>
     </center>
     \n\n
     and the version with load balancing
     \n\n
     <center>
     <A HREF="../../../../demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity_load_balance.cc">
     demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity_load_balance.cc
     </A>
     </center>
     \n\n
   - The original version of the distributed code for the
     doubly-adaptive solution of the unsteady
     heat equation, \n\n
     <center>
     <A HREF="../../../../demo_drivers/unsteady_heat/two_d_unsteady_heat_2adapt/two_d_unsteady_heat_2adapt.cc">
     demo_drivers/unsteady_heat/two_d_unsteady_heat_2adapt/two_d_unsteady_heat_2adapt.cc
     </A>
     </center>
     \n\n
     and the version with load balancing 
     \n\n
     <center>
     <A HREF="../../../../demo_drivers/mpi/distribution/restart/two_d_unsteady_heat_2adapt_load_balance.cc">
     demo_drivers/mpi/distribution/restart/two_d_unsteady_heat_2adapt_load_balance.cc
     </A>
     </center>
     \n\n
   - The original version of the distributed code for the
     solution of Turek & Hron's FSI benchmark problem
     heat equation, \n\n
     <center>
     <A HREF="../../../../demo_drivers/mpi/multi_domain/turek_flag/turek_flag.cc">
     demo_drivers/mpi/multi_domain/turek_flag/turek_flag.cc
     </A>
     </center>
     \n\n
     and the version with load balancing 
     \n\n
     <center>
     <A HREF="../../../../demo_drivers/mpi/multi_domain/turek_flag/turek_flag_load_balance.cc">
     demo_drivers/mpi/multi_domain/turek_flag/turek_flag_load_balance.cc
     </A>
     </center>
     \n\n
   .
   Comparing the codes with sdiff is particularly instructive.
.
<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/distribution/adaptive_driven_cavity">
demo_drivers/mpi/distribution/adaptive_driven_cavity
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity_load_balance.cc">
demo_drivers/mpi/distribution/adaptive_driven_cavity/adaptive_driven_cavity_load_balance.cc
</A>
</CENTER>




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

